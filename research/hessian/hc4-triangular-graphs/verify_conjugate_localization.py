"""Exact controls for LD0/1;not a whole-proof verifier or a JC decision test."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def check(name,condition):
 assert condition,name
 checks.append(name)
def zero(name,expr):check(name,s.cancel(expr)==0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
t,z,u=s.symbols('t z u');pol=t**3-3*t-z;h=3*t**2-3;Delta=4-z**2
zero('cubic derivative is source unit after prescribed localization',s.diff(pol,t)-h)
zero('cubic discriminant exact reduced branch polynomial',s.discriminant(pol,t)-27*Delta)
zero('norm of source Jacobian',s.resultant(pol,h,t)+27*Delta)
zero('branch minus2 has double critical and simple retained root',pol.subs(z,-2)-(t-1)**2*(t+2))
zero('branch plus2 has double critical and simple retained root',pol.subs(z,2)-(t+1)**2*(t-2))
check('simple retained root minus2 is etale',h.subs(t,-2)==9)
check('simple retained root plus2 is etale',h.subs(t,2)==9)
zero('all-conjugate localization additionally removes both retained roots',Delta.subs(z,t**3-3*t)+(t-2)*(t+2)*(t-1)**2*(t+1)**2)
M=s.Matrix([[0,0,z],[1,0,3],[0,1,0]]);H=3*M**2-3*s.eye(3)
check('companion cubic polynomial',M**3-3*M-z*s.eye(3)==s.zeros(3))
zero('companion Jacobian norm equals root resultant',H.det()+27*Delta)
check('companion source unit inverse exact',s.simplify(H*H.inv())==s.eye(3))
zero('sum of reciprocal source Jacobians cancels',s.trace(H.inv()))
zero('square reciprocal source Jacobians has nonzero branch poles',s.trace(H.inv()**2)-s.Rational(2,3)/Delta)
zero('invariant reciprocal is product of all conjugate units',-27/H.det()-1/Delta)
for n in [0,1,2,5]:
 zero(f'all-conjugate invariant trace includes pole depth{n}',s.trace(s.eye(6)/Delta**n)-6/Delta**n)

# Finite tame orders independently check trace-dual finite stages.
# The infinite-stage assertion is proved using finite-module valuations in the note.
for e in [1,2,3,5]:
 C=s.zeros(e)
 for i in range(e-1):C[i+1,i]=1
 C[0,e-1]=u
 Gram=s.Matrix(e,e,lambda i,j:s.trace(C**(i+j)))
 check(f'tame trace Gram nondegenerate e{e}',Gram.det()!=0)
 Inv=Gram.inv()
 check(f'tame trace inverse basis e{e}',s.simplify(Gram*Inv)==s.eye(e))
 for n in [0,2]:
  stage=u**n*Inv
  check(f'finite-stage trace dual scaling e{e} n{n}',s.simplify(Gram*stage/u**n)==s.eye(e))
 for a in [-3,0,4]:
  j=(-a)%e;n=(a+j)//e+2
  value=s.trace(C**(a+j))/u**n
  zero(f'fixed multiplier monomial develops deeper trace pole e{e} a{a}',value-e/u**2)
  check(f'finite multiplier cannot clear this deeper stage e{e} a{a}',s.denom(s.cancel(value))==u**2)
ALG=len(checks)
context=json.loads((HERE/'sources/conjugate-localization-context.json').read_text())
for p,digest in context['retained_proof_and_source_sha256'].items():check('retained bytes '+p,sha(ROOT/p)==digest)
for p,digest in context['new_primary_sha256'].items():check('new primary bytes '+p,sha(HERE/p)==digest)
paper=(HERE/'sources/collision-geometry-paper2-eae2206.tex').read_text()
lean=(HERE/'sources/collision-geometry-secant-frame-denominator-eae2206.lean').read_text()
check('primary manuscript title and source limitation retained','Planar Secant Projectors' in paper and 'not necessarily' in paper and 'E^G=B' in paper)
check('Lean finite coefficient theorem and uniform limitation retained','finiteCoefficientDenominatorIdeal_ne_bot' in lean and 'uniform landing assertion' in lean)
check('Lean B-linear auxiliary theorem not called Keller frame','is not `B`-linear' in lean and 'linearOperatorDenominatorIdeal_ne_bot' in lean)
proof=HERE/'ALL_CONJUGATE_LOCALIZATION_AND_UNIFORM_TRACE_LANDING.md';flat=' '.join(proof.read_text().split())
for phrase in ['CURRENT LAST UNSOLVED IMPLICATION','E=S_M[1/f]','E^G=B[1/f]','OPEN AND CLOSED','FREE','proper locus','ONLY AFTER localizing','T_(R\') != 0','UNIFORM','ENTIRE source','intersection_(n>=0)','UNRAMIFIED source point','NOT the whole affine plane','Conjugate-sheet fixed loci','No Lean build','DEFERRED']:
 check('written scope '+phrase,phrase in flat)
check('no new external theorem imports',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=ALG,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-conjugate-localization.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),algebraic_controls=ALG,elapsed_seconds=result['elapsed_seconds']),indent=2))
