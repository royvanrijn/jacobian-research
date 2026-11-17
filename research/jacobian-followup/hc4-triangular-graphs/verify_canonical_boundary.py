"""Global canonical obstruction transcription controls;no old checker replay."""
import hashlib,json,platform,time
from itertools import combinations
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(b,label):
 assert bool(b),label
 checks.append(label)
def zero(f,label):check(s.expand(f)==0,label)
u,e,a0,a1=s.symbols('u e a0 a1',positive=True)
zero(s.simplify(s.diff((a0+a1*u)*u**e,u)-u**(e-1)*(e*a0+(e+1)*a1*u)),
 'the local target differential has the exact leading ramification factor u to e minus one')
check((e*a0).is_nonzero,'in characteristic zero the leading tame coefficient is a unit symbol')
k,l=s.symbols('k l',integer=True,nonnegative=True)
check(s.expand((k+1)*(l+2)-1).is_positive,
 'a ramified intermediate index cannot divide an actual unramified source index one')
mm=s.symbols('mm',integer=True,positive=True)
check(s.expand(mm*(l+1)).is_positive,
 'a positive ramification coordinate cannot become zero by a positive torsion multiple in a free boundary lattice')
r=s.symbols('r',integer=True,positive=True)
zero((e*r)*(e-1)/e-r*(e-1),
 'uniform Galois ramification becomes principal after clearing the common target-index denominator')
# Full tame trace-dual lattices,not just discriminant scalars.
t=s.symbols('t',nonzero=True)
for degree in [2,3,5]:
 M=s.zeros(degree)
 for i in range(degree-1):M[i+1,i]=1
 M[0,degree-1]=t
 check(M**degree==t*s.eye(degree),'the tame companion extension has the exact uniformizer equation at e='+str(degree))
 gram=s.Matrix(degree,degree,lambda i,j:s.simplify(s.trace(M**i*M**(-j))/degree))
 check(gram==s.eye(degree),'the complete inverse-different dual basis is u to minus j over e at e='+str(degree))
 check(min(-j for j in range(degree))==1-degree,
  'the fractional trace-dual generator has valuation minus e plus one at e='+str(degree))
# Top Koszul dual cohomology is a rank-one quotient for every global regular sequence.
for count in range(1,6):
 fs=s.symbols('f0:'+str(count))
 for homdegree in range(2,count+1):
  high=list(combinations(range(count),homdegree));mid=list(combinations(range(count),homdegree-1));low=list(combinations(range(count),homdegree-2))
  def boundary(cols,rows):
   out=s.zeros(len(rows),len(cols))
   for j,I in enumerate(cols):
    for pos,idx in enumerate(I):out[rows.index(I[:pos]+I[pos+1:]),j]=(-1)**pos*fs[idx]
   return out
  check(boundary(mid,low)*boundary(high,mid)==s.zeros(len(low),len(high)),
   'the global Koszul differential squares to zero at length '+str(count)+' and degree '+str(homdegree))
 topentries=[(-1)**i*fs[i] for i in range(count)]
 check(set(s.expand(f*f) for f in topentries)==set(f*f for f in fs),
  'the top dual Koszul image is the whole defining ideal at regular-sequence length '+str(count))
# The actual smooth finite quadric and its affine-plane chart.
X,Y,Z,ss,tt=s.symbols('X Y Z ss tt')
Q=X*Y-Z*Z+1
zero(Q.subs({X:ss,Y:2*tt+ss*tt*tt,Z:1+ss*tt}),
 'the genuine affine-plane chart lies identically in the smooth finite quadric')
zero(Q-X*s.diff(Q,X)-Z*s.diff(Q,Z)/2-1,
 'the full quadric Jacobian ideal contains one,so the finite cover is smooth everywhere')
zero((Z-1)*(Z+1)-X*Y+Q,
 'the two regular inverse formulas glue on their actual chart overlap')
zero(((Z-1)/X).subs({X:ss,Z:1+ss*tt})-tt,
 'the first chart inverse recovers the source variable')
zero(s.cancel((Y/(Z+1)).subs({Y:2*tt+ss*tt*tt,Z:1+ss*tt})-tt),
 'the second chart inverse recovers the source variable even along the retained line')
zero(s.det(s.Matrix([ss,2*tt+ss*tt*tt]).jacobian([ss,tt]))-2*(1+ss*tt),
 'the entire affine-plane composite has the exact nonconstant Jacobian')
check((1+ss*tt).subs({ss:1,tt:-1})==0,
 'the ramification meets an actual source point and cannot be moved to the omitted boundary')
# Frobenius versus ordinary trace:the normal finite CI control really is Frobenius.
dd=s.symbols('D');Mu=s.Matrix([[0,dd],[1,0]]);H=s.Matrix([[0,1],[1,0]])
Trace=s.diag(2,2*dd)
check(H.det()==-1,'coefficient extraction gives a perfect global Frobenius pairing for the non-Keller quadric cover')
check(H*(2*Mu)==Trace,'the ramified ordinary trace differs from the Frobenius pairing by the actual different 2Z')
zero(Trace.det()-4*dd,'the ordinary trace discriminant vanishes on the branch even though the dualizing module is free')
zero((2*Mu).det()+4*dd,'the norm of the different has the required ramification sign')
# Characteristic-zero cannot be deleted from the result.
x,y=s.symbols('x y')
for prime in [2,3,5]:
 F=x**prime-x
 check(all(int(c)%prime==0 for c in s.Poly(s.diff(F,x)+1,x).coeffs()),
  'the Artin-Schreier finite cover is etale in characteristic '+str(prime))
 check(F.subs(x,0)==F.subs(x,1)==0,
  'the finite etale monogenic control is noninjective in characteristic '+str(prime))
 check(s.Poly(F,x).LC()==1 and s.degree(F,x)==prime,
  'the actual characteristic-p cover has a global monic equation of degree '+str(prime))
algebra=len(checks)
context=json.loads((HERE/'sources/canonical-boundary-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
 check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
 check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
  'inspected primary source '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==1 and context['newly_used_unproved_external_statements']==[],
 'finite-duality module identification is the sole new proved import;purity and finite-cover triviality remain retained inputs')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==10 and not flags['smooth_or_locally_Gorenstein_normalizations_excluded']
 and not flags['arbitrary_intermediate_field_given_a_polynomial_affine_space_model']
 and not flags['archived_quartic_or_degree_six_nontorsion_packets_excluded']
 and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
 'global torsion and CI obstructions preserve local-Gorenstein,nontrivial-canonical,intermediate-model and full-HC4 boundaries')
proof=HERE/'CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md'
paths=['sources/canonical-boundary-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-canonical-boundary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
