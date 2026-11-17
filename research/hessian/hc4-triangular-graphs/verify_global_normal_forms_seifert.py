"""Exact controls for the weighted global normal form and Seifert consequences."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
p,c,l,a,X,P,Q,z,rho,s=sp.symbols('p c l a X P Q z rho s')
M=sp.Matrix([[1-p*p,2*p-p**3],[-p,1-p*p]])
zero(M.det()-1,'the weighted root-twist matrix is globally unimodular in every degree')
zeros(list(M*M.adjugate()-sp.eye(2)),'the weighted twist matrix has its exact polynomial inverse')
old=M*sp.Matrix([X,P]);e=1-2*p;q=sp.expand(old[0]*old[1]-e*X*P)
fac=sp.div(q,p*(p-1),p)[0]
zero(q-p*(p-1)*fac,'the weighted product discrepancy has its universal p(p-1) factor')
zero(sp.rem(q,p*(p-1),p),'CRT root values0/1 suffice in arbitrary multiplier degree')
forward=[old[0],old[1],Q-a*l*fac]
back=M.adjugate()*sp.Matrix([X,P])
inverse=[back[0],back[1],Q+a*l*fac.subs({X:back[0],P:back[1]},simultaneous=True)]
zeros([v.subs(dict(zip([X,P,Q],forward)),simultaneous=True)-r for v,r in zip(inverse,[X,P,Q])],
      'the full weighted twist inverse cancels its forward fibre map')
zeros([v.subs(dict(zip([X,P,Q],inverse)),simultaneous=True)-r for v,r in zip(forward,[X,P,Q])],
      'the forward weighted fibre map cancels its complete inverse')
zero(a*forward[0]*forward[1]+c*forward[2]-(a*e*X*P+c*Q)
     -a*fac*(p*(p-1)-c*l),
     'the actual potential identity holds universally modulo p(p-1)=C_N*l')
zero(sp.Matrix(forward).jacobian([X,P,Q]).det()-1,'the full fibre twist has determinant one')
for choice in [0,1]:
    zero((a*e).subs(p,choice)-a*(1-2*choice),'a chosen root sign changes the weighted coefficient orientation '+str(choice))
    zero(q.subs(p,choice),'the product correction vanishes on each CRT root choice '+str(choice))
# Entire four-coordinate determinant:arbitrary coefficient dependence on z.
pf=sp.Function('p')(z);af=sp.Function('a')(z);lf=sp.Function('l')(z)
fwdz=[v.subs({p:pf,a:af,l:lf},simultaneous=True) for v in forward]
zero(sp.Matrix([fwdz[0],z,fwdz[1],fwdz[2]]).jacobian([X,z,P,Q]).det()-1,
     'derivatives of arbitrary interpolating coefficient polynomials do not change the source Jacobian')
# General cotangent block formula and removal of base-only differences in the coefficient ideal.
f1,f2,d,dd,alpha,beta=sp.symbols('f1 f2 d dd alpha beta')
zero((dd+alpha*f1+beta*f2)+(P-alpha)*f1+(Q-beta)*f2-(dd+P*f1+Q*f2),
     'an arbitrary base difference in the coefficient ideal is removed by an actual polynomial fibre translation')
j11,j12,j21,j22,t11,t12,t22=sp.symbols('j11 j12 j21 j22 t11 t12 t22')
HH=sp.Matrix([[t11,t12,j11,j21],[t12,t22,j12,j22],[j11,j12,0,0],[j21,j22,0,0]])
zero(HH.det()-(j11*j22-j12*j21)**2,'the whole cotangent Hessian is the square of its ordered coefficient Jacobian')
# The explicit three-point source-alignment target is checked as one required family,not a degree sweep.
C3=(z-1)*(z-2)*(z-3);a3=(3*z*z-12*z+10)/2
D3=((rho+3)*z*z-(3*rho+13)*z+2*rho+12)/2
R=a3*X*P+C3*Q+D3
bez=sp.cancel((1-a3*sp.diff(C3,z))/C3)
check(sp.Poly(bez,z).is_univariate,'the normalizing coefficient inverse has an actual polynomial Bezout complement')
zero(a3*sp.diff(C3,z)+C3*bez-1,'a3 C3prime equals one modulo the full multiplier')
zero(sp.diff(C3,z)*(a3*X)+bez*X*C3-X,'the whole explicit coefficient ideal is (X,C3),not just its zero set')
zero(sp.Matrix([a3*X,C3]).jacobian([X,z]).det()-a3*sp.diff(C3,z),'the explicit target coefficient Jacobian')
zero(sp.hessian(R,[X,z,P,Q]).det()-(a3*sp.diff(C3,z))**2,
     'the native target Hessian is not a unit;unit existence remains the separate E3 gate')
check(sp.Poly(a3*sp.diff(C3,z),z).degree()>0,'the explicit target is not incorrectly declared native unit Hessian')
for i,value in [(1,1),(2,-1),(3,rho)]:
    zero(D3.subs(z,i)-value,'the target interpolates the exact critical value at base point '+str(i))
    zero((a3*sp.diff(C3,z)).subs(z,i)-1,'all target point orientations are literally one '+str(i))
    point={X:0,z:i,P:0,Q:-sp.diff(D3,z).subs(z,i)/sp.diff(C3,z).subs(z,i)}
    zeros([sp.diff(R,v).subs(point,simultaneous=True) for v in [X,z,P,Q]],'the complete target critical equations at point '+str(i))
    zero(sp.hessian(R,[X,z,P,Q]).det().subs(point,simultaneous=True)-1,
         'the target full critical determinant is one for every rho at point '+str(i))
    zero(R.subs(point,simultaneous=True)-value,'the entire target critical point has its stated value '+str(i))
for collision in [-1,1]:
    vals=[sp.cancel(D3.subs({z:i,rho:collision},simultaneous=True)) for i in [1,2,3]]
    check(len(set(vals))==2,'critical-value collision retains three distinct simple base points '+str(collision))
# An actual fibre matrix,not mere point data,is required by the orientation repair.
b11,b12,b21,b22,rr,ss=sp.symbols('b11 b12 b21 b22 rr ss')
B=sp.Matrix([[b11,b12],[b21,b22]]);ff=sp.Matrix([f1,f2]);kk=B*ff
repair=B+sp.Matrix([rr,ss])*sp.Matrix([[-f2,f1]])
zeros(list(repair*ff-kk),'the universal GL2 determinant repair preserves the whole ordered pair')
zero(repair.det()-B.det()-ss*kk[0]+rr*kk[1],'the global repair changes the determinant by the required ideal combination')
# Degree-three transpose is ordinary;identity monodromy and block triangularity force epsilon I.
check((-1)**((3+1)*(3+1))==1,'the graded transpose has positive sign in the actual middle degree')
eps,A,Bb,Cc=sp.symbols('eps A B C')
L=sp.Matrix([[eps,0,0],[A,eps,0],[Bb,Cc,eps]])
solutions=sp.solve(list(L-L.T),[A,Bb,Cc],dict=True)
check(solutions==[{A:0,Bb:0,Cc:0}],
      'triangular Seifert form with identity full monodromy has no cross blocks')
zero((L-L.T).subs({A:0,Bb:0,Cc:0}).norm(),'the corresponding ordinary intersection form vanishes')
cycle=sp.Matrix([[1,0,1],[0,1,1],[0,0,1]])
puncture=sp.Matrix([[2,1,-1],[1,2,-1],[-1,-1,1]])
check(cycle.det()==1 and puncture.det()==1,'both the collision-uniform cycle matrix and puncture Seifert matrix are unimodular')
zeros(list(cycle*puncture*cycle.T-sp.eye(3)),'the S_rho puncture Seifert matrix is integrally congruent to identity')
zeros(list(cycle.inv()*sp.eye(3)*cycle.inv().T-puncture),'the inverse basis formula gives the exact native puncture matrix')
algebra=len(checks);context=json.loads((HERE/'sources/global-normal-forms-seifert-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected full nonproper Seifert primary '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==3
      and all(e['hypotheses_checked'] for e in context['newly_imported_proved_external_statements']),
      'three explicit proved imports use only isolated finite singularities and permit unrestricted infinity behaviour')
check(context['newly_used_unproved_external_statements']==[],'no global alignment,Keller fibre or tameness conjecture is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['K2_h1_source_equivalence_constructed']
      and not flags['S_rho_normal_form_source_alignment_constructed']
      and not flags['all_four_variable_polynomials_classified']
      and not flags['infinity_Milnor_pieces_proved_absent_or_contractible']
      and not flags['original_goal_complete'],
      'explicit normal-form gates and homological Seifert closure do not claim their unresolved geometric conditions')
proof=HERE/'GLOBAL_ORIENTED_NORMAL_FORMS_AND_SEIFERT_BARRIER.md'
paths=['sources/global-normal-forms-seifert-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-global-normal-forms-seifert.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
