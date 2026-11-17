"""Exact universal identity controls;the canonical note supplies all-degree proofs."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
X,P,Q,a,c,A,B,d,at,ct,At,Bt,dt,au,cu,Au,Bu,du=sp.symbols('X P Q a c A B d at ct At Bt dt au cu Au Bu du')
relations={B:(1-A*a)/c,Bt:(1-At*at)/ct,Bu:(1-Au*au)/cu}
def bz(value):return sp.cancel(value.subs(relations,simultaneous=True))
def bz0(value,label):zero(bz(value),label)
def bzs(values,label):zeros([bz(v) for v in values],label)
h=a*X*P+c*Q+d;H=h-d;K=-X*B*P+A*Q
bz0(A*H-c*K-X*P,'the first chart inverse is exact under the universal Bezout relation')
bz0(B*H+a*K-Q,'the second chart inverse is exact under the universal Bezout relation')
bz0(sp.Matrix([H,K]).jacobian([P,Q]).det()-X,'the entire chart determinant is exactly X')
Pt=(At*(h-dt)-ct*K)/X;Qt=Bt*(h-dt)+at*K
ht=at*X*Pt+ct*Qt+dt;Kt=-X*Bt*Pt+At*Qt
bz0(ht-h,'universal finite transport preserves the potential')
bz0(Kt-K,'universal finite transport preserves the lifted invariant')
bzs([(A*(ht-d)-c*Kt)/X-P,B*(ht-d)+a*Kt-Q],
    'the universal reverse transport is a two-sided inverse')
bzs([(Au*(ht-du)-cu*Kt)/X-(Au*(h-du)-cu*K)/X,
     Bu*(ht-du)+au*Kt-(Bu*(h-du)+au*K)],
    'three arbitrary endpoint frames satisfy the universal composition law')
# Independent endpoint increments model all polynomial divided differences at once.
da,dc,dA,dB,dd=sp.symbols('da dc dA dB dd')
increment={at:a+X*da,ct:c+X*dc,At:A+X*dA,Bt:B+X*dB,dt:d+X*dd}
Ppoly=P+dA*H-(A+X*dA)*dd-dc*K
Qpoly=Q+X*(dB*H-(B+X*dB)*dd+da*K)
bz0(Pt.subs(increment,simultaneous=True)-Ppoly,'the universal apparent X denominator cancels in P')
bz0(Qt.subs(increment,simultaneous=True)-Qpoly,'the universal Q increment is divisible by X')
check(sp.Poly(Ppoly,X,P,Q,da,dc,dA,dB,dd).is_multivariate,
      'the lifted expression is polynomial without localizing any coefficient')
bz0(a*K+B*(h-d)-Q,'the exact Derksen recovery identity supplies Q')
# First-order endpoint increments satisfy the differentiated Bezout relation.
LP=dA*H-A*dd-dc*K;LQ=X*(dB*H-B*dd+da*K)
derivative_B=-(dA*a+A*da+B*dc)/c
def dzbz(value):return bz(value.subs(dB,derivative_B))
zero(dzbz(da*X*X*P+a*X*LP+dc*X*Q+c*LQ+dd*X),
     'the lifted infinitesimal vector field preserves h universally')
zero(dzbz(-X*(dB*X*P+B*LP)+dA*X*Q+A*LQ),
     'the lifted infinitesimal vector field preserves K universally')
EX=sp.Matrix([0,c,0,-a*X]);EP=sp.Matrix([c,0,0,-a*P]);LX=sp.Matrix([0,LP,X,LQ])
# Rows X,P,z give a nonzero minor exactly on C!=0,X!=0.
minor=sp.Matrix.hstack(EX,EP,LX).extract([0,1,2],[0,1,2]).det()
zero(minor+c*c*X,'the three global fields span the fibre tangent on the stated open')
bz0(c*sp.diff(K,P)-a*X*sp.diff(K,Q)+X,
     'E_X acts as minus X partial_K in the exact chart')
# At a smooth residual point the lift escapes the doubly-zero hyperbolic coordinates.
root={X:0,P:0,c:0,A:1/a}
zero(LP.subs(root,simultaneous=True)+(dc*Q+dd)/a,
     'at C=0 and X=P=0 the P velocity is minus kappa divided by a')
zero(LQ.subs(root,simultaneous=True),'the escaping root flow keeps Q fixed')
zero((c*Q+a*X*P).subs(root,simultaneous=True),'the residual point has H=0')
# One explicit family is retained because it is the actual alignment target,not a degree sample.
z,t,rho=sp.symbols('z t rho');C3=(z-1)*(z-2)*(z-3);a3=(3*z*z-12*z+10)/2
A3=sp.diff(C3,z);B3=sp.cancel((1-a3*A3)/C3)
D3=((rho+3)*z*z-(3*rho+13)*z+2*rho+12)/2
check(sp.Poly(B3,z).is_univariate,'the actual three-point family has a polynomial Bezout complement')
zero(A3*a3+B3*C3-1,'the actual three-point Bezout certificate is exact')
R3=a3*X*P+C3*Q+D3;K3=-X*B3*P+A3*Q
def shift(f):return sp.expand(f.subs(z,z+t*X))
def delta(f):return sp.cancel((shift(f)-f)/X)
for f,label in [(a3,'a3'),(C3,'C3'),(A3,'A3'),(B3,'B3'),(D3,'D_rho')]:
    check(sp.Poly(delta(f),X,z,t,rho).is_multivariate,'actual divided difference polynomial '+label)
    zero(shift(f)-f-X*delta(f),'actual divided difference identity '+label)
p3=P+delta(A3)*(R3-D3)-shift(A3)*delta(D3)-delta(C3)*K3
q3=Q+X*(delta(B3)*(R3-D3)-shift(B3)*delta(D3)+delta(a3)*K3)
zero(shift(a3)*X*p3+shift(C3)*q3+shift(D3)-R3,'the full actual R_rho flow preserves its potential for every rho')
zero(-X*shift(B3)*p3+shift(A3)*q3-K3,'the full actual R_rho flow preserves K for every rho')
zero(p3.subs(t,0)-P,'the actual P flow has identity time zero')
zero(q3.subs(t,0)-Q,'the actual Q flow has identity time zero')
zero(sp.diff(p3,t).subs(t,0)-(sp.diff(A3,z)*(R3-D3)-A3*sp.diff(D3,z)-sp.diff(C3,z)*K3),
     'the explicit flow differentiates to the displayed P derivation')
zero(sp.diff(q3,t).subs(t,0)-X*(sp.diff(B3,z)*(R3-D3)-B3*sp.diff(D3,z)+sp.diff(a3,z)*K3),
     'the explicit flow differentiates to the displayed Q derivation')
for r in [1,2,3]:
    point={X:0,P:0,z:r,Q:-sp.diff(D3,z).subs(z,r)/sp.diff(C3,z).subs(z,r)}
    zero((sp.diff(A3,z)*(R3-D3)-A3*sp.diff(D3,z)-sp.diff(C3,z)*K3).subs(point,simultaneous=True),
         'the actual critical point is fixed by the lift '+str(r))
    zero(sp.diff(R3,z).subs(point,simultaneous=True),'the actual full critical residual equation '+str(r))
# Universal ordinary-suspension tangent calculation,including the whole-value K2 model.
x,y,b1,b2=sp.symbols('x y b1 b2')
Du=sp.Matrix([0,-b1,2*x,0]);Dv=sp.Matrix([0,-b2,0,2*x]);Mu=sp.Matrix([-b1,0,2*y,0])
grad=sp.Matrix([2*y,2*x,b1,b2])
zeros([grad.dot(v) for v in [Du,Dv,Mu]],'the ordinary suspension fields are tangent to every value fibre')
zero(sp.Matrix.hstack(Du,Dv,Mu).extract([0,2,3],[0,1,2]).det()+4*b1*x*x,
     'the ordinary suspension has a spanning triple wherever x and B_u are nonzero')
u,v,s=sp.symbols('u v s');F=(s+u*v)**2+v-1
zero(sp.diff(F,u)-2*v*(s+u*v),'the global K2 u translation coefficient is polynomial in the whole value')
zero(sp.diff(F,v)-(2*u*(s+u*v)+1),'the global K2 v translation coefficient is polynomial in the whole value')
check(sp.Poly(F,u,v).coeff_monomial(u*u*v*v)==1,'every K2 specialization has nonconstant base polynomial')
algebra=len(checks);context=json.loads((HERE/'sources/relative-flexibility-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'primary terminology context '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no transitivity theorem or conjectural action classification is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==4 and not flags['one_required_free_action_constructed']
      and not flags['exact_global_kernel_of_lifted_derivation_classified']
      and not flags['S_rho_normal_form_source_alignment_constructed']
      and not flags['original_goal_complete'],
      'collective orbit closure does not assert an individual free quotient,source alignment or completed goal')
proof=HERE/'GLOBAL_FINITE_DIFFERENCE_LIFTS_AND_RELATIVE_FLEXIBILITY.md'
paths=['sources/relative-flexibility-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-relative-flexibility.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
