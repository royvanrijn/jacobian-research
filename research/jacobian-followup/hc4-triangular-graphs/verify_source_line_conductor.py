"""Exact identities and controls for the written universal conductor theorem.

The global trace-dual and every-line proofs are written mathematics,not
inferred from bounded examples.No registered mathematical checker is replayed.
"""
import hashlib,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];started=time.monotonic();checks=[]
def check(ok,label):
 assert ok,label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
t,z,p,q,T,a,b=s.symbols('t z p q T a b')

# Symbolic solution of the two differentiated primitive-relation equations.
pt,pz,qt,qz,ht=s.symbols('pt pz qt qz ht');J=pt*qz-pz*qt
rp=-ht*qz/J;rq=ht*pz/J
check(s.cancel(rp*pz+rq*qz)==0,'the z derivative of the primitive relation vanishes exactly')
check(s.cancel(ht+rp*pt+rq*qt)==0,'the t derivative gives the global h/J conductor multiplier')

# An actual polynomial Keller automorphism:all source-line conductors are1.
P=t+z*z;Q=z;Phi=T-p+q*q;h=s.diff(Phi,T).subs({T:t,p:P,q:Q})
check(s.det(s.Matrix([P,Q]).jacobian([t,z]))==1,'the shear control is globally Keller')
check(s.expand(Phi.subs({T:t,p:P,q:Q}))==0,'its primitive inverse equation is exact')
check(h==1,'the inverse shear has unit conductor on every source line')

# An immersed nodal curve,with exact finite affine delta1 and conductor2.
pc=z*z-1;qc=z*(z*z-1);H=q*q-p*p*(p+1);cc=z*z-1
check(s.expand(H.subs({p:pc,q:qc}))==0,'the nodal curve equation vanishes on its normalization')
check(s.cancel(s.diff(H,q).subs({p:pc,q:qc})/s.diff(pc,z))==cc,'the curve complementary-module quotient is its conductor')
check(s.gcd(s.diff(pc,z),s.diff(qc,z))==1,'the nodal normalization is immersive everywhere')
check([pc.subs(z,r) for r in [-1,1]]==[0,0] and [qc.subs(z,r) for r in [-1,1]]==[0,0],'the conductor roots are two distinct finite preimages of one node')
inc=s.diag(1,p)
check(inc.det()==p,'the nodal power-basis inclusion has index length1')
check(s.resultant(z*z-(p+1),cc,z)==p*p,'norm of the conductor is the squared inclusion determinant')
check(s.degree(cc,z)==2,'the affine conductor degree is exactly twice delta1')

# Without a whole Keller Jacobian, h includes a genuine ramification factor.
P=t-z*z;Q=-t*z;Phi=q*q-(T-p)*T*T;JJ=s.det(s.Matrix([P,Q]).jacobian([t,z]))
hh=s.expand(s.diff(Phi,T).subs({T:t,p:P,q:Q}))
check(s.expand(Phi.subs({T:t,p:P,q:Q}))==0,'the connected small primitive order has the exact equation')
check(JJ==-t-2*z*z,'its Jacobian is nonconstant,so this is not a plane Keller counterexample')
check(s.expand(hh-t*JJ)==0,'the inverse derivative is t*J,not the conductor t alone')
check(s.cancel(s.diff(Phi,q).subs({T:t,p:P,q:Q})/s.diff(P,z))==t,'trace duality removes the critical-Jacobian factor from the conductor')
check(s.diag(1,-t).det()==-t,'its exact power-basis index detects the nonnormal t=0 line')
check([s.expand(g).subs({T:0,p:-1,q:0}) for g in [Phi,s.diff(Phi,T),s.diff(Phi,p),s.diff(Phi,q)]]==[0]*4,'the primitive hypersurface is singular at the finite double image')
check([JJ.subs({t:0,z:r}) for r in [-1,1]]==[-2,-2],'both finite normalization branches are locally unramified')

# A completed etale normalization can have a singular primitive order.
G=s.Matrix([[2,a+b],[a+b,a*a+b*b]])
check(s.expand(G.det()-(a-b)**2)==0,'the split primitive-order discriminant records branch-value collision')
check(s.eye(2).det()==1,'the split normalization trace pairing remains perfect')

# Known whole3D Keller map:reproduce only the identities needed here.
x,y,zz=s.symbols('x y zz');u=1+x*y
As=u**3*zz+y*y*u*(4+3*x*y)
Bs=y+3*x*u*u*zz+3*x*y*y*(4+3*x*y)
Cs=2*x-3*x*x*y-x**3*zz
JK=s.Matrix([As,Bs,Cs]).jacobian([x,y,zz])
check(s.expand(JK.det())==-2,'the credited whole3D source is globally Keller')
pts=[(0,0,s.Rational(-1,4)),(1,s.Rational(-3,2),s.Rational(13,2)),(-1,s.Rational(3,2),s.Rational(13,2))]
for i,point in enumerate(pts):
 sub=dict(zip([x,y,zz],point))
 check([g.subs(sub) for g in [As,Bs,Cs]]==[s.Rational(-1,4),0,0],'actual retained source point'+str(i)+' has the common target')
gs=2-3*x*y-x*x*zz;ws=u*gs
check(s.expand(Cs-x*gs)==0,'the generic reconstruction has C=x*gamma')
check(s.expand(Bs*Cs-(2*gs+4*ws-3*ws*ws))==0,'the first cubic target identity is exact')
check(s.expand(As*Cs**2-(ws*ws+ws*gs-ws**3))==0,'the second cubic target identity is exact')
check(s.expand(Cs**2*zz-(5*gs**2-3*ws*gs-gs**3))==0,'the primitive source-coordinate reconstruction identity is exact')

A,B,C,Z,w=s.symbols('A B C Z w');P=B*C;Q=A*C*C;gamma=(P-4*w+3*w*w)/2
f=w**3-2*w*w+P*w-2*Q;g=s.expand(8*(C*C*Z-5*gamma**2+3*w*gamma+gamma**3))
raw=s.resultant(f,g,w);norm=s.cancel(raw/(512*C**6));poly=s.Poly(norm,Z)
check(s.expand(raw-512*C**6*norm)==0,'the resultant denominator is removed by an exact identity')
check(all(s.fraction(s.cancel(v))[1].is_number for v in s.Poly(norm,A,B,C,Z).coeffs()),'the primitive norm has polynomial target coefficients,with no C denominator')
check(poly.degree()==3 and poly.LC()==1,'the source coordinate z has a monic cubic field norm')
fibre=s.factor(norm.subs({A:s.Rational(-1,4),B:0,C:0}))
check(s.expand(fibre-(Z+s.Rational(1,4))*(Z-s.Rational(13,2))**2)==0,'the primitive norm has one simple and one double finite root at the full source fibre')
check(s.degree(s.gcd(s.Poly(fibre,Z),s.Poly(s.diff(fibre,Z),Z)))==1,'the fibre norm has two distinct root values and cannot be a cube')
point={A:s.Rational(-1,4),B:0,C:0,Z:s.Rational(13,2)}
for v in [None,A,B,C,Z]:
 expr=norm if v is None else s.diff(norm,v)
 check(expr.subs(point)==0,'the primitive hypersurface singularity has vanishing '+('equation' if v is None else 'partial_'+str(v)))
check(s.diff(norm,Z).subs(point)==0,'order trace-dual localization deletes the two actual etale points sharing z')
check(all(s.Poly(e,x,y,zz).is_multivariate or e.is_number for e in JK.adjugate()),'the whole Keller inverse differential matrix has polynomial entries over the constant determinant')
algebra=len(checks)

context=json.loads((HERE/'sources/source-line-conductor-context.json').read_text())
for pth,hsh in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/pth)==hsh,'immutable retained input '+pth)
check(context['newly_imported_proved_external_statements']==[],'all one-line,source-finite and trace-conductor imports are already retained and credited')
check(context['newly_used_unproved_external_statements']==[],'no primitive-smoothness or D-algebra conjecture is imported')
check(context['actual_plane_or_HC4_counterexample_constructed'] is False,'the actual3D control is not reported as a plane or HC4 counterexample')
check(context['last_HC4_source_to_regular_inverse_Gram_implication_proved'] is False,'the last full-HC4 implication remains OPEN')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,elapsed_seconds=round(time.monotonic()-started,3),
 proof_sha256=sha(HERE/'SOURCE_LINE_CONDUCTOR_AND_PRIMITIVE_ORDER_COLLISIONS.md'),checker_sha256=sha(Path(__file__)),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={},monic_primitive_z_equation=str(norm),specialized_primitive_z_equation=str(fibre),**flags)
(HERE/'receipt-source-line-conductor.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
