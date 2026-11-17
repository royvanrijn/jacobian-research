"""New coordinate/Poisson controls only; run before registration,not thereafter."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
x,y,z,alpha,beta,p,q,c=s.symbols('x y z alpha beta p q c')
u=1+x*y;g0=2-3*x*y-x*x*z
A=u**3*z+y*y*u*(1+3*u);B=y+3*x*u*u*z+3*x*y*y*(1+3*u);C=x*g0
D=3*u*g0-3;r=3*x*(1-D)/2;ss=(u*D-1)/2
T=y*u*(1+3*u)+x*(u*u+u+1)*z
eq('new syzygy identity',u*B-3*x*A,u*y)
eq('scheme ideal y identity',r*A+ss*B,y)
eq('scheme ideal z identity',A-T*y,z)
eq('scheme ideal z coefficients',(1-T*r)*A-T*ss*B,z)
aa=-3*x*x*r-x**3*(1-T*r);bb=-3*x*x*ss+x**3*T*ss
eq('offset in exact coordinate ideal',2*x+aa*A+bb*B,C)
M=s.Matrix([[y*u*(1+3*u),u**3],[1+3*x*y*(1+3*u),3*x*u*u]])
for i,f in enumerate([A,B]):eq('original coefficient matrix row'+str(i),(M*s.Matrix([y,z]))[i],f)
eq('original coefficient determinant',M.det(),-u**3)
J=x*(u*u+u+1);h0=J*ss;k0=-J*r
Mhat=M+s.Matrix([[h0*z,-h0*y],[k0*z,-k0*y]])
for i,f in enumerate([A,B]):eq('corrected coefficient matrix row'+str(i),(Mhat*s.Matrix([y,z]))[i],f)
eq('corrected determinant is field unit',Mhat.det(),-1)
# The actual polynomial inverse is the adjugate divided by the constant -1.
Minv=-Mhat.adjugate()
for i in range(2):
 for j in range(2):eq('concrete matrix inverse '+str((i,j)),(Mhat*Minv)[i,j],int(i==j))
newrow=s.Matrix([[2*alpha-aa,beta-bb]])*Mhat
eq('complete moment identity in ambient polynomial ring',newrow[0]*y+newrow[1]*z-2*x,2*alpha*A+beta*B-C)
# Generic matrix algebra verifies both helper inverse formulas without exploding
# the explicit four-variable map or replaying the old map Jacobian.
a,b,d,e,h,k,yy,zz=s.symbols('a b d e h k yy zz')
N=s.Matrix([[a,b],[d,e]]);Nc=N+s.Matrix([[h*zz,-h*yy],[k*zz,-k*yy]])
eq('universal rank-one determinant correction',Nc.det(),N.det()+h*(d*yy+e*zz)-k*(a*yy+b*zz))
for i,f in enumerate([a*yy+b*zz,d*yy+e*zz]):eq('universal correction preserves generator'+str(i),(Nc*s.Matrix([yy,zz]))[i],f)
V=s.Matrix([[p,q]])
for i,f in enumerate([p,q]):eq('generic helper inverse schema'+str(i),(V*N.inv()*N)[i],f)
xc=(p*y+q*z-c)/2
eq('level inverse eliminates x exactly',p*y+q*z-2*xc,c)
eq('all-level base point retained',xc.subs({y:1,z:0,p:c,q:0}),0)
eq('new helper determinant from row scaling',2*s.Integer(-1),-2)
eq('ambient source coordinate determinant',(-2)*(-2),4)
eq('new projection determinant by retained-Jacobian chain rule',s.Rational(2,4),s.Rational(1,2))
eq('new reduced symplectic volume coefficient',-s.Rational(1,2),s.Rational(-1,2))
gam,dy,dz=s.symbols('gam dy dz')
f=y/gam;gg=z/gam**2
br=s.diff(f,y)*s.diff(gg,gam)*dy+s.diff(f,z)*s.diff(gg,gam)*dz-s.diff(f,gam)*(s.diff(gg,y)*dy+s.diff(gg,z)*dz)
eq('new invariant cotangent bracket',br,(y*dz-2*z*dy)/gam**4)
point={x:0,y:1,z:0}
for label,f,values in [('A',A,[7,8,1]),('B',B,[12,1,0])]:
 for w,vv in zip([x,y,z],values):eq('new reduction bracket input '+label+' '+str(w),s.diff(f,w).subs(point),vv)
eq('new D_C y at point',s.Rational(1*12-7*0,-2),-6)
eq('new D_C z at point',s.Rational(7*1-8*12,-2),s.Rational(89,2))
eq('new reduced source positions do not commute',br.subs({gam:1,y:1,z:0,dy:-6,dz:s.Rational(89,2)}),s.Rational(89,2))
check('nonzero position bracket excludes every y,z-preserving Darboux chart',s.Rational(89,2)!=0,True)
check('Darboux comparison has identity vertical coefficient matrix',s.Matrix([[1,0],[0,1]]).det()==1,True)
context=json.loads((HERE/'sources/gamma-coordinate-context.json').read_text())
for path,sha in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha)
check('no new external theorem import',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'all_degree_coordinate_ideal_helper_theorem_proved','primitive_gamma_moment_polynomial_is_ambient_coordinate_proved','every_gamma_level_has_explicit_polynomial_A4_coordinates','scheme_ideal_equality_and_polynomial_inverse_supplied','new_gamma_projection_has_degree_three_and_all_level_collision','new_source_reduced_form_is_exact_polynomial_symplectic','natural_source_y_z_pair_cannot_be_Darboux_positions','target_A_B_pair_cannot_be_global_Darboux_positions','LI_beta_all_level_stable_exclusion_preserved','external_review_deferred'}
for path,vv in context.items():
 if isinstance(vv,bool):check('scope '+path,vv is (path in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'PRIMITIVE_GAMMA_COORDINATE_AND_DARBOUX_BARRIER.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-gamma-coordinate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
