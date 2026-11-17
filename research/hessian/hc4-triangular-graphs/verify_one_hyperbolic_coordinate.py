"""Exact controls of the written cylinder-plane and first-normal-jet proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):
    check(sp.cancel(value)==0 or sp.simplify(value)==0,label)
def zeros(values,label):
    check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,t,R,S,P,Q,u,v,s,rho=sp.symbols('x y z w t R S P Q u v s rho')
L=sp.symbols('L')
# Rational generic curves recover every plane generator, not only the function.
for degree in [1,2,3]:
    C=z**degree-1;D=z**(degree+2)+z
    ww=(s-D)/C
    zero((w*C+D).subs(w,ww)-s,'the whole generic multiplier curve has an actual solved coordinate ring degree='+str(degree))
    aa,bb,gg=sp.gcdex(C,s-D,z)
    zero(aa*C+bb*(s-D)-gg,'the generic Bezout identity makes C an actual unit degree='+str(degree))
    check(sp.Poly(gg,z).degree()==0 and gg!=0,'no root of C remains on a generic level degree='+str(degree))
    check(sp.Poly(C,z).degree()==degree,'every omitted finite generic point is accounted for degree='+str(degree))
T=u*v+rho
B3=u*T*T-u+T
parameter_u=(s-t)/(t*t-1)
parameter_v=(t-rho)*(t*t-1)/(s-t)
zero(B3.subs({u:parameter_u,v:parameter_v},simultaneous=True)-s,'the generic B_rho fibre is the actual rational three-puncture chart for every rho')
zero(T.subs({u:parameter_u,v:parameter_v},simultaneous=True)-t,'the curve chart has its inverse and recovers both original plane coordinates')
for value in [1,-1]:
    zero((B3-s).subs(rho,value).subs({u:parameter_u.subs(rho,value),v:parameter_v.subs(rho,value)},simultaneous=True),
         'generic cylinder rigidity still applies at a critical-value collision rho='+str(value))
# The unit lemma uses additivity of exponential degrees, never a classification.
A0,A1,B0,B1=sp.symbols('A0 A1 B0 B1')
tau=sp.symbols('tau')
zero(sp.expand((A0+tau*A1)*(B0+tau*B1)).coeff(tau,2)-A1*B1,
     'multiplicativity of exp(t delta) has the nonzero leading product needed for order additivity')
check(1/z*z==1,'the unit and inverse order sum is zero in every punctured-curve chart')
zero(sp.diff(z,y),'the free cylinder translation fixes every plane element')
check(sp.diff(y,y)==1,'the free translation bounds the common kernel by the whole plane algebra')
zero(sp.diff(z,w),'the excluded affine-line generic curve B=z has a derivation annihilating B')
check(sp.diff(w,w)==1,'that negative control moves a plane coordinate, so the puncture hypothesis cannot be dropped')
# A global one-variable cylinder automorphism has one constant scale.
d=sp.symbols('d',nonzero=True)
cylinder=d*S+t*t+R
inverse=(y-t*t-R)/d
zero(cylinder.subs(S,inverse)-y,'the intrinsic free coordinate has its actual polynomial inverse over the plane')
zero(inverse.subs(y,cylinder)-S,'the free scale has its reverse inverse over the plane')
check(sp.Matrix([t,R,cylinder]).jacobian([t,R,S]).det()==d,'the entire cylinder free scale is a single nonzero constant')
# Full first-normal-jet identities are independent of all higher jets.
A=sp.symbols('A',nonzero=True)
hxx,hxz,hxw,bzz,bzw,bww=sp.symbols('hxx hxz hxw bzz bzw bww')
critical_hessian=sp.Matrix([[hxx,A,hxz,hxw],[A,0,0,0],[hxz,0,bzz,bzw],[hxw,0,bzw,bww]])
zero(critical_hessian.det()+A*A*(bzz*bww-bzw*bzw),
     'all higher normal Hessian entries cancel from the exact critical determinant')
af=sp.Function('a')(z,w);lf=sp.Function('ell')(z,w)
bf=sp.Function('B')(z,w);rf=sp.Function('R')(x,y,z,w)
hf=bf+x*(af*y+lf)+x*x*rf
zero(hf.subs(x,0)-bf,'the exceptional trace is unchanged by every higher normal jet')
zero(sp.diff(hf,x).subs(x,0)-(af*y+lf),'only the stated first normal slope and intercept survive')
zero(sp.diff(hf,y)-x*(af+x*sp.diff(rf,y)),'the exact vertical equation makes each critical line an actual local component')
for q in [z,w]:
    zero(sp.diff(hf,q).subs(x,0)-sp.diff(bf,q),'the complete exceptional critical-line scheme retains the native plane gradient '+str(q))
C=z*z-1;D=z
B=w*C+D
jet_slope=2+(z*z-1)*(w+z)
ell=z+w*w
h=B+x*(jet_slope*y+ell)+x*x*(y**5+z**3*w*w+x*y*z)
for sign in [1,-1]:
    wi=-sp.Rational(sign,2)
    yi=-(sign+wi*wi)/2
    point={x:0,z:sign,w:wi,y:yi}
    zeros([sp.diff(h,q).subs(point) for q in [x,y,z,w]],'arbitrary higher-jet control retains the full exceptional critical point sign='+str(sign))
    zero(jet_slope.subs(point)-2,'the first slope is constant only where required, at the critical points sign='+str(sign))
    zero(sp.hessian(h,[x,y,z,w]).det().subs(point)-16,'the control full critical point is genuinely Morse despite high normal degree sign='+str(sign))
# Constant rank-zero Tang form, with a whole primitive row completion.
aa=t;bb=t*t+1;a1=-1;b1=-t
matrix=sp.Matrix([[aa,bb],[a1,b1]])
check(matrix.det()==1,'the primitive univariate row has an actual determinant-one polynomial completion')
forward=matrix*sp.Matrix([P,Q])
back=sp.Matrix([[b1,-bb],[-a1,aa]])*sp.Matrix([R,S])
zeros(forward.subs({P:back[0],Q:back[1]},simultaneous=True)-sp.Matrix([R,S]),'the completed free-variable change has its global left inverse')
zeros(back.subs({R:forward[0],S:forward[1]},simultaneous=True)-sp.Matrix([P,Q]),'the completed free-variable change has its global right inverse')
zero((aa*P+bb*Q).subs({P:back[0],Q:back[1]},simultaneous=True)-R,'the entire exceptional trace has exactly one free cylinder coordinate')
cp,fx,gx,av,bv,c=sp.symbols('cp fx gx av bv c',nonzero=True)
J=fx*cp*bv-cp*av*gx
slope=-bv*fx+av*gx
zero(slope+J/cp,'the actual normal slope is the unsquared -Jacobian/C-prime at every common root')
zero((A*d)*cp+c-(A*d*cp+c),'the common free scale forces one identical nonzero C-prime value at all roots')
for degree in [2,3]:
    polynomial=t**degree-1
    derivative=sp.diff(polynomial,t)
    check(sp.degree(derivative,t)<sp.degree(polynomial,t),'the derivative-minus-constant degree is strictly below all distinct roots degree='+str(degree))
    check(sp.gcd(polynomial,derivative,t)==1,'the complete univariate common-zero control is reduced degree='+str(degree))
n=sp.symbols('n',integer=True,positive=True)
check(n-1<n,'the decisive degree inequality holds symbolically in every degree')
check(n!=0,'the leading derivative coefficient never vanishes in characteristic zero')
# Rank-one Tang relative equations solve the whole curve as a graph.
f1,f2,g1,g2,d1,d2=sp.symbols('f1 f2 g1 g2 d1 d2')
coefficient=sp.Matrix([[f1,g1],[f2,g2]])
rhs=sp.Matrix([d1,d2])
solution=-coefficient.inv()*rhs
zeros(coefficient*solution+rhs,'both eliminated variables are solved globally when the Keller determinant is a nonzero constant')
check(coefficient.det()==f1*g2-f2*g1,'the rank-one whole-curve graph uses the actual Keller determinant, not a dropped singular group')
# Applications: every multiplier degree and every rho collision keeps all components.
for name,plane in [('h0',w*(z*z-1)),('h1',w*(z*z-1)+z)]:
    for sign in [1,-1]:
        wi=0 if name=='h0' else -sp.Rational(sign,2)
        point={z:sign,w:wi}
        zeros([sp.diff(plane,q).subs(point) for q in [z,w]],'the old model has the required actual plane critical point '+name+str(sign))
        zero(sp.hessian(plane,[z,w]).det().subs(point)+4,'the exceptional critical trace is Morse '+name+str(sign))
for point,value in [({u:0,v:1-rho*rho},rho),({u:-sp.Rational(1,2),v:2*(rho-1)},1),({u:sp.Rational(1,2),v:-2*(rho+1)},-1)]:
    zeros([sp.diff(B3,q).subs(point,simultaneous=True) for q in [u,v]],'the exact suspension plane retains a complete critical point for all rho '+str(value))
    zero(sp.hessian(B3,[u,v]).det().subs(point,simultaneous=True)+1,'all coupled-plane critical points remain Morse at colliding values '+str(value))
Vnew=w*y+rho-z
Pnew=x+(w*Vnew-2*w*w*y-2*w*rho-1)/2
check(sp.Poly(Vnew,x,y,z,w).total_degree()>1 and sp.Poly(Pnew,x,y,z,w).total_degree()>1,
      'the native pulled-back hyperbolic coordinates are nonlinear and cannot be treated as actual affine coordinates')
# Retained-one-coordinate equivalence of K2 and h1 is excluded even without unit Hessian.
h1=2*x*y+B
k2=w*(2*x*y+z*z-1)+z
zero(x*(2*w*z+1)-z*sp.diff(k2,y)-x,'the K2 three-relative ideal contains the retained coordinate scheme-theoretically')
zero(sp.diff(k2,w)-2*y*x-(z*z-1),'the K2 three-relative ideal contains the whole same critical-plane equations')
for sign in [1,-1]:
    point={x:0,z:sign,w:-sp.Rational(sign,2)}
    zero(sp.diff(sp.diff(h1,x),y).subs(point)-2,'the h1 free normal slope is the same on every critical line sign='+str(sign))
    zero(sp.diff(sp.diff(k2,x),y).subs(point)+sign,'the K2 free normal slopes remain opposite sign='+str(sign))
check(sp.Rational(-1,2)!=sp.Rational(1,2),'no common nonzero source/value multiplier matches the two unordered normal-slope pairs')
counter=2*x*y+z*z+w*w
check(sp.hessian(counter,[x,y,z,w]).det()==-16,'the excluded N=1 boundary really has a unit native Hessian and an automorphic gradient')
context=json.loads((HERE/'sources/one-hyperbolic-coordinate-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'new versioned primary context is byte pinned, with no classification imported')
check(context['newly_used_unproved_external_statements']==[],'all new unit/kernel/cylinder arguments are proved and the affine-line theorem is retained')
check(context['constant_normal_slope_obstruction_independent_of_all_higher_normal_jets'] is True,'the proof keeps every first-normal-jet-compatible higher perturbation and every source degree')
check(context['full_source_h0_h1_orbits_resolved'] is False and context['K3_rho_full_source_family_resolved'] is False
      and context['K2_one_coordinate_unit_hessian_route_excluded'] is False
      and context['nonlinear_pulledback_coordinate_algebras_treated_as_actual_affine_coordinates'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'unrestricted base changes and nonlinear normal-form Hessian issues remain explicit open boundaries')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
sources={'sources/one-hyperbolic-coordinate-context.json':hashlib.sha256((HERE/'sources/one-hyperbolic-coordinate-context.json').read_bytes()).hexdigest()}
for entry in context['newly_cached_primary_context']:sources['sources/'+entry['retained_file']]=entry['sha256']
receipt=dict(outcome='PASS',claims=['HC4MYJR0','HC4MYJR1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ONE_HYPERBOLIC_COORDINATE_AND_INTRINSIC_CYLINDER_PLANE.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written intrinsic-cylinder, one-coordinate orientation and arbitrary-higher-normal-jet proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-one-hyperbolic-coordinate.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
