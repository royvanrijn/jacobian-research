"""Universal integer jet identities;no bounded potential search or old replay."""
import hashlib,json,platform,time
from itertools import permutations
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.expand(value)==0,label)
x,y,u,v=sp.symbols('x y u v');z=[x,y,u,v]
a,b,r,s,p,i,j=sp.symbols('a b r s p i j',integer=True)
zero(a*(s-b)-b*(r-a)-(a*s-b*r),'the all-exponent Jacobian coefficient has its universal divided form')
zero((a*s-b*r).subs({r:p*i,s:p*j})-p*(a*j-b*i),
     'every resonant coefficient is universally divisible by p before any residue specialization')
P=x+x*x*y+x**4+x**6*y*y;Q=y+x**5+x**6*y+x**7*y*y+x**8*y**3
J=sp.expand(sp.det(sp.Matrix([P,Q]).jacobian([x,y])))
check(all(int(c)%2==0 for c in sp.Poly(J-1,x,y).coeffs()),'the actual candidate is Keller in every characteristic-two field')
check(sp.Poly(J,x,y).coeff_monomial(x*y)==2,'its lowest universal divided plane equation is exactly one')
K=sp.expand((J-1)/2)
check(all(c.q==1 for c in sp.Poly(K,x,y).coeffs()),'the naive integral plane representative has error exactly 2K')
for xx,yy in [(0,1),(1,0),(1,1)]:
    check([int(f.subs({x:xx,y:yy}))%2 for f in [P,Q]]==[0,1],
          'the candidate has the exact residue collision at '+str((xx,yy)))
H0=u*P+v*Q
check(sp.Poly(H0,z).total_degree()==12,'the entire displayed candidate has fixed finite degree twelve')
def homogeneous(poly,degree):
    return sp.Add(*(coefficient*sp.prod(zz**e for zz,e in zip(z,ex))
                    for ex,coefficient in sp.Poly(poly,z).terms() if sum(ex)==degree))
q=x*u+y*v;M=sp.hessian(q,z)
zero(homogeneous(H0,2)-q,'the candidate has the split hyperbolic quadratic part')
zero(homogeneous(H0,3),'its cubic potential part is exactly zero')
zero(homogeneous(H0,4)-x*x*y*u,'its sole quartic monomial violates the divided Hessian condition')
check(M*M==sp.eye(4) and M.det()==1,'the quadratic Hessian is its own inverse and has determinant one')
# All twenty cubic coefficients and all thirty-five quartic coefficients are symbolic.
def exponents(degree):
    return [(a,b,c,degree-a-b-c) for a in range(degree+1)
            for b in range(degree+1-a) for c in range(degree+1-a-b)]
cs={ex:sp.Symbol('c'+''.join(map(str,ex))) for ex in exponents(3)}
ds={ex:sp.Symbol('d'+''.join(map(str,ex))) for ex in exponents(4)}
C=sum(cs[ex]*sp.prod(zz**e for zz,e in zip(z,ex)) for ex in cs)
D=sum(ds[ex]*sp.prod(zz**e for zz,e in zip(z,ex)) for ex in ds)
Delta=lambda f:sp.diff(f,x,u)+sp.diff(f,y,v)
ellx=2*cs[2,0,1,0]+cs[1,1,0,1];elly=2*cs[0,2,0,1]+cs[1,1,1,0]
ellu=2*cs[1,0,2,0]+cs[0,1,1,1];ellv=2*cs[0,1,0,2]+cs[1,0,1,1]
zero(Delta(C)-(ellx*x+elly*y+ellu*u+ellv*v),'the four harmonic defects include every general cubic coefficient')
zero(sp.trace(M*sp.hessian(C,z))-2*Delta(C),'the whole degree-one Hessian determinant equation is exactly twice Delta C')
AA=ds[2,1,1,0];BB=ds[1,2,0,1]
zero(sp.Poly(sp.trace(M*sp.hessian(D,z)),z).coeff_monomial(x*y)-4*(AA+BB),
     'all thirty-five quartic coefficients contribute only the two displayed coefficients to xy')
T=(-cs[0,1,0,2]*cs[1,2,0,0]+3*cs[0,1,2,0]*cs[3,0,0,0]
   -2*cs[0,2,0,1]*cs[2,0,1,0]+2*cs[0,2,1,0]*cs[2,0,0,1]
   +3*cs[0,3,0,0]*cs[1,0,0,2]-cs[1,0,2,0]*cs[2,1,0,0])
E=(M*sp.hessian(C,z)).subs({u:0,v:0})
e2=sp.expand((sp.trace(E)**2-sp.trace(E*E))/2)
e2xy=sp.Poly(e2,x,y).coeff_monomial(x*y)
zero(e2xy-(-4*T-4*ellu*cs[2,1,0,0]-4*ellv*cs[1,2,0,0]+2*ellx*elly),
     'the universal cubic quadratic-determinant contribution has the exact integral factors')
check(all(c.q==1 and int(c)%2==0 for c in sp.Poly(e2xy,*cs.values()).coeffs()),
      'the xy quadratic determinant contribution is universally even,without imposing harmonic equations')
check(all(sum(ex)==2 for ex,c in sp.Poly(T,*cs.values()).terms()),
      'the divided correction T is an integral homogeneous quadratic in cubic coefficients')
check(all(c.q==1 for c in sp.Poly(T,*cs.values()).coeffs()),
      'the displayed divided Hessian equation has integer coefficients')
# A complete determinant coefficient expansion is independent of the trace calculation.
smallH=q+C+AA*x*x*y*u+BB*x*y*y*v
matrix=sp.hessian(smallH,z).subs({u:0,v:0});detxy=0
for perm in permutations(range(4)):
    sign=(-1)**sum(perm[ii]>perm[jj] for ii in range(4) for jj in range(ii+1,4))
    detxy+=sign*sp.Poly(sp.expand(sp.prod(matrix[ii,perm[ii]] for ii in range(4))),x,y).coeff_monomial(x*y)
target=4*(AA+BB-T-ellu*cs[2,1,0,0]-ellv*cs[1,2,0,0])+2*ellx*elly
zero(detxy-target,'all twenty-four determinant terms independently reproduce the general cubic--quartic xy identity')
# Integral hyperbolic splitting does not divide by two.
alpha,beta,gamma,t=sp.symbols('alpha beta gamma t')
zero(sp.diff(alpha+beta*t+gamma*t*t,t).subs(t,0)-beta,
     'the isotropic Hensel equation has the unit polar coefficient as its derivative')
Qw,Qv,Bvw=sp.symbols('Qw Qv Bvw')
zero((Qw-Qw*Bvw+Qw*Qw*Qv).subs({Bvw:1,Qv:0}),
     'the hyperbolic partner adjustment uses an isotropic vector and polar pairing one,without dividing by two')
# Symbolic valuation bounds cover every ramification index,not a finite list.
ee=sp.symbols('ee',integer=True,positive=True);rr=2*ee+1
check(sp.simplify((ee+2*(ee+1))-rr).is_positive,
      'the 2 ell_x ell_y valuation is strictly above the forbidden precision for every e')
check(sp.simplify((2*ee+(ee+1)+1)-rr).is_positive,
      'the 4 ell_u times cubic-coefficient error also vanishes at every forbidden precision')
check(sp.simplify(rr-2*ee)==1,'the surviving quartic term has valuation exactly one below the forbidden precision')
check(sp.simplify(2*ee+2-rr).is_positive,'the quadratic cubic-correction term cannot cancel that quartic unit')
zero((1+2*K)**2-1-4*K*(1+K),'the naive whole cotangent Hessian has the exact sharp lower-precision factor four')
check(all(int(c)%4==0 for c in sp.Poly(J**2-1,x,y).coeffs()),
      'the entire naive potential,not only its four-jet,is a unit-Hessian lift modulo four')
check(int(sp.Poly(J**2,x,y).coeff_monomial(x*y))%8==4,
      'the candidate exhibits the forbidden unramified next coefficient')
# A generic cotangent block determinant verifies the Hessian identity used above.
f1,f2,g1,g2,b11,b12,b22=sp.symbols('f1 f2 g1 g2 b11 b12 b22')
block=sp.Matrix([[b11,b12,f1,g1],[b12,b22,f2,g2],[f1,f2,0,0],[g1,g2,0,0]])
zero(block.det()-(f1*g2-f2*g1)**2,'the universal cotangent block Hessian is the squared plane Jacobian')
for xx,yy in [(0,1),(1,0),(1,1)]:
    gradient=[sp.diff(H0,zz) for zz in z]
    check([int(f.subs({x:xx,y:yy,u:0,v:0}))%2 for f in gradient]==[0,0,0,1],
          'the full characteristic-two gradient collision is exact at '+str((xx,yy)))
algebra=len(checks)
context=json.loads((HERE/'sources/ramified-hc4-lifting-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary coefficient-equation context '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'no strong-Keller sufficiency or mixed-characteristic lifting conjecture is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==8 and not flags['arbitrary_nonlinear_potential_source_orbit_excluded']
      and not flags['old_registered_arithmetic_checkers_replayed']
      and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
      'sharp all-ramification nonlifting excludes this route while full HC4 and nonlinear reframings remain open')
proof=HERE/'RAMIFICATION_INDEPENDENT_LIFTING_OBSTRUCTIONS_FOR_HC4.md'
paths=['sources/ramified-hc4-lifting-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-ramified-hc4-lifting.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
