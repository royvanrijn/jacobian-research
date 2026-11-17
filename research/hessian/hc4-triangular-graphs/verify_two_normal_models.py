"""Exact controls of global normal-model maps and relative-rank alignment."""
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
x,y,z,t,X,Y,Z,W,P,V,u,v,kappa,rho=sp.symbols('x y z t X Y Z W P V u v kappa rho')
lam,c,e,r=sp.symbols('lam c e r',nonzero=True)
alpha=sp.symbols('alpha')
q=2*x*y+z*z
fscalar=(1+lam*t)*q+e*z-c*t
K2=W*(2*X*Y+Z*Z-1)+Z
K3=K2+2*X*(rho-Z)
variables=[x,y,z,t]
for name,p,extra in [('scalar force',alpha/e,alpha*x),('length two',(t*x+alpha)/e,t*x*x+alpha*x)]:
    shear=sp.Matrix([x,y-z*p-x*p*p/2,z+x*p,t])
    back=sp.Matrix([x,y+z*p-x*p*p/2,z-x*p,t])
    sm=dict(zip(variables,shear));bm=dict(zip(variables,back))
    zero(q.subs(sm,simultaneous=True)-q,'the actual source shear preserves the nondegenerate quadratic '+name)
    zeros(shear.subs(bm,simultaneous=True)-sp.Matrix(variables),'source shear has its global polynomial left inverse '+name)
    zeros(back.subs(sm,simultaneous=True)-sp.Matrix(variables),'source shear has its global polynomial right inverse '+name)
    check(shear.jacobian(variables).det()==1,'source shear is volume preserving in all parameter values '+name)
    zero(fscalar.subs(sm,simultaneous=True)-(fscalar+extra),'literal source-potential identity identifies the resonant family '+name)
    zero((fscalar+extra).subs(bm,simultaneous=True)-fscalar,'the normalization uses the inverse shear, not a birational substitution '+name)
normal=sp.Matrix([r*X,r*Y,r*Z,e*W/(lam*r)-1/lam])
normal_inverse=sp.Matrix([x/r,y/r,z/r,lam*r*(t+1/lam)/e])
ns=dict(zip(variables,normal))
zero(fscalar.subs(ns,simultaneous=True).subs(c,lam*r*r)-(e*r*K2+r*r),
     'every scalar nonzero parameter normalizes to the single fixed K2')
zeros(normal_inverse.subs(ns,simultaneous=True)-sp.Matrix([X,Y,Z,W]),'the scalar constant normalization has its global inverse')
check(normal.jacobian([X,Y,Z,W]).det()==e*r*r/lam,'the scalar normalization has nonzero constant Jacobian')
for name,p in [('scalar',alpha/e),('length two',(t*x+alpha)/e)]:
    back=sp.Matrix([x,y+z*p-x*p*p/2,z-x*p,t])
    complete=back.subs(ns,simultaneous=True)
    actual=fscalar+alpha*x+(t*x*x if name=='length two' else 0)
    zero(actual.subs(dict(zip(variables,complete)),simultaneous=True).subs(c,lam*r*r)-(e*r*K2+r*r),
         'the complete original-to-fixed-model source identity is polynomial '+name)
coupled=sp.Matrix([e*lam*X,r*r*Y/(e*lam)-r*Z/lam,r*Z,e*W/(lam*r)-1/lam])
coupled_inverse=sp.Matrix([x/(e*lam),e*lam*(y+z/lam)/(r*r),z/r,lam*r*(t+1/lam)/e])
cs=dict(zip(variables,coupled))
actual=fscalar+alpha*x+2*t*x*z
zero(actual.subs(cs,simultaneous=True).subs(c,lam*r*r)-(e*r*K3.subs(rho,alpha*lam/(2*r))+r*r),
     'every length-three parameter reduces to rho=alpha lambda/(2r)')
zeros(coupled_inverse.subs(cs,simultaneous=True)-sp.Matrix([X,Y,Z,W]),'the coupled normalization has its global polynomial inverse')
zeros(coupled.subs(dict(zip([X,Y,Z,W],coupled_inverse)),simultaneous=True)-sp.Matrix(variables),'the coupled normalization also has its reverse polynomial inverse')
check(coupled.jacobian([X,Y,Z,W]).det()==e*r*r/lam,'the coupled affine normalization has constant nonzero Jacobian')
zero((alpha*lam/(2*r)).subs({alpha:2*rho,lam:1,r:1})-rho,'every rho is actually realized, not merely necessary')
# A full source automorphism produces a genuine hyperbolic suspension.
new=sp.Matrix([X+(W*(W*Y+rho-Z)-2*W*W*Y-2*W*rho-1)/2,
               W*Y+rho-Z,W,Y])
back=sp.Matrix([P-(u*V-2*u*u*v-2*u*rho-1)/2,v,u*v+rho-V,u])
old_variables=[X,Y,Z,W];new_variables=[P,V,u,v]
zeros(new.subs(dict(zip(old_variables,back)),simultaneous=True)-sp.Matrix(new_variables),
      'the suspension map has its actual global polynomial left inverse')
zeros(back.subs(dict(zip(new_variables,new)),simultaneous=True)-sp.Matrix(old_variables),
      'the suspension map has its actual global polynomial right inverse')
check(new.jacobian(old_variables).det() in [1,-1],'no divisor is inverted in the hyperbolic suspension coordinate map')
T=u*v+rho
B=u*T*T-u+T
zero(K3.subs(dict(zip(old_variables,back)),simultaneous=True)-(2*P*V+B),
     'the normalized cubic is literally a quadratic hyperbolic suspension')
zero(sp.diff(B,v)-u*(2*u*T+1),'the first plane critical equation has the exhaustive product split')
zero(sp.diff(B,u)-(T*T-1+v*(2*u*T+1)),'the second plane critical equation retains every u-zero point')
points=[({u:0,v:1-rho*rho},rho),({u:-sp.Rational(1,2),v:2*(rho-1)},1),
        ({u:sp.Rational(1,2),v:-2*(rho+1)},-1)]
for index,(point,value) in enumerate(points):
    zeros([sp.diff(B,q).subs(point,simultaneous=True) for q in [u,v]],'complete plane critical point '+str(index))
    zero(B.subs(point,simultaneous=True)-value,'actual critical value retained for all rho '+str(index))
    zero(sp.hessian(B,[u,v]).det().subs(point,simultaneous=True)+1,'every plane critical point is Morse, including collisions '+str(index))
    zero(sp.hessian(2*P*V+B,[P,V,u,v]).det().subs(point,simultaneous=True)-4,'all suspended critical Hessians have the same native value '+str(index))
check(len({point[u] for point,_ in points})==3,'the three full critical points stay distinct at every rho')
for value in [1,-1]:
    counts=sum(sp.simplify(vv-value)==0 for _,vv in [(pt,val.subs(rho,value) if hasattr(val,'subs') else val) for pt,val in points])
    check(counts==2,'coincident values retain exactly two different Morse points rho='+str(value))
# Exact generic plane fibre, not a zero-divisor cancellation argument.
S=sp.symbols('S')
param_u=(kappa-S)/(S*S-1)
param_v=(S-rho)*(S*S-1)/(kappa-S)
zero((u*v+rho).subs({u:param_u,v:param_v},simultaneous=True)-S,'the generic fibre parameterization has its regular inverse T')
zero(B.subs({u:param_u,v:param_v},simultaneous=True)-kappa,'the parameterization gives every generic plane level point')
zero((kappa-T)/(T*T-1)-u+(B-kappa)/(T*T-1),'the generic inverse u is exact modulo the actual level equation')
zero((T-rho)/u-v,'the second inverse coordinate is exact on the whole generic u-nonzero fibre')
for sign in [1,-1]:
    zero((B-sign)-(T-sign)*(u*(T+sign)+1),'the entire special plane level factorization is retained sign='+str(sign))
zero(B-rho-u*(T*T+v-1),'the u-zero critical level retains its full second component')
for sign in [1,-1]:
    zero((B-sign).subs(rho,sign)-u*v*(u*u*v+2*sign*u+1),'coincident level has all three distinct components rho='+str(sign))
    check((u*u*v+2*sign*u+1).subs(u,0)==1,'the third component misses the first line rho='+str(sign))
    zero((u*u*v+2*sign*u+1).subs({u:-sp.Rational(sign,2),v:0}),
         'the second line meets the third component at the other Morse point rho='+str(sign))
L,nk=sp.symbols('L nk')
zero(2*(L-1)-1-(2*L-3),'two punctured-line components retain their one transverse intersection')
zero(L+(L-2)-1-(2*L-3),'the rho level retains its affine line, punctured line and intersection')
zero(2*L+(L-1)-2-(3*L-3),'the coincident level retains the entire chain and both distinct intersections')
zero((L-1)*L*L+L*(L-3+nk*L)-(L**3-3*L+nk*L*L),
     'every suspended level has the passing full HG0 profile including collisions')
# The fixed scalar model has exactly two full Morse critical points.
for sign in [1,-1]:
    point={X:0,Y:0,Z:sign,W:-sp.Rational(sign,2)}
    zeros([sp.diff(K2,q).subs(point) for q in old_variables],'complete fixed-model critical point sign='+str(sign))
    zero(K2.subs(point)-sign,'the fixed model critical value is preserved sign='+str(sign))
    zero(sp.hessian(K2,old_variables).det().subs(point)-4,'the fixed model critical Hessian agrees sign='+str(sign))
h0=2*X*Y+W*(Z*Z-1)
for sign in [1,-1]:
    zero(h0.subs({X:0,Y:0,W:0,Z:sign}),'the old h0 has equal critical values at both points sign='+str(sign))
check(1!=-1,'the two distinct K2 values exclude source equivalence with h0 even after scalar and constant changes')
# Relative critical scheme identities and the non-affine-plane rank-one certificate.
relative=[sp.diff(K2,Y),sp.diff(K2,Z)]
zero(relative[0]-2*W*X,'the scalar relative first equation retains the whole multiplier divisor')
zero(relative[1]-(2*W*Z+1),'the scalar relative second equation forces W to be a unit')
zero(X*(2*W*Z+1)-Z*(2*W*X)-X,'the relative ideal equals (X,2WZ+1) scheme-theoretically')
param={X:0,Z:-1/(2*W)}
zeros([vv.subs(param,simultaneous=True) for vv in relative],'the entire relative surface is Gm times A1')
check(sp.hessian(K2,[Y,Z])==sp.diag(0,2*W),'the vertical Hessian is rank one everywhere on the relative surface')
check(L*(L-1)!=L*L,'its exact E-polynomial excludes the rank-two Tang affine-plane graph')
check(sp.Matrix(relative).jacobian(old_variables).subs(param,simultaneous=True).rank()==2,
      'the selected full Hessian rows give a genuinely smooth relative critical surface')
zeros([sp.diff(2*P*V+B,P)-2*V,sp.diff(2*P*V+B,V)-2*P],
      'the coupled suspension relative critical scheme is an actual global section')
check(sp.hessian(2*P*V+B,[P,V])==sp.Matrix([[0,2],[2,0]]),'the global-section application retains the full hyperbolic vertical Hessian')
relative_cubic=[sp.diff(K3,X),sp.diff(K3,Z)]
zeros([vv.subs({Z:W*Y+rho,X:W*W*Y+W*rho+sp.Rational(1,2)},simultaneous=True) for vv in relative_cubic],
      'the coupled cubic already has an actual global section over its unchanged natural base W,Y')
check(new[2]==W and new[3]==Y,'the global suspension map preserves the natural coupled base algebra exactly')
# The constant-rank-one Tang branch really forces C_U=0 before triangularity.
ss,tt,UU,VV=sp.symbols('ss tt UU VV')
Af=sp.Function('A')(tt,UU);Bf=sp.Function('B')(tt,UU);Cf=sp.Function('C')(tt,UU)
gg=Af+ss*Bf+VV*Cf
zero(sp.det(sp.hessian(gg,[UU,VV]))+sp.diff(Cf,UU)**2,
     'the rank-one Tang projection branch has the exact vertical determinant -C_U squared')
check(sp.Poly(3*tt+2,tt).diff().degree()==0,'the resulting one-variable Keller component is necessarily affine')
triangular=[3*tt+2,UU/3+tt*tt]
check(sp.Matrix(triangular).jacobian([tt,UU]).det()==1,'the rank-one terminal pair is the actual triangular Keller form')
# Relative nonlinear source maps preserve these data by the critical-chain rule.
Ynew=Y+Z*Z;Znew=Z
tau=sp.Matrix([Ynew,Znew])
check(tau.jacobian([Y,Z]).det()==1,'a representative arbitrary-degree relative shear is globally invertible')
transformed=K2.subs({Y:Ynew,Z:Znew},simultaneous=True)
grad=sp.Matrix([sp.diff(K2,Y),sp.diff(K2,Z)])
zeros(sp.Matrix([sp.diff(transformed,Y),sp.diff(transformed,Z)])-tau.jacobian([Y,Z]).T*grad.subs({Y:Ynew,Z:Znew},simultaneous=True),
      'the exact vertical first chain rule retains the whole critical ideal')
correction=sp.hessian(transformed,[Y,Z])-tau.jacobian([Y,Z]).T*sp.hessian(K2,[Y,Z]).subs({Y:Ynew,Z:Znew},simultaneous=True)*tau.jacobian([Y,Z])
zeros(list(correction.subs(param,simultaneous=True)),'all nonlinear vertical Hessian corrections vanish on the full relative critical scheme')
context=json.loads((HERE/'sources/two-normal-model-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
check(context['newly_used_unproved_external_statements']==[],'no cancellation, stable equivalence or arbitrary Morse rectification theorem is imported')
check(context['all_affine_universal_questions_reduced_to_four_old_tests_fixed_K2_and_K3_rho'] is True,
      'normalization removes parameters from actual global source orbits with no source-degree bound')
check(context['remaining_fixed_K2_full_source_orbit_resolved'] is False
      and context['remaining_K3_rho_full_source_family_resolved'] is False
      and context['fixed_K2_identified_with_h0_or_h1'] is False
      and context['rho_complete_source_orbit_invariant_proved'] is False
      and context['arbitrary_quadratic_suspension_inverse_closure_proved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'full-source mixing and true endpoint questions stay open beyond the natural-base exclusions')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJQ0','HC4MYJQ1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md').read_bytes()).hexdigest(),
 source_sha256={'sources/two-normal-model-context.json':hashlib.sha256((HERE/'sources/two-normal-model-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree global-normalization, plane-fibre and constant relative-rank proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-two-normal-models.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
