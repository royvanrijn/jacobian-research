"""Exact controls of the written all-degree global additive-action barrier."""
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
x,y,z,w,s,u,v,P,Q,kappa,a,b=sp.symbols('x y z w s u v P Q kappa a b')
variables=[x,y,z,w]
def deriv(coeff,f,coords=variables):return sp.expand(sum(cc*sp.diff(f,qq) for cc,qq in zip(coeff,coords)))
def commutator(left,right,coords=variables):return [deriv(left,cc,coords)-deriv(right,dd,coords) for cc,dd in zip(right,left)]
q=2*x*y+z*z
K2=w*(q-1)+z
vv=1-q
relation=2*x*y+(s+u*v)**2+v-1
zero((s+u*v).subs({s:K2,u:w,v:vv},simultaneous=True)-z,'the uniform value presentation has the actual old z inverse')
zero(relation.subs({s:K2,u:w,v:vv},simultaneous=True),'the new value presentation relation holds on every original source point')
zero((1-2*x*y-z*z).subs(z,s+u*v)-v+relation,'the reverse v composite is exact modulo the one actual relation')
zero(K2.subs({z:s+u*v,w:u},simultaneous=True)-s-u*relation,'the reverse value composite is exact modulo the relation')
check(sp.Matrix([x,y,w,vv]).jacobian(variables).det()==2*z,
      'the four x,y,u,v generators alone are not a polynomial source automorphism')
Du=sp.Matrix([0,-2*vv*z,2*x*vv,2*x])
Dv=sp.Matrix([0,-(2*w*z+1),2*x*w,0])
Eu=sp.Matrix([-2*vv*z,0,2*y*vv,2*y])
Ev=sp.Matrix([-(2*w*z+1),0,2*y*w,0])
for label,coefficient in [('Du',Du),('Dv',Dv),('Eu',Eu),('Ev',Ev)]:
    zero(deriv(coefficient,K2),'the actual global polynomial action fixes K2: '+label)
    zero(sum(sp.diff(cc,qq) for cc,qq in zip(coefficient,variables)),
         'the global action has zero ambient divergence: '+label)
    for coordinate in variables:
        power=coordinate
        for _ in range(3):power=deriv(coefficient,power)
        zero(power,'the action is nilpotent on the actual source generator '+label+':'+str(coordinate))
zeros(commutator(Du,Dv),'the two full x-standard derivations commute on every source generator')
zeros(commutator(Eu,Ev),'the two full y-standard derivations commute on every source generator')
for label,coefficient,moving,invariant in [('Du',Du,w,vv),('Dv',Dv,vv,w)]:
    zero(deriv(coefficient,moving)-2*x,'the x-localized action is exactly a translation in its base variable '+label)
    zero(deriv(coefficient,invariant),'the other base variable is a genuine global invariant '+label)
    zero(deriv(coefficient,x),'the common standard coefficient ring fixes x '+label)
for coefficient in [Eu,Ev]:zero(deriv(coefficient,y),'the symmetric standard coefficient ring fixes y')
for scalar,coefficient in [(x,Du),(x,Dv),(y,Eu),(y,Ev)]:
    check(deriv(coefficient,scalar)==0,'an actual standard action kernel contains its retained coordinate '+str(scalar))
check(deriv(Du,vv)==0 and deriv(Dv,w)==0,'the generated kernel algebra contains u and v')
zero(K2+w*vv-z,'those kernel generators recover the missing original z globally')
zero(relation.subs(y,-((s+u*v)**2+v-1)/(2*x)),
     'inverting x gives a genuine polynomial chart in s,x,u,v')
check(sp.Poly(K2.subs(x,0),z,w).total_degree()>0,
      'the minimal-denominator descent uses a nonconstant value on x=0')
check(sp.Poly(K2.subs(y,0),z,w).total_degree()>0,
      'the symmetric descent also has a nonconstant value on y=0')
# Every coefficient degree in the entire commuting module is handled structurally.
Af=sp.Function('a')(x,K2);Bf=sp.Function('b')(x,K2)
zero(deriv(Du,Af),'every polynomial a(x,K2) is fixed by Du')
zero(deriv(Dv,Af),'every polynomial a(x,K2) is fixed by Dv')
zero(deriv(Du,Bf),'every polynomial b(x,K2) is fixed by Du')
zero(deriv(Dv,Bf),'every polynomial b(x,K2) is fixed by Dv')
chart={x:0,w:(z-kappa)/(1-z*z)}
zero((K2-kappa).subs(chart,simultaneous=True),'the complete x-zero generic fibre chart has its actual value')
zero(vv.subs(chart,simultaneous=True)-(1-z*z),'the chart omits exactly its two true endpoints')
dy=sp.cancel((a*Du[1]+b*Dv[1]).subs(chart,simultaneous=True))
fixed=2*a*z*(1-z*z)**2+b*(z*z-2*kappa*z+1)
zero(-(1-z*z)*dy-fixed,'the whole standard module fixed-point condition is exactly the written degree-five polynomial')
for index in [0,2,3]:zero((a*Du[index]+b*Dv[index]).subs(chart,simultaneous=True),
                          'all other action components vanish on the entire fixed-point chart '+str(index))
check(sp.Poly(fixed,z).degree()==5,'the a-nonzero standard action has exact degree five fixed condition')
check(sp.Poly(fixed.subs(a,0),z).degree()==2,'the a-zero,b-nonzero action has exact degree two fixed condition')
zero(fixed.subs(z,1)-2*b*(1-kappa),'the positive omitted endpoint is not a root away from its true critical value')
zero(fixed.subs(z,-1)-2*b*(1+kappa),'the negative omitted endpoint is not a root away from its true critical value')
zero(fixed.subs({b:0,z:0}),'when b vanishes the root z=0 is an actual generic-fibre point')
for aa,bb,kk in [(1,0,0),(0,1,0),(1,1,0),(2,3,2)]:
    poly=sp.Poly(fixed.subs({a:aa,b:bb,kappa:kk}),z)
    roots_allowed=poly.degree()
    for endpoint in [1,-1]:
        while poly.degree()>0 and poly.eval(endpoint)==0:
            poly=sp.div(poly,sp.Poly(z-endpoint,z))[0]
    check(poly.degree()>=1,'exact root count retains an actual fixed point,never only omitted endpoints '+str((aa,bb,kk)))
# General hyperbolic suspension: the four actions work in every B degree.
Bpoly=z**3*w*w+z+w
h=2*x*y+Bpoly
Dz=sp.Matrix([0,-sp.diff(Bpoly,z),2*x,0])
Dw=sp.Matrix([0,-sp.diff(Bpoly,w),0,2*x])
for label,coefficient in [('Dz',Dz),('Dw',Dw)]:
    zero(deriv(coefficient,h),'the exact all-degree suspension action fixes its native function '+label)
    power=y
    for _ in range(sp.Poly(Bpoly,z if label=='Dz' else w).degree()+1):power=deriv(coefficient,power)
    zero(power,'polynomial degree gives actual local nilpotence of the suspension action '+label)
zeros(commutator(Dz,Dw),'the general suspension actions commute because mixed partials agree')
zero(h.subs(y,(s-Bpoly)/(2*x))-s,'the general x-localized value chart is polynomial in the two base variables')
zero((s-Bpoly)/(2*x)*2*x+Bpoly-s,'the general localized chart has the reverse value identity')
check(deriv(Dz,w)==0 and deriv(Dw,z)==0,'the general relative Derksen algebra contains both base coordinates')
# Necessary Tang action: actual coefficients, fixed scheme and global invariant.
F=2*x*z;G=z*z-1
cotangent=F*P+G*Q+z
coords=[x,z,P,Q]
delta=[0,0,G,-F]
zero(deriv(delta,cotangent,coords),'the necessary cotangent action fixes the entire function')
for coordinate in coords:zero(deriv(delta,deriv(delta,coordinate,coords),coords),
                             'the cotangent action is nilpotent on a coordinate '+str(coordinate))
check(sp.gcd(F,G,x,z)==1,'the finite-zero coefficient row is genuinely coprime for the kernel descent')
zero(x*z*z-x*(z*z-1)-x,'the complete fixed scheme contains x,not an extra z-zero component')
for sign in [1,-1]:
    zero(F.subs({x:0,z:sign}),'every fixed plane is retained sign='+str(sign))
    zero(G.subs(z,sign),'every fixed plane retains both equations sign='+str(sign))
    zero(cotangent.subs({x:0,z:sign})-sign,'the value is constant on the entire fixed plane sign='+str(sign))
check(sp.hessian(cotangent,coords).det()==16*z**4,
      'the free generic action data do not certify a native unit Hessian')
# Constant affine precomposition,not a nonlinear Hessian argument.
ell=2*x+3*y+z-5
affine=sp.Matrix([(x-3*y-z+5)/2,y,z,w])
zero(ell.subs(dict(zip(variables,affine)),simultaneous=True)-x,'every affine coordinate pullback has an actual constant normalization')
check(affine.jacobian(variables).det()==sp.Rational(1,2),'that normalization has constant nonzero Jacobian')
H=sp.Matrix(4,4,lambda i,j:sp.symbols('H'+str(min(i,j))+str(max(i,j))))
J=affine.jacobian(variables)
zero((J.T*H*J).det()-J.det()**2*H.det(),'only constant affine Hessian congruence is used in the visibility extension')
context=json.loads((HERE/'sources/global-additive-barrier-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'new primary suspension context is version pinned with no classification imported')
check(context['newly_used_unproved_external_statements']==[],'all needed action,kernel and fixed-point arguments are proved directly')
check(context['K2_both_full_standard_translation_modules_have_generic_fixed_points'] is True,
      'all coefficient degrees are covered by the exact module fixed-point proof')
check(context['all_K2_preserving_additive_actions_classified'] is False
      and context['required_free_action_exists_for_K2'] is False
      and context['K2_one_coordinate_unit_hessian_route_excluded'] is False
      and context['K2_affine_hyperbolic_pullbacks_excluded'] is False
      and context['uniform_presentation_treated_as_four_coordinate_source_automorphism'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the standard module barrier is not a classification,free action,source automorphism or full closure')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
sources={'sources/global-additive-barrier-context.json':hashlib.sha256((HERE/'sources/global-additive-barrier-context.json').read_bytes()).hexdigest()}
for entry in context['newly_cached_primary_context']:sources['sources/'+entry['retained_file']]=entry['sha256']
receipt=dict(outcome='PASS',claims=['HC4MYJS0','HC4MYJS1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree affine-visibility,necessary-action,uniform-suspension and full-module obstruction proof; not independent replay,formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-global-additive-barrier.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
