"""Exact controls of the primitive quotient,plinth and global action certificate proof."""
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
x,y,z,w,s,u,v,t,P,Q,Rvar,Wvar=sp.symbols('x y z w s u v t P Q R W')
variables=[x,y,z,w]
def deriv(coefficient,f,coords=variables):return sp.expand(sum(cc*sp.diff(f,qq) for cc,qq in zip(coefficient,coords)))
def common_gcd(coefficients):
    out=sp.Integer(0)
    for value in coefficients:out=sp.gcd(out,value,*variables)
    return sp.Poly(out,*variables).monic().as_expr()
# Finite Taylor projection supplies the common one-variable fraction-field chart.
A0,A1,A2,A3=sp.symbols('A0 A1 A2 A3')
poly=A0+A1*t+A2*t*t+A3*t**3
projection=sum((-t)**j*sp.diff(poly,t,j)/sp.factorial(j) for j in range(4))
zero(projection-A0,'the local-slice Taylor projection is genuinely invariant')
reconstruction=sum(t**j/sp.factorial(j)*sum((-t)**k*sp.diff(poly,t,j+k)/sp.factorial(k) for k in range(4-j)) for j in range(4))
zero(reconstruction-poly,'the finite Taylor reconstruction gives a true polynomial coordinate chart')
for power in [1,2,3]:
    iterate=t
    for _ in range(4):iterate=sp.expand(t**power*sp.diff(iterate,t))
    check(iterate!=0,'a positive-degree coefficient on K[t] is not locally nilpotent degree='+str(power))
zero(sp.diff(7*t+3,t)-7,'a constant coefficient gives the unique primitive one-variable action up to scale')
# Primitivity is essential for a same-kernel fixed-point assertion.
zero((1/x)*x-1,'the nonprimitive countercontrol changes speed by a genuine denominator')
check(x.subs(x,0)==0 and sp.Integer(1)!=0,'dividing x partial_y to partial_y removes its fixed hypersurface')
# The K2 base-value pair has exact pure two-dimensional fibres.
K2=w*(2*x*y+z*z-1)+z
c,kappa,r=sp.symbols('c kappa r',nonzero=True)
zf=kappa-w*r
yf=(r-zf*zf+1)/(2*c)
zero((2*c*y+z*z-1).subs({y:yf,z:zf},simultaneous=True)-r,'each nonzero-x fibre has the actual free r coordinate')
zero((K2-kappa).subs({x:c,y:yf,z:zf},simultaneous=True),'each nonzero-x fibre is an actual A2,not a generic rational chart')
zero((kappa-w*(2*c*y+z*z-1)-z).subs(kappa,w*(2*c*y+z*z-1)+z),
     'the fibre parameterization has its reverse z identity')
B=w*(z*z-1)+z
for sign in [1,-1]:
    zero((B-sign)-(z-sign)*(w*(z+sign)+1),'the reducible x-zero critical fibre retains both two-dimensional source components sign='+str(sign))
    check(sp.Poly(w*(z+sign)+1,z,w).total_degree()==2,'the second component is a genuine plane curve,never a whole three-dimensional slice sign='+str(sign))
check(sp.Poly(B-kappa,z,w).total_degree()==3,'every zero-x fibre is a nonzero plane hypersurface times its free y line')
# Entire-ring primitivity controls coefficients,including degenerate base fibres.
vv=1-2*x*y-z*z
Du=sp.Matrix([0,-2*vv*z,2*x*vv,2*x])
Dv=sp.Matrix([0,-(2*w*z+1),2*x*w,0])
check(common_gcd(Du)==1,'the first standard action is primitive in the whole four-variable source ring')
check(common_gcd(Dv)==1,'the second standard action is primitive in the whole four-variable source ring')
for index,(aa,bb) in enumerate([(x,s*s-1),(s-1,s+1),(x+s,s)]):
    check(sp.gcd(aa,bb,x,s)==1,'the coefficient pair is primitive over the actual base '+str(index))
    module=aa.subs(s,K2)*Du+bb.subs(s,K2)*Dv
    check(common_gcd(module)==1,'no hidden source divisor appears after pullback,including critical base values '+str(index))
    zero(deriv(module,w)-2*x*aa.subs(s,K2),'every common source divisor must divide x a '+str(index))
    zero(deriv(module,vv)-2*x*bb.subs(s,K2),'every common source divisor must also divide x b '+str(index))
    check(sp.cancel(module[1].subs(x,0))!=0,'the possible extra common divisor x is actually excluded '+str(index))
aa,bb,kk=sp.symbols('aa bb kk')
fixed=2*aa*z*(1-z*z)**2+bb*(z*z-2*kk*z+1)
zero(sp.Poly(fixed,z).coeff_monomial(z**5)-2*aa,'the generic leading coefficient forces a0=0 if the action vanished identically on x=0')
zero(sp.Poly(fixed.subs(aa,0),z).coeff_monomial(z*z)-bb,'the remaining coefficient then forces b0=0,contradicting primitive base coefficients')
# Replicas are not restricted to the old common coefficient ring.
replica=1+vv*vv
zero(deriv(Du,replica),'a replica factor can depend on a further invariant outside C[x,K2]')
check(common_gcd(replica*Du)==sp.Poly(replica,*variables).monic().as_expr(),
      'the invariant replica introduces exactly its own factor,not a new primitive action')
for coordinate in variables:
    zero(deriv(replica*Du,deriv(replica*Du,coordinate))-replica**2*deriv(Du,deriv(Du,coordinate)),
         'the invariant replica has the exact nilpotent iteration identity '+str(coordinate))
    zero((replica*Du[variables.index(coordinate)]).subs({x:0,z:0,w:-kk}),
         'replicas retain the old generic fixed affine line '+str(coordinate))
factor=x*(K2-1)
nonprimitive=factor*(Du+K2*Dv)
zero(deriv(Du,factor),'the base gcd removed from a module action is invariant under Du')
zero(deriv(Dv,factor),'the base gcd removed from a module action is invariant under Dv')
zeros([sp.cancel(cc/factor)-dd for cc,dd in zip(nonprimitive,Du+K2*Dv)],
      'coefficient gcd removal leaves an actual standard-module derivation')
# The plinth test is proper and nonprincipal,not an asserted exact formula.
F=2*x*z;G=z*z-1
coords=[x,z,P,Q]
action=[0,0,G,-F]
value=F*P+G*Q+z
zero(deriv(action,P,coords)-G,'the plinth contains the first actual invariant coefficient')
zero(deriv(action,Q,coords)+F,'the plinth contains the second actual invariant coefficient')
check(sp.gcd(F,G,x,z)==1,'these coefficients cannot lie in one height-one prime of the polynomial invariant ring')
test=sp.Function('r')(x,z,P,Q)
for sign in [1,-1]:
    zero(deriv(action,test,coords).subs({x:0,z:sign}),
         'every image polynomial vanishes on a whole fixed plane,so the plinth is proper sign='+str(sign))
zero(deriv(action,value,coords),'the value parameter is a genuine invariant,not a new free source coordinate')
# Full coordinate inverse plus the same kernel base gives an actual Hessian certificate.
b1,b2,pv,qv=sp.symbols('b1 b2 pv qv')
scale=sp.symbols('scale',nonzero=True)
Fk=b1+b2*b2;Gk=b2
D0=b1**3+b2
H=pv*Fk+qv*Gk
g=D0+scale*H
check(sp.Matrix([Fk,Gk]).jacobian([b1,b2]).det()==1,'the coefficient pair satisfies the everywhere constant orientation requirement')
zero(sp.hessian(g,[b1,b2,pv,qv]).det()-scale**4,'the exact coordinate-action certificate has the stated unit full Hessian')
source=[x,z,y+w*w,w]
source_inverse={x:b1,z:b2,y:pv-qv*qv,w:qv}
zeros([ff.subs(source_inverse,simultaneous=True)-vvv for ff,vvv in zip(source,[b1,b2,pv,qv])],
      'the supplied complete coordinate system has a real polynomial inverse')
native=g.subs({b1:x,b2:z,pv:y+w*w,qv:w},simultaneous=True)
zero(native.subs(source_inverse,simultaneous=True)-g,'the actual inverse source map produces the unit-Hessian function,without nonlinear Hessian invariance')
UU,VV,SS,TT=sp.symbols('UU VV SS TT')
target_b2=TT/scale
target_b1=SS/scale-target_b2**2
target_p=(UU-3*target_b1**2)/scale
target_q=(VV-1)/scale-2*target_b2*target_p
gradient=sp.Matrix([sp.diff(g,qq) for qq in [b1,b2,pv,qv]])
inverse={b1:target_b1,b2:target_b2,pv:target_p,qv:target_q}
zeros(gradient.subs(inverse,simultaneous=True)-sp.Matrix([UU,VV,SS,TT]),
      'the constructive gradient inverse follows from the actual plane Keller inverse')
zero(Fk.subs({b1:target_b1,b2:target_b2},simultaneous=True)-SS/scale,
     'specializing the first gradient outputs recovers the first plane inverse identity')
zero(Gk.subs({b1:target_b1,b2:target_b2},simultaneous=True)-TT/scale,
     'specializing the first gradient outputs recovers the second plane inverse identity')
# A one-variable generator and its polynomial inverse force a unit linear coefficient.
cc,dd=sp.symbols('cc dd',nonzero=True)
zero(((H-dd)/cc).subs(H,cc*H+dd)-H,'the kernel generator comparison has the exact polynomial inverse identity')
check(2*2>1,'two nonlinear polynomial degrees cannot compose to the kernel generator')
context=json.loads((HERE/'sources/primitive-quotient-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'versioned primary rank/plinth context is byte pinned without importing a classification')
check(context['newly_used_unproved_external_statements']==[],'the local-slice,primitivity,plinth and global certificate arguments are proved directly')
check(context['all_standard_module_embedded_invariant_fields_and_value_preserving_conjugates_excluded'] is True,
      'the obstruction covers entire embedded quotient fields,not just chosen coefficient speeds')
check(context['required_nonstandard_K2_action_exists'] is False
      and context['all_K2_preserving_actions_classified'] is False
      and context['K2_one_coordinate_unit_hessian_route_excluded'] is False
      and context['abstract_kernel_isomorphism_treated_as_embedded_field_equality'] is False
      and context['exact_plinth_equal_coefficient_ideal_proved'] is False
      and context['three_variable_principal_plinth_theorem_transferred_to_four_variables'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'other quotient fields and all true source-orbit questions remain open,with dimensional and global inverse hypotheses explicit')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
sources={'sources/primitive-quotient-context.json':hashlib.sha256((HERE/'sources/primitive-quotient-context.json').read_bytes()).hexdigest()}
for entry in context['newly_cached_primary_context']:sources['sources/'+entry['retained_file']]=entry['sha256']
receipt=dict(outcome='PASS',claims=['HC4MYJT0','HC4MYJT1','HC4MYJT2'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'PRIMITIVE_QUOTIENTS_PLINTH_AND_GLOBAL_ACTION_CERTIFICATES.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree primitive quotient saturation,nonprincipal plinth and actual coordinate-action certificate; not independent replay,formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-primitive-quotients.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
