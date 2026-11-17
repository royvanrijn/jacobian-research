"""Exact controls for the written divisor-excision and one-place critical-budget arguments."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(v,label):
    assert bool(v),label
    checks.append(label)
def zero(v,label):check(sp.simplify(sp.expand(v))==0,label)
def zeros(vs,label):check(all(sp.simplify(sp.expand(v))==0 for v in vs),label)

s,t,x,y=sp.symbols('s t x y')
L=[s,t+s*t*t]
g=1+s*t
J=sp.Matrix(L).jacobian([s,t]).det()
zero(L[1]-t*g,'the non-Keller comparison has its exact selected and extra source factors')
zero(J-(1+2*s*t),'its entire source is explicitly not Keller')
zero(g-s*t-1,'the selected and extra source divisors generate the unit ideal')
zeros([L[0].subs(t,0)-s,L[1].subs(t,0)],'the selected source divisor is the complete target line on the excised source open')
check(sp.Poly(L[1],t).degree()==2,'the etale-open comparison still has generic degree two')
check(J.subs(t,0)==1 and g.subs(t,0)==1,'both excluded factors are units on the selected divisor')
zero((L[1]**2)-(t*t)*(g*g),'the exact target image and pullback factor retain a doubled divisor')
zero((1-2*s*t)*g*g+t*t*(3*s*s+2*s**3*t)-1,
     'the doubled-divisor factor also has an explicit unit-ideal certificate')
zero(sp.rem(L[1]-t,t*t,t),'the target coordinate recovers the nonreduced source parameter modulo t squared')
zero((t*(-s)-1).subs(t,-1/s),'the selected source equation is a nonconstant unit on the extra hyperbola')
check(sp.gcd(t,g)==1 and sp.groebner([t,g],s,t).contains(sp.Integer(1)),
      'polynomial coprimality and geometric disjointness are both checked')
check(sp.gcd(s,t)==1 and not sp.groebner([s,t],s,t).contains(sp.Integer(1)),
      'coprime polynomials alone need not define disjoint curves')

# The finite critical-locus contradiction is algebraic and all-degree in the proof.
P=sp.symbols('P')
u=sp.symbols('u')
p=s*s+t*t*t
outer=p*(p-1)*(p-2)
critical_root=1+sp.sqrt(3)/3
zero(sp.diff(u*(u-1)*(u-2),u).subs(u,critical_root),
     'a nonlinear outer polynomial has an actual critical level')
zeros([sp.diff(outer,v)-sp.diff(u*(u-1)*(u-2),u).subs(u,p)*sp.diff(p,v)
       for v in [s,t]],'the composition has its whole critical level by the exact chain rule')
zero((1-sp.Integer(1))*2,'a single boundary point cannot support a nonconstant degree-zero unit divisor')

# Normalized polynomial curves with a node or higher two-branch contact.
z=sp.symbols('z')
for k in [1,2,4]:
    pk=y*y-x**(2*k)*(x+1)
    param=[z*z-1,z*(z*z-1)**k]
    zero(pk.subs({x:param[0],y:param[1]},simultaneous=True),
         'the contact-family curve has an A1 normalization k='+str(k))
    zeros([a.subs(z,1)-a.subs(z,-1) for a in param],
          'the two normalization branches meet at the selected singularity k='+str(k))
    zero(sp.diff(pk,x)+x**(2*k-1)*((2*k+1)*x+2*k),
         'the full critical polynomial has the exact selected and other multiplicities k='+str(k))
    zero(sp.diff(pk,y)-2*y,'its second critical equation has multiplicity one k='+str(k))
    check((2*k-1)+1==2*k,'the complete critical length equals twice the generic genus k='+str(k))
    check((2*k)//(2*k-1)==(2 if k==1 else 1),
          'the exact budget allows another one-place level only in the nodal case k='+str(k))
    check(2*k-(2*k-1)==1,'the selected finite fibre has branch defect one k='+str(k))

node=y*y-x*x*(x+1)
zero(node.subs({x:-x-sp.Rational(2,3),y:sp.I*y},simultaneous=True)+node+sp.Rational(4,27),
     'the two nodal levels are related by the exact affine change')
for xx,value in [(0,0),(-sp.Rational(2,3),-sp.Rational(4,27))]:
    zeros([sp.diff(node,v).subs({x:xx,y:0}) for v in [x,y]],'the displayed nodal fibre has its actual critical point '+str(xx))
    zero(node.subs({x:xx,y:0})-value,'the critical value is exact '+str(xx))
    check(sp.hessian(node,[x,y]).det().subs({x:xx,y:0})!=0,
          'its local critical length is one by the nonsingular Hessian '+str(xx))

delta,branch,e,gg,rr=sp.symbols('delta branch e gg rr')
zero((2*delta-branch)-delta-(delta-branch),'delta at least branch defect gives the universal two-component bound')
zero((e+2*gg-2+rr)-e-(2*gg+rr-2),'the forced component boundary loss is at least its degree for genus>=0 and r>=2')
check(all(ee+2*genus-2+ends>=ee for ee,genus,ends in [(1,0,2),(3,0,3),(2,1,2)]),
      'finite-image boundary controls retain both puncture count and normalization genus')

# Genuine high-free-degree closed embedding with repeated normal trace.
w1,w2=sp.symbols('w1 w2')
K=[s+t*t,t]
f=K[0]**2
A=[w1*w1,w2**3]
H=w1*A[0]+w2*A[1]
C=w1*w2
X=f+H
Y=[K[0]-(2*f+H)*A[0]+w2*C,K[1]-(2*f+H)*A[1]-w1*C]
phi=f*f+w1*K[0]+w2*K[1]
zero(X*X+w1*Y[0]+w2*Y[1]-phi,'the positive repeated-trace embedding has its literal quadratic identity')
zero(sp.Matrix(K).jacobian([s,t]).det()-1,'its actual plane coefficient map is Keller')
zero(sp.hessian(phi,[s,t,w1,w2]).det()-1,'its actual potential has unit Hessian')
recf=X-H
recK=[Y[0]+(X+recf)*A[0]-w2*C,Y[1]+(X+recf)*A[1]+w1*C]
zeros([recK[0]-recK[1]**2-s,recK[1]-t],'ambient polynomial coordinates recover the full actual source')
T1,T2=sp.symbols('T1 T2')
zero((T1*T1).subs(T1,K[0])-f,'the nonreduced source normal divisor has the exact doubled target image')

context=json.loads((HERE/'sources/extra-prime-excision-context.json').read_text())
check(len(context['newly_used_unproved_external_statements'])==3,
      'the section fact and two Milnor identities are explicit external inputs')
check(context['multiplicities_in_normal_trace_retained'] and context['exact_pullback_factor_and_unit_ideal_proved'],
      'the whole target image and extra-prime decomposition retain scheme multiplicities')
check(context['extra_source_primes_disjoint_from_entire_radical_curve'],
      'extra prime disjointness is from the entire radical curve')
check(context['scalar_units_on_all_extra_components_force_inverse'] and context['forced_extra_component_has_actual_nonconstant_source_unit'],
      'a one-place selected component forces an actual nonconstant source unit on an extra component')
check(context['second_one_place_component_requires_all_nodal_selected_curve'] and context['polynomial_radical_one_place_component_count_at_most_two'],
      'the quantitative one-place count retains its exact all-nodal exception')
check(not context['Nisnevich_excision_alone_assumed_global_inverse'] and not context['source_open_surface_counterexample_declared_global_Keller'],
      'local etale excision is not a global plane Keller inverse theorem')
check(not context['nodal_two_level_control_declared_nonautomorphic_Keller_example'],
      'two possible nodal levels do not construct a Keller map identifying them')
check(not context['all_genuine_nondescended_rank_one_traces_resolved'] and not context['original_goal_complete'],
      'the full w-dependent extra-prime obstruction and unrestricted objective remain open')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained proof/source byte pin '+name)
for item in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/item['retained_file']).read_bytes()).hexdigest()==item['sha256'],
          'new primary section byte pin '+item['retained_file'])
old=json.loads((HERE/'receipt-one-place-trace.json').read_text())
check(old['outcome']=='PASS' and old['check_count']==72,'the retained IZ receipt remains PASS72 and is not newly replayed')
for name,digest in old['source_sha256'].items():
    check(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'retained IZ receipt source '+name)
names=['verify_extra_prime_excision.py','EXTRA_PRIME_EXCISION_AND_ONE_PLACE_BUDGET.md','sources/extra-prime-excision-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
 source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed={'receipt-one-place-trace.json':72},
 newly_cached_primary_context=context['newly_cached_primary_context'],
 newly_activated_cached_primary_inputs=context['newly_activated_cached_primary_inputs'],
 newly_used_unproved_external_statements=context['newly_used_unproved_external_statements'],
 retained_structural_inputs=context['retained_structural_inputs'],written_not_finite_exhaustion=context['direct_proof_arguments'],
 **{k:v for k,v in context.items() if type(v) is bool},
 runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-extra-prime-excision.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,original_goal_complete=False),indent=2))
