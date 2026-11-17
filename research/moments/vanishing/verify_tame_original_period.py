"""Selected exact controls; the universal proof is the companion note."""
from pathlib import Path
import hashlib,json
import sympy as s

t,z,a,b,q=s.symbols('t z a b q')
rows=[]
for j in range(8):
    # delta_c=C[u]: partial_s is multiplication by u; s-c=-d/du.
    u=s.symbols('u')
    v=u**j
    assert s.expand(-s.diff(u*v,u)+(j+1)*v)==0
    rows.append({'derivative_order':j,'Euler_index':-j-1})
# A genuine CRT projector, acting on two polynomial-exponential summands.
projector=s.invert(q**2,(q-1)**3)*q**2
assert s.rem(projector,q**2,q)==0
assert s.rem(projector-1,(q-1)**3,q)==0
f=(t**2+2*t+3)*s.exp(t)+(t+4)
poly=s.Poly(s.expand(projector),q)
projected=sum(coeff*s.diff(f,t,k[0]) for k,coeff in poly.terms())
assert s.simplify(projected-(t**2+2*t+3)*s.exp(t))==0
# Both boundary charts of the proper quadric family are product charts.
g=1-z**2
eq=a*b-g
assert s.expand(eq.subs(b,g/a))==0
boundary_eq=b-g*a
assert s.diff(boundary_eq,b)==1
assert s.expand(boundary_eq.subs(b,g*a))==0
# Compact sphere moments of P=z-1, with z uniform on [-1,1].
moments=[]
for m in range(10):
    moment=s.integrate((z-1)**m,(z,-1,1))/2
    assert moment==(-2)**m/s.Integer(m+1)
    moments.append(str(moment))
F=(1-s.exp(-2*t))/(2*t)
assert s.simplify((2*t+1)*F+t*s.diff(F,t)-1)==0
assert s.limit(F,t,0)==1
assert s.limit(s.diff(F,t),t,0)==-1
out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'point_module_Euler_indices':rows,'CRT_original_marker_projection_exact':True,
 'quadric_boundary_product_chart_identities_exact':True,
 'quadric_original_height_moments':moments,
 'quadric_t_dependent_amplitude_has_constant_period':True,
 'quadric_original_unit_pure':False,
 'universal_tame_lattice_and_strictness_proofs_machine_verified':False,
 'wildness_sufficient_for_purity_or_failure':False,
 'exceptional_counterexample_constructed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('tame_original_period_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
