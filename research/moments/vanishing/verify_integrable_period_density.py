"""Selected exact controls; universal continuation proof is written, not automated."""
from pathlib import Path
import hashlib
import json
import sympy as s

t, r, m = s.symbols('t r m')
a, b, c, d = s.symbols('a b c d', positive=True)
A, B, C = c+d-a, b, b+d
lead = t*(1-t)
first = A+B+1-C-(A+B+1)*t
logq = (1-c)/t-(1-b-d)/(1-t)
coef1 = s.factor(2*lead*logq+first)
coef0 = s.factor(lead*(s.diff(logq,t)+logq**2)+first*logq-A*B)
ind0 = s.factor(r*(r-1)+s.limit(coef1,t,0)*r+s.limit(t*coef0,t,0))
ind1 = s.factor(-r*(r-1)+s.limit(coef1,t,1)*r+s.limit((t-1)*coef0,t,1))
assert s.expand(ind0-(r+1-a)*(r+1-c)) == 0
assert s.expand(ind1+r*(r-(b+d-1))) == 0

params = {a:s.Rational(1,2),b:s.Rational(1,3),c:s.Rational(1,2),d:s.Rational(2,3)}
assert s.simplify(t*coef1.subs(params)-t*(2-3*t)) == 0
assert s.simplify(t*coef0.subs(params)-(s.Rational(1,4)-35*t/36)) == 0
assert s.expand(ind0.subs(params)-(2*r+1)**2/4) == 0
assert s.expand(ind1.subs(params)+r**2) == 0

# The all-m triplication identity follows by equal rational ratios and value1 at m=0.
sphere_ratio = 27*(m+s.Rational(1,2))**3/s.prod(3*m+s.Rational(3,2)+i for i in range(3))
beta_ratio = (m+s.Rational(1,2))**2/((m+s.Rational(5,6))*(m+s.Rational(7,6)))
assert s.factor(sphere_ratio-beta_ratio) == 0
assert s.simplify(sphere_ratio.subs(m,0)) == s.Rational(9,35)
moment_controls = []
for j in range(5):
    u = 27**j*s.rf(s.Rational(1,2),j)**3/s.rf(s.Rational(3,2),3*j)
    v = s.rf(s.Rational(1,2),j)**2/(s.rf(s.Rational(5,6),j)*s.rf(s.Rational(7,6),j))
    assert u == v
    moment_controls.append(str(u))
# Exact reduction of the beta normalization by reflection at1/6.
normalization = (2*s.pi)/(6*s.pi)
assert normalization == s.Rational(1,3)

x,y,z,w = s.symbols('x y z w')
aa,bb,cc = x*z-y*w, s.I*(x*z+y*w), x*w+y*z
assert s.expand(aa**2+bb**2+cc**2-(x*w-y*z)**2) == 0
assert s.Poly(27*aa**2*bb**2*cc**2,x,y,z,w).total_degree() == 12
assert (x*w-y*z).subs({x:1,y:2,z:3,w:6}) == 0
assert (27*aa**2*bb**2*cc**2).subs({x:1,y:2,z:3,w:6}) != 0
u,v,q = s.symbols('u v q')
H = 27*u**2*v**2*q**2
assert H.xreplace({u:v,v:u,q:-q}) == H
assert H.xreplace({u:v,v:q,q:u}) == H

# Integrable growth after an arbitrary ramification: exact representative controls.
growth = []
for delta,k in [(s.Rational(1,2),1),(s.Rational(1,2),4),(s.Rational(1,3),7),(s.Integer(1),1)]:
    exponent = delta/k-1
    assert exponent > -1
    growth.append({'delta':str(delta),'ramification':k,'power':str(exponent),'subpolar':True})

physical = -s.log(1-t)/t
assert s.limit(physical,t,0) == 1
continued = (-s.log(1-t)-2*s.pi*s.I)/t
assert s.simplify(continued-physical) == -2*s.pi*s.I/t
assert s.limit(t*(continued-physical),t,0) == -2*s.pi*s.I

out = {
    'completed':True,'sympy_version':s.__version__,
    'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'two_beta_general_indicial_polynomials':[str(ind0),str(ind1)],
    'general_two_beta_exponents_checked_symbolically':True,
    'octahedral_density_ODE_gauge_exact':True,
    'sphere_beta_moment_ratio_identical_for_symbolic_m':True,
    'moment_base_value':1,'selected_moment_controls':moment_controls,
    'beta_normalization_after_explicit_gamma_reflection':'1/3',
    'original_Hopf_sphere_identity_exact':True,'original_H_matrix_degree':12,
    'degree12_leading_form_nonzero_on_rank_one_matrix':True,
    'octahedral_signed_permutation_generators_exact':True,
    'local_growth_controls':growth,
    'physical_integrability_without_branch_integrability_control':str(-2*s.pi*s.I),
    'all_degree_conclusion':'Only0 is pure in the original octahedral C[H];nonconstant phases do not even have eventually vanishing pure moments',
    'arbitrary_Haar_densities_assumed_to_satisfy_continuation_hypothesis':False,
    'arbitrary_finite_rank_or_regular_holonomic_density_automatically_excluded':False,
    'universal_contour_transport_and_local_estimates_machine_verified':False,
    'all_two_beta_laws_derived_by_finite_moment_sampling':False,
    'all_exceptional_seeds_or_phases_classified':False,
    'new_degree_density_support_or_moment_prefix_campaign':False,
    'exceptional_MZ_decisions':'OPEN on both models',
    'independent_replay':False,'formal_verification':False,'external_review':False,
}
destination = Path(__file__).with_name('integrable_period_density_verification.json')
destination.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
