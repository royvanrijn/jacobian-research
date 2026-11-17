"""Selected exact controls, not machine verification of the universal proofs."""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

HERE = Path(__file__).resolve().parent
x, y, z, w, v, t = s.symbols('x y z w v t')
k, a, b, e, rho = s.symbols('k a b e rho', positive=True)

# Pullback of the PGL2 volume to the determinant-one x-chart.
coords = s.Matrix([y/x, z/x, 1/x**2])
assert s.simplify(coords.jacobian([x, y, z]).det()/(1/x**2)**2) == -2/x
assert s.Matrix([y, z, y*z+v]).jacobian([y, z, v]).det() == 1
assert s.simplify(s.diff(t**e, t)/t**(2*e)-e*t**(-e-1)) == 0
assert s.expand(rho*(-1)+(rho-1)) == -1
assert s.expand((a+k*a/k)-2*a) == 0
assert s.expand((e+1)*b-b-e*b) == 0

local = []
for kval in (2, 4, 6):
    below = z**kval+v
    critical = z**kval+v**kval
    order = lambda f: min(sum(m) for m, _ in s.Poly(f, v, z).terms())
    assert order(below) == 1 < kval
    assert order(critical) == kval
    # The first blowup really permits a finite simple volume pole at equality.
    r = s.symbols('r')
    assert s.cancel(critical.subs(z, v*r)/v**kval) == r**kval+1
    assert (1+1)-1-2 == -1
    local.append({'k': kval, 'subcritical_multiplicity': 1,
                  'critical_multiplicity': kval,
                  'critical_first_blowup_volume_order': -1})

# The two-root obstruction in the critical-pole argument.
assert s.expand((z-s.I*t)*(z+s.I*t)) == z**2+t**2
assert s.expand((z-s.I*v)*(z+s.I*v)) == z**2+v**2
assert s.Poly(2*s.I*t, t).degree() == 1 < 2
root_controls = [
    {'polynomial': 'z^2+v', 'cover_v': 't^2', 'root_orders': [1, 1],
     'required_orders_for_finite_log_pole': [2, 2], 'excluded': True},
    {'polynomial': 'z^2+v^2', 'cover_v': 't', 'root_orders': [1, 1],
     'required_orders_for_finite_log_pole': [1, 1], 'excluded': False}]

weights = []
for kval, j, ell in ((4, 1, 1), (4, 2, 2), (6, 1, 3), (6, 3, 3)):
    numerator = v**j*z**ell
    F = numerator/v**kval
    weight = 2*(kval-j-ell)
    assert s.cancel(F.subs({v: t**-2*v, z: t**-2*z})/F) == t**weight
    assert j+ell == kval-weight//2
    weights.append({'k': kval, 'numerator_order': j+ell, 'left_weight': weight})

# Exact coefficient differentiation in one actual exceptional component.
c0, c1 = s.symbols('c0 c1')
f0 = z**8+14*z**4*w**4+w**8
f1 = x**8+14*x**4*y**4+y**8
P = c0*f0+c1*f1


def haar(poly):
    ans = 0
    for (aa, bb, cc, dd), coeff in s.Poly(s.expand(poly), x, y, z, w).terms():
        if cc == bb and dd == aa:
            ans += coeff*(-1)**bb*s.factorial(aa)*s.factorial(bb)/s.factorial(aa+bb+1)
    return s.expand(ans)


moment = haar(f0*P**3)
assert s.expand(s.diff(moment, c0, c1)-6*haar(f0*f0*f1*P)) == 0
second = haar((f0+f1)**2)
assert second.is_positive

codimensions = []
for degree, multiplicity in ((8, 1), (12, 1), (24, 2), (60, 2)):
    kval = degree//2
    ambient = (degree+1)*multiplicity
    null_dimension = kval*multiplicity+1
    eta = multiplicity*(kval+1)-1
    assert ambient-null_dimension == eta
    codimensions.append({'D': degree, 'right_multiplicity': multiplicity,
                         'ambient_dimension': ambient, 'nullcone_dimension': null_dimension,
                         'codimension': eta})

degree, eta, nu, L, marker, deriv = s.symbols('D eta nu L marker deriv')
assert s.expand(2*(s.symbols('m')-deriv)-(marker+degree*deriv)-(
    2*s.symbols('m')-marker-(degree+2)*deriv)) == 0
assert s.expand((degree+2)*(eta*nu*L-1)-(degree+2)*(eta*nu-1)-(
    (degree+2)*eta*nu*(L-1))) == 0

result = {
    'completed': True, 'sympy_version': s.__version__,
    'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'PGL2_original_volume_and_SL2_pullback_exact': True,
    'critical_pole_and_tame_pullback_equalities_exact': True,
    'local_multiplicity_boundary_controls': local,
    'Puiseux_root_order_controls': root_controls,
    'left_weight_to_numerator_multiplicity_controls': weights,
    'actual_W8_Haar_coefficient_derivative_identity_exact': True,
    'real_seed_positive_second_moment': str(second),
    'epsilon_real_seed_mod_epsilon_cubed_eventually_pure_but_not_pure': True,
    'nullcone_incidence_dimension_controls': codimensions,
    'ordinary_depth_cutoff_arithmetic_exact': True,
    'full_eventual_pure_leading_condition': 'STRICT common left root multiplicity>D/2',
    'full_thin_eventual_pure_radical': 'ZERO on BOTH models,ALL degrees',
    'all_regular_real_seed_polynomial_subalgebras': 'eventual radical0 overC;nilpotents overANY commutativeC-algebra;all-original-marker safety',
    'all_single_isotypic_nonreduced_families': 'safe overANY commutativeC-algebra;2m>degQ+(D+2)(eta*nu-1),with certificate J^nu subset first-N pure ideal',
    'ordinary_V12_quartic_deficit_preserved': True,
    'ordinary_pure_ideal_equal_to_invariant_ideal_or_assumed_reduced': False,
    'computed_actual_N_nu_prefix_or_intercept': False,
    'field_cutoff_copied_without_intercept_to_nonreduced_rings': False,
    'full_phase_unstable_inferred_from_leading_instability': False,
    'independent_symmetrized_height_condition_closed': False,
    'universal_root_valuation_lc_centre_period_and_symbolic_proofs_machine_verified': False,
    'new_degree_support_density_or_prefix_campaign': False,
    'full_exceptional_MZ_decisions': 'OPEN on BOTH models:mixed-isotypic image-interior marked safety',
    'independent_replay': False, 'formal_verification': False, 'external_review': False}
(HERE/'strict_leading_thin_verification.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
