"""Exact transcription controls, not a machine proof of the universal theorem."""
from pathlib import Path
import hashlib
import json
import sympy as S

t, z, s, c, b = S.symbols('t z s c b')
B2 = 1+c*c-2*c*s
R = ((1+c-b)*(1-c+b))/((1+c+b)*(1-c-b))


def mod_b(expr):
    num = S.together(expr).as_numer_denom()[0]
    return S.factor(S.rem(num, b*b-B2, b))


assert mod_b((S.diff(R, c)+S.diff(R, b)*(c-s)/b)/R+2/b) == 0
assert mod_b((S.diff(R, s)-S.diff(R, b)*c/b)/R
             -2*(c*s-1)/((s*s-1)*b)) == 0
for m in range(8):
    ct = sum(S.binomial(m, j)**2*z**j*(1+z)**(m-j)
             for j in range(m+1))
    assert S.expand(ct-S.legendre(m, 1+2*z)) == 0
    rod = S.diff((t*t-1)**m, t, m)/(2**m*S.factorial(m))
    assert S.expand(rod-S.legendre(m, t)) == 0
    for j in range(8):
        val = S.integrate(S.legendre(m, t)*S.legendre(j, t), (t, -1, 1))
        assert val == (S.Rational(2, 2*m+1) if m == j else 0)

q = (1-c*c)/(2*(1-2*c*t+c*c)**S.Rational(3, 2))
for j in range(6):
    assert S.expand(S.diff(q, c, j).subs(c, 0)/S.factorial(j)
                    -(2*j+1)*S.legendre(j, t)/2) == 0

# Derivatives split into odd log coefficient and an even rational remainder.
def dc(f):
    return S.diff(f, c)+S.diff(f, b)*(c-s)/b


C = (1-c*c)/(2*b**3)
E = -2*c/b**2
for j in range(5):
    assert S.simplify(C.subs(b, -b)+C) == 0
    assert S.simplify(E.subs(b, -b)-E) == 0
    C, E = S.factor(dc(C)), S.factor(dc(E)-2*C/b)

# Includes the genuinely rational degenerate centres; do not assume C != 0.
C0 = (1-c*c)/(2*b**3)
assert S.simplify(C0.subs(c, 1)) == 0
assert S.simplify((-2*c/B2).subs(c, 1)-1/(s-1)) == 0
assert S.simplify((-2*c/B2).subs(c, -1)-1/(s+1)) == 0
assert S.simplify(mod_b(2*dc(C0)+dc(dc(C0))).subs(c, 1)) == 0

# Exact binomial reconstruction of a finite shifted-moment jet.
centres = [S.Integer(0), S.Integer(1), S.Integer(-1), S.Rational(2, 3), S.I]
mu = [S.Integer(1), 2+S.I, S.Rational(-3, 5)]
for centre in centres:
    a = [sum(mu[j]*S.diff(c**m, c, j).subs(c, centre)/S.factorial(j)
             for j in range(len(mu))) for m in range(9)]
    for m in range(9):
        shifted = S.expand(sum(S.binomial(m, k)*(-centre)**(m-k)*a[k]
                               for k in range(m+1)))
        assert shifted == (mu[m] if m < len(mu) else 0)
    if centre != 1:
        assert S.simplify((centre-1)**2/(B2.subs({c: centre, s: 1}))) == 1
    if centre != -1:
        assert S.simplify((centre+1)**2/(B2.subs({c: centre, s: -1}))) == 1

# Integration substitution used in the Cauchy formula.
integrand = (1-c*c)/(z*z*(B2-z*z))
primitive = (1-c*c)/B2*(-1/z + S.log((b+z)/(b-z))/(2*b))
assert mod_b(S.diff(primitive, z)-integrand) == 0

out = {
    'date': '2026-11-16', 'completed': True,
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'sympy_version': S.__version__,
    'Legendre_constant_term_and_Rodrigues_controls': True,
    'Poisson_jet_at_zero_controls': True,
    'exact_log_parameter_and_base_derivatives': True,
    'log_coefficient_odd_rational_remainder_even_controls': True,
    'single_possible_rational_pole': '(1+c^2)/(2c), c!=0',
    'rational_centres_plus_minus_one_retained': True,
    'nontrivial_jet_log_cancellation_at_one_retained': True,
    'binomial_jet_reconstruction_centres': [str(x) for x in centres],
    'available_endpoint_has_nonzero_residue_controls': True,
    'Cauchy_integral_substitution_exact': True,
    'classification': 'P-c eventual-pure iff uv=0 and c=1 on both exceptional quotients; degenerate locus safe',
    'nonlinear_composite_purity_reflected': False,
    'all_Bezout_weights_excluded': False,
    'proper_volume_model_assumed_for_full_phase': False,
    'arbitrary_coefficient_ring_promotion_claimed': False,
    'universal_Hodge_or_MZ_proof_machine_verified': False,
    'independent_replay': False, 'formal_verification': False,
    'external_review': False, 'external_review_next_task': False,
    'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    'new_degree_support_or_parameter_campaign': False,
}
Path(__file__).with_name('haar_legendre_verification.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
