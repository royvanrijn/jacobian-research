"""Finite exact transcription controls, not an analytic proof of rigidity."""
from pathlib import Path
import hashlib
import json
import sympy as s
import moment_bridges as b

ROOT = Path(__file__).resolve().parent
checks = {}


def check(group, assertion):
    if not assertion:
        raise AssertionError(group)
    checks[group] = checks.get(group, 0) + 1


def inverse_by_composition(a, z, order):
    """Recursive inversion using composition, independently of Lagrange's sum."""
    ans = z
    coeffs = {j: a.coeff(z, j) for j in range(2, s.degree(a, z)+1)}
    for k in range(2, order+1):
        error = sum(c * s.Poly(ans**j, z).nth(k) for j, c in coeffs.items())
        ans -= s.expand(error) * z**k
    return s.expand(ans)


z, d = s.symbols('z delta')
r = 2-3*z+z*z
for k in range(1, 7):
    F = z+5*z**(k+1)+z**(k+2)*r
    G = (1+d**k)*F
    rhs = 1-k*d**k+(1+d**k)*z**(k+1)*(r+z*s.diff(r, z))+(k+1)*(G-z)/z
    check('perturbation_multiplier', s.expand(s.diff(G, z)-rhs) == 0)

for mu in [(2,), (1, -1), (1, 2, -3), (0, 1, -1), (1, 1, 1, 1)]:
    ys = s.symbols('y:' + str(len(mu)))
    a = s.Poly(z * s.prod(1-c*z for c in mu), z).as_expr()
    inv = inverse_by_composition(a, z, 8)
    f = b.rigidity_phase(mu, ys)
    for k in range(1, 8):
        moment = b.factorial(s.expand(f**k), ys)
        rhs = (k+1)*s.factorial(k)**(len(mu)+1)*inv.coeff(z, k+1)
        check('inverse_factorial', s.expand(moment-rhs) == 0)
        simple = b.dirichlet(s.expand(f**k), ys, [s.Integer(1)]*len(mu))
        check('inverse_simplex', s.expand(moment-s.rf(len(mu), (len(mu)+1)*k)*simple) == 0)

# A sharp r-term gap: roots of unity give a=z(1-z^r).
for roots in [(1, -1), (1, s.I, -1, -s.I)]:
    ys = s.symbols('y:' + str(len(roots)))
    f = b.rigidity_phase(roots, ys)
    moments = [b.factorial(s.expand(f**k), ys) for k in range(1, len(roots)+1)]
    check('sharp_root_unity_gap', all(t == 0 for t in moments[:-1]) and moments[-1] != 0)

# Critical-value collisions: a=((1+z)^D-1)/D has just one critical value.
for degree in (3, 5, 7):
    a = s.expand(((1+z)**degree-1)/degree)
    inv = inverse_by_composition(a, z, 9)
    check('one_critical_value_inverse', all(inv.coeff(z, k) != 0 for k in range(2, 10)))
    check('one_critical_value_binomial', all(inv.coeff(z, k) == s.binomial(s.Rational(1, degree), k)*degree**k for k in range(1, 10)))

y = s.symbols('y:2')
f = y[0]**2-y[0]*y[1]/3+2*y[1]+1
for sizes in ([2, 2], [1, 3]):
    alpha = [s.Rational(k, 2) for k in sizes]
    F, xs, blocks = b.gaussian_lift(f, y, sizes)
    for marker in (y[0], y[1]**2, y[0]*y[1]+1):
        Q, _, _ = b.gaussian_lift(marker, y, sizes)
        for k in range(4):
            check('gamma_gaussian_markers', s.expand(b.gamma(s.expand(marker*f**k), y, alpha)-b.gaussian(s.expand(Q*F**k), xs)) == 0)
    for exponents in ((2, 0, 0, 0), (2, 0, 2, 0), (4, 2, 0, 0), (1, 0, 0, 0)):
        powers, factor = b.radial_monomial(exponents, sizes)
        Q = s.prod(x**a for x, a in zip(xs, exponents))
        for k in range(3):
            rhs = 0 if powers is None else factor*b.gamma(s.expand(s.prod(v**a for v, a in zip(y, powers))*f**k), y, alpha)
            check('arbitrary_gaussian_monomial', s.expand(b.gaussian(s.expand(Q*F**k), xs)-rhs) == 0)

ts = s.symbols('t:3')
ys = s.symbols('y:3')
p = ts[0]+2*ts[1]**2-ts[2]/3
F = b.homogeneous_lift(p, ts, ys, 2)
q = ts[0]*ts[2]+1
H = b.homogeneous_lift(q, ts, ys, 2)
for k in range(4):
    check('inhomogeneous_simplex', s.expand(b.gamma(H*F**k, ys, [1, 1, 1])-s.rf(3, 2+2*k)*b.dirichlet(q*p**k, ts, [1, 1, 1])) == 0)

for prime in (2, 3, 5):
    composed = s.Poly((z-z**prime)+(z-z**prime)**prime-z, z, modulus=prime)
    check('positive_characteristic_boundary', composed == s.Poly(-z**(prime**2), z, modulus=prime))

# Replay the retained sharper Image counterexample by direct balanced contraction.
u, v, x, yvar = s.symbols('xi1 xi2 z1 z2')
R, Z, W, T = u*x+v*yvar, u*yvar, 2*v*x, u*x-v*yvar
image_phase = s.expand((R+Z)*(R**2*W-(2*R+Z)*T**2/2))
def balanced_contraction(poly):
    return sum(c*s.factorial(exps[0])*s.factorial(exps[1])
               for exps, c in s.Poly(poly, u, v, x, yvar).terms()
               if exps[:2] == exps[2:])
for k in range(1, 5):
    power = s.Poly(image_phase, u, v, x, yvar)**k
    check('retained_image_pure', balanced_contraction(power) == 0)
    check('retained_image_marker', balanced_contraction(s.Poly(Z, u, v, x, yvar)*power)
          == s.factorial(4*k+2)*s.factorial(k)/s.factorial2(2*k+1))

# Check the unchanged existing producer's actual schema, not a guessed adapter.
import sys
sys.path.insert(0, str(ROOT.parent / 'moment-vanishing-audit'))
from gaussian_source_hom_producer import normalize
lift, xs, _ = b.gaussian_lift(y[0]-2*y[1], y)
encoded = b.producer_input(lift, xs)
normalize(encoded)
check('existing_producer_schema', encoded['dimension'] == 4)
parameter = s.Symbol('a')
normalize(b.producer_input(parameter*xs[0], xs, {'kind': 'number_field', 'minimal_polynomial': 'a^2+1'}))
check('existing_producer_number_field', True)
for invalid in (s.Float('0.5')*xs[0], s.I*xs[0]):
    try:
        b.producer_input(invalid, xs)
    except ValueError:
        check('reject_inexact_or_unrepresented_field', True)
    else:
        raise AssertionError('invalid producer coefficient accepted')
try:
    b.producer_input(s.I*xs[0], xs, {'kind': 'number_field', 'minimal_polynomial': 'a^2+1'})
except ValueError:
    check('reject_inexact_or_unrepresented_field', True)
else:
    raise AssertionError('an unrepresented algebraic constant was accepted')

receipt = {
    'date': '2026-11-16',
    'status': 'pass',
    'checks': checks,
    'total': sum(checks.values()),
    'scope': 'Finite exact bridge, inverse-index, perturbation-algebra and adapter checks. Universal rigidity uses the written complex-analytic proof.',
    'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), ROOT/'moment_bridges.py', ROOT/'RIGIDITY_CRITICAL_VALUES.tex']},
    'analytic_theorem_formally_verified': False,
    'independent_replay': False,
    'external_review': False,
}
(ROOT / 'verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
