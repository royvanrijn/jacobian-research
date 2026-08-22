#!/usr/bin/env sage -python
"""Construct and screen a component-nef chord level over the old base.

The generic component-nef chord is the pullback of the standard marked chord
through ``tau_-P0``.  For a rational level ``c``, use the translated
coordinates (x',y'); there the equation is the line
``y'=c*(x'-x(S))-y(S)``.  Its fixed intersection at ``-S`` divides out,
leaving a quadratic in x'.  This quadratic is the function field of the
degree-two curve in each old elliptic fibre, so its discriminant is the exact
quadratic-extension squareclass over QQ(T).  Pullback by ``tau_-P0`` gives
the physical equation without expanding a prohibitively large rational
group-law substitution.  The branch degree is then checked: only degree two
can be a rational bisection eligible for the collision protocol.

This is a generic-fibre construction.  It does not identify the resolved
pencil member or its complete branch divisor on a regular model.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
CHORD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-chord.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-branch.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(polynomial(ring, data[denominator]))


def point(field, ring, curve, data):
    return curve(
        rational(field, ring, data, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high"),
        rational(field, ring, data, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high"),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--chord", type=Path, default=CHORD)
parser.add_argument("--level", type=QQ, default=QQ(0))
parser.add_argument(
    "--skip-canonical-squareclass", action="store_false", dest="canonicalize_squareclass",
    help="skip the exact squarefree reduction (diagnostic only; not a collision candidate)",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.set_defaults(canonicalize_squareclass=True)
args = parser.parse_args()
for name in ("child", "chord", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
chord = json.loads(args.chord.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert chord["status"] == "PASS_EXACT_COMPONENT_NEF_OLD_ZERO_CHORD"

base = PolynomialRing(QQ, "T")
T = base.gen()
field = base.fraction_field()
A = field(polynomial(base, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]))
B = field(polynomial(base, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]))
curve = EllipticCurve(field, [0, 0, 0, A, B])
points = chord["points_in_standard_weierstrass_group"]
s = point(field, base, curve, points["translated_marked_point_S"])
sx, sy = s.xy()

line_ring = PolynomialRing(field, "x_prime")
x_prime = line_ring.gen()
# Substituting y'=c*(x'-x(S))-y(S) into the Weierstrass equation gives a
# cubic.  It vanishes at x'=x(S), the fixed base intersection -S, and its
# quotient is the required degree-two moving intersection.
line_cubic = line_ring((args.level * (x_prime - sx) - sy)**2 - (x_prime**3 + A*x_prime + B))
quadratic, remainder = line_cubic.quo_rem(x_prime - sx)
assert remainder == 0
quadratic = quadratic.monic()
assert quadratic.degree() == 2
coefficients = quadratic.list()
discriminant = field(coefficients[1]**2 - 4 * coefficients[0] * coefficients[2])
print(
    "H92Q6CHILDCOMPNEFBRANCH|quadratic_ready|level={}|disc_num_deg={}|disc_den_deg={}".format(
        args.level, base(discriminant.numerator()).degree(), base(discriminant.denominator()).degree()
    ),
    flush=True,
)

numerator = base(discriminant.numerator())
denominator = base(discriminant.denominator())


def odd_squarefree_part(value):
    """Return the product of the odd-multiplicity squarefree factors.

    This deliberately uses gcd-based squarefree decomposition rather than a
    full factorization over QQ: the latter is unnecessary for the squareclass
    and is expensive for the degree-192 numerator.
    """
    decomposition = value.squarefree_decomposition()
    result = base.one()
    records = []
    for factor, exponent in decomposition:
        factor = base(factor).monic()
        if exponent % 2:
            result *= factor
            records.append([str(factor), int(exponent)])
    return result, records


squareclass_data = {"canonicalized": False}
geometric_branch_degree = None
if args.canonicalize_squareclass:
    # At level zero the line is y'=-y(S), so the quadratic discriminant is
    # -3*x(S)^2-4*A.  The x(S) denominator is a square h^2, hence the
    # discriminant denominator is itself a square and no full factorization
    # is needed to obtain the extension squareclass.
    assert args.level == 0
    sx_numerator = base(sx.numerator())
    sx_denominator = base(sx.denominator())
    assert A.denominator() in QQ
    branch_numerator = -3 * sx_numerator**2 - 4 * base(A) * sx_denominator**2
    assert discriminant == field(branch_numerator) / field(sx_denominator**2)
    assert branch_numerator.gcd(sx_denominator) in QQ
    print(
        "H92Q6CHILDCOMPNEFBRANCH|squarefree_gcd|branch_num_deg={}".format(branch_numerator.degree()),
        flush=True,
    )
    repeated_part = branch_numerator.gcd(branch_numerator.derivative())
    assert repeated_part in QQ, "the level-zero branch numerator is not squarefree"
    squareclass = field(branch_numerator.monic())
    geometric_branch_degree = int(branch_numerator.degree())
    squareclass_data = {
        "canonicalized": True,
        "squareclass": str(squareclass),
        "representative": "monic(-3*Nx^2-4*A*Dx^2)",
        "branch_numerator_degree": int(branch_numerator.degree()),
        "denominator_square": "Dx^2",
        "branch_numerator_squarefree": True,
        "geometric_branch_degree": geometric_branch_degree,
        "rational_bisection_eligible": geometric_branch_degree == 2,
    }

payload = {
    "schema": "elkies-k3.h92-q6-child-component-nef-bisection-branch.v1",
    "status": "PASS_EXACT_COMPONENT_NEF_GENERIC_BISECTION_BRANCH",
    "inputs": {
        "child": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "component_nef_chord": {"path": str(args.chord.relative_to(ROOT)), "sha256": digest(args.chord)},
    },
    "level": str(args.level),
    "generic_equation": {
        "description": "m_component_nef=level, written after tau_-P0 in (x_prime,y_prime)",
        "translated_line": "y_prime=level*(x_prime-x(S))-y(S)",
        "fixed_intersection": "-S at x_prime=x(S)",
        "physical_equation": "pull back the translated line through tau_-P0",
        "quadratic_in_x_coefficients_low_to_high": [str(value) for value in coefficients],
    },
    "quadratic_extension": {
        "model": "x satisfies the displayed monic quadratic over QQ(T)",
        "discriminant": str(discriminant),
        "squareclass_data": squareclass_data,
    },
    "boundary": (
        "This is the generic quadratic function field of one level of the transported "
        "component-nef chord. It does not identify the resolved closure, infinity "
        "branch contribution, a global pencil certificate, an extension collision, "
        "or a rank result."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDCOMPNEFBRANCH|level={}|quadratic=1|canonicalized={}|branch_degree={}|status=PASS_EXACT_COMPONENT_NEF_GENERIC_LEVEL_BRANCH".format(
        args.level, args.canonicalize_squareclass, geometric_branch_degree
    ),
    flush=True,
)
