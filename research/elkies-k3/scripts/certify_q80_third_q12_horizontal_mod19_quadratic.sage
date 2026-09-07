#!/usr/bin/env sage
"""Recover the Q80 third-q12 horizontal over a quadratic p=19 field.

The exact MW5 embedding certificate reduces the target to ``H=P4+Q``, where
``P4`` is the fourth ordered rational polynomial section and ``Q`` is a
polynomial section over the algebraic closure.  This verifier rebuilds the
polynomial-section scheme, takes its radical, proves the complete factorization
of its degree-twelve lexicographic eliminant, reconstructs every quadratic
solution, and filters ``P4+Q`` by the pinned target fingerprint.

The result is an exact finite-field marked horizontal.  It is a seed for
deformation/CRT lifting, not a characteristic-zero section or q12 child.
"""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, prod
from sage.repl.load import load


ROOT = Path(__file__).resolve().parents[2]
EXPORTER = ROOT / "elkies-k3/scripts/export_q80_third_q12_polynomial_closure_scheme.sage"
INPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-height-shell-complete.json"
LATTICE = ROOT / "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-quadratic-horizontal.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Keep the exporter's INPUT/OUTPUT/ROOT constants out of this verifier's
# namespace while retaining its exact scheme objects.
export_namespace = dict(globals())
load(str(EXPORTER), export_namespace)
scheme = export_namespace["scheme"]
residual = export_namespace["residual"]
finite = export_namespace["finite"]
names = export_namespace["names"]
A = export_namespace["A"]
B = export_namespace["B"]
payload = json.loads(INPUT.read_text())
lattice_payload = json.loads(LATTICE.read_text())
modular = payload["parameters"][0]["modular"][0]

ideal = scheme.ideal(residual)
if ideal.dimension() != 0 or ideal.vector_space_dimension() != 16:
    raise ArithmeticError("unexpected polynomial-section scheme dimension")
radical = ideal.radical()
if radical.vector_space_dimension() != 12:
    raise ArithmeticError("unexpected radical degree")

lex_ring = PolynomialRing(finite, names=names, order="lex")
lex_variables = lex_ring.gens()
lex_ideal = lex_ring.ideal([lex_ring(str(value)) for value in radical.gens()])
lex_basis = tuple(lex_ideal.groebner_basis())
if len(lex_basis) != 6:
    raise ArithmeticError("unexpected lexicographic basis length")
univariate = lex_basis[-1]
expected_univariate = (
    lex_variables[-1]**12 - 2*lex_variables[-1]**10
    + 5*lex_variables[-1]**8 + 3*lex_variables[-1]**6
    - lex_variables[-1]**4 + lex_variables[-1]**2 - 7
)
if univariate != expected_univariate:
    raise ArithmeticError("unexpected radical eliminant")
factorization = tuple(univariate.factor())
if [factor.degree() for factor, exponent in factorization] != [1]*6+[2]*3:
    raise ArithmeticError("unexpected eliminant factor degrees")
if any(exponent != 1 for factor, exponent in factorization):
    raise ArithmeticError("radical eliminant is not squarefree")

base_ring = PolynomialRing(finite, "W")
W = base_ring.gen()
A_base = base_ring(A)
B_base = base_ring(B)
star_factor = next(
    factor.monic() for factor, exponent in (4*A_base**3+27*B_base**2).factor()
    if int(exponent) == 7
)
star_root = -star_factor[0]/star_factor[1]


def polynomial_point(record, curve, ring, field):
    x_value = ring(record["x_coefficients_low_to_high"])
    y_value = ring(record["y_coefficients_low_to_high"])
    point = curve(field(x_value), field(y_value))
    if point.is_zero():
        raise ArithmeticError("stored polynomial point is zero")
    return point


def polynomial_square_root(poly):
    if not poly:
        return poly.parent().zero()
    factorization = tuple(poly.factor())
    if any(int(exponent) % 2 for factor, exponent in factorization):
        return None
    unit = poly.factor().unit()
    if not unit.is_square():
        return None
    answer = poly.parent()(unit.sqrt())
    for factor, exponent in factorization:
        answer *= factor**(int(exponent)//2)
    return answer if answer**2 == poly else None


def po_from_x(point):
    if point.is_zero():
        return None
    x_value = point[0]
    root = polynomial_square_root(x_value.denominator())
    if root is None:
        return None
    finite_poles = root.degree()
    numerator_degree = x_value.numerator().degree()
    twice = max(2*finite_poles, numerator_degree-4)
    if twice < 0 or twice % 2:
        return None
    return int(twice//2)


def relation_value(relation, variable, extension, alpha):
    remainder = relation-variable
    answer = extension.zero()
    for monomial, coefficient in remainder.dict().items():
        if any(monomial[index] for index in range(len(monomial)-1)):
            raise ArithmeticError("lex relation is not triangular in sat")
        answer += extension(coefficient)*alpha**monomial[len(monomial)-1]
    return -answer


def finite_coefficients(poly):
    return [str(value) for value in poly.list()]


exact_fingerprints = {
    tuple(embedding["ordered_P_dot_O_of_H_minus_eight_polynomial_points"])
    for embedding in lattice_payload["polynomial_rank3_subgroup"]["embeddings"]
}
quadratic_records = []
accepted = []
modulus_ring = PolynomialRing(finite, "s_modulus")
for factor_index, (factor, exponent) in enumerate(factorization, 1):
    if factor.degree() != 2:
        continue
    modulus = modulus_ring(str(factor).replace("sat", "s_modulus"))
    extension = finite.extension(modulus, f"a{factor_index}")
    extension_ring = PolynomialRing(extension, "W")
    extension_W = extension_ring.gen()
    extension_field = extension_ring.fraction_field()
    extension_A = extension_ring([extension(value) for value in A_base.list()])
    extension_B = extension_ring([extension(value) for value in B_base.list()])
    curve = EllipticCurve(
        extension_field,
        [0, 0, 0, extension_field(extension_A), extension_field(extension_B)],
    )
    rational_points = tuple(
        polynomial_point(record, curve, extension_ring, extension_field)
        for record in modular["polynomial_shell"]
    )
    P4 = rational_points[3]
    first_alpha = extension.gen()
    # If s^2+b*s+c is the modulus, the conjugate root is -b-first_alpha.
    conjugate_alpha = -extension(modulus[1])-first_alpha
    for conjugate_index, alpha in enumerate((first_alpha, conjugate_alpha), 1):
        coordinate_values = [
            relation_value(lex_basis[index], lex_variables[index], extension, alpha)
            for index in range(5)
        ]
        l_value, x0_value, x1_value, x2_value, x3_value = coordinate_values
        if l_value*alpha != 1:
            raise ArithmeticError("sat*l=1 failed on a quadratic solution")
        X = (
            l_value**2*extension_W**4+x3_value*extension_W**3
            + x2_value*extension_W**2+x1_value*extension_W+x0_value
        )
        square_value = X**3+extension_A*X+extension_B
        Y = polynomial_square_root(square_value)
        if Y is None or Y.degree() != 6:
            raise ArithmeticError("quadratic polynomial-section square root failed")
        if Y[6] != l_value**3:
            Y = -Y
        if Y[6] != l_value**3 or Y**2 != square_value:
            raise ArithmeticError("quadratic polynomial-section orientation failed")
        Q = curve(extension_field(X), extension_field(Y))
        H = P4+Q
        if H.is_zero() or po_from_x(H) != 2:
            profile_ok = False
            fingerprint = None
        else:
            fingerprint = tuple(po_from_x(H-point) for point in rational_points)
            profile_ok = fingerprint in exact_fingerprints

            # The target is on the identity component at the finite I1* place.
            x_value = H[0]
            if not x_value.denominator()(extension(star_root)):
                profile_ok = False
            else:
                node_ring = PolynomialRing(extension, "x_node")
                x_node = node_ring.gen()
                cubic = (
                    x_node**3+extension_A(extension(star_root))*x_node
                    + extension_B(extension(star_root))
                )
                singular = cubic.gcd(cubic.derivative()).roots(multiplicities=False)
                if len(singular) != 1 or x_value(extension(star_root)) == singular[0]:
                    profile_ok = False
            if (
                H[0].numerator().degree()-H[0].denominator().degree() != 4
                or H[1].numerator().degree()-H[1].denominator().degree() != 6
            ):
                profile_ok = False

        record = {
            "factor": str(factor),
            "extension_modulus": str(extension.modulus()),
            "conjugate_index_one_based": conjugate_index,
            "alpha": str(alpha),
            "Q": {
                "x_coefficients_low_to_high": finite_coefficients(X),
                "y_coefficients_low_to_high": finite_coefficients(Y),
            },
            "H": {
                "P_dot_O": po_from_x(H) if not H.is_zero() else None,
                "x_numerator_coefficients_low_to_high": finite_coefficients(H[0].numerator()),
                "x_denominator_coefficients_low_to_high": finite_coefficients(H[0].denominator()),
                "y_numerator_coefficients_low_to_high": finite_coefficients(H[1].numerator()),
                "y_denominator_coefficients_low_to_high": finite_coefficients(H[1].denominator()),
                "ordered_intersection_fingerprint": list(fingerprint) if fingerprint else None,
            },
            "target_profile": bool(profile_ok),
        }
        quadratic_records.append(record)
        if profile_ok:
            accepted.append(record)

accepted_factors = {record["factor"] for record in accepted}
if len(accepted) != 2 or len(accepted_factors) != 1:
    raise ArithmeticError(
        f"expected one conjugate quadratic target orbit, found {len(accepted)} records"
    )

output = {
    "schema": "elkies-k3.q80-third-q12-quadratic-horizontal-modp.v1",
    "status": "PASS_EXACT_QUADRATIC_THIRD_Q12_HORIZONTAL_MOD19",
    "prime": 19,
    "specialization": "u=-2",
    "scheme": {
        "quotient_dimension": 16,
        "radical_degree": 12,
        "lex_basis": [str(value) for value in lex_basis],
        "eliminant_factorization": [[str(factor), int(exponent)] for factor, exponent in factorization],
    },
    "quadratic_solutions": quadratic_records,
    "accepted_horizontal_orbit": accepted,
    "target": {
        "edge": "D7+D5/MW5 --q12--> A5+A3+3A1/MW6",
        "height": "8",
        "P_dot_O": 2,
        "component_profile": "identity at finite I1* and infinity I3*",
    },
    "inputs": {
        "polynomial_scheme_exporter": {"path": str(EXPORTER.relative_to(ROOT)), "sha256": sha256(EXPORTER)},
        "modular_shell": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
        "height_lattice": {"path": str(LATTICE.relative_to(ROOT)), "sha256": sha256(LATTICE)},
    },
    "claim_boundary": {
        "proved": [
            "complete radical decomposition of the exact-degree-four polynomial-section chart",
            "exact quadratic polynomial section identity",
            "exact group-law construction H=P4+Q",
            "target P.O=2 and ordered intersection fingerprint",
            "identity-component equation profile at both reducible fibres",
        ],
        "not_proved": [
            "the connected D7+D5 resolved q12 quotient",
            "the A5+A3+3A1 child equation",
            "alignment across another prime",
            "characteristic-zero lifting",
        ],
    },
    "reproduce": "sage elkies-k3/scripts/certify_q80_third_q12_horizontal_mod19_quadratic.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True, default=int)+"\n")
print(
    "Q80QUADRATICQ12|prime=19|radical_degree=12|quadratic_solutions=6|"
    "accepted_orbit_size=2|P.O=2|height=8|"
    "status=PASS_EXACT_QUADRATIC_THIRD_Q12_HORIZONTAL_MOD19",
    flush=True,
)
