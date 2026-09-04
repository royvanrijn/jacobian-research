#!/usr/bin/env sage-python
"""Certify the Frobenius factor of an R17 singleton or product twist.

This verifier consumes the exact output of the open-source
ToricControlledReduction ``readfile`` driver.  The regular twist model

    d(t)*y^2 = x^3 + A(t)*x + B(t)

has a 36-dimensional primitive toric factor with Hodge vector ``[3,30,3]``.
Eight dimensions are a boundary permutation motive.  If ``D={d=0}`` and
``Z={d=0, x^3+A*x+B=0}``, its characteristic polynomial is

    det(T-p*F | H^0(Z)) / det(T-p*F | H^0(D)).

The quotient is the degree-24 or degree-28 elliptic L-polynomial.  The verifier derives
the boundary factor independently from finite-field factor degrees, checks
the two already-audited power sums, the weight-two functional equation, and
every cyclotomic factor allowed by its degree.  In the product case, absence
of such a factor gives ``rho<=18`` and hence geometric Mordell--Weil rank zero.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import sys

from sage.all import AA, GF, PolynomialRing, QQ, ZZ, cyclotomic_polynomial, euler_phi

from parse_toric_controlled_reduction_output import parse_readfile_output


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
DEFAULT_AUDIT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json"
)
DEFAULT_SINGLETON_AUDIT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-singleton-twist-finite-field-bound-audit-v1.json"
)
EXPORTER = ROOT / "elkies-k3/scripts/export_r17_product_toric_frobenius_input.sage"
SINGLETON_EXPORTER = ROOT / "elkies-k3/scripts/export_r17_singleton_toric_frobenius_input.sage"
PARSER = ROOT / "elkies-k3/scripts/parse_toric_controlled_reduction_output.py"
VERIFIER = Path(__file__).resolve()
RUNNER = ROOT / "elkies-k3/scripts/run_r17_product_toric_frobenius.sh"
SINGLETON_RUNNER = ROOT / "elkies-k3/scripts/run_r17_singleton_toric_frobenius.sh"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pair-key", required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--toric-output", type=Path, required=True)
parser.add_argument("--toric-commit", required=True)
parser.add_argument("--toric-executable", type=Path, required=True)
parser.add_argument("--sage-prefix", type=Path, required=True)
parser.add_argument("--bisections", type=Path, default=DEFAULT_BISECTIONS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
parser.add_argument("--singleton-audit", type=Path, default=DEFAULT_SINGLETON_AUDIT)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if not args.toric_executable.is_file():
    raise FileNotFoundError(args.toric_executable)

prime = int(args.prime)
field = GF(prime)
if prime < 5 or not field.is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")
labels = tuple(args.pair_key.split(":"))
if len(labels) not in (1, 2) or len(set(labels)) != len(labels):
    raise ValueError("--pair-key must contain one label or two distinct labels")
twist_degree = 2 * len(labels)
elliptic_degree = 20 + 2 * twist_degree
boundary_degree = 2 * twist_degree
toric_degree = elliptic_degree + boundary_degree
arithmetic_genus = 2 + twist_degree // 2
expected_hodge_numbers = [
    arithmetic_genus - 1,
    toric_degree - 2 * (arithmetic_genus - 1),
    arithmetic_genus - 1,
]
trivial_lattice_rank = 2 + 4 * twist_degree
ambient_degree = 2 + 2 * twist_degree
full_h2_degree = 22 + 6 * twist_degree
active_exporter = EXPORTER if len(labels) == 2 else SINGLETON_EXPORTER
active_runner = RUNNER if len(labels) == 2 else SINGLETON_RUNNER

bisections = json.loads(args.bisections.read_text())
if bisections.get("schema") != "elkies-k3.bisection-extension-input.v1":
    raise ValueError("unexpected bisection schema")
by_label = {record["label"]: record for record in bisections["bisections"]}
if any(label not in by_label for label in labels):
    raise ValueError("twist label absent from bisection input")
Rt = PolynomialRing(field, "t")
t = Rt.gen()
quadratics = tuple(
    Rt([field(QQ(value)) for value in by_label[label]["branch"]["numerator_coefficients"]])
    for label in labels
)
d = quadratics[0]
for quadratic in quadratics[1:]:
    d *= quadratic

model = json.loads(args.model.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("unexpected direct-model status")
weierstrass = model["weierstrass_model"]
A = Rt([field(QQ(value)) for value in weierstrass["A_coefficients_low_to_high"]])
B = Rt([field(QQ(value)) for value in weierstrass["B_coefficients_low_to_high"]])
base_discriminant = -field(16) * (field(4) * A**3 + field(27) * B**2)
if (
    any(q.degree() != 2 or not q.is_squarefree() for q in quadratics)
    or d.degree() != twist_degree
    or not d.is_squarefree()
    or A.degree() != 8
    or B.degree() != 12
    or base_discriminant.degree() != 24
    or not base_discriminant.is_squarefree()
    or d.gcd(base_discriminant).degree()
):
    raise ArithmeticError("good-reduction gate failed")

parsed_output = parse_readfile_output(args.toric_output.read_text())
output_label = parsed_output["label"]
output_monomials = parsed_output["monomials"]
output_coefficients = parsed_output["coefficients"]
output_halfspace_A = parsed_output["halfspace_A"]
output_halfspace_b = parsed_output["halfspace_b"]
output_prime = parsed_output["prime"]
hodge_numbers = parsed_output["hodge_numbers"]
full_coefficients = [ZZ(value) for value in parsed_output["frobenius_coefficients"]]
if output_prime != prime or hodge_numbers != expected_hodge_numbers:
    raise ArithmeticError("unexpected toric prime or Hodge vector")

expected_terms = {(0, 3, 0): -field.one()}
for index, coefficient in enumerate(d):
    if coefficient:
        expected_terms[(index, 0, 2)] = coefficient
for index, coefficient in enumerate(-A):
    if coefficient:
        expected_terms[(index, 1, 0)] = coefficient
for index, coefficient in enumerate(-B):
    if coefficient:
        expected_terms[(index, 0, 0)] = coefficient
expected_monomials = sorted(expected_terms)
expected_coefficients = [int(expected_terms[monomial]) for monomial in expected_monomials]
first_facet_y_coefficient = 6 - twist_degree // 2
expected_halfspace_A = [
    (0, 0, 1),
    (0, 1, 0),
    (-1, -4, -first_facet_y_coefficient),
    (0, -2, -3),
    (1, 0, 0),
]
expected_halfspace_b = [0, 0, 12, 6, 0]
if (
    output_monomials != expected_monomials
    or output_coefficients != expected_coefficients
    or output_halfspace_A != expected_halfspace_A
    or output_halfspace_b != expected_halfspace_b
):
    raise ArithmeticError("toric output input data does not match the exact product twist")

Rpoly = PolynomialRing(ZZ, "T")
T = Rpoly.gen()
full_polynomial = Rpoly(full_coefficients)
if full_polynomial.degree() != toric_degree or full_polynomial.leading_coefficient() != 1:
    raise ArithmeticError("unexpected primitive toric Frobenius degree")

# Determine the closed-point degrees in D and in the finite etale degree-three
# cover Z -> D.  A degree-e point of D followed by a degree-f factor of the
# cubic over its residue field gives a degree-e*f point of Z.
branch_degrees = []
cover_degrees = []
branch_factor_records = []
for branch_index, (branch_factor, multiplicity) in enumerate(d.factor()):
    if multiplicity != 1:
        raise ArithmeticError("branch divisor is not reduced")
    extension_degree = int(branch_factor.degree())
    branch_degrees.append(extension_degree)
    if extension_degree == 1:
        residue_field = field
        branch_root = -branch_factor[0] / branch_factor[1]
    else:
        residue_field = GF(
            prime**extension_degree,
            name=f"a{branch_index}",
            modulus=branch_factor,
        )
        branch_root = residue_field.gen()
    Rx = PolynomialRing(residue_field, "x")
    x = Rx.gen()
    cubic = x**3 + residue_field(A(branch_root)) * x + residue_field(B(branch_root))
    cubic_factors = [
        int(factor.degree())
        for factor, factor_multiplicity in cubic.factor()
        for _unused in range(int(factor_multiplicity))
    ]
    if sum(cubic_factors) != 3 or not cubic.is_squarefree():
        raise ArithmeticError("two-torsion cover is not finite etale above the branch divisor")
    point_degrees = [extension_degree * value for value in cubic_factors]
    cover_degrees.extend(point_degrees)
    branch_factor_records.append(
        {
            "branch_degree": extension_degree,
            "cubic_factor_degrees_over_residue_field": cubic_factors,
            "closed_point_degrees_in_cover": point_degrees,
        }
    )

boundary_numerator = Rpoly.one()
for degree in cover_degrees:
    boundary_numerator *= T**degree - prime**degree
boundary_denominator = Rpoly.one()
for degree in branch_degrees:
    boundary_denominator *= T**degree - prime**degree
if boundary_numerator % boundary_denominator:
    raise ArithmeticError("boundary permutation quotient is not polynomial")
boundary_polynomial = boundary_numerator // boundary_denominator
if boundary_polynomial.degree() != boundary_degree or full_polynomial % boundary_polynomial:
    raise ArithmeticError("toric boundary factor does not divide Frobenius")

elliptic_polynomial = full_polynomial // boundary_polynomial
if elliptic_polynomial.degree() != elliptic_degree or elliptic_polynomial.leading_coefficient() != 1:
    raise ArithmeticError("elliptic Frobenius quotient has the wrong degree")
power_sum_1 = -elliptic_polynomial[elliptic_degree - 1]
power_sum_2 = (
    elliptic_polynomial[elliptic_degree - 1] ** 2
    - 2 * elliptic_polynomial[elliptic_degree - 2]
)

if len(labels) == 2:
    audit = json.loads(args.audit.read_text())
    audit_target = next(
        record for record in audit["targets"] if record["pair_key"] == args.pair_key
    )
    audit_reduction = next(
        record for record in audit_target["reductions"] if int(record["prime"]) == prime
    )
    if [int(power_sum_1), int(power_sum_2)] != [
        int(value)
        for value in audit_reduction["elliptic_L_frobenius_power_sums_n1_n2"]
    ]:
        raise ArithmeticError("complete toric Frobenius disagrees with audited power sums")
    independent_moment_check = "PASS_AGAINST_STORED_FIBREWISE_N1_N2_AUDIT"
else:
    singleton_audit = json.loads(args.singleton_audit.read_text())
    singleton_target = next(
        record for record in singleton_audit["targets"] if record["label"] == labels[0]
    )
    singleton_reduction = next(
        record
        for record in singleton_target["reductions"]
        if int(record["prime"]) == prime
    )
    if [int(power_sum_1), int(power_sum_2)] != [
        int(value)
        for value in singleton_reduction["elliptic_L_frobenius_power_sums_n1_n2"]
    ]:
        raise ArithmeticError("complete toric Frobenius disagrees with singleton audit")
    independent_moment_check = "PASS_AGAINST_STORED_FIBREWISE_N1_N2_AUDIT"

Rnorm = PolynomialRing(QQ, "Z")
Z = Rnorm.gen()
normalized = Rnorm(elliptic_polynomial(prime * Z) / prime**elliptic_degree)
functional_equation_sign = ZZ(normalized[0])
if functional_equation_sign not in (-1, 1):
    raise ArithmeticError("normalized constant term is not a functional-equation sign")
if normalized != functional_equation_sign * Z**elliptic_degree * normalized(1 / Z):
    raise ArithmeticError("weight-two functional equation failed")

# Certify the Weil absolute-value condition exactly.  Strip the endpoint roots
# forced by an anti-reciprocal sign.  The remaining reciprocal polynomial Q0
# has all roots on the unit circle precisely when Q0(Z)/Z^m=R(Z+1/Z) and every
# root of R lies in [-2,2].  Sage's AA root isolation and comparisons are
# certified.
weil_core = normalized
endpoint_multiplicities = {}
for endpoint in (-1, 1):
    divisor = Z - endpoint
    multiplicity = 0
    while not weil_core(endpoint):
        quotient, remainder = weil_core.quo_rem(divisor)
        if remainder:
            raise ArithmeticError("endpoint root division failed")
        weil_core = quotient
        multiplicity += 1
    endpoint_multiplicities[str(endpoint)] = multiplicity
if weil_core.degree() % 2:
    raise ArithmeticError("Weil core has odd degree")
core_half_degree = weil_core.degree() // 2
core_sign = QQ(weil_core[0] / weil_core.leading_coefficient())
if core_sign != 1 or weil_core != Z ** weil_core.degree() * weil_core(1 / Z):
    raise ArithmeticError("endpoint-stripped Weil core is not reciprocal")
Rtrace = PolynomialRing(QQ, "W")
W = Rtrace.gen()
chebyshev_traces = [Rtrace(2), W]
for index in range(2, core_half_degree + 1):
    chebyshev_traces.append(W * chebyshev_traces[-1] - chebyshev_traces[-2])
trace_polynomial = Rtrace(weil_core[core_half_degree])
for index in range(1, core_half_degree + 1):
    trace_polynomial += weil_core[core_half_degree + index] * chebyshev_traces[index]
if weil_core != Z**core_half_degree * trace_polynomial(Z + 1 / Z):
    raise ArithmeticError("reciprocal trace-polynomial reconstruction failed")
trace_roots = trace_polynomial.roots(ring=AA)
trace_root_count = sum(int(multiplicity) for _root, multiplicity in trace_roots)
if trace_root_count != core_half_degree or any(
    root < -2 or root > 2 for root, _multiplicity in trace_roots
):
    raise ArithmeticError("elliptic quotient fails the exact Weil-circle test")

# Reconstruct the missing ambient toric divisor classes, the complete trivial
# lattice, and hence all of H^2.  Here D is the twist branch divisor
# and Z is its degree-three two-torsion cover.  The I0* root lattice contributes
# H^0(D) (central components) plus H^0(Z) (the three outer components).
ambient_polynomial = (T - prime) ** 2 * boundary_denominator**2
trivial_lattice_polynomial = (
    (T - prime) ** 2 * boundary_denominator * boundary_numerator
)
full_h2_polynomial = ambient_polynomial * full_polynomial
if (
    ambient_polynomial.degree() != ambient_degree
    or trivial_lattice_polynomial.degree() != trivial_lattice_rank
    or full_h2_polynomial.degree() != full_h2_degree
    or full_h2_polynomial != trivial_lattice_polynomial * elliptic_polynomial
):
    raise ArithmeticError("full H^2/trivial-lattice reconstruction failed")

# The elementary bound phi(m)>=sqrt(m/2) gives a finite exhaustive loop.
cyclotomic_search_bound = 2 * elliptic_degree**2
cyclotomic_hits = []
for order in range(1, cyclotomic_search_bound + 1):
    if euler_phi(order) > elliptic_degree:
        continue
    cyclotomic = Rnorm(cyclotomic_polynomial(order))
    common = normalized.gcd(cyclotomic)
    if common.degree():
        multiplicity = 0
        remaining = normalized
        while True:
            quotient, remainder = remaining.quo_rem(cyclotomic)
            if remainder:
                break
            remaining = quotient
            multiplicity += 1
        cyclotomic_hits.append(
            {
                "order": order,
                "degree": int(common.degree()),
                "multiplicity": multiplicity,
                "total_degree": int(common.degree()) * multiplicity,
                "factor_coefficients_low_to_high": [str(value) for value in common],
            }
        )
cyclotomic_degree = sum(record["total_degree"] for record in cyclotomic_hits)
picard_upper_bound = trivial_lattice_rank + cyclotomic_degree
geometric_mw_rank_upper_bound = cyclotomic_degree
known_mw_rank_lower_bound = 0 if len(labels) == 2 else 1
closed_rank_zero = len(labels) == 2 and cyclotomic_degree == 0
closed_singleton_rank_one = len(labels) == 1 and cyclotomic_degree == 1

if args.output is None:
    tag = args.pair_key.replace(":", "--")
    args.output = (
        ROOT
        / "artifacts/generated-results"
        / (
            f"elkies-k3-r17-product-{tag}-p{prime}-toric-frobenius-v1.json"
            if len(labels) == 2
            else f"elkies-k3-r17-singleton-{tag}-p{prime}-toric-frobenius-v1.json"
        )
    )
args.output.parent.mkdir(parents=True, exist_ok=True)
record = {
    "schema": (
        "elkies-k3.r17-product-toric-frobenius.v1"
        if len(labels) == 2
        else "elkies-k3.r17-singleton-toric-frobenius.v1"
    ),
    "status": (
        "PASS_GEOMETRIC_PRODUCT_TWIST_RANK_ZERO"
        if closed_rank_zero
        else (
            "PASS_GEOMETRIC_SINGLETON_TWIST_RANK_ONE"
            if closed_singleton_rank_one
            else "PASS_COMPLETE_FROBENIUS_PICARD_BOUND"
        )
    ),
    "claim": (
        "The product twist has geometric Mordell-Weil rank zero."
        if closed_rank_zero
        else (
            "The singleton twist has geometric Mordell-Weil rank one."
            if closed_singleton_rank_one
            else "The full Frobenius polynomial is certified and gives the displayed rank interval."
        )
    ),
    "character_kind": "product" if len(labels) == 2 else "singleton",
    "pair_key": args.pair_key,
    "labels": list(labels),
    "prime": prime,
    "good_reduction": {
        "status": "PASS",
        "checks": [
            "p is an odd prime at least five",
            "every quadratic character retains degree two and is squarefree",
            f"their product retains degree {twist_degree} and is squarefree",
            "A, B, and Delta retain degrees 8, 12, and 24",
            "Delta is squarefree and coprime to the twist divisor",
            "the exact toric Jacobian nondegeneracy calculation passes",
        ],
        "geometric_fibre_configuration": f"{twist_degree}I0*+24I1",
    },
    "toric_model": "d(t)*y^2=x^3+A(t)*x+B(t)",
    "toric_output_label": output_label,
    "toric_hodge_numbers": hodge_numbers,
    "toric_nondegeneracy_jacobian_reduction": "PASS",
    "toric_primitive_frobenius_coefficients_low_to_high": [
        str(value) for value in full_polynomial
    ],
    "boundary": {
        "description": "p*(H^0(Z)-H^0(D)), D={d=0}, Z={d=0,x^3+A*x+B=0}",
        "branch_factor_records": branch_factor_records,
        "H0_D_frobenius_coefficients_low_to_high": [
            str(value) for value in boundary_denominator
        ],
        "H0_Z_frobenius_coefficients_low_to_high": [
            str(value) for value in boundary_numerator
        ],
        "frobenius_coefficients_low_to_high": [str(value) for value in boundary_polynomial],
    },
    "H2": {
        "degree": full_h2_degree,
        "ambient_toric_divisor_factor_coefficients_low_to_high": [
            str(value) for value in ambient_polynomial
        ],
        "trivial_lattice_factor_description": (
            "(T-p)^2*det(T-pF|H^0(D))*det(T-pF|H^0(Z))"
        ),
        "trivial_lattice_frobenius_coefficients_low_to_high": [
            str(value) for value in trivial_lattice_polynomial
        ],
        "full_frobenius_coefficients_low_to_high": [
            str(value) for value in full_h2_polynomial
        ],
        "factorization_identity": "P_H2=P_trivial*P_elliptic",
    },
    "elliptic_L": {
        "degree": elliptic_degree,
        "frobenius_characteristic_coefficients_low_to_high": [
            str(value) for value in elliptic_polynomial
        ],
        "power_sums_n1_n2": [str(power_sum_1), str(power_sum_2)],
        "independent_power_sum_check": independent_moment_check,
        "functional_equation_sign": int(functional_equation_sign),
        "weil_circle_check": {
            "status": "PASS_EXACT_REAL_ROOT_ISOLATION",
            "trace_polynomial_degree": int(trace_polynomial.degree()),
            "stripped_endpoint_multiplicities": endpoint_multiplicities,
            "trace_polynomial_coefficients_low_to_high": [
                str(value) for value in trace_polynomial
            ],
            "roots_in_closed_interval_minus2_2_counted_with_multiplicity": trace_root_count,
        },
        "cyclotomic_search_bound": cyclotomic_search_bound,
        "cyclotomic_hits_after_T_equals_pZ": cyclotomic_hits,
    },
    "bounds": {
        "trivial_lattice_rank": trivial_lattice_rank,
        "picard_number_upper_bound": picard_upper_bound,
        "geometric_twist_mw_rank_lower_bound": known_mw_rank_lower_bound,
        "geometric_twist_mw_rank_upper_bound": geometric_mw_rank_upper_bound,
        "derivation": (
            f"rho(Fpbar)<={trivial_lattice_rank}+tate_degree; Shioda-Tate gives "
            f"rank(MW)<=rho-2-{twist_degree}*rank(D4)=tate_degree; the explicit "
            f"U+{twist_degree}D4 lattice gives rho>={trivial_lattice_rank}"
        ),
        "characteristic_zero_specialization": (
            f"NS(Qbar) injects into NS(Fpbar), while U+{twist_degree}D4 exists in characteristic zero"
        ),
    },
    "software": {
        "ToricControlledReduction_repository": "https://github.com/edgarcosta/ToricControlledReduction",
        "ToricControlledReduction_commit": args.toric_commit,
        "ToricControlledReduction_output_sha256": digest(args.toric_output),
        "ToricControlledReduction_input_sha256": digest(
            args.toric_output.with_name("toric-controlled-reduction.input")
        ),
        "ToricControlledReduction_executable_sha256": digest(args.toric_executable),
        "exporter_sha256": digest(active_exporter),
        "raw_output_parser_sha256": digest(PARSER),
        "independent_verifier_sha256": digest(VERIFIER),
        "runner_sha256": digest(active_runner),
        "sage_prefix": str(args.sage_prefix.resolve()),
        "build_invocations": [
            shlex.join(
                [
                    "./configure",
                    "--disable-gdb",
                    f"--with-ntl={args.sage_prefix.resolve()}",
                    f"--with-gmp={args.sage_prefix.resolve()}",
                ]
            ),
            "make",
        ],
        "frobenius_invocation": shlex.join(
            [
                str(args.toric_executable.resolve()),
                str(args.toric_output.with_name("toric-controlled-reduction.input").resolve()),
                str(args.toric_output.resolve()),
            ]
        ),
        "certificate_invocation": shlex.join(sys.argv),
    },
    "inputs": {
        str(path.resolve().relative_to(ROOT)): digest(path)
        for path in (
            (args.bisections, args.model, args.audit)
            if len(labels) == 2
            else (args.bisections, args.model, args.singleton_audit)
        )
    },
}
args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
key_field = "pair" if len(labels) == 2 else "label"
print(
    f"R17TORICFROB|{key_field}={args.pair_key}|p={prime}|degree={elliptic_degree}"
    f"|tate_degree={cyclotomic_degree}|rho_upper={picard_upper_bound}"
    f"|mw_upper={geometric_mw_rank_upper_bound}|output={args.output}"
    f"|status={record['status']}",
    flush=True,
)
