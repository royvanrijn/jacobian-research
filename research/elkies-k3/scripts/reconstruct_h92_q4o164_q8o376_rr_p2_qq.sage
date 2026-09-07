#!/usr/bin/env sage -python
"""Reconstruct and exactly compile the q8/orbit376 RR direction over QQ.

Inputs are complete finite-field P^2 scans produced by
``probe_h92_q4o164_q8o376_rr_p2_modp.sage``. The unique semistable-4A1
candidate in each good-prime scan determines one projective quadratic

    BB=b0+b1*t+b2*t^2.

The normalized projective coordinates are reconstructed by coefficientwise
CRT/rational reconstruction and, when necessary, one tiny projective LLL.
The candidate is then accepted only after rebuilding over QQ:

* the smooth-collision recurrence AA=BB*Y/X mod Z^2;
* exact Z^6 collision removal;
* exact vertical square content;
* a degree-four binary quartic;
* its invariant Jacobian with semistable 4A1 fibres;
* reduction back to every construction and held-out prime.

This produces an exact unpointed equation candidate. Attaching the selected
P1229 origin and the marked NS edge remains a separate certification gate.
"""

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import CRT_list, GF, PolynomialRing, QQ, ZZ, gcd, inverse_mod, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument(
    "--held-out",
    action="append",
    type=Path,
    default=[],
    help="complete scan excluded from CRT and used only for replay",
)
parser.add_argument(
    "--selector",
    choices=("unique-4a1", "unique-candidate"),
    default="unique-4a1",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q4o164-q8o376-rr-p2-qq.json",
)
args = parser.parse_args()

inputs = args.inputs or sorted(LOCAL.glob("q4o164-q8o376-rr-p2-scan-mod*.json"))
inputs = [path if path.is_absolute() else ROOT / path for path in inputs]
held_out = [path if path.is_absolute() else ROOT / path for path in args.held_out]
output = args.output if args.output.is_absolute() else ROOT / args.output
if not inputs:
    raise SystemExit("no modular P2 scan inputs supplied")
if len(inputs) < 2:
    raise SystemExit("use at least two construction primes")
if set(inputs) & set(held_out):
    raise SystemExit("construction and held-out input sets must be disjoint")
for path in inputs + held_out + [MODEL, HORIZONTAL]:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q8O376RRQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


def scan_complete(record):
    search = record["search_space"]
    processed = int(search.get("processed", search.get("scanned", 0)))
    return bool(
        search.get("complete", False)
        or (
            search.get("stop_after") is None
            and processed == int(search["expected_size"])
        )
    )


def select_candidate(record):
    if not scan_complete(record):
        raise ArithmeticError(f"prime {record['prime']} scan is incomplete")
    candidates = record["candidates"]
    if args.selector == "unique-4a1":
        selected = [
            item
            for item in candidates
            if item.get("child")
            and item["child"].get("semistable_4A1_fingerprint")
        ]
    else:
        selected = candidates
    if len(selected) != 1:
        raise ArithmeticError(
            f"prime {record['prime']} has {len(selected)} candidates under selector {args.selector}"
        )
    return selected[0]


records = [json.loads(path.read_text()) for path in inputs]
held_records = [json.loads(path.read_text()) for path in held_out]
selected = [select_candidate(record) for record in records]
held_selected = [select_candidate(record) for record in held_records]
primes = [ZZ(record["prime"]) for record in records]
held_primes = [ZZ(record["prime"]) for record in held_records]
if len(set(primes + held_primes)) != len(primes) + len(held_primes):
    raise ArithmeticError("duplicate construction or held-out prime")

# Every scan must reduce the same exact inputs.
input_hash_sets = [record["inputs"]["sha256"] for record in records + held_records]
if input_hash_sets and any(item != input_hash_sets[0] for item in input_hash_sets[1:]):
    raise ArithmeticError("modular scans do not share one exact parent/horizontal input set")


def chart_index(values):
    return next(index for index, value in enumerate(values) if value)


projective_rows = [item["projective_BB_coefficients_low_to_high"] for item in selected]
charts = [
    int(item.get("projective_chart_index", chart_index(row)))
    for item, row in zip(selected, projective_rows)
]
if len(set(charts)) != 1:
    raise ArithmeticError(f"projective chart changes across construction primes: {charts}")
pivot = charts[0]
if any(row[pivot] % prime != 1 for row, prime in zip(projective_rows, primes)):
    raise ArithmeticError("selected modular BB rows are not canonically normalized")
free_indices = [index for index in range(3) if index != pivot]
modulus = math.prod(primes)
crt_residues = [
    ZZ(CRT_list([row[index] for row in projective_rows], primes))
    for index in free_indices
]


def vector_bits(values):
    return max(
        max(abs(QQ(value).numerator()).nbits(), QQ(value).denominator().nbits())
        for value in values
    )


def insert_pivot(values):
    result = []
    iterator = iter(values)
    for index in range(3):
        result.append(QQ.one() if index == pivot else QQ(next(iterator)))
    return vector(QQ, result)


def canonical_mod(values, prime):
    field = GF(prime)
    reduced = []
    for value in values:
        value = QQ(value)
        if value.denominator() % prime == 0:
            raise ArithmeticError(f"bad reconstructed denominator modulo {prime}")
        reduced.append(field(value.numerator()) / field(value.denominator()))
    local_pivot = next(index for index, value in enumerate(reduced) if value)
    scale = reduced[local_pivot]
    return tuple(int(value / scale) for value in reduced), local_pivot


def matches_all(values, rows, row_primes):
    for row, prime in zip(rows, row_primes):
        reduced, local_pivot = canonical_mod(values, prime)
        if local_pivot != pivot or reduced != tuple(int(value) for value in row):
            return False
    return True


projective_candidates = []
try:
    coefficientwise = insert_pivot(
        [residue.rational_reconstruction(modulus) for residue in crt_residues]
    )
    projective_candidates.append({
        "method": "coefficientwise rational reconstruction",
        "BB": coefficientwise,
        "primitive_vector_max_bits": vector_bits(coefficientwise),
        "projective_scale": "1",
    })
except ArithmeticError:
    pass

# Simultaneous projective reconstruction in dimension three. A short row is
# scale*(r0,r1,1) modulo the CRT modulus.
basis = matrix(ZZ, 3, 3)
basis[0, 0] = modulus
basis[1, 1] = modulus
basis[2] = vector(ZZ, list(crt_residues) + [1])
for row in sorted(basis.LLL(delta=0.99).rows(), key=lambda item: item * item):
    scale = ZZ(row[-1])
    if not scale or gcd(scale, modulus) != 1:
        continue
    if any((row[index] - scale * crt_residues[index]) % modulus for index in range(2)):
        continue
    candidate = insert_pivot([QQ(row[index]) / scale for index in range(2)])
    projective_candidates.append({
        "method": "simultaneous projective LLL",
        "BB": candidate,
        "primitive_vector_max_bits": max(abs(ZZ(value)).nbits() for value in row),
        "projective_scale": str(scale),
        "lll_norm_squared": str(row * row),
    })

unique_candidates = []
seen = set()
for record in sorted(
    projective_candidates,
    key=lambda item: (item["primitive_vector_max_bits"], item["method"]),
):
    key = tuple(record["BB"])
    if key in seen or not matches_all(record["BB"], projective_rows, primes):
        continue
    seen.add(key)
    unique_candidates.append(record)
if not unique_candidates:
    raise ArithmeticError("CRT/LLL produced no projective BB vector replaying every input prime")
log(
    "RECONSTRUCT",
    primes=",".join(map(str, primes)),
    modulus_bits=ZZ(modulus).nbits(),
    candidates=len(unique_candidates),
    chart=pivot,
)

model = json.loads(MODEL.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"

TQ = PolynomialRing(QQ, "t")
t = TQ.gen()
UQ = PolynomialRing(QQ, "u")
u = UQ.gen()
RT = PolynomialRing(UQ, "t")

A = TQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = TQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
section = horizontal["section"]


def polynomial_from_record(record, key):
    return TQ([QQ(value) for value in record[key]])


def exact_power_root(polynomial, exponent):
    polynomial = polynomial.monic()
    answer = polynomial.parent().one()
    for factor, multiplicity in polynomial.factor():
        multiplicity = int(multiplicity)
        if multiplicity % exponent:
            raise ArithmeticError(
                f"QQ polynomial is not an exact {exponent}-th power: multiplicity {multiplicity}"
            )
        answer *= factor.monic() ** (multiplicity // exponent)
    if answer**exponent != polynomial:
        raise ArithmeticError("exact QQ polynomial power-root verification failed")
    return answer.monic()


X = polynomial_from_record(section["x"], "numerator_coefficients_low_to_high")
x_den = polynomial_from_record(section["x"], "denominator_coefficients_low_to_high").monic()
Y = polynomial_from_record(section["y"], "numerator_coefficients_low_to_high")
y_den = polynomial_from_record(section["y"], "denominator_coefficients_low_to_high").monic()
Zx = exact_power_root(x_den, 2)
Zy = exact_power_root(y_den, 3)
if Zx != Zy:
    raise ArithmeticError("exact x/y denominators do not share one projective Z")
Z = Zx
if (X.degree(), Z.degree(), Y.degree()) != (12, 4, 18):
    raise ArithmeticError("exact q8 horizontal degree profile changed")
if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
    raise ArithmeticError("exact q8 horizontal misses compact parent")
if gcd(X, Z) != 1:
    raise ArithmeticError("exact q8 horizontal pole divisor is not smooth")

Z2 = Z**2
X_inverse = inverse_mod(X, Z2)
residue = TQ((Y * X_inverse) % Z2)


def lift(poly):
    return RT([UQ(value) for value in TQ(poly).list()])


def coefficients_in_u(poly):
    if not poly:
        return [TQ.zero()]
    u_degree = max(
        (coefficient.degree() for coefficient in poly.list() if coefficient),
        default=0,
    )
    return [
        TQ([QQ(coefficient[index]) for coefficient in poly.list()])
        for index in range(u_degree + 1)
    ]


def content_in_t(poly):
    coefficients = [item for item in coefficients_in_u(poly) if item]
    if not coefficients:
        return TQ.zero()
    answer = coefficients[0]
    for item in coefficients[1:]:
        answer = answer.gcd(item)
        if answer.degree() == 0:
            break
    return answer.monic()


def exact_square_root(polynomial):
    polynomial = polynomial.monic()
    answer = polynomial.parent().one()
    for factor, multiplicity in polynomial.factor():
        multiplicity = int(multiplicity)
        if multiplicity % 2:
            raise ArithmeticError("vertical content is not a square over QQ")
        answer *= factor.monic() ** (multiplicity // 2)
    if answer**2 != polynomial:
        raise ArithmeticError("vertical square-root verification failed")
    return answer.monic()


def polynomial_record(poly):
    return [str(value) for value in poly.list()]


def nested_record(poly):
    return [polynomial_record(value) for value in poly.list()]


def fibre_profile(A_child, B_child):
    Delta = UQ(-16 * (4 * A_child**3 + 27 * B_child**2))
    if not Delta or A_child.degree() > 8 or B_child.degree() > 12 or Delta.degree() > 24:
        raise ArithmeticError("child invariant model is not an elliptic K3 degree profile")
    degree_A = int(A_child.degree())
    degree_B = int(B_child.degree())
    degree_delta = int(Delta.degree())
    infinity_order = 24 - degree_delta
    infinity_multiplicative = bool(
        infinity_order == 0
        or (infinity_order in (1, 2) and degree_A == 8 and degree_B == 12)
    )
    if not infinity_multiplicative:
        raise ArithmeticError("child infinity fibre is not semistable")
    if A_child.gcd(Delta) != 1 or B_child.gcd(Delta) != 1:
        raise ArithmeticError("child has a finite additive fibre")

    repeated_gcd = Delta.gcd(Delta.derivative()).monic()
    repeated_support = (
        repeated_gcd.squarefree_part().monic()
        if repeated_gcd.degree() > 0
        else UQ.one()
    )
    nodal, remainder = Delta.quo_rem(repeated_support**2)
    if remainder:
        raise ArithmeticError("finite repeated support does not divide Delta exactly twice")
    if not nodal.is_squarefree() or nodal.gcd(repeated_support) != 1:
        raise ArithmeticError("child finite fibres have multiplicity greater than two")
    finite_double_degree = int(repeated_support.degree()) if repeated_support != 1 else 0
    finite_simple_degree = int(nodal.degree())
    geometric_root_rank = finite_double_degree + (1 if infinity_order == 2 else 0)
    geometric_I1_count = finite_simple_degree + (1 if infinity_order == 1 else 0)
    if geometric_root_rank != 4 or geometric_I1_count != 16:
        raise ArithmeticError(
            f"child is not semistable 4A1: roots={geometric_root_rank}, I1={geometric_I1_count}"
        )
    return {
        "degrees_A_B_Delta": [degree_A, degree_B, degree_delta],
        "Delta_coefficients_low_to_high": polynomial_record(Delta),
        "finite_I2_support_polynomial_coefficients_low_to_high": polynomial_record(repeated_support),
        "finite_I1_support_polynomial_coefficients_low_to_high": polynomial_record(nodal),
        "finite_I2_count_geometric": finite_double_degree,
        "finite_I1_count_geometric": finite_simple_degree,
        "infinity_delta_order": infinity_order,
        "infinity_kodaira": "smooth" if infinity_order == 0 else f"I{infinity_order}",
        "ADE": "4A1",
        "root_rank": 4,
        "geometric_I1_count": 16,
        "euler_number": degree_delta + infinity_order,
    }


def reduce_list(poly, prime):
    field = GF(prime)
    values = []
    for coefficient in poly.list():
        coefficient = QQ(coefficient)
        if coefficient.denominator() % prime == 0:
            raise ArithmeticError(f"bad exact child denominator modulo {prime}")
        values.append(int(field(coefficient.numerator()) / field(coefficient.denominator())))
    return values


def verify_scan_replay(compiled, scan_records, scan_candidates):
    for scan, candidate in zip(scan_records, scan_candidates):
        prime = ZZ(scan["prime"])
        bb_list = list(compiled["BB"].list())
        bb_values = bb_list + [QQ.zero()] * (3 - len(bb_list))
        reduced_bb, local_pivot = canonical_mod(bb_values, prime)
        if local_pivot != pivot or reduced_bb != tuple(candidate["projective_BB_coefficients_low_to_high"]):
            return False
        child = candidate.get("child")
        if child:
            if reduce_list(compiled["A_child"], prime) != child["A_coefficients_low_to_high"]:
                return False
            if reduce_list(compiled["B_child"], prime) != child["B_coefficients_low_to_high"]:
                return False
    return True


def compile_direction(projective_record):
    BB = TQ(projective_record["BB"])
    if not BB or BB.degree() > 2:
        raise ArithmeticError("reconstructed BB is not a nonzero quadratic")
    AA = TQ((BB * residue) % Z2)
    if AA.degree() > 8 or (AA * X - BB * Y) % Z2:
        raise ArithmeticError("reconstructed direction fails the collision recurrence")

    A_nested, X_nested, Y_nested, Z_nested, AA_nested, BB_nested = map(
        lift, (A, X, Y, Z, AA, BB)
    )
    N = AA_nested - RT(u) * Z_nested**2
    Db = -BB_nested
    raw = (
        N**4
        - 6 * X_nested * N**2 * Db**2
        - 8 * Y_nested * N * Db**3
        - 3 * X_nested**2 * Db**4
        - 4 * A_nested * Z_nested**4 * Db**4
    )
    after_collision, remainder = raw.quo_rem(Z_nested**6)
    if remainder or after_collision.degree() > 8:
        raise ArithmeticError("exact raw chord class does not remove Z^6")
    content = content_in_t(after_collision)
    square_root = exact_square_root(content)
    if content.degree() > 4:
        raise ArithmeticError("exact finite vertical square has degree greater than four")
    infinity_square_degree = 4 - int(content.degree())
    if infinity_square_degree < 0 or infinity_square_degree % 2:
        raise ArithmeticError("exact infinity vertical square degree is inconsistent")
    quartic, remainder = after_collision.quo_rem(lift(content))
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("exact vertical square removal does not leave a quartic")

    values = list(quartic.list()) + [UQ.zero()] * 5
    e, d, c, b, a = values[:5]
    I = UQ(12 * a * e - 3 * b * d + c**2)
    J = UQ(72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3)
    A_child = UQ(-27 * I)
    B_child = UQ(-27 * J)
    profile = fibre_profile(A_child, B_child)
    compiled = {
        **projective_record,
        "BB": BB,
        "AA": AA,
        "finite_vertical_square_root": square_root,
        "finite_vertical_square_degree": int(content.degree()),
        "infinity_vertical_square_degree": infinity_square_degree,
        "quartic": quartic,
        "I": I,
        "J": J,
        "A_child": A_child,
        "B_child": B_child,
        "fibre_profile": profile,
    }
    if not verify_scan_replay(compiled, records, selected):
        raise ArithmeticError("exact child does not replay every construction-prime scan")
    if held_records and not verify_scan_replay(compiled, held_records, held_selected):
        raise ArithmeticError("exact child misses a held-out scan")
    return compiled


compiled = None
failures = []
for projective_record in unique_candidates:
    try:
        compiled = compile_direction(projective_record)
        break
    except ArithmeticError as error:
        failures.append({
            "method": projective_record["method"],
            "BB": [str(value) for value in projective_record["BB"]],
            "reason": str(error),
        })
if compiled is None:
    raise ArithmeticError(f"no reconstructed direction compiled exactly: {failures}")

log(
    "COMPILE",
    BB=",".join(str(value) for value in compiled["BB"].list()),
    finite_square=compiled["finite_vertical_square_degree"],
    infinity_square=compiled["infinity_vertical_square_degree"],
    child="4A1",
)

payload = {
    "schema": "elkies-k3.q4o164-q8o376-rr-p2-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_Q8O376_UNPOINTED_RR_AND_4A1_JACOBIAN",
    "reconstruction": {
        "construction_primes": list(map(int, primes)),
        "held_out_primes": list(map(int, held_primes)),
        "CRT_modulus": str(modulus),
        "CRT_modulus_bits": int(ZZ(modulus).nbits()),
        "projective_chart_index": pivot,
        "method": compiled["method"],
        "projective_scale": compiled["projective_scale"],
        "primitive_vector_max_bits": compiled["primitive_vector_max_bits"],
        "failed_exact_candidates": failures,
    },
    "resolved_RR": {
        "divisor_target": "q8/orbit376 from the C8-pointed q4/orbit164 model",
        "ambient_dimension": 12,
        "smooth_collision_rank": 8,
        "post_collision_dimension": 4,
        "vertical_condition_rank": 2,
        "kernel_dimension": 2,
        "basis": [
            {"AA_coefficients_low_to_high": polynomial_record(Z2), "BB_coefficients_low_to_high": []},
            {
                "AA_coefficients_low_to_high": polynomial_record(compiled["AA"]),
                "BB_coefficients_low_to_high": polynomial_record(compiled["BB"]),
            },
        ],
        "finite_vertical_square_root_coefficients_low_to_high": polynomial_record(
            compiled["finite_vertical_square_root"]
        ),
        "finite_vertical_square_degree": compiled["finite_vertical_square_degree"],
        "infinity_vertical_square_degree": compiled["infinity_vertical_square_degree"],
        "exact_Z6_collision_removal": True,
        "exact_vertical_square_removal": True,
    },
    "quartic": {
        "degree_in_old_t": int(compiled["quartic"].degree()),
        "coefficients_in_old_t_low_to_high": nested_record(compiled["quartic"]),
        "I_coefficients_low_to_high": polynomial_record(compiled["I"]),
        "J_coefficients_low_to_high": polynomial_record(compiled["J"]),
    },
    "child": {
        "A_coefficients_low_to_high": polynomial_record(compiled["A_child"]),
        "B_coefficients_low_to_high": polynomial_record(compiled["B_child"]),
        **compiled["fibre_profile"],
    },
    "checks": {
        "reduction_to_every_construction_prime": True,
        "reduction_to_every_held_out_prime": True,
        "exact_parent_horizontal_identity": True,
        "exact_binary_quartic_invariants": True,
        "exact_semistable_4A1_fibre_audit": True,
        "large_Groebner_required": False,
    },
    "proof_boundary": (
        "This is an exact characteristic-zero two-plane, quartic and semistable 4A1 invariant "
        "Jacobian reconstructed from complete finite-field scans and checked at every supplied "
        "prime. It is still unpointed: the degree-one P1229 curve must be restricted to the "
        "pencil, inverted, and used to attach the origin and the certified q8/orbit376 marked "
        "NS transport before promoting the equation edge."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs + held_out + [MODEL, HORIZONTAL]],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in inputs + held_out + [MODEL, HORIZONTAL]
        },
    },
    "runtime_seconds": time.monotonic() - started,
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("DONE", status=payload["status"], output=output)
