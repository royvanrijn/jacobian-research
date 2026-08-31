#!/usr/bin/env sage -python
"""Classify the complete p=89 polynomial section shell by lattice Abel words.

Invert each filtered polynomial P.O=0 section on the P1229-pointed q8 child
to a point (T(u),W(u)) on the resolved quartic and then to the corresponding
point (x(u),y(u)) on the old q4/o164 Weierstrass surface by the exact chord
formula.  At several ordinary old-base values, split T(u)-t0 over GF(p^6),
sum the conjugate q4 points, and compare literal two/three-fibre signatures
with all 938 physical lattice P.O=0 classes.  Then repeat the exact four-word
MW-quotient search using only classes actually realized by the modular shell.
No Groebner basis, multivariate elimination, or QQ lift is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
HEIGHT = LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json"
FRAME = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt"
FRONTIER = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json"
INPUTS_BASE = (Q8, MODEL, BASIS, C8, HORIZONTAL, Q3, HEIGHT, MARKING, FRAME, FRONTIER)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=89)
parser.add_argument("--shell", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument("--fibres", type=int, default=3)
args = parser.parse_args()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
SHELL = (
    args.shell.resolve() if args.shell else
    LOCAL / f"q12o5867-p0-shell-all-records-mod{prime}.json"
)
OUTPUT = (
    args.output.resolve() if args.output else
    LOCAL / f"q12o5867-p0-shell-lattice-classification-mod{prime}.json"
)
started = time.monotonic()
F = GF(prime)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
c8 = json.loads(C8.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
q3_exact = json.loads(Q3.read_text())
height = json.loads(HEIGHT.read_text())
marking = json.loads(MARKING.read_text())
frontier = json.loads(FRONTIER.read_text())
shell = json.loads(SHELL.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert shell["prime"] == int(prime)
assert shell["status"].startswith("PASS_MODP_Q12O5867_P0_SHELL")
assert "all_records" in shell and len(shell["all_records"]) == 300

PT = PolynomialRing(F, "T")
T = PT.gen()
PK = PT.fraction_field()
U = PolynomialRing(F, "u")
u = U.gen()
KU = U.fraction_field()
TU = PolynomialRing(KU, "T")
Tvar = TU.gen()
TUP = PolynomialRing(U, "T")


def reduce_qq(value):
    value = QQ(value)
    assert value.denominator() % prime
    return F(value.numerator())/F(value.denominator())


def parent_poly(values):
    return PT([reduce_qq(value) for value in values])


def child_poly(values):
    return U([reduce_qq(value) for value in values])


def parent_function(record):
    return PK(parent_poly(record["numerator_coefficients_low_to_high"])) / PK(
        parent_poly(record["denominator_coefficients_low_to_high"])
    )


def child_function(record):
    return KU(child_poly(record["numerator_coefficients_low_to_high"])) / KU(
        child_poly(record["denominator_coefficients_low_to_high"])
    )


parent_A = parent_poly(model["compact_model"]["A_coefficients_low_to_high"])
parent_B = parent_poly(model["compact_model"]["B_coefficients_low_to_high"])
child_A = child_poly(q8["child"]["minimal_A_coefficients_low_to_high"])
child_B = child_poly(q8["child"]["minimal_B_coefficients_low_to_high"])


def checked_parent_point(x_coordinate, y_coordinate):
    point = PK(x_coordinate), PK(y_coordinate)
    assert point[1]**2 == point[0]**3+PK(parent_A)*point[0]+PK(parent_B)
    return point


def parent_neg(point):
    return None if point is None else (point[0], -point[1])


def parent_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+PK(parent_A))/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return checked_parent_point(x3, slope*(x1-x3)-y1)


def parent_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return parent_mul(-coefficient, parent_neg(point))
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = parent_add(answer, addend)
        addend = parent_add(addend, addend)
        coefficient >>= 1
    return answer


parent_basis = [
    checked_parent_point(
        parent_poly(section["x_coefficients_low_to_high"]),
        parent_poly(section["y_coefficients_low_to_high"]),
    )
    for section in basis["resolved_hensel"]["sections"]
]
H = checked_parent_point(
    parent_function(horizontal["section"]["x"]),
    parent_function(horizontal["section"]["y"]),
)
c8_record = c8["opposite_constant_support_section"]
c8_old = parent_function(c8_record["x"]), parent_function(c8_record["y"])
base_scale = reduce_qq(model["exact_coordinate_change"]["c"])
xy_scale = reduce_qq(model["exact_coordinate_change"]["s"])


def substitute_scaled(function):
    return PK(PT(function.numerator()(base_scale*T))) / PK(
        PT(function.denominator()(base_scale*T))
    )


C8opposite = checked_parent_point(
    substitute_scaled(c8_old[0])/xy_scale**2,
    substitute_scaled(c8_old[1])/xy_scale**3,
)


def word_point(terms):
    answer = None
    for coefficient, point in terms:
        answer = parent_add(answer, parent_mul(coefficient, point))
    return answer


target_words = {
    "Q1": word_point(((1, H), (-1, parent_basis[3]))),
    "Q2": word_point((
        (1, H), (1, parent_basis[0]), (1, parent_basis[2]),
        (-1, parent_basis[3]), (1, parent_basis[4]),
    )),
    "Q4": word_point((
        (1, H), (-1, C8opposite), (-1, parent_basis[1]),
        (-1, parent_basis[2]), (1, parent_basis[3]),
        (1, parent_basis[5]), (-1, parent_basis[6]), (-1, parent_basis[7]),
    )),
    "Q3": word_point((
        (2, C8opposite), (1, parent_basis[1]), (2, parent_basis[2]),
        (-2, parent_basis[3]), (1, parent_basis[4]), (-2, parent_basis[5]),
        (1, parent_basis[6]), (1, parent_basis[7]),
    )),
}

# Reconstruct the normalized q8 quartic square factor.
XH = parent_poly(horizontal["section"]["x"]["numerator_coefficients_low_to_high"])
denominator_x = parent_poly(horizontal["section"]["x"]["denominator_coefficients_low_to_high"])
YH = parent_poly(horizontal["section"]["y"]["numerator_coefficients_low_to_high"])
denominator_y = parent_poly(horizontal["section"]["y"]["denominator_coefficients_low_to_high"])
Z = PT.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic()**(int(exponent)//2)
assert Z**2 == denominator_x and Z**3 == denominator_y
resolved_pairs = [
    (
        parent_poly(item["AA_coefficients_low_to_high"]),
        parent_poly(item["BB_coefficients_low_to_high"]),
    )
    for item in q8["resolved_RR"]["resolved_basis_pairs"]
]
(AA0, BB0), (AA1, BB1) = resolved_pairs


def lift_TU(poly):
    return TU([KU(value) for value in PT(poly).list()])


aa = lift_TU(AA0)+KU(u)*lift_TU(AA1)
bb = lift_TU(BB0)+KU(u)*lift_TU(BB1)
X_u, Y_u, Z_u, A_u = map(lift_TU, (XH, YH, Z, parent_A))
raw = (
    aa**4-6*X_u*aa**2*bb**2+8*Y_u*aa*bb**3
    -3*X_u**2*bb**4-4*A_u*bb**4*Z_u**4
)
after_collision, remainder = raw.quo_rem(Z_u**4)
assert not remainder
quartic = TUP([
    U([reduce_qq(value) for value in coefficient])
    for coefficient in q8["quartic"]["coefficients_in_old_T_low_to_high"]
])
square_quotient, remainder = TUP(after_collision).quo_rem(quartic)
assert not remainder
factorization = square_quotient.factor()
assert all(int(exponent) % 2 == 0 for unused, exponent in factorization)
unit = F(factorization.unit())
assert unit.is_square()
square_factor = TUP(unit.sqrt())
for factor, exponent in factorization:
    square_factor *= factor**(int(exponent)//2)
assert TUP(after_collision) == quartic*square_factor**2


def evaluate_bivariate(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(TUP(poly).list()):
        answer = answer*old_t+KU(U(coefficient))
    return answer


def evaluate_parent_poly(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(PT(poly).list()):
        answer = answer*old_t+KU(coefficient)
    return answer


def compose_parent(function, old_t):
    function = PK(function)
    return evaluate_parent_poly(function.numerator(), old_t) / evaluate_parent_poly(
        function.denominator(), old_t
    )


# Pointed-quartic inverse on the child.
pointing = q8["preferred_pointed_zero"]
q_origin = child_function(pointing["quartic_ordinate"])
p1229_support = reduce_qq(pointing["old_base_coordinate"])
quartic_coefficients = [KU(U(coefficient)) for coefficient in TUP(quartic).list()]
translated = []
from math import comb
for new_degree in range(5):
    translated.append(sum(
        quartic_coefficients[old_degree]*F(comb(old_degree, new_degree))
        * KU(p1229_support)**(old_degree-new_degree)
        for old_degree in range(new_degree, 5)
    ))
e, d, c, b, a = translated
assert e == q_origin**2
a1 = d/q_origin
a2 = c-d**2/(4*q_origin**2)
a3 = 2*q_origin*b
b2 = a1**2+4*a2


def invert_child_section(X, Y):
    X = KU(X)
    Y = KU(Y)
    x_general = X/9-b2/12
    y_general = Y/27-(a1*x_general+a3)/2
    local_t = 2*q_origin*(x_general+a2)/y_general
    old_t = KU(p1229_support)+local_t
    ordinate = (x_general*local_t**2-d*local_t)/(2*q_origin)-q_origin
    assert ordinate**2 == evaluate_bivariate(quartic, old_t)
    AA = evaluate_parent_poly(AA0, old_t)+KU(u)*evaluate_parent_poly(AA1, old_t)
    BB = evaluate_parent_poly(BB0, old_t)+KU(u)*evaluate_parent_poly(BB1, old_t)
    ZZ = evaluate_parent_poly(Z, old_t)
    slope = -AA/(BB*ZZ)
    square = evaluate_bivariate(square_factor, old_t)
    hx = compose_parent(H[0], old_t)
    hy = compose_parent(H[1], old_t)
    old_x = (ordinate*square/BB**2-hx+slope**2)/2
    old_y = slope*(old_x-hx)-hy
    assert old_y**2 == compose_parent(parent_A, old_t)*old_x + old_x**3 + compose_parent(parent_B, old_t)
    restrictions = [
        evaluate_parent_poly(AAi, old_t)+evaluate_parent_poly(BBi, old_t)*ZZ*slope
        for AAi, BBi in resolved_pairs
    ]
    assert KU(u) == -restrictions[0]/restrictions[1]
    return old_t, ordinate, old_x, old_y


def rational_degree(value):
    return max(value.numerator().degree(), value.denominator().degree())


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def jacobian_rank(X, Y):
    dx = -3*X**2-child_A
    dy = 2*Y
    columns = [padded(u**power*dx, 13) for power in range(5)]
    columns += [padded(u**power*dy, 13) for power in range(7)]
    return int(matrix(F, columns).transpose().rank())


# Prove the inverse formulas on the exact degree-one Q3 branch.
q3_section = q3_exact["section"]
Q3_child = (
    child_poly(q3_section["x_coefficients_low_to_high"]),
    child_poly(q3_section["y_coefficients_low_to_high"]),
)
q3_old_t, unused_q3_w, q3_old_x, q3_old_y = invert_child_section(*Q3_child)
assert rational_degree(q3_old_t) == 1
assert q3_old_x == compose_parent(target_words["Q3"][0], q3_old_t)
assert q3_old_y == compose_parent(target_words["Q3"][1], q3_old_t)

# Reconstruct the complete exact physical P.O=0 shell and its q4 MW words.
def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


U2 = matrix(ZZ, ((0, 1), (1, 0)))
frame = load_matrix(FRAME)
gram = block_diagonal_matrix(U2, -frame)
basis_in_parent = matrix(ZZ, marking["basis_in_source"])
parent_in_basis = matrix(ZZ, marking["source_in_basis"])
parent_fibre = vector(ZZ, [1, 0] + [0]*17)*parent_in_basis
explicit = {
    name: vector(ZZ, value)
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
known_sections = [
    (name, vector(ZZ, curve[-13:]))
    for name, curve in explicit.items()
    if curve[1] == 1 and any(curve[-13:])
]
known_matrix = matrix(ZZ, [list(value) for unused, value in known_sections])
smith, unused_left, smith_right = known_matrix.smith_form()
known_rank = int(known_matrix.rank())
smith_diagonal = [abs(ZZ(smith[index, index])) for index in range(known_rank)]


def quotient_key(value):
    transformed = vector(ZZ, value)*smith_right
    return tuple(
        [int(transformed[index] % smith_diagonal[index]) for index in range(known_rank)]
        + [int(transformed[index]) for index in range(known_rank, 13)]
    )


def quotient_subtract(left, right):
    return (
        tuple((left[index]-right[index]) % smith_diagonal[index] for index in range(known_rank))
        + tuple(left[index]-right[index] for index in range(known_rank, 13))
    )


embedding = next(
    item for item in height["marked_embedding_enumeration"]["embeddings"]
    if item["embedding_index"] == 15
)
B_rows = matrix(QQ, embedding["rows_B0_through_B7_in_marked_MW9"])
H_tail = vector(QQ, height["marked_embedding_enumeration"]["q8_horizontal_marked_MW9_tail"])
equation_rows = B_rows.stack(matrix(QQ, [H_tail]))
assert equation_rows.det() == -3
c8_relation = vector(ZZ, [-2, -3, -4, 3, -2, 2, -1, -2])
c8_basis_coordinates = vector(QQ, [QQ(value)/3 for value in c8_relation])


def equation_word(parent_curve):
    target_tail = vector(QQ, parent_curve[-9:])
    rational_word = equation_rows.transpose().solve_right(target_tail)
    assert rational_word[-1] in ZZ
    choices = []
    for c8_coefficient in range(3):
        integral_B = vector(QQ, rational_word[:8])-c8_coefficient*c8_basis_coordinates
        if all(value in ZZ for value in integral_B):
            choices.append((c8_coefficient, vector(ZZ, integral_B), ZZ(rational_word[-1])))
    assert len(choices) == 1
    c8_coefficient, integral_B, h_coefficient = choices[0]
    return {
        "B0_through_B7": [int(value) for value in integral_B],
        "C8opposite": int(c8_coefficient),
        "H": int(h_coefficient),
    }


physical_classes = []
norm_shell = pari(frame).qfminim(4)
half = [vector(ZZ, column) for column in matrix(ZZ, norm_shell[2]).columns()]
seen_mw = set()
for tail in half+[-value for value in half]:
    if tail*frame*tail != 4:
        continue
    component_pairings = vector(ZZ, tail)*frame[:, :4]
    affine_pairings = vector(ZZ, [1]*4)-component_pairings
    if min(tuple(component_pairings)+tuple(affine_pairings)) < 0:
        continue
    section = vector(ZZ, [1, 1]+list(tail))
    parent_curve = section*basis_in_parent
    mw = vector(ZZ, tail[-13:])
    assert tuple(mw) not in seen_mw
    seen_mw.add(tuple(mw))
    word = equation_word(parent_curve)
    word_sha = hashlib.sha256(json.dumps(word, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    physical_classes.append({
        "class_index": len(physical_classes),
        "current_4A1_section": section,
        "current_4A1_mw": mw,
        "current_component_pairings": component_pairings,
        "current_affine_pairings": affine_pairings,
        "q4o164_parent_curve": parent_curve,
        "q4o164_parent_degree": int(section*gram*parent_fibre),
        "q4o164_parent_a_minus_b": int(parent_curve[0]-parent_curve[1]),
        "q4_equation_word": word,
        "q4_equation_word_sha256": word_sha,
    })
assert len(physical_classes) == 938


def fibre_add(left, right, coefficient_A):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2 or not y1:
            return None
        slope = (3*x1**2+coefficient_A)/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return x3, slope*(x1-x3)-y1


def fibre_mul(coefficient, point, coefficient_A):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return fibre_mul(-coefficient, None if point is None else (point[0], -point[1]), coefficient_A)
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = fibre_add(answer, addend, coefficient_A)
        addend = fibre_add(addend, addend, coefficient_A)
        coefficient >>= 1
    return answer


parent_generators = parent_basis+[C8opposite, H]


def evaluate_parent_F(function, t0):
    function = PK(function)
    denominator = function.denominator()(t0)
    if not denominator:
        raise ZeroDivisionError
    return F(function.numerator()(t0))/F(denominator)


safe_t0 = []
generator_fibres = {}
for integer_t0 in range(int(prime)):
    t0 = F(integer_t0)
    if 4*parent_A(t0)**3+27*parent_B(t0)**2 == 0:
        continue
    try:
        points = [
            (evaluate_parent_F(point[0], t0), evaluate_parent_F(point[1], t0))
            for point in parent_generators
        ]
    except ZeroDivisionError:
        continue
    safe_t0.append(integer_t0)
    generator_fibres[integer_t0] = points


def expected_point(item, integer_t0):
    points = generator_fibres[integer_t0]
    coefficient_A = F(parent_A(F(integer_t0)))
    word = item["q4_equation_word"]
    answer = None
    for coefficient, point in zip(word["B0_through_B7"], points[:8]):
        answer = fibre_add(answer, fibre_mul(coefficient, point, coefficient_A), coefficient_A)
    answer = fibre_add(answer, fibre_mul(word["C8opposite"], points[8], coefficient_A), coefficient_A)
    answer = fibre_add(answer, fibre_mul(word["H"], points[9], coefficient_A), coefficient_A)
    return None if answer is None else (int(answer[0]), int(answer[1]))


from math import lcm
extensions = {}


def splitting_field(factor_degrees):
    extension_degree = 1
    for degree in factor_degrees:
        extension_degree = lcm(extension_degree, degree)
    if extension_degree not in extensions:
        E = GF(prime**extension_degree, f"z{extension_degree}")
        extensions[extension_degree] = (
            E, PolynomialRing(E, "v"), {E(F(value)): value for value in range(int(prime))}
        )
    return extensions[extension_degree]


def evaluate_over_extension(function, value, E):
    function = KU(function)
    def evaluate(poly):
        answer = E.zero()
        for coefficient in reversed(U(poly).list()):
            answer = answer*value+E(coefficient)
        return answer
    denominator = evaluate(function.denominator())
    if not denominator:
        raise ZeroDivisionError
    return evaluate(function.numerator())/denominator


def trace_signature(old_t, old_x, old_y, degree):
    trials = []
    for integer_t0 in safe_t0:
        if len(trials) >= args.fibres:
            break
        t0 = F(integer_t0)
        polynomial = U(old_t.numerator()-t0*old_t.denominator())
        if polynomial.degree() != degree or not polynomial.is_squarefree():
            continue
        factor_degrees = [int(factor.degree()) for factor, unused in polynomial.factor()]
        E, EU, descent = splitting_field(factor_degrees)
        roots = EU([E(value) for value in polynomial.list()]).roots(multiplicities=False)
        if len(roots) != degree:
            continue
        try:
            points = [
                (evaluate_over_extension(old_x, root, E), evaluate_over_extension(old_y, root, E))
                for root in roots
            ]
        except ZeroDivisionError:
            continue
        coefficient_A = E(parent_A(t0))
        coefficient_B = E(parent_B(t0))
        assert all(y**2 == x**3+coefficient_A*x+coefficient_B for x, y in points)
        trace = None
        for point in points:
            trace = fibre_add(trace, point, coefficient_A)
        if trace is None:
            descended = None
        else:
            assert trace[0] in descent and trace[1] in descent
            descended = (int(descent[trace[0]]), int(descent[trace[1]]))
        trials.append({"old_base_value": integer_t0, "root_count": len(roots), "trace": descended})
    return trials


expected_cache = {}
def expected_index(trials):
    supports = tuple(trial["old_base_value"] for trial in trials)
    if supports not in expected_cache:
        table = {}
        for item in physical_classes:
            signature = tuple(expected_point(item, t0) for t0 in supports)
            table.setdefault(signature, []).append(item["class_index"])
        expected_cache[supports] = table
    signature = tuple(trial["trace"] for trial in trials)
    return expected_cache[supports].get(signature, [])


def expected_equation_profile(item):
    profile = list(map(int, item["current_component_pairings"]))
    if profile == [0, 0, 0, 0]:
        return [0, 0, 0, 0]
    if profile == [0, 1, 0, 0]:
        return [1, 0, 0, 0]
    return None


classifications = []
realizations = {}
for shell_index, record in enumerate(shell["all_records"]):
    X = U(record["x_coefficients_low_to_high"])
    Y = U(record["y_coefficients_low_to_high"])
    old_t, unused_W, old_x, old_y = invert_child_section(X, Y)
    degree = rational_degree(old_t)
    assert degree == record["inverse_parent_degree"]
    if degree == 0:
        trials = []
        trace_matches = []
    else:
        trials = trace_signature(old_t, old_x, old_y, degree)
        assert len(trials) == args.fibres
        trace_matches = expected_index(trials)
    compatible = [
        index for index in trace_matches
        if expected_equation_profile(physical_classes[index]) == record["equation_component_profile"]
        and physical_classes[index]["q4o164_parent_degree"] == degree
    ]
    row = {
        "shell_index": shell_index,
        "equation_component_profile": record["equation_component_profile"],
        "inverse_parent_degree": int(degree),
        "constant_parent_base_modp": int(F(old_t)) if degree == 0 else None,
        "ordinary_coefficient_jacobian_rank": record["ordinary_coefficient_jacobian_rank"],
        "ordinary_fibre_trials": trials,
        "trace_matching_lattice_class_indices": trace_matches,
        "profile_compatible_lattice_class_indices": compatible,
        "classification": (
            "NONDOMINANT_PARENT_MAP" if degree == 0 else
            "UNIQUE_PHYSICAL_LATTICE_CLASS" if len(compatible) == 1 else
            "AMBIGUOUS_PHYSICAL_LATTICE_CLASS" if compatible else
            "NO_PHYSICAL_LATTICE_CLASS"
        ),
    }
    classifications.append(row)
    if len(compatible) == 1:
        realizations.setdefault(compatible[0], []).append(shell_index)


# Repeat the exact four-section quotient search first on uniquely identified
# classes, then on the larger set-valued classification supplied by the same
# exact invariants.  The latter is explicitly reported as conditional.
target_record = next(
    item for item in frontier["targets"]
    if item["fibre_fingerprint_sha256"] == "d676cab5918a08add2f743081cf76932d02adc9e18ab9c808dc05760fe0157dd"
)
target_mw = vector(ZZ, target_record["mw"])
target_key = quotient_key(target_mw)
possible_realizations = {}
for row in classifications:
    for class_index in row["profile_compatible_lattice_class_indices"]:
        possible_realizations.setdefault(class_index, []).append(row["shell_index"])


def word_order(items, pole_first):
    degrees = [item["q4o164_parent_degree"] for item in items]
    poles = [item["q4o164_parent_a_minus_b"] for item in items]
    if pole_first:
        return (sum(poles), sum(degrees), max(poles), max(degrees), tuple(item["class_index"] for item in items))
    return (sum(degrees), sum(poles), max(degrees), max(poles), tuple(item["class_index"] for item in items))


def shell_assignment(items, realization_map):
    class_indices = []
    for item in items:
        if item["class_index"] not in class_indices:
            class_indices.append(item["class_index"])
    choices = {
        class_index: sorted(realization_map[class_index], key=lambda index: (
            -classifications[index]["ordinary_coefficient_jacobian_rank"], index
        ))
        for class_index in class_indices
    }
    selected = {}
    def recurse(position, used):
        if position == len(class_indices):
            return True
        class_index = class_indices[position]
        for shell_index in choices[class_index]:
            if shell_index in used:
                continue
            selected[class_index] = shell_index
            if recurse(position+1, used | {shell_index}):
                return True
        selected.pop(class_index, None)
        return False
    if not recurse(0, set()):
        return None
    return [selected[item["class_index"]] for item in items]


def best_four(realization_map, pole_first):
    pool = [physical_classes[index] for index in sorted(realization_map)]
    pairs = {}
    for left in range(len(pool)):
        for right in range(left, len(pool)):
            key = quotient_key(pool[left]["current_4A1_mw"]+pool[right]["current_4A1_mw"])
            indices = (left, right)
            ordering = word_order([pool[left], pool[right]], pole_first)
            pairs.setdefault(key, []).append((ordering, indices))
    for key in pairs:
        pairs[key] = sorted(pairs[key])[:20]
    words = []
    for left_key, left_pairs in pairs.items():
        needed = quotient_subtract(target_key, left_key)
        if needed not in pairs or left_key > needed:
            continue
        for unused_left_order, left_pair in left_pairs:
            for unused_right_order, right_pair in pairs[needed]:
                indices = left_pair+right_pair
                items = [pool[index] for index in indices]
                assignment = shell_assignment(items, realization_map)
                if assignment is not None:
                    words.append((word_order(items, pole_first), indices, assignment))
    if not words:
        return None
    unused_order, indices, assignment = min(words)
    items = [pool[index] for index in indices]
    residual = target_mw-sum((item["current_4A1_mw"] for item in items), vector(ZZ, 13))
    coefficients = vector(QQ, known_matrix.transpose().solve_right(residual.column()).column(0))
    assert all(value in ZZ for value in coefficients)
    coefficients = vector(ZZ, coefficients)
    assert sum((item["current_4A1_mw"] for item in items), vector(ZZ, 13))+coefficients*known_matrix == target_mw
    return {
        "lattice_class_indices": [item["class_index"] for item in items],
        "representative_shell_indices": assignment,
        "all_realizing_shell_indices": [realization_map[item["class_index"]] for item in items],
        "q4o164_parent_degrees": [item["q4o164_parent_degree"] for item in items],
        "q4o164_parent_a_minus_b": [item["q4o164_parent_a_minus_b"] for item in items],
        "q4o164_parent_degree_sum": sum(item["q4o164_parent_degree"] for item in items),
        "q4o164_parent_a_minus_b_sum": sum(item["q4o164_parent_a_minus_b"] for item in items),
        "known_section_correction": [
            {"name": known_sections[index][0], "coefficient": int(value)}
            for index, value in enumerate(coefficients) if value
        ],
        "exact_mw_identity_pass": True,
    }


best_degree = best_four(realizations, False)
best_pole = best_four(realizations, True)
possible_best_degree = best_four(possible_realizations, False)
possible_best_pole = best_four(possible_realizations, True)


def polynomial_section_intersection(left_record, right_record):
    x_difference = U(left_record["x_coefficients_low_to_high"])-U(right_record["x_coefficients_low_to_high"])
    y_difference = U(left_record["y_coefficients_low_to_high"])-U(right_record["y_coefficients_low_to_high"])
    if not x_difference and not y_difference:
        return None
    if not x_difference:
        common = y_difference.monic()
    elif not y_difference:
        common = x_difference.monic()
    else:
        common = x_difference.gcd(y_difference).monic()
    # Reject a raw Weierstrass collision at a singular finite fibre: it needs
    # a resolved local chart rather than this smooth-section intersection.
    discriminant = 4*child_A**3+27*child_B**2
    if common.gcd(discriminant).degree() > 0:
        return None
    infinity_orders = []
    if x_difference:
        infinity_orders.append(4-int(x_difference.degree()))
    if y_difference:
        infinity_orders.append(6-int(y_difference.degree()))
    infinity_intersection = max(0, min(infinity_orders))
    return int(common.degree()+infinity_intersection)


intersection_disambiguation = []
if possible_best_pole is not None:
    selected_shells = sorted(set(sum(possible_best_pole["all_realizing_shell_indices"], [])))
    anchors = [
        row for row in classifications
        if len(row["profile_compatible_lattice_class_indices"]) == 1
        and row["equation_component_profile"] == [0, 0, 0, 0]
        and row["shell_index"] not in selected_shells
    ]
    for shell_index in selected_shells:
        row = classifications[shell_index]
        alternatives = row["profile_compatible_lattice_class_indices"]
        if len(alternatives) <= 1:
            continue
        comparisons = []
        surviving = set(alternatives)
        for anchor in anchors:
            observed = polynomial_section_intersection(
                shell["all_records"][shell_index], shell["all_records"][anchor["shell_index"]]
            )
            if observed is None:
                continue
            anchor_class = physical_classes[anchor["profile_compatible_lattice_class_indices"][0]]
            expected = {
                index: int(
                    physical_classes[index]["current_4A1_section"]
                    * gram * anchor_class["current_4A1_section"]
                )
                for index in alternatives
            }
            comparisons.append({
                "anchor_shell_index": anchor["shell_index"],
                "anchor_lattice_class_index": anchor_class["class_index"],
                "observed_smooth_section_intersection": observed,
                "expected_lattice_intersections": {str(index): value for index, value in expected.items()},
            })
            surviving &= {index for index, value in expected.items() if value == observed}
            if len(surviving) <= 1:
                break
        intersection_disambiguation.append({
            "shell_index": shell_index,
            "initial_class_alternatives": alternatives,
            "surviving_class_alternatives": sorted(surviving),
            "comparisons": comparisons,
            "resolved_uniquely": len(surviving) == 1,
        })

q3_lattice_class = next(
    item for item in physical_classes
    if list(item["current_4A1_mw"]) == [6, 1, 0, 0, 1, -3, -4, 1, 1, 0, -1, 1, 0]
)
q2_lattice_class = physical_classes[500]
q3_record = {
    "x_coefficients_low_to_high": list(map(int, Q3_child[0].list())),
    "y_coefficients_low_to_high": list(map(int, Q3_child[1].list())),
}
targeted_anchor_disambiguation = []
for shell_index in sorted(set(
    index
    for group in (possible_best_pole["all_realizing_shell_indices"] if possible_best_pole else [])
    for index in group
)):
    alternatives = classifications[shell_index]["profile_compatible_lattice_class_indices"]
    if len(alternatives) <= 1:
        continue
    comparisons = []
    for anchor_name, anchor_record, anchor_class in (
        ("exact_Q3", q3_record, q3_lattice_class),
        ("named_Q2_shell116", shell["all_records"][116], q2_lattice_class),
    ):
        observed = polynomial_section_intersection(shell["all_records"][shell_index], anchor_record)
        comparisons.append({
            "anchor": anchor_name,
            "anchor_lattice_class_index": anchor_class["class_index"],
            "observed_smooth_section_intersection": observed,
            "expected_lattice_intersections": {
                str(index): int(
                    physical_classes[index]["current_4A1_section"]
                    * gram * anchor_class["current_4A1_section"]
                ) for index in alternatives
            },
        })
    targeted_anchor_disambiguation.append({
        "shell_index": shell_index,
        "class_alternatives": alternatives,
        "Q2_Q3_comparisons": comparisons,
    })

# Feed the smooth-anchor resolutions back into the possible realization map
# and rerun the bounded word search.  Empty survivors reject that modular
# record as a lift of any of its initial generic alternatives.
refined_possible_realizations = {
    class_index: list(indices) for class_index, indices in possible_realizations.items()
}
for result in intersection_disambiguation:
    shell_index = result["shell_index"]
    survivors = set(result["surviving_class_alternatives"])
    for class_index in result["initial_class_alternatives"]:
        if class_index in survivors:
            continue
        if class_index in refined_possible_realizations:
            refined_possible_realizations[class_index] = [
                index for index in refined_possible_realizations[class_index]
                if index != shell_index
            ]
            if not refined_possible_realizations[class_index]:
                del refined_possible_realizations[class_index]
refined_possible_best_degree = best_four(refined_possible_realizations, False)
refined_possible_best_pole = best_four(refined_possible_realizations, True)


def disambiguate_word_by_smooth_anchors(word, realization_map):
    if word is None:
        return []
    selected_shells = sorted(set(
        index for group in word["all_realizing_shell_indices"] for index in group
    ))
    anchors = [
        row for row in classifications
        if len(row["profile_compatible_lattice_class_indices"]) == 1
        and row["equation_component_profile"] == [0, 0, 0, 0]
        and row["shell_index"] not in selected_shells
    ]
    results = []
    for shell_index in selected_shells:
        alternatives = sorted(
            class_index for class_index, indices in realization_map.items()
            if shell_index in indices
        )
        if len(alternatives) <= 1:
            continue
        surviving = set(alternatives)
        comparisons = []
        for anchor in anchors:
            observed = polynomial_section_intersection(
                shell["all_records"][shell_index], shell["all_records"][anchor["shell_index"]]
            )
            if observed is None:
                continue
            anchor_class = physical_classes[anchor["profile_compatible_lattice_class_indices"][0]]
            expected = {
                index: int(
                    physical_classes[index]["current_4A1_section"]
                    * gram * anchor_class["current_4A1_section"]
                ) for index in alternatives
            }
            comparisons.append({
                "anchor_shell_index": anchor["shell_index"],
                "anchor_lattice_class_index": anchor_class["class_index"],
                "observed_smooth_section_intersection": observed,
                "expected_lattice_intersections": {str(index): value for index, value in expected.items()},
            })
            surviving &= {index for index, value in expected.items() if value == observed}
            if len(surviving) <= 1:
                break
        results.append({
            "shell_index": shell_index,
            "initial_class_alternatives": alternatives,
            "surviving_class_alternatives": sorted(surviving),
            "comparisons": comparisons,
            "resolved_uniquely": len(surviving) == 1,
        })
    return results


refined_word_disambiguation = disambiguate_word_by_smooth_anchors(
    refined_possible_best_pole, refined_possible_realizations
)
refined_Q2_Q3_disambiguation = []
for result in refined_word_disambiguation:
    shell_index = result["shell_index"]
    alternatives = result["initial_class_alternatives"]
    comparisons = []
    for anchor_name, anchor_record, anchor_class in (
        ("exact_Q3", q3_record, q3_lattice_class),
        ("named_Q2_shell116", shell["all_records"][116], q2_lattice_class),
    ):
        observed = polynomial_section_intersection(shell["all_records"][shell_index], anchor_record)
        comparisons.append({
            "anchor": anchor_name,
            "anchor_lattice_class_index": anchor_class["class_index"],
            "observed_smooth_section_intersection": observed,
            "expected_lattice_intersections": {
                str(index): int(
                    physical_classes[index]["current_4A1_section"]
                    * gram * anchor_class["current_4A1_section"]
                ) for index in alternatives
            },
        })
    refined_Q2_Q3_disambiguation.append({
        "shell_index": shell_index,
        "class_alternatives": alternatives,
        "Q2_Q3_comparisons": comparisons,
    })

q2_refined_possible_realizations = {
    class_index: list(indices) for class_index, indices in refined_possible_realizations.items()
}
for result in refined_Q2_Q3_disambiguation:
    shell_index = result["shell_index"]
    q2_comparison = next(
        comparison for comparison in result["Q2_Q3_comparisons"]
        if comparison["anchor"] == "named_Q2_shell116"
    )
    observed = q2_comparison["observed_smooth_section_intersection"]
    if observed is None:
        continue
    survivors = {
        int(index) for index, expected in q2_comparison["expected_lattice_intersections"].items()
        if expected == observed
    }
    for class_index in result["class_alternatives"]:
        if class_index in survivors or class_index not in q2_refined_possible_realizations:
            continue
        q2_refined_possible_realizations[class_index] = [
            index for index in q2_refined_possible_realizations[class_index]
            if index != shell_index
        ]
        if not q2_refined_possible_realizations[class_index]:
            del q2_refined_possible_realizations[class_index]
q2_refined_best_degree = best_four(q2_refined_possible_realizations, False)
q2_refined_best_pole = best_four(q2_refined_possible_realizations, True)


def integer_list(value):
    return [int(item) for item in value]


lattice_payload = []
for item in physical_classes:
    lattice_payload.append({
        "class_index": item["class_index"],
        "current_4A1_mw": integer_list(item["current_4A1_mw"]),
        "current_component_pairings": integer_list(item["current_component_pairings"]),
        "q4o164_parent_curve": integer_list(item["q4o164_parent_curve"]),
        "q4o164_parent_degree": item["q4o164_parent_degree"],
        "q4o164_parent_a_minus_b": item["q4o164_parent_a_minus_b"],
        "q4_equation_word": item["q4_equation_word"],
        "q4_equation_word_sha256": item["q4_equation_word_sha256"],
        "realizing_shell_indices": realizations.get(item["class_index"], []),
        "possible_realizing_shell_indices": possible_realizations.get(item["class_index"], []),
    })
word_hash = hashlib.sha256(json.dumps(
    [item["q4_equation_word"] for item in lattice_payload],
    sort_keys=True, separators=(",", ":")
).encode()).hexdigest()


def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


counts = {}
for row in classifications:
    counts[row["classification"]] = counts.get(row["classification"], 0)+1
payload = {
    "schema": "elkies-k3.h92-q12o5867-p0-shell-lattice-classification-modp.v1",
    "status": (
        "PASS_MODP_Q12O5867_REALIZED_SHELL_FOUR_WORD_FOUND"
        if best_pole is not None else
        "PARTIAL_MODP_Q12O5867_AMBIGUOUS_SHELL_FOUR_WORD_FOUND"
        if q2_refined_best_pole is not None else
        "REJECTED_MODP_Q12O5867_COMPATIBLE_SHELL_NO_FOUR_WORD"
    ),
    "prime": int(prime),
    "ordinary_fibres_per_positive_degree_record": args.fibres,
    "exact_Q3_inverse_regression": True,
    "lattice_shell": {
        "physical_P_dot_O_zero_class_count": len(physical_classes),
        "equation_basis_determinant": int(equation_rows.det()),
        "expected_q4_word_list_sha256": word_hash,
        "classes": lattice_payload,
    },
    "polynomial_shell": {
        "signed_record_count": len(classifications),
        "classification_counts": counts,
        "uniquely_realized_lattice_class_count": len(realizations),
        "profile_compatible_possible_lattice_class_count": len(possible_realizations),
        "records": classifications,
    },
    "q12o5867_realized_four_section_search": {
        "target_fibre_fingerprint_sha256": target_record["fibre_fingerprint_sha256"],
        "target_mw": integer_list(target_mw),
        "known_section_subgroup_rank": known_rank,
        "known_section_subgroup_smith_diagonal": integer_list(smith_diagonal),
        "realized_class_count": len(realizations),
        "best_parent_degree_first_word": best_degree,
        "best_parent_a_minus_b_first_word": best_pole,
        "conditional_possible_class_count": len(possible_realizations),
        "conditional_best_parent_degree_first_word": possible_best_degree,
        "conditional_best_parent_a_minus_b_first_word": possible_best_pole,
        "selected_word_intersection_disambiguation": intersection_disambiguation,
        "selected_word_Q2_Q3_intersection_disambiguation": targeted_anchor_disambiguation,
        "refined_conditional_possible_class_count": len(refined_possible_realizations),
        "refined_conditional_best_parent_degree_first_word": refined_possible_best_degree,
        "refined_conditional_best_parent_a_minus_b_first_word": refined_possible_best_pole,
        "refined_selected_word_intersection_disambiguation": refined_word_disambiguation,
        "refined_selected_word_Q2_Q3_intersection_disambiguation": refined_Q2_Q3_disambiguation,
        "Q2_intersection_refined_possible_class_count": len(q2_refined_possible_realizations),
        "Q2_intersection_refined_best_parent_degree_first_word": q2_refined_best_degree,
        "Q2_intersection_refined_best_parent_a_minus_b_first_word": q2_refined_best_pole,
    },
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "QQ_lift_attempted": False,
        "construction": "exact pointed-quartic/chord inverse, finite-field Abel trace, exact marked MW quotient meet-in-the-middle",
        "splitting_extension_degrees_used": sorted(extensions),
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "This is an exact bounded classification over the displayed finite field. Unique labels require both the Abel signature and the transported component profile. The four-word search is exact modulo the displayed known subgroup among uniquely realized modular classes. It is not a QQ lift or endpoint equation certificate."
    ),
    "inputs": {
        "paths": [display_path(path) for path in INPUTS_BASE+(SHELL,)],
        "sha256": {display_path(path): sha256(path) for path in INPUTS_BASE+(SHELL,)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867LATTICECLASS|prime={}|records={}|classes={}|unique_realized={}|"
    "classification={}|four_word={}|runtime={:.3f}|status={}|output={}".format(
        prime, len(classifications), len(physical_classes), len(realizations), counts,
        int((best_pole or q2_refined_best_pole) is not None), payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
