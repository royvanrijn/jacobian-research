#!/usr/bin/env python3
"""Replay the exact algebra in THREE_MERIDIAN_FRONTIER_AUDIT.md.

Run with Python 3 and no third-party packages:
    python3 research/meridional-rank/representations/verify_three_meridian_candidate.py

Candidate: closure(sigma1^4 (sigma1 sigma2 sigma3)^5) = T(4,5;2,2).
The default replay treats the presentation and meridian as software inputs.
The optional --verify-exterior flag, run with Sage and SnapPy 3.3.2,
reconstructs the literal braid exterior, compares its presentation/meridian
and verifies hyperbolic gluing equations with interval arithmetic. The
diagram/cell-group and peripheral identification remains a SOFTWARE IMPORT.
Bridge number and the (1,1)-presentation are separate literature inputs.

The checks prove exact identities for the displayed presentation and the
trefoil permutations and the odd-cone-order countercheck. They exhibit no generating meridian triple, do not
replay the monodromy/coset searches, and prove no meridional-rank bound.
Five finite quotients replay the negative part of the bounded 595-triple
screen. All words use signed generator indices; permutation products act
rightmost first. No external inputs or files are read except this
file itself for the reported SHA-256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import platform


BETA = (1, 1, 1, 1) + (1, 2, 3) * 5
RELATOR_AB = "aabbbaabbbbaBBBAABBBabbbb"
MU_AB = (-2, -1)  # BA; m=ab=mu^-1.
# Indices 1,2 denote x,m in XM words and a,b in AB words.
AB_TO_XM = {1: [1, 2, 2, 2, 2, 2], 2: [-2, -2, -2, -2, -2, -1, 2]}
XM_TO_AB = {
    1: [1, -2, -1, -2, -1, -2, -1, -2, -1, -2, -1],
    2: [1, 2],
}
MAGNUS = [
    (0, 1), (1, -1), (-3, -1), (-2, 1), (-1, -1),
    (-5, -1), (-9, -1), (-8, 1), (-4, 1), (0, 1),
    (-1, -1), (-2, 1), (2, 1), (6, 1), (7, -1),
    (3, -1), (-1, -1),
]
# x_j=m^j*x*m^-j; e_i=x_(i-10), 1<=i<=16; phi(c)=m*c*m^-1.
PHI = {i: [i + 1] for i in range(1, 16)}
PHI[16] = [-13, -9, 10, -11, -7, 8, -9, -5, -1, 2, 6, 10, -9, 8, 12, 16]
PHI_INVERSE = {i: [i - 1] for i in range(2, 17)}
PHI_INVERSE[1] = [1, 5, 9, -8, 7, 11, 15, -16, -12, -8, 9, -10, -6, 7, -8, -4]

# Five literal finite quotients suffice for the retained 286 negative witnesses.
# (a image, b image, full image order); each satisfies the original relator.
FINITE_MODELS = (
    ((1, 0, 3, 4, 5, 2), (0, 2, 3, 5, 1, 4), 360),
    ((1, 0, 3, 4, 5, 2), (1, 2, 4, 0, 3, 5), 360),
    ((0, 1, 3, 4, 2), (2, 3, 1, 4, 0), 60),
    ((1, 2, 3, 4, 0), (0, 3, 4, 1, 2), 60),
    ((0, 2, 3, 4, 1), (1, 4, 3, 0, 2), 20),
)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def reduce_word(word):
    result = []
    for letter in word:
        if result and result[-1] == -letter:
            result.pop()
        else:
            result.append(letter)
    return result


def inverse_word(word):
    return [-letter for letter in reversed(word)]


def substitute(word, images):
    return reduce_word(
        letter
        for index in word
        for letter in (
            images[index] if index > 0 else inverse_word(images[-index])
        )
    )


def word_equal(left, right, label):
    require(reduce_word(left) == reduce_word(right), label)


def conjugate_x(index, sign=1):
    stable = [2 if index > 0 else -2] * abs(index)
    return stable + [sign] + inverse_word(stable)


def check_inverse_maps(phi, psi):
    for index in range(1, 17):
        word_equal(substitute(phi[index], psi), [index], "psi(phi(e_%d))" % index)
        word_equal(substitute(psi[index], phi), [index], "phi(psi(e_%d))" % index)


def action(word, exponent):
    for _ in range(abs(exponent)):
        word = substitute(word, PHI if exponent > 0 else PHI_INVERSE)
    return word


def multiply_nf(left, right):
    """Normal form (c,k)=c*m^k; multiplication uses phi^k."""
    return (
        reduce_word(left[0] + action(right[0], left[1])),
        left[1] + right[1],
    )


def normal_form(word, images):
    result = ([], 0)
    for letter in word:
        result = multiply_nf(result, images[letter])
    return result


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse_permutation(permutation):
    return tuple(permutation.index(index) for index in range(len(permutation)))


def power_permutation(permutation, exponent):
    result = tuple(range(len(permutation)))
    if exponent < 0:
        permutation = inverse_permutation(permutation)
    for _ in range(abs(exponent)):
        result = compose(result, permutation)
    return result


def conjugate_permutation(left, right):
    return compose(compose(left, right), inverse_permutation(left))


def subgroup_order(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = [identity]
    for element in queue:
        for generator in generators:
            product = compose(element, generator)
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return len(seen)


def rejects(callback):
    try:
        callback()
    except AssertionError:
        return True
    return False


def permutation_word(word, a, b):
    images = {1: a, 2: b, -1: inverse_permutation(a), -2: inverse_permutation(b)}
    result = tuple(range(len(a)))
    for letter in word:
        result = compose(result, images[letter])
    return result


def determinant(matrix):
    """Bareiss determinant with exact integer divisions and row pivoting."""
    matrix = [row[:] for row in matrix]
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "square matrix")
    previous, sign = 1, 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            alternate = next(
                (j for j in range(pivot_index + 1, size) if matrix[j][pivot_index]),
                None,
            )
            if alternate is None:
                return 0
            matrix[pivot_index], matrix[alternate] = matrix[alternate], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                value = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index] * matrix[pivot_index][column]
                )
                require(value % previous == 0, "exact Bareiss division")
                matrix[row][column] = value // previous
            matrix[row][pivot_index] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def exponent_vector(word):
    return [word.count(index) - word.count(-index) for index in range(1, 17)]


def check_bounded_negative_certificates(relator):
    # Enumerate ALL freely reduced conjugators of lengths 0..3. Deduplicate
    # conjugate words in the free group; no quotient equality is assumed.
    conjugators = [()]
    for length in range(1, 4):
        for word in list(conjugators):
            if len(word) != length - 1:
                continue
            for letter in (1, 2, -1, -2):
                if not word or letter != -word[-1]:
                    conjugators.append(word + (letter,))
    require(len(conjugators) == 53, "all reduced conjugators of length <=3")
    conjugate_words = {
        tuple(reduce_word(list(word) + list(MU_AB) + inverse_word(word)))
        for word in conjugators
    }
    conjugate_words.remove(MU_AB)
    require(len(conjugate_words) == 35, "35 nontrivial free conjugate words")
    trials = [(MU_AB, left, right) for left, right in combinations(sorted(conjugate_words), 2)]
    require(len(trials) == 595, "595 unordered pairs")
    rejected = set()
    for a, b, full_order in FINITE_MODELS:
        identity = tuple(range(len(a)))
        require(permutation_word(relator, a, b) == identity, "finite quotient relator")
        require(subgroup_order([a, b]) == full_order, "finite quotient full image order")
        images = {
            word: permutation_word(word, a, b)
            for triple in trials
            for word in triple
        }
        for index, trial in enumerate(trials):
            if subgroup_order([images[word] for word in trial]) < full_order:
                rejected.add(index)
    require(len(rejected) == 286, "286 finite negative generation certificates")
    return {
        "conjugators": 53, "distinct_nontrivial_conjugates": 35,
        "triples": 595, "finite_quotients": 5,
        "generation_disproved": len(rejected), "unresolved": len(trials) - len(rejected),
    }


def check_fiber_nilpotent_barrier(relator):
    # On fiber abelianization Q=phi^-1_ab. Each column of I-Q is the
    # exponent vector of e_i*phi^-1(e_i^-1).
    q_columns = [exponent_vector(PHI_INVERSE[index]) for index in range(1, 17)]
    difference = [
        [(1 if row == column else 0) - q_columns[column][row] for column in range(16)]
        for row in range(16)
    ]
    require(determinant(difference) == 1, "det(I-phi^-1_ab)=1")
    orbit_columns = [exponent_vector(action([10, -9], power)) for power in range(-9, 7)]
    orbit_matrix = [[column[row] for column in orbit_columns] for row in range(16)]
    require(orbit_matrix == difference, "sixteen seed-orbit columns equal I-Q")
    require(determinant(orbit_matrix) == 1, "unimodular seed-orbit columns")

    # Thus the phi-orbit subgroup L of d=e10*e9^-1 surjects onto F16_ab.
    # The group-theoretic argument in the audit gives L*gamma_j(F16)=F16
    # for every finite j. This does not imply L=F16.
    # A literal finite quotient detects properness of H=<m,a*m*a^-1>.
    a, b, full_order = FINITE_MODELS[1]
    require(permutation_word(relator, a, b) == tuple(range(len(a))), "A6 barrier relator")
    m = compose(a, b)
    pair = [m, conjugate_permutation(a, m)]
    require(subgroup_order([a, b]) == full_order == 360, "A6 barrier full order")
    require(subgroup_order(pair) == 8, "A6 barrier proper pair order")
    return {
        "det_I_minus_phi_inverse_ab": 1,
        "sixteen_seed_orbit_vector_determinant": 1,
        "full_quotient_order": 360, "two_meridian_subgroup_image_order": 8,
        "scope": "nilpotent-density countercheck, not an MRC counterexample",
    }


def check_odd_cone_lifting_barrier():
    """Exact controls for the written all-odd-order trefoil countercheck."""
    identity = ((1, 0), (0, 1))
    t = ((1, 1), (0, 1))
    b = ((1, 0), (-1, 1))
    v4 = ((1, 0), (4, 1))

    def matrix_multiply(left, right):
        return tuple(
            tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
            for i in range(2)
        )

    def matrix_inverse(matrix):
        a, b0 = matrix[0]
        c, d = matrix[1]
        require(a * d - b0 * c == 1, "determinant-one matrix")
        return ((d, -b0), (-c, a))

    # q=(b^-1*a)*b*(b^-1*a)^-1 and q*a=-V4.
    v = matrix_multiply(matrix_inverse(b), t)
    q = matrix_multiply(matrix_multiply(v, b), matrix_inverse(v))
    require(q == ((-1, 1), (-4, 3)), "odd-cone second meridian matrix")
    require(matrix_multiply(q, t) == tuple(tuple(-x for x in row) for row in v4),
            "Gamma0(4) generators")
    require(matrix_multiply(matrix_multiply(t, b), t)
            == matrix_multiply(matrix_multiply(b, t), b), "matrix braid relation")
    require(matrix_multiply(t, matrix_inverse(t)) == identity, "matrix inverse control")

    lines = ((1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (1, 2))

    def normalize_line(vector):
        vector = tuple(x % 4 for x in vector)
        require(any(x % 2 for x in vector), "primitive projective line")
        return min(vector, tuple(3 * x % 4 for x in vector))

    line_ids = {normalize_line(vector): i for i, vector in enumerate(lines)}
    primitive_lines = {
        normalize_line((x, y))
        for x in range(4) for y in range(4) if x % 2 or y % 2
    }
    require(set(line_ids) == primitive_lines and len(line_ids) == 6,
            "all six primitive projective lines modulo four")

    def projective_action(matrix):
        return tuple(
            line_ids[normalize_line(tuple(sum(matrix[r][k] * vector[k] for k in range(2))
                                         for r in range(2)))]
            for vector in lines
        )

    ta, tb, tq = (projective_action(matrix) for matrix in (t, b, q))
    full = {tuple(range(6))}
    queue = list(full)
    for element in queue:
        for generator in (ta, tb):
            product = compose(element, generator)
            if product not in full:
                full.add(product)
                queue.append(product)
    require(len(full) == 24 and subgroup_order([ta, tq]) == 4,
            "six-line full/stabilizer image orders")
    require(ta[0] == tq[0] == 0, "subgroup fixes the distinguished line")
    require(power_permutation(ta, 4) == tuple(range(6))
            and power_permutation(ta, 2) != tuple(range(6)), "meridian image order four")
    for residue in (1, 3):
        element = power_permutation(ta, residue)
        require(subgroup_order([conjugate_permutation(g, element) for g in full]) == 24,
                "normal closure at odd residue modulo four")

    # Properness is detected in the original trefoil's S3 quotient.
    a3, b3 = (1, 0, 2), (0, 2, 1)
    v3 = compose(inverse_permutation(b3), a3)
    q3 = conjugate_permutation(v3, b3)
    require(q3 == a3 and subgroup_order([a3, b3]) == 6
            and subgroup_order([a3, q3]) == 2, "proper two-meridian subgroup in S3")

    # Literal single braid rewrites certify Delta*a=b*Delta,
    # Delta*b=a*Delta and (ab)^3=Delta^2, hence centrality.
    def braid_rewrite(word, index, old, new):
        require(word[index:index + 3] == old and {old, new} == {"aba", "bab"},
                "valid single application of the original braid relator")
        return word[:index] + new + word[index + 3:]

    require(braid_rewrite("baba", 0, "bab", "aba") == "abaa", "Delta*a rewrite")
    require(braid_rewrite("abab", 1, "bab", "aba") == "aaba", "Delta*b rewrite")
    require(braid_rewrite("ababab", 3, "bab", "aba") == "abaaba", "central element rewrite")
    return {
        "proper_S3_full_and_subgroup_orders": [6, 2],
        "mod_four_projective_full_and_subgroup_orders": [24, 4],
        "meridian_projective_order": 4,
        "odd_residues_checked": [1, 3],
        "scope": "finite controls of the written all-odd-order proof; no MRC counterexample",
    }


def check_candidate_exterior():
    """Optional software-import replay; requires Sage with SnapPy 3.3.2."""
    import snappy

    require(snappy.__version__ == "3.3.2", "pinned SnapPy version")
    snappy.set_rand_seed(0)
    exterior = snappy.Link(braid_closure=list(BETA)).exterior()
    group = exterior.fundamental_group()
    require(group.generators() == ["a", "b"], "original presentation generators")
    require(group.relators() == [RELATOR_AB], "original presentation relator")
    require(group.meridian() == "BA", "original peripheral meridian")
    verified, _shapes = exterior.verify_hyperbolicity()
    require(verified, "interval hyperbolicity of the literal exterior")
    return {
        "snappy_version": snappy.__version__,
        "random_seed": 0,
        "triangulation_isosig": exterior.triangulation_isosig(),
        "presentation_and_meridian_match": True,
        "interval_hyperbolicity": True,
        "scope": "exterior only; diagram/peripheral software import; no closed-cover claim",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-exterior", action="store_true",
                        help="also replay SnapPy 3.3.2 imports and interval hyperbolicity in Sage")
    args = parser.parse_args()
    counts = {}
    require(len(BETA) == 19 and set(BETA) == {1, 2, 3}, "literal braid input")
    relator = [{"a": 1, "b": 2, "A": -1, "B": -2}[c] for c in RELATOR_AB]
    require(len(relator) == 25, "literal group relator")

    for index in (1, 2):
        word_equal(substitute(AB_TO_XM[index], XM_TO_AB), [index], "AB inverse")
        word_equal(substitute(XM_TO_AB[index], AB_TO_XM), [index], "XM inverse")
    counts["coordinate_inverse_identities"] = 4

    rewritten = substitute(relator, AB_TO_XM)
    require(len(rewritten) == 61, "rewritten relator length")
    word_equal(substitute(rewritten, XM_TO_AB), relator, "relator back to AB")
    height = 0
    magnus = []
    for letter in rewritten:
        if abs(letter) == 2:
            height += 1 if letter > 0 else -1
        else:
            magnus.append((height, letter))
    require(height == 0 and magnus == MAGNUS, "Magnus reconstruction")
    word_equal(
        [letter for index, sign in MAGNUS for letter in conjugate_x(index, sign)],
        rewritten,
        "expanded Magnus relator",
    )
    indices = [index for index, _ in magnus]
    require(min(indices) == -9 and max(indices) == 7, "Magnus extreme indices")
    require(indices.count(-9) == indices.count(7) == 1, "unique extreme letters")
    counts["magnus_reconstruction_and_expansion"] = 2
    counts["magnus_extrema"] = {"minimum": -9, "maximum": 7, "multiplicities": [1, 1]}

    check_inverse_maps(PHI, PHI_INVERSE)
    counts["automorphism_inverse_identities"] = 32

    # Mapping-torus generator 17 denotes m. Under substitution into <x,m>,
    # the first 15 relators vanish freely; the last is a conjugate of W^-1.
    basis_xm = {index: conjugate_x(index - 10) for index in range(1, 17)}
    mapping_to_xm = dict(basis_xm)
    mapping_to_xm[17] = [2]
    top_position = MAGNUS.index((7, -1))
    q_word = reduce_word(
        letter
        for index, sign in MAGNUS[top_position + 1:]
        for letter in conjugate_x(index, sign)
    )
    for index in range(1, 17):
        mapping_relator = [17, index, -17] + inverse_word(PHI[index])
        expected = (
            []
            if index < 16
            else q_word + inverse_word(rewritten) + inverse_word(q_word)
        )
        word_equal(
            substitute(mapping_relator, mapping_to_xm), expected,
            "mapping-torus relator transport %d" % index,
        )
    counts["mapping_torus_relator_transports"] = 16

    nf_ab = {
        1: ([10], 5), 2: ([-5], -4),
        -1: ([-5], -5), -2: ([9], 4),
    }
    nf_xm = {1: ([10], 0), -1: ([-10], 0), 2: ([], 1), -2: ([], -1)}
    require(normal_form(relator, nf_ab) == ([], 0), "original relator normal form")
    require(normal_form(MU_AB, nf_ab) == ([], -1), "meridian normal form")
    for index in range(1, 17):
        require(
            normal_form(basis_xm[index], nf_xm) == ([index], 0),
            "basis normal form %d" % index,
        )
        basis_ab = substitute(basis_xm[index], XM_TO_AB)
        require(normal_form(basis_ab, nf_ab) == ([index], 0), "AB basis normal form")
    counts["original_relator_normal_forms"] = 1
    counts["meridian_normal_forms"] = 1
    counts["basis_normal_forms"] = 32

    # For arbitrary w=c*m^q, (w*mu*w^-1)*m=c*phi^-1(c^-1).
    # These are exact finite controls, not a decision of the orbit criterion.
    # The criterion itself follows since <m,d1,d2>=L semidirect <m>, where
    # L=<phi^j(di): j in Z, i=1,2>; hence the triple generates iff L=F16.
    for index in range(1, 17):
        c, exponent = [index], index - 8
        w_inverse = (action(inverse_word(c), -exponent), -exponent)
        lhs = multiply_nf(
            multiply_nf(multiply_nf((c, exponent), ([], -1)), w_inverse),
            ([], 1),
        )
        rhs = (reduce_word(c + action(inverse_word(c), -1)), 0)
        require(lhs == rhs, "conjugate-meridian seed formula %d" % index)
    counts["seed_formula_checks"] = 16

    orders = {}
    # Zero-based tuple notation; A4 inputs are a=(123), b=(134).
    for name, a, b, full_order, triple_order in (
        ("S3", (1, 0, 2), (0, 2, 1), 6, 6),
        ("A4", (1, 2, 0, 3), (2, 1, 3, 0), 12, 3),
    ):
        require(compose(compose(a, b), a) == compose(compose(b, a), b),
                "trefoil relation in " + name)
        n = compose(
            compose(compose(power_permutation(a, 2), b), power_permutation(a, -2)),
            inverse_permutation(b),
        )
        b_prime = conjugate_permutation(n, b)
        c = conjugate_permutation(a, b_prime)
        observed = [subgroup_order([a, b]), subgroup_order([a, b_prime, c])]
        require(observed == [full_order, triple_order], "trefoil orders in " + name)
        if name == "S3":
            require(len({a, b_prime, c}) == 3, "three distinct transpositions")
            require(all(power_permutation(p, 2) == (0, 1, 2) for p in (a, b_prime, c)),
                    "meridian squares vanish in S3")
        else:
            require(n == (2, 3, 0, 1), "n=(13)(24)")
            require(a == b_prime == c, "A4 triple collapses")
        orders[name] = observed
    counts["trefoil_relation_checks"] = 2
    counts["trefoil_full_and_triple_orders"] = orders

    corrupt_phi = {index: word[:] for index, word in PHI.items()}
    corrupt_phi[16][0] *= -1
    require(rejects(lambda: check_inverse_maps(corrupt_phi, PHI_INVERSE)),
            "corrupt automorphism rejected")
    corrupt_coordinate = AB_TO_XM[1] + [1]
    require(rejects(lambda: word_equal(
        substitute(corrupt_coordinate, XM_TO_AB), [1], "corrupt coordinate")),
        "corrupt coordinate rejected")
    counts["corruption_rejection_controls"] = 2
    counts["bounded_negative_generation_certificates"] = check_bounded_negative_certificates(relator)
    counts["fiber_nilpotent_density_barrier"] = check_fiber_nilpotent_barrier(relator)
    counts["odd_cone_order_lifting_barrier"] = check_odd_cone_lifting_barrier()
    if args.verify_exterior:
        counts["candidate_exterior"] = check_candidate_exterior()

    print(json.dumps({
        "status": "pass_exact_algebra_only_no_MRC_resolution",
        "python": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "literal_braid": list(BETA),
        "counts": counts,
        "assurance": (
            "Exact standard-library replay for literal algebraic inputs. "
            "Diagram/peripheral identification remains a SnapPy 3.3.2 software "
            "import. No three-meridian generation or universal lower bound."
        ),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
