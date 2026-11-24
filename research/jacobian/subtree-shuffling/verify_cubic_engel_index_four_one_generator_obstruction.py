#!/usr/bin/env python3
"""Exact direct controls for the classical four-dimensional cubic example.

Standard library only. Symbolic polynomial identities and bounded tree-word
controls are separate from the written all-grade induction/universal argument.
"""
import argparse
import hashlib
import itertools
import json
import platform
from fractions import Fraction as Q
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = Path("research/jacobian/subtree-shuffling")
NOTE = "CUBIC_ENGEL_INDEX_FOUR_ONE_GENERATOR_OBSTRUCTION.md"
CHECKER = "verify_cubic_engel_index_four_one_generator_obstruction.py"
LEDGER = "cubic_engel_index_four_one_generator_obstruction_source_ledger.json"
RECEIPT = "cubic_engel_index_four_one_generator_obstruction_verification.json"
SOURCES = {
    "direct_proof": (NOTE, "Self-contained Jacobian, polarization, cyclic generation, all-grade words and inverse formulas."),
    "exact_checker": (CHECKER, "Literal five-monomial input and exact arithmetic controls; no asymptotic or finite-prefix inference."),
    "primary_context": ("sources/umirbaev-polarization-algebras-2014.pdf",
                        "Section 5, printed p12/PDFp16: exact Gorni-Zampieri map and one-generator polarization context; bibliography [19] on printed p13. Context only, no theorem imported."),
    "grading_convention": ("SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md",
                           "Free symmetric ternary algebra and internal-vertex grading convention only; no inverse-degree or index-three theorem imported."),
}
ASSURANCE = {
    "formal_verification": False, "external_review": False,
    "independent_replay": False, "worldwide_priority": False,
    "bounded_tree_tests_prove_all_grades": False,
    "unbounded_inverse_degree_claim": False,
}
PRIMARY = {
    "registry_id": "EXT-PAPER-UMIRBAEV-POLARIZATION-ALGEBRAS-2014",
    "author": "Ualbai Umirbaev", "title": "Polarization algebras and their relations",
    "edition": "MPIM Preprint 2014-60",
    "url": "https://archive.mpim-bonn.mpg.de/id/eprint/2035/1/preprint_2014_60.pdf",
    "locator": "Section 5, printed p12/PDFp16; bibliography [19], printed p13",
    "scope": "Exact displayed map and its attribution inspected. Classical context, not a logical premise.",
    "attribution": "G. Gorni and G. Zampieri, Yagzhev polynomial mappings: on the structure of the Taylor expansion of their local inverse, Ann. Polon. Math.64(1996),285-290.",
    "original_1996_paper_inspected": False,
}
# Coordinate indices x,y,s,t are 0,1,2,3. Monomial input slots are symmetric.
TERMS = ((0, (0, 2, 3), Q(1)), (0, (1, 2, 2), Q(-1)),
         (1, (0, 3, 3), Q(1)), (1, (1, 2, 3), Q(-1)),
         (2, (3, 3, 3), Q(1)))
PERM3 = tuple(itertools.permutations(range(3)))
ZERO = (Q(0),) * 4
BASIS = tuple(tuple(Q(i == j) for i in range(4)) for j in range(4))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    ledger = json.loads((HERE / LEDGER).read_text())
    require(set(ledger) == {"schema", "sources", "primary_context", "assurance"}, "Ledger keys")
    require(type(ledger["schema"]) is int and ledger["schema"] == 1, "Ledger schema")
    require(ledger["assurance"] == ASSURANCE and ledger["primary_context"] == PRIMARY,
            "Changed source/assurance scope")
    rows = ledger["sources"]
    require(len(rows) == len(SOURCES) and {r["role"] for r in rows} == set(SOURCES),
            "Exact source roles")
    for row in rows:
        require(set(row) == {"role", "path", "scope", "sha256"}, "Source schema")
        name, scope = SOURCES[row["role"]]
        require(row["path"] == str(PREFIX / name) and row["scope"] == scope, "Source scope")
        path = repo / row["path"]
        require(path.is_file() and sha(path) == row["sha256"], "Stale source: " + row["role"])
    return rows


def clean(p):
    return {m: Q(c) for m, c in p.items() if c}


def add(*polys):
    out = {}
    for p in polys:
        for m, c in p.items():
            out[m] = out.get(m, Q(0)) + c
    return clean(out)


def scale(p, c):
    return clean({m: v * c for m, v in p.items()})


def mul(p, q):
    out = {}
    for m, a in p.items():
        for n, b in q.items():
            e = tuple(x + y for x, y in zip(m, n))
            out[e] = out.get(e, Q(0)) + a * b
    return clean(out)


def power(p, exponent):
    out = {(0, 0, 0, 0): Q(1)}
    for _ in range(exponent):
        out = mul(out, p)
    return out


def monomial(slots):
    return {tuple(slots.count(i) for i in range(4)): Q(1)}


VARIABLES = tuple(monomial((i,)) for i in range(4))
X, Y, S, T = VARIABLES


def literal_map(terms=TERMS):
    return [add(*(scale(monomial(slots), c) for out, slots, c in terms if out == i))
            for i in range(4)]


def derivative(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            n = list(m)
            n[i] -= 1
            out[tuple(n)] = c * m[i]
    return clean(out)


def matrix_mul(a, b):
    return [[add(*(mul(a[i][k], b[k][j]) for k in range(4)))
             for j in range(4)] for i in range(4)]


def matrix_nonzero(a):
    return any(p for row in a for p in row)


def evaluate(p, vector):
    return sum((c * prod(v ** e for v, e in zip(vector, mon))
                for mon, c in p.items()), Q(0))


def substitute(p, images):
    result = {}
    for mon, c in p.items():
        term = {(0, 0, 0, 0): c}
        for e, image in zip(mon, images):
            term = mul(term, power(image, e))
        result = add(result, term)
    return result


def vadd(*vectors):
    return tuple(sum((v[i] for v in vectors), Q(0)) for i in range(4))


def vscale(vector, scalar):
    return tuple(scalar * v for v in vector)


def mu(a, b, c):
    out = [Q(0)] * 4
    for output, slots, coefficient in TERMS:
        out[output] += coefficient / 6 * sum(
            (a[slots[p[0]]] * b[slots[p[1]]] * c[slots[p[2]]] for p in PERM3), Q(0))
    return tuple(out)


def mu_inclusion_exclusion(a, b, c, homogeneous):
    def hvalue(v):
        return tuple(evaluate(p, v) for p in homogeneous)
    return vscale(vadd(hvalue(vadd(a, b, c)), vscale(hvalue(vadd(a, b)), -1),
                      vscale(hvalue(vadd(a, c)), -1), vscale(hvalue(vadd(b, c)), -1),
                      hvalue(a), hvalue(b), hvalue(c)), Q(1, 6))


def determinant_columns(columns):
    total = Q(0)
    for perm in itertools.permutations(range(4)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
        total += sign * prod(columns[j][perm[j]] for j in range(4))
    return total


def serialized_poly(p):
    return [{"exponents": list(m), "coefficient": str(c)} for m, c in sorted(p.items())]


def serialized_vector(v):
    return [str(c) for c in v]


def run_checks():
    homogeneous = literal_map()
    jacobian = [[derivative(p, i) for i in range(4)] for p in homogeneous]
    j2 = matrix_mul(jacobian, jacobian)
    j3 = matrix_mul(j2, jacobian)
    j4 = matrix_mul(j3, jacobian)
    w = add(mul(X, T), scale(mul(Y, S), -1))
    expected_cube = [[{} for _ in range(4)] for _ in range(4)]
    expected_cube[0][3] = scale(mul(mul(S, power(T, 3)), w), 3)
    expected_cube[1][3] = scale(mul(power(T, 4), w), 3)
    require(j3 == expected_cube and matrix_nonzero(j3), "Exact Jacobian cube")
    require(not matrix_nonzero(j4), "Jacobian fourth power is not zero")
    a0 = (Q(1), Q(0), Q(0), Q(1))
    cube_at_a0_et = tuple(evaluate(j3[i][3], a0) for i in range(4))
    require(cube_at_a0_et == vscale(BASIS[1], 3), "Exact-index witness")
    tensor = {}
    for i, j, k in itertools.product(range(4), repeat=3):
        value = mu(BASIS[i], BASIS[j], BASIS[k])
        require(value == mu_inclusion_exclusion(BASIS[i], BASIS[j], BASIS[k], homogeneous),
                "Polarization implementations disagree")
        tensor[i, j, k] = value
    diagonal = [add(*(scale(monomial((i, j, k)), value[out])
                      for (i, j, k), value in tensor.items())) for out in range(4)]
    require(diagonal == homogeneous, "Polarization diagonal polynomial identity")
    for out, coordinate in itertools.product(range(4), repeat=2):
        predicted = add(*(scale(mul(VARIABLES[i], VARIABLES[j]), 3 * tensor[i, j, coordinate][out])
                          for i, j in itertools.product(range(4), repeat=2)))
        require(predicted == jacobian[out][coordinate], "JH(a)=3L_a polynomial identity")
    g = (Q(0), Q(1), Q(0), Q(1))
    b = mu(g, g, g)
    c = mu(g, g, b)
    e = mu(g, b, b)
    d = mu(b, b, c)
    require(b == BASIS[2] and c == vscale(BASIS[1], Q(-1, 3))
            and e == vscale(BASIS[0], Q(-1, 3)) and d == vscale(BASIS[0], Q(1, 9)),
            "Cyclic word signs")
    det = determinant_columns((g, b, c, e))
    require(det == Q(1, 9), "Cyclic generation determinant")
    require(vadd(g, vscale(c, 3)) == BASIS[3], "Recover e_t")
    mixed = [mu(g, b, v) for v in BASIS]
    expected_mixed = [(Q(1, 6), 0, 0, 0), (0, Q(-1, 6), 0, 0),
                      (Q(-1, 3), 0, 0, 0), (0, Q(-1, 6), 0, 0)]
    require(mixed == expected_mixed, "Mixed multiplier formula")
    require(mu(g, b, c) == vscale(c, Q(-1, 6)), "Even-family eigenvalue")
    require(mu(g, b, e) == vscale(e, Q(1, 6)), "Odd-family eigenvalue")
    # Tree leaves are 0; a tuple is a ternary internal vertex.
    def grade(tree):
        return 0 if tree == 0 else 1 + sum(grade(child) for child in tree)
    def leaves(tree):
        return 1 if tree == 0 else sum(leaves(child) for child in tree)
    def tree_value(tree):
        return g if tree == 0 else mu(*(tree_value(child) for child in tree))
    bt = (0, 0, 0)
    ct, et = (0, 0, bt), (0, bt, bt)
    word_rows = []
    for k in range(13):
        for label, tree, base, sign, start in (("even", ct, c, -1, 2),
                                              ("odd", et, e, 1, 3)):
            value = tree_value(tree)
            require(value == vscale(base, Q(sign, 6) ** k) and value != ZERO,
                    "Bounded persistent tree evaluation")
            require(grade(tree) == start + 2 * k and leaves(tree) == 2 * grade(tree) + 1,
                    "Tree grade or leaves")
            word_rows.append({"family": label, "k": k, "grade": grade(tree),
                              "leaves": leaves(tree), "value": serialized_vector(value)})
        ct, et = (0, bt, ct), (0, bt, et)
    require({0, 1} | {row["grade"] for row in word_rows} == set(range(28)),
            "Displayed bounded grades have a gap")
    inverse_rows = []
    for epsilon in (-1, 1):
        forward = [add(v, scale(p, epsilon)) for v, p in zip(VARIABLES, homogeneous)]
        u = add(S, scale(power(T, 3), -epsilon))
        inverse = [add(X, scale(mul(mul(u, T), X), -epsilon), scale(mul(power(u, 2), Y), epsilon)),
                   add(Y, scale(mul(power(T, 2), X), -epsilon), scale(mul(mul(u, T), Y), epsilon)),
                   u, T]
        require([substitute(p, inverse) for p in forward] == list(VARIABLES),
                "Forward after inverse")
        require([substitute(p, forward) for p in inverse] == list(VARIABLES),
                "Inverse after forward")
        degrees = [max(sum(mon) for mon in p) for p in inverse]
        require(degrees == [7, 5, 3, 1] and inverse[0][(0, 1, 0, 6)] == epsilon,
                "Exact inverse degree seven")
        # Parameter z is represented by the first polynomial variable.
        line = [scale(X, v) for v in g]
        expected_line = [add(scale(X, g[i]), scale(power(X, 3), -epsilon * b[i]),
                             scale(power(X, 5), 3 * c[i]), scale(power(X, 7), -3 * epsilon * e[i]))
                         for i in range(4)]
        require([substitute(p, line) for p in inverse] == expected_line,
                "Exact weighted inverse line cancellation")
        inverse_rows.append({"epsilon": epsilon, "coordinate_degrees": degrees,
                             "both_compositions_identity": True, "generator_line_identity": True})
    rejected = []
    def reject(ok, label):
        require(not ok, "Corrupted control accepted: " + label)
        rejected.append(label)
    corrupted = list(TERMS)
    corrupted[1] = (0, (1, 2, 2), Q(1))
    bad_j = [[derivative(p, i) for i in range(4)] for p in literal_map(corrupted)]
    bad_j2 = matrix_mul(bad_j, bad_j)
    bad_j4 = matrix_mul(bad_j2, bad_j2)
    reject(not matrix_nonzero(bad_j4), "Reverse the y*s^2 cubic sign and retain J^4=0")
    reject(vscale(mu(g, g, g), 2) == tuple(evaluate(p, g) for p in homogeneous),
           "Use polarization denominator three instead of six")
    reject(mu(g, b, c) == vscale(c, Q(1, 6)), "Reverse the even-family eigenvalue sign")
    reject(grade((0, bt, (0, 0, bt))) == grade((0, 0, bt)) + 1,
           "Count the fixed b word as a generator leaf")
    return {
        "literal_cubic_terms": [{"output": i, "inputs": list(slots), "coefficient": str(c)}
                               for i, slots, c in TERMS],
        "jacobian_cube": [{"row": i, "column": j, "polynomial": serialized_poly(p)}
                          for i, row in enumerate(j3) for j, p in enumerate(row) if p],
        "jacobian_fourth_power_nonzero_entries": 0,
        "exact_index_witness": {"a": serialized_vector(a0),
                               "L_a_cubed_e_t": serialized_vector(vscale(cube_at_a0_et, Q(1, 27)))},
        "polarization_basis_triples": len(tensor),
        "polarization_diagonal_and_derivative_polynomial_identities": True,
        "cyclic_words": {name: serialized_vector(value)
                         for name, value in (("g", g), ("b", b), ("c", c), ("e", e), ("d", d))},
        "cyclic_basis_determinant": str(det),
        "mixed_multiplier_columns": [serialized_vector(v) for v in mixed],
        "bounded_words": word_rows, "internal_grades_checked_inclusive": [0, 27],
        "inverse_controls": inverse_rows, "negative_controls_rejected": rejected,
        "corrupted_jacobian_fourth_power_nonzero_entries": sum(bool(p) for row in bad_j4 for p in row),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parents[2])
    parser.add_argument("--output", type=Path, default=HERE / RECEIPT)
    args = parser.parse_args()
    sources = check_sources(args.repo_root.resolve())
    results = run_checks()
    receipt = {"schema": 1, "status": "PASS", "python": platform.python_version(),
               "arithmetic": "Exact rational sparse polynomials; standard library only.",
               "source_ledger_sha256": sha(HERE / LEDGER), "sources": sources,
               "checks": results, "assurance": ASSURANCE,
               "scope": "Exact polynomial and linear identities plus bounded word controls. All grades and all characteristic-zero extensions follow from the written induction and universal evaluation argument; no inverse-growth or index-three conclusion."}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: exact J^3/J^4; 64 polarization triples; cyclic determinant 1/9; "
          "26 persistent word controls; both inverse compositions for each sign; "
          "4 corrupted controls rejected.")
    print(args.output)


if __name__ == "__main__":
    main()
