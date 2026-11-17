#!/usr/bin/env python3
"""QTC1 exact controls; universal composition and boundary results are written proofs."""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import sympy as s
import check_matrix_carrier as pbw

HERE = Path(__file__).resolve().parent
x, y, p, z, h = pbw.x, pbw.y, pbw.p, pbw.z, pbw.h
t = s.Symbol("t")
variables = (x, y, p, z)


@lru_cache(maxsize=None)
def product(f, g):
    return pbw.mul(f, g)


def mm(a, b):
    return [[s.expand(sum(product(s.sympify(a[i][k]), s.sympify(b[k][j]))
                         for k in range(len(b))))
             for j in range(len(b[0]))] for i in range(len(a))]


def plus(a, b, factor=1):
    return [[s.expand(v+factor*b[i][j]) for j, v in enumerate(row)]
            for i, row in enumerate(a)]


def scale(c, a):
    return [[s.expand(c*v) for v in row] for row in a]


def zero(a, label):
    assert all(s.expand(v) == 0 for row in a for v in row), label


def corner_action(a, q):
    return [[s.sympify(a), product(s.sympify(a), q)], [0, 0]]


def test_relations(E, images):
    zero(plus(mm(E, E), E, -1), "idempotent")
    for image in images:
        zero(plus(mm(mm(E, image), E), image, -1), "corner membership")
    for i in range(4):
        for j in range(i+1, 4):
            kappa = -1 if (i, j) in ((0, 2), (1, 3)) else 0
            comm = plus(mm(images[i], images[j]), mm(images[j], images[i]), -1)
            zero(plus(comm, scale(h*kappa, E), -1), "full Weyl relation")


def degrees(matrices):
    return max(s.Poly(s.expand(v), *variables).total_degree()
               for a in matrices for row in a for v in row if v != 0)


def main():
    retained = json.loads((HERE/"identity_lift_controls.json").read_text())
    assert retained["claim"] == "QNI1" and retained["status"] == "PASS"
    assert retained["uniform_forward_bound"] == 12
    for name, expected in retained["hashes"].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, name
    retained_path = HERE.parent.parent/retained["retained_input"]["path"]
    assert hashlib.sha256(retained_path.read_bytes()).hexdigest() == retained["retained_input"]["sha256"]
    psi = (x, y, p+x**2, z)
    chi = (x+p**2, y, p, z)
    E = corner_action(1, t*x)
    F = corner_action(1, t*y)
    rho = [corner_action(a, t*x) for a in psi]
    sigma = [corner_action(a, t*y) for a in chi]
    test_relations(E, rho)
    test_relations(F, sigma)
    assert E[0][0] == F[0][0] == 1  # positive trace-ideal generators.

    @lru_cache(maxsize=None)
    def evaluate(f):
        out = [[s.Integer(0), s.Integer(0)], [s.Integer(0), s.Integer(0)]]
        for powers, coeff in s.Poly(s.expand(f), *variables).terms():
            term = E
            for i, exponent in enumerate(powers):
                for _ in range(exponent):
                    term = mm(term, rho[i])
            out = plus(out, scale(coeff, term))
        return out

    def block(a):
        blocks = [[evaluate(s.sympify(v)) for v in row] for row in a]
        return [[blocks[i//2][j//2][i % 2][j % 2]
                 for j in range(2*len(a[0]))] for i in range(2*len(a))]

    H = block(F)
    composed = [block(image) for image in sigma]
    test_relations(H, composed)
    expected = (x+(p+x**2)**2, y, p+x**2, z)
    classical = lambda a: [[s.expand(v.subs({t: 0, h: 0})) for v in row] for row in a]
    e11 = [[int(i == j == 0) for j in range(4)] for i in range(4)]
    zero(plus(classical(H), e11, -1), "normalized composite special idempotent")
    for image, target in zip(composed, expected):
        zero(plus(classical(image), scale(target, e11), -1),
             "classical composition, not reversed")
    dM, dN, dP = degrees([E, *rho]), degrees([F, *sigma]), degrees([H, *composed])
    assert dP <= dM*dN

    # Long's actual source chart and admitted second point.
    r, v, a, b, d, w, X = s.symbols("r v a b d w X")
    A = 3*r*a-2
    Q = A**2-8*r*v
    ss = 2*v+4*a-3*r*a**2
    tt = (r*a**3-a**2-v*a)/2
    f = r*X**3-2*X**2+ss*X+4*tt
    g = r*X**2+(r*a-2)*X+(r*a**2-2*a+ss)
    assert s.expand(f-(X-a)*g) == 0
    assert s.expand(g.subs(X, a)-2*v) == 0
    diagonal = g/(2*v)
    identity = s.cancel((diagonal**2-diagonal) -
                        f*(r*X+2*r*a-2)/(4*v**2))
    assert identity == 0

    b_expr = (w+2-r*a)/(2*r)
    reduce_w = lambda e: s.cancel(s.rem(s.together(e).as_numer_denom()[0],
                                     w**2-Q, w) /
                                  s.together(e).as_numer_denom()[1])
    assert reduce_w(g.subs(X, b_expr)) == 0
    assert s.cancel(b_expr-a-(w-A)/(2*r)) == 0
    v_b = w*(w-A)/(4*r)
    assert reduce_w(s.diff(f, X).subs(X, b_expr)/2-v_b) == 0
    assert s.expand((w-A)*(w+A)+8*r*v-(w**2-Q)) == 0
    # w-A is already a unit: its inverse is -(w+A)/(8*r*v).
    assert reduce_w((w-A)*(-(w+A)/(8*r*v))-1) == 0
    # Inverting v_b forces w^-1; conversely w^-1 and (w-A)^-1 give v_b^-1.
    assert s.cancel((w-A)/(4*r*v_b)-1/w) == 0
    assert reduce_w(1/(w**2)-1/Q) == 0
    assert s.diff(Q, v) == -8*r
    assert Q.subs({r: 1, a: 1, v: s.Rational(1, 8)}) == 0
    assert 2*s.Rational(1, 8) != 0

    # These are index arithmetic examples, not newly constructed actions.
    index_cases = []
    for deltaM, NM, deltaN, NN in [(1, 3, 1, 3), (3, 3, 1, 3),
                                   (3, 4, 3, 5), (3, 3, 3, 3)]:
        gapM, gapN = NM-deltaM, NN-deltaN
        gap = NM*NN-deltaM*deltaN
        assert gap == gapM*deltaN+deltaM*gapN+gapM*gapN
        assert (gap == 0) == (gapM == gapN == 0)
        index_cases.append({"factor_completed_degrees": [deltaM, deltaN],
                            "factor_generic_indices": [NM, NN],
                            "composite_excess": gap})

    inputs = [
        "COMPOSITION_DUALITY_AND_BOUNDARY_OF_QUANTUM_GRAPHS.md",
        "check_quantum_graph_composition.py", "check_matrix_carrier.py",
        "INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "QUANTUM_DIVISION_AND_COMPLETION.md",
        "IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md", "identity_lift_controls.json",
        "RAMIFICATION_AND_FINITE_CLOSURE_CARRIERS.md",
        "ARITHMETIC_COTANGENT_EXCLUSION.md", "BOUNDED_REALIZATION.md",
    ]
    out = {
        "claim": "QTC1", "status": "PASS",
        "forward_composition": {
            "factor_matrix_sizes": [2, 2], "composite_matrix_size": 4,
            "factor_degree_bounds": [int(dM), int(dN)],
            "composite_degree": int(dP), "product_degree_bound": int(dM*dN),
            "all_six_relations": "PASS", "classical_composition_order": "chi(psi)",
            "parameters_independent": ["h", "t"],
            "special_idempotent": "e11", "actual_Long_carrier": False},
        "Long_off_diagonal": {
            "cubic_factorization": "PASS", "diagonal_idempotent": "PASS",
            "quadratic_relation": "w^2=Q",
            "second_source_v": "w*(w-A)/(4*r)",
            "other_factor_already_unit": "(w-A)*(w+A)=-8*r*v",
            "admitted_ring": "S0[w,w^-1]/(w^2-Q)",
            "contains_Q_inverse": True,
            "interior_branch_point": {"r": 1, "a": 1, "v": "1/8"},
            "finite_over_original_source": False,
            "finite_etale_after_inverting_Q": True,
            "universal_nonfiniteness_uses_written_integrality_argument": True},
        "generic_index_arithmetic": index_cases,
        "actual_QNI_iterates": {"source": "QNI1", "ranks": "n>=3",
            "generic_left_and_right_indices": "3^k",
            "completed_degree": 1, "forward_degree_bound": "12^k",
            "arbitrary_identity_contact": True, "k_range": "all k>=1",
            "universal_proof": "forward composition and index tower law",
            "new_rank_two_example": False},
        "research_context": [
            {"source": "https://stacks.math.columbia.edu/tag/02LA",
             "used": "Finiteness descends under the fpqc map X->f(X)."},
            {"source": "https://stacks.math.columbia.edu/tag/0BMB",
             "used": "Finite etale covers extend across a codimension>=2 "
                     "complement in regular affine space by purity."},
            {"source": "https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27",
             "used": "Characteristic-zero affine space has no nontrivial "
                     "geometrically connected finite etale covers."},
            {"source": "https://arxiv.org/html/2604.16636v1",
             "used": "Theorem4.6 is a first-order Azumaya lifting criterion; "
                     "it gives no uniform finite realization theorem."},
            {"source": "https://arxiv.org/html/2601.23110v1",
             "used": "The Weyl degree<prime result concerns injectivity and "
                     "first-order Witt lifting, not all-order bounded support."},
        ],
        "coverage_boundary": (
            "Integral progenerator, block composition, multiplicative "
            "division indices, self-correspondence nonfiniteness and "
            "evaluation-kernel reduction are written proofs. Symbolic "
            "positive controls are actual forward full Weyl actions. "
            "The Long root algebra is a classical correspondence; its "
            "nonfiniteness does not exclude a hypothetical original finite "
            "carrier. No generic left nonfiniteness is inferred solely "
            "from a nonfinite special fiber. First-order arithmetic "
            "lifting is not prescribed constant arithmetic support."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs},
    }
    (HERE/"quantum_graph_composition_controls.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps({"claim": "QTC1", "status": "PASS",
                      "corner_relations": 6,
                      "degrees": [int(dM), int(dN), int(dP)],
                      "Long_localization": "PASS", "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
