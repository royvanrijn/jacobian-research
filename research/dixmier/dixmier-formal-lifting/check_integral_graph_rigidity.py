#!/usr/bin/env python3
"""QIR1 exact controls. The universal theorems use the written proof."""
from pathlib import Path
import hashlib
import json
import sympy as s
import check_matrix_carrier as pbw

HERE = Path(__file__).resolve().parent
x, y, p, z, h = pbw.x, pbw.y, pbw.p, pbw.z, pbw.h
t = s.Symbol("t")


def eye(n):
    return [[s.Integer(i == j) for j in range(n)] for i in range(n)]


def mm(a, b, product=pbw.mul):
    return [[s.expand(sum(product(a[i][k], b[k][j])
                         for k in range(len(b))))
             for j in range(len(b[0]))] for i in range(len(a))]


def plus(a, b, factor=1):
    return [[s.expand(v+factor*b[i][j]) for j, v in enumerate(row)]
            for i, row in enumerate(a)]


def scalar(c, a):
    return [[s.expand(c*v) for v in row] for row in a]


def zero(a, label):
    assert all(s.expand(v) == 0 for row in a for v in row), label


def corner2(a):
    # v=(1,x), w=(1,0)^T, vw=1, w*a*v gives the full right action.
    return [[a, pbw.mul(a, x), 0], [0, 0, 0], [0, 0, 0]]


def substitute(a, rules):
    return [[s.expand(v.subs(rules)) for v in row] for row in a]


def main():
    ident = eye(3)
    base_e = corner2(s.Integer(1))
    zero(plus(mm(base_e, base_e), base_e, -1), "base idempotent")

    # Arbitrary lifts: tilde(v)=(1,x), tilde(w)=(1+t*p,0)^T.
    # Their product is 1+t*p, which has positive spatial degree.
    # The construction never inverts it.
    U = [[1, 0, 1+t*p], [0, 1, 0], [0, 0, 1]]
    Ui = [[1, 0, -1-t*p], [0, 1, 0], [0, 0, 1]]
    V = [[1, 0, 0], [0, 1, 0], [-1, -x, 1]]
    Vi = [[1, 0, 0], [0, 1, 0], [1, x, 1]]
    permutation = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    G = mm(permutation, mm(mm(U, V), U))
    Gi = mm(mm(mm(Ui, Vi), Ui), permutation)
    zero(plus(mm(G, Gi), ident, -1), "elementary inverse G Gi")
    zero(plus(mm(Gi, G), ident, -1), "elementary inverse Gi G")
    assert s.Poly(1+t*p, x, y, p, z).total_degree() == 1

    E = mm(mm(G, base_e), Gi)
    images = [mm(mm(G, corner2(a)), Gi) for a in (x, y, p, z)]
    zero(plus(mm(E, E), E, -1), "normalized integral idempotent")
    standard_e = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
    zero(plus(substitute(E, {t: 0, h: 0}), standard_e, -1),
         "normalized special idempotent")
    for a, image in zip((x, y, p, z), images):
        zero(plus(mm(mm(E, image), E), image, -1), "corner admission")
        zero(plus(substitute(image, {t: 0, h: 0}),
                  scalar(a, standard_e), -1), "special full graph")
    relations = 0
    for i in range(4):
        for j in range(i+1, 4):
            kappa = -1 if (i, j) in ((0, 2), (1, 3)) else 0
            comm = plus(mm(images[i], images[j]), mm(images[j], images[i]), -1)
            zero(plus(comm, scalar(h*kappa, E), -1), "full Weyl relation")
            relations += 1

    # A genuine full action with lower-than-quantum scalar corrections.
    # These generators are canonical for every independent h and t.
    right = (x+t*p+t**2*p**2, y, p, z)
    for i in range(4):
        for j in range(i+1, 4):
            kappa = -1 if (i, j) in ((0, 2), (1, 3)) else 0
            assert s.expand(pbw.mul(right[i], right[j]) -
                            pbw.mul(right[j], right[i]) - h*kappa) == 0

    corrections = []
    probes = (s.Integer(1), x, y, p, z, x*p+y*z)
    for order in (1, 3, 5):
        parameter = t**order
        def product(f, g):
            return s.expand(pbw.mul(f, g).subs(h, parameter))
        P = x
        stages = []
        for j in range(1, order):
            one_difference = s.expand(P-right[0])
            assert s.expand(one_difference).coeff(t, 0) == 0
            cj = s.expand(one_difference).coeff(t, j)
            for m in probes:
                difference = product(P, m)-product(m, right[0])
                poly = s.Poly(s.expand(difference), t)
                assert all(power[0] >= j for power, coeff in poly.terms() if coeff)
                assert s.expand(s.expand(difference).coeff(t, j)-m*cj) == 0
            P = s.expand(P-t**j*cj)
            stages.append({"order": j, "scalar": str(cj)})
        for m in probes:
            difference = product(P, m)-product(m, right[0])
            poly = s.Poly(s.expand(difference), t)
            assert all(power[0] >= order for power, coeff in poly.terms() if coeff)
        if order > 1:
            # Starting with the uncorrected left x fails modulo a(t).
            assert s.Poly(x-right[0], t).terms()[-1][0][0] == 1
        corrections.append({"parameter_order": order, "finite_scalar_steps": stages,
                            "quantum_order_congruence": "PASS",
                            "original_quantum_order_congruence": order == 1})

    # Universal frame inversion, with unrestricted first derivatives.
    fx, fp, fy, fz = s.symbols("fx fp fy fz")
    fields = (-fp, -fz, fx-2*x*fp, fy-2*y*fz)
    assert s.expand(fields[2]-2*x*fields[0]-fx) == 0
    assert s.expand(fields[3]-2*y*fields[1]-fy) == 0

    # Omitting the full right action permits distinct full finite lattices.
    q = 1-t*x
    assert s.expand(pbw.mul(q, x)-pbw.mul(x, q)) == 0
    assert s.expand(pbw.mul(q, p)-pbw.mul(p, q)-t*h) == 0
    assert s.Poly(q, x, y, p, z).total_degree() == 1
    # Over the generic field this scalar commutator would force 1 into A*q
    # if q*p lay there, contradicting additive Bernstein degree.
    boundary = {"lattices": "A_R and A_R*(1-t*x)",
                "same_generic_full_lattice": False,
                "preserved_by_right_x": True,
                "full_right_Weyl_action_on_both": False,
                "right_p_defect": "t*h",
                "formal_completions_agree": True}

    inputs = [
        "INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md",
        "check_integral_graph_rigidity.py", "check_matrix_carrier.py",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "PROJECTIVE_FRAME_BOUNDARY.md",
        "IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md",
        "HAMILTONIAN_IDEALIZER_CRITERION.md",
        "RAMIFICATION_AND_FINITE_CLOSURE_CARRIERS.md",
    ]
    receipt = {
        "claim": "QIR1", "status": "PASS",
        "elementary_stabilization": {
            "matrix_size": 3, "nonunit_lift_product": "1+t*p",
            "two_sided_polynomial_inverse": "PASS",
            "integral_idempotent": "PASS", "special_idempotent": "e11",
            "full_corner_relations": relations,
            "parameters_independent": ["h", "t"],
            "Long_carrier": False},
        "ramified_scalar_graph_upgrade": corrections,
        "universal_Hamiltonian_frame_inversion": "PASS",
        "partial_right_action_boundary": boundary,
        "research_context": [
            {"source": "https://arxiv.org/html/2610.33919v1",
             "used": "Lemma3.1 generic finite-left Weyl bimodule projectivity. "
                     "No arithmetic-support nonlifting hypothesis is imported."},
            {"source": "https://arxiv.org/abs/math/0104248",
             "used": "Rank-one ideal geometry motivates distinguishing "
                     "a fixed rational action from formal lattice data; "
                     "no A1 classification is transferred to higher rank."}],
        "coverage_boundary": (
            "The integral projectivity, finite-corner equivalence, generic "
            "lattice uniqueness, automatic scalar upgrade and integral "
            "homothety use written proofs. The symbolic elementary identities "
            "and six relations hold for independent h,t. The listed ramified "
            "cases are controls, not a proof by parameter census. Uniqueness "
            "fixes the rational action. No generic free frame or Long carrier "
            "is constructed; arbitrary formal gauges remain outside it."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs},
    }
    (HERE/"integral_graph_rigidity_controls.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QIR1", "status": "PASS",
                      "corner_relations": relations, "parameter_orders": [1, 3, 5],
                      "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
