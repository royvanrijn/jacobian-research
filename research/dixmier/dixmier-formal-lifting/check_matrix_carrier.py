#!/usr/bin/env python3
"""QMC1 finite matrix controls. Universal carrier results use written proofs."""
from pathlib import Path
import hashlib
import json
import math
import sympy as s

HERE = Path(__file__).resolve().parent
x, y, p, z, h = s.symbols("x y p z h")
positions, momenta = (x, y), (p, z)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mul(f, g):
    f, g = s.expand(f), s.expand(g)
    if f == 0 or g == 0:
        return s.Integer(0)
    degree = s.Poly(f, *momenta)
    out = 0
    for a in range(degree.degree(p)+1):
        for b in range(degree.degree(z)+1):
            left = s.diff(f, p, a, z, b)
            if left == 0:
                continue
            out += h**(a+b)*left*s.diff(g, x, a, y, b) / (
                math.factorial(a)*math.factorial(b))
    return s.expand(out)


def mm(a, b):
    return [[s.expand(sum(mul(a[i][k], b[k][j]) for k in range(2)))
             for j in range(2)] for i in range(2)]


def add(a, b, factor=1):
    return [[s.expand(a[i][j]+factor*b[i][j]) for j in range(2)]
            for i in range(2)]


def scale(c, a):
    return [[s.expand(c*v) for v in row] for row in a]


def zero(a, label):
    assert all(s.expand(v) == 0 for row in a for v in row), label


def comm(a, b):
    return add(mm(a, b), mm(b, a), -1)


def diag(a):
    return [[a, s.Integer(0)], [s.Integer(0), a]]


def corner(a):
    return [[a, s.Integer(0)], [s.Integer(0), s.Integer(0)]]


def bind_receipt(name):
    receipt = json.loads((HERE/name).read_text())
    assert receipt["status"] == "PASS", name
    for rel, expected in receipt["hashes"].items():
        path = HERE/rel
        if not path.exists():
            path = HERE.parent/"dixmier-audit"/rel
        assert path.is_file(), (name, rel)
        assert sha(path) == expected, (name, rel)
    return receipt


def main():
    for i in range(2):
        for j in range(2):
            assert s.expand(mul(momenta[i], positions[j]) -
                            mul(positions[j], momenta[i]) -
                            h*int(i == j)) == 0
    retained_projective = bind_receipt("projective_frame_controls.json")
    retained_degree = bind_receipt("division_index_controls.json")
    assert retained_projective["polynomial_splitting"]
    assert retained_degree["Long_classical_field_degree"]["value"] == 3

    one = diag(s.Integer(1))
    e0 = corner(s.Integer(1))
    # Product of two elementary matrices and its exact polynomial inverse.
    U = [[1, h*x], [h*p, 1+h**2*mul(p, x)]]
    V = [[1+h**2*mul(x, p), -h*x], [-h*p, 1]]
    zero(add(mm(U, V), one, -1), "U V inverse")
    zero(add(mm(V, U), one, -1), "V U inverse")
    e = mm(mm(V, e0), U)
    zero(add(mm(e, e), e, -1), "positive idempotent")
    images = [mm(mm(V, corner(a)), U) for a in (x, y, p, z)]
    assert any(v.has(x, p) for row in e for v in row)
    for i in range(2):
        for j in range(2):
            assert s.expand(e[i][j].subs(h, 0)-e0[i][j]) == 0
    for a, image in zip((x, y, p, z), images):
        zero(add(mm(mm(e, image), e), image, -1), "corner admission")
        zero([[s.expand(image[i][j].subs(h, 0)-corner(a)[i][j])
               for j in range(2)] for i in range(2)],
             "classical identity graph")
    relations = 0
    for i in range(4):
        for j in range(i+1, 4):
            kappa = -1 if (i, j) in ((0, 2), (1, 3)) else 0
            zero(add(comm(images[i], images[j]), scale(h*kappa, e), -1),
                 f"full positive relation {i},{j}")
            relations += 1

    # Actual QPF1 idempotent at m=2, before any h inversion.
    q = 1-h**2*x
    a0 = h*mul(q, x**2)
    b0 = 1+h**2*x+h*mul(p, x**2)
    row, col = (-p, q), (a0, b0)
    assert s.expand(mul(row[0], col[0])+mul(row[1], col[1])-1) == 0
    bad = [[s.expand(int(i == j)-mul(col[i], row[j]))
            for j in range(2)] for i in range(2)]
    zero(add(mm(bad, bad), bad, -1), "QPF1 idempotent")
    complement = add(one, bad, -1)
    A, B = diag(p), diag(x)
    compressed_A, compressed_B = mm(mm(bad, A), bad), mm(mm(bad, B), bad)
    defect = add(comm(compressed_A, compressed_B), scale(h, bad), -1)
    curvature = add(mm(mm(mm(bad, B), complement), mm(A, bad)),
                    mm(mm(mm(bad, A), complement), mm(B, bad)), -1)
    zero(add(defect, curvature, -1), "exact compression defect")
    assert any(v != 0 for rr in defect for v in rr), "nonzero compression defect"
    witnesses = []
    for i in range(2):
        for j in range(2):
            if defect[i][j] != 0:
                poly = s.Poly(defect[i][j], h, x, p)
                witnesses.append({"entry": [i, j],
                                  "terms": len(poly.terms()),
                                  "minimum_h_order": min(m[0] for m, c in poly.terms()),
                                  "factored": str(s.factor(defect[i][j]))})
    # Spectator variables commute with the projector exactly.
    for spectator in (y, z):
        zero(comm(bad, diag(spectator)), "spectator commutes with projector")

    names = ["FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
             "check_matrix_carrier.py", "PROJECTIVE_FRAME_BOUNDARY.md",
             "projective_frame_controls.json", "QUANTUM_DIVISION_AND_COMPLETION.md",
             "division_index_controls.json"]
    receipt = {
        "claim": "QMC1", "status": "PASS",
        "positive_polynomial_change_and_inverse": "PASS",
        "positive_nonconstant_idempotent": "PASS",
        "positive_full_rank_two_corner_relations": relations,
        "positive_classical_graph": "identity; not Long",
        "actual_QPF1_idempotent_specialization": "m=2",
        "compression_curvature_identity": "PASS",
        "compression_failure_witnesses": witnesses,
        "spectator_projector_commutators": 2,
        "retained_Long_classical_field_degree": 3,
        "hashes": {name: sha(HERE/name) for name in names},
        "external_dependencies": [
            {"source": "https://arxiv.org/html/2610.33919v1",
             "used": "Lemma3.1 generic projectivity; Lemma5.2 matrix Morita construction",
             "independently_reviewed": False},
            {"source": "https://doi.org/10.1112/jlms/s2-18.3.429",
             "used": "Stafford Theorem3.6: finite projective Weyl modules of rank>=2 are free",
             "independently_reviewed": False}],
        "coverage_boundary": (
            "Exact finite matrix identities and retained input hashes only. "
            "The universal carrier criterion, generic rank, division bounds, "
            "and non-surjectivity use written arguments and stated external theorems. "
            "No full right action on the nonfree QPF1 module and no finite Long "
            "carrier or scalar/matrix counterexample are constructed."),
        "independent_replay": False,
        "formal_verification": False, "external_review": False}
    (HERE/"matrix_carrier_controls.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QMC1", "status": "PASS",
                      "full_corner_relations": relations,
                      "nonzero_compression_entries": len(witnesses)}, indent=2))


if __name__ == "__main__":
    main()
