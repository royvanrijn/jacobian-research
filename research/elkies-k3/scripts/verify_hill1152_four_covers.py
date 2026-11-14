#!/usr/bin/env sage -python
"""Load and certify four generic rank>=18 covers from the frozen Hill1152 bank.

No specialization search. The maps are determined by exact rational conic
points; full section coordinates are reproducible from the retained loader.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from sage.all import EllipticCurve, PolynomialRing, QQ, matrix

ROOT = Path(__file__).resolve().parents[2]
PILOT = ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2"


def load_pilot():
    path = Path(__file__).with_name("pilot_hill1152_transfer.py")
    spec = importlib.util.spec_from_file_location("hill1152_generic_loader", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_family(index):
    """Return (curve over QQ(U), 18 ordered points, source Gram, map metadata)."""
    pilot = load_pilot()
    result = json.loads((PILOT/"result.json").read_text())
    row = result["quadratic_covers"][index]
    assert row.get("rational_conic") and row.get("smooth_branch")
    data = json.loads(pilot.SOURCE.read_text())
    R, K, A, B, E0, old_basis, G = pilot.build(data)
    v = R.gen()
    c = QQ(result["common_chart"]["c"])
    a, b = R(result["common_chart"]["A"]), R(result["common_chart"]["B"])
    q = R(row["q"])
    T, Z, S = map(QQ, row["conic_point"])
    assert Z != 0  # All four frozen witnesses satisfy this chart.
    t0, s0 = T/Z, S/Z
    N = R(t0*v*v-2*s0*v+q[2]*t0+q[1])
    D = R(v*v-q[2])
    W = R(-s0*v*v+(2*q[2]*t0+q[1])*v-q[2]*s0)
    assert q[2]*N*N+q[1]*N*D+q[0]*D*D == W*W
    assert N.gcd(D).degree() == 0
    u, s = K(N)/D, K(W)/D
    old_t = c+1/u
    curve = EllipticCurve(K, [R(D**8*a(u)), R(D**12*b(u))])
    # The first 17 points are inherited. The 18th is the residual chord.
    points = [curve(D**4*u**4*pilot.compose(P[0], old_t),
                    D**6*u**6*pilot.compose(P[1], old_t)) for P in old_basis]
    x0, x1, y0, y1 = [R(row[k]) for k in ("x0", "x1", "y0", "y1")]
    points.append(curve(D**4*(x0(u)+x1(u)*s), D**6*(y0(u)+y1(u)*s)))
    assert len(points) == 18
    return curve, points, G, {
        "index": index, "source_label": row["label"],
        "N": pilot.coeffs(N), "D": pilot.coeffs(D), "W": pilot.coeffs(W),
        "parent_old_parameter": str(old_t),
        "constant_squareclass_preserved": True,
    }


def height(point, chi=4):
    if not point:
        return QQ(0)
    x = point[0]
    return QQ(max(x.numerator().degree(), x.denominator().degree()+2*chi))


def verify_all():
    records = []
    for index in [0, 1, 2, 4]:
        curve, points, G, mapping = load_family(index)
        R = curve.a4().parent().ring()
        delta = R(curve.discriminant())
        assert delta.degree() == 48 and delta.is_squarefree()
        H = matrix(QQ, 18, 18)
        for i in range(18):
            H[i, i] = height(points[i])
            for j in range(i):
                H[i, j] = H[j, i] = (height(points[i]+points[j])-H[i, i]-H[j, j])/2
        assert H.is_positive_definite()
        assert H[:17, :17] == 2*G
        assert H.det() == 452984832
        schur = H.det()/H[:17, :17].det()
        assert schur == 3
        records.append({**mapping,
            "height_matrix": [[str(x) for x in row] for row in H.rows()],
            "determinant": str(H.det()), "extra_schur_complement": str(schur),
            "discriminant_degree": 48, "all_singular_fibres": "I1",
            "generic_rank_lower_bound": 18,
            "generic_exact_rank_or_saturation": "UNKNOWN"})
    return {"schema": "hill1152-four-generic-covers-v1", "status": "PASS",
        "pilot_input_sha256": sha256((PILOT/"input.json").read_bytes()).hexdigest(),
        "pilot_result_sha256": sha256((PILOT/"result.json").read_bytes()).hexdigest(),
        "families": records,
        "boundary": "Four certified rank>=18 function-field subgroups; no exact rank, saturation, curve-record or priority claim. Specialization must verify independence again."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    result = verify_all()
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "PASS", "families": 4, "rank_lower_bound": 18,
                      "determinant": "452984832", "schur_complement": 3}))
