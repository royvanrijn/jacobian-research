#!/usr/bin/env python3
"""One genus-two Jacobian for each large inherited Sha block."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
from cubic_bridge import Cubic
from cover_experiment import mul, sub, trim

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "JACOBIAN_SHA_BLOCK_PROTOCOL.json"
OUTPUT = r.OUT / "rank_jump_jacobian_sha_blocks_v1.json"


def scale(poly, scalar):
    return trim([scalar * x for x in poly])


def add(a, b):
    return sub(a, scale(b, -1))


def symplectic_basis(B):
    """Return paired coordinates and the leftover radical, with exact checks."""
    n = len(B)
    assert all(B[i][i] == 0 for i in range(n))
    assert all(B[i][j] == B[j][i] for i in range(n) for j in range(n))
    active = [1 << i for i in range(n)]
    pairs = []
    while True:
        pair = next(((a, b) for i, a in enumerate(active)
                     for b in active[i+1:] if lc.pairing(a, b, B)), None)
        if pair is None:
            break
        a, b = pair
        pairs.append([a, b])
        active = lc.canonical([
            x ^ (a if lc.pairing(x, b, B) else 0)
              ^ (b if lc.pairing(x, a, B) else 0)
            for x in active
        ])
    flat = [x for pair in pairs for x in pair]
    assert r.rank(flat + active) == n
    assert len(flat) == r.rank(map(r.pack, B))
    assert active == lc.orthogonal(map(r.pack, B), n)
    for i, x in enumerate(flat):
        for j, y in enumerate(flat):
            assert lc.pairing(x, y, B) == int((i ^ 1) == j)
    return pairs, active


def geometry(A, B, u):
    assert u
    D = 1 + A*u*u + B*u**3
    disc = -4*A**3 - 27*B**2
    assert D and disc
    # Y=u^2 b, a=(1-z^2)/u.
    a = [r.F(1, u), 0, r.F(-1, u)]
    sextic = scale(add(add(mul(mul(a, a), a), scale(a, A)), [B]), u**4)
    expected = [u*D, 0, -3*u-A*u**3, 0, 3*u, 0, -u]
    assert sextic == expected
    numerator = add(scale(a, 1+A*u*u), [B*u*u])
    F = [B+A*B*u*u-B*B*u**3, A+3*B*u+A*A*u*u, 2*A*u, 1]
    image_equation = add(
        add(mul(mul(numerator, numerator), numerator),
            scale(mul(mul(numerator, numerator), [0, 0, 1]), F[2])),
        add(scale(mul(numerator, [0, 0, 0, 0, 1]), F[1]),
            [0, 0, 0, 0, 0, 0, F[0]]))
    assert image_equation == scale(sextic, D*D/r.F(u**4))
    K = Cubic(A, B)
    th = K.theta
    gamma = K.sub(K.one, K.scale(th, u))
    branch_numerator = K.add(K.scale(th, 1+A*u*u), K.scalar(B*u*u))
    branch_image = K.mul(branch_numerator, K.inverse(gamma))
    alpha = K.add(th, K.scale(K.square(th), u))
    assert branch_image == alpha
    return {
        "sextic_coefficients_low_to_high": list(map(str, sextic)),
        "sextic_discriminant": str(64*u**22*D*disc**2),
        "anchor_x_numerator": list(map(str, a)),
        "anchor_y_scale": str(r.F(1, u*u)),
        "specialized_x_numerator": list(map(str, numerator)),
        "specialized_x_denominator": "z^2",
        "specialized_y_multiplier": str(D/r.F(u*u)),
        "specialized_y_denominator": "z^3",
        "branch_image_in_cubic_basis": list(map(str, branch_image)),
        "both_quotient_equations_verified": True,
        "labelled_two_torsion_identification_verified": True,
    }


def nonisogeny_witness(A, B, u, primes):
    anchor = [B, A, 0, 1]
    target = [B+A*B*u*u-B*B*u**3, A+3*B*u+A*A*u*u, 2*A*u, 1]
    disc = -4*A**3 - 27*B**2
    D = 1+A*u*u+B*u**3
    checked = []
    def residue(x, p):
        return x.numerator * pow(x.denominator, -1, p) % p
    def count(coefficients, p):
        coefficients = [residue(r.F(c), p) for c in coefficients]
        total = 1
        for x in range(p):
            value = sum(c*pow(x, i, p) for i, c in enumerate(coefficients)) % p
            total += 1 if value == 0 else 2 if pow(value, (p-1)//2, p) == 1 else 0
        return total
    for p in primes:
        if residue(disc*D, p) == 0:
            continue
        old, new = count(anchor, p), count(target, p)
        checked.append({"prime": p, "anchor_point_count": old, "target_point_count": new})
        if old != new:
            return {"status": "PROVED_NOT_Q_ISOGENOUS", "checked_good_primes": checked,
                    "extension_splitting_status": "PROVED_NONSPLIT"}
    return {"status": "UNKNOWN", "checked_good_primes": checked,
            "extension_splitting_status": "UNKNOWN"}


def build(check=False):
    inp = r.read(lc.INPUT)
    A, B = map(r.F, inp["anchor"]["short_model_ainvariants"][3:])
    rows = []
    for u in r.read(PROTOCOL)["parameters"]:
        source = next(x for x in inp["rows"] if int(x["parameter_u"]) == u)
        assert source["all_local_kummer_images_complete"]
        w = source["W_u_basis"]
        matrix = next(x["matrix"] for x in inp["ct"] if x["u"] == u)
        pairs, radical = symplectic_basis(matrix)
        masks = [[lc.lift(x, w) for x in pair] for pair in pairs]
        rows.append({
            "u": u, "geometry": geometry(A, B, u),
            "isogeny_test": nonisogeny_witness(A, B, u, r.read(PROTOCOL)["isogeny_witness_primes"]),
            "inherited_selmer_dimension": len(w),
            "symplectic_pair_coordinate_masks": pairs,
            "independent_Sha_pair_anchor_masks": masks,
            "inherited_radical_anchor_masks": [lc.lift(x, w) for x in radical],
            "independent_Sha_block_dimension_lower_bound": 2*len(pairs),
            "projection_liftable_subspace_dimension_upper_bound": len(radical),
            "projection_index_in_locally_admissible_anchor_space_lower_bound": 2**(2*len(pairs)),
            "common_abelian_surface": "Jac(Y^2 = retained_sextic(z))",
            "Sha_block_killed_in_common_surface": True,
            "layer": "solubility obstruction",
        })
    paths = (PROTOCOL, Path(__file__), lc.INPUT, HERE/"cubic_bridge.py",
             HERE/"local_collision.py", HERE/"retrospective.py", HERE/"cover_experiment.py")
    out = {
        "schema": "rank-jump.jacobian-sha-block.v1",
        "bindings": {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths},
        "rows": rows,
        "scope": "A certified independent Sha[2] subspace killed in one explicit genus-two Jacobian per fibre. No full Sha dimension, numerical whole-curve rank upper bound, or new rational point is claimed.",
    }
    if check:
        assert r.read(OUTPUT) == out
        print("PASS common Jacobian and inherited Sha block")
    else:
        r.write_new(OUTPUT, out)
    for row in rows:
        print(row["u"], "Sha dimension >=", row["independent_Sha_block_dimension_lower_bound"],
              "projection index >=", row["projection_index_in_locally_admissible_anchor_space_lower_bound"])


def verify():
    from sage.all import QQ, PolynomialRing, EllipticCurve, GF
    inp = r.read(lc.INPUT)
    A, B = map(QQ, inp["anchor"]["short_model_ainvariants"][3:])
    R = PolynomialRing(QQ, "z")
    z = R.gen()
    for row in r.read(OUTPUT)["rows"]:
        u = row["u"]
        g = R(list(map(QQ, row["geometry"]["sextic_coefficients_low_to_high"])))
        assert g.discriminant() == QQ(row["geometry"]["sextic_discriminant"]) != 0
        a = (1-z*z)/u
        D = 1+A*u*u+B*u**3
        assert g == u**4*(a**3+A*a+B)
        x = (a+(A*a+B)*u*u)/(z*z)
        assert x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3 == g*D*D/(u**4*z**6)
        assert x(-z) == x
        for witness in row["isogeny_test"]["checked_good_primes"]:
            field = GF(witness["prime"])
            E0 = EllipticCurve(field, [A, B])
            Eu = EllipticCurve(field, [0, 2*A*u, 0, A+3*B*u+A*A*u*u,
                                      B+A*B*u*u-B*B*u**3])
            assert E0.cardinality() == witness["anchor_point_count"]
            assert Eu.cardinality() == witness["target_point_count"]
        print("PASS independent sextic discriminant and quotient maps", u)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "check", "verify"))
    mode = parser.parse_args().mode
    if mode == "verify":
        verify()
    else:
        build(mode == "check")
