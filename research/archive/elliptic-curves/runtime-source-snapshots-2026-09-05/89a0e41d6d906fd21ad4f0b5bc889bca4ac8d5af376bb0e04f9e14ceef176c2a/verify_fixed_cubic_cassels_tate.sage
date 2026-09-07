#!/usr/bin/env sage
"""Replay the u=-1 Cassels--Tate certificate without descent or point search.

Fisher, On binary quartics and the Cassels--Tate pairing (2022), Theorem 3.1
and Remark 3.3. All covering maps, square roots, local square witnesses,
Hilbert symbols and the complete finite support are checked exactly.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector, prod, pari

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from fixed_cubic_field_curve_family import field_multiply, field_product

RESULTS = ROOT / "artifacts/generated-results/elliptic-curves"
SOURCE = RESULTS / "fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
EVIDENCE = RESULTS / "fixed_cubic_u_minus1_cassels_tate_evidence_v1.json.gz"
SUMMARY = RESULTS / "fixed_cubic_u_minus1_cassels_tate_v1.json"
R = PolynomialRing(QQ, "x")
x = R.gen()


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def invariants(q):
    e, d, c, b, a = [q[i] for i in range(5)]
    return (12*a*e - 3*b*d + c*c,
            72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3)


def rational_unit(value, prime):
    value, prime = QQ(value), ZZ(prime)
    require(value != 0, "zero has no unit squareclass")
    valuation = value.valuation(prime)
    return valuation, value / prime**valuation


def residue(value, modulus):
    value = QQ(value)
    return int(value.numerator() * value.denominator().inverse_mod(modulus) % modulus)


def legendre_unit(value, prime):
    return 1 if pow(residue(value, prime), (int(prime)-1)//2, int(prime)) == 1 else -1


def local_square(value, prime):
    if value == 0:
        return True
    valuation, unit = rational_unit(value, prime)
    return valuation % 2 == 0 and (residue(unit, 8) == 1 if prime == 2
                                 else legendre_unit(unit, prime) == 1)


def hilbert_symbol(left, right, prime):
    """Elementary rational Hilbert formula; independent of PARI hilbert."""
    require(left != 0 and right != 0, "Hilbert arguments must be nonzero")
    if prime == "infinity":
        return -1 if left < 0 and right < 0 else 1
    va, ua = rational_unit(left, prime)
    vb, ub = rational_unit(right, prime)
    if prime == 2:
        a, b = residue(ua, 8), residue(ub, 8)
        exponent = (a-1)*(b-1)//4 + va*(b*b-1)//8 + vb*(a*a-1)//8
        return -1 if exponent % 2 else 1
    answer = -1 if (va*vb*((int(prime)-1)//2)) % 2 else 1
    if vb % 2:
        answer *= legendre_unit(ua, prime)
    if va % 2:
        answer *= legendre_unit(ub, prime)
    return answer


def primitive(q):
    denominator = q.denominator()
    integral = q * denominator
    content = ZZ(0)
    for coefficient in integral:
        content = content.gcd(ZZ(coefficient))
    require(content != 0, "zero pairing quadratic")
    answer = integral / content
    return -answer if answer.leading_coefficient() < 0 else answer


def fisher_gamma(quartics, square_root, I, J):
    """Check m^2=z1*z2*z3 and derive the primitive gamma of Theorem 3.1."""
    T = PolynomialRing(QQ, "phi")
    phi = T.gen()
    L = T.quotient(phi**3 - 3*I*phi + J, "phi_bar")
    phi = L.gen()
    S = PolynomialRing(L, "t")
    cubic = []
    for q in quartics:
        require(q.degree() == 4 and invariants(q) == (I, J), "quartic invariant mismatch")
        a, b, c = q[4], q[3], q[2]
        cubic.append((4*a*phi + 3*b*b - 8*a*c)/3)
    m = sum(QQ(c)*phi**i for i, c in enumerate(square_root))
    require(m != 0 and m*m == prod(cubic), "invalid cubic square-root identity")
    q = quartics[0]
    e, d, c, b, a = [q[i] for i in range(5)]
    hessian = R([3*d*d-8*c*e, 4*(c*d-6*b*e),
                 2*(2*c*c-24*a*e-3*b*d), 4*(b*c-6*a*d), 3*b*b-8*a*c])
    H = (4*phi*S(q.derivative(2)) + S(hessian.derivative(2)))/36
    H += QQ(2)/9*(I-phi*phi)
    require(H[2] == cubic[0], "H normalization mismatch")
    product = (m/cubic[0])*H
    return primitive(R([product[k].lift()[2] for k in range(3)]))


def source_data():
    source = json.loads(SOURCE.read_text())
    require(source["status"] == "PASS_EXACT_FULL_SPAN_LOCAL_INTERSECTIONS_NO_CLASS_GROUP",
            "uncertified local intersection input")
    run = next(row for row in source["runs"] if row["parameter_u"] == "-1")
    require(run["all_local_kummer_images_complete"] is True and run["W_u_dimension"] == 18,
            "wrong surviving subspace")
    B, A, _, _ = map(QQ, source["anchor"]["base_polynomial_ascending"])
    _, a2, _, a4, a6 = map(QQ, run["raw_curve_ainvariants"])
    return source, run, A, B, a2*a2-3*a4, -2*a2**3+9*a2*a4-27*a6


def verify_cover(record, source, A, B, I, J):
    mask = record["anchor_mask"]
    require(type(mask) is int and 0 < mask < 2**20, "invalid anchor mask")
    beta = field_product([list(map(QQ, row)) for i, row in enumerate(
        source["anchor"]["known_kummer_basis_beta_power_coordinates"]) if mask >> i & 1], A, B)
    require(tuple(map(QQ, record["beta"])) == beta, "cover class mismatch")
    q = R(list(map(QQ, record["quartic"])))
    require(q.degree() == 4 and q.denominator() == 1, "nonintegral or degenerate quartic")
    require(invariants(q) == (I, J), "cover has the wrong Jacobian invariants")
    parametrization = matrix(QQ, record["parametrization"])
    change = matrix(QQ, record["parameter_transform"])
    require(parametrization.dimensions() == (3, 3) and parametrization.det() != 0,
            "degenerate conic parametrization")
    require(change.dimensions() == (2, 2) and change.det() != 0, "degenerate parameter change")
    s, t = change * vector([x, 1])
    gamma = parametrization * vector([s*s, s*t, t*t])
    scale = QQ(record["d_over_quartic_y"])
    require(scale != 0, "zero covering scale")
    coefficients = field_multiply(beta, field_multiply(gamma, gamma, A, B), A, B)
    require(coefficients[1] + scale**2*q == 0 and coefficients[2] - scale**2*q == 0,
            "quartic does not map to the declared two-quadric cover")
    return q


def verify_pair(record, covers, support, I, J, proved_primes):
    masks = record["basis_anchor_masks"]
    require(len(masks) == 2 and masks[0] != masks[1], "invalid pair masks")
    sum_mask = masks[0] ^ masks[1]
    require(record["sum_anchor_mask"] == sum_mask, "wrong sum cover")
    quartics = [covers[mask] for mask in masks + [sum_mask]]
    gamma = fisher_gamma(quartics, record["square_root_phi_coefficients"], I, J)
    require(gamma == R(list(map(QQ, record["gamma"]))), "wrong pairing quadratic")
    leading = quartics[1][4]
    factors = record["leading_coefficient_factorization"]
    primes = [int(prime) for prime, exponent in factors]
    require(len(primes) == len(set(primes)), "duplicate prime in factorization")
    require(all(type(exponent) is int and exponent > 0 for prime, exponent in factors),
            "invalid factor exponent")
    require(prod(ZZ(prime)**exponent for prime, exponent in factors) == abs(leading),
            "incomplete leading coefficient factorization")
    expected_places = sorted(set(support + [2, 3, 5, 7] + primes))
    for prime in expected_places:
        if prime not in proved_primes:
            require(ZZ(prime).is_prime(proof=True), f"unproved prime {prime}")
            proved_primes.add(prime)
    remaining = abs(ZZ(quartics[0].discriminant()))
    require(remaining != 0, "singular quartic")
    for prime in expected_places:
        remaining //= ZZ(prime)**remaining.valuation(prime)
    require(remaining == 1, "missing discriminant prime")
    # gamma is primitive integral. Remark 3.3 therefore proves that every
    # place outside expected_places + infinity contributes +1.
    terms = record["local_terms"]
    require([term["place"] for term in terms] == expected_places + ["infinity"],
            "incomplete local pairing support")
    answer = 1
    for term in terms:
        prime, xx = term["place"], QQ(term["x"])
        qvalue, gvalue = quartics[0](xx), gamma(xx)
        require(qvalue != 0 and gvalue != 0, "zero local evaluation")
        require(qvalue == QQ(term["q_value"]) and gvalue == QQ(term["gamma_value"]),
                "incorrect local evaluation")
        require(qvalue > 0 if prime == "infinity" else local_square(qvalue, prime),
                "local witness is not on the cover")
        symbol = hilbert_symbol(leading, gvalue, prime)
        require(type(term["hilbert_symbol"]) is int and symbol == term["hilbert_symbol"],
                "incorrect local Hilbert symbol")
        answer *= symbol
    value = int(answer == -1)
    require(type(record["value"]) is int and record["value"] == value, "incorrect pairing entry")
    return value


def coordinates_record(v, basis):
    anchor_mask = 0
    for bit, row in zip(v, basis):
        if bit:
            anchor_mask ^= row["mask"]
    return {"W_coordinates": list(map(int, v)),
            "one_based_W_indices": [i+1 for i, bit in enumerate(v) if bit],
            "anchor_mask": anchor_mask,
            "one_based_anchor_indices": [i+1 for i in range(20) if anchor_mask >> i & 1]}


def decompose(M, basis):
    """Construct a symplectic complement and the full restricted radical."""
    vectors = list(M.row_ambient_module().basis())
    pairs = []
    while vectors:
        pair = next(((i, j) for i in range(len(vectors)) for j in range(i+1, len(vectors))
                     if vectors[i]*M*vectors[j]), None)
        if pair is None:
            break
        i, j = pair
        a, b = vectors[i], vectors[j]
        pairs.append([a, b])
        vectors = [v + (v*M*b)*a + (v*M*a)*b for k, v in enumerate(vectors) if k not in pair]
    change = matrix(GF(2), [v for pair in pairs for v in pair] + vectors)
    require(change.rank() == M.nrows(), "symplectic change is not invertible")
    normal = change*M*change.transpose()
    target = matrix(GF(2), M.nrows())
    for i in range(len(pairs)):
        target[2*i, 2*i+1] = target[2*i+1, 2*i] = 1
    require(normal == target, "incorrect symplectic decomposition")
    radical = list(M.right_kernel().basis())
    return {"pairing_rank": int(M.rank()), "restricted_radical_dimension": len(radical),
            "restricted_radical_basis": [coordinates_record(v, basis) for v in radical],
            "symplectic_pairs": [[coordinates_record(v, basis) for v in pair] for pair in pairs],
            "obstructed_class_count": 2**M.nrows()-2**len(radical),
            "nonzero_compatible_class_count": 2**len(radical)-1,
            "rational_kummer_intersection_dimension_upper": len(radical),
            "sha_two_torsion_dimension_lower": int(M.rank()),
            "full_curve_rank_upper": None,
            "full_selmer_radical_computed": False}


def verify_published_control(record):
    quartics = [R(list(map(QQ, q))) for q in record["quartics"]]
    require(quartics == [R([-64, -164, -52, 68, -11]),
                         R([-3, -52, -232, -60, -4]),
                         R([-53, 102, 32, -78, -31])], "wrong published control")
    gamma = fisher_gamma(quartics, record["square_root_phi_coefficients"], 44608, 18842960)
    require(gamma == 5*x*x - 16*x - 12, "Fisher Example 3.4 gamma mismatch")
    require([t["place"] for t in record["local_terms"]] == [2, 3, 5, 7, 571, "infinity"],
            "incomplete published control")
    answer = 1
    for term in record["local_terms"]:
        prime, xx = term["place"], QQ(term["x"])
        qvalue, gvalue = quartics[0](xx), gamma(xx)
        require(qvalue != 0 and (qvalue > 0 if prime == "infinity" else local_square(qvalue, prime)),
                "invalid control point")
        require(qvalue == QQ(term["q_value"]) and gvalue == QQ(term["gamma_value"]),
                "invalid control evaluation")
        symbol = hilbert_symbol(-4, gvalue, prime)
        require(symbol == term["hilbert_symbol"], "wrong control Hilbert symbol")
        answer *= symbol
    require(answer == -1 and record["value"] == 1, "published nonzero control failed")


def verify_radical_search(rows, M, basis, covers):
    require(M.right_kernel().dimension() == 2 and len(rows) == 3, "wrong radical search scope")
    seen, result = set(), []
    for row in rows:
        v = vector(GF(2), row["W_coordinates"])
        require(len(v) == 18 and v != 0 and v*M == 0, "search class is outside the radical")
        coordinates = coordinates_record(v, basis)
        mask = coordinates["anchor_mask"]
        require(mask not in seen, "duplicate radical search class")
        seen.add(mask)
        search = row["search"]
        require(search["bound"] == 10000 and search["status"] == "BOUNDED_SEARCH_COMPLETE",
                "unrecognized bounded point search")
        require(search["rank_claim"] is None, "bounded search cannot assert rank")
        q = R(list(map(QQ, search["quartic"])))
        require(q*QQ(search["quartic_y_rescaling"])**2 == covers[mask], "wrong search cover")
        # The fixed height-10000 replay is cheap on these three quartics.
        points = [[str(z) for z in pt] for pt in pari.hyperellratpoints(q, 10000)]
        require(points == search["raw_points"], "bounded point-search replay mismatch")
        infinity = [] if not q[4].is_square() else [str(q[4].sqrt())]
        require(infinity == search["rational_points_at_infinity"], "wrong infinity-point check")
        result.append({**coordinates, "quartic_ascending": search["quartic"],
                       "coefficient_height": str(max(abs(z) for z in q)),
                       "point_height_bound": 10000, "affine_point_count": len(points),
                       "infinity_point_count": len(infinity),
                       "point_or_sha_status": "UNKNOWN"})
    return sorted(result, key=lambda row: (ZZ(row["coefficient_height"]), row["anchor_mask"]))


def verify(evidence):
    source, run, A, B, I, J = source_data()
    require(evidence["source_sha256"] == digest(SOURCE), "source hash mismatch")
    require(evidence["status"] == "COMPLETE", "incomplete arithmetic evidence")
    covers = {}
    for record in evidence["covers"]:
        mask = record["anchor_mask"]
        require(mask not in covers, "duplicate covering")
        covers[mask] = verify_cover(record, source, A, B, I, J)
    M = matrix(GF(2), 18)
    seen, proved_primes = set(), set()
    for row in evidence["pairings"]:
        i, j = row["i"], row["j"]
        require(type(i) is int and type(j) is int and 0 <= i < j < 18 and (i, j) not in seen,
                "duplicate or invalid pairing indices")
        require(row["basis_anchor_masks"] == [run["W_u_basis"][k]["mask"] for k in (i, j)],
                "pairing is not bound to the declared basis")
        seen.add((i, j))
        M[i, j] = M[j, i] = verify_pair(row, covers, run["complete_finite_place_support"], I, J, proved_primes)
    require(len(seen) == 153, "missing pairing entries; no radical claim permitted")
    for row in evidence.get("cross_checks", []):
        value = verify_pair(row, covers, run["complete_finite_place_support"], I, J, proved_primes)
        left, right = [vector(GF(2), v) for v in row["W_vectors"]]
        require(row["basis_anchor_masks"] == [coordinates_record(v, run["W_u_basis"])["anchor_mask"]
                                             for v in (left, right)], "wrong cross-check classes")
        require(value == int(left*M*right), "symmetry or bilinearity check failed")
    if "published_control" in evidence:
        verify_published_control(evidence["published_control"])
    searches = (verify_radical_search(evidence["radical_search"], M, run["W_u_basis"], covers)
                if "radical_search" in evidence else [])
    return {"matrix": [list(map(int, row)) for row in M.rows()],
            "verified_cover_count": len(covers), "verified_pairing_entry_count": len(seen),
            "verified_cross_check_count": len(evidence.get("cross_checks", [])),
            "published_nonzero_control_verified": "published_control" in evidence,
            "radical_search_candidates_by_coefficient_height": searches,
            "preferred_radical_generator_anchor_masks": [row["anchor_mask"] for row in searches[:2]],
            "verified_distinct_primes": len(proved_primes), **decompose(M, run["W_u_basis"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    evidence = json.loads(gzip.decompress(args.evidence.read_bytes()))
    result = verify(evidence)
    if args.check:
        summary = json.loads(SUMMARY.read_text())
        require(summary["evidence_sha256"] == digest(args.evidence), "evidence hash mismatch")
        for path, expected_hash in summary["source_hashes"].items():
            require(digest(ROOT / path) == expected_hash, f"stale source hash: {path}")
        require(summary["arithmetic"] == result, "stale arithmetic summary")
    print(json.dumps(result, indent=2) if not args.check else
          f"FIXEDCUBICCT|status=PASS|entries=153|rank={result['pairing_rank']}"
          f"|radical_dimension={result['restricted_radical_dimension']}")


if __name__ == "__main__":
    main()
