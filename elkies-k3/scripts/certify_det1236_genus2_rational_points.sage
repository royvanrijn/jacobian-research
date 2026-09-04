#!/usr/bin/env sage
"""Certify all rational points on the determinant-1236 genus-two quotient.

This is a quadratic-Chabauty plus Mordell--Weil-sieve certificate for

    B: y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9.

It deliberately proves only B(QQ).  It does not construct or decide rational
points on the marked double cover X_0^6(103)/<w_618> -> B.
"""

import argparse
import hashlib
import json
import sys
import tempfile
import urllib.request
from pathlib import Path


ROOT = Path.cwd().resolve()
if not (ROOT / "elkies-k3/AGENTS.md").is_file():
    raise RuntimeError("run this certificate from the repository root")
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / "elkies-k3-det1236-genus2-rational-points-v1.json"

UPSTREAM_COMMIT = "84af22e9cd1244c3d44e3c083073b44b8d728159"
UPSTREAM_URL = (
    "https://raw.githubusercontent.com/jbalakrishnan/QC_bielliptic/"
    + UPSTREAM_COMMIT
    + "/qc_g2_bielliptic.sage"
)
UPSTREAM_SHA256 = "8ed7c1d61282d5da3c46cf83ca5b315d4a0a6ad2db8674f3623a54d3b41d3210"
PATCHED_SHA256 = "ac350585ecdd795c1ff24ea6e257d4f49cdc43f84482099527e852bd518a423c"
WORKING_PRECISION = 20
COEFFICIENT_EXPONENT = 4
QC_PRIMES = [7, 11]
SIEVE_BOUND = 7000


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    return sha256_bytes(path.read_bytes())


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def prepare_upstream(source_path, cache_dir):
    if source_path is None:
        original_path = cache_dir / (
            "qc_g2_bielliptic_" + UPSTREAM_COMMIT + ".sage"
        )
        if not original_path.is_file():
            with urllib.request.urlopen(UPSTREAM_URL, timeout=60) as response:
                data = response.read()
            if sha256_bytes(data) != UPSTREAM_SHA256:
                raise AssertionError("downloaded quadratic-Chabauty source hash changed")
            original_path.write_bytes(data)
    else:
        original_path = source_path.resolve()
    data = original_path.read_bytes()
    if sha256_bytes(data) != UPSTREAM_SHA256:
        raise AssertionError(
            "quadratic-Chabauty source is not the pinned upstream revision"
        )

    text = data.decode("utf-8")
    old_infinity_removal = """    if GF(p)(a6).is_square() == False:\n        D.remove(Hp(0,1,0))\n"""
    new_infinity_removal = """    if GF(p)(a6).is_square() == False:\n        # Sage 10.9 no longer constructs the formal weighted-projective\n        # placeholder Hp(0,1,0) when the two points at infinity are not\n        # rational.  Remove infinity by coordinates instead.\n        D = [P for P in D if P[2] != 0]\n"""
    old_infinity_test = (
        "            if P == H(0, 1, 0) or P == HK(0, 1, 0):\n"
    )
    new_infinity_test = "            if P[2] == 0:\n"
    if text.count(old_infinity_removal) != 1 or text.count(old_infinity_test) != 1:
        raise AssertionError("pinned compatibility-patch context changed")
    text = text.replace(old_infinity_removal, new_infinity_removal)
    text = text.replace(old_infinity_test, new_infinity_test)
    patched_path = cache_dir / (
        "qc_g2_bielliptic_" + UPSTREAM_COMMIT + "_sage10_9.sage"
    )
    patched_path.write_text(text)
    if sha256_file(patched_path) != PATCHED_SHA256:
        raise AssertionError("compatibility-patched source hash changed")
    return patched_path


parser = argparse.ArgumentParser()
parser.add_argument("--check", action="store_true")
parser.add_argument("--fresh", action="store_true")
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--upstream-source", type=Path)
parser.add_argument(
    "--cache-dir",
    type=Path,
    default=Path(tempfile.gettempdir()) / "det1236-qc-bielliptic-v1",
)
cli_args = sys.argv[1:]
if cli_args and cli_args[0].endswith(
    "certify_det1236_genus2_rational_points.sage"
):
    cli_args = cli_args[1:]
args = parser.parse_args(cli_args)
args.cache_dir.mkdir(parents=True, exist_ok=True)
patched_source = prepare_upstream(args.upstream_source, args.cache_dir)
load(str(patched_source))


R.<x> = PolynomialRing(QQ)
f = 1944*x^6 + 441*x^4 - 90*x^2 + 9
H = HyperellipticCurve(f)
a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6^2])
E2 = EllipticCurve([0, a2, 0, a0*a4, a0^2*a6])
G1 = E1.gens(proof=True)[0]
G2 = E2.gens(proof=True)[0]
assert E1.gens_certain() and E2.gens_certain()
assert E1.rank_bounds() == (1, 1) and E2.rank_bounds() == (1, 1)
assert E1.torsion_subgroup().order() == 1
assert E2.torsion_subgroup().order() == 1
assert E1.minimal_model().cremona_label() == "618f1"
assert E2.minimal_model().cremona_label() == "618e1"

P0 = H(0, 3)
assert E1(0, a6*3) == 3*G1
im_divisors = [[G1, E2(0)], [E1(0), G2]]

known_affine_points = []
for xx, yy in [
    (0, 3),
    (1, 48),
    (QQ(1)/3, QQ(8)/3),
    (QQ(1)/5, QQ(312)/125),
]:
    x_values = [QQ(xx)] if xx == 0 else [QQ(xx), -QQ(xx)]
    for xxx in x_values:
        for yyy in [QQ(yy), -QQ(yy)]:
            known_affine_points.append(H(xxx, yyy))
assert len(set(known_affine_points)) == 14
assert not QQ(a6).is_square()


def quotient_images(P):
    xx, yy = QQ(P[0]), QQ(P[1])
    image1 = E1(a6*xx^2, a6*yy) - E1(0, a6*3)
    image2 = (
        E2(0)
        if xx == 0
        else E2(a0/xx^2, a0*yy/xx^3)
    )
    return image1, image2


# Check the coefficient normalization and both automorphism formulas exactly,
# rather than relying on the later p-adic recovery to infer them.
for xx, yy, n1, n2 in [
    (0, 3, 0, 0),
    (1, 48, 1, 1),
    (QQ(1)/3, QQ(8)/3, -4, -2),
    (QQ(1)/5, QQ(312)/125, -13, 3),
]:
    P = H(xx, yy)
    image1, image2 = quotient_images(P)
    assert image1 == n1*G1 and image2 == n2*G2
    x_image1, x_image2 = quotient_images(H(-xx, yy))
    y_image1, y_image2 = quotient_images(H(xx, -yy))
    assert x_image1 == image1 and x_image2 == -image2
    assert y_image1 == (-6-n1)*G1 and y_image2 == -image2


def compute_qc_coefficients(p):
    cache = args.cache_dir / (
        "coeffs_%s_p%s_prec%s_exp%s.sobj"
        % (UPSTREAM_COMMIT[:12], p, WORKING_PRECISION, COEFFICIENT_EXPONENT)
    )
    if cache.is_file() and not args.fresh:
        record = load(str(cache))
        if record.get("schema") != 1 or record.get("prime") != p:
            raise AssertionError("invalid quadratic-Chabauty cache")
        return record

    rational, other = quadratic_chabauty_bielliptic(
        f,
        p,
        WORKING_PRECISION,
        up_to_auto=True,
        omega_info=True,
    )
    assert len(rational) == len(other) == 621
    flat_points = [P for points in other for P in points]
    flat_coeffs = coefficients_mod_pN_v2(
        f,
        flat_points,
        im_divisors,
        P0,
        p,
        COEFFICIENT_EXPONENT,
    )
    assert len(flat_coeffs) == len(flat_points)
    for coeffs in flat_coeffs:
        assert coeffs[0].precision_absolute() >= COEFFICIENT_EXPONENT
        assert coeffs[1].precision_absolute() >= COEFFICIENT_EXPONENT

    compact = []
    pos = 0
    modulus = p^COEFFICIENT_EXPONENT
    for points in other:
        coeffs = flat_coeffs[pos:pos + len(points)]
        pos += len(points)
        compact.append([
            (ZZ(c[0].lift()) % modulus, ZZ(c[1].lift()) % modulus)
            for c in coeffs
        ])
    assert pos == len(flat_coeffs)

    rational_strings = sorted(
        str(P) for points in rational for P in points
    )
    record = {
        "schema": 1,
        "prime": p,
        "omega_count": len(compact),
        "orbit_fake_count": sum(map(len, compact)),
        "nonempty_omega_count": sum(bool(v) for v in compact),
        "recognized_rational_orbit_representatives": rational_strings,
        "coefficients": compact,
    }
    save(record, str(cache))
    return record


qc_records = {p: compute_qc_coefficients(p) for p in QC_PRIMES}
expected_recognized_orbits = {
    7: [
        "(-1/3 : -8/3 : 1)",
        "(1 : -48 : 1)",
        "(1/5 : 312/125 : 1)",
    ],
    11: [
        "(-1/5 : 312/125 : 1)",
        "(1 : 48 : 1)",
        "(1/3 : -8/3 : 1)",
    ],
}
for p in QC_PRIMES:
    # Up to (x,y) -> (+/-x,+/-y), these are exactly the three known
    # nonzero-x orbits.  The x=0 orbit is an exceptional coordinate of the
    # second elliptic quotient and is accounted for separately above.
    assert qc_records[p]["recognized_rational_orbit_representatives"] == (
        expected_recognized_orbits[p]
    )
fake_coeffs = {}
for p in QC_PRIMES:
    compact = qc_records[p]["coefficients"]
    modulus = p^COEFFICIENT_EXPONENT
    full_compact = []
    for coeffs in compact:
        expanded = set()
        for n1, n2 in coeffs:
            # x -> -x: (n1,n2) -> (n1,-n2).
            # y -> -y: (n1,n2) -> (-6-n1,-n2), since phi_1(P0)=3G1.
            expanded.update([
                (n1 % modulus, n2 % modulus),
                (n1 % modulus, (-n2) % modulus),
                ((-6-n1) % modulus, (-n2) % modulus),
                ((-6-n1) % modulus, n2 % modulus),
            ])
        full_compact.append(sorted(expanded))
    fake_coeffs[p] = full_compact

m7 = 7^COEFFICIENT_EXPONENT
m11 = 11^COEFFICIENT_EXPONENT
M = m7*m11
candidates = set()
for omega in range(621):
    for c7 in fake_coeffs[7][omega]:
        for c11 in fake_coeffs[11][omega]:
            candidates.add((
                omega,
                CRT(c7[0], c11[0], m7, m11) % M,
                CRT(c7[1], c11[1], m7, m11) % M,
            ))
initial_candidate_count = len(candidates)
initial_nonempty_omega = len(set(c[0] for c in candidates))
assert initial_candidate_count == 231760
assert initial_nonempty_omega == 405


def point_key(P):
    if P.is_zero():
        return ("O",)
    return (int(ZZ(P[0])), int(ZZ(P[1])))


def subgroup_logs(E, G):
    order = ZZ(G.order())
    logs = {}
    P = E(0)
    for n in range(order):
        logs[point_key(P)] = ZZ(n)
        P += G
    assert P.is_zero()
    return order, logs


def allowed_pairs(q):
    F = GF(q)
    E1q = E1.change_ring(F)
    E2q = E2.change_ring(F)
    G1q = E1q(G1)
    G2q = E2q(G2)
    order1, logs1 = subgroup_logs(E1q, G1q)
    order2, logs2 = subgroup_logs(E2q, G2q)
    d1, d2 = gcd(M, order1), gcd(M, order2)
    if d1 == d2 == 1:
        return d1, d2, {(0, 0)}, order1, order2

    base1 = E1q(F(0), F(a6*3))
    allowed = set()
    for xx in F:
        rhs = F(a6)*xx^6 + F(a4)*xx^4 + F(a2)*xx^2 + F(a0)
        if not rhs.is_square():
            continue
        y_values = [F(0)] if rhs == 0 else [rhs.sqrt(), -rhs.sqrt()]
        for yy in y_values:
            image1 = E1q(F(a6)*xx^2, F(a6)*yy) - base1
            image2 = (
                E2q(0)
                if xx == 0
                else E2q(F(a0)/xx^2, F(a0)*yy/xx^3)
            )
            n1 = logs1.get(point_key(image1))
            n2 = logs2.get(point_key(image2))
            if n1 is not None and n2 is not None:
                allowed.add((n1 % d1, n2 % d2))

    if F(a6).is_square():
        root_a6 = F(a6).sqrt()
        for sign in [1, -1]:
            image1 = -base1
            image2 = E2q(F(0), F(sign*a0)*root_a6)
            n1 = logs1.get(point_key(image1))
            n2 = logs2.get(point_key(image2))
            if n1 is not None and n2 is not None:
                allowed.add((n1 % d1, n2 % d2))
    return d1, d2, allowed, order1, order2


bad_primes = set(ZZ(f.discriminant()).prime_divisors()) | set(QC_PRIMES)
sieve_trace = []
for q in prime_range(5, SIEVE_BOUND + 1):
    if q in bad_primes:
        continue
    d1, d2, allowed, order1, order2 = allowed_pairs(q)
    if d1 == d2 == 1:
        continue
    before = len(candidates)
    candidates = {
        c for c in candidates if (c[1] % d1, c[2] % d2) in allowed
    }
    if len(candidates) < before:
        sieve_trace.append({
            "prime": int(q),
            "generator_orders": [int(order1), int(order2)],
            "gcds_with_crt_modulus": [int(d1), int(d2)],
            "allowed_reduction_pairs": len(allowed),
            "before": before,
            "after": len(candidates),
        })
    if not candidates:
        break

assert not candidates
assert sieve_trace[-1]["prime"] == 6599

payload = {
    "schema_version": 1,
    "status": "PASS_DET1236_GENUS2_RATIONAL_POINTS",
    "scope_boundary": (
        "This certifies B(QQ), not rational lifts to the marked genus-six "
        "curve C_1236. The determinant-1236 arithmetic gate remains unresolved."
    ),
    "curve": {
        "label": "B = X_0^6(103)/<w_2,w_309>",
        "model": "y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9",
        "rational_points_complete": True,
        "rational_point_count": 14,
        "rational_points": [
            "(0,+/-3)",
            "(+/-1,+/-48)",
            "(+/-1/3,+/-8/3)",
            "(+/-1/5,+/-312/125)",
        ],
        "rational_points_at_infinity": 0,
    },
    "elliptic_quotients": [
        {
            "cremona_label": "618f1",
            "rank_bounds": [1, 1],
            "torsion_order": 1,
            "generator_on_working_model": [216, -5184],
        },
        {
            "cremona_label": "618e1",
            "rank_bounds": [1, 1],
            "torsion_order": 1,
            "generator_on_working_model": [9, 432],
        },
    ],
    "quadratic_chabauty": {
        "method_source": "Balakrishnan QC_bielliptic",
        "source_url": UPSTREAM_URL,
        "source_commit": UPSTREAM_COMMIT,
        "source_sha256": UPSTREAM_SHA256,
        "compatibility_patched_sha256": PATCHED_SHA256,
        "compatibility_patch_boundary": (
            "Only two Sage-10.9 point-at-infinity representation tests are "
            "changed; the p-adic height and root computations are unchanged."
        ),
        "working_precision": WORKING_PRECISION,
        "coefficient_modulus_exponent": COEFFICIENT_EXPONENT,
        "rational_orbit_accounting": (
            "At each prime the recognized rational roots are exactly the "
            "three known nonzero-x orbits. The exceptional x=0 orbit "
            "(0,+/-3) is checked directly; the nonsquare leading coefficient "
            "excludes rational points at infinity."
        ),
        "primes": [
            {
                "prime": p,
                "omega_count": qc_records[p]["omega_count"],
                "orbit_mock_points": qc_records[p]["orbit_fake_count"],
                "full_mock_points_after_symmetry": sum(
                    map(len, fake_coeffs[p])
                ),
                "nonempty_omega_count": qc_records[p]["nonempty_omega_count"],
                "recognized_rational_orbit_representatives": qc_records[p][
                    "recognized_rational_orbit_representatives"
                ],
            }
            for p in QC_PRIMES
        ],
    },
    "mordell_weil_sieve": {
        "coordinate_lattice": "Z*G_618f1 x Z*G_618e1",
        "automorphism_actions": [
            "x -> -x: (n1,n2) -> (n1,-n2)",
            "y -> -y: (n1,n2) -> (-6-n1,-n2)",
        ],
        "crt_modulus": int(M),
        "initial_candidate_cosets": initial_candidate_count,
        "initial_nonempty_omega_count": initial_nonempty_omega,
        "sieve_prime_search_bound": SIEVE_BOUND,
        "trace": sieve_trace,
        "remaining_candidate_cosets": 0,
    },
    "theorem_inputs": [
        "Bianchi's bielliptic quadratic-Chabauty finite-height-value theorem and implementation",
        "the formal logarithm recovery of Mordell--Weil coefficients modulo p^4",
        "the Mordell--Weil sieve on reductions of the two saturated rank-one elliptic quotients",
    ],
    "reproduce": (
        "sage -- elkies-k3/scripts/"
        "certify_det1236_genus2_rational_points.sage --fresh"
    ),
}

rendered = json.dumps(payload, indent=2, sort_keys=True, default=int) + "\n"
if args.check:
    if not args.output.is_file():
        raise FileNotFoundError(args.output)
    if args.output.read_text() != rendered:
        raise AssertionError("generated rational-point artifact changed")
    print(json.dumps({
        "status": "PASS_DET1236_GENUS2_RATIONAL_POINTS_CHECK",
        "output": relative(args.output),
        "sha256": sha256_file(args.output),
        "rational_point_count": 14,
        "remaining_candidate_cosets": 0,
    }, sort_keys=True, default=int))
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(json.dumps({
        "status": "WROTE_DET1236_GENUS2_RATIONAL_POINTS_CERTIFICATE",
        "output": relative(args.output),
        "sha256": sha256_file(args.output),
        "rational_point_count": 14,
        "remaining_candidate_cosets": 0,
    }, sort_keys=True, default=int))
