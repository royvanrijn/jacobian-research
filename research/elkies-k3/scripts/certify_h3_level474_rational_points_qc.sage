#!/usr/bin/env sage
"""Certify all rational points on the level-474 H3 base without Magma.

For

    C: y^2 = -3*x^6 + 22*x^4 - 19*x^2 + 64,

the two degree-two elliptic quotients both have unconditional rank one and
trivial rational torsion.  This verifier uses the pinned Sage-10.9-compatible
revision of Bianchi--Padurariu bielliptic quadratic Chabauty at 11 and 41.
It translates every unrecognised p-adic root into quotient Mordell--Weil
coordinates modulo p^4 and eliminates their CRT classes with exact good
reduction at the listed finite primes.

The script is a current, independently implemented global proof route.  It
does not alter the historic Magma certificate or pretend to replay it.
"""

import argparse
import hashlib
import json
import sys
import tempfile
import urllib.request
from pathlib import Path

RESEARCH_ROOT = Path.cwd().resolve()
if not (RESEARCH_ROOT / "elkies-k3/AGENTS.md").is_file():
    raise RuntimeError("run this verifier from the research directory")
OUTPUT = (
    RESEARCH_ROOT
    / "artifacts/generated-results/elkies-k3-h3-level474-rational-points-qc-sage109.json"
)
UPSTREAM_COMMIT = "84af22e9cd1244c3d44e3c083073b44b8d728159"
UPSTREAM_URL = (
    "https://raw.githubusercontent.com/jbalakrishnan/QC_bielliptic/"
    + UPSTREAM_COMMIT + "/qc_g2_bielliptic.sage"
)
UPSTREAM_SHA256 = "8ed7c1d61282d5da3c46cf83ca5b315d4a0a6ad2db8674f3623a54d3b41d3210"
PATCHED_SHA256 = "ac350585ecdd795c1ff24ea6e257d4f49cdc43f84482099527e852bd518a423c"
SIEVE_BOUND = 1987


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    return sha256_bytes(path.read_bytes())


def prepare_upstream(source_path, cache_dir):
    if source_path is None:
        original_path = cache_dir / ("qc_g2_bielliptic_" + UPSTREAM_COMMIT + ".sage")
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
        raise AssertionError("quadratic-Chabauty source is not the pinned upstream revision")

    source = data.decode("utf-8")
    old_infinity_removal = """    if GF(p)(a6).is_square() == False:\n        D.remove(Hp(0,1,0))\n"""
    new_infinity_removal = """    if GF(p)(a6).is_square() == False:\n        # Sage 10.9 no longer constructs the formal weighted-projective\n        # placeholder Hp(0,1,0) when the two points at infinity are not\n        # rational.  Remove infinity by coordinates instead.\n        D = [P for P in D if P[2] != 0]\n"""
    old_infinity_test = "            if P == H(0, 1, 0) or P == HK(0, 1, 0):\n"
    new_infinity_test = "            if P[2] == 0:\n"
    if source.count(old_infinity_removal) != 1 or source.count(old_infinity_test) != 1:
        raise AssertionError("pinned compatibility-patch context changed")
    source = source.replace(old_infinity_removal, new_infinity_removal)
    source = source.replace(old_infinity_test, new_infinity_test)
    patched_path = cache_dir / ("qc_g2_bielliptic_" + UPSTREAM_COMMIT + "_sage10_9.sage")
    patched_path.write_text(source)
    if sha256_file(patched_path) != PATCHED_SHA256:
        raise AssertionError("compatibility-patched source hash changed")
    return patched_path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--upstream-source", type=Path)
parser.add_argument(
    "--cache-dir",
    type=Path,
    default=RESEARCH_ROOT / "artifacts/local/h3-provenance-20260913/qc-bielliptic",
)
cli_args = sys.argv[1:]
if cli_args and cli_args[0].endswith("certify_h3_level474_rational_points_qc.sage"):
    cli_args = cli_args[1:]
args = parser.parse_args(cli_args)
args.cache_dir.mkdir(parents=True, exist_ok=True)
patched_source = prepare_upstream(args.upstream_source, args.cache_dir)
load(str(patched_source))

R.<x> = PolynomialRing(QQ)
f = -3*x^6 + 22*x^4 - 19*x^2 + 64
H = HyperellipticCurve(f)
a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6^2])
E2 = EllipticCurve([0, a2, 0, a0*a4, a0^2*a6])
for E in [E1, E2]:
    assert E.rank_bounds() == (1, 1)
    assert E.gens(proof=True)
    assert E.gens_certain()
    assert E.torsion_subgroup().order() == 1
G1 = E1.gens(proof=True)[0]
G2 = E2.gens(proof=True)[0]
P0 = H(0, 8)


def quotient_images(P):
    xx, yy = QQ(P[0]), QQ(P[1])
    image1 = E1(a6*xx^2, a6*yy)
    image2 = E2(0) if xx == 0 else E2(a0/xx^2, a0*yy/xx^3)
    return image1, image2


# These exact images fix the quotient maps, generator orientation and the
# base-point offset used below.  They are not a claim that the displayed
# product-coordinate normalizers form a saturated basis of J(QQ).
for point, n1, n2 in [
    (H(0, 8), 6, 0),
    (H(1, 8), -1, 1),
    (H(-1, 8), -1, -1),
    (H(QQ(13)/7, QQ(4016)/343), 9, -3),
    (H(-QQ(13)/7, QQ(4016)/343), 9, 3),
]:
    assert quotient_images(point) == (n1*G1, n2*G2)

# ``coefficients_mod_pN_v2`` inverts formal quotient logarithms.  Supplying
# these product-coordinate normalizers yields necessary E1 x E2 coordinates
# for every rational C-point; the sieve never assumes they lift to a basis of
# the genus-two Jacobian.
coordinate_normalizers = [[G1, E2(0)], [E1(0), G2]]

PRIMES = [11, 41]
PRECISION = 20
COEFFICIENT_EXPONENT = 4
expected_orbits = {
    (QQ(0), QQ(8)),
    (QQ(1), QQ(8)),
    (QQ(13)/7, QQ(4016)/343),
}

def orbit_key(P):
    return (abs(QQ(P[0])), abs(QQ(P[1])))

def compact_qc(p):
    assert E1.has_good_reduction(p) and E2.has_good_reduction(p)
    assert E1.is_ordinary(p) and E2.is_ordinary(p)
    rational, other = quadratic_chabauty_bielliptic(
        f, p, PRECISION, up_to_auto=True, omega_info=True
    )
    assert len(rational) == len(other) == 72
    found = {orbit_key(P) for row in rational for P in row}
    assert found == expected_orbits
    flat = [P for row in other for P in row]
    raw = coefficients_mod_pN_v2(
        f, flat, coordinate_normalizers, P0, p, COEFFICIENT_EXPONENT
    )
    assert len(raw) == len(flat)
    modulus = p^COEFFICIENT_EXPONENT
    compact, index = [], 0
    for row in other:
        coefficients = raw[index:index + len(row)]
        index += len(row)
        packed = []
        for coefficient in coefficients:
            assert coefficient[0].precision_absolute() >= COEFFICIENT_EXPONENT
            assert coefficient[1].precision_absolute() >= COEFFICIENT_EXPONENT
            packed.append((
                ZZ(coefficient[0].lift()) % modulus,
                ZZ(coefficient[1].lift()) % modulus,
            ))
        compact.append(packed)
    assert index == len(raw)
    return {
        'modulus': modulus,
        'rational_orbit_count': len(found),
        'other_counts': [len(row) for row in other],
        'compact': compact,
    }

def close_under_automorphisms(rows, modulus):
    out = []
    for row in rows:
        expanded = set()
        for a, b in row:
            # With coefficients relative to P0, x -> -x fixes a and negates b;
            # y -> -y sends (a,b) to (-12-a,-b), since phi1(P0)=6G1.
            expanded.update([
                (a % modulus, b % modulus),
                (a % modulus, (-b) % modulus),
                ((-12-a) % modulus, (-b) % modulus),
                ((-12-a) % modulus, b % modulus),
            ])
        out.append(sorted(expanded))
    return out

def point_key(P):
    if P.is_zero():
        return ('O',)
    return (int(ZZ(P[0])), int(ZZ(P[1])))

def subgroup_logs(E, G):
    order = ZZ(G.order())
    logs, P = {}, E(0)
    for n in range(order):
        logs[point_key(P)] = ZZ(n)
        P += G
    assert P.is_zero()
    return order, logs

def allowed_pairs(q, modulus):
    F = GF(q)
    E1q, E2q = E1.change_ring(F), E2.change_ring(F)
    G1q, G2q = E1q(G1), E2q(G2)
    order1, logs1 = subgroup_logs(E1q, G1q)
    order2, logs2 = subgroup_logs(E2q, G2q)
    d1, d2 = gcd(modulus, order1), gcd(modulus, order2)
    if d1 == d2 == 1:
        return d1, d2, {(0, 0)}, order1, order2
    allowed = set()
    for xx in F:
        rhs = F(a6)*xx^6 + F(a4)*xx^4 + F(a2)*xx^2 + F(a0)
        if not rhs.is_square():
            continue
        y_values = [F(0)] if rhs == 0 else [rhs.sqrt(), -rhs.sqrt()]
        for yy in y_values:
            image1 = E1q(F(a6)*xx^2, F(a6)*yy)
            image2 = E2q(0) if xx == 0 else E2q(F(a0)/xx^2, F(a0)*yy/xx^3)
            n1, n2 = logs1.get(point_key(image1)), logs2.get(point_key(image2))
            if n1 is not None and n2 is not None:
                allowed.add(((n1 - 6) % d1, n2 % d2))
    if F(a6).is_square():
        root = F(a6).sqrt()
        for sign in [1, -1]:
            image1 = E1q(0)
            image2 = E2q(F(0), F(sign*a0)*root)
            n1, n2 = logs1.get(point_key(image1)), logs2.get(point_key(image2))
            if n1 is not None and n2 is not None:
                allowed.add(((n1 - 6) % d1, n2 % d2))
    return d1, d2, allowed, order1, order2

records = {p: compact_qc(p) for p in PRIMES}
assert records[11]['other_counts'] == [5,5,5,5,6,6,6,6,1,1,1,1,8,6,6,6,1,1,1,1,6,6,6,6,4,4,4,4,1,1,1,1,8,6,6,6,4,4,4,4,7,7,7,7,6,6,6,6,9,7,9,9,7,7,7,7,6,6,6,6,4,4,4,4,9,9,9,9,5,5,5,5]
assert records[41]['other_counts'] == [16,9,16,12,19,14,13,7,11,10,7,14,10,15,19,11,9,14,16,5,10,16,7,14,11,12,11,16,16,11,11,10,7,9,9,18,10,12,11,16,9,5,8,16,11,14,15,7,11,5,12,13,8,10,15,11,15,18,14,9,11,7,16,16,11,11,11,19,16,11,16,7]
rows11 = close_under_automorphisms(records[11]['compact'], records[11]['modulus'])
rows41 = close_under_automorphisms(records[41]['compact'], records[41]['modulus'])
M11, M41 = records[11]['modulus'], records[41]['modulus']
M = M11*M41
candidates = set()
for omega in range(72):
    for c11 in rows11[omega]:
        for c41 in rows41[omega]:
            candidates.add((
                omega,
                CRT(c11[0], c41[0], M11, M41) % M,
                CRT(c11[1], c41[1], M11, M41) % M,
            ))
initial_candidate_count = len(candidates)

bad_primes = set(ZZ(f.discriminant()).prime_divisors()) | set(PRIMES)
sieve_trace = []
for q in prime_range(5, SIEVE_BOUND + 1):
    if q in bad_primes:
        continue
    d1, d2, allowed, order1, order2 = allowed_pairs(q, M)
    if d1 == d2 == 1:
        continue
    before = len(candidates)
    candidates = {
        c for c in candidates
        if (c[1] % d1, c[2] % d2) in allowed
    }
    if len(candidates) < before:
        sieve_trace.append({
            'prime': int(q), 'generator_orders': [int(order1), int(order2)],
            'gcds_with_crt_modulus': [int(d1), int(d2)],
            'allowed_reduction_pairs': len(allowed),
            'before': before, 'after': len(candidates),
        })
    if not candidates:
        break
assert not candidates
assert sieve_trace[-1]["prime"] == SIEVE_BOUND
payload = {
    "schema_version": 1,
    "artifact_kind": "exact_genus_two_rational_points",
    "status": "PASS_H3_LEVEL474_RATIONAL_POINTS_QC_SAGE109",
    "scope_boundary": (
        "This certifies C(QQ) for the displayed level-474 H3 base. It is an "
        "independent current Sage quadratic-Chabauty and finite "
        "Mordell--Weil-sieve proof, not a replay of the historic Magma route "
        "or a construction of a downstream rootless MW17 equation."
    ),
    "curve": {
        "model": "y^2 = -3*x^6 + 22*x^4 - 19*x^2 + 64",
        "published_model": "Y^2 = -27*X^6 + 198*X^4 - 171*X^2 + 576; Y=3*y",
        "rational_points_complete": True,
        "rational_points": [
            "(0,+/-8)", "(+/-1,+/-8)",
            "(+/-13/7,+/-4016/343)"
        ],
        "rational_points_at_infinity": 0,
        "infinity_exclusion": "the leading coefficient -3 is not a rational square",
    },
    "elliptic_quotients": [
        {
            "model": "y^2=x^3+22*x^2+57*x+576",
            "rank_bounds": [1, 1],
            "torsion_order": 1,
            "generator": str(G1),
            "conductor": int(E1.conductor()),
        },
        {
            "model": "y^2=x^3-19*x^2+1408*x-12288",
            "rank_bounds": [1, 1],
            "torsion_order": 1,
            "generator": str(G2),
            "conductor": int(E2.conductor()),
        },
    ],
    "quadratic_chabauty": {
        "method_source": "Bianchi--Padurariu QC_bielliptic",
        "source_url": UPSTREAM_URL,
        "source_commit": UPSTREAM_COMMIT,
        "source_sha256": UPSTREAM_SHA256,
        "compatibility_patched_sha256": PATCHED_SHA256,
        "compatibility_patch_boundary": (
            "Only the Sage-10.9 point-at-infinity representation tests are "
            "adapted; the p-adic height and root computations are unchanged."
        ),
        "working_precision": int(PRECISION),
        "coefficient_modulus_exponent": int(COEFFICIENT_EXPONENT),
        "primes": [
            {
                "prime": int(p),
                "omega_count": 72,
                "recognized_rational_orbits": int(records[p]["rational_orbit_count"]),
                "unrecognised_padic_orbit_counts": [int(v) for v in records[p]["other_counts"]],
            }
            for p in PRIMES
        ],
    },
    "mordell_weil_sieve": {
        "quotient_coordinate_origin": "P0=(0,8), with first quotient image 6*G1 and second quotient image O",
        "coordinate_lattice": "Z*G1 x Z*G2; every rational C-point maps here because both quotients are rank-one and torsion-free",
        "automorphism_actions": [
            "x -> -x: (a,b) -> (a,-b)",
            "y -> -y: (a,b) -> (-12-a,-b)",
        ],
        "crt_modulus": int(M),
        "initial_candidate_cosets": int(initial_candidate_count),
        "sieve_prime_search_bound": SIEVE_BOUND,
        "trace": sieve_trace,
        "remaining_candidate_cosets": 0,
    },
    "theorem_inputs": [
        "Bianchi--Padurariu bielliptic quadratic-Chabauty finite-height-value theorem and pinned implementation",
        "unconditional rank-one, torsion-free certificates for both elliptic quotients",
        "exact formal-log quotient coordinates modulo 11^4 and 41^4",
        "finite good-reduction Mordell--Weil sieve through the final eliminating prime 1987",
    ],
    "reproduce": "cd research && sage elkies-k3/scripts/certify_h3_level474_rational_points_qc.sage --check",
}
rendered = json.dumps(payload, indent=2, sort_keys=True, default=int) + "\n"
args.output.parent.mkdir(parents=True, exist_ok=True)
if args.check:
    if not args.output.is_file():
        raise FileNotFoundError(args.output)
    if args.output.read_text() != rendered:
        raise AssertionError("generated rational-point artifact changed")
    print(json.dumps({
        "status": "PASS_H3_LEVEL474_RATIONAL_POINTS_QC_SAGE109_CHECK",
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "remaining_candidate_cosets": 0,
    }, sort_keys=True, default=int))
else:
    args.output.write_text(rendered)
    print(json.dumps({
        "status": "WROTE_H3_LEVEL474_RATIONAL_POINTS_QC_SAGE109_CERTIFICATE",
        "output": str(args.output),
        "sha256": sha256_file(args.output),
        "remaining_candidate_cosets": 0,
    }, sort_keys=True, default=int))
