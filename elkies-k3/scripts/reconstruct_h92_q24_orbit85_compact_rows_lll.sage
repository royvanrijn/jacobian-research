#!/usr/bin/env sage -python
"""
Projective LLL recovery of the q24/orbit85 resolved-RR kernel rows.

The modular signatures store a canonical 2x56 RREF plane in the coefficient
frame

    A(U) / Z(U)^2  +  B(U) / Z(U) * m

with deg(A)<=40 and deg(B)<=14.  The smooth q24 collision condition gives

    A*X - B*Y = C*Z^2.

Following the q32 compact-row precedent, this script reconstructs the same
2-plane in the compact intrinsic coordinates

    (coeffs(A), coeffs(B*Z), coeffs(C))

of dimensions 41 + 39 + 45 = 125.  It uses all but the final prime to build a
CRT lattice for each projective row and reserves the final prime as a held-out
check.

A PASS output is a held-out validated QQ row basis for the q24/orbit85
resolved-RR kernel.  It is still not the final characteristic-zero D12 equation
certificate; the exact RR replay and neighbor compilation remain separate.
"""

import argparse
import json
import math
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

ACCEPTED_SIGNATURE_STATUSES = {
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
}


def available_signature_primes():
    primes = []
    for path in sorted(LOCAL.glob("q24-orbit85-d12-signature-mod-*.json")):
        try:
            primes.append(ZZ(path.stem.split("-")[-1]))
        except (TypeError, ValueError):
            continue
    return primes


def parse_primes(raw):
    if not raw:
        return available_signature_primes()
    return [ZZ(token) for token in str(raw).replace(",", " ").split()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--primes",
    help="Comma/space separated prime list. Defaults to all q24/orbit85 signatures.",
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

selected_primes = parse_primes(args.primes)
records = []
for p in selected_primes:
    spath = LOCAL / f"q24-orbit85-d12-signature-mod-{p}.json"
    qpath = LOCAL / f"q24-degree46-direct-global-mod-{p}.json"
    if not spath.exists():
        raise SystemExit(f"missing signature artifact: {spath}")
    if not qpath.exists():
        raise SystemExit(f"missing q24 direct bridge artifact: {qpath}")
    sig = json.loads(spath.read_text())
    q24 = json.loads(qpath.read_text())
    assert sig["status"] in ACCEPTED_SIGNATURE_STATUSES
    assert ZZ(sig["prime"]) == p
    assert sig["source_neighbor"] == {
        "q": 24,
        "orbit": 85,
        "source": "D13/MW4",
        "child": "D12/MW5",
    }
    assert sig["rr"] == {
        "ambient": 56,
        "collision_rank": 48,
        "post_collision": 8,
        "resolved_rank": 6,
        "kernel": 2,
        "geometric_fibre_twist": -8,
    }
    assert sig["plane_pivots"] == [0, 1]
    assert q24["status"] == "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
    records.append((p, sig, q24))

records.sort(key=lambda item: int(item[0]))
if len(records) < 3:
    raise SystemExit("need at least three complete modular signatures")


DA, DY, D0 = 40, 38, 44
NC = 125


def compact_plane(p, sig, q24):
    F = GF(p)
    R = PolynomialRing(F, "U")
    U = R.gen()
    sec = q24["section_mod_p"]
    Z = R([F(v) for v in sec["Z_coefficients_low_to_high"]])
    X = R([F(v) for v in sec["X_coefficients_low_to_high"]])
    Y = R([F(v) for v in sec["Y_coefficients_low_to_high"]])

    assert Z.degree() == 24
    P = matrix(F, sig["plane_rref_2x56"])
    rows = []
    for row in P.rows():
        A = sum(row[i] * U**i for i in range(41))
        B = sum(row[41 + i] * U**i for i in range(15))
        C, remainder = (A * X - B * Y).quo_rem(Z**2)
        assert remainder == 0
        BY = B * Z
        assert A.degree() <= DA
        assert BY.degree() <= DY
        assert C.degree() <= D0
        rows.append(
            [A[i] for i in range(DA + 1)]
            + [BY[i] for i in range(DY + 1)]
            + [C[i] for i in range(D0 + 1)]
        )
    M = matrix(F, rows).echelon_form()
    assert M.dimensions() == (2, NC)
    assert tuple(M.pivots()) == (0, 1)
    return M


planes = [(p, compact_plane(p, sig, q24)) for p, sig, q24 in records]
print(
    "Q24O85ROWLLL_INPUT|"
    f"primes={','.join(str(p) for p, _ in planes)}|count={len(planes)}|"
    "coords=125|pivots=0,1|status=PASS",
    flush=True,
)

train = planes[:-1]
holdp, H = planes[-1]
mods = [p for p, _ in train]


def crt_scalar(vals, mods):
    x = ZZ(0)
    M = ZZ(1)
    for rr, p in zip(vals, mods):
        p = ZZ(p)
        rr = ZZ(rr) % p
        t = ((rr - x) % p) * ((M % p).inverse_mod(p)) % p
        x = (x + M * t) % (M * p)
        M *= p
    if x > M // 2:
        x -= M
    return x, M


_, MOD = crt_scalar([0] * len(mods), mods)
print(
    "Q24O85ROWLLL_MODULUS|"
    f"train={len(train)}|holdout={holdp}|"
    f"modulus_bits={MOD.nbits()}|status=PASS",
    flush=True,
)

nonpivot = list(range(2, NC))
N = len(nonpivot)
assert N == 123

results = []
for ri in range(2):
    residues = []
    for j in nonpivot:
        x, M = crt_scalar([int(P[ri, j]) for _, P in train], mods)
        assert M == MOD
        residues.append(x)

    dim = N + 1
    B = matrix(ZZ, dim, dim)
    for j in range(N):
        B[j, j] = MOD
    for j, residue in enumerate(residues):
        B[N, j] = residue
    B[N, N] = 1

    print(
        "Q24O85ROWLLL_START|"
        f"row={ri}|dimension={dim}|modulus_bits={MOD.nbits()}|status=START",
        flush=True,
    )
    R = B.LLL(delta=0.99)
    print(
        "Q24O85ROWLLL_REDUCED|"
        f"row={ri}|dimension={dim}|status=PASS",
        flush=True,
    )

    candidates = []
    for bi, v in enumerate(R.rows()):
        d = ZZ(v[N])
        if d == 0:
            continue
        nums = [ZZ(v[j]) for j in range(N)]
        g = abs(d)
        for a in nums:
            g = math.gcd(int(g), abs(int(a)))
        g = ZZ(g)
        if g > 1:
            d //= g
            nums = [a // g for a in nums]
        if d < 0:
            d = -d
            nums = [-a for a in nums]

        if d % holdp == 0:
            matches = -1
        else:
            dinv = (d % holdp).inverse_mod(holdp)
            matches = sum(
                int((a % holdp) * dinv % holdp) == int(H[ri, j])
                for a, j in zip(nums, nonpivot)
            )

        bits = max([abs(d).nbits()] + [abs(a).nbits() for a in nums])
        norm2 = sum(a * a for a in nums) + d * d
        candidates.append((matches, bits, norm2, bi, d, nums))

    if not candidates:
        raise RuntimeError(f"LLL row {ri}: no nonzero-denominator candidate")

    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    matches, bits, norm2, bi, d, nums = candidates[0]
    print(
        "Q24O85ROWLLL_BEST|"
        f"row={ri}|basis_index={bi}|heldout={matches}/{N}|"
        f"height_bits={bits}|den_bits={abs(d).nbits()}|"
        f"status={'PASS_HELDOUT' if matches == N else 'PARTIAL'}",
        flush=True,
    )

    for rank, item in enumerate(sorted(candidates, key=lambda x: (x[2], x[1]))[:5]):
        mm, bb, nn, ii, dd, aa = item
        print(
            "Q24O85ROWLLL_SHORT|"
            f"row={ri}|rank={rank}|basis_index={ii}|heldout={mm}/{N}|"
            f"height_bits={bb}|den_bits={abs(dd).nbits()}",
            flush=True,
        )

    results.append({
        "row": int(ri),
        "heldout_matches": int(matches),
        "height_bits": int(bits),
        "denominator": str(d),
        "numerators": [str(a) for a in nums],
        "pass_heldout": bool(matches == N),
    })

complete = all(row["pass_heldout"] for row in results)

payload = {
    "schema": "elkies-k3.h3-q24-orbit85-compact-row-lll.v1",
    "status": (
        "PASS_Q24_ORBIT85_COMPACT_ROWS_HELDOUT"
        if complete
        else "PARTIAL_Q24_ORBIT85_COMPACT_ROWS_LLL"
    ),
    "source_neighbor": {
        "q": 24,
        "orbit": 85,
        "source": "D13/MW4",
        "child": "D12/MW5",
    },
    "training_primes": [int(p) for p, _ in train],
    "heldout_prime": int(holdp),
    "modulus_bits": int(MOD.nbits()),
    "coordinate_degrees": [DA, DY, D0],
    "coordinate_blocks": ["A", "B_times_Z", "C_from_A_X_minus_B_Y_over_Z2"],
    "rows": results,
    "proof_boundary": (
        "Held-out validated projective recovery of compact q24/orbit85 RR "
        "kernel rows. The exact QQ RR replay and D12 neighbor compilation "
        "remain separate."
    ),
}

if complete:
    qrows = []
    for ri, result in enumerate(results):
        d = ZZ(result["denominator"])
        nums = [ZZ(value) for value in result["numerators"]]
        row = ["0"] * NC
        row[ri] = "1"
        for a, j in zip(nums, nonpivot):
            row[j] = str(a / d)
        qrows.append(row)
    payload["qq_compact_plane_rref_2x125"] = qrows

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit85-compact-row-lll.json"
)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24O85ROWLLL_RESULT|"
    f"rows_passed={sum(row['pass_heldout'] for row in results)}/2|"
    f"status={payload['status']}",
    flush=True,
)
