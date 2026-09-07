#!/usr/bin/env sage -python
"""
CRT/rational-reconstruct the q24/orbit85 resolved-RR kernel candidate.

Inputs are compact modular signatures produced by
extract_h92_q24_d12_modp_signature.sage.  The signatures store the canonical
RREF basis of the 2-dimensional resolved RR kernel in the fixed 56-column
ambient, plus the resulting binary quartic/Jacobian data.

This script reconstructs from all but the final selected prime and reserves
the final prime as an independent held-out check.  A passing artifact is still
a QQ candidate for the q24/orbit85 construction, not the final equation-level
certificate: the next step must replay the exact RR functions and compile the
D12 child equation over QQ.
"""

import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, Zmod


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

ACCEPTED_SIGNATURE_STATUSES = {
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
}
PARTIAL_STATUS = "NEED_MORE_Q24_ORBIT85_CRT_PRIMES"
COMPLETE_STATUS = "PASS_Q24_ORBIT85_RR_KERNEL_CRT_HELDOUT_VALIDATED"


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
    primes = []
    for token in str(raw).replace(",", " ").split():
        primes.append(ZZ(token))
    return primes


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--primes",
    help="Comma/space separated prime list. Defaults to all q24/orbit85 signatures.",
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

primes = parse_primes(args.primes)
if len(primes) < 3:
    raise ValueError("need at least three signatures: two train primes plus one held-out")

sigs = []
for p in primes:
    path = LOCAL / f"q24-orbit85-d12-signature-mod-{p}.json"
    if not path.exists():
        raise SystemExit(f"missing signature artifact: {path}")
    d = json.loads(path.read_text())
    assert d["status"] in ACCEPTED_SIGNATURE_STATUSES
    assert ZZ(d["prime"]) == p
    sigs.append(d)

ref = sigs[0]
for d in sigs[1:]:
    assert d["source_neighbor"] == ref["source_neighbor"] == {
        "q": 24,
        "orbit": 85,
        "source": "D13/MW4",
        "child": "D12/MW5",
    }
    assert d["rr"] == ref["rr"] == {
        "ambient": 56,
        "collision_rank": 48,
        "post_collision": 8,
        "resolved_rank": 6,
        "kernel": 2,
        "geometric_fibre_twist": -8,
    }
    assert d["plane_pivots"] == ref["plane_pivots"] == [0, 1]
    assert d["quartic_degree"] == ref["quartic_degree"] == 4
    assert (d["child_root_rank"], d["child_root_det"], d["child_euler"]) == (
        12,
        4,
        24,
    )
    assert d["jacobian_A"]["num_degree"] == ref["jacobian_A"]["num_degree"]
    assert d["jacobian_A"]["den_degree"] == ref["jacobian_A"]["den_degree"]
    assert d["jacobian_B"]["num_degree"] == ref["jacobian_B"]["num_degree"]
    assert d["jacobian_B"]["den_degree"] == ref["jacobian_B"]["den_degree"]
    assert len(d["quartic_coefficients"]) == len(ref["quartic_coefficients"]) == 5
    for got, expected in zip(d["quartic_coefficients"], ref["quartic_coefficients"]):
        assert got["num_degree"] == expected["num_degree"]
        assert got["den_degree"] == expected["den_degree"]


def crt_scalar(residues, mods):
    x = ZZ(0)
    M = ZZ(1)
    for rr, pp in zip(residues, mods):
        pp = ZZ(pp)
        rr = ZZ(rr) % pp
        t = ((rr - x) % pp) * ((M % pp).inverse_mod(pp)) % pp
        x += M * t
        M *= pp
        x %= M
    return x, M


def rr_scalar(residues, mods):
    x, M = crt_scalar(residues, mods)
    try:
        q = QQ(Zmod(M)(x).rational_reconstruction())
    except (ValueError, ArithmeticError, ZeroDivisionError):
        return None, M
    return q, M


def reduce_q(q, p):
    q = QQ(q)
    p = ZZ(p)
    den = ZZ(q.denominator()) % p
    if den == 0:
        return None
    return int((ZZ(q.numerator()) % p) * den.inverse_mod(p) % p)


train_primes = primes[:-1]
hold_p = primes[-1]
train = sigs[:-1]
hold = sigs[-1]
_, crt_modulus = crt_scalar([0] * len(train_primes), train_primes)


def reconstruct_array(name, getter):
    shape0 = getter(train[0])
    vals = []
    total = 0
    recovered = 0
    heldout_valid = 0
    heldout_mismatch = 0
    heldout_undefined = 0
    max_num_bits = 0
    max_den_bits = 0

    for j in range(len(shape0)):
        residues = [getter(d)[j] for d in train]
        q, M = rr_scalar(residues, train_primes)
        total += 1
        if q is None:
            vals.append(None)
            continue
        recovered += 1
        max_num_bits = max(max_num_bits, abs(ZZ(q.numerator())).nbits())
        max_den_bits = max(max_den_bits, abs(ZZ(q.denominator())).nbits())
        hv = reduce_q(q, hold_p)
        if hv is None:
            heldout_undefined += 1
        elif hv == int(getter(hold)[j]) % int(hold_p):
            heldout_valid += 1
        else:
            heldout_mismatch += 1
        vals.append(q)

    status = "PASS_HELDOUT" if heldout_valid == total else "PARTIAL"
    print(
        "Q24O85CRT_PROGRESS|"
        f"object={name}|train={len(train_primes)}|holdout={hold_p}|"
        f"recovered={recovered}/{total}|heldout={heldout_valid}/{total}|"
        f"mismatch={heldout_mismatch}|undefined={heldout_undefined}|"
        f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
        f"status={status}",
        flush=True,
    )
    return {
        "values": vals,
        "total": total,
        "recovered": recovered,
        "heldout_valid": heldout_valid,
        "heldout_mismatch": heldout_mismatch,
        "heldout_undefined": heldout_undefined,
        "max_num_bits": max_num_bits,
        "max_den_bits": max_den_bits,
    }


def flat_plane(d):
    return [int(v) for row in d["plane_rref_2x56"] for v in row]


def flat_norm_rf(item):
    return [int(v) for v in item["num"]] + [int(v) for v in item["den"]]


def flat_jacobian_A(d):
    return flat_norm_rf(d["jacobian_A"])


def flat_jacobian_B(d):
    return flat_norm_rf(d["jacobian_B"])


def flat_quartic(d):
    values = []
    for item in d["quartic_coefficients"]:
        values.extend(flat_norm_rf(item))
    return values


rp = reconstruct_array("plane", flat_plane)
rq = reconstruct_array("quartic", flat_quartic)
ra = reconstruct_array("jacobian_A", flat_jacobian_A)
rb = reconstruct_array("jacobian_B", flat_jacobian_B)

objects = (rp, rq, ra, rb)
complete = all(r["heldout_valid"] == r["total"] for r in objects)


def progress_without_values(result):
    return {key: value for key, value in result.items() if key != "values"}


payload = {
    "schema": "elkies-k3.h3-q24-orbit85-rr-kernel-crt-qq.v1",
    "status": COMPLETE_STATUS if complete else PARTIAL_STATUS,
    "source_neighbor": ref["source_neighbor"],
    "rr": ref["rr"],
    "training_primes": [int(p) for p in train_primes],
    "heldout_prime": int(hold_p),
    "crt_modulus_bits": int(crt_modulus.nbits()),
    "plane_pivots": ref["plane_pivots"],
    "progress": {
        "plane": progress_without_values(rp),
        "quartic": progress_without_values(rq),
        "jacobian_A": progress_without_values(ra),
        "jacobian_B": progress_without_values(rb),
    },
    "proof_boundary": (
        "CRT/rational reconstruction of the canonical modular q24/orbit85 "
        "resolved-RR kernel and compact quartic/Jacobian data. A PASS status "
        "means the reconstructed rational coefficients reduce to the held-out "
        "signature prime. The exact QQ RR replay and D12 equation certificate "
        "remain separate."
    ),
}

if complete:
    pv = rp["values"]
    payload["qq_candidate"] = {
        "plane_rref_2x56": [
            [str(pv[56 * i + j]) for j in range(56)]
            for i in range(2)
        ],
        "quartic_coefficients": [],
        "jacobian_A": {
            "num_degree": ref["jacobian_A"]["num_degree"],
            "den_degree": ref["jacobian_A"]["den_degree"],
        },
        "jacobian_B": {
            "num_degree": ref["jacobian_B"]["num_degree"],
            "den_degree": ref["jacobian_B"]["den_degree"],
        },
    }

    offset = 0
    qvals = rq["values"]
    for item in ref["quartic_coefficients"]:
        n = len(item["num"])
        d = len(item["den"])
        payload["qq_candidate"]["quartic_coefficients"].append({
            "num_degree": item["num_degree"],
            "den_degree": item["den_degree"],
            "num": [str(q) for q in qvals[offset:offset + n]],
            "den": [str(q) for q in qvals[offset + n:offset + n + d]],
        })
        offset += n + d

    avals = ra["values"]
    An = len(ref["jacobian_A"]["num"])
    payload["qq_candidate"]["jacobian_A"]["num"] = [
        str(q) for q in avals[:An]
    ]
    payload["qq_candidate"]["jacobian_A"]["den"] = [
        str(q) for q in avals[An:]
    ]

    bvals = rb["values"]
    Bn = len(ref["jacobian_B"]["num"])
    payload["qq_candidate"]["jacobian_B"]["num"] = [
        str(q) for q in bvals[:Bn]
    ]
    payload["qq_candidate"]["jacobian_B"]["den"] = [
        str(q) for q in bvals[Bn:]
    ]

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit85-rr-kernel-crt-qq.json"
)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24O85CRT_RESULT|"
    f"primes={','.join(map(str, primes))}|train={len(train_primes)}|"
    f"holdout={hold_p}|bits={crt_modulus.nbits()}|status={payload['status']}",
    flush=True,
)
