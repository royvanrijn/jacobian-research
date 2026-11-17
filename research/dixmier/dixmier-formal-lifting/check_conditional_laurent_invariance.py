#!/usr/bin/env python3
"""Exact finite controls for the conditional Laurent-invariance lemma.

Independent Python integers/Fraction Laurent arithmetic. The all-N result
is a written theorem with explicit hypotheses, not a finite extrapolation.
No DC-pair or ordinary differential-operator realization is asserted.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile

HERE = Path(__file__).resolve().parent


def clean(p):
    return {e: c for e, c in p.items() if c}


def add(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            out[e + f] = out.get(e + f, Q(0)) + c * d
    return clean(out)


def power(p, n):
    out = {0: Q(1)}
    for _ in range(n):
        out = mul(out, p)
    return out


def binom(a, n):
    out = Q(1)
    for j in range(n):
        out *= (a - j) / (j + 1)
    return out


def frames(Y, D, q1):
    # p=4, q0=6, P_N=t^(2D), epsilon0=epsilon1=1.
    F1 = add(power(Y, 4), {12 * D: Q(-1)})
    F2 = add(power(F1, 4), {2 * D * q1: Q(-1)})
    return F1, F2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    assert binom(Q(1, 4), 5) == Q(231, 8192)
    nonresonant = []
    resonant = []
    for N in [3, 5, 9]:
        D = 2 * N + 1
        pN = 2 * D
        gamma = -pN + N + 1
        # Original orders 4,6,22,81; distances 2+7=9=p+q0-1.
        # The compressed root is z^3(1+z^-1)^(1/4). Omit its fifth
        # binomial coefficient to enforce the exact trace coefficient0.
        Y = {(3 - j) * D: binom(Q(1, 4), j) for j in range(5)}
        Y = add(Y, {gamma: Q(1)})
        F1, F2 = frames(Y, D, 22)
        qI = 40 * D + N + 1
        delta = 4 * 11 * D - qI
        C = (4 * 22 - 3 * (6 + 22) + 4) // 2
        assert C == 4 and delta == D * C - N - 1
        assert delta == 7 * N + 3 and delta > N + 1
        assert max(F1) == 11 * D and F1[11 * D] == 1
        assert max(F2) == qI and F2[qI] == 16
        assert Y.get(-pN, Q(0)) == 0
        assert F2.get(qI - 1, Q(0)) == 0
        assert F1[7 * D] == Q(-231, 2048)
        assert F2[40 * D] == Q(-231, 512)
        # c_N[2D]=c_original[4]=1/4 is the compressed convention.
        # The literal uncompressed c_N[4D]=c_original[4] fails.
        assert Y[2 * D] == Q(1, 4) and Y.get(4 * D, Q(0)) == 0
        nonresonant.append({
            "N": N, "D": D, "p_N": pN, "C": C, "last_gap": delta,
            "final_head": qI, "final_leading_coefficient": "16",
            "trace_coefficient": "0", "head_minus_one": "0",
            "prefix_distinguished": str(F1[7 * D]),
            "final_distinguished": str(F2[40 * D]),
            "compressed_index_value": "1/4",
            "literal_raw_index_countercontrol": "0 != 1/4",
        })
        # Original orders 4,6,16,63; distances 8+1=9.
        # The last N-gap makes the missing head-minus-one condition
        # necessary. Retain v!=0 as deliberately invalid controls.
        for v in [Q(0), Q(1), Q(-2)]:
            Y = add({3 * D: Q(1), -D: Q(1, 4)},
                    {gamma: Q(1)}, {gamma - 1: v})
            F1, F2 = frames(Y, D, 16)
            qI = 31 * D + N + 1
            delta = 4 * 8 * D - qI
            C = (4 * 16 - 3 * (6 + 16) + 4) // 2
            assert C == 1 and delta == D * C - N - 1
            assert delta == N
            assert max(F1) == 8 * D and F1[8 * D] == 1
            assert max(F2) == qI and F2[qI] == 16
            assert Y.get(-pN, Q(0)) == 0
            assert F1.get(7 * D, Q(0)) == 0
            assert F2.get(qI - 1, Q(0)) == 16 * v
            assert F2.get(31 * D, Q(0)) == 192 * v
            resonant.append({
                "N": N, "D": D, "v": str(v), "C": C, "last_gap": delta,
                "final_head": qI, "final_leading_coefficient": "16",
                "trace_coefficient": "0",
                "head_minus_one": str(16 * v),
                "final_distinguished": str(F2.get(31 * D, Q(0))),
                "satisfies_required_missing_head": not bool(v),
            })
    # A direct finite cyclic-vector trace control for arbitrary q mod s.
    # The universal statement is the written permutation-sum proof.
    cyclic_controls = []
    for s in [2, 3, 5]:
        vector = [Q(j * j - 2 * j + 3, j + 1) for j in range(s)]
        for q in range(2 * s):
            difference = [vector[(j + q) % s] - vector[j] for j in range(s)]
            assert sum(difference, Q(0)) == 0
        cyclic_controls.append({"components": s, "tested_q": list(range(2 * s))})

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    receipt = {
        "status": "PASS",
        "claim": "Exact finite controls for an independently written conditional Laurent-coefficient invariance lemma",
        "environment": {"python": platform.python_version(), "arithmetic": "Python integers and Fraction; no algebra or nullspace solver"},
        "primary_version": "arXiv:2410.06959v5, 19 January 2026",
        "nonresonant_data": {"p": 4, "q0": 6, "q1": 22, "d0": 2, "l0": 1, "original_final_order": 81},
        "nonresonant_controls": nonresonant,
        "resonant_data": {"p": 4, "q0": 6, "q1": 16, "d0": 2, "l0": 1, "original_final_order": 63},
        "resonant_controls": resonant,
        "written_gap_dichotomy": "delta_N=D*C-N-1 with integral C>=1; only C=1,d0=2,l0=1 gives delta_N=N, and every other case gives delta_N>N+1. No separate gap assumption or N+1 boundary remains.",
        "cyclic_trace_controls": cyclic_controls,
        "hashes": {
            "CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md": sha(HERE / "CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md"),
            "check_conditional_laurent_invariance.py": sha(Path(__file__)),
            "sources/zheglov-2410.06959v5.pdf": sha(HERE / "sources/zheglov-2410.06959v5.pdf"),
        },
        "coverage_boundary": "Finite commutative Laurent controls only. The all-N lemma is proved in the note under explicit head, final-order, trace and exceptional normalization hypotheses; its gap dichotomy follows from integrality. No DC-pair, differential-operator realization, earlier-source theorem, integer-polynomial coefficient assertion or main DC1 theorem is certified.",
        "formal_verification": False, "external_review": False,
    }
    out = args.output_dir or Path(tempfile.mkdtemp(prefix="conditional-laurent-invariance-"))
    out.mkdir(parents=True, exist_ok=True)
    destination = out / "conditional_laurent_invariance_controls.json"
    destination.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: 3 nonzero invariance controls, 9 resonant/invalid controls, cyclic trace controls")
    print("Receipt:", destination)


if __name__ == "__main__":
    main()
