#!/usr/bin/env python3
"""Independent finite audit of the uniform cubic resultant valuation.

The unrestricted statement is proved in the accompanying research note.
This script uses FLINT integer resultants, not the root-by-root proof.
"""
import argparse
from hashlib import sha256
import json
import platform
import math
from fractions import Fraction
from pathlib import Path

import flint
from flint import fmpz_poly, nmod_poly


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def v3(value):
    value = abs(int(value))
    require(value != 0, "zero resultant")
    answer = 0
    while value % 3 == 0:
        value //= 3
        answer += 1
    return answer


def model(n):
    out = [0] * (n + 1)
    for k in range(n + 1):
        out[n - k] = (-1) ** k * math.comb(n + 2 * k, 3 * k)
    return fmpz_poly(out)


EXPECTED_NOTE_SHA256 = "0dff58f64cd64b3e324e71c17d7314de849eaeccf797e4177414a89c7684b57a"
NOTE_PATH = "research/factorial/image-rigidity/SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md"


def run(repo_root):
    note_hash = sha256((repo_root / NOTE_PATH).read_bytes()).hexdigest()
    require(note_hash == EXPECTED_NOTE_SHA256, "reviewed proof-source hash mismatch")
    records = []
    for n in list(range(1, 33)) + [39, 52]:
        p, q = model(n), model(n + 1)
        residual = nmod_poly(list(p), 3).gcd(nmod_poly(list(q), 3))
        expected = residual.degree()
        actual = v3(p.resultant(q))
        require(actual == expected, (n, actual, expected))
        records.append({"n": n, "residual_gcd_degree": expected,
                        "resultant_valuation": actual})

    deformations = []
    for n in range(1, 17):
        p, q = model(n), model(n + 1)
        expected = nmod_poly(list(p), 3).gcd(nmod_poly(list(q), 3)).degree()
        for variant in range(3):
            dp = fmpz_poly([9 * (((k + 3) * (variant + 2) + n) % 11 - 5)
                           for k in range(n)])
            dq = fmpz_poly([9 * (((k + 1) * (variant + 3) + n) % 13 - 6)
                           for k in range(n + 1)])
            actual = v3((p + dp).resultant(q + dq))
            require(actual == expected, (n, variant, actual, expected))
            deformations.append({"n": n, "variant": variant, "valuation": actual})

    raw = []
    for n in range(1, 13):
        def weighted(m):
            return fmpz_poly([math.comb(m, k) * math.factorial(m + 2 * k)
                              * (1 + 9 * ((m + 2 * k) ** 2 + 1))
                              for k in range(m + 1)])
        p, q = model(n), model(n + 1)
        deg = nmod_poly(list(p), 3).gcd(nmod_poly(list(q), 3)).degree()
        expected = (n * (n + 1) + (n + 1) * v3(math.factorial(n))
                    + n * v3(math.factorial(n + 1)) + deg)
        actual = v3(weighted(n).resultant(weighted(n + 1)))
        require(actual == expected, (n, actual, expected))
        raw.append({"n": n, "valuation": actual})

    # The weaker modulus-three hypothesis really fails.
    units = {3: 220, 4: -12860}
    for n in (3, 4):
        value = sum(Fraction(math.comb(n, k) * math.factorial(n + 2 * k)
                             * units.get(n + 2 * k, 1), (-6) ** k)
                    for k in range(n + 1))
        require(value == 0, ("mod3 countercontrol", n, value))
    require(all(u % 3 == 1 for u in units.values()), "incorrect countercontrol units")

    result = {"status": "PASS", "arithmetic": "FLINT exact integer resultants",
              "model_pairs": records, "deformed_pairs": deformations,
              "raw_weighted_pairs": raw,
              "mod3_common_root_countercontrol": {"orders": [3, 4], "phase": "-1/6"},
              "scope": "Finite audit; unrestricted quantifiers rely on the written proof."}
    own_path = "research/factorial/image-rigidity/verify_sfc_cubic_resultant_independent.py"
    result.update({
        "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
        "python_version": platform.python_version(),
        "python_flint_version": flint.__version__,
        "flint_version": getattr(flint, "__FLINT_VERSION__", "not exposed"),
        "source_pins": [
            {"path": NOTE_PATH, "sha256": note_hash,
             "expected_sha256": EXPECTED_NOTE_SHA256, "reviewed_pin_verified": True},
            {"path": own_path, "sha256": sha256((repo_root / own_path).read_bytes()).hexdigest(),
             "role": "current checker bytes; external manifests may pin this hash"}
        ],
        "unrestricted_quantifiers_proved_by_finite_sampling": False,
        "universal_proof": "SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md, Sections 1-3"
    })
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.repo_root.resolve())
    target = args.output or Path(__file__).with_name("sfc_cubic_resultant_independent_verification.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "model_pairs": len(result["model_pairs"]),
                      "deformed_pairs": len(result["deformed_pairs"]),
                      "raw_pairs": len(result["raw_weighted_pairs"]), "output": str(target)}))
