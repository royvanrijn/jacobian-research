#!/usr/bin/env python3
"""Build a closed conductor-first pilot around the three family anchors.

The population is deliberately small and exact: offsets -4 through 4 around
the Fermigier numerator 28917 (denominator 20), Mestre family 2 at u=483,
and Mestre family 3 at u=660.  Every fibre receives a cheap partial
discriminant factorization.  In each lane, the anchor and the two best other
sieve records receive global minimalization and a full local Tate replay.  A
full-dimensional known Kummer subgroup is then certified on those survivors.
No standalone point search or Selmer computation is run here.  The `u=481`
record pins two exact rational points returned incidentally by a separate
provisional descent diagnostic, while withholding that diagnostic's upper
bound.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence

from sympy import factorint, isprime


ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))

from analyze_icarm_7fff_zip_sequence import kodaira_symbol, valuation  # noqa: E402
from build_conductor_first_near_miss_targets import (  # noqa: E402
    discriminant,
    mod2_certificate,
    no_rational_two_torsion_witness,
    pari_small_prime_saturation,
    point_record,
)
from conductor_first_pipeline import work_queues  # noqa: E402
from fermigier_mestre import FermigierMestreFamily  # noqa: E402
from search_mestre_dsquare_four import (  # noqa: E402
    FAMILIES,
    base_parameter,
    known_jacobian_points,
)
from search_mestre_root_tuple_scale_max200 import (  # noqa: E402
    mod3_independence_certificate,
)


Q = Fraction
OFFSETS = tuple(range(-4, 5))
FACTOR_LIMIT = 10_000
TATE_KEEP_PER_LANE = 3
TATE_TIMEOUT_SECONDS = 45.0
GP_STACK_BYTES = 1_000_000_000
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_family_anchor_pilot_v1.json"
)
FIXED_TARGET_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_near_miss_descent_targets_v1.json"
)
FIXED_TARGET_MANIFEST_SHA256 = (
    "cfd2ac0d0ec995c102df534fc216fbb849a21970ed68ce6785ce000f12cba0ff"
)
ANCHOR_TARGET_IDS = {
    "fermigier-u28917-20": "fermigier-u28917-20",
    "family2-u483": "family2-u483",
    "family3-u660": "family3-u660",
}
FAMILY2_U481_DESCENT_POINTS = (
    (
        Q(1_751_764_273_991_163, 114_244),
        Q(316_701_696_272_026_942_800, 371_293),
    ),
    (
        Q(4_656_599_310_055_995, 114_244),
        Q(160_497_779_458_185_494_424, 28_561),
    ),
)


def fraction_text(value: Fraction | int) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def stable_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def model_record(model: Sequence[Fraction]) -> list[str]:
    return [fraction_text(value) for value in model]


def population() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for offset in OFFSETS:
        numerator = 28_917 + offset
        parameter_u = Q(numerator, 20)
        parameter_t = 2 * parameter_u
        rows.append(
            {
                "id": f"fermigier-u{numerator}-20",
                "lane": "fermigier",
                "offset": offset,
                "anchor": offset == 0,
                "parameter_u": fraction_text(parameter_u),
                "family_parameter": fraction_text(parameter_t),
                "source_model": model_record(
                    FermigierMestreFamily.coefficients(parameter_t)
                ),
            }
        )
    for family_index, anchor in ((2, 483), (3, 660)):
        family = FAMILIES[family_index]
        for offset in OFFSETS:
            parameter_u = Q(anchor + offset)
            parameter_t = base_parameter(family, parameter_u)
            rows.append(
                {
                    "id": f"family{family_index}-u{anchor + offset}",
                    "lane": f"family{family_index}",
                    "offset": offset,
                    "anchor": offset == 0,
                    "parameter_u": fraction_text(parameter_u),
                    "family_parameter": fraction_text(parameter_t),
                    "roots": list(family.roots),
                    "source_model": model_record(
                        family.construction.primitive_jacobian_coefficients(
                            parameter_t
                        )
                    ),
                }
            )
    return rows


def cheap_sieve(model: Sequence[str]) -> dict[str, Any]:
    delta = Q(discriminant(tuple(Q(value) for value in model)))
    numerator = abs(delta.numerator)
    if numerator == 0:
        return {"status": "reject", "reason": "singular discriminant"}
    factors = {int(base): int(exponent) for base, exponent in factorint(
        numerator, limit=FACTOR_LIMIT
    ).items()}
    complete = all(isprime(base) for base in factors)
    radical_proxy = prod(factors)
    squareful_proxy = numerator // radical_proxy
    composite_cofactors = sorted(base for base in factors if not isprime(base))
    return {
        "status": "pass",
        "discriminant_numerator": str(numerator),
        "discriminant_denominator": str(delta.denominator),
        "factor_limit": FACTOR_LIMIT,
        "factorization": [[str(base), exponent] for base, exponent in sorted(factors.items())],
        "factorization_complete": complete,
        "composite_cofactors": [str(value) for value in composite_cofactors],
        "radical": str(radical_proxy) if complete else None,
        "radical_proxy": str(radical_proxy),
        "squareful_quotient_proxy": str(squareful_proxy),
        "priority_key": [
            radical_proxy.bit_length(),
            radical_proxy,
            -squareful_proxy.bit_length(),
        ],
        "interpretation": (
            "cheap discriminant-numerator ordering only; incomplete rows use a "
            "composite-cofactor proxy and are not conductor claims"
        ),
    }


def gp_rational(value: str) -> str:
    rational = Q(value)
    return f"({rational.numerator}/{rational.denominator})"


def exact_tate(model: Sequence[str]) -> dict[str, Any]:
    gp = shutil.which("gp")
    if gp is None:
        raise FileNotFoundError("PARI/GP executable 'gp' is required")
    vector = ",".join(gp_rational(value) for value in model)
    program = f'''
default(realprecision,80);
E=ellinit([{vector}]);v=0;M=ellminimalmodel(E,&v);G=ellglobalred(M);
print("PARI_VERSION\\t",version());
print("MIN_MODEL\\t",M.a1,"\\t",M.a2,"\\t",M.a3,"\\t",M.a4,"\\t",M.a6);
print("MIN_CHANGE\\t",v[1],"\\t",v[2],"\\t",v[3],"\\t",v[4]);
print("MIN_DISC\\t",M.disc);
print("CONDUCTOR\\t",G[1]);
print("LOG_CONDUCTOR\\t",log(G[1]));
print("TAMAGAWA_PRODUCT\\t",G[3]);
print("ROOT_NUMBER\\t",ellrootno(M));
for(i=1,matsize(G[4])[1],p=G[4][i,1];L=elllocalred(M,p);print("LOCAL\\t",p,"\\t",L[1],"\\t",L[2],"\\t",L[3][1],"\\t",L[3][2],"\\t",L[3][3],"\\t",L[3][4],"\\t",L[4]));
quit
'''
    completed = subprocess.run(
        [gp, "-q", "-s", str(GP_STACK_BYTES), "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=TATE_TIMEOUT_SECONDS,
    )
    combined = completed.stdout + completed.stderr
    if "***" in combined:
        raise RuntimeError(" ".join(combined.split())[:1000])
    scalar: dict[str, list[str]] = {}
    local_rows: list[list[str]] = []
    for line in completed.stdout.splitlines():
        fields = line.strip().split("\t")
        if not fields or len(fields) == 1:
            continue
        if fields[0] == "LOCAL":
            local_rows.append(fields[1:])
        else:
            scalar[fields[0]] = fields[1:]
    minimal_discriminant = int(scalar["MIN_DISC"][0])
    conductor = int(scalar["CONDUCTOR"][0])
    local_reductions = []
    for row in local_rows:
        prime, exponent, code, u, r, s, t, tamagawa = map(int, row)
        local_reductions.append(
            {
                "prime": str(prime),
                "minimal_discriminant_valuation": valuation(
                    abs(minimal_discriminant), prime
                ),
                "conductor_exponent": exponent,
                "kodaira_code": code,
                "kodaira_symbol": kodaira_symbol(code),
                "local_minimal_change": [u, r, s, t],
                "tamagawa_number": tamagawa,
            }
        )
    reconstructed = prod(
        int(row["prime"]) ** int(row["conductor_exponent"])
        for row in local_reductions
    )
    if reconstructed != conductor:
        raise ArithmeticError("local conductor exponents did not reconstruct N")
    tamagawa_product = prod(
        int(row["tamagawa_number"]) for row in local_reductions
    )
    if tamagawa_product != int(scalar["TAMAGAWA_PRODUCT"][0]):
        raise ArithmeticError("local Tamagawa numbers did not reconstruct the product")
    return {
        "status": "complete",
        "engine": "PARI/GP " + scalar["PARI_VERSION"][0].strip("[]").replace(", ", "."),
        "global_minimal": True,
        "global_minimal_model": scalar["MIN_MODEL"],
        "source_to_minimal_change": scalar["MIN_CHANGE"],
        "minimal_discriminant": str(minimal_discriminant),
        "conductor": str(conductor),
        "log_conductor": scalar["LOG_CONDUCTOR"][0],
        "root_number": int(scalar["ROOT_NUMBER"][0]),
        "tamagawa_product": tamagawa_product,
        "local_reductions": local_reductions,
        "conductor_reconstructed_from_local_exponents": str(reconstructed),
    }


def selected_ids(records: Sequence[dict[str, Any]]) -> set[str]:
    selected: set[str] = set()
    for lane in ("fermigier", "family2", "family3"):
        rows = [record for record in records if record["lane"] == lane]
        anchor = next(record for record in rows if record["anchor"])
        selected.add(anchor["id"])
        others = sorted(
            (record for record in rows if not record["anchor"]),
            key=lambda record: (
                record["cheap_sieve"]["priority_key"], record["id"]
            ),
        )
        selected.update(record["id"] for record in others[: TATE_KEEP_PER_LANE - 1])
    return selected


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixed_anchor_subgroups() -> dict[str, dict[str, Any]]:
    actual = sha256_file(FIXED_TARGET_MANIFEST)
    if actual != FIXED_TARGET_MANIFEST_SHA256:
        raise AssertionError(
            f"changed fixed-target manifest: {actual} != {FIXED_TARGET_MANIFEST_SHA256}"
        )
    manifest = json.loads(FIXED_TARGET_MANIFEST.read_text())
    targets = {target["id"]: target for target in manifest["targets"]}
    answer = {}
    for pilot_id, target_id in ANCHOR_TARGET_IDS.items():
        target = targets[target_id]
        answer[pilot_id] = {
            "status": "complete",
            "method": "imported source-pinned fixed-target subgroup certificate",
            "source_manifest": str(FIXED_TARGET_MANIFEST.relative_to(ROOT)),
            "source_manifest_sha256": actual,
            "source_target_id": target_id,
            "certified_rank_lower_bound": target["certified_known_rank"],
            "kummer_image_dimension": target["known_kummer_image_dimension"],
            "rational_two_torsion_dimension": target[
                "rational_two_torsion_dimension"
            ],
            "exact_point_membership_checked": True,
            "known_basis": target["known_basis"],
            "known_basis_sha256": target["known_basis_sha256"],
            "mod2_certificate": target["known_basis_mod2_certificate"],
        }
    return answer


def new_subgroup(record: dict[str, Any]) -> dict[str, Any]:
    model = tuple(Q(value) for value in record["source_model"])
    if record["lane"] == "fermigier":
        parameter_t = Q(record["family_parameter"])
        points = FermigierMestreFamily.known_jacobian_points(parameter_t)
        point_source = "twelve visible plus one extra Fermigier--Mestre section"
    else:
        family = FAMILIES[int(record["lane"].removeprefix("family"))]
        parameter_u = Q(record["parameter_u"])
        points = known_jacobian_points(family, parameter_u)
        point_source = "twelve visible plus split-infinity Mestre section"
    mod3 = mod3_independence_certificate(model, points, prime_bound=499)
    baseline_rank = int(mod3["combined_exact_rank_over_F3"])
    indices = [int(index) for index in mod3["independent_subset_indices_one_based"]]
    independent = tuple(points[index - 1] for index in indices)
    saturated = pari_small_prime_saturation(model, independent)
    basis_points = saturated
    descent_discovery = None
    if record["id"] == "family2-u481":
        basis_points += FAMILY2_U481_DESCENT_POINTS
        descent_discovery = {
            "engine": "temporary PARI 2.19.0-development commit 6af5b91",
            "command": "ellrank(E,0,known_basis) with an 8-GB PARI stack",
            "provisional_output": {
                "rank_interval": [14, 14],
                "sha_indicator": 0,
                "elapsed_milliseconds": 330449,
            },
            "claim_boundary": (
                "the upper bound is withheld because the BNF was not certified; "
                "only the displayed rational points are used below"
            ),
            "new_points": [
                point_record(point) for point in FAMILY2_U481_DESCENT_POINTS
            ],
        }
    rank = len(basis_points)
    mod2 = mod2_certificate(model, basis_points)
    if int(mod2["combined_binary_rank"]) != rank:
        raise ArithmeticError("bounded saturation did not yield a full mod-2 image")
    torsion = no_rational_two_torsion_witness(model)
    basis = [point_record(point) for point in basis_points]
    return {
        "status": "complete",
        "method": (
            "exact mod-3 independence, bounded PARI saturation through 3 used "
            "only for discovery, then independent exact mod-2 finite quotients"
        ),
        "point_source": point_source,
        "visible_point_count": len(points),
        "visible_baseline_rank_lower_bound": baseline_rank,
        "mod3_independent_subset_indices_one_based": indices,
        "mod3_certificate": mod3,
        "saturation_prime_bound": 3,
        "global_saturation_claim": False,
        "descent_point_discovery": descent_discovery,
        "certified_rank_lower_bound": rank,
        "kummer_image_dimension": rank,
        "rational_two_torsion_dimension": 0,
        "rational_two_torsion_exclusion": torsion,
        "exact_point_membership_checked": True,
        "known_basis": basis,
        "known_basis_sha256": stable_sha256(basis),
        "mod2_certificate": mod2,
    }


def build() -> dict[str, Any]:
    records = population()
    for record in records:
        record["cheap_sieve"] = cheap_sieve(record["source_model"])
    keep = selected_ids(records)
    selected = [record for record in records if record["id"] in keep]
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            record["id"]: executor.submit(exact_tate, record["source_model"])
            for record in selected
        }
        for record in selected:
            try:
                record["tate"] = futures[record["id"]].result()
            except subprocess.TimeoutExpired:
                record["tate"] = {
                    "status": "timeout",
                    "timeout_seconds": TATE_TIMEOUT_SECONDS,
                }
            except Exception as error:
                record["tate"] = {"status": "error", "error": str(error)[:1000]}
    anchors = fixed_anchor_subgroups()
    for record in selected:
        if record.get("tate", {}).get("status") != "complete":
            continue
        record["known_subgroup"] = (
            anchors[record["id"]]
            if record["id"] in anchors
            else new_subgroup(record)
        )
    records.sort(key=lambda record: (record["lane"], record["offset"]))
    payload: dict[str, Any] = {
        "schema": "elliptic-curves.conductor-first-family-anchor-pilot.v1",
        "reproducing_command": (
            ".venv/bin/python elliptic-curves/cas/"
            "build_conductor_first_family_pilot.py --check"
        ),
        "population": {
            "offsets": list(OFFSETS),
            "lanes": {
                "fermigier": "u=(28917+offset)/20",
                "family2": "u=483+offset",
                "family3": "u=660+offset",
            },
            "candidate_count": len(records),
            "tate_selection_per_lane": TATE_KEEP_PER_LANE,
            "selection_rule": (
                "anchor plus two non-anchor rows with lexicographically smallest "
                "cheap discriminant-radical proxy priority"
            ),
        },
        "candidates": records,
        "queues": work_queues(records),
        "claim_boundary": (
            "cheap factorization is heuristic ordering; exact Tate records "
            "and known-subgroup records are certificates; the u=481 point lower "
            "bound is exact, but no complete Selmer or cover classification is claimed"
        ),
    }
    payload["result_sha256"] = stable_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != text:
            raise AssertionError(f"stale or missing pilot artifact: {args.output}")
        print(f"PASS {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
        print(f"WROTE {args.output}")
    for stage, ids in payload["queues"].items():
        print(f"{stage}: {len(ids)}")


if __name__ == "__main__":
    main()
