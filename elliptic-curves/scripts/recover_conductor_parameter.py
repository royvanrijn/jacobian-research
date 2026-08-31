#!/usr/bin/env python3
"""Recover bounded family parameters from local discriminant fingerprints."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.conductor_engineering import (  # noqa: E402
    discover_projective_root_balls,
    factor_over_known_primes,
    recover_parameters,
    select_projective_ball,
    select_repeated_prime_constraints,
    small_prime_valuations,
    weierstrass_invariant_polynomials,
    weierstrass_invariants,
)
from ecsearch.crt_lattice import p_adic_valuation  # noqa: E402


def _resolve(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else REPOSITORY_ROOT / path


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(path)


def _nested(mapping: dict[str, Any], keys: Sequence[str]) -> Any:
    value: Any = mapping
    for key in keys:
        value = value[key]
    return value


def load_family(specification: dict[str, Any]) -> dict[str, Any]:
    family_spec = specification["family"]
    source_path = _resolve(family_spec["file"])
    source = json.loads(source_path.read_text())
    model = _nested(
        source,
        family_spec.get("model_path", ["canonical_weierstrass_model"]),
    )
    polynomial_model = model["coefficients_low_to_high"]
    invariants = weierstrass_invariant_polynomials(polynomial_model)
    declared_discriminant = tuple(
        int(value) for value in model["discriminant_coefficients_low_to_high"]
    )
    if invariants["discriminant"] != declared_discriminant:
        raise AssertionError("the family discriminant does not match its a-invariants")
    return {
        "id": source.get("id", family_spec.get("id", source_path.stem)),
        "source": _display_path(source_path),
        "model": polynomial_model,
        "invariants": invariants,
    }


def load_target(specification: dict[str, Any]) -> dict[str, Any]:
    target_spec = specification["target"]
    source_path = _resolve(target_spec["record_file"])
    source = json.loads(source_path.read_text())
    records = _nested(source, target_spec.get("records_path", ["curves"]))
    match = target_spec["match"]
    matching = [record for record in records if record[match["field"]] == match["value"]]
    if len(matching) != 1:
        raise ValueError(f"target selector found {len(matching)} records")
    record = matching[0]
    coefficients = tuple(int(value) for value in record[target_spec.get("ainvs_field", "ainvs")])
    invariants = weierstrass_invariants(coefficients)
    reported_field = target_spec.get("discriminant_field", "discriminant")
    if reported_field in record and invariants["discriminant"] != int(record[reported_field]):
        raise AssertionError("the target's reported discriminant is inconsistent")
    return {
        "label": target_spec.get("label", str(match["value"])),
        "source": _display_path(source_path),
        "selector": match,
        "ainvs": coefficients,
        "submitted_invariants": invariants,
        "bad_primes": [int(value) for value in record.get("bad_primes", [])],
    }


def pari_minimal_and_local_data(
    coefficients: Sequence[int], primes: Sequence[int]
) -> dict[str, Any]:
    """Minimalize a target and replay Tate data at the selected primes."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP executable 'gp' is required")
    model = ",".join(map(str, coefficients))
    local_lines = "\n".join(
        f'L=elllocalred(M,{prime});print("LOCAL\\t{prime}\\t",L[1],"\\t",L[2],"\\t",L[4],"\\t",ellrootno(M,{prime}));'
        for prime in primes
    )
    program = f"""
default(parisizemax,1000000000);
E=ellinit([{model}]);v=0;M=ellminimalmodel(E,&v);
print("MODEL\\t",M.a1,"\\t",M.a2,"\\t",M.a3,"\\t",M.a4,"\\t",M.a6);
print("CHANGE\\t",v[1],"\\t",v[2],"\\t",v[3],"\\t",v[4]);
print("DISCRIMINANT\\t",M.disc);
{local_lines}
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    combined = completed.stdout + completed.stderr
    errors = [
        line
        for line in combined.splitlines()
        if "***" in line and "Warning:" not in line
    ]
    if errors:
        raise RuntimeError(combined)
    minimal_model: list[int] | None = None
    change: list[int] | None = None
    discriminant: int | None = None
    local: list[dict[str, Any]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if fields[0] == "MODEL":
            minimal_model = [int(value) for value in fields[1:]]
        elif fields[0] == "CHANGE":
            change = [int(value) for value in fields[1:]]
        elif fields[0] == "DISCRIMINANT":
            discriminant = int(fields[1])
        elif fields[0] == "LOCAL":
            prime, conductor_exponent, kodaira_code, tamagawa, local_root_number = map(
                int, fields[1:]
            )
            multiplicative_n = kodaira_code - 4 if kodaira_code >= 5 else None
            local.append(
                {
                    "prime": prime,
                    "conductor_exponent": conductor_exponent,
                    "kodaira_code": kodaira_code,
                    "kodaira_symbol": (
                        f"I{multiplicative_n}" if multiplicative_n is not None else None
                    ),
                    "tamagawa_number": tamagawa,
                    "local_root_number": local_root_number,
                    "split_multiplicative": (
                        multiplicative_n is not None and local_root_number == -1
                    ),
                }
            )
    if minimal_model is None or change is None or discriminant is None:
        raise RuntimeError(f"incomplete PARI output: {combined}")
    version = subprocess.run(
        [gp, "--version-short"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    return {
        "engine": f"PARI/GP {version}",
        "minimal_model": minimal_model,
        "change_u_r_s_t": change,
        "submitted_model_is_global_minimal": change == [1, 0, 0, 0],
        "minimal_discriminant": str(discriminant),
        "selected_local_reductions": local,
    }


def build(specification: dict[str, Any], raw_specification: bytes) -> dict[str, Any]:
    family = load_family(specification)
    target = load_target(specification)
    declared_constraints = specification.get("constraints")
    if declared_constraints:
        constraints = declared_constraints
        constraint_selection = {
            "method": "explicit",
            "reason": (
                "the specification pins local branches or family-specific "
                "prime exclusions"
            ),
        }
        primes = [int(constraint["prime"]) for constraint in constraints]
        pari = pari_minimal_and_local_data(target["ainvs"], primes)
    else:
        fingerprint = specification.get("automatic_constraints", {})
        maximum_prime = int(fingerprint.get("maximum_prime", 97))
        minimum_valuation = int(fingerprint.get("minimum_valuation", 2))
        excluded_primes = {
            int(prime) for prime in fingerprint.get("excluded_primes", [])
        }
        preliminary = pari_minimal_and_local_data(target["ainvs"], ())
        preliminary_discriminant = int(preliminary["minimal_discriminant"])
        selected_constraints = select_repeated_prime_constraints(
            preliminary_discriminant,
            maximum_prime=maximum_prime,
            minimum_valuation=minimum_valuation,
            excluded_primes=excluded_primes,
        )
        constraints = [
            {"prime": prime, "exponent": exponent}
            for prime, exponent in selected_constraints
        ]
        if not constraints:
            raise ValueError("automatic selection found no repeated small-prime powers")
        primes = [int(constraint["prime"]) for constraint in constraints]
        pari = pari_minimal_and_local_data(target["ainvs"], primes)
        constraint_selection = {
            "method": "automatic repeated-small-prime threshold",
            "maximum_prime": maximum_prime,
            "minimum_valuation": minimum_valuation,
            "excluded_primes": sorted(excluded_primes),
            "selected_prime_exponents": [
                [constraint["prime"], constraint["exponent"]]
                for constraint in constraints
            ],
        }
    if len(set(primes)) != len(primes):
        raise ValueError("constraints must use distinct primes")
    minimal_discriminant = int(pari["minimal_discriminant"])
    target_invariants = weierstrass_invariants(pari["minimal_model"])
    if target_invariants["discriminant"] != minimal_discriminant:
        raise AssertionError("minimal target invariants are inconsistent")

    discriminant_coefficients = family["invariants"]["discriminant"]
    local_profiles = []
    constraint_groups = []
    for constraint in constraints:
        prime = int(constraint["prime"])
        exponent = int(constraint.get("exponent", 0))
        target_valuation = p_adic_valuation(minimal_discriminant, prime)
        if exponent == 0:
            exponent = target_valuation
        if exponent > target_valuation:
            raise ValueError(
                f"requested v_{prime}>={exponent}, but target has {target_valuation}"
            )
        balls, profile = discover_projective_root_balls(
            discriminant_coefficients,
            prime,
            exponent,
            maximum_roots=int(specification["search"].get("maximum_roots", 100_000)),
        )
        if not balls:
            raise ValueError(
                f"the family discriminant has no projective root modulo {prime}^{exponent}"
            )
        branch = constraint.get("branch")
        selected = (select_projective_ball(balls, branch),) if branch else balls
        profile["target_minimal_discriminant_valuation"] = target_valuation
        profile["selected_branches"] = [ball.label for ball in selected]
        profile["selection"] = "declared branch" if branch else "all maximal balls"
        local_profiles.append(profile)
        constraint_groups.append(selected)

    search = specification["search"]
    recovery = recover_parameters(
        discriminant_coefficients,
        constraint_groups,
        coefficient_radius=int(search["coefficient_radius"]),
        height_cap=(int(search["height_cap"]) if search.get("height_cap") else None),
        weights=tuple(int(value) for value in search.get("weights", [1, 1])),
        family_c4=family["invariants"]["c4"],
        target_c4=target_invariants["c4"],
        target_discriminant=minimal_discriminant,
        candidate_limit=int(search.get("candidate_limit", 20)),
    )
    maximum_fingerprint_prime = int(
        specification.get("maximum_fingerprint_prime", max(primes))
    )
    declared_factorization, unfactored_cofactor = factor_over_known_primes(
        minimal_discriminant, target["bad_primes"]
    )
    return {
        "schema": "elliptic-curves.conductor-parameter-recovery.v1",
        "claim_level": "exact bounded computation",
        "specification_sha256": hashlib.sha256(raw_specification).hexdigest(),
        "family": {
            "id": family["id"],
            "source": family["source"],
            "discriminant_degree": len(discriminant_coefficients) - 1,
            "a_invariants_and_discriminant_identity_verified": True,
        },
        "target": {
            "label": target["label"],
            "source": target["source"],
            "selector": target["selector"],
            "submitted_ainvs": list(target["ainvs"]),
            "pari_minimalization": pari,
            "minimal_discriminant_small_prime_valuations": [
                list(pair)
                for pair in small_prime_valuations(
                    minimal_discriminant, maximum_fingerprint_prime
                )
            ],
            "factorization_over_record_bad_primes": [
                [str(prime), exponent]
                for prime, exponent in declared_factorization
            ],
            "unfactored_cofactor_after_record_bad_primes": str(
                unfactored_cofactor
            ),
        },
        "constraint_selection": constraint_selection,
        "local_fingerprints": local_profiles,
        "recovery": recovery,
        "interpretation": {
            "proved": (
                "Every returned candidate lies in the declared bounded CRT/Gauss "
                "search and has all displayed forced family-discriminant valuations. "
                "Every exact match has the target j-invariant."
            ),
            "boundary": (
                "Local valuations can select multiple p-adic balls. Declared branch "
                "choices are replay inputs, not consequences of the discriminant "
                "valuation alone. A j-match proves Qbar-isomorphism of nonsingular "
                "fibres, not the submitter's construction method."
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("specification", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw = args.specification.read_bytes()
    result = build(json.loads(raw), raw)
    if args.check:
        if args.output is None:
            raise SystemExit("--check requires --output")
        if args.output.exists():
            expected = json.loads(args.output.read_text())
            # The release string is provenance, not mathematical output.
            result["target"]["pari_minimalization"]["engine"] = expected[
                "target"
            ]["pari_minimalization"]["engine"]
        rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing artifact: {args.output}")
        print(f"PASS {args.output}")
    elif args.output is not None:
        rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output.exists():
            raise SystemExit(f"refusing to overwrite existing output: {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"WROTE {args.output}")
    matches = result["recovery"]["exact_j_matches"]
    print(
        "RECOVERY "
        f"family={result['family']['id']} target={result['target']['label']} "
        f"bounded_parameters={result['recovery']['unique_bounded_parameters']} "
        f"exact_j_matches={[row['parameter'] for row in matches]}"
    )


if __name__ == "__main__":
    main()
