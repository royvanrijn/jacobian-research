#!/usr/bin/env python3
"""Discover target-curve families in declared generated construction spaces."""

from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Sequence


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
CAS_ROOT = PROGRAM_ROOT / "cas"
sys.path[:0] = [str(PROGRAM_ROOT), str(CAS_ROOT)]

from ecsearch.family_discovery import (  # noqa: E402
    DiscoveryTarget,
    PolynomialWeierstrassFamily,
    SixRootMestreFamily,
    discover_target_families,
)
from icarm_curve245_mestre import fermigier_roots  # noqa: E402
from mestre_root_tuples import (  # noqa: E402
    SixRootMestreConstruction,
    normalize_integer_root_tuple,
)


def _resolve(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else REPOSITORY_ROOT / path


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(path)


def _load_json(path: Path) -> Any:
    if path.name.endswith(".json.gz"):
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def _nested(value: Any, path: Sequence[str | int]) -> Any:
    for key in path:
        value = value[key]
    return value


def load_targets(specification: dict[str, Any]) -> tuple[DiscoveryTarget, ...]:
    targets = []
    for target_spec in specification["targets"]:
        path = _resolve(target_spec["record_file"])
        source = _load_json(path)
        if "records_path" in target_spec:
            records = _nested(source, target_spec["records_path"])
            match = target_spec["match"]
            selected = [
                record
                for record in records
                if record[match["field"]] == match["value"]
            ]
            if len(selected) != 1:
                raise ValueError(
                    f"target {target_spec['label']} selector found {len(selected)} rows"
                )
            record = selected[0]
            coefficients = _nested(record, target_spec.get("ainvs_path", ["ainvs"]))
        else:
            coefficients = _nested(source, target_spec["ainvs_path"])
        targets.append(
            DiscoveryTarget(
                label=target_spec["label"],
                coefficients=tuple(Fraction(str(value)) for value in coefficients),
                metadata={
                    "source": _display_path(path),
                    "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                },
            )
        )
    labels = [target.label for target in targets]
    if len(labels) != len(set(labels)):
        raise ValueError("target labels must be unique")
    return tuple(targets)


def load_polynomial_families(
    specifications: Sequence[dict[str, Any]],
) -> list[PolynomialWeierstrassFamily]:
    families = []
    for family_spec in specifications:
        path = _resolve(family_spec["file"])
        source = _load_json(path)
        model = _nested(
            source,
            family_spec.get("model_path", ["canonical_weierstrass_model"]),
        )
        coefficients = {
            name: tuple(int(value) for value in polynomial)
            for name, polynomial in model["coefficients_low_to_high"].items()
        }
        family = PolynomialWeierstrassFamily(
            identifier=source.get("id", family_spec.get("id", path.stem)),
            coefficient_polynomials=coefficients,
            parameter_name=family_spec.get("parameter_name", "t"),
            metadata={"source": _display_path(path)},
        )
        declared_discriminant = model.get("discriminant_coefficients_low_to_high")
        if declared_discriminant is not None and tuple(
            map(int, declared_discriminant)
        ) != family.invariant_polynomials["discriminant"]:
            raise ValueError(
                f"declared discriminant does not match model in {_display_path(path)}"
            )
        families.append(family)
    return families


def load_mestre_roots(
    specification: dict[str, Any]
) -> tuple[tuple[tuple[tuple[int, ...], dict[str, Any]], ...], dict[str, int]]:
    origins_by_roots: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    emission_count = 0
    generator_parameter_pairs_tested = 0
    generator_degenerate_pairs = 0

    def add(roots: tuple[int, ...], origin: dict[str, Any]) -> None:
        nonlocal emission_count
        emission_count += 1
        origins_by_roots.setdefault(roots, []).append(origin)

    for source_spec in specification.get("census_sources", []):
        path = _resolve(source_spec["file"])
        payload = _load_json(path)
        roots = _nested(payload, source_spec["roots_path"])
        for values in roots:
            add(
                tuple(map(int, values)),
                {
                    "source": _display_path(path),
                    "population": "/".join(map(str, source_spec["roots_path"])),
                },
            )

    for generator_spec in specification.get("generators", []):
        if generator_spec["kind"] != "fermigier-six-root":
            raise ValueError(f"unsupported six-root generator {generator_spec['kind']}")
        for u_text in generator_spec["u_values"]:
            for v_text in generator_spec["v_values"]:
                generator_parameter_pairs_tested += 1
                u, v = Fraction(str(u_text)), Fraction(str(v_text))
                generated = fermigier_roots(u, v)
                if any(value.denominator != 1 for value in generated):
                    raise ValueError("integer Fermigier generator produced nonintegral roots")
                try:
                    roots = normalize_integer_root_tuple(
                        value.numerator for value in generated
                    )
                except ValueError:
                    generator_degenerate_pairs += 1
                    continue
                add(
                    roots,
                    {
                        "source": _display_path(
                            CAS_ROOT / "icarm_curve245_mestre.py"
                        ),
                        "construction": "fermigier_roots",
                        "generator_kind": generator_spec["kind"],
                        "generator_parameters": {
                            "u": str(u_text),
                            "v": str(v_text),
                        },
                    },
                )

    for explicit in specification.get("explicit_families", []):
        add(
            tuple(map(int, explicit["roots"])),
            {
                "source": "explicit control in discovery specification",
                "label": explicit["label"],
            },
        )

    result = tuple(
        (roots, {"origins": origins})
        for roots, origins in origins_by_roots.items()
    )
    return result, {
        "emitted_family_count_before_deduplication": emission_count,
        "duplicate_family_emission_count": emission_count - len(result),
        "generator_parameter_pairs_tested": generator_parameter_pairs_tested,
        "generator_degenerate_parameter_pairs": generator_degenerate_pairs,
    }


def load_mestre_families(
    specification: dict[str, Any]
) -> tuple[list[SixRootMestreFamily], dict[str, int]]:
    families = []
    roots_with_metadata, generation_summary = load_mestre_roots(specification)
    for roots, metadata in roots_with_metadata:
        construction = SixRootMestreConstruction(tuple(roots))
        if not construction.is_quartic_family:
            raise AssertionError(f"census tuple is not a quartic family: {roots}")
        families.append(
            SixRootMestreFamily(
                roots=roots,
                exact_model=construction.primitive_jacobian_coefficients,
                metadata=metadata,
            )
        )
    return families, generation_summary


def build(specification: dict[str, Any], raw_specification: bytes) -> dict[str, Any]:
    targets = load_targets(specification)
    families = []
    families.extend(
        load_polynomial_families(specification.get("polynomial_families", []))
    )
    mestre_specification = specification.get("six_root_mestre")
    mestre_count = 0
    mestre_generation_summary: dict[str, int] = {}
    if mestre_specification:
        mestre_families, mestre_generation_summary = load_mestre_families(
            mestre_specification
        )
        mestre_count = len(mestre_families)
        families.extend(mestre_families)
    identifiers = [family.identifier for family in families]
    if not families or len(identifiers) != len(set(identifiers)):
        raise ValueError("family identifiers must be nonempty and unique")
    modular_primes = tuple(int(value) for value in specification["modular_primes"])
    records = [
        discover_target_families(
            target,
            families,
            modular_primes=modular_primes,
            fingerprint_prime_bound=int(
                specification.get("fingerprint_prime_bound", 97)
            ),
        )
        for target in targets
    ]
    return {
        "schema": "elliptic-curves.generated-family-discovery.v1",
        "claim_level": "complete bounded exact computation",
        "specification_sha256": hashlib.sha256(raw_specification).hexdigest(),
        "construction_space": {
            "polynomial_weierstrass_family_count": len(
                specification.get("polynomial_families", [])
            ),
            "six_root_mestre_family_count": mestre_count,
            "six_root_generation_summary": mestre_generation_summary,
            "total_family_count": len(families),
            "family_identifier_sha256": hashlib.sha256(
                "\n".join(identifiers).encode()
            ).hexdigest(),
        },
        "targets": records,
        "interpretation": {
            "proved": (
                "Every family in the declared generated construction space was "
                "screened. Modular no-root witnesses are exact exclusions; every "
                "survivor was factored over QQ and re-specialized exactly."
            ),
            "boundary": (
                "A no-match result excludes only the declared bounded generators. "
                "It does not exclude a larger root census, another construction "
                "template, a different fibration, an isogeny, or a private family."
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("specification", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    raw = arguments.specification.read_bytes()
    result = build(json.loads(raw), raw)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output is None:
            raise SystemExit("--check requires --output")
        if not arguments.output.exists() or arguments.output.read_text() != rendered:
            raise SystemExit(f"stale or missing artifact: {arguments.output}")
        print(f"PASS {arguments.output}")
    elif arguments.output is not None:
        if arguments.output.exists():
            raise SystemExit(f"refusing to overwrite existing output: {arguments.output}")
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        print(f"WROTE {arguments.output}")
    summary = {
        target["target"]: [
            (match["family_id"], match["parameter"])
            for match in target["q_isomorphism_matches"]
        ]
        for target in result["targets"]
    }
    print(f"DISCOVERY {summary}")


if __name__ == "__main__":
    main()
