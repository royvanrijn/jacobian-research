#!/usr/bin/env python3
"""Audit q8 representative consistency before finite/infinity gluing.

This compares the dominant-D13, abstract-nef, and actual component-nef
physical root targets.  For the component-nef target it also converts the
E8/II* and both E6/IV* chart orientations into valuation cycles, computes the
minimal common fibre shift making each cycle nonnegative, and prints the
resulting monomial valuation envelopes.  These envelopes are diagnostics, not
complete local modules when a branch cancellation is possible.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
DEFAULT_OUTPUT = LOCAL / "q8-finite-infinity-target-audit.json"

II_PHYSICAL_TO_CHART = (2, 7, 6, 4, 8, 3, 5, 1)
II_U = (2, 2, 4, 6, 3, 4, 3, 5)
II_X = (2, 4, 6, 10, 4, 7, 5, 8)
II_Y = (3, 5, 9, 15, 6, 10, 8, 12)

IV_PHYSICAL_TO_CHART = {
    "minus": (2, 4, 6, 5, 3, 1),
    "plus": (3, 5, 6, 4, 2, 1),
}
IV_U = (2, 1, 1, 2, 2, 3)
IV_X = (2, 2, 2, 3, 3, 4)
IV_Y = (3, 2, 2, 4, 4, 6)


def load(path: Path, expected_status: str | None = None) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if expected_status is not None and payload.get("status") != expected_status:
        raise ValueError(f"{path}: expected status {expected_status}, got {payload.get('status')}")
    return payload


def chart_cycle(physical: Iterable[int], mapping: tuple[int, ...]) -> list[int]:
    physical = list(map(int, physical))
    if len(physical) != len(mapping):
        raise ValueError("cycle/mapping dimension mismatch")
    result: list[int | None] = [None] * len(mapping)
    for physical_index, chart_index in enumerate(mapping):
        result[chart_index - 1] = physical[physical_index]
    if any(value is None for value in result):
        raise ValueError("mapping is not a permutation")
    return [int(value) for value in result]


def minimal_fibre_shift(cycle: list[int], fibre_values: tuple[int, ...]) -> tuple[int, list[int]]:
    if len(cycle) != len(fibre_values) or any(value <= 0 for value in fibre_values):
        raise ValueError("invalid fibre valuation cycle")
    shift = max(
        [0]
        + [(-coefficient + value - 1) // value for coefficient, value in zip(cycle, fibre_values) if coefficient < 0]
    )
    shifted = [coefficient + shift * value for coefficient, value in zip(cycle, fibre_values)]
    if any(value < 0 for value in shifted):
        raise AssertionError("fibre shift did not clear the negative cycle")
    return shift, shifted


def dominates(left: tuple[int, int, int], right: tuple[int, int, int]) -> bool:
    return all(a <= b for a, b in zip(left, right))


def monomial_envelope(
    target: list[int],
    u_values: tuple[int, ...],
    x_values: tuple[int, ...],
    y_values: tuple[int, ...],
) -> list[tuple[int, int, int]]:
    bound = max(8, max(target, default=0) + 4)
    generators: list[tuple[int, int, int]] = []
    for y_exponent in (0, 1):
        for u_exponent in range(bound + 1):
            for x_exponent in range(bound + 1):
                valuation = [
                    u_exponent * u + x_exponent * x + y_exponent * y
                    for u, x, y in zip(u_values, x_values, y_values)
                ]
                if any(value < required for value, required in zip(valuation, target)):
                    continue
                candidate = (u_exponent, x_exponent, y_exponent)
                if any(dominates(existing, candidate) for existing in generators):
                    continue
                generators = [existing for existing in generators if not dominates(candidate, existing)]
                generators.append(candidate)
    if not generators:
        raise RuntimeError(f"no monomial envelope found through exponent bound {bound}")
    return sorted(generators, key=lambda item: (sum(item), item[2], item[1], item[0]))


def monomial_label(exponents: tuple[int, int, int]) -> str:
    names = ("u", "X", "Y")
    pieces = []
    for name, exponent in zip(names, exponents):
        if exponent == 1:
            pieces.append(name)
        elif exponent > 1:
            pieces.append(f"{name}^{exponent}")
    return "1" if not pieces else "*".join(pieces)


def summarize_target(payload: dict[str, Any]) -> dict[str, Any]:
    selected = payload["selected_q8"]
    normalization = payload["normalization"]
    return {
        "representative": normalization["representative"],
        "effective_reflection_count": len(normalization.get("effective_component_reflections", [])),
        "vertical_fibre_coefficient": int(selected["vertical_fibre_coefficient"]),
        "E6_cycle": list(map(int, selected["E6"]["vertical_cycle"])),
        "E6_degrees": list(map(int, selected["E6"]["component_degrees"])),
        "E6_affine_degree": int(selected["E6"]["affine_component_degree"]),
        "E8_cycle": list(map(int, selected["E8"]["vertical_cycle"])),
        "E8_degrees": list(map(int, selected["E8"]["component_degrees"])),
        "E8_affine_degree": int(selected["E8"]["affine_component_degree"]),
    }


def local_signature(summary: dict[str, Any]) -> tuple[Any, ...]:
    return (
        summary["vertical_fibre_coefficient"],
        tuple(summary["E6_cycle"]),
        tuple(summary["E6_degrees"]),
        summary["E6_affine_degree"],
        tuple(summary["E8_cycle"]),
        tuple(summary["E8_degrees"]),
        summary["E8_affine_degree"],
    )


def finite_mode(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
    payload = load(path, "PASS_EXACT_FINITE_Q_CONDITIONS")
    first = payload["module"]["exact_basis"][0][0]
    mode = "nef" if first == "1" else "dominant-d13"
    return {
        "path": str(path),
        "mode": mode,
        "basis": payload["module"]["exact_basis"],
        "smith_degrees": payload["module"]["smith_degrees"],
    }


def fractional_target(
    cycle: list[int],
    mapping: tuple[int, ...],
    u_values: tuple[int, ...],
    x_values: tuple[int, ...],
    y_values: tuple[int, ...],
) -> dict[str, Any]:
    chart = chart_cycle(cycle, mapping)
    shift, shifted = minimal_fibre_shift(chart, u_values)
    generators = monomial_envelope(shifted, u_values, x_values, y_values)
    return {
        "chart_cycle": chart,
        "common_fibre_denominator_power": shift,
        "shifted_nonnegative_cycle": shifted,
        "monomial_valuation_envelope": [monomial_label(item) for item in generators],
        "fractional_form": f"I/({('u^' + str(shift)) if shift != 1 else 'u'})" if shift else "I",
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--dominant", type=Path, default=LOCAL / "q8-target-dominant-d13.json")
parser.add_argument("--nef", type=Path, default=LOCAL / "q8-target-nef.json")
parser.add_argument("--component-nef", type=Path, default=LOCAL / "q8-target-component-nef.json")
parser.add_argument("--finite-dominant", type=Path, default=LOCAL / "q8-finite-dominant.json")
parser.add_argument("--finite-nef", type=Path, default=LOCAL / "q8-finite-nef.json")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

for name in ("dominant", "nef", "component_nef", "finite_dominant", "finite_nef", "output"):
    value = getattr(args, name)
    setattr(args, name, value.resolve() if value is not None else None)

targets = {
    "dominant-d13": summarize_target(load(args.dominant, "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET")),
    "nef": summarize_target(load(args.nef, "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET")),
    "component-nef": summarize_target(load(args.component_nef, "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET")),
}
component = targets["component-nef"]

same_as_dominant = local_signature(component) == local_signature(targets["dominant-d13"])
same_as_nef = local_signature(component) == local_signature(targets["nef"])

finite_records = [record for record in (
    finite_mode(args.finite_dominant), finite_mode(args.finite_nef)
) if record is not None]
for record in finite_records:
    record["compatible_with_component_nef"] = (
        same_as_nef if record["mode"] == "nef" else same_as_dominant
    )

component_ii = fractional_target(
    component["E8_cycle"], II_PHYSICAL_TO_CHART, II_U, II_X, II_Y
)
component_iv = {
    sign: fractional_target(component["E6_cycle"], mapping, IV_U, IV_X, IV_Y)
    for sign, mapping in IV_PHYSICAL_TO_CHART.items()
}
iv_common_shift = (
    component_iv["minus"]["common_fibre_denominator_power"]
    if component_iv["minus"]["common_fibre_denominator_power"]
    == component_iv["plus"]["common_fibre_denominator_power"]
    else None
)

# At the smooth old-base infinity, finite E6/E8 components are absent.  In the
# unshifted generic frame the signed requirement is -a, where a is the fibre
# coefficient: positive means vanishing, negative means an allowed pole.
raw_infinity_order = -component["vertical_fibre_coefficient"]
common_finite_denominator_shift = (
    component_ii["common_fibre_denominator_power"] + iv_common_shift
    if iv_common_shift is not None else None
)
shifted_frame_infinity_order = (
    raw_infinity_order - common_finite_denominator_shift
    if common_finite_denominator_shift is not None else None
)

compatible_finite = any(record["compatible_with_component_nef"] for record in finite_records)
status = (
    "READY_TO_DERIVE_EXACT_INFINITY_INTERSECTION"
    if compatible_finite else "NEED_COMPONENT_NEF_FINITE_MODULE"
)

payload = {
    "schema": "elkies-k3.q8-finite-infinity-target-audit.v1",
    "status": status,
    "targets": targets,
    "comparison": {
        "component_nef_same_local_signature_as_dominant": same_as_dominant,
        "component_nef_same_local_signature_as_nef": same_as_nef,
        "existing_finite_modules": finite_records,
    },
    "component_nef_fractional_local_targets": {
        "II_star": component_ii,
        "IV_star": component_iv,
        "warning": (
            "The monomial envelopes are exact valuation envelopes.  If the two IV* "
            "orientations differ, a branch binomial such as Y+-c*u^2 may be needed "
            "before this becomes the complete local module."
        ),
    },
    "infinity": {
        "vertical_fibre_coefficient": component["vertical_fibre_coefficient"],
        "raw_signed_required_order": raw_infinity_order,
        "interpretation": (
            "positive=must vanish; zero=regular; negative=may have a pole"
        ),
        "common_finite_denominator_shift": common_finite_denominator_shift,
        "signed_order_in_common_shifted_numerator_frame": shifted_frame_infinity_order,
    },
    "next_gate": (
        "Derive the complete component-nef II*/IV* coefficient module in the "
        "same translated chord frame, then intersect it with the certified smooth "
        "infinity lattice."
        if not compatible_finite else
        "Use the compatible finite module and the signed infinity order to compute "
        "the exact polynomial-module intersection; no degree-window search is needed."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

for name, summary in targets.items():
    print(
        "Q8FININF|rep={}|reflections={}|vf={}|E6_cycle={}|E6_deg={}|E6_aff={}|"
        "E8_cycle={}|E8_deg={}|E8_aff={}".format(
            name,
            summary["effective_reflection_count"],
            summary["vertical_fibre_coefficient"],
            ",".join(map(str, summary["E6_cycle"])),
            ",".join(map(str, summary["E6_degrees"])),
            summary["E6_affine_degree"],
            ",".join(map(str, summary["E8_cycle"])),
            ",".join(map(str, summary["E8_degrees"])),
            summary["E8_affine_degree"],
        )
    )
print(
    "Q8FININF|component_same_dominant={}|component_same_nef={}|finite_compatible={}".format(
        int(same_as_dominant), int(same_as_nef), int(compatible_finite)
    )
)
print(
    "Q8FININF|II_shift={}|II_cycle={}|II_shifted={}|II_envelope={}".format(
        component_ii["common_fibre_denominator_power"],
        ",".join(map(str, component_ii["chart_cycle"])),
        ",".join(map(str, component_ii["shifted_nonnegative_cycle"])),
        ",".join(component_ii["monomial_valuation_envelope"]),
    )
)
for sign in ("minus", "plus"):
    record = component_iv[sign]
    print(
        "Q8FININF|IV_orientation={}|shift={}|cycle={}|shifted={}|envelope={}".format(
            sign,
            record["common_fibre_denominator_power"],
            ",".join(map(str, record["chart_cycle"])),
            ",".join(map(str, record["shifted_nonnegative_cycle"])),
            ",".join(record["monomial_valuation_envelope"]),
        )
    )
print(
    "Q8FININF|raw_infinity_order={}|finite_denominator_shift={}|"
    "shifted_frame_order={}|status={}".format(
        raw_infinity_order,
        "ambiguous" if common_finite_denominator_shift is None else common_finite_denominator_shift,
        "ambiguous" if shifted_frame_infinity_order is None else shifted_frame_infinity_order,
        status,
    )
)
