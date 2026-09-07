#!/usr/bin/env sage-python
"""M27 wrapper for the exhaustive ten-new-bit residual strict fibre diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
M27 = ART / "curve302_m26_fibre_orbit58145_extension472_mod2_v1.json"
OUTPUT = ART / "curve302_m27_remaining_strict_fibre_scan_v1.json"
IMPLEMENTATION = CAS / "audit_curve302_m26_remaining_strict_fibre_scan.sage"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    impl = importlib.machinery.SourceFileLoader("curve302_m27_fibre_impl", str(IMPLEMENTATION)).load_module()
    impl.M26 = M27
    impl.RANK = 27
    impl.EXTENSION_DIMENSION = 10
    # The retained M26 implementation asserts a five-row panel.  Scan the
    # already recovered strict-02 row as a sealed diagnostic control, then
    # remove it from the M27 residual report below.
    impl.EXCLUDED = {"residual-strict-04", "residual-strict-06"}
    payload = impl.build()
    payload["schema"] = "elliptic-curves.curve302-m27-remaining-strict-fibre-scan.v1"
    payload["status"] = "PASS_RETROSPECTIVE_M27_PARITY_FIBRE_SCAN"
    payload["inputs"][str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    payload["protocol"].update(
        directions=4,
        parities_per_direction=1024,
        fixed_M17_shell_projection="hold each original M17 chart parity fixed and enumerate every newly available M27 parity bit",
        centre_and_target_reduction="DD/MPFR-agreeing CVP in a rank-27 rounded 384-bit canonical-height metric",
    )
    payload["reproducing_command"] = "sage -python elliptic-curves/cas/audit_curve302_m27_remaining_strict_fibre_scan.sage --check"
    payload["directions"] = [row for row in payload["directions"] if row["id"] != "residual-strict-02"]
    for direction in payload["directions"]:
        direction["M27_extension_count"] = direction.pop("M26_extension_count")
        direction["M27_extension_winner"] = direction.pop("M26_extension_winner")
        direction["M27_extension_coordinate_trials"] = direction.pop("M26_extension_coordinate_trials")
        winner = direction["M27_extension_winner"]
        trials = direction["M27_extension_coordinate_trials"]
        winner["extension_mask_in_new_M27_coordinates"] = winner.pop("extension_mask_in_new_M26_coordinates")
        winner["M27_parity_mask"] = winner.pop("M26_parity_mask")
        winner["centre_M27_word"] = winner.pop("centre_M26_word")
        winner["target_translation_M27_word"] = winner.pop("target_translation_M26_word")
        for trial in trials:
            trial["extension_mask_in_new_M27_coordinates"] = trial.pop("extension_mask_in_new_M26_coordinates")
            trial["M27_parity_mask"] = trial.pop("M26_parity_mask")
            trial["centre_M27_word"] = trial.pop("centre_M26_word")
            trial["target_translation_M27_word"] = trial.pop("target_translation_M26_word")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.build == arguments.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve immutable M27 parity-fibre diagnostic")
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("stored M27 parity-fibre diagnostic did not replay")
    print("CURVE302M27FIBRE|directions=4|parities_per_direction=1024|status=PASS", flush=True)


if __name__ == "__main__":
    main()
