#!/usr/bin/env sage-python
"""M29 wrapper for the 4096-parity scan of one remaining strict direction."""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
M29 = ART / "curve302_m28_fibre_orbit4761_extension1772_mod2_v1.json"
IMPLEMENTATION = CAS / "audit_curve302_m28_one_residual_fibre.sage"
ALLOWED = {"residual-strict-03", "residual-strict-07"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output(direction: str) -> Path:
    return ART / ("curve302_m29_" + direction.replace("-", "_") + "_fibre_v1.json")


def convert_keys(value):
    if isinstance(value, list):
        return [convert_keys(item) for item in value]
    if isinstance(value, dict):
        return {key.replace("M28", "M29"): convert_keys(item) for key, item in value.items()}
    return value


def build(direction: str):
    impl = importlib.machinery.SourceFileLoader("curve302_m29_one_impl", str(IMPLEMENTATION)).load_module()
    impl.M28, impl.RANK, impl.EXTENSION_DIMENSION = M29, 29, 12
    payload = convert_keys(impl.build(direction))
    payload["schema"] = "elliptic-curves.curve302-m29-one-residual-fibre.v1"
    payload["status"] = "PASS_RETROSPECTIVE_M29_ONE_RESIDUAL_FIBRE_SCAN"
    payload["inputs"][str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    payload["protocol"]["parities"] = 4096
    payload["reproducing_command"] = "sage -python elliptic-curves/cas/audit_curve302_m29_one_residual_fibre.sage --check --direction " + direction
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direction", choices=sorted(ALLOWED), required=True)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    payload = build(args.direction)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    path = output(args.direction)
    if args.build:
        if path.exists():
            raise FileExistsError("preserve immutable M29 one-direction diagnostic")
        path.write_text(rendered)
    elif path.read_text() != rendered:
        raise ArithmeticError("stored M29 diagnostic did not replay")
    print("CURVE302M29ONE|{}|parities=4096|winner_digits={}".format(args.direction, payload["winner"]["reduced_coordinate_decimal_digits"]), flush=True)


if __name__ == "__main__":
    main()
