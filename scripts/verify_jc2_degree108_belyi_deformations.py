#!/usr/bin/env python3
"""Verify the pinned dessin-first degree-(72,108) deformation certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CAS = ROOT / "plane-jc" / "cas"
CORE = CAS / "jc2_degree108_belyi_deformations.py"
EXPECTED_CORE_SHA256 = "2930cfba18577b9986ca224d2dce8b40f65effae16cedd0b945932d6db2b36d6"
assert hashlib.sha256(CORE.read_bytes()).hexdigest() == EXPECTED_CORE_SHA256
sys.path.insert(0, str(CAS))

from jc2_degree108_belyi_deformations import OUT, compile_report  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="skip replaying the terminal five-parameter Singular unit-ideal check",
    )
    parser.add_argument("--certificate", type=Path, default=OUT)
    parser.add_argument("--singular-timeout", type=int, default=300)
    arguments = parser.parse_args()

    pinned = json.loads(arguments.certificate.read_text(encoding="utf-8"))
    run_full = not arguments.quick
    replay = compile_report(run_full, arguments.singular_timeout)
    if not run_full:
        replay["terminal_open_check"] = pinned["terminal_open_check"]
    if replay != pinned:
        raise RuntimeError("pinned Belyi-deformation certificate differs from exact replay")

    assert pinned["coefficient_field"]["galois_group"] == "S5"
    assert pinned["dessins"]["rotation_orbits"] == 5
    assert pinned["top_reconstruction"]["passport"] == [
        [2] * 10 + [1],
        [3] * 7,
        [17, 1, 1, 1, 1],
    ]
    assert pinned["deformations"]["linear_maps"]["BE"]["kernel_dimension"] == 2
    assert pinned["deformations"]["linear_maps"]["CF"]["kernel_dimension"] == 3
    assert pinned["deformations"]["G_constant_kernel_dimension"] == 1
    assert pinned["terminal_open_check"]["unit_ideal"] is True
    print("JC2_DEGREE108_BELYI_DEFORMATIONS_PASS")
    print(f"MODE={'full' if run_full else 'quick'}")
    print(f"CERTIFICATE={arguments.certificate.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
