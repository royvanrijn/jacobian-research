#!/usr/bin/env python3
"""Thin unchanged-method adapter for E29 and ICARM 398--400."""

from pathlib import Path

import run_frozen_wgxli_latent_lattice as frozen


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"

frozen.FREEZE = ARTIFACTS / "latent_lattice_secondary_method_freeze_v1.json"
frozen.OUTPUT = ARTIFACTS / "latent_lattice_secondary_frozen_dimension_v1.json"
frozen.EXPECTED_FREEZE_SHA256 = (
    "8795cdd203ba1c698e0f0534c14a45a91c37a6f4d405e795a9f2f295f86bfcba"
)
frozen.EXPECTED_TAG = "LATENT-LATTICE-E29-398-400-FROZEN-2026-09-01-v1"
frozen.OUTPUT_SCHEMA = (
    "elliptic-curves.latent-lattice-secondary-frozen-dimension.v1"
)
frozen.ADAPTER_PATH = Path(__file__).resolve()
frozen.TARGET_CURVE_IDS = (12, 398, 399, 400)
frozen.TARGET_SOURCES = {
    12: "206ef6992f433155d349618d55c289a00cf9014eb222da059c57e0db76131c0e",
    398: "5e09b5ed49cde24d20fcf300794e58a47f7f75ac7bba98c92b68ff3654df49f4",
    399: "92125a3aafd44ff45028ade826ef62338667164a402b3264507816e3c2009ead",
    400: "704f2292b4395923e2887a4ddd7a35f03d46baed42f5c75aac2dfe8519dd4275",
}


if __name__ == "__main__":
    frozen.main()
