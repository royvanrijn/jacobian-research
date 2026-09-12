#!/usr/bin/env python3
"""Replay the characteristic-zero mixed quartic--sextic obstruction."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
LOW_LAYER_CHECKER = (
    ROOT / "verify_hc4_meng_four_low_layer_quartics.py"
)
LINE_CHECKER = ROOT / "explore_hc4_meng_mixed_quartic_sextic.py"
PLANE_CHECKER = (
    ROOT / "verify_hc4_meng_mixed_quartic_sextic_planes_qq.py"
)
EXPECTED_HASHES = {
    LOW_LAYER_CHECKER: (
        "1e874888b9db32d58c4de418ed79b86d3942157dde5b747e47256671b3bbfb85"
    ),
    LINE_CHECKER: (
        "d458900769afd3ed564be082bbe673ff8bad0043cdf0911039a3821b173241e0"
    ),
    PLANE_CHECKER: (
        "2b80517d8a1d9386e0b70055c8456d0f6ab8c1ec0d4fe25db9d711ebf402e324"
    ),
}


def replay(path: Path, *arguments: str) -> None:
    result = subprocess.run(
        [sys.executable, str(path), *arguments],
        check=False,
    )
    assert result.returncode == 0, (
        f"{path.name} failed with exit code {result.returncode}"
    )


for path, expected in EXPECTED_HASHES.items():
    assert path.is_file()
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual == expected, (
        f"stale component hash for {path.name}: {actual}"
    )

replay(LOW_LAYER_CHECKER)
replay(LINE_CHECKER, "--lines-only")
replay(PLANE_CHECKER)

print(
    "PASS: mixed sparse quartic--sextic line and plane strata are "
    "excluded over QQ"
)
print(
    "SCOPE: the 234 quartic principal parts from HC4MQ1, with no cubic "
    "and a zero-gradient sextic supported on at most four monomials"
)
