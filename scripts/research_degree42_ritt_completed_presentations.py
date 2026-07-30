#!/usr/bin/env python3
"""Extract completed two-variable presentations for rotated Ritt conormals."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree42_ritt_rotated_conormal_flags import (  # noqa: E402
    graph_normal_map,
    rotated_source_ideal_data,
    serialize_ideal,
)


def presentation_cache(word: tuple[int, int, int]) -> Path:
    """Return the compressed completed-presentation cache."""

    label = "".join(map(str, word))
    return (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"degree42_ritt_completed_presentations_{label}.txt.gz"
    )


def completed_presentation_output(word: tuple[int, int, int]) -> str:
    """Return two seven-generator presentations over ``Q[tau,zeta]``."""

    cache = presentation_cache(word)
    if cache.is_file():
        with gzip.open(cache, "rt") as source:
            return source.read()
    (
        parameters,
        factor_variables,
        thick,
        _thin,
        boundary,
        _thick_omission,
        _thin_omission,
    ) = rotated_source_ideal_data(word)
    normals, base_coordinates, images = graph_normal_map(
        word, factor_variables
    )
    map_images = ",".join(
        str(images[parameter]).replace("**", "^")
        for parameter in parameters
    )
    normal_module = ",".join(f"[{normal}]" for normal in normals)
    base_map = ",".join(
        ["0"] * len(normals) + [str(value) for value in base_coordinates]
    )
    program = f"""
ring src=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal IBsource={serialize_ideal(boundary)};
ring q=0,({",".join(map(str, normals + base_coordinates))}),(dp({len(normals)}),dp(2));
map phi=src,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IB=phi(IBsource);
module MK={normal_module};
module RT=modulo(MK,IT);
module RS=modulo(MK,IB);
ideal normalIdeal={",".join(map(str, normals))};
module NT=RT+normalIdeal*freemodule({len(normals)});
module NS=RS+normalIdeal*freemodule({len(normals)});
ring b=0,(tau,zeta),dp;
map psi=q,{base_map};
option(redSB);
module PT=std(psi(NT));
module PS=std(psi(NS));
print("PRESENTATION_TOTAL");
for(int i=1;i<=ncols(PT);i++){{print(PT[i]);}}
print("PRESENTATION_SPECTATOR");
for(int j=1;j<=ncols(PS);j++){{print(PS[j]);}}
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    markers = ("PRESENTATION_TOTAL", "PRESENTATION_SPECTATOR")
    missing = tuple(
        marker for marker in markers if marker not in result.stdout
    )
    if missing or "? error" in result.stdout:
        raise RuntimeError(
            "Singular returned an invalid completed presentation; "
            f"missing={missing}, stderr={result.stderr.strip()}, "
            f"stdout_tail={result.stdout[-2000:]}"
        )
    with gzip.open(cache, "wt", compresslevel=9) as target:
        target.write(result.stdout)
    return result.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--word", choices=("237", "327"), required=True)
    arguments = parser.parse_args()
    word = tuple(int(character) for character in arguments.word)
    print(completed_presentation_output(word))


if __name__ == "__main__":
    main()
