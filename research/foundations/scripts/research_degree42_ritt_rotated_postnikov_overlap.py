#!/usr/bin/env python3
"""Research completed Postnikov overlap on a rotated degree-42 half-braid.

For the chart selected by ``--word``, let

    I = thick-path ideal, J = boundary ideal, L = thin-path ideal,
    K = graph ideal.

The exact containment ``J <= L + mJ`` proves ``J^ = L^`` by Nakayama,
because ``L <= J`` is tautological from the residual construction.

The quadratic overlap

    (J intersect (I + K^2)) / (I + KJ)

is tested modulo the square of the Dickson-base ideal.  A zero remainder
after an explicit maximal-adic cutoff is an Artin--Rees certificate, not a
stabilization heuristic.  This script is exploratory until word-specific
cutoffs and outputs are pinned by a verifier.
"""

from __future__ import annotations

import argparse
import re
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


def rotated_postnikov_overlap_audit(
    word: tuple[int, int, int],
    sector_cutoff: int = 4,
    overlap_cutoff: int = 5,
) -> dict[str, object]:
    """Run exact Nakayama and Artin--Rees containment tests."""

    (
        parameters,
        factor_variables,
        thick,
        thin,
        boundary,
        thick_omission,
        thin_omission,
    ) = rotated_source_ideal_data(word)
    normals, base_coordinates, images = graph_normal_map(
        word, factor_variables
    )
    local_variables = normals + base_coordinates
    map_images = ",".join(
        str(images[parameter]).replace("**", "^")
        for parameter in parameters
    )
    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal ILsource={serialize_ideal(thin)};
ideal IBsource={serialize_ideal(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IL=phi(ILsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
ideal baseIdeal={",".join(map(str, base_coordinates))};
ideal maximalIdeal={",".join(map(str, local_variables))};

ideal thinNakayamaDenominatorStd=std(IL+maximalIdeal*IB);
ideal thinNakayamaRemainder=simplify(
  reduce(IB,thinNakayamaDenominatorStd),2
);
ideal tautologicalRemainder=simplify(reduce(IL,std(IB)),2);
print("THIN_BOUNDARY_NAKAYAMA");
print(size(thinNakayamaRemainder));
print(size(tautologicalRemainder));

ideal sectorDenominator=IT+K*IB+baseIdeal^2*IB;
ideal sectorDenominatorStd=std(sectorDenominator);
ideal sectorCutoffPower=maximalIdeal^{sector_cutoff};
ideal sectorCutoffIntersection=intersect(IB,sectorCutoffPower);
ideal sectorCutoffRemainder=simplify(
  reduce(sectorCutoffIntersection,sectorDenominatorStd),2
);
print("SECTOR_SOURCE");
print(size(sectorCutoffIntersection));
print(size(sectorCutoffRemainder));
print(vdim(std(sectorDenominator+sectorCutoffPower)));
print(vdim(std(IB+sectorCutoffPower)));

ideal overlapNumerator=intersect(IB,IT+K^2);
ideal overlapDenominator=IT+K*IB+baseIdeal^2*overlapNumerator;
ideal overlapDenominatorStd=std(overlapDenominator);
ideal overlapCutoffPower=maximalIdeal^{overlap_cutoff};
ideal overlapCutoffIntersection=intersect(
  overlapNumerator,overlapCutoffPower
);
ideal overlapCutoffRemainder=simplify(
  reduce(overlapCutoffIntersection,overlapDenominatorStd),2
);
print("QUADRATIC_OVERLAP");
print(size(overlapNumerator));
print(size(overlapCutoffIntersection));
print(size(overlapCutoffRemainder));
print(vdim(std(overlapDenominator+overlapCutoffPower)));
print(vdim(std(overlapNumerator+overlapCutoffPower)));
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
    compact = " ".join(
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip() and not line.strip().startswith("//")
    )
    thin_match = re.search(
        r"THIN_BOUNDARY_NAKAYAMA ([0-9]+) ([0-9]+)", compact
    )
    sector_match = re.search(
        r"SECTOR_SOURCE ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)",
        compact,
    )
    overlap_match = re.search(
        r"QUADRATIC_OVERLAP "
        r"([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)",
        compact,
    )
    assert thin_match is not None, result.stdout
    assert sector_match is not None, result.stdout
    assert overlap_match is not None, result.stdout
    thin_values = tuple(map(int, thin_match.groups()))
    sector_values = tuple(map(int, sector_match.groups()))
    overlap_values = tuple(map(int, overlap_match.groups()))
    return {
        "word": word,
        "thick_composite_omission": thick_omission,
        "thin_prime_omission": thin_omission,
        "thin_boundary_nakayama": {
            "boundary_remainder_generators": thin_values[0],
            "tautological_reverse_remainder_generators": thin_values[1],
            "completed_equality": thin_values == (0, 0),
        },
        "sector_source_mod_base_square": {
            "artin_rees_cutoff": sector_cutoff,
            "cutoff_intersection_generators": sector_values[0],
            "cutoff_remainder_generators": sector_values[1],
            "denominator_quotient_length": sector_values[2],
            "numerator_quotient_length": sector_values[3],
            "dimension": sector_values[2] - sector_values[3],
        },
        "quadratic_overlap_mod_base_square": {
            "numerator_generators": overlap_values[0],
            "artin_rees_cutoff": overlap_cutoff,
            "cutoff_intersection_generators": overlap_values[1],
            "cutoff_remainder_generators": overlap_values[2],
            "denominator_quotient_length": overlap_values[3],
            "numerator_quotient_length": overlap_values[4],
            "dimension": overlap_values[3] - overlap_values[4],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--word", choices=("237", "327"), required=True)
    parser.add_argument("--sector-cutoff", type=int, default=4)
    parser.add_argument("--overlap-cutoff", type=int, default=5)
    arguments = parser.parse_args()
    word = tuple(int(character) for character in arguments.word)
    print(
        rotated_postnikov_overlap_audit(
            word,
            sector_cutoff=arguments.sector_cutoff,
            overlap_cutoff=arguments.overlap_cutoff,
        )
    )


if __name__ == "__main__":
    main()
