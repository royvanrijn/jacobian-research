#!/usr/bin/env sage-python
"""Certify the 2-primary parity reduction for product twists.

For a product character ``d`` put ``A=E(K(sqrt(d)))``, and let ``M`` and
``N`` be the invariant and anti-invariant Mordell--Weil lattices.  The exact
identity certified here is

    Hhat^-1(<sigma>,A) = N / (1-sigma)A = N / (2 B_-),

where ``B_-=((1-sigma)/2)A``.  Thus the quotient is ``(N/2N)/Gamma_-``,
with ``Gamma_-`` the anti-invariant projection of the integral character
glue.  The Kummer form of the same test is

    [T]=0 iff delta_L(T)=delta_L(tau) for some tau in E(K)
          iff T+tau is divisible by two in E(L).

The product-twist anti-invariant lattices are not known, so this script does
not assert their Tate cohomology.  For any height-eight ``T`` whose Tate class
is zero, it exhausts the invariant R17(2) trace parities and isolates the deep
trace parities not covered by the norm-eight bisection inversion.  Separately,
it records the two glue possibilities under the genuinely additional
hypothesis that ``N=<16>`` is primitive of rank one.  This keeps the remaining
search finite without changing UNKNOWN to a rank or vanishing theorem.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import runpy

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
ALTERNATE = (
    ROOT
    / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
)
ORBIT_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.json"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
NORM8_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json"
)
NORM8_TABLE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
NORM10_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.json"
)
NORM10_TABLE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.tsv"
)
INVERSION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.json"
)
RANKS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json"
)
DIMENSION = 17


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def mask_vector(mask: int) -> tuple[int, ...]:
    return tuple((mask >> index) & 1 for index in range(DIMENSION))


def exact_norm(entries, gram) -> int:
    value = vector(ZZ, entries)
    return int(value * gram * value)


def parse_masks(path: Path) -> set[int]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    masks = {int(row["orbit_mask"]) for row in rows}
    if len(masks) != len(rows):
        raise ArithmeticError(f"duplicate parity mask in {path}")
    return masks


class CosetOracle:
    """fplll CVP decisions with exact integral norm checks."""

    def __init__(self, gram, *, float_type="dd", precision=160):
        if float_type == "mpfr":
            FPLLL.set_precision(precision)
        self.gram = gram
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix([[int(entry) for entry in row] for row in gram.rows()]),
            gram=True,
            float_type=float_type,
            update=True,
        )
        self.mu = [
            [self.gso.get_mu(i, j) if i > j else 0.0 for j in range(DIMENSION)]
            for i in range(DIMENSION)
        ]
        # The digit vector itself is a representative.  This deliberately
        # generous bound prevents an incorrect empty CVP ball.
        self.distance_bound = (
            sum(abs(int(entry)) for row in gram.rows() for entry in row) / 4 + 1
        )

    def solve(self, mask: int) -> tuple[int, tuple[int, ...], float]:
        residue = mask_vector(mask)
        target = [
            -(
                residue[index]
                + sum(
                    residue[row] * self.mu[row][index]
                    for row in range(index + 1, DIMENSION)
                )
            )
            / 2
            for index in range(DIMENSION)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0, DIMENSION, self.distance_bound, 0, target=target
        )
        if not solutions:
            raise ArithmeticError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(entry)) for entry in coordinates)
        if any(
            abs(float(entry) - integer) > 1.0e-7
            for entry, integer in zip(coordinates, closest)
        ):
            raise ArithmeticError("CVP returned nonintegral lattice coordinates")
        representative = tuple(
            residue[index] + 2 * closest[index] for index in range(DIMENSION)
        )
        if any(
            (representative[index] - residue[index]) % 2
            for index in range(DIMENSION)
        ):
            raise ArithmeticError("CVP representative is in the wrong parity coset")
        norm = exact_norm(representative, self.gram)
        error = abs(4 * float(reported_distance) - norm)
        if error > 1.0e-7 or norm < 0 or norm % 2:
            raise ArithmeticError(
                f"invalid CVP output norm={norm}, distance error={error}"
            )
        return norm, representative, error


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    alternate = json.loads(ALTERNATE.read_text())
    orbit_certificate = json.loads(ORBIT_CERTIFICATE.read_text())
    direct = json.loads(DIRECT.read_text())
    norm8_certificate = json.loads(NORM8_CERTIFICATE.read_text())
    norm10_certificate = json.loads(NORM10_CERTIFICATE.read_text())
    inversion = json.loads(INVERSION.read_text())
    ranks = json.loads(RANKS.read_text())

    if orbit_certificate.get("status") != "PASS_ALTERNATE_ROOTLESS_LATTICE_BISECTION_ORBITS":
        raise ArithmeticError("alternate rootless orbit certificate is not exact")
    if norm8_certificate.get("status") != "PASS_EXACT_COMPLETE_ALTERNATE_NORM8_PENCIL_PRIORITY":
        raise ArithmeticError("norm-eight trace certificate is not complete")
    if norm10_certificate.get("status") != "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATION_PRIORITY":
        raise ArithmeticError("norm-ten trace certificate is not complete")
    if inversion.get("status") != "PASS_EXACT_COMPLETE_NORM8_BISECTION_INVERSION":
        raise ArithmeticError("product bisection inversion is not complete")
    if inversion["search"]["squareclass_hit_count"] != 0:
        raise ArithmeticError("product bisection inversion now has a hit")
    if direct["sections"].get("status") != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("direct invariant section basis is not saturated")

    historical_gram = matrix(ZZ, alternate["rootless_frame"])
    short_change = matrix(
        ZZ, orbit_certificate["input"]["row_short_basis_change"]
    )
    short_gram = short_change * historical_gram * short_change.transpose()
    direct_gram = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    historical_to_direct = matrix(
        ZZ, direct["frame_certificate"]["integral_isometry_to_alternate_Q80"]
    )
    section_coordinates = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    if historical_gram.nrows() != DIMENSION or historical_gram.det() != 948:
        raise ArithmeticError("alternate frame is not the determinant-948 rank-17 lattice")
    if abs(short_change.det()) != 1 or abs(section_coordinates.det()) != 1:
        raise ArithmeticError("coordinate chain is not unimodular")
    if (
        historical_to_direct
        * historical_gram
        * historical_to_direct.transpose()
        != direct_gram
    ):
        raise ArithmeticError("historical/direct frame isometry changed")
    short_to_direct = short_change * historical_to_direct.inverse()
    short_to_section = short_to_direct * section_coordinates.inverse()
    if any(entry.denominator() != 1 for entry in short_to_section):
        raise ArithmeticError("short-to-section transport is not integral")
    short_to_direct = matrix(ZZ, short_to_direct)
    short_to_section = matrix(ZZ, short_to_section)

    # Completeness of the 49-class complement is exact and does not depend on
    # floating CVP decisions.  Re-enumerate the cheap norm <=6 shells with the
    # repository's Fincke--Pohst traversal (whose signed count is checked by
    # PARI), and consume the already-certified complete norm-8 and norm-10
    # tables for the two large shells.
    enumerator = runpy.run_path(str(ENUMERATOR))
    short_shells = enumerator["streaming_short_vectors"](short_gram, bound=6)
    norm4_masks = set(short_shells["masks_by_norm"][4])
    norm6_masks = set(short_shells["masks_by_norm"][6])
    if len(norm4_masks) != 1313 or len(norm6_masks) != 26645:
        raise ArithmeticError("unexpected exact norm-4/norm-6 parity shells")
    norm8_masks = parse_masks(NORM8_TABLE)
    norm10_masks = parse_masks(NORM10_TABLE)
    if len(norm8_masks) != 63917 or len(norm10_masks) != 39147:
        raise ArithmeticError("unexpected exact norm-eight/norm-ten table sizes")
    classified_sets = [
        {0},
        norm4_masks,
        norm6_masks,
        norm8_masks,
        norm10_masks,
    ]
    if any(
        classified_sets[left] & classified_sets[right]
        for left in range(len(classified_sets))
        for right in range(left)
    ):
        raise ArithmeticError("certified minimum-norm parity shells overlap")
    classified = set().union(*classified_sets)
    if len(classified) != 131023:
        raise ArithmeticError("exact shells do not classify the expected 131023 cosets")
    deep_masks = sorted(set(range(1 << DIMENSION)) - classified)
    if len(deep_masks) != 49:
        raise ArithmeticError("exact shell complement does not contain 49 parities")

    # CVP is used only to discover a norm-12 witness in each exact complement
    # class.  Exact integral recomputation proves the upper bound 12; exclusion
    # from all exact shells through 10 proves that this is the minimum.
    oracle = CosetOracle(short_gram)
    representatives: dict[int, tuple[int, ...]] = {}
    maximum_distance_error = 0.0
    for mask in deep_masks:
        norm, representative, error = oracle.solve(mask)
        if norm != 12:
            raise ArithmeticError("deep parity has no returned norm-12 witness")
        representatives[mask] = representative
        maximum_distance_error = max(maximum_distance_error, error)
    histogram = Counter(
        {
            0: 1,
            4: len(norm4_masks),
            6: len(norm6_masks),
            8: len(norm8_masks),
            10: len(norm10_masks),
            12: len(deep_masks),
        }
    )

    quadratic_values = []
    for mask in range(1 << DIMENSION):
        residue_norm = exact_norm(mask_vector(mask), short_gram)
        if residue_norm % 2:
            raise ArithmeticError("alternate lattice stopped being even")
        quadratic_value = (residue_norm // 2) % 2
        quadratic_values.append(quadratic_value)
    quadratic_histogram = Counter(quadratic_values)
    if quadratic_histogram != Counter({0: 65280, 1: 65792}):
        raise ArithmeticError("unexpected F2 quadratic-form value distribution")

    category_by_mask = {}
    for norm, masks in (
        (0, {0}),
        (4, norm4_masks),
        (6, norm6_masks),
        (8, norm8_masks),
        (10, norm10_masks),
        (12, set(deep_masks)),
    ):
        for mask in masks:
            category_by_mask[mask] = norm
            if quadratic_values[mask] != (norm // 2) % 2:
                raise ArithmeticError("quadratic parity disagrees with a minimum shell")
    audit = CosetOracle(short_gram, float_type="mpfr", precision=256)
    audit_masks = deep_masks
    maximum_audit_error = 0.0
    for mask in audit_masks:
        norm, unused_representative, error = audit.solve(mask)
        if norm != 12:
            raise ArithmeticError("MPFR audit disagrees with double-double CVP")
        maximum_audit_error = max(maximum_audit_error, error)

    deep_records = []
    for mask in deep_masks:
        short_vector = vector(ZZ, representatives[mask])
        # Sign does not change a parity class.  Normalize the displayed CVP
        # witness without making any uniqueness claim for shortest lifts.
        if tuple(map(int, -short_vector)) < tuple(map(int, short_vector)):
            short_vector = -short_vector
        historical_vector = short_vector * short_change
        direct_vector = short_vector * short_to_direct
        section_vector = short_vector * short_to_section
        if any(
            exact_norm(value, gram) != 12
            for value, gram in (
                (short_vector, short_gram),
                (historical_vector, historical_gram),
                (direct_vector, direct_gram),
                (section_vector, matrix(ZZ, direct["sections"]["height_gram"])),
            )
        ):
            raise ArithmeticError("deep trace coordinate transport changed its norm")
        deep_records.append(
            {
                "orbit_mask": mask,
                "orbit_hex": f"0x{mask:05x}",
                "minimum_norm": 12,
                "quadratic_value_mod_2": 0,
                "binary_short_basis_parity": list(mask_vector(mask)),
                "short_basis_minimum_witness": list(map(int, short_vector)),
                "historical_alternate_w": list(map(int, historical_vector)),
                "direct_alternate_w": list(map(int, direct_vector)),
                "section_basis_w": list(map(int, section_vector)),
            }
        )

    exact_rank_one_keys = {
        row["pair_key"]
        for row in ranks["results"]
        if row.get("status") == "completed"
        and int(row.get("rank_lower_bound", -1)) == 1
        and int(row.get("rank_upper_bound", -1)) == 1
    }
    target_keys = {row["pair_key"] for row in inversion["targets"]}
    if len(target_keys) != 17 or target_keys != exact_rank_one_keys:
        raise ArithmeticError("the product targets are not the seventeen exact rank-one bases")

    payload = {
        "schema": "elkies-k3.r17-norm12-11952-product-tate-parity.v1",
        "status": "PASS_EXACT_PRODUCT_TATE_PARITY_REDUCTION",
        "mathematical_status": "UNKNOWN_PRODUCT_TWIST_TATE_COHOMOLOGY",
        "theorem": {
            "setup": (
                "For A=E(L), M=A^+, N=A^-, and B_-=((1-sigma)/2)A, the "
                "integral character glue projects isomorphically to B_-/N."
            ),
            "tate_cohomology": (
                "Hhat^-1(<sigma>,A)=N/(1-sigma)A=N/(2B_-)="
                "(N/2N)/Gamma_-."
            ),
            "kummer_criterion": (
                "For T in N, [T]=0 iff there is tau in M with T+tau in 2A; "
                "equivalently delta_L(T)=delta_L(tau)."
            ),
            "exponent": 2,
            "torsion_hypothesis": (
                "The pulled-back rootless 48I1 elliptic surface has no nonzero "
                "torsion section by the height formula, so the lattice calculation "
                "has no hidden 2-torsion term."
            ),
        },
        "conditional_rank_one_anti_lattice": {
            "hypothesis": (
                "The anti-invariant lattice generated in the displayed rational span "
                "is primitive rank one N=<16>, generated by a twist-height-eight T "
                "after degree-two pullback."
            ),
            "possibilities": [
                {
                    "anti_glue_projection": "0",
                    "tate_cohomology": "Z/2",
                    "unique_nonzero_class_representative": "T modulo 2N",
                },
                {
                    "anti_glue_projection": "N/2N",
                    "tate_cohomology": "0",
                    "glue_identity": "2R=T+tau for an invariant trace tau",
                },
            ],
            "not_an_actual_rank_statement": True,
        },
        "invariant_trace_parity": {
            "lattice": "alternate-Q80 Mordell-Weil lattice M, pulled back as M(2)",
            "rank": DIMENSION,
            "determinant": int(historical_gram.det()),
            "total_classes": 1 << DIMENSION,
            "quadratic_form": "q(tau mod 2M)=height_M(tau)/2 mod 2",
            "quadratic_value_histogram": {
                str(key): value for key, value in sorted(quadratic_histogram.items())
            },
            "minimum_norm_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "nonzero_isotropic_classes": quadratic_histogram[0] - 1,
            "isotropic_minimum_norm_partition": {
                "4": histogram[4],
                "8": histogram[8],
                "12": histogram[12],
            },
            "deep_norm12_class_count": len(deep_records),
            "deep_norm12_classes": deep_records,
            "coordinate_chain": (
                "short -> historical alternate -> direct compiled frame -> "
                "saturated equation section basis"
            ),
        },
        "interaction_with_completed_inversion": {
            "product_target_count": len(target_keys),
            "product_pair_keys": sorted(target_keys),
            "norm8_trace_classes_tested": len(norm8_masks),
            "norm8_squareclass_hits": inversion["search"]["squareclass_hit_count"],
            "rootless_height_lower_bound": 8,
            "height_identity": (
                "Writing r=R.O on the rootless chi=4 pullback, orthogonality gives "
                "4*height_L(R)=height_L(T)+height_L(tau)=16+2*height_K(tau), "
                "while height_L(R)=8+2r; hence height_K(tau)=8+4r."
            ),
            "norm4_glue_exclusion": (
                "If height_M(tau)=4 and pulled height(T)=16, then "
                "height((T+tau)/2)=6, below the rootless chi=4 minimum 8."
            ),
            "residual_zero_class_trace_parities": len(deep_records),
            "residual_trace_target_pairs": len(deep_records) * len(target_keys),
            "residual_geometry": (
                "For a norm-12 trace, R=(T+tau)/2 has pulled-back height 10 and "
                "R.O=1. These are deep genus-one bisection carriers, not members "
                "of the completed P.O=0/norm-eight inversion."
            ),
            "no_other_trace_minima": (
                "The exact M/2M covering spectrum has maximum minimum norm 12. "
                "Together with height_K(tau)=8+4(R.O), only minima 8 and 12 can "
                "occur for a half-point carrier; minimum 8 is already inverted and "
                "minimum 12 is the 49-class residual."
            ),
        },
        "direct_polynomial_parity_encoding": {
            "height_eight_candidate": (
                "A direct chi=4 polynomial solution T must first satisfy the exact "
                "twist Weierstrass equations and the local height-eight component gates."
            ),
            "zero_class_disjunction": (
                "Its Tate class is zero only if at least one of the 49 stored trace "
                "parities admits half-point variables R over L satisfying 2R=T+tau."
            ),
            "nonzero_class_conjunction": (
                "Its Tate class is nonzero exactly when all 49 Kummer equalities "
                "delta_L(T)=delta_L(tau) fail."
            ),
            "finite_field_use": (
                "Finite fields may be used to reject individual half-point systems. "
                "Only surviving systems may be lifted and verified over QQ; a local "
                "failure is not a standalone computation of the global Tate quotient."
            ),
            "case_index": "Cartesian product of product_pair_keys and deep_norm12_classes",
        },
        "numerical_certificate": {
            "completeness_algorithm": (
                "exact Fincke-Pohst/PARI shells through norm 6 plus the pinned "
                "complete norm-8 and norm-10 tables; their complement has 49 classes"
            ),
            "exact_classified_through_norm10": len(classified),
            "witness_algorithm": (
                "fplll CVP followed by exact integral norm recomputation; floating "
                "decisions are not used for complement completeness"
            ),
            "double_double_norm12_witnesses": len(deep_masks),
            "maximum_double_double_distance_error": maximum_distance_error,
            "mpfr_precision_bits": 256,
            "mpfr_audited_cosets": len(audit_masks),
            "mpfr_audit_includes_all_deep_norm12_classes": True,
            "maximum_mpfr_distance_error": maximum_audit_error,
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                ENUMERATOR,
                ALTERNATE,
                ORBIT_CERTIFICATE,
                DIRECT,
                NORM8_CERTIFICATE,
                NORM8_TABLE,
                NORM10_CERTIFICATE,
                NORM10_TABLE,
                INVERSION,
                RANKS,
            )
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "fpylll_float_types": ["dd", "mpfr-256"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_norm12_11952_product_tate_parity.sage --check"
        ),
        "proof_boundary": {
            "proved": (
                "The Tate/graph-glue identity, its Kummer parity criterion, the "
                "conditional primitive anti-rank-one classification, the complete "
                "alternate-M/2M minimum spectrum, and the rank-independent reduction "
                "of any possible zero-class height-eight carrier to 49 deep trace "
                "parities per target."
            ),
            "not_proved": (
                "No anti-invariant product-twist Mordell-Weil lattice N, integral "
                "base-change glue Gamma_-, Tate quotient H_d, height-eight product "
                "section, finite-field survivor classification, or characteristic-zero "
                "lift is supplied. All seventeen product-twist questions remain UNKNOWN."
            ),
        },
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "PRODUCTTATEPARITY|"
        f"targets={len(target_keys)}|deep_trace_parities={len(deep_records)}|"
        f"residual_cases={len(target_keys) * len(deep_records)}|"
        f"status={payload['status']}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
