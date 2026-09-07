"""Staged JSON compiler for the global low-degree support census."""

from __future__ import annotations

import gc
import shutil
import subprocess
from collections import Counter, defaultdict
from typing import Any, Callable, Sequence

import sympy as sp
import z3

from .low_degree_census import (
    Support,
    SupportEnumerator,
    bucket_summary,
    collision_support_counts_by_size,
    collision_support_space_counts,
    degree_profiles_below,
    dense_quadratic_collision_status,
    enumerate_valuation_faces,
    groebner_status,
    profile_rank_function,
    sha256_json,
    sign_smt,
    singular_batch_status,
)


STAGE_FILENAMES = (
    "global_low_degree_census_01_profiles.json",
    "global_low_degree_census_02_supports.json",
    "global_low_degree_census_03_buckets.json",
    "global_low_degree_census_04_valuations.json",
    "global_low_degree_census_05_smt.json",
    "global_low_degree_census_06_modular.json",
    "global_low_degree_census_07_exact.json",
    "global_low_degree_census_08_boundary.json",
)
MANIFEST_FILENAME = "global_low_degree_census_manifest.json"


def _software_versions() -> dict[str, str | None]:
    singular = shutil.which("Singular")
    singular_version = None
    if singular is not None:
        completed = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode == 0 and completed.stdout:
            singular_version = completed.stdout.splitlines()[0]
    return {
        "sympy": sp.__version__,
        "z3": z3.get_version_string(),
        "singular": singular_version,
    }


def _representative_orbits(
    supports_by_size: dict[int, tuple[Support, ...]],
) -> tuple[
    tuple[Support, ...],
    dict[int, tuple[dict[str, Any], ...]],
]:
    representatives: list[Support] = []
    by_size: dict[int, tuple[dict[str, Any], ...]] = {}
    for size, supports in supports_by_size.items():
        grouped: dict[Support, set[Support]] = defaultdict(set)
        for support in supports:
            grouped[support.canonical_under_collision_stabilizer()].add(support)
        rows = []
        for representative, members in sorted(grouped.items()):
            representatives.append(representative)
            rows.append(
                {
                    "support_id": representative.identifier,
                    "support": representative.to_json(),
                    "orbit_size": len(members),
                    "member_ids": sorted(member.identifier for member in members),
                    "second_member_rule": (
                        None if len(members) == 1 else "simultaneous x2<->x3, F2<->F3"
                    ),
                }
            )
        by_size[size] = tuple(rows)
    return tuple(representatives), by_size


def _common_header(
    schema: str,
    *,
    max_degree: int,
    target_profile: Sequence[int],
    support_bound: int,
) -> dict[str, Any]:
    return {
        "schema": schema,
        "parameters": {
            "dimension": 3,
            "maximum_coordinate_degree": max_degree,
            "target_invariant_degree_profile": list(target_profile),
            "maximum_nonlinear_support": support_bound,
            "collision_axis": 1,
            "normalization": ["F(0)=0", "JF(0)=I", "F(e1)=0"],
        },
    }


def build_low_degree_census(
    *,
    max_degree: int = 7,
    target_profile: Sequence[int] = (7, 6, 4),
    support_bound: int = 6,
    primes: Sequence[int] = (11, 13, 17),
    progress: Callable[[str], None] | None = None,
) -> dict[str, dict[str, Any]]:
    """Build all eight deterministic stage artifacts in memory."""

    target = tuple(map(int, target_profile))
    if len(target) != 3 or tuple(sorted(target, reverse=True)) != target:
        raise ValueError("target_profile must be a nonincreasing triple")
    if target[0] > max_degree:
        raise ValueError(
            "max_degree must cover the largest degree occurring below target_profile"
        )

    def announce(message: str) -> None:
        if progress is not None:
            progress(message)

    announce("stage 1/8: invariant degree profiles")
    profiles = degree_profiles_below(target)
    stage1 = {
        **_common_header(
            "global-low-degree-census.profiles.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "profile_notation": {
            "name": "filtered output-span degree flag",
            "definition": (
                "For V=span(F1,F2,F3), rank(C_{>k}) equals the number of "
                "profile entries greater than k. This is invariant under affine "
                "source changes and invertible affine target recombination."
            ),
            "warning": (
                "Sorted degrees of the three coordinates in a fixed collision frame "
                "are not affine-invariant. Exact support alone cannot decide top-degree "
                "cancellation between output rows; rank minors are an algebra-stage gate."
            ),
        },
        "profile_count": len(profiles),
        "profiles": [list(profile) for profile in profiles],
        "profile_rank_gates": [
            {
                "profile": list(profile),
                "rank_above_degree": {
                    str(threshold): profile_rank_function(profile, threshold)
                    for threshold in range(max_degree + 1)
                },
            }
            for profile in profiles
        ],
        "complete": True,
    }

    announce("stage 2/8: exact sparse support closure")
    enumerator = SupportEnumerator(max_degree, support_bound)
    supports_by_size, search_statistics = enumerator.enumerate()
    representatives, orbit_rows = _representative_orbits(supports_by_size)
    labelled_counts = collision_support_counts_by_size(max_degree, support_bound)
    balanced_counts = {
        size: len(supports_by_size.get(size, ()))
        for size in range(1, support_bound + 1)
    }
    stage2 = {
        **_common_header(
            "global-low-degree-census.supports.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "full_boolean_support_space": collision_support_space_counts(max_degree),
        "collision_admissible_supports_by_size": {
            str(size): count for size, count in labelled_counts.items()
        },
        "determinant_balanced_supports_by_size": {
            str(size): count for size, count in balanced_counts.items()
        },
        "determinant_balanced_orbits_by_size": {
            str(size): len(orbit_rows.get(size, ()))
            for size in range(1, support_bound + 1)
        },
        "orbits": {
            str(size): list(orbit_rows.get(size, ()))
            for size in range(1, support_bound + 1)
        },
        "search_statistics": search_statistics,
        "completeness": {
            "status": "complete",
            "scope": (
                f"Every exact normalized collision support of coordinate degree at most "
                f"{max_degree} and at most {support_bound} nonlinear monomial occurrences."
            ),
            "proof_invariant": (
                "A singleton active determinant bucket must acquire a second contributing "
                "triple. All such triples are branched. Once balanced, every optional next "
                "monomial is branched. Temporary one-axis collision supports in F2/F3 are "
                "repaired by every possible second pure-axis monomial."
            ),
            "not_claimed": "No support-cardinality-unbounded enumeration is claimed.",
        },
    }

    announce("stage 3/8: determinant bucket ledger")
    stage3 = {
        **_common_header(
            "global-low-degree-census.buckets.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "universal_bucket_hypergraph": enumerator.universe_statistics,
        "formula": {
            "bucket_exponent": "alpha+beta+gamma-(1,1,1)",
            "triple_multiplier": "det(rows alpha,beta,gamma)",
            "singleton_rule": (
                "A nonconstant bucket with one nonzero determinant triple is a nonzero "
                "coefficient monomial on the exact-support torus and cannot vanish."
            ),
        },
        "representative_count": len(representatives),
        "representatives": [bucket_summary(support) for support in representatives],
        "all_representatives_have_no_singleton": all(
            bucket_summary(support)["singleton_bucket_count"] == 0
            for support in representatives
        ),
    }

    # The closure cache contains tens of thousands of transient Support objects.
    del enumerator
    gc.collect()

    announce("stage 4/8: finite exposed-face valuation census")
    valuation_rows = []
    valuation_counts: Counter[int] = Counter()
    for support in representatives:
        valuations = enumerate_valuation_faces(support)
        valuation_counts[len(valuations)] += 1
        valuation_rows.append(
            {
                "support_id": support.identifier,
                "valuation_class_count": len(valuations),
                "classes": list(valuations),
                "classes_sha256": sha256_json(valuations),
            }
        )
    stage4 = {
        **_common_header(
            "global-low-degree-census.valuations.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "finite_equivalence": (
            "Weights are equivalent when they have the same active coordinate stratum "
            "and expose the same face in each Newton support. Primitive integer "
            "representatives are minimized in L1 norm."
        ),
        "necessary_bounded_image_rule": (
            "If m_i(w)>0, the exposed support of F_i has at least two monomials. "
            "All seven coordinate strata are included. This is necessary, not sufficient."
        ),
        "representative_count_by_valuation_class_count": {
            str(count): multiplicity for count, multiplicity in sorted(valuation_counts.items())
        },
        "supports_without_candidate_valuation": [
            row["support_id"] for row in valuation_rows if not row["classes"]
        ],
        "representatives": valuation_rows,
        "completeness": "complete modulo exposed-face equivalence within the sparse scope",
    }

    announce("stage 5/8: coefficient-sign SMT")
    sign_rows = [sign_smt(support) for support in representatives]
    sign_counts = Counter(row["status"] for row in sign_rows)
    stage5 = {
        **_common_header(
            "global-low-degree-census.smt.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "solver": f"Z3 {z3.get_version_string()} exact Boolean arithmetic",
        "decision_problem": (
            "Can coefficient signs make every determinant bucket and every nontrivial "
            "collision sum contain both signs?"
        ),
        "scope_warning": (
            "This is a necessary test over ordered fields. It is not a complex-field "
            "exclusion and it does not replace polynomial coefficient algebra."
        ),
        "status_counts": dict(sorted(sign_counts.items())),
        "representatives": sign_rows,
    }

    announce("stage 6/8: modular coefficient-torus algebra")
    modular_by_prime: dict[str, list[dict[str, Any]]] = {}
    modular_summary: dict[str, dict[str, int]] = {}
    for prime in primes:
        rows = list(singular_batch_status(representatives, int(prime)))
        modular_by_prime[str(prime)] = rows
        modular_summary[str(prime)] = {
            "unit_ideal": sum(row["unit_ideal"] for row in rows),
            "isolated_or_zero_dimensional": sum(
                not row["unit_ideal"] and row["dimension"] == 0 for row in rows
            ),
            "positive_dimensional": sum(
                not row["unit_ideal"]
                and row["dimension"] is not None
                and row["dimension"] > 0
                for row in rows
            ),
        }
    stage6 = {
        **_common_header(
            "global-low-degree-census.modular.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "backend": _software_versions()["singular"],
        "primes": list(map(int, primes)),
        "summary": modular_summary,
        "representatives_by_prime": modular_by_prime,
        "logical_status": (
            "routing evidence only: an empty special fibre at a fixed prime does not by "
            "itself prove that the characteristic-zero generic fibre is empty"
        ),
    }

    announce("stage 7/8: exact rational Groebner replay")
    singular_exact = list(singular_batch_status(representatives, 0))
    sympy_exact = [groebner_status(support) for support in representatives]
    if [row["unit_ideal"] for row in singular_exact] != [
        row["unit_ideal"] for row in sympy_exact
    ]:
        raise AssertionError("Singular and SymPy disagree on an exact support ideal")
    dense_quadratic = dense_quadratic_collision_status()
    profile_set = set(profiles)
    eliminated_profiles = [
        row
        for row in dense_quadratic["eliminated_degree_profiles"]
        if tuple(row) in profile_set
    ]
    exact_survivors = [
        row["support_id"] for row in singular_exact if not row["unit_ideal"]
    ]
    if exact_survivors:
        sparse_floor: int | None = None
        sparse_scope = (
            f"Exact coefficient systems survive within support at most {support_bound}; "
            "this run emits no global nonlinear-support lower bound."
        )
    else:
        sparse_floor = support_bound + 1
        sparse_scope = (
            f"Every normalized characteristic-zero Keller collision of raw coordinate "
            f"degree at most {max_degree} has at least {sparse_floor} nonlinear "
            "monomial occurrences. Hence the same lower bound holds for every invariant "
            f"degree profile below {target}. No support of size {sparse_floor} is "
            "constructed or asserted to exist."
        )
    stage7 = {
        **_common_header(
            "global-low-degree-census.exact.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "backends": {
            "primary": _software_versions()["singular"],
            "independent_replay": f"SymPy {sp.__version__} exact Groebner basis",
        },
        "representative_count": len(representatives),
        "unit_ideal_count": sum(row["unit_ideal"] for row in singular_exact),
        "surviving_support_ids": exact_survivors,
        "singular_results": singular_exact,
        "sympy_results_sha256": sha256_json(sympy_exact),
        "dense_quadratic_collision_ideal": dense_quadratic,
        "completely_eliminated_profiles": eliminated_profiles,
        "sparse_frontier_theorem": {
            "lower_bound_nonlinear_support": sparse_floor,
            "attainment_proved": False,
            "scope": sparse_scope,
        },
    }

    announce("stage 8/8: source/coefficient boundary ledger")
    stage8 = {
        **_common_header(
            "global-low-degree-census.boundary.v1",
            max_degree=max_degree,
            target_profile=target,
            support_bound=support_bound,
        ),
        "source_infinity": {
            "audited_by": STAGE_FILENAMES[3],
            "coordinate_strata": 7,
            "status": "all exposed-face classes compiled before coefficient algebra",
        },
        "coefficient_boundary": {
            "compactification_used": False,
            "reason": (
                "No coefficient was normalized to one. Each exact support was checked on "
                "its full coefficient torus using rho*product(c)-1. Setting a coefficient "
                "to zero is exactly a smaller support, and every smaller support through "
                "the sparse bound is a separate census row."
            ),
            "downward_support_boundary_complete": True,
        },
        "projective_charts_requiring_audit_after_exact_lifting": (
            0 if not exact_survivors else None
        ),
        "supports_requiring_component_boundary_audit": exact_survivors,
        "projective_boundary_status": (
            "not required: exact coefficient-torus ideals are empty"
            if not exact_survivors
            else "pending for exact-algebra survivors"
        ),
        "survives_projective_boundary_analysis": (
            False if not exact_survivors else None
        ),
        "scope_warning": (
            f"This closes only the support-at-most-{support_bound} stratum and the dense "
            f"quadratic profiles. It is not a complete census of all supports below {target}."
        ),
    }

    payloads = {
        STAGE_FILENAMES[0]: stage1,
        STAGE_FILENAMES[1]: stage2,
        STAGE_FILENAMES[2]: stage3,
        STAGE_FILENAMES[3]: stage4,
        STAGE_FILENAMES[4]: stage5,
        STAGE_FILENAMES[5]: stage6,
        STAGE_FILENAMES[6]: stage7,
        STAGE_FILENAMES[7]: stage8,
    }
    return payloads
