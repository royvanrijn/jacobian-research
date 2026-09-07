#!/usr/bin/env sage-python
"""Certify the projective scaled-shell description of p-neighbour births.

For an even positive lattice K and a good odd prime p, a new norm-two
vector v of the p-neighbour K_ell is replaced by z=p*v in K.  The vector z
has norm 2*p^2 and its nonzero reduction modulo p spans ell.  Conversely
every such z gives a new root of K_ell.  Thus the rootless neighbour lines
are the isotropic quadric with the old-root hyperplane sections and the
projected scaled shell removed.

This checker exhausts every line for every state/prime pair in the three
mass-closed ternary control genera.  It compares the roots predicted from
the parent roots and the single scaled shell with an independent root
enumeration in the materialized child lattice.  The general graph-glue
formula is proved in RANK_MUTATION_AND_LIFT_THEOREMS.md; these controls use
the trivial bridge so that the complete projective strata remain small.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import product
import json
from pathlib import Path
import runpy

from sage.all import (
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
SMALL_SCRIPT = ROOT / "elkies-k3/scripts/analyze_small_genus_defect_graphs.sage"
PLANNER_SCRIPT = ROOT / "elkies-k3/scripts/plan_inverse_ade_targets.sage"
SMALL_ARTIFACT = (
    ROOT / "artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v2.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-inverse-ade-projective-birth-strata-v1.json"
)


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_signed_shell(gram, target_norm):
    """Return every signed integral vector of one exact norm."""

    enumeration = pari(gram).qfminim(target_norm)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    answer = set()
    for representative in representatives:
        if representative * gram * representative != target_norm:
            continue
        answer.add(tuple(map(int, representative)))
        answer.add(tuple(map(int, -representative)))
    return tuple(sorted(answer))


def projective_reduction(value, prime, canonical_projective):
    """Reduce a rational K-coordinate vector modulo a good prime."""

    reduced = []
    for entry in map(QQ, value):
        denominator = int(entry.denominator())
        assert denominator % prime
        reduced.append(
            (int(entry.numerator()) * pow(denominator, -1, prime)) % prime
        )
    return canonical_projective(reduced, prime)


def rational_lattice_basis(generators):
    """Return a square basis for the Z-span of rational row generators."""

    generator_matrix = matrix(QQ, generators)
    denominator = lcm(entry.denominator() for entry in generator_matrix.list())
    integral_generators = matrix(ZZ, denominator * generator_matrix)
    basis = integral_generators.row_module(ZZ).basis_matrix() / denominator
    assert basis.nrows() == basis.ncols() == generator_matrix.ncols()
    return matrix(QQ, basis)


def child_roots_in_parent(gram, prime, line, small):
    """Independently materialize K_ell and return its roots in K tensor QQ."""

    form = small["quadratic_form"](gram)
    child_to_parent = form.find_p_neighbor_from_vec(
        prime, line, return_matrix=True
    ).transpose()
    child_gram = child_to_parent * gram * child_to_parent.transpose()
    roots = small["physical_roots"](child_gram)
    return {
        tuple(vector(QQ, root) * child_to_parent)
        for root in roots
    }


def compile_core_strata(gram, prime, small, canonical_projective):
    """Compile old-root hyperplanes and projected norm-2*p^2 birth points."""

    parent_roots = small["physical_roots"](gram)
    scaled_shell = exact_signed_shell(gram, 2 * prime**2)
    births = defaultdict(set)
    divisible_scaled_roots = 0
    for scaled in scaled_shell:
        key = projective_reduction(scaled, prime, canonical_projective)
        if key is None:
            divisible_scaled_roots += 1
            continue
        line = vector(ZZ, key)
        assert (line * gram * line) % prime == 0
        births[key].add(tuple(QQ(entry) / prime for entry in scaled))
    return parent_roots, scaled_shell, births, divisible_scaled_roots


def nonzero_graph_glue_control(small, canonical_projective):
    """Exercise births in the nonzero glue coset of an index-two completion."""

    prime = 5
    core = 2 * identity_matrix(ZZ, 3)
    bridge = matrix(ZZ, [[2]])
    ambient_gram = block_diagonal_matrix(core, bridge)
    half = QQ(1) / 2
    standard = identity_matrix(ZZ, 4)
    parent_basis = rational_lattice_basis(
        [list(row) for row in standard.rows()] + [[half, half, half, half]]
    )
    parent_gram = parent_basis * ambient_gram * parent_basis.transpose()
    parent_roots = {
        tuple(vector(QQ, root) * parent_basis)
        for root in small["physical_roots"](parent_gram)
    }
    assert len(parent_roots) == 24

    births = defaultdict(set)
    shell_counts = {}
    bridge_cells = ((0, (QQ(-1), QQ(0), QQ(1))), (1, (-half, half)))
    for glue_class, bridge_vectors in bridge_cells:
        for bridge_vector in bridge_vectors:
            core_norm = QQ(2) - 2 * bridge_vector**2
            shell = set()
            for numerators in product(range(-2 * prime, 2 * prime + 1), repeat=3):
                if glue_class:
                    if not all(numerator % 2 for numerator in numerators):
                        continue
                    scaled = vector(QQ, [QQ(numerator) / 2 for numerator in numerators])
                else:
                    scaled = vector(ZZ, numerators)
                if scaled * core * scaled != prime**2 * core_norm:
                    continue
                shell.add(tuple(scaled))
                key = projective_reduction(
                    scaled, prime, canonical_projective
                )
                if key is not None:
                    births[key].add(tuple(list(scaled / prime) + [bridge_vector]))
            shell_counts[f"class_{glue_class}_bridge_{bridge_vector}"] = len(shell)

    lines = small["projective_isotropic_lines"](core, prime)
    assert set(births) == {tuple(line) for line in lines}
    line_records = []
    nonzero_glue_births = 0
    for line in lines:
        key = tuple(map(int, line))
        surviving = set()
        for root in parent_roots:
            core_part = vector(QQ, root[:3])
            pairing = core_part * core * vector(ZZ, line)
            residue = (
                int(pairing.numerator())
                * pow(int(pairing.denominator()), -1, prime)
            ) % prime
            if residue == 0:
                surviving.add(root)

        born = births[key]
        line_nonzero_glue_births = sum(
            1 for root in born if QQ(root[3]).denominator() == 2
        )
        assert len(surviving) == 4
        assert len(born) == 20
        assert line_nonzero_glue_births == 16

        form = small["quadratic_form"](core)
        child_core_basis = form.find_p_neighbor_from_vec(
            prime, line, return_matrix=True
        ).transpose()
        parent_half_class = vector(QQ, [half, half, half])
        child_half_classes = []
        for correction in product(range(prime), repeat=3):
            representative = parent_half_class + vector(QQ, correction) / prime
            pairings = representative * core * child_core_basis.transpose()
            if all(pairing in ZZ for pairing in pairings):
                child_half_classes.append(representative)
        assert len(child_half_classes) == 1
        child_glue = child_half_classes[0]
        child_generators = [
            list(row) + [0] for row in child_core_basis.rows()
        ]
        child_generators.append([0, 0, 0, 1])
        child_generators.append(list(child_glue) + [half])
        child_basis = rational_lattice_basis(child_generators)
        child_gram = child_basis * ambient_gram * child_basis.transpose()
        assert all(entry in ZZ for entry in child_gram.list())
        assert all(entry % 2 == 0 for entry in child_gram.diagonal())
        actual = {
            tuple(vector(QQ, root) * child_basis)
            for root in small["physical_roots"](child_gram)
        }
        assert surviving | born == actual
        assert not surviving & born
        nonzero_glue_births += line_nonzero_glue_births
        line_records.append(
            {
                "line": list(key),
                "old_surviving_signed_roots": len(surviving),
                "born_core_signed_roots": len(born) - line_nonzero_glue_births,
                "born_nonzero_glue_signed_roots": line_nonzero_glue_births,
                "child_signed_roots": len(actual),
            }
        )

    assert nonzero_glue_births == 96
    return {
        "name": "A1^3 plus A1 with diagonal half-class glue",
        "prime": prime,
        "parent_completion_root_system": "D4",
        "isotropic_quadric_points": len(lines),
        "scaled_shell_signed_vectors_by_cell": shell_counts,
        "born_nonzero_glue_signed_roots": nonzero_glue_births,
        "exact_set_comparisons_pass": len(lines),
        "line_root_counts": line_records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    small = runpy.run_path(str(SMALL_SCRIPT))
    planner = runpy.run_path(str(PLANNER_SCRIPT))
    canonical_projective = planner["canonical_projective"]
    source = json.loads(SMALL_ARTIFACT.read_text())

    cases = []
    totals = {
        "state_prime_cases": 0,
        "isotropic_lines": 0,
        "scaled_shell_signed_vectors": 0,
        "projected_birth_stratum_points": 0,
        "predicted_rootless_lines": 0,
        "actual_rootless_lines": 0,
    }
    for genus in source["genera"]:
        edge_index = {
            (row["source"], int(row["prime"]), tuple(row["isotropic_line"])): row
            for row in genus["line_edges"]
        }
        for state in genus["states"]:
            label = state["state"]
            gram = matrix(ZZ, state["gram"])
            for prime in genus["prime_set_optimization"]["prime_set_summaries"][-1][
                "primes"
            ]:
                prime = int(prime)
                assert gram.det() % prime
                parent_roots, scaled_shell, births, divisible = compile_core_strata(
                    gram, prime, small, canonical_projective
                )
                assert divisible == len(parent_roots)
                lines = small["projective_isotropic_lines"](gram, prime)
                predicted_rootless = []
                actual_rootless = []
                line_records = []
                for line in lines:
                    key = tuple(map(int, line))
                    surviving = {
                        tuple(map(QQ, root))
                        for root in parent_roots
                        if int(vector(ZZ, root) * gram * line) % prime == 0
                    }
                    born = births.get(key, set())
                    predicted = surviving | born
                    actual = child_roots_in_parent(gram, prime, line, small)
                    assert predicted == actual
                    if not predicted:
                        predicted_rootless.append(list(key))
                    if not actual:
                        actual_rootless.append(list(key))
                    truth = edge_index[(label, prime, key)]
                    assert len(surviving) == int(
                        truth["old_survivor_signed_root_count"]
                    )
                    assert len(born) == int(truth["new_birth_signed_root_count"])
                    assert len(actual) == int(truth["child_signed_root_count"])
                    line_records.append(
                        {
                            "line": list(key),
                            "old_surviving_signed_roots": len(surviving),
                            "born_signed_roots": len(born),
                            "child_signed_roots": len(actual),
                            "rootless": not actual,
                        }
                    )
                assert predicted_rootless == actual_rootless
                case = {
                    "genus": genus["name"],
                    "source_state": label,
                    "source_gram_sha256": state["gram_sha256"],
                    "prime": prime,
                    "isotropic_quadric_points": len(lines),
                    "old_root_hyperplanes": len(parent_roots) // 2,
                    "scaled_shell_norm": 2 * prime**2,
                    "scaled_shell_signed_vectors": len(scaled_shell),
                    "scaled_shell_vectors_reducing_to_zero": divisible,
                    "projected_birth_stratum_points": len(births),
                    "rootless_lines_from_quadric_complement": predicted_rootless,
                    "line_root_counts": line_records,
                    "exact_set_comparisons_pass": len(lines),
                }
                cases.append(case)
                totals["state_prime_cases"] += 1
                totals["isotropic_lines"] += len(lines)
                totals["scaled_shell_signed_vectors"] += len(scaled_shell)
                totals["projected_birth_stratum_points"] += len(births)
                totals["predicted_rootless_lines"] += len(predicted_rootless)
                totals["actual_rootless_lines"] += len(actual_rootless)

    assert totals["state_prime_cases"] == sum(
        len(genus["states"])
        * len(
            genus["prime_set_optimization"]["prime_set_summaries"][-1]["primes"]
        )
        for genus in source["genera"]
    )
    assert totals["predicted_rootless_lines"] == totals["actual_rootless_lines"]
    graph_glue_control = nonzero_graph_glue_control(small, canonical_projective)

    payload = {
        "schema": "elkies-k3.inverse-ade-projective-birth-strata.v1",
        "status": "PASS_EXHAUSTIVE_PROJECTIVE_BIRTH_STRATA_CONTROLS",
        "theorem": {
            "core_birth_shell": (
                "A born norm-two vector v of K_ell is equivalent to z=p*v in K "
                "with z^2=2*p^2 and nonzero projective reduction [z mod p]=ell."
            ),
            "rootless_locus": (
                "Q_p^iso minus the union of old-root hyperplane sections and "
                "the projective reduction of the norm-2*p^2 shell."
            ),
            "graph_glue_extension": (
                "For bridge vector c and transported graph class a, replace K by "
                "the dual coset [z]=p*a and 2*p^2 by p^2*(2-c^2)."
            ),
        },
        "inputs": {
            relative(SMALL_ARTIFACT): digest(SMALL_ARTIFACT),
            relative(SMALL_SCRIPT): digest(SMALL_SCRIPT),
            relative(PLANNER_SCRIPT): digest(PLANNER_SCRIPT),
        },
        "control_scope": {
            "mass_closed_ternary_genera": len(source["genera"]),
            "uses_marked_target_core": False,
            "uses_target_isometry_class": False,
            "target_predicate": "rootless",
            "comparison": (
                "Every predicted physical root set is compared in the parent "
                "rational space with an independent materialized-child enumeration."
            ),
        },
        "totals": totals,
        "cases": cases,
        "nonzero_graph_glue_control": graph_glue_control,
        "proof_boundary": (
            "The theorem is exact in all ranks, but explicitly expanding the scaled "
            "shell may be more expensive than candidate-wise affine CVP in rank 15. "
            "These exhaustive controls certify correctness, not a uniform complexity "
            "bound or readiness of the 936 foundry rows, which still lack bridge/glue "
            "and compatible source markings."
        ),
        "software": {"sage": SAGE_VERSION},
        "reproduce": (
            "sage -python "
            "elkies-k3/scripts/certify_inverse_ade_projective_birth_strata.sage"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
