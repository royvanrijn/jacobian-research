#!/usr/bin/env sage
"""Find the CM24-only q4 lattice marking of the pair14 fifth pencil.

The unit-preserving pair14 Jacobian has geometric root data (15,82,360),
namely A5+A5+A4+A1.  It does not occur among the generic fourth-frame q4
neighbors with horizontal class +/-(1,0), so this search works directly in
the rank-18 CM24 fourth frame.  It scans a bounded indexed window of its
norm-eight shell, applies one q4 factorization, and retains exact children
with the pair14 root data and horizontal marking.

This is a finite lattice-marking calculation; the equation-level fiber
signature is certified separately by audit_q80_fifth_q4_pair14_twist_gf73.sage.
The default target is pair14.  ``--target-root-data`` permits the same bounded
marking search for another equation-level CM24 fiber signature.
"""

import argparse
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, lcm, matrix, pari, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--start-half", type=int, default=0)
parser.add_argument("--stop-half", type=int, default=None)
parser.add_argument("--a", type=int, default=2)
parser.add_argument("--b", type=int, default=2)
parser.add_argument("--mw-coordinates", default="1,0")
parser.add_argument("--mw-up-to-sign", action="store_true")
parser.add_argument(
    "--target-root-data",
    default="15,82,360",
    help="comma-separated root rank,count,determinant (default: pair14)",
)
parser.add_argument("--stop-after", type=int, default=1)
parser.add_argument("--output", type=Path, default=None)
search_arguments = parser.parse_args()

a = ZZ(search_arguments.a)
b = ZZ(search_arguments.b)
if a*b != 4:
    parser.error("the norm-eight q4 shell requires --a * --b = 4")
target_coordinates = tuple(
    QQ(value) for value in search_arguments.mw_coordinates.split(",")
)
if len(target_coordinates) != 2:
    parser.error("--mw-coordinates requires two comma-separated values")
target_roots = tuple(
    ZZ(value) for value in search_arguments.target_root_data.split(",")
)
if len(target_roots) != 3:
    parser.error("--target-root-data requires rank,count,determinant")

# Hide this wrapper's options from the readiness checker's loader chain.
saved_argv = list(sys.argv)
try:
    sys.argv = [sys.argv[0]]
    load(str(HERE / "analyze_q80_fifth_q4_cm24_readiness.sage"))
finally:
    sys.argv = saved_argv

source = special_fourth
coordinate_map = (
    source
    * new_projected_basis.transpose()
    * new_optimal_height.inverse()
)

# Do not enumerate the full rank-18 norm-eight shell: it exceeds PARI's 1GB
# stack.  First impose the two exact horizontal-coordinate equations.  Their
# saturated integer kernel has rank 16; completing the square then asks for a
# small sphere in one affine root/glue coset.
denominator = lcm(value.denominator() for value in coordinate_map.list())
constraints = matrix(ZZ, denominator*coordinate_map.transpose())
kernel = constraints.right_kernel_matrix()
assert kernel.nrows() == kernel.rank() == 16
assert kernel*coordinate_map == 0
# The saturated kernel returned by echelon reduction has very large entries.
# Reduce it in the source metric before invoking any closest-vector routine;
# otherwise the generic enumerator can consume a gigabyte on this tiny-radius
# problem despite the underlying ADE decomposition.
kernel = kernel.LLL(gram=source)
assert kernel.nrows() == kernel.rank() == 16
assert kernel*coordinate_map == 0

particular = None
particular_pair = None
for left in range(len(candidate_rows)):
    for right in range(left, len(candidate_rows)):
        candidate = vector(
            ZZ,
            vector(ZZ, candidate_rows[left][2][2:])
            + vector(ZZ, candidate_rows[right][2][2:]),
        )
        if tuple(candidate*coordinate_map) == target_coordinates:
            particular = candidate
            particular_pair = (left, right)
            break
    if particular is not None:
        break
assert particular is not None

center = vector(QQ, target_coordinates)*new_projected_basis
assert new_project_mw(particular) == center
center_norm = center*source*center
residual_norm = QQ(8)-center_norm
assert residual_norm > 0

# The kernel is the order-three primitive closure of the displayed
# A5+A3+4A2 root lattice.  Split it into its three root-lattice cosets and
# enumerate each ADE component separately.  The largest enumeration then has
# rank five, instead of a memory-hungry rank-16 generic CVP.
root_module = new_simple.row_module()
kernel_module = kernel.row_module()
assert all(row in kernel_module for row in new_simple.rows())
glue = next(vector(ZZ, row) for row in kernel.rows() if row not in root_module)
assert glue not in root_module and 2*glue not in root_module
assert 3*glue in root_module
assert all(
    any(vector(ZZ, row)-multiple*glue in root_module for multiple in range(3))
    for row in kernel.rows()
)

root_components = []
unseen = set(range(new_root_gram.nrows()))
while unseen:
    seed = min(unseen)
    component = {seed}
    frontier = [seed]
    unseen.remove(seed)
    while frontier:
        current = frontier.pop()
        neighbors = {
            index for index in tuple(unseen)
            if new_root_gram[current, index] != 0
        }
        component.update(neighbors)
        unseen.difference_update(neighbors)
        frontier.extend(neighbors)
    root_components.append(tuple(sorted(component)))
root_components = tuple(sorted(root_components, key=lambda row: (-len(row), row)))
assert tuple(sorted(map(len, root_components), reverse=True)) == (5, 3, 2, 2, 2, 2)


def component_close_vectors(gram, offset):
    lattice = IntegralLattice(gram)
    iterator = lattice.enumerate_close_vectors(-offset)
    rows = []
    for _index in range(100_000):
        coefficients = vector(ZZ, next(iterator))
        distance = (coefficients+offset)*gram*(coefficients+offset)
        if distance > residual_norm:
            break
        rows.append((distance, coefficients))
    else:
        raise RuntimeError("component close-vector enumeration exceeded hard cap")
    return tuple(rows)


positive_candidates = []
distance_counts = {}
for glue_multiple in range(3):
    offset_ambient = particular+glue_multiple*glue-center
    offset_coordinates = (
        offset_ambient
        * source
        * new_simple.transpose()
        * new_root_gram.inverse()
    )
    assert offset_coordinates*new_simple == offset_ambient
    component_choices = []
    for component in root_components:
        gram = new_root_gram.matrix_from_rows_and_columns(component, component)
        offset = vector(QQ, [offset_coordinates[index] for index in component])
        component_choices.append(component_close_vectors(gram, offset))

    def assemble(component_index, distance, coefficients):
        if distance > residual_norm:
            return
        if component_index == len(root_components):
            distance_counts[str(distance)] = distance_counts.get(str(distance), 0)+1
            if distance != residual_norm:
                return
            root_coefficients = vector(ZZ, [0]*new_simple.nrows())
            for component, values in zip(root_components, coefficients):
                for index, value in zip(component, values):
                    root_coefficients[index] = value
            witness = (
                particular+glue_multiple*glue+root_coefficients*new_simple
            )
            assert witness*source*witness == 8
            assert tuple(witness*coordinate_map) == target_coordinates
            positive_candidates.append(witness)
            return
        for component_distance, values in component_choices[component_index]:
            assemble(
                component_index+1,
                distance+component_distance,
                coefficients+(values,),
            )

    assemble(0, QQ(0), ())

affine_candidates = list(positive_candidates)
if search_arguments.mw_up_to_sign:
    affine_candidates.extend(-witness for witness in positive_candidates)
affine_candidates = tuple(
    vector(ZZ, row) for row in sorted(set(map(tuple, affine_candidates)))
)
stop_half = search_arguments.stop_half or len(affine_candidates)
window = affine_candidates[search_arguments.start_half:stop_half]
print(
    "Q80FIFTHCM24MARK|stage=affine_cvp|"
    f"kernel_rank={kernel.rank()}|particular_pair={particular_pair}|"
    f"center_norm={center_norm}|residual_norm={residual_norm}|"
    f"positive_candidates={len(positive_candidates)}|"
    f"oriented_candidates={len(affine_candidates)}|"
    f"window={search_arguments.start_half}:{stop_half}|a={a}|b={b}",
    flush=True,
)

candidates = []
oriented_horizontal = 0
primitive = 0
root_count_matches = 0
for witness in window:
    for _single_orientation in (0,):
        coordinates = tuple(witness*coordinate_map)
        coordinate_match = coordinates == target_coordinates
        if search_arguments.mw_up_to_sign:
            coordinate_match = coordinate_match or coordinates == tuple(
                -value for value in target_coordinates
            )
        if not coordinate_match:
            continue
        oriented_horizontal += 1
        try:
            child, transport = neighbor(source, ZZ(4), a, b, witness)
        except AssertionError:
            continue
        primitive += 1
        minimum = pari(child).qfminim(2)
        root_count = int(minimum[0])
        if root_count != target_roots[1]:
            continue
        root_count_matches += 1
        roots = matrix(ZZ, minimum[2]).transpose()
        root_basis = roots.row_module().basis_matrix()
        root_gram = root_basis*child*root_basis.transpose()
        root_data = (
            int(root_basis.rank()),
            root_count,
            int(abs(root_gram.det())),
        )
        if root_data != target_roots:
            continue

        raw_fiber = vector(ZZ, [a, b] + list(witness))
        reduced_fiber, reflections = chamber_reduce(raw_fiber, new_curves, new_ns)
        assert reduced_fiber*new_ns*reduced_fiber == 0
        assert all(
            intersection(reduced_fiber, curve, new_ns) >= 0
            for _name, curve in new_curves
        )
        old_pullback = vector(ZZ, reduced_fiber*special_fourth_basis)
        row = {
            "a": int(a),
            "b": int(b),
            "v": list(map(int, witness)),
            "horizontal_coordinates": [str(value) for value in coordinates],
            "child_root_data": list(root_data),
            "child_frame": [list(map(int, item)) for item in child.rows()],
            "transport": [list(map(int, item)) for item in transport.rows()],
            "reflection_count": len(reflections),
            "reduced_fiber": list(map(int, reduced_fiber)),
            "reduced_D_dot_old_F": int(
                intersection(reduced_fiber, new_fiber, new_ns)
            ),
            "reduced_D_dot_old_O": int(
                intersection(reduced_fiber, new_zero, new_ns)
            ),
            "old_pullback": list(map(int, old_pullback)),
            "old_fiber_degree": int(
                intersection(old_pullback, old_fiber, old_special_ns)
            ),
            "old_zero_pairing": int(
                intersection(old_pullback, old_zero, old_special_ns)
            ),
        }
        candidates.append(row)
        print(
            "Q80FIFTHCM24MARK|stage=hit|"
            f"hit={len(candidates)}|coordinates={coordinates}|"
            f"root_data={root_data}|D.oldF={row['reduced_D_dot_old_F']}|"
            f"D.oldO={row['reduced_D_dot_old_O']}|"
            f"v={tuple(map(int, witness))}",
            flush=True,
        )
        if search_arguments.stop_after and len(candidates) >= search_arguments.stop_after:
            break
    if search_arguments.stop_after and len(candidates) >= search_arguments.stop_after:
        break

output = search_arguments.output or (
    ROOT / "artifacts/local/q80-fifth-q4-cm24-pair14-marking.json"
)
output.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": (
        "q80-fifth-q4-cm24-pair14-lattice-marking-v1"
        if target_roots == (15, 82, 360)
        else "q80-fifth-q4-cm24-target-lattice-marking-v1"
    ),
    "status": "bounded_indexed_norm_eight_shell_experiment",
    "prime": 73,
    "neighbor_factorization": [int(a), int(b)],
    "target_horizontal_coordinates": [str(value) for value in target_coordinates],
    "horizontal_up_to_sign": bool(search_arguments.mw_up_to_sign),
    "target_root_data": list(target_roots),
    "affine_kernel_rank": int(kernel.rank()),
    "affine_particular_pair": list(particular_pair),
    "affine_center_norm": str(center_norm),
    "affine_residual_norm": str(residual_norm),
    "affine_first_distance_counts": distance_counts,
    "positive_affine_candidates": len(positive_candidates),
    "oriented_affine_candidates": len(affine_candidates),
    "half_shell_window": [search_arguments.start_half, stop_half],
    "oriented_horizontal": oriented_horizontal,
    "primitive": primitive,
    "root_count_matches": root_count_matches,
    "candidates": candidates,
    "rank_claim": None,
    "reproduce": (
        "sage elkies-k3/scripts/"
        "search_q80_fifth_q4_cm24_pair14_marking.sage "
        f"--start-half {search_arguments.start_half} --stop-half {stop_half} "
        f"--a {a} --b {b} --mw-coordinates "
        f"{','.join(map(str, target_coordinates))}"
        f" --target-root-data {','.join(map(str, target_roots))}"
        + (" --mw-up-to-sign" if search_arguments.mw_up_to_sign else "")
        + f" --stop-after {search_arguments.stop_after}"
    ),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n")
print(
    "Q80FIFTHCM24MARK|"
    f"window={search_arguments.start_half}:{stop_half}|"
    f"oriented_horizontal={oriented_horizontal}|primitive={primitive}|"
    f"root_count_matches={root_count_matches}|hits={len(candidates)}|"
    f"output={output}|status=PASS_BOUNDED_CM24_SHELL",
    flush=True,
)
