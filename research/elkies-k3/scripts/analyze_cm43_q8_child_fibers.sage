#!/usr/bin/env sage
"""Lift the marked q=8 child's E7+E8 roots to explicit CM-43 divisors.

This is an equation-construction aid.  It computes a deterministic simple
root system in the q=8 child frame, lifts every component through the exact
neighbor basis, and reports its class in
[F,O,E7(7),E8(8),P1,P2,P3].  Horizontal degrees are intersections with the
old Kumar fiber.  The output is not by itself an effectivity/chamber proof.
"""

from pathlib import Path
from collections import deque

from sage.all import *


SCRIPT_DIR = Path(__file__).resolve().parent
load(str(SCRIPT_DIR / "verify_kumar_cm43_q8_q9_factor.sage"))
load(str(SCRIPT_DIR / "verify_cm43_marked_divisor_transport.sage"))


def marked_to_geometric(marked_divisor):
    marked_divisor = vector(ZZ, marked_divisor)
    split_coordinates = vector(
        ZZ,
        list(marked_divisor[:2])
        + list(marked_divisor[2:]*frame_isometry.transpose()),
    )
    return vector(ZZ, split_coordinates*geometric_split)


root_result = pari(q8_child).qfminim(2)
half_roots = [
    vector(ZZ, column) for column in matrix(ZZ, root_result[2]).columns()
]
roots = half_roots+[-root for root in half_roots]
old_fiber_geometric = vector(ZZ, [1]+[0]*19)


def lifted_geometric(root):
    child_divisor = vector(ZZ, [0, 0]+list(root))
    return marked_to_geometric(vector(ZZ, child_divisor*q8_basis))


def effectivity_sign(root):
    """Orient roots by old horizontal degree, then old vertical chamber."""
    divisor = lifted_geometric(root)
    old_degree = divisor*geometric_ns*old_fiber_geometric
    if old_degree:
        return sign(old_degree)
    vertical_weight = sum(divisor[2:17])
    assert vertical_weight != 0
    return sign(vertical_weight)


positive_roots = [root for root in roots if effectivity_sign(root) > 0]
positive_set = {tuple(root) for root in positive_roots}
simple_roots = [
    root for root in positive_roots
    if not any(tuple(root-left) in positive_set for left in positive_roots)
]
assert len(simple_roots) == 15
simple = matrix(ZZ, [list(root) for root in simple_roots])
cartan = simple*q8_child*simple.transpose()

parents = list(range(15))


def find(index):
    while parents[index] != index:
        parents[index] = parents[parents[index]]
        index = parents[index]
    return index


def union(left, right):
    left = find(left)
    right = find(right)
    if left != right:
        parents[right] = left


for row in range(15):
    for column in range(row):
        if cartan[row, column]:
            union(row, column)
components = {}
for index in range(15):
    components.setdefault(find(index), []).append(index)
components = sorted(components.values(), key=len)
assert tuple(map(len, components)) == (7, 8)
assert tuple(
    cartan.matrix_from_rows_and_columns(indices, indices).det()
    for indices in components
) == (2, 1)

for component_number, indices in enumerate(components, 1):
    print(
        f"CM43Q8CHILD|component={component_number}|rank={len(indices)}"
        f"|det={cartan.matrix_from_rows_and_columns(indices, indices).det()}",
        flush=True,
    )

    component_positive_roots = []
    for root in positive_roots:
        coordinates = vector(
            QQ, simple.transpose().solve_right(root.column()).column(0)
        )
        if any(coordinates[index] for index in range(15) if index not in indices):
            continue
        assert all(coordinates[index] in ZZ and coordinates[index] >= 0
                   for index in indices)
        component_positive_roots.append((sum(coordinates), coordinates, root))
    _, highest_coordinates, _ = max(
        component_positive_roots, key=lambda item: item[0]
    )

    # Search the finite Weyl chambers for an affine simple system whose
    # representatives all have nonnegative degree over the old Kumar base.
    # Adding the new fiber changes old degree by four, so choose residues in
    # {0,1,2,3}.  The affine Kac-weighted degree sum must then be exactly four.
    if len(indices) == 7:
        initial_basis = matrix(
            ZZ, [list(simple_roots[index]) for index in indices]
        )
        component_cartan = cartan.matrix_from_rows_and_columns(indices, indices)
        multiplicities = vector(
            ZZ, [ZZ(highest_coordinates[index]) for index in indices]
        )

        def effective_old_curve(divisor):
            degree = divisor*geometric_ns*old_fiber_geometric
            if degree > 0:
                return True
            if degree < 0:
                return False
            if any(divisor[index] for index in (0, 1, 17, 18, 19)):
                return False
            vertical = list(divisor[2:17])
            return all(value >= 0 for value in vertical) and any(vertical)

        queue = deque([initial_basis])
        seen = {tuple(initial_basis.list())}
        effective_affine_basis = None
        while queue and len(seen) <= 300000:
            basis = queue.popleft()
            raw_divisors = [lifted_geometric(row) for row in basis.rows()]
            degrees = [
                ZZ(divisor*geometric_ns*old_fiber_geometric)
                for divisor in raw_divisors
            ]
            residues = [degree % 4 for degree in degrees]
            affine_residue = (-sum(
                multiplicities[index]*residues[index]
                for index in range(7)
            )) % 4
            weighted_degree = affine_residue+sum(
                multiplicities[index]*residues[index]
                for index in range(7)
            )
            if weighted_degree == 4:
                adjusted = [
                    raw_divisors[index]
                    + ((residues[index]-degrees[index])//4)*q8_geometric
                    for index in range(7)
                ]
                affine = q8_geometric-sum(
                    multiplicities[index]*adjusted[index]
                    for index in range(7)
                )
                if (
                    affine*geometric_ns*old_fiber_geometric == affine_residue
                    and all(effective_old_curve(divisor) for divisor in adjusted)
                    and effective_old_curve(affine)
                ):
                    effective_affine_basis = (adjusted, affine, residues)
                    break
            for reflection in range(7):
                transform = matrix.identity(ZZ, 7)
                transform[reflection, reflection] = -1
                for adjacent in range(7):
                    if component_cartan[adjacent, reflection] == -1:
                        transform[adjacent, reflection] = 1
                reflected = transform*basis
                key = tuple(reflected.list())
                if key not in seen:
                    seen.add(key)
                    queue.append(reflected)
        assert effective_affine_basis is not None
        adjusted, affine, residues = effective_affine_basis
        print(
            f"CM43Q8CHILD|effective_affine=E7|weyl_chambers={len(seen)}"
            f"|degrees={tuple(residues)+(affine_residue,)}",
            flush=True,
        )
        for index, divisor in enumerate(adjusted):
            print(
                f"CM43Q8CHILD|effective_E7_simple={index}"
                f"|divisor={tuple(divisor)}",
                flush=True,
            )
        print(
            f"CM43Q8CHILD|effective_E7_affine={tuple(affine)}",
            flush=True,
        )
    for simple_index in indices:
        geometric_divisor = lifted_geometric(simple_roots[simple_index])
        assert geometric_divisor*geometric_ns*geometric_divisor == -2
        old_degree = geometric_divisor*geometric_ns*old_fiber_geometric
        print(
            f"CM43Q8CHILD|simple={simple_index}|old_degree={old_degree}"
            f"|divisor={tuple(geometric_divisor)}",
            flush=True,
        )

    highest_divisor = sum(
        ZZ(highest_coordinates[index])*lifted_geometric(simple_roots[index])
        for index in indices
    )
    affine_divisor = q8_geometric-highest_divisor
    assert affine_divisor*geometric_ns*affine_divisor == -2
    assert all(
        affine_divisor*geometric_ns*lifted_geometric(simple_roots[index])
        in (0, 1)
        for index in indices
    )
    affine_degree = affine_divisor*geometric_ns*old_fiber_geometric
    print(
        f"CM43Q8CHILD|component={component_number}|affine"
        f"|multiplicities={tuple(ZZ(highest_coordinates[index]) for index in indices)}"
        f"|old_degree={affine_degree}|divisor={tuple(affine_divisor)}",
        flush=True,
    )

print("CM43Q8CHILD|status=PASS", flush=True)
