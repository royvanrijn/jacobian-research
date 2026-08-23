#!/usr/bin/env sage -python
"""
Certify the geometric affine-D13 component graph for the H3 q24 parent fibre
and match it to the deterministic D13 root basis used by the exact q24
Neron--Severi transport.

Status: ACTIVE_SEARCH, promoted to a graph/vertical-profile proof when all
assertions pass.

This script is the bridge between

  * the explicit 12-centre / 13-component resolution of the actual I9* germ,
  * the deterministic D13 simple-root basis used by
    close_h92_q8_q24_by_q6_translation.sage, and
  * the resolved Riemann--Roch module still to be constructed.

The crucial point is that blow-up centres are not components.  A rank-two
quadratic tangent cone contributes two geometric exceptional curves.  When a
later blow-up centre lies on an aggregate split label, all geometric branches
of that label are expanded before updating the dual graph.

The graph update is the standard one for a point blow-up on the resolved
configuration: remove pairwise incidences among all components through the
centre, then connect the new exceptional component to each of them.  For a
rank-two split tangent cone, first introduce the two exceptional branches and
their common incidence; a subsequent centre can then separate the resulting
multiple intersection.

Outputs:
  artifacts/local/elkies-k3/q24-i9star-component-graph-mod-<p>.json

The output contains:
  * the affine D13 graph including the old identity component F0;
  * the D13 root graph after removing F0;
  * the primitive Kodaira multiplicity vector from the affine intersection
    matrix kernel;
  * the two possible graph isomorphisms to the deterministic lattice simple
    roots (they differ only by the D13 spinor-arm swap);
  * the corresponding full vertical coefficients of -7F + V_D13 on every
    geometric component.

This does NOT yet impose the resolved marked-chord quotient.  If the two
spinor coefficients differ, the next module compiler must carry both
orientations until actual branch incidence/trivialization selects one.
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, gcd, lcm, matrix, vector


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
RESOLUTION = LOCAL / f"q24-i9star-resolution-mod-{args.prime}.json"
TRANSLATION = LOCAL / "q8-q24-physical-to-equation-translation.json"
CLOSE = ROOT / "elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"

for path in (RESOLUTION, TRANSLATION, CLOSE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")


def canonicalize_resolution_schema(resolution):
    centers = resolution.get("centers", [])
    labels = [str(record["label"]) for record in centers]
    old_e_labels = labels and all(
        label.startswith("E") and label[1:].isdigit()
        for label in labels
    )
    label_map = (
        {label: "C" + label[1:] for label in labels}
        if old_e_labels
        else {}
    )

    def relabel(value):
        value = str(value)
        for old, new in sorted(
            label_map.items(),
            key=lambda item: len(item[0]),
            reverse=True,
        ):
            value = value.replace(old, new)
        return value

    for record in centers:
        old_label = str(record["label"])
        record["label"] = label_map.get(old_label, old_label)
        record["active_components"] = [
            label_map.get(str(component), str(component))
            for component in record.get("active_components", [])
        ]
        if "path" in record:
            record["path"] = relabel(record["path"])

    split_centers = resolution.get(
        "split_tangent_cone_centers",
        resolution.get("split_center_labels", []),
    )
    split_centers = [label_map.get(str(label), str(label)) for label in split_centers]
    resolution["split_tangent_cone_centers"] = split_centers

    if "actual_geometric_exceptional_components" not in resolution:
        resolution["actual_geometric_exceptional_components"] = resolution.get(
            "actual_exceptional_irreducible_components"
        )

    for record in centers:
        tangent_cone = record.get("tangent_cone")
        if isinstance(tangent_cone, dict):
            continue
        factors = [
            {
                "factor": str(item.get("factor")),
                "multiplicity": int(item.get("multiplicity", 1)),
            }
            for item in record.get("tangent_factors", [])
        ]
        component_count = int(
            record.get(
                "new_exceptional_irreducible_components",
                2 if record["label"] in split_centers else 1,
            )
        )
        record["tangent_cone"] = {
            "polynomial": str(tangent_cone),
            "factorization_over_base_field": factors,
            "geometric_exceptional_components": component_count,
        }

    return resolution


resolution = json.loads(RESOLUTION.read_text())
resolution = canonicalize_resolution_schema(resolution)
translation = json.loads(TRANSLATION.read_text())
resolution_pass_statuses = {
    "PASS_EXPLICIT_MODP_I9STAR_D13_COMPONENT_RESOLUTION",
    "PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS",
}
assert resolution["status"] in resolution_pass_statuses
assert int(resolution["actual_blowup_centers"]) == 12
assert int(resolution["actual_geometric_exceptional_components"]) == 13
assert resolution["expected_root_lattice"] == "D13"
assert translation["status"] == "PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"
assert translation["q24_equation"]["child_root_lattice"] == "D12"

# Re-run the current exact transport in an embedded scope so the component
# graph is matched to exactly the same deterministic D13 root basis that
# produced the stored vertical_root_coefficients.  This is intentionally a
# replay, not a second independently chosen root ordering.
saved_argv = list(sys.argv)
scope = {"__name__": "__embedded_q24_close__"}
try:
    sys.argv = [str(CLOSE)]
    exec(compile(CLOSE.read_text(), str(CLOSE), "exec"), scope)
finally:
    sys.argv = saved_argv

for name in ("adapted", "vr", "vf"):
    if name not in scope:
        raise SystemExit(f"embedded q24 closeout did not expose {name}")

adapted = matrix(ZZ, scope["adapted"])
root_gram = matrix(ZZ, adapted[:13, :13])
vr = vector(ZZ, scope["vr"])
vf = ZZ(scope["vf"])
assert root_gram.nrows() == root_gram.ncols() == 13
assert len(vr) == 13
assert list(map(int, vr)) == list(
    map(int, translation["q24_equation"]["vertical_root_coefficients"])
)
assert int(vf) == int(translation["q24_equation"]["vertical_fibre_coefficient"])
assert vf == -7


def edge(left, right):
    if left == right:
        raise ValueError("loop in component graph")
    return tuple(sorted((str(left), str(right))))


def adjacency(vertices, edges):
    answer = {str(v): set() for v in vertices}
    for left, right in edges:
        answer[left].add(right)
        answer[right].add(left)
    return answer


def connected(vertices, edges):
    vertices = {str(v) for v in vertices}
    if not vertices:
        return True
    adj = adjacency(vertices, edges)
    start = next(iter(vertices))
    seen = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for nxt in adj[current]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return seen == vertices


# ---------------------------------------------------------------------------
# 1. Reconstruct the geometric dual graph from the actual blow-up chronology.
# ---------------------------------------------------------------------------
vertices = {"F0"}
edges = set()
center_components = {}
chronology = []

for record in resolution["centers"]:
    label = str(record["label"])
    component_count = int(
        record["tangent_cone"]["geometric_exceptional_components"]
    )
    if component_count not in (1, 2):
        raise ArithmeticError(f"unexpected component count at {label}")
    new_components = (
        [label]
        if component_count == 1
        else [label + "a", label + "b"]
    )

    active = []
    for active_label in record["active_components"]:
        active_label = str(active_label)
        if active_label == "F0":
            active.append("F0")
        else:
            if active_label not in center_components:
                raise ArithmeticError(
                    f"centre {label} refers to unknown active component {active_label}"
                )
            active.extend(center_components[active_label])
    if not active:
        raise ArithmeticError(f"centre {label} has no active old component")
    if len(set(active)) != len(active):
        raise ArithmeticError(f"centre {label} repeats an active geometric component")

    # All old branches listed at this point meet at the centre.  Their proper
    # transforms are separated by the blow-up, so remove every old incidence
    # among them before attaching the new exceptional locus.
    removed = []
    for left, right in itertools.combinations(active, 2):
        e = edge(left, right)
        if e in edges:
            edges.remove(e)
            removed.append(e)

    vertices.update(new_components)
    added = []
    for new_component in new_components:
        for old_component in active:
            e = edge(new_component, old_component)
            edges.add(e)
            added.append(e)

    # A rank-two tangent cone is the union of two lines.  Before any later
    # centre separates their intersection, those two exceptional branches
    # meet each other on the exceptional P^2.
    if len(new_components) == 2:
        e = edge(new_components[0], new_components[1])
        edges.add(e)
        added.append(e)

    center_components[label] = list(new_components)
    chronology.append({
        "center": label,
        "active_geometric_components": list(active),
        "new_geometric_components": list(new_components),
        "removed_edges": [list(e) for e in removed],
        "added_edges": [list(e) for e in added],
    })

assert len(vertices) == 14
assert len(vertices - {"F0"}) == 13
assert connected(vertices, edges)

root_vertices = set(vertices) - {"F0"}
root_edges = {
    e for e in edges if "F0" not in e
}
assert connected(root_vertices, root_edges)
assert len(root_edges) == 12  # tree on 13 vertices

root_adj = adjacency(root_vertices, root_edges)
root_degrees = {name: len(root_adj[name]) for name in root_vertices}
root_forks = sorted(name for name, degree in root_degrees.items() if degree == 3)
root_leaves = sorted(name for name, degree in root_degrees.items() if degree == 1)
assert len(root_forks) == 1
assert len(root_leaves) == 3
assert all(degree in (1, 2, 3) for degree in root_degrees.values())


def arm_data(adj):
    forks = [name for name, neighbors in adj.items() if len(neighbors) == 3]
    if len(forks) != 1:
        raise ArithmeticError(f"expected one trivalent D13 root, got {forks}")
    fork = forks[0]
    arms = []
    for first in sorted(adj[fork]):
        path = [first]
        previous = fork
        current = first
        while len(adj[current]) == 2:
            nxt = next(value for value in adj[current] if value != previous)
            previous, current = current, nxt
            path.append(current)
        if len(adj[current]) != 1:
            raise ArithmeticError("D13 arm did not terminate in a leaf")
        arms.append(path)
    arms.sort(key=lambda path: (len(path), tuple(map(str, path))))
    return fork, arms


geometric_fork, geometric_arms = arm_data(root_adj)
geometric_arm_lengths = sorted(len(path) for path in geometric_arms)
assert geometric_arm_lengths == [1, 1, 10]

print(
    "Q24D13GRAPH_GEOM|"
    f"vertices={len(root_vertices)}|edges={len(root_edges)}|"
    f"fork={geometric_fork}|leaves={','.join(root_leaves)}|"
    f"arms={','.join(map(str,geometric_arm_lengths))}|status=PASS_D13_GRAPH",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Recover Kodaira multiplicities from the affine intersection matrix.
# ---------------------------------------------------------------------------
affine_vertices = ["F0"] + sorted(root_vertices)
affine_index = {name: i for i, name in enumerate(affine_vertices)}
affine_matrix = matrix(ZZ, len(affine_vertices), len(affine_vertices))
for i in range(len(affine_vertices)):
    affine_matrix[i, i] = -2
for left, right in edges:
    i = affine_index[left]
    j = affine_index[right]
    affine_matrix[i, j] = affine_matrix[j, i] = 1

kernel = affine_matrix.right_kernel().basis_matrix()
assert kernel.nrows() == 1
raw = vector(QQ, kernel[0])
scale = ZZ(1)
for value in raw:
    scale = lcm(scale, ZZ(value.denominator()))
mult = vector(ZZ, [ZZ(scale * value) for value in raw])
content = gcd(tuple(mult))
if content:
    mult = vector(ZZ, [value // content for value in mult])
if any(value < 0 for value in mult):
    mult = -mult
assert all(value > 0 for value in mult)
assert affine_matrix * mult == vector(ZZ, [0] * len(mult))
assert mult[affine_index["F0"]] == 1
multiplicities = {
    name: int(mult[affine_index[name]]) for name in affine_vertices
}

print(
    "Q24D13GRAPH_FIBRE|"
    f"multiplicities={','.join(name+':'+str(multiplicities[name]) for name in affine_vertices)}|"
    "status=PASS_AFFINE_KERNEL",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Build the deterministic lattice D13 graph and match its arms.
# ---------------------------------------------------------------------------
for i in range(13):
    assert root_gram[i, i] == 2
    for j in range(13):
        if i == j:
            continue
        if root_gram[i, j] not in (0, -1):
            raise ArithmeticError(
                f"deterministic D13 simple-root Gram has unexpected entry {i},{j}={root_gram[i,j]}"
            )

lattice_vertices = {str(i) for i in range(13)}
lattice_edges = {
    edge(str(i), str(j))
    for i in range(13)
    for j in range(i + 1, 13)
    if root_gram[i, j] == -1
}
assert len(lattice_edges) == 12
assert connected(lattice_vertices, lattice_edges)
lattice_adj = adjacency(lattice_vertices, lattice_edges)
lattice_fork, lattice_arms = arm_data(lattice_adj)
lattice_arm_lengths = sorted(len(path) for path in lattice_arms)
assert lattice_arm_lengths == [1, 1, 10]

lat_short = [path for path in lattice_arms if len(path) == 1]
lat_long = next(path for path in lattice_arms if len(path) == 10)
geo_short = [path for path in geometric_arms if len(path) == 1]
geo_long = next(path for path in geometric_arms if len(path) == 10)
assert len(lat_short) == len(geo_short) == 2


def verify_mapping(mapping):
    if set(mapping) != lattice_vertices:
        return False
    if set(mapping.values()) != root_vertices:
        return False
    image_edges = {
        edge(mapping[left], mapping[right]) for left, right in lattice_edges
    }
    return image_edges == root_edges


mappings = []
for swap in (0, 1):
    mapping = {lattice_fork: geometric_fork}
    for left, right in zip(lat_long, geo_long):
        mapping[left] = right
    mapping[lat_short[0][0]] = geo_short[swap][0]
    mapping[lat_short[1][0]] = geo_short[1 - swap][0]
    assert verify_mapping(mapping)
    mappings.append(mapping)

# The two maps must differ only by the two short spinor arms.
assert all(
    mappings[0][vertex] == mappings[1][vertex]
    for vertex in lattice_vertices
    if vertex not in {lat_short[0][0], lat_short[1][0]}
)

# ---------------------------------------------------------------------------
# 4. Push the exact -7F + V_D13 coefficients to geometric components.
# ---------------------------------------------------------------------------
profiles = []
for orientation, mapping in enumerate(mappings):
    root_coefficients = {
        mapping[str(i)]: int(vr[i]) for i in range(13)
    }
    full_vertical = {"F0": int(vf * multiplicities["F0"])}
    for name in sorted(root_vertices):
        full_vertical[name] = int(
            vf * multiplicities[name] + ZZ(root_coefficients[name])
        )
    profiles.append({
        "orientation": orientation,
        "lattice_to_geometric": {
            key: mapping[key] for key in sorted(mapping, key=int)
        },
        "geometric_root_coefficients": root_coefficients,
        "full_vertical_coefficients": full_vertical,
    })

short_lattice_indices = [int(lat_short[0][0]), int(lat_short[1][0])]
short_coefficients = [int(vr[index]) for index in short_lattice_indices]
spinor_swap_irrelevant = short_coefficients[0] == short_coefficients[1]

split_centers = list(resolution.get("split_tangent_cone_centers", []))
assert len(split_centers) == 1
split_center = split_centers[0]
split_record = next(
    item for item in resolution["centers"] if item["label"] == split_center
)
base_factorization = split_record["tangent_cone"]["factorization_over_base_field"]
split_over_base_field = len(base_factorization) == 2

print(
    "Q24D13GRAPH_MATCH|"
    f"lattice_fork={lattice_fork}|geometric_fork={geometric_fork}|"
    f"spinor_indices={','.join(map(str,short_lattice_indices))}|"
    f"spinor_coefficients={','.join(map(str,short_coefficients))}|"
    f"swap_irrelevant={int(spinor_swap_irrelevant)}|"
    f"split_center={split_center}|split_over_Fp={int(split_over_base_field)}|"
    "status=PASS_TWO_D13_ORIENTATIONS",
    flush=True,
)
for profile in profiles:
    print(
        "Q24D13GRAPH_VERTICAL|"
        f"orientation={profile['orientation']}|"
        + "|".join(
            f"{name}={profile['full_vertical_coefficients'][name]}"
            for name in affine_vertices
        )
        + "|status=PASS",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h3-q24-i9star-component-graph-modp.v1",
    "status": "PASS_H3_Q24_AFFINE_D13_COMPONENT_GRAPH",
    "prime": int(args.prime),
    "resolution_input": str(RESOLUTION.relative_to(ROOT)),
    "translation_input": str(TRANSLATION.relative_to(ROOT)),
    "geometric_graph": {
        "affine_vertices": affine_vertices,
        "affine_edges": [list(value) for value in sorted(edges)],
        "root_vertices": sorted(root_vertices),
        "root_edges": [list(value) for value in sorted(root_edges)],
        "root_fork": geometric_fork,
        "root_arm_lengths": geometric_arm_lengths,
        "chronology": chronology,
    },
    "fibre_multiplicities": multiplicities,
    "lattice_graph": {
        "root_gram": [[int(value) for value in row] for row in root_gram.rows()],
        "edges": [list(value) for value in sorted(lattice_edges)],
        "fork": lattice_fork,
        "arm_lengths": lattice_arm_lengths,
        "vertical_root_coefficients": list(map(int, vr)),
        "vertical_fibre_coefficient": int(vf),
    },
    "spinor_orientation": {
        "short_lattice_indices": short_lattice_indices,
        "short_root_coefficients": short_coefficients,
        "swap_irrelevant": bool(spinor_swap_irrelevant),
        "split_center": split_center,
        "split_tangent_factorization_over_base_field": base_factorization,
        "split_over_base_field": bool(split_over_base_field),
        "candidate_profiles": profiles,
    },
    "boundary": (
        "The affine D13 graph, fibre multiplicities, and the two graph-theoretic "
        "matches to the deterministic lattice root basis are certified.  If the "
        "two spinor coefficients differ, graph data alone does not orient the "
        "two short arms.  The resolved marked-chord compiler must carry both "
        "profiles until actual branch/trivialization incidence selects one.  "
        "No ordinary-jet quotient is promoted here."
    ),
    "next": (
        "Pull the ten post-collision marked-chord basis functions through the "
        "actual resolution charts, apply each candidate full vertical profile, "
        "derive the connected resolved quotient, and require rank 8 / h0=2 / "
        "squarefree quartic degree 4 before certifying the D12 child."
    ),
}

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-i9star-component-graph-mod-{args.prime}.json"
)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24D13GRAPH_RESULT|"
    f"root_vertices=13|root_edges=12|arms=1,1,10|"
    f"orientations=2|swap_irrelevant={int(spinor_swap_irrelevant)}|"
    "status=PASS_H3_Q24_AFFINE_D13_COMPONENT_GRAPH",
    flush=True,
)
