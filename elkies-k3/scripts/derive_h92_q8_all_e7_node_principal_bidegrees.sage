#!/usr/bin/env sage -python
"""Derive actual two-parameter q8 principal-bidegree templates at all E7 nodes.

Each actual H92 E7 edge node now has q6 frame t*R, so the q8 local condition
is ``g*f/t^9 in R``.  This script records the leading bidegree in the two
completed regular parameters for every endpoint ambient monomial
``u^i*x^a*m^b/h^k``:

    ord(g) + (4k-i-9) ord(t) + a ord(x) + b ord(m).

The vectors come from actual chart pullbacks and implicit surface equations;
they are not Kodaira-label data.  Negative bidegrees are necessary local
principal-part conditions.  Equal or incomparable leading bidegrees still
require the subsequent exact two-variable quotient/overlap calculation.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
UNMARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-node-frames.json"
MARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-2-5-node-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-e7-node-principal-bidegrees.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--unmarked", type=Path, default=UNMARKED)
parser.add_argument("--marked", type=Path, default=MARKED)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
unmarked = json.loads(args.unmarked.read_text())
marked = json.loads(args.marked.read_text())
assert ambient["status"] in (
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT", "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
)
assert unmarked["status"] == "PASS_EXACT_Q8_UNMARKED_E7_NODE_FRAMES"
assert marked["status"] == "PASS_EXACT_Q8_E7_2_5_NODE_FRAME"

# Coordinates are in the actual completed regular parameters listed below.
# For example, E7_4--E7_3 has U=Y^2*unit and parameters (Z,Y), giving
# t=(3,4), x=(4,6), m=(2,3), g=(4,6).  Every row is independently tied to
# an emitted actual node frame and its Cartier factor.
frame_data = {
    "E7_1--E7_4": {"parameters": ("U", "Y"), "t": (3, 2), "x": (4, 2), "m": (2, 1), "g": (4, 2), "g_text": "U^4*Y^2"},
    "E7_4--E7_3": {"parameters": ("Z", "Y"), "t": (3, 4), "x": (4, 6), "m": (2, 3), "g": (4, 6), "g_text": "Z^4*Y^6"},
    "E7_3--E7_7": {"parameters": ("U", "Y"), "t": (3, 4), "x": (5, 6), "m": (2, 3), "g": (5, 6), "g_text": "U^5*Y^6"},
    "E7_7--E7_2": {"parameters": ("Z", "Y"), "t": (3, 2), "x": (5, 4), "m": (2, 1), "g": (5, 5), "g_text": "Z^5*Y^5"},
    "E7_3--E7_6": {"parameters": ("Z", "Y"), "t": (2, 4), "x": (3, 6), "m": (2, 3), "g": (3, 6), "g_text": "Z^3*Y^6"},
    "E7_2--E7_5": {"parameters": ("Z", "Y"), "t": (1, 2), "x": (2, 4), "m": (0, 1), "g": (6, 5), "g_text": "Z^6*Y^5"},
}
unmarked_by_name = {entry["chart"]: entry for entry in unmarked["unmarked_nodes"]}
assert set(unmarked_by_name) == set(frame_data)-{"E7_2--E7_5"}
for name, record in unmarked_by_name.items():
    assert record["q6_module"] == "t*R"
    assert record["q8_cartier_factor"] == frame_data[name]["g_text"]
assert marked["chart"]["name"] == "E7_2--E7_5"
assert marked["module_frame"]["q6_fractional_module"].endswith("=t*(1,m)=t*R")
assert marked["module_frame"]["q8_twist_cartier_equation"] == "g=Z^6*Y^5"


def bidegree(entry, frame):
    delta = 4*int(entry["h_power"])-int(entry["u_power"])-9
    return tuple(
        int(frame["g"][index] + delta*frame["t"][index]
            + int(entry["x_power"])*frame["x"][index]
            + int(entry["m_power"])*frame["m"][index])
        for index in range(2)
    )


def pareto_minimal(groups):
    """Return bidegrees no other negative group can precede componentwise."""
    degrees = tuple(groups)
    return tuple(
        degree for degree in degrees
        if not any(
            other != degree and other[0] <= degree[0] and other[1] <= degree[1]
            for other in degrees
        )
    )


node_records = []
for name, frame in frame_data.items():
    groups = {}
    for index, entry in enumerate(ambient["ambient_basis"]):
        degree = bidegree(entry, frame)
        if min(degree) < 0:
            groups.setdefault(degree, []).append(index)
    groups_payload = [
        {"bidegree": list(degree), "basis_indices": labels,
         "singleton": len(labels) == 1}
        for degree, labels in sorted(groups.items())
    ]
    singleton_indices = [
        item["basis_indices"][0] for item in groups_payload if item["singleton"]
    ]
    minima = set(pareto_minimal(groups))
    minimal_payload = [item for item in groups_payload if tuple(item["bidegree"]) in minima]
    unique_minimal_indices = [
        item["basis_indices"][0] for item in minimal_payload if item["singleton"]
    ]
    node_records.append({
        "chart": name,
        "parameters": list(frame["parameters"]),
        "leading_orders": {key: list(frame[key]) for key in ("t", "x", "m", "g")},
        "formula": "g+(4*k-i-9)*t+a*x+b*m",
        "negative_bidegree_groups": groups_payload,
        "negative_group_count": len(groups_payload),
        "unique_negative_leading_terms": singleton_indices,
        "unique_negative_leading_count": len(singleton_indices),
        "pareto_minimal_negative_groups": minimal_payload,
        "pareto_minimal_negative_count": len(minimal_payload),
        "exact_unique_pareto_leading_terms": unique_minimal_indices,
        "exact_unique_pareto_leading_count": len(unique_minimal_indices),
    })

payload = {
    "schema": "elkies-k3.h92-q8-all-e7-node-principal-bidegrees.v1",
    "status": "PASS_EXACT_Q8_ALL_E7_NODE_PRINCIPAL_BIDEGREE_TEMPLATE",
    "inputs": {
        "endpoint_ambient": {"path": path_label(args.ambient), "sha256": digest(args.ambient)},
        "unmarked_node_frames": {"path": path_label(args.unmarked), "sha256": digest(args.unmarked)},
        "cancellation_node_frame": {"path": path_label(args.marked), "sha256": digest(args.marked)},
    },
    "ambient_dimension": int(ambient["ambient_dimension"]),
    "node_conditions": node_records,
    "exact_initial_constraint_rule": (
        "A singleton Pareto-minimal negative bidegree has a nonzero unit "
        "leading coefficient and cannot receive a contribution from another "
        "term after multiplication by a regular unit. Its ambient coefficient "
        "therefore vanishes as an exact necessary node condition."
    ),
    "boundary": (
        "These are exact leading two-parameter node templates. They do not "
        "turn equal/incomparable groups into sufficient quotient conditions, "
        "handle the marked smooth point, prove overlap compatibility, or "
        "produce a complete q8 matrix, kernel, or pencil."
    ),
}
exact_initial_indices = sorted({
    index for item in node_records for index in item["exact_unique_pareto_leading_terms"]
})
payload["exact_initial_unique_basis_indices"] = exact_initial_indices
payload["exact_initial_coordinate_rank"] = len(exact_initial_indices)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8ALLNODEBIDEGREES|nodes=6|ambient={}|negative_groups={}|singletons={}|pareto={}|exact_rows={}|"
    "status=PASS_EXACT_Q8_ALL_E7_NODE_PRINCIPAL_BIDEGREE_TEMPLATE".format(
        payload["ambient_dimension"],
        sum(item["negative_group_count"] for item in node_records),
        sum(item["unique_negative_leading_count"] for item in node_records),
        sum(item["pareto_minimal_negative_count"] for item in node_records),
        len(exact_initial_indices),
    ),
    flush=True,
)
