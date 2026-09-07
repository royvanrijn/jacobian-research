#!/usr/bin/env sage
"""Regression checks for the reusable exact neighbor engine and certificates."""

from pathlib import Path

from sage.all import CartanMatrix, QQ, ZZ, identity_matrix, matrix, vector


HERE = Path(__file__).resolve().parent
load(str(HERE / "exact_neighbor_examples.sage"))


for factory, expected_fixed_components, expected_mw_rank in (
    (q80_first_q4_example, 17, 4),
    (h3_d13_q24_example, 0, 5),
):
    example = factory()
    result, certificate = run_example(example)
    assert result["nonnegative_on_supplied_curves"]
    assert result["child_root_lattice_primitive"]
    assert result["minimization_status"] == "ROOT_MW_MINIMIZED"
    assert len(result["fixed_component_sequence"]) == expected_fixed_components
    assert result["child_mw_height"].nrows() == expected_mw_rank
    assert certificate["schema"] == CERTIFICATE_SCHEMA
    assert certificate["child"]["root_data"] == list(example["expected_root_data"])
    assert certificate["certificate_sha256"] == json_digest({
        key: value for key, value in certificate.items()
        if key != "certificate_sha256"
    })
    assert len(certificate["child"]["lifted_simple_components"]) == example["expected_root_data"][0]
    serialized_input = neighbor_input(
        example["ns"], example["old_fiber"], example["divisor"], example["curves"],
        example["proof_metadata"],
    )
    replay_result, replay_certificate = run_neighbor_input(serialized_input)
    assert replay_certificate == certificate
    first_root = result["child_simple_roots"].row(0)
    lifted_root = lift_child_frame_vector(result, first_root)
    assert transport_parent_vector_to_child(result, lifted_root) == vector(
        ZZ, [0, 0] + list(first_root)
    )
    marked = transport_marked_parent_vectors(
        result, {"old_fiber": example["old_fiber"]}
    )
    assert len(marked["old_fiber"]) == 19
    assert replay_result["child_frame"] == result["child_frame"]
    print(
        "EXACTNEIGHBOR|example={}|degree=2|fixed_components={}|"
        "child_root_data={}|mw_rank={}|child_det={}|status=PASS".format(
            example["name"],
            len(result["fixed_component_sequence"]),
            result["child_root_data"],
            result["child_mw_height"].nrows(),
            result["child_frame"].det(),
        ),
        flush=True,
    )

# Keep the terminal rootless branch and the nonprimitive-root fallback covered.
rootless_fixture = matrix(ZZ, ((4, 0), (0, 4)))
rootless_result = minimize_child_frame(rootless_fixture)
assert rootless_result["root_data"] == (0, 0, 1)
assert rootless_result["minimization_status"] == "ROOTLESS_LLL_MINIMIZED"
assert (
    rootless_result["basis"] * rootless_fixture
    * rootless_result["basis"].transpose() == rootless_result["frame"]
)

# The even unimodular overlattice D16+ has root system D16 of index two.  The
# engine retains its exact child and reports that torsion/glue data are needed.
d16_simple = matrix(ZZ, [
    [1 if column == row else -1 if column == row + 1 else 0
     for column in range(16)]
    for row in range(15)
] + [[0] * 14 + [1, 1]])
d16_glue = vector(QQ, [QQ(1) / 2] * 16) * d16_simple.inverse()
nonprimitive_basis = matrix(
    QQ, list(identity_matrix(ZZ, 16).rows()) + [d16_glue]
).row_module(ZZ).basis_matrix()
nonprimitive_fixture = (
    nonprimitive_basis * CartanMatrix(["D", 16])
    * nonprimitive_basis.transpose()
).change_ring(ZZ)
nonprimitive_result = minimize_child_frame(nonprimitive_fixture)
assert nonprimitive_result["root_data"] == (16, 480, 4)
assert not nonprimitive_result["root_lattice_primitive"]
assert nonprimitive_result["root_smith_invariants"][-1] == 2
assert nonprimitive_result["mw_height"] is None
assert nonprimitive_result["minimization_status"] == "PARTIAL_NONPRIMITIVE_ROOT_LATTICE"
print("EXACTNEIGHBOR|fixtures=rootless,nonprimitive|status=PASS", flush=True)

# The U-split and parent/child transport are degree-independent.  This
# primitive degree-three isotropic class in U plus negative A2 is the compact
# regression for a future q=3 equation-level RR compiler.
degree_three_ns = matrix(ZZ, [
    [0, 1, 0, 0], [1, 0, 0, 0],
    [0, 0, -2, 1], [0, 0, 1, -2],
])
degree_three_old_fiber = vector(ZZ, [1, 0, 0, 0])
degree_three_divisor = vector(ZZ, [1, 3, 2, 1])
assert intersection(degree_three_divisor, degree_three_divisor, degree_three_ns) == 0
degree_three_result = degree_q_neighbor(
    degree_three_ns, degree_three_divisor, degree_three_old_fiber,
    expected_old_fiber_degree=3,
)
assert degree_three_result["old_fiber_degree"] == 3
assert degree_three_result["reduced_divisor"] == degree_three_divisor
assert degree_three_result["transport"] * degree_three_ns * degree_three_result["transport"].transpose() == matrix(ZZ, [
    [0, 1, 0, 0], [1, 0, 0, 0],
    [0, 0, -degree_three_result["child_frame"][0, 0], -degree_three_result["child_frame"][0, 1]],
    [0, 0, -degree_three_result["child_frame"][1, 0], -degree_three_result["child_frame"][1, 1]],
])
serialized_degree_three = neighbor_input(
    degree_three_ns, degree_three_old_fiber, degree_three_divisor, (),
    expected_old_fiber_degree=3,
)
assert serialized_degree_three["expected_old_fiber_degree"] == 3
assert run_neighbor_input(serialized_degree_three)[0]["old_fiber_degree"] == 3
try:
    degree_q_neighbor(
        degree_three_ns, degree_three_divisor, degree_three_old_fiber,
        expected_old_fiber_degree=2,
    )
except ValueError as error:
    assert "expected" in str(error)
else:
    raise AssertionError("wrong declared old-fiber degree was accepted")
print("EXACTNEIGHBOR|fixture=primitive-degree-three|degree=3|status=PASS", flush=True)
