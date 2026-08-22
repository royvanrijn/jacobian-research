#!/usr/bin/env sage
"""Regression checks for the strict resolved-chart neighbour compiler core."""

import json

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector

CORE = "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
exec(compile(open(CORE).read(), CORE, "exec"))

ambient = ("1", "z", "z2", "z3")
first = quotient_condition(
    "chart_A", ambient,
    lambda value: {
        "1": (1, 0), "z": (0, 1), "z2": (0, 0), "z3": (0, 0),
    }[value],
    ("1", "z"), "synthetic resolved chart",
)
compiled = compile_resolved_conditions(ambient, (first,), complete=True)
assert compiled["ambient_dimension"] == 4
assert compiled["condition_rows"] == 2
assert compiled["rank"] == 2
assert compiled["kernel_dimension"] == 2
assert compiled["kernel_materialization"] == "zero_columns"
assert compiled["kernel_basis"] == matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]])
assert require_h0_two(compiled).nrows() == 2

# A finite ambient image can be recorded without falsely asserting that its
# target local quotient is finite.  This is the modular normal-form hand-off
# used by the high-degree node compiler.
finite_image = finite_ambient_image_condition(
    "synthetic infinite quotient image",
    ("a", "b", "c"),
    lambda value: {
        "a": {(0, 0): 1, (1, 0): 2},
        "b": {(1, 0): 3},
        "c": {},
    }[value],
    lambda key: key,
    GF(7),
    "synthetic local normal forms",
)
assert finite_image["coordinate_keys"] == ((0, 0), (1, 0))
assert finite_image["matrix"] == matrix(GF(7), [[1, 0, 0], [2, 3, 0]])
assert finite_image["matrix"].rank() == 2

# Transition identities make actual chart gluing auditable before a finite
# overlap quotient is chosen.  The frame ratio need not be a unit globally.
transition_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
tz, tu, ty = transition_ring.gens()
transition = verify_resolved_chart_transition(
    "synthetic sibling blow-up charts",
    transition_ring,
    ty**2-tu*tz**2,
    ty**2-tz*tu,
    (tu*tz, 1/tz, ty/tz),
    {"old_t": tz*tu},
    {"old_t": tz},
    tz*tu,
    tz,
    "synthetic U/Z chart overlap",
)
assert transition["surface_ratio"] == 1/tz**2
assert transition["transported_pullbacks"]["old_t"] == tz*tu
assert transition["frame_ratio"] == 1

# An overlap condition compares only caller-supplied resolved-chart maps in a
# common finite quotient.  No component-label transition is inferred.
overlap = resolved_chart_overlap_condition(
    "synthetic_overlap",
    ambient,
    lambda label: {
        "1": (1, 0), "z": (0, 1), "z2": (1, 1), "z3": (0, 0),
    }[label],
    lambda label: {
        "1": (0, 1), "z": (1, 0), "z2": (1, 1), "z3": (0, 0),
    }[label],
    matrix(QQ, [[1, 0], [0, 1]]),
    matrix(QQ, [[0, 1], [1, 0]]),
    ("overlap_0", "overlap_1"),
    "synthetic resolved-chart transition",
)
assert overlap["matrix"] == matrix(QQ, [[0, 0, 0, 0], [0, 0, 0, 0]])
try:
    resolved_chart_overlap_condition(
        "bad_overlap", ambient, lambda unused: (0,), lambda unused: (0,),
        matrix(QQ, [[1, 0]]), matrix(QQ, [[1]]), ("q",), "bad dimensions",
    )
except ValueError as error:
    assert "incompatible" in str(error)
else:
    raise AssertionError("incompatible overlap maps were accepted")

# Dimension-only compilation remains exact while avoiding eager rational
# nullspace materialization for large coefficient matrices.
dimension_only = compile_resolved_conditions(
    ambient, (first,), complete=True, compute_kernel=False
)
assert dimension_only["rank"] == 2
assert dimension_only["kernel_dimension"] == 2
assert dimension_only["kernel_basis"] is None
assert dimension_only["kernel_materialization"] == "not_requested"
assert dimension_only["h0_certified"]
try:
    require_h0_two(dimension_only)
except ValueError as error:
    assert "not materialized" in str(error)
else:
    raise AssertionError("dimension-only compilation returned a pencil basis")
assert certify_explicit_pencil_basis(
    dimension_only, matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]])
) == matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]])
try:
    certify_explicit_pencil_basis(dimension_only, matrix(QQ, [[1, 0, 0, 0], [0, 0, 0, 1]]))
except ValueError as error:
    assert "violates" in str(error)
else:
    raise AssertionError("invalid explicit pencil basis was accepted")

# If zero columns do not exhaust the nullity, retain them and compute only the
# residual active kernel.  This is the fallback used after sparse reduction.
mixed = quotient_condition(
    "mixed-sparse-chart", ambient,
    lambda value: {"1": (1,), "z": (1,), "z2": (0,), "z3": (0,)}[value],
    ("constant",), "synthetic chart with a residual active relation",
)
mixed_compilation = compile_resolved_conditions(ambient, (mixed,), complete=True)
assert mixed_compilation["rank"] == 1
assert mixed_compilation["kernel_dimension"] == 3
assert mixed_compilation["kernel_materialization"] == "zero_columns_plus_reduced_right_kernel"
assert mixed_compilation["condition_matrix"] * mixed_compilation["kernel_basis"].transpose() == matrix(QQ, 1, 3)

# The chart adapter starts from a stated local quotient and uses exact normal
# forms; it is intentionally not a mere valuation table.  This miniature
# chart produces the three standard quotient functionals on the ambient rows.
local_ring = PolynomialRing(QQ, names=("z", "u"), order="degrevlex")
z, u = local_ring.gens()
chart_block = resolved_chart_quotient_condition(
    "synthetic-resolved-chart",
    ambient,
    local_ring,
    lambda value: local_ring(value == "1") + local_ring(value == "z") * z
    + local_ring(value == "z2") * u,
    local_ring.ideal((z**2, z * u, u**2)),
    (1, z, u),
    "synthetic blow-up chart with an explicitly declared trivialization",
)
assert chart_block["matrix"] == matrix(QQ, [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
])

# A marked degree-two ambient uses caller-derived coefficients and the chord
# residue in the same resolved quotient machinery.  The result is not a
# valuation table: both the ideal and its normal-form basis are explicit.
chord_block = resolved_marked_chord_condition(
    "synthetic-marked-chord-chart",
    ("a0", "a1", "b0", "b1"),
    lambda label: {
        "a0": (1, 0), "a1": (u, 0),
        "b0": (0, 1), "b1": (0, u),
    }[label],
    local_ring,
    local_ring.ideal((u**2, z)),
    (1, u),
    1+u,
    "synthetic resolved marked chord",
)
assert chord_block["matrix"] == matrix(QQ, [
    [1, 0, 1, 0],
    [0, 1, 1, 1],
])

# A complete cover is declared as an ordered collection of actual quotient
# charts (and, when needed, caller-derived overlap rows).  It does not infer
# completeness from the resulting nullity.
cover = compile_resolved_chart_cover(
    ambient,
    ({
        "name": "synthetic-cover-chart",
        "local_ring": local_ring,
        "trivialized_pullback": lambda value: local_ring(value == "1")
        + local_ring(value == "z") * z,
        "quotient_ideal": local_ring.ideal((z**2, u)),
        "quotient_basis": (1, z),
        "provenance": "synthetic resolved-chart cover",
    },),
    complete=True,
)
assert len(cover["chart_blocks"]) == 1 and not cover["overlap_blocks"]
assert cover["compilation"]["condition_matrix"] == matrix(QQ, [
    [1, 0, 0, 0], [0, 1, 0, 0],
])
assert cover["compilation"]["h0_certified"]
try:
    compile_resolved_chart_cover(
        ambient,
        ({
            "name": "missing-provenance",
            "local_ring": local_ring,
            "trivialized_pullback": lambda unused: 0,
            "quotient_ideal": local_ring.ideal((z, u)),
            "quotient_basis": (1,),
        },),
    )
except ValueError as error:
    assert "provenance" in str(error)
else:
    raise AssertionError("incomplete chart record was accepted")

# One-variable jet quotients are finite too, although Sage does not expose
# vector_space_dimension or multivariate ``reduce`` on their ideals.
jet_ring = PolynomialRing(QQ, "v")
v = jet_ring.gen()
jet_block = resolved_marked_chord_condition(
    "synthetic-univariate-jet",
    ("a0", "b0"),
    lambda label: (1, 0) if label == "a0" else (0, 1),
    jet_ring,
    jet_ring.ideal(v**2),
    (1, v),
    1+v,
    "synthetic finite jet quotient",
)
assert jet_block["matrix"] == matrix(QQ, [[1, 1], [0, 1]])

incomplete = compile_resolved_conditions(ambient, (first,), complete=False)
assert not incomplete["h0_certified"]
try:
    require_h0_two(incomplete)
except ValueError as error:
    assert "incomplete" in str(error)
else:
    raise AssertionError("incomplete chart cover incorrectly certified h0")

assert bounded_weierstrass_monomials(6, (0,)) == (
    {"t_power": 0, "x_power": 0, "y_power": 0},
    {"t_power": 0, "x_power": 1, "y_power": 0},
    {"t_power": 0, "x_power": 2, "y_power": 0},
    {"t_power": 0, "x_power": 0, "y_power": 1},
    {"t_power": 0, "x_power": 1, "y_power": 1},
)
assert bounded_weierstrass_monomials(2, (0,)) == (
    {"t_power": 0, "x_power": 0, "y_power": 0},
)
assert marked_section_generic_fibre_basis(2, "m_-P1") == (
    {"kind": "weierstrass_monomial", "x_power": 0, "y_power": 0},
    {"kind": "marked_chord", "symbol": "m_-P1"},
)
assert marked_section_generic_fibre_basis(1, "unused") == (
    {"kind": "weierstrass_monomial", "x_power": 0, "y_power": 0},
)
assert balanced_marked_chord_power_basis(1) == (
    {"kind": "m_power", "m_power": 0, "x_power": 0, "chord_symbol": "m",
     "pole_order_at_O": 0, "pole_order_at_marked_section": 0},
    {"kind": "m_power", "m_power": 1, "x_power": 0, "chord_symbol": "m",
     "pole_order_at_O": 1, "pole_order_at_marked_section": 1},
)
balanced_q8_basis = balanced_marked_chord_power_basis(9, "m_-P1", "x")
assert len(balanced_q8_basis) == 18
assert balanced_q8_basis[-1] == {
    "kind": "x_m_power", "m_power": 7, "x_power": 1,
    "chord_symbol": "m_-P1", "x_symbol": "x",
    "pole_order_at_O": 9, "pole_order_at_marked_section": 7,
}
# Before any monomial space is constructed, the generic-fibre marked point
# and vertical correction must reconstruct the actual Neron--Severi class.
# This tiny hyperbolic test also rejects a purported vertical curve of degree
# one against the old fibre.
decomposition_gram = matrix(ZZ, ((0, 1), (1, 0)))
decomposition_fibre = vector(ZZ, (1, 0))
decomposition_zero = vector(ZZ, (-1, 1))
decomposition_divisor = 2*decomposition_zero-decomposition_fibre
decomposition = certify_generic_fibre_divisor_decomposition(
    decomposition_gram, decomposition_divisor, decomposition_fibre,
    decomposition_zero, decomposition_zero, (), fiber_twist=-1,
    expected_old_fiber_degree=2,
)
assert decomposition["generic_restriction"] == {
    "zero_section_coefficient": 1, "marked_section_coefficient": 1,
}
assert decomposition["reconstructed_divisor"] == tuple(decomposition_divisor)
try:
    certify_generic_fibre_divisor_decomposition(
        decomposition_gram, decomposition_divisor, decomposition_fibre,
        decomposition_zero, decomposition_zero,
        (("not_vertical", 1, decomposition_zero),), fiber_twist=-1,
    )
except ValueError as error:
    assert "vertical" in str(error)
else:
    raise AssertionError("positive-degree vertical support was accepted")
# Higher degree inputs retain a literal repeated-section representative until
# the caller has separately certified its group-law/chord conversion.
supported_divisor = 3*decomposition_zero-2*decomposition_fibre
support_certificate = certify_generic_fibre_horizontal_support(
    decomposition_gram, supported_divisor, decomposition_fibre,
    (("O", 2, decomposition_zero), ("Q", 1, decomposition_zero)),
    (), fiber_twist=-2, expected_old_fiber_degree=3,
)
assert support_certificate["old_fiber_degree"] == 3
assert [row["multiplicity"] for row in support_certificate["horizontal_support"]] == [2, 1]
try:
    certify_generic_fibre_horizontal_support(
        decomposition_gram, supported_divisor, decomposition_fibre,
        (("not_a_section", 1, decomposition_fibre),), (), fiber_twist=-2,
    )
except ValueError as error:
    assert "section" in str(error)
else:
    raise AssertionError("degree-zero horizontal curve was accepted")
assert endpoint_coefficient_interval(11, 6, 4) == {
    "denominator_power": 2, "u_power_lower": 11, "u_power_upper": 14,
}
assert endpoint_coefficient_interval(17, 9, 4) == {
    "denominator_power": 2, "u_power_lower": 17, "u_power_upper": 17,
}
try:
    endpoint_coefficient_interval(1, 0, 0)
except ValueError as error:
    assert "denominator degree" in str(error)
else:
    raise AssertionError("nonpositive denominator degree was accepted")

# High-height marked sections are evaluated in the caller's local arithmetic
# from an exact group-law DAG, rather than globally normalized first.
point_dag = {
    "operation": "add",
    "left": {"operation": "scalar", "scalar": 22,
             "point": {"operation": "negate", "point": "P1"}},
    "right": "reconstructed_-P2",
}
assert evaluate_marked_point_dag(
    point_dag, {"P1": 3, "reconstructed_-P2": 5},
    lambda left, right: left + right, lambda value: -value,
    lambda multiplier, value: multiplier * value,
) == -61
try:
    evaluate_marked_point_dag("missing", {}, lambda a, b: a + b, lambda a: -a,
                              lambda a, b: a * b)
except ValueError as error:
    assert "unknown" in str(error)
else:
    raise AssertionError("unknown DAG leaf was accepted")

# Chamber moves must preserve an auditable ordered pairing record, rather than
# merely exhibit an isometric target vector.
root_gram = matrix(QQ, [[-2]])
reflected, reflection_record = replay_weyl_reflections(
    vector([1]), root_gram, (("A1", vector([1])),), expected_pairings=(-2,)
)
assert reflected == vector([-1])
assert reflection_record == ({"root": "A1", "pairing_before": -2},)

# The reusable orchestration layer binds the ordered chamber move, declared
# nef check, complete resolved condition cover, and displayed pencil to one
# exact object.  It cannot infer a chart condition from the root label.
neighbor_gram = matrix(ZZ, [[0, 1, 0], [1, 0, 0], [0, 0, -2]])
neighbor_raw = vector(ZZ, [1, 1, 1])
neighbor_fiber = vector(ZZ, [1, 0, 0])
neighbor_root = vector(ZZ, [0, 0, 1])
neighbor_pipeline = compile_elliptic_neighbor_rr_pencil(
    neighbor_gram,
    neighbor_raw,
    neighbor_fiber,
    (("A1", neighbor_root),),
    (("A1", neighbor_root),),
    ambient,
    (first,),
    complete_resolved_cover=True,
    expected_reflection_pairings=(-2,),
    expected_nef_divisor=(1, 1, -1),
    pencil_basis=matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]]),
)
assert neighbor_pipeline["nef_divisor"] == vector(ZZ, [1, 1, -1])
assert neighbor_pipeline["preflight"]["old_fiber_degree"] == 1
assert neighbor_pipeline["rr"]["h0_certified"]
assert neighbor_pipeline["pencil_basis"].nrows() == 2
try:
    compile_elliptic_neighbor_rr_pencil(
        neighbor_gram, neighbor_raw, neighbor_fiber,
        (("A1", neighbor_root),), (("A1", neighbor_root),), ambient, (first,),
        complete_resolved_cover=True, expected_reflection_pairings=(-2,),
        expected_nef_divisor=(1, 1, 1),
    )
except ValueError as error:
    assert "expected nef" in str(error)
else:
    raise AssertionError("wrong Weyl endpoint was accepted")

# A degree-independent exact-neighbour lattice certificate can enter the
# equation layer without relabelling its old-fibre degree or promoting its
# supplied-wall statement to a global nef proof.
with open("artifacts/generated-results/elkies-k3-exact-neighbor-h3-d13-q24.json") as handle:
    q24_lattice_certificate = json.load(handle)
certificate_handoff = exact_neighbor_certificate_handoff(
    q24_lattice_certificate, expected_old_fiber_degree=2
)
assert certificate_handoff["child_root_data"] == (12, 264, 4)
assert certificate_handoff["child_mw_height"].nrows() == 5
assert certificate_handoff["proof_boundary"].startswith("The engine proves only")
try:
    exact_neighbor_certificate_handoff(q24_lattice_certificate, expected_old_fiber_degree=24)
except ValueError as error:
    assert "old-fibre degree" in str(error)
else:
    raise AssertionError("search-shell label was accepted as an old-fibre degree")

# A degree-two marked chord has an exact quartic/Jacobian route.  Retaining
# the square factor check is essential: throwing away a non-square scalar
# would instead compile a quadratic twist of the requested child.
base_ring = PolynomialRing(QQ, "t")
t = base_ring.gen()
base_field = base_ring.fraction_field()
x_point, y_point, old_a, old_b = (
    base_field(0), base_field(1), base_field(1), base_field(1)
)
assert y_point**2 == x_point**3 + old_a * x_point + old_b
assert chord_tangent_slope(x_point, y_point, old_a) == QQ(1) / 2
parameter = base_field.gen()
chord = pencil_chord_solution(1, 1, 0, 1, parameter)
assert pencil_value_on_marked_section(1, 1, 0, 1, chord) == parameter
assert marked_chord_value(x_point, y_point, x_point, y_point, old_a, True) == QQ(1) / 2
assert rational_map_degree(parameter) == 1
radicand = chord_discriminant(x_point, y_point, old_a, chord)
quartic, square_factor = squarefree_binary_quartic(radicand, base_ring)
assert radicand == square_factor**2 * quartic
assert quartic.degree() in (3, 4)
coefficient_a, coefficient_b, discriminant = binary_quartic_jacobian_coefficients(quartic)
assert discriminant and coefficient_a and coefficient_b
compiled_chord_hop = compile_degree_two_chord_hop(
    base_ring,
    parameter,
    1, 1, 0, 1,
    x_point, y_point, old_a, old_b,
    marked_chords=(("P", QQ(1) / 2),),
)
assert compiled_chord_hop["chord"] == chord
assert compiled_chord_hop["radicand"] == radicand
assert compiled_chord_hop["binary_quartic"] == quartic
assert compiled_chord_hop["jacobian_a"] == coefficient_a
assert compiled_chord_hop["jacobian_b"] == coefficient_b
assert compiled_chord_hop["transported_parameter_values"]["P"] == QQ(1) / 3
try:
    compile_degree_two_chord_hop(
        base_ring, parameter, 1, 1, 0, 1,
        x_point, y_point, old_a, old_b + 1,
    )
except ValueError as error:
    assert "does not lie" in str(error)
else:
    raise AssertionError("off-curve marked point was accepted")

# A resolved binary-quartic point transports to the exact Jacobian through
# covariants; this is a coordinate transport, never a search for a child
# section.  The scaling branch is checked too.
quartic_ring = PolynomialRing(QQ, "u")
u = quartic_ring.gen()
transport_quartic = u**4-u**2+1
covariant_transport = transport_binary_quartic_point_to_jacobian(
    transport_quartic, 0, 1, 1, QQ(2),
)
assert covariant_transport["standard_y"]**2 == (
    covariant_transport["standard_x"]**3
    + covariant_transport["standard_a"]*covariant_transport["standard_x"]
    + covariant_transport["standard_b"]
)
assert covariant_transport["child_y"]**2 == (
    covariant_transport["child_x"]**3
    + covariant_transport["child_a"]*covariant_transport["child_x"]
    + covariant_transport["child_b"]
)
try:
    transport_binary_quartic_point_to_jacobian(transport_quartic, 0, 1, 2)
except ValueError as error:
    assert "does not satisfy" in str(error)
else:
    raise AssertionError("off-quartic transport point was accepted")

# The degree-two conversion must be linked to a complete resolved RR matrix,
# not merely called with a displayed pair of chord functions.
rr_ambient = ("c0", "c1", "c2", "c3")
rr_block = quotient_condition(
    "synthetic_resolved_block", rr_ambient,
    lambda label: {
        "c0": (1, 0), "c1": (0, 1), "c2": (0, 0), "c3": (0, 0),
    }[label],
    ("q0", "q1"), "synthetic resolved quotient",
)
rr_compilation = compile_resolved_conditions(
    rr_ambient, (rr_block,), complete=True,
)
resolved_hop = compile_resolved_degree_two_chord_hop(
    rr_compilation,
    matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]]),
    ((0, 0), (0, 0), (1, 1), (0, 1)),
    base_ring, parameter, x_point, y_point, old_a, old_b,
    marked_chords=(("P", QQ(1) / 2),),
)
assert resolved_hop["chord_coefficients"] == ((1, 1), (0, 1))
assert resolved_hop["conversion"]["chord"] == chord
try:
    compile_resolved_degree_two_chord_hop(
        rr_compilation, matrix(QQ, [[1, 0, 0, 0], [0, 1, 0, 0]]),
        ((0, 0), (0, 0), (1, 1), (0, 1)),
        base_ring, parameter, x_point, y_point, old_a, old_b,
    )
except ValueError as error:
    assert "violates" in str(error)
else:
    raise AssertionError("non-kernel pencil was accepted for degree-two conversion")

# The child Jacobian is compiled directly from that same resolved pencil.
# Its finite-fibre data may be used for root bookkeeping, but its infinity
# record is deliberately not promoted to a resolved additive-fibre chart.
resolved_child = compile_resolved_degree_two_child_jacobian(
    rr_compilation,
    matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]]),
    ((0, 0), (0, 0), (1, 1), (0, 1)),
    base_ring, base_ring, parameter, x_point, y_point, old_a, old_b,
    marked_chords=(("P", QQ(1) / 2),),
)
assert resolved_child["resolved_hop"]["conversion"]["binary_quartic"] == quartic
assert resolved_child["jacobian_a"] == coefficient_a
assert resolved_child["jacobian_b"] == coefficient_b
assert "infinity_boundary" in resolved_child["finite_classification"]
try:
    compile_resolved_degree_two_child_jacobian(
        rr_compilation,
        matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]]),
        ((0, 0), (0, 0), (1, 1), (0, 1)),
        base_ring, PolynomialRing(QQ, names=("s", "t")), parameter,
        x_point, y_point, old_a, old_b,
    )
except ValueError as error:
    assert "new base ring" in str(error)
else:
    raise AssertionError("multivariate new base was accepted for a curve fibration")

# A higher-degree compiler can retain an exact raw genus-one relation by
# eliminating only the selected old Weierstrass variable.  This is not yet a
# minimization or a Jacobian claim.
elimination_ring = PolynomialRing(QQ, names=("T", "u", "x", "y"), order="lex")
new_t, old_u, old_x, old_y = elimination_ring.gens()
elimination_data = eliminate_cleared_weierstrass_pencil(
    elimination_ring,
    old_y**2-old_x**3-old_u*old_x-1,
    1,
    old_x,
    new_t,
    (old_x,),
)
assert elimination_data["pencil_relation"] == old_x-new_t
expected_eliminated_relation = old_y**2-new_t**3-old_u*new_t-1
assert len(elimination_data["relations"]) == 1
assert expected_eliminated_relation in elimination_data["elimination_ideal"]

# The higher-degree bridge must derive that same relation from a pencil that
# is certified against the resolved RR matrix, rather than taking two free
# cleared functions as an independent input.
resolved_elimination = compile_resolved_genus_one_elimination(
    rr_compilation,
    matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]]),
    (old_y, old_u, 1, old_x),
    elimination_ring,
    old_y**2-old_x**3-old_u*old_x-1,
    new_t,
    (old_x,),
)
assert resolved_elimination["pencil_zero"] == 1
assert resolved_elimination["pencil_one"] == old_x
assert expected_eliminated_relation in resolved_elimination["elimination"]["elimination_ideal"]
try:
    compile_resolved_genus_one_elimination(
        rr_compilation, matrix(QQ, [[1, 0, 0, 0], [0, 1, 0, 0]]),
        (old_y, old_u, 1, old_x), elimination_ring,
        old_y**2-old_x**3-old_u*old_x-1, new_t, (old_x,),
    )
except ValueError as error:
    assert "violates" in str(error)
else:
    raise AssertionError("non-kernel pencil was accepted for genus-one elimination")
try:
    eliminate_cleared_weierstrass_pencil(
        elimination_ring, old_y**2-old_x**3, 1, old_x, new_t, (),
    )
except ValueError as error:
    assert "nonempty" in str(error)
else:
    raise AssertionError("empty variable elimination was accepted")

# A cleared denominator can introduce an entire affine component.  The core
# removes it only when the exact chart denominator is explicitly declared.
spurious_data = eliminate_cleared_weierstrass_pencil(
    elimination_ring,
    old_x*(old_y**2-old_x**3-old_u*old_x-1),
    old_x,
    old_x**2,
    new_t,
    (old_x,),
    saturate_by=(old_x,),
)
assert spurious_data["saturation_product"] == old_x
assert spurious_data["saturated_ideal"] != spurious_data["source_ideal"]
assert expected_eliminated_relation in spurious_data["elimination_ideal"]
try:
    eliminate_cleared_weierstrass_pencil(
        elimination_ring, old_y**2-old_x**3, 1, old_x, new_t, (old_x,),
        saturate_by=(0,),
    )
except ValueError as error:
    assert "zero" in str(error)
else:
    raise AssertionError("zero saturation factor was accepted")

# The higher-degree plane-cubic adapter requires the exact tangent frame; it
# does not search for a flex or choose an origin implicitly.
plane_ring = PolynomialRing(QQ, names=("X", "Y", "Z"))
plane_x, plane_y, plane_z = plane_ring.gens()
plane_cubic = plane_y**2*plane_z-plane_x**3-plane_x*plane_z**2-plane_z**3

def tangent_third(point):
    px, py = (QQ(value) for value in point)
    affine = PolynomialRing(QQ, names=("x", "y"))
    ax, ay = affine.gens()
    relation = ay**2-ax**3-ax-1
    derivative_x = relation.derivative(ax)(px, py)
    derivative_y = relation.derivative(ay)(px, py)
    assert derivative_y
    tangent = py-derivative_x/derivative_y*(ax-px)
    intersection = PolynomialRing(QQ, "x")(relation(ax, tangent))
    third_x = -intersection[2]/intersection[3]-2*px
    return third_x, QQ(tangent.subs({ax: third_x}))

frame_p = (QQ(0), QQ(1))
frame_q = tangent_third(frame_p)
frame_r = tangent_third(frame_q)
plane_data = pointed_plane_cubic_to_weierstrass(
    plane_cubic,
    ((frame_p[0], frame_p[1], 1),
     (frame_q[0], frame_q[1], 1),
     (frame_r[0], frame_r[1], 1)),
)
assert plane_data["curve"].discriminant()
assert len(plane_data["source_to_weierstrass"]) == 3
try:
    pointed_plane_cubic_to_weierstrass(
        plane_cubic,
        ((0, 1, 1), (0, 1, 1), (1, 0, 1)),
    )
except (ValueError, ArithmeticError) as error:
    assert "frame" in str(error) or "point" in str(error)
else:
    raise AssertionError("dependent tangent frame was accepted")

# Finite Weierstrass minimization records the exact fourth/sixth-power scale
# and reports the independent infinity chart instead of overclaiming global
# minimality from affine valuations alone.
min_ring = PolynomialRing(QQ, "s")
s = min_ring.gen()
finite_minimization = finite_minimize_short_weierstrass(min_ring, s**4, s**6)
assert finite_minimization["scaling_unit"] == 1/s
assert finite_minimization["minimal_a"] == 1
assert finite_minimization["minimal_b"] == 1
assert finite_minimization["finite_places"][0]["minimal_orders"] == (0, 0, 0)
try:
    finite_minimize_short_weierstrass(min_ring, 0, 0)
except ValueError as error:
    assert "singular" in str(error)
else:
    raise AssertionError("singular short model was accepted")
assert kodaira_data_from_short_orders(3, 4, 8) == (6, 8, 3, "IV*")
assert kodaira_data_from_short_orders(4, 5, 10) == (8, 10, 1, "II*")
assert kodaira_data_from_short_orders(0, 0, 1) == (0, 1, 1, "I1")
finite_classification = classify_finite_short_weierstrass_fibres(
    min_ring, s**3, s**4
)
assert [(item["degree"], item["minimal_orders"], item["kodaira"])
        for item in finite_classification["finite_fibres"]] == [
    (1, (3, 4, 8), "IV*"), (1, (0, 0, 1), "I1"),
]
assert finite_classification["finite_root_rank"] == 6
assert finite_classification["finite_euler_number"] == 9
assert finite_classification["finite_root_determinant"] == 3
assert "normalized_orders" in finite_classification["infinity_boundary"]
try:
    kodaira_data_from_short_orders(1, 1, 1)
except ValueError as error:
    assert "unrecognized" in str(error)
else:
    raise AssertionError("invalid Kodaira valuation triple was accepted")
ns_discriminant = certify_shioda_tate_discriminant(
    3,
    [[QQ(8)/3, QQ(1)/3, -1], [QQ(1)/3, QQ(8)/3, 1], [-1, 1, 46]],
    expected_ns_discriminant=948,
)
assert ns_discriminant["height_determinant"] == 316
assert ns_discriminant["absolute_ns_discriminant"] == 948
try:
    certify_shioda_tate_discriminant(3, [[1]], torsion_order=2)
except ArithmeticError as error:
    assert "integral" in str(error)
else:
    raise AssertionError("nonintegral Shioda--Tate discriminant was accepted")

# Resolved component transport is checked against explicit chart pairings and
# degree balance; names or a Kodaira label cannot manufacture this evidence.
component_transport = certify_component_pairing_transport(
    {"zero": "old O", "first": "old E7_7 exceptional component"},
    {"zero": 1, "first": 1},
    ("old_E7_1", "old_E7_2"),
    {"horizontal": (0, 0), "corrected": (0, 22)},
    {"horizontal_degree": 23, "correction_degree": -22,
     "transported_section_degree": 1},
    {
        "section_sources": {"zero": "old O", "first": "old E7_7 exceptional component"},
        "source_fiber_degrees": {"zero": 1, "first": 1},
        "pairing_basis": ("old_E7_1", "old_E7_2"),
        "resolved_pairings": {"horizontal": (0, 0), "corrected": (0, 22)},
        "vertical_correction": {"horizontal_degree": 23, "correction_degree": -22,
                                "transported_section_degree": 1},
    },
)
assert component_transport["vertical_correction"]["transported_section_degree"] == 1
try:
    certify_component_pairing_transport(
        {"zero": "old O"}, {"zero": 2}, ("old_E7_1",), {"row": (0,)},
        {"horizontal_degree": 1, "correction_degree": 0, "transported_section_degree": 1},
        {"section_sources": {"zero": "old O"}, "source_fiber_degrees": {"zero": 2},
         "pairing_basis": ("old_E7_1",), "resolved_pairings": {"row": (0,)},
         "vertical_correction": {"horizontal_degree": 1, "correction_degree": 0,
                                 "transported_section_degree": 1}},
    )
except ValueError as error:
    assert "not a section" in str(error)
else:
    raise AssertionError("non-section component transport was accepted")

# The end-to-end adapter binds the independently certified lattice, RR,
# equation, and transport layers.  It insists on the exact target Gram and
# section words, not merely matching ranks.
hop = certify_exact_neighbor_hop(
    {
        "square": 0, "primitive": True, "old_fiber_degree": 2,
        "nef_on_declared_walls": True, "weyl_reflection_count": 3,
    },
    {
        "complete_resolved_chart_cover": True, "ambient_dimension": 4,
        "condition_rank": 2, "condition_codimension": 2,
        "kernel_dimension": 2, "h0": 2,
    },
    {
        "root_lattice": "A1", "root_rank": 1,
        "root_determinant": 2, "mordell_weil_rank": 1,
    },
    {
        "height_gram": [[QQ(2)]], "section_words": ["P"],
        "section_new_fiber_degrees": [1],
    },
    {
        "rr": {"ambient_dimension": 4, "condition_rank": 2,
               "condition_codimension": 2, "kernel_dimension": 2, "h0": 2},
        "child": {"root_lattice": "A1", "root_rank": 1,
                  "root_determinant": 2, "mordell_weil_rank": 1},
        "height_gram": [[QQ(2)]], "section_words": ["P"],
    },
)
assert hop["rr"]["h0"] == 2 and hop["transport"]["height_gram"] == [["2"]]
try:
    certify_exact_neighbor_hop(
        {
            "square": 0, "primitive": True, "old_fiber_degree": 2,
            "nef_on_declared_walls": True, "weyl_reflection_count": 0,
        },
        {
            "complete_resolved_chart_cover": True, "ambient_dimension": 4,
            "condition_rank": 2, "condition_codimension": 2,
            "kernel_dimension": 2, "h0": 2,
        },
        {"root_lattice": "A1", "root_rank": 1,
         "root_determinant": 2, "mordell_weil_rank": 1},
        {"height_gram": [[QQ(3)]], "section_words": ["P"],
         "section_new_fiber_degrees": [1]},
        {"rr": {"ambient_dimension": 4, "condition_rank": 2,
                "condition_codimension": 2, "kernel_dimension": 2, "h0": 2},
         "child": {"root_lattice": "A1", "root_rank": 1,
                   "root_determinant": 2, "mordell_weil_rank": 1},
         "height_gram": [[QQ(2)]], "section_words": ["P"]},
    )
except ValueError as error:
    assert "height" in str(error)
else:
    raise AssertionError("wrong transported Gram was accepted")

# A small kernel restriction gives the rank of a stacked modular condition
# matrix without forming that stack's full row reduction.
finite = GF(43)
overlay = modular_condition_overlay_rank(
    matrix(finite, [[1, 0, 0], [0, 1, 0]]),
    matrix(finite, [[0, 0, 1], [1, 0, 0]]),
)
assert overlay["base_rank"] == 2
assert overlay["base_kernel_dimension"] == 1
assert overlay["overlay_rank_on_base_kernel"] == 1
assert overlay["stacked_rank"] == 3 and overlay["stacked_kernel_dimension"] == 0

print("ELLIPTICNEIGHBOR|ambient=4|conditions=2|rank=2|kernel=2|status=PASS")
