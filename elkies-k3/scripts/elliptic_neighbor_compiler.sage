#!/usr/bin/env sage
"""Exact linear-algebra core for equation-level elliptic neighbours.

This module intentionally separates a lattice-neighbour search label from the
old-fibre degree of the resulting divisor.  The latter is the ``q`` in the
Riemann--Roch calculation.  The first H3 ``q6`` shell, for example, has
old-fibre degree two after fixed components are removed.

``compile_resolved_conditions`` accepts an ambient basis and maps to finite
quotients of actual resolved charts.  It performs no inference from a Kodaira
symbol.  An incomplete chart list can yield a diagnostic kernel, but only a
complete list can certify ``h0(D)=2``.
"""

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, gcd, matrix, prod, vector


def intersection(left, right, gram):
    """Exact intersection of row vectors in an integral Neron--Severi frame."""
    return ZZ(vector(ZZ, left) * matrix(ZZ, gram) * vector(ZZ, right))


def primitive(vector_value):
    """Return whether an integral vector has content one."""
    return gcd([ZZ(entry) for entry in vector(ZZ, vector_value)]) == 1


def weyl_reflection(class_value, root, gram):
    """Reflect an NS class in one declared ``(-2)`` root exactly."""
    class_value = vector(ZZ, class_value)
    root = vector(ZZ, root)
    if intersection(root, root, gram) != -2:
        raise ValueError("Weyl reflection root does not have square -2")
    return class_value + intersection(class_value, root, gram) * root


def replay_weyl_reflections(class_value, gram, roots, expected_pairings=None):
    """Apply and record an ordered Weyl reflection sequence.

    ``roots`` is an ordered sequence of ``(name, root)`` pairs.  Supplying
    ``expected_pairings`` makes the chamber move a certificate rather than a
    bare lattice isometry; each value is checked before its reflection.
    """
    roots = tuple((str(name), vector(ZZ, root)) for name, root in roots)
    if expected_pairings is not None and len(expected_pairings) != len(roots):
        raise ValueError("Weyl pairing record has the wrong length")
    result = vector(ZZ, class_value)
    record = []
    for index, (name, root) in enumerate(roots):
        pairing = intersection(result, root, gram)
        if expected_pairings is not None and pairing != ZZ(expected_pairings[index]):
            raise ValueError(
                "Weyl reflection {} has pairing {}, expected {}".format(
                    name, pairing, expected_pairings[index]
                )
            )
        record.append({"root": name, "pairing_before": int(pairing)})
        result = weyl_reflection(result, root, gram)
    return result, tuple(record)


def isotropic_nef_preflight(gram, divisor, old_fiber, wall_curves):
    """Check the finite, declared part of a primitive-nef certificate.

    ``wall_curves`` is an ordered sequence of ``(name, curve)`` pairs.  This
    does not assert that it is the full effective cone; that assertion belongs
    to the caller's independent chamber certificate.
    """
    divisor = vector(ZZ, divisor)
    old_fiber = vector(ZZ, old_fiber)
    pairings = []
    for name, curve in wall_curves:
        curve = vector(ZZ, curve)
        if intersection(curve, curve, gram) != -2:
            raise ValueError("declared wall {} is not a (-2)-curve".format(name))
        pairings.append((str(name), int(intersection(divisor, curve, gram))))
    return {
        "square": int(intersection(divisor, divisor, gram)),
        "primitive": bool(primitive(divisor)),
        "old_fiber_degree": int(intersection(divisor, old_fiber, gram)),
        "declared_wall_pairings": pairings,
        "nonnegative_on_declared_walls": all(value >= 0 for _, value in pairings),
    }


def bounded_weierstrass_monomials(old_fiber_degree, base_powers):
    """Return ``t^i*x^a*y^b`` labels for a standard ambient RR space.

    The generic-fibre pole bound is ``old_fiber_degree-1``.  ``base_powers``
    must be explicit: the vertical twist is geometry, not a guess.
    """
    degree = ZZ(old_fiber_degree)
    if degree < 1:
        raise ValueError("old-fibre degree must be positive")
    answer = []
    for i in base_powers:
        i = ZZ(i)
        for b in (ZZ(0), ZZ(1)):
            for a in range((degree - 1 - 3 * b) // 2 + 1):
                answer.append({"t_power": int(i), "x_power": int(a), "y_power": int(b)})
    return tuple(answer)


def certify_generic_fibre_divisor_decomposition(
    gram,
    divisor,
    old_fiber,
    zero_section,
    marked_section,
    vertical_support=(),
    fiber_twist=0,
    expected_old_fiber_degree=None,
):
    """Certify ``D=(q-1)O+P+V+kF`` before building an RR ambient.

    A lattice neighbour vector alone does not say which marked point occurs
    on the old generic fibre.  This adapter makes that indispensable choice
    explicit: the caller supplies the zero section ``O``, a marked section
    ``P``, the resolved vertical support of ``V``, and the integer fibre
    twist ``k``.  It then checks the equality in the declared Neron--Severi
    frame and the fibre degrees of every displayed term.

    ``vertical_support`` is an ordered sequence of ``(name, coefficient,
    component)`` records.  Each component must have degree zero against the
    old fibre; this verifies it is vertical *relative to the supplied old
    fibration*, but does not claim that the caller's list is the full set of
    fibre components.  The returned generic restriction is therefore a
    certificate of the form ``(q-1)O+P``, not a procedure that discovers P.
    """
    gram = matrix(ZZ, gram)
    divisor = vector(ZZ, divisor)
    old_fiber = vector(ZZ, old_fiber)
    zero_section = vector(ZZ, zero_section)
    marked_section = vector(ZZ, marked_section)
    dimension = len(divisor)
    if gram.nrows() != gram.ncols() or gram.nrows() != dimension:
        raise ValueError("Neron--Severi Gram matrix and divisor have incompatible sizes")
    if any(len(value) != dimension for value in (old_fiber, zero_section, marked_section)):
        raise ValueError("generic-fibre decomposition has incompatible Neron--Severi sizes")
    q = intersection(divisor, old_fiber, gram)
    if q <= 0:
        raise ValueError("target divisor has nonpositive old-fibre degree")
    if expected_old_fiber_degree is not None and q != ZZ(expected_old_fiber_degree):
        raise ValueError("target divisor has the wrong old-fibre degree")
    if intersection(zero_section, old_fiber, gram) != 1:
        raise ValueError("declared zero section does not have old-fibre degree one")
    if intersection(marked_section, old_fiber, gram) != 1:
        raise ValueError("declared marked section does not have old-fibre degree one")
    vertical = vector(ZZ, [0]*dimension)
    records = []
    for item in vertical_support:
        if len(item) != 3:
            raise ValueError("vertical support record must be (name, coefficient, component)")
        name, coefficient, component = item
        coefficient = ZZ(coefficient)
        component = vector(ZZ, component)
        if len(component) != dimension:
            raise ValueError("vertical component has incompatible Neron--Severi size")
        degree = intersection(component, old_fiber, gram)
        if degree != 0:
            raise ValueError("declared vertical component {} has old-fibre degree {}".format(name, degree))
        vertical += coefficient*component
        records.append({
            "name": str(name), "coefficient": int(coefficient),
            "old_fiber_degree": int(degree),
        })
    fiber_twist = ZZ(fiber_twist)
    generic_restriction = (q-1)*zero_section+marked_section
    reconstructed = generic_restriction+vertical+fiber_twist*old_fiber
    if reconstructed != divisor:
        raise ValueError("declared generic-fibre decomposition does not reconstruct target divisor")
    return {
        "old_fiber_degree": int(q),
        "zero_section_degree": int(intersection(zero_section, old_fiber, gram)),
        "marked_section_degree": int(intersection(marked_section, old_fiber, gram)),
        "generic_restriction": {
            "zero_section_coefficient": int(q-1),
            "marked_section_coefficient": 1,
        },
        "vertical_support": tuple(records),
        "fiber_twist": int(fiber_twist),
        "reconstructed_divisor": tuple(int(value) for value in reconstructed),
    }


def endpoint_coefficient_interval(e8_floor, e7_pole, denominator_degree):
    """Return the least base-coefficient interval satisfying two endpoint bounds.

    A coefficient ``u^i/h(u)^k`` has order ``i`` at ``u=0`` when ``h(0)``
    is a unit, and has t-order ``denominator_degree*k-i`` at ``u=infinity``.
    Given a required floor ``i >= e8_floor`` and an allowed t-pole of order
    ``e7_pole``, choose the least nonnegative ``k`` with

        e8_floor <= denominator_degree*k + e7_pole.

    The returned inclusive interval of allowed ``i`` is exact for those two
    endpoint inequalities.  Finite collision and other resolved-chart
    conditions remain separate compiler blocks.
    """
    e8_floor = ZZ(e8_floor)
    e7_pole = ZZ(e7_pole)
    denominator_degree = ZZ(denominator_degree)
    if e8_floor < 0 or e7_pole < 0 or denominator_degree <= 0:
        raise ValueError("endpoint bounds must be nonnegative and denominator degree positive")
    gap = e8_floor-e7_pole
    denominator_power = max(ZZ(0), -((-gap) // denominator_degree))
    upper = denominator_degree*denominator_power+e7_pole
    if upper < e8_floor:
        raise ArithmeticError("endpoint interval construction failed")
    return {
        "denominator_power": int(denominator_power),
        "u_power_lower": int(e8_floor),
        "u_power_upper": int(upper),
    }


def evaluate_marked_point_dag(expression, named_points, add, negate, scalar):
    """Evaluate an exact marked-point expression in a caller's local model.

    A high-height marked point should not be forced into one globally reduced
    rational-function pair before its divisor conditions are imposed.  The
    compiler therefore accepts a small expression DAG whose leaves name
    certified points and whose operations are supplied by the caller.  A
    resolved chart can use its own local-ring, projective, or quotient-point
    arithmetic; the same DAG then has identical marked-divisor semantics in
    every chart.

    Supported nodes are a leaf string and dictionaries with ``operation``
    equal to ``add``, ``negate``, or ``scalar``.  ``scalar`` takes an integer
    ``scalar`` and a child ``point``.
    """
    if isinstance(expression, str):
        try:
            return named_points[expression]
        except KeyError as error:
            raise ValueError("unknown marked-point leaf {}".format(expression)) from error
    if not isinstance(expression, dict) or "operation" not in expression:
        raise ValueError("marked-point expression must be a leaf or operation dictionary")
    operation = expression["operation"]
    if operation == "add":
        return add(
            evaluate_marked_point_dag(expression["left"], named_points, add, negate, scalar),
            evaluate_marked_point_dag(expression["right"], named_points, add, negate, scalar),
        )
    if operation == "negate":
        return negate(evaluate_marked_point_dag(
            expression["point"], named_points, add, negate, scalar
        ))
    if operation == "scalar":
        return scalar(
            ZZ(expression["scalar"]),
            evaluate_marked_point_dag(expression["point"], named_points, add, negate, scalar),
        )
    raise ValueError("unsupported marked-point operation {}".format(operation))


def quotient_condition(name, ambient_basis, evaluator, quotient_basis, provenance):
    """Build one exact chart-condition block.

    ``evaluator`` maps one ambient element to coordinates in the specified
    finite quotient.  Keeping this map chart-specific makes its blow-up and
    non-Cartier module provenance auditable.
    """
    quotient_basis = tuple(quotient_basis)
    columns = []
    for basis_element in ambient_basis:
        residue = vector(QQ, evaluator(basis_element))
        if len(residue) != len(quotient_basis):
            raise ValueError("{} returned a residue of the wrong length".format(name))
        columns.append(residue)
    return {
        "name": str(name),
        "matrix": matrix(QQ, columns).transpose(),
        "quotient_basis": tuple(map(str, quotient_basis)),
        "provenance": str(provenance),
    }


def finite_ambient_image_condition(
    name, ambient_basis, evaluator, coordinate_order, coefficient_field, provenance
):
    """Represent a finite ambient image inside a possibly infinite quotient.

    ``evaluator`` returns a sparse dictionary from coordinate keys to
    coefficients for one ambient element.  The ambient is finite, so the
    union of its supports is finite even when the target quotient is not.
    This is the correct interface for a local principal-divisibility problem
    such as ``R/(t^T)``: it records the finite image of the chosen ambient,
    rather than claiming that the entire quotient has a finite monomial basis.

    ``coordinate_order`` must be an explicit key function, so the emitted
    matrix is reproducible.  The helper is field-generic; callers may use QQ
    for an exact block or a finite field for a modular regression.
    """
    ambient_basis = tuple(ambient_basis)
    images = []
    coordinates = set()
    for basis_element in ambient_basis:
        image = {
            key: coefficient_field(value)
            for key, value in dict(evaluator(basis_element)).items()
            if coefficient_field(value)
        }
        images.append(image)
        coordinates.update(image)
    coordinates = tuple(sorted(coordinates, key=coordinate_order))
    coordinate_index = {key: index for index, key in enumerate(coordinates)}
    entries = {
        (coordinate_index[key], column): value
        for column, image in enumerate(images)
        for key, value in image.items()
    }
    return {
        "name": str(name),
        "matrix": matrix(
            coefficient_field, len(coordinates), len(ambient_basis), entries,
            sparse=True,
        ),
        "coordinate_keys": coordinates,
        "coordinate_basis": tuple(map(str, coordinates)),
        "provenance": str(provenance),
    }


def verify_resolved_chart_transition(
    name,
    source_ring,
    source_surface,
    target_surface,
    target_coordinates_in_source,
    source_pullbacks,
    target_pullbacks,
    source_frame=None,
    target_frame=None,
    provenance="",
):
    """Verify an explicit resolved-chart overlap map exactly.

    Both charts use the variables of ``source_ring`` for their displayed
    equations.  ``target_coordinates_in_source`` gives the target variables
    as rational functions on the source overlap, in that variable order.
    The helper checks that the target strict transform is a nonzero rational
    multiple of the source strict transform and that every supplied old-model
    pullback agrees after transport.  Optional frame functions return the
    exact transition ratio needed by a later Čech or residue block.

    No component labels, fibre types, or guessed normal forms enter this
    verification; all expressions are caller-derived from actual charts.
    """
    source_ring = source_ring
    field = source_ring.fraction_field()
    target_coordinates_in_source = tuple(
        field(value) for value in target_coordinates_in_source
    )
    if len(target_coordinates_in_source) != source_ring.ngens():
        raise ValueError("{} has the wrong number of overlap coordinates".format(name))
    source_surface = field(source_surface)
    target_surface = field(target_surface)
    transported_surface = field(target_surface(*target_coordinates_in_source))
    if not source_surface or not transported_surface:
        raise ValueError("{} has a zero strict transform on the overlap".format(name))
    surface_ratio = transported_surface/source_surface
    if not surface_ratio:
        raise ValueError("{} has a zero surface transition ratio".format(name))
    if set(source_pullbacks) != set(target_pullbacks):
        raise ValueError("{} has mismatched old-coordinate pullback labels".format(name))
    transported_pullbacks = {}
    for label in sorted(source_pullbacks):
        source_value = field(source_pullbacks[label])
        target_value = field(target_pullbacks[label])
        transported = field(target_value(*target_coordinates_in_source))
        if transported != source_value:
            raise ValueError("{} does not transport old coordinate {}".format(name, label))
        transported_pullbacks[str(label)] = source_value
    if (source_frame is None) != (target_frame is None):
        raise ValueError("{} must provide both frame functions or neither".format(name))
    frame_ratio = None
    if source_frame is not None:
        source_frame = field(source_frame)
        target_frame = field(target_frame)
        if not source_frame or not target_frame:
            raise ValueError("{} has a zero frame function".format(name))
        frame_ratio = field(target_frame(*target_coordinates_in_source))/source_frame
        if not frame_ratio:
            raise ValueError("{} has a zero frame transition ratio".format(name))
    return {
        "name": str(name),
        "surface_ratio": surface_ratio,
        "transported_pullbacks": transported_pullbacks,
        "frame_ratio": frame_ratio,
        "provenance": str(provenance),
    }


def resolved_chart_quotient_condition(
    name,
    ambient_basis,
    local_ring,
    trivialized_pullback,
    quotient_ideal,
    quotient_basis,
    provenance,
):
    """Build a condition block by reducing in an actual resolved chart.

    ``trivialized_pullback`` must map an ambient numerator into the local
    coordinate ring *after* the caller has supplied the line-bundle
    trivialization on this blow-up chart.  The explicit trivialization is
    essential: raw pullback alone can introduce exceptional factors and is
    not a divisorial condition.  ``quotient_ideal`` must be a finite-colength
    ideal in ``local_ring`` representing the required chart/module quotient.

    This adapter deliberately has no Kodaira-type argument.  A caller must
    derive the chart map and quotient from the resolution, then this function
    turns that derivation into an auditable exact matrix block.
    """
    quotient_basis = tuple(local_ring(value) for value in quotient_basis)
    groebner_basis = quotient_ideal.groebner_basis()
    if not quotient_basis:
        raise ValueError("{} has an empty quotient basis".format(name))
    quotient_dimension = quotient_ideal.vector_space_dimension()
    if quotient_dimension != len(quotient_basis):
        raise ValueError(
            "{} quotient basis has length {}, but quotient dimension is {}".format(
                name, len(quotient_basis), quotient_dimension
            )
        )
    if any(
        local_ring(value).reduce(groebner_basis) != local_ring(value)
        for value in quotient_basis
    ):
        raise ValueError("{} quotient basis is not reduced".format(name))

    def evaluator(basis_element):
        remainder = local_ring(trivialized_pullback(basis_element)).reduce(
            groebner_basis
        )
        coordinates = vector(
            QQ,
            tuple(
                remainder.monomial_coefficient(monomial)
                for monomial in quotient_basis
            ),
        )
        reconstructed = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coordinates, quotient_basis)
        )
        if remainder != reconstructed:
            raise ValueError(
                "{} quotient basis does not span the pullback remainder".format(name)
            )
        return coordinates

    return quotient_condition(
        name,
        ambient_basis,
        evaluator,
        quotient_basis,
        provenance,
    )


def resolved_chart_overlap_condition(
    name,
    ambient_basis,
    left_evaluator,
    right_evaluator,
    left_to_overlap,
    right_to_overlap,
    overlap_basis,
    provenance,
):
    """Compile equality on an explicit overlap quotient of two resolved charts.

    The evaluators return coordinates in their own finite chart quotients.
    ``left_to_overlap`` and ``right_to_overlap`` are caller-derived exact maps
    to the same finite overlap quotient.  The resulting rows impose equality
    there.  This does not infer a transition from a fibre type or a component
    label: all chart maps and trivialization ratios must be supplied by the
    resolved geometry.
    """
    ambient_basis = tuple(ambient_basis)
    overlap_basis = tuple(overlap_basis)
    left_to_overlap = matrix(QQ, left_to_overlap)
    right_to_overlap = matrix(QQ, right_to_overlap)
    if not overlap_basis:
        raise ValueError("overlap quotient basis is empty")
    if left_to_overlap.nrows() != len(overlap_basis):
        raise ValueError("left overlap map has incompatible codomain")
    if right_to_overlap.nrows() != len(overlap_basis):
        raise ValueError("right overlap map has incompatible codomain")

    def evaluator(basis_element):
        left = vector(QQ, left_evaluator(basis_element))
        right = vector(QQ, right_evaluator(basis_element))
        if len(left) != left_to_overlap.ncols():
            raise ValueError("left evaluator has incompatible quotient dimension")
        if len(right) != right_to_overlap.ncols():
            raise ValueError("right evaluator has incompatible quotient dimension")
        return left_to_overlap*left-right_to_overlap*right

    return quotient_condition(
        name,
        ambient_basis,
        evaluator,
        overlap_basis,
        provenance,
    )


def compile_resolved_conditions(
    ambient_basis, condition_blocks, complete=False, compute_kernel=True
):
    """Stack exact local conditions and return rank/codimension/kernel data.

    ``h0_certified`` is true only when a complete resolved-chart cover was
    explicitly supplied and the exact nullity is two.  Set
    ``compute_kernel=False`` when a caller needs only the exact dimension:
    this avoids asking a rational-nullspace backend to materialize an
    impractically large basis before the condition matrix has been simplified.
    When materializing, individually zero ambient columns are split off first;
    if they account for the full nullity, this returns their coordinate basis
    without a dense rational-nullspace call.  ``kernel_basis`` is ``None`` in
    dimension-only mode and cannot then be used to form a pencil.
    """
    ambient_basis = tuple(ambient_basis)
    width = len(ambient_basis)
    blocks = tuple(condition_blocks)
    for block in blocks:
        if block["matrix"].ncols() != width:
            raise ValueError("condition block {} has incompatible width".format(block["name"]))
    condition_matrix = matrix(QQ, 0, width)
    for block in blocks:
        condition_matrix = condition_matrix.stack(block["matrix"])
    rank = condition_matrix.rank()
    nullity = width-rank
    kernel_basis = None
    kernel_materialization = "not_requested"
    if compute_kernel:
        zero_columns = tuple(
            column for column in range(width)
            if all(condition_matrix[row, column] == 0
                   for row in range(condition_matrix.nrows()))
        )
        active_columns = tuple(
            column for column in range(width) if column not in zero_columns
        )
        coordinate_kernel = matrix(
            QQ, len(zero_columns), width,
            lambda row, column: QQ(1) if column == zero_columns[row] else QQ(0),
        )
        active_nullity = len(active_columns)-rank
        if active_nullity < 0:
            raise ArithmeticError("condition rank exceeds its active columns")
        if active_nullity == 0:
            kernel_basis = coordinate_kernel
            kernel_materialization = "zero_columns"
        else:
            active_matrix = condition_matrix.matrix_from_columns(active_columns)
            active_kernel = active_matrix.right_kernel()
            if active_kernel.dimension() != active_nullity:
                raise ArithmeticError("active rank and right-kernel dimensions disagree")
            lifted_kernel = matrix(
                QQ, active_kernel.dimension(), width,
                lambda row, column: (
                    active_kernel.basis_matrix()[row, active_columns.index(column)]
                    if column in active_columns else QQ(0)
                ),
            )
            kernel_basis = coordinate_kernel.stack(lifted_kernel)
            kernel_materialization = "zero_columns_plus_reduced_right_kernel"
        if kernel_basis.nrows() != nullity or condition_matrix*kernel_basis.transpose() != matrix(QQ, condition_matrix.nrows(), nullity):
            raise ArithmeticError("materialized kernel does not match the exact condition matrix")
    return {
        "ambient_dimension": width,
        "condition_rows": condition_matrix.nrows(),
        "rank": rank,
        "codimension": rank,
        "kernel_dimension": nullity,
        "kernel_basis": kernel_basis,
        "kernel_materialization": kernel_materialization,
        "complete_resolved_chart_cover": bool(complete),
        "h0_certified": bool(complete and nullity == 2),
        "condition_matrix": condition_matrix,
    }


def require_h0_two(compilation):
    """Return only after an equation-level two-dimensional pencil is certified."""
    if not compilation["complete_resolved_chart_cover"]:
        raise ValueError("resolved-chart cover is incomplete; h0(D) is not certified")
    if compilation["kernel_dimension"] != 2:
        raise ValueError("resolved condition kernel has dimension {}, not 2".format(
            compilation["kernel_dimension"]
        ))
    if compilation["kernel_basis"] is None:
        raise ValueError("kernel dimension is certified but its basis was not materialized")
    return compilation["kernel_basis"]


def certify_explicit_pencil_basis(compilation, pencil_basis):
    """Validate an independently constructed basis of a complete RR pencil.

    Exact Riemann--Roch compilers often obtain a small kernel by a structured
    quotient solve (for example a unit inversion in ``QQ[u]/(h^2)``), rather
    than by asking a dense rational linear-algebra backend to materialize a
    nullspace with very large coefficients.  This adapter binds that efficient
    construction back to the *same* complete resolved-chart condition matrix:
    it accepts a basis only when it has the certified nullity, full row rank,
    and is annihilated by every registered condition block.

    It therefore provides a reusable hand-off from local-module compilation
    to pencil elimination without treating an independently displayed pair of
    functions as evidence for ``h0(D)=2``.
    """
    if not compilation["complete_resolved_chart_cover"]:
        raise ValueError("cannot certify a pencil from an incomplete resolved-chart cover")
    if compilation["kernel_dimension"] != 2:
        raise ValueError(
            "resolved condition kernel has dimension {}, not 2".format(
                compilation["kernel_dimension"]
            )
        )
    condition_matrix = matrix(QQ, compilation["condition_matrix"])
    pencil_basis = matrix(QQ, pencil_basis)
    if pencil_basis.ncols() != condition_matrix.ncols():
        raise ValueError("pencil basis has incompatible ambient dimension")
    if pencil_basis.nrows() != compilation["kernel_dimension"]:
        raise ValueError("pencil basis has the wrong number of generators")
    if pencil_basis.rank() != pencil_basis.nrows():
        raise ValueError("pencil basis is linearly dependent")
    if condition_matrix * pencil_basis.transpose() != matrix(
        QQ, condition_matrix.nrows(), pencil_basis.nrows()
    ):
        raise ValueError("pencil basis violates a resolved-chart condition")
    return pencil_basis


def compile_elliptic_neighbor_rr_pencil(
    gram,
    raw_divisor,
    old_fiber,
    reflection_roots,
    wall_curves,
    ambient_basis,
    condition_blocks,
    complete_resolved_cover,
    expected_reflection_pairings=None,
    expected_nef_divisor=None,
    pencil_basis=None,
):
    """Compile the lattice-to-pencil half of one exact elliptic neighbour.

    This is the reusable orchestration boundary for the first five stages of
    a fibration hop.  It deliberately accepts *chart-derived* condition
    blocks rather than a Kodaira type or a component-label table:

    1. replay the ordered Weyl move from ``raw_divisor``;
    2. verify the primitive isotropic target and its declared nef walls;
    3. stack the exact resolved-chart conditions on the stated RR ambient;
    4. certify that the completed cover has a two-dimensional kernel; and
    5. bind an optional displayed basis to that same matrix.

    The caller owns the actual Weierstrass functions represented by
    ``ambient_basis`` and the trivialized chart maps in ``condition_blocks``.
    This function therefore cannot manufacture vertical corrections from the
    lattice and cannot claim that the declared walls exhaust the effective
    cone.  Its result is the strict hand-off to a degree-two chord compiler or
    to ``eliminate_cleared_weierstrass_pencil`` for a higher-degree genus-one
    presentation.
    """
    gram = matrix(ZZ, gram)
    raw_divisor = vector(ZZ, raw_divisor)
    old_fiber = vector(ZZ, old_fiber)
    if gram.nrows() != gram.ncols() or gram.nrows() != len(raw_divisor):
        raise ValueError("Neron--Severi Gram matrix and raw divisor have incompatible sizes")
    if len(old_fiber) != len(raw_divisor):
        raise ValueError("old fibre has incompatible Neron--Severi size")
    reflection_roots = tuple((str(name), vector(ZZ, root)) for name, root in reflection_roots)
    wall_curves = tuple((str(name), vector(ZZ, curve)) for name, curve in wall_curves)
    if any(len(root) != len(raw_divisor) for unused_name, root in reflection_roots):
        raise ValueError("Weyl root has incompatible Neron--Severi size")
    if any(len(curve) != len(raw_divisor) for unused_name, curve in wall_curves):
        raise ValueError("declared wall has incompatible Neron--Severi size")

    nef_divisor, reflection_record = replay_weyl_reflections(
        raw_divisor, gram, reflection_roots, expected_reflection_pairings
    )
    if expected_nef_divisor is not None and nef_divisor != vector(ZZ, expected_nef_divisor):
        raise ValueError("recorded Weyl sequence did not reach the expected nef divisor")
    preflight = isotropic_nef_preflight(gram, nef_divisor, old_fiber, wall_curves)
    if preflight["square"] != 0 or not preflight["primitive"]:
        raise ValueError("Weyl-reduced target is not primitive isotropic")
    if preflight["old_fiber_degree"] <= 0:
        raise ValueError("Weyl-reduced target has nonpositive old-fibre degree")
    if not preflight["nonnegative_on_declared_walls"]:
        raise ValueError("Weyl-reduced target is negative on a declared wall")

    compilation = compile_resolved_conditions(
        ambient_basis,
        condition_blocks,
        complete=complete_resolved_cover,
        compute_kernel=pencil_basis is None,
    )
    certified_pencil = None
    if pencil_basis is not None:
        certified_pencil = certify_explicit_pencil_basis(compilation, pencil_basis)
    elif compilation["h0_certified"]:
        certified_pencil = require_h0_two(compilation)
    return {
        "raw_divisor": raw_divisor,
        "nef_divisor": nef_divisor,
        "weyl_reflections": reflection_record,
        "preflight": preflight,
        "ambient_basis": tuple(ambient_basis),
        "rr": compilation,
        "pencil_basis": certified_pencil,
    }


def exact_neighbor_certificate_handoff(certificate, expected_old_fiber_degree=None):
    """Normalize a serialized exact-neighbour lattice certificate for equations.

    ``exact_neighbor_engine.sage`` is the source of the repository's long
    Neron--Severi neighbour chains.  This adapter makes its versioned output a
    first-class input to the equation compiler without pretending that a
    lattice certificate contains chart functions or a global nef proof.  It
    validates the fixed-component endpoint, positive old-fibre degree, and
    root/MW child data, then returns only JSON-safe data which a caller can
    pair with an independently constructed RR ambient and resolved charts.
    """
    if not isinstance(certificate, dict):
        raise TypeError("exact-neighbour certificate must be a dictionary")
    required = {"schema", "status", "proof_boundary", "divisor", "child"}
    missing = required.difference(certificate)
    if missing:
        raise ValueError("exact-neighbour certificate is missing {}".format(sorted(missing)))
    if certificate["schema"] != "elkies-k3.exact-neighbor-certificate.v1":
        raise ValueError("unrecognized exact-neighbour certificate schema")
    if certificate["status"] != "PASS_EXACT_NEIGHBOR":
        raise ValueError("exact-neighbour certificate is not passing")
    divisor = certificate["divisor"]
    child = certificate["child"]
    divisor_required = {
        "raw", "reduced", "old_fiber_degree", "nonnegative_on_supplied_curves",
        "supplied_curve_count", "fixed_component_sequence",
    }
    child_required = {
        "root_data", "root_lattice_primitive", "minimization_status", "mw_height",
    }
    if divisor_required.difference(divisor):
        raise ValueError("exact-neighbour divisor record is incomplete")
    if child_required.difference(child):
        raise ValueError("exact-neighbour child record is incomplete")
    degree = ZZ(divisor["old_fiber_degree"])
    if degree <= 0:
        raise ValueError("exact-neighbour certificate has nonpositive old-fibre degree")
    if expected_old_fiber_degree is not None and degree != ZZ(expected_old_fiber_degree):
        raise ValueError("exact-neighbour certificate old-fibre degree does not match target")
    if not bool(divisor["nonnegative_on_supplied_curves"]):
        raise ValueError("exact-neighbour certificate is negative on a supplied curve")
    root_data = tuple(ZZ(value) for value in child["root_data"])
    if len(root_data) != 3 or min(root_data) < 0 or root_data[2] <= 0:
        raise ValueError("exact-neighbour child root data are invalid")
    height = child["mw_height"]
    if child["root_lattice_primitive"]:
        if height is None:
            raise ValueError("primitive root-lattice certificate lacks MW height data")
        height = matrix(QQ, [[QQ(value) for value in row] for row in height])
        if height.nrows() != height.ncols() or height != height.transpose():
            raise ValueError("exact-neighbour child MW height is not symmetric")
    elif height is not None:
        raise ValueError("nonprimitive root-lattice certificate cannot claim saturated MW height")
    return {
        "raw_divisor": tuple(ZZ(value) for value in divisor["raw"]),
        "reduced_divisor": tuple(ZZ(value) for value in divisor["reduced"]),
        "fixed_component_sequence": tuple(
            {"name": str(item["name"]), "pairing": int(ZZ(item["pairing"]))}
            for item in divisor["fixed_component_sequence"]
        ),
        "old_fiber_degree": int(degree),
        "supplied_curve_count": int(ZZ(divisor["supplied_curve_count"])),
        "proof_boundary": str(certificate["proof_boundary"]),
        "child_root_data": tuple(int(value) for value in root_data),
        "child_root_lattice_primitive": bool(child["root_lattice_primitive"]),
        "child_minimization_status": str(child["minimization_status"]),
        "child_mw_height": height,
    }


def certify_exact_neighbor_hop(divisor, rr, child, transport, expected):
    """Check the normalized end-to-end evidence for one equation-level hop.

    The lattice engine and resolved-chart compiler deliberately have different
    inputs: the former owns the integral Neron--Severi transport, while the
    latter owns the function-space calculation.  This adapter is their strict
    hand-off.  It accepts *normalized* evidence dictionaries so individual
    surfaces need not share artifact schemas, but refuses to turn a declared
    wall calculation into an unqualified nef claim or a dimension-only RR
    calculation into a displayed pencil.

    Required dictionaries::

        divisor = {"square", "primitive", "old_fiber_degree",
                   "nef_on_declared_walls", "weyl_reflection_count"}
        rr = {"complete_resolved_chart_cover", "ambient_dimension",
              "condition_rank", "condition_codimension",
              "kernel_dimension", "h0"}
        child = {"root_lattice", "root_rank", "root_determinant",
                 "mordell_weil_rank"}
        transport = {"height_gram", "section_words",
                     "section_new_fiber_degrees"}

    ``expected`` must pin the corresponding target data, including the exact
    Gram matrix and section identities.  The returned dictionary is JSON-ready
    apart from callers' optional provenance fields.  It is intentionally a
    verifier, not a substitute for elimination or a resolved-chart cover.
    """
    required = {
        "divisor": {
            "square", "primitive", "old_fiber_degree",
            "nef_on_declared_walls", "weyl_reflection_count",
        },
        "rr": {
            "complete_resolved_chart_cover", "ambient_dimension",
            "condition_rank", "condition_codimension", "kernel_dimension", "h0",
        },
        "child": {"root_lattice", "root_rank", "root_determinant", "mordell_weil_rank"},
        "transport": {"height_gram", "section_words", "section_new_fiber_degrees"},
        "expected": {"rr", "child", "height_gram", "section_words"},
    }
    for label, values in (("divisor", divisor), ("rr", rr), ("child", child),
                          ("transport", transport), ("expected", expected)):
        missing = required[label].difference(values)
        if missing:
            raise ValueError("{} evidence is missing {}".format(label, sorted(missing)))

    if ZZ(divisor["square"]) != 0 or not bool(divisor["primitive"]):
        raise ValueError("target divisor is not primitive isotropic")
    if ZZ(divisor["old_fiber_degree"]) <= 0:
        raise ValueError("target divisor has nonpositive old-fibre degree")
    if not bool(divisor["nef_on_declared_walls"]):
        raise ValueError("target divisor is negative on a declared wall")
    if ZZ(divisor["weyl_reflection_count"]) < 0:
        raise ValueError("Weyl reflection count is negative")

    if not bool(rr["complete_resolved_chart_cover"]):
        raise ValueError("RR evidence does not declare a complete resolved-chart cover")
    if (ZZ(rr["condition_rank"]) != ZZ(rr["condition_codimension"]) or
            ZZ(rr["ambient_dimension"])-ZZ(rr["condition_rank"]) != ZZ(rr["kernel_dimension"])):
        raise ValueError("RR rank-nullity data are inconsistent")
    if ZZ(rr["h0"]) != ZZ(rr["kernel_dimension"]):
        raise ValueError("RR h0 does not equal its certified kernel dimension")
    if ZZ(rr["h0"]) != 2:
        raise ValueError("RR evidence does not certify a two-dimensional pencil")

    for key, value in expected["rr"].items():
        if key not in rr or ZZ(rr[key]) != ZZ(value):
            raise ValueError("RR target {} does not match".format(key))
    for key, value in expected["child"].items():
        if key not in child or child[key] != value:
            raise ValueError("child target {} does not match".format(key))

    height = matrix(QQ, [[QQ(entry) for entry in row] for row in transport["height_gram"]])
    expected_height = matrix(QQ, [[QQ(entry) for entry in row] for row in expected["height_gram"]])
    if height.nrows() != height.ncols() or height != height.transpose():
        raise ValueError("transported height matrix is not symmetric")
    if height != expected_height:
        raise ValueError("transported height matrix does not match the target")
    if tuple(map(str, transport["section_words"])) != tuple(map(str, expected["section_words"])):
        raise ValueError("transported section identities do not match the target")
    section_degrees = tuple(ZZ(value) for value in transport["section_new_fiber_degrees"])
    if len(section_degrees) != height.nrows():
        raise ValueError("section-degree list has the wrong length")
    if any(value <= 0 for value in section_degrees):
        raise ValueError("a transported section does not meet the new fibre positively")

    return {
        "divisor": {
            "square": int(ZZ(divisor["square"])),
            "primitive": True,
            "old_fiber_degree": int(ZZ(divisor["old_fiber_degree"])),
            "nef_on_declared_walls": True,
            "weyl_reflection_count": int(ZZ(divisor["weyl_reflection_count"])),
        },
        "rr": {
            "ambient_dimension": int(ZZ(rr["ambient_dimension"])),
            "condition_rank": int(ZZ(rr["condition_rank"])),
            "condition_codimension": int(ZZ(rr["condition_codimension"])),
            "kernel_dimension": int(ZZ(rr["kernel_dimension"])),
            "h0": int(ZZ(rr["h0"])),
        },
        "child": {
            "root_lattice": str(child["root_lattice"]),
            "root_rank": int(ZZ(child["root_rank"])),
            "root_determinant": int(ZZ(child["root_determinant"])),
            "mordell_weil_rank": int(ZZ(child["mordell_weil_rank"])),
        },
        "transport": {
            "section_words": [str(value) for value in transport["section_words"]],
            "height_gram": [[str(value) for value in row] for row in height.rows()],
            "section_new_fiber_degrees": [int(value) for value in section_degrees],
        },
    }


def certify_component_pairing_transport(
    section_sources, source_fiber_degrees, pairing_basis, resolved_pairings,
    vertical_correction, expected,
):
    """Certify transported sections against declared component pairing data.

    Component names alone do not identify a chart, an intersection point, or
    even a branch of an additive fibre.  This deliberately small adapter only
    accepts declared named source curves, their new-fibre degrees, a named
    ordered component basis, and exact pairing rows.  It compares every entry
    against caller-pinned targets and verifies the displayed
    horizontal-plus-vertical degree balance.  It does *not* promote lattice
    component data to a resolved-chart trace; callers must record that
    provenance separately.

    ``vertical_correction`` has the keys ``horizontal_degree``,
    ``correction_degree``, and ``transported_section_degree``.  It records an
    actual divisor correction; it is not inferred from a Kodaira symbol.
    """
    required_expected = {
        "section_sources", "source_fiber_degrees", "pairing_basis",
        "resolved_pairings", "vertical_correction",
    }
    missing = required_expected.difference(expected)
    if missing:
        raise ValueError("component-transport target is missing {}".format(sorted(missing)))
    if not isinstance(section_sources, dict) or not isinstance(source_fiber_degrees, dict):
        raise TypeError("component sources and fibre degrees must be dictionaries")
    if set(section_sources) != set(source_fiber_degrees):
        raise ValueError("component source labels and fibre-degree labels differ")
    if section_sources != expected["section_sources"]:
        raise ValueError("transported source-curve identities do not match the target")
    expected_degrees = expected["source_fiber_degrees"]
    if set(source_fiber_degrees) != set(expected_degrees):
        raise ValueError("component fibre-degree labels do not match the target")
    normalized_degrees = {key: ZZ(value) for key, value in source_fiber_degrees.items()}
    if normalized_degrees != {key: ZZ(value) for key, value in expected_degrees.items()}:
        raise ValueError("component fibre degrees do not match the target")
    if any(value != 1 for value in normalized_degrees.values()):
        raise ValueError("a named transported source curve is not a section")

    pairing_basis = tuple(map(str, pairing_basis))
    if pairing_basis != tuple(map(str, expected["pairing_basis"])):
        raise ValueError("component pairing basis does not match the target")
    if not isinstance(resolved_pairings, dict):
        raise TypeError("resolved pairings must be a dictionary")
    if set(resolved_pairings) != set(expected["resolved_pairings"]):
        raise ValueError("component pairing rows do not match the target")
    normalized_pairings = {}
    for label, row in resolved_pairings.items():
        row = tuple(ZZ(value) for value in row)
        expected_row = tuple(ZZ(value) for value in expected["resolved_pairings"][label])
        if len(row) != len(pairing_basis):
            raise ValueError("component pairing row {} has the wrong length".format(label))
        if row != expected_row:
            raise ValueError("component pairing row {} does not match the target".format(label))
        normalized_pairings[label] = row

    required_vertical = {"horizontal_degree", "correction_degree", "transported_section_degree"}
    if required_vertical.difference(vertical_correction):
        raise ValueError("vertical correction is missing a degree field")
    normalized_vertical = {key: ZZ(vertical_correction[key]) for key in required_vertical}
    expected_vertical = {key: ZZ(expected["vertical_correction"][key]) for key in required_vertical}
    if normalized_vertical != expected_vertical:
        raise ValueError("vertical correction degrees do not match the target")
    if (normalized_vertical["horizontal_degree"] + normalized_vertical["correction_degree"]
            != normalized_vertical["transported_section_degree"]):
        raise ValueError("horizontal and vertical correction degrees do not balance")
    if normalized_vertical["transported_section_degree"] != 1:
        raise ValueError("corrected transported divisor is not a section")

    return {
        "section_sources": {key: str(value) for key, value in section_sources.items()},
        "source_fiber_degrees": {key: int(value) for key, value in normalized_degrees.items()},
        "pairing_basis": list(pairing_basis),
        "resolved_pairings": {
            key: [int(value) for value in row] for key, row in normalized_pairings.items()
        },
        "vertical_correction": {key: int(value) for key, value in normalized_vertical.items()},
    }


def modular_condition_overlay_rank(base_matrix, overlay_matrix):
    """Rank an additional condition block only on a prior modular kernel.

    For two matrices over the same field, the stacked rank is

        rank(base) + rank(overlay restricted to ker(base)).

    The identity is useful when a resolved-chart block has a small modular
    kernel but its full rational matrix is too large to materialize densely.
    A full-column result is a good-reduction certificate for the corresponding
    characteristic-zero stacked matrix, provided the caller has separately
    checked that every displayed rational coefficient has good reduction.
    """
    base_matrix = matrix(base_matrix)
    overlay_matrix = matrix(overlay_matrix)
    if base_matrix.base_ring() != overlay_matrix.base_ring():
        raise ValueError("modular condition matrices have different base fields")
    if base_matrix.ncols() != overlay_matrix.ncols():
        raise ValueError("modular condition matrices have incompatible widths")
    kernel_basis = base_matrix.right_kernel().basis_matrix()
    restricted = overlay_matrix * kernel_basis.transpose()
    base_rank = base_matrix.rank()
    overlay_rank = restricted.rank()
    stacked_rank = base_rank+overlay_rank
    if stacked_rank > base_matrix.ncols():
        raise ArithmeticError("modular overlay rank exceeds ambient dimension")
    return {
        "ambient_dimension": base_matrix.ncols(),
        "base_rank": base_rank,
        "base_kernel_dimension": kernel_basis.nrows(),
        "overlay_rank_on_base_kernel": overlay_rank,
        "stacked_rank": stacked_rank,
        "stacked_kernel_dimension": base_matrix.ncols()-stacked_rank,
        "kernel_basis": kernel_basis,
        "restricted_overlay": restricted,
    }


def chord_tangent_slope(x_point, y_point, old_a):
    """Return the chord slope at a marked smooth point of ``y^2=x^3+a*x+b``."""
    if y_point == 0:
        raise ValueError("marked point is two-torsion; its tangent is vertical")
    return (3 * x_point**2 + old_a) / (2 * y_point)


def chord_discriminant(x_point, y_point, old_a, chord_slope):
    """Return the residual-quadratic discriminant for the marked chord.

    If ``m=(y-y(P))/(x-x(P))``, eliminating ``x`` gives a quadratic whose
    discriminant is the returned function.  This is the exact genus-one
    presentation used by a degree-two neighbour; it is independent of any
    guessed Kodaira normal form.
    """
    return (
        chord_slope**4 - 6 * x_point * chord_slope**2
        + 8 * y_point * chord_slope - 3 * x_point**2 - 4 * old_a
    )


def pencil_chord_solution(a0, b0, a1, b1, parameter):
    """Solve ``parameter=(a1+b1*m)/(a0+b0*m)`` for the marked chord ``m``."""
    denominator = parameter * b0 - b1
    if not denominator:
        raise ValueError("pencil basis is dependent in the chord direction")
    return (a1 - parameter * a0) / denominator


def pencil_value_on_marked_section(a0, b0, a1, b1, chord_value):
    """Evaluate the pencil ratio on a marked old-model section or divisor."""
    denominator = a0 + b0 * chord_value
    if not denominator:
        raise ValueError("pencil denominator vanishes identically on marked section")
    return (a1 + b1 * chord_value) / denominator


def marked_chord_value(
    x_section, y_section, x_marked, y_marked, old_a, at_marked_point=False
):
    """Evaluate the marked chord on an old-model section.

    At the marked point itself the quotient is resolved by the tangent slope;
    callers must request that branch explicitly, preventing a silent 0/0
    substitution during section transport.
    """
    if at_marked_point:
        return chord_tangent_slope(x_marked, y_marked, old_a)
    denominator = x_section - x_marked
    if not denominator:
        raise ValueError("marked chord has 0/0 value; request tangent branch explicitly")
    return (y_section - y_marked) / denominator


def rational_map_degree(value):
    """Return the degree of a rational map P1 -> P1 in a fixed base coordinate."""
    numerator = value.numerator()
    denominator = value.denominator()
    return max(numerator.degree(), denominator.degree())


def squarefree_binary_quartic(radicand, old_base_ring):
    """Extract and validate the squarefree binary quartic of a chord pencil.

    The factors of the numerator and denominator with odd valuation determine
    the binary quartic.  The omitted factor must be an *exact square* in the
    function field; checking this prevents silently replacing the child by a
    quadratic twist.  The return value is ``(quartic, square_factor)`` with
    ``radicand = square_factor^2 * quartic``.
    """
    numerator = old_base_ring(radicand.numerator())
    denominator = old_base_ring(radicand.denominator())
    numerator_factorization = numerator.factor()
    denominator_factorization = denominator.factor()
    odd_factors = tuple(
        factor
        for factorization in (numerator_factorization, denominator_factorization)
        for factor, exponent in factorization
        if exponent % 2
    )
    # Since 1/g and g agree up to the square g^-2, an odd denominator factor
    # enters the quartic with positive exponent.  Its factorization unit does
    # not disappear, however: c_num/c_den is the required scalar class.
    quartic = old_base_ring(
        numerator_factorization.unit() / denominator_factorization.unit()
    ) * old_base_ring(prod(odd_factors))
    if not quartic:
        raise ValueError("chord discriminant is a square, not a genus-one model")
    quotient = radicand / old_base_ring.fraction_field()(quartic)
    if not quotient.is_square():
        raise ValueError("odd-factor quartic lost a non-square scalar (quadratic twist)")
    square_factor = quotient.sqrt()
    if radicand != square_factor**2 * quartic:
        raise ArithmeticError("squarefree binary-quartic reconstruction failed")
    return quartic, square_factor


def binary_quartic_invariants(quartic):
    """Return classical ``I,J`` for ``a*x^4+b*x^3*z+c*x^2*z^2+d*x*z^3+e*z^4``."""
    coefficients = list(quartic.list()) + [quartic.base_ring()(0)] * 5
    e, d, c, b, a = coefficients[:5]
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e + 9 * b * c * d - 27 * a * d**2
        - 27 * b**2 * e - 2 * c**3
    )
    return invariant_i, invariant_j


def binary_quartic_jacobian_coefficients(quartic):
    """Return ``(A,B,Delta)`` for the short Jacobian ``Y^2=X^3+A*X+B``."""
    invariant_i, invariant_j = binary_quartic_invariants(quartic)
    coefficient_a = -27 * invariant_i
    coefficient_b = -27 * invariant_j
    discriminant = 4 * invariant_i**3 - invariant_j**2
    if not discriminant:
        raise ValueError("binary quartic is singular")
    return coefficient_a, coefficient_b, discriminant


def binary_quartic_covariants(quartic):
    """Return the classical ``(H,G)`` covariants of an exact binary quartic.

    For ``f(x,z)`` associated to the univariate quartic, set
    ``H=(f_xx*f_zz-f_xz^2)/3`` and ``G=f_x*H_z-f_z*H_x``.  The returned
    binary polynomial ring and covariants satisfy the checked identity used
    by the Jacobian map.  Keeping this in the core lets section transport use
    an actual resolved divisor point, rather than search the child equation.
    """
    base_ring = quartic.base_ring()
    binary_ring = PolynomialRing(base_ring, names=("x", "z"))
    x, z = binary_ring.gens()
    binary_quartic = sum(
        binary_ring(quartic[index]) * x**index * z**(4-index)
        for index in range(5)
    )
    hessian = (
        binary_quartic.derivative(x, 2) * binary_quartic.derivative(z, 2)
        - binary_quartic.derivative(x).derivative(z)**2
    ) / 3
    jacobian_covariant = (
        binary_quartic.derivative(x) * hessian.derivative(z)
        - binary_quartic.derivative(z) * hessian.derivative(x)
    )
    invariant_i, invariant_j = binary_quartic_invariants(quartic)
    if jacobian_covariant**2 != (
        -QQ(16)/3 * hessian**3
        + 256 * invariant_i * hessian * binary_quartic**2
        - QQ(1024)/3 * invariant_j * binary_quartic**3
    ):
        raise ArithmeticError("binary-quartic covariant identity failed")
    return binary_ring, binary_quartic, hessian, jacobian_covariant


def transport_binary_quartic_point_to_jacobian(
    quartic, x_value, z_value, w_value, minimalizing_unit=1
):
    """Transport an exact point ``[x:z:w]`` on ``w^2=f(x,z)`` to its Jacobian.

    The map is the covariant map
    ``X=-3H/(4f)``, ``Y=9G*w/(32f^2)`` to the standard short Jacobian.  An
    explicitly supplied fourth/sixth-power ``minimalizing_unit`` then yields
    coordinates on a caller's minimal child model.  The function verifies the
    source point and both Weierstrass equations exactly.
    """
    base_field = quartic.base_ring().fraction_field()
    binary_ring, binary_quartic, hessian, jacobian_covariant = binary_quartic_covariants(
        quartic
    )
    x_value, z_value, w_value = tuple(
        base_field(value) for value in (x_value, z_value, w_value)
    )
    quartic_value = base_field(binary_quartic(x=x_value, z=z_value))
    if not quartic_value or w_value**2 != quartic_value:
        raise ValueError("binary-quartic point does not satisfy w^2=f or lies at f=0")
    hessian_value = base_field(hessian(x=x_value, z=z_value))
    covariant_value = base_field(jacobian_covariant(x=x_value, z=z_value))
    raw_x = -QQ(3)/4 * hessian_value/quartic_value
    raw_y = QQ(9)/32 * covariant_value*w_value/quartic_value**2
    coefficient_a, coefficient_b, unused_discriminant = binary_quartic_jacobian_coefficients(
        quartic
    )
    if raw_y**2 != raw_x**3 + coefficient_a*raw_x + coefficient_b:
        raise ArithmeticError("binary-quartic covariant map missed the standard Jacobian")
    minimalizing_unit = base_field(minimalizing_unit)
    child_x = minimalizing_unit**2 * raw_x
    child_y = minimalizing_unit**3 * raw_y
    if child_y**2 != child_x**3 + (minimalizing_unit**4*coefficient_a)*child_x + minimalizing_unit**6*coefficient_b:
        raise ArithmeticError("minimalized binary-quartic transport missed its child equation")
    return {
        "standard_x": raw_x,
        "standard_y": raw_y,
        "standard_a": coefficient_a,
        "standard_b": coefficient_b,
        "child_x": child_x,
        "child_y": child_y,
        "child_a": minimalizing_unit**4*coefficient_a,
        "child_b": minimalizing_unit**6*coefficient_b,
    }


def compile_degree_two_chord_hop(
    old_base_ring,
    parameter,
    a0,
    b0,
    a1,
    b1,
    x_marked,
    y_marked,
    old_a,
    old_b=None,
    marked_chords=(),
):
    """Compile a degree-two chord pencil to its exact quartic and Jacobian.

    The two pencil generators are represented in the marked-chord frame as
    ``s_i=a_i+b_i*m``.  Their ratio is the new base parameter.  This routine
    performs the complete algebraic hand-off shared by every degree-two
    neighbour:

    * validate the marked point on the old short Weierstrass fibre when
      ``old_b`` is supplied;
    * solve the pencil relation for ``m`` over the new parameter field;
    * form the exact chord discriminant and remove *only* a certified square;
    * return the binary quartic and its short Jacobian invariants; and
    * evaluate the new parameter on explicitly transported marked chords.

    ``old_base_ring`` is the polynomial ring in the old base coordinate over
    the new parameter field.  The returned values remain exact elements of
    its fraction field, so minimization and local fibre classification stay
    surface-specific downstream steps.  ``marked_chords`` is an iterable of
    ``(name, chord_value)`` pairs; no old-model section is rediscovered.
    """
    old_base_field = old_base_ring.fraction_field()
    parameter = old_base_field(parameter)
    a0, b0, a1, b1, x_marked, y_marked, old_a = tuple(
        old_base_field(value)
        for value in (a0, b0, a1, b1, x_marked, y_marked, old_a)
    )
    if old_b is not None:
        old_b = old_base_field(old_b)
        if y_marked**2 != x_marked**3 + old_a*x_marked + old_b:
            raise ValueError("marked point does not lie on the old Weierstrass fibre")
    determinant = a0*b1-a1*b0
    if not determinant:
        raise ValueError("degree-two pencil generators are dependent in the chord frame")
    chord = pencil_chord_solution(a0, b0, a1, b1, parameter)
    radicand = chord_discriminant(x_marked, y_marked, old_a, chord)
    quartic, square_factor = squarefree_binary_quartic(radicand, old_base_ring)
    coefficient_a, coefficient_b, discriminant = binary_quartic_jacobian_coefficients(quartic)

    transported_values = {}
    for name, chord_value in marked_chords:
        name = str(name)
        if name in transported_values:
            raise ValueError("duplicate transported marked divisor {}".format(name))
        transported_values[name] = pencil_value_on_marked_section(
            a0, b0, a1, b1, old_base_field(chord_value)
        )
    return {
        "parameter": parameter,
        "pencil_determinant": determinant,
        "chord": chord,
        "radicand": radicand,
        "binary_quartic": quartic,
        "square_factor": square_factor,
        "jacobian_a": coefficient_a,
        "jacobian_b": coefficient_b,
        "jacobian_discriminant": discriminant,
        "transported_parameter_values": transported_values,
    }


def compile_resolved_degree_two_chord_hop(
    compilation,
    pencil_basis,
    chord_expansions,
    old_base_ring,
    parameter,
    x_marked,
    y_marked,
    old_a,
    old_b=None,
    marked_chords=(),
):
    """Bind a complete resolved RR pencil to the degree-two chord compiler.

    ``chord_expansions`` gives, in the exact ambient-basis order, pairs
    ``(a,b)`` representing each ambient function as ``a+b*m``.  The supplied
    two-row ``pencil_basis`` is first checked against the *same* complete
    resolved condition matrix.  Applying it to these pairs produces the two
    chord-frame generators, so the new parameter and binary quartic cannot be
    detached from the vertical-condition calculation.

    This is the equation-level hand-off for a degree-two neighbour.  For
    higher old-fibre degree the resolved RR interface remains applicable, but
    a surface-specific genus-one presentation/Jacobian adapter is required.
    """
    certified_basis = certify_explicit_pencil_basis(compilation, pencil_basis)
    if len(chord_expansions) != certified_basis.ncols():
        raise ValueError("chord expansion list has incompatible ambient dimension")
    old_base_field = old_base_ring.fraction_field()
    expansions = tuple(
        (old_base_field(entry[0]), old_base_field(entry[1]))
        for entry in chord_expansions
    )
    coefficients = []
    for row in range(certified_basis.nrows()):
        coefficients.append(tuple(
            sum(
                old_base_field(certified_basis[row, column]) * expansions[column][coordinate]
                for column in range(certified_basis.ncols())
            )
            for coordinate in range(2)
        ))
    (a0, b0), (a1, b1) = coefficients
    conversion = compile_degree_two_chord_hop(
        old_base_ring, parameter, a0, b0, a1, b1,
        x_marked, y_marked, old_a, old_b, marked_chords,
    )
    return {
        "pencil_basis": certified_basis,
        "chord_coefficients": ((a0, b0), (a1, b1)),
        "conversion": conversion,
    }


def eliminate_cleared_weierstrass_pencil(
    ambient_ring,
    old_weierstrass_relation,
    pencil_zero,
    pencil_one,
    new_parameter,
    eliminate_variables,
    saturate_by=(),
):
    """Eliminate selected old variables from an exact cleared pencil ideal.

    This is the degree-independent algebraic first stage after a certified RR
    pencil.  All inputs must already lie in one polynomial ring on a declared
    affine chart: ``pencil_zero`` and ``pencil_one`` are common-denominator
    representatives of the two sections, and the new-pencil equation is
    ``pencil_one-new_parameter*pencil_zero``.  The caller explicitly chooses
    which *old* variables to eliminate; retained variables describe the raw
    genus-one presentation over the new parameter.

    Saturation is permitted only through the explicit ``saturate_by`` list of
    declared chart denominators.  The unsaturated ideal, product, and
    saturated ideal are all returned.  Projective closure, Jacobian, and
    minimization remain separate certificate stages.  Returning the full
    elimination ideal prevents a caller from silently discarding an
    additional relation.
    """
    old_weierstrass_relation = ambient_ring(old_weierstrass_relation)
    pencil_zero = ambient_ring(pencil_zero)
    pencil_one = ambient_ring(pencil_one)
    new_parameter = ambient_ring(new_parameter)
    variables = tuple(ambient_ring(value) for value in eliminate_variables)
    if not old_weierstrass_relation:
        raise ValueError("old Weierstrass relation is zero")
    if not pencil_zero or not pencil_one:
        raise ValueError("cleared pencil generators must both be nonzero")
    if not variables or len(set(variables)) != len(variables):
        raise ValueError("elimination variables must be a nonempty distinct list")
    if any(value not in ambient_ring.gens() for value in variables):
        raise ValueError("an elimination variable is not a generator of the ambient ring")
    pencil_relation = pencil_one-new_parameter*pencil_zero
    source_ideal = ambient_ring.ideal((old_weierstrass_relation, pencil_relation))
    saturation_factors = tuple(ambient_ring(value) for value in saturate_by)
    if any(not value for value in saturation_factors):
        raise ValueError("declared saturation factor is zero")
    saturation_product = ambient_ring(prod(saturation_factors))
    saturated_ideal = source_ideal
    if saturation_factors:
        saturated_ideal = source_ideal.saturation(saturation_product)[0]
    elimination = saturated_ideal.elimination_ideal(list(variables))
    if elimination.is_one():
        raise ArithmeticError("cleared pencil ideal has empty elimination locus")
    return {
        "pencil_relation": pencil_relation,
        "source_ideal": source_ideal,
        "saturation_factors": saturation_factors,
        "saturation_product": saturation_product,
        "saturated_ideal": saturated_ideal,
        "eliminated_variables": variables,
        "elimination_ideal": elimination,
        "relations": tuple(elimination.gens()),
    }


def compile_resolved_genus_one_elimination(
    compilation,
    pencil_basis,
    ambient_representatives,
    ambient_ring,
    old_weierstrass_relation,
    new_parameter,
    eliminate_variables,
    saturate_by=(),
):
    """Eliminate a higher-degree neighbour from its certified RR pencil.

    This is the degree-independent counterpart of
    ``compile_resolved_degree_two_chord_hop``.  Each entry of
    ``ambient_representatives`` is a caller-derived, common-denominator
    representative of the corresponding RR basis element in ``ambient_ring``.
    The two certified kernel rows are combined there before any elimination.
    Thus a raw genus-one relation cannot be substituted from a separately
    displayed pencil after the resolved-chart calculation has certified a
    different one.

    The result is an exact elimination ideal, not a claim that a selected
    generator is a plane cubic, that it has a rational point, or that its
    Jacobian is minimal.  Those operations must use the appropriate explicit
    genus-one and section-transport adapters downstream.
    """
    certified_basis = certify_explicit_pencil_basis(compilation, pencil_basis)
    representatives = tuple(ambient_ring(value) for value in ambient_representatives)
    if len(representatives) != certified_basis.ncols():
        raise ValueError("ambient representative list has incompatible dimension")
    if any(not value for value in representatives):
        raise ValueError("an ambient representative is zero")
    pencil_zero, pencil_one = tuple(
        sum(
            ambient_ring(certified_basis[row, column])*representatives[column]
            for column in range(certified_basis.ncols())
        )
        for row in range(2)
    )
    if not pencil_zero or not pencil_one:
        raise ArithmeticError("certified pencil has a zero cleared representative")
    elimination = eliminate_cleared_weierstrass_pencil(
        ambient_ring,
        old_weierstrass_relation,
        pencil_zero,
        pencil_one,
        new_parameter,
        eliminate_variables,
        saturate_by=saturate_by,
    )
    return {
        "pencil_basis": certified_basis,
        "pencil_zero": pencil_zero,
        "pencil_one": pencil_one,
        "elimination": elimination,
    }


def pointed_plane_cubic_to_weierstrass(cubic, tangent_frame):
    """Normalize a pointed plane cubic using an explicit non-flex tangent frame.

    ``cubic`` is a homogeneous ternary cubic over a field.  ``tangent_frame``
    is the ordered triple ``(P, Q, R)`` of the point and its two successive
    third tangent intersections.  Requiring that frame prevents this routine
    from invoking a flex search or guessing a rational origin.  It returns a
    Sage elliptic curve together with the rational source-to-Weierstrass map.

    The map may have a base point at a tangent-frame point.  A resolved-chart
    caller must evaluate its limit there; this routine records the map but
    does not silently substitute a 0/0 projective triple.  Minimalization and
    identification with a target Weierstrass model are separate steps.
    """
    ring = cubic.parent()
    if ring.ngens() != 3 or cubic.total_degree() != 3:
        raise ValueError("plane-cubic normalizer requires a homogeneous ternary cubic")
    field = ring.base_ring()
    x, y, z = ring.gens()
    if len(tangent_frame) != 3:
        raise ValueError("non-flex tangent frame must contain three points")
    points = tuple(tuple(field(value) for value in point) for point in tangent_frame)
    if any(len(point) != 3 or not any(point) for point in points):
        raise ValueError("tangent-frame point has invalid projective coordinates")
    if any(cubic(*point) for point in points):
        raise ValueError("tangent-frame point does not lie on the plane cubic")
    transform = matrix(field, points).transpose()
    if not transform.det():
        raise ValueError("non-flex tangent frame is projectively dependent")
    transformed = ring(transform.act_on_polynomial(cubic))
    first_substitution = (x*x, y*z, x*z)
    first_denominator = x*x*z
    numerator = ring(transformed(first_substitution))
    quartic, remainder = numerator.quo_rem(first_denominator)
    if remainder:
        raise ArithmeticError("tangent-frame quadratic transform has residual base factor")
    leading_x = field(quartic.coefficient(x**3))
    leading_y = field(quartic.coefficient(y*y*z))
    if not leading_x or not leading_y:
        raise ArithmeticError("non-flex tangent frame gives degenerate Weierstrass normalization")
    second_substitution = (-x, y/leading_y, leading_x*leading_y*z)
    normalized = ring(quartic(second_substitution) / leading_x)
    curve = EllipticCurve(normalized([x, y, 1]))
    if curve.discriminant() == 0:
        raise ArithmeticError("plane cubic normalized to a singular Weierstrass curve")

    transformed_coordinates = tuple(
        sum(transform[row, column] * (x, y, z)[column] for column in range(3))
        for row in range(3)
    )
    after_quadratic = tuple(coordinate(first_substitution) for coordinate in transformed_coordinates)
    source_to_weierstrass = tuple(
        ring(coordinate(second_substitution)) for coordinate in after_quadratic
    )
    return {
        "curve": curve,
        "normalized_cubic": normalized,
        "tangent_frame": points,
        "source_to_weierstrass": source_to_weierstrass,
    }


def finite_minimize_short_weierstrass(base_ring, coefficient_a, coefficient_b):
    """Finite-place minimize ``Y^2=X^3+A*X+B`` over a rational base.

    The return value records the exact product ``u`` for the change
    ``X=u^2 X_min, Y=u^3 Y_min`` and the valuation data at every irreducible
    finite base factor.  It does not claim global minimality at infinity:
    that chart's raw and normalized valuations are reported separately because
    a K3 (or other elliptic-surface) line-bundle scaling is extra geometry.
    """
    base_field = base_ring.fraction_field()
    coefficient_a = base_field(coefficient_a)
    coefficient_b = base_field(coefficient_b)
    discriminant = -16*(4*coefficient_a**3 + 27*coefficient_b**2)
    if not discriminant:
        raise ValueError("short Weierstrass model is singular")

    factors = set()
    for value in (coefficient_a, coefficient_b, discriminant):
        numerator = base_ring(value.numerator())
        denominator = base_ring(value.denominator())
        factors.update(factor for factor, unused in numerator.factor())
        factors.update(factor for factor, unused in denominator.factor())
    scaling_unit = base_field(1)
    finite_data = []
    for factor in sorted(factors, key=str):
        raw_a = int(coefficient_a.valuation(factor))
        raw_b = int(coefficient_b.valuation(factor))
        raw_delta = int(discriminant.valuation(factor))
        scaling = min(raw_a // 4, raw_b // 6)
        scaling_unit *= base_field(factor)**(-scaling)
        finite_data.append({
            "factor": factor,
            "raw_orders": (raw_a, raw_b, raw_delta),
            "scaling": scaling,
        })
    minimal_a = coefficient_a*scaling_unit**4
    minimal_b = coefficient_b*scaling_unit**6
    minimal_delta = discriminant*scaling_unit**12
    if any(value.denominator() != 1 for value in (minimal_a, minimal_b, minimal_delta)):
        raise ArithmeticError("finite minimization did not clear all base denominators")
    minimal_a = base_ring(minimal_a.numerator())
    minimal_b = base_ring(minimal_b.numerator())
    minimal_delta = base_ring(minimal_delta.numerator())
    for item in finite_data:
        factor = item["factor"]
        item["minimal_orders"] = (
            int(minimal_a.valuation(factor)),
            int(minimal_b.valuation(factor)),
            int(minimal_delta.valuation(factor)),
        )
        if min(item["minimal_orders"][0] // 4, item["minimal_orders"][1] // 6) != 0:
            raise ArithmeticError("finite minimization left a fourth/sixth-power factor")
    infinity_raw = tuple(-value.degree() for value in (minimal_a, minimal_b, minimal_delta))
    infinity_scaling = min(infinity_raw[0] // 4, infinity_raw[1] // 6)
    infinity_normalized = tuple(
        infinity_raw[index] - (4, 6, 12)[index]*infinity_scaling
        for index in range(3)
    )
    return {
        "minimal_a": minimal_a,
        "minimal_b": minimal_b,
        "minimal_discriminant": minimal_delta,
        "scaling_unit": scaling_unit,
        "finite_places": tuple(finite_data),
        "infinity": {
            "raw_orders": infinity_raw,
            "scaling": infinity_scaling,
            "normalized_orders": infinity_normalized,
        },
    }


def kodaira_data_from_short_orders(order_a, order_b, order_discriminant):
    """Classify a minimal characteristic-zero short Weierstrass valuation triple.

    The result is ``(root_rank, euler_number, root_determinant, symbol)``.
    The input must already be minimal at the place and have positive
    discriminant order.  This is a local classifier only; component labels
    and section intersections still require the actual resolved fibre chart.
    """
    order_a, order_b, order_discriminant = tuple(
        ZZ(value) for value in (order_a, order_b, order_discriminant)
    )
    if min(order_a, order_b, order_discriminant) < 0 or order_discriminant == 0:
        raise ValueError("Kodaira classifier requires a minimal singular valuation triple")
    if order_a == 0 or order_b == 0:
        return int(order_discriminant-1), int(order_discriminant), int(order_discriminant), "I{}".format(order_discriminant)
    if order_discriminant == 2:
        return 0, 2, 1, "II"
    if order_discriminant == 3:
        return 1, 3, 2, "III"
    if order_discriminant == 4:
        return 2, 4, 3, "IV"
    if order_discriminant == 6 and order_a >= 2 and order_b >= 3:
        return 4, 6, 4, "I0*"
    if order_discriminant >= 7 and order_a == 2 and order_b == 3:
        multiplicative_index = order_discriminant-6
        return int(multiplicative_index+4), int(multiplicative_index+6), 4, "I{}*".format(multiplicative_index)
    if order_discriminant == 8:
        return 6, 8, 3, "IV*"
    if order_discriminant == 9:
        return 7, 9, 2, "III*"
    if order_discriminant == 10:
        return 8, 10, 1, "II*"
    raise ValueError("unrecognized minimal short-Weierstrass valuation triple")


def classify_finite_short_weierstrass_fibres(base_ring, coefficient_a, coefficient_b):
    """Aggregate exact *finite* Kodaira data of a short Weierstrass model.

    This is deliberately a two-stage operation: first perform the exact
    fourth/sixth-power finite minimization, then classify every finite place
    with positive minimal discriminant valuation.  The result contains no
    assertion about the fibre at infinity.  That chart has to be supplied
    with its own line-bundle trivialization and resolution, so its valuation
    record is returned only as an explicit boundary diagnostic.

    In particular, the Kodaira symbols below certify root ranks and
    determinants.  They do not assign a marked section to a fibre component;
    that still requires a resolved-chart transport calculation.
    """
    minimization = finite_minimize_short_weierstrass(
        base_ring, coefficient_a, coefficient_b
    )
    finite_fibres = []
    root_rank = ZZ(0)
    euler_number = ZZ(0)
    root_determinant = ZZ(1)
    for place in minimization["finite_places"]:
        orders = tuple(ZZ(value) for value in place["minimal_orders"])
        if orders[2] == 0:
            continue
        rank, euler, determinant, symbol = kodaira_data_from_short_orders(*orders)
        degree = ZZ(place["factor"].degree())
        finite_fibres.append({
            "factor": place["factor"],
            "degree": int(degree),
            "raw_orders": tuple(int(value) for value in place["raw_orders"]),
            "scaling": int(place["scaling"]),
            "minimal_orders": tuple(int(value) for value in orders),
            "kodaira": symbol,
            "root_rank": int(rank),
            "euler_number": int(euler),
            "root_determinant": int(determinant),
        })
        root_rank += degree*rank
        euler_number += degree*euler
        root_determinant *= ZZ(determinant)**degree
    return {
        "finite_minimization": minimization,
        "finite_fibres": tuple(finite_fibres),
        "finite_root_rank": int(root_rank),
        "finite_euler_number": int(euler_number),
        "finite_root_determinant": int(root_determinant),
        "infinity_boundary": minimization["infinity"],
    }


def certify_shioda_tate_discriminant(
    root_determinant, height_gram, torsion_order=1, expected_ns_discriminant=None
):
    """Certify the absolute Neron--Severi discriminant from Shioda--Tate data.

    For a fibration with trivial-lattice root determinant ``d_R``, Mordell--
    Weil height Gram ``H``, and torsion/glue index ``n``, the absolute
    discriminant is ``d_R*det(H)/n^2``.  The torsion order is kept explicit:
    setting it to one is a certified torsion-free assertion, not a hidden
    convention.  This arithmetic check complements, but never replaces,
    resolved component/glue maps.
    """
    root_determinant = ZZ(root_determinant)
    torsion_order = ZZ(torsion_order)
    if root_determinant <= 0 or torsion_order <= 0:
        raise ValueError("root determinant and torsion/glue index must be positive")
    height = matrix(QQ, [[QQ(value) for value in row] for row in height_gram])
    if height.nrows() != height.ncols() or height != height.transpose():
        raise ValueError("height Gram must be square and symmetric")
    height_determinant = height.det()
    if height_determinant <= 0:
        raise ValueError("height Gram is not positive determinant")
    discriminant = QQ(root_determinant)*height_determinant/(torsion_order**2)
    if discriminant.denominator() != 1:
        raise ArithmeticError("Shioda--Tate discriminant is not integral")
    discriminant = ZZ(discriminant)
    if expected_ns_discriminant is not None and discriminant != ZZ(expected_ns_discriminant):
        raise ValueError("Neron--Severi discriminant does not match its target")
    return {
        "root_determinant": root_determinant,
        "height_determinant": height_determinant,
        "torsion_glue_index": torsion_order,
        "absolute_ns_discriminant": discriminant,
    }
