#!/usr/bin/env sage
"""Reusable exact lattice engine for degree-two elliptic neighbors.

The input is an integral Neron--Severi Gram matrix, an isotropic divisor,
the old fiber, and named effective (-2)-curves (sections and components).
The engine performs only exact integer/rational arithmetic:

* remove fixed components by deterministic reflection/reduction;
* split the resulting primitive isotropic vector together with a hyperbolic
  mate, producing an integral child frame and transport; and
* make the child convenient for the next neighbor by choosing deterministic
  simple roots and LLL-reducing its Mordell--Weil quotient.

``nonnegative_on_supplied_curves`` means exactly that.  A full nef proof still
needs the appropriate section/multisection argument for the given chamber;
this module never promotes a finite wall list to such a proof.
"""

from sage.all import (
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
    xgcd,
)

import hashlib
import json
from pathlib import Path

from sage.version import version as SAGE_VERSION


CERTIFICATE_SCHEMA = "elkies-k3.exact-neighbor-certificate.v1"
INPUT_SCHEMA = "elkies-k3.exact-neighbor-input.v1"


def matrix_rows(value):
    """Serialize an integral Sage matrix with no implementation-dependent form."""
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def rational_rows(value):
    """Serialize rational matrices exactly as strings."""
    return [[str(entry) for entry in row] for row in matrix(QQ, value).rows()]


def vector_entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def matrix_digest(value):
    """SHA-256 of a canonical integral row-major matrix representation."""
    payload = ";".join(
        ",".join(str(entry) for entry in row)
        for row in matrix(ZZ, value).rows()
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def json_digest(value):
    """SHA-256 of a canonical JSON payload."""
    text = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def intersection(left, right, ns):
    """Return the exact intersection of two Neron--Severi row vectors."""
    return ZZ(vector(ZZ, left) * ns * vector(ZZ, right))


def reduce_fixed_components(divisor, curves, ns, max_steps=10000):
    """Remove fixed supplied (-2)-curves in the supplied deterministic order.

    Each item of ``curves`` is ``(name, class)``.  When ``D.C<0``, the update
    ``D <- D + (D.C) C`` removes ``-D.C`` copies of ``C``.  The return value is
    ``(reduced_divisor, sequence)`` where each sequence entry is
    ``(name, negative_pairing)``.
    """
    divisor = vector(ZZ, divisor)
    curves = tuple((str(name), vector(ZZ, curve)) for name, curve in curves)
    square = intersection(divisor, divisor, ns)
    for name, curve in curves:
        if intersection(curve, curve, ns) != -2:
            raise ValueError(f"wall {name} is not a (-2)-curve")

    sequence = []
    for _ in range(max_steps):
        for name, curve in curves:
            pairing = intersection(divisor, curve, ns)
            if pairing < 0:
                divisor += pairing * curve
                if intersection(divisor, divisor, ns) != square:
                    raise AssertionError("fixed-component reduction changed D^2")
                sequence.append((name, pairing))
                break
        else:
            return divisor, tuple(sequence)
    raise RuntimeError("fixed-component reduction did not terminate")


def primitive_hyperbolic_split(ns, fiber):
    """Split a primitive isotropic vector and return child frame and transport."""
    ns = matrix(ZZ, ns)
    fiber = vector(ZZ, fiber)
    if ns.nrows() != ns.ncols() or ns != ns.transpose():
        raise ValueError("Neron--Severi Gram matrix must be square symmetric")
    if len(fiber) != ns.nrows() or intersection(fiber, fiber, ns) != 0:
        raise ValueError("fiber must be an isotropic Neron--Severi vector")

    pairings = tuple(ns * fiber)
    current = ZZ(0)
    mate_entries = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, old_scale, new_scale = xgcd(current, ZZ(value))
        mate_entries = [old_scale * entry for entry in mate_entries]
        mate_entries[index] += new_scale
        current = divisor
    if abs(current) != 1:
        raise ValueError("fiber is not primitive in the Neron--Severi lattice")
    if current == -1:
        mate_entries = [-entry for entry in mate_entries]
    mate = vector(ZZ, mate_entries)
    if intersection(fiber, mate, ns) != 1:
        raise AssertionError("Bezout mate has incorrect fiber pairing")
    mate_square = intersection(mate, mate, ns)
    if mate_square % 2:
        raise AssertionError("even Neron--Severi lattice expected")
    mate -= (mate_square // 2) * fiber
    if intersection(mate, mate, ns) != 0 or intersection(fiber, mate, ns) != 1:
        raise AssertionError("failed to normalize hyperbolic mate")

    complement = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in complement.rows()]
    )
    hyperbolic = matrix(ZZ, ((0, 1), (1, 0)))
    if abs(transport.det()) != 1:
        raise AssertionError("neighbor transport is not unimodular")
    if transport * ns * transport.transpose() != block_diagonal_matrix(
        hyperbolic, -child
    ):
        raise AssertionError("neighbor transport does not split U")
    if not child.is_positive_definite():
        raise AssertionError("child frame is not positive definite")
    return {
        "fiber": fiber,
        "mate": mate,
        "child_frame": child,
        "transport": transport,
    }


def transport_parent_vector_to_child(neighbor, parent_vector):
    """Express a parent Neron--Severi vector in the child's U-plus-frame basis."""
    transport = matrix(ZZ, neighbor["transport"])
    parent_vector = vector(ZZ, parent_vector)
    if len(parent_vector) != transport.ncols():
        raise ValueError("parent vector has the wrong Neron--Severi rank")
    return vector(ZZ, parent_vector * transport.inverse())


def lift_child_frame_vector(neighbor, child_vector):
    """Lift a positive-child-frame vector to the parent Neron--Severi lattice."""
    transport = matrix(ZZ, neighbor["transport"])
    child_vector = vector(ZZ, child_vector)
    if len(child_vector) != transport.nrows() - 2:
        raise ValueError("child vector has the wrong frame rank")
    return vector(ZZ, vector(ZZ, [0, 0] + list(child_vector)) * transport)


def lift_child_component_data(neighbor):
    """Lift deterministic child simple roots to named parent divisor classes."""
    simple = matrix(ZZ, neighbor["child_simple_roots"])
    return tuple(
        (
            f"child_R{index:02d}",
            lift_child_frame_vector(neighbor, root),
        )
        for index, root in enumerate(simple.rows(), 1)
    )


def transport_marked_parent_vectors(neighbor, markings):
    """Transport named parent divisors to full child Neron--Severi coordinates."""
    return {
        str(name): transport_parent_vector_to_child(neighbor, divisor)
        for name, divisor in markings.items()
    }


def roots_and_data(frame):
    """Return all norm-two roots, their saturated span, and its invariants."""
    frame = matrix(ZZ, frame)
    result = pari(frame).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, frame.ncols()), (0, 0, 1)
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2]).columns()
    ]
    roots = tuple(sorted(half + [-root for root in half], key=tuple))
    if len(roots) != count:
        raise AssertionError("PARI root count disagrees with the root list")
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    gram = basis * frame * basis.transpose()
    return roots, basis, (basis.rank(), count, abs(ZZ(gram.det())))


def deterministic_simple_roots(frame):
    """Choose lexicographically positive simple roots of the root subsystem."""
    roots, _, data = roots_and_data(frame)
    rank = data[0]
    if rank == 0:
        return matrix(ZZ, 0, matrix(ZZ, frame).ncols()), (), ()
    positive = tuple(
        root for root in roots
        if next(value for value in root if value) > 0
    )
    positive_set = {tuple(root) for root in positive}
    simple = tuple(
        root for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    )
    simple_matrix = matrix(ZZ, [list(root) for root in simple])
    cartan = simple_matrix * matrix(ZZ, frame) * simple_matrix.transpose()
    if simple_matrix.nrows() != rank or simple_matrix.rank() != rank:
        raise AssertionError("failed to recover a simple root basis")
    if any(cartan[i, i] != 2 for i in range(rank)):
        raise AssertionError("simple roots do not have norm two")
    if any(
        cartan[i, j] not in (0, -1)
        for i in range(rank) for j in range(rank) if i != j
    ):
        raise AssertionError("simple roots are not ADE")
    return simple_matrix, positive, tuple(tuple(row) for row in cartan.rows())


def component_walls(frame, old_fiber, include_zero=None):
    """Build named simple and affine old-fiber components from a child frame.

    ``old_fiber`` and the optional ``include_zero`` are Neron--Severi vectors;
    the positive frame coordinates are embedded after their first two entries.
    """
    frame = matrix(ZZ, frame)
    old_fiber = vector(ZZ, old_fiber)
    simple, positive, cartan_rows = deterministic_simple_roots(frame)
    walls = []
    if include_zero is not None:
        walls.append(("O", vector(ZZ, include_zero)))
    walls.extend(
        (f"R{index:02d}", vector(ZZ, [0, 0] + list(root)))
        for index, root in enumerate(simple.rows(), 1)
    )

    cartan = matrix(ZZ, cartan_rows)
    unseen = set(range(cartan.nrows()))
    component_number = 0
    inverse_simple = simple.pseudoinverse() if simple.nrows() else None
    while unseen:
        component_number += 1
        pending = [min(unseen)]
        unseen.remove(pending[0])
        component = []
        while pending:
            index = pending.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other] != 0]
            for other in adjacent:
                unseen.remove(other)
                pending.append(other)
        component = tuple(sorted(component))
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root) * inverse_simple
            if not all(value in ZZ and value >= 0 for value in coordinates):
                continue
            support = tuple(index for index, value in enumerate(coordinates) if value)
            if support and all(index in component for index in support):
                candidates.append((sum(coordinates), tuple(root), root))
        if not candidates:
            raise AssertionError("could not find a highest root")
        _, _, highest = max(candidates)
        walls.append(
            (
                f"Theta0_{component_number}",
                old_fiber - vector(ZZ, [0, 0] + list(highest)),
            )
        )
    return tuple(walls), simple


def minimize_child_frame(frame):
    """Root-adapt a child and LLL-reduce its saturated MW quotient."""
    frame = matrix(ZZ, frame)
    roots, root_basis, root_data = roots_and_data(frame)
    root_rank = root_data[0]
    if root_rank == 0:
        lll_columns = frame.LLL_gram()
        basis = lll_columns.transpose()
        adapted = basis * frame * basis.transpose()
        if abs(basis.det()) != 1:
            raise AssertionError("rootless LLL change is not unimodular")
        return {
            "basis": basis,
            "frame": adapted,
            "root_data": root_data,
            "simple_roots": matrix(ZZ, 0, frame.ncols()),
            "mw_height": adapted.change_ring(QQ),
            "root_smith_invariants": (),
            "root_lattice_primitive": True,
            "minimization_status": "ROOTLESS_LLL_MINIMIZED",
        }

    smith, smith_left, smith_right = root_basis.smith_form()
    if smith != smith_left * root_basis * smith_right:
        raise AssertionError("unexpected Smith normal form convention")
    smith_invariants = tuple(
        abs(ZZ(smith[index, index])) for index in range(root_rank)
    )
    if smith_invariants != (1,) * root_rank:
        # The child itself and its root data are still exact.  What is not
        # available without additional torsion/glue data is an integral
        # root-plus-MW coordinate system, so return a structured partial
        # result rather than discarding the child or pretending saturation.
        simple, _, _ = deterministic_simple_roots(frame)
        return {
            "basis": identity_matrix(ZZ, frame.nrows()),
            "frame": frame,
            "root_data": root_data,
            "simple_roots": simple,
            "mw_height": None,
            "root_smith_invariants": smith_invariants,
            "root_lattice_primitive": False,
            "minimization_status": "PARTIAL_NONPRIMITIVE_ROOT_LATTICE",
        }
    completion = smith_right.inverse()
    simple, _, _ = deterministic_simple_roots(frame)
    basis = simple.stack(completion[root_rank:])
    if abs(basis.det()) != 1:
        raise AssertionError("root adaptation did not give a unimodular basis")
    adapted = basis * frame * basis.transpose()
    cartan = adapted[:root_rank, :root_rank]
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling

    if height.nrows():
        scale = lcm(value.denominator() for value in height.list())
        lll_columns = (scale * height).change_ring(ZZ).LLL_gram()
        quotient_change = lll_columns.transpose()
        if abs(quotient_change.det()) != 1:
            raise AssertionError("MW LLL change is not unimodular")
        full_change = block_diagonal_matrix(
            identity_matrix(ZZ, root_rank), quotient_change
        )
        basis = full_change * basis
        adapted = basis * frame * basis.transpose()
        coupling = adapted[:root_rank, root_rank:]
        tail = adapted[root_rank:, root_rank:]
        height = tail - coupling.transpose() * cartan.inverse() * coupling
    return {
        "basis": basis,
        "frame": adapted,
        "root_data": root_data,
        "simple_roots": simple,
        "mw_height": height,
        "root_smith_invariants": smith_invariants,
        "root_lattice_primitive": True,
        "minimization_status": "ROOT_MW_MINIMIZED",
    }


def degree_q_neighbor(
    ns, divisor, old_fiber, curves=(), expected_old_fiber_degree=None,
    max_reduction_steps=10000,
):
    """Reduce, construct, and minimize a positive-degree child fibration.

    The result contains the reduced divisor, fixed-component sequence, exact
    child transport, raw child frame, and root/MW-minimized child frame.

    The lattice U-split is independent of ``q=D.F``.  Equation-level callers
    will choose a degree-``q`` generic-fibre Riemann--Roch ambient separately;
    retaining the actual degree here prevents the old degree-two neighbour
    label from being mistaken for a restriction of the Neron--Severi engine.
    Supplying ``expected_old_fiber_degree`` makes a pinned chain step fail
    closed if its declared degree changes under fixed-component reduction.
    """
    ns = matrix(ZZ, ns)
    old_fiber = vector(ZZ, old_fiber)
    divisor = vector(ZZ, divisor)
    if intersection(old_fiber, old_fiber, ns) != 0:
        raise ValueError("old_fiber must be isotropic")
    if intersection(divisor, divisor, ns) != 0:
        raise ValueError("divisor must be isotropic")

    reduced, reduction = reduce_fixed_components(
        divisor, curves, ns, max_steps=max_reduction_steps
    )
    old_fiber_degree = intersection(reduced, old_fiber, ns)
    if old_fiber_degree <= 0:
        raise ValueError("reduced divisor does not have positive degree over old_fiber")
    if (expected_old_fiber_degree is not None
            and old_fiber_degree != ZZ(expected_old_fiber_degree)):
        raise ValueError(
            "reduced divisor has old-fiber degree {}, expected {}".format(
                old_fiber_degree, expected_old_fiber_degree
            )
        )
    if any(intersection(reduced, curve, ns) < 0 for _, curve in curves):
        raise AssertionError("reduction left a negative supplied curve")
    split = primitive_hyperbolic_split(ns, reduced)
    minimized = minimize_child_frame(split["child_frame"])
    return {
        "raw_divisor": divisor,
        "reduced_divisor": reduced,
        "fixed_component_sequence": reduction,
        "old_fiber_degree": old_fiber_degree,
        "nonnegative_on_supplied_curves": True,
        "supplied_curve_count": len(curves),
        **split,
        "minimized_child_frame": minimized["frame"],
        "minimized_child_basis": minimized["basis"],
        "child_root_data": minimized["root_data"],
        "child_mw_height": minimized["mw_height"],
        "child_simple_roots": minimized["simple_roots"],
        "child_root_smith_invariants": minimized["root_smith_invariants"],
        "child_root_lattice_primitive": minimized["root_lattice_primitive"],
        "minimization_status": minimized["minimization_status"],
    }


def degree_two_neighbor(ns, divisor, old_fiber, curves=(), max_reduction_steps=10000):
    """Compatibility wrapper for pinned degree-two neighbour certificates."""
    return degree_q_neighbor(
        ns, divisor, old_fiber, curves,
        expected_old_fiber_degree=2, max_reduction_steps=max_reduction_steps,
    )


def neighbor_certificate(ns, old_fiber, curves, result, proof_metadata=None):
    """Build a versioned, JSON-ready exact neighbor certificate.

    ``proof_metadata`` is caller-supplied contextual evidence (for example a
    link to a global chamber proof).  It is retained separately from the
    engine's finite supplied-wall check and is never interpreted by this code.
    """
    ns = matrix(ZZ, ns)
    old_fiber = vector(ZZ, old_fiber)
    curves = tuple((str(name), vector(ZZ, curve)) for name, curve in curves)
    transport = matrix(ZZ, result["transport"])
    child = matrix(ZZ, result["child_frame"])
    minimized = matrix(ZZ, result["minimized_child_frame"])
    wall_records = [
        {
            "name": name,
            "class": vector_entries(curve),
            "self_intersection": int(intersection(curve, curve, ns)),
            "reduced_pairing": int(
                intersection(result["reduced_divisor"], curve, ns)
            ),
        }
        for name, curve in curves
    ]
    certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "status": "PASS_EXACT_NEIGHBOR",
        "proof_boundary": (
            "The engine proves only nonnegativity against the supplied "
            "curves.  Global nefness requires the caller's separate "
            "section/multisection certificate."
        ),
        "software": {"sage_version": str(SAGE_VERSION)},
        "source": {
            "ns_gram": matrix_rows(ns),
            "ns_gram_sha256": matrix_digest(ns),
            "old_fiber": vector_entries(old_fiber),
        },
        "divisor": {
            "raw": vector_entries(result["raw_divisor"]),
            "reduced": vector_entries(result["reduced_divisor"]),
            "old_fiber_degree": int(
                intersection(result["reduced_divisor"], old_fiber, ns)
            ),
            "nonnegative_on_supplied_curves": result[
                "nonnegative_on_supplied_curves"
            ],
            "supplied_curve_count": result["supplied_curve_count"],
            "fixed_component_sequence": [
                {"name": name, "pairing": int(pairing)}
                for name, pairing in result["fixed_component_sequence"]
            ],
        },
        "supplied_curves": wall_records,
        "child": {
            "frame": matrix_rows(child),
            "frame_sha256": matrix_digest(child),
            "transport": matrix_rows(transport),
            "transport_sha256": matrix_digest(transport),
            "minimized_frame": matrix_rows(minimized),
            "minimized_frame_sha256": matrix_digest(minimized),
            "minimized_basis": matrix_rows(result["minimized_child_basis"]),
            "root_data": [int(value) for value in result["child_root_data"]],
            "root_smith_invariants": [
                int(value) for value in result["child_root_smith_invariants"]
            ],
            "root_lattice_primitive": result["child_root_lattice_primitive"],
            "minimization_status": result["minimization_status"],
            "mw_height": (
                None if result["child_mw_height"] is None
                else rational_rows(result["child_mw_height"])
            ),
            "lifted_simple_components": [
                {"name": name, "class": vector_entries(divisor)}
                for name, divisor in lift_child_component_data(result)
            ],
        },
        "proof_metadata": proof_metadata or {},
    }
    certificate["certificate_sha256"] = json_digest(certificate)
    return certificate


def write_neighbor_certificate(path, certificate):
    """Write a canonical certificate and return its SHA-256 payload digest."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    return certificate["certificate_sha256"]


def neighbor_input(
    ns, old_fiber, divisor, curves, proof_metadata=None,
    expected_old_fiber_degree=None,
):
    """Build a versioned, JSON-ready input for the standalone engine runner."""
    ns = matrix(ZZ, ns)
    payload = {
        "schema": INPUT_SCHEMA,
        "ns_gram": matrix_rows(ns),
        "old_fiber": vector_entries(old_fiber),
        "divisor": vector_entries(divisor),
        "supplied_curves": [
            {"name": str(name), "class": vector_entries(curve)}
            for name, curve in curves
        ],
        "proof_metadata": proof_metadata or {},
    }
    if expected_old_fiber_degree is not None:
        expected_old_fiber_degree = ZZ(expected_old_fiber_degree)
        if expected_old_fiber_degree <= 0:
            raise ValueError("expected_old_fiber_degree must be positive")
        payload["expected_old_fiber_degree"] = int(expected_old_fiber_degree)
    return payload


def run_neighbor_input(data):
    """Validate and execute a versioned serialized exact-neighbor input."""
    if data.get("schema") != INPUT_SCHEMA:
        raise ValueError(f"expected schema {INPUT_SCHEMA}")
    required = {"ns_gram", "old_fiber", "divisor", "supplied_curves"}
    missing = required.difference(data)
    if missing:
        raise ValueError(f"neighbor input is missing {sorted(missing)}")
    ns = matrix(ZZ, data["ns_gram"])
    old_fiber = vector(ZZ, data["old_fiber"])
    divisor = vector(ZZ, data["divisor"])
    curves = tuple(
        (str(record["name"]), vector(ZZ, record["class"]))
        for record in data["supplied_curves"]
    )
    result = degree_q_neighbor(
        ns, divisor, old_fiber, curves,
        expected_old_fiber_degree=data.get("expected_old_fiber_degree"),
    )
    certificate = neighbor_certificate(
        ns, old_fiber, curves, result, data.get("proof_metadata", {})
    )
    return result, certificate
