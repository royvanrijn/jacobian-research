#!/usr/bin/env sage-python
"""Plan literal marked-U realizations in a fixed Neron--Severi lattice.

status: ACTIVE_COMPILER
claim: bounded literal-U enumeration with fail-closed physical gates
inputs: one explicit marked (NS,U,W) JSON configuration
outputs: a marked-U plan JSON record

This is deliberately separate from ``plan_inverse_ade_targets.sage``.  The
inverse-ADE planner mutates a positive-definite core.  This planner starts with
an explicit marked splitting

    NS = <F,F+O> + W(-1)

and enumerates target splittings by their physical cross intersections.  For
each tuple ``(F.F', F.O', O.F', O.O')`` it forms

    G_A = A^t J A - J,

enumerates ordered representations of ``G_A`` in ``W``, constructs the
literal primitive ``U'``, and computes its orthogonal frame.  Exact target
frame/root gates are intrinsic.  Nefness and effective-zero claims are
fail-closed and require hash-pinned external evidence.

The representation source is either a complete PARI shell enumeration or a
declared replay catalog.  Catalog mode is never reported as a complete search
of all representations unless its own scope says so explicitly.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
J = matrix(ZZ, [[0, 1], [1, 0]])


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def relative(path):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matrix_digest(value):
    encoded = json.dumps(rows(value), separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def load_matrix(specification, base=ROOT):
    """Load an integral matrix from rows or a whitespace-delimited file."""

    if isinstance(specification, list):
        return matrix(ZZ, specification)
    if not isinstance(specification, dict):
        raise TypeError("matrix specification must be rows or an object")
    if "rows" in specification:
        return matrix(ZZ, specification["rows"])
    if "path" not in specification:
        raise ValueError("matrix object needs rows or path")
    path = Path(specification["path"])
    path = path if path.is_absolute() else Path(base) / path
    if specification.get("sha256") and digest(path) != specification["sha256"]:
        raise ValueError(f"matrix input hash mismatch: {path}")
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def explicit_source(source, base=ROOT):
    """Validate and normalize an explicit ``(NS,U,W)`` source marking."""

    required = ("ns_gram", "u_basis_in_ns", "frame_basis_in_ns")
    missing = [key for key in required if key not in source]
    if missing:
        raise ValueError(f"explicit marked source is missing {missing}")
    ns = load_matrix(source["ns_gram"], base)
    source_u = load_matrix(source["u_basis_in_ns"], base)
    source_w = load_matrix(source["frame_basis_in_ns"], base)
    if not ns.is_symmetric() or ns.nrows() != ns.ncols():
        raise ValueError("NS Gram matrix must be square and symmetric")
    if source_u.nrows() != 2 or source_u.ncols() != ns.nrows():
        raise ValueError("source U must have two ambient NS rows")
    if source_w.nrows() != ns.nrows() - 2 or source_w.ncols() != ns.nrows():
        raise ValueError("source frame basis has the wrong shape")
    if source_u * ns * source_u.transpose() != J:
        raise ValueError("source U basis must be ordered as (F,F+O)")
    if source_u * ns * source_w.transpose():
        raise ValueError("source U and frame are not orthogonal")
    transport = source_u.stack(source_w)
    if abs(int(transport.det())) != 1:
        raise ValueError("the declared U+W splitting is not unimodular")
    frame = -(source_w * ns * source_w.transpose())
    if not frame.is_positive_definite():
        raise ValueError("the source frame is not positive definite")
    split = transport * ns * transport.transpose()
    if split != block_diagonal_matrix(J, -frame):
        raise ArithmeticError("marked-source splitting identity failed")
    source_zero = source_u.row(1) - source_u.row(0)
    if source_zero * ns * source_zero != -2:
        raise ArithmeticError("source zero does not have square -2")
    return {
        "label": source.get("label", "unnamed-marked-source"),
        "ns": ns,
        "source_u": source_u,
        "source_w": source_w,
        "split_transport": transport,
        "frame": frame,
        "source_zero": source_zero,
    }


def physical_cross_matrix(coordinates):
    d, s, t, z = map(ZZ, coordinates)
    return matrix(
        ZZ,
        [[d, d + s], [d + t, d + s + t + z]],
    )


def physical_coordinates(cross):
    d = ZZ(cross[0, 0])
    s = ZZ(cross[0, 1] - d)
    t = ZZ(cross[1, 0] - d)
    z = ZZ(cross[1, 1] - d - s - t)
    return tuple(map(int, (d, s, t, z)))


def declared_values(specification):
    if isinstance(specification, list):
        return sorted(set(map(int, specification)))
    if isinstance(specification, int):
        return [int(specification)]
    if not isinstance(specification, dict):
        raise TypeError("bound must be an integer, list, or {min,max}")
    minimum = int(specification["min"])
    maximum = int(specification["max"])
    step = int(specification.get("step", 1))
    if step <= 0 or maximum < minimum:
        raise ValueError("invalid finite bound")
    return list(range(minimum, maximum + 1, step))


def intersection_tuples(bounds):
    names = (
        "F_dot_F_prime",
        "F_dot_O_prime",
        "O_dot_F_prime",
        "O_dot_O_prime",
    )
    missing = [name for name in names if name not in bounds]
    if missing:
        raise ValueError(f"intersection box is missing {missing}")
    values = [declared_values(bounds[name]) for name in names]
    tuples = sorted(
        (d, s, t, z)
        for d in values[0]
        for s in values[1]
        for t in values[2]
        for z in values[3]
    )
    if bounds.get("physical_intersections_only", True):
        for d, s, t, z in tuples:
            if d <= 0 or s < 0 or t < 0 or (z < 0 and z != -2):
                raise ValueError(
                    "physical tuples require d>0, s,t>=0, and z=-2 or z>=0"
                )
    return tuples


def smith_invariants(gram):
    if gram.nrows() == 0:
        return []
    return [
        abs(int(entry))
        for entry in gram.elementary_divisors()
        if abs(int(entry)) > 1
    ]


def p_primary_invariants(invariants, prime):
    result = []
    for invariant in invariants:
        value = int(invariant)
        power = 1
        while value % prime == 0:
            power *= prime
            value //= prime
        if power > 1:
            result.append(power)
    return result


def prime_support(value):
    value = abs(ZZ(value))
    if value in (0, 1):
        return []
    return [int(prime) for prime, _ in value.factor()]


def prime_local_screen(cross, gram, constraints):
    """Apply representation-independent bridge constraints first."""

    rank = int(gram.rank())
    semidefinite = (
        gram[0, 0] >= 0 and gram[1, 1] >= 0 and gram.det() >= 0
    )
    checks = [{"name": "positive_semidefinite_G_A", "pass": bool(semidefinite)}]
    reasons = []
    if not semidefinite:
        reasons.append("G_A_NOT_POSITIVE_SEMIDEFINITE")
    expected_rank = constraints.get("relative_rank")
    if expected_rank is not None:
        passed = rank == int(expected_rank)
        checks.append({"name": "relative_rank", "pass": passed})
        if not passed:
            reasons.append("RELATIVE_RANK_MISMATCH")

    raw_determinant = abs(int(gram.det()))
    expected_raw = constraints.get("raw_bridge_determinant")
    if expected_raw is not None:
        passed = raw_determinant == int(expected_raw)
        checks.append({"name": "raw_bridge_determinant", "pass": passed})
        if not passed:
            reasons.append("RAW_BRIDGE_DETERMINANT_MISMATCH")

    saturated_determinant = constraints.get("saturated_bridge_determinant")
    saturation_index = constraints.get("saturation_index")
    if rank == 2 and saturated_determinant is not None and saturation_index is not None:
        expected = int(saturation_index) ** 2 * int(saturated_determinant)
        passed = raw_determinant == expected
        checks.append(
            {
                "name": "det_G_equals_index_squared_det_C",
                "pass": passed,
                "expected_raw_determinant": expected,
            }
        )
        if not passed:
            reasons.append("BRIDGE_SQUARE_INDEX_OBSTRUCTION")

    parity = constraints.get("det_A_parity")
    if parity is not None:
        observed = "even" if int(cross.det()) % 2 == 0 else "odd"
        passed = observed == parity
        checks.append({"name": "det_A_parity", "pass": passed, "observed": observed})
        if not passed:
            reasons.append("DET_A_PARITY_MISMATCH")

    declared_gram = constraints.get("saturated_bridge_gram")
    declared_invariants = None
    if declared_gram is not None:
        target_bridge = load_matrix(declared_gram)
        declared_invariants = smith_invariants(target_bridge)
        if rank == 2 and int(constraints.get("saturation_index", 1)) == 1:
            passed = pari(gram).qfisom(pari(target_bridge)) != 0
            checks.append({"name": "index_one_raw_bridge_isometry", "pass": passed})
            if not passed:
                reasons.append("INDEX_ONE_BRIDGE_LOCAL_CLASS_MISMATCH")

    required_primary = constraints.get("p_primary_invariants", {})
    if required_primary and declared_invariants is None:
        raise ValueError("p_primary_invariants requires saturated_bridge_gram")
    for prime_text, expected in sorted(required_primary.items(), key=lambda row: int(row[0])):
        prime = int(prime_text)
        observed = p_primary_invariants(declared_invariants, prime)
        passed = observed == list(map(int, expected))
        checks.append(
            {
                "name": f"p_primary_bridge_invariants_at_{prime}",
                "pass": passed,
                "observed": observed,
            }
        )
        if not passed:
            reasons.append(f"P_PRIMARY_BRIDGE_MISMATCH_AT_{prime}")

    return {
        "pass": not reasons,
        "rejection_reasons": reasons,
        "checks": checks,
        "relative_rank": rank,
        "raw_bridge_determinant": raw_determinant,
        "det_A": int(cross.det()),
        "prime_support": sorted(
            set(prime_support(raw_determinant))
            | set(prime_support(constraints.get("saturated_bridge_determinant", 1)))
        ),
        "ran_before_representation_enumeration": True,
    }


def signed_vectors_of_norm(gram, norm, cache):
    norm = int(norm)
    if norm < 0:
        return []
    if norm == 0:
        return [vector(ZZ, [0] * gram.nrows())]
    key = (matrix_digest(gram), norm)
    if key in cache:
        return cache[key]
    result = pari(gram).qfminim(norm)
    representatives = matrix(ZZ, result[2].sage()).columns()
    answer = set()
    for representative in representatives:
        value = vector(ZZ, representative)
        if value * gram * value != norm:
            continue
        positive = tuple(map(int, value))
        answer.add(positive)
        answer.add(tuple(-entry for entry in positive))
    cache[key] = [vector(ZZ, value) for value in sorted(answer)]
    return cache[key]


def binary_reduction(gram):
    rank = int(gram.rank())
    if rank == 2:
        change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    elif rank == 1:
        kernel = gram.right_kernel_matrix().row(0)
        k0, k1 = map(ZZ, kernel)
        gcd_value, bezout_0, bezout_1 = k0.xgcd(k1)
        if gcd_value != 1:
            raise ArithmeticError("rank-one binary kernel is not primitive")
        change = matrix(ZZ, [[bezout_1, -bezout_0], [k0, k1]])
    else:
        change = matrix.identity(ZZ, 2)
    if abs(int(change.det())) != 1:
        raise ArithmeticError("binary reduction is not unimodular")
    reduced = change * gram * change.transpose()
    return change, reduced


def exact_representations(frame, gram, cache, maximum=None):
    """Enumerate all ordered representations, unless ``maximum`` truncates."""

    change, reduced = binary_reduction(gram)
    first_shell = signed_vectors_of_norm(frame, reduced[0, 0], cache)
    second_shell = signed_vectors_of_norm(frame, reduced[1, 1], cache)
    inverse = change.inverse().change_ring(ZZ)
    answer = []
    truncated = False
    for first in first_shell:
        for second in second_shell:
            if first * frame * second != reduced[0, 1]:
                continue
            reduced_rows = matrix(ZZ, [list(first), list(second)])
            original = inverse * reduced_rows
            if original * frame * original.transpose() != gram:
                raise ArithmeticError("representation back-transform failed")
            answer.append(
                {
                    "id": "exact-" + matrix_digest(original)[:16],
                    "projected_vectors": rows(original),
                    "metadata": {},
                    "gate_evidence": {},
                }
            )
            if maximum is not None and len(answer) >= int(maximum):
                truncated = True
                return answer, {
                    "complete": False,
                    "truncated": truncated,
                    "first_shell_vectors": len(first_shell),
                    "second_shell_vectors": len(second_shell),
                    "binary_reduction": {
                        "basis_change": rows(change),
                        "reduced_gram": rows(reduced),
                    },
                }
    return answer, {
        "complete": True,
        "truncated": truncated,
        "first_shell_vectors": len(first_shell),
        "second_shell_vectors": len(second_shell),
        "binary_reduction": {
            "basis_change": rows(change),
            "reduced_gram": rows(reduced),
        },
    }


def catalog_representations(catalog, coordinates, gram, frame):
    answer = []
    for record in catalog.get("candidates", []):
        if tuple(map(int, record["intersection_coordinates"])) != tuple(coordinates):
            continue
        projected = matrix(ZZ, record["projected_vectors"])
        if projected.nrows() != 2 or projected.ncols() != frame.nrows():
            raise ValueError(f"catalog candidate {record['id']} has the wrong shape")
        if projected * frame * projected.transpose() != gram:
            raise ValueError(f"catalog candidate {record['id']} does not represent G_A")
        answer.append(record)
    answer.sort(key=lambda row: (row.get("id", ""), row["projected_vectors"]))
    return answer, {
        "complete": bool(catalog.get("complete_for_declared_box", False)),
        "truncated": False,
        "catalog_scope": catalog.get(
            "scope", "unspecified replay catalog; no global completeness claim"
        ),
        "catalog_candidates_for_cross_matrix": len(answer),
    }


def construct_literal_u(source, cross, projected):
    split_ns = block_diagonal_matrix(J, -source["frame"])
    u_part = cross.transpose() * J
    target_split = u_part.augment(projected)
    if target_split * split_ns * target_split.transpose() != J:
        raise ArithmeticError("constructed target basis is not U")
    complement_split = (target_split * split_ns).right_kernel_matrix()
    full_split = target_split.stack(complement_split)
    if abs(int(full_split.det())) != 1:
        raise ArithmeticError("constructed target U is not a primitive direct summand")
    child = -(complement_split * split_ns * complement_split.transpose())
    target_ambient = target_split * source["split_transport"]
    complement_ambient = complement_split * source["split_transport"]
    if target_ambient * source["ns"] * target_ambient.transpose() != J:
        raise ArithmeticError("ambient target U identity failed")
    if complement_ambient * source["ns"] * target_ambient.transpose():
        raise ArithmeticError("ambient target frame is not orthogonal to U'")
    return {
        "target_split": target_split,
        "target_ambient": target_ambient,
        "target_fibre": target_ambient.row(0),
        "target_zero": target_ambient.row(1) - target_ambient.row(0),
        "complement_split": complement_split,
        "complement_ambient": complement_ambient,
        "child": child,
        "transport_determinant": int(full_split.det()) * int(source["split_transport"].det()),
    }


def irreducible_ade_type(rank, signed_count):
    if signed_count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and signed_count == 2 * rank * (rank - 1):
        return f"D{rank}"
    exceptional = {(6, 72): "E6", (7, 126): "E7", (8, 240): "E8"}
    if (rank, signed_count) in exceptional:
        return exceptional[(rank, signed_count)]
    raise ArithmeticError(
        f"unclassified simply-laced component rank={rank}, roots={signed_count}"
    )


def root_signature(frame):
    enumeration = pari(frame).qfminim(2)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    lines = []
    seen = set()
    for representative in representatives:
        root = vector(ZZ, representative)
        if root * frame * root != 2:
            continue
        positive = tuple(map(int, root))
        negative = tuple(-entry for entry in positive)
        key = min(positive, negative)
        if key not in seen:
            seen.add(key)
            lines.append(vector(ZZ, key))
    lines.sort(key=lambda item: tuple(map(int, item)))
    if not lines:
        return {
            "ade_type": "rootless",
            "root_rank": 0,
            "signed_root_count": 0,
            "mordell_weil_rank_from_shioda_tate": frame.nrows(),
            "components": [],
        }
    root_matrix = matrix(ZZ, [list(item) for item in lines])
    pairings = root_matrix * frame * root_matrix.transpose()
    adjacency = [set() for _ in lines]
    for left in range(len(lines)):
        for right in range(left + 1, len(lines)):
            if pairings[left, right]:
                adjacency[left].add(right)
                adjacency[right].add(left)
    unseen = set(range(len(lines)))
    components = []
    labels = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        indices = []
        while stack:
            index = stack.pop()
            indices.append(index)
            for neighbor in sorted(adjacency[index]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        vectors = matrix(ZZ, [list(lines[index]) for index in indices])
        rank = int(vectors.rank())
        signed_count = 2 * len(indices)
        label = irreducible_ade_type(rank, signed_count)
        labels.append(label)
        components.append(
            {
                "type": label,
                "rank": rank,
                "signed_root_count": signed_count,
            }
        )
    counts = Counter(labels)
    family_order = {"A": 0, "D": 1, "E": 2}
    ordered = sorted(counts, key=lambda label: (family_order[label[0]], int(label[1:])))
    ade = "+".join(
        (str(counts[label]) if counts[label] > 1 else "") + label
        for label in ordered
    )
    return {
        "ade_type": ade,
        "root_rank": int(root_matrix.rank()),
        "signed_root_count": 2 * len(lines),
        "mordell_weil_rank_from_shioda_tate": (
            frame.nrows() - int(root_matrix.rank())
        ),
        "components": sorted(
            components, key=lambda row: (family_order[row["type"][0]], int(row["type"][1:]))
        ),
    }


def json_pointer(document, pointer):
    if pointer in ("", "/"):
        return document
    current = document
    for raw in pointer.lstrip("/").split("/"):
        token = raw.replace("~1", "/").replace("~0", "~")
        current = current[int(token)] if isinstance(current, list) else current[token]
    return current


def candidate_field(candidate, name):
    fields = {
        "candidate_id": candidate["candidate_id"],
        "primitive_u_basis_in_ambient_ns": candidate[
            "primitive_u_basis_in_ambient_ns"
        ],
        "target_fibre_in_ambient_ns": candidate["target_fibre_in_ambient_ns"],
        "target_zero_in_ambient_ns": candidate["target_zero_in_ambient_ns"],
        "child_frame_gram": candidate["child_frame_gram"],
        "child_frame_rank": len(candidate["child_frame_gram"]),
        "child_frame_determinant": candidate["child_frame_determinant"],
    }
    if name not in fields:
        raise ValueError(f"unsupported evidence candidate field {name}")
    return fields[name]


def verify_gate_evidence(name, evidence, candidate, cache):
    if not evidence:
        return {
            "gate": name,
            "pass": False,
            "reason": "MISSING_EXACT_EVIDENCE",
        }
    if evidence.get("verdict") != "PASS_EXACT":
        return {
            "gate": name,
            "pass": False,
            "reason": "EVIDENCE_VERDICT_IS_NOT_PASS_EXACT",
        }
    artifact = evidence.get("artifact")
    if not artifact:
        return {"gate": name, "pass": False, "reason": "MISSING_ARTIFACT"}
    path = Path(artifact["path"])
    path = path if path.is_absolute() else ROOT / path
    observed_hash = digest(path)
    if observed_hash != artifact.get("sha256"):
        return {"gate": name, "pass": False, "reason": "ARTIFACT_HASH_MISMATCH"}
    cache_key = str(path.resolve())
    if cache_key not in cache:
        cache[cache_key] = json.loads(path.read_text())
    payload = cache[cache_key]
    for field in ("schema", "status"):
        if artifact.get(field) is not None and payload.get(field) != artifact[field]:
            return {
                "gate": name,
                "pass": False,
                "reason": f"ARTIFACT_{field.upper()}_MISMATCH",
            }
    bindings = []
    for binding in evidence.get("bindings", []):
        observed = json_pointer(payload, binding["pointer"])
        if "candidate_field" in binding:
            expected = candidate_field(candidate, binding["candidate_field"])
        else:
            expected = binding["equals"]
        passed = observed == expected
        bindings.append({"pointer": binding["pointer"], "pass": passed})
        if not passed:
            return {
                "gate": name,
                "pass": False,
                "reason": "ARTIFACT_BINDING_MISMATCH",
                "bindings": bindings,
            }
    for binding in evidence.get("literal_bindings", []):
        observed = candidate_field(candidate, binding["candidate_field"])
        passed = observed == binding["equals"]
        bindings.append(
            {
                "candidate_field": binding["candidate_field"],
                "literal": True,
                "pass": passed,
            }
        )
        if not passed:
            return {
                "gate": name,
                "pass": False,
                "reason": "LITERAL_CANDIDATE_BINDING_MISMATCH",
                "bindings": bindings,
            }
    if not bindings:
        return {
            "gate": name,
            "pass": False,
            "reason": "EVIDENCE_NOT_BOUND_TO_LITERAL_CANDIDATE",
        }
    return {
        "gate": name,
        "pass": True,
        "reason": "HASH_PINNED_EXACT_EVIDENCE",
        "artifact": relative(path),
        "bindings": bindings,
    }


def target_gate(frame, signature, source_frame, target):
    checks = []
    if target.get("root_rank") is not None:
        checks.append(
            {
                "name": "root_rank",
                "pass": signature["root_rank"] == int(target["root_rank"]),
            }
        )
    if target.get("ade_type") is not None:
        checks.append(
            {
                "name": "ade_type",
                "pass": signature["ade_type"] == target["ade_type"],
            }
        )
    if target.get("frame_gram") is not None:
        desired = load_matrix(target["frame_gram"])
        checks.append(
            {
                "name": "integral_target_frame_isometry",
                "pass": pari(frame).qfisom(pari(desired)) != 0,
            }
        )
    if target.get("exclude_source_frame_isometry", False):
        checks.append(
            {
                "name": "exclude_source_frame_isometry",
                "pass": pari(frame).qfisom(pari(source_frame)) == 0,
            }
        )
    return {"pass": all(row["pass"] for row in checks), "checks": checks}


def saturated_bridge(projected, frame):
    raw_module = projected.row_module(ZZ)
    saturation = raw_module.saturation()
    raw_basis = raw_module.basis_matrix()
    saturated_basis = saturation.basis_matrix()
    index = abs(int(raw_module.index_in(saturation)))
    saturated_gram = saturated_basis * frame * saturated_basis.transpose()
    return {
        "rank": int(raw_module.rank()),
        "saturation_index": index,
        "raw_basis": rows(raw_basis),
        "saturated_basis": rows(saturated_basis),
        "saturated_gram": rows(saturated_gram),
        "saturated_discriminant_invariants": smith_invariants(saturated_gram),
        "saturated_determinant": abs(int(saturated_gram.det())),
    }


def bridge_representation_gate(record, constraints):
    checks = []
    expected_index = constraints.get("saturation_index")
    if expected_index is not None:
        checks.append(
            {
                "name": "saturation_index",
                "pass": record["saturation_index"] == int(expected_index),
            }
        )
    expected_det = constraints.get("saturated_bridge_determinant")
    if expected_det is not None:
        checks.append(
            {
                "name": "saturated_bridge_determinant",
                "pass": record["saturated_determinant"] == int(expected_det),
            }
        )
    if constraints.get("saturated_bridge_gram") is not None:
        desired = load_matrix(constraints["saturated_bridge_gram"])
        observed = matrix(ZZ, record["saturated_gram"])
        checks.append(
            {
                "name": "saturated_bridge_integral_isometry",
                "pass": pari(observed).qfisom(pari(desired)) != 0,
            }
        )
    return {"pass": all(row["pass"] for row in checks), "checks": checks}


def equation_cost(coordinates, gram, projected, metadata):
    d, s, t, z = map(int, coordinates)
    declared = metadata.get("equation_complexity", {})
    rr_ambient = metadata.get("expected_riemann_roch_ambient")
    values = {
        "old_fibre_degree": d,
        "old_zero_degree": t,
        "new_zero_old_fibre_degree": s,
        "zero_zero_intersection": z,
        "maximum_projection_norm_half": max(map(int, gram.diagonal())) // 2,
        "projected_coordinate_l1": sum(abs(int(entry)) for entry in projected.list()),
        "expected_riemann_roch_ambient": rr_ambient,
        "coefficient_l1": declared.get("coefficient_l1"),
        "coordinate_input_bits": declared.get("coordinate_input_bits"),
    }
    fields = list(values)
    return {
        "fields": fields,
        "values": values,
        "vector": [values[field] for field in fields],
        "declared_metrics": declared,
    }


def within_cost_bound(cost, bound):
    values = cost["values"]
    reasons = []
    aliases = {
        "maximum_old_fibre_degree": "old_fibre_degree",
        "maximum_old_zero_degree": "old_zero_degree",
        "maximum_new_zero_old_fibre_degree": "new_zero_old_fibre_degree",
        "maximum_coefficient_l1": "coefficient_l1",
        "maximum_coordinate_input_bits": "coordinate_input_bits",
        "maximum_projected_coordinate_l1": "projected_coordinate_l1",
    }
    for key, field in aliases.items():
        if key not in bound:
            continue
        observed = values[field]
        if observed is None or int(observed) > int(bound[key]):
            reasons.append({"bound": key, "observed": observed, "maximum": bound[key]})
    return not reasons, reasons


def cost_sort_key(candidate):
    values = candidate["equation_facing_cost"]["values"]
    large = 10**100
    return (
        values["old_fibre_degree"],
        values["old_zero_degree"],
        values["new_zero_old_fibre_degree"],
        values["zero_zero_intersection"],
        values["expected_riemann_roch_ambient"] if values["expected_riemann_roch_ambient"] is not None else large,
        values["coefficient_l1"] if values["coefficient_l1"] is not None else large,
        values["coordinate_input_bits"] if values["coordinate_input_bits"] is not None else large,
        values["projected_coordinate_l1"],
        candidate["candidate_id"],
    )


def plan(configuration, base=ROOT):
    """Run the bounded marked-U planner and return a JSON-safe certificate."""

    source = explicit_source(configuration["source"], base)
    target = dict(configuration.get("target", {}))
    if target.get("frame_gram") is not None:
        target_frame = load_matrix(target["frame_gram"], base)
        if target_frame.nrows() != source["frame"].nrows():
            raise ValueError("target frame rank differs from source frame rank")
        if abs(int(target_frame.det())) != abs(int(source["frame"].det())):
            raise ValueError("target frame determinant differs from the fixed NS determinant")
        # All later gates operate on the normalized in-memory matrix.  This
        # also makes paths relative to the input manifest behave consistently.
        target["frame_gram"] = rows(target_frame)
    tuples = intersection_tuples(configuration["intersection_box"])
    representation_source = configuration.get("representations", {"mode": "exact"})
    mode = representation_source.get("mode", "exact")
    if mode not in ("exact", "catalog"):
        raise ValueError("representations.mode must be exact or catalog")
    constraints = dict(configuration.get("prime_local_bridge_constraints", {}))
    if constraints.get("saturated_bridge_gram") is not None:
        constraints["saturated_bridge_gram"] = rows(
            load_matrix(constraints["saturated_bridge_gram"], base)
        )
    physical_required = configuration.get("physical_gates", {}).get("required", True)
    cost_bound = configuration.get("equation_cost_bound", {})
    shell_cache = {}
    artifact_cache = {}
    counters = Counter()
    matrix_summaries = []
    retained = []
    pending = []
    all_representation_searches_complete = True

    for coordinates in tuples:
        counters["cross_matrices_enumerated"] += 1
        cross = physical_cross_matrix(coordinates)
        gram = cross.transpose() * J * cross - J
        local = prime_local_screen(cross, gram, constraints)
        summary = {
            "intersection_coordinates": list(coordinates),
            "cross_pairing_A": rows(cross),
            "positive_projection_gram_G_A": rows(gram),
            "prime_local_screen": local,
        }
        if not local["pass"]:
            counters["cross_matrices_rejected_before_representations"] += 1
            summary["status"] = "REJECTED_PRIME_LOCAL_BEFORE_REPRESENTATIONS"
            summary["representations_enumerated"] = 0
            matrix_summaries.append(summary)
            continue
        counters["cross_matrices_reaching_representations"] += 1
        if mode == "exact":
            records, enumeration = exact_representations(
                source["frame"],
                gram,
                shell_cache,
                representation_source.get("maximum_representations_per_matrix"),
            )
        else:
            records, enumeration = catalog_representations(
                representation_source, coordinates, gram, source["frame"]
            )
        all_representation_searches_complete &= enumeration["complete"]
        summary["representation_enumeration"] = enumeration
        summary["representations_enumerated"] = len(records)
        summary["status"] = "REPRESENTATIONS_ENUMERATED"
        counters["representations_enumerated"] += len(records)

        for record in records:
            projected = matrix(ZZ, record["projected_vectors"])
            bridge = saturated_bridge(projected, source["frame"])
            bridge_gate = bridge_representation_gate(bridge, constraints)
            if not bridge_gate["pass"]:
                counters["representations_rejected_by_exact_bridge_gate"] += 1
                continue
            literal = construct_literal_u(source, cross, projected)
            counters["literal_primitive_u_constructed"] += 1
            observed_coordinates = physical_coordinates(
                source["source_u"]
                * source["ns"]
                * literal["target_ambient"].transpose()
            )
            if observed_coordinates != tuple(coordinates):
                raise ArithmeticError("literal U' does not recover its physical cross tuple")
            if literal["target_zero"] * source["ns"] * literal["target_zero"] != -2:
                raise ArithmeticError("target pseudo-zero does not have square -2")
            if literal["target_fibre"] * source["ns"] * literal["target_zero"] != 1:
                raise ArithmeticError("target pseudo-zero is not a section class")
            signature = root_signature(literal["child"])
            target_result = target_gate(
                literal["child"], signature, source["frame"], target
            )
            if not target_result["pass"]:
                counters["literal_u_rejected_by_target_gate"] += 1
                continue

            candidate = {
                "candidate_id": record.get(
                    "id", "candidate-" + matrix_digest(projected)[:16]
                ),
                "intersection_coordinates": {
                    "F_dot_F_prime": coordinates[0],
                    "F_dot_O_prime": coordinates[1],
                    "O_dot_F_prime": coordinates[2],
                    "O_dot_O_prime": coordinates[3],
                },
                "cross_pairing_A": rows(cross),
                "positive_projection_gram_G_A": rows(gram),
                "projected_vectors_in_source_frame": rows(projected),
                "primitive_u_basis_in_source_split_coordinates": rows(
                    literal["target_split"]
                ),
                "primitive_u_basis_in_ambient_ns": rows(literal["target_ambient"]),
                "target_fibre_in_ambient_ns": list(map(int, literal["target_fibre"])),
                "target_zero_in_ambient_ns": list(map(int, literal["target_zero"])),
                "child_frame_basis_in_ambient_ns": rows(
                    literal["complement_ambient"]
                ),
                "child_frame_gram": rows(literal["child"]),
                "child_frame_determinant": abs(int(literal["child"].det())),
                "child_frame_gram_sha256": matrix_digest(literal["child"]),
                "unimodular_transport_determinant": literal[
                    "transport_determinant"
                ],
                "root_ADE_gate": signature,
                "target_gate": target_result,
                "bridge": bridge,
                "bridge_gate": bridge_gate,
                "metadata": record.get("metadata", {}),
            }
            cost = equation_cost(
                coordinates, gram, projected, candidate["metadata"]
            )
            candidate["equation_facing_cost"] = cost
            inside, cost_rejections = within_cost_bound(cost, cost_bound)
            candidate["equation_cost_bound_gate"] = {
                "pass": inside,
                "rejections": cost_rejections,
            }
            if not inside:
                counters["literal_u_rejected_by_equation_cost_bound"] += 1
                continue

            evidence = record.get("gate_evidence", {})
            nefness = verify_gate_evidence(
                "nefness", evidence.get("nefness"), candidate, artifact_cache
            )
            effective_zero = verify_gate_evidence(
                "effective_zero",
                evidence.get("effective_zero"),
                candidate,
                artifact_cache,
            )
            shared_zero_consistency = {
                "required": coordinates[3] == -2,
                "pass": (
                    coordinates[3] != -2
                    or literal["target_zero"] == source["source_zero"]
                ),
            }
            candidate["physical_gates"] = {
                "nefness": nefness,
                "effective_zero": effective_zero,
                "shared_zero_consistency": shared_zero_consistency,
                "pass": (
                    nefness["pass"]
                    and effective_zero["pass"]
                    and shared_zero_consistency["pass"]
                ),
                "required": physical_required,
            }
            if physical_required and not candidate["physical_gates"]["pass"]:
                counters["lattice_hits_pending_physical_gates"] += 1
                pending.append(candidate)
                continue
            counters["retained_realizations"] += 1
            retained.append(candidate)
        matrix_summaries.append(summary)

    retained.sort(key=cost_sort_key)
    pending.sort(key=cost_sort_key)
    complete = mode == "exact" and all_representation_searches_complete
    if retained:
        status = "PASS_MARKED_U_REALIZATIONS_FOUND"
    elif pending:
        status = "LATTICE_HITS_PENDING_PHYSICAL_GATES"
    else:
        status = "NO_REALIZATION_IN_DECLARED_SEARCH"
    return {
        "schema": "elkies-k3.marked-u-realization-plan.v1",
        "status": status,
        "planner_kind": "marked-U elliptic-fibration realization planner",
        "separate_from_core_planner": True,
        "source": {
            "label": source["label"],
            "ns_rank": source["ns"].nrows(),
            "frame_rank": source["frame"].nrows(),
            "frame_determinant": abs(int(source["frame"].det())),
            "explicit_NS_U_W_validated": True,
        },
        "target_request": {
            "root_rank": target.get("root_rank"),
            "ade_type": target.get("ade_type"),
            "target_frame_supplied": target.get("frame_gram") is not None,
            "exclude_source_frame_isometry": target.get(
                "exclude_source_frame_isometry", False
            ),
        },
        "search_order": (
            "lexicographic (F.F',F.O',O.F',O.O'); representations deterministic "
            "inside each cross matrix; retained candidates ranked by equation cost"
        ),
        "intersection_tuples": [list(value) for value in tuples],
        "prime_local_bridge_constraints": constraints,
        "representation_mode": mode,
        "search_completeness": {
            "complete_over_all_integral_representations_in_declared_box": complete,
            "complete_within_declared_catalog": mode == "catalog",
            "catalog_scope": representation_source.get("scope") if mode == "catalog" else None,
            "no_global_claim_from_catalog_mode": mode == "catalog",
        },
        "equation_cost_bound": cost_bound,
        "counters": dict(sorted(counters.items())),
        "cross_matrix_summaries": matrix_summaries,
        "selected_realization": retained[0] if retained else None,
        "retained_realizations": retained,
        "pending_physical_realizations": pending,
        "proof_boundary": {
            "intrinsic": (
                "Exact integral construction of every enumerated U', automatic primitivity, "
                "orthogonal frame, root/ADE classification, target isometry, saturation, "
                "and declared prime-local bridge checks."
            ),
            "external": (
                "Nefness and effective-zero acceptance comes only from hash-pinned exact "
                "evidence bound to the literal candidate. Equation compilation and arithmetic "
                "rank are downstream gates, not consequences of a lattice hit."
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    input_path = arguments.input.resolve()
    configuration = json.loads(input_path.read_text())
    result = plan(configuration, input_path.parent)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(encoded, end="")
        return
    output = arguments.output
    output = output if output.is_absolute() else ROOT / output
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale planner output: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(f"{result['status']}|output={relative(output)}")


if __name__ == "__main__":
    main()
