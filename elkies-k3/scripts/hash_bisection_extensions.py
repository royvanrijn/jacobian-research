#!/usr/bin/env python3
"""Canonicalize quadratic bisection covers and inspect exact collisions.

The equation-level stage supplies records for covers ``z^2=f(t)`` with
``f in QQ(t)^*``.  This utility keeps the rational constant squareclass as
well as odd irreducible factors, so it hashes extensions of ``QQ(t)`` rather
than only their geometric branch divisors.  Records with equal keys are a
collision bucket.  If exact anti-invariant-section coordinates in a declared
twist height lattice are present, the script returns their height matrix and
its rational rank; it never invents those coordinates from a branch collision.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import sympy as sp


SCHEMA = "elkies-k3.bisection-extension-input.v1"
OUTPUT_SCHEMA = "elkies-k3.bisection-extension-collisions.v1"
ROOT = Path(__file__).resolve().parents[2]
RANK17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"


def rational_text(value: sp.Rational | Fraction | int) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def parse_rational(value: Any) -> sp.Rational:
    return sp.Rational(str(value))


def parse_nonnegative_integer(value: Any, *, name: str) -> int:
    text = str(value)
    if not text.isdigit():
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(text)


def parse_integer_vector(value: Any, *, name: str, length: int) -> tuple[int, ...]:
    """Read a pinned lattice vector without accepting rational coordinates."""

    entries = value if isinstance(value, (list, tuple)) else str(value).split()
    if len(entries) != length:
        raise ValueError(f"{name}: expected {length} integral coordinates")
    parsed: list[int] = []
    for entry in entries:
        text = str(entry)
        if text.startswith("+"):
            text = text[1:]
        if not text or (text.startswith("-") and not text[1:].isdigit()) or (
            not text.startswith("-") and not text.isdigit()
        ):
            raise ValueError(f"{name}: coordinates must be integers")
        parsed.append(int(text))
    return tuple(parsed)


def load_integral_gram(path: Path, *, rank: int) -> tuple[tuple[int, ...], ...]:
    rows = tuple(
        tuple(int(entry) for entry in line.split())
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )
    if len(rows) != rank or any(len(row) != rank for row in rows):
        raise ValueError(f"{path}: expected a {rank} by {rank} integral Gram matrix")
    if any(rows[left][right] != rows[right][left] for left in range(rank) for right in range(rank)):
        raise ValueError(f"{path}: Gram matrix is not symmetric")
    return rows


def integral_norm(vector: tuple[int, ...], gram: tuple[tuple[int, ...], ...]) -> int:
    return sum(
        vector[left] * gram[left][right] * vector[right]
        for left in range(len(vector))
        for right in range(len(vector))
    )


def polynomial_from_ascending(
    coefficients: list[Any], variable: sp.Symbol, *, require_nonzero: bool = True
) -> sp.Poly:
    if not coefficients:
        raise ValueError("polynomial coefficients must not be empty")
    expression = sum(parse_rational(value) * variable**index for index, value in enumerate(coefficients))
    polynomial = sp.Poly(expression, variable, domain=sp.QQ)
    if require_nonzero and polynomial.is_zero:
        raise ValueError("branch numerator and denominator must be nonzero")
    return polynomial


def polynomial_key(polynomial: sp.Poly) -> tuple[str, ...]:
    monic = polynomial.monic()
    return tuple(rational_text(value) for value in reversed(monic.all_coeffs()))


def rational_squareclass(value: sp.Rational) -> tuple[int, tuple[int, ...]]:
    """Represent the class of a nonzero rational modulo rational squares."""

    value = sp.Rational(value)
    if not value:
        raise ValueError("the zero function has no quadratic squareclass")
    sign = -1 if value < 0 else 1
    parity_primes: set[int] = set()
    for integer in (abs(int(value.p)), int(value.q)):
        for prime, exponent in sp.factorint(integer).items():
            if exponent % 2:
                if prime in parity_primes:
                    parity_primes.remove(prime)
                else:
                    parity_primes.add(prime)
    return sign, tuple(sorted(parity_primes))


def reduced_branch_fraction(
    branch: dict[str, Any], variable: sp.Symbol
) -> tuple[sp.Poly, sp.Poly]:
    """Return the coprime numerator and denominator of a declared branch."""

    numerator = polynomial_from_ascending(branch["numerator_coefficients"], variable)
    denominator = polynomial_from_ascending(branch["denominator_coefficients"], variable)
    reduced_numerator, reduced_denominator = sp.fraction(
        sp.cancel(numerator.as_expr() / denominator.as_expr())
    )
    return (
        sp.Poly(reduced_numerator, variable, domain=sp.QQ),
        sp.Poly(reduced_denominator, variable, domain=sp.QQ),
    )


def geometric_branch_divisor(branch: dict[str, Any], variable: sp.Symbol) -> dict[str, Any]:
    """Compute the geometric branch divisor of ``z^2=f(t)`` exactly.

    A smooth rational bisection has a connected degree-two map to the base
    ``P1``.  Riemann--Hurwitz therefore forces precisely two geometric branch
    points.  Finite points are the irreducible factors occurring to odd order;
    infinity is included when the degree difference is odd.  We retain the
    factor labels only for auditability: the extension hash below is the
    stronger rational squareclass key.
    """

    numerator, denominator = reduced_branch_fraction(branch, variable)
    finite_factors: list[dict[str, Any]] = []
    finite_degree = 0
    for source, polynomial in (("numerator", numerator), ("denominator", denominator)):
        _, factors = sp.factor_list(polynomial)
        for factor, exponent in factors:
            if int(exponent) % 2:
                degree = int(factor.degree())
                finite_degree += degree
                finite_factors.append({
                    "source": source,
                    "monic_ascending_coefficients": list(polynomial_key(factor)),
                    "degree": degree,
                })
    infinity = (numerator.degree() - denominator.degree()) % 2 == 1
    return {
        "finite_odd_factors": sorted(
            finite_factors,
            key=lambda item: (item["source"], item["monic_ascending_coefficients"]),
        ),
        "infinity": infinity,
        "geometric_degree": finite_degree + int(infinity),
        "reduced_numerator_ascending_coefficients": polynomial_coefficients_ascending(numerator),
        "reduced_denominator_ascending_coefficients": polynomial_coefficients_ascending(denominator),
    }


def extension_key(branch: dict[str, Any], variable: sp.Symbol) -> dict[str, Any]:
    """Return a canonical key for ``QQ(t)^*/QQ(t)^{*2}``."""

    numerator, denominator = reduced_branch_fraction(branch, variable)
    numerator_constant, numerator_factors = sp.factor_list(numerator)
    denominator_constant, denominator_factors = sp.factor_list(denominator)
    constant_key = rational_squareclass(sp.Rational(numerator_constant) / sp.Rational(denominator_constant))
    parity: dict[tuple[str, ...], int] = {}
    for factor, exponent in (*numerator_factors, *denominator_factors):
        key = polynomial_key(factor)
        parity[key] = (parity.get(key, 0) + int(exponent)) % 2
    factors = [list(key) for key in sorted(key for key, exponent in parity.items() if exponent)]
    return {
        "constant_squareclass": {
            "sign": constant_key[0],
            "odd_primes": list(constant_key[1]),
        },
        "monic_odd_factors_ascending_coefficients": factors,
    }


def polynomial_coefficients_ascending(polynomial: sp.Poly) -> list[str]:
    """Serialize an exact polynomial in the input coefficient convention."""

    return [rational_text(polynomial.nth(index)) for index in range(polynomial.degree() + 1)]


def branch_from_quadratic_cover(
    cover: dict[str, Any], variable: sp.Symbol
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return the branch squareclass for ``a*z^2+b*z+c=0``.

    Completing the square gives ``(2*a*z+b)^2=b^2-4*a*c``.  The
    discriminant therefore describes exactly the same quadratic extension,
    even when the affine quadratic model is not monic.  A zero leading
    coefficient would not define a quadratic bisection and is rejected.
    """

    for name in ("leading_coefficients", "linear_coefficients", "constant_coefficients"):
        if name not in cover:
            raise ValueError(f"quadratic_cover: missing {name}")
    leading = polynomial_from_ascending(
        cover["leading_coefficients"], variable, require_nonzero=False
    )
    linear = polynomial_from_ascending(
        cover["linear_coefficients"], variable, require_nonzero=False
    )
    constant = polynomial_from_ascending(
        cover["constant_coefficients"], variable, require_nonzero=False
    )
    if leading.is_zero:
        raise ValueError("quadratic_cover: leading coefficient must be nonzero")
    discriminant = linear * linear - 4 * leading * constant
    if discriminant.is_zero:
        raise ValueError("quadratic_cover: discriminant must be nonzero")
    branch = {
        "numerator_coefficients": polynomial_coefficients_ascending(discriminant),
        "denominator_coefficients": ["1"],
    }
    provenance = {
        "kind": "quadratic_cover_discriminant",
        "relation": "leading*z^2 + linear*z + constant = 0",
        "discriminant_ascending_coefficients": branch["numerator_coefficients"],
    }
    return branch, provenance


def branch_from_record(record: dict[str, Any], variable: sp.Symbol) -> tuple[dict[str, Any], dict[str, Any]]:
    """Read one of the two exact branch-model forms accepted by the schema."""

    has_branch = "branch" in record
    has_cover = "quadratic_cover" in record
    if has_branch == has_cover:
        raise ValueError(
            f"{record['label']}: supply exactly one of branch or quadratic_cover"
        )
    if has_branch:
        branch = record["branch"]
        if not isinstance(branch, dict):
            raise ValueError(f"{record['label']}: branch must be an object")
        if set(branch) != {"numerator_coefficients", "denominator_coefficients"}:
            raise ValueError(
                f"{record['label']}: branch requires numerator_coefficients and denominator_coefficients"
            )
        return branch, {"kind": "declared_branch_rational_function"}
    if not isinstance(record["quadratic_cover"], dict):
        raise ValueError(f"{record['label']}: quadratic_cover must be an object")
    branch, provenance = branch_from_quadratic_cover(record["quadratic_cover"], variable)
    # The production equation compiler supplies the resolved chord denominator
    # h and the reduced quadratic q together with the displayed bisection.
    # Verify discriminant=h^2*q exactly, then factor only q.  Factoring the
    # unreduced degree-eight discriminant with enormous rational coefficients
    # is needlessly expensive on the complete 39,120-orbit batch.
    chord = record.get("residual_chord")
    trace = record.get("trace_section")
    if isinstance(chord, dict) and isinstance(trace, dict) and (
        "q_coefficients" in chord and "h_coefficients" in trace
    ):
        discriminant = polynomial_from_ascending(
            branch["numerator_coefficients"], variable
        )
        h = polynomial_from_ascending(trace["h_coefficients"], variable)
        q = polynomial_from_ascending(chord["q_coefficients"], variable)
        if discriminant != h * h * q:
            raise ValueError(
                f"{record['label']}: declared chord reduction does not satisfy "
                "quadratic discriminant=h^2*q"
            )
        provenance = {
            **provenance,
            "kind": "quadratic_cover_discriminant_with_exact_square_reduction",
            "square_multiplier_ascending_coefficients": polynomial_coefficients_ascending(h),
            "reduced_squareclass_ascending_coefficients": polynomial_coefficients_ascending(q),
            "square_reduction_identity": "discriminant=h^2*q",
        }
        return {
            "numerator_coefficients": polynomial_coefficients_ascending(q),
            "denominator_coefficients": ["1"],
        }, provenance
    return branch, provenance


def is_trivial_squareclass(key: dict[str, Any]) -> bool:
    return (
        key["constant_squareclass"] == {"sign": 1, "odd_primes": []}
        and not key["monic_odd_factors_ascending_coefficients"]
    )


def key_digest(key: dict[str, Any]) -> str:
    payload = json.dumps(key, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def parse_matrix(rows: list[list[Any]]) -> sp.Matrix:
    matrix = sp.Matrix([[parse_rational(value) for value in row] for row in rows])
    if matrix.rows != matrix.cols or matrix != matrix.transpose():
        raise ValueError("twist_height_gram must be square and symmetric")
    # The free Mordell--Weil height pairing is positive definite.  Merely
    # accepting a symmetric matrix would allow an indefinite formal matrix to
    # manufacture a full-rank collision contribution that cannot be a height
    # lattice.  Sylvester's criterion is exact over QQ and avoids a numerical
    # eigenvalue test.
    if matrix.rows == 0 or any(matrix[:size, :size].det() <= 0 for size in range(1, matrix.rows + 1)):
        raise ValueError("twist_height_gram must be positive definite")
    return matrix


def matrix_payload(labels: list[str], heights: sp.Matrix, *, source: str) -> dict[str, Any]:
    """Serialize one exact positive-definite anti-invariant height matrix."""

    if heights != heights.transpose() or heights.rows != len(labels):
        raise ValueError("anti-invariant height matrix has incompatible dimensions")
    if any(heights[:size, :size].det() <= 0 for size in range(1, heights.rows + 1)):
        raise ValueError("anti-invariant height matrix must be positive definite")
    return {
        "status": "COMPUTED_FROM_" + source,
        "labels": labels,
        "matrix": [[rational_text(value) for value in row] for row in heights.tolist()],
        "rank": int(heights.rank()),
        "determinant": rational_text(heights.det()),
    }


def declared_twist_height_data(
    records: list[dict[str, Any]], height_gram: sp.Matrix | None
) -> dict[str, Any] | None:
    if height_gram is None or any("anti_invariant_coordinates" not in record for record in records):
        return None
    coordinates = sp.Matrix([
        [parse_rational(value) for value in record["anti_invariant_coordinates"]]
        for record in records
    ])
    if coordinates.cols != height_gram.cols:
        raise ValueError("anti-invariant coordinate dimension does not match twist_height_gram")
    heights = coordinates * height_gram * coordinates.transpose()
    return matrix_payload(
        [str(record["label"]) for record in records], heights,
        source="DECLARED_TWIST_HEIGHT_LATTICE",
    )


def pair_key(left: str, right: str) -> str:
    """Canonical key for an unordered pair of distinct or equal labels."""

    return "|".join(sorted((left, right)))


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((left, right)))


def parse_declared_pairings(
    geometry: dict[str, Any], key: str, labels: list[str]
) -> dict[tuple[str, str], sp.Rational]:
    source = geometry.get(key)
    if not isinstance(source, dict):
        raise ValueError(f"rootless_smooth_base_change: missing {key}")
    expected = {pair_key(left, right) for offset, left in enumerate(labels) for right in labels[offset:]}
    if set(source) != expected:
        missing = len(expected-set(source))
        unexpected = len(set(source)-expected)
        raise ValueError(
            f"rootless_smooth_base_change.{key}: pair coverage mismatch "
            f"(missing={missing}, unexpected={unexpected})"
        )
    return {
        canonical_pair(left, right): parse_rational(source[pair_key(left, right)])
        for offset, left in enumerate(labels) for right in labels[offset:]
    }


def geometric_lift_height_data(
    records: list[dict[str, Any]], geometry: dict[str, Any] | None
) -> dict[str, Any] | None:
    """Compute the anti-invariant height form from lifted section intersections.

    This route applies only after an equation-level certificate has established
    that the common quadratic cover branches at smooth fibres and the pulled
    back fibration remains rootless.  If P_i is a lift of B_i and tau is the
    deck involution, then chi(Y)=4, P_i^2=-4, P_i.(tau P_i)=2, and

        <P_i-tau(P_i), P_j-tau(P_j)>
          = 2*(P_i.tau(P_j)-P_i.P_j).

    The diagonal values are checked rather than inferred from user text.
    """

    if geometry is None:
        return None
    if not isinstance(geometry, dict):
        raise ValueError("rootless_smooth_base_change must be an object")
    required = {
        "base_change_chi", "branch_fibres_smooth", "base_change_rootless",
        "same_sheet_intersections", "conjugate_sheet_intersections",
    }
    if set(geometry) != required:
        raise ValueError("rootless_smooth_base_change has an invalid field set")
    if parse_rational(geometry["base_change_chi"]) != 4:
        raise ValueError("rootless_smooth_base_change requires base_change_chi=4")
    if geometry["branch_fibres_smooth"] is not True or geometry["base_change_rootless"] is not True:
        raise ValueError(
            "rootless_smooth_base_change requires certified smooth branch fibres "
            "and a rootless pullback"
        )
    labels = [str(record["label"]) for record in records]
    same = parse_declared_pairings(geometry, "same_sheet_intersections", labels)
    conjugate = parse_declared_pairings(geometry, "conjugate_sheet_intersections", labels)
    for label in labels:
        if same[canonical_pair(label, label)] != -4 or conjugate[canonical_pair(label, label)] != 2:
            raise ValueError(
                "rootless_smooth_base_change diagonal intersections must be "
                "P^2=-4 and P.tau(P)=2"
            )
    heights = sp.Matrix([
        [2*(conjugate[canonical_pair(left, right)]-same[canonical_pair(left, right)]) for right in labels]
        for left in labels
    ])
    result = matrix_payload(labels, heights, source="GEOMETRIC_LIFT_INTERSECTIONS")
    result["formula"] = "2*(P_i.tau(P_j)-P_i.P_j)"
    result["base_change_chi"] = "4"
    return result


def collision_height_data(
    records: list[dict[str, Any]], height_gram: sp.Matrix | None,
    geometry: dict[str, Any] | None,
) -> dict[str, Any]:
    declared = declared_twist_height_data(records, height_gram)
    geometric = geometric_lift_height_data(records, geometry)
    if declared is not None and geometric is not None:
        if declared["matrix"] != geometric["matrix"]:
            raise ValueError("declared twist heights disagree with geometric lift intersections")
        return geometric
    if geometric is not None:
        return geometric
    if declared is not None:
        return declared
    if height_gram is None:
        return {"status": "MISSING_TWIST_HEIGHT_GRAM_OR_GEOMETRIC_LIFTS"}
    return {"status": "MISSING_ANTI_INVARIANT_COORDINATES"}


def verify_lattice_orbit_coverage(
    payload: dict[str, Any], records: list[dict[str, Any]]
) -> dict[str, Any]:
    """Require exactly one equation-level record per enumerated lattice orbit.

    This optional gate binds a future large collision run to the exact TSV
    emitted by ``enumerate_rootless_bisection_orbits.sage``.  It deliberately
    identifies an orbit by its class in M/2M, not by a chosen norm-ten vector:
    an equation may use any section-translate of that representative.
    """

    specification = payload.get("required_lattice_orbits")
    if specification is None:
        return {"status": "NOT_REQUESTED"}
    allowed = {
        "table", "sha256", "frame_artifact", "frame_sha256", "vector_key",
    }
    if not isinstance(specification, dict) or set(specification) - allowed:
        raise ValueError(
            "required_lattice_orbits supports table, optional sha256, and optional "
            "frame_artifact/frame_sha256/vector_key"
        )
    if "table" not in specification:
        raise ValueError("required_lattice_orbits: missing table")
    table = Path(str(specification["table"]))
    if not table.is_file():
        raise ValueError(f"required_lattice_orbits: missing table {table}")
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    if "sha256" in specification and str(specification["sha256"]) != digest:
        raise ValueError("required_lattice_orbits: table SHA-256 does not match")
    vector_key = str(specification.get("vector_key", "pinned_rank17_w"))
    gram = load_integral_gram(RANK17_GRAM, rank=17)
    frame_artifact = specification.get("frame_artifact")
    frame_digest = None
    if frame_artifact is not None:
        frame_path = Path(str(frame_artifact))
        if not frame_path.is_file():
            raise ValueError(f"required_lattice_orbits: missing frame artifact {frame_path}")
        frame_digest = hashlib.sha256(frame_path.read_bytes()).hexdigest()
        if "frame_sha256" in specification and str(specification["frame_sha256"]) != frame_digest:
            raise ValueError("required_lattice_orbits: frame artifact SHA-256 does not match")
        try:
            frame_rows = json.loads(frame_path.read_text())["rootless_frame"]
        except (json.JSONDecodeError, KeyError, TypeError) as error:
            raise ValueError("required_lattice_orbits: invalid rootless frame artifact") from error
        gram = tuple(tuple(int(value) for value in row) for row in frame_rows)
        if len(gram) != 17 or any(len(row) != 17 for row in gram):
            raise ValueError("required_lattice_orbits: frame artifact must contain a 17 by 17 frame")
        if any(gram[left][right] != gram[right][left] for left in range(17) for right in range(17)):
            raise ValueError("required_lattice_orbits: frame artifact is not symmetric")
        if "vector_key" not in specification:
            vector_key = "alternate_rank17_w"
    with table.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or {"orbit_mask", vector_key} - set(rows[0]):
        raise ValueError("required_lattice_orbits: invalid orbit TSV")
    expected = {int(row["orbit_mask"], 0) for row in rows}
    if len(expected) != len(rows):
        raise ValueError("required_lattice_orbits: duplicate orbit_mask in TSV")
    representatives = {
        int(row["orbit_mask"], 0): parse_integer_vector(
            row[vector_key],
            name=f"required_lattice_orbits[{row['orbit_mask']}].{vector_key}",
            length=17,
        )
        for row in rows
    }
    supplied: list[int] = []
    for record in records:
        if "lattice_orbit_mask" not in record:
            raise ValueError(f"{record['label']}: missing lattice_orbit_mask")
        orbit = int(str(record["lattice_orbit_mask"]), 0)
        supplied.append(orbit)
        if vector_key not in record:
            raise ValueError(f"{record['label']}: missing {vector_key}")
        vector = parse_integer_vector(
            record[vector_key],
            name=f"{record['label']}.{vector_key}",
            length=17,
        )
        vector_norm = integral_norm(vector, gram)
        if vector_norm < 10 or vector_norm % 4 != 2:
            raise ValueError(
                f"{record['label']}: {vector_key} is not a section-nonnegative "
                "degree-two (-2)-class vector"
            )
        if orbit in representatives and any(
            (value - representative) % 2
            for value, representative in zip(vector, representatives[orbit])
        ):
            raise ValueError(
                f"{record['label']}: {vector_key} is not in its declared "
                "section-translation orbit modulo 2M"
            )
    supplied_set = set(supplied)
    if len(supplied_set) != len(supplied):
        raise ValueError("required_lattice_orbits: duplicate lattice_orbit_mask in bisections")
    missing = sorted(expected - supplied_set)
    unexpected = sorted(supplied_set - expected)
    if missing or unexpected:
        raise ValueError(
            "required_lattice_orbits: coverage mismatch "
            f"(missing={len(missing)}, unexpected={len(unexpected)})"
        )
    return {
        "status": "COMPLETE_EXACT_ORBIT_COVERAGE",
        "table": str(table),
        "table_sha256": digest,
        "frame_artifact": None if frame_artifact is None else str(frame_artifact),
        "frame_artifact_sha256": frame_digest,
        "vector_key": vector_key,
        "orbit_count": len(expected),
        "record_lattice_vectors": (
            "Every equation-level record supplied the declared rank-17 vector of "
            "norm 2 modulo 4, at least 10, congruent modulo 2M to its "
            "declared canonical orbit."
        ),
    }


def analyze(
    payload: dict[str, Any], *, require_collision_heights: bool = False,
    require_rank_at_least: int | None = None, compact: bool = False,
) -> dict[str, Any]:
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"expected schema {SCHEMA}")
    variable_name = str(payload.get("base_parameter", "t"))
    variable = sp.Symbol(variable_name)
    records = payload.get("bisections")
    if not isinstance(records, list) or not records:
        raise ValueError("bisections must be a nonempty list")
    labels = [str(record.get("label", "")) for record in records]
    if not all(labels) or len(labels) != len(set(labels)):
        raise ValueError("bisection labels must be nonempty and distinct")
    height_gram = parse_matrix(payload["twist_height_gram"]) if "twist_height_gram" in payload else None
    geometric_lifts = payload.get("rootless_smooth_base_change")
    invariant_mw_rank = (
        parse_nonnegative_integer(payload["invariant_mw_rank"], name="invariant_mw_rank")
        if "invariant_mw_rank" in payload
        else None
    )
    orbit_coverage = verify_lattice_orbit_coverage(payload, records)

    # A section translation is an automorphism over the old base.  Thus two
    # supplied representatives of one M/2M class are the same quadratic cover
    # even before their equations are normalized; retaining both can never
    # give two anti-invariant directions.  Reject this globally, rather than
    # only after an equal-squareclass bucket happens to expose the duplicate.
    declared_orbits: dict[int, str] = {}
    for record in records:
        if "lattice_orbit_mask" not in record:
            continue
        orbit = int(str(record["lattice_orbit_mask"]), 0)
        label = str(record["label"])
        if orbit in declared_orbits:
            raise ValueError(
                "multiple bisection records use section-translation orbit {} "
                "({} and {})".format(orbit, declared_orbits[orbit], label)
            )
        declared_orbits[orbit] = label

    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    keys: dict[str, dict[str, Any]] = {}
    normalized_records: dict[str, dict[str, Any]] = {}
    record_extension_digests: dict[str, str] = {}
    for record in records:
        branch, provenance = branch_from_record(record, variable)
        key = extension_key(branch, variable)
        if is_trivial_squareclass(key):
            raise ValueError(f"{record['label']}: the quadratic cover is split over QQ({variable_name})")
        branch_divisor = geometric_branch_divisor(branch, variable)
        if branch_divisor["geometric_degree"] != 2:
            raise ValueError(
                f"{record['label']}: a smooth rational bisection must have "
                "geometric quadratic branch degree two, got {}".format(
                    branch_divisor["geometric_degree"]
                )
            )
        digest = key_digest(key)
        keys[digest] = key
        buckets[digest].append(record)
        record_extension_digests[str(record["label"])] = digest
        normalized_records[str(record["label"])] = {
            **provenance,
            "geometric_branch_divisor": branch_divisor,
        }

    groups = []
    for digest in sorted(buckets):
        group = sorted(buckets[digest], key=lambda record: str(record["label"]))
        declared_orbits = [record.get("lattice_orbit_mask") for record in group]
        if any(orbit is not None for orbit in declared_orbits):
            if any(orbit is None for orbit in declared_orbits):
                raise ValueError(
                    "an extension bucket mixes records with and without lattice_orbit_mask"
                )
            orbit_masks = [int(str(orbit), 0) for orbit in declared_orbits]
            if len(set(orbit_masks)) != len(orbit_masks):
                raise ValueError(
                    "an extension bucket contains multiple representatives of the same "
                    "section-translation orbit"
                )
        else:
            orbit_masks = None
        groups.append({
            "extension_sha256": digest,
            "extension_squareclass": keys[digest],
            "bisection_labels": [str(record["label"]) for record in group],
            "translation_orbit_masks": orbit_masks,
            "branch_provenance": {
                str(record["label"]): normalized_records[str(record["label"])]
                for record in group
            },
            "bisection_count": len(group),
            "collision": len(group) >= 2,
            "anti_invariant_height": (
                collision_height_data(group, height_gram, geometric_lifts)
                if len(group) >= 2 else None
            ),
        })
    collisions = [group for group in groups if group["collision"]]
    if require_collision_heights:
        missing = [
            group["extension_sha256"]
            for group in collisions
            if not group["anti_invariant_height"]["status"].startswith("COMPUTED_FROM_")
        ]
        if missing:
            raise ValueError(
                "collision heights are required, but {} collision bucket(s) lack "
                "a declared twist height lattice/coordinates or certified geometric "
                "lift intersections".format(
                    len(missing)
                )
            )
    for group in collisions:
        height = group["anti_invariant_height"]
        group["anti_invariant_rank"] = (
            height["rank"]
            if height["status"].startswith("COMPUTED_FROM_")
            else None
        )
        group["base_changed_rank_lower_bound"] = (
            invariant_mw_rank + group["anti_invariant_rank"]
            if invariant_mw_rank is not None and group["anti_invariant_rank"] is not None
            else None
        )
    if require_rank_at_least is not None:
        if require_rank_at_least < 0:
            raise ValueError("require_rank_at_least must be nonnegative")
        achieved = [
            group["base_changed_rank_lower_bound"]
            for group in collisions
            if group["base_changed_rank_lower_bound"] is not None
        ]
        if not achieved or max(achieved) < require_rank_at_least:
            raise ValueError(
                "no collision has the declared invariant rank plus certified "
                f"anti-invariant height rank required to reach {require_rank_at_least}"
            )
    extension_manifest = [
        {
            "label": str(record["label"]),
            "lattice_orbit_mask": record.get("lattice_orbit_mask"),
            "extension_sha256": record_extension_digests[str(record["label"])],
        }
        for record in records
    ]
    manifest_digest = hashlib.sha256(
        json.dumps(extension_manifest, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_EXTENSION_CANONICALIZATION",
        "base_parameter": variable_name,
        "input_bisection_count": len(records),
        "distinct_quadratic_extensions": len(groups),
        "collision_count": len(collisions),
        "collision_height_requirement": bool(require_collision_heights),
        "invariant_mw_rank": invariant_mw_rank,
        "required_base_changed_rank_lower_bound": require_rank_at_least,
        "lattice_orbit_coverage": orbit_coverage,
        "collisions": collisions,
        "all_extension_groups": collisions if compact else groups,
        "compact_output": compact,
        "extension_manifest": extension_manifest if compact else None,
        "extension_manifest_sha256": manifest_digest,
        "proof_boundary": (
            "An equal squareclass is an exact equality of quadratic extensions of QQ(t). "
            "Every accepted record has exactly two geometric branch points, as required "
            "for a smooth rational bisection. "
            "A rank conclusion requires equation-level anti-invariant sections and a "
            "height matrix; any such matrix here is computed only from declared "
            "twist coordinates or certified smooth-rootless lift intersections. "
            "When invariant_mw_rank is declared, the lower bound adds it to the "
            "certified anti-invariant height rank, using the rational +/- eigenspace "
            "decomposition for the quadratic deck involution."
        ),
    }


def self_test_payload() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "base_parameter": "t",
        "invariant_mw_rank": 17,
        "twist_height_gram": [["2", "0"], ["0", "2"]],
        "bisections": [
            {
                "label": "B1",
                "branch": {"numerator_coefficients": ["-2", "2"], "denominator_coefficients": ["1"]},
                "anti_invariant_coordinates": ["1", "0"],
            },
            {
                "label": "B2",
                "quadratic_cover": {
                    "leading_coefficients": ["1"],
                    "linear_coefficients": ["0"],
                    "constant_coefficients": ["2", "-2"],
                },
                "anti_invariant_coordinates": ["0", "1"],
            },
            {
                "label": "B3",
                "branch": {"numerator_coefficients": ["0", "1"], "denominator_coefficients": ["1"]},
            },
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--compact", action="store_true",
        help=(
            "emit the complete label/orbit/digest manifest but omit singleton "
            "extension-group provenance; collision groups remain complete"
        ),
    )
    parser.add_argument(
        "--require-collision-heights",
        action="store_true",
        help="reject every hash collision without declared anti-invariant height data",
    )
    parser.add_argument(
        "--require-rank-at-least",
        type=int,
        help=(
            "require a collision whose declared invariant Mordell--Weil rank plus "
            "certified anti-invariant height rank reaches this lower bound"
        ),
    )
    arguments = parser.parse_args()
    if arguments.self_test == (arguments.input is not None):
        parser.error("supply exactly one of --input or --self-test")
    payload = self_test_payload() if arguments.self_test else json.loads(arguments.input.read_text())
    result = analyze(
        payload,
        require_collision_heights=arguments.require_collision_heights,
        require_rank_at_least=arguments.require_rank_at_least,
        compact=arguments.compact,
    )
    if arguments.input is not None:
        result["input"] = {
            "path": str(arguments.input),
            "sha256": hashlib.sha256(arguments.input.read_bytes()).hexdigest(),
        }
    if arguments.self_test:
        assert result["collision_count"] == 1
        collision = result["collisions"][0]
        assert collision["bisection_labels"] == ["B1", "B2"]
        assert collision["branch_provenance"]["B2"]["kind"] == "quadratic_cover_discriminant"
        assert collision["branch_provenance"]["B1"]["geometric_branch_divisor"]["geometric_degree"] == 2
        assert collision["branch_provenance"]["B1"]["geometric_branch_divisor"]["infinity"]
        # Common factors of a declared rational function do not alter either
        # its squareclass or its branch divisor.
        t = sp.Symbol("t")
        unreduced = {
            "numerator_coefficients": ["-2", "0", "2"],
            "denominator_coefficients": ["1", "1"],
        }
        assert extension_key(unreduced, t) == extension_key(
            {"numerator_coefficients": ["-2", "2"], "denominator_coefficients": ["1"]}, t
        )
        assert geometric_branch_divisor(unreduced, t)["geometric_degree"] == 2
        assert collision["anti_invariant_height"]["matrix"] == [["2", "0"], ["0", "2"]]
        assert collision["anti_invariant_height"]["rank"] == 2
        assert collision["anti_invariant_rank"] == 2
        assert collision["base_changed_rank_lower_bound"] == 19
        split_payload = {
            "schema": SCHEMA,
            "bisections": [{
                "label": "split",
                "quadratic_cover": {
                    "leading_coefficients": ["1"],
                    "linear_coefficients": ["0"],
                    "constant_coefficients": ["-1"],
                },
            }],
        }
        try:
            analyze(split_payload)
        except ValueError as error:
            assert "split over QQ(t)" in str(error)
        else:
            raise AssertionError("split quadratic cover was accepted")
        reduced_chord_payload = {
            "schema": SCHEMA,
            "bisections": [{
                "label": "resolved-chord",
                "quadratic_cover": {
                    "leading_coefficients": ["1"],
                    "linear_coefficients": ["0"],
                    "constant_coefficients": ["-1/4", "-1/2", "-1/2", "-1/2", "-1/4"],
                },
                "trace_section": {"h_coefficients": ["1", "1"]},
                "residual_chord": {"q_coefficients": ["1", "0", "1"]},
            }],
        }
        reduced_chord = analyze(reduced_chord_payload, compact=True)
        assert reduced_chord["distinct_quadratic_extensions"] == 1
        assert reduced_chord["compact_output"] is True
        assert len(reduced_chord["extension_manifest"]) == 1
        bad_reduction = json.loads(json.dumps(reduced_chord_payload))
        bad_reduction["bisections"][0]["residual_chord"]["q_coefficients"][0] = "2"
        try:
            analyze(bad_reduction)
        except ValueError as error:
            assert "discriminant=h^2*q" in str(error)
        else:
            raise AssertionError("an incorrect resolved-chord square reduction was accepted")
        genus_one_payload = {
            "schema": SCHEMA,
            "bisections": [{
                "label": "four-branch-points",
                "branch": {
                    "numerator_coefficients": ["0", "2", "-3", "1"],
                    "denominator_coefficients": ["1"],
                },
            }],
        }
        try:
            analyze(genus_one_payload)
        except ValueError as error:
            assert "branch degree two, got 4" in str(error)
        else:
            raise AssertionError("a genus-one double cover was accepted as a bisection")
        duplicate_orbit_payload = self_test_payload()
        duplicate_orbit_payload["bisections"][0]["lattice_orbit_mask"] = "0x0e"
        duplicate_orbit_payload["bisections"][1]["lattice_orbit_mask"] = 14
        try:
            analyze(duplicate_orbit_payload)
        except ValueError as error:
            assert "section-translation orbit" in str(error)
        else:
            raise AssertionError("a duplicate translation orbit was accepted as a collision")
        distinct_branch_duplicate_orbit = self_test_payload()
        distinct_branch_duplicate_orbit["bisections"][0]["lattice_orbit_mask"] = 14
        distinct_branch_duplicate_orbit["bisections"][1]["lattice_orbit_mask"] = 14
        distinct_branch_duplicate_orbit["bisections"][1]["branch"] = {
            "numerator_coefficients": ["0", "1"], "denominator_coefficients": ["1"],
        }
        distinct_branch_duplicate_orbit["bisections"][1].pop("quadratic_cover")
        try:
            analyze(distinct_branch_duplicate_orbit)
        except ValueError as error:
            assert "multiple bisection records use section-translation orbit" in str(error)
        else:
            raise AssertionError("a duplicate translation orbit with inconsistent branch data was accepted")
        with TemporaryDirectory() as directory:
            table = Path(directory) / "orbits.tsv"
            vector_14 = "0 0 -1 0 -1 0 -1 0 0 0 -1 0 2 -1 -1 0 1"
            vector_53 = "0 -1 0 2 -1 0 0 -1 2 -1 0 1 0 1 0 -2 -1"
            table.write_text(
                "orbit_mask\thex\tpinned_rank17_w\n"
                f"14\t0x0000e\t{vector_14}\n"
                f"53\t0x00035\t{vector_53}\n"
            )
            covered_payload = {
                "schema": SCHEMA,
                "required_lattice_orbits": {
                    "table": str(table),
                    "sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                },
                "bisections": [
                    {
                        "label": "orbit-14",
                        "lattice_orbit_mask": 14,
                        "pinned_rank17_w": vector_14,
                        "branch": {
                            "numerator_coefficients": ["-2", "2"],
                            "denominator_coefficients": ["1"],
                        },
                    },
                    {
                        "label": "orbit-53",
                        "lattice_orbit_mask": "0x35",
                        "pinned_rank17_w": vector_53.split(),
                        "branch": {
                            "numerator_coefficients": ["0", "1"],
                            "denominator_coefficients": ["1"],
                        },
                    },
                ],
            }
            covered = analyze(covered_payload)
            assert covered["lattice_orbit_coverage"]["status"] == "COMPLETE_EXACT_ORBIT_COVERAGE"
            shifted_payload = dict(covered_payload)
            shifted_payload["bisections"] = [dict(record) for record in covered_payload["bisections"]]
            shifted_payload["bisections"][0]["pinned_rank17_w"] = vector_53
            try:
                analyze(shifted_payload)
            except ValueError as error:
                assert "section-translation orbit" in str(error)
            else:
                raise AssertionError("a vector from a different translation orbit was accepted")
        with TemporaryDirectory() as directory:
            # Exercise the alternate-frame schema without depending on large,
            # ignored generated q80 artifacts.  The matrix need only be a
            # declared exact rank-17 frame for this parser/attachment test;
            # the actual alternate q80 enumeration has its own Sage replay.
            alternate_frame = Path(directory) / "alternate-frame.json"
            pinned_gram = load_integral_gram(RANK17_GRAM, rank=17)
            alternate_frame.write_text(json.dumps({
                "rootless_frame": [list(row) for row in pinned_gram],
            }))
            alternate_rows = [
                {"orbit_mask": "14", "hex": "0x0000e", "alternate_rank17_w": vector_14},
                {"orbit_mask": "53", "hex": "0x00035", "alternate_rank17_w": vector_53},
            ]
            table = Path(directory) / "alternate-orbits.tsv"
            table.write_text(
                "orbit_mask\thex\talternate_rank17_w\n"
                +"\n".join(
                    f"{row['orbit_mask']}\t{row['hex']}\t{row['alternate_rank17_w']}"
                    for row in alternate_rows
                )+"\n"
            )
            alternate_payload = {
                "schema": SCHEMA,
                "required_lattice_orbits": {
                    "table": str(table),
                    "sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                    "frame_artifact": str(alternate_frame),
                    "frame_sha256": hashlib.sha256(alternate_frame.read_bytes()).hexdigest(),
                    "vector_key": "alternate_rank17_w",
                },
                "bisections": [
                    {
                        "label": f"alternate-{index}",
                        "lattice_orbit_mask": row["orbit_mask"],
                        "alternate_rank17_w": row["alternate_rank17_w"],
                        "branch": {
                            "numerator_coefficients": ["-2", "2"] if index == 0 else ["0", "1"],
                            "denominator_coefficients": ["1"],
                        },
                    }
                    for index, row in enumerate(alternate_rows)
                ],
            }
            alternate_covered = analyze(alternate_payload)
            assert alternate_covered["lattice_orbit_coverage"] == {
                "status": "COMPLETE_EXACT_ORBIT_COVERAGE",
                "table": str(table),
                "table_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                "frame_artifact": str(alternate_frame),
                "frame_artifact_sha256": hashlib.sha256(alternate_frame.read_bytes()).hexdigest(),
                "vector_key": "alternate_rank17_w",
                "orbit_count": 2,
                "record_lattice_vectors": (
                    "Every equation-level record supplied the declared rank-17 vector of "
                    "norm 2 modulo 4, at least 10, congruent modulo 2M to its "
                    "declared canonical orbit."
                ),
            }
            shifted_payload = dict(alternate_payload)
            shifted_payload["bisections"] = [
                dict(record) for record in alternate_payload["bisections"]
            ]
            shifted_payload["bisections"][0]["alternate_rank17_w"] = alternate_rows[1]["alternate_rank17_w"]
            try:
                analyze(shifted_payload)
            except ValueError as error:
                assert "section-translation orbit" in str(error)
            else:
                raise AssertionError("an alternate-frame vector from a different orbit was accepted")
        missing_height_payload = {
            "schema": SCHEMA,
            "bisections": [
                {
                    "label": "missing-height-1",
                    "branch": {"numerator_coefficients": ["0", "1"], "denominator_coefficients": ["1"]},
                },
                {
                    "label": "missing-height-2",
                    "branch": {"numerator_coefficients": ["0", "1"], "denominator_coefficients": ["1"]},
                },
            ],
        }
        try:
            analyze(missing_height_payload, require_collision_heights=True)
        except ValueError as error:
            assert "collision heights are required" in str(error)
        else:
            raise AssertionError("heightless collision was accepted in strict mode")
        try:
            analyze(self_test_payload(), require_rank_at_least=20)
        except ValueError as error:
            assert "required to reach 20" in str(error)
        else:
            raise AssertionError("insufficient rank collision was accepted")
        indefinite_height_payload = self_test_payload()
        indefinite_height_payload["twist_height_gram"] = [["2", "0"], ["0", "-2"]]
        try:
            analyze(indefinite_height_payload)
        except ValueError as error:
            assert "positive definite" in str(error)
        else:
            raise AssertionError("an indefinite declared twist height lattice was accepted")
        geometric_height_payload = self_test_payload()
        geometric_height_payload.pop("twist_height_gram")
        for record in geometric_height_payload["bisections"]:
            record.pop("anti_invariant_coordinates", None)
        geometric_height_payload["rootless_smooth_base_change"] = {
            "base_change_chi": "4",
            "branch_fibres_smooth": True,
            "base_change_rootless": True,
            "same_sheet_intersections": {
                "B1|B1": "-4", "B1|B2": "0", "B2|B2": "-4",
            },
            "conjugate_sheet_intersections": {
                "B1|B1": "2", "B1|B2": "2", "B2|B2": "2",
            },
        }
        geometric = analyze(
            geometric_height_payload, require_collision_heights=True,
            require_rank_at_least=19,
        )
        geometric_height = geometric["collisions"][0]["anti_invariant_height"]
        assert geometric_height["status"] == "COMPUTED_FROM_GEOMETRIC_LIFT_INTERSECTIONS"
        assert geometric_height["matrix"] == [["12", "4"], ["4", "12"]]
        assert geometric_height["rank"] == 2
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "BISECTIONEXT|bisections={}|extensions={}|collisions={}|status={}".format(
            result["input_bisection_count"], result["distinct_quadratic_extensions"],
            result["collision_count"], result["status"],
        )
    )


if __name__ == "__main__":
    main()
