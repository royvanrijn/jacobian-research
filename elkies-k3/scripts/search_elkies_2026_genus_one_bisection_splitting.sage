#!/usr/bin/env sage-python
"""High-throughput simultaneous splitting search for R17 genus-one bisections.

The search uses a diversity-greedy sample of equation-cheap norm-eight trace
classes and every one of the 43 norm-twelve deep classes.  Norm-eight traces
carry the pencil ``M=M0+lambda*h^2``.  A neutral specialization of one of the
generic R17 sections supplies a point on a member; subsequent parameter
samples come from multiples on that pointed binary quartic's Jacobian.  The
norm-twelve deep classes have ``deg(h)=4`` and their unique regular member
already has a quartic branch polynomial.

The C++ hot loop intersects dynamically sized Legendre-symbol bitmasks in
three or more disjoint prime blocks.  It ranks parameters by the weakest-block
surviving cover count and F2-rank of the trace masks.  Only parameters with at
least two all-block survivors receive exact square tests.  Every exact
simultaneous split is mapped to points on the specialized R17 curve and then
reduced against the specialized generic MW17 subgroup by exact finite-group
certificates.

This is a bounded search.  Modular rankings are heuristics; positive exact
square identities and positive finite-quotient escape are certificates, while
failure to escape the chosen finite quotients is not a dependence proof.
"""

from __future__ import annotations

import argparse
from array import array
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, prime_range, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3/scripts"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
IDENTIFICATION = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
RANK_SCRIPT = SCRIPTS / "rank_elkies_2026_bisection_orbits.sage"
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"
PILOT_SCRIPT = SCRIPTS / "search_elkies_2026_rank28_genus_one_bisections.sage"
DIVERSITY_SCRIPT = SCRIPTS / "analyze_r17_multisection_diversity.py"
SCANNER = SCRIPTS / "scan_elkies_2026_genus_one_bisection_splitting.cpp"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
LOCAL = ROOT / "artifacts/local/elkies-k3/r17-genus-one-bisection-splitting"
MAGIC = 0x47315332


def load_script(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def polynomial_record(polynomial, width: int | None = None) -> list[str]:
    size = width if width is not None else polynomial.degree() + 1
    return [rational_text(polynomial[index]) for index in range(size)]


def f2_rank(values) -> int:
    pivots = {}
    for raw in values:
        value = int(raw)
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def select_diverse_norm_eight(ranked: list[dict], count: int, pool_size: int) -> list[dict]:
    """Greedy maximin selection inside a hard equation-cost prefix."""

    pool = ranked[: min(pool_size, len(ranked))]
    if count > len(pool):
        raise ValueError("norm-eight count exceeds the equation-cheap pool")
    chosen = [pool[0]]
    remaining = pool[1:]
    span_rank = f2_rank([chosen[0]["orbit_mask"]])
    support_signatures = {
        sum(1 << (index - 1) for index in chosen[0]["support_one_based"])
    }
    while len(chosen) < count:
        chosen_masks = [row["orbit_mask"] for row in chosen]
        best = None
        best_key = None
        for row in remaining:
            mask = int(row["orbit_mask"])
            candidate_rank = f2_rank([*chosen_masks, mask])
            support = sum(1 << (index - 1) for index in row["support_one_based"])
            key = (
                candidate_rank > span_rank,
                support not in support_signatures,
                min((mask ^ prior).bit_count() for prior in chosen_masks),
                -int(row["equation_rank"]),
            )
            if best_key is None or key > best_key:
                best, best_key = row, key
        assert best is not None
        chosen.append(best)
        remaining.remove(best)
        span_rank = f2_rank(row["orbit_mask"] for row in chosen)
        support_signatures.add(
            sum(1 << (index - 1) for index in best["support_one_based"])
        )
    return chosen


def enumerate_deep_classes(diversity_module, rank_module) -> tuple[list[dict], dict]:
    gram = diversity_module.load_gram(SHORT_GRAM)
    oracle = diversity_module.CosetOracle(gram, 2)
    short_coordinates = rank_module.load_matrix(SHORT_COORDS)
    pinned = rank_module.load_matrix(PINNED)
    identification = json.loads(IDENTIFICATION.read_text())
    basis_change = matrix(ZZ, identification["pinned_identification"]["basis_change_matrix"])
    short_to_published = short_coordinates * basis_change.transpose()
    section_data = json.loads(SECTIONS.read_text())
    direct_costs, closures = rank_module.section_input_costs(section_data)
    records = []
    maximum_error = 0.0
    histogram = Counter()
    for mask in range(1 << 17):
        norm, representative, error = oracle.solve(diversity_module.decode_residue(mask, 2))
        histogram[norm] += 1
        maximum_error = max(maximum_error, error)
        if norm != 12:
            continue
        short = vector(ZZ, representative)
        key, score, published, sign = rank_module.canonical_scored_orientation(
            short, short_to_published, direct_costs, closures
        )
        short *= sign
        pinned_vector = short * short_coordinates
        if pinned_vector * pinned * pinned_vector != 12:
            raise ArithmeticError("deep representative has the wrong pinned norm")
        records.append(
            {
                "equation_rank": None,
                "orbit_mask": mask,
                "orbit_hex": f"0x{mask:05x}",
                "published_basis_w": tuple(int(value) for value in published),
                "pinned_rank17_w": tuple(int(value) for value in pinned_vector),
                "short_basis_w": tuple(int(value) for value in short),
                **score,
            }
        )
    if len(records) != 43:
        raise ArithmeticError(f"expected 43 norm-twelve deep classes, got {len(records)}")
    records.sort(key=lambda row: (row["group_addition_upper_bound"], row["coordinate_input_bits"], row["orbit_mask"]))
    for index, row in enumerate(records, start=1):
        row["deep_equation_rank"] = index
    return records, {
        "complete_coset_count": 1 << 17,
        "minimum_norm_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "deep_class_count": len(records),
        "maximum_double_double_integrality_error": maximum_error,
    }


def rational_square_root(value):
    value = QQ(value)
    if value < 0:
        return None
    numerator = math.isqrt(int(value.numerator()))
    denominator = math.isqrt(int(value.denominator()))
    if numerator * numerator != value.numerator() or denominator * denominator != value.denominator():
        return None
    return QQ(numerator) / denominator


def pointed_quartic_jacobian(q, t0, s0, multiple_bound: int):
    z_ring = PolynomialRing(QQ, "z")
    z = z_ring.gen()
    shifted = z_ring(q(z + t0))
    e, d, c, b, a = (QQ(shifted[index]) for index in range(5))
    v0 = QQ(s0)
    if not v0 or v0**2 != e:
        raise ArithmeticError("invalid pointed-quartic base point")
    curve = EllipticCurve(QQ, [
        d / v0,
        c - d**2 / (4 * v0**2),
        2 * v0 * b,
        -4 * v0**2 * a,
        a * (d**2 - 4 * v0**2 * c),
    ])
    opposite_x = d**2 / (4 * v0**2) - c
    points = curve.lift_x(opposite_x, all=True)
    generator = next((point for point in points if not point.has_finite_order()), None)
    if generator is None:
        return curve, None, []
    samples = []
    current = generator
    for multiple in range(2, multiple_bound + 1):
        current += generator
        if current.is_zero() or current[1] == 0:
            continue
        x_value, y_value = current[:2]
        local_t = (4 * v0**2 * (x_value + c) - d**2) / (2 * v0 * y_value)
        if not local_t:
            continue
        quartic_s = (x_value * local_t**2 - d * local_t) / (2 * v0) - v0
        parameter = local_t + t0
        if quartic_s**2 != q(parameter):
            raise ArithmeticError("pointed quartic inverse map failed")
        samples.append((multiple, parameter, quartic_s))
    return curve, generator, samples


def fit_neutral_norm_eight_member(
    frame, trace_point, basis, A, B, delta, ring, field, seed_height,
    jacobian_multiples, seed_salt, forbidden_seed_parameters,
):
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("selected norm-eight trace is not in the finite-pole chart")
    candidates = []
    for denominator in range(1, seed_height + 1):
        for numerator in range(-seed_height, seed_height + 1):
            if math.gcd(abs(numerator), denominator) == 1:
                candidates.append(QQ(numerator) / denominator)
    # Avoid manufacturing a high-multiplicity collision at one common seed.
    # The rotation is deterministic and depends only on the lattice class.
    offset = seed_salt % len(candidates)
    candidates = candidates[offset:] + candidates[:offset]
    for seed_index, t0 in enumerate(candidates):
        if t0 in forbidden_seed_parameters:
            continue
        if not h(t0) or not delta(t0):
            continue
        target_index = (seed_index + seed_salt) % len(basis)
        target = basis[target_index]
        target_x, target_y = QQ(target[0](t0)), QQ(target[1](t0))
        trace_x, trace_y = QQ(trace_point[0](t0)), QQ(trace_point[1](t0))
        if target_x == trace_x:
            continue
        target_slope = (target_y + trace_y) / (target_x - trace_x)
        pencil_parameter = (target_slope * h(t0) - M0(t0)) / h(t0)**2
        M = M0 + pencil_parameter * h**2
        try:
            data = CHORD_MODULE.chord_data_from_slope_numerator(
                h, Nx, Ny, M, A, B, delta, ring, field, expected_q_degree=4
            )
        except ArithmeticError:
            continue
        q = data["q"]
        if not q.is_irreducible() or q.gcd(delta).degree() or q.gcd(h).degree():
            continue
        s0 = (2 * target_x - data["sum_x"](t0)) / h(t0)
        if not s0 or s0**2 != q(t0):
            continue
        curve, generator, samples = pointed_quartic_jacobian(q, t0, s0, jacobian_multiples)
        if generator is None:
            continue
        forbidden_seed_parameters.add(t0)
        return data, {
            "method": "generic-R17-section seed followed by pointed-quartic Jacobian multiples",
            "generic_section_index_one_based": target_index + 1,
            "t": rational_text(t0),
            "s": rational_text(s0),
            "lambda": rational_text(pencil_parameter),
        }, curve, generator, samples
    raise ArithmeticError("no neutral pointed norm-eight member found in the declared seed box")


def search_deep_point(q, height: int):
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if math.gcd(abs(numerator), denominator) != 1:
                continue
            t0 = QQ(numerator) / denominator
            s0 = rational_square_root(q(t0))
            if s0:
                return t0, s0
    return None


def trace_record(
    row, norm, basis, generic_curve, A, B, delta, ring, field, seed_height,
    deep_point_height, jacobian_multiples, forbidden_seed_parameters,
):
    published = vector(ZZ, row["published_basis_w"])
    trace = sum((coefficient * point for coefficient, point in zip(published, basis)), generic_curve(0))
    frame = CHORD_MODULE.trace_chord_frame(trace[0], trace[1], ring)
    if norm == 8:
        data, seed, jacobian, generator, samples = fit_neutral_norm_eight_member(
            frame, trace, basis, A, B, delta, ring, field, seed_height,
            jacobian_multiples, int(row["orbit_mask"]), forbidden_seed_parameters,
        )
        member = "one rational member of M=M0+lambda*h^2"
    else:
        if frame["h"].degree() != 4:
            raise ArithmeticError("norm-twelve deep trace does not have four finite poles")
        data = CHORD_MODULE.chord_data_from_slope_numerator(
            frame["h"], frame["Nx"], frame["Ny"], frame["M0"],
            A, B, delta, ring, field, expected_q_degree=4,
        )
        q = data["q"]
        if not q.is_irreducible() or q.gcd(delta).degree() or q.gcd(frame["h"]).degree():
            raise ArithmeticError("invalid norm-twelve deep quartic")
        found = search_deep_point(q, deep_point_height)
        if found is None:
            seed, jacobian, generator, samples = None, None, None, []
        else:
            t0, s0 = found
            jacobian, generator, samples = pointed_quartic_jacobian(q, t0, s0, jacobian_multiples)
            seed = {"method": "bounded primitive point search", "t": rational_text(t0), "s": rational_text(s0)}
        member = "unique regular M0 member"
    q = data["q"]
    record = {
        "label": f"norm{norm}-orbit-{int(row['orbit_mask']):05x}",
        "trace_norm": norm,
        "lattice_orbit_mask": int(row["orbit_mask"]),
        "published_basis_w": list(map(int, row["published_basis_w"])),
        "pinned_rank17_w": list(map(int, row["pinned_rank17_w"])),
        "equation_complexity": {
            key: row.get(key) for key in (
                "equation_rank", "deep_equation_rank", "group_addition_upper_bound",
                "support_count", "dependency_count", "coordinate_input_bits",
                "maximum_absolute_coefficient", "coefficient_l1",
            ) if row.get(key) is not None
        },
        "trace_section": {
            "h_coefficients_low_to_high": polynomial_record(frame["h"]),
            "Nx_coefficients_low_to_high": polynomial_record(frame["Nx"]),
            "Ny_coefficients_low_to_high": polynomial_record(frame["Ny"]),
            "M0_coefficients_low_to_high": polynomial_record(frame["M0"]),
        },
        "member_selection": member,
        "branch_polynomial_q_coefficients_low_to_high": polynomial_record(q, 5),
        "branch_polynomial_irreducible_over_Q": True,
        "branch_polynomial_squarefree": True,
        "branch_polynomial_coprime_to_surface_discriminant": True,
        "branch_polynomial_coprime_to_trace_denominator": True,
        "lifted_section": {
            "cover": "s^2=q(t)",
            "x0_coefficients_low_to_high": polynomial_record(data["x0"]),
            "x1_coefficients_low_to_high": polynomial_record(data["x1"]),
            "y0_coefficients_low_to_high": polynomial_record(data["y0"]),
            "y1_coefficients_low_to_high": polynomial_record(data["y1"]),
            "trace": "P(t,s)+P(t,-s)=tau(t)",
        },
        "rational_seed": seed,
        "pointed_quartic_jacobian": None if jacobian is None else {
            "weierstrass_coefficients": [rational_text(value) for value in jacobian.a_invariants()],
            "generator": [rational_text(value) for value in generator[:2]],
            "generator_nontorsion": True,
            "sample_count": len(samples),
        },
        "jacobian_parameter_samples": [
            {"multiple": multiple, "t": rational_text(t), "s": rational_text(s)}
            for multiple, t, s in samples
        ],
    }
    record["_data"] = data
    return record


def normalized_integer_quartic(record):
    values = [Fraction(value) for value in record["branch_polynomial_q_coefficients_low_to_high"]]
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    coefficients = tuple(value.numerator * (denominator // value.denominator) for value in values)
    return denominator, coefficients


def homogeneous_value(curve, a: int, b: int) -> int:
    unused, coefficients = curve
    return sum(coefficients[index] * a**index * b**(4-index) for index in range(5))


def exact_square_root(curve, a: int, b: int):
    scalar, coefficients = curve
    radicand = scalar * homogeneous_value((1, coefficients), a, b)
    if radicand < 0:
        return None
    root = math.isqrt(radicand)
    return root if root * root == radicand else None


def choose_primes(curves, count: int):
    # If a quartic's clearing scalar vanishes modulo p, that curve receives an
    # all-one local mask at p.  This loses filtering power but never rejects a
    # rational square, and keeps the simultaneous tables compact.
    candidates = []
    for prime in prime_range(19, 500):
        prime = int(prime)
        blind = sum(scalar % prime == 0 for scalar, unused in curves)
        # A blind curve is harmless but contributes no filtering at this
        # place.  Balance that loss against the p^2 table size.
        candidates.append((blind * 1000 + prime * prime // 10, blind, prime))
    selected = sorted(candidates)[:count]
    return [prime for unused, blind, prime in selected]


def export_tables(path: Path, curves, lattice_masks, primes, block_count: int):
    if sys.byteorder != "little":
        raise RuntimeError("the packed C++ table format currently requires little endian")
    word_count = (len(curves) + 63) // 64
    with path.open("wb") as stream:
        stream.write(struct.pack("<IIIII", MAGIC, len(curves), word_count, len(primes), block_count))
        stream.write(array("I", lattice_masks).tobytes())
        for prime_index, prime in enumerate(primes):
            block = prime_index % block_count
            stream.write(struct.pack("<II", prime, block))
            squares = {value * value % prime for value in range(prime)}
            reduced = [(scalar % prime, tuple(value % prime for value in coefficients)) for scalar, coefficients in curves]
            masks = array("Q")
            for a in range(prime):
                for b in range(prime):
                    words = [0] * word_count
                    for index, curve in enumerate(reduced):
                        value = curve[0] * homogeneous_value((1, curve[1]), a, b)
                        if value % prime in squares:
                            words[index // 64] |= 1 << (index % 64)
                    masks.extend(words)
            stream.write(masks.tobytes())
    return word_count


def compile_and_scan(table, scan_output, binary, args):
    compiler = shutil.which("g++")
    if compiler is None:
        raise SystemExit("g++ is required")
    compile_command = [compiler, "-O3", "-std=c++17", "-Wall", "-Wextra", "-pedantic", str(SCANNER), "-o", str(binary)]
    subprocess.run(compile_command, check=True)
    scan_command = [
        str(binary), str(table), str(args.random_count), str(args.coordinate_bound),
        str(args.box_height), str(args.top_k), str(args.random_seed), str(scan_output),
    ]
    completed = subprocess.run(scan_command, check=True, text=True, capture_output=True)
    if "|status=PASS" not in completed.stdout:
        raise ArithmeticError("scanner did not report PASS")
    return compile_command, scan_command, completed.stdout.strip()


def parse_scan(path: Path, word_count: int):
    extremes, rankings = [], []
    for line in path.read_text().splitlines():
        fields = line.split()
        row = {
            "projective_pair": [int(fields[1]), int(fields[2])],
            "weakest_block_lattice_rank": int(fields[3]),
            "weakest_block_cover_count": int(fields[4]),
            "block_rank_sum": int(fields[5]),
            "block_cover_count_sum": int(fields[6]),
            "all_block_survivor_count": int(fields[7]),
        }
        if fields[0] == "E":
            words = [int(value, 16) for value in fields[8:]]
            if len(words) != word_count:
                raise ArithmeticError("scanner returned the wrong mask width")
            row["survivor_words_hex_low_to_high"] = [f"{word:016x}" for word in words]
            row["_words"] = words
            extremes.append(row)
        elif fields[0] == "R":
            rankings.append(row)
        else:
            raise ArithmeticError("unknown scanner record")
    return extremes, rankings


def set_indices(words, count):
    for index in range(count):
        if (words[index // 64] >> (index % 64)) & 1:
            yield index


def finite_quotient_profile(model, generic_points, new_points, labels, relation_primes, prime_bound):
    cas = ROOT / "elliptic-curves/cas"
    if str(cas) not in sys.path:
        sys.path.insert(0, str(cas))
    from elliptic_candidate_record import build_finite_quotient_certificate, verify_finite_quotient_certificate
    from finite_quotient_escape import QuotientBlock, analyze_escape

    def fraction(value):
        return Fraction(rational_text(value))

    rational_model = tuple(fraction(value) for value in model)
    rational_points = tuple((fraction(x), fraction(y)) for x, y in (*generic_points, *new_points))
    attempts = []
    maximum = 0
    basis_labels = []
    torsion_witness_orders = []
    for relation_prime in relation_primes:
        try:
            certificate = build_finite_quotient_certificate(
                rational_model, rational_points, relation_prime=relation_prime, prime_bound=prime_bound
            )
            verify_finite_quotient_certificate(rational_model, rational_points, certificate)
            blocks = tuple(
                QuotientBlock.build(
                    modulus=relation_prime,
                    rows=signature["rows"],
                    column_count=len(rational_points),
                    source=f"good-reduction-p={signature['prime']}",
                ) for signature in certificate["signatures"]
            )
            profile = analyze_escape(blocks, known_column_count=17, candidate_labels=labels)
            record = profile.to_record()
            record["relation_prime"] = relation_prime
            record["status"] = "PASS_EXACT_FINITE_QUOTIENT_PROFILE"
            record["finite_quotient_certificate"] = certificate
            attempts.append(record)
            torsion_witness_orders.append(int(certificate["torsion_witness"]["group_order"]))
            if profile.marginal_dimension > maximum:
                maximum = profile.marginal_dimension
                basis_labels = list(profile.independent_escape_basis_labels)
        except (RuntimeError, ArithmeticError, AssertionError, ValueError) as error:
            attempts.append({"relation_prime": relation_prime, "status": "bounded-certificate-failure", "detail": str(error)})
    torsion_order_divides = 0
    for order in torsion_witness_orders:
        torsion_order_divides = math.gcd(torsion_order_divides, order)
    return {
        "quotient": "specialized points modulo the specialized generic MW17 subgroup",
        "known_column_count": 17,
        "candidate_count": len(new_points),
        "maximum_certified_marginal_dimension": maximum,
        "independent_escape_basis_labels": basis_labels,
        "torsion_order_divides_gcd": torsion_order_divides,
        "candidate_escape_is_nontorsion": bool(maximum and torsion_order_divides == 1),
        "certified_rank_lower_bound": (
            17 + maximum if maximum and torsion_order_divides == 1 else None
        ),
        "relation_prime_profiles": attempts,
        "claim_boundary": "positive finite-quotient escape is exact; zero observed escape is not rational dependence",
    }


def exact_test_extremes(extremes, curves, records, A, B, basis, exact_limit, relation_primes, prime_bound):
    ordering = sorted(
        extremes,
        key=lambda row: (
            row["all_block_survivor_count"], row["weakest_block_lattice_rank"],
            row["weakest_block_cover_count"], row["block_rank_sum"],
        ), reverse=True,
    )
    total_extremes = len(ordering)
    if exact_limit:
        ordering = ordering[:exact_limit]
    exact_tests = 0
    tested_parameters = 0
    hits = []
    seen_pairs = set()
    for row in ordering:
        a, b = row["projective_pair"]
        if (a, b) in seen_pairs:
            continue
        seen_pairs.add((a, b))
        tested_parameters += 1
        t = QQ(a) / b
        split = []
        new_points = []
        labels = []
        for index in set_indices(row["_words"], len(curves)):
            exact_tests += 1
            root = exact_square_root(curves[index], a, b)
            if root is None:
                continue
            scalar, unused = curves[index]
            s = QQ(root) / (scalar * b**2)
            data = records[index]["_data"]
            x = data["x0"](t) + data["x1"](t) * s
            y = data["y0"](t) + data["y1"](t) * s
            if y**2 != x**3 + A(t) * x + B(t):
                raise ArithmeticError("constructed split-cover point missed the R17 fibre")
            labels.append(records[index]["label"])
            new_points.append((x, y))
            split.append({
                "cover_index_zero_based": index,
                "cover_label": records[index]["label"],
                "trace_norm": records[index]["trace_norm"],
                "cover_coordinate_s": rational_text(s),
                "point_on_R17_fibre": [rational_text(x), rational_text(y)],
            })
        if len(split) < 2:
            continue
        generic_points = [(QQ(point[0](t)), QQ(point[1](t))) for point in basis]
        model = (QQ(0), QQ(0), QQ(0), QQ(A(t)), QQ(B(t)))
        quotient = finite_quotient_profile(
            model, generic_points, new_points, labels, relation_primes, prime_bound
        )
        clean_rank = {key: value for key, value in row.items() if not key.startswith("_")}
        hits.append({
            **clean_rank,
            "t": rational_text(t),
            "exact_split_cover_count": len(split),
            "exact_trace_mask_rank": f2_rank(records[item["cover_index_zero_based"]]["lattice_orbit_mask"] for item in split),
            "split_covers": split,
            "generic_MW17_points": [[rational_text(x), rational_text(y)] for x, y in generic_points],
            "quotient_by_generic_MW17": quotient,
        })
    return {
        "modular_extreme_count_available": total_extremes,
        "modular_extremes_exact_tested": tested_parameters,
        "exact_stage_truncated": bool(exact_limit and total_extremes > exact_limit),
        "exact_square_tests": exact_tests,
        "simultaneous_split_hits": hits,
    }


def jacobian_cross_scan(records, curves, primes):
    """Test Jacobian-generated parameters against every other compiled cover."""
    squares = {prime: {value * value % prime for value in range(prime)} for prime in primes}
    seen = set()
    extremes = []
    for source_index, record in enumerate(records):
        for sample in record["jacobian_parameter_samples"]:
            t = Fraction(sample["t"])
            pair = (t.numerator, t.denominator)
            if pair in seen:
                continue
            seen.add(pair)
            survivors = []
            for index, curve in enumerate(curves):
                scalar, coefficients = curve
                good = True
                for prime in primes:
                    value = scalar * homogeneous_value(
                        (1, tuple(coefficient % prime for coefficient in coefficients)),
                        pair[0] % prime, pair[1] % prime,
                    )
                    if value % prime not in squares[prime]:
                        good = False
                        break
                if good:
                    survivors.append(index)
            if len(survivors) >= 2:
                words = [0] * ((len(curves) + 63) // 64)
                for index in survivors:
                    words[index // 64] |= 1 << (index % 64)
                extremes.append({
                    "projective_pair": list(pair),
                    "weakest_block_lattice_rank": f2_rank(records[index]["lattice_orbit_mask"] for index in survivors),
                    "weakest_block_cover_count": len(survivors),
                    "block_rank_sum": 0,
                    "block_cover_count_sum": 0,
                    "all_block_survivor_count": len(survivors),
                    "source": "pointed-quartic-Jacobian population",
                    "_words": words,
                })
    return extremes, len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--norm8-count", type=int, default=100)
    parser.add_argument("--equation-pool-size", type=int, default=1000)
    parser.add_argument("--pari-stack-gb", type=int, default=2)
    parser.add_argument("--seed-height", type=int, default=32)
    parser.add_argument("--deep-point-height", type=int, default=12)
    parser.add_argument("--jacobian-multiples", type=int, default=10)
    parser.add_argument("--prime-count", type=int, default=12)
    parser.add_argument("--prime-block-count", type=int, default=3)
    parser.add_argument("--random-count", type=int, default=1_000_000)
    parser.add_argument("--coordinate-bound", type=int, default=1_000_000_000)
    parser.add_argument("--box-height", type=int, default=2000)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--random-seed", type=int, default=0x17B15EC7)
    parser.add_argument("--exact-limit", type=int, default=0, help="zero exact-tests every modular extreme")
    parser.add_argument("--relation-primes", default="2,3,5")
    parser.add_argument("--reduction-prime-bound", type=int, default=250)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--local-directory", type=Path, default=LOCAL)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not 100 <= args.norm8_count <= 1000:
        parser.error("--norm8-count must lie in the requested interval 100..1000")
    if args.equation_pool_size < args.norm8_count or args.prime_count < args.prime_block_count * 2:
        parser.error("the equation pool or prime ensemble is too small")
    if args.prime_block_count < 3:
        parser.error("at least three disjoint prime blocks are required")
    relation_primes = tuple(int(value) for value in args.relation_primes.split(","))

    global CHORD_MODULE
    rank_module = load_script("r17_rank_bisections_search", RANK_SCRIPT)
    CHORD_MODULE = load_script("r17_chord_search", CHORD_SCRIPT)
    pilot_module = load_script("r17_genus_one_pilot_search", PILOT_SCRIPT)
    diversity_module = load_script("r17_diversity_search", DIVERSITY_SCRIPT)
    ranked, norm8_enumeration = pilot_module.equation_ranked_norm_eight(rank_module, args.pari_stack_gb)
    norm8 = select_diverse_norm_eight(ranked, args.norm8_count, args.equation_pool_size)
    deep, deep_enumeration = enumerate_deep_classes(diversity_module, rank_module)

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    ring = PolynomialRing(QQ, "t")
    field = ring.fraction_field()
    A = ring([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = ring([QQ(value) for value in model["B_coefficients_low_to_high"]])
    delta = ring(-16 * (4 * A**3 + 27 * B**2))
    basis_coordinates = CHORD_MODULE.reconstruct_basis(ring, A, B, section_data)
    generic_curve = EllipticCurve(field, [A, B])
    basis = [generic_curve(field(x), field(y)) for x, y in basis_coordinates]

    records = []
    used_seed_parameters = set()
    for row in norm8:
        records.append(trace_record(
            row, 8, basis, generic_curve, A, B, delta, ring, field,
            args.seed_height, args.deep_point_height, args.jacobian_multiples,
            used_seed_parameters,
        ))
    for row in deep:
        records.append(trace_record(
            row, 12, basis, generic_curve, A, B, delta, ring, field,
            args.seed_height, args.deep_point_height, args.jacobian_multiples,
            used_seed_parameters,
        ))
    curves = [normalized_integer_quartic(record) for record in records]
    primes = choose_primes(curves, args.prime_count)

    args.local_directory.mkdir(parents=True, exist_ok=True)
    table_path = args.local_directory / "legendre-mask-tables.bin"
    scan_path = args.local_directory / "scan-output.txt"
    binary_path = args.local_directory / "scan-genus-one-bisections"
    word_count = export_tables(
        table_path, curves, [record["lattice_orbit_mask"] for record in records],
        primes, args.prime_block_count,
    )
    compile_command, scan_command, scanner_summary = compile_and_scan(
        table_path, scan_path, binary_path, args
    )
    random_extremes, rankings = parse_scan(scan_path, word_count)
    jacobian_extremes, jacobian_population = jacobian_cross_scan(records, curves, primes)
    exact = exact_test_extremes(
        [*random_extremes, *jacobian_extremes], curves, records, A, B,
        basis_coordinates, args.exact_limit, relation_primes, args.reduction_prime_bound,
    )

    for record in records:
        samples = record.pop("jacobian_parameter_samples")
        manifest = "".join(
            f"{sample['multiple']}\t{sample['t']}\t{sample['s']}\n"
            for sample in samples
        )
        record["jacobian_parameter_sample_summary"] = {
            "count": len(samples),
            "multiple_range": (
                None if not samples else [samples[0]["multiple"], samples[-1]["multiple"]]
            ),
            "manifest_sha256": sha256(manifest.encode()).hexdigest(),
            "reconstruction": "replay the stored pointed Jacobian generator and declared multiple range",
        }
        del record["_data"]
    for row in random_extremes:
        row.pop("_words", None)
    for row in jacobian_extremes:
        row.pop("_words", None)
    result = {
        "schema": "elkies-k3.r17-genus-one-bisection-splitting-search.v1",
        "status": (
            "PASS_EXACT_SIMULTANEOUS_SPLIT_HITS_QUOTIENTED"
            if exact["simultaneous_split_hits"]
            else "PASS_BOUNDED_SEARCH_NO_EXACT_SIMULTANEOUS_SPLIT"
        ),
        "inputs": {relative(path): digest(path) for path in (
            MODEL, SECTIONS, PINNED, SHORT_COORDS, SHORT_GRAM, IDENTIFICATION,
            RANK_SCRIPT, CHORD_SCRIPT, PILOT_SCRIPT, DIVERSITY_SCRIPT, SCANNER,
            Path(__file__),
        )},
        "trace_selection": {
            "norm8": {
                "selected_count": len(norm8),
                "equation_pool_size": args.equation_pool_size,
                "selection": "greedy F2-rank/support-novelty/maximin-Hamming inside equation-cost prefix",
                "selected_span_dimension_over_F2": f2_rank(row["orbit_mask"] for row in norm8),
                "selected_equation_rank_range": [min(row["equation_rank"] for row in norm8), max(row["equation_rank"] for row in norm8)],
                "complete_frontier": norm8_enumeration,
            },
            "norm12_deep": {"selected_count": len(deep), "selection": "all deep classes", **deep_enumeration},
            "total_trace_and_cover_count": len(records),
        },
        "construction": {
            "norm8_member_method": "neutral generic-section seed, then rational parameter sampling on the pointed quartic Jacobian",
            "known_exceptional_points_used": False,
            "norm12_member_method": "unique regular M0 member; all 43 exact quartics compiled",
            "pointed_cover_count": sum(record["rational_seed"] is not None for record in records),
            "jacobian_generated_parameter_count": jacobian_population,
            "all_branch_quartics_irreducible_squarefree_and_smoothly_branched": True,
            "records": records,
        },
        "simultaneous_legendre_sieve": {
            "dynamic_mask_word_count": word_count,
            "primes": primes,
            "denominator_blind_curve_counts_by_prime": {
                str(prime): sum(scalar % prime == 0 for scalar, unused in curves)
                for prime in primes
            },
            "prime_blocks": [[prime for index, prime in enumerate(primes) if index % args.prime_block_count == block] for block in range(args.prime_block_count)],
            "ranking": "lexicographic weakest-block F2 trace-mask rank, cover count, then block sums",
            "random_population_requested": args.random_count,
            "random_coordinate_bound": args.coordinate_bound,
            "exhaustive_box_height": args.box_height,
            "scanner_summary": scanner_summary,
            "modular_extreme_count_random_or_box": len(random_extremes),
            "modular_extreme_count_jacobian_population": len(jacobian_extremes),
            "top_ranked_parameters": rankings,
            "table_sha256": digest(table_path),
        },
        "exact_stage": exact,
        "proof_boundary": (
            "The declared compact box is exhaustive and the pseudorandom/Jacobian populations are deterministic bounded samples. "
            "Legendre masks and their diversity ranks are only necessary local screens. Every reported hit has exact square, point, "
            "and finite-quotient records; no-hit output is not a global nonexistence theorem."
        ),
        "reproducing_command": "sage -python elkies-k3/scripts/search_elkies_2026_genus_one_bisection_splitting.sage",
        "implementation": {
            "compile_command": [
                "g++", "-O3", "-std=c++17", "-Wall", "-Wextra", "-pedantic",
                relative(SCANNER), "-o", "<local>/scan-genus-one-bisections",
            ],
            "scan_command": [
                "<local>/scan-genus-one-bisections", "<local>/legendre-mask-tables.bin",
                str(args.random_count), str(args.coordinate_bound), str(args.box_height),
                str(args.top_k), str(args.random_seed), "<local>/scan-output.txt",
            ],
            "scanner_source_sha256": digest(SCANNER),
        },
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != serialized:
            raise ArithmeticError("stored search artifact differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"R17G1SPLITV2|norm8={len(norm8)}|deep={len(deep)}|curves={len(records)}|"
        f"jacobian_population={jacobian_population}|modular_extremes={len(random_extremes)+len(jacobian_extremes)}|"
        f"hits={len(exact['simultaneous_split_hits'])}|status={result['status']}|output={args.output}"
    )


if __name__ == "__main__":
    main()
