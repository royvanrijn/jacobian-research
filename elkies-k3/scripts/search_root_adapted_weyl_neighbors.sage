#!/usr/bin/env sage
"""Exhaustively classify low-q neighbors of a root-adapted frame.

The input is a positive-definite frame Gram matrix of rank 17 whose first
``--root-rank`` basis vectors are simple roots spanning the full primitive
root lattice.  Instead of enumerating an enormous shell vector by vector,
the script enumerates it modulo the Weyl group: a representative is given by
its Mordell--Weil quotient coordinate and nonnegative Dynkin labels.  The
integrality congruence then recovers the unique root coordinates.

For every primitive horizontal orbit the checker constructs the exact
neighbor frame, classifies its roots, and root-adapts the child whenever its
root lattice is primitive.  Optional frame exports make exact breadth-first
continuation possible without reverting to a capped raw-shell search.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
U = matrix(ZZ, ((0, 1), (1, 0)))


def display_path(path):
    """Use a repository-relative path whenever the input lives in ROOT."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)

parser = argparse.ArgumentParser()
parser.add_argument("--frame", type=Path, required=True)
parser.add_argument("--root-rank", type=int, required=True)
parser.add_argument(
    "--q",
    type=int,
    action="append",
    required=True,
    help="neighbor norm q=ab; repeat for several q values",
)
parser.add_argument(
    "--degree",
    type=int,
    default=2,
    help="intersection with the old fiber (the b factor; default: 2)",
)
parser.add_argument("--output", type=Path)
parser.add_argument("--frames-dir", type=Path)
parser.add_argument(
    "--summary-only",
    action="store_true",
    help=(
        "retain exact per-q counts and root histograms but omit the potentially "
        "large neighbor-witness list from JSON output"
    ),
)
parser.add_argument(
    "--pari-gb",
    type=int,
    default=1,
    help="PARI stack size in GiB for large Mordell--Weil shells (default: 1)",
)
parser.add_argument(
    "--filter-marking",
    type=Path,
    help="optional exact source marking used to filter by one marked target degree",
)
parser.add_argument("--filter-target", help="target key in --filter-marking")
parser.add_argument(
    "--filter-max-degree",
    type=int,
    help="retain only fibres meeting --filter-target in at most this degree",
)
parser.add_argument(
    "--filter-linear-vector",
    help="optional comma-separated source coordinates for a second linear pairing filter",
)
parser.add_argument("--filter-linear-max", type=int)
parser.add_argument(
    "--adapt-mw-at-least",
    type=int,
    help=(
        "root-adapt and export only children of at least this MW rank; "
        "the default is the input MW rank"
    ),
)
parser.add_argument(
    "--rank-growth-only",
    action="store_true",
    help=(
        "screen non-growing children by exact root rank/count only; compute "
        "the full root lattice and determinant only for adapted MW hits"
    ),
)
parser.add_argument(
    "--stop-after-first-growth",
    action="store_true",
    help=(
        "stop each q pass after the first primitive root-adaptable child "
        "meeting --adapt-mw-at-least; this is a deterministic first-hit "
        "search, not an exhaustive classification"
    ),
)
parser.add_argument(
    "--stream-first-growth",
    action="store_true",
    help=(
        "with --stop-after-first-growth, test dominant witnesses in their "
        "deterministic generation order instead of materializing and sorting "
        "the full orbit set"
    ),
)
parser.add_argument(
    "--mw-vectors-cache",
    type=Path,
    help=(
        "optional JSON checkpoint for the Mordell--Weil quotient vectors; "
        "currently requires exactly one --q value"
    ),
)
parser.add_argument(
    "--mw-vector-cap",
    type=int,
    help=(
        "pass an exact storage cap m to PARI qfminim(B,m); the resulting "
        "quotient-vector sample is bounded and must not be called exhaustive"
    ),
)
parser.add_argument(
    "--include-zero-mw",
    action="store_true",
    help=(
        "also enumerate the zero Mordell--Weil projection; this is useful "
        "for auditing root-supported neighbours, which are omitted by default"
    ),
)
parser.add_argument(
    "--stream-skip",
    type=int,
    default=0,
    help="skip this many deterministic dominant witnesses before streaming tests",
)
parser.add_argument(
    "--stream-limit",
    type=int,
    help="test at most this many streamed witnesses after --stream-skip",
)
parser.add_argument(
    "--stream-progress-every",
    type=int,
    default=500,
    help="emit streamed-test progress every N witnesses (default: 500)",
)
parser.add_argument(
    "--select-orbit-index",
    type=int,
    action="append",
    help=(
        "process only these deterministic dominant-orbit indices; repeat for "
        "several indices. This is an exact selected-candidate reconstruction, "
        "not an exhaustive shell classification"
    ),
)
args = parser.parse_args()
if args.summary_only and args.frames_dir is not None:
    parser.error("--summary-only cannot be combined with --frames-dir")
if args.pari_gb < 1:
    parser.error("--pari-gb must be positive")
if args.pari_gb != 1:
    pari.allocatemem(args.pari_gb * 1024**3)
if args.stream_first_growth and not args.stop_after_first_growth:
    parser.error("--stream-first-growth requires --stop-after-first-growth")
if args.mw_vectors_cache is not None and len(set(args.q)) != 1:
    parser.error("--mw-vectors-cache requires exactly one distinct --q value")
if args.stream_skip < 0 or (args.stream_limit is not None and args.stream_limit <= 0):
    parser.error("stream skip must be nonnegative and stream limit must be positive")
if args.stream_progress_every <= 0:
    parser.error("--stream-progress-every must be positive")
if args.mw_vector_cap is not None and args.mw_vector_cap <= 0:
    parser.error("--mw-vector-cap must be positive")
if any(value is not None for value in
       (args.filter_marking, args.filter_target, args.filter_max_degree)):
    if any(value is None for value in
           (args.filter_marking, args.filter_target, args.filter_max_degree)):
        parser.error("the three --filter-* arguments must be supplied together")
    if args.filter_max_degree < 0:
        parser.error("--filter-max-degree must be nonnegative")
if (args.filter_linear_vector is None) != (args.filter_linear_max is None):
    parser.error("--filter-linear-vector and --filter-linear-max must be supplied together")


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def matrix_rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def bezout_vector_for_pairing(ns, fiber):
    pairings = list(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    if abs(current) != 1:
        return None
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_frame(ns, fiber, determinant):
    assert fiber * ns * fiber == 0
    mate = bezout_vector_for_pairing(ns, fiber)
    if mate is None:
        return None
    mate_square = ZZ(mate * ns * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    kernel = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    assert child.is_positive_definite() and child.det() == determinant
    neighbor_basis = matrix(
        ZZ,
        [list(fiber), list(mate)] + [list(row) for row in kernel.rows()],
    )
    assert abs(neighbor_basis.det()) == 1
    assert (
        neighbor_basis * ns * neighbor_basis.transpose()
        == block_diagonal_matrix(U, -child)
    )
    return child, neighbor_basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2]).columns()
    ]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(
        ZZ, [list(root) for root in roots]
    ).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (
        root_basis.rank(),
        count,
        abs(ZZ(root_gram.det())),
    )


def root_rank_and_count(gram):
    """Fast exact screen that avoids an integral root-lattice HNF."""
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return 0, 0
    half = matrix(ZZ, result[2])
    return half.rank(), count


def deterministic_simple_roots(gram):
    roots, _, data = roots_and_data(gram)
    root_rank = data[0]
    # Lexicographic sign is an additive total ordering, hence defines a
    # positive root system without relying on a bounded search for a regular
    # linear functional in a possibly very large coordinate basis.
    positive = [
        root
        for root in roots
        if next(value for value in root if value != 0) > 0
    ]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root
        for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == root_rank
    cartan = simple * gram * simple.transpose()
    assert set(cartan.diagonal()) == {2}
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(root_rank)
        for column in range(root_rank)
        if row != column
    )
    return simple, cartan


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other]]
            for other in adjacent:
                unseen.remove(other)
                todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(sorted(result, key=lambda component: (len(component), component)))


def component_name(cartan, component):
    block = cartan.matrix_from_rows_and_columns(component, component)
    rank = block.nrows()
    determinant = abs(ZZ(block.det()))
    root_count = ZZ(pari(block).qfminim(2)[0])
    if determinant == rank + 1 and root_count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and determinant == 4 and root_count == 2 * rank * (rank - 1):
        return f"D{rank}"
    exceptional = {
        (6, 3, 72): "E6",
        (7, 2, 126): "E7",
        (8, 1, 240): "E8",
    }
    return exceptional.get((rank, determinant, root_count), f"R{rank}d{determinant}n{root_count}")


def ade_name(cartan):
    return "+".join(
        component_name(cartan, component)
        for component in connected_components(cartan)
    )


def dominant_weights(cartan, component, bound):
    block = cartan.matrix_from_rows_and_columns(component, component)
    inverse = block.inverse()
    assert all(value >= 0 for value in inverse.list())
    weights = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == len(component):
            weights.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index] * value**2
            added += 2 * value * sum(
                inverse[index, previous] * prefix[previous]
                for previous in range(index)
            )
            new_norm = norm + added
            if new_norm > bound:
                break
            recurse(prefix + [value], new_norm)
            value += 1

    recurse([], QQ(0))
    return tuple(weights)


def root_adaptation(child):
    _, root_basis, invariants = roots_and_data(child)
    root_rank = invariants[0]
    if root_rank == 0:
        lll = matrix(ZZ, pari(child).qflllgram())
        adapted_basis = lll.transpose()
        adapted = adapted_basis * child * adapted_basis.transpose()
        assert abs(adapted_basis.det()) == 1
        return adapted, adapted_basis, adapted.change_ring(QQ)
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    smith_diagonal = tuple(
        abs(smith[index, index]) for index in range(root_rank)
    )
    if smith_diagonal != (1,) * root_rank:
        return None
    simple, cartan = deterministic_simple_roots(child)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    assert abs(adapted_basis.det()) == 1
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank), lll.transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    assert abs(adapted_basis.det()) == 1
    assert adapted[:root_rank, :root_rank] == cartan
    return adapted, adapted_basis, height


def json_default(value):
    if isinstance(value, Integer):
        return int(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


frame_path = args.frame.resolve()
frame_display_path = display_path(frame_path)
frame = load_gram(frame_path)
assert frame.nrows() == frame.ncols() == 17
assert frame.is_positive_definite()
determinant = abs(ZZ(frame.det()))
root_rank = args.root_rank
assert 0 <= root_rank < 17
cartan = frame[:root_rank, :root_rank]
if root_rank:
    assert set(cartan.diagonal()) == {2}
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(root_rank)
        for column in range(root_rank)
        if row != column
    )
    assert cartan.rank() == root_rank
input_roots, input_root_basis, input_root_data = roots_and_data(frame)
assert input_root_data[0] == root_rank
if root_rank:
    assert matrix(ZZ, [list(row) for row in input_roots]).row_module() == input_root_basis.row_module()
components = connected_components(cartan)
input_ade = ade_name(cartan)
adapt_mw_at_least = (
    17 - root_rank
    if args.adapt_mw_at_least is None
    else args.adapt_mw_at_least
)
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail if root_rank == 0 else tail - coupling.transpose() * cartan.inverse() * coupling
assert height.is_positive_definite()

all_records = []
summaries = []
ns = block_diagonal_matrix(U, -frame)
marked_filter = None
linear_filter = None
if args.filter_marking is not None:
    filter_payload = json.loads(args.filter_marking.resolve().read_text())
    target_container = filter_payload.get(
        "target_fibres_in_root_adapted_hub",
        filter_payload.get("target_fibres_in_child"),
    )
    assert target_container is not None and args.filter_target in target_container
    marked_filter = vector(ZZ, target_container[args.filter_target])
    assert marked_filter * ns * marked_filter == 0
if args.filter_linear_vector is not None:
    linear_filter = vector(ZZ, [ZZ(value) for value in args.filter_linear_vector.split(",")])
    assert len(linear_filter) == 19
for q in sorted(set(args.q)):
    assert q > 0 and q % args.degree == 0
    target = ZZ(2 * q)
    factor_a = ZZ(q // args.degree)
    factor_b = ZZ(args.degree)

    height_scale = lcm(entry.denominator() for entry in height.list())
    scaled_height = (height_scale * height).change_ring(ZZ)
    cache_key = {
        "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
        "q": int(q),
        "target": int(target),
        "height": rational_rows(height),
    }
    if args.mw_vector_cap is not None:
        cache_key["mw_vector_cap"] = int(args.mw_vector_cap)
    if args.include_zero_mw:
        cache_key["include_zero_mw"] = True
    if args.mw_vectors_cache is not None and args.mw_vectors_cache.exists():
        cache_payload = json.loads(args.mw_vectors_cache.read_text())
        assert cache_payload["schema"] == "root-adapted-mw-vectors.v1"
        assert cache_payload["key"] == cache_key
        mw_vectors = tuple(
            vector(ZZ, values) for values in cache_payload["vectors"]
        )
        mw_pari_count = int(
            cache_payload.get("pari_vector_count", 2 * len(mw_vectors))
        )
        assert all(
            (args.include_zero_mw or value != 0)
            and value * height * value <= target
            for value in mw_vectors
        )
        print(
            "ROOTWEYL_MW_CACHE|q={}|vectors={}|path={}|status=LOADED".format(
                q, len(mw_vectors), display_path(args.mw_vectors_cache)
            ),
            flush=True,
        )
    else:
        mw_result = (
            pari(scaled_height).qfminim(height_scale * target)
            if args.mw_vector_cap is None
            else pari(scaled_height).qfminim(
                height_scale * target, args.mw_vector_cap
            )
        )
        mw_pari_count = int(mw_result[0])
        mw_vector_map = {}
        for column in matrix(ZZ, mw_result[2]).columns():
            for sign in (1, -1):
                value = sign * vector(ZZ, column)
                if value == 0 or value * height * value > target:
                    continue
                canonical = min(tuple(value), tuple(-value))
                mw_vector_map[canonical] = vector(ZZ, canonical)
        mw_vectors = tuple(
            sorted(
                mw_vector_map.values(),
                key=lambda value: (value * height * value, tuple(value)),
            )
        )
        if args.include_zero_mw:
            mw_vectors = (vector(ZZ, [0] * height.nrows()),) + mw_vectors
        if args.mw_vectors_cache is not None:
            args.mw_vectors_cache.parent.mkdir(parents=True, exist_ok=True)
            args.mw_vectors_cache.write_text(
                json.dumps(
                    {
                        "schema": "root-adapted-mw-vectors.v1",
                        "key": cache_key,
                        "pari_vector_count": mw_pari_count,
                        "vectors": [list(map(int, value)) for value in mw_vectors],
                    },
                    separators=(",", ":"),
                )
                + "\n"
            )
            print(
                "ROOTWEYL_MW_CACHE|q={}|vectors={}|path={}|status=SAVED".format(
                    q, len(mw_vectors), display_path(args.mw_vectors_cache)
                ),
                flush=True,
            )

    component_weights = tuple(
        dominant_weights(cartan, component, QQ(target))
        for component in components
    )
    # Convolve the component weight lists once and index them by exact norm.
    # Repeating their full Cartesian product for every MW vector is the main
    # avoidable cost for split types such as A9+D4.
    combined_weights_by_norm = {}

    def combine_component_weights(index, choices, norm):
        if index == len(component_weights):
            combined_weights_by_norm.setdefault(norm, []).append(tuple(choices))
            return
        for values, weight_norm in component_weights[index]:
            new_norm = norm + weight_norm
            if new_norm <= target:
                combine_component_weights(
                    index + 1,
                    choices + [(values, weight_norm)],
                    new_norm,
                )

    combine_component_weights(0, [], QQ(0))
    cartan_inverse = cartan.inverse() if root_rank else matrix(QQ, 0, 0)
    dominant_orbits = {}
    stream_seen = set()
    stream_hit = None
    stream_nonprimitive = 0
    stream_tested = 0
    stream_limit_reached = False
    for mw in mw_vectors:
        mw_norm = mw * height * mw
        for choices in combined_weights_by_norm.get(target - mw_norm, ()):
            labels = vector(ZZ, [0] * root_rank)
            for component, (values, _) in zip(components, choices):
                for index, value in zip(component, values):
                    labels[index] = value
            root_coordinates = (
                vector(ZZ, [])
                if root_rank == 0
                else cartan_inverse * (labels - coupling * mw)
            )
            if not all(value in ZZ for value in root_coordinates):
                continue
            witness = vector(ZZ, list(root_coordinates) + list(mw))
            assert witness * frame * witness == target
            if root_rank:
                assert cartan * root_coordinates + coupling * mw == labels
            witness_key = tuple(witness)
            orbit_data = (tuple(mw), tuple(labels))
            if not args.stream_first_growth:
                dominant_orbits[witness_key] = orbit_data
                continue
            if witness_key in stream_seen:
                continue
            stream_seen.add(witness_key)
            stream_index = len(stream_seen)
            if stream_index <= args.stream_skip:
                continue
            if args.stream_limit is not None and stream_tested >= args.stream_limit:
                stream_limit_reached = True
                break
            stream_tested += 1
            if stream_tested % args.stream_progress_every == 0:
                print(
                    "ROOTWEYL_STREAM_PROGRESS|q={}|index={}|tested={}|status=RUNNING".format(
                        q, stream_index, stream_tested
                    ),
                    flush=True,
                )
            fiber = vector(ZZ, [factor_a, factor_b] + list(witness))
            result = child_frame(ns, fiber, determinant)
            if result is None:
                stream_nonprimitive += 1
                continue
            child, _ = result
            child_root_rank, _ = root_rank_and_count(child)
            if 17 - int(child_root_rank) < adapt_mw_at_least:
                continue
            if root_adaptation(child) is None:
                continue
            stream_hit = (witness_key, orbit_data, len(stream_seen))
            break
        if stream_hit is not None or stream_limit_reached:
            break

    orbit_count = (
        args.stream_skip + stream_tested
        if args.stream_first_growth
        else len(dominant_orbits)
    )

    print(
        "ROOTWEYL_ENUM|q={}|ab={},{}|mw_vectors={}|component_weight_counts={}|"
        "combined_weights={}|orbits={}|status=PASS".format(
            q,
            factor_a,
            factor_b,
            len(mw_vectors),
            tuple(map(len, component_weights)),
            sum(map(len, combined_weights_by_norm.values())),
            orbit_count,
        ),
        flush=True,
    )

    histogram = Counter()
    nonprimitive = 0
    q_records = []
    primitive_neighbors = 0
    maximum_child_mw_rank = 17 - root_rank
    screened_orbits = orbit_count if args.stream_first_growth else 0
    stopped_early = False
    if args.stream_first_growth:
        orbit_entries = (
            []
            if stream_hit is None
            else [(stream_hit[2], (stream_hit[0], stream_hit[1]))]
        )
    else:
        orbit_entries = enumerate(sorted(dominant_orbits.items()), start=1)
    for orbit_index, (witness_tuple, orbit_data) in orbit_entries:
        screened_orbits = orbit_index
        if (args.select_orbit_index is not None
                and orbit_index not in set(args.select_orbit_index)):
            continue
        mw_tuple, labels_tuple = orbit_data
        witness = vector(ZZ, witness_tuple)
        fiber = vector(ZZ, [factor_a, factor_b] + list(witness))
        marked_filter_degree = None
        if marked_filter is not None:
            marked_filter_degree = ZZ(fiber * ns * marked_filter)
            assert marked_filter_degree >= 0
            if marked_filter_degree > args.filter_max_degree:
                continue
        linear_filter_value = None
        if linear_filter is not None:
            linear_filter_value = ZZ(fiber * ns * linear_filter)
            if linear_filter_value < 0 or linear_filter_value > args.filter_linear_max:
                continue
        result = child_frame(ns, fiber, determinant)
        if result is None:
            nonprimitive += 1
            continue
        child, neighbor_basis = result
        if args.rank_growth_only:
            child_root_rank, child_root_count = root_rank_and_count(child)
            child_mw_rank = 17 - int(child_root_rank)
            if child_mw_rank >= adapt_mw_at_least:
                _, _, child_root_data = roots_and_data(child)
                assert child_root_data[:2] == (
                    child_root_rank,
                    child_root_count,
                )
            else:
                child_root_data = (child_root_rank, child_root_count, None)
        else:
            _, _, child_root_data = roots_and_data(child)
            child_mw_rank = 17 - int(child_root_data[0])
        histogram[child_root_data] += 1
        adapted_data = (
            root_adaptation(child)
            if child_mw_rank >= adapt_mw_at_least
            else None
        )
        if adapted_data is None:
            child_adapted = adapted_basis = child_height = None
            child_ade = "not-adapted"
        else:
            child_adapted, adapted_basis, child_height = adapted_data
            child_cartan = child_adapted[: child_root_data[0], : child_root_data[0]]
            child_ade = ade_name(child_cartan) if child_root_data[0] else "rootless"
        record = {
            "q": int(q),
            "factor_order": [int(factor_a), int(factor_b)],
            "old_fiber_degree": int(factor_b),
            "orbit_index": int(orbit_index),
            "mw_projection": list(map(int, mw_tuple)),
            "dominant_labels": list(map(int, labels_tuple)),
            "witness": list(map(int, witness)),
            "fiber": list(map(int, fiber)),
            "child_root_data": [
                None if value is None else int(value)
                for value in child_root_data
            ],
            "child_mw_rank": child_mw_rank,
            "child_ade": child_ade,
        }
        if marked_filter_degree is not None:
            record["marked_filter"] = {
                "target": args.filter_target,
                "degree": int(marked_filter_degree),
                "maximum": args.filter_max_degree,
            }
        if linear_filter_value is not None:
            record["marked_linear_filter"] = {
                "value": int(linear_filter_value),
                "maximum": args.filter_linear_max,
                "vector": list(map(int, linear_filter)),
            }
        if adapted_data is not None:
            record.update(
                {
                    "child_frame": matrix_rows(child),
                    "neighbor_basis": matrix_rows(neighbor_basis),
                    "child_root_adapted_frame": matrix_rows(child_adapted),
                    "child_root_adapted_basis": matrix_rows(adapted_basis),
                    "child_mw_height": rational_rows(child_height),
                }
            )
        primitive_neighbors += 1
        maximum_child_mw_rank = max(maximum_child_mw_rank, child_mw_rank)
        if not args.summary_only:
            q_records.append(record)
            all_records.append(record)
        if args.stop_after_first_growth and adapted_data is not None:
            stopped_early = True
            break

    summary = {
        "q": int(q),
        "factor_order": [int(factor_a), int(factor_b)],
        "mw_projection_representatives": len(mw_vectors),
        "mw_pari_vector_count": mw_pari_count,
        "mw_vector_cap": args.mw_vector_cap,
        "mw_enumeration_complete": args.mw_vector_cap is None,
        "dominant_orbits": orbit_count,
        "dominant_orbits_complete": not args.stream_first_growth,
        "stream_precheck_nonprimitive": int(stream_nonprimitive),
        "stream_skip": int(args.stream_skip),
        "stream_tested": int(stream_tested),
        "stream_limit_reached": stream_limit_reached,
        "screened_orbits": int(screened_orbits),
        "search_stopped_early": stopped_early,
        "primitive_neighbors": primitive_neighbors,
        "nonprimitive_orbits": int(nonprimitive),
        "root_histogram": [
            {
                "root_rank": int(invariants[0]),
                "root_count": int(invariants[1]),
                "root_determinant": (
                    None if invariants[2] is None else int(invariants[2])
                ),
                "mw_rank": 17 - int(invariants[0]),
                "orbit_count": int(count),
            }
            for invariants, count in sorted(histogram.items())
        ],
    }
    summaries.append(summary)
    print(
        "ROOTWEYL|q={}|ab={},{}|mw_vectors={}|orbits={}|primitive={}|"
        "nonprimitive={}|screened={}|max_mw={}|status={}".format(
            q,
            factor_a,
            factor_b,
            len(mw_vectors),
            orbit_count,
            primitive_neighbors,
            nonprimitive,
            screened_orbits,
            maximum_child_mw_rank,
            (
                "PASS_FIRST_HIT"
                if stopped_early
                else "PASS_STREAM_CHUNK_NO_HIT"
                if stream_limit_reached
                else "PASS"
            ),
        ),
        flush=True,
    )
    for item in summary["root_histogram"]:
        print(
            "ROOTWEYL_HIST|q={}|root_data={},{},{}|MW={}|orbits={}".format(
                q,
                item["root_rank"],
                item["root_count"],
                "NA" if item["root_determinant"] is None else item["root_determinant"],
                item["mw_rank"],
                item["orbit_count"],
            ),
            flush=True,
        )

if args.frames_dir is not None:
    args.frames_dir.mkdir(parents=True, exist_ok=True)
    for record in all_records:
        if "child_root_adapted_frame" not in record:
            continue
        digest = hashlib.sha256(
            repr(record["child_root_adapted_frame"]).encode()
        ).hexdigest()[:12]
        data = record["child_root_data"]
        path = args.frames_dir / (
            f"q{record['q']}-o{record['orbit_index']:04d}-"
            f"r{data[0]}-n{data[1]}-d{data[2]}-{digest}.txt"
        )
        lines = [
            f"# source = {frame_display_path}",
            f"# q = {record['q']}",
            f"# factor_order = {tuple(record['factor_order'])}",
            f"# orbit_index = {record['orbit_index']}",
            f"# witness = {tuple(record['witness'])}",
            f"# root_data = {tuple(data)}",
            f"# ADE = {record['child_ade']}",
            f"# MW_rank = {record['child_mw_rank']}",
        ]
        lines.extend(
            " ".join(map(str, row))
            for row in record["child_root_adapted_frame"]
        )
        path.write_text("\n".join(lines) + "\n")

payload = {
    "status": (
        "PASS_ROOT_ADAPTED_WEYL_SELECTED_NEIGHBORS"
        if args.select_orbit_index is not None
        else "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS_TARGET_FILTERED"
        if marked_filter is not None else "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    ),
    "frame": frame_display_path,
    "determinant": int(determinant),
    "input_root_data": list(map(int, input_root_data)),
    "input_ade": input_ade,
    "input_mw_rank": 17 - int(input_root_data[0]),
    "input_mw_height": rational_rows(height),
    "summaries": summaries,
    "neighbors": all_records,
}
if args.summary_only:
    payload["summary_only"] = True
    payload["witness_retention"] = (
        "Exact enumeration and classification were performed, but individual "
        "neighbor witnesses were deliberately omitted from this artifact."
    )
if args.select_orbit_index is not None:
    payload["selected_orbit_indices"] = sorted(set(args.select_orbit_index))
    payload["selection_scope"] = (
        "Exact reconstruction and root adaptation of the listed deterministic "
        "dominant-orbit indices; no completeness claim for the surrounding shell."
    )
if marked_filter is not None:
    payload["marked_target_filter"] = {
        "marking": display_path(args.filter_marking),
        "target": args.filter_target,
        "maximum_degree": args.filter_max_degree,
        "scope": "complete dominant-orbit enumeration at the listed q values",
    }
if linear_filter is not None:
    payload["marked_linear_filter"] = {
        "vector": list(map(int, linear_filter)),
        "maximum": args.filter_linear_max,
        "scope": "complete dominant-orbit enumeration at the listed q values",
    }
if args.output is not None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n"
    )
