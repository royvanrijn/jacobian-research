from pathlib import Path
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--h", type=float, default=22.94)
parser.add_argument("--norm-tol", type=float, default=0.015)
parser.add_argument("--pair-tol", type=float, default=0.06)

parser.add_argument("--anchors", type=int, default=4000)
parser.add_argument("--beam", type=int, default=64)
parser.add_argument("--branch", type=int, default=6)
parser.add_argument("--seed", type=int, default=210017)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "results"

Q = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

H = np.loadtxt(args.gram, dtype=float)

Vall = np.load(
    args.vectors,
    mmap_mode="r"
)

Nall = np.load(
    args.norms,
    mmap_mode="r"
)

assert Q.shape == (17, 17)
assert np.all(np.diag(Q) == 4)
assert H.shape == (Vall.shape[1], Vall.shape[1])

rng = np.random.default_rng(args.seed)

h = args.h
target_norm = 4.0 * h

# ============================================================
# Narrow fixed shell.
# ============================================================

mask = (
    np.abs(Nall / target_norm - 1.0)
    <= args.norm_tol
)

global_ids = np.flatnonzero(mask)

print("TARGET")
print("h =", h)
print("target_norm =", target_norm)
print("norm_tol =", args.norm_tol)
print("pair_tol =", args.pair_tol)
print()
print("ambient dimension =", Vall.shape[1])
print("full pool =", len(Vall))
print("shell =", len(global_ids))

if len(global_ids) < 17:
    raise SystemExit("shell too small")

# Copy shell into RAM. For rank21 this is manageable.
W = np.asarray(
    Vall[global_ids],
    dtype=np.int64
)

norms = np.asarray(
    Nall[global_ids],
    dtype=float
)

m = len(W)

# WH lets us compute pairings with any candidate vector:
#
#       <W_i,v> = (W H) v
#
WH = W @ H

print("WH ready", flush=True)

# ============================================================
# Target ordering.
#
# Start from restrictive row 3, matching the successful probes.
# Then greedily select the row having strongest constraints
# against already selected rows.
# ============================================================

anchor = 3

order = [anchor]
remaining = set(range(17))
remaining.remove(anchor)

while remaining:
    ti = max(
        remaining,
        key=lambda x: (
            sum(Q[x, j] != 0 for j in order),
            sum(abs(int(Q[x, j])) for j in order),
            np.count_nonzero(Q[x])
        )
    )

    order.append(ti)
    remaining.remove(ti)

print("target order =", order)

# ============================================================
# Candidate quality.
# ============================================================

norm_error = np.abs(
    norms / target_norm - 1.0
)

# Prefer vectors closest to exact target norm for anchors.
anchor_order = np.argsort(norm_error)

if len(anchor_order) > args.anchors:
    # Half deterministic closest vectors, half spread randomly
    # throughout the valid shell.
    nclose = args.anchors // 2

    close = anchor_order[:nclose]

    rest = anchor_order[nclose:]

    random_part = rng.choice(
        rest,
        size=args.anchors - nclose,
        replace=False
    )

    anchor_candidates = np.concatenate(
        [close, random_part]
    )
else:
    anchor_candidates = anchor_order

print("anchors =", len(anchor_candidates))

# ============================================================
# A state is:
#
#   score
#   assignment {target_index: shell_index}
#   signs      {target_index: +/-1}
#
# We explicitly permit sign orientation because our pool stores
# unoriented vectors.
# ============================================================

best_depth = 0
best_state = None
started = time.time()

# Pairing cache:
#
# shell candidate index -> pairings against every vector in shell
#
# Each cached item is ~m*8 bytes. Keep bounded.
pair_cache = {}
pair_cache_order = []
CACHE_MAX = 512


def pairing_column(ci):
    if ci in pair_cache:
        return pair_cache[ci]

    p = WH @ W[ci]

    pair_cache[ci] = p
    pair_cache_order.append(ci)

    if len(pair_cache_order) > CACHE_MAX:
        old = pair_cache_order.pop(0)
        pair_cache.pop(old, None)

    return p


def extend_state(state, ti):
    score, assignment, signs = state

    used = set(assignment.values())

    possible = np.ones(
        m,
        dtype=bool
    )

    total_error = np.zeros(
        m,
        dtype=float
    )

    # Norm condition already enforced by shell, but include
    # continuous norm error in ranking.
    total_error += norm_error

    for tj, cj in assignment.items():

        sj = signs[tj]

        p = pairing_column(cj) * sj

        wanted = h * Q[ti, tj]

        # candidate itself may be used as +v or -v.
        ep = np.abs(p - wanted) / h
        em = np.abs(-p - wanted) / h

        e = np.minimum(ep, em)

        possible &= (
            e <= args.pair_tol
        )

        total_error += e

    if used:
        possible[
            np.fromiter(
                used,
                dtype=np.int64
            )
        ] = False

    ids = np.flatnonzero(possible)

    if not len(ids):
        return []

    # Take strongest local matches only.
    take = min(
        args.branch,
        len(ids)
    )

    if len(ids) > take:
        local = np.argpartition(
            total_error[ids],
            take - 1
        )[:take]

        ids = ids[local]

    ids = ids[
        np.argsort(total_error[ids])
    ]

    out = []

    for ci in ids:

        # Determine best orientation against all previous vectors.
        plus_error = 0.0
        minus_error = 0.0

        for tj, cj in assignment.items():

            sj = signs[tj]

            actual = (
                pairing_column(cj)[ci]
                * sj
            )

            wanted = (
                h * Q[ti, tj]
            )

            plus_error += abs(
                actual - wanted
            ) / h

            minus_error += abs(
                -actual - wanted
            ) / h

        sign = (
            1
            if plus_error <= minus_error
            else -1
        )

        na = dict(assignment)
        ns = dict(signs)

        na[ti] = int(ci)
        ns[ti] = sign

        new_score = (
            score
            + float(total_error[ci])
        )

        out.append(
            (
                new_score,
                na,
                ns
            )
        )

    return out


# ============================================================
# Search anchors.
# ============================================================

states = []

for ci in anchor_candidates:

    states.append(
        (
            float(norm_error[ci]),
            {anchor: int(ci)},
            {anchor: 1}
        )
    )

# Global sign symmetry lets anchor always be +.
states.sort(key=lambda x:x[0])
states = states[:args.beam]

best_depth = 1
best_state = states[0]

print(
    f"BEST|depth=1"
    f"|score={best_state[0]:.8g}",
    flush=True
)

for depth in range(1, 17):

    ti = order[depth]

    print()
    print(
        f"LEVEL|depth={depth+1}"
        f"|target={ti}"
        f"|states={len(states)}"
        f"|cache={len(pair_cache)}",
        flush=True
    )

    candidates = []

    for state in states:
        candidates.extend(
            extend_state(
                state,
                ti
            )
        )

    if not candidates:
        print("DEAD END")
        break

    candidates.sort(
        key=lambda x:x[0]
    )

    # Deduplicate states by chosen ambient vectors.
    seen = set()
    states = []

    for state in candidates:

        key = tuple(
            sorted(
                state[1].values()
            )
        )

        if key in seen:
            continue

        seen.add(key)
        states.append(state)

        if len(states) >= args.beam:
            break

    if not states:
        break

    best_depth = depth + 1
    best_state = states[0]

    print(
        f"BEST|depth={best_depth}"
        f"|score={best_state[0]:.8g}"
        f"|states={len(states)}"
        f"|seconds={time.time()-started:.3f}",
        flush=True
    )

    print(
        " assignment=",
        best_state[1],
        flush=True
    )

    if best_depth == 17:
        break


# ============================================================
# Verify best result.
# ============================================================

print()
print("FINAL")
print("best_depth =", best_depth)
print("h =", h)
print("seconds =", time.time() - started)

if best_state is not None:

    score, assignment, signs = best_state

    print("score =", score)
    print("assignment =", assignment)
    print("signs =", signs)

    selected_targets = sorted(
        assignment
    )

    A = np.array(
        [
            signs[ti]
            * W[assignment[ti]]
            for ti in selected_targets
        ],
        dtype=np.int64
    )

    G = (
        A @ H
    ) @ A.T

    QT = Q[
        np.ix_(
            selected_targets,
            selected_targets
        )
    ]

    residual = (
        G - h * QT
    )

    rel = (
        np.linalg.norm(residual)
        /
        (
            h
            * np.linalg.norm(QT)
        )
    )

    print("relative residual =", rel)

    print()
    print("selected global pool IDs:")

    for ti in selected_targets:

        si = assignment[ti]

        print(
            f"target={ti}"
            f" shell={si}"
            f" global={int(global_ids[si])}"
            f" sign={signs[ti]}"
            f" norm={norms[si]:.12g}"
        )

    np.savetxt(
        OUT / "rank17-rank21-targeted-partial-A.txt",
        A,
        fmt="%d"
    )

    np.savetxt(
        OUT / "rank17-rank21-targeted-partial-residual.txt",
        residual,
        fmt="%.17g"
    )

