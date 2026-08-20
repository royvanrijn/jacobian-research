from pathlib import Path
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--h-min", type=float, default=22.84)
parser.add_argument("--h-max", type=float, default=23.04)
parser.add_argument("--h-step", type=float, default=0.02)

parser.add_argument("--norm-tol", type=float, default=0.015)
parser.add_argument("--pair-tol", type=float, default=0.06)

parser.add_argument("--anchors", type=int, default=8000)
parser.add_argument("--beam", type=int, default=256)
parser.add_argument("--branch", type=int, default=10)

parser.add_argument("--anchor-target", type=int, default=3)
parser.add_argument("--seed", type=int, default=210017)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]

J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "results"

Q = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

H = np.loadtxt(
    args.gram,
    dtype=float
)

VALL = np.load(
    args.vectors,
    mmap_mode="r"
)

NALL = np.load(
    args.norms,
    mmap_mode="r"
)

assert Q.shape == (17, 17)
assert np.all(np.diag(Q) == 4)

assert VALL.ndim == 2
assert H.shape == (VALL.shape[1], VALL.shape[1])
assert len(VALL) == len(NALL)

rng = np.random.default_rng(args.seed)

# Generator writes these sorted. Check rather than assume.
norms_sorted = bool(
    np.all(
        NALL[1:] >= NALL[:-1]
    )
)

print("ambient_dim =", VALL.shape[1])
print("full_pool =", len(VALL))
print("norms_sorted =", norms_sorted)
print()

summary = []

GLOBAL_BEST = None


def shell_indices(target_norm):

    lo = target_norm * (1.0 - args.norm_tol)
    hi = target_norm * (1.0 + args.norm_tol)

    if norms_sorted:

        a = int(
            np.searchsorted(
                NALL,
                lo,
                side="left"
            )
        )

        b = int(
            np.searchsorted(
                NALL,
                hi,
                side="right"
            )
        )

        return np.arange(
            a,
            b,
            dtype=np.int64
        )

    return np.flatnonzero(
        (NALL >= lo)
        &
        (NALL <= hi)
    )


def run_h(h):

    global GLOBAL_BEST

    started = time.time()

    target_norm = 4.0 * h
    pair_abs_tol = args.pair_tol * h

    gids = shell_indices(
        target_norm
    )

    print("=" * 76)
    print(
        f"H_START|h={h:.8f}"
        f"|target_norm={target_norm:.8f}"
        f"|shell={len(gids)}"
    )

    if len(gids) < 17:
        print("SHELL TOO SMALL")
        return None

    # --------------------------------------------------------
    # Load this shell into RAM.
    # --------------------------------------------------------

    W = np.asarray(
        VALL[gids],
        dtype=np.int64
    )

    norms = np.asarray(
        NALL[gids],
        dtype=float
    )

    m = len(W)

    norm_error = np.abs(
        norms / target_norm - 1.0
    )

    # Pairing with arbitrary v:
    #
    #     <W_i,v> = (W H) v
    #
    WH = W @ H

    print(
        f"SHELL_READY|m={m}",
        flush=True
    )

    # --------------------------------------------------------
    # Anchor selection.
    #
    # We know target 3 behaved well previously, and global
    # overall sign symmetry lets anchor sign be +1.
    # --------------------------------------------------------

    by_norm = np.argsort(
        norm_error
    )

    if len(by_norm) <= args.anchors:

        anchor_candidates = by_norm

    else:

        close_count = args.anchors // 2

        close = by_norm[:close_count]

        remainder = by_norm[close_count:]

        random_part = rng.choice(
            remainder,
            size=args.anchors - close_count,
            replace=False
        )

        anchor_candidates = np.concatenate(
            [
                close,
                random_part
            ]
        )

    anchor_target = args.anchor_target

    # --------------------------------------------------------
    # State:
    #
    # (
    #   score,
    #   assignment : {target -> shell candidate},
    #   signs      : {target -> +/-1},
    #   history    : [target order actually chosen]
    # )
    # --------------------------------------------------------

    states = []

    for ci in anchor_candidates:

        states.append(
            (
                float(norm_error[ci]),
                {
                    anchor_target:
                    int(ci)
                },
                {
                    anchor_target:
                    1
                },
                [
                    anchor_target
                ]
            )
        )

    states.sort(
        key=lambda x: x[0]
    )

    states = states[
        :args.beam
    ]

    best_depth = 1
    best_state = states[0]

    print(
        f"BEST|h={h:.8f}"
        f"|depth=1"
        f"|score={best_state[0]:.10g}",
        flush=True
    )

    # --------------------------------------------------------
    # Small LRU-like pairing cache.
    #
    # Do NOT use the old 512-column cache with a ~650k shell;
    # that can consume several GB.
    # --------------------------------------------------------

    pair_cache = {}
    pair_age = []
    CACHE_MAX = 64

    def pairing_column(ci):

        if ci in pair_cache:
            return pair_cache[ci]

        p = WH @ W[ci]

        pair_cache[ci] = p
        pair_age.append(ci)

        if len(pair_age) > CACHE_MAX:

            old = pair_age.pop(0)

            pair_cache.pop(
                old,
                None
            )

        return p

    # --------------------------------------------------------
    # Evaluate one potential next target.
    #
    # Critical correctness point:
    #
    # Candidate +v must satisfy ALL existing constraints,
    # OR candidate -v must satisfy ALL existing constraints.
    #
    # We may not choose sign independently for each pairing.
    # --------------------------------------------------------

    def domain_for_target(
        ti,
        assignment,
        signs
    ):

        plus_ok = np.ones(
            m,
            dtype=bool
        )

        minus_ok = np.ones(
            m,
            dtype=bool
        )

        plus_error = np.zeros(
            m,
            dtype=np.float64
        )

        minus_error = np.zeros(
            m,
            dtype=np.float64
        )

        for tj, cj in assignment.items():

            sj = signs[tj]

            # Pairing of +candidate with the already oriented
            # selected vector.
            actual = (
                pairing_column(cj)
                * sj
            )

            wanted = (
                h * Q[ti, tj]
            )

            ep = np.abs(
                actual - wanted
            )

            em = np.abs(
                -actual - wanted
            )

            plus_ok &= (
                ep <= pair_abs_tol
            )

            minus_ok &= (
                em <= pair_abs_tol
            )

            plus_error += (
                ep / h
            )

            minus_error += (
                em / h
            )

        # Same unoriented lattice line cannot be reused.
        for ci in assignment.values():

            plus_ok[ci] = False
            minus_ok[ci] = False

        line_ok = (
            plus_ok
            |
            minus_ok
        )

        count = int(
            np.count_nonzero(
                line_ok
            )
        )

        return (
            count,
            plus_ok,
            minus_ok,
            plus_error,
            minus_error
        )

    # --------------------------------------------------------
    # MRV search.
    #
    # For every state:
    #   - evaluate EVERY remaining target;
    #   - if any target has zero domain, the state is impossible;
    #   - choose the target with the smallest nonzero domain;
    #   - branch on its best oriented candidates.
    # --------------------------------------------------------

    for depth in range(
        1,
        17
    ):

        children = []

        zero_domain_states = 0

        mrv_hist = {}

        for (
            parent_score,
            assignment,
            signs,
            history
        ) in states:

            remaining = [
                ti
                for ti in range(17)
                if ti not in assignment
            ]

            if not remaining:
                continue

            best_ti = None
            best_domain = None
            impossible = False

            # Find MRV target.
            for ti in remaining:

                d = domain_for_target(
                    ti,
                    assignment,
                    signs
                )

                count = d[0]

                if count == 0:

                    # Every target must eventually be assigned.
                    # If even one currently has no candidate,
                    # this entire state cannot become full.
                    impossible = True
                    break

                if (
                    best_domain is None
                    or count < best_domain[0]
                ):

                    best_ti = ti
                    best_domain = d

            if impossible:

                zero_domain_states += 1
                continue

            (
                count,
                plus_ok,
                minus_ok,
                plus_error,
                minus_error
            ) = best_domain

            mrv_hist[best_ti] = (
                mrv_hist.get(
                    best_ti,
                    0
                )
                + 1
            )

            # -----------------------------------------------
            # Construct orientation-specific candidates.
            # -----------------------------------------------

            plus_ids = np.flatnonzero(
                plus_ok
            )

            minus_ids = np.flatnonzero(
                minus_ok
            )

            ids = np.concatenate(
                [
                    plus_ids,
                    minus_ids
                ]
            )

            orientations = np.concatenate(
                [
                    np.ones(
                        len(plus_ids),
                        dtype=np.int8
                    ),
                    -np.ones(
                        len(minus_ids),
                        dtype=np.int8
                    )
                ]
            )

            local_scores = np.concatenate(
                [
                    (
                        plus_error[plus_ids]
                        +
                        norm_error[plus_ids]
                    ),
                    (
                        minus_error[minus_ids]
                        +
                        norm_error[minus_ids]
                    )
                ]
            )

            if not len(ids):
                continue

            take = min(
                args.branch,
                len(ids)
            )

            if len(ids) > take:

                pick = np.argpartition(
                    local_scores,
                    take - 1
                )[:take]

                ids = ids[pick]
                orientations = orientations[pick]
                local_scores = local_scores[pick]

            order = np.argsort(
                local_scores
            )

            ids = ids[order]
            orientations = orientations[order]
            local_scores = local_scores[order]

            for (
                ci,
                si,
                local_score
            ) in zip(
                ids,
                orientations,
                local_scores
            ):

                na = dict(
                    assignment
                )

                ns = dict(
                    signs
                )

                na[best_ti] = int(ci)
                ns[best_ti] = int(si)

                nh = list(
                    history
                )

                nh.append(
                    best_ti
                )

                children.append(
                    (
                        parent_score
                        +
                        float(local_score),
                        na,
                        ns,
                        nh
                    )
                )

        if not children:

            print(
                f"DEAD_END"
                f"|h={h:.8f}"
                f"|attempt_depth={depth+1}"
                f"|states={len(states)}"
                f"|zero_domain_states={zero_domain_states}",
                flush=True
            )

            break

        children.sort(
            key=lambda x: x[0]
        )

        # Deduplicate exact partial embeddings.
        seen = set()
        states = []

        for child in children:

            (
                score,
                assignment,
                signs,
                history
            ) = child

            key = tuple(
                sorted(
                    (
                        ti,
                        assignment[ti],
                        signs[ti]
                    )
                    for ti in assignment
                )
            )

            if key in seen:
                continue

            seen.add(key)
            states.append(
                child
            )

            if len(states) >= args.beam:
                break

        best_depth = depth + 1
        best_state = states[0]

        print(
            f"BEST|h={h:.8f}"
            f"|depth={best_depth}"
            f"|score={best_state[0]:.10g}"
            f"|states={len(states)}"
            f"|zero_domain_states={zero_domain_states}"
            f"|mrv={mrv_hist}"
            f"|order={best_state[3]}"
            f"|seconds={time.time()-started:.3f}",
            flush=True
        )

        if best_depth == 17:
            break

    # --------------------------------------------------------
    # Independent residual calculation.
    # --------------------------------------------------------

    (
        score,
        assignment,
        signs,
        history
    ) = best_state

    targets = sorted(
        assignment
    )

    A = np.array(
        [
            signs[ti]
            *
            W[
                assignment[ti]
            ]
            for ti in targets
        ],
        dtype=np.int64
    )

    G = (
        A @ H
    ) @ A.T

    QT = Q[
        np.ix_(
            targets,
            targets
        )
    ]

    residual = (
        G
        -
        h * QT
    )

    rel = (
        np.linalg.norm(
            residual
        )
        /
        (
            h
            *
            np.linalg.norm(
                QT
            )
        )
    )

    elapsed = (
        time.time()
        -
        started
    )

    print()
    print(
        f"H_FINAL"
        f"|h={h:.8f}"
        f"|depth={best_depth}"
        f"|score={score:.10g}"
        f"|rel={rel:.12g}"
        f"|order={history}"
        f"|seconds={elapsed:.3f}"
    )

    print(
        "assignment =",
        assignment
    )

    print(
        "signs =",
        signs
    )

    print(
        "selected:"
    )

    for ti in targets:

        ci = assignment[ti]

        print(
            f" target={ti}"
            f" shell={ci}"
            f" global={int(gids[ci])}"
            f" sign={signs[ti]}"
            f" norm={norms[ci]:.12g}"
        )

    result = {
        "h": h,
        "depth": best_depth,
        "score": score,
        "rel": rel,
        "order": history,
        "assignment": assignment,
        "signs": signs,
        "A": A,
        "residual": residual,
        "elapsed": elapsed
    }

    if (
        GLOBAL_BEST is None
        or best_depth > GLOBAL_BEST["depth"]
        or (
            best_depth == GLOBAL_BEST["depth"]
            and rel < GLOBAL_BEST["rel"]
        )
    ):

        GLOBAL_BEST = result

        np.savetxt(
            OUT /
            "rank17-rank21-mrv-best-A.txt",
            A,
            fmt="%d"
        )

        np.savetxt(
            OUT /
            "rank17-rank21-mrv-best-residual.txt",
            residual,
            fmt="%.17g"
        )

    return result


# ============================================================
# Sweep h.
# ============================================================

hs = []

h = args.h_min

while h <= args.h_max + 1e-12:

    hs.append(
        round(
            h,
            12
        )
    )

    h += args.h_step


for h in hs:

    result = run_h(
        h
    )

    if result is not None:

        summary.append(
            result
        )


print()
print("=" * 76)
print("SWEEP SUMMARY")

for r in sorted(
    summary,
    key=lambda x: (
        -x["depth"],
        x["rel"]
    )
):

    print(
        f"h={r['h']:.8f}"
        f" depth={r['depth']}"
        f" rel={r['rel']:.10g}"
        f" score={r['score']:.10g}"
        f" order={r['order']}"
    )


if GLOBAL_BEST is not None:

    print()
    print("GLOBAL BEST")

    print(
        f"h={GLOBAL_BEST['h']:.8f}"
    )

    print(
        f"depth={GLOBAL_BEST['depth']}"
    )

    print(
        f"rel={GLOBAL_BEST['rel']:.12g}"
    )

    print(
        "order=",
        GLOBAL_BEST["order"]
    )

