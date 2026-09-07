from pathlib import Path
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("vectors")
parser.add_argument("norms")

parser.add_argument("--trials", type=int, default=500)
parser.add_argument("--slice", type=int, default=4000)
parser.add_argument("--norm-tol", type=float, default=0.035)
parser.add_argument("--pair-tol", type=float, default=0.10)
parser.add_argument("--seed", type=int, default=170021)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

Q = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

H = np.loadtxt(
    args.gram,
    dtype=float
)

V = np.load(
    args.vectors,
    mmap_mode="r"
)

N = np.load(
    args.norms,
    mmap_mode="r"
)

assert Q.shape == (17,17)
assert np.all(np.diag(Q) == 4)

assert V.ndim == 2
assert H.shape == (V.shape[1], V.shape[1])
assert len(V) == len(N)

rng = np.random.default_rng(args.seed)

print("ambient_dim =", V.shape[1])
print("pool =", len(V))
print("norm min =", float(N.min()))
print("norm max =", float(N.max()))

# ------------------------------------------------------------
# Target ordering.
#
# Start with a very constrained row, then greedily add the
# target having the most nonzero constraints to what is already
# placed.
# ------------------------------------------------------------

anchor = max(
    range(17),
    key=lambda i: (
        np.count_nonzero(Q[i]),
        np.sum(np.abs(Q[i]))
    )
)

order = [anchor]
remaining = set(range(17))
remaining.remove(anchor)

while remaining:
    j = max(
        remaining,
        key=lambda x: (
            sum(Q[x,k] != 0 for k in order),
            sum(abs(int(Q[x,k])) for k in order),
            np.count_nonzero(Q[x])
        )
    )

    order.append(j)
    remaining.remove(j)

print("target order =", order)

# Sorted norms let us pull narrow shells cheaply.
perm = np.argsort(N)
SN = np.asarray(N[perm])

best_depth = 0
best_trial = None
best_h = None
best_assignment = None

start = time.time()


def window(center):
    lo = center * (1.0 - args.norm_tol)
    hi = center * (1.0 + args.norm_tol)

    a = np.searchsorted(SN, lo, side="left")
    b = np.searchsorted(SN, hi, side="right")

    return perm[a:b]


for trial in range(args.trials):

    # Choose a random candidate line. Its norm implies h.
    anchor_global = int(
        rng.integers(len(V))
    )

    center = float(
        N[anchor_global]
    )

    h = center / 4.0

    shell = window(center)

    if len(shell) < 17:
        continue

    # Don't always take the same part of a giant shell.
    if len(shell) > args.slice:
        chosen = rng.choice(
            shell,
            size=args.slice,
            replace=False
        )

        # Ensure our anchor itself is present.
        chosen[0] = anchor_global
        shell = chosen

    else:
        shell = np.asarray(shell)

    W = np.asarray(
        V[shell],
        dtype=np.int64
    )

    # Put anchor at index zero if possible.
    where = np.flatnonzero(
        shell == anchor_global
    )

    if not len(where):
        continue

    ai = int(where[0])

    # Precompute W*H once.
    WH = W @ H

    assignment = {
        order[0]: (
            ai,
            1
        )
    }

    used = {ai}

    # --------------------------------------------------------
    # Greedy + limited branching.
    #
    # State:
    #   assignment target -> (ambient row, sign)
    # --------------------------------------------------------

    states = [
        (
            0.0,
            assignment,
            used
        )
    ]

    for depth in range(1,17):

        ti = order[depth]

        next_states = []

        # Cap beam: enough diversity without exploding.
        states = sorted(
            states,
            key=lambda x:x[0]
        )[:24]

        for score, ass, usedset in states:

            prev_targets = list(
                ass.keys()
            )

            A = np.array(
                [
                    sign * W[idx]
                    for idx,sign in (
                        ass[t]
                        for t in prev_targets
                    )
                ],
                dtype=np.int64
            )

            # Pair every shell vector with already selected
            # oriented vectors.
            P = WH @ A.T

            target = np.array(
                [
                    h * Q[ti,t]
                    for t in prev_targets
                ],
                dtype=float
            )

            # Try each candidate with both possible signs.
            err_plus = np.max(
                np.abs(P - target),
                axis=1
            )

            err_minus = np.max(
                np.abs(-P - target),
                axis=1
            )

            e = np.minimum(
                err_plus,
                err_minus
            )

            signs = np.where(
                err_plus <= err_minus,
                1,
                -1
            )

            valid = np.flatnonzero(
                e <= args.pair_tol * h
            )

            if not len(valid):
                continue

            # Don't reuse the same unoriented lattice line.
            valid = [
                int(x)
                for x in valid
                if int(x) not in usedset
            ]

            if not valid:
                continue

            valid.sort(
                key=lambda x:e[x]
            )

            # Branch only on best few.
            for x in valid[:8]:

                na = dict(ass)
                na[ti] = (
                    x,
                    int(signs[x])
                )

                nu = set(usedset)
                nu.add(x)

                ns = (
                    score
                    + float(
                        e[x] / max(h,1e-12)
                    )
                )

                next_states.append(
                    (
                        ns,
                        na,
                        nu
                    )
                )

        states = next_states

        reached = depth + 1

        if states and reached > best_depth:

            best_depth = reached
            winner = min(
                states,
                key=lambda x:x[0]
            )

            best_trial = trial
            best_h = h
            best_assignment = winner[1]

            print(
                f"BEST|trial={trial}"
                f"|depth={best_depth}"
                f"|h={h:.12g}"
                f"|shell={len(shell)}"
                f"|score={winner[0]:.8g}",
                flush=True
            )

        if not states:
            break

    if best_depth == 17:

        winner = min(
            states,
            key=lambda x:x[0]
        )

        ass = winner[1]

        A = np.array(
            [
                ass[i][1] * W[ass[i][0]]
                for i in range(17)
            ],
            dtype=np.int64
        )

        G = (A @ H) @ A.T
        R = G - h * Q

        rel = (
            np.linalg.norm(R)
            /
            (
                h * np.linalg.norm(Q)
            )
        )

        print()
        print("FOUND")
        print("h =",h)
        print("relative residual =",rel)
        print("A =")
        print(A)

        break


print()
print("FINAL")
print("best_depth =", best_depth)
print("best_trial =", best_trial)
print("best_h =", best_h)
print("seconds =", time.time()-start)

if best_assignment is not None:
    print("best_assignment =",best_assignment)
