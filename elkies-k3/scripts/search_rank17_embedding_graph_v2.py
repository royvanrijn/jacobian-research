from pathlib import Path
import argparse
import math
import time
from collections import defaultdict, Counter

import numpy as np
from sage.all import Matrix, RDF


parser = argparse.ArgumentParser()

parser.add_argument("gram")

parser.add_argument(
    "--limit",
    type=int,
    default=15000,
    help="maximum unoriented short-vector lines"
)

parser.add_argument(
    "--norm-tol",
    type=float,
    default=0.02,
    help="relative tolerance for equal target norms"
)

parser.add_argument(
    "--pair-tol",
    type=float,
    default=0.06,
    help="pairing error tolerance in units of h"
)

parser.add_argument(
    "--max-shells",
    type=int,
    default=100
)

parser.add_argument(
    "--max-shell-lines",
    type=int,
    default=1400,
    help="skip shells larger than this before +/- expansion"
)

parser.add_argument(
    "--node-limit",
    type=int,
    default=2_000_000,
    help="maximum DFS nodes per shell"
)

parser.add_argument(
    "--seconds-per-shell",
    type=float,
    default=180.0
)

args = parser.parse_args()


BASE = Path(__file__).resolve().parents[1]

J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "results"

OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# Target rank-17 lattice
# ============================================================

Q = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

assert Q.shape == (17, 17)
assert np.all(np.diag(Q) == 4)

target_values = sorted(
    set(Q[np.triu_indices(17, 1)].tolist())
)

print(
    "target offdiag values =",
    target_values
)

assert set(target_values).issubset(
    {-2, -1, 0, 1, 2}
)


# ============================================================
# Choose a good target anchor.
#
# Prefer a row with a distinctive pairing profile.
# ============================================================

profiles = []

for i in range(17):
    vals = [
        int(Q[i, j])
        for j in range(17)
        if j != i
    ]

    c = Counter(vals)

    profiles.append(c)


global_counts = Counter()

for i in range(17):
    for j in range(i):
        global_counts[int(Q[i, j])] += 1


def anchor_score(i):
    # Larger = more restrictive.
    c = profiles[i]

    score = 0.0

    for q, n in c.items():
        rarity = 1.0 / max(global_counts[q], 1)
        score += n * rarity

    # Non-zero constraints generally prune more.
    score += 0.1 * sum(
        n
        for q, n in c.items()
        if q != 0
    )

    return score


anchor = max(
    range(17),
    key=anchor_score
)

print(
    "anchor target =",
    anchor,
    "profile =",
    dict(sorted(profiles[anchor].items()))
)


# ============================================================
# Ambient lattice
# ============================================================

HA_np = np.loadtxt(
    args.gram,
    dtype=float
)

r = HA_np.shape[0]

assert HA_np.shape == (r, r)

print()
print(f"ambient_rank={r}")

HA = Matrix(
    RDF,
    HA_np.tolist()
)

print("LLL start", flush=True)

U = HA.LLL_gram()

print("LLL done", flush=True)

R = U.transpose() * HA * U

R_np = np.array(
    R,
    dtype=float
)

U_np = np.array(
    U,
    dtype=np.int64
)


# ============================================================
# Generate short ambient vector LINES.
#
# We canonicalize +/- here only to avoid duplicate lines.
# Later every selected line is expanded to BOTH orientations.
# ============================================================

cand = {}


def canonical(z):
    z = np.asarray(
        z,
        dtype=np.int64
    )

    nz = np.flatnonzero(z)

    if not len(nz):
        return None

    if z[nz[0]] < 0:
        z = -z

    return tuple(map(int, z))


def add(z):
    key = canonical(z)

    if key is None:
        return

    if key in cand:
        return

    z = np.array(
        key,
        dtype=np.int64
    )

    norm = float(
        z @ R_np @ z
    )

    if norm <= 0:
        return

    cand[key] = norm


# Single basis vectors.
for i in range(r):

    z = np.zeros(
        r,
        dtype=np.int64
    )

    z[i] = 1
    add(z)

    z = np.zeros(
        r,
        dtype=np.int64
    )

    z[i] = 2
    add(z)


# Pairs.
for i in range(r):

    for j in range(i + 1, r):

        for s in (-1, 1):

            z = np.zeros(
                r,
                dtype=np.int64
            )

            z[i] = 1
            z[j] = s

            add(z)


# Triples.
for i in range(r):

    for j in range(i + 1, r):

        for k in range(j + 1, r):

            for sj in (-1, 1):

                for sk in (-1, 1):

                    z = np.zeros(
                        r,
                        dtype=np.int64
                    )

                    z[i] = 1
                    z[j] = sj
                    z[k] = sk

                    add(z)


items = sorted(
    cand.items(),
    key=lambda x: x[1]
)[:args.limit]

Z = np.array(
    [
        np.array(
            key,
            dtype=np.int64
        )
        for key, norm in items
    ]
)

norms = np.array(
    [
        norm
        for key, norm in items
    ],
    dtype=float
)

# Reduced -> original ambient coordinates.
V = Z @ U_np.T


print()
print(
    "unoriented candidate lines =",
    len(V)
)

print(
    "norm range =",
    float(norms.min()),
    float(norms.max())
)


# ============================================================
# Build narrow norm shells.
#
# Use logarithmic bins and two offsets so candidates close to
# a bin boundary are not systematically separated.
# ============================================================

width = math.log1p(
    args.norm_tol
)

shell_map = {}


for shift in (0.0, 0.5):

    bins = defaultdict(list)

    for i, norm in enumerate(norms):

        b = math.floor(
            math.log(norm) / width
            + shift
        )

        bins[b].append(i)

    for inds in bins.values():

        # Need at least 17 distinct lines.
        if len(inds) < 17:
            continue

        inds = tuple(sorted(inds))

        shell_map[inds] = True


shells = list(
    shell_map.keys()
)


# Small shells are much cheaper and more informative.
shells.sort(
    key=lambda inds: (
        len(inds),
        np.median(norms[list(inds)])
    )
)


eligible = []

skipped_large = 0

for inds in shells:

    if len(inds) > args.max_shell_lines:
        skipped_large += 1
        continue

    eligible.append(inds)


eligible = eligible[
    :args.max_shells
]


print()
print(
    "shells total =",
    len(shells)
)

print(
    "shells skipped_large =",
    skipped_large
)

print(
    "shells searched =",
    len(eligible)
)


# ============================================================
# Bitset helpers
# ============================================================

def bool_row_to_mask(row):
    packed = np.packbits(
        row,
        bitorder="little"
    )

    return int.from_bytes(
        packed.tobytes(),
        byteorder="little"
    )


def iter_bits(mask):

    while mask:

        lsb = mask & -mask

        i = lsb.bit_length() - 1

        yield i

        mask ^= lsb


# ============================================================
# Global best / checkpointing
# ============================================================

best_depth = 0
best_payload = None

BEST_PATH = (
    OUT
    / "rank17-E29-best-partial.txt"
)


def save_best(
    depth,
    h,
    assignment,
    W,
    shell_index,
    nodes
):

    global best_depth
    global best_payload

    if depth <= best_depth:
        return

    best_depth = depth

    selected = []

    for ti, ci in sorted(
        assignment.items()
    ):
        selected.append(
            (
                ti,
                ci,
                W[ci].copy()
            )
        )

    best_payload = (
        h,
        selected
    )

    print(
        f"BEST|shell={shell_index}"
        f"|depth={depth}"
        f"|h={h:.12g}"
        f"|nodes={nodes}"
        f"|assignment="
        + str({
            ti: ci
            for ti, ci in assignment.items()
        }),
        flush=True
    )

    with BEST_PATH.open("w") as f:

        f.write(
            f"depth {depth}\n"
        )

        f.write(
            f"h {h:.17g}\n"
        )

        f.write(
            f"shell {shell_index}\n"
        )

        f.write(
            f"nodes {nodes}\n"
        )

        for ti, ci, v in selected:

            f.write(
                f"target {ti}"
                f" candidate {ci}"
                f" vector "
                + " ".join(
                    map(str, v)
                )
                + "\n"
            )


# ============================================================
# Search shells
# ============================================================

total_nodes = 0


for shell_index, line_inds in enumerate(
    eligible
):

    line_inds = np.array(
        line_inds,
        dtype=np.int64
    )

    line_norms = norms[
        line_inds
    ]

    center = float(
        np.median(
            line_norms
        )
    )

    h = center / 4.0

    # Tight norm filter around the shell center.
    keep = np.abs(
        line_norms / center - 1.0
    ) <= args.norm_tol

    line_inds = line_inds[
        keep
    ]

    if len(line_inds) < 17:
        continue


    # --------------------------------------------------------
    # BOTH orientations.
    # --------------------------------------------------------

    W0 = V[
        line_inds
    ]

    W = np.vstack(
        [
            W0,
            -W0
        ]
    )

    m = len(W)

    print()
    print(
        f"SHELL|i={shell_index}"
        f"|norm={center:.12g}"
        f"|h={h:.12g}"
        f"|lines={len(W0)}"
        f"|oriented={m}",
        flush=True
    )


    # --------------------------------------------------------
    # Full pairing matrix for this shell.
    # --------------------------------------------------------

    P = (
        W @ HA_np
    ) @ W.T


    # Norm-valid candidate bitmask.
    diag = np.diag(P)

    norm_ok = (
        np.abs(
            diag - 4.0 * h
        )
        <= args.norm_tol
        * 4.0
        * h
    )

    base_mask = bool_row_to_mask(
        norm_ok
    )

    if base_mask.bit_count() < 17:

        print(
            "SKIP|too_few_norm_candidates",
            flush=True
        )

        continue


    # --------------------------------------------------------
    # Precompute compatibility masks:
    #
    # compat[q][candidate_i]
    #
    # is a Python integer bitset of every candidate j having
    #
    #     <i,j> ~= h*q
    #
    # --------------------------------------------------------

    compat = {}

    for q in target_values:

        wanted = h * q

        A = (
            np.abs(
                P - wanted
            )
            <= args.pair_tol * h
        )

        # Never map two target vertices to same candidate.
        np.fill_diagonal(
            A,
            False
        )

        compat[q] = [
            bool_row_to_mask(
                A[i]
            )
            for i in range(m)
        ]


    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    assignment = {}

    used_mask = 0

    shell_nodes = 0

    shell_start = time.monotonic()

    aborted = False


    # Global sign symmetry:
    #
    # W[:len(W0)] is our canonical orientation.
    # Force anchor image into that half.
    #
    positive_mask = (
        (1 << len(W0)) - 1
    )


    def domain_for(ti):

        mask = (
            base_mask
            & ~used_mask
        )

        for tj, cj in assignment.items():

            q = int(
                Q[ti, tj]
            )

            mask &= compat[q][cj]

            if not mask:
                break

        return mask


    def search():

        nonlocal_vars = None

        # Python nested mutation wrappers.
        global total_nodes

        # Mutable values via outer list.
        state[0] += 1
        total_nodes += 1

        nodes = state[0]

        if nodes > args.node_limit:
            state[1] = True
            return

        if (
            time.monotonic()
            - shell_start
            > args.seconds_per_shell
        ):
            state[1] = True
            return

        depth = len(
            assignment
        )

        save_best(
            depth,
            h,
            assignment,
            W,
            shell_index,
            nodes
        )

        if depth == 17:

            Arows = np.array(
                [
                    W[
                        assignment[i]
                    ]
                    for i in range(17)
                ],
                dtype=np.int64
            )

            G = (
                Arows @ HA_np
            ) @ Arows.T

            residual = (
                G - h * Q
            )

            rel = (
                np.linalg.norm(
                    residual
                )
                /
                (
                    h
                    * np.linalg.norm(Q)
                )
            )

            max_rel = (
                np.max(
                    np.abs(
                        residual
                    )
                )
                / h
            )

            print()
            print(
                "FOUND FULL EMBEDDING",
                flush=True
            )

            print(
                f"h={h:.17g}"
            )

            print(
                f"relative_residual={rel:.17g}"
            )

            print(
                f"max_residual_over_h={max_rel:.17g}"
            )

            print("A=")
            print(Arows)

            np.savetxt(
                OUT
                / "rank17-E29-embedding-A.txt",
                Arows,
                fmt="%d"
            )

            np.savetxt(
                OUT
                / "rank17-E29-embedding-gram.txt",
                G,
                fmt="%.17g"
            )

            np.savetxt(
                OUT
                / "rank17-E29-embedding-residual.txt",
                residual,
                fmt="%.17g"
            )

            raise SystemExit(0)


        # ----------------------------------------------------
        # MRV:
        #
        # choose whichever unassigned target vertex currently
        # has the smallest domain.
        # ----------------------------------------------------

        best_t = None
        best_domain = None
        best_count = None

        for ti in range(17):

            if ti in assignment:
                continue

            D = domain_for(
                ti
            )

            # Fix global sign using the anchor.
            if (
                depth == 0
                and ti == anchor
            ):
                D &= positive_mask

            count = (
                D.bit_count()
            )

            if count == 0:
                return

            if (
                best_count is None
                or count < best_count
            ):
                best_t = ti
                best_domain = D
                best_count = count


        ti = best_t


        # ----------------------------------------------------
        # Candidate ordering:
        #
        # prefer candidates that leave the smallest number of
        # possibilities for future target vertices.
        #
        # This is a cheap "fail-first" heuristic.
        # ----------------------------------------------------

        choices = []

        for ci in iter_bits(
            best_domain
        ):

            score = 0

            bit = (
                1 << ci
            )

            # Temporarily evaluate future compatibility.
            for tj in range(17):

                if (
                    tj == ti
                    or tj in assignment
                ):
                    continue

                q = int(
                    Q[ti, tj]
                )

                score += (
                    compat[q][ci]
                    & base_mask
                    & ~used_mask
                ).bit_count()

            choices.append(
                (
                    score,
                    ci
                )
            )


        # Smaller continuation space first.
        choices.sort()


        for _, ci in choices:

            if state[1]:
                return

            assignment[
                ti
            ] = ci

            old_used = (
                state[2]
            )

            state[2] |= (
                1 << ci
            )

            # used_mask is accessed through state below.
            search()

            state[2] = old_used

            del assignment[
                ti
            ]


    # Work around Python nonlocal assignment to integer mask:
    #
    # state[0] = nodes
    # state[1] = aborted
    # state[2] = used_mask
    #
    state = [
        0,
        False,
        0
    ]


    # Override helper to read state mask.
    def domain_for(ti):

        mask = (
            base_mask
            & ~state[2]
        )

        for tj, cj in assignment.items():

            q = int(
                Q[ti, tj]
            )

            mask &= compat[q][cj]

            if not mask:
                break

        return mask


    search()


    elapsed = (
        time.monotonic()
        - shell_start
    )

    print(
        f"SHELL_DONE|i={shell_index}"
        f"|nodes={state[0]}"
        f"|seconds={elapsed:.3f}"
        f"|aborted={state[1]}"
        f"|global_best_depth={best_depth}",
        flush=True
    )


print()
print("NO FULL EMBEDDING")
print(
    "best_depth =",
    best_depth
)
print(
    "total_nodes =",
    total_nodes
)
print(
    "best partial saved to",
    BEST_PATH
)
