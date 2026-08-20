from pathlib import Path
import argparse
import random
import time
import numpy as np
from sage.all import Matrix, RDF

parser = argparse.ArgumentParser()

parser.add_argument("gram")
parser.add_argument("--samples", type=int, default=100_000_000)
parser.add_argument("--keep", type=int, default=1_000_000)
parser.add_argument("--walkers", type=int, default=32)
parser.add_argument("--seed", type=int, default=290017)
parser.add_argument("--cutoff", type=float, default=350.0)
parser.add_argument("--max-coeff", type=int, default=20)
parser.add_argument("--checkpoint", type=int, default=5_000_000)

args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results"
OUT.mkdir(parents=True, exist_ok=True)

HA = np.loadtxt(args.gram, dtype=float)
r = HA.shape[0]

G = Matrix(RDF, HA.tolist())

print("LLL start", flush=True)
U = G.LLL_gram()
print("LLL done", flush=True)

R = np.array(
    U.transpose() * G * U,
    dtype=float
)

U_np = np.array(U, dtype=np.int64)

rng = random.Random(args.seed)

# Adaptive cutoff starts at requested maximum, but tightens when
# the candidate dictionary becomes too large.
active_cutoff = args.cutoff

# canonical reduced-coordinate vector -> norm
best = {}

# Cheap restart reservoir.
reservoir = []
reservoir_limit = 200_000


def canonical(z):
    nz = np.flatnonzero(z)

    if not len(nz):
        return None

    z = z.copy()

    if z[nz[0]] < 0:
        z = -z

    return tuple(map(int, z))


def norm_of(z):
    return float(z @ R @ z)


def consider(z):
    key = canonical(z)

    if key is None:
        return False

    if key in best:
        return False

    a = np.asarray(key, dtype=np.int64)
    n = norm_of(a)

    if n > active_cutoff:
        return False

    best[key] = n

    # Reservoir sampling-ish restart set.
    if len(reservoir) < reservoir_limit:
        reservoir.append(key)
    else:
        # Slowly replace old entries.
        if rng.random() < reservoir_limit / len(best):
            reservoir[rng.randrange(reservoir_limit)] = key

    return True


# ------------------------------------------------------------
# Seed pool
# ------------------------------------------------------------

for i in range(r):
    z = np.zeros(r, dtype=np.int64)
    z[i] = 1
    consider(z)

for i in range(r):
    for j in range(i + 1, r):
        for s in (-1, 1):
            z = np.zeros(r, dtype=np.int64)
            z[i] = 1
            z[j] = s
            consider(z)

for i in range(r):
    for j in range(i + 1, r):
        for k in range(j + 1, r):
            for sj in (-1, 1):
                for sk in (-1, 1):
                    z = np.zeros(r, dtype=np.int64)
                    z[i] = 1
                    z[j] = sj
                    z[k] = sk
                    consider(z)


# ------------------------------------------------------------
# Independent walkers
# ------------------------------------------------------------

walkers = []
walker_norms = []

for _ in range(args.walkers):
    key = reservoir[rng.randrange(len(reservoir))]
    z = np.array(key, dtype=np.int64)

    walkers.append(z)
    walker_norms.append(norm_of(z))


# ------------------------------------------------------------
# Checkpoint
# ------------------------------------------------------------

def save_checkpoint(step, final=False):
    count = min(args.keep, len(best))

    items = sorted(
        best.items(),
        key=lambda x: x[1]
    )[:count]

    Z = np.array(
        [np.array(k, dtype=np.int64) for k, n in items],
        dtype=np.int64
    )

    N = np.array(
        [n for k, n in items],
        dtype=float
    )

    V = Z @ U_np.T

    suffix = "final" if final else f"step{step}"

    np.save(
        OUT / f"E29-short-pool-v2-vectors-{suffix}.npy",
        V
    )

    np.save(
        OUT / f"E29-short-pool-v2-norms-{suffix}.npy",
        N
    )

    # Stable aliases always point to latest checkpoint.
    np.save(
        OUT / "E29-short-pool-v2-vectors.npy",
        V
    )

    np.save(
        OUT / "E29-short-pool-v2-norms.npy",
        N
    )

    print(
        f"CHECKPOINT|step={step}"
        f"|unique={len(best)}"
        f"|saved={len(V)}"
        f"|min={N.min():.6g}"
        f"|max={N.max():.6g}",
        flush=True
    )


# ------------------------------------------------------------
# Exploration
# ------------------------------------------------------------

start = time.time()
accepted = 0

for step in range(1, args.samples + 1):

    wi = step % args.walkers

    current = walkers[wi]
    current_norm = walker_norms[wi]

    # Independent restart.
    if rng.random() < 0.003 and reservoir:
        key = reservoir[rng.randrange(len(reservoir))]
        current = np.array(key, dtype=np.int64)
        current_norm = norm_of(current)

    p = current.copy()

    move = rng.randrange(100)

    if move < 55:
        # local coordinate move
        i = rng.randrange(r)
        p[i] += rng.choice((-1, 1))

    elif move < 70:
        i = rng.randrange(r)
        p[i] += rng.choice((-2, 2))

    elif move < 90:
        # shear using another coordinate
        i, j = rng.sample(range(r), 2)
        p[i] += rng.choice((-1, 1)) * p[j]

    elif move < 97:
        # two simultaneous mutations
        i, j = rng.sample(range(r), 2)
        p[i] += rng.choice((-1, 1))
        p[j] += rng.choice((-1, 1))

    else:
        # bigger rare jump
        for _ in range(rng.randint(2, 5)):
            i = rng.randrange(r)
            p[i] += rng.choice((-1, 1))

    if np.max(np.abs(p)) > args.max_coeff:
        continue

    if not np.any(p):
        continue

    pn = norm_of(p)

    if pn <= active_cutoff:
        if consider(p):
            accepted += 1

        # Annealed preference for short vectors.
        #
        # Temperature cycles so walkers alternate between
        # exploitation and exploration.
        phase = (step % 2_000_000) / 2_000_000
        temperature = 8.0 + 40.0 * abs(2.0 * phase - 1.0)

        delta = pn - current_norm

        if (
            delta <= 0
            or rng.random() < np.exp(-delta / temperature)
        ):
            walkers[wi] = p
            walker_norms[wi] = pn

    if step % 100_000 == 0:
        elapsed = time.time() - start

        rate = step / max(elapsed, 1e-9)

        vals = walker_norms

        print(
            f"STEP|n={step}"
            f"|unique={len(best)}"
            f"|accepted={accepted}"
            f"|rate={rate:.1f}/s"
            f"|walker_min={min(vals):.3f}"
            f"|walker_med={np.median(vals):.3f}"
            f"|walker_max={max(vals):.3f}",
            flush=True
        )

    # --------------------------------------------------------
    # Bound memory.
    #
    # Final output keeps args.keep vectors (3M in our run).
    # Let exploration grow beyond that, but if the dictionary
    # exceeds 2*keep, retain approximately the shortest
    # 4/3*keep and lower the active norm cutoff accordingly.
    #
    # np.partition avoids sorting millions of Python tuples.
    # --------------------------------------------------------

    prune_trigger = 2 * args.keep
    prune_target = max(
        args.keep,
        (4 * args.keep) // 3
    )

    if (
        step % 5_000_000 == 0
        and len(best) > prune_trigger
    ):
        print(
            f"PRUNE_START|step={step}"
            f"|unique={len(best)}"
            f"|active_cutoff={active_cutoff:.9g}",
            flush=True
        )

        values = np.fromiter(
            best.values(),
            dtype=np.float64,
            count=len(best)
        )

        kth = min(
            prune_target - 1,
            len(values) - 1
        )

        threshold = float(
            np.partition(values, kth)[kth]
        )

        del values

        # Tighten monotonically.
        active_cutoff = min(
            active_cutoff,
            threshold
        )

        best = {
            k: n
            for k, n in best.items()
            if n <= active_cutoff
        }

        # Rebuild the fixed-size restart reservoir from survivors.
        reservoir.clear()

        for k in best:
            reservoir.append(k)

            if len(reservoir) >= reservoir_limit:
                break

        # Any walker that now lies outside the retained region
        # is restarted from a surviving candidate.
        for wi in range(args.walkers):
            if walker_norms[wi] > active_cutoff:
                key = reservoir[
                    rng.randrange(len(reservoir))
                ]

                z = np.array(
                    key,
                    dtype=np.int64
                )

                walkers[wi] = z
                walker_norms[wi] = norm_of(z)

        print(
            f"PRUNE_DONE|step={step}"
            f"|unique={len(best)}"
            f"|active_cutoff={active_cutoff:.9g}"
            f"|reservoir={len(reservoir)}",
            flush=True
        )

    if (
        args.checkpoint > 0
        and step % args.checkpoint == 0
    ):
        save_checkpoint(step)


save_checkpoint(args.samples, final=True)

print()
print("DONE")
print("samples =", args.samples)
print("unique =", len(best))
print("accepted =", accepted)
