from sage.all import *
from sage.quadratic_forms.genera.genus import Genus
from pathlib import Path
import argparse
import hashlib
import math
import random


def read_hessian(path):
    return matrix(ZZ, [
        [ZZ(x) for x in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def form_from_hessian(H):
    # Q(x) = 1/2 x^T H x.
    n = H.nrows()
    coeffs = []

    for i in range(n):
        assert H[i,i] % 2 == 0
        for j in range(i, n):
            coeffs.append(
                H[i,i] // 2 if i == j
                else H[i,j]
            )

    Q = QuadraticForm(ZZ, n, coeffs)

    assert Q.Hessian_matrix() == H
    return Q


def gram_key(Q):
    H = Q.Hessian_matrix()
    return hashlib.sha256(
        repr(tuple(tuple(map(int,r)) for r in H.rows())).encode()
    ).hexdigest()


def fingerprint(Q):
    # Q=1 corresponds to lattice norm 2.
    # Q=2 corresponds to lattice norm 4.
    S = Q.short_vector_list_up_to_length(3, True)
    return len(S[1]), len(S[2])


def energy(fp):
    roots, n4 = fp

    # Rootlessness is primary, shell distance secondary.
    #
    # r=0,n4=1320 -> 9
    # r=1,n4=1311 -> 1000
    #
    # So once we find rootless target-genus territory we explore it
    # preferentially.
    return 8 * roots + abs(n4 - 1311)


def save_form(Q, path, **metadata):
    H = Q.Hessian_matrix()

    lines = []
    for k,v in metadata.items():
        lines.append(f"# {k} = {v}")

    lines.extend(
        " ".join(map(str,row))
        for row in H.rows()
    )

    Path(path).write_text("\n".join(lines) + "\n")


def genus_key(Q):
    g = Genus(Q.Hessian_matrix())
    syms = {
        int(s.prime()): str(s)
        for s in g.local_symbols()
    }
    return (
        syms.get(2,""),
        syms.get(3,""),
        syms.get(79,""),
    )


ap = argparse.ArgumentParser()

ap.add_argument("--gram", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--steps", type=int, default=10000)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--primes", default="5,7,11")
ap.add_argument("--restart-every", type=int, default=1000)
ap.add_argument("--temperature", type=float, default=20.0)
ap.add_argument("--genus-check-every", type=int, default=250)

args = ap.parse_args()

random.seed(args.seed)
set_random_seed(args.seed)

out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)

primes = [ZZ(x) for x in args.primes.split(",")]

H0 = read_hessian(args.gram)

assert H0.nrows() == 17
assert H0.det() == 948
assert all(H0[i,i] % 2 == 0 for i in range(17))

Q0 = form_from_hessian(H0)

start_genus = genus_key(Q0)
fp0 = fingerprint(Q0)

print(
    "PNEIGHBOR|stage=start"
    f"|seed={args.seed}"
    f"|roots={fp0[0]}"
    f"|n4={fp0[1]}"
    f"|det={H0.det()}"
    f"|primes={','.join(map(str,primes))}",
    flush=True
)

assert fp0 == (1,1311), fp0

current = Q0
current_fp = fp0
current_E = energy(fp0)

best = Q0
best_fp = fp0
best_E = current_E

visited = {gram_key(Q0)}

accepted = 0
duplicates = 0
failures = 0
evaluated = 0
rootless_seen = 0

save_form(
    Q0,
    out / "start-gram.txt",
    seed=args.seed,
    roots=fp0[0],
    n4=fp0[1],
)


for step in range(1, args.steps + 1):

    if (
        args.restart_every > 0
        and step > 1
        and step % args.restart_every == 0
    ):
        current = best
        current_fp = best_fp
        current_E = best_E

        print(
            "PNEIGHBOR|stage=restart"
            f"|step={step}"
            f"|best_roots={best_fp[0]}"
            f"|best_n4={best_fp[1]}",
            flush=True
        )

    p = random.choice(primes)

    try:
        v = current.find_primitive_p_divisible_vector__random(p)
        neighbor = current.find_p_neighbor_from_vec(p, v)

        HN = neighbor.Hessian_matrix()

        if HN.det() != 948:
            raise RuntimeError(
                f"determinant changed to {HN.det()}"
            )

        if not all(HN[i,i] % 2 == 0 for i in range(17)):
            raise RuntimeError("neighbor is not even")

    except Exception as e:
        failures += 1

        if failures <= 10:
            print(
                "PNEIGHBOR|stage=failure"
                f"|step={step}|p={p}|error={e}",
                flush=True
            )
        continue

    key = gram_key(neighbor)

    if key in visited:
        duplicates += 1
        continue

    visited.add(key)

    try:
        fp = fingerprint(neighbor)
    except Exception as e:
        failures += 1
        continue

    evaluated += 1

    roots,n4 = fp
    E = energy(fp)

    # Periodic hard assertion that the intrinsic neighbor walk really
    # remains in our starting genus.
    if (
        args.genus_check_every > 0
        and evaluated % args.genus_check_every == 0
    ):
        g = genus_key(neighbor)

        if g != start_genus:
            print(
                "PNEIGHBOR|stage=GENUS_ERROR"
                f"|step={step}|p={p}"
                f"|roots={roots}|n4={n4}",
                flush=True
            )
            raise RuntimeError("p-neighbor changed genus")

        print(
            "PNEIGHBOR|stage=genus_check"
            f"|step={step}|status=OK",
            flush=True
        )

    # --------------------------------------------------------
    # JACKPOT
    # --------------------------------------------------------

    if roots == 0:
        rootless_seen += 1

        save_form(
            neighbor,
            out / f"rootless-step{step}-p{p}-n{n4}.txt",
            seed=args.seed,
            step=step,
            prime=p,
            roots=roots,
            n4=n4,
        )

        print(
            "PNEIGHBOR|stage=rootless"
            f"|step={step}|p={p}|n4={n4}",
            flush=True
        )

    if fp == (0,1311):

        jackpot = out / f"JACKPOT-step{step}-p{p}.txt"

        save_form(
            neighbor,
            jackpot,
            seed=args.seed,
            step=step,
            prime=p,
            roots=0,
            n4=1311,
        )

        # Final genus assertion.
        assert genus_key(neighbor) == start_genus

        print(
            "PNEIGHBOR|stage=JACKPOT"
            f"|step={step}|p={p}"
            f"|file={jackpot}",
            flush=True
        )

        raise SystemExit(0)

    # Save useful near-target rooted states.
    if roots <= 1 and abs(n4 - 1311) <= 6:
        save_form(
            neighbor,
            out / f"near-step{step}-p{p}-r{roots}-n{n4}.txt",
            seed=args.seed,
            step=step,
            prime=p,
            roots=roots,
            n4=n4,
        )

    # --------------------------------------------------------
    # BEST
    # --------------------------------------------------------

    if E < best_E:
        best = neighbor
        best_fp = fp
        best_E = E

        save_form(
            best,
            out / "best-gram.txt",
            seed=args.seed,
            step=step,
            prime=p,
            roots=roots,
            n4=n4,
        )

        print(
            "PNEIGHBOR|stage=best"
            f"|step={step}|p={p}"
            f"|roots={roots}|n4={n4}"
            f"|energy={E}",
            flush=True
        )

    # --------------------------------------------------------
    # SA acceptance
    # --------------------------------------------------------

    dE = E - current_E

    if args.restart_every > 0:
        phase = step % args.restart_every
        frac = phase / args.restart_every
    else:
        frac = step / max(1,args.steps)

    T = max(
        0.1,
        args.temperature * (1.0-frac)
    )

    if dE <= 0:
        accept = True
    else:
        accept = (
            random.random()
            < math.exp(-float(dE)/T)
        )

    if accept:
        current = neighbor
        current_fp = fp
        current_E = E
        accepted += 1

    if step % 100 == 0:
        print(
            "PNEIGHBOR|stage=progress"
            f"|step={step}"
            f"|current={current_fp}"
            f"|best={best_fp}"
            f"|visited={len(visited)}"
            f"|evaluated={evaluated}"
            f"|accepted={accepted}"
            f"|rootless={rootless_seen}"
            f"|duplicates={duplicates}"
            f"|failures={failures}",
            flush=True
        )


print(
    "PNEIGHBOR|stage=done"
    f"|best={best_fp}"
    f"|visited={len(visited)}"
    f"|evaluated={evaluated}"
    f"|rootless={rootless_seen}"
    f"|accepted={accepted}"
    f"|duplicates={duplicates}"
    f"|failures={failures}",
    flush=True
)
