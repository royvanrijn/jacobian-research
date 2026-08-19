from sage.all import *
from pathlib import Path
import argparse
import hashlib
import random


def read_hessian(path):
    return matrix(ZZ, [
        [ZZ(x) for x in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def form_from_hessian(H):
    n = H.nrows()
    coeffs = []

    for i in range(n):
        assert H[i,i] % 2 == 0
        for j in range(i,n):
            coeffs.append(
                H[i,i]//2 if i == j
                else H[i,j]
            )

    Q = QuadraticForm(ZZ,n,coeffs)
    assert Q.Hessian_matrix() == H
    return Q


def fingerprint(Q):
    S = Q.short_vector_list_up_to_length(3, True)
    return (len(S[1]), len(S[2]))


def raw_key(H):
    return hashlib.sha256(
        repr(tuple(tuple(map(int,r)) for r in H.rows())).encode()
    ).hexdigest()


def reduced_key(Q):
    """
    Hash a reduced equivalent quadratic form where possible.

    Reduction is not being used as a proof of isometry uniqueness;
    it is a strong practical deduplication heuristic.
    """
    try:
        R = Q.reduced_form()
        H = R.Hessian_matrix()
    except Exception:
        H = Q.Hessian_matrix()

    return raw_key(H)


def target_score(fp):
    r,n4 = fp

    if fp == (0,1311):
        return -10**9

    # Mild preference only.
    return 6*r + abs(n4-1311)


def save_form(Q,path,**meta):
    H = Q.Hessian_matrix()

    lines = [f"# {k} = {v}" for k,v in meta.items()]
    lines += [
        " ".join(map(str,row))
        for row in H.rows()
    ]

    Path(path).write_text("\n".join(lines)+"\n")


ap = argparse.ArgumentParser()

ap.add_argument("--gram-dir",required=True)
ap.add_argument("--out",required=True)
ap.add_argument("--generations",type=int,default=30)
ap.add_argument("--beam",type=int,default=80)
ap.add_argument("--samples-per-parent",type=int,default=120)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--primes",default="5,7,11,13,17,19,23")

# How much of the beam is reserved for target-ish states.
ap.add_argument("--elite",type=int,default=20)

# How much is purely random diversity.
ap.add_argument("--random-slots",type=int,default=30)

args = ap.parse_args()

random.seed(args.seed)
set_random_seed(args.seed)

out = Path(args.out)
out.mkdir(parents=True,exist_ok=True)

primes = [ZZ(x) for x in args.primes.split(",")]

frontier = []
seen = set()

for p in sorted(Path(args.gram_dir).glob("*.txt")):
    try:
        Q = form_from_hessian(read_hessian(p))
        H = Q.Hessian_matrix()

        assert H.det() == 948

        fp = fingerprint(Q)
        k = reduced_key(Q)

    except Exception as e:
        print("DIVERSITY|start_skip",p,e,flush=True)
        continue

    if k in seen:
        continue

    seen.add(k)
    frontier.append((Q,fp,k,str(p)))

print(
    f"DIVERSITY|start"
    f"|states={len(frontier)}"
    f"|primes={','.join(map(str,primes))}",
    flush=True
)

if not frontier:
    raise SystemExit("no starting states")


for gen in range(1,args.generations+1):

    generated = []
    rootless_hits = []

    for parent_idx,(Q,parent_fp,parent_key,parent_src) in enumerate(frontier):

        for sample in range(args.samples_per_parent):

            p = random.choice(primes)

            try:
                v = Q.find_primitive_p_divisible_vector__random(p)
                N = Q.find_p_neighbor_from_vec(p,v)

                H = N.Hessian_matrix()

                if H.det() != 948:
                    continue

                if not all(H[i,i] % 2 == 0 for i in range(17)):
                    continue

                k = reduced_key(N)

                if k in seen:
                    continue

                seen.add(k)

                fp = fingerprint(N)

            except Exception:
                continue

            r,n4 = fp

            if fp == (0,1311):
                path = out / (
                    f"JACKPOT-g{gen}-p{p}-{k[:12]}.txt"
                )

                save_form(
                    N,path,
                    generation=gen,
                    prime=p,
                    roots=0,
                    n4=1311,
                    parent=parent_fp,
                )

                print(
                    f"DIVERSITY|JACKPOT"
                    f"|gen={gen}|p={p}"
                    f"|parent={parent_fp}"
                    f"|file={path}",
                    flush=True
                )

                raise SystemExit(0)

            if r == 0:
                rootless_hits.append((n4,parent_fp,p,k))

                # Persist all rootless discoveries.
                save_form(
                    N,
                    out / (
                        f"rootless-g{gen}-p{p}"
                        f"-n{n4}-{k[:12]}.txt"
                    ),
                    generation=gen,
                    prime=p,
                    roots=0,
                    n4=n4,
                    parent=parent_fp,
                )

                print(
                    f"DIVERSITY|rootless"
                    f"|gen={gen}|p={p}"
                    f"|n4={n4}"
                    f"|parent={parent_fp}",
                    flush=True
                )

            generated.append(
                (
                    N,
                    fp,
                    k,
                    target_score(fp),
                )
            )

    if not generated:
        print(f"DIVERSITY|dead|gen={gen}",flush=True)
        break

    # --------------------------------------------------------
    # Reduce redundancy by shell fingerprint:
    # keep at most ONE state per (roots,n4) initially.
    # --------------------------------------------------------

    by_fp = {}

    random.shuffle(generated)

    for item in generated:
        Q,fp,k,score = item

        if fp not in by_fp:
            by_fp[fp] = item

    representatives = list(by_fp.values())

    # --------------------------------------------------------
    # 1. Elite near-target slice.
    # --------------------------------------------------------

    elite = sorted(
        representatives,
        key=lambda x:(
            x[3],
            x[1][0],
            abs(x[1][1]-1311)
        )
    )[:args.elite]

    chosen = []
    chosen_keys = set()

    def add(item):
        Q,fp,k,score = item
        if k in chosen_keys:
            return
        chosen_keys.add(k)
        chosen.append(item)

    for x in elite:
        add(x)

    # --------------------------------------------------------
    # 2. Root-count diversity.
    #
    # Preserve best states for several root-count buckets,
    # including deliberately ugly regions.
    # --------------------------------------------------------

    roots_buckets = [
        0,1,2,3,4,5,
        (6,10),
        (11,20),
        (21,1000),
    ]

    per_bucket = max(
        1,
        (args.beam - args.elite - args.random_slots)
        // len(roots_buckets)
    )

    for bucket in roots_buckets:

        if isinstance(bucket,tuple):
            lo,hi = bucket
            pool = [
                x for x in representatives
                if lo <= x[1][0] <= hi
            ]
        else:
            pool = [
                x for x in representatives
                if x[1][0] == bucket
            ]

        # Within an ugly bucket, prefer shell diversity rather than
        # closeness alone.
        random.shuffle(pool)

        # Sort only weakly by target distance.
        pool.sort(
            key=lambda x:abs(x[1][1]-1311)
        )

        for x in pool[:per_bucket]:
            add(x)

    # --------------------------------------------------------
    # 3. Pure random novelty.
    # --------------------------------------------------------

    remaining = [
        x for x in generated
        if x[2] not in chosen_keys
    ]

    random.shuffle(remaining)

    for x in remaining[:args.random_slots]:
        add(x)

    # --------------------------------------------------------
    # 4. Fill remaining slots with arbitrary novel states.
    # --------------------------------------------------------

    if len(chosen) < args.beam:
        leftovers = [
            x for x in generated
            if x[2] not in chosen_keys
        ]

        random.shuffle(leftovers)

        for x in leftovers:
            add(x)

            if len(chosen) >= args.beam:
                break

    frontier = [
        (Q,fp,k,f"g{gen}")
        for Q,fp,k,score in chosen[:args.beam]
    ]

    # Persist frontier.
    gdir = out / f"frontier-g{gen:02d}"
    gdir.mkdir(exist_ok=True)

    for i,(Q,fp,k,src) in enumerate(frontier):
        save_form(
            Q,
            gdir/f"{i:03d}-r{fp[0]}-n{fp[1]}.txt",
            generation=gen,
            roots=fp[0],
            n4=fp[1],
        )

    root_counts = {}

    for _,fp,_,_ in frontier:
        root_counts[fp[0]] = root_counts.get(fp[0],0)+1

    best = min(
        (fp for _,fp,_,_ in frontier),
        key=lambda f:target_score(f)
    )

    floor = min(
        (n4 for n4,_,_,_ in rootless_hits),
        default=-1
    )

    print(
        f"DIVERSITY|generation={gen}"
        f"|generated={len(generated)}"
        f"|fingerprints={len(representatives)}"
        f"|frontier={len(frontier)}"
        f"|seen={len(seen)}"
        f"|best={best}"
        f"|root_counts={root_counts}"
        f"|rootless_floor={floor}",
        flush=True
    )

print("DIVERSITY|done",flush=True)
