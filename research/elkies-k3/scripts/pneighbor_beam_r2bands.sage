from sage.all import *
from sage.quadratic_forms.genera.genus import Genus
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
            coeffs.append(H[i,i]//2 if i == j else H[i,j])

    Q = QuadraticForm(ZZ, n, coeffs)
    assert Q.Hessian_matrix() == H
    return Q


def key(Q):
    H = Q.Hessian_matrix()
    return hashlib.sha256(
        repr(tuple(tuple(map(int,r)) for r in H.rows())).encode()
    ).hexdigest()


def fp(Q):
    S = Q.short_vector_list_up_to_length(3, True)
    return len(S[1]), len(S[2])


def score(f):
    r,n4 = f

    # Absolute target first.
    if r == 0 and n4 == 1311:
        return -1000000

    # Prefer rootless only when close.
    if r == 0:
        return 4*abs(n4-1311)

    # Keep rooted ridge around 1311 alive.
    if r == 1:
        return 20 + abs(n4-1311)

    return 100 + 20*r + abs(n4-1311)


def save(Q,path,**meta):
    H=Q.Hessian_matrix()
    lines=[f"# {k} = {v}" for k,v in meta.items()]
    lines += [" ".join(map(str,row)) for row in H.rows()]
    Path(path).write_text("\n".join(lines)+"\n")


ap=argparse.ArgumentParser()
ap.add_argument("--gram-dir",required=True)
ap.add_argument("--out",required=True)
ap.add_argument("--generations",type=int,default=20)
ap.add_argument("--beam",type=int,default=30)
ap.add_argument("--samples-per-parent",type=int,default=150)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--primes",default="5,7,11")

args=ap.parse_args()

random.seed(args.seed)
set_random_seed(args.seed)

out=Path(args.out)
out.mkdir(parents=True,exist_ok=True)

primes=[ZZ(x) for x in args.primes.split(",")]

frontier=[]
seen=set()

for p in sorted(Path(args.gram_dir).glob("*.txt")):
    H=read_hessian(p)
    Q=form_from_hessian(H)

    assert H.det()==948

    f=fp(Q)
    assert f==(1,1311), (p,f)

    k=key(Q)
    if k in seen:
        continue
    seen.add(k)

    frontier.append((score(f),f,Q,str(p)))

print("BEAM|start|states=",len(frontier),flush=True)

for gen in range(1,args.generations+1):

    candidates=[]
    rootless=[]

    for parent_idx,(_,parent_fp,Q,parent_src) in enumerate(frontier):

        for sample in range(args.samples_per_parent):
            p=random.choice(primes)

            try:
                v=Q.find_primitive_p_divisible_vector__random(p)
                N=Q.find_p_neighbor_from_vec(p,v)
                H=N.Hessian_matrix()

                if H.det()!=948:
                    continue

                k=key(N)

                if k in seen:
                    continue

                seen.add(k)

                f=fp(N)
                s=score(f)

            except Exception:
                continue

            roots,n4=f

            if roots==0:
                rootless.append((n4,p,N,parent_fp))

                save(
                    N,
                    out/f"rootless-g{gen}-p{p}-n{n4}-{k[:10]}.txt",
                    generation=gen,
                    prime=p,
                    roots=roots,
                    n4=n4,
                    parent=parent_fp,
                )

                print(
                    f"BEAM|rootless|gen={gen}|p={p}|n4={n4}|parent={parent_fp}",
                    flush=True
                )

            if f==(0,1311):
                jackpot=out/f"JACKPOT-g{gen}-p{p}-{k[:10]}.txt"

                save(
                    N,
                    jackpot,
                    generation=gen,
                    prime=p,
                    roots=0,
                    n4=1311,
                    parent=parent_fp,
                )

                print(
                    f"BEAM|JACKPOT|gen={gen}|p={p}|file={jackpot}",
                    flush=True
                )
                raise SystemExit(0)

            # Keep only reasonably interesting candidates.
            if roots <= 2 and abs(n4-1311) <= 40:
                candidates.append((s,f,N,k))

    candidates.sort(
        key=lambda x:(x[0], x[1][0], abs(x[1][1]-1311))
    )

    # --------------------------------------------------------
    # STRATIFIED BEAM
    #
    # Do not let the dense (1,1311) ridge consume the beam.
    # Explicitly preserve:
    #   - one-root states near 1311
    #   - two-root states near 1311
    #   - rootless states
    #   - diverse leftovers
    # --------------------------------------------------------

    one_root = sorted(
        [x for x in candidates if x[1][0] == 1],
        key=lambda x:(abs(x[1][1]-1311), x[0])
    )

    rootless_pool = sorted(
        [x for x in candidates if x[1][0] == 0],
        key=lambda x:abs(x[1][1]-1311)
    )

    # Deliberately preserve different regions of the r=2 corridor.
    r2bands = [
        (1295,1304),
        (1305,1308),
        (1309,1311),
        (1312,1315),
        (1316,1322),
    ]

    two_bands = []

    for lo,hi in r2bands:
        pool = sorted(
            [
                x for x in candidates
                if x[1][0] == 2
                and lo <= x[1][1] <= hi
            ],
            key=lambda x:(
                abs(x[1][1]-1311),
                x[0]
            )
        )

        two_bands.append((lo,hi,pool))

    next_frontier = []
    selected_keys = set()
    per_fp = {}

    def add_from(pool, quota, per_fingerprint=3):
        added = 0

        for item in pool:
            s0,f,Q,k = item

            if k in selected_keys:
                continue

            if per_fp.get(f,0) >= per_fingerprint:
                continue

            per_fp[f] = per_fp.get(f,0) + 1
            selected_keys.add(k)
            next_frontier.append(item)
            added += 1

            if added >= quota:
                break

    # 10 good one-root states.
    add_from(one_root, 10, 2)

    # 8 states from EACH r=2 band = up to 40 r=2 states.
    for lo,hi,pool in two_bands:
        add_from(pool, 8, 2)

    # Keep up to 4 rootless states as boundary anchors.
    add_from(rootless_pool, 4, 2)

    # Fill any remaining beam slots globally.
    for item in candidates:
        if len(next_frontier) >= args.beam:
            break

        s0,f,Q,k = item

        if k in selected_keys:
            continue

        if per_fp.get(f,0) >= 2:
            continue

        per_fp[f] = per_fp.get(f,0) + 1
        selected_keys.add(k)
        next_frontier.append(item)

    if not next_frontier:
        print("BEAM|dead|gen=",gen,flush=True)
        break

    frontier=next_frontier

    print(
        f"BEAM|generation={gen}"
        f"|candidates={len(candidates)}"
        f"|frontier={len(frontier)}"
        f"|seen={len(seen)}"
        f"|best={frontier[0][1]}"
        f"|beam_r0={sum(1 for _,f,_,_ in frontier if f[0] == 0)}"
        f"|beam_r1={sum(1 for _,f,_,_ in frontier if f[0] == 1)}"
        f"|beam_r2={sum(1 for _,f,_,_ in frontier if f[0] == 2)}"
        f"|rootless={sorted(set(n for n,_,_,_ in rootless))}",
        flush=True
    )

    # Persist frontier.
    gdir=out/f"frontier-g{gen:02d}"
    gdir.mkdir(exist_ok=True)

    for i,(s,f,Q,k) in enumerate(frontier):
        save(
            Q,
            gdir/f"{i:03d}-r{f[0]}-n{f[1]}.txt",
            generation=gen,
            score=s,
            roots=f[0],
            n4=f[1],
        )

print("BEAM|done",flush=True)
