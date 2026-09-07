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
                H[i,i] // 2 if i == j
                else H[i,j]
            )

    Q = QuadraticForm(ZZ,n,coeffs)
    assert Q.Hessian_matrix() == H
    return Q


def key(Q):
    H = Q.Hessian_matrix()
    return hashlib.sha256(
        repr(tuple(tuple(map(int,r)) for r in H.rows())).encode()
    ).hexdigest()


def fp(Q):
    S = Q.short_vector_list_up_to_length(3, True)
    return len(S[1]),len(S[2])


def save(Q,path,**meta):
    H = Q.Hessian_matrix()

    lines = [
        f"# {k} = {v}"
        for k,v in meta.items()
    ]

    lines += [
        " ".join(map(str,row))
        for row in H.rows()
    ]

    Path(path).write_text("\n".join(lines)+"\n")


ap = argparse.ArgumentParser()

ap.add_argument("--gram-dir",required=True)
ap.add_argument("--out",required=True)
ap.add_argument("--rounds",type=int,default=20)

# Number of first-step neighbors sampled from each rootless anchor.
ap.add_argument("--bridges-per-anchor",type=int,default=100)

# How many good rooted intermediates to retain globally.
ap.add_argument("--bridge-beam",type=int,default=40)

# Number of second-step neighbors per rooted intermediate.
ap.add_argument("--second-per-bridge",type=int,default=150)

ap.add_argument("--anchors",type=int,default=20)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--primes",default="5,7,11,13")

args = ap.parse_args()

random.seed(args.seed)
set_random_seed(args.seed)

out = Path(args.out)
out.mkdir(parents=True,exist_ok=True)

primes = [ZZ(x) for x in args.primes.split(",")]

# ------------------------------------------------------------
# Load rootless anchors.
# ------------------------------------------------------------

anchors = []
seen = set()

for p in Path(args.gram_dir).glob("*.txt"):
    try:
        Q = form_from_hessian(read_hessian(p))
        f = fp(Q)
    except Exception:
        continue

    if f[0] != 0:
        continue

    k = key(Q)

    if k in seen:
        continue

    seen.add(k)
    anchors.append((f[1],Q,k,str(p)))

anchors.sort(key=lambda x:abs(x[0]-1311))
anchors = anchors[:args.anchors]

print(
    f"TWOSTEP|start|anchors={len(anchors)}"
    f"|floor={min((x[0] for x in anchors),default=-1)}",
    flush=True
)

global_seen = set(x[2] for x in anchors)
best_floor = min(x[0] for x in anchors)

for round_no in range(1,args.rounds+1):

    # ========================================================
    # STEP 1:
    # rootless -> rooted bridge
    # ========================================================

    bridges = []

    for anchor_n4,Q,anchor_key,src in anchors:

        for sample in range(args.bridges_per_anchor):

            p = random.choice(primes)

            try:
                v = Q.find_primitive_p_divisible_vector__random(p)
                N = Q.find_p_neighbor_from_vec(p,v)

                if N.Hessian_matrix().det() != 948:
                    continue

                k = key(N)

                if k in global_seen:
                    continue

                global_seen.add(k)

                f = fp(N)

            except Exception:
                continue

            r,n4 = f

            # We specifically want rooted bridges.
            if r not in (1,2):
                continue

            # Keep a wide shell corridor.
            if not (1280 <= n4 <= 1340):
                continue

            # Ranking deliberately allows r=2 and large shell excursions.
            bridge_score = (
                0 if r == 1 else 3,
                abs(n4-1311)
            )

            bridges.append(
                (bridge_score,f,N,k,anchor_n4,p)
            )

    bridges.sort(
        key=lambda x:x[0]
    )

    # Diversity by fingerprint.
    retained = []
    per_fp = {}

    for item in bridges:
        f = item[1]

        if per_fp.get(f,0) >= 3:
            continue

        per_fp[f] = per_fp.get(f,0)+1
        retained.append(item)

        if len(retained) >= args.bridge_beam:
            break

    print(
        f"TWOSTEP|round={round_no}"
        f"|bridges={len(bridges)}"
        f"|retained={len(retained)}",
        flush=True
    )

    # ========================================================
    # STEP 2:
    # rooted bridge -> rootless endpoint
    # ========================================================

    new_rootless = []

    for idx,(bridge_score,bfp,Q,k,anchor_n4,p1) in enumerate(retained):

        for sample in range(args.second_per_bridge):

            p2 = random.choice(primes)

            try:
                v = Q.find_primitive_p_divisible_vector__random(p2)
                N = Q.find_p_neighbor_from_vec(p2,v)

                if N.Hessian_matrix().det() != 948:
                    continue

                nk = key(N)

                if nk in global_seen:
                    continue

                global_seen.add(nk)

                f = fp(N)

            except Exception:
                continue

            r,n4 = f

            if r != 0:
                continue

            save(
                N,
                out / (
                    f"rootless-r{round_no}"
                    f"-from-r{bfp[0]}n{bfp[1]}"
                    f"-p{p1}-{p2}"
                    f"-n{n4}-{nk[:10]}.txt"
                ),
                round=round_no,
                roots=0,
                n4=n4,
                bridge_roots=bfp[0],
                bridge_n4=bfp[1],
                prime1=p1,
                prime2=p2,
                anchor_n4=anchor_n4,
            )

            print(
                f"TWOSTEP|rootless"
                f"|round={round_no}"
                f"|n4={n4}"
                f"|bridge={bfp}"
                f"|p1={p1}|p2={p2}"
                f"|anchor={anchor_n4}",
                flush=True
            )

            new_rootless.append(
                (n4,N,nk)
            )

            if n4 < best_floor:
                best_floor = n4

                print(
                    f"TWOSTEP|NEW_FLOOR"
                    f"|round={round_no}"
                    f"|n4={n4}"
                    f"|bridge={bfp}"
                    f"|p1={p1}|p2={p2}",
                    flush=True
                )

            if n4 == 1311:

                jackpot = out / (
                    f"JACKPOT-round{round_no}"
                    f"-p{p1}-{p2}-{nk[:10]}.txt"
                )

                save(
                    N,
                    jackpot,
                    round=round_no,
                    roots=0,
                    n4=1311,
                    bridge=bfp,
                    prime1=p1,
                    prime2=p2,
                )

                print(
                    f"TWOSTEP|JACKPOT"
                    f"|round={round_no}"
                    f"|file={jackpot}",
                    flush=True
                )

                raise SystemExit(0)

    # ========================================================
    # Refresh rootless anchors:
    # best old anchors + all interesting new rootless endpoints.
    # ========================================================

    pool = [
        (n4,Q,k,"old")
        for n4,Q,k,_ in anchors
    ]

    pool += [
        (n4,Q,k,"new")
        for n4,Q,k in new_rootless
    ]

    pool.sort(
        key=lambda x:abs(x[0]-1311)
    )

    anchors = []
    used = set()

    for item in pool:
        n4,Q,k,src = item

        if k in used:
            continue

        used.add(k)
        anchors.append(item)

        if len(anchors) >= args.anchors:
            break

    print(
        f"TWOSTEP|round_done={round_no}"
        f"|new_rootless={len(new_rootless)}"
        f"|anchors={len(anchors)}"
        f"|floor={best_floor}"
        f"|seen={len(global_seen)}",
        flush=True
    )

print(
    f"TWOSTEP|done|floor={best_floor}",
    flush=True
)
