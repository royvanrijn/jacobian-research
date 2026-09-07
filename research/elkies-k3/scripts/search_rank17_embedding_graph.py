from pathlib import Path
import argparse
import numpy as np
from sage.all import Matrix, RDF

parser = argparse.ArgumentParser()
parser.add_argument("gram")
parser.add_argument("--tol", type=float, default=0.08)
parser.add_argument("--limit", type=int, default=10000)
parser.add_argument("--max-shells", type=int, default=50)
args = parser.parse_args()

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

# ------------------------------------------------------------
# Target rank-17 Gram.
# ------------------------------------------------------------

Q = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

assert Q.shape == (17,17)
assert np.all(np.diag(Q) == 4)

print("target offdiag values =",
      sorted(set(Q[np.triu_indices(17,1)].tolist())))

# ------------------------------------------------------------
# Ambient Gram + LLL.
# ------------------------------------------------------------

HA_np = np.loadtxt(args.gram, dtype=float)
r = HA_np.shape[0]

assert HA_np.shape == (r,r)

HA = Matrix(RDF, HA_np.tolist())

print(f"ambient_rank={r}")
print("LLL start", flush=True)

U = HA.LLL_gram()

print("LLL done", flush=True)

R = U.transpose() * HA * U

R_np = np.array(R, dtype=float)
U_np = np.array(U, dtype=np.int64)

# ------------------------------------------------------------
# Cheap short-vector pool:
#   e_i
#   2e_i
#   e_i +/- e_j
#   e_i +/- e_j +/- e_k (bounded)
# ------------------------------------------------------------

cand = {}

def canonical(z):
    z = np.asarray(z, dtype=np.int64)

    nz = np.flatnonzero(z)
    if not len(nz):
        return None

    if z[nz[0]] < 0:
        z = -z

    return tuple(map(int,z))

def add(z):
    key = canonical(z)

    if key is None or key in cand:
        return

    z = np.array(key, dtype=np.int64)
    norm = float(z @ R_np @ z)

    cand[key] = norm


for i in range(r):

    z=np.zeros(r,dtype=np.int64)
    z[i]=1
    add(z)

    z=np.zeros(r,dtype=np.int64)
    z[i]=2
    add(z)


for i in range(r):
    for j in range(i+1,r):
        for s in (-1,1):

            z=np.zeros(r,dtype=np.int64)
            z[i]=1
            z[j]=s

            add(z)


# triples: only +/-1 combinations
for i in range(r):
    for j in range(i+1,r):
        for k in range(j+1,r):

            for sj in (-1,1):
                for sk in (-1,1):

                    z=np.zeros(r,dtype=np.int64)
                    z[i]=1
                    z[j]=sj
                    z[k]=sk

                    add(z)


items=sorted(
    cand.items(),
    key=lambda x:x[1]
)[:args.limit]

Z=np.array(
    [np.array(k,dtype=np.int64) for k,n in items]
)

norms=np.array([n for k,n in items])

# Original E29 coordinates.
V=Z @ U_np.T

print("candidate vectors =",len(V))
print(
    "norm range =",
    float(norms.min()),
    float(norms.max())
)

# Complete ambient pairing matrix.
P = (V @ HA_np) @ V.T

# ------------------------------------------------------------
# Target basis ordering.
#
# Start with target vertex having the most restrictive pairing
# profile.
# ------------------------------------------------------------

profiles=[]

for i in range(17):

    vals=Q[i].copy()
    vals=np.delete(vals,i)

    counts=tuple(
        (x,int(np.sum(vals==x)))
        for x in sorted(set(vals))
    )

    profiles.append((counts,i))

# Score restrictive rows by rarity / nonzero count.
order=sorted(
    range(17),
    key=lambda i:(
        -np.sum(Q[i]!=0),
        tuple(Q[i])
    )
)

print("target search order =",order)

# ------------------------------------------------------------
# Shells.
#
# All 17 target basis vectors have norm 4, so an embedded
# copy must use vectors from essentially one ambient norm shell.
# ------------------------------------------------------------

shells=[]

for center in norms:

    mask=np.abs(norms/center-1.0) <= args.tol

    inds=np.flatnonzero(mask)

    if len(inds) < 17:
        continue

    # Deduplicate near-identical shells.
    key=tuple(inds.tolist())

    if shells and shells[-1][1]==key:
        continue

    shells.append((center,key))


# Prefer shells with lots of available vectors.
shells.sort(
    key=lambda x:-len(x[1])
)

shells=shells[:args.max_shells]

print("shells =",len(shells))

best_depth=0
best_assignment=None
best_h=None
nodes=0


def pairing_ok(actual, target, h):

    wanted=h*target

    # Scale-aware absolute tolerance.
    err=abs(actual-wanted)

    return err <= args.tol * max(
        h,
        abs(wanted),
        1.0
    )


# ------------------------------------------------------------
# Backtracking embedding.
# ------------------------------------------------------------

for si,(center,key) in enumerate(shells):

    C=np.array(key,dtype=np.int64)

    h=center/4.0

    print(
        f"SHELL|i={si}"
        f"|norm={center:.12g}"
        f"|h={h:.12g}"
        f"|size={len(C)}",
        flush=True
    )

    assignment={}
    used=set()

    # Pairings restricted to shell.
    PS=P[np.ix_(C,C)]

    # Candidate lists initially all shell vertices.
    poss={
        ti:list(range(len(C)))
        for ti in range(17)
    }

    # Enforce norm more tightly.
    for ti in range(17):
        poss[ti]=[
            x for x in poss[ti]
            if abs(
                P[C[x],C[x]]/(4*h)-1
            ) <= args.tol
        ]


    def search(depth):
        global best_depth,best_assignment,best_h,nodes

        nodes += 1

        if depth > best_depth:
            best_depth=depth

            best_assignment=dict(assignment)
            best_h=h

            print(
                f"DEPTH|shell={si}"
                f"|depth={depth}"
                f"|h={h:.12g}"
                f"|nodes={nodes}"
                f"|assignment={best_assignment}",
                flush=True
            )

        if depth==17:

            print()
            print("FOUND FULL EMBEDDING")
            print("h =",h)

            rows=[
                V[C[assignment[i]]]
                for i in range(17)
            ]

            A=np.array(rows,dtype=np.int64)

            G=A @ HA_np @ A.T
            residual=G-h*Q

            rel=np.linalg.norm(residual)/(
                h*np.linalg.norm(Q)
            )

            print("relative residual =",rel)
            print("A =")
            print(A)

            np.savetxt(
                BASE/"results/rank17-E29-embedding-A.txt",
                A,
                fmt="%d"
            )

            np.savetxt(
                BASE/"results/rank17-E29-embedding-residual.txt",
                residual,
                fmt="%.17g"
            )

            raise SystemExit(0)

        # Choose next target vertex by order.
        ti=order[depth]

        candidates=[]

        for x in poss[ti]:

            if x in used:
                continue

            ok=True

            for tj,y in assignment.items():

                actual=PS[x,y]
                target=Q[ti,tj]

                if not pairing_ok(
                    actual,target,h
                ):
                    ok=False
                    break

            if ok:
                candidates.append(x)

        # Fail-fast.
        if not candidates:
            return

        # Prefer candidates that leave many compatible neighbors.
        scored=[]

        for x in candidates:

            compatibility=0

            for future in order[depth+1:]:

                for y in poss[future]:

                    if y in used or y==x:
                        continue

                    if pairing_ok(
                        PS[x,y],
                        Q[ti,future],
                        h
                    ):
                        compatibility+=1

            scored.append(
                (-compatibility,x)
            )

        scored.sort()

        for _,x in scored:

            assignment[ti]=x
            used.add(x)

            search(depth+1)

            used.remove(x)
            del assignment[ti]


    search(0)


print()
print("NO FULL EMBEDDING")

print("best_depth =",best_depth)
print("best_h =",best_h)
print("best_assignment =",best_assignment)
print("nodes =",nodes)

if best_assignment:

    print()
    print("BEST PARTIAL AMBIENT VECTORS")

    for ti,x in sorted(best_assignment.items()):
        # Need recover shell corresponding to best_h.
        print(
            "target",ti,
            "candidate-index",x
        )
