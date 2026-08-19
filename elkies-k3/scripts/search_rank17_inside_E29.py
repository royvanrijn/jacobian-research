from pathlib import Path
import numpy as np
import random
import math

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/rank29"

HG = np.loadtxt(
    J/"short_vector_basis_gram.txt",
    dtype=float
)

HS = np.loadtxt(
    D/"E29-height-gram.txt",
    dtype=float
)

assert HG.shape == (17,17)
assert HS.shape == (29,29)

# Generic shape, made scale-invariant.
eigG = np.linalg.eigvalsh(HG)
detG = np.linalg.det(HG)

# Normalize determinant to 1 so scale disappears.
GN = HG / (detG ** (1.0/17.0))

target_eigs = np.linalg.eigvalsh(GN)

def score_subset(indices):
    M = HS[np.ix_(indices,indices)]

    sign,logdet=np.linalg.slogdet(M)
    if sign <= 0:
        return 1e100,None

    scale=np.exp(logdet/17.0)
    MN=M/scale

    eig=np.linalg.eigvalsh(MN)

    # Compare logarithmic spectra; basis-dependent Gram entries
    # themselves are not meaningful under GL(17,Z).
    s=np.mean(
        (np.log(eig)-np.log(target_eigs))**2
    )

    return float(s),(scale,eig)

rng=random.Random(29017)

current=sorted(rng.sample(range(29),17))
current_score,_=score_subset(current)

best=(current_score,list(current))

print("initial",best,flush=True)

T=0.1

for step in range(500000):
    pos=rng.randrange(17)

    outside=[
        i for i in range(29)
        if i not in current
    ]

    newv=rng.choice(outside)

    proposal=list(current)
    proposal[pos]=newv
    proposal=sorted(proposal)

    score,_=score_subset(proposal)

    delta=score-current_score

    if delta <= 0 or rng.random() < math.exp(-delta/max(T,1e-12)):
        current=proposal
        current_score=score

    if score < best[0]:
        best=(score,list(proposal))

        print(
            f"BEST|step={step}"
            f"|score={score:.12g}"
            f"|indices={best[1]}",
            flush=True
        )

    T *= 0.99998

    if T < 1e-5:
        T=0.05

print()
print("FINAL BEST")
print(best)
