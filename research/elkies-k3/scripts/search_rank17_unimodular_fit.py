from pathlib import Path

import numpy as np
import random
import math

BASE = Path(__file__).resolve().parents[1]

J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/rank29"

HG = np.loadtxt(J / "short_vector_basis_gram.txt", dtype=float)
HS = np.loadtxt(D / "E29-height-gram.txt", dtype=float)

assert HG.shape == (17, 17)
assert HS.shape == (29, 29)

SUBSET = [
    0, 1, 4, 5, 6, 8, 9, 11, 15,
    16, 17, 19, 22, 23, 24, 25, 26
]

M = HS[np.ix_(SUBSET, SUBSET)]


def fit_score(U):
    X = U.T @ M @ U

    # Best scalar h in Frobenius norm:
    # minimize ||X - h HG||_F^2
    h = np.sum(X * HG) / np.sum(HG * HG)

    if h <= 0:
        return 1e100, None

    R = X - h * HG

    rel_frob = np.linalg.norm(R, "fro") / (
        h * np.linalg.norm(HG, "fro")
    )

    max_rel = np.max(np.abs(R)) / (
        h * np.max(np.abs(HG))
    )

    sign_x, logdet_x = np.linalg.slogdet(X)
    sign_g, logdet_g = np.linalg.slogdet(HG)

    if sign_x <= 0 or sign_g <= 0:
        det_h = float("nan")
    else:
        det_h = math.exp((logdet_x - logdet_g) / 17.0)

    return rel_frob, {
        "h": h,
        "rel_frob": rel_frob,
        "max_rel": max_rel,
        "det_h": det_h,
        "X": X,
        "R": R,
    }


rng = random.Random(170029)

U = np.eye(17, dtype=np.int64)

current_score, current_info = fit_score(U)
best_score = current_score
best_U = U.copy()
best_info = current_info

print(
    "INITIAL"
    f"|score={current_score:.12g}"
    f"|h={current_info['h']:.12g}"
    f"|det_h={current_info['det_h']:.12g}"
    f"|max_rel={current_info['max_rel']:.12g}",
    flush=True,
)

T = 0.02

for step in range(2_000_000):
    P = U.copy()

    move = rng.randrange(100)

    if move < 10:
        # sign flip
        i = rng.randrange(17)
        P[:, i] *= -1

    elif move < 20:
        # swap columns
        i, j = rng.sample(range(17), 2)
        P[:, [i, j]] = P[:, [j, i]]

    else:
        # elementary unimodular shear
        i, j = rng.sample(range(17), 2)

        # Mostly ±1, occasionally ±2/±3.
        r = rng.random()
        if r < 0.80:
            k = rng.choice([-1, 1])
        elif r < 0.95:
            k = rng.choice([-2, 2])
        else:
            k = rng.choice([-3, 3])

        P[:, i] += k * P[:, j]

    score, info = fit_score(P)

    delta = score - current_score

    if delta <= 0 or rng.random() < math.exp(
        -delta / max(T, 1e-12)
    ):
        U = P
        current_score = score
        current_info = info

    if score < best_score:
        best_score = score
        best_U = P.copy()
        best_info = info

        print(
            f"BEST|step={step}"
            f"|score={score:.12g}"
            f"|h={info['h']:.12g}"
            f"|det_h={info['det_h']:.12g}"
            f"|h_ratio={info['det_h']/info['h']:.12g}"
            f"|max_rel={info['max_rel']:.12g}"
            f"|maxU={np.max(np.abs(best_U))}",
            flush=True,
        )

    T *= 0.999995

    if T < 1e-5:
        T = 0.01


print()
print("FINAL BEST")
print("score =", best_score)
print("h =", best_info["h"])
print("det_h =", best_info["det_h"])
print("h_ratio =", best_info["det_h"] / best_info["h"])
print("max_rel =", best_info["max_rel"])
print()
print("U =")
print(best_U)
print()
print("residual =")
print(best_info["R"])

np.savetxt(
    BASE / "results/rank17-best-U.txt",
    best_U,
    fmt="%d",
)

np.savetxt(
    BASE / "results/rank17-best-residual.txt",
    best_info["R"],
    fmt="%.17g",
)

