#!/usr/bin/env python3
"""Directly test whether the refined Q_10 is genuinely outside the Coxeter rank-9 lattice.

The rank-10 reconstruction only enforced the nine target |height pairing|=2
relations.  A possible false-positive mode is that Q is actually a norm-4
section already lying in the saturated rank-9 Coxeter lattice, while matching
those nine required edges and also having *additional* +/-2 edges that were not
part of the target fingerprint.

This script avoids tangent-dimension inference and tests the section itself.
For every one of the 45 explicit Coxeter minimal sections (9 V_i and 36
D_ij=V_i-V_j), fit both possible quadratic chord identities

    y_Q + y_P = m(T) (x_Q-x_P)   [pairing +2]
    y_Q - y_P = m(T) (x_Q-x_P)   [pairing -2]

and compare the observed +/-2 edge set with the target fingerprint

    (0,0,2,2,2,1,1,1,1).

If the observed edge set is exactly the intended nine edges, there is no
rank-9 norm-4 pairing completion with the same +/-2 incidence pattern for the
raw Coxeter Gram 2(I+J).  This rules out the most important saturation
false-positive, subject to the usual height-4/no-reducible-fiber assumptions.

The script also reports effective degrees, distance from all 45 known sections,
and a numerical squarefreeness diagnostic for the discriminant polynomial.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import itertools
import numpy as np

BASE = Path(__file__).resolve().parents[1]
PAIRS = [(i, j) for i in range(9) for j in range(i + 1, 9)]
PI = {p: n for n, p in enumerate(PAIRS)}
TRIPLES = [(i, j, k) for i in range(9) for j in range(i + 1, 9) for k in range(j + 1, 9)]

TINC = np.zeros((84, 9), dtype=float)
EINC = np.zeros((36, 9), dtype=float)
for r, t in enumerate(TRIPLES):
    TINC[r, list(t)] = 1.0
for r, (i, j) in enumerate(PAIRS):
    EINC[r, [i, j]] = 1.0
TP = np.linalg.inv(TINC.T @ TINC) @ TINC.T
EP = np.linalg.inv(EINC.T @ EINC) @ EINC.T


def pmul(a, b):
    return np.convolve(a, b)


def p3(a):
    return pmul(pmul(a, a), a)


def padd(a, b, sa=1.0, sb=1.0):
    z = np.zeros(max(len(a), len(b)), dtype=float)
    z[: len(a)] += sa * a
    z[: len(b)] += sb * b
    return z


def pslope(s, i, j):
    return -s[PI[(j, i)]] if i > j else s[PI[(i, j)]]


def reconstruct(s):
    f = np.empty((84, 5), dtype=float)
    for r, (i, j, k) in enumerate(TRIPLES):
        mij, mik, mjk = pslope(s, i, j), pslope(s, i, k), pslope(s, j, k)
        f[r] = pmul(mik, mij) + pmul(mik, mjk) - pmul(mij, mjk)
    x = TP @ f

    q = np.empty((36, 7), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        q[r] = pmul(s[r], x[i] - x[j])
    y = EP @ q

    ae = np.empty((36, 9), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        ae[r] = (
            pmul(s[r], y[i] - y[j])
            - pmul(x[i], x[i])
            - pmul(x[i], x[j])
            - pmul(x[j], x[j])
        )
    A = ae.mean(axis=0)

    bn = np.empty((9, 13), dtype=float)
    for i in range(9):
        bn[i] = pmul(y[i], y[i]) - p3(x[i]) - pmul(A, x[i])
    B = bn.mean(axis=0)
    return x, y, A, B


def derived45(s):
    x, y, A, B = reconstruct(s)
    xs = [x[i].copy() for i in range(9)]
    ys = [y[i].copy() for i in range(9)]
    labels = [f"V{i}" for i in range(9)]
    for r, (i, j) in enumerate(PAIRS):
        m = s[r]
        xd = pmul(m, m) - x[i] - x[j]
        yd = -y[i] - pmul(m, xd - x[i])
        xs.append(xd)
        ys.append(yd)
        labels.append(f"D{i}{j}")
    return np.asarray(xs), np.asarray(ys), labels, A, B


def conv_matrix(dx):
    M = np.zeros((7, 3), dtype=float)
    for j in range(3):
        M[j : j + 5, j] = dx
    return M


def line_fit(xq, yq, xp, yp, sign):
    dx = xq - xp
    lhs = yq + sign * yp
    M = conv_matrix(dx)
    m, *_ = np.linalg.lstsq(M, lhs, rcond=None)
    err = lhs - M @ m
    abs_err = float(np.max(np.abs(err)))
    scale = max(float(np.max(np.abs(lhs))), float(np.max(np.abs(M @ m))), 1e-300)
    rel_err = abs_err / scale
    return abs_err, rel_err, m


def fingerprint45(pv):
    fp = list(map(int, pv))
    fp += [int(pv[i] - pv[j]) for i, j in PAIRS]
    return tuple(fp)


def effective_degree(a, rel=1e-9):
    a = np.asarray(a, dtype=float)
    scale = float(np.max(np.abs(a)))
    if scale == 0.0:
        return -1
    nz = np.flatnonzero(np.abs(a) > rel * scale)
    return int(nz[-1]) if len(nz) else -1


def discriminant_root_separation(A, B):
    core = padd(p3(A), pmul(B, B), 4.0, 27.0)
    deg = effective_degree(core)
    if deg <= 0:
        return deg, float("nan"), float("nan")
    coeff = core[: deg + 1]
    coeff = coeff / max(float(np.max(np.abs(coeff))), 1e-300)
    roots = np.roots(coeff[::-1])
    if len(roots) < 2:
        return deg, float("inf"), float("inf")
    mind = float("inf")
    relmind = float("inf")
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            d = abs(roots[i] - roots[j])
            mind = min(mind, float(d))
            relmind = min(relmind, float(d / max(1.0, abs(roots[i]), abs(roots[j]))))
    return deg, mind, relmind


def rank9_completion_count(target_pv):
    """Count norm-4 vectors in the rational raw-Coxeter span with exactly the
    same +/-2 incidence pattern against all 45 Coxeter nodes.

    For distinct norm-4 vectors in a min-height-4 lattice, pairings lie in
    [-2,2].  The Gram is G=2(I+J), so G^-1=(1/2)(I-J/10).
    """
    target_fp = fingerprint45(target_pv)
    target_edges = tuple((abs(p) == 2, int(np.sign(p))) for p in target_fp)
    count = 0
    examples = []
    # The target edges force b2=b3=b4=2 and b0=b1=0; only enumerate the rest.
    for rest in itertools.product(range(-2, 3), repeat=4):
        b = (0, 0, 2, 2, 2) + tuple(rest)
        fp = fingerprint45(b)
        observed = tuple((abs(p) == 2, int(np.sign(p))) for p in fp)
        if observed != target_edges:
            continue
        norm = 0.5 * (sum(x * x for x in b) - (sum(b) ** 2) / 10.0)
        if abs(norm - 4.0) < 1e-12:
            count += 1
            examples.append(b)
    return count, examples


ap = argparse.ArgumentParser()
ap.add_argument(
    "--refined",
    type=Path,
    default=BASE / "results" / "rank10-winning-refine-v1",
)
ap.add_argument("--v-pairings", default="0,0,2,2,2,1,1,1,1")
ap.add_argument("--abs-tol", type=float, default=1e-8)
ap.add_argument("--rel-tol", type=float, default=1e-8)
args = ap.parse_args()

refined = args.refined.resolve()
pv = tuple(map(int, args.v_pairings.split(",")))
if len(pv) != 9:
    raise SystemExit("--v-pairings must contain nine integers")

required = [refined / "slopes.txt", refined / "x-q.txt", refined / "y-q.txt"]
for p in required:
    if not p.exists():
        raise SystemExit(f"missing required input: {p}")

slopes = np.loadtxt(refined / "slopes.txt", dtype=float).reshape(36, 3)
xq = np.loadtxt(refined / "x-q.txt", dtype=float).reshape(-1)
yq = np.loadtxt(refined / "y-q.txt", dtype=float).reshape(-1)
xs, ys, labels, A, B = derived45(slopes)
target_fp = fingerprint45(pv)
expected = {i: int(p) for i, p in enumerate(target_fp) if abs(int(p)) == 2}

observed = {}
rows = []
for i, label in enumerate(labels):
    plus_abs, plus_rel, _ = line_fit(xq, yq, xs[i], ys[i], +1.0)
    minus_abs, minus_rel, _ = line_fit(xq, yq, xs[i], ys[i], -1.0)
    if (plus_rel, plus_abs) <= (minus_rel, minus_abs):
        best_pairing, best_abs, best_rel = +2, plus_abs, plus_rel
    else:
        best_pairing, best_abs, best_rel = -2, minus_abs, minus_rel
    hit = best_abs <= args.abs_tol and best_rel <= args.rel_tol
    if hit:
        observed[i] = best_pairing
    rows.append((i, label, plus_abs, plus_rel, minus_abs, minus_rel, hit, best_pairing, best_abs, best_rel))

missing = [(i, labels[i], p) for i, p in expected.items() if observed.get(i) != p]
extra = [(i, labels[i], p) for i, p in observed.items() if i not in expected]
wrong_sign = [(i, labels[i], expected[i], observed[i]) for i in expected if i in observed and observed[i] != expected[i]]

print("TARGET +/-2 EDGES")
print(" ".join(f"{labels[i]}:{p:+d}" for i, p in expected.items()))
print()
print("OBSERVED QUADRATIC-LINE EDGES")
for i, label, pa, pr, ma, mr, hit, bp, ba, br in rows:
    marker = "EDGE" if hit else "----"
    exp = expected.get(i)
    print(
        f"{marker}|{label}|best={bp:+d}|abs={ba:.3e}|rel={br:.3e}"
        f"|plus_abs={pa:.3e}|plus_rel={pr:.3e}|minus_abs={ma:.3e}|minus_rel={mr:.3e}"
        + (f"|expected={exp:+d}" if exp is not None else "")
    )

# Distinctness from the 45 explicitly reconstructed sections, including sign.
xd = np.linalg.norm(xs - xq[None, :], axis=1)
yplus = np.linalg.norm(ys - yq[None, :], axis=1)
yminus = np.linalg.norm(ys + yq[None, :], axis=1)
nearest = int(np.argmin(xd))
print()
print("DISTINCTNESS")
print(f"nearest_x={labels[nearest]}|x_distance={xd[nearest]:.6e}|y_same={yplus[nearest]:.6e}|y_opposite={yminus[nearest]:.6e}")
print(f"xQ_degree={effective_degree(xq)}|yQ_degree={effective_degree(yq)}")

# Discriminant squarefreeness heuristic.
ddeg, dsep, drelsep = discriminant_root_separation(A, B)
print()
print("DISCRIMINANT")
print(f"degree={ddeg}|min_root_separation={dsep:.6e}|min_relative_root_separation={drelsep:.6e}")

completion_count, examples = rank9_completion_count(pv)
exact_edges = not missing and not extra and not wrong_sign
print()
print("SUMMARY")
print("expected_edge_count =", len(expected))
print("observed_edge_count =", len(observed))
print("missing_edges =", len(missing))
print("extra_edges =", len(extra))
print("wrong_sign_edges =", len(wrong_sign))
print("exact_target_abs2_pattern =", exact_edges)
print("rank9_norm4_completion_count_with_same_abs2_pattern =", completion_count)
if examples:
    print("rank9_completion_examples =", examples[:12])
print("independence_pattern_pass =", bool(exact_edges and completion_count == 0))

if missing:
    print("MISSING", " ".join(f"{lab}:{p:+d}" for _, lab, p in missing))
if extra:
    print("EXTRA", " ".join(f"{lab}:{p:+d}" for _, lab, p in extra))
if wrong_sign:
    print("WRONG_SIGN", wrong_sign)
