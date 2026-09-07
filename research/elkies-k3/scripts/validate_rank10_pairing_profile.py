#!/usr/bin/env python3
"""Infer the Q10/Coxeter height-pairing profile from rational group-law denominators.

The target numerical generator fingerprint is

    (0,0,2,2,2,1,1,1,1).

Only the three +2 entries were enforced directly in the first rank-10 solve.
A norm-4 point already in the rational Coxeter rank-9 span with those same
+2 anchors must instead have the last four pairings equal to a permutation of

    (0,1,1,2).

This script distinguishes those possibilities without tangent-dimension
heuristics.  For two height-4 sections P,Q, use the group-law operation with
smaller height:

    <Q,P> >= 0 : R = Q-P, slope numerator y_Q+y_P
    <Q,P> <  0 : R = Q+P, slope numerator y_Q-y_P.

With no reducible fibers,

    h(R) = 8 - 2*abs(<Q,P>) = 4 + 2(R.O).

Writing x(R)=N/D^2 in lowest terms therefore gives

    deg D = 2-abs(<Q,P>).

Before cancellation the denominator is (x_Q-x_P)^2.  Thus, when
x_Q-x_P is quartic, the expected polynomial gcd degree between the raw
numerator N and raw denominator is

    pairing 0   -> gcd degree 4
    pairing +/-1-> gcd degree 6
    pairing +/-2-> gcd degree 8.

We estimate the gcd degree from the bottom singular-value gap of the Sylvester
matrix.  As an internal calibration, the script first measures the same
quantity on all pairs among the 45 explicit Coxeter sections; their exact
pairings are known from G=2(I+J).

This is a numerical diagnostic, not an exact proof.  Stability of the inferred
degrees and a clean separation of the last four entries as 1,1,1,1 are the
important signals.
"""
from __future__ import annotations

from collections import Counter, defaultdict
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
    coords = [np.eye(9, dtype=int)[i] for i in range(9)]
    for r, (i, j) in enumerate(PAIRS):
        m = s[r]
        xd = pmul(m, m) - x[i] - x[j]
        yd = -y[i] - pmul(m, xd - x[i])
        xs.append(xd)
        ys.append(yd)
        labels.append(f"D{i}{j}")
        v = np.zeros(9, dtype=int)
        v[i] = 1
        v[j] = -1
        coords.append(v)
    return np.asarray(xs), np.asarray(ys), labels, np.asarray(coords), A, B


def trim(a, rel=1e-11):
    a = np.asarray(a, dtype=float).copy()
    scale = float(np.max(np.abs(a))) if len(a) else 0.0
    if scale == 0.0:
        return np.zeros(1)
    nz = np.flatnonzero(np.abs(a) > rel * scale)
    if not len(nz):
        return np.zeros(1)
    return a[: int(nz[-1]) + 1]


def sylvester(f_low, g_low):
    f = trim(f_low)[::-1]
    g = trim(g_low)[::-1]
    m = len(f) - 1
    n = len(g) - 1
    S = np.zeros((m + n, m + n), dtype=float)
    for r in range(n):
        S[r, r : r + m + 1] = f
    for r in range(m):
        S[n + r, r : r + n + 1] = g
    return S, m, n


def gcd_spectrum(num, den):
    nf = trim(num)
    df = trim(den)
    ns = max(float(np.max(np.abs(nf))), 1e-300)
    ds = max(float(np.max(np.abs(df))), 1e-300)
    S, m, n = sylvester(nf / ns, df / ds)
    sv = np.linalg.svd(S, compute_uv=False)
    rel = sv / max(float(sv[0]), 1e-300)
    maxd = min(m, n)
    gaps = []
    for d in range(1, maxd + 1):
        if d >= len(sv):
            break
        hi = float(sv[-d - 1])
        lo = float(sv[-d])
        gaps.append((hi / max(lo, 1e-300), d, hi, lo))
    if gaps:
        gap, d, hi, lo = max(gaps)
    else:
        gap, d, hi, lo = float("nan"), 0, float("nan"), float("nan")
    return {
        "gcd_degree": d,
        "gap": gap,
        "smax": float(sv[0]),
        "smin": float(sv[-1]),
        "tail_rel": rel[-min(10, len(rel)) :],
        "num_degree": len(nf) - 1,
        "den_degree": len(df) - 1,
    }


def addition_gcd_profile(xa, ya, xb, yb, pairing):
    # For p>=0 use A-B, whose chord is through A and -B: numerator yA+yB.
    # For p<0 use A+B: numerator yA-yB.
    dx = xa - xb
    lhs = ya + yb if pairing >= 0 else ya - yb
    den = pmul(dx, dx)
    num = pmul(lhs, lhs) - pmul(xa + xb, den)
    return gcd_spectrum(num, den)


def expected_gcd_degree(pairing):
    return 4 + 2 * abs(int(pairing))


def infer_pairing_from_gcd(d):
    # nearest expected class among |p|=0,1,2
    choices = [(abs(d - (4 + 2 * p)), p) for p in (0, 1, 2)]
    return min(choices)[1]


def fake_rank9_completions():
    # Given b0=b1=0, b2=b3=b4=2 and norm 4 in G=2(I+J).
    out = []
    for rest in itertools.product(range(-2, 3), repeat=4):
        b = (0, 0, 2, 2, 2) + tuple(rest)
        norm = 0.5 * (sum(x * x for x in b) - (sum(b) ** 2) / 10.0)
        if abs(norm - 4.0) < 1e-12:
            out.append(b)
    return out


ap = argparse.ArgumentParser()
ap.add_argument(
    "--refined",
    type=Path,
    default=BASE / "results/rank10-winning-refine-v1",
)
ap.add_argument("--v-pairings", default="0,0,2,2,2,1,1,1,1")
args = ap.parse_args()

refined = args.refined.resolve()
pv = tuple(map(int, args.v_pairings.split(",")))
if len(pv) != 9:
    raise SystemExit("--v-pairings must contain nine integers")

slopes = np.loadtxt(refined / "slopes.txt", dtype=float).reshape(36, 3)
xq = np.loadtxt(refined / "x-q.txt", dtype=float).reshape(-1)
yq = np.loadtxt(refined / "y-q.txt", dtype=float).reshape(-1)
xs, ys, labels, coords, A, B = derived45(slopes)
G = 2 * (np.eye(9, dtype=int) + np.ones((9, 9), dtype=int))

print("EXACT RANK-9 FALSE-POSITIVE COMPLETIONS")
fakes = fake_rank9_completions()
print("count =", len(fakes))
for b in fakes:
    print("FAKE", " ".join(f"{x:+d}" for x in b))
print("target", " ".join(f"{x:+d}" for x in pv))
print()

# Internal controls from the 45 known sections.
control = defaultdict(list)
for a in range(45):
    for b in range(a + 1, 45):
        p = int(coords[a] @ G @ coords[b])
        if abs(p) > 2:
            continue
        prof = addition_gcd_profile(xs[a], ys[a], xs[b], ys[b], p)
        control[p].append((prof["gcd_degree"], prof["gap"]))

print("COXETER CONTROL PROFILES")
for p in sorted(control):
    vals = control[p]
    counts = Counter(d for d, _ in vals)
    expected = expected_gcd_degree(p)
    gaps_good = [g for d, g in vals if d == expected and np.isfinite(g)]
    print(
        f"CONTROL|pairing={p:+d}|expected_gcd={expected}|n={len(vals)}"
        f"|degree_counts={dict(sorted(counts.items()))}"
        + (f"|median_good_gap={np.median(gaps_good):.3e}" if gaps_good else "")
    )
print()

print("Q AGAINST COXETER GENERATORS")
matches = 0
last4 = []
for i in range(9):
    p = int(pv[i])
    prof = addition_gcd_profile(xq, yq, xs[i], ys[i], p)
    expected = expected_gcd_degree(p)
    inferred_abs = infer_pairing_from_gcd(prof["gcd_degree"])
    ok = prof["gcd_degree"] == expected
    matches += int(ok)
    if i >= 5:
        last4.append(inferred_abs)
    tail = " ".join(f"{x:.2e}" for x in prof["tail_rel"])
    print(
        f"PAIR|V{i}|target={p:+d}|expected_gcd={expected}"
        f"|observed_gcd={prof['gcd_degree']}|gap={prof['gap']:.3e}"
        f"|inferred_abs_pairing={inferred_abs}|match={int(ok)}"
        f"|numdeg={prof['num_degree']}|dendeg={prof['den_degree']}|tail={tail}"
    )

print()
print("SUMMARY")
print("generator_profile_matches =", f"{matches}/9")
print("last4_inferred_abs_pairings =", " ".join(map(str, last4)))
print("target_last4 = 1 1 1 1")
print("rank9_fake_last4_multiset = 0 1 1 2")
print("target_profile_supported =", last4 == [1, 1, 1, 1] and matches == 9)
