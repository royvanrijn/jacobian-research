#!/usr/bin/env python3
"""Calibrate the numerical tangent-dimension drop from Coxeter rank 9 to rank 10.

The equilibrated explicit rank-10 Jacobian has shown a stable strongest bottom
spectral gap at d=15, while geometric K3 moduli counting predicts dimension 8
for a generic MW-rank-10 elliptic K3.  This script separates *parameter-space*
tangent redundancy from physical surface moduli by making two apples-to-apples
comparisons at the saved refined rank-10 solution:

1. Compute the Coxeter-only (rank-9) tangent spectrum at the *same refined base
   point* and compare its strongest bottom gap with the explicit rank-10 one.
   A stable 16 -> 15 drop means the added section contributes exactly one new
   independent condition even if the slope parameterization has extra tangent
   redundancy.

2. Map each numerical tangent space through the differential of the actual
   Weierstrass coefficients (A_0..A_8, B_0..B_12).  The rank of that image is a
   much closer proxy for physical surface-moduli dimension.  The ideal pattern
   is parameter dimensions 16 -> 15 but surface-image dimensions 9 -> 8.

All rank decisions are based on equilibrated singular spectra and reported over
multiple finite-difference steps.  Stable gaps across steps matter more than any
single hard threshold.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import json
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
    rx = f - TINC @ x

    q = np.empty((36, 7), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        q[r] = pmul(s[r], x[i] - x[j])
    y = EP @ q
    ry = q - EINC @ y

    ae = np.empty((36, 9), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        ae[r] = (
            pmul(s[r], y[i] - y[j])
            - pmul(x[i], x[i])
            - pmul(x[i], x[j])
            - pmul(x[j], x[j])
        )
    A = ae.mean(axis=0)
    rA = ae - A

    bn = np.empty((9, 13), dtype=float)
    for i in range(9):
        bn[i] = pmul(y[i], y[i]) - p3(x[i]) - pmul(A, x[i])
    B = bn.mean(axis=0)
    rB = bn - B
    return x, y, A, B, (rx, ry, rA, rB)


def derived45(s):
    x, y, A, B, blocks = reconstruct(s)
    xs = [x[i].copy() for i in range(9)]
    ys = [y[i].copy() for i in range(9)]
    for r, (i, j) in enumerate(PAIRS):
        m = s[r]
        xd = pmul(m, m) - x[i] - x[j]
        yd = -y[i] - pmul(m, xd - x[i])
        xs.append(xd)
        ys.append(yd)
    return np.asarray(xs), np.asarray(ys), A, B, blocks


def fixed_from_meta(meta):
    sign = int(meta["gauge_sign"])
    sp = tuple(sorted(map(int, meta.get("second_pair", [0, 2]))))
    return {
        3 * PI[(0, 1)]: float(sign),
        3 * PI[(0, 1)] + 1: 0.0,
        3 * PI[(0, 1)] + 2: 1.0,
        3 * PI[sp] + 1: 1.0,
    }


def freepos(fixed):
    return [i for i in range(108) if i not in fixed]


def unpack(free, fixed):
    flat = np.empty(108, dtype=float)
    pos = freepos(fixed)
    flat[pos] = free
    for i, v in fixed.items():
        flat[i] = v
    return flat.reshape(36, 3)


def base_raw(free, fixed):
    *_, blocks = reconstruct(unpack(free, fixed))
    return np.concatenate(
        [
            blocks[0].ravel(),
            blocks[1].ravel(),
            0.5 * blocks[2].ravel(),
            0.25 * blocks[3].ravel(),
        ]
    )


def surface_observable_from_base(free, fixed):
    _, _, A, B, _ = reconstruct(unpack(free, fixed))
    return np.concatenate([A, B])


def fingerprint45(pv):
    fp = list(map(int, pv))
    fp += [int(pv[i] - pv[j]) for i, j in PAIRS]
    return tuple(fp)


def anchors_from_pv(pv):
    return [(i, int(p)) for i, p in enumerate(fingerprint45(pv)) if abs(int(p)) == 2]


def choose_primary(anchors):
    vertices = [a for a in anchors if a[0] < 9]
    return vertices[0] if vertices else anchors[0]


def conv_matrix(dx):
    M = np.zeros((7, 3), dtype=float)
    for j in range(3):
        M[j : j + 5, j] = dx
    return M


def explicit_state(saved_state, fixed, anchors):
    """Expand [104 base, 5 xQ, 3 primary slope] to all anchor slopes."""
    base = saved_state[:104]
    xq = saved_state[104:109]
    m0 = saved_state[109:112]
    s = unpack(base, fixed)
    xs, ys, A, B, _ = derived45(s)
    primary = choose_primary(anchors)
    ordered = [primary] + [a for a in anchors if a != primary]
    ai, p = ordered[0]
    yq = pmul(m0, xq - xs[ai]) - (p / 2.0) * ys[ai]
    ms = [m0]
    for aj, pj in ordered[1:]:
        lhs = yq + (pj / 2.0) * ys[aj]
        M = conv_matrix(xq - xs[aj])
        m, *_ = np.linalg.lstsq(M, lhs, rcond=None)
        ms.append(m)
    return np.concatenate([base, xq, np.asarray(ms).reshape(-1)])


def rank10_raw_explicit(z, fixed, anchors):
    base = z[:104]
    xq = z[104:109]
    ms = z[109:].reshape(len(anchors), 3)
    s = unpack(base, fixed)
    xs, ys, A, B, blocks = derived45(s)
    primary = choose_primary(anchors)
    ordered = [primary] + [a for a in anchors if a != primary]
    ai, p = ordered[0]
    yq = pmul(ms[0], xq - xs[ai]) - (p / 2.0) * ys[ai]

    pieces = [
        blocks[0].ravel(),
        blocks[1].ravel(),
        0.5 * blocks[2].ravel(),
        0.25 * blocks[3].ravel(),
        0.5 * (pmul(yq, yq) - p3(xq) - pmul(A, xq) - B),
    ]
    for k, (aj, pj) in enumerate(ordered[1:], start=1):
        pieces.append(yq + (pj / 2.0) * ys[aj] - pmul(ms[k], xq - xs[aj]))
    return np.concatenate(pieces)


def central_jacobian(fun, z, step):
    f0 = fun(z)
    J = np.empty((len(f0), len(z)), dtype=float)
    for c in range(len(z)):
        h = step * max(1.0, abs(float(z[c])))
        a = z.copy()
        b = z.copy()
        a[c] += h
        b[c] -= h
        J[:, c] = (fun(a) - fun(b)) / (2.0 * h)
    return J


def equilibrate(J, iterations=8, floor=1e-300):
    """Ruiz row/column equilibration: A = diag(rs) J diag(cs)."""
    A = np.asarray(J, dtype=float).copy()
    rs = np.ones(A.shape[0], dtype=float)
    cs = np.ones(A.shape[1], dtype=float)
    for _ in range(iterations):
        rn = np.linalg.norm(A, axis=1)
        rf = 1.0 / np.sqrt(np.maximum(rn, floor))
        A *= rf[:, None]
        rs *= rf
        cn = np.linalg.norm(A, axis=0)
        cf = 1.0 / np.sqrt(np.maximum(cn, floor))
        A *= cf[None, :]
        cs *= cf
    return A, rs, cs


def bottom_gaps(sv, max_d=24):
    out = []
    for d in range(1, min(max_d, len(sv) - 1) + 1):
        hi = float(sv[-d - 1])
        lo = float(sv[-d])
        out.append((hi / max(lo, 1e-300), d, hi, lo))
    out.sort(reverse=True)
    return out


def best_gap(sv, max_d=24):
    return bottom_gaps(sv, max_d=max_d)[0]


def observable_jacobian(base, fixed, step):
    return central_jacobian(lambda q: surface_observable_from_base(q, fixed), base, step)


def physical_image_spectrum(Jeq, cs, sv, Vt, d, base_dim, Jobs):
    """Map the d-dimensional equilibrated tangent space to d(A,B).

    Jeq = R J C.  If v is a right-null vector of Jeq, then C v is a
    corresponding tangent vector in the original coordinates.  We QR the
    resulting original-coordinate tangent basis before applying Jobs.
    """
    Veq = Vt[-d:].T                       # nvars x d
    Torig = cs[:, None] * Veq             # C v
    Q, _ = np.linalg.qr(Torig, mode="reduced")
    O = Jobs @ Q[:base_dim, :]
    Oeq, _, _ = equilibrate(O)
    so = np.linalg.svd(Oeq, compute_uv=False)
    return so, O


def image_rank_report(so, tangent_dim):
    if len(so) == 0:
        return 0, float("nan"), []
    # Here the desired quantity is rank, so inspect gaps in descending singular
    # values between retained and near-zero image directions.
    candidates = []
    for r in range(1, min(tangent_dim, len(so))):
        gap = float(so[r - 1]) / max(float(so[r]), 1e-300)
        candidates.append((gap, r))
    candidates.sort(reverse=True)
    if not candidates:
        return min(tangent_dim, len(so)), float("inf"), []
    gap, rank = candidates[0]
    return rank, gap, candidates[:6]


ap = argparse.ArgumentParser()
ap.add_argument(
    "--root",
    type=Path,
    default=BASE / "results/coxeter9-slope-numeric-v1/root-000029",
)
ap.add_argument(
    "--refined",
    type=Path,
    default=BASE / "results/rank10-winning-refine-v1",
)
ap.add_argument("--v-pairings", default="0,0,2,2,2,1,1,1,1")
ap.add_argument("--steps", default="1e-4,3e-5,1e-5,3e-6,1e-6")
ap.add_argument("--max-gap-d", type=int, default=24)
ap.add_argument(
    "--out",
    type=Path,
    default=BASE / "results/rank9-rank10-tangent-calibration-v1",
)
args = ap.parse_args()

root = args.root.resolve()
refined = args.refined.resolve()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)

meta = json.loads((root / "candidate.json").read_text())
fixed = fixed_from_meta(meta)
pv = tuple(map(int, args.v_pairings.split(",")))
if len(pv) != 9:
    raise SystemExit("--v-pairings must contain nine integers")
anchors = anchors_from_pv(pv)
saved = np.load(refined / "state.npy")
if len(saved) != 112:
    raise SystemExit(f"expected rank10 state length 112, got {len(saved)}")
base = saved[:104].copy()
z10 = explicit_state(saved, fixed, anchors)

print("rank9_variables =", len(base))
print("rank10_explicit_variables =", len(z10))
print("rank10_anchors =", len(anchors))
print("rank9_raw_max = %.6e" % float(np.max(np.abs(base_raw(base, fixed)))))
print("rank10_raw_max = %.6e" % float(np.max(np.abs(rank10_raw_explicit(z10, fixed, anchors)))))
print("expected physical dimensions: rank9=9 rank10=8")
print()

steps = [float(x) for x in args.steps.split(",") if x.strip()]
rows = []
for step in steps:
    J9 = central_jacobian(lambda q: base_raw(q, fixed), base, step)
    E9, _, c9 = equilibrate(J9)
    U9, s9, Vt9 = np.linalg.svd(E9, full_matrices=False)
    g9, d9, hi9, lo9 = best_gap(s9, args.max_gap_d)

    J10 = central_jacobian(lambda q: rank10_raw_explicit(q, fixed, anchors), z10, step)
    E10, _, c10 = equilibrate(J10)
    U10, s10, Vt10 = np.linalg.svd(E10, full_matrices=False)
    g10, d10, hi10, lo10 = best_gap(s10, args.max_gap_d)

    Jobs = observable_jacobian(base, fixed, step)
    so9, _ = physical_image_spectrum(E9, c9, s9, Vt9, d9, 104, Jobs)
    so10, _ = physical_image_spectrum(E10, c10, s10, Vt10, d10, 104, Jobs)
    r9, og9, top9 = image_rank_report(so9, d9)
    r10, og10, top10 = image_rank_report(so10, d10)

    print(
        "CALIBRATE|step=%.1e|d9=%d|gap9=%.3e|d10=%d|gap10=%.3e|drop=%d|surface9=%d|surface10=%d|surface_drop=%d"
        % (step, d9, g9, d10, g10, d9 - d10, r9, r10, r9 - r10),
        flush=True,
    )
    print(
        "SURFACE_GAPS|step=%.1e|rank9_gap=%.3e|rank10_gap=%.3e|rank9_tail=%s|rank10_tail=%s"
        % (
            step,
            og9,
            og10,
            " ".join("%.3e" % x for x in so9[-min(12, len(so9)) :]),
            " ".join("%.3e" % x for x in so10[-min(12, len(so10)) :]),
        )
    )

    np.savetxt(out / ("sv-rank9-%.0e.txt" % step), s9.reshape(1, -1), fmt="%.17g")
    np.savetxt(out / ("sv-rank10-%.0e.txt" % step), s10.reshape(1, -1), fmt="%.17g")
    np.savetxt(out / ("surface-image-sv-rank9-%.0e.txt" % step), so9.reshape(1, -1), fmt="%.17g")
    np.savetxt(out / ("surface-image-sv-rank10-%.0e.txt" % step), so10.reshape(1, -1), fmt="%.17g")

    rows.append((step, d9, g9, d10, g10, d9 - d10, r9, r10, r9 - r10))

print()
print("CROSS-STEP CALIBRATION SUMMARY")
for step, d9, g9, d10, g10, drop, r9, r10, rdrop in rows:
    print(
        "DROP|step=%.1e|parameter=%d->%d|parameter_drop=%d|surface=%d->%d|surface_drop=%d|gap9=%.3e|gap10=%.3e"
        % (step, d9, d10, drop, r9, r10, rdrop, g9, g10)
    )

(out / "summary.txt").write_text(
    "\n".join(
        "step=%.1e parameter=%d->%d drop=%d surface=%d->%d surface_drop=%d gap9=%.6e gap10=%.6e"
        % row
        for row in rows
    )
    + "\n"
)
print("saved =", out)
