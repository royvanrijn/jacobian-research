#!/usr/bin/env python3
"""Measure the *incremental* geometric condition imposed by Q_10.

Previous tangent diagnostics compared strongest SVD gaps in two different
parameterizations:

  rank 9  : 104 Coxeter variables
  rank 10 : 104 Coxeter + 5 x_Q + 27 explicit anchor-slope variables

Those raw nullity counts are not directly comparable.  The invariant local
question is instead:

    inside the rank-9 tangent space, how many independent conditions remain
    after the coordinates of Q_10 are allowed to adjust?

This script answers that by linearizing the reduced rank-10 equations at the
saved healthy refinement.  The reduced Q variables are only

    5 coefficients of x_Q + 3 coefficients of the primary quadratic slope.

All eight secondary anchor slopes are eliminated by least-squares projection,
exactly as in refine_rank10_winning_fingerprint.py.

Let T_9 be the numerical Coxeter tangent space, A the derivative of the Q
constraints with respect to base variables, and B their derivative with
respect to Q variables.  A base tangent u survives iff there exists dq with

    A T_9 u + B dq = 0.

Equivalently, after projecting to the left nullspace of B, surviving base
motions satisfy

    N_B^T A T_9 u = 0.

The rank of this projected matrix is the incremental codimension.  For one
new independent Mordell--Weil section we expect

    incremental_codimension = 1.

The script repeats the calculation across finite-difference step sizes and
also reports the tangent dimension of the 112-variable reduced rank-10 system.
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
    for i, value in fixed.items():
        flat[i] = value
    return flat.reshape(36, 3)


def base_raw(base, fixed):
    *_, blocks = reconstruct(unpack(base, fixed))
    return np.concatenate(
        [
            blocks[0].ravel(),
            blocks[1].ravel(),
            0.5 * blocks[2].ravel(),
            0.25 * blocks[3].ravel(),
        ]
    )


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


def projection_error(lhs, dx):
    M = conv_matrix(dx)
    m, *_ = np.linalg.lstsq(M, lhs, rcond=None)
    return lhs - M @ m


def q_raw(base, qstate, fixed, anchors):
    """Only the new-section equations; qstate=[x_Q(5), primary slope(3)]."""
    s = unpack(base, fixed)
    xs, ys, A, B, _ = derived45(s)
    xq = qstate[:5]
    m0 = qstate[5:8]
    primary = choose_primary(anchors)
    ordered = [primary] + [a for a in anchors if a != primary]
    ai, p = ordered[0]
    yq = pmul(m0, xq - xs[ai]) - (p / 2.0) * ys[ai]

    pieces = [0.5 * (pmul(yq, yq) - p3(xq) - pmul(A, xq) - B)]
    for aj, pj in ordered[1:]:
        pieces.append(projection_error(yq + (pj / 2.0) * ys[aj], xq - xs[aj]))
    return np.concatenate(pieces)


def reduced_rank10_raw(z, fixed, anchors):
    base = z[:104]
    qstate = z[104:112]
    return np.concatenate([base_raw(base, fixed), q_raw(base, qstate, fixed, anchors)])


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
    """Ruiz row/column equilibration: E = diag(rs) J diag(cs)."""
    E = np.asarray(J, dtype=float).copy()
    rs = np.ones(E.shape[0], dtype=float)
    cs = np.ones(E.shape[1], dtype=float)
    for _ in range(iterations):
        rn = np.linalg.norm(E, axis=1)
        rf = 1.0 / np.sqrt(np.maximum(rn, floor))
        E *= rf[:, None]
        rs *= rf
        cn = np.linalg.norm(E, axis=0)
        cf = 1.0 / np.sqrt(np.maximum(cn, floor))
        E *= cf[None, :]
        cs *= cf
    return E, rs, cs


def strongest_bottom_gap(sv, max_d=24):
    candidates = []
    for d in range(1, min(max_d, len(sv) - 1) + 1):
        hi = float(sv[-d - 1])
        lo = float(sv[-d])
        candidates.append((hi / max(lo, 1e-300), d, hi, lo))
    return max(candidates)


def original_tangent_basis(J, d):
    E, _, cs = equilibrate(J)
    _, sv, Vt = np.linalg.svd(E, full_matrices=False)
    Veq = Vt[-d:].T
    T = cs[:, None] * Veq
    T, _ = np.linalg.qr(T, mode="reduced")
    return T, sv


def normalize_columns(M, floor=1e-300):
    n = np.linalg.norm(M, axis=0)
    scale = 1.0 / np.maximum(n, floor)
    return M * scale[None, :], scale


def rank_at(sv, reltol):
    if len(sv) == 0 or float(sv[0]) == 0.0:
        return 0
    return int(np.count_nonzero(sv > reltol * float(sv[0])))


def projected_increment(A, B, reltol):
    """Return rank(B), rank(N_B^T A), spectra, using common row scaling."""
    M = np.hstack([A, B])
    rn = np.linalg.norm(M, axis=1)
    rs = 1.0 / np.maximum(rn, 1e-300)
    Ar = A * rs[:, None]
    Br = B * rs[:, None]
    Ar, _ = normalize_columns(Ar)
    Br, _ = normalize_columns(Br)

    U, sb, _ = np.linalg.svd(Br, full_matrices=True)
    rank_b = rank_at(sb, reltol)
    left_null = U[:, rank_b:]
    C = left_null.T @ Ar
    sc = np.linalg.svd(C, compute_uv=False)
    inc = rank_at(sc, reltol)
    return rank_b, inc, sb, sc


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
ap.add_argument("--steps", default="1e-4,3e-5,1e-5,3e-6,1e-6,3e-7")
ap.add_argument("--rank-tols", default="1e-5,1e-6,1e-7,1e-8,1e-9")
ap.add_argument("--max-gap-d", type=int, default=24)
ap.add_argument(
    "--out",
    type=Path,
    default=BASE / "results/rank10-incremental-codimension-v1",
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
z = np.load(refined / "state.npy")
if len(z) != 112:
    raise SystemExit(f"expected saved rank10 state length 112, got {len(z)}")
base = z[:104].copy()
qstate = z[104:112].copy()

print("base_variables = 104")
print("q_variables_reduced = 8")
print("anchors =", len(anchors))
print("base_raw_max = %.6e" % float(np.max(np.abs(base_raw(base, fixed)))))
print("q_raw_max = %.6e" % float(np.max(np.abs(q_raw(base, qstate, fixed, anchors)))))
print("expected_incremental_codimension = 1")
print()

steps = [float(x) for x in args.steps.split(",") if x.strip()]
tols = [float(x) for x in args.rank_tols.split(",") if x.strip()]
summary = []
for step in steps:
    J9 = central_jacobian(lambda b: base_raw(b, fixed), base, step)
    E9, _, _ = equilibrate(J9)
    _, s9, _ = np.linalg.svd(E9, full_matrices=False)
    gap9, d9, hi9, lo9 = strongest_bottom_gap(s9, args.max_gap_d)
    T9, _ = original_tangent_basis(J9, d9)

    # Derivatives of only the new-section constraints.
    A = central_jacobian(lambda b: q_raw(b, qstate, fixed, anchors), base, step)
    B = central_jacobian(lambda q: q_raw(base, q, fixed, anchors), qstate, step)
    AT = A @ T9

    # Reduced full system, with no explicit secondary-slope variables.
    Jr = central_jacobian(lambda zz: reduced_rank10_raw(zz, fixed, anchors), z, step)
    Er, _, _ = equilibrate(Jr)
    sr = np.linalg.svd(Er, compute_uv=False)
    gapr, dr, hir, lor = strongest_bottom_gap(sr, args.max_gap_d)

    print(
        "STEP|h=%.1e|d9=%d|gap9=%.3e|reduced_d10=%d|reduced_gap=%.3e"
        % (step, d9, gap9, dr, gapr),
        flush=True,
    )

    chosen = None
    for tol in tols:
        rank_b, inc, sb, sc = projected_increment(AT, B, tol)
        qfiber = 8 - rank_b
        predicted = d9 - inc + qfiber
        print(
            "INC|h=%.1e|tol=%.0e|rankB=%d|qfiber=%d|incremental_codim=%d|predicted_reduced_d10=%d|observed_reduced_d10=%d|B_relmin=%.3e|C_s1=%.3e|C_s2=%.3e"
            % (
                step,
                tol,
                rank_b,
                qfiber,
                inc,
                predicted,
                dr,
                float(sb[-1] / max(sb[0], 1e-300)) if len(sb) else 0.0,
                float(sc[0]) if len(sc) else 0.0,
                float(sc[1]) if len(sc) > 1 else 0.0,
            ),
            flush=True,
        )
        if abs(tol - 1e-7) < 1e-20:
            chosen = (rank_b, inc, predicted)

    if chosen is None:
        rank_b, inc, _, _ = projected_increment(AT, B, 1e-7)
        chosen = (rank_b, inc, d9 - inc + (8 - rank_b))
    summary.append((step, d9, gap9, dr, gapr, *chosen))

    np.savetxt(out / ("sv-rank9-%.0e.txt" % step), s9.reshape(1, -1), fmt="%.17g")
    np.savetxt(out / ("sv-rank10-reduced-%.0e.txt" % step), sr.reshape(1, -1), fmt="%.17g")

print()
print("CROSS-STEP INCREMENTAL SUMMARY (tol=1e-7)")
for step, d9, gap9, dr, gapr, rank_b, inc, predicted in summary:
    print(
        "CODIM|h=%.1e|d9=%d|rankB=%d|incremental_codim=%d|predicted_d10=%d|observed_reduced_d10=%d|gap9=%.3e|gap10=%.3e"
        % (step, d9, rank_b, inc, predicted, dr, gap9, gapr)
    )

(out / "summary.txt").write_text(
    "\n".join(
        "h=%.1e d9=%d gap9=%.6e reduced_d10=%d gap10=%.6e rankB=%d incremental_codim=%d predicted_d10=%d"
        % row
        for row in summary
    )
    + "\n"
)
print("saved =", out)
