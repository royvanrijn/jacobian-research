#!/usr/bin/env python3
"""Rank Coxeter-9 numerical roots by geometric usefulness.

The slope solver deliberately records every small-residual root.  Many of those
belong to the singular/cuspidal component A=B=0, so residual size alone is not
an appropriate ranking criterion.

For every ``root-*`` directory this script computes:

* scale-aware A/B and discriminant strengths;
* the polynomial non-isotriviality invariant 3*A'*B - 2*A*B';
* effective degrees of A, B and Delta;
* finite-difference Jacobian singular values for the full gauge-fixed residual;
* numerical tangent nullity at several relative SVD thresholds.

For a genuine Coxeter-9 elliptic K3 scaffold, the expected post-gauge moduli
nullity is about 9.  The full rank-17 locus should later cut this to dimension 1.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json
import math

import numpy as np


BASE = Path(__file__).resolve().parents[1]
PAIRS = [(i, j) for i in range(9) for j in range(i + 1, 9)]
PAIR_INDEX = {p: n for n, p in enumerate(PAIRS)}
TRIPLES = [(i, j, k) for i in range(9) for j in range(i + 1, 9) for k in range(j + 1, 9)]

TINC = np.zeros((len(TRIPLES), 9), dtype=float)
for r, (i, j, k) in enumerate(TRIPLES):
    TINC[r, [i, j, k]] = 1.0
TINC_PINV = np.linalg.inv(TINC.T @ TINC) @ TINC.T

EINC = np.zeros((len(PAIRS), 9), dtype=float)
for r, (i, j) in enumerate(PAIRS):
    EINC[r, [i, j]] = 1.0
EINC_PINV = np.linalg.inv(EINC.T @ EINC) @ EINC.T


def pmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.convolve(a, b)


def p3(a: np.ndarray) -> np.ndarray:
    return pmul(pmul(a, a), a)


def pder(a: np.ndarray) -> np.ndarray:
    if len(a) <= 1:
        return np.zeros(1, dtype=float)
    return np.asarray([k * a[k] for k in range(1, len(a))], dtype=float)


def padd(a: np.ndarray, b: np.ndarray, sa: float = 1.0, sb: float = 1.0) -> np.ndarray:
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=float)
    out[: len(a)] += sa * a
    out[: len(b)] += sb * b
    return out


def pmax(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if len(a) else 0.0


def effective_degree(a: np.ndarray, rel: float = 1e-8) -> int:
    scale = pmax(a)
    if scale == 0.0:
        return -1
    nz = np.flatnonzero(np.abs(a) > rel * scale)
    return int(nz[-1]) if len(nz) else -1


def pair_slope(slopes: np.ndarray, i: int, j: int) -> np.ndarray:
    if i > j:
        return -slopes[PAIR_INDEX[(j, i)]]
    return slopes[PAIR_INDEX[(i, j)]]


def reconstruct(slopes: np.ndarray):
    f = np.empty((len(TRIPLES), 5), dtype=float)
    for r, (i, j, k) in enumerate(TRIPLES):
        mij = pair_slope(slopes, i, j)
        mik = pair_slope(slopes, i, k)
        mjk = pair_slope(slopes, j, k)
        f[r] = pmul(mik, mij) + pmul(mik, mjk) - pmul(mij, mjk)

    x = TINC_PINV @ f
    rx = f - TINC @ x

    rhs_y = np.empty((len(PAIRS), 7), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        rhs_y[r] = pmul(slopes[r], x[i] - x[j])
    y = EINC_PINV @ rhs_y
    ry = rhs_y - EINC @ y

    a_edges = np.empty((len(PAIRS), 9), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        a_edges[r] = (
            pmul(slopes[r], y[i] - y[j])
            - pmul(x[i], x[i])
            - pmul(x[i], x[j])
            - pmul(x[j], x[j])
        )
    A = a_edges.mean(axis=0)
    rA = a_edges - A

    b_nodes = np.empty((9, 13), dtype=float)
    for i in range(9):
        b_nodes[i] = pmul(y[i], y[i]) - p3(x[i]) - pmul(A, x[i])
    B = b_nodes.mean(axis=0)
    rB = b_nodes - B

    return x, y, A, B, (rx, ry, rA, rB)


def raw_residual(slopes: np.ndarray) -> np.ndarray:
    _, _, _, _, blocks = reconstruct(slopes)
    return np.concatenate([b.reshape(-1) for b in blocks])


def gauge_fixed_positions(sign: int, second_pair: tuple[int, int]) -> dict[int, float]:
    p01 = PAIR_INDEX[(0, 1)]
    p2 = PAIR_INDEX[second_pair]
    return {
        3 * p01 + 0: float(sign),
        3 * p01 + 1: 0.0,
        3 * p01 + 2: 1.0,
        3 * p2 + 1: 1.0,
    }


def free_positions(fixed: dict[int, float]) -> list[int]:
    return [i for i in range(108) if i not in fixed]


def jacobian_singular_values(slopes: np.ndarray, fixed: dict[int, float], step: float) -> np.ndarray:
    flat = np.asarray(slopes, dtype=float).reshape(-1)
    positions = free_positions(fixed)
    base_shape = raw_residual(slopes).shape[0]
    J = np.empty((base_shape, len(positions)), dtype=float)

    for col, pos in enumerate(positions):
        h = step * max(1.0, abs(float(flat[pos])))
        plus = flat.copy()
        minus = flat.copy()
        plus[pos] += h
        minus[pos] -= h
        fp = raw_residual(plus.reshape(len(PAIRS), 3))
        fm = raw_residual(minus.reshape(len(PAIRS), 3))
        J[:, col] = (fp - fm) / (2.0 * h)

    return np.linalg.svd(J, compute_uv=False)


def geom_metrics(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> dict[str, float | int | bool]:
    a3 = p3(A)
    b2 = pmul(B, B)
    disc_core = padd(a3, b2, 4.0, 27.0)

    Ap = pder(A)
    Bp = pder(B)
    term1 = pmul(Ap, B)
    term2 = pmul(A, Bp)
    jvar = padd(term1, term2, 3.0, -2.0)

    x2_scales = np.asarray([pmax(pmul(v, v)) for v in x])
    x3_scales = np.asarray([pmax(p3(v)) for v in x])
    x2_scale = float(np.median(x2_scales))
    x3_scale = float(np.median(x3_scales))

    eps = 1e-300
    a_strength = pmax(A) / max(x2_scale, eps)
    b_strength = pmax(B) / max(x3_scale, eps)
    surface_strength = max(a_strength, b_strength)
    disc_vs_x6 = pmax(disc_core) / max(x3_scale * x3_scale, eps)

    jden = 3.0 * pmax(term1) + 2.0 * pmax(term2)
    jvar_rel = pmax(jvar) / max(jden, eps)

    # A=B=0 is the ubiquitous cuspidal component.  The fixed slope gauge makes
    # x-scale meaningful, so this threshold is scale-aware rather than absolute.
    cusp_like = surface_strength < 1e-7 or disc_vs_x6 < 1e-12
    nonisotrivial = (not cusp_like) and jvar_rel > 1e-6

    return {
        "A_max": pmax(A),
        "B_max": pmax(B),
        "A_strength": a_strength,
        "B_strength": b_strength,
        "surface_strength": surface_strength,
        "delta_max": pmax(disc_core),
        "delta_vs_x6": disc_vs_x6,
        "jvar_max": pmax(jvar),
        "jvar_rel": jvar_rel,
        "A_degree": effective_degree(A),
        "B_degree": effective_degree(B),
        "delta_degree": effective_degree(disc_core),
        "cusp_like": cusp_like,
        "nonisotrivial": nonisotrivial,
    }


parser = argparse.ArgumentParser()
parser.add_argument(
    "--roots",
    type=Path,
    default=BASE / "results" / "coxeter9-slope-numeric-v1",
    help="directory containing root-* subdirectories",
)
parser.add_argument("--jacobian-step", type=float, default=2e-6)
parser.add_argument("--skip-jacobian", action="store_true")
parser.add_argument("--top", type=int, default=20)
args = parser.parse_args()

root_dir = args.roots.resolve()
roots = sorted(p for p in root_dir.glob("root-*") if p.is_dir())
if not roots:
    raise SystemExit(f"no root-* directories under {root_dir}")

rows = []
for n, d in enumerate(roots, start=1):
    meta_path = d / "candidate.json"
    slopes_path = d / "slopes.txt"
    if not meta_path.exists() or not slopes_path.exists():
        print(f"SKIP|root={d.name}|reason=missing_files", flush=True)
        continue

    meta = json.loads(meta_path.read_text())
    slopes = np.loadtxt(slopes_path, dtype=float).reshape(len(PAIRS), 3)
    x, y, A, B, blocks = reconstruct(slopes)
    raw_max = max(float(np.max(np.abs(b))) for b in blocks)
    gm = geom_metrics(x, A, B)

    sign = int(meta["gauge_sign"])
    second_pair = tuple(map(int, meta.get("second_pair", [0, 2])))
    second_pair = tuple(sorted(second_pair))
    fixed = gauge_fixed_positions(sign, second_pair)

    nullities = {}
    smin = float("nan")
    smax = float("nan")
    if not args.skip_jacobian:
        sv = jacobian_singular_values(slopes, fixed, args.jacobian_step)
        smax = float(sv[0]) if len(sv) else 0.0
        smin = float(sv[-1]) if len(sv) else 0.0
        for exponent in (6, 8, 10, 12):
            tol = (10.0 ** (-exponent)) * max(smax, 1e-300)
            rank = int(np.count_nonzero(sv > tol))
            nullities[f"nullity_1e-{exponent}"] = len(free_positions(fixed)) - rank
        np.savetxt(d / "jacobian-singular-values.txt", sv.reshape(1, -1), fmt="%.17g")

    row = {
        "root": d.name,
        "run_index": int(meta.get("run_index", -1)),
        "gauge_sign": sign,
        "restart": int(meta.get("restart", -1)),
        "raw_max": raw_max,
        "slope_rank": int(np.linalg.matrix_rank(slopes, tol=1e-9)),
        **gm,
        **nullities,
        "jacobian_smax": smax,
        "jacobian_smin": smin,
    }
    rows.append(row)

    print(
        f"ROOT|{d.name}"
        f"|raw={raw_max:.3e}"
        f"|cusp={int(bool(gm['cusp_like']))}"
        f"|noniso={int(bool(gm['nonisotrivial']))}"
        f"|strength={float(gm['surface_strength']):.3e}"
        f"|delta={float(gm['delta_vs_x6']):.3e}"
        f"|jvar={float(gm['jvar_rel']):.3e}"
        + (f"|null8={nullities.get('nullity_1e-8', -1)}" if nullities else ""),
        flush=True,
    )


def rank_key(r):
    # Non-cuspidal, non-isotrivial roots first.  Among them prefer the expected
    # rank-9 scaffold tangent dimension, then stronger discriminant and residual.
    n8 = int(r.get("nullity_1e-8", 999))
    strength = max(float(r["surface_strength"]), 1e-300)
    delta = max(float(r["delta_vs_x6"]), 1e-300)
    return (
        bool(r["cusp_like"]),
        not bool(r["nonisotrivial"]),
        abs(n8 - 9),
        -math.log10(strength),
        -math.log10(delta),
        float(r["raw_max"]),
    )

rows.sort(key=rank_key)
out_tsv = root_dir / "root-diagnostics.tsv"
fields = list(rows[0].keys())
with out_tsv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

print()
print("TOP ROOTS")
for r in rows[: args.top]:
    print(
        f"RANK|root={r['root']}"
        f"|run={r['run_index']}"
        f"|raw={float(r['raw_max']):.3e}"
        f"|cusp={int(bool(r['cusp_like']))}"
        f"|noniso={int(bool(r['nonisotrivial']))}"
        f"|strength={float(r['surface_strength']):.3e}"
        f"|delta={float(r['delta_vs_x6']):.3e}"
        f"|jvar={float(r['jvar_rel']):.3e}"
        f"|deg=({r['A_degree']},{r['B_degree']},{r['delta_degree']})"
        f"|null8={r.get('nullity_1e-8', 'NA')}"
        f"|null10={r.get('nullity_1e-10', 'NA')}"
    )

print("saved =", out_tsv)
