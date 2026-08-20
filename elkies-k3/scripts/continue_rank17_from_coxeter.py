#!/usr/bin/env python3
"""Continue a numerical Coxeter-9 K3 root toward the recovered rank-17 lattice.

The Coxeter slope reconstruction gives nine explicit height-4 sections on a
9-dimensional K3 locus (after coordinate gauge).  The exact lattice extension
chain selects eight further height-4 sections.  Existence of each section cuts
the moduli locus, while its |height pairing|=2 links identify the desired
Mordell--Weil component.

For an added section Q with quartic x_Q and an earlier anchor P satisfying
< Q, P > = +/-2, the corresponding sum/difference is again minimal.  Hence
the chord slope is quadratic.  With p=<Q,P> and s=p/2 in {+/-1},

    y_Q + s*y_P = m(T) * (x_Q - x_P),       deg m <= 2.

We use one such relation to *define* y_Q from x_Q and a primary quadratic
slope, impose the curve identity

    y_Q^2 = x_Q^3 + A*x_Q + B,

and use additional +/-2 anchors as polynomial line constraints.  Previous
extensions remain constrained while the Coxeter slope variables are allowed to
move along their numerical moduli directions.

Run rank 10 first.  If it converges cleanly, increase --target-rank one stage at
a time; the solver resumes from the previous rank directory when present.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json
import math

import numpy as np

try:
    from scipy.optimize import least_squares
except ImportError as exc:
    raise SystemExit("scipy is required for this numerical continuation") from exc


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


def pmul(a, b):
    return np.convolve(a, b)


def p3(a):
    return pmul(pmul(a, a), a)


def pder(a):
    if len(a) <= 1:
        return np.zeros(1)
    return np.asarray([k * a[k] for k in range(1, len(a))], dtype=float)


def padd(a, b, sa=1.0, sb=1.0):
    n = max(len(a), len(b))
    out = np.zeros(n)
    out[: len(a)] += sa * a
    out[: len(b)] += sb * b
    return out


def pmax(a):
    return float(np.max(np.abs(a))) if len(a) else 0.0


def pair_slope(slopes, i, j):
    if i > j:
        return -slopes[PAIR_INDEX[(j, i)]]
    return slopes[PAIR_INDEX[(i, j)]]


def reconstruct_base(slopes):
    f = np.empty((len(TRIPLES), 5))
    for r, (i, j, k) in enumerate(TRIPLES):
        mij = pair_slope(slopes, i, j)
        mik = pair_slope(slopes, i, k)
        mjk = pair_slope(slopes, j, k)
        f[r] = pmul(mik, mij) + pmul(mik, mjk) - pmul(mij, mjk)
    x = TINC_PINV @ f
    rx = f - TINC @ x

    rhs_y = np.empty((len(PAIRS), 7))
    for r, (i, j) in enumerate(PAIRS):
        rhs_y[r] = pmul(slopes[r], x[i] - x[j])
    y = EINC_PINV @ rhs_y
    ry = rhs_y - EINC @ y

    a_edges = np.empty((len(PAIRS), 9))
    for r, (i, j) in enumerate(PAIRS):
        a_edges[r] = (
            pmul(slopes[r], y[i] - y[j])
            - pmul(x[i], x[i])
            - pmul(x[i], x[j])
            - pmul(x[j], x[j])
        )
    A = a_edges.mean(axis=0)
    rA = a_edges - A

    b_nodes = np.empty((9, 13))
    for i in range(9):
        b_nodes[i] = pmul(y[i], y[i]) - p3(x[i]) - pmul(A, x[i])
    B = b_nodes.mean(axis=0)
    rB = b_nodes - B
    return x, y, A, B, (rx, ry, rA, rB)


def geom_metrics(x, A, B):
    a3 = p3(A)
    b2 = pmul(B, B)
    disc = padd(a3, b2, 4.0, 27.0)
    Ap = pder(A)
    Bp = pder(B)
    t1 = pmul(Ap, B)
    t2 = pmul(A, Bp)
    jvar = padd(t1, t2, 3.0, -2.0)
    x2 = np.median([pmax(pmul(v, v)) for v in x])
    x3 = np.median([pmax(p3(v)) for v in x])
    eps = 1e-300
    strength = max(pmax(A) / max(x2, eps), pmax(B) / max(x3, eps))
    delta_rel = pmax(disc) / max(x3 * x3, eps)
    jden = 3.0 * pmax(t1) + 2.0 * pmax(t2)
    jvar_rel = pmax(jvar) / max(jden, eps)
    return strength, delta_rel, jvar_rel


def fixed_from_meta(meta):
    sign = int(meta["gauge_sign"])
    second_pair = tuple(sorted(map(int, meta.get("second_pair", [0, 2]))))
    p01 = PAIR_INDEX[(0, 1)]
    p2 = PAIR_INDEX[second_pair]
    return {
        3 * p01 + 0: float(sign),
        3 * p01 + 1: 0.0,
        3 * p01 + 2: 1.0,
        3 * p2 + 1: 1.0,
    }


def base_free_positions(fixed):
    return [i for i in range(108) if i not in fixed]


def pack_base(slopes, fixed):
    flat = np.asarray(slopes).reshape(-1)
    return np.asarray([flat[i] for i in base_free_positions(fixed)])


def unpack_base(free, fixed):
    flat = np.empty(108)
    positions = base_free_positions(fixed)
    flat[positions] = free
    for i, value in fixed.items():
        flat[i] = value
    return flat.reshape(36, 3)


def load_anchor_plan(path: Path):
    groups = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            qpos = int(row["section_position"])
            groups.setdefault(qpos, []).append(
                {
                    "anchor_position": int(row["anchor_position"]),
                    "pairing": int(row["pairing"]),
                    "role": row["role"],
                    "anchor_signed_index": int(row["anchor_signed_index"]),
                    "section_signed_index": int(row["section_signed_index"]),
                }
            )
    for qpos, rows in groups.items():
        rows.sort(key=lambda r: (0 if r["role"] == "primary" else 1, r["anchor_position"]))
        if rows[0]["role"] != "primary":
            raise RuntimeError(f"qpos {qpos} has no primary anchor")
    return groups


def variable_layout(groups, target_rank):
    layout = {}
    pos = 104
    for qpos in range(9, target_rank):
        anchors = groups[qpos]
        entry = {"x": slice(pos, pos + 5)}
        pos += 5
        entry["primary_m"] = slice(pos, pos + 3)
        pos += 3
        secondary = []
        for _ in anchors[1:]:
            secondary.append(slice(pos, pos + 3))
            pos += 3
        entry["secondary_m"] = secondary
        layout[qpos] = entry
    return layout, pos


def decode(z, fixed, groups, target_rank):
    slopes = unpack_base(z[:104], fixed)
    bx, by, A, B, base_blocks = reconstruct_base(slopes)
    xs = [bx[i].copy() for i in range(9)]
    ys = [by[i].copy() for i in range(9)]
    layout, nvars = variable_layout(groups, target_rank)
    if len(z) != nvars:
        raise RuntimeError(f"state length {len(z)} != expected {nvars}")
    extras = {}

    for qpos in range(9, target_rank):
        anchors = groups[qpos]
        L = layout[qpos]
        xq = np.asarray(z[L["x"]])
        m0 = np.asarray(z[L["primary_m"]])
        primary = anchors[0]
        apos = primary["anchor_position"]
        s = primary["pairing"] / 2.0
        dx = xq - xs[apos]
        yq = pmul(m0, dx) - s * ys[apos]
        sec_ms = [np.asarray(z[slice_]) for slice_ in L["secondary_m"]]
        extras[qpos] = {
            "x": xq,
            "y": yq,
            "primary_m": m0,
            "secondary_m": sec_ms,
            "anchors": anchors,
        }
        xs.append(xq)
        ys.append(yq)

    return slopes, xs, ys, A, B, base_blocks, extras


def raw_blocks(z, fixed, groups, target_rank, new_secondary_limit=None):
    slopes, xs, ys, A, B, base_blocks, extras = decode(
        z, fixed, groups, target_rank
    )
    blocks = {
        "base_x": base_blocks[0],
        "base_y": base_blocks[1],
        "base_A": base_blocks[2],
        "base_B": base_blocks[3],
    }

    for qpos in range(9, target_rank):
        q = extras[qpos]
        curve = pmul(q["y"], q["y"]) - p3(q["x"]) - pmul(A, q["x"]) - B
        blocks[f"q{qpos}_curve"] = curve

        secondaries = q["anchors"][1:]
        ms = q["secondary_m"]
        use = len(secondaries)
        if qpos == target_rank - 1 and new_secondary_limit is not None:
            use = min(use, new_secondary_limit)
        for n in range(use):
            anchor = secondaries[n]
            apos = anchor["anchor_position"]
            s = anchor["pairing"] / 2.0
            lhs = q["y"] + s * ys[apos]
            rhs = pmul(ms[n], q["x"] - xs[apos])
            blocks[f"q{qpos}_line{n}"] = lhs - rhs

    return slopes, xs, ys, A, B, blocks


def residual(z, fixed, groups, target_rank, args, new_secondary_limit=None):
    slopes, xs, ys, A, B, blocks = raw_blocks(
        z, fixed, groups, target_rank, new_secondary_limit
    )
    pieces = []
    for name, block in blocks.items():
        if name == "base_A":
            w = 0.5
        elif name == "base_B":
            w = 0.25
        elif name.endswith("_curve"):
            w = args.curve_weight
        else:
            w = 1.0
        pieces.append(w * np.asarray(block).reshape(-1))

    strength, delta_rel, jvar_rel = geom_metrics(np.asarray(xs[:9]), A, B)
    # Numerical branch guards only.  They vanish identically once the surface
    # stays away from the cuspidal/isotrivial components.
    pieces.append(np.asarray([
        args.guard_weight * max(0.0, args.min_surface_strength - strength),
        args.guard_weight * max(0.0, args.min_delta_rel - delta_rel),
        args.guard_weight * max(0.0, args.min_jvar_rel - jvar_rel),
    ]))

    # Prevent the newly introduced x-coordinate from collapsing exactly onto
    # its primary anchor during branch discovery.
    if target_rank > 9:
        qpos = target_rank - 1
        q = decode(z, fixed, groups, target_rank)[-1][qpos]
        primary = q["anchors"][0]
        apos = primary["anchor_position"]
        dxn = float(np.linalg.norm(q["x"] - xs[apos]))
        scale = max(1.0, float(np.median([np.linalg.norm(x) for x in xs[:9]])))
        pieces.append(np.asarray([
            args.guard_weight * max(0.0, args.min_distinctness * scale - dxn)
        ]))

    return np.concatenate(pieces)


def all_abs2_validation(xs, ys, Gsel, target_rank):
    worst = 0.0
    count = 0
    for qpos in range(9, target_rank):
        for apos in range(qpos):
            pairing = int(Gsel[qpos, apos])
            if abs(pairing) != 2:
                continue
            dx = xs[qpos] - xs[apos]
            lhs = ys[qpos] + (pairing / 2.0) * ys[apos]
            # convolution matrix mapping quadratic m to m*dx
            M = np.zeros((7, 3))
            for j in range(3):
                M[j : j + 5, j] = dx
            m, *_ = np.linalg.lstsq(M, lhs, rcond=None)
            err = float(np.max(np.abs(lhs - M @ m)))
            worst = max(worst, err)
            count += 1
    return worst, count


def diagnostics(z, fixed, groups, Gsel, target_rank):
    slopes, xs, ys, A, B, blocks = raw_blocks(z, fixed, groups, target_rank)
    maxima = {name: float(np.max(np.abs(block))) for name, block in blocks.items()}
    raw_max = max(maxima.values()) if maxima else 0.0
    base_max = max(v for k, v in maxima.items() if k.startswith("base_"))
    curve_max = max([v for k, v in maxima.items() if k.endswith("_curve")] or [0.0])
    line_max = max([v for k, v in maxima.items() if "_line" in k] or [0.0])
    strength, delta_rel, jvar_rel = geom_metrics(np.asarray(xs[:9]), A, B)
    abs2_max, abs2_count = all_abs2_validation(xs, ys, Gsel, target_rank)
    return {
        "raw_max": raw_max,
        "base_max": base_max,
        "curve_max": curve_max,
        "line_max": line_max,
        "all_abs2_max": abs2_max,
        "all_abs2_count": abs2_count,
        "surface_strength": strength,
        "delta_rel": delta_rel,
        "jvar_rel": jvar_rel,
        "slope_rank": int(np.linalg.matrix_rank(slopes, tol=1e-9)),
    }, (slopes, xs, ys, A, B)


def save_stage(directory, z, fixed, groups, Gsel, target_rank, metadata):
    directory.mkdir(parents=True, exist_ok=True)
    stats, data = diagnostics(z, fixed, groups, Gsel, target_rank)
    slopes, xs, ys, A, B = data
    np.save(directory / "state.npy", z)
    np.savetxt(directory / "slopes.txt", slopes, fmt="%.17g")
    np.savetxt(directory / "x-selected.txt", np.asarray(xs), fmt="%.17g")
    np.savetxt(directory / "y-selected.txt", np.asarray(ys), fmt="%.17g")
    np.savetxt(directory / "A.txt", A.reshape(1, -1), fmt="%.17g")
    np.savetxt(directory / "B.txt", B.reshape(1, -1), fmt="%.17g")
    payload = dict(metadata)
    payload.update(stats)
    payload["target_rank"] = target_rank
    payload["variables"] = len(z)
    (directory / "candidate.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return stats


def initialize_rank10(root_slopes, fixed, groups, target_rank, rng, scale, previous=None):
    layout, nvars = variable_layout(groups, target_rank)
    z = np.zeros(nvars)
    z[:104] = pack_base(root_slopes, fixed) if previous is None else previous[:104]

    if previous is not None:
        z[: len(previous)] = previous

    # Decode the already-present part to get realistic x coefficient scales.
    previous_rank = target_rank - 1
    if previous_rank == 9:
        bx, _, _, _, _ = reconstruct_base(unpack_base(z[:104], fixed))
        existing_x = bx
    else:
        _, xs, _, _, _, _, _ = decode(previous, fixed, groups, previous_rank)
        existing_x = np.asarray(xs)

    qpos = target_rank - 1
    L = layout[qpos]
    coeff = rng.normal(size=len(existing_x))
    coeff /= max(np.linalg.norm(coeff), 1e-12)
    xguess = coeff @ existing_x
    xguess += rng.normal(0.0, 0.15 * scale, size=5)
    z[L["x"]] = xguess

    slope_scale = max(0.25, float(np.median(np.abs(root_slopes))))
    z[L["primary_m"]] = rng.normal(0.0, slope_scale, size=3)
    for sl in L["secondary_m"]:
        z[sl] = rng.normal(0.0, slope_scale, size=3)
    return z


parser = argparse.ArgumentParser()
parser.add_argument(
    "--root",
    type=Path,
    default=BASE / "results" / "coxeter9-slope-numeric-v1" / "root-000001",
)
parser.add_argument(
    "--chain",
    type=Path,
    default=BASE / "results" / "rank17-extension-chain-v1",
)
parser.add_argument(
    "--audit",
    type=Path,
    default=BASE / "results" / "rank17-extension-audit-v1",
)
parser.add_argument("--target-rank", type=int, default=10)
parser.add_argument("--restarts", type=int, default=12)
parser.add_argument("--max-nfev", type=int, default=1800)
parser.add_argument("--seed", type=int, default=20260820)
parser.add_argument("--init-scale", type=float, default=0.5)
parser.add_argument("--curve-weight", type=float, default=0.5)
parser.add_argument("--guard-weight", type=float, default=10.0)
parser.add_argument("--min-surface-strength", type=float, default=1e-3)
parser.add_argument("--min-delta-rel", type=float, default=1e-6)
parser.add_argument("--min-jvar-rel", type=float, default=1e-8)
parser.add_argument("--min-distinctness", type=float, default=1e-3)
parser.add_argument(
    "--out",
    type=Path,
    default=BASE / "results" / "rank17-continuation-v1",
)
args = parser.parse_args()

if not 10 <= args.target_rank <= 17:
    raise SystemExit("--target-rank must be between 10 and 17")

root = args.root.resolve()
chain = args.chain.resolve()
audit = args.audit.resolve()
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)

for p in (
    root / "slopes.txt",
    root / "candidate.json",
    chain / "selected-gram.txt",
    audit / "continuation-anchors.tsv",
):
    if not p.exists():
        raise SystemExit(f"missing required input: {p}")

root_slopes = np.loadtxt(root / "slopes.txt", dtype=float).reshape(36, 3)
root_meta = json.loads((root / "candidate.json").read_text())
fixed = fixed_from_meta(root_meta)
groups = load_anchor_plan(audit / "continuation-anchors.tsv")
Gsel = np.loadtxt(chain / "selected-gram.txt", dtype=np.int64)

rng = np.random.default_rng(args.seed)

print("root =", root)
print("root raw_max =", root_meta.get("raw_max"))
print("target_rank =", args.target_rank)
print("base free variables = 104")
print()

previous = None
start_rank = 10
for r in range(10, args.target_rank + 1):
    stage_dir = out / f"rank{r}"
    state_path = stage_dir / "state.npy"
    if state_path.exists():
        previous = np.load(state_path)
        print(f"RESUME|rank={r}|state={state_path}")
        continue

    if r > 10:
        prev_path = out / f"rank{r-1}" / "state.npy"
        if not prev_path.exists():
            raise SystemExit(f"missing previous continuation state: {prev_path}")
        previous = np.load(prev_path)

    best = None
    for restart in range(args.restarts):
        z0 = initialize_rank10(
            root_slopes,
            fixed,
            groups,
            r,
            rng,
            args.init_scale,
            previous=previous,
        )
        if restart:
            # Slightly perturb the already solved locus as well, but much less
            # than the newly introduced section variables.
            z0[:104] += rng.normal(0.0, 2e-3, size=104)

        # First use one secondary anchor for the newly introduced section.
        stage_results = []
        limits = [1, None]
        z = z0
        for limit in limits:
            fun = lambda zz, lim=limit: residual(
                zz, fixed, groups, r, args, new_secondary_limit=lim
            )
            res = least_squares(
                fun,
                z,
                method="trf",
                jac="2-point",
                max_nfev=args.max_nfev,
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                x_scale="jac",
                verbose=0,
            )
            z = res.x
            stage_results.append(
                {
                    "secondary_limit": limit,
                    "nfev": int(res.nfev),
                    "cost_norm": float(np.linalg.norm(res.fun)),
                }
            )

        stats, _ = diagnostics(z, fixed, groups, Gsel, r)
        usable = (
            stats["surface_strength"] >= args.min_surface_strength
            and stats["delta_rel"] >= args.min_delta_rel
            and stats["jvar_rel"] >= args.min_jvar_rel
        )
        score = max(stats["raw_max"], stats["all_abs2_max"])
        is_best = usable and (best is None or score < best[0])

        print(
            f"CONT|rank={r}"
            f"|restart={restart}"
            f"|score={score:.3e}"
            f"|raw={stats['raw_max']:.3e}"
            f"|base={stats['base_max']:.3e}"
            f"|curve={stats['curve_max']:.3e}"
            f"|line={stats['line_max']:.3e}"
            f"|allabs2={stats['all_abs2_max']:.3e}"
            f"|strength={stats['surface_strength']:.3e}"
            f"|delta={stats['delta_rel']:.3e}"
            f"|jvar={stats['jvar_rel']:.3e}"
            f"|usable={int(usable)}"
            + ("|BEST" if is_best else ""),
            flush=True,
        )

        if is_best:
            best = (score, z.copy(), stats, restart, stage_results)

    if best is None:
        raise SystemExit(f"no usable continuation root found for rank {r}")

    score, zbest, stats, restart, stage_results = best
    metadata = {
        "version": 1,
        "kind": "rank17-coxeter-continuation",
        "source_root": str(root),
        "restart": restart,
        "seed": args.seed,
        "stage_results": stage_results,
    }
    save_stage(stage_dir, zbest, fixed, groups, Gsel, r, metadata)
    previous = zbest
    print(
        f"ACCEPT|rank={r}|score={score:.3e}|state={stage_dir / 'state.npy'}",
        flush=True,
    )
    print()

print("DONE")
print("saved =", out)
