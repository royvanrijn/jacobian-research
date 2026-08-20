#!/usr/bin/env python3
"""Numerically reconstruct the rank-17 K3 from the Coxeter-9 slope system.

This is the second reconstruction stage after
``build_coxeter9_x_reconstruction.py``.

Let V_0,...,V_8 be the nine height-4 clique sections and define, for i<j,

    D_ij = V_i - V_j.

Write m_ij(T) for the slope of the line through V_i and -V_j.  Each m_ij is
quadratic.  Associativity of the chord-and-tangent law gives the important
coherent-slope identity

    slope(D_ij,D_jk) = m_ik - m_ij - m_jk.

Consequently all 84 triangle slopes are determined by the 36 pair slopes.
Eliminating the D_ij x-coordinates from

    x(V_i)+x(V_j)+x(D_ij) = m_ij^2

and

    x(D_ij)+x(D_jk)+x(D_ik)
        = (m_ik-m_ij-m_jk)^2

gives, for every i<j<k,

    x_i+x_j+x_k
      = m_ik*m_ij + m_ik*m_jk - m_ij*m_jk.              (1)

Thus 36 quadratic slopes (108 scalar coefficients) determine the nine quartic
x_i linearly when the 84 right hand sides of (1) lie in the column space of
the 84x9 triple-incidence matrix.

Once x is known, the chosen orientation gives

    y_i+y_j = m_ij * (x_i-x_j),                            (2)

so the nine sextic y_i are also recovered linearly from the 36 pair equations.
Finally every edge must produce one common degree-8 A(T),

    A = m_ij*(y_i-y_j) - (x_i^2+x_i*x_j+x_j^2),           (3)

and every vertex must produce one common degree-12 B(T),

    B = y_i^2 - x_i^3 - A*x_i.                            (4)

The optimizer therefore has only the slope coefficients as nonlinear
variables.  x, y, A and B are reconstructed at every residual evaluation.

There is a four-dimensional coordinate gauge (PGL_2 on the base plus
Weierstrass scaling) and one genuine K3-moduli dimension.  For numerical
search we fix

    m_01(T) = T^2 + s,  s in {-1,+1}

and the T coefficient of m_02 to 1.  Both signs can be searched.  This gauge
is only a numerical discovery device; any candidate is subsequently meant to
be exactified and transformed to a rational gauge.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json
import math

import numpy as np

try:
    from scipy.optimize import least_squares
except ImportError as exc:
    raise SystemExit(
        "scipy is required; run this with ordinary python or sage -python "
        "from an environment containing scipy"
    ) from exc


BASE = Path(__file__).resolve().parents[1]

PAIRS = [(i, j) for i in range(9) for j in range(i + 1, 9)]
PAIR_INDEX = {p: n for n, p in enumerate(PAIRS)}
TRIPLES = [
    (i, j, k)
    for i in range(9)
    for j in range(i + 1, 9)
    for k in range(j + 1, 9)
]

# Triple-to-vertex incidence: each row represents x_i+x_j+x_k.
TINC = np.zeros((len(TRIPLES), 9), dtype=float)
for r, (i, j, k) in enumerate(TRIPLES):
    TINC[r, [i, j, k]] = 1.0
TINC_PINV = np.linalg.inv(TINC.T @ TINC) @ TINC.T

# Edge-to-vertex plus-incidence: each row represents y_i+y_j.
EINC = np.zeros((len(PAIRS), 9), dtype=float)
for r, (i, j) in enumerate(PAIRS):
    EINC[r, [i, j]] = 1.0
EINC_PINV = np.linalg.inv(EINC.T @ EINC) @ EINC.T


def pmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Multiply low-to-high coefficient vectors."""
    return np.convolve(a, b)


def p3(a: np.ndarray) -> np.ndarray:
    return pmul(pmul(a, a), a)


def peval(a: np.ndarray, t: float) -> float:
    z = 0.0
    for c in a[::-1]:
        z = z * t + float(c)
    return z


def pair_slope(slopes: np.ndarray, i: int, j: int) -> np.ndarray:
    if i > j:
        # m_ji is the slope through V_j and -V_i.  From the definition
        # (y_j+y_i)/(x_j-x_i) it is -m_ij.
        return -slopes[PAIR_INDEX[(j, i)]]
    return slopes[PAIR_INDEX[(i, j)]]


def reconstruct(slopes: np.ndarray):
    """Recover x,y,A,B and raw consistency residual blocks from slopes.

    slopes has shape (36,3), with coefficients [constant,T,T^2].
    """
    # ------------------------------------------------------------
    # x from the 84 triangle identities (1).
    # ------------------------------------------------------------
    f = np.empty((len(TRIPLES), 5), dtype=float)

    for r, (i, j, k) in enumerate(TRIPLES):
        mij = pair_slope(slopes, i, j)
        mik = pair_slope(slopes, i, k)
        mjk = pair_slope(slopes, j, k)

        f[r] = (
            pmul(mik, mij)
            + pmul(mik, mjk)
            - pmul(mij, mjk)
        )

    # For each quartic coefficient, solve TINC*x_k=f_k.
    x = TINC_PINV @ f  # shape 9x5
    rx = f - TINC @ x

    # ------------------------------------------------------------
    # y from pair equations (2).
    # ------------------------------------------------------------
    rhs_y = np.empty((len(PAIRS), 7), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        rhs_y[r] = pmul(slopes[r], x[i] - x[j])

    y = EINC_PINV @ rhs_y  # shape 9x7
    ry = rhs_y - EINC @ y

    # ------------------------------------------------------------
    # Every pair must give the same A via (3).
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Every vertex must give the same B via (4).
    # ------------------------------------------------------------
    b_nodes = np.empty((9, 13), dtype=float)
    for i in range(9):
        b_nodes[i] = (
            pmul(y[i], y[i])
            - p3(x[i])
            - pmul(A, x[i])
        )

    B = b_nodes.mean(axis=0)
    rB = b_nodes - B

    return x, y, A, B, (rx, ry, rA, rB)


def block_stats(blocks) -> dict[str, float]:
    names = ("x", "y", "A", "B")
    out = {}
    for name, block in zip(names, blocks):
        out[f"{name}_rms"] = float(np.sqrt(np.mean(block * block)))
        out[f"{name}_max"] = float(np.max(np.abs(block)))
    out["raw_max"] = max(out[f"{name}_max"] for name in names)
    return out


def residual_from_slopes(slopes: np.ndarray, stage: str) -> np.ndarray:
    _, _, _, _, blocks = reconstruct(slopes)
    rx, ry, rA, rB = blocks

    if stage == "x":
        selected = ((1.0, rx),)
    elif stage == "xy":
        selected = ((1.0, rx), (1.0, ry))
    elif stage == "xya":
        selected = ((1.0, rx), (1.0, ry), (0.5, rA))
    elif stage == "full":
        selected = (
            (1.0, rx),
            (1.0, ry),
            (0.5, rA),
            (0.25, rB),
        )
    else:
        raise ValueError(stage)

    return np.concatenate([
        weight * block.reshape(-1)
        for weight, block in selected
    ])


def make_gauge(sign: int, second_pair: tuple[int, int]):
    """Return fixed flat-coordinate dictionary for the numerical gauge."""
    if sign not in (-1, 1):
        raise ValueError(sign)
    if second_pair not in PAIR_INDEX:
        raise ValueError(second_pair)
    if second_pair == (0, 1):
        raise ValueError("second gauge pair must differ from (0,1)")

    p01 = PAIR_INDEX[(0, 1)]
    p2 = PAIR_INDEX[second_pair]

    # Low-to-high: m01 = sign + 0*T + 1*T^2.
    fixed = {
        3 * p01 + 0: float(sign),
        3 * p01 + 1: 0.0,
        3 * p01 + 2: 1.0,
        # Kill the remaining one-dimensional stabilizer numerically.
        3 * p2 + 1: 1.0,
    }
    return fixed


def unpack(free: np.ndarray, fixed: dict[int, float]) -> np.ndarray:
    flat = np.empty(108, dtype=float)
    free_positions = [i for i in range(108) if i not in fixed]
    flat[free_positions] = free
    for i, value in fixed.items():
        flat[i] = value
    return flat.reshape(len(PAIRS), 3)


def pack(slopes: np.ndarray, fixed: dict[int, float]) -> np.ndarray:
    flat = np.asarray(slopes, dtype=float).reshape(-1)
    return np.asarray([flat[i] for i in range(108) if i not in fixed])


def j_samples(A: np.ndarray, B: np.ndarray):
    values = []
    for t in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0):
        a = peval(A, t)
        b = peval(B, t)
        den = 4.0 * a**3 + 27.0 * b**2
        if not math.isfinite(den) or abs(den) < 1e-20:
            values.append(float("nan"))
        else:
            values.append(1728.0 * 4.0 * a**3 / den)
    finite = np.asarray([v for v in values if math.isfinite(v)], dtype=float)
    spread = float(np.std(finite)) if len(finite) else float("nan")
    return values, spread


def self_test() -> None:
    """Check all sign/orientation formulas on a constant real elliptic curve."""
    A0 = -1.0
    B0 = 1.0
    xs = np.arange(9, dtype=float)
    ys = np.sqrt(xs**3 + A0 * xs + B0)

    slopes = np.zeros((len(PAIRS), 3), dtype=float)
    for r, (i, j) in enumerate(PAIRS):
        # slope through P_i and -P_j
        slopes[r, 0] = (ys[i] + ys[j]) / (xs[i] - xs[j])

    x, y, A, B, blocks = reconstruct(slopes)
    stats = block_stats(blocks)

    err_x = float(np.max(np.abs(x[:, 0] - xs)))
    err_y = float(np.max(np.abs(y[:, 0] - ys)))
    err_A = abs(float(A[0]) - A0) + float(np.max(np.abs(A[1:])))
    err_B = abs(float(B[0]) - B0) + float(np.max(np.abs(B[1:])))

    print("SELFTEST")
    print("x_error =", err_x)
    print("y_error =", err_y)
    print("A_error =", err_A)
    print("B_error =", err_B)
    print("raw_max =", stats["raw_max"])

    if max(err_x, err_y, err_A, err_B, stats["raw_max"]) > 1e-8:
        raise SystemExit("self-test failed")
    print("PASS")


def save_candidate(out: Path, tag: str, slopes, x, y, A, B, stats, metadata):
    d = out / tag
    d.mkdir(parents=True, exist_ok=True)

    np.savetxt(d / "slopes.txt", slopes, fmt="%.17g")
    np.savetxt(d / "x-vertices.txt", x, fmt="%.17g")
    np.savetxt(d / "y-vertices.txt", y, fmt="%.17g")
    np.savetxt(d / "A.txt", A.reshape(1, -1), fmt="%.17g")
    np.savetxt(d / "B.txt", B.reshape(1, -1), fmt="%.17g")

    with (d / "slopes.tsv").open("w") as handle:
        handle.write("i\tj\tc0\tc1\tc2\n")
        for (i, j), m in zip(PAIRS, slopes):
            handle.write(
                f"{i}\t{j}\t{m[0]:.17g}\t{m[1]:.17g}\t{m[2]:.17g}\n"
            )

    payload = dict(metadata)
    payload.update(stats)
    payload["slope_coefficient_rank"] = int(np.linalg.matrix_rank(slopes))
    js, spread = j_samples(A, B)
    payload["j_samples"] = js
    payload["j_sample_std"] = spread

    (d / "candidate.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )


parser = argparse.ArgumentParser()
parser.add_argument("--restarts", type=int, default=16)
parser.add_argument("--seed", type=int, default=20260820)
parser.add_argument("--max-nfev", type=int, default=1500)
parser.add_argument("--init-scale", type=float, default=0.75)
parser.add_argument(
    "--gauge-sign",
    choices=("split", "nonsplit", "both"),
    default="both",
    help="m01=T^2-1, T^2+1, or try both",
)
parser.add_argument(
    "--second-pair",
    default="0,2",
    help="pair whose T slope coefficient is fixed to 1, e.g. 0,2",
)
parser.add_argument(
    "--out",
    type=Path,
    default=BASE / "results" / "coxeter9-slope-numeric-v1",
)
parser.add_argument("--self-test", action="store_true")
args = parser.parse_args()

if args.self_test:
    self_test()
    raise SystemExit(0)

try:
    second_pair = tuple(map(int, args.second_pair.split(",")))
except Exception as exc:
    raise SystemExit("--second-pair must look like 0,2") from exc
if len(second_pair) != 2:
    raise SystemExit("--second-pair must contain two indices")
second_pair = tuple(sorted(second_pair))

if args.gauge_sign == "split":
    signs = (-1,)
elif args.gauge_sign == "nonsplit":
    signs = (1,)
else:
    signs = (-1, 1)

out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(args.seed)
best = None
run_index = 0

print("pairs =", len(PAIRS))
print("triples =", len(TRIPLES))
print("nonlinear scalar slope coefficients =", 108)
print("gauge fixed coefficients =", 4)
print("free variables =", 104)
print("second gauge pair =", second_pair)
print()

for sign in signs:
    fixed = make_gauge(sign, second_pair)

    def fun_for(stage):
        def f(free):
            return residual_from_slopes(unpack(free, fixed), stage)
        return f

    for restart in range(args.restarts):
        run_index += 1

        # Random free coefficients; fixed coordinates are inserted by unpack.
        z = rng.normal(0.0, args.init_scale, size=104)

        # Stage the increasingly strong geometric conditions.  The first stage
        # discovers a square-compatible x configuration; later stages enforce
        # coherent y coordinates and finally a single common Weierstrass model.
        stage_results = []
        for stage in ("x", "xy", "xya", "full"):
            res = least_squares(
                fun_for(stage),
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
                (stage, int(res.nfev), float(np.linalg.norm(res.fun)))
            )

        slopes = unpack(z, fixed)
        x, y, A, B, blocks = reconstruct(slopes)
        stats = block_stats(blocks)
        slope_rank = int(np.linalg.matrix_rank(slopes, tol=1e-9))
        js, jspread = j_samples(A, B)

        score = stats["raw_max"]
        is_best = best is None or score < best[0]

        print(
            f"RUN|n={run_index}"
            f"|sign={sign:+d}"
            f"|restart={restart}"
            f"|raw_max={score:.6e}"
            f"|x={stats['x_max']:.3e}"
            f"|y={stats['y_max']:.3e}"
            f"|A={stats['A_max']:.3e}"
            f"|B={stats['B_max']:.3e}"
            f"|slope_rank={slope_rank}"
            f"|j_std={jspread:.6g}"
            + ("|BEST" if is_best else ""),
            flush=True,
        )

        if is_best:
            metadata = {
                "version": 1,
                "kind": "coxeter9-slope-numeric",
                "run_index": run_index,
                "restart": restart,
                "seed": args.seed,
                "gauge_sign": sign,
                "second_pair": list(second_pair),
                "gauge": "m01=T^2+sign; second-pair T coefficient=1",
                "stage_results": stage_results,
            }
            save_candidate(
                out,
                "best",
                slopes,
                x,
                y,
                A,
                B,
                stats,
                metadata,
            )
            best = (score, sign, restart, stats, slope_rank, jspread)

        # Preserve genuinely small roots separately, including degenerate ones;
        # the rank/j-variation diagnostics tell us whether they are interesting.
        if score < 1e-7:
            metadata = {
                "version": 1,
                "kind": "coxeter9-slope-numeric-root",
                "run_index": run_index,
                "restart": restart,
                "seed": args.seed,
                "gauge_sign": sign,
                "second_pair": list(second_pair),
                "stage_results": stage_results,
            }
            save_candidate(
                out,
                f"root-{run_index:06d}",
                slopes,
                x,
                y,
                A,
                B,
                stats,
                metadata,
            )

print()
if best is None:
    raise SystemExit("no runs completed")

score, sign, restart, stats, slope_rank, jspread = best
print("BEST")
print("raw_max =", score)
print("gauge_sign =", sign)
print("restart =", restart)
print("slope_coefficient_rank =", slope_rank)
print("j_sample_std =", jspread)
for key in sorted(stats):
    print(key, "=", stats[key])
print("saved =", out / "best")
