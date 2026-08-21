from sage.all import *
from pathlib import Path
import argparse
import ast
import json
import re
import time


ap = argparse.ArgumentParser(
    description=(
        "Reconstruct many normalized A10/P1 hits and apply the canonical P2 "
        "component-6 and target-intersection gates in one Sage process."
    )
)
ap.add_argument("--seed-start", type=int, required=True)
ap.add_argument("--seed-end", type=int, required=True)
ap.add_argument("--dir", default="artifacts/local/elkies-k3/mw3-a10-p1")
ap.add_argument(
    "--formula-meta",
    default="artifacts/local/elkies-k3/mw3-a10-p1/p31-component2-valid.meta.txt",
)
ap.add_argument("--target-jsonl", default=None)
args = ap.parse_args()

root = Path(args.dir)
hit_re = re.compile(r"MW3A10SCAN_HIT\|(.+)")
jobs = []
for seed in range(args.seed_start, args.seed_end + 1):
    meta = root / f"component2-valid-seed{seed}.meta.txt"
    scan = root / f"component2-valid-seed{seed}.scan.log"
    if not meta.exists() or not scan.exists():
        continue
    values_line = next(
        line.split("=", 1)[1]
        for line in meta.read_text().splitlines()
        if line.startswith("values=")
    )
    fixed = ast.literal_eval(values_line)
    hit_index = 0
    for line in scan.read_text().splitlines():
        match = hit_re.fullmatch(line)
        if not match:
            continue
        hit_index += 1
        hit = {}
        for item in re.split(r"[|,]", match.group(1)):
            name, value = item.split("=", 1)
            hit[name] = int(value)
        jobs.append((seed, hit_index, {**fixed, **hit}))

# Parse the reconstruction formulas once.  This is the dominant repeated
# startup cost in the older one-process-per-hit batch driver.
meta_lines = Path(args.formula_meta).read_text().splitlines()
p = int(next(line.split("=", 1)[1] for line in meta_lines if line.startswith("prime=")))
K = GF(p)
base_names = [
    "rho", "r1", "s1", "lam", "x2", "x3", "y3", "y4", "y5",
    "a4", "a5", "a6",
]
R = PolynomialRing(K, base_names, order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)
active_names = ["rho", "r1", "s1", "lam", "x2", "x3", "y4"]

derived_text = {}
inside = False
for line in meta_lines:
    if line == "DERIVED":
        inside = True
        continue
    if line == "EQUATIONS":
        break
    if not inside or " <- " not in line:
        continue
    name, rhs = line.split(" <- ", 1)
    derived_text[name.strip()] = rhs.split("    #", 1)[0].strip()
derived_expressions = {name: RF(rhs) for name, rhs in derived_text.items()}

Kt = PolynomialRing(K, "t")
t0 = Kt.gen()
FF = FractionField(Kt)


def valuation_at(poly, point):
    factor = t0 - point
    valuation = 0
    while poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


def reconstruct(point):
    assignments = {name: K(point[name]) for name in active_names}
    values = {d[name]: value for name, value in assignments.items()}
    values[d["y5"]] = K(0)
    derived_values = {}
    pending = dict(derived_expressions)
    for _ in range(len(pending) + 2):
        progress = False
        for name, expression in list(pending.items()):
            variables = set(expression.numerator().variables()) | set(
                expression.denominator().variables()
            )
            if any(variable not in values for variable in variables):
                continue
            value = K(expression.subs(values))
            derived_values[name] = value
            if name in d:
                values[d[name]] = value
            del pending[name]
            progress = True
        if not pending or not progress:
            break
    if pending:
        raise RuntimeError(f"could not resolve derived coordinates: {sorted(pending)}")

    A = Kt([derived_values[f"a{i}"] for i in range(9)])
    X1 = Kt([
        derived_values["s0"], derived_values["x1"], assignments["x2"],
        assignments["x3"], derived_values["sinf"],
    ])
    Y1 = Kt([
        K(0), derived_values["y1"], derived_values["y2"],
        derived_values["y3"], assignments["y4"], K(0),
    ])
    B = Y1**2 - X1**3 - A * X1
    Delta = -16 * (4 * A**3 + 27 * B**2)
    lam = assignments["lam"]
    valuations = [
        valuation_at(Delta, K(0)), valuation_at(Delta, K(1)),
        valuation_at(Delta, lam), 24 - Delta.degree(),
    ]
    fixed_discriminant = t0**3 * (t0 - 1)**2 * (t0 - lam)**2
    residual = Delta // fixed_discriminant if valuations == [3, 2, 2, 11] else Kt(0)
    valid = (
        valuations == [3, 2, 2, 11]
        and residual.degree() == 6
        and gcd(residual, residual.derivative()).degree() == 0
    )
    return {
        "point": {name: int(value) for name, value in assignments.items()},
        "valid_semistable": bool(valid),
        "valuations": [int(value) for value in valuations],
        "A": A,
        "B": B,
        "X1": X1,
        "Y1": Y1,
        "nodes": [K(3), assignments["s1"], derived_values["sl"]],
        "sinf": derived_values["sinf"],
    }


def p2_profile_orientation(P1, P2):
    """Return raw (I11,I3) labels and the sign giving target (6,2), if any."""
    translated = P1 + P2
    x = FF(translated[0])
    y = FF(translated[1])
    if x.denominator()(0) == 0 or y.denominator()(0) == 0:
        c3 = 2
    else:
        c3 = 1 if K(x(0)) == K(3) and K(y(0)) == K(0) else 2

    degree_at_infinity = y.numerator().degree() - y.denominator().degree()
    distance = 6 - degree_at_infinity
    if distance not in (3, 4):
        raise RuntimeError(f"unexpected translated P2 I11 distance {distance}")
    c11 = 6 if distance == 3 else 5
    raw = (c11, c3)
    if raw == (6, 2):
        return raw, 1
    if ((11 - c11) % 11, (3 - c3) % 3) == (6, 2):
        return raw, -1
    return raw, 0


# Reuse the P2 coefficient ring across surfaces.
Q = PolynomialRing(K, ("r", "q0", "q1", "q2"), order="degrevlex")
r, q0, q1, q2 = Q.gens()
Qt = PolynomialRing(Q, "t")
t = Qt.gen()


def canonical_p2(record):
    A = Qt(record["A"])
    B = Qt(record["B"])
    lam = K(record["point"]["lam"])
    sinf = K(record["sinf"])
    fiber_points = [K(0), K(1), lam]
    F = t * (t - 1) * (t - lam)
    C = Qt(0)
    for i, (point, node) in enumerate(zip(fiber_points, record["nodes"])):
        basis = Qt(1)
        denominator = K(1)
        for j, other in enumerate(fiber_points):
            if i == j:
                continue
            basis *= t - other
            denominator *= point - other
        C += Q(node) * (Q(point) - r)**2 * basis / Q(denominator)

    X = C + F * (q0 + q1 * t + q2 * t**2 + Q(sinf) * t**3)
    z = t - r
    H = X**3 + A * X * z**4 + B * z**6
    if H[18] != 0 or H[17] != 0:
        raise RuntimeError("P2 infinity incidence failed")
    equations = [Q(H[k]) for k in range(16, 8, -1)]
    solutions = Q.ideal(equations).variety()
    E = EllipticCurve(FF, [0, 0, 0, FF(record["A"]), FF(record["B"])])
    P1_point = E(FF(record["X1"]), FF(record["Y1"]))

    hits = []
    for solution in solutions:
        pole = K(solution[r])
        if pole in fiber_points:
            continue
        X_special = Kt([K(coefficient.subs(solution)) for coefficient in X.list()])
        H_special = Kt([K(coefficient.subs(solution)) for coefficient in H.list()])
        factorization = H_special.factor()
        roots = K(factorization.unit()).sqrt(all=True)
        if not roots or any(exponent % 2 for _, exponent in factorization):
            continue
        Y_special = Kt(roots[0])
        for factor, exponent in factorization:
            Y_special *= factor ** (exponent // 2)
        if any(Y_special(point) != 0 for point in fiber_points):
            continue

        z_special = t0 - pole
        P2_point = E(FF(X_special) / FF(z_special)**2, FF(Y_special) / FF(z_special)**3)
        raw_components, profile_sign = p2_profile_orientation(P1_point, P2_point)

        x2_function = FractionField(Kt)(X_special) / FractionField(Kt)(z_special)**2
        y2_function = FractionField(Kt)(Y_special) / FractionField(Kt)(z_special)**3
        intersections = []
        for sign in (K(1), K(-1)):
            dx = Kt((FractionField(Kt)(record["X1"]) - x2_function).numerator())
            dy = Kt((FractionField(Kt)(record["Y1"]) - sign * y2_function).numerator())
            smooth = gcd(dx, dy)
            for point in fiber_points:
                while smooth(point) == 0:
                    smooth //= t0 - point
            tangent_equal = (
                record["X1"].derivative()(lam) == x2_function.derivative()(lam)
                and record["Y1"].derivative()(lam)
                == (sign * y2_function).derivative()(lam)
            )
            intersections.append(int(smooth.degree()) + int(tangent_equal))
        canonical_intersection = (
            intersections[0] if profile_sign == 1
            else intersections[1] if profile_sign == -1
            else None
        )
        hits.append(
            {
                "r": int(pole),
                "q0": int(solution[q0]),
                "q1": int(solution[q1]),
                "q2": int(solution[q2]),
                "X2": [int(c) for c in X_special.list()],
                "Y2": [int(c) for c in Y_special.list()],
                "intersections": intersections,
                "raw_components": [int(c) for c in raw_components],
                "profile_sign": int(profile_sign),
                "target_profile": profile_sign != 0,
                "target": profile_sign != 0 and canonical_intersection == 1,
            }
        )
    return hits


counts = {"boundary": 0, "no-p2": 0, "P2-HIT": 0, "TARGET-P2-HIT": 0}
targets = []
started = time.monotonic()
for job_number, (seed, hit_index, point) in enumerate(jobs, 1):
    job_started = time.monotonic()
    record = reconstruct(point)
    if not record["valid_semistable"]:
        status = "boundary"
        p2_hits = []
    else:
        p2_hits = canonical_p2(record)
        status = "TARGET-P2-HIT" if any(hit["target"] for hit in p2_hits) else (
            "P2-HIT" if p2_hits else "no-p2"
        )
    counts[status] += 1
    print(
        f"MW3A10MULTIP2|job={job_number}/{len(jobs)}|seed={seed}|hit={hit_index}"
        f"|status={status}|p2_hits={len(p2_hits)}"
        f"|seconds={time.monotonic()-job_started:.3f}",
        flush=True,
    )
    for p2_hit in p2_hits:
        if not p2_hit["target"]:
            continue
        output = {
            "seed": int(seed),
            "hit": int(hit_index),
            "p": int(p),
            "point": record["point"],
            "A": [int(c) for c in record["A"].list()],
            "B": [int(c) for c in record["B"].list()],
            "X1": [int(c) for c in record["X1"].list()],
            "Y1": [int(c) for c in record["Y1"].list()],
            "nodes": [int(c) for c in record["nodes"]],
            "sinf": int(record["sinf"]),
            "P2": p2_hit,
        }
        targets.append(output)
        print("MW3A10MULTIP2_TARGET|" + json.dumps(output, sort_keys=True), flush=True)

if args.target_jsonl:
    target_path = Path(args.target_jsonl)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text("".join(json.dumps(item, sort_keys=True) + "\n" for item in targets))

print(
    "MW3A10MULTIP2|done=1|jobs=" + str(len(jobs))
    + "|targets=" + str(len(targets))
    + "|counts=" + json.dumps({name: int(value) for name, value in counts.items()}, sort_keys=True)
    + f"|seconds={time.monotonic()-started:.3f}",
    flush=True,
)
