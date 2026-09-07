from sage.all import *
from pathlib import Path
import argparse
import json


ap = argparse.ArgumentParser(
    description="Fast reconstruction of a scanned A10/MW3 P1 chart hit from pinned formulas."
)
ap.add_argument("--meta", required=True)
ap.add_argument("--point", required=True)
ap.add_argument("--out-json", default=None)
args = ap.parse_args()

meta_lines = Path(args.meta).read_text().splitlines()
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
assignments = {}
for item in args.point.split(","):
    name, value = item.split("=", 1)
    assignments[name.strip()] = K(int(value.strip()))
if set(assignments) != set(active_names):
    raise RuntimeError(f"point must give exactly {active_names}")
values = {d[name]: value for name, value in assignments.items()}
values[d["y5"]] = K(0)

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
    rhs = rhs.split("    #", 1)[0]
    derived_text[name.strip()] = rhs.strip()

derived_expressions = {name: RF(rhs) for name, rhs in derived_text.items()}
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

Kt = PolynomialRing(K, "t")
t = Kt.gen()
A = Kt([derived_values[f"a{i}"] for i in range(9)])
X = Kt([
    derived_values["s0"], derived_values["x1"], assignments["x2"],
    assignments["x3"], derived_values["sinf"],
])
Y = Kt([
    K(0), derived_values["y1"], derived_values["y2"],
    derived_values["y3"], assignments["y4"], K(0),
])
B = Y**2 - X**3 - A * X
Delta = -16 * (4 * A**3 + 27 * B**2)


def valuation_at(poly, point):
    factor = t - point
    valuation = 0
    while poly(point) == 0:
        poly //= factor
        valuation += 1
    return valuation


lam = assignments["lam"]
valuations = [
    valuation_at(Delta, K(0)), valuation_at(Delta, K(1)),
    valuation_at(Delta, lam), 24 - Delta.degree(),
]
fixed = t**3 * (t - 1)**2 * (t - lam)**2
residual = Delta // fixed if valuations == [3, 2, 2, 11] else Kt(0)
squarefree = residual.degree() == 6 and gcd(residual, residual.derivative()).degree() == 0
valid = valuations == [3, 2, 2, 11] and squarefree

record = {
    "p": p,
    "point": {name: int(value) for name, value in assignments.items()},
    "valid_semistable": bool(valid),
    "valuations": [int(value) for value in valuations],
    "residual_squarefree": bool(squarefree),
    "A": [int(c) for c in A.list()],
    "B": [int(c) for c in B.list()],
    "X1": [int(c) for c in X.list()],
    "Y1": [int(c) for c in Y.list()],
    "nodes": [int(3), int(assignments["s1"]), int(derived_values["sl"])],
    "sinf": int(derived_values["sinf"]),
}
print("MW3A10FAST|" + json.dumps(record, sort_keys=True), flush=True)
if args.out_json:
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
