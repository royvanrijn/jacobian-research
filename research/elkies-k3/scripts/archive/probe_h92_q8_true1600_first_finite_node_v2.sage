#!/usr/bin/env sage -python
"""
H3 q=8 true-1600 global survivor + first exact finite E7 node quotient test.

This reconstructs the documented p=43 global matrix on the true uniform h^-18
envelope, with the E7 helper normalization u_power = d - 9, and refuses to
continue unless it reproduces the known ranks

    E8 + H:        rank 1446, kernel 154
    E8 + H + E7:  rank 1582, kernel 18.

It then runs the repository's exact node-principal clearing and Singular local
normal-form machinery at E7_4--E7_3, restricts that finite quotient image to
the 18 global survivors, and reports the new rank/codimension.

Run:
    sage -python ~/Downloads/probe_h92_q8_true1600_first_finite_node.sage

Optional:
    sage -python ... --repo /path/to/jacobian-research --prime 43
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, binomial, gcd, matrix


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit(
        "Could not locate jacobian-research. Re-run with "
        "--repo /path/to/jacobian-research"
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, default=None)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--node", default="E7_4--E7_3")
args = parser.parse_args()

ROOT = locate_repo(args.repo)
P = int(args.prime)
if P != 43:
    print(
        "NOTE: regression ranks are pinned for p=43; for another prime the "
        "script still computes but does not enforce the p43 rank tuple.",
        flush=True,
    )

GEN = ROOT / "artifacts" / "generated-results"
SCRIPTS = ROOT / "elkies-k3" / "scripts"
P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
GENERIC_PATH = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"

# Keep these inside ROOT because several audited repository scripts deliberately
# record paths relative to ROOT.
PREFIX = "zz-h92-q8-true1600"
AMBIENT_PATH = GEN / f"{PREFIX}-ambient.json"
COND_PATH = GEN / f"{PREFIX}-generic-conditions.json"
R13_PATH = GEN / f"{PREFIX}-e7-1-3-residue-rows.json"
R56_PATH = GEN / f"{PREFIX}-e7-5-6-residue-rows.json"
R47_PATH = GEN / f"{PREFIX}-e7-4-7-residue-rows.json"
KERNEL_PATH = GEN / f"{PREFIX}-global-kernel-mod-{P}.json"
CLEAR_PATH = GEN / f"{PREFIX}-node-principal-clearings.json"
NODE_PATH = GEN / (
    f"{PREFIX}-{args.node.replace('--','-')}-local-normal-form-mod-{P}.json"
)
RESULT_PATH = GEN / (
    f"{PREFIX}-{args.node.replace('--','-')}-restricted-mod-{P}.json"
)

sage = shutil.which("sage")
if not sage:
    raise SystemExit("sage executable not found")


def run_repo_script(name, *extra):
    command = [sage, "-python", str(SCRIPTS / name), *map(str, extra)]
    print("RUN|" + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def qq_to_finite(value, finite):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ZeroDivisionError(
            f"prime {P} divides denominator of exact coefficient {value}"
        )
    return finite(ZZ(value.numerator())) / den


# ---------------------------------------------------------------------------
# 1. True 1600-column ambient.
#
# Actual reduced q8 coefficients are u^d/h^18.
# The E7 q6^9 helper compiler must see i=d-9 exactly once.
# ---------------------------------------------------------------------------

generic = json.loads(GENERIC_PATH.read_text())
assert generic["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert len(generic["basis"]) == 18

ambient_basis = []
for family in generic["basis"]:
    a = int(family["x_power"])
    b = int(family["m_power"])
    if a == 0:
        assert 0 <= b <= 9
        max_d = 87
    else:
        assert a == 1 and 0 <= b <= 7
        max_d = 89
    for d in range(max_d + 1):
        ambient_basis.append(
            {
                "kind": family["kind"],
                "x_power": a,
                "m_power": b,
                # CRITICAL: helper-normalized exponent, not actual d.
                "u_power": d - 9,
                "h_power": 18,
                "actual_u_power": d,
                "coefficient": f"u^{d}/h(u)^18",
            }
        )

assert len(ambient_basis) == 1600
ambient = {
    "schema": "elkies-k3.h92-q8-true1600-helper-ambient.v1",
    "status": "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
    "ambient_dimension": 1600,
    "ambient_basis": ambient_basis,
    "normalization": {
        "actual_coefficient": "u^d/h(u)^18*x^a*m^b",
        "stored_u_power": "d-9",
        "reason": (
            "Repository E7 compilers include the q6^9 helper -9 shift; "
            "true reduced q8 d is emulated by helper i=d-9."
        ),
        "global_minus_11F": "applied exactly once in the E8 block",
    },
}
AMBIENT_PATH.write_text(json.dumps(ambient, indent=2, sort_keys=True) + "\n")
print(f"TRUEAMBIENT|path={AMBIENT_PATH}|columns=1600|h=18", flush=True)

# Exact generic E7 valuation template + exact function-field residue rows.
run_repo_script(
    "derive_h92_q8_all_component_generic_conditions.sage",
    "--ambient", AMBIENT_PATH,
    "--output", COND_PATH,
)
run_repo_script(
    "derive_h92_q8_e7_1_3_generic_residue_rows.sage",
    "--conditions", COND_PATH,
    "--output", R13_PATH,
)
run_repo_script(
    "derive_h92_q8_e7_5_6_generic_residue_rows.sage",
    "--conditions", COND_PATH,
    "--output", R56_PATH,
)
run_repo_script(
    "derive_h92_q8_e7_4_7_generic_residue_rows.sage",
    "--conditions", COND_PATH,
    "--output", R47_PATH,
)

# ---------------------------------------------------------------------------
# 2. Direct p-adic/global condition matrices on the same 1600 coefficients.
# ---------------------------------------------------------------------------

finite = GF(P)
N = len(ambient_basis)
u_ring = PolynomialRing(finite, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()

p1 = json.loads(P1_PATH.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"


def finite_poly(coefficients):
    return u_ring([qq_to_finite(v, finite) for v in coefficients])


h = finite_poly(p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and h(0)
assert gcd(h, h.derivative()) == 1

x_p = u_field(finite_poly(
    p1["x_entrance_base"]["numerator_coefficients"]
)) / u_field(finite_poly(
    p1["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(finite_poly(
    p1["y_entrance_base"]["numerator_coefficients"]
)) / u_field(finite_poly(
    p1["y_entrance_base"]["denominator_coefficients"]
))


def residue_mod(value, modulus):
    value = u_field(value)
    num = u_ring(value.numerator())
    den = u_ring(value.denominator())
    if gcd(den, modulus) != 1:
        raise ZeroDivisionError(
            f"non-unit denominator modulo {modulus}: gcd={gcd(den, modulus)}"
        )
    return u_ring((num * den.inverse_mod(modulus)) % modulus)


# ---- Smooth h-collision block.
#
# q=(m-p)/h, X=h^2*x.  For k=18 the worst h-pole is 27, hence
# 18 * deg(h^27) = 18 * 108 = 1944 rows.
H_POLE = 27
h_modulus = h ** H_POLE
h_residue_dimension = h_modulus.degree()
assert h_residue_dimension == 108
rho = u_field(h) * y_p / x_p
assert gcd(h, u_ring(rho.numerator())) == 1
assert gcd(h, u_ring(rho.denominator())) == 1

frame_coordinates = (
    [(0, j) for j in range(10)] + [(1, j) for j in range(8)]
)
coord_index = {v: i for i, v in enumerate(frame_coordinates)}
H_ROWS = len(frame_coordinates) * h_residue_dimension
assert H_ROWS == 1944
H = matrix(finite, H_ROWS, N)

for col, entry in enumerate(ambient_basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    for j in range(b + 1):
        exponent = 2 * j - b - 18 - 2 * a
        if exponent >= 0:
            continue
        value = (
            finite(binomial(b, j))
            * u ** d
            * rho ** (b - j)
            * h ** (H_POLE + exponent)
        )
        remainder = residue_mod(value, h_modulus)
        offset = coord_index[(a, j)] * h_residue_dimension
        for degree, coefficient in enumerate(remainder.list()):
            if coefficient:
                H[offset + degree, col] += coefficient

# ---- E8 block.
#
# C=[1,Q,...,Q^9,X,...,XQ^7], Q=u^2*m, X=u^4*x.
# Reduced E8 module has ideal (u^2,X,Y), colength 2.
# Apply the literal -11F exactly once here:
#
#   m^b floor     d >= 11 + 2b
#   x*m^b floor   d >= 15 + 2b
#
# These give 178 + 11*18 = 376 one-coordinate conditions, plus the two
# quotient jets from evaluation Q=s=Y_P/X_P.
E8 = matrix(finite, 378, N)
floor_row = 0
for col, entry in enumerate(ambient_basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    floor = (11 + 2 * b) if a == 0 else (15 + 2 * b)
    if d < floor:
        E8[floor_row, col] = 1
        floor_row += 1
assert floor_row == 376

u2 = u ** 2
s_e8 = u_field(u ** 2) * y_p / x_p
# s=Y_P/X_P must be a local E8 unit.
s0 = residue_mod(s_e8, u2)
assert s0[0]

for col, entry in enumerate(ambient_basis):
    if int(entry["x_power"]) != 0:
        continue
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    floor = 11 + 2 * b
    if d < floor or d - floor >= 2:
        continue
    value = u_field(u ** (d - floor)) * s_e8 ** b / u_field(h ** 18)
    remainder = residue_mod(value, u2)
    for jet in range(2):
        E8[376 + jet, col] = remainder[jet]

# ---- Exact generic E7 rows from repository compilers.
conditions = json.loads(COND_PATH.read_text())
singleton_indices = list(
    map(int, conditions["singleton_coordinate_block"]["basis_indices"])
)
assert len(singleton_indices) == 6

e7_row_specs = [("singleton", idx, None) for idx in singleton_indices]
residue_payloads = [
    json.loads(R13_PATH.read_text()),
    json.loads(R56_PATH.read_text()),
    json.loads(R47_PATH.read_text()),
]
residue_count_by_file = []
for payload in residue_payloads:
    before = len(e7_row_specs)
    for component in payload["components"]:
        for row in component["non_singleton_residue_rows"]:
            e7_row_specs.append(
                (
                    component["component"],
                    int(row["residual_order"]),
                    row["entries"],
                )
            )
    residue_count_by_file.append(len(e7_row_specs) - before)

assert residue_count_by_file == [1435, 791, 255], residue_count_by_file
assert len(e7_row_specs) == 2487
E7 = matrix(finite, len(e7_row_specs), N)

for row_index, spec in enumerate(e7_row_specs):
    if spec[0] == "singleton":
        E7[row_index, int(spec[1])] = 1
        continue
    for item in spec[2]:
        E7[row_index, int(item["basis_index"])] += qq_to_finite(
            item["coefficient"], finite
        )

# ---------------------------------------------------------------------------
# 3. Regression guard: these are the current frontier invariants at p=43.
# ---------------------------------------------------------------------------

rank_h = H.rank()
rank_e8 = E8.rank()
rank_e7 = E7.rank()
EH = E8.stack(H)
rank_eh = EH.rank()
GLOBAL = EH.stack(E7)
rank_global = GLOBAL.rank()
kernel_dim = N - rank_global

print(
    "Q8TRUEGENV3|prime={}|ambient={}|"
    "H={}:{}:{}|E8={}:{}:{}|E7GEN={}:{}:{}|"
    "E8H={}:{}:{}|E8H_E7GEN={}:{}:{}".format(
        P, N,
        H.nrows(), rank_h, N-rank_h,
        E8.nrows(), rank_e8, N-rank_e8,
        E7.nrows(), rank_e7, N-rank_e7,
        EH.nrows(), rank_eh, N-rank_eh,
        GLOBAL.nrows(), rank_global, kernel_dim,
    ),
    flush=True,
)

if P == 43:
    expected = {
        "H": (1944, 1068, 532),
        "E8": (378, 378, 1222),
        "E7": (2487, 218, 1382),
        "EH": (2322, 1446, 154),
        "GLOBAL": (4809, 1582, 18),
    }
    actual = {
        "H": (H.nrows(), rank_h, N-rank_h),
        "E8": (E8.nrows(), rank_e8, N-rank_e8),
        "E7": (E7.nrows(), rank_e7, N-rank_e7),
        "EH": (EH.nrows(), rank_eh, N-rank_eh),
        "GLOBAL": (GLOBAL.nrows(), rank_global, kernel_dim),
    }
    if actual != expected:
        print("REGRESSION_MISMATCH", flush=True)
        print("expected=" + json.dumps(expected, sort_keys=True), flush=True)
        print("actual=" + json.dumps(actual, sort_keys=True), flush=True)
        raise SystemExit(2)

K = GLOBAL.right_kernel().basis_matrix()
assert K.nrows() == kernel_dim and K.ncols() == N

kernel_payload = {
    "schema": "elkies-k3.h92-q8-true1600-global-kernel-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_GLOBAL_GENERIC_E7_KERNEL",
    "prime": P,
    "ambient_path": str(AMBIENT_PATH.relative_to(ROOT)),
    "dimensions": {
        "ambient": int(N),
        "H_rows": int(H.nrows()),
        "E8_rows": int(E8.nrows()),
        "E7_generic_rows": int(E7.nrows()),
        "global_rank": int(rank_global),
        "kernel": int(kernel_dim),
    },
    "kernel_basis_rows": [
        [int(value) for value in row] for row in K.rows()
    ],
}
KERNEL_PATH.write_text(json.dumps(kernel_payload, indent=2) + "\n")
print(f"GLOBAL_KERNEL|path={KERNEL_PATH}|dimension={kernel_dim}", flush=True)

# ---------------------------------------------------------------------------
# 4. First genuine finite local quotient: E7_4--E7_3.
# ---------------------------------------------------------------------------

run_repo_script(
    "derive_h92_q8_e7_node_principal_clearings.sage",
    "--ambient", AMBIENT_PATH,
    "--output", CLEAR_PATH,
)
clearings = json.loads(CLEAR_PATH.read_text())
params = clearings["common_parameters"]
assert int(params["ambient_dimension"]) == 1600
assert int(params["K"]) == 18
assert int(params["T"]) == 17
print("NODECLEARINGCHECK|K=18|T=17|status=PASS", flush=True)

run_repo_script(
    "probe_h92_q8_e7_node_principal_local_normal_form_modp.sage",
    "--ambient", AMBIENT_PATH,
    "--clearings", CLEAR_PATH,
    "--chart", args.node,
    "--prime", str(P),
    "--mode", "local-normal-form",
    "--output", NODE_PATH,
)

node = json.loads(NODE_PATH.read_text())
image = node["finite_ambient_image"]
assert int(image["ambient_dimension"]) == N
node_rows = int(image["rows"])
sparse_columns = image["sparse_columns"]
assert len(sparse_columns) == N

# Build N*K^T directly from sparse node columns, avoiding materializing the
# usually much larger node quotient matrix N itself.
restricted = matrix(finite, node_rows, K.nrows(), sparse=True)
for col, sparse_column in enumerate(sparse_columns):
    if not sparse_column:
        continue
    survivor_coefficients = [K[j, col] for j in range(K.nrows())]
    if not any(survivor_coefficients):
        continue
    for row, coefficient in sparse_column:
        coefficient = finite(int(coefficient))
        if not coefficient:
            continue
        for j, survivor_value in enumerate(survivor_coefficients):
            if survivor_value:
                restricted[int(row), j] += coefficient * survivor_value

node_rank_on_survivors = restricted.rank()
remaining = K.nrows() - node_rank_on_survivors

result = {
    "schema": "elkies-k3.h92-q8-true1600-first-finite-node-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_FINITE_NODE_RESTRICTION",
    "prime": P,
    "node": args.node,
    "global_survivor_dimension": int(K.nrows()),
    "node_full_ambient_image": {
        "rows": node_rows,
        "rank": int(image["rank"]),
        "kernel_dimension": int(image["kernel_dimension"]),
    },
    "restricted_node_rank": int(node_rank_on_survivors),
    "remaining_dimension": int(remaining),
    "paths": {
        "ambient": str(AMBIENT_PATH.relative_to(ROOT)),
        "global_kernel": str(KERNEL_PATH.relative_to(ROOT)),
        "clearings": str(CLEAR_PATH.relative_to(ROOT)),
        "node_local_normal_form": str(NODE_PATH.relative_to(ROOT)),
    },
}
RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUEFINITENODE|prime={}|node={}|global_survivor={}|"
    "node_rows={}|node_ambient_rank={}|restricted_rank={}|remaining={}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_FINITE_NODE_RESTRICTION".format(
        P,
        args.node,
        K.nrows(),
        node_rows,
        int(image["rank"]),
        node_rank_on_survivors,
        remaining,
    ),
    flush=True,
)
print(f"RESULT_JSON|{RESULT_PATH}", flush=True)
