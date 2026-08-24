#!/usr/bin/env sage -python
"""
Exact H92 q24/orbit85 D13 -> D12 compiler using the CURRENT component-valuation
construction contract.

This intentionally supersedes the older 58 -> 10 -> point-cluster QQ probe.
The certified modular construction is:

    geometric fibre twist -8
    ambient 56
    smooth collision rank 48 -> dimension 8
    direct resolved-component valuation rank 6 -> dimension 2
    binary quartic degree 4
    D12/MW5 child

The expensive component conditions are imposed directly at the generic points
of the exceptional components, exactly as in
probe_h92_q24_d12_component_valuation_rr_modp.sage.  Mod 100003 is used only
for pivot/rank selection and final regression; all construction arithmetic and
all redundant valuation replays are exact over QQ.
"""

import json
import os
import hashlib
import gzip
import sys
import time
from pathlib import Path

from sage.all import (
    QQ, ZZ, GF, PolynomialRing, matrix, vector, gcd, lcm, identity_matrix
)

# ---------------------------------------------------------------------------
# Locate repo and replay ONLY the trusted exact prefix of the existing QQ
# compiler: exact parent/q24 input + exact I9* resolution.  Stop before the
# obsolete cluster interpretation.
# ---------------------------------------------------------------------------
def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h / "Documents/jacobian-research",
        h / "jacobian-research",
        h / "src/jacobian-research",
        h / "git/jacobian-research",
    ]
    seen = set()
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


ROOT = locate_repo()
SCRIPTS = ROOT / "elkies-k3/scripts"
LOCAL = ROOT / "artifacts/local/elkies-k3"
LEGACY = SCRIPTS / "lift_q24_d13_to_d12_resolved_rr_qq.sage"
CORE = SCRIPTS / "elliptic_neighbor_compiler.sage"

MANIFEST = LOCAL / "q24-orbit85-exact-construction-manifest.json"
COMPONENT_MOD = LOCAL / "q24-d12-component-valuation-rr-mod-100003.json"
SIGNATURE = LOCAL / "q24-orbit85-d12-signature-mod-100003.json"
MOD_RESOLUTION = LOCAL / "q24-i9star-resolution-mod-100003.json"
OUTPUT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"

CHECKPOINT_ROOT = Path(
    os.environ.get(
        "Q24DIVVAL_CHECKPOINT_DIR",
        str(LOCAL / "q24-divval-qq-checkpoints"),
    )
).expanduser().resolve()
PRIMARY_CHECKPOINT_DIR = CHECKPOINT_ROOT / "primary"
CHECKPOINT_MANIFEST = CHECKPOINT_ROOT / "checkpoint_manifest.json"
PRIMARY_CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

STOP_BEFORE_MOD_COMPONENT = os.environ.get(
    "Q24DIVVAL_STOP_BEFORE", ""
).strip()

FORCE_PRIMARY_RECOMPUTE = os.environ.get(
    "Q24DIVVAL_FORCE_PRIMARY_RECOMPUTE", ""
).strip().lower() in ("1", "true", "yes")


def atomic_write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp-{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def read_checkpoint_manifest():
    if not CHECKPOINT_MANIFEST.exists():
        return {
            "schema": "elkies-k3.q24-divval-qq-checkpoint-manifest.v1",
            "primary": {},
            "components": {},
        }
    try:
        data = json.loads(CHECKPOINT_MANIFEST.read_text())
    except Exception:
        return {
            "schema": "elkies-k3.q24-divval-qq-checkpoint-manifest.v1",
            "primary": {},
            "components": {},
        }
    data.setdefault("primary", {})
    data.setdefault("components", {})
    return data


def update_checkpoint_manifest(section, key, payload):
    data = read_checkpoint_manifest()
    data.setdefault(section, {})[str(key)] = payload
    atomic_write_json(CHECKPOINT_MANIFEST, data)


def primary_checkpoint_path(record_label, kind, threshold, branch_index, factor):
    digest = hashlib.sha256(str(factor).encode("utf-8")).hexdigest()[:16]
    return PRIMARY_CHECKPOINT_DIR / (
        f"{record_label}__{kind}__t{int(threshold)}"
        f"__b{int(branch_index)}__{digest}.json"
    )


def load_primary_checkpoint(
    path, record_label, kind, threshold, branch_index, factor
):
    if FORCE_PRIMARY_RECOMPUTE or not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except Exception as exc:
        print(
            "Q24DIVVALQQ_CHECKPOINT|"
            f"kind=PRIMARY|file={path.name}|status=INVALID_JSON|"
            f"reason={type(exc).__name__}:{exc}",
            flush=True,
        )
        return None

    expected = {
        "status": "PASS_EXACT_PRIMARY_GROEBNER_BASIS",
        "exact_center": str(record_label),
        "chart": str(kind),
        "threshold": int(threshold),
        "branch_index": int(branch_index),
        "factor": str(factor),
    }
    for field, value in expected.items():
        if data.get(field) != value:
            return None

    raw = data.get("groebner_basis")
    if not isinstance(raw, list) or not raw:
        return None

    gb = [S(str(poly)) for poly in raw]
    print(
        "Q24DIVVALQQ_CHECKPOINT|"
        f"kind=PRIMARY|center={record_label}|chart={kind}|"
        f"threshold={threshold}|branch={branch_index}|"
        f"basis={len(gb)}|status=HIT",
        flush=True,
    )
    return gb



def discovery_state_path(component):
    return CHECKPOINT_ROOT / f"state_after_{component}.json.gz"


def save_discovery_state(
    component,
    current_exact,
    raw_surface,
    raw_basis,
    current_transform,
    component_ledger,
):
    path = discovery_state_path(component)
    payload = {
        "schema": "elkies-k3.q24-divval-qq-discovery-state.v2",
        "status": "PASS_EXACT_DISCOVERY_STATE",
        "after_component": str(component),
        "current_exact": str(current_exact),
        "raw_surface": str(raw_surface),
        "raw_basis": [str(poly) for poly in raw_basis],
        "current_transform": [
            [str(v) for v in row] for row in current_transform.rows()
        ],
        "component_ledger": component_ledger,
    }

    tmp = path.with_name(path.name + f".tmp-{os.getpid()}")
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, sort_keys=True)
    tmp.replace(path)

    atomic_write_json(
        CHECKPOINT_ROOT / "latest_state.json",
        {
            "schema": "elkies-k3.q24-divval-qq-latest-state.v1",
            "status": "PASS",
            "after_component": str(component),
            "file": str(path),
        },
    )

    print(
        "Q24DIVVALQQ_CHECKPOINT|"
        f"kind=DISCOVERY_STATE|component={component}|"
        f"file={path.name}|basis={len(raw_basis)}|status=SAVED",
        flush=True,
    )


def load_latest_discovery_state(discovery_labels):
    pointer = CHECKPOINT_ROOT / "latest_state.json"
    if not pointer.exists():
        return None

    try:
        info = json.loads(pointer.read_text())
        path = Path(info["file"])
        if not path.exists():
            return None

        with gzip.open(path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)

        if data.get("status") != "PASS_EXACT_DISCOVERY_STATE":
            return None

        after = str(data["after_component"])
        if after not in discovery_labels:
            return None

        transform = matrix(
            QQ,
            [[QQ(v) for v in row] for row in data["current_transform"]],
        )
        surface = S(str(data["raw_surface"]))
        basis = [S(str(poly)) for poly in data["raw_basis"]]

        expected_rows = 8 - (discovery_labels.index(after) + 1)
        if transform.dimensions() != (expected_rows, 8):
            raise ArithmeticError(
                f"saved transform has dimensions {transform.dimensions()}, "
                f"expected {expected_rows}x8"
            )

        return {
            "after_component": after,
            "index": discovery_labels.index(after) + 1,
            "current_exact": str(data["current_exact"]),
            "raw_surface": surface,
            "raw_basis": basis,
            "current_transform": transform,
            "component_ledger": list(data.get("component_ledger", [])),
            "file": str(path),
        }
    except Exception as exc:
        print(
            "Q24DIVVALQQ_RESUME|"
            f"status=STATE_LOAD_FAILED|reason={type(exc).__name__}:{exc}",
            flush=True,
        )
        return None



def save_primary_checkpoint(
    path, record_label, kind, threshold, branch_index,
    factor, scheme_multiplicity, gb
):
    payload = {
        "schema": "elkies-k3.q24-divval-qq-primary-gb.v1",
        "status": "PASS_EXACT_PRIMARY_GROEBNER_BASIS",
        "exact_center": str(record_label),
        "chart": str(kind),
        "threshold": int(threshold),
        "branch_index": int(branch_index),
        "factor": str(factor),
        "scheme_multiplicity": int(scheme_multiplicity),
        "ring": "QQ[u,x,y]",
        "monomial_order": "degrevlex",
        "groebner_basis": [str(poly) for poly in gb],
    }
    atomic_write_json(path, payload)
    key = (
        f"{record_label}:{kind}:t{int(threshold)}:"
        f"b{int(branch_index)}:{path.stem.split('__')[-1]}"
    )
    update_checkpoint_manifest(
        "primary",
        key,
        {
            "file": str(path),
            "exact_center": str(record_label),
            "chart": str(kind),
            "threshold": int(threshold),
            "branch_index": int(branch_index),
            "groebner_basis_length": len(gb),
            "status": "PASS",
        },
    )
    print(
        "Q24DIVVALQQ_CHECKPOINT|"
        f"kind=PRIMARY|center={record_label}|chart={kind}|"
        f"threshold={threshold}|branch={branch_index}|"
        f"file={path.name}|basis={len(gb)}|status=SAVED",
        flush=True,
    )


for path in (LEGACY, CORE, MANIFEST, COMPONENT_MOD, SIGNATURE, MOD_RESOLUTION):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

started = time.monotonic()

print(
    "Q24DIVVALQQ_CHECKPOINT_INIT|"
    f"root={CHECKPOINT_ROOT}|"
    f"stop_before={STOP_BEFORE_MOD_COMPONENT or 'NONE'}|"
    f"force_recompute={int(FORCE_PRIMARY_RECOMPUTE)}|status=READY",
    flush=True,
)

def log(stage, **kw):
    extra = "|".join(f"{k}={v}" for k, v in kw.items())
    print(
        f"Q24DIVVALQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{extra}" if extra else ""),
        flush=True,
    )


legacy_text = LEGACY.read_text()
cut_marker = 'records = {r["label"]: r for r in center_records}'
cut = legacy_text.find(cut_marker)
if cut < 0:
    raise SystemExit(
        "Could not locate trusted-prefix boundary in legacy QQ compiler"
    )
prefix = legacy_text[:cut]

scope = {"__name__": "__q24_divval_exact_prefix__"}
saved_argv = list(sys.argv)
try:
    sys.argv = [str(LEGACY)]
    exec(compile(prefix, str(LEGACY), "exec"), scope)
finally:
    sys.argv = saved_argv

need = (
    "A", "B", "Z", "X", "Y", "R", "U", "K",
    "S", "u", "x", "y", "alpha", "surface0",
    "center_records", "chart_substitutions",
)
missing = [name for name in need if name not in scope]
if missing:
    raise SystemExit("exact prefix missing: " + ",".join(missing))

A = scope["A"]; B = scope["B"]
Z = scope["Z"]; X = scope["X"]; Y = scope["Y"]
R = scope["R"]; U = scope["U"]; K = scope["K"]
S = scope["S"]; u = scope["u"]; x = scope["x"]; y = scope["y"]
alpha = QQ(scope["alpha"])
surface0 = S(scope["surface0"])
center_records = list(scope["center_records"])
chart_substitutions = scope["chart_substitutions"]

records = {str(r["label"]): r for r in center_records}
exact_by_path = {str(r["path"]): r for r in center_records}

core = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), core)
squarefree_binary_quartic = core["squarefree_binary_quartic"]
binary_quartic_invariants = core["binary_quartic_invariants"]
classify_finite_short_weierstrass_fibres = core[
    "classify_finite_short_weierstrass_fibres"
]
kodaira_data_from_short_orders = core["kodaira_data_from_short_orders"]
pencil_chord_solution = core["pencil_chord_solution"]
chord_discriminant = core["chord_discriminant"]

manifest = json.loads(MANIFEST.read_text())
component_mod = json.loads(COMPONENT_MOD.read_text())
signature = json.loads(SIGNATURE.read_text())

assert manifest["status"] == "PASS_Q24_ORBIT85_EXACT_CONSTRUCTION_MANIFEST"
contract = manifest["resolved_rr_contract"]
assert contract["global_rr"] == {
    "ambient_dimension": 56,
    "smooth_collision_rank": 48,
    "post_collision_dimension": 8,
}
assert contract["resolved_rank"] == 6
assert contract["kernel_dimension"] == 2
assert component_mod["status"] == "CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP"
assert signature["rr"] == {
    "ambient": 56,
    "collision_rank": 48,
    "geometric_fibre_twist": -8,
    "kernel": 2,
    "post_collision": 8,
    "resolved_rank": 6,
}

log("CONTRACT", ambient=56, collision=48, post=8, resolved=6, kernel=2, status="PASS")

# ---------------------------------------------------------------------------
# Exact QQ resolution tree -> modular geometric labels.
# ---------------------------------------------------------------------------
PIVOT_PRIME = ZZ(100003)
Fp = GF(PIVOT_PRIME)
mod_resolution = json.loads(MOD_RESOLUTION.read_text())

def red_q(q):
    q = QQ(q)
    d = ZZ(q.denominator())
    if d % PIVOT_PRIME == 0:
        raise ZeroDivisionError(
            f"denominator divisible by pivot prime {PIVOT_PRIME}"
        )
    return Fp(ZZ(q.numerator())) / Fp(d)

def canon_mod_direction(direction):
    vals = [red_q(v) for v in direction]
    pivot = next((i for i, v in enumerate(vals) if v), None)
    if pivot is None:
        raise ArithmeticError("zero projective direction")
    inv = vals[pivot] ** -1
    return tuple(int(v * inv) for v in vals)

mod_centers = list(mod_resolution["centers"])
mod_by_path = {str(r["path"]): r for r in mod_centers}
exact_root = exact_by_path["root"]
mod_root = mod_by_path["root"]

exact_to_mod = {str(exact_root["label"]): str(mod_root["label"])}
queue = [(exact_root, mod_root)]

while queue:
    ep, mp = queue.pop(0)
    mod_edges = {}
    for edge in mp.get("children", []):
        kind = str(edge["selected_chart"])
        direction = tuple(map(int, edge["direction"]))
        child_path = (
            f"{mp['path']}/{mp['label']}:{kind}:"
            + ",".join(map(str, direction))
        )
        mod_edges[(kind, direction)] = mod_by_path[child_path]

    for edge in ep.get("children", []):
        kind = str(edge["selected_chart"])
        direction = canon_mod_direction(edge["direction"])
        key = (kind, direction)
        if key not in mod_edges:
            raise ArithmeticError(
                f"no modular edge for exact {ep['label']}:{key}"
            )
        ec = exact_by_path[str(edge["path"])]
        mc = mod_edges[key]
        exact_to_mod[str(ec["label"])] = str(mc["label"])
        queue.append((ec, mc))

if len(exact_to_mod) != len(center_records):
    raise ArithmeticError(
        f"incomplete exact/modular centre map "
        f"{len(exact_to_mod)}/{len(center_records)}"
    )

mod_to_exact = {v: k for k, v in exact_to_mod.items()}
if len(mod_to_exact) != len(exact_to_mod):
    raise ArithmeticError("noninjective exact/modular centre map")

log(
    "GEOMETRY_MAP",
    mapping=",".join(f"{m}->{mod_to_exact[m]}" for m in sorted(mod_to_exact)),
    status="PASS",
)

# ---------------------------------------------------------------------------
# Small exact linear algebra helpers.
# ---------------------------------------------------------------------------
def primitive_integer_basis(M):
    rows = []
    for row in M.rows():
        den = ZZ.one()
        for value in row:
            den = lcm(den, ZZ(QQ(value).denominator()))
        ints = [ZZ(QQ(value) * den) for value in row]
        content = ZZ.zero()
        for value in ints:
            content = gcd(content, abs(value))
        if content > 1:
            ints = [value // content for value in ints]
        pivot = next((v for v in ints if v), ZZ.zero())
        if pivot < 0:
            ints = [-v for v in ints]
        rows.append([QQ(v) for v in ints])
    return matrix(QQ, rows)

def fast_right_kernel_QQ(M, label):
    Mmod = matrix(
        Fp, M.nrows(), M.ncols(),
        [red_q(v) for v in M.list()]
    )
    rank = int(Mmod.rank())
    if rank == 0:
        return 0, matrix.identity(QQ, M.ncols())

    pivot_cols = tuple(map(int, Mmod.pivots()))
    Pmod = Mmod.matrix_from_columns(pivot_cols)
    pivot_rows = tuple(map(int, Pmod.transpose().pivots()))
    free_cols = tuple(i for i in range(M.ncols()) if i not in pivot_cols)

    P = M.matrix_from_rows_and_columns(pivot_rows, pivot_cols)
    Ffree = M.matrix_from_rows_and_columns(pivot_rows, free_cols)
    values = P.solve_right(-Ffree) if free_cols else matrix(QQ, rank, 0)

    rows = []
    for j, free_col in enumerate(free_cols):
        row = [QQ.zero()] * M.ncols()
        row[free_col] = QQ.one()
        for i, pivot_col in enumerate(pivot_cols):
            row[pivot_col] = values[i, j]
        rows.append(row)

    Kmat = primitive_integer_basis(matrix(QQ, rows))
    if M * Kmat.transpose() != matrix(QQ, M.nrows(), Kmat.nrows()):
        raise ArithmeticError(f"{label}: exact kernel replay failed")

    log(
        "EXACT_KERNEL",
        label=label,
        rows=M.nrows(),
        cols=M.ncols(),
        rank=rank,
        kernel=Kmat.nrows(),
        status="PASS",
    )
    return rank, Kmat

# ---------------------------------------------------------------------------
# 56 -> 8 smooth collision, using the reduced B-coordinate formulation.
#
# Ambient:
#   deg A <= 40 : 41 coefficients
#   deg B <= 14 : 15 coefficients
# total 56.
#
# Collision says A*X == B*Y mod Z^2.  Since X is a unit mod Z^2,
# A == B*Y/X mod Z^2.  The only conditions are that coefficients U^41..U^47
# vanish: a 7x15 matrix, not a 48x56 QQ solve.
# ---------------------------------------------------------------------------
modulus = Z**2
assert modulus.degree() == 48 and modulus.leading_coefficient() == 1

cache_path = LOCAL / "q24-xinv-mod-z2-qq.json"
Xinv = None
if cache_path.exists():
    try:
        data = json.loads(cache_path.read_text())
        for key in (
            "coefficients_low_to_high",
            "Xinv_coefficients_low_to_high",
            "xinv_coefficients_low_to_high",
        ):
            values = data.get(key)
            if isinstance(values, list) and len(values) <= 48:
                candidate = R([QQ(v) for v in values])
                if (X * candidate) % modulus == 1:
                    Xinv = candidate
                    break
    except Exception:
        Xinv = None

if Xinv is None:
    log("COLLISION_INVERSE", method="MOD_Z_NEWTON_LIFT", status="BEGIN")
    invZ = X.inverse_mod(Z)
    Xinv = (invZ * (2 - X * invZ)) % modulus
    if (X * Xinv) % modulus != 1:
        raise ArithmeticError("X inverse modulo Z^2 failed")
    cache_path.write_text(json.dumps({
        "schema": "elkies-k3.q24-xinv-mod-z2-qq.v2",
        "coefficients_low_to_high": [str(v) for v in Xinv.list()],
    }, indent=2, sort_keys=True) + "\n")
else:
    log("COLLISION_INVERSE", method="CACHE", status="PASS")

YXinv = (Y * Xinv) % modulus
A_columns = [YXinv]
for i in range(14):
    q = R(U * A_columns[-1])
    if q.degree() == 48:
        q -= q[48] * modulus
    if q.degree() >= 48:
        q %= modulus
    A_columns.append(R(q))

H = matrix(
    QQ, 7, 15,
    lambda row, col: QQ(A_columns[col][41 + row])
)
rank7, K8 = fast_right_kernel_QQ(H, "SMOOTH_COLLISION_56")
if rank7 != 7 or K8.nrows() != 8:
    raise ArithmeticError(
        f"56D collision gave rank/kernel {rank7}/{K8.nrows()}, expected 7/8"
    )

post_pairs = []
ambient_rows = []
for row in K8.rows():
    BB = R(sum(row[i] * U**i for i in range(15)))
    AA = R(sum(row[i] * A_columns[i] for i in range(15)))
    if AA.degree() > 40 or BB.degree() > 14:
        raise ArithmeticError(
            f"geometric ambient bounds failed: degA={AA.degree()} degB={BB.degree()}"
        )
    if (AA * X - BB * Y) % modulus:
        raise ArithmeticError("exact smooth collision replay failed")
    post_pairs.append((AA, BB))
    ambient_rows.append(
        [QQ(AA[i]) for i in range(41)]
        + [QQ(BB[i]) for i in range(15)]
    )

plane8_56 = matrix(QQ, ambient_rows)
assert plane8_56.dimensions() == (8, 56)

log(
    "GEOM_COLLISION",
    ambient=56,
    full_rank=48,
    reduced_rank=7,
    post=8,
    status="PASS",
)

# ---------------------------------------------------------------------------
# Build the eight exact local numerator functions efficiently.
#
# Elementary B_i=U^i columns satisfy
#   A_{i+1}=U*A_i-c_i*Z^2.
# Every numerator is linear in x,y, so do the expensive work in QQ[w].
# ---------------------------------------------------------------------------
LU = PolynomialRing(QQ, "w")
w = LU.gen()

def shift_to_LU(poly):
    result = LU.zero()
    aw = LU(alpha) + w
    for c in reversed(R(poly).list()):
        result = result * aw + QQ(c)
    return LU(result)

def LU_to_S(poly):
    poly = LU(poly)
    return S(sum(
        QQ(c) * u**i for i, c in enumerate(poly.list()) if c
    ))

log("LOCAL_NUMERATORS", method="UNIVARIATE_TRIPLE_RECURRENCE", status="BEGIN")

z = shift_to_LU(Z)
xs = shift_to_LU(X)
ys = shift_to_LU(Y)
a0 = shift_to_LU(A_columns[0])

z2 = z*z
z3 = z2*z
z4 = z3*z
z5 = z4*z
xz = xs*z
yz = ys*z

N_columns = [(
    yz - a0*xz,   # constant
    a0*z3,        # x
    z4,            # y
)]
correction = (
    -xs*z3,
    z5,
    LU.zero(),
)
aw = LU(alpha) + w

for i in range(14):
    c = QQ(A_columns[i][47])
    C, D, E = N_columns[-1]
    N_columns.append((
        aw*C - c*correction[0],
        aw*D - c*correction[1],
        aw*E,
    ))

local_numerators = []
for row in K8.rows():
    C = LU.zero(); D = LU.zero(); E = LU.zero()
    for i, q in enumerate(row):
        if q:
            C += QQ(q) * N_columns[i][0]
            D += QQ(q) * N_columns[i][1]
            E += QQ(q) * N_columns[i][2]
    local_numerators.append(
        S(LU_to_S(C) + x*LU_to_S(D) + y*LU_to_S(E))
    )

assert len(local_numerators) == 8
log("LOCAL_NUMERATORS", basis=8, status="PASS")

# ---------------------------------------------------------------------------
# Direct resolved-component valuations.
#
# IMPORTANT: basis functions are TOTAL pullbacks; do not divide them by
# exceptional factors.  Only the surface is strict-transformed.
# ---------------------------------------------------------------------------
def divide_power_raw(poly, exceptional, power):
    q = S(poly)
    for unused in range(int(power)):
        q, rem = q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("required exceptional surface power does not divide")
    return S(q)

parent_link = {}
for parent in center_records:
    for edge in parent.get("children", []):
        child_record = exact_by_path[str(edge["path"])]
        parent_link[str(child_record["label"])] = (
            str(parent["label"]),
            str(edge["selected_chart"]),
        )

root_exact = str(exact_root["label"])
state_cache = {root_exact: (surface0, list(local_numerators))}

def state_before(exact_label):
    if exact_label in state_cache:
        return state_cache[exact_label]
    parent_label, kind = parent_link[exact_label]
    parent_surface, parent_basis = state_before(parent_label)
    parent_record = records[parent_label]
    point = tuple(QQ(v) for v in parent_record["point"])
    subs, e = chart_substitutions(point, kind)

    child_surface = divide_power_raw(
        S(parent_surface(*subs)),
        e,
        int(parent_record["multiplicity"]),
    )
    child_basis = [S(poly(*subs)) for poly in parent_basis]
    state_cache[exact_label] = (child_surface, child_basis)
    return state_cache[exact_label]

thresholds = {
    str(k): int(v)
    for k, v in contract["centre_thresholds"].items()
}
expected_thresholds = {
    "C01":2, "C02":4, "C03":3, "C04":6,
    "C05":5, "C06":8, "C07":7, "C08":10,
    "C09":9, "C10":6, "C11":11, "C12":1,
}
if thresholds != expected_thresholds:
    raise ArithmeticError(f"unexpected manifest thresholds {thresholds}")

mod_ledger = {
    str(row["component"]): row
    for row in component_mod["resolved_cluster"]["condition_ledger"]
}

PRIMARY_GB_CACHE = {}

def quotient_rows_for_chart(
    surface_before, basis_before, record, kind, threshold
):
    point = tuple(QQ(v) for v in record["point"])
    subs, e = chart_substitutions(point, kind)

    strict = divide_power_raw(
        S(surface_before(*subs)),
        e,
        int(record["multiplicity"]),
    )
    exceptional_restriction = S(strict.subs({e: 0}))

    # Do not pull huge section polynomials into a chart that misses E.
    if exceptional_restriction.is_constant() and exceptional_restriction != 0:
        return None

    factors = [
        (S(factor), int(mult))
        for factor, mult in exceptional_restriction.factor()
    ]
    if not factors:
        raise ArithmeticError(
            f"empty exceptional factorization at {record['label']}:{kind}"
        )

    # This is the expensive operation.  Perform it only for a chart that
    # actually contains the exceptional component.
    log(
        "CHART_PULLBACK",
        center=record["label"],
        chart=kind,
        basis=len(basis_before),
        status="BEGIN",
    )
    pulled = [S(poly(*subs)) for poly in basis_before]
    log(
        "CHART_PULLBACK",
        center=record["label"],
        chart=kind,
        basis=len(basis_before),
        status="PASS",
    )

    matrices = []
    branch_info = []

    for branch_index, (h, scheme_mult) in enumerate(factors):
        P = S.ideal([e, h])
        J = S.ideal([strict]) + P**int(threshold)

        cache_key = (
            str(record["label"]),
            str(kind),
            int(threshold),
            str(h),
        )
        checkpoint_path = primary_checkpoint_path(
            record["label"], kind, threshold, branch_index, h
        )

        gb = None
        if cache_key in PRIMARY_GB_CACHE:
            gb = PRIMARY_GB_CACHE[cache_key]
            log(
                "PRIMARY",
                center=record["label"],
                chart=kind,
                branch=branch_index,
                threshold=threshold,
                method="MEMORY_CACHE",
                status="PASS",
            )
        else:
            gb = load_primary_checkpoint(
                checkpoint_path,
                record["label"],
                kind,
                threshold,
                branch_index,
                h,
            )
            if gb is not None:
                PRIMARY_GB_CACHE[cache_key] = gb
                log(
                    "PRIMARY",
                    center=record["label"],
                    chart=kind,
                    branch=branch_index,
                    threshold=threshold,
                    method="DISK_CHECKPOINT",
                    status="PASS",
                )

        if gb is None:
            log(
                "PRIMARY",
                center=record["label"],
                chart=kind,
                branch=branch_index,
                threshold=threshold,
                method="COMPUTE",
                status="BEGIN",
            )
            pieces = J.primary_decomposition()
            primary = []
            for Q in pieces:
                radical = Q.radical()
                if radical <= P and P <= radical:
                    primary.append(Q)
            if len(primary) != 1:
                raise ArithmeticError(
                    f"expected unique P-primary component at "
                    f"{record['label']}:{kind}:{h}; got {len(primary)}"
                )
            gb = primary[0].groebner_basis()
            save_primary_checkpoint(
                checkpoint_path,
                record["label"],
                kind,
                threshold,
                branch_index,
                h,
                scheme_mult,
                gb,
            )
            PRIMARY_GB_CACHE[cache_key] = gb
        remainders = [poly.reduce(gb) for poly in pulled]
        exponents = sorted({
            exp
            for rem in remainders
            for exp, coef in rem.dict().items()
            if coef
        })
        M = matrix(
            QQ,
            len(exponents),
            len(basis_before),
            lambda i, j: QQ(
                remainders[j].dict().get(exponents[i], QQ.zero())
            ),
        )
        matrices.append(M)
        branch_info.append({
            "factor": str(h),
            "rows": int(M.nrows()),
            "scheme_multiplicity": scheme_mult,
        })
        log(
            "PRIMARY",
            center=record["label"],
            chart=kind,
            branch=branch_index,
            rows=M.nrows(),
            status="PASS",
        )

    total = matrix(QQ, 0, len(basis_before))
    for M in matrices:
        if M.nrows():
            total = total.stack(M)

    return {
        "matrix": total,
        "chart": kind,
        "branches": branch_info,
        "strict": strict,
        "pulled": pulled,
    }

def mod_rank(M):
    return int(matrix(
        Fp, M.nrows(), M.ncols(),
        [red_q(v) for v in M.list()]
    ).rank())

# ---------------------------------------------------------------------------
# Two-pass exact divisorial solve.
#
# Pass 1 is only a discovery accelerator.  The modular ledger tells us which
# components increased the global rank and whether one affine chart already
# realizes the full component-cover rank.  We compute only that subset over
# QQ and require its exact reduction/rank growth to match.
#
# Pass 2 is the proof: after obtaining the exact 2-plane, replay EVERY
# component and EVERY chart over QQ on those two functions.  Therefore no
# modularly-pruned condition is trusted as an exact vanishing statement.
# ---------------------------------------------------------------------------
def ledger_new_rank(entry):
    value = entry.get("new_global_rank", entry.get("new_rank"))
    if value is None:
        raise KeyError(
            "component ledger has neither new_global_rank nor new_rank"
        )
    return int(value)


def discovery_chart_kinds(entry):
    cover = int(entry["component_cover_rank"])
    charts = list(entry.get("charts", []))

    # Prefer a single chart which already realizes the full component rank.
    full = [
        str(row["chart"])
        for row in charts
        if int(row.get("rank", -1)) == cover
    ]
    for preferred in ("u", "x", "y"):
        if preferred in full:
            return (preferred,)
    if full:
        return (full[0],)

    # Otherwise retain every chart with a nonzero modular quotient rank.
    nonzero = [
        str(row["chart"])
        for row in charts
        if int(row.get("rank", 0)) > 0
    ]
    return tuple(nonzero) if nonzero else ("u", "x", "y")


discovery_labels = [
    label
    for label in sorted(thresholds, key=lambda name: int(name[1:]))
    if ledger_new_rank(mod_ledger[label]) > 0
]

print(
    "Q24DIVVALQQ_DISCOVERY_PLAN|"
    f"components={','.join(discovery_labels)}|"
    f"count={len(discovery_labels)}|"
    "method=LAZY_8D_GEOMETRY_PLUS_SMALL_TRANSFORM|status=PASS",
    flush=True,
)

active_exact = [mod_to_exact[label] for label in discovery_labels]

def ancestor_chain_to(descendant, ancestor):
    rev = []
    cur = descendant
    while cur != ancestor:
        if cur not in parent_link:
            raise ArithmeticError(
                f"{ancestor} is not an ancestor of {descendant}"
            )
        parent, kind = parent_link[cur]
        rev.append((parent, cur, kind))
        cur = parent
    return list(reversed(rev))

def path_first_chart(parent, descendant):
    chain = ancestor_chain_to(descendant, parent)
    return chain[0][2] if chain else None

# Keep the eight geometric basis functions untouched by kernel combinations.
# Only their chart coordinates change.  The current RR subspace is represented
# by current_transform (k x 8).
raw_basis = list(local_numerators)
raw_surface = surface0
current_exact = root_exact
current_transform = identity_matrix(QQ, 8)

discovery_rows_original = matrix(QQ, 0, 8)
component_ledger = []

RESUME_START_INDEX = 0
RESUME_RANK_OFFSET = 0

full_state = load_latest_discovery_state(discovery_labels)

if full_state is not None:
    RESUME_START_INDEX = int(full_state["index"])
    RESUME_RANK_OFFSET = RESUME_START_INDEX
    current_exact = str(full_state["current_exact"])
    raw_surface = S(full_state["raw_surface"])
    raw_basis = list(full_state["raw_basis"])
    current_transform = matrix(QQ, full_state["current_transform"])
    component_ledger = list(full_state["component_ledger"])

    print(
        "Q24DIVVALQQ_RESUME|"
        f"after_component={full_state['after_component']}|"
        f"completed={RESUME_START_INDEX}|"
        f"dimension={current_transform.nrows()}|"
        f"state={full_state['file']}|"
        "status=PASS_FULL_STATE_CHECKPOINT",
        flush=True,
    )

else:
    manifest = read_checkpoint_manifest()
    saved = manifest.get("components", {})
    last_label = None

    for idx, label in enumerate(discovery_labels):
        entry = saved.get(label)
        if not isinstance(entry, dict):
            break
        if entry.get("status") != "PASS_EXACT_DISCOVERY_COMPONENT":
            break

        raw_transform = entry.get("current_transform")
        if not isinstance(raw_transform, list):
            break

        candidate = matrix(
            QQ,
            [[QQ(v) for v in row] for row in raw_transform],
        )
        expected_rows = 8 - (idx + 1)
        if candidate.dimensions() != (expected_rows, 8):
            break

        current_transform = candidate
        RESUME_START_INDEX = idx + 1
        RESUME_RANK_OFFSET = idx + 1
        last_label = label

        component_ledger.append({
            "component": label,
            "exact_center": str(entry.get("exact_center", "")),
            "threshold": int(entry.get("threshold", thresholds[label])),
            "discovery_charts": list(entry.get("discovery_charts", [])),
            "dimension_before": int(
                entry.get("dimension_before", expected_rows + 1)
            ),
            "local_rank": int(entry.get("local_rank", 1)),
            "dimension_after": int(
                entry.get("dimension_after", expected_rows)
            ),
            "cumulative_rank": int(
                entry.get("cumulative_rank", idx + 1)
            ),
            "resumed_from_v1_manifest": True,
        })

    if RESUME_START_INDEX:
        next_label = (
            discovery_labels[RESUME_START_INDEX]
            if RESUME_START_INDEX < len(discovery_labels)
            else None
        )

        print(
            "Q24DIVVALQQ_RESUME|"
            f"after_component={last_label}|"
            f"completed={RESUME_START_INDEX}|"
            f"dimension={current_transform.nrows()}|"
            f"next={next_label or 'RESOLVED'}|"
            "status=PASS_V1_COMPONENT_CHECKPOINT",
            flush=True,
        )

        if next_label is not None:
            next_exact = mod_to_exact[next_label]
            log(
                "RESUME_GEOMETRY",
                target_component=next_label,
                qq_center=next_exact,
                basis=8,
                status="BEGIN",
            )

            raw_surface, raw_basis = state_before(next_exact)
            raw_surface = S(raw_surface)
            raw_basis = list(raw_basis)
            current_exact = next_exact

            log(
                "RESUME_GEOMETRY",
                target_component=next_label,
                qq_center=next_exact,
                basis=8,
                status="PASS",
            )
    else:
        print("Q24DIVVALQQ_RESUME|status=FRESH_DISCOVERY", flush=True)


for active_index, mod_label in enumerate(discovery_labels):
    if active_index < RESUME_START_INDEX:
        print(
            "Q24DIVVALQQ_DISCOVERY|"
            f"component={mod_label}|status=SKIP_COMPONENT_CHECKPOINT",
            flush=True,
        )
        continue

    if STOP_BEFORE_MOD_COMPONENT and mod_label == STOP_BEFORE_MOD_COMPONENT:
        update_checkpoint_manifest(
            "components",
            f"STOP_BEFORE_{mod_label}",
            {
                "status": "CLEAN_STOP",
                "next_component": mod_label,
                "checkpoint_root": str(CHECKPOINT_ROOT),
            },
        )
        print(
            "Q24DIVVALQQ_CHECKPOINT_STOP|"
            f"before_component={mod_label}|"
            f"checkpoint_root={CHECKPOINT_ROOT}|"
            "status=CLEAN_EXIT",
            flush=True,
        )
        raise SystemExit(0)

    exact_label = mod_to_exact[mod_label]

    # Walk the raw 8-function basis down to the active centre.
    if current_exact != exact_label:
        chain = ancestor_chain_to(exact_label, current_exact)
        for parent_label, child_label, kind in chain:
            parent_record = records[parent_label]
            point = tuple(QQ(v) for v in parent_record["point"])
            subs, e = chart_substitutions(point, kind)

            raw_surface = divide_power_raw(
                S(raw_surface(*subs)),
                e,
                int(parent_record["multiplicity"]),
            )

            log(
                "DISCOVERY_RAW_PULLBACK",
                parent=parent_label,
                child=child_label,
                chart=kind,
                basis=8,
                current_dim=current_transform.nrows(),
                status="BEGIN",
            )
            raw_basis = [S(poly(*subs)) for poly in raw_basis]
            log(
                "DISCOVERY_RAW_PULLBACK",
                parent=parent_label,
                child=child_label,
                chart=kind,
                basis=8,
                current_dim=current_transform.nrows(),
                status="PASS",
            )
            current_exact = child_label

    record = records[exact_label]
    threshold = thresholds[mod_label]
    expected = mod_ledger[mod_label]

    next_exact = (
        active_exact[active_index + 1]
        if active_index + 1 < len(active_exact)
        else None
    )
    path_kind = (
        path_first_chart(exact_label, next_exact)
        if next_exact is not None
        else None
    )

    candidate_kinds = list(discovery_chart_kinds(expected))
    kinds = (
        (path_kind,)
        if path_kind in candidate_kinds
        else tuple(candidate_kinds)
    )

    centre_raw = matrix(QQ, 0, 8)
    results = {}
    used = []

    for kind in kinds:
        result = quotient_rows_for_chart(
            raw_surface,
            raw_basis,
            record,
            kind,
            threshold,
        )
        if result is None:
            continue
        results[kind] = result
        used.append(kind)
        if result["matrix"].nrows():
            centre_raw = centre_raw.stack(result["matrix"])

    # Restrict the raw quotient map to the current k-dimensional RR subspace.
    centre_current = centre_raw * current_transform.transpose()

    expected_new = ledger_new_rank(expected)
    local_rank_mod = mod_rank(centre_current)
    if expected_new != 1 or local_rank_mod != 1:
        raise ArithmeticError(
            f"lazy discovery rank at {mod_label} is {local_rank_mod}; "
            f"expected {expected_new}"
        )

    before_dim = current_transform.nrows()
    local_rank_exact, Klocal = fast_right_kernel_QQ(
        centre_current,
        f"DISCOVERY_{mod_label}",
    )
    if local_rank_exact != 1 or Klocal.nrows() != before_dim - 1:
        raise ArithmeticError(
            f"unexpected lazy exact kernel at {mod_label}: "
            f"rank={local_rank_exact}, kernel={Klocal.nrows()}, "
            f"before={before_dim}"
        )

    # The actual condition in the original post-collision coordinates.
    M_original = centre_current * current_transform
    if M_original.nrows():
        discovery_rows_original = discovery_rows_original.stack(M_original)

    # Pure small-matrix update.  NO giant polynomial linear combination.
    current_transform = primitive_integer_basis(
        Klocal * current_transform
    )

    new_cuts_done = active_index - RESUME_START_INDEX + 1
    cumulative_rank = RESUME_RANK_OFFSET + new_cuts_done
    exp_cumulative = int(expected["cumulative_rank"])
    print(
        "Q24DIVVALQQ_DISCOVERY|"
        f"component={mod_label}|qq_center={exact_label}|"
        f"threshold={threshold}|charts={','.join(used)}|"
        f"before={before_dim}|local_rank=1|"
        f"after={current_transform.nrows()}|"
        f"cumulative_rank={cumulative_rank}|"
        f"mod_cumulative={exp_cumulative}|"
        f"status={'PASS' if cumulative_rank==exp_cumulative else 'MISMATCH'}",
        flush=True,
    )
    if cumulative_rank != exp_cumulative:
        raise ArithmeticError(
            f"lazy cumulative mismatch at {mod_label}: "
            f"QQmod={cumulative_rank}, mod={exp_cumulative}"
        )

    component_record = {
        "component": mod_label,
        "exact_center": exact_label,
        "threshold": threshold,
        "discovery_charts": list(used),
        "dimension_before": before_dim,
        "local_rank": 1,
        "dimension_after": int(current_transform.nrows()),
        "cumulative_rank": cumulative_rank,
    }
    component_ledger.append(component_record)
    update_checkpoint_manifest(
        "components",
        mod_label,
        {
            **component_record,
            "status": "PASS_EXACT_DISCOVERY_COMPONENT",
            "current_transform": [
                [str(v) for v in row]
                for row in current_transform.rows()
            ],
        },
    )
    print(
        "Q24DIVVALQQ_CHECKPOINT|"
        f"kind=COMPONENT|component={mod_label}|"
        f"after={current_transform.nrows()}|"
        f"cumulative_rank={cumulative_rank}|status=SAVED",
        flush=True,
    )

    # If this is the continuation chart, reuse its already-computed raw
    # pullback state.  Again, no kernel combination is applied to polynomials.
    if next_exact is not None and path_kind in results:
        result = results[path_kind]
        chain = ancestor_chain_to(next_exact, exact_label)
        immediate_child = chain[0][1]
        raw_surface = result["strict"]
        raw_basis = list(result["pulled"])
        current_exact = immediate_child

    save_discovery_state(
        mod_label,
        current_exact,
        raw_surface,
        raw_basis,
        current_transform,
        component_ledger,
    )

if current_transform.dimensions() != (2, 8):
    raise ArithmeticError(
        f"lazy discovery ended with transform "
        f"{current_transform.dimensions()}, expected 2x8"
    )
# On a resumed run, discovery_rows_original contains only cuts computed
# in THIS process. Restored cuts are already encoded in current_transform.
# The full all-component QQ replay below independently re-verifies them.
expected_new_discovery_rank = 6 - RESUME_RANK_OFFSET
actual_new_discovery_rank = mod_rank(discovery_rows_original)

if actual_new_discovery_rank != expected_new_discovery_rank:
    raise ArithmeticError(
        f"resume-local discovery rank is {actual_new_discovery_rank}, "
        f"expected {expected_new_discovery_rank} "
        f"(resolved_total=6, restored={RESUME_RANK_OFFSET})"
    )

if current_transform.dimensions() != (2, 8):
    raise ArithmeticError(
        f"resolved transform has dimensions {current_transform.dimensions()}, "
        "expected 2x8"
    )

print(
    "Q24DIVVALQQ_RESUME_RANK|"
    f"restored_rank={RESUME_RANK_OFFSET}|"
    f"new_rank={actual_new_discovery_rank}|"
    f"total_rank={RESUME_RANK_OFFSET + actual_new_discovery_rank}|"
    f"kernel={current_transform.nrows()}|status=PASS",
    flush=True,
)

K2 = primitive_integer_basis(current_transform)

atomic_write_json(
    CHECKPOINT_ROOT / "resolved_2plane.json",
    {
        "schema": "elkies-k3.q24-divval-qq-resolved-plane.v1",
        "status": "PASS_EXACT_RESOLVED_2PLANE",
        "K2_post_collision_2x8": [
            [str(v) for v in row] for row in K2.rows()
        ],
        "discovery_components": list(discovery_labels),
        "resolved_rank": 6,
        "kernel_dimension": 2,
    },
)
update_checkpoint_manifest(
    "components",
    "RESOLVED_2PLANE",
    {
        "file": str(CHECKPOINT_ROOT / "resolved_2plane.json"),
        "status": "PASS",
    },
)

if discovery_rows_original.nrows():
    if discovery_rows_original * K2.transpose() != matrix(
        QQ, discovery_rows_original.nrows(), 2
    ):
        raise ArithmeticError(
            "exact 2-plane fails newly-computed discovery rows"
        )
log(
    "EXACT_KERNEL",
    label="RESOLVED_DIVISORIAL",
    rows=discovery_rows_original.nrows(),
    cols=8,
    rank=6,
    kernel=2,
    method="LAZY_SMALL_TRANSFORM",
    status="PASS",
)

# Construct only the final two local numerator functions from the elementary
# univariate numerator basis.
B2_verify = primitive_integer_basis(K2 * K8)
assert B2_verify.dimensions() == (2, 15)

final_local_numerators = []
for row in B2_verify.rows():
    C = LU.zero()
    D = LU.zero()
    E = LU.zero()
    for i, coefficient in enumerate(row):
        if coefficient:
            q = QQ(coefficient)
            C += q * N_columns[i][0]
            D += q * N_columns[i][1]
            E += q * N_columns[i][2]
    final_local_numerators.append(
        S(LU_to_S(C) + x*LU_to_S(D) + y*LU_to_S(E))
    )

assert len(final_local_numerators) == 2

final_state_cache = {
    root_exact: (surface0, list(final_local_numerators))
}

def final_state_before(exact_label):
    if exact_label in final_state_cache:
        return final_state_cache[exact_label]

    parent_label, kind = parent_link[exact_label]
    parent_surface, parent_basis = final_state_before(parent_label)
    parent_record = records[parent_label]
    point = tuple(QQ(v) for v in parent_record["point"])
    subs, e = chart_substitutions(point, kind)

    child_surface = divide_power_raw(
        S(parent_surface(*subs)),
        e,
        int(parent_record["multiplicity"]),
    )

    log(
        "FINAL_STATE_PULLBACK",
        parent=parent_label,
        child=exact_label,
        chart=kind,
        basis=2,
        status="BEGIN",
    )
    child_basis = [S(poly(*subs)) for poly in parent_basis]
    log(
        "FINAL_STATE_PULLBACK",
        parent=parent_label,
        child=exact_label,
        chart=kind,
        basis=2,
        status="PASS",
    )

    final_state_cache[exact_label] = (
        child_surface,
        child_basis,
    )
    return final_state_cache[exact_label]


# Full exact proof replay, optimized at the DIVISORIAL level.
#
# A divisorial valuation is a condition at the generic point of an
# irreducible exceptional component.  It is therefore sufficient to certify
# membership on one affine chart containing a nonempty open subset of that
# component.  Rechecking every overlapping chart is redundant.
#
# Safety rules:
#   * the chosen exact chart must meet the exceptional divisor;
#   * it must see every irreducible branch created at this blow-up centre;
#   * the final two functions must reduce to zero in the exact P-primary
#     quotient;
#   * if no single chart sees all branches, fall back to the old all-chart
#     audit for that component.
#
# Each successful generic-component audit is checkpointed against a
# fingerprint of the exact final 2x8 post-collision plane.

def exact_plane_fingerprint(Kmat):
    text = "\\n".join(
        ",".join(str(QQ(v)) for v in row)
        for row in Kmat.rows()
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


FINAL_PLANE_FINGERPRINT = exact_plane_fingerprint(K2)

discovery_chart_hint = {}
for row in component_ledger:
    label = str(row.get("component", ""))
    hints = list(row.get("discovery_charts", []))
    if label and hints:
        discovery_chart_hint[label] = [str(v) for v in hints]


def generic_audit_chart_order(mod_label):
    expected = mod_ledger[mod_label]
    cover_rank = int(expected.get("component_cover_rank", 0))
    chart_rows = list(expected.get("charts", []))
    order_index = {"u": 0, "x": 1, "y": 2}

    result = []

    # First preference: the exact discovery chart.  For active components
    # this normally gives us an immediate disk-checkpoint hit.
    for kind in discovery_chart_hint.get(mod_label, []):
        if kind not in result:
            result.append(kind)

    # Next: modular charts which individually realize the whole component
    # cover rank, then remaining charts by decreasing rank.
    ranked = sorted(
        chart_rows,
        key=lambda row: (
            0 if int(row.get("rank", -1)) == cover_rank else 1,
            -int(row.get("rank", 0)),
            order_index.get(str(row.get("chart")), 99),
        ),
    )
    for row in ranked:
        kind = str(row.get("chart"))
        if kind in ("u", "x", "y") and kind not in result:
            result.append(kind)

    for kind in ("u", "x", "y"):
        if kind not in result:
            result.append(kind)

    return tuple(result)


def audit_checkpoint_valid(mod_label):
    manifest = read_checkpoint_manifest()
    entry = manifest.get("components", {}).get(f"AUDIT_{mod_label}")
    return (
        isinstance(entry, dict)
        and entry.get("status") == "PASS_EXACT_GENERIC_COMPONENT"
        and entry.get("plane_fingerprint") == FINAL_PLANE_FINGERPRINT
        and int(entry.get("threshold", -1)) == int(thresholds[mod_label])
    )


print(
    "Q24DIVVALQQ_AUDIT_PLAN|"
    "method=ONE_GENERIC_CHART_PER_IRREDUCIBLE_COMPONENT|"
    f"plane={FINAL_PLANE_FINGERPRINT[:16]}|"
    "fallback=ALL_CHARTS_IF_BRANCH_COVER_INCOMPLETE|status=PASS",
    flush=True,
)

verified_charts = 0
verified_rows = 0
verified_components = 0
audit_fallbacks = 0

for mod_label in sorted(
    thresholds, key=lambda name: int(name[1:])
):
    exact_label = mod_to_exact[mod_label]
    record = records[exact_label]
    threshold = thresholds[mod_label]

    if audit_checkpoint_valid(mod_label):
        manifest = read_checkpoint_manifest()
        entry = manifest["components"][f"AUDIT_{mod_label}"]
        verified_components += 1
        verified_charts += int(entry.get("charts_checked", 1))
        verified_rows += int(entry.get("rows_checked", 0))
        print(
            "Q24DIVVALQQ_GENERIC_AUDIT|"
            f"component={mod_label}|qq_center={exact_label}|"
            f"threshold={threshold}|chart={entry.get('chart','?')}|"
            f"branches={entry.get('branches','?')}|"
            "status=SKIP_AUDIT_CHECKPOINT",
            flush=True,
        )
        continue

    surface_before, basis_before = final_state_before(exact_label)
    expected_branches = int(
        record.get("geometric_exceptional_components", 1)
    )

    chosen = None
    incomplete = []

    for kind in generic_audit_chart_order(mod_label):
        result = quotient_rows_for_chart(
            surface_before,
            basis_before,
            record,
            kind,
            threshold,
        )
        if result is None:
            continue

        branch_count = len(result.get("branches", []))

        if branch_count < expected_branches:
            incomplete.append((kind, branch_count))
            print(
                "Q24DIVVALQQ_GENERIC_AUDIT_CHART|"
                f"component={mod_label}|qq_center={exact_label}|"
                f"chart={kind}|branches={branch_count}|"
                f"expected_branches={expected_branches}|"
                "status=INCOMPLETE_BRANCH_COVER",
                flush=True,
            )
            continue

        M2 = result["matrix"]
        if M2 and any(M2.list()):
            raise ArithmeticError(
                f"final exact 2-plane violates generic component valuation "
                f"{mod_label}:{kind}; rank={mod_rank(M2)}"
            )

        chosen = (kind, branch_count, int(M2.nrows()))
        break

    if chosen is None:
        # Rare safety fallback, principally for a split exceptional divisor
        # whose irreducible branches are not simultaneously visible in one
        # affine chart.  This restores the old stronger audit for this one
        # component only.
        audit_fallbacks += 1
        rows_here = 0
        charts_here = 0
        max_branches = 0

        for kind in ("u", "x", "y"):
            result = quotient_rows_for_chart(
                surface_before,
                basis_before,
                record,
                kind,
                threshold,
            )
            if result is None:
                continue

            M2 = result["matrix"]
            if M2 and any(M2.list()):
                raise ArithmeticError(
                    f"final exact 2-plane violates fallback component "
                    f"valuation {mod_label}:{kind}; rank={mod_rank(M2)}"
                )

            charts_here += 1
            rows_here += int(M2.nrows())
            max_branches = max(
                max_branches,
                len(result.get("branches", [])),
            )

        if charts_here == 0 or max_branches < expected_branches:
            raise ArithmeticError(
                f"could not certify all irreducible branches of "
                f"{mod_label}->{exact_label}; expected={expected_branches}, "
                f"max_seen={max_branches}, incomplete={incomplete}"
            )

        chosen = (
            "ALL",
            max_branches,
            rows_here,
        )
        verified_charts += charts_here
        verified_rows += rows_here

    else:
        verified_charts += 1
        verified_rows += chosen[2]

    verified_components += 1

    audit_record = {
        "status": "PASS_EXACT_GENERIC_COMPONENT",
        "component": mod_label,
        "exact_center": exact_label,
        "threshold": int(threshold),
        "chart": chosen[0],
        "branches": int(chosen[1]),
        "expected_branches": int(expected_branches),
        "rows_checked": int(chosen[2]),
        "charts_checked": 1 if chosen[0] != "ALL" else 3,
        "plane_fingerprint": FINAL_PLANE_FINGERPRINT,
    }
    update_checkpoint_manifest(
        "components",
        f"AUDIT_{mod_label}",
        audit_record,
    )

    print(
        "Q24DIVVALQQ_GENERIC_AUDIT|"
        f"component={mod_label}|qq_center={exact_label}|"
        f"threshold={threshold}|chart={chosen[0]}|"
        f"branches={chosen[1]}|rows={chosen[2]}|"
        f"expected_branches={expected_branches}|"
        "status=PASS_EXACT_GENERIC_COMPONENT",
        flush=True,
    )

if verified_components != len(thresholds):
    raise ArithmeticError(
        f"generic component audit certified {verified_components}/"
        f"{len(thresholds)} components"
    )

print(
    "Q24DIVVALQQ_GENERIC_AUDIT_SUMMARY|"
    f"components={verified_components}|charts={verified_charts}|"
    f"rows={verified_rows}|fallbacks={audit_fallbacks}|"
    f"plane={FINAL_PLANE_FINGERPRINT[:16]}|"
    "status=PASS_EXACT_ALL_COMPONENTS",
    flush=True,
)

print(
    "Q24DIVVALQQ_REDUNDANT_REPLAY|"
    f"components={len(thresholds)}|charts={verified_charts}|"
    f"rows={verified_rows}|basis=2|status=PASS_EXACT_ALL_COMPONENTS",
    flush=True,
)

log(
    "DIVISORIAL_RR",
    post=8,
    resolved_rank=6,
    kernel=2,
    redundant_rows=verified_rows,
    exact_replay=1,
    discovery_components=len(discovery_labels),
    status="PASS",
)

# ---------------------------------------------------------------------------
# Recover exact 2x56 ambient plane and regress it mod 100003.
# ---------------------------------------------------------------------------
B2 = primitive_integer_basis(K2 * K8)
assert B2.dimensions() == (2, 15)

def ambient56_from_B2(Bmat):
    rows = []
    pairs = []
    for row in Bmat.rows():
        BB = R(sum(row[i] * U**i for i in range(15)))
        AA = R(sum(row[i] * A_columns[i] for i in range(15)))
        if AA.degree() > 40 or BB.degree() > 14:
            raise ArithmeticError(
                "final exact plane violates 56D ambient bounds"
            )
        if (AA*X - BB*Y) % modulus:
            raise ArithmeticError("final exact plane fails collision")
        rows.append(
            [QQ(AA[i]) for i in range(41)]
            + [QQ(BB[i]) for i in range(15)]
        )
        pairs.append((AA, BB, K(AA)/K(Z**2), K(BB)/K(Z)))
    return matrix(QQ, rows), pairs

# The modular signature fixes RREF pivots (0,1).  Canonicalize the exact
# basis with the same two columns BEFORE forming the pencil.  This removes
# arbitrary GL_2(Q) basis changes from quartic/Jacobian regression.
raw56, unused_pairs = ambient56_from_B2(B2)
signature_pivots = tuple(map(int, signature["plane_pivots"]))
if signature_pivots != (0, 1):
    raise ArithmeticError(
        f"unexpected modular plane pivots {signature_pivots}"
    )

pivot_minor = raw56.matrix_from_columns(signature_pivots)
if pivot_minor.det() == 0:
    raise ArithmeticError(
        "exact final plane has singular modular-signature pivot minor"
    )

basis_change = pivot_minor.inverse()
B2 = basis_change * B2
final56, final_pairs = ambient56_from_B2(B2)

# Exact canonical pivot replay.
if final56.matrix_from_columns(signature_pivots) != identity_matrix(QQ, 2):
    raise ArithmeticError("exact final plane canonicalization failed")

final56_mod = matrix(
    Fp, 2, 56, [red_q(v) for v in final56.list()]
)
sig_plane = matrix(Fp, signature["plane_rref_2x56"])
if final56_mod != sig_plane:
    raise ArithmeticError(
        "canonical exact 2x56 plane does not reduce to modular signature"
    )

log(
    "PLANE_REGRESSION",
    plane="2x56",
    pivots="0,1",
    canonical_basis=1,
    prime=100003,
    status="PASS",
)

# ---------------------------------------------------------------------------
# Degree-two chord pencil -> exact quartic -> exact D12 Jacobian.
# ---------------------------------------------------------------------------
VR = PolynomialRing(QQ, "V")
V = VR.gen()
VF = VR.fraction_field()
UR = PolynomialRing(VF, "U")
UK = UR.fraction_field()

def lift_poly(poly):
    poly = R(poly)
    return UR([VF(c) for c in poly.list()])

def lift_rf(value):
    value = K(value)
    return (
        UK(lift_poly(R(value.numerator())))
        / UK(lift_poly(R(value.denominator())))
    )

a0, b0 = lift_rf(final_pairs[0][2]), lift_rf(final_pairs[0][3])
a1, b1 = lift_rf(final_pairs[1][2]), lift_rf(final_pairs[1][3])
xPV, yPV = lift_rf(K(X)/K(Z**2)), lift_rf(K(Y)/K(Z**3))
AV = lift_poly(A)

# -----------------------------------------------------------------------
# Fast exact chord discriminant.
#
# m = N/(Z*D), where
#   N = A1 - V*A0
#   D = V*B0 - B1.
#
# For marked point (-P)=(X/Z^2,-Y/Z^3),
#
# disc =
# [N^4 - 6 X N^2 D^2 - 8 Y N D^3
#  - 3 X^2 D^4 - 4 A Z^4 D^4] / (Z^4 D^4).
#
# The denominator is (Z*D)^4, hence already a square.  Extract the odd
# multiplicity part of the numerator with Yun's gcd-only squarefree
# decomposition rather than irreducibly factoring numerator+denominator.
# -----------------------------------------------------------------------
QUARTIC_CHECKPOINT = CHECKPOINT_ROOT / "quartic_exact_yun.json.gz"

def exact_poly_quo(left, right, label):
    q, r = UR(left).quo_rem(UR(right))
    if r:
        raise ArithmeticError(
            f"{label}: polynomial quotient has nonzero remainder"
        )
    return UR(q)


def load_quartic_checkpoint():
    if not QUARTIC_CHECKPOINT.exists():
        return None
    try:
        with gzip.open(QUARTIC_CHECKPOINT, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
        if data.get("status") != "PASS_EXACT_YUN_QUARTIC":
            return None
        if data.get("plane_fingerprint") != FINAL_PLANE_FINGERPRINT:
            return None
        q = UR([VF(v) for v in data["quartic_coefficients"]])
        if int(q.degree()) != 4:
            return None
        print(
            "Q24DIVVALQQ_QUARTIC_CHECKPOINT|"
            f"file={QUARTIC_CHECKPOINT.name}|degree=4|status=HIT",
            flush=True,
        )
        return q
    except Exception as exc:
        print(
            "Q24DIVVALQQ_QUARTIC_CHECKPOINT|"
            f"status=INVALID|reason={type(exc).__name__}:{exc}",
            flush=True,
        )
        return None


def save_quartic_checkpoint(quartic, multiplicities):
    payload = {
        "schema": "elkies-k3.q24-divval-qq-quartic-yun.v1",
        "status": "PASS_EXACT_YUN_QUARTIC",
        "plane_fingerprint": FINAL_PLANE_FINGERPRINT,
        "degree": int(quartic.degree()),
        "quartic_coefficients": [str(v) for v in quartic.list()],
        "yun_multiplicities": list(map(int, multiplicities)),
    }
    tmp = QUARTIC_CHECKPOINT.with_name(
        QUARTIC_CHECKPOINT.name + f".tmp-{os.getpid()}"
    )
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, sort_keys=True)
    tmp.replace(QUARTIC_CHECKPOINT)
    print(
        "Q24DIVVALQQ_QUARTIC_CHECKPOINT|"
        f"file={QUARTIC_CHECKPOINT.name}|degree={quartic.degree()}|status=SAVED",
        flush=True,
    )


log("QUARTIC", method="DIRECT_NUMERATOR_YUN", status="BEGIN")

quartic = load_quartic_checkpoint()
square_factor = None

if quartic is None:
    AA0, BB0 = final_pairs[0][0], final_pairs[0][1]
    AA1, BB1 = final_pairs[1][0], final_pairs[1][1]

    A0U = lift_poly(AA0)
    A1U = lift_poly(AA1)
    B0U = lift_poly(BB0)
    B1U = lift_poly(BB1)
    ZU = lift_poly(Z)
    XU = lift_poly(X)
    YU = lift_poly(Y)
    AU = lift_poly(A)

    N = UR(A1U - VF(V) * A0U)
    D = UR(VF(V) * B0U - B1U)

    if not N or not D:
        raise ArithmeticError("degenerate exact chord pencil")

    log(
        "QUARTIC_INPUT",
        degN=N.degree(),
        degD=D.degree(),
        degZ=ZU.degree(),
        method="CLEARED_DIRECT_FORMULA_NESTED_QQU_V",
        status="PASS",
    )

    # Build the cleared numerator in (QQ[U])[V], NOT QQ(V)[U].
    #
    # N and D are linear in V, hence Qdisc has V-degree <= 4.  Doing the
    # degree-160 U arithmetic over the polynomial coefficient ring QQ[U]
    # avoids fraction-field normalization after every multiplication.
    log(
        "QUARTIC_NUMERATOR",
        method="NESTED_QQU_V",
        status="BEGIN",
    )

    RV = PolynomialRing(R, "W")
    W = RV.gen()

    Nv = RV(AA1) - W * RV(AA0)
    Dv = W * RV(BB0) - RV(BB1)

    if Nv.degree() != 1 or Dv.degree() != 1:
        raise ArithmeticError(
            f"unexpected V-degrees N/D={Nv.degree()}/{Dv.degree()}"
        )

    Nv2 = Nv * Nv
    Dv2 = Dv * Dv
    Nv4 = Nv2 * Nv2
    Dv3 = Dv2 * Dv
    Dv4 = Dv2 * Dv2

    ZR2 = Z * Z
    ZR4 = ZR2 * ZR2
    X2R = X * X

    log(
        "QUARTIC_NUMERATOR_PRODUCTS",
        degN4=Nv4.degree(),
        degD4=Dv4.degree(),
        status="PASS",
    )

    Qdisc_v = RV(
        Nv4
        - 6 * RV(X) * Nv2 * Dv2
        - 8 * RV(Y) * Nv * Dv3
        - 3 * RV(X2R) * Dv4
        - 4 * RV(A * ZR4) * Dv4
    )

    if not Qdisc_v:
        raise ArithmeticError("nested cleared chord discriminant vanished")
    if Qdisc_v.degree() > 4:
        raise ArithmeticError(
            f"cleared discriminant V-degree={Qdisc_v.degree()}, expected <=4"
        )

    log(
        "QUARTIC_NUMERATOR_NESTED",
        V_degree=Qdisc_v.degree(),
        U_degrees=",".join(
            str(R(c).degree()) if c else "-1"
            for c in Qdisc_v.list()
        ),
        status="PASS",
    )

    # Transpose (QQ[U])[V] -> QQ(V)[U] once.
    #
    # Only five V coefficient polynomials exist, so this performs a handful
    # of coefficient-field coercions instead of thousands during products.
    log(
        "QUARTIC_TRANSPOSE",
        source="QQ[U][V]",
        target="QQ(V)[U]",
        status="BEGIN",
    )

    Qdisc = UR.zero()
    for j, coeff_u in enumerate(Qdisc_v.list()):
        coeff_u = R(coeff_u)
        if coeff_u:
            Qdisc += lift_poly(coeff_u) * VF(V)**j
    Qdisc = UR(Qdisc)

    if not Qdisc:
        raise ArithmeticError("transposed chord discriminant vanished")

    expected_u_degree = max(
        (R(c).degree() for c in Qdisc_v.list() if c),
        default=-1,
    )
    if int(Qdisc.degree()) != int(expected_u_degree):
        raise ArithmeticError(
            f"transpose changed U degree "
            f"{expected_u_degree}->{Qdisc.degree()}"
        )

    log(
        "QUARTIC_TRANSPOSE",
        U_degree=Qdisc.degree(),
        V_degree=Qdisc_v.degree(),
        status="PASS",
    )

    log(
        "QUARTIC_NUMERATOR",
        method="NESTED_QQU_V",
        degree=Qdisc.degree(),
        status="PASS",
    )

    # denominator = (Z*D)^4; retain its factorization unit so our quartic
    # normalization exactly matches squarefree_binary_quartic().
    denom_base = UR(ZU * D)
    denom_lc = VF(denom_base.leading_coefficient())
    denom_unit = denom_lc**4
    denom_monic = UR(denom_base / denom_lc)

    numerator_lc = VF(Qdisc.leading_coefficient())

    log(
        "QUARTIC_GCD",
        degree=Qdisc.degree(),
        method="BIVARIATE_SINGULAR_YUN",
        status="BEGIN",
    )

    # Primitive/content split with respect to U, but entirely over QQ[V].
    PVFF = PolynomialRing(QQ, "Vff")
    Vff = PVFF.gen()

    qv_coeffs = [R(c) for c in Qdisc_v.list()]
    max_u_degree = max((c.degree() for c in qv_coeffs if c), default=-1)

    u_coeffs = []
    for i_u in range(max_u_degree + 1):
        u_coeffs.append(
            PVFF([
                QQ(c[i_u]) if i_u <= c.degree() else QQ(0)
                for c in qv_coeffs
            ])
        )

    nonzero_u_coeffs = [c for c in u_coeffs if c]
    if not nonzero_u_coeffs:
        raise ArithmeticError("cleared discriminant has no U coefficients")

    u_content = nonzero_u_coeffs[0]
    for c in nonzero_u_coeffs[1:]:
        u_content = u_content.gcd(c)
        if u_content.degree() == 0:
            break
    if u_content:
        u_content = u_content / u_content.leading_coefficient()

    primitive_u_coeffs = []
    for c in u_coeffs:
        if not c:
            primitive_u_coeffs.append(PVFF.zero())
            continue
        q, r = c.quo_rem(u_content)
        if r:
            raise ArithmeticError("U-content did not divide a coefficient")
        primitive_u_coeffs.append(PVFF(q))

    # Singular-backed bivariate gcd over QQ[U,V].
    BUV = PolynomialRing(QQ, names=("Uff", "Vff2"), order="degrevlex")
    Uff, Vff2 = BUV.gens()

    Fbi = BUV.zero()
    for i_u, coeff_v in enumerate(primitive_u_coeffs):
        for j_v, coeff in enumerate(coeff_v.list()):
            if coeff:
                Fbi += QQ(coeff) * Uff**i_u * Vff2**j_v

    if Fbi.degree(Uff) != max_u_degree:
        raise ArithmeticError(
            f"bivariate transpose changed U degree "
            f"{max_u_degree}->{Fbi.degree(Uff)}"
        )

    dFbi = Fbi.derivative(Uff)
    Gbi = Fbi.gcd(dFbi)

    def bivar_quo(left, right, label):
        q, r = BUV(left).quo_rem(BUV(right))
        if r:
            raise ArithmeticError(
                f"{label}: bivariate quotient has nonzero remainder"
            )
        return BUV(q)

    Wbi = bivar_quo(Fbi, Gbi, "yun-initial-bivariate")

    log(
        "QUARTIC_GCD",
        gcd_U_degree=Gbi.degree(Uff),
        gcd_V_degree=Gbi.degree(Vff2),
        squarefree_head_U_degree=Wbi.degree(Uff),
        method="BIVARIATE_SINGULAR_YUN",
        status="PASS",
    )

    odd_bi = BUV.one()
    square_bi = BUV.one()
    multiplicities = []
    i_mult = 1

    while Wbi.degree(Uff) > 0:
        Ybi = Wbi.gcd(Gbi)
        Zbi = bivar_quo(Wbi, Ybi, f"yun-z-bivariate-{i_mult}")

        if Zbi.degree(Uff) > 0:
            multiplicities.append(i_mult)
            if i_mult % 2:
                odd_bi *= Zbi
            half = i_mult // 2
            if half:
                square_bi *= Zbi**half

        Wbi = Ybi
        Gbi = bivar_quo(Gbi, Ybi, f"yun-g-bivariate-{i_mult}")
        i_mult += 1

    if Gbi.degree(Uff) > 0:
        raise ArithmeticError(
            f"bivariate Yun ended with residual U-degree {Gbi.degree(Uff)}"
        )

    recon_bi = odd_bi * square_bi**2
    recon_unit, recon_rem = Fbi.quo_rem(recon_bi)
    if recon_rem:
        raise ArithmeticError("bivariate Yun reconstruction has remainder")
    if recon_unit.degree(Uff) > 0 or recon_unit.degree(Vff2) > 0:
        raise ArithmeticError(
            "bivariate Yun reconstruction left a nonconstant unit"
        )

    def bivar_to_ur(poly):
        by_u = {}
        for exponent, coeff in poly.dict().items():
            i_u, j_v = exponent
            by_u.setdefault(i_u, {})[j_v] = QQ(coeff)
        if not by_u:
            return UR.zero()
        max_i = max(by_u)
        coeffs = []
        for i_u in range(max_i + 1):
            terms = by_u.get(i_u, {})
            coeffs.append(
                VF(sum(
                    coeff * V**j_v
                    for j_v, coeff in terms.items()
                ))
            )
        return UR(coeffs)

    odd_ur = bivar_to_ur(odd_bi)
    square_ur = bivar_to_ur(square_bi)
    if not odd_ur or not square_ur:
        raise ArithmeticError("empty Yun odd/square factor")

    odd_monic = UR(odd_ur / odd_ur.leading_coefficient())
    square_monic = UR(square_ur / square_ur.leading_coefficient())

    # Exact QQ[U,V] reconstruction above already proves the square-class.
    # Avoid re-expanding a degree-78 square in QQ(V)[U].
    log(
        "QUARTIC_VERIFY_SQUARE",
        odd_degree=odd_monic.degree(),
        square_degree=square_monic.degree(),
        method="LEADING_COEFFICIENT_CERTIFICATE",
        status="BEGIN",
    )

    content_vf = VF(sum(
        QQ(c) * V**j
        for j, c in enumerate(u_content.list())
        if c
    ))
    recon_scalar = QQ(recon_unit.constant_coefficient())

    expected_numerator_lc = VF(
        content_vf
        * recon_scalar
        * odd_ur.leading_coefficient()
        * square_ur.leading_coefficient()**2
    )

    if expected_numerator_lc != numerator_lc:
        raise ArithmeticError(
            "fraction-free Yun leading-coefficient certificate failed"
        )

    if (
        int(odd_monic.degree()) + 2 * int(square_monic.degree())
        != int(Qdisc.degree())
    ):
        raise ArithmeticError(
            "fraction-free Yun degree reconstruction failed"
        )

    log(
        "QUARTIC_VERIFY_SQUARE",
        odd_degree=odd_monic.degree(),
        square_degree=square_monic.degree(),
        method="LEADING_COEFFICIENT_CERTIFICATE",
        status="PASS",
    )

    quartic = UR((numerator_lc / denom_unit) * odd_monic)
    square_factor = UK(square_monic) / UK(denom_monic**2)

    log(
        "QUARTIC_YUN",
        method="BIVARIATE_SINGULAR_YUN",
        pieces=len(multiplicities),
        multiplicities=",".join(map(str, multiplicities)),
        square_degree=square_monic.degree(),
        quartic_degree=quartic.degree(),
        status="PASS",
    )

    save_quartic_checkpoint(quartic, multiplicities)

quartic_degree = int(quartic.degree())
if quartic_degree != 4:
    raise ArithmeticError(
        f"exact q24 quartic degree {quartic_degree}, expected modular degree 4"
    )
log("QUARTIC", method="BIVARIATE_SINGULAR_YUN", degree=quartic_degree, status="PASS")

log("INVARIANTS", status="BEGIN")
I, J = binary_quartic_invariants(quartic)
jacA = VF(-27) * VF(I)
jacB = VF(-27) * VF(J)
log("INVARIANTS", status="PASS")

def fast_classify_finite_short_weierstrass(base_ring, coefficient_a, coefficient_b):
    base_field = base_ring.fraction_field()
    coefficient_a = base_field(coefficient_a)
    coefficient_b = base_field(coefficient_b)
    discriminant = -16 * (
        4 * coefficient_a**3 + 27 * coefficient_b**2
    )
    if not discriminant:
        raise ValueError("short Weierstrass model is singular")

    # One support polynomial is enough:
    # poles of A/B plus zero/pole support of Delta.
    support_inputs = (
        base_ring(coefficient_a.denominator()),
        base_ring(coefficient_b.denominator()),
        base_ring(discriminant.numerator()),
        base_ring(discriminant.denominator()),
    )

    support = base_ring.one()
    for poly in support_inputs:
        if poly.degree() <= 0:
            continue
        deriv = poly.derivative()
        sf = poly
        if deriv:
            repeated = poly.gcd(deriv)
            q, r = poly.quo_rem(repeated)
            if r:
                raise ArithmeticError("squarefree support quotient failed")
            sf = base_ring(q)
        common = support.gcd(sf)
        q, r = sf.quo_rem(common)
        if r:
            raise ArithmeticError("support lcm quotient failed")
        support *= q

    if support.degree() > 0:
        support /= support.leading_coefficient()

    log(
        "CHILD_SUPPORT",
        degree=support.degree(),
        method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
        status="BEGIN",
    )
    factorization = support.factor() if support.degree() > 0 else ()
    factors = tuple(f for f, unused_e in factorization)
    log(
        "CHILD_SUPPORT",
        degree=support.degree(),
        factors=len(factors),
        method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
        status="PASS",
    )

    def poly_order(poly, factor):
        poly = base_ring(poly)
        if not poly:
            raise ArithmeticError("valuation requested on zero polynomial")
        n = 0
        while poly.degree() >= factor.degree():
            q, r = poly.quo_rem(factor)
            if r:
                break
            n += 1
            poly = base_ring(q)
        return n

    def rf_order(value, factor):
        value = base_field(value)
        return (
            poly_order(value.numerator(), factor)
            - poly_order(value.denominator(), factor)
        )

    scaling_unit = base_field.one()
    finite_data = []
    for factor in sorted(factors, key=str):
        raw_a = rf_order(coefficient_a, factor)
        raw_b = rf_order(coefficient_b, factor)
        raw_delta = rf_order(discriminant, factor)
        scaling = min(raw_a // 4, raw_b // 6)
        scaling_unit *= base_field(factor)**(-scaling)
        finite_data.append({
            "factor": factor,
            "raw_orders": (raw_a, raw_b, raw_delta),
            "scaling": scaling,
        })

    minimal_a_rf = coefficient_a * scaling_unit**4
    minimal_b_rf = coefficient_b * scaling_unit**6
    minimal_delta_rf = discriminant * scaling_unit**12

    if any(
        value.denominator() != 1
        for value in (minimal_a_rf, minimal_b_rf, minimal_delta_rf)
    ):
        raise ArithmeticError(
            "fast finite minimization did not clear all base denominators"
        )

    minimal_a = base_ring(minimal_a_rf.numerator())
    minimal_b = base_ring(minimal_b_rf.numerator())
    minimal_delta = base_ring(minimal_delta_rf.numerator())

    finite_fibres = []
    root_rank = ZZ(0)
    euler_number = ZZ(0)
    root_determinant = ZZ(1)

    for item in finite_data:
        factor = item["factor"]
        orders = (
            poly_order(minimal_a, factor),
            poly_order(minimal_b, factor),
            poly_order(minimal_delta, factor),
        )
        item["minimal_orders"] = orders

        if min(orders[0] // 4, orders[1] // 6) != 0:
            raise ArithmeticError(
                "fast finite minimization left a fourth/sixth-power factor"
            )
        if orders[2] == 0:
            continue

        rank, euler, determinant, symbol = kodaira_data_from_short_orders(
            *orders
        )
        degree = ZZ(factor.degree())
        finite_fibres.append({
            "factor": factor,
            "degree": int(degree),
            "raw_orders": tuple(map(int, item["raw_orders"])),
            "scaling": int(item["scaling"]),
            "minimal_orders": tuple(map(int, orders)),
            "kodaira": symbol,
            "root_rank": int(rank),
            "euler_number": int(euler),
            "root_determinant": int(determinant),
        })
        root_rank += degree * rank
        euler_number += degree * euler
        root_determinant *= ZZ(determinant)**degree

    infinity_raw = tuple(
        -value.degree()
        for value in (minimal_a, minimal_b, minimal_delta)
    )
    infinity_scaling = min(
        infinity_raw[0] // 4,
        infinity_raw[1] // 6,
    )
    infinity_normalized = tuple(
        infinity_raw[i] - (4, 6, 12)[i] * infinity_scaling
        for i in range(3)
    )

    return {
        "finite_minimization": {
            "minimal_a": minimal_a,
            "minimal_b": minimal_b,
            "minimal_discriminant": minimal_delta,
            "scaling_unit": scaling_unit,
            "finite_places": tuple(finite_data),
            "infinity": {
                "raw_orders": infinity_raw,
                "scaling": infinity_scaling,
                "normalized_orders": infinity_normalized,
            },
        },
        "finite_fibres": tuple(finite_fibres),
        "finite_root_rank": int(root_rank),
        "finite_euler_number": int(euler_number),
        "finite_root_determinant": int(root_determinant),
        "infinity_boundary": {
            "raw_orders": infinity_raw,
            "scaling": infinity_scaling,
            "normalized_orders": infinity_normalized,
        },
    }

log(
    "CHILD_CLASSIFICATION",
    method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
    status="BEGIN",
)
classification = fast_classify_finite_short_weierstrass(VR, jacA, jacB)
log(
    "CHILD_CLASSIFICATION",
    method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
    status="PASS",
)
root_rank = int(classification["finite_root_rank"])
root_euler = int(classification["finite_euler_number"])
root_det = int(classification["finite_root_determinant"])

infinity = classification["infinity_boundary"]
infinity_orders = tuple(map(int, infinity["normalized_orders"]))
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    ir, ie, idet, infinity_kind = kodaira_data_from_short_orders(
        *infinity_orders
    )
    root_rank += int(ir)
    root_euler += int(ie)
    root_det *= int(idet)

if (root_rank, root_det, root_euler) != (12, 4, 24):
    raise ArithmeticError(
        "exact child is not D12: "
        f"rank={root_rank}, det={root_det}, euler={root_euler}"
    )

minimal_A = classification["finite_minimization"]["minimal_a"]
minimal_B = classification["finite_minimization"]["minimal_b"]
minimal_Delta = classification["finite_minimization"]["minimal_discriminant"]

log(
    "CHILD",
    root_rank=root_rank,
    root_det=root_det,
    euler=root_euler,
    MW=5,
    status="PASS_D12",
)

# ---------------------------------------------------------------------------
# Modular quartic/Jacobian regression.
# ---------------------------------------------------------------------------
def normalized_rf_mod(value):
    value = VF(value)
    n = VR(value.numerator())
    d = VR(value.denominator())
    lc = d.leading_coefficient()
    n /= lc
    d /= lc
    return {
        "num": [int(red_q(c)) for c in n.list()],
        "den": [int(red_q(c)) for c in d.list()],
    }

def trimmed(values):
    values = list(values)
    while values and values[-1] == 0:
        values.pop()
    return values

def assert_sig_rf(value, expected, label):
    got = normalized_rf_mod(value)
    if (
        trimmed(got["num"]) != trimmed(expected["num"])
        or trimmed(got["den"]) != trimmed(expected["den"])
    ):
        raise ArithmeticError(f"{label} does not reduce to modular signature")

for i in range(5):
    assert_sig_rf(
        VF(quartic[i]),
        signature["quartic_coefficients"][i],
        f"quartic[{i}]",
    )
assert_sig_rf(jacA, signature["jacobian_A"], "jacobian_A")
assert_sig_rf(jacB, signature["jacobian_B"], "jacobian_B")

log("MODULAR_SIGNATURE", plane=1, quartic=1, jacobian=1, status="PASS")

def qlist(poly):
    return [str(v) for v in poly.list()]

finite_data = [
    {
        "factor": str(item["factor"]),
        "degree": int(item["degree"]),
        "minimal_orders": list(map(int, item["minimal_orders"])),
        "kodaira": item["kodaira"],
    }
    for item in classification["finite_fibres"]
]

payload = {
    "schema": "elkies-k3.h92-q24-d13-to-d12-component-valuation-qq.v1",
    "status": "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR",
    "inputs": {
        "exact_construction_manifest": str(MANIFEST.relative_to(ROOT)),
        "component_rr_modp": str(COMPONENT_MOD.relative_to(ROOT)),
        "signature_modp": str(SIGNATURE.relative_to(ROOT)),
    },
    "rr": {
        "ambient_dimension": 56,
        "smooth_collision_rank": 48,
        "post_collision_dimension": 8,
        "resolved_component_rank": 6,
        "kernel_dimension": 2,
        "geometric_fibre_twist": -8,
        "component_ledger": component_ledger,
        "redundant_exact_condition_rows": int(verified_rows),
        "generic_audit_components": int(verified_components),
        "generic_audit_charts": int(verified_charts),
        "generic_audit_fallbacks": int(audit_fallbacks),
        "plane_2x56": [
            [str(v) for v in row] for row in final56.rows()
        ],
    },
    "quartic": {
        "degree": quartic_degree,
        "coefficients_in_U_low_to_high": [str(v) for v in quartic.list()],
        "binary_quartic_I": str(I),
        "binary_quartic_J": str(J),
    },
    "jacobian_raw": {
        "A": str(jacA),
        "B": str(jacB),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": qlist(minimal_A),
        "minimal_B_coefficients_low_to_high": qlist(minimal_B),
        "minimal_discriminant_coefficients_low_to_high": qlist(minimal_Delta),
        "finite_fibres": finite_data,
        "infinity_orders": list(infinity_orders),
        "infinity_kind": infinity_kind,
        "root_rank": root_rank,
        "root_determinant": root_det,
        "euler_number": root_euler,
        "MW_rank_if_rho19": 5,
    },
    "verification": {
        "exact_q24_input": True,
        "exact_collision": True,
        "all_redundant_component_valuations": True,
        "mod_100003_plane_regression": True,
        "mod_100003_quartic_regression": True,
        "mod_100003_jacobian_regression": True,
        "D12_root_data": True,
    },
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q24DIVVALQQ_RESULT|"
    "ambient=56|collision=48|post=8|resolved=6|kernel=2|"
    "quartic=4|root_rank=12|root_det=4|euler=24|MW=5|"
    "status=PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR",
    flush=True,
)
