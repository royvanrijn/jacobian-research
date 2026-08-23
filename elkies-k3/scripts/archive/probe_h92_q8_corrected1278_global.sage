#!/usr/bin/env sage -python
"""
First global rank probe for the corrected H92 q8 movable class.

The previous degree-18 class was isotropic but not fully nef.  Complete chamber
reduction removes one horizontal bisection and then the vertical chain

    E7_2, E7_4, E7_3, E7_1,

leaving a primitive isotropic movable class of old-fibre degree 16 with

    generic support: 8*O + 8*(-P1)
    vertical difference:
      (-10,0,2,4,5,8,7,6,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0).

This probe rebuilds the generic/global compiler rather than recycling q6^9
rows:

  * generic basis: 1,m,...,m^8, x,xm,...,xm^6  (rank 16);
  * common smooth-collision denominator h^16;
  * endpoint caps d<=78 (m-family), d<=80 (x-family), 1278 columns;
  * E8 global fibre twist -10F and unchanged complete ideal (u^2,X,Y);
  * E7 target
        c_mov = 8*c6 + (2,6,8,5,6,4,7),
    so exact generic residues use q6^8, not q6^9.

The existing exact component-function-field evaluators are reused only after
being mechanically patched in a temporary directory from q6^9 to q6^8 and to
the new Cartier exponents.  The repository files are not modified.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_global.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
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
        home / "Documents" / "jacobian-research",
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
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
SCRIPTS = ROOT / "elkies-k3" / "scripts"
P = int(args.prime)
if not ZZ(P).is_prime() or P in (2, 3):
    raise ValueError("prime must be odd and != 3")
finite = GF(P)

sage = shutil.which("sage")
if not sage:
    raise SystemExit("sage executable not found")

CHAMBER = GEN / "zz-h92-q8-complete-chamber-reduction.json"
P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
FRAMES = GEN / "elkies-k3-h92-q8-generic-component-chart-frames.json"

PREFIX = "zz-h92-q8-corrected1278"
AMBIENT_PATH = GEN / f"{PREFIX}-ambient.json"
COND_PATH = GEN / f"{PREFIX}-generic-conditions.json"
R13_PATH = GEN / f"{PREFIX}-e7-1-3-residue-rows.json"
R56_PATH = GEN / f"{PREFIX}-e7-5-6-residue-rows.json"
R47_PATH = GEN / f"{PREFIX}-e7-4-7-residue-rows.json"
KERNEL_PATH = GEN / f"{PREFIX}-global-kernel-mod-{P}.json"

for path in (CHAMBER, P1_PATH, FRAMES):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

chamber = json.loads(CHAMBER.read_text())
assert chamber["status"] == "PASS_BOUNDED_COMPLETE_Q8_CHAMBER_REDUCTION"
assert chamber["final_old_fibre_degree"] == 16
assert chamber["final_square"] == 0
assert chamber["final_primitive_pairing_gcd"] == 1
assert chamber["generic_fibre_support"] == "8*O + 8*(-P1)"
EXPECTED_VERTICAL = [
    -10,0,2,4,5,8,7,6,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0
]
assert chamber["vertical_difference_from_generic_support"] == EXPECTED_VERTICAL
assert chamber["final_fiber_pairings"][:7] == [1,0,0,0,0,1,2]
assert chamber["final_fiber_pairings"][7:15] == [1,0,0,0,0,0,0,0]

# Actual resolved E7 component order is pinned in the repository.
# source E7_i -> resolved E7_j:
SOURCE_TO_RESOLVED = (1, 6, 4, 3, 7, 2, 5)
source_deg = chamber["final_fiber_pairings"][:7]
resolved_deg = [0] * 7
for source_index, resolved_index in enumerate(SOURCE_TO_RESOLVED):
    resolved_deg[resolved_index - 1] = source_deg[source_index]
assert resolved_deg == [1,1,0,0,2,0,0]

# Exact reduced q6 E7 cycle and corrected q6^8 integral twist.
# c_mov=(-6,-10,-16,-11,-6,-8,-13)
# c6=(-1,-2,-3,-2,-3/2,-3/2,-5/2)
E7_TWIST = (2, 6, 8, 5, 6, 4, 7)
assert E7_TWIST == (2, 6, 8, 5, 6, 4, 7)

# Actual E7 generic valuations from the certified atlas.
V_T = (2, 2, 4, 3, 1, 2, 3)
V_X = (2, 4, 6, 4, 2, 3, 5)
V_M = (1, 1, 3, 2, 0, 2, 2)
COMPONENTS = tuple(f"E7_{i}" for i in range(1, 8))

# ---------------------------------------------------------------------------
# 1. Corrected rank-16 generic basis and sharp 1278-column endpoint envelope.
# ---------------------------------------------------------------------------

families = (
    [{"kind": "m_power", "x_power": 0, "m_power": b} for b in range(9)]
    + [{"kind": "x_m_power", "x_power": 1, "m_power": b} for b in range(7)]
)
assert len(families) == 16

ambient_basis = []
for family in families:
    a = family["x_power"]
    b = family["m_power"]
    max_d = 78 if a == 0 else 80
    for d in range(max_d + 1):
        ambient_basis.append({
            "kind": family["kind"],
            "x_power": a,
            "m_power": b,
            # q6^8 helper normalization: patched residue evaluator uses
            # (4*k-i-8), hence i=d-8 makes this exactly 4*k-d.
            "u_power": d - 8,
            "h_power": 16,
            "actual_u_power": d,
            "coefficient": f"u^{d}/h(u)^16",
        })

N = len(ambient_basis)
assert N == 9*79 + 7*81 == 1278

ambient = {
    "schema": "elkies-k3.h92-q8-corrected1278-helper-ambient.v1",
    "status": "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
    "ambient_dimension": N,
    "ambient_basis": ambient_basis,
    "source_movable_class": chamber["final_class"],
    "generic_fibre_divisor": "8*O + 8*(-P1)",
    "generic_fibre_degree": 16,
    "vertical_difference": EXPECTED_VERTICAL,
    "normalization": {
        "actual_coefficient": "u^d/h(u)^16*x^a*m^b",
        "stored_u_power": "d-8",
        "e7_tensor_power": 8,
        "e7_integral_twist": list(E7_TWIST),
        "e8_global_fibre_twist": -10,
        "marked_e7_frame": [
            "m^b/t^14",
            "x*m^b/t^16",
        ],
        "degree_caps": {"m_family": 78, "x_family": 80},
    },
}
AMBIENT_PATH.write_text(json.dumps(ambient, indent=2, sort_keys=True) + "\n")
print(
    f"Q8CORRECTEDAMBIENT|columns={N}|basis=16|h=16|"
    "m_cap=78|xm_cap=80|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Corrected q6^8 E7 valuation condition template.
# ---------------------------------------------------------------------------

component_conditions = []
for component_index, component in enumerate(COMPONENTS):
    grouped = {}
    for basis_index, entry in enumerate(ambient_basis):
        a = ZZ(entry["x_power"])
        b = ZZ(entry["m_power"])
        i = ZZ(entry["u_power"])
        k = ZZ(entry["h_power"])
        # q6^8 helper: w+(4k-i-8)*ord(t)+...
        residual = (
            ZZ(E7_TWIST[component_index])
            + (4*k - i - 8) * ZZ(V_T[component_index])
            + a * ZZ(V_X[component_index])
            + b * ZZ(V_M[component_index])
        )
        # With i=d-8 this must equal the reduced true coefficient formula.
        d = ZZ(entry["actual_u_power"])
        direct = (
            ZZ(E7_TWIST[component_index])
            + (64-d) * ZZ(V_T[component_index])
            + a * ZZ(V_X[component_index])
            + b * ZZ(V_M[component_index])
        )
        assert residual == direct
        if residual < 0:
            grouped.setdefault(int(residual), []).append(basis_index)

    records = [
        {
            "residual_order": order,
            "basis_indices": indices,
            "basis_labels": [ambient_basis[i] for i in indices],
            "singleton_exact_linear_condition": len(indices) == 1,
        }
        for order, indices in sorted(grouped.items())
    ]
    component_conditions.append({
        "component": component,
        "orders": {
            "t": V_T[component_index],
            "x": V_X[component_index],
            "m": V_M[component_index],
            "q8_twist": E7_TWIST[component_index],
            "q6_eighth_generator": 8*V_T[component_index],
        },
        "residual_formula": "w+(4*k-i-8)*ord(t)+a*ord(x)+b*ord(m)",
        "negative_order_groups": records,
        "singleton_condition_count": sum(
            row["singleton_exact_linear_condition"] for row in records
        ),
    })

singleton_indices = sorted({
    row["basis_indices"][0]
    for condition in component_conditions
    for row in condition["negative_order_groups"]
    if row["singleton_exact_linear_condition"]
})

conditions = {
    "schema": "elkies-k3.h92-q8-corrected-e7-generic-conditions.v1",
    "status": "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE",
    "inputs": {
        "endpoint_ambient": {
            "path": str(AMBIENT_PATH.relative_to(ROOT)),
            "sha256": digest(AMBIENT_PATH),
        }
    },
    "component_order": list(COMPONENTS),
    "ambient_basis_sha256": canonical_digest([
        [
            e["kind"], e["x_power"], e["m_power"],
            e["u_power"], e["h_power"]
        ]
        for e in ambient_basis
    ]),
    "q8_twist": list(E7_TWIST),
    "coefficient_convention": (
        "corrected movable q8: u^i/h^k*x^a*m^b with helper i=d-8, "
        "q6 tensor power 8"
    ),
    "component_conditions": component_conditions,
    "singleton_coordinate_block": {
        "basis_indices": singleton_indices,
        "matrix": [
            [1 if col == index else 0 for col in range(N)]
            for index in singleton_indices
        ],
        "rank": len(singleton_indices),
    },
}
COND_PATH.write_text(json.dumps(conditions, indent=2, sort_keys=True) + "\n")

negative_groups = sum(
    len(c["negative_order_groups"]) for c in component_conditions
)
print(
    "Q8CORRECTEDE7CONDITIONS|"
    f"negative_groups={negative_groups}|singletons={len(singleton_indices)}|"
    f"twist={','.join(map(str,E7_TWIST))}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Reuse exact component-function-field evaluators, patched q6^9 -> q6^8.
# ---------------------------------------------------------------------------

def patched(source_name, replacements, target_dir):
    source = (SCRIPTS / source_name).read_text()
    for old, new in replacements:
        if old not in source:
            raise RuntimeError(
                f"patch anchor not found in {source_name}: {old!r}"
            )
        source = source.replace(old, new)
    path = target_dir / source_name
    path.write_text(source)
    return path


with tempfile.TemporaryDirectory(prefix="h92-q8-q6pow8-") as tmp:
    tmp = Path(tmp)

    p13 = patched(
        "derive_h92_q8_e7_1_3_generic_residue_rows.sage",
        [
            (
                'for component, g_power in (("E7_1", (2, 4)), '
                '("E7_2", (5, 6)), ("E7_3", (6, 4))):',
                'for component, g_power in (("E7_1", (2, 5)), '
                '("E7_2", (6, 6)), ("E7_3", (8, 5))):',
            ),
            ("degree-9", "degree-8"),
        ],
        tmp,
    )
    p56 = patched(
        "derive_h92_q8_e7_5_6_generic_residue_rows.sage",
        [
            (
                'g_y_power = {"E7_5": 5, "E7_6": 6}[component]',
                'g_y_power = {"E7_5": 6, "E7_6": 8}[component]',
            ),
            ("degree-9", "degree-8"),
            (
                '9 - (4*int(basis[index]["h_power"])-int(basis[index]["u_power"]))',
                '8 - (4*int(basis[index]["h_power"])-int(basis[index]["u_power"]))',
            ),
            ("component_ring(lx)**9", "component_ring(lx)**8"),
        ],
        tmp,
    )
    p47 = patched(
        "derive_h92_q8_e7_4_7_generic_residue_rows.sage",
        [
            ("g = series_field(e**4 * s**6)", "g = series_field(e**5 * s**8)"),
            ("g = series_field(e**5 * s**6)", "g = series_field(e**7 * s**8)"),
            ("degree-9", "degree-8"),
        ],
        tmp,
    )

    def run(path, output):
        cmd = [
            sage, "-python", str(path),
            "--frames", str(FRAMES),
            "--conditions", str(COND_PATH),
            "--p1", str(P1_PATH),
            "--output", str(output),
        ]
        print("RUN|" + " ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=ROOT, check=True)

    run(p13, R13_PATH)
    run(p56, R56_PATH)
    run(p47, R47_PATH)

# ---------------------------------------------------------------------------
# 4. Smooth h block on the corrected h^-16 frame.
# ---------------------------------------------------------------------------

p1 = json.loads(P1_PATH.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"

u_ring = PolynomialRing(finite, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()


def qq_to_finite(value):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ZeroDivisionError(f"prime {P} divides denominator {value}")
    return finite(ZZ(value.numerator())) / den


def finite_poly(values):
    return u_ring([qq_to_finite(v) for v in values])


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
        raise ZeroDivisionError("non-unit denominator in residue computation")
    return u_ring((num * den.inverse_mod(modulus)) % modulus)


H_POLE = 24
h_modulus = h**H_POLE
h_dim = h_modulus.degree()
assert h_dim == 96

rho = u_field(h) * y_p / x_p
assert gcd(h, u_ring(rho.numerator())) == 1
assert gcd(h, u_ring(rho.denominator())) == 1

frame_coordinates = [(0,j) for j in range(9)] + [(1,j) for j in range(7)]
coord_index = {v:i for i,v in enumerate(frame_coordinates)}
H = matrix(finite, len(frame_coordinates)*h_dim, N)
assert H.nrows() == 1536

for col, entry in enumerate(ambient_basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    for j in range(b+1):
        exponent = 2*j - b - 16 - 2*a
        if exponent >= 0:
            continue
        value = (
            finite(binomial(b,j))
            * u**d
            * rho**(b-j)
            * h**(H_POLE + exponent)
        )
        remainder = residue_mod(value, h_modulus)
        offset = coord_index[(a,j)] * h_dim
        for degree, coefficient in enumerate(remainder.list()):
            if coefficient:
                H[offset+degree, col] += coefficient

# ---------------------------------------------------------------------------
# 5. E8: same complete ideal (u^2,X,Y), but global -10F.
# ---------------------------------------------------------------------------

floor_count = 0
for entry in ambient_basis:
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    floor = (10 + 2*b) if a == 0 else (14 + 2*b)
    if d < floor:
        floor_count += 1
assert floor_count == 302

E8 = matrix(finite, floor_count + 2, N)
row = 0
for col, entry in enumerate(ambient_basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    floor = (10 + 2*b) if a == 0 else (14 + 2*b)
    if d < floor:
        E8[row,col] = 1
        row += 1
assert row == floor_count

u2 = u**2
s_e8 = u_field(u**2) * y_p / x_p
assert residue_mod(s_e8, u2)[0]

for col, entry in enumerate(ambient_basis):
    if int(entry["x_power"]) != 0:
        continue
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    floor = 10 + 2*b
    if d < floor or d-floor >= 2:
        continue
    value = u_field(u**(d-floor)) * s_e8**b / u_field(h**16)
    remainder = residue_mod(value, u2)
    for jet in range(2):
        E8[floor_count+jet,col] = remainder[jet]

# ---------------------------------------------------------------------------
# 6. Stack corrected exact generic E7 residue rows.
# ---------------------------------------------------------------------------

e7_specs = [("singleton", idx, None) for idx in singleton_indices]
row_counts = []
for path in (R13_PATH, R56_PATH, R47_PATH):
    payload = json.loads(path.read_text())
    before = len(e7_specs)
    for component in payload["components"]:
        for r in component["non_singleton_residue_rows"]:
            e7_specs.append((
                component["component"],
                int(r["residual_order"]),
                r["entries"],
            ))
    row_counts.append(len(e7_specs)-before)

E7 = matrix(finite, len(e7_specs), N)
for ri, spec in enumerate(e7_specs):
    if spec[0] == "singleton":
        E7[ri, int(spec[1])] = 1
    else:
        for item in spec[2]:
            E7[ri, int(item["basis_index"])] += qq_to_finite(
                item["coefficient"]
            )

# ---------------------------------------------------------------------------
# 7. First corrected global ranks.
# ---------------------------------------------------------------------------

rank_h = int(H.rank())
rank_e8 = int(E8.rank())
rank_e7 = int(E7.rank())
E8H = E8.stack(H)
rank_e8h = int(E8H.rank())
GLOBAL = E8H.stack(E7)
rank_global = int(GLOBAL.rank())
kernel = GLOBAL.right_kernel().basis_matrix()
kernel_dim = int(kernel.nrows())

print(
    "Q8CORRECTEDGLOBAL|"
    f"prime={P}|ambient={N}|"
    f"H={H.nrows()}:{rank_h}:{N-rank_h}|"
    f"E8={E8.nrows()}:{rank_e8}:{N-rank_e8}|"
    f"E7GEN={E7.nrows()}:{rank_e7}:{N-rank_e7}|"
    f"E8H={E8H.nrows()}:{rank_e8h}:{N-rank_e8h}|"
    f"ALL={GLOBAL.nrows()}:{rank_global}:{kernel_dim}|"
    f"e7_rows_by_evaluator={','.join(map(str,row_counts))}",
    flush=True,
)

KERNEL_PATH.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-corrected1278-global-kernel-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED_Q8_GLOBAL_GENERIC_KERNEL",
    "prime": P,
    "ambient": str(AMBIENT_PATH.relative_to(ROOT)),
    "conditions": str(COND_PATH.relative_to(ROOT)),
    "dimensions": {
        "ambient": N,
        "H_rows": H.nrows(),
        "H_rank": rank_h,
        "E8_rows": E8.nrows(),
        "E8_rank": rank_e8,
        "E7_rows": E7.nrows(),
        "E7_rank": rank_e7,
        "E8H_rank": rank_e8h,
        "global_rank": rank_global,
        "kernel": kernel_dim,
    },
    "kernel_basis_rows": [
        [int(v) for v in row] for row in kernel.rows()
    ],
    "boundary": (
        "This is the corrected movable degree-16 global generic-component "
        "compiler. Finite E7 translated/node/marked overlap conditions have "
        "not yet been added."
    ),
}, indent=2, sort_keys=True) + "\n")

print(f"GLOBAL_KERNEL|path={KERNEL_PATH}|dimension={kernel_dim}", flush=True)
