#!/usr/bin/env sage -python
"""
Inspect the first D13 q24 neighbour as a divisor on the already-exact q6 parent.

This follows the repaired H3 neighbour architecture:

    E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4

For the next hop, do NOT first reconstruct the q24 horizontal point on D13.
Instead pull the WHOLE q24 isotropic divisor back to the q6 simple frame and
write its restriction to the generic q6 fibre as

    D24 ~ (d-1) O6 + Q + V,

where Q is the effective q6 section representing the Pic^0 class and V is
vertical (F6 plus E8/E6 simple components).

If d=4, H^0 on the generic q6 fibre has the explicit four-function basis

    1, x, y, m_Q=(y+y_Q)/(x-x_Q),

before vertical/resolved-component conditions are imposed.
"""

import json
import sys
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, lcm, matrix, vector
)


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h / "Documents" / "jacobian-research",
        h / "jacobian-research",
        h / "src" / "jacobian-research",
        h / "git" / "jacobian-research",
        h / "projects" / "jacobian-research",
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
AUDIT = SCRIPTS / "audit_h3_d13_q24_pullback_invariants.sage"
TMP = LOCAL / "q24-q6-rr-route-native-audit.json"
OUT = LOCAL / "q24-q6-rr-route.json"

if not AUDIT.exists():
    raise SystemExit(f"missing {AUDIT}")


def run_script(path, argv):
    saved = list(sys.argv)
    scope = {"__name__": "__embedded__"}
    try:
        sys.argv = [str(path)] + list(argv)
        exec(compile(path.read_text(), str(path), "exec"), scope)
    finally:
        sys.argv = saved
    return scope


print("Q24Q6RR|stage=pullback_audit", flush=True)
s = run_script(AUDIT, ["--output", str(TMP)])

needed = (
    "q6_simple_ns", "simple_frame", "F_q6_simple", "O_q6_abstract",
    "P_simple", "O8_simple", "D_simple", "F8_simple", "Bsimple",
)
missing = [name for name in needed if name not in s]
if missing:
    raise SystemExit(f"pullback audit did not expose {missing}")

ns = s["q6_simple_ns"]
G = s["simple_frame"]
F6 = vector(ZZ, s["F_q6_simple"])
O6 = vector(ZZ, s["O_q6_abstract"])
P24 = vector(ZZ, s["P_simple"])
O8 = vector(ZZ, s["O8_simple"])
D24 = vector(ZZ, s["D_simple"])
F8 = vector(ZZ, s["F8_simple"])
Bsimple = s["Bsimple"]

assert G.nrows() == G.ncols() == 17
assert ns.nrows() == ns.ncols() == 19
assert F6 * ns * F6 == 0
assert O6 * ns * O6 == -2
assert O6 * ns * F6 == 1
assert D24 * ns * D24 == 0
assert D24 * ns * F8 == 2
assert F8 * ns * F6 == 2

dP = ZZ(P24 * ns * F6)
dO8 = ZZ(O8 * ns * F6)
dD = ZZ(D24 * ns * F6)

print(
    "Q24Q6RR_DEGREE|"
    f"P24={dP}|O8={dO8}|F8={F8*ns*F6}|D24={dD}|"
    f"identity={int(dD == dO8 + dP - 20)}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# q6 root/MW decomposition.
# The q6 simple positive frame is E8+E6 roots (first 14) + MW rank 3.
# ---------------------------------------------------------------------------

root = G[:14, :14]
coupling = G[:14, 14:]
tail = G[14:, 14:]
H = tail - coupling.transpose() * root.inverse() * coupling

assert root.rank() == 14
assert H.nrows() == H.ncols() == 3

zD = vector(ZZ, D24[-3:])
zP = vector(ZZ, P24[-3:])
zO8 = vector(ZZ, O8[-3:])

hD = QQ(zD * H * zD)
hP = QQ(zP * H * zP)
hO8 = QQ(zO8 * H * zO8)

print(
    "Q24Q6RR_PIC0|"
    f"D_mw={','.join(map(str,zD))}|D_height={hD}|"
    f"P_mw={','.join(map(str,zP))}|P_height={hP}|"
    f"O8_mw={','.join(map(str,zO8))}|O8_height={hO8}|"
    "status=PASS",
    flush=True,
)


def coset_order(dual):
    result = ZZ(1)
    for value in dual:
        result = lcm(result, ZZ(QQ(value).denominator()))
    return result


def effective_q6_section(z, cap=200000):
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * 14 + list(z))
    pairing = vector(QQ, base * G[:, :14])
    dual = pairing * root.inverse()
    order = coset_order(dual)

    # E8 is unimodular; the only nontrivial discriminant contribution is E6.
    if order == 1:
        correction = QQ(0)
    elif order == 3:
        correction = QQ(4) / 3
    else:
        raise ArithmeticError(f"unexpected E8+E6 discriminant order {order}")

    height = QQ(z * H * z)
    target = height + correction
    if target not in ZZ:
        raise ArithmeticError(
            f"section norm not integral: h={height}, corr={correction}"
        )
    target = ZZ(target)
    if target < 4 or target % 2:
        raise ArithmeticError(
            f"section norm incompatible with K3 (-2)-section: {target}"
        )

    lattice = IntegralLattice(root)
    it = lattice.enumerate_close_vectors(-dual)
    chosen = None
    examined = 0
    for _ in range(cap):
        shift = vector(ZZ, next(it))
        examined += 1
        w = base + vector(ZZ, list(shift) + [0] * 3)
        norm = ZZ(w * G * w)
        if norm == target:
            chosen = w
            break
        if norm > target:
            break

    if chosen is None:
        raise ArithmeticError(
            f"CVP missed q6 section z={tuple(z)}, target={target}, "
            f"examined={examined}"
        )

    a = ZZ((target - 2) // 2)
    pole = ZZ((target - 4) // 2)
    Q = vector(ZZ, [a, 1] + list(chosen))
    assert Q * ns * Q == -2
    assert Q * ns * F6 == 1
    assert Q * ns * O6 == pole
    assert vector(ZZ, Q[-3:]) == z

    return {
        "section": Q,
        "height": height,
        "correction": correction,
        "class_order": order,
        "pole": pole,
        "target_norm": target,
        "examined": examined,
    }


q = effective_q6_section(zD)
Q = q["section"]

print(
    "Q24Q6RR_SECTION|"
    f"mw={','.join(map(str,zD))}|height={q['height']}|"
    f"corr={q['correction']}|class_order={q['class_order']}|"
    f"PdotO={q['pole']}|cvp={q['examined']}|"
    f"source_degree={Q*ns*F6}|status=PASS_EFFECTIVE_Q6_SECTION",
    flush=True,
)

# ---------------------------------------------------------------------------
# Whole q24 divisor = generic horizontal representative + vertical correction.
# ---------------------------------------------------------------------------

if dD <= 0:
    raise ArithmeticError(f"q24 divisor has nonpositive q6 degree {dD}")

generic = (dD - 1) * O6 + Q
V = D24 - generic

assert V * ns * F6 == 0
assert vector(ZZ, V[-3:]) == vector(ZZ, (0, 0, 0))
assert V[1] == 0

fibre_coeff = ZZ(V[0])
root_coeffs = vector(ZZ, V[2:16])
assert len(root_coeffs) == 14

# Identify E8/E6 connected components from the root Cartan graph.
adj = {
    i: {j for j in range(14) if i != j and root[i, j] != 0}
    for i in range(14)
}
components = []
unused = set(range(14))
while unused:
    seed = min(unused)
    stack = [seed]
    comp = set()
    while stack:
        i = stack.pop()
        if i in comp:
            continue
        comp.add(i)
        stack.extend(adj[i] - comp)
    unused -= comp
    components.append(tuple(sorted(comp)))

component_records = []
for comp in components:
    sub = root.matrix_from_rows_and_columns(comp, comp)
    det = abs(ZZ(sub.det()))
    if len(comp) == 8 and det == 1:
        name = "E8"
    elif len(comp) == 6 and det == 3:
        name = "E6"
    else:
        name = f"rank{len(comp)}_det{det}"
    coeff = [int(root_coeffs[i]) for i in comp]
    component_records.append((name, comp, coeff))

print(
    "Q24Q6RR_VERTICAL|"
    f"fibre={fibre_coeff}|"
    f"roots={','.join(map(str,root_coeffs))}|"
    f"support={sum(bool(x) for x in root_coeffs)}|"
    f"V2={V*ns*V}|status=PASS_EXACT_VERTICAL",
    flush=True,
)

for name, comp, coeff in component_records:
    print(
        "Q24Q6RR_COMPONENT|"
        f"type={name}|nodes={','.join(str(i+1) for i in comp)}|"
        f"coeff={','.join(map(str,coeff))}|status=PASS",
        flush=True,
    )

assert D24 == generic + V

# For degree four, the generic elliptic RR basis is immediate.
generic_basis = None
if dD == 4:
    generic_basis = [
        "1",
        "x",
        "y",
        "m_Q=(y+y_Q)/(x-x_Q)",
    ]
    print(
        "Q24Q6RR_GENERIC|degree=4|dimension=4|"
        "basis=1,x,y,m_Q|m_Q=(y+y_Q)/(x-x_Q)|"
        "status=PASS_SMALL_GENERIC_RR",
        flush=True,
    )
else:
    print(
        "Q24Q6RR_GENERIC|"
        f"degree={dD}|dimension={dD}|"
        "status=NEEDS_GENERIC_RR_BASIS",
        flush=True,
    )

Q_source = Q * Bsimple
D_source = D24 * Bsimple
V_source = V * Bsimple

payload = {
    "schema": "elkies-k3.h92-q24-q6-rr-route.v1",
    "status": (
        "PASS_Q24_DEGREE4_Q6_RR_ROUTE"
        if dD == 4
        else "PASS_Q24_Q6_PULLBACK_ROUTE"
    ),
    "parent": "exact E8+E6/MW3 q6 fibration",
    "degrees_on_q6": {
        "q24_horizontal_P": int(dP),
        "q8_zero": int(dO8),
        "q8_fibre": int(F8 * ns * F6),
        "q24_divisor": int(dD),
    },
    "q24_generic_picard": {
        "degree": int(dD),
        "mw_coordinates_in_q6_simple_basis": list(map(int, zD)),
        "height": str(hD),
        "effective_section": {
            "source_h3_ns": list(map(int, Q_source)),
            "q6_simple_ns": list(map(int, Q)),
            "height": str(q["height"]),
            "local_correction": str(q["correction"]),
            "class_order": int(q["class_order"]),
            "P_dot_O": int(q["pole"]),
        },
        "equivalent_generic_divisor": f"{int(dD-1)}*O6 + Q",
        "generic_rr_basis": generic_basis,
    },
    "vertical_completion_in_q6_simple_frame": {
        "fibre_coefficient": int(fibre_coeff),
        "root_coefficients": list(map(int, root_coeffs)),
        "components": [
            {
                "type": name,
                "nodes_1_based": [i + 1 for i in comp],
                "coefficients": coeff,
            }
            for name, comp, coeff in component_records
        ],
        "square": int(V * ns * V),
        "source_h3_ns": list(map(int, V_source)),
    },
    "q24_divisor": {
        "q6_simple_ns": list(map(int, D24)),
        "source_h3_ns": list(map(int, D_source)),
        "square": int(D24 * ns * D24),
    },
    "next": (
        "Build the global two-dimensional L(D24) on the exact q6 equation: "
        "start from the generic degree-d basis above, impose the E8/E6 "
        "resolved-component ideals encoded by the vertical coefficients, then "
        "derive infinity degree bounds. Eliminate the resulting two generators "
        "to a genus-one pencil and take binary-quartic invariants."
    ),
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24Q6RR_RESULT|"
    f"degree={dD}|mw={','.join(map(str,zD))}|"
    f"vertical_F={fibre_coeff}|vertical_support={sum(bool(x) for x in root_coeffs)}|"
    f"status={payload['status']}",
    flush=True,
)
