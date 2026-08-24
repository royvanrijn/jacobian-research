#!/usr/bin/env sage -python
"""
Exact H3 q6/orbit42 lift:
    q24/orbit85 D12/MW5  ->  q6/orbit42 A11/MW6.

This consumes the newly certified exact q24 D12 model and the already
certified mod-100003 orbit42 marked section.

The historical "q6" neighbor has actual old-fibre degree two.  Therefore,
once P42 is recovered exactly, the new base is simply the marked chord slope.

Stages:
  1. reproduce the modular D12 canonical/actual-twist normalization over QQ;
  2. Hensel-lift the P.O=3 orbit42 section (deg Z = 3) from mod 100003;
  3. rationally reconstruct and verify the section exactly;
  4. compile the chord discriminant fraction-free in QQ[u,m];
  5. extract its squarefree degree-3/4 model by bivariate Yun;
  6. compute binary-quartic invariants/Jacobian;
  7. classify the exact child and require A11 / MW6;
  8. regress section, quartic and Jacobian mod 100003.

The p-adic section residue is checkpointed.  Increase --precision and rerun;
the script automatically resumes the saved p-adic residue when possible.
"""

import argparse
import json
import subprocess
import time
from pathlib import Path

from sage.all import (
    GF, QQ, ZZ, PolynomialRing, Zp, lcm, matrix, sage_eval, vector
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
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
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/local/elkies-k3").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--precision", type=int, default=64, help="p-adic digits")
parser.add_argument("--seed", type=Path, help="optional previous Hensel checkpoint")
parser.add_argument("--rank-only", action="store_true")
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
SCRIPTS = ROOT / "elkies-k3/scripts"
CORE = SCRIPTS / "elliptic_neighbor_compiler.sage"

EXACT_D12 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
BRIDGE = LOCAL / "q24-orbit42-current-equation-bridge.json"
MOD_SECTION = LOCAL / f"q24-orbit42-current-equation-section-mod-{args.prime}.json"
MOD_CHORD = LOCAL / f"q24-orbit42-a11-chord-mod-{args.prime}.json"
DEFAULT_OUTPUT = LOCAL / "q24-d12-to-a11-orbit42-qq.json"
CHECKPOINT = LOCAL / f"q24-orbit42-section-hensel-p{args.prime}.json"

if not EXACT_D12.exists():
    raise SystemExit(f"Missing exact D12 prerequisite: {EXACT_D12}")
if not CORE.exists():
    raise SystemExit(f"Missing compiler core: {CORE}")

# AUTO_MODULAR_PREREQUISITES
def run_seed_producer(script_name, *extra):
    script = SCRIPTS / script_name
    if not script.exists():
        raise SystemExit(f"Missing modular seed producer: {script}")
    command = ["sage", "-python", str(script), *map(str, extra)]
    print(
        "Q24O42QQ_PREREQ|"
        f"script={script_name}|command={' '.join(command)}|status=BEGIN",
        flush=True,
    )
    subprocess.run(command, cwd=str(ROOT), check=True)
    print(
        "Q24O42QQ_PREREQ|"
        f"script={script_name}|status=PASS",
        flush=True,
    )

if not BRIDGE.exists():
    run_seed_producer(
        "extract_h92_q24_orbit42_current_equation_bridge.sage"
    )
    if not BRIDGE.exists():
        raise SystemExit(f"Bridge producer did not create {BRIDGE}")

if not MOD_SECTION.exists():
    run_seed_producer(
        "recover_h92_q24_orbit42_current_equation_section_modp.sage",
        "--prime", args.prime,
    )
    if not MOD_SECTION.exists():
        raise SystemExit(
            f"Orbit42 section producer did not create {MOD_SECTION}"
        )

if not MOD_CHORD.exists():
    run_seed_producer(
        "compile_h92_q24_orbit42_a11_chord_modp.sage",
        "--prime", args.prime,
    )
    if not MOD_CHORD.exists():
        raise SystemExit(
            f"Orbit42 chord producer did not create {MOD_CHORD}"
        )

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{k}={v}" for k, v in fields.items())
    print(
        f"Q24O42QQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


exact = json.loads(EXACT_D12.read_text())
mod_section = json.loads(MOD_SECTION.read_text())
mod_chord = json.loads(MOD_CHORD.read_text())

assert exact["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert mod_section["status"] == "PASS_Q24_ORBIT42_CURRENT_EQUATION_SECTION_MODP"
assert mod_chord["status"] == "PASS_Q24_ORBIT42_A11_CHORD_MODP"

p = ZZ(args.prime)
F = GF(p)

VR = PolynomialRing(QQ, "V")
V = VR.gen()
VF = VR.fraction_field()

RQ = PolynomialRing(QQ, "u")
u = RQ.gen()
KQ = RQ.fraction_field()

RF = PolynomialRing(F, "u")
uf = RF.gen()


def parse_vf(text):
    return VF(sage_eval(str(text), locals={"V": V}))


def red_q(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % p == 0:
        raise ArithmeticError(f"bad reduction denominator at p={p}")
    return F(ZZ(value.numerator())) / F(den)


def red_poly(poly, target_ring=RF):
    poly = RQ(poly)
    return target_ring([red_q(c) for c in poly.list()])


def max_q_bits(polys):
    bits = 0
    for poly in polys:
        for c in poly:
            q = QQ(c)
            bits = max(
                bits,
                abs(ZZ(q.numerator())).nbits(),
                abs(ZZ(q.denominator())).nbits(),
            )
    return bits


# ---------------------------------------------------------------------------
# 1. Reconstruct the exact actual-twist D12 model in the same canonical
#    coordinate used by the modular orbit42 section recovery.
# ---------------------------------------------------------------------------
raw = exact["jacobian_raw"]
Aorig = parse_vf(raw["A"])
Borig = parse_vf(raw["B"])

Jmap = VF(QQ(6912)) * Aorig**3 / (
    VF(QQ(4)) * Aorig**3 + VF(QQ(27)) * Borig**2
)
N = VR(Jmap.numerator())
Den = VR(Jmap.denominator())
g = N.gcd(Den)
if g.degree() > 0:
    N //= g
    Den //= g

if (N.degree(), Den.degree()) != (18, 18):
    raise ArithmeticError(
        f"unexpected exact D12 j degrees {(N.degree(), Den.degree())}, expected (18,18)"
    )

# Unique multiplicity-8 linear factor, found without factorization.
h = VR(Den)
for _ in range(7):
    h = h.gcd(h.derivative())
if h.degree() != 1:
    # Diagnostic fallback only.
    candidates = [
        fac for fac, exponent in Den.factor()
        if int(exponent) == 8 and fac.degree() == 1
    ]
    if len(candidates) != 1:
        raise ArithmeticError(
            f"could not isolate unique exact I8* factor; repeated gcd degree={h.degree()}"
        )
    h = candidates[0]
h = h.monic()
r = -QQ(h[0]) / QQ(h[1])

RT = PolynomialRing(QQ, "T")
T = RT.gen()


def invpoly(poly):
    poly = VR(poly)
    return RT(sum(
        QQ(poly[i]) * (QQ(r) * T + 1)**i * T**(18-i)
        for i in range(poly.degree() + 1)
    ))


Pj = invpoly(N)
Qj = invpoly(Den)
g = Pj.gcd(Qj)
if g.degree() > 0:
    Pj //= g
    Qj //= g

if (Pj.degree(), Qj.degree()) != (18, 10):
    raise ArithmeticError(
        f"unexpected inverted j degrees {(Pj.degree(), Qj.degree())}"
    )

lc = Qj.leading_coefficient()
Pj /= lc
Qj /= lc
if not Qj.is_monic():
    raise ArithmeticError("exact inverted j denominator is not monic")

center = -QQ(Qj[9]) / QQ(10)

RS = PolynomialRing(QQ, "S")
S = RS.gen()
P1 = RS(Pj(S + center))
Q1 = RS(Qj(S + center))
if not (Q1.is_monic() and Q1[9] == 0 and Q1[8] and Q1[7]):
    raise ArithmeticError("exact D12 centered j normalization failed")

base_scale = QQ(Q1[7]) / QQ(Q1[8])
P2 = RQ(P1(base_scale * u))
Q2 = RQ(Q1(base_scale * u))
lc2 = Q2.leading_coefficient()
P2 /= lc2
Q2 /= lc2
if not (
    Q2.is_monic()
    and Q2.degree() == 10
    and Q2[9] == 0
    and Q2[8] == Q2[7]
):
    raise ArithmeticError("exact D12 scaled j normalization failed")


def monic_power_root(poly, exponent):
    """Recover the monic exact e-th root without factorization."""
    poly = RQ(poly)
    if not poly:
        raise ArithmeticError("zero has no power root")
    lc = QQ(poly.leading_coefficient())
    f = RQ(poly / lc)
    if f.degree() == 0:
        if f != 1:
            raise ArithmeticError("constant monic normalization failed")
        return RQ.one()
    if f.degree() % exponent:
        raise ArithmeticError(
            f"degree {f.degree()} not divisible by exponent {exponent}"
        )
    degree = f.degree() // exponent
    coeffs = [QQ(0)] * degree + [QQ(1)]
    for j in range(degree - 1, -1, -1):
        current = RQ(coeffs)
        target_degree = (exponent - 1) * degree + j
        known = (current**exponent)[target_degree]
        coeffs[j] = QQ(f[target_degree] - known) / QQ(exponent)
    root = RQ(coeffs)
    if root**exponent != f:
        raise ArithmeticError(f"exact monic {exponent}-th root reconstruction failed")
    return root


a = monic_power_root(P2, 3)
H2 = RQ(P2 - QQ(1728) * Q2)
b = monic_power_root(H2, 2)

if (a.degree(), b.degree()) != (6, 9):
    raise ArithmeticError(
        f"unexpected exact canonical D12 degrees {(a.degree(), b.degree())}"
    )
if (a**3 - b**2).degree() != 10:
    raise ArithmeticError("exact canonical D12 discriminant profile failed")

vmap = KQ(r) + KQ(1) / (KQ(base_scale) * KQ(u) + KQ(center))


def eval_v_rational(value, argument):
    value = VF(value)
    return KQ(value.numerator()(argument)) / KQ(value.denominator()(argument))


Aeval = eval_v_rational(Aorig, vmap)
Beval = eval_v_rational(Borig, vmap)
Acan = KQ(-QQ(3) * a)
Bcan = KQ(QQ(2) * b)

cA = Aeval / Acan
cB = Beval / Bcan
wfun = cB / cA
if cA != wfun**2 or cB != wfun**3:
    raise ArithmeticError("exact D12 Weierstrass scaling relation failed")

wn = RQ(wfun.numerator())
wd = RQ(wfun.denominator())
wn_lc = QQ(wn.leading_coefficient())
wd_lc = QQ(wd.leading_coefficient())
wn_monic = RQ(wn / wn_lc)
wd_monic = RQ(wd / wd_lc)
sn = monic_power_root(wn_monic, 2)
sd = monic_power_root(wd_monic, 2)
square_part = KQ(sn) / KQ(sd)

Dfun = wfun / square_part**2
Dnum = RQ(Dfun.numerator())
Dden = RQ(Dfun.denominator())
if Dnum.degree() > 0 or Dden.degree() > 0:
    raise ArithmeticError("exact D12 twist did not reduce to a constant")
twist = QQ(Dnum[0]) / QQ(Dden[0])

At = RQ(-QQ(3) * twist**2 * a)
Bt = RQ(QQ(2) * twist**3 * b)

if Aeval != KQ(square_part**4) * KQ(At):
    raise ArithmeticError("exact D12 A twist reconstruction failed")
if Beval != KQ(square_part**6) * KQ(Bt):
    raise ArithmeticError("exact D12 B twist reconstruction failed")

mod_model = mod_section["actual_twist_model"]
if red_q(r) != F(mod_model["I8star_root"]):
    raise ArithmeticError("exact I8* root does not reduce to modular normalization")
if red_q(center) != F(mod_model["center"]):
    raise ArithmeticError("exact center does not reduce to modular normalization")
if red_q(base_scale) != F(mod_model["base_scale"]):
    raise ArithmeticError("exact base scale does not reduce to modular normalization")
if red_q(twist) != F(mod_model["twist"]):
    raise ArithmeticError("exact twist does not reduce to modular normalization")

At_mod = red_poly(At)
Bt_mod = red_poly(Bt)
if At_mod != RF(mod_model["A"]) or Bt_mod != RF(mod_model["B"]):
    raise ArithmeticError("exact actual-twist D12 model misses modular model")

log(
    "MODEL",
    Adeg=At.degree(),
    Bdeg=Bt.degree(),
    max_bits=max_q_bits((At, Bt)),
    I8star=r,
    status="PASS_EXACT_ACTUAL_TWIST_D12",
)


# ---------------------------------------------------------------------------
# 2. Select the known modular orbit42 section that produces the A11 chord.
# ---------------------------------------------------------------------------
a11_results = [row for row in mod_chord["results"] if row.get("has_I12_A11")]
if not a11_results:
    raise ArithmeticError("modular orbit42 chord artifact contains no A11 target")

preferred_indices = [int(row["target_index"]) for row in a11_results]
targets_by_index = {
    int(row["index"]): row
    for row in mod_section["orbit42_section_candidates"]
}
target_index = preferred_indices[0]
if target_index not in targets_by_index:
    raise ArithmeticError(f"modular A11 target index {target_index} is missing")
target = targets_by_index[target_index]

Xf = RF(target["X"])
Yf = RF(target["Y"])
Zf = RF(target["Z"])
dx, dy, dz = tuple(map(int, target["degrees"]))

if Zf.degree() != dz or Zf.leading_coefficient() != 1:
    raise ArithmeticError("modular P42 denominator is not monic")
if Yf**2 != Xf**3 + At_mod * Xf * Zf**4 + Bt_mod * Zf**6:
    raise ArithmeticError("selected modular P42 is not on exact model reduction")

log(
    "MODULAR_SECTION",
    target=target_index,
    Xdeg=dx,
    Ydeg=dy,
    Zdeg=dz,
    PdotO=3,
    status="PASS",
)


# ---------------------------------------------------------------------------
# 3. Nonsingular p-adic Hensel lift of the section coefficients.
# ---------------------------------------------------------------------------
unknown_count = dz + (dx + 1) + (dy + 1)
equation_degree = max(
    2 * dy,
    3 * dx,
    At.degree() + dx + 4 * dz,
    Bt.degree() + 6 * dz,
)


def coeffs(poly, degree):
    return [
        poly[i] if i <= poly.degree() else poly.base_ring()(0)
        for i in range(degree + 1)
    ]


def identity_and_jacobian(ring, z_values, x_values, y_values):
    uu = ring.gen()
    base = ring.base_ring()
    z = sum(base(z_values[i]) * uu**i for i in range(dz)) + uu**dz
    xpoly = sum(base(x_values[i]) * uu**i for i in range(dx + 1))
    ypoly = sum(base(y_values[i]) * uu**i for i in range(dy + 1))
    aa = ring(At)
    bb = ring(Bt)

    identity = ypoly**2 - xpoly**3 - aa * xpoly * z**4 - bb * z**6

    derivatives = []
    for i in range(dz):
        derivatives.append(
            (-4 * aa * xpoly * z**3 - 6 * bb * z**5) * uu**i
        )
    for i in range(dx + 1):
        derivatives.append(
            (-3 * xpoly**2 - aa * z**4) * uu**i
        )
    for i in range(dy + 1):
        derivatives.append(2 * ypoly * uu**i)

    residue = vector(base, coeffs(identity, equation_degree))
    jacobian = matrix(base, [
        [coeffs(derivative, equation_degree)[row] for derivative in derivatives]
        for row in range(equation_degree + 1)
    ])
    return residue, jacobian


z0 = coeffs(Zf, dz - 1)
x0 = coeffs(Xf, dx)
y0 = coeffs(Yf, dy)

residue_f, jacobian_f = identity_and_jacobian(RF, z0, x0, y0)
if residue_f:
    raise ArithmeticError("modular Hensel seed has nonzero residual")
rank = jacobian_f.rank()
if rank != unknown_count:
    raise ArithmeticError(
        f"orbit42 section Jacobian rank={rank}, unknowns={unknown_count}"
    )
rows = tuple(jacobian_f.transpose().pivots())
if len(rows) != unknown_count:
    raise ArithmeticError("failed to select a square Hensel Jacobian")
square_f = jacobian_f.matrix_from_rows(rows)
if not square_f.is_invertible():
    raise ArithmeticError("selected Hensel Jacobian is singular")

log(
    "HENSEL_JACOBIAN",
    equations=equation_degree + 1,
    unknowns=unknown_count,
    rank=rank,
    selected_rows=len(rows),
    status="PASS_NONSINGULAR",
)

if args.rank_only:
    raise SystemExit(0)

target_precision = int(args.precision)
if target_precision < 2:
    raise ValueError("--precision must be at least 2")

padic = Zp(p, prec=target_precision)
RP = PolynomialRing(padic, "u")

values = [padic(ZZ(v)) for v in z0 + x0 + y0]
seed_precision = 1

seed_path = args.seed.resolve() if args.seed else CHECKPOINT
if seed_path.exists():
    try:
        seed = json.loads(seed_path.read_text())
        if (
            ZZ(seed.get("prime")) == p
            and int(seed.get("target_index", -1)) == target_index
            and isinstance(seed.get("residues"), list)
            and len(seed["residues"]) == len(values)
        ):
            sp = int(seed.get("precision", 1))
            if 1 < sp < target_precision:
                values = [padic(ZZ(v)) for v in seed["residues"]]
                seed_precision = sp
                log(
                    "HENSEL_RESUME",
                    seed=seed_path,
                    seed_precision=sp,
                    target_precision=target_precision,
                    status="PASS",
                )
    except Exception as exc:
        log(
            "HENSEL_RESUME",
            seed=seed_path,
            status="IGNORED_INVALID",
            reason=f"{type(exc).__name__}:{exc}",
        )


def valuation_floor(vector_values):
    vals = [v.valuation() for v in vector_values if v]
    return target_precision if not vals else min(vals)


for iteration in range(1, 2 * target_precision + 4):
    residual_p, derivative_p = identity_and_jacobian(
        RP,
        values[:dz],
        values[dz:dz + dx + 1],
        values[dz + dx + 1:],
    )
    valuation = valuation_floor(residual_p)
    log(
        "HENSEL",
        iteration=iteration,
        residual_valuation=valuation,
        target=target_precision,
    )
    if valuation >= target_precision:
        break

    correction = derivative_p.matrix_from_rows(rows).solve_right(
        -vector(padic, [residual_p[row] for row in rows])
    )
    values = [v + delta for v, delta in zip(values, correction)]
else:
    raise ArithmeticError("orbit42 p-adic Newton iteration did not converge")

modulus = p**target_precision


def reconstruct_one(value):
    try:
        return QQ(ZZ(value.lift()).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError):
        return None


reconstructed = [reconstruct_one(v) for v in values]

checkpoint_payload = {
    "schema": "elkies-k3.h3-q24-orbit42-section-hensel.v1",
    "status": "PASS_PADIC_RESIDUE",
    "prime": int(p),
    "precision": target_precision,
    "seed_precision": int(seed_precision),
    "target_index": target_index,
    "degrees": [dx, dy, dz],
    "jacobian_rank": int(rank),
    "selected_rows": [int(row) for row in rows],
    "residues": [str(ZZ(v.lift())) for v in values],
    "reconstruction_complete": not any(v is None for v in reconstructed),
}
CHECKPOINT.write_text(json.dumps(checkpoint_payload, indent=2, sort_keys=True) + "\n")
log(
    "HENSEL_CHECKPOINT",
    file=CHECKPOINT,
    precision=target_precision,
    complete=int(checkpoint_payload["reconstruction_complete"]),
    status="SAVED",
)

if any(v is None for v in reconstructed):
    print(
        "Q24O42QQ_NEEDS_MORE_PRECISION|"
        f"current={target_precision}|checkpoint={CHECKPOINT}|"
        f"suggested={target_precision*2}|status=NEEDS_MORE_PRECISION",
        flush=True,
    )
    raise SystemExit(2)

Z = RQ(reconstructed[:dz] + [QQ(1)])
X = RQ(reconstructed[dz:dz + dx + 1])
Y = RQ(reconstructed[dz + dx + 1:])

if Y**2 != X**3 + At * X * Z**4 + Bt * Z**6:
    print(
        "Q24O42QQ_NEEDS_MORE_PRECISION|"
        f"current={target_precision}|checkpoint={CHECKPOINT}|"
        f"suggested={target_precision*2}|reason=EXACT_IDENTITY_FAILED|"
        "status=NEEDS_MORE_PRECISION",
        flush=True,
    )
    raise SystemExit(2)

if red_poly(X) != Xf or red_poly(Y) != Yf or red_poly(Z) != Zf:
    raise ArithmeticError("exact orbit42 section does not reduce to modular seed")

section_bits = max_q_bits((X, Y, Z))
log(
    "SECTION",
    Xdeg=X.degree(),
    Ydeg=Y.degree(),
    Zdeg=Z.degree(),
    max_bits=section_bits,
    identity=1,
    modp=1,
    status="PASS_EXACT_ORBIT42_SECTION",
)


# ---------------------------------------------------------------------------
# 4. Exact chord discriminant in QQ[u,m], then fraction-free Yun in u.
# ---------------------------------------------------------------------------
Pm = PolynomialRing(QQ, "m")
m = Pm.gen()
Km = Pm.fraction_field()
RUm = PolynomialRing(Km, "u")
uk = RUm.gen()

BUV = PolynomialRing(QQ, names=("ub", "mb"), order="degrevlex")
ub, mb = BUV.gens()


def uq_to_bi(poly):
    poly = RQ(poly)
    return sum(QQ(c) * ub**i for i, c in enumerate(poly.list()) if c)


Xb = uq_to_bi(X)
Yb = uq_to_bi(Y)
Zb = uq_to_bi(Z)
Ab = uq_to_bi(At)

raw_bi = BUV(
    mb**4 * Zb**4
    - QQ(6) * Xb * mb**2 * Zb**2
    - QQ(8) * Yb * mb * Zb
    - QQ(3) * Xb**2
    - QQ(4) * Ab * Zb**4
)

raw_u_degree = int(raw_bi.degree(ub))
raw_m_degree = int(raw_bi.degree(mb))
log(
    "CHORD_RAW",
    Udeg=raw_u_degree,
    Mdeg=raw_m_degree,
    status="PASS",
)

# U-content over QQ[m].
coeff_by_u = []
for iu in range(raw_u_degree + 1):
    terms = {}
    for (eu, em), c in raw_bi.dict().items():
        if eu == iu and c:
            terms[em] = QQ(c)
    coeff_by_u.append(
        Pm([
            terms.get(j, QQ(0))
            for j in range(max(terms.keys(), default=-1) + 1)
        ]) if terms else Pm.zero()
    )

nonzero = [c for c in coeff_by_u if c]
content = nonzero[0]
for c in nonzero[1:]:
    content = content.gcd(c)
    if content.degree() == 0:
        break
if content:
    content /= content.leading_coefficient()

primitive = BUV.zero()
for iu, c in enumerate(coeff_by_u):
    if not c:
        continue
    q, rem = c.quo_rem(content)
    if rem:
        raise ArithmeticError("chord U-content division failed")
    for em, aij in enumerate(q.list()):
        if aij:
            primitive += QQ(aij) * ub**iu * mb**em

Gbi = primitive.gcd(primitive.derivative(ub))


def bi_quo(left, right, label):
    q, rem = BUV(left).quo_rem(BUV(right))
    if rem:
        raise ArithmeticError(f"{label}: bivariate quotient has remainder")
    return BUV(q)


Wbi = bi_quo(primitive, Gbi, "chord-yun-initial")
odd_bi = BUV.one()
square_bi = BUV.one()
multiplicities = []
mult = 1

while Wbi.degree(ub) > 0:
    Yg = Wbi.gcd(Gbi)
    Zi = bi_quo(Wbi, Yg, f"chord-yun-z-{mult}")
    if Zi.degree(ub) > 0:
        multiplicities.append(mult)
        if mult % 2:
            odd_bi *= Zi
        if mult // 2:
            square_bi *= Zi**(mult // 2)
    Wbi = Yg
    Gbi = bi_quo(Gbi, Yg, f"chord-yun-g-{mult}")
    mult += 1

if Gbi.degree(ub) > 0:
    raise ArithmeticError("chord Yun left a residual U factor")


def bi_to_rum(poly):
    by_u = {}
    for (eu, em), c in poly.dict().items():
        by_u.setdefault(eu, {})[em] = QQ(c)
    if not by_u:
        return RUm.zero()
    out = []
    for iu in range(max(by_u) + 1):
        terms = by_u.get(iu, {})
        out.append(
            Km(sum(QQ(c) * m**em for em, c in terms.items()))
        )
    return RUm(out)


odd_u = bi_to_rum(odd_bi)
square_u = bi_to_rum(square_bi)
odd_u = RUm(odd_u / odd_u.leading_coefficient())
square_u = RUm(square_u / square_u.leading_coefficient())


def bi_to_rum_full(poly):
    return bi_to_rum(poly)


raw_k = bi_to_rum_full(raw_bi)
product = RUm(odd_u * square_u**2)
q, rem = raw_k.quo_rem(product)
if rem or q.degree() > 0:
    raise ArithmeticError(
        f"chord square-class quotient is not a base-field unit; "
        f"qdeg={q.degree()}, rem={bool(rem)}"
    )
unit = Km(q[0])
quartic = RUm(unit * odd_u)
square_factor = RUm(square_u)

if raw_k != quartic * square_factor**2:
    raise ArithmeticError("exact orbit42 chord square-class reconstruction failed")

quartic_degree = int(quartic.degree())
if quartic_degree not in (3, 4):
    raise ArithmeticError(
        f"orbit42 squarefree chord degree={quartic_degree}, expected 3/4"
    )

log(
    "QUARTIC",
    raw_degree=raw_u_degree,
    square_degree=square_factor.degree(),
    degree=quartic_degree,
    multiplicities=",".join(map(str, multiplicities)),
    status="PASS_EXACT",
)


# ---------------------------------------------------------------------------
# 5. Binary quartic invariants and exact A11 child classification.
# ---------------------------------------------------------------------------
core = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), core)
binary_quartic_invariants = core["binary_quartic_invariants"]
kodaira_data_from_short_orders = core["kodaira_data_from_short_orders"]

I, J = binary_quartic_invariants(quartic)
jacA = Km(-27) * Km(I)
jacB = Km(-27) * Km(J)

Delta = Km(-16) * (Km(4) * jacA**3 + Km(27) * jacB**2)
if not Delta:
    raise ArithmeticError("orbit42 child Jacobian is singular")

log("JACOBIAN", status="PASS")


def squarefree_part(poly):
    poly = Pm(poly)
    if poly.degree() <= 0:
        return Pm.one()
    g = poly.gcd(poly.derivative())
    q, rem = poly.quo_rem(g)
    if rem:
        raise ArithmeticError("squarefree support quotient failed")
    return Pm(q)


support = Pm.one()
for poly in (
    Pm(jacA.denominator()),
    Pm(jacB.denominator()),
    Pm(Delta.numerator()),
    Pm(Delta.denominator()),
):
    sf = squarefree_part(poly)
    common = support.gcd(sf)
    q, rem = sf.quo_rem(common)
    if rem:
        raise ArithmeticError("support LCM quotient failed")
    support *= q
if support.degree() > 0:
    support /= support.leading_coefficient()

log(
    "CHILD_SUPPORT",
    degree=support.degree(),
    method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
    status="BEGIN",
)
support_factorization = support.factor() if support.degree() > 0 else ()
factors = tuple(f for f, unused in support_factorization)
log(
    "CHILD_SUPPORT",
    degree=support.degree(),
    factors=len(factors),
    method="SINGLE_SQUAREFREE_SUPPORT_FACTOR",
    status="PASS",
)


def poly_order(poly, factor):
    poly = Pm(poly)
    if not poly:
        raise ArithmeticError("valuation of zero polynomial")
    order = 0
    while poly.degree() >= factor.degree():
        q, rem = poly.quo_rem(factor)
        if rem:
            break
        order += 1
        poly = Pm(q)
    return order


def rf_order(value, factor):
    value = Km(value)
    return (
        poly_order(value.numerator(), factor)
        - poly_order(value.denominator(), factor)
    )


scaling_unit = Km.one()
finite_places = []
for factor in sorted(factors, key=str):
    va = rf_order(jacA, factor)
    vb = rf_order(jacB, factor)
    vd = rf_order(Delta, factor)
    scaling = min(va // 4, vb // 6)
    scaling_unit *= Km(factor)**(-scaling)
    finite_places.append({
        "factor": factor,
        "raw_orders": (va, vb, vd),
        "scaling": scaling,
    })

minimal_A_rf = jacA * scaling_unit**4
minimal_B_rf = jacB * scaling_unit**6
minimal_D_rf = Delta * scaling_unit**12

if any(
    v.denominator() != 1
    for v in (minimal_A_rf, minimal_B_rf, minimal_D_rf)
):
    raise ArithmeticError("orbit42 finite minimization left denominators")

minimal_A = Pm(minimal_A_rf.numerator())
minimal_B = Pm(minimal_B_rf.numerator())
minimal_D = Pm(minimal_D_rf.numerator())

root_rank = ZZ(0)
root_det = ZZ(1)
euler = ZZ(0)
finite_fibres = []

for place in finite_places:
    factor = place["factor"]
    orders = (
        poly_order(minimal_A, factor),
        poly_order(minimal_B, factor),
        poly_order(minimal_D, factor),
    )
    place["minimal_orders"] = orders
    if orders[2] == 0:
        continue
    rank_i, euler_i, det_i, symbol = kodaira_data_from_short_orders(*orders)
    degree_i = ZZ(factor.degree())
    root_rank += degree_i * ZZ(rank_i)
    euler += degree_i * ZZ(euler_i)
    root_det *= ZZ(det_i)**degree_i
    finite_fibres.append({
        "factor": str(factor),
        "degree": int(degree_i),
        "minimal_orders": list(map(int, orders)),
        "kodaira": symbol,
    })

infinity_raw = tuple(
    -poly.degree() for poly in (minimal_A, minimal_B, minimal_D)
)
infinity_scaling = min(infinity_raw[0] // 4, infinity_raw[1] // 6)
infinity_orders = tuple(
    infinity_raw[i] - (4, 6, 12)[i] * infinity_scaling
    for i in range(3)
)
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    rr, ee, dd, infinity_kind = kodaira_data_from_short_orders(*infinity_orders)
    root_rank += ZZ(rr)
    euler += ZZ(ee)
    root_det *= ZZ(dd)

if (root_rank, root_det, euler) != (11, 12, 24):
    raise ArithmeticError(
        f"exact orbit42 child is not A11: "
        f"rank={root_rank}, det={root_det}, euler={euler}"
    )

log(
    "CHILD",
    root_rank=root_rank,
    root_det=root_det,
    euler=euler,
    MW=6,
    infinity=infinity_kind,
    status="PASS_A11",
)


# ---------------------------------------------------------------------------
# 6. Regression to the already certified modular orbit42 chord.
# ---------------------------------------------------------------------------
Fm = PolynomialRing(F, "m")
mf = Fm.gen()
Kmf = Fm.fraction_field()


def red_km(value):
    value = Km(value)
    num = Pm(value.numerator())
    den = Pm(value.denominator())
    nmod = Fm([red_q(c) for c in num.list()])
    dmod = Fm([red_q(c) for c in den.list()])
    return Kmf(nmod) / Kmf(dmod)


selected_mod = next(
    row for row in a11_results
    if int(row["target_index"]) == target_index
)

expected_quartic = [
    Kmf(sage_eval(str(text), locals={"m": mf}))
    for text in selected_mod["reduced_quartic_coefficients_low_to_high"]
]
actual_quartic = [red_km(quartic[i]) for i in range(quartic_degree + 1)]

if actual_quartic != expected_quartic:
    raise ArithmeticError("exact orbit42 quartic misses modular A11 chord")
if red_km(jacA) != Kmf(sage_eval(str(selected_mod["jacobian_A"]), locals={"m": mf})):
    raise ArithmeticError("exact orbit42 jacobian A misses modular chord")
if red_km(jacB) != Kmf(sage_eval(str(selected_mod["jacobian_B"]), locals={"m": mf})):
    raise ArithmeticError("exact orbit42 jacobian B misses modular chord")

log(
    "MODULAR_REGRESSION",
    section=1,
    quartic=1,
    jacobian=1,
    prime=p,
    status="PASS",
)


# ---------------------------------------------------------------------------
# 7. Artifact.
# ---------------------------------------------------------------------------
def qlist(poly):
    return [str(v) for v in poly.list()]


payload = {
    "schema": "elkies-k3.h3-q24-d12-to-a11-orbit42-qq.v1",
    "status": "PASS_EXACT_Q24_D12_TO_A11_ORBIT42",
    "inputs": {
        "exact_d12": str(EXACT_D12.relative_to(ROOT)),
        "modular_section": str(MOD_SECTION.relative_to(ROOT)),
        "modular_chord": str(MOD_CHORD.relative_to(ROOT)),
    },
    "normalization": {
        "I8star_root": str(r),
        "center": str(center),
        "base_scale": str(base_scale),
        "twist": str(twist),
        "A_coefficients_low_to_high": qlist(At),
        "B_coefficients_low_to_high": qlist(Bt),
    },
    "orbit42_section": {
        "target_index_mod_100003": target_index,
        "P_dot_O": 3,
        "X_coefficients_low_to_high": qlist(X),
        "Y_coefficients_low_to_high": qlist(Y),
        "Z_coefficients_low_to_high": qlist(Z),
        "degrees": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "max_coefficient_bits": int(section_bits),
        "hensel_precision": target_precision,
    },
    "quartic": {
        "degree": quartic_degree,
        "raw_degree_in_old_base": raw_u_degree,
        "square_factor_degree": int(square_factor.degree()),
        "coefficients_in_u_low_to_high": [str(v) for v in quartic.list()],
        "I": str(I),
        "J": str(J),
    },
    "jacobian_raw": {
        "A": str(jacA),
        "B": str(jacB),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": qlist(minimal_A),
        "minimal_B_coefficients_low_to_high": qlist(minimal_B),
        "minimal_discriminant_coefficients_low_to_high": qlist(minimal_D),
        "finite_fibres": finite_fibres,
        "infinity_orders": list(map(int, infinity_orders)),
        "infinity_kind": infinity_kind,
        "root_rank": int(root_rank),
        "root_determinant": int(root_det),
        "euler_number": int(euler),
        "MW_rank_if_rho19": 6,
    },
    "verification": {
        "exact_d12_normalization": True,
        "exact_section_identity": True,
        "mod_100003_section_regression": True,
        "exact_chord_square_class": True,
        "mod_100003_quartic_regression": True,
        "mod_100003_jacobian_regression": True,
        "exact_A11_classification": True,
    },
}

OUT = args.output.resolve() if args.output else DEFAULT_OUTPUT
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24O42QQ_RESULT|"
    f"section_deg={X.degree()},{Y.degree()},{Z.degree()}|"
    f"quartic={quartic_degree}|"
    f"root_rank={root_rank}|root_det={root_det}|euler={euler}|MW=6|"
    "status=PASS_EXACT_Q24_D12_TO_A11_ORBIT42",
    flush=True,
)
