#!/usr/bin/env sage -python
"""
CRT/rational-reconstruct the q32 D12 short Weierstrass model and the canonical
spinor x-coordinate over QQ.

Inputs:
  artifacts/local/elkies-k3/q32-signature-mod-*.json
  artifacts/local/elkies-k3/q32-pointed-spinor-weierstrass-anchor.json

Strategy:
  1. Reconstruct A(V), B(V), x_P(V) coefficientwise from canonical monic-
     denominator modular rational functions.
  2. Try a 10-prime reconstruction and use the 11th as a genuine holdout.
  3. Reconstruct again from all primes.
  4. Recover y_P over QQ(V) from y^2=x^3+A*x+B instead of CRT'ing y, because
     the modular square-root sign is intentionally arbitrary prime-to-prime.
  5. Verify the exact Weierstrass identity and reduction to every modular prime.

No expensive q32 RR compiler is rerun.
"""

import json, math
from pathlib import Path
from sage.all import GF, Integers, PolynomialRing, QQ, ZZ

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

POINTED = LOCAL / "q32-pointed-spinor-weierstrass-anchor.json"
if not POINTED.exists():
    raise SystemExit(f"missing {POINTED}")
pointed = json.loads(POINTED.read_text())
if pointed.get("status") != "PASS_POINTED_Q32_D12_SPINOR_MARKING":
    raise SystemExit("pointed spinor anchor is not passing")

pointed_by_p = {int(r["prime"]): r for r in pointed["primes"]}

records = []
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig = json.loads(path.read_text())
        p = int(sig["prime"])
    except Exception:
        continue
    if sig.get("status") != "PASS_Q32_MODP_SIGNATURE":
        continue
    if p not in pointed_by_p:
        continue
    records.append((p, sig, pointed_by_p[p]))
records.sort(key=lambda z: z[0])

if len(records) < 3:
    raise SystemExit(f"need at least 3 compatible primes, have {len(records)}")

print(
    "Q32QQ_INPUT|"
    f"primes={','.join(str(p) for p,_,_ in records)}|"
    f"count={len(records)}|status=PASS",
    flush=True,
)


def rf_record(rec, key):
    if key == "A":
        return rec[1]["jacobian_A"]
    if key == "B":
        return rec[1]["jacobian_B"]
    if key == "x":
        return rec[2]["marked_section_x"]
    if key == "y_mod":
        return rec[2]["marked_section_y"]
    raise KeyError(key)


def degree_signature(key, recs):
    return sorted(set(
        (int(rf_record(r, key)["num_degree"]), int(rf_record(r, key)["den_degree"]))
        for r in recs
    ))


for key in ("A", "B", "x", "y_mod"):
    ds = degree_signature(key, records)
    print(f"Q32QQ_DEGREES|field={key}|patterns={ds}|status={'PASS' if len(ds)==1 else 'MISMATCH'}", flush=True)
    if len(ds) != 1:
        raise SystemExit(f"inconsistent modular degree pattern for {key}: {ds}")


def crt_residues(values):
    r = ZZ(0)
    m = ZZ(1)
    for a, p in values:
        p = ZZ(p)
        a = ZZ(a) % p
        t = ((a - r) % p) * (m % p).inverse_mod(p) % p
        r += m * t
        m *= p
    return r % m, m


def rational_reconstruct(values):
    residue, modulus = crt_residues(values)
    try:
        q = QQ(Integers(modulus)(residue).rational_reconstruction())
    except Exception:
        return None, modulus
    return q, modulus


def reconstruct_rf(key, recs):
    nr, dr = degree_signature(key, recs)[0]
    nums = []
    dens = []
    mod = None
    for i in range(nr + 1):
        vals = []
        for rec in recs:
            p = rec[0]
            arr = rf_record(rec, key)["num"]
            vals.append((arr[i] if i < len(arr) else 0, p))
        q, mod = rational_reconstruct(vals)
        if q is None:
            return None, mod, f"num[{i}]"
        nums.append(q)
    for i in range(dr + 1):
        vals = []
        for rec in recs:
            p = rec[0]
            arr = rf_record(rec, key)["den"]
            vals.append((arr[i] if i < len(arr) else 0, p))
        q, mod = rational_reconstruct(vals)
        if q is None:
            return None, mod, f"den[{i}]"
        dens.append(q)

    R = PolynomialRing(QQ, "V")
    n = R(nums)
    d = R(dens)
    if not d:
        return None, mod, "zero-denominator"

    # Canonical exact normalization matching the modular representation.
    lc = d.leading_coefficient()
    n /= lc
    d /= lc
    K = R.fraction_field()
    return K(n) / K(d), mod, None


def reduce_q(q, p, F):
    q = QQ(q)
    den = ZZ(q.denominator()) % p
    if not den:
        raise ZeroDivisionError
    return F(ZZ(q.numerator()) % p) / F(den)


def norm_mod_rf(value, p):
    RQ = value.parent().ring()
    F = GF(p)
    RF = PolynomialRing(F, "V")
    K = RF.fraction_field()
    nQ = RQ(value.numerator())
    dQ = RQ(value.denominator())
    n = RF([reduce_q(c, p, F) for c in nQ.list()])
    d = RF([reduce_q(c, p, F) for c in dQ.list()])
    if not d:
        raise ZeroDivisionError
    lc = d.leading_coefficient()
    n /= lc
    d /= lc
    return K(n) / K(d)


def record_to_mod_rf(rec, key):
    p = rec[0]
    F = GF(p)
    R = PolynomialRing(F, "V")
    K = R.fraction_field()
    rr = rf_record(rec, key)
    n = R([F(int(v)) for v in rr["num"]])
    d = R([F(int(v)) for v in rr["den"]])
    return K(n) / K(d)


def matches_record(value, rec, key):
    try:
        return norm_mod_rf(value, rec[0]) == record_to_mod_rf(rec, key)
    except ZeroDivisionError:
        return False


def bits_of_rf(v):
    R = v.parent().ring()
    vals = list(R(v.numerator())) + list(R(v.denominator()))
    nb = max(
        (abs(ZZ(QQ(c).numerator())).nbits() for c in vals if c),
        default=0,
    )
    db = max(
        (abs(ZZ(QQ(c).denominator())).nbits() for c in vals if c),
        default=0,
    )
    return int(nb), int(db)


def reconstruct_bundle(recs, label):
    out = {}
    modulus = None
    for key in ("A", "B", "x"):
        v, modulus, failure = reconstruct_rf(key, recs)
        if v is None:
            print(
                f"Q32QQ_RECON|set={label}|field={key}|"
                f"modulus_bits={ZZ(modulus).nbits()}|failure={failure}|status=NEED_MORE_PRIMES",
                flush=True,
            )
            return None, modulus
        nb, db = bits_of_rf(v)
        print(
            f"Q32QQ_RECON|set={label}|field={key}|"
            f"modulus_bits={ZZ(modulus).nbits()}|"
            f"max_num_bits={nb}|max_den_bits={db}|status=PASS",
            flush=True,
        )
        out[key] = v
    return out, modulus


# First: genuine holdout if possible.
holdout_pass = False
holdout_prime = None
if len(records) >= 4:
    train = records[:-1]
    test = records[-1]
    candidate10, mod10 = reconstruct_bundle(train, "holdout-train")
    holdout_prime = test[0]
    if candidate10 is not None:
        checks = {
            key: matches_record(candidate10[key], test, key)
            for key in ("A", "B", "x")
        }
        holdout_pass = all(checks.values())
        print(
            "Q32QQ_HOLDOUT|"
            f"prime={holdout_prime}|"
            + "|".join(f"{k}={int(v)}" for k,v in checks.items())
            + f"|status={'PASS_INDEPENDENT_PRIME' if holdout_pass else 'FAIL_HOLDOUT'}",
            flush=True,
        )

# Full reconstruction.
exact, modulus = reconstruct_bundle(records, "all")
if exact is None:
    print(
        "Q32QQ_RESULT|"
        f"primes={len(records)}|modulus_bits={ZZ(modulus).nbits()}|"
        "status=NEED_MORE_MODULAR_PRECISION",
        flush=True,
    )
    raise SystemExit(0)

A, B, x = exact["A"], exact["B"], exact["x"]
K = A.parent()
R = K.ring()

# Every modular record must match exactly.
all_reduce = True
for rec in records:
    checks = {key: matches_record(exact[key], rec, key) for key in ("A", "B", "x")}
    ok = all(checks.values())
    all_reduce &= ok
    print(
        "Q32QQ_REDUCTION|"
        f"prime={rec[0]}|"
        + "|".join(f"{k}={int(v)}" for k,v in checks.items())
        + f"|status={'PASS' if ok else 'FAIL'}",
        flush=True,
    )
if not all_reduce:
    raise ArithmeticError("full reconstructed A/B/x does not reduce to all inputs")


def qq_square_root(q):
    q = QQ(q)
    if q < 0:
        return None
    n = ZZ(q.numerator())
    d = ZZ(q.denominator())
    if not n.is_square() or not d.is_square():
        return None
    return QQ(n.sqrt()) / QQ(d.sqrt())


def poly_sqrt_QQ(P):
    P = R(P)
    if not P:
        return R.zero()
    fac = P.factor()
    u = qq_square_root(QQ(fac.unit()))
    if u is None:
        return None
    out = R(u)
    for f, e in fac:
        e = int(e)
        if e % 2:
            return None
        out *= f ** (e // 2)
    assert out**2 == P
    return out


def rat_sqrt_QQ(value):
    value = K(value)
    if not value:
        return K.zero()
    n = R(value.numerator())
    d = R(value.denominator())
    nr = poly_sqrt_QQ(n)
    dr = poly_sqrt_QQ(d)
    if nr is None or dr is None or not dr:
        return None
    y = K(nr) / K(dr)
    assert y**2 == value
    return y


rhs = K(x**3 + A*x + B)
y = rat_sqrt_QQ(rhs)
if y is None:
    print(
        "Q32QQ_EXACT_SECTION|square=0|"
        "status=RECONSTRUCTED_ABX_BUT_Y_NOT_RATIONAL_SQUARE",
        flush=True,
    )
    print(
        "Q32QQ_RESULT|"
        f"primes={len(records)}|modulus_bits={ZZ(modulus).nbits()}|"
        f"holdout={int(holdout_pass)}|status=PARTIAL_ABX_ONLY",
        flush=True,
    )
    raise SystemExit(0)

assert y**2 == x**3 + A*x + B

# Pick the sign of y that agrees with as many saved modular choices as possible.
# The saved choices are arbitrary, so this is cosmetic; +/-y represent the two spinors.
plus = 0
minus = 0
for rec in records:
    try:
        ym = norm_mod_rf(y, rec[0])
        target = record_to_mod_rf(rec, "y_mod")
        plus += int(ym == target)
        minus += int(-ym == target)
    except ZeroDivisionError:
        pass
if minus > plus:
    y = -y
    plus, minus = minus, plus

ynb, ydb = bits_of_rf(y)
print(
    "Q32QQ_EXACT_SECTION|"
    f"square=1|max_num_bits={ynb}|max_den_bits={ydb}|"
    f"saved_sign_matches={plus}|saved_opposite_matches={minus}|"
    "status=PASS_EXACT_SECTION",
    flush=True,
)

# Fibre discriminant is an additional exact sanity check.
Delta = K(-16 * (4*A**3 + 27*B**2))
if not Delta:
    raise ArithmeticError("reconstructed Jacobian has zero discriminant")

def serial_rf(v):
    n = R(v.numerator())
    d = R(v.denominator())
    return {
        "num_degree": int(n.degree()),
        "den_degree": int(d.degree()),
        "num": [str(c) for c in n.list()],
        "den": [str(c) for c in d.list()],
    }

payload = {
    "schema": "elkies-k3.h3-q32-pointed-d12-qq-reconstruction.v1",
    "status": (
        "PASS_EXACT_Q32_POINTED_D12_QQ_WITH_HOLDOUT"
        if holdout_pass else
        "PASS_EXACT_Q32_POINTED_D12_QQ"
    ),
    "primes": [p for p,_,_ in records],
    "crt_modulus_bits": int(ZZ(modulus).nbits()),
    "independent_holdout": {
        "prime": holdout_prime,
        "passed": bool(holdout_pass),
    },
    "jacobian_A": serial_rf(A),
    "jacobian_B": serial_rf(B),
    "spinor_section_x": serial_rf(x),
    "spinor_section_y": serial_rf(y),
    "exact_weierstrass_identity": True,
    "interpretation": (
        "Characteristic-zero short D12 Jacobian candidate and transported spinor "
        "section reconstructed from modular q32 signatures. The section sign only "
        "records the arbitrary choice E10a versus E10b as zero."
    ),
}

out = LOCAL / "q32-pointed-d12-qq-reconstruction.json"
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{out}", flush=True)
print(
    "Q32QQ_RESULT|"
    f"primes={len(records)}|modulus_bits={ZZ(modulus).nbits()}|"
    f"holdout={int(holdout_pass)}|exact_identity=1|"
    f"status={payload['status']}",
    flush=True,
)
