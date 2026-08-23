#!/usr/bin/env sage -python
"""
Projective/LLL reconstruction of the q32 D12 model and spinor section.

Scalar rational reconstruction wastes roughly half the CRT bits on every
coefficient.  But each modular rational function

    N(V) / D(V),  D monic

is the reduction of ONE projective QQ coefficient vector.  Clear all QQ
denominators simultaneously and reconstruct the primitive integer vector

    [coeffs(N), coeffs(D)]

from the congruences a_j == d*r_j (mod M), where d is the final (leading
denominator) projective coordinate.

Use all but the final prime for LLL and reserve the final prime as an
independent holdout.  If A,B,x reconstruct, recover y exactly from the
Weierstrass equation, avoiding modular square-root sign ambiguity.
"""

import json, math
from pathlib import Path
from itertools import combinations
from sage.all import GF, PolynomialRing, QQ, ZZ, matrix

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

POINTED = LOCAL / "q32-pointed-spinor-weierstrass-anchor.json"
if not POINTED.exists():
    raise SystemExit(f"missing {POINTED}")
pdata = json.loads(POINTED.read_text())
if pdata.get("status") != "PASS_POINTED_Q32_D12_SPINOR_MARKING":
    raise SystemExit("pointed anchor not passing")
pby = {int(r["prime"]): r for r in pdata["primes"]}

records = []
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig = json.loads(path.read_text())
        p = int(sig["prime"])
    except Exception:
        continue
    if sig.get("status") != "PASS_Q32_MODP_SIGNATURE" or p not in pby:
        continue
    records.append((ZZ(p), sig, pby[p]))
records.sort(key=lambda t: int(t[0]))

if len(records) < 4:
    raise SystemExit("need at least four compatible primes")

train = records[:-1]
hold = records[-1]
mods = [r[0] for r in train]
holdp = hold[0]

print(
    "Q32PROJLLL_INPUT|"
    f"train={','.join(map(str,mods))}|holdout={holdp}|"
    f"count={len(records)}|status=PASS",
    flush=True,
)


def rec_for(row, key):
    if key == "A":
        return row[1]["jacobian_A"]
    if key == "B":
        return row[1]["jacobian_B"]
    if key == "x":
        return row[2]["marked_section_x"]
    if key == "y":
        return row[2]["marked_section_y"]
    raise KeyError(key)


def flat(rec):
    return [ZZ(v) for v in rec["num"]] + [ZZ(v) for v in rec["den"]]


def crt_scalar(vals, mods):
    x = ZZ(0)
    M = ZZ(1)
    for rr, p in zip(vals, mods):
        rr = ZZ(rr) % p
        t = ((rr - x) % p) * ((M % p).inverse_mod(p)) % p
        x = (x + M*t) % (M*p)
        M *= p
    if x > M//2:
        x -= M
    return x, M


_, MOD = crt_scalar([0]*len(mods), mods)
print(
    f"Q32PROJLLL_MODULUS|bits={MOD.nbits()}|train_count={len(train)}|status=PASS",
    flush=True,
)


def primitive(v):
    vals = [ZZ(a) for a in v]
    g = ZZ(0)
    for a in vals:
        g = ZZ(math.gcd(int(g), abs(int(a))))
    if g > 1:
        vals = [a//g for a in vals]
    if vals[-1] < 0:
        vals = [-a for a in vals]
    return vals


def holdout_matches(vec, hvals):
    d = ZZ(vec[-1])
    if d % holdp == 0:
        return -1
    inv = (d % holdp).inverse_mod(holdp)
    return sum(
        int((ZZ(a) % holdp) * inv % holdp) == int(ZZ(r) % holdp)
        for a, r in zip(vec[:-1], hvals[:-1])
    )


def candidate_vectors(R):
    rows = [list(v) for v in R.rows()]
    seen = set()

    def emit(v):
        pv = primitive(v)
        key = tuple(pv)
        if key not in seen and pv[-1] != 0:
            seen.add(key)
            return pv
        return None

    for v in rows:
        q = emit(v)
        if q is not None:
            yield q

    # Sometimes the sought projective vector is a tiny combination of the
    # first reduced rows rather than literally a reduced basis row.
    k = min(12, len(rows))
    for i, j in combinations(range(k), 2):
        for s in (1, -1):
            v = [rows[i][c] + s*rows[j][c] for c in range(len(rows[i]))]
            q = emit(v)
            if q is not None:
                yield q


def reconstruct_projective(key):
    tref = rec_for(train[0], key)
    dims = (int(tref["num_degree"]), int(tref["den_degree"]))
    all_train = [flat(rec_for(r, key)) for r in train]
    hvals = flat(rec_for(hold, key))
    L = len(all_train[0])

    assert all(len(v) == L for v in all_train)
    assert len(hvals) == L
    # norm_rf() makes the leading denominator coefficient 1 at every prime.
    assert all((vals[-1] % p) == 1 for vals, (p,_,_) in zip(all_train, train))
    assert hvals[-1] % holdp == 1

    N = L - 1
    residues = []
    for j in range(N):
        x, M = crt_scalar([vals[j] for vals in all_train], mods)
        assert M == MOD
        residues.append(x)

    B = matrix(ZZ, N+1, N+1)
    for j in range(N):
        B[j,j] = MOD
    for j, r in enumerate(residues):
        B[N,j] = r
    B[N,N] = 1

    print(
        f"Q32PROJLLL_START|field={key}|dimension={N+1}|"
        f"modulus_bits={MOD.nbits()}|status=START",
        flush=True,
    )
    Rred = B.LLL(delta=0.99)
    print(f"Q32PROJLLL_REDUCED|field={key}|status=PASS", flush=True)

    scored = []
    for v in candidate_vectors(Rred):
        m = holdout_matches(v, hvals)
        bits = max(abs(a).nbits() for a in v)
        norm2 = sum(a*a for a in v)
        scored.append((m, bits, norm2, v))

    scored.sort(key=lambda z: (-z[0], z[1], z[2]))
    best = scored[0]
    m, bits, norm2, vec = best
    print(
        f"Q32PROJLLL_BEST|field={key}|heldout={m}/{N}|"
        f"height_bits={bits}|projective_den_bits={abs(vec[-1]).nbits()}|"
        f"status={'PASS_HELDOUT' if m == N else 'PARTIAL'}",
        flush=True,
    )
    for rank, (mm, bb, nn, vv) in enumerate(scored[:5]):
        print(
            f"Q32PROJLLL_SHORT|field={key}|rank={rank}|heldout={mm}/{N}|"
            f"height_bits={bb}|den_bits={abs(vv[-1]).nbits()}",
            flush=True,
        )

    if m != N:
        return None

    # Convert primitive integer projective vector into exact rational function.
    dproj = QQ(vec[-1])
    qvals = [QQ(a)/dproj for a in vec]
    nlen = len(tref["num"])
    num = qvals[:nlen]
    den = qvals[nlen:]
    RQ = PolynomialRing(QQ, "V")
    KQ = RQ.fraction_field()
    n = RQ(num)
    d = RQ(den)
    assert d and d.leading_coefficient() == 1
    value = KQ(n)/KQ(d)
    return value, vec


exact = {}
rawvec = {}
for key in ("A", "B", "x"):
    ans = reconstruct_projective(key)
    if ans is None:
        print(
            f"Q32PROJLLL_RESULT|failed={key}|modulus_bits={MOD.nbits()}|"
            "status=NEED_MORE_PRIMES_OR_STRONGER_RECONSTRUCTION",
            flush=True,
        )
        raise SystemExit(0)
    exact[key], rawvec[key] = ans

A, B, x = exact["A"], exact["B"], exact["x"]
KQ = A.parent()
RQ = KQ.ring()


def reduce_q(q, p, F):
    q = QQ(q)
    den = ZZ(q.denominator()) % p
    if not den:
        raise ZeroDivisionError
    return F(ZZ(q.numerator()) % p) / F(den)


def reduce_rf(v, p):
    F = GF(p)
    RF = PolynomialRing(F, "V")
    KF = RF.fraction_field()
    n = RQ(v.numerator())
    d = RQ(v.denominator())
    nn = RF([reduce_q(c, p, F) for c in n.list()])
    dd = RF([reduce_q(c, p, F) for c in d.list()])
    lc = dd.leading_coefficient()
    return KF(nn/lc) / KF(dd/lc)


def mod_record(row, key):
    p = row[0]
    F = GF(p)
    RF = PolynomialRing(F, "V")
    KF = RF.fraction_field()
    rr = rec_for(row, key)
    return KF(RF([F(v) for v in rr["num"]])) / KF(RF([F(v) for v in rr["den"]]))


# Full reduction audit over all 11 primes.
for row in records:
    p = row[0]
    checks = {
        key: reduce_rf(exact[key], p) == mod_record(row, key)
        for key in ("A", "B", "x")
    }
    print(
        "Q32PROJLLL_REDUCTION|"
        f"prime={p}|A={int(checks['A'])}|B={int(checks['B'])}|x={int(checks['x'])}|"
        f"status={'PASS' if all(checks.values()) else 'FAIL'}",
        flush=True,
    )
    if not all(checks.values()):
        raise ArithmeticError(f"reduction mismatch at {p}")


def qq_sqrt(q):
    q = QQ(q)
    if q < 0:
        return None
    n, d = ZZ(q.numerator()), ZZ(q.denominator())
    if not n.is_square() or not d.is_square():
        return None
    return QQ(n.sqrt()) / QQ(d.sqrt())


def poly_sqrt(P):
    P = RQ(P)
    if not P:
        return RQ.zero()
    fac = P.factor()
    u = qq_sqrt(QQ(fac.unit()))
    if u is None:
        return None
    out = RQ(u)
    for f, e in fac:
        e = int(e)
        if e % 2:
            return None
        out *= f**(e//2)
    assert out**2 == P
    return out


def rat_sqrt(v):
    v = KQ(v)
    n = poly_sqrt(v.numerator())
    d = poly_sqrt(v.denominator())
    if n is None or d is None or not d:
        return None
    y = KQ(n)/KQ(d)
    assert y**2 == v
    return y


rhs = KQ(x**3 + A*x + B)
y = rat_sqrt(rhs)
if y is None:
    print(
        "Q32PROJLLL_SECTION|rhs_square=0|"
        "status=ABX_RECONSTRUCTED_BUT_SECTION_SQRT_FAILED",
        flush=True,
    )
    raise SystemExit(0)

assert y**2 == x**3 + A*x + B

# Choose cosmetic sign agreeing with the majority of stored modular choices.
plus = minus = 0
for row in records:
    p = row[0]
    ym = reduce_rf(y, p)
    target = mod_record(row, "y")
    plus += int(ym == target)
    minus += int(-ym == target)
if minus > plus:
    y = -y
    plus, minus = minus, plus

print(
    f"Q32PROJLLL_SECTION|rhs_square=1|saved_sign_matches={plus}|"
    f"opposite={minus}|status=PASS_EXACT_SECTION",
    flush=True,
)


def serial(v):
    n = RQ(v.numerator())
    d = RQ(v.denominator())
    return {
        "num_degree": int(n.degree()),
        "den_degree": int(d.degree()),
        "num": [str(c) for c in n.list()],
        "den": [str(c) for c in d.list()],
    }


payload = {
    "schema": "elkies-k3.h3-q32-pointed-d12-projective-lll-qq.v1",
    "status": "PASS_EXACT_Q32_POINTED_D12_QQ_PROJECTIVE_LLL_HELDOUT",
    "training_primes": [int(r[0]) for r in train],
    "heldout_prime": int(holdp),
    "crt_modulus_bits": int(MOD.nbits()),
    "jacobian_A": serial(A),
    "jacobian_B": serial(B),
    "spinor_section_x": serial(x),
    "spinor_section_y": serial(y),
    "exact_weierstrass_identity": True,
    "projective_integer_vectors": {
        key: [str(a) for a in rawvec[key]]
        for key in ("A", "B", "x")
    },
}
out = LOCAL / "q32-pointed-d12-qq-projective-lll.json"
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{out}", flush=True)
print(
    "Q32PROJLLL_RESULT|"
    f"modulus_bits={MOD.nbits()}|holdout={holdp}|exact_identity=1|"
    f"status={payload['status']}",
    flush=True,
)
