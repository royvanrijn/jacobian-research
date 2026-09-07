#!/usr/bin/env sage -python
"""
Reconstruct the q32 D12 polynomial Weierstrass model over QQ.

Modular evidence gives, at every good prime,
    A = A0 / h^4
    B = B0 / h^6
    xP = X0 / h^2
with h monic quadratic and deg(A0,B0,X0)=(8,12,4).

This removes the redundant rational-function denominator coordinates that made
the earlier projective LLL problem unnecessarily large.

We reconstruct h, X0, A0, B0 as projective rational coefficient vectors using
all but the last prime and reserve the last prime as an independent holdout.
Then recover Y0 exactly from
    Y0^2 = X0^3 + A0*X0 + B0.
"""

import json, math
from itertools import combinations
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, matrix

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

pointed_path = LOCAL / "q32-pointed-spinor-weierstrass-anchor.json"
if not pointed_path.exists():
    raise SystemExit(f"missing {pointed_path}")
pointed = json.loads(pointed_path.read_text())
if pointed.get("status") != "PASS_POINTED_Q32_D12_SPINOR_MARKING":
    raise SystemExit("pointed anchor not passing")
pby = {int(r["prime"]): r for r in pointed["primes"]}

records = []
for spath in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig = json.loads(spath.read_text())
        p = int(sig["prime"])
    except Exception:
        continue
    if sig.get("status") != "PASS_Q32_MODP_SIGNATURE" or p not in pby:
        continue
    records.append((ZZ(p), sig, pby[p]))
records.sort(key=lambda r: int(r[0]))

if len(records) < 5:
    raise SystemExit("need at least five good signatures")

train = records[:-1]
hold = records[-1]
mods = [r[0] for r in train]
holdp = hold[0]

print(
    "Q32POLYQQ_INPUT|"
    f"train={','.join(map(str,mods))}|holdout={holdp}|"
    f"count={len(records)}|status=PASS",
    flush=True,
)


def mod_objects(row):
    p, sig, pt = row
    F = GF(p)
    R = PolynomialRing(F, "V")

    def rf(rec):
        n = R([F(v) for v in rec["num"]])
        d = R([F(v) for v in rec["den"]])
        lc = d.leading_coefficient()
        return n/lc, d/lc

    An, Ad = rf(sig["jacobian_A"])
    Bn, Bd = rf(sig["jacobian_B"])
    Xn, Xd = rf(pt["marked_section_x"])

    def exact_root(P, e):
        P = R(P)
        out = R.one()
        for f, m in P.factor():
            m = int(m)
            if m % e:
                return None
            out *= f.monic() ** (m // e)
        return out.monic()

    hA = exact_root(Ad, 4)
    hB = exact_root(Bd, 6)
    hX = exact_root(Xd, 2)
    if hA is None or hA != hB or hA != hX:
        raise ArithmeticError(f"prime {p}: common denominator structure failed")
    h = hA

    assert Ad == h**4 and Bd == h**6 and Xd == h**2
    assert An.degree() == 8 and Bn.degree() == 12 and Xn.degree() == 4
    assert h.degree() == 2

    return {
        "h": [ZZ(int(h[i])) for i in range(2)],  # leading coeff is fixed 1
        "X0": [ZZ(int(Xn[i])) for i in range(5)],
        "A0": [ZZ(int(An[i])) for i in range(9)],
        "B0": [ZZ(int(Bn[i])) for i in range(13)],
    }


moddata = [(r[0], mod_objects(r)) for r in records]
train_data = moddata[:-1]
hold_data = moddata[-1][1]


def crt_scalar(vals, primes):
    x = ZZ(0)
    M = ZZ(1)
    for a, p in zip(vals, primes):
        a = ZZ(a) % p
        t = ((a - x) % p) * ((M % p).inverse_mod(p)) % p
        x = (x + M*t) % (M*p)
        M *= p
    if x > M//2:
        x -= M
    return x, M


_, MOD = crt_scalar([0]*len(mods), mods)
print(
    f"Q32POLYQQ_MODULUS|bits={MOD.nbits()}|train_count={len(train)}|status=PASS",
    flush=True,
)


def primitive(vals):
    vals = [ZZ(x) for x in vals]
    g = ZZ(0)
    for x in vals:
        g = ZZ(math.gcd(int(g), abs(int(x))))
    if g > 1:
        vals = [x//g for x in vals]
    if vals[-1] < 0:
        vals = [-x for x in vals]
    return vals


def candidate_vectors(Rred):
    rows = [list(v) for v in Rred.rows()]
    seen = set()

    def add(v):
        v = primitive(v)
        k = tuple(v)
        if v[-1] and k not in seen:
            seen.add(k)
            return v
        return None

    for v in rows:
        q = add(v)
        if q is not None:
            yield q

    k = min(10, len(rows))
    for i, j in combinations(range(k), 2):
        for s in (1, -1):
            q = add([rows[i][c] + s*rows[j][c] for c in range(len(rows[i]))])
            if q is not None:
                yield q

    # Low dimension: also try combinations of three reduced rows.
    if len(rows) <= 10:
        k = min(7, len(rows))
        for i, j, l in combinations(range(k), 3):
            for s1 in (1, -1):
                for s2 in (1, -1):
                    q = add([
                        rows[i][c] + s1*rows[j][c] + s2*rows[l][c]
                        for c in range(len(rows[i]))
                    ])
                    if q is not None:
                        yield q


def reconstruct(name):
    arrs = [d[name] for _, d in train_data]
    hvals = hold_data[name]
    N = len(arrs[0])
    assert all(len(a) == N for a in arrs)

    residues = []
    for j in range(N):
        x, M = crt_scalar([a[j] for a in arrs], mods)
        assert M == MOD
        residues.append(x)

    # Lattice vectors [n_0,...,n_{N-1},d] satisfying n_j=d*r_j mod M.
    B = matrix(ZZ, N+1, N+1)
    for j in range(N):
        B[j,j] = MOD
    for j, r in enumerate(residues):
        B[N,j] = r
    B[N,N] = 1

    print(
        f"Q32POLYQQ_LLL_START|object={name}|dimension={N+1}|"
        f"modulus_bits={MOD.nbits()}|status=START",
        flush=True,
    )
    Rred = B.LLL(delta=0.99)
    print(f"Q32POLYQQ_LLL_REDUCED|object={name}|status=PASS", flush=True)

    scored = []
    for v in candidate_vectors(Rred):
        d = v[-1]
        if d % holdp == 0:
            matches = -1
        else:
            inv = (d % holdp).inverse_mod(holdp)
            matches = sum(
                int((v[j] % holdp) * inv % holdp) == int(hvals[j] % holdp)
                for j in range(N)
            )
        bits = max(abs(x).nbits() for x in v)
        norm2 = sum(x*x for x in v)
        scored.append((matches, bits, norm2, v))

    scored.sort(key=lambda z: (-z[0], z[1], z[2]))
    best = scored[0]
    matches, bits, _, vec = best

    print(
        f"Q32POLYQQ_LLL_BEST|object={name}|heldout={matches}/{N}|"
        f"height_bits={bits}|common_den_bits={abs(vec[-1]).nbits()}|"
        f"status={'PASS_HELDOUT' if matches == N else 'PARTIAL'}",
        flush=True,
    )
    for rank, (mm, bb, nn, vv) in enumerate(scored[:5]):
        print(
            f"Q32POLYQQ_LLL_SHORT|object={name}|rank={rank}|"
            f"heldout={mm}/{N}|height_bits={bb}|"
            f"common_den_bits={abs(vv[-1]).nbits()}",
            flush=True,
        )

    if matches != N:
        return None

    d = QQ(vec[-1])
    coeffs = [QQ(x)/d for x in vec[:-1]]
    return coeffs, vec


# Smallest/most diagnostic first.
exact = {}
vectors = {}
for name in ("h", "X0", "A0", "B0"):
    ans = reconstruct(name)
    if ans is None:
        print(
            f"Q32POLYQQ_RESULT|failed={name}|modulus_bits={MOD.nbits()}|"
            "status=NEED_MORE_PRIMES",
            flush=True,
        )
        raise SystemExit(0)
    exact[name], vectors[name] = ans


RQ = PolynomialRing(QQ, "V")
V = RQ.gen()
h = RQ(exact["h"] + [QQ(1)])
X0 = RQ(exact["X0"])
A0 = RQ(exact["A0"])
B0 = RQ(exact["B0"])

assert h.degree() == 2
assert X0.degree() == 4
assert A0.degree() == 8
assert B0.degree() == 12

print(
    "Q32POLYQQ_MODEL|"
    f"hdeg={h.degree()}|X0deg={X0.degree()}|A0deg={A0.degree()}|B0deg={B0.degree()}|"
    "status=PASS_RECONSTRUCTED",
    flush=True,
)


def reduce_q(q, p, F):
    q = QQ(q)
    den = ZZ(q.denominator()) % p
    if not den:
        raise ZeroDivisionError
    return F(ZZ(q.numerator()) % p) / F(den)


def reduce_poly(P, p):
    F = GF(p)
    R = PolynomialRing(F, "V")
    return R([reduce_q(c, p, F) for c in RQ(P).list()])


# Audit against every prime including the holdout.
for p, md in moddata:
    F = GF(p)
    R = PolynomialRing(F, "V")
    checks = (
        reduce_poly(h, p) == R(md["h"] + [1]),
        reduce_poly(X0, p) == R(md["X0"]),
        reduce_poly(A0, p) == R(md["A0"]),
        reduce_poly(B0, p) == R(md["B0"]),
    )
    print(
        f"Q32POLYQQ_REDUCTION|prime={p}|h={int(checks[0])}|X0={int(checks[1])}|"
        f"A0={int(checks[2])}|B0={int(checks[3])}|"
        f"status={'PASS' if all(checks) else 'FAIL'}",
        flush=True,
    )
    if not all(checks):
        raise ArithmeticError(f"reduction mismatch at {p}")


def qq_sqrt(q):
    q = QQ(q)
    if q < 0:
        return None
    n = ZZ(q.numerator())
    d = ZZ(q.denominator())
    if not n.is_square() or not d.is_square():
        return None
    return QQ(n.sqrt())/QQ(d.sqrt())


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


rhs = X0**3 + A0*X0 + B0
Y0 = poly_sqrt(rhs)
if Y0 is None:
    print(
        "Q32POLYQQ_SECTION|square=0|"
        "status=MODEL_RECONSTRUCTED_SECTION_NOT_SQUARE",
        flush=True,
    )
    raise SystemExit(0)

assert Y0**2 == rhs
print(
    f"Q32POLYQQ_SECTION|square=1|Y0deg={Y0.degree()}|"
    "status=PASS_EXACT_POLYNOMIAL_SECTION",
    flush=True,
)

Delta = -16*(4*A0**3 + 27*B0**2)
assert Delta
fac = Delta.factor()
print(
    "Q32POLYQQ_DISCRIMINANT|"
    f"degree={Delta.degree()}|"
    f"factors={';'.join(str(f.degree())+'^'+str(int(e)) for f,e in fac)}|"
    "status=PASS",
    flush=True,
)


def serial_poly(P):
    P = RQ(P)
    return {
        "degree": int(P.degree()),
        "coefficients_low_to_high": [str(c) for c in P.list()],
    }


payload = {
    "schema": "elkies-k3.h3-q32-pointed-d12-polynomial-qq.v1",
    "status": "PASS_EXACT_Q32_POINTED_D12_POLYNOMIAL_QQ_HELDOUT",
    "training_primes": [int(p) for p in mods],
    "heldout_prime": int(holdp),
    "crt_modulus_bits": int(MOD.nbits()),
    "h": serial_poly(h),
    "A0": serial_poly(A0),
    "B0": serial_poly(B0),
    "spinor_X0": serial_poly(X0),
    "spinor_Y0": serial_poly(Y0),
    "exact_weierstrass_identity": True,
    "discriminant": serial_poly(Delta),
    "projective_integer_vectors": {
        k: [str(x) for x in vectors[k]]
        for k in ("h","X0","A0","B0")
    },
    "interpretation": (
        "The rational q32 Jacobian signatures are one polynomial K3 Weierstrass "
        "model y^2=x^3+A0*x+B0 viewed through x=X0/h^2, y=Y0/h^3 and "
        "A=A0/h^4, B=B0/h^6. The displayed spinor section is exact over QQ(V)."
    ),
}

out = LOCAL / "q32-pointed-d12-polynomial-qq.json"
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{out}", flush=True)
print(
    "Q32POLYQQ_RESULT|"
    f"modulus_bits={MOD.nbits()}|holdout={holdp}|exact_identity=1|"
    f"status={payload['status']}",
    flush=True,
)
