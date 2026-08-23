#!/usr/bin/env python3
"""
Reconstruct exact mu = a + b*sqrt(-3) for Q80 orbit 1222.

Sources are auto-discovered in the permanent modular artifact directory:

1. Full modular Jacobian embedding files:
     q80_modj_p<P>_s+1_j+1.json, etc.
   These are converted to the monic-invariant mu residue.

2. Fast four-fiber files:
     q80_mu_fast_p<P>.json
   These already contain mu_a_mod_p and mu_j_mod_p.

If both exist for the same prime, they must agree.

The two rational components a,b are CRT/rational-reconstructed independently.
Every candidate is tested against all unused cached primes and against the
independent pinned p=73 original/conjugate kernels.

Writes on success to the permanent orbit data directory:
  q80_char0_orbit1222_mu_exact.json
  q80_char0_orbit1222_mu_exact.sage

Usage (from any working directory):
  python3 elkies-k3/scripts/q80-orbit1222-char0/mu_reconstruct_mixed.py
"""

import json
import math
import re
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
EMBS = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def inv(a, p):
    return pow(int(a) % p, -1, p)


def ef(p, ss, sj):
    return MODULAR_DATA / f"q80_modj_p{p}_s{ss:+d}_j{sj:+d}.json"


def load_json(path):
    return json.loads(path.read_text())


def pad(values, n):
    values = list(values)
    return values + [0] * (n - len(values))


def pmul(a, b, p, n=25):
    out = [0] * n
    for i, x in enumerate(a):
        x %= p
        if not x:
            continue
        for j, y in enumerate(b):
            if i + j >= n:
                break
            out[i + j] = (out[i + j] + x * (y % p)) % p
    return out


def ppow(a, e, p, n=25):
    out = [1] + [0] * (n - 1)
    base = pad(a, n)[:n]
    while e:
        if e & 1:
            out = pmul(out, base, p, n)
        e //= 2
        if e:
            base = pmul(base, base, p, n)
    return out


def proot(target, exponent, p, degree):
    target = [int(x) % p for x in pad(target, 25)[:25]]
    if target[0] != 1:
        raise ArithmeticError(f"root constant {target[0]} != 1")

    root = [0] * 25
    root[0] = 1
    ie = inv(exponent, p)

    for k in range(1, 25):
        known = ppow(root, exponent, p, 25)[k]
        root[k] = (target[k] - known) * ie % p

    if ppow(root, exponent, p, 25) != target:
        raise ArithmeticError("formal root check failed")
    if any(root[k] % p for k in range(degree + 1, 25)):
        raise ArithmeticError("formal root degree overflow")
    return root[: degree + 1]


def embedding_mu(kernel, p):
    kernel = [int(x) % p for x in kernel]
    if len(kernel) != 50:
        raise ArithmeticError(f"kernel length {len(kernel)} != 50")

    N = kernel[:25]
    D = kernel[25:]
    if N[0] != 1:
        raise ArithmeticError(f"N[0]={N[0]} != 1")

    C = proot(N, 3, p, 8)

    M = [(N[i] - (1728 % p) * D[i]) % p for i in range(25)]
    m0 = M[0]
    if m0 == 0:
        raise ArithmeticError("m0=0")

    S = proot([x * inv(m0, p) % p for x in M], 2, p, 12)

    c8 = C[8] % p
    s12 = S[12] % p
    if not c8 or not s12:
        raise ArithmeticError("degree drop in C/S")

    # mu for the monic normalization C/C8, S/S12.
    return m0 * s12 * s12 * inv(c8 * c8 * c8, p) % p


def full_packet_mu(p):
    rec = {e: load_json(ef(p, *e)) for e in EMBS}
    rs = int(rec[(1, 1)]["s_root"]) % p
    rj = int(rec[(1, 1)]["j_root"]) % p

    if rs * rs % p != (-6) % p or rj * rj % p != (-3) % p:
        raise ArithmeticError(f"p={p}: invalid roots")

    vals = {e: embedding_mu(rec[e]["kernel"], p) for e in EMBS}
    fpp = vals[(1, 1)]
    fpm = vals[(1, -1)]
    fmp = vals[(-1, 1)]
    fmm = vals[(-1, -1)]

    a = (fpp + fpm + fmp + fmm) * inv(4, p) % p
    s = (fpp + fpm - fmp - fmm) * inv(4 * rs, p) % p
    b = (fpp - fpm + fmp - fmm) * inv(4 * rj, p) % p
    sj = (fpp - fpm - fmp + fmm) * inv(4 * rs * rj, p) % p

    if s or sj:
        raise ArithmeticError(
            f"p={p}: unexpected sqrt(-6)/sqrt(18) support {s},{sj}"
        )

    return a, b


def discover_full():
    result = set()
    pattern = re.compile(r"q80_modj_p(\d+)_s\+1_j\+1\.json$")

    for path in MODULAR_DATA.glob("q80_modj_p*_s+1_j+1.json"):
        m = pattern.match(path.name)
        if not m:
            continue
        p = int(m.group(1))
        if p == 73:
            continue
        try:
            ok = all(
                ef(p, *e).exists()
                and len(load_json(ef(p, *e)).get("kernel", ())) == 50
                for e in EMBS
            )
        except Exception:
            ok = False
        if ok:
            result.add(p)

    return sorted(result)


def discover_fast():
    result = {}
    pattern = re.compile(r"q80_mu_fast_p(\d+)\.json$")

    for path in MODULAR_DATA.glob("q80_mu_fast_p*.json"):
        m = pattern.match(path.name)
        if not m:
            continue
        p = int(m.group(1))
        try:
            r = load_json(path)
            a = int(r["mu_a_mod_p"]) % p
            b = int(r["mu_j_mod_p"]) % p
        except Exception as exc:
            print(
                f"Q80MUMIX|p={p}|source=fast|status=SKIP|"
                f"type={type(exc).__name__}|message={exc}",
                flush=True,
            )
            continue
        result[p] = (a, b)

    return result


def parse_kernel(path):
    text = path.read_text()
    m = re.search(
        r"expected_kernel\s*=\s*vector\(\s*finite\s*,\s*\[(.*?)\]\s*,?\s*\)",
        text,
        re.S,
    )
    if not m:
        raise RuntimeError(f"cannot parse expected_kernel from {path}")

    values = [int(x) for x in re.findall(r"-?\d+", m.group(1))]
    if len(values) != 50:
        raise RuntimeError(f"{path}: parsed {len(values)} values")
    return values


def p73_mu():
    plus = parse_kernel(
        REPO_ROOT
        / "elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage"
    )
    minus = parse_kernel(
        REPO_ROOT
        / "elkies-k3/scripts/analyze_q80_third_q12_galois_descent_gf73.sage"
    )

    up = embedding_mu(plus, 73)
    um = embedding_mu(minus, 73)

    a = (up + um) * inv(2, 73) % 73
    b = (up - um) * inv(2 * 17, 73) % 73
    return a, b


def crt_pair(x, m, y, p):
    t = ((y - x) % p) * inv(m, p) % p
    return (x + m * t) % (m * p), m * p


def crt(residues):
    x, m = 0, 1
    for p, y in residues:
        x, m = crt_pair(x, m, y, p)
    return x, m


def ratrec(x, m):
    x %= m
    if x == 0:
        return Fraction(0, 1)

    bound = math.isqrt(m // 2)
    r0, r1 = m, x
    t0, t1 = 0, 1

    while abs(r1) > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1

    n, d = r1, t1
    if d == 0:
        return None
    if d < 0:
        n, d = -n, -d

    g = math.gcd(abs(n), d)
    n //= g
    d //= g

    if d <= 0 or abs(n) > bound or d > bound:
        return None
    if (n - x * d) % m:
        return None

    return Fraction(n, d)


def qmod(q, p):
    if q.denominator % p == 0:
        raise ZeroDivisionError
    return q.numerator % p * inv(q.denominator, p) % p


def main():
    residues = {}
    sources = {}

    full_primes = discover_full()
    fast = discover_fast()

    for p in full_primes:
        try:
            ab = full_packet_mu(p)
        except Exception as exc:
            print(
                f"Q80MUMIX|p={p}|source=full|status=SKIP|"
                f"type={type(exc).__name__}|message={exc}",
                flush=True,
            )
            continue
        residues[p] = ab
        sources[p] = ["full"]

    for p, ab in fast.items():
        if p in residues and residues[p] != ab:
            raise SystemExit(
                f"p={p}: FAST/FULL DISAGREEMENT: "
                f"full={residues[p]} fast={ab}"
            )
        residues[p] = ab
        sources.setdefault(p, []).append("fast")

    primes = sorted(residues)
    if len(primes) < 3:
        raise SystemExit(f"need >=3 residues; found {primes}")

    a73, b73 = p73_mu()

    print(
        f"Q80MUMIX|full_primes={len(full_primes)}|"
        f"fast_primes={len(fast)}|unique_primes={len(primes)}|"
        f"p73={a73},{b73}",
        flush=True,
    )
    print(
        "Q80MUMIX|primes="
        + ",".join(
            f"{p}[{'+'.join(sources[p])}]" for p in primes
        ),
        flush=True,
    )

    accepted = None

    for k in range(3, len(primes) + 1):
        train = primes[:k]
        holdout = primes[k:]

        xa, M = crt([(p, residues[p][0]) for p in train])
        xb, M2 = crt([(p, residues[p][1]) for p in train])
        assert M == M2

        a = ratrec(xa, M)
        b = ratrec(xb, M)

        if a is None or b is None:
            print(
                f"Q80MUMIX|k={k}|modulus_digits={len(str(M))}|"
                f"a={'?' if a is None else a}|"
                f"b={'?' if b is None else b}|status=UNRESOLVED",
                flush=True,
            )
            continue

        try:
            train_ok = all(
                qmod(a, p) == residues[p][0]
                and qmod(b, p) == residues[p][1]
                for p in train
            )
            holdout_ok = all(
                qmod(a, p) == residues[p][0]
                and qmod(b, p) == residues[p][1]
                for p in holdout
            )
            p73_ok = qmod(a, 73) == a73 and qmod(b, 73) == b73
        except ZeroDivisionError:
            train_ok = holdout_ok = p73_ok = False

        print(
            f"Q80MUMIX|k={k}|modulus_digits={len(str(M))}|"
            f"a_num_digits={len(str(abs(a.numerator)))}|"
            f"a_den_digits={len(str(a.denominator))}|"
            f"b_num_digits={len(str(abs(b.numerator)))}|"
            f"b_den_digits={len(str(b.denominator))}|"
            f"train={train_ok}|holdout={holdout_ok}|p73={p73_ok}",
            flush=True,
        )

        if train_ok and holdout_ok and p73_ok:
            accepted = (a, b, train, holdout, M)
            break

    if accepted is None:
        print("Q80MUMIX|status=NEED_MORE_RESIDUES", flush=True)
        raise SystemExit(2)

    a, b, train, holdout, M = accepted

    out_json = ORBIT_DATA / "q80_char0_orbit1222_mu_exact.json"
    out_sage = ORBIT_DATA / "q80_char0_orbit1222_mu_exact.sage"

    out_json.write_text(
        json.dumps(
            {
                "version": 1,
                "field": "QQ(sqrt(-3))",
                "mu_a": str(a),
                "mu_j": str(b),
                "mu": f"({a}) + ({b})*sqrt(-3)",
                "reconstruction_primes": train,
                "holdout_primes": holdout,
                "crt_modulus": str(M),
                "p73_validation": True,
                "sources": {
                    str(p): sources[p] for p in primes
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    def qtext(q):
        return str(q.numerator) if q.denominator == 1 else (
            f"({q.numerator}/{q.denominator})"
        )

    out_sage.write_text(
        "#!/usr/bin/env sage\n"
        "from sage.all import QuadraticField\n"
        'K = QuadraticField(-3, "j")\n'
        "j = K.gen()\n"
        f"mu = {qtext(a)} + {qtext(b)}*j\n"
        'print(f"Q80MUCHAR0|mu={mu}|status=PASS_EXACT_MU")\n'
    )

    print(
        f"Q80MUMIX|a={a}|b={b}|"
        f"train={len(train)}|holdout={len(holdout)}|p73=PASS|"
        "status=PASS_EXACT_MU",
        flush=True,
    )
    print(f"Q80MUMIX|json={out_json}", flush=True)
    print(f"Q80MUMIX|sage={out_sage}", flush=True)


if __name__ == "__main__":
    main()
