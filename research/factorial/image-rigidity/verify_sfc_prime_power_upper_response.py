#!/usr/bin/env python3
"""Exact finite controls for the higher-response proof; not a parameter proof."""
if not __debug__:
    raise RuntimeError("Assertions are required; do not run Python with -O.")

from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research/factorial/image-rigidity").is_dir())
PROOFDIR = ROOT / "research/factorial/image-rigidity"
BASE_PATH = PROOFDIR / "verify_sfc_prime_power_scalar_transfer.py"
EXPECTED_SOURCE_SHA256 = {
    "SFC_PRIME_POWER_UPPER_RESPONSE.md": "c32d3c76f705aa5fbb5e046b23b1ca4c18fdc91568b5f3e5e36b3fd4774dbdf2",
    "SFC_PRIME_POWER_SCALAR_TRANSFER.md": "9fa457f230397c670c7b4152c0510df93da24df03fc08112f1c5bf5c73bc0349",
    "verify_sfc_prime_power_scalar_transfer.py": "cc64f7e35237b6b2fe674c5da8b7b54bdf3ceab17c6494e391ebeeb1c8bca18c",
    "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex": "9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
    "LOWER_TWO_FOURTH_POWER_GAP_SFC.tex": "e238569c21a0ccc7fb507d910608a734461fefaab4d96d5df0124afe73407631",
    "DIAGONAL_PRIME_GAP_SFC.tex": "db034b5fb29ece9c3a0cb8cbce4c82f692b933e3679b27ee823db5ca5e7f98eb"
}


def verify_sources():
    found = {}
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        actual = sha256((PROOFDIR/name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("Source SHA256 mismatch: "+name)
        found[name] = actual
    return found


spec = spec_from_file_location("scalar_transfer", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)
vp, ev, mul = base.vp, base.ev, base.mul


def valuation_json(x, p):
    return vp(x, p) if x else None


def response_series(p, s, n, gamma, h, order):
    """Solve actual nonresonant scalar leaves over exact rational series."""
    gamma = F(gamma)
    q = p**s
    assert n >= 1 and 1 <= h <= q and gamma > 0 and vp(gamma, p) >= 0
    c = gamma + q - 1
    b = (q - 1) // (p - 1)
    B = b - s
    states, nu = base.leaves(p, s, c, order)
    leaves = [z for z in states if z["e"]]
    nodes = [c - z["j"] for z in leaves]
    Q0 = [F(1)]
    for a in nodes:
        Q0 = mul(Q0, [-a, 1])
    Gs = []
    for i in range(q - 1):
        g = [F(1)]
        for j, a in enumerate(nodes):
            if j != i:
                g = mul(g, [-a, 1])
        Gs.append(g)
    U = [[[F(0)] * (q - 1) for _ in range(q - 1)]
         for _ in range(order + 1)]
    forcing = [[F(0)] * (q - 1) for _ in range(order + 1)]
    for i, state in enumerate(leaves):
        a = nodes[i]
        H = [x / state["h"][0] for x in state["h"]]
        assert H[0] == 1
        assert all(vp(x, p) >= m for m, x in enumerate(H))
        for ell in range(n):
            H = [x / (a + q * ell) for x in base.apply_operator(
                H, a + q * ell, state["F"], q, order)]
            assert H[0] == 1
            assert all(vp(x, p) >= m for m, x in enumerate(H))
        # Check the new slope filtration directly on the actual inner terms.
        TH = H[:]
        for k in range(1, q):
            TH = base.apply_operator(TH, F(0), state["F"], 1, order)
            tau = base.series_div(TH, H, order)
            assert tau[0] == 0
            assert all(vp(x, p) >= max(0, m-k)
                       for m, x in enumerate(tau))
        w = ev(Gs[i], a)
        for j, g in enumerate(Gs):
            row = base.series_div(base.poly_operator(
                g, H, a, state["F"], q, order), H, order)
            for m in range(order + 1):
                U[m][i][j] = row[m] * nu**m / w
            assert U[0][i][j] == int(i == j)
            assert all(vp(U[m][i][j], p) >= m
                       for m in range(1, order + 1))
        row = base.series_div(base.poly_operator(
            Q0, H, a, state["F"], q, order), H, order)
        for m in range(order + 1):
            forcing[m][i] = row[m] * nu**m / w
        assert forcing[0][i] == 0
        assert all(vp(forcing[m][i], p) >= s+m-1
                   for m in range(1, order + 1))
    bs = [[F(0)] * (q - 1) for _ in range(order + 1)]
    Qv = [Q0] + [[F(0)] * q for _ in range(order)]
    for m in range(1, order + 1):
        for i in range(q - 1):
            bs[m][i] = -forcing[m][i] - sum(
                U[r][i][j] * bs[m-r][j]
                for r in range(1, m+1) for j in range(q-1))
            assert vp(bs[m][i], p) >= s+m-1
            for k, coeff in enumerate(Gs[i]):
                Qv[m][k] += bs[m][i] * coeff
    Tv = [base.upper_functional(P, n, gamma, q, h) for P in Qv]
    Av = [ev(P, -q*n) for P in Qv]
    Ts = [base.upper_functional(g, n, gamma, q, h) for g in Gs]
    D = min(vp(x, p) for x in Ts)
    assert vp(Av[0], p) == B
    for m in range(1, order+1):
        assert vp(Av[m], p) >= B+m
        if D < 10**9:
            assert vp(Tv[m], p) >= s+m-1+D
        else:
            assert Tv[m] == 0
    C = gamma + q*n + q
    kappa = F(q*(n+1), 1) / base.rising(C, h)
    kval = vp(kappa, p)
    Fv = [Av[0]] + [Av[m] - (-1)**q * kappa * p**b * Tv[m-1]
                       for m in range(1, order+1)]
    if D < 10**9:
        e = kval + 2*s + D
        for m in range(2, order+1):
            assert vp(Fv[m], p) >= min(B+m, B+e+m-2)
            if e > 0:
                assert vp(Fv[m], p) > B
        R = max(1, 2-e)
        assert all(vp(Fv[m], p) > B
                   for m in range(R+1, order+1))
    else:
        e, R = None, 1
    return dict(p=p, s=s, q=q, n=n, gamma=gamma, h=h, order=order,
                b=b, B=B, nu=nu, Q0=Q0, Gs=Gs, Ts=Ts, D=D,
                Qv=Qv, Tv=Tv, Av=Av, Fv=Fv, kappa=kappa,
                kval=kval, e=e, R=R)


def formal_translation(r):
    """Independent polynomial substitution in F_p[[w]][X]."""
    p, order, Q0 = r["p"], r["order"], r["Q0"]
    G = [F(0)] * (order+1)
    power = 1
    while power <= order:
        G[power] = -1
        power *= p
    minus_G_powers = [[F(1)] + [F(0)]*order]
    for k in range(1, len(Q0)):
        minus_G_powers.append(base.series_mul(minus_G_powers[-1], G, order))
    expected = [[F(0)] * len(Q0) for _ in range(order+1)]
    for j, a in enumerate(Q0):
        for k in range(j+1):
            for m, coeff in enumerate(minus_G_powers[k]):
                expected[m][j-k] += a*comb(j, k)*coeff
    actual = [[a/(p*r["nu"])**m for a in row]
              for m, row in enumerate(r["Qv"])]
    assert all(vp(a, p) >= 0 for row in actual for a in row)
    if r["s"] == 1:
        assert all(base.residue(actual[m][j]-expected[m][j], p) == 0
                   for m in range(order+1) for j in range(len(Q0)))
        h = r["h"]
        if h >= 3:
            d = h-2
            assert order >= d+1
            assert vp(r["Tv"][d], p) == d
            expected_leading = (-1)**p * comb(p-1, d)
            assert base.residue(r["Tv"][d]/(p*r["nu"])**d, p) == (
                expected_leading % p)
            if r["kval"] < 0:
                threshold = r["kval"] + h-1
                assert vp(r["Fv"][h-1], p) == threshold
                assert all(vp(r["Fv"][m], p) > threshold
                           for m in range(h, order+1))
        return "prime translation identity checked coefficientwise"
    assert all(vp(a, p) >= r["s"]-1 for a in actual[1])
    # The actual w coefficient vanishes mod p, but -Q0' has unit leading term.
    assert base.residue(expected[1][-2], p) == (-(r["q"]-1)) % p != 0
    assert base.residue(actual[1][-2], p) == 0
    return "proper-power naive translation refuted at the w coefficient"


def poly_divrem(a, b):
    a, b = base.trim(a), base.trim(b)
    assert any(b)
    quotient = [F(0)] * max(1, len(a)-len(b)+1)
    while any(a) and len(a) >= len(b):
        d, x = len(a)-len(b), a[-1]/b[-1]
        quotient[d] += x
        for j, y in enumerate(b):
            a[d+j] -= x*y
        a = base.trim(a)
    return base.trim(quotient), a


def exact_moment_gcd(r):
    n, gamma, q, h = r["n"], r["gamma"], r["q"], r["h"]
    lower = [F(comb(n, k))*base.rising(gamma, q*k) for k in range(n+1)]
    upper = [F(comb(n+1, k))*base.rising(gamma+h, q*k) for k in range(n+2)]
    a, b = lower, upper
    while any(b):
        a, b = b, poly_divrem(a, b)[1]
    degree = len(base.trim(a))-1
    assert degree == 0
    return degree


def record(r, translation):
    p = r["p"]
    return dict(p=p, s=r["s"], q=r["q"], n=r["n"], gamma=str(r["gamma"]),
                h=r["h"], checked_outer_degree=r["order"],
                actual_complete_leaves=r["q"],
                kappa=str(r["kappa"]), kappa_valuation=r["kval"],
                D=None if r["D"] == 10**9 else r["D"],
                e=r["e"], denominator_depth_bound=r["R"],
                T0=str(r["Tv"][0]), Tj=[str(x) for x in r["Ts"]],
                upper_response_valuations=[valuation_json(x, p) for x in r["Tv"]],
                mixed_auxiliary_valuations=[valuation_json(x, p) for x in r["Fv"]],
                translation_control=translation,
                exact_moment_gcd_degree=exact_moment_gcd(r))


def main():
    sources = verify_sources()
    base.verify_sources()
    cases = [
        (2, 1, 1, F(1, 3), 2, 4),
        (3, 1, 2, F(1, 2), 3, 5),
        (5, 1, 1, F(11), 5, 7),
        (5, 1, 1, F(3111), 5, 7),
        (5, 1, 2, F(1, 2), 4, 6),
        (7, 1, 1, F(29), 7, 8),
        (2, 3, 2, F(5, 3), 7, 5),
        (3, 2, 1, F(1, 2), 8, 5),
    ]
    controls = []
    for args in cases:
        r = response_series(*args)
        t = formal_translation(r)
        if args[:5] == (5, 1, 1, F(11), 5):
            assert r["Tv"][0] == -32400
            assert r["Ts"] == [1146, 1184, 1224, 1266]
            assert sum(r["Ts"]) == 4820
            assert r["kappa"] == F(1, 637560)
            assert r["kval"] == -1 and r["e"] == 1
            assert vp(r["Tv"][1], 5) == 2
            assert vp(r["Tv"][3], 5) == 3
        if args[:5] == (5, 1, 1, F(3111), 5):
            assert r["kval"] == -4
        controls.append(record(r, t))
        print("PASS", args, flush=True)
    receipt = dict(
        status="PASS", evidence_type="exact finite algebra controls; not a quantified proof",
        imported_source_sha256=sources,
        imported_checker_sha256=sha256(BASE_PATH.read_bytes()).hexdigest(),
        own_note_sha256=sha256((HERE / "SFC_PRIME_POWER_UPPER_RESPONSE.md").read_bytes()).hexdigest(),
        own_checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        finite_cases=controls,
        independent_algebra="actual scalar-leaf reconstruction versus formal polynomial translation; direct rational Euclidean moment gcds",
        universal_claims="proved in companion note from the pinned complete scalar transfer; finite cases are controls only")
    output = HERE / "sfc_prime_power_upper_response_verification.json"
    output.write_text(json.dumps(receipt, indent=2)+"\n")
    print("PASS:", len(controls), "exact finite controls; receipt", output)


if __name__ == "__main__":
    main()
