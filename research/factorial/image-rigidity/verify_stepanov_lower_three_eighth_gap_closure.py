#!/usr/bin/env python3
"""Exact finite controls for the all-index gap-eight proof; not an all-index scan."""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def scale(a, k):
    return trim([k * x for x in a])


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x * y
    return trim(out)


def ev(a, x):
    s = 0
    for v in reversed(a):
        s = s*x + v
    return s


def vp(x, p):
    if x == 0:
        return None
    x = Fraction(x)
    v = 0
    a, b = abs(x.numerator), x.denominator
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def rank_mod(cols, p):
    a = [[x % p for x in row] for row in zip(*cols)]
    r = 0
    for j in range(len(cols)):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][j], -1, p)
        a[r] = [x*inv % p for x in a[r]]
        for i in range(len(a)):
            if i != r:
                k = a[i][j]
                a[i] = [(x-k*y) % p for x, y in zip(a[i], a[r])]
        r += 1
    return r


def reduced_A(v, d, p, c=0):
    out = [0] * d
    for j, x in enumerate(v):
        out[j] += (c-j)*x
        if j+1 < d:
            out[j+1] += x
    return [x % p for x in out]


def actual_Bs(n, d):
    gamma = 3*n+1
    out = [[1]]
    for j in range(n+d):
        f = out[-1]
        g = [0]*(len(f)+1)
        for k, a in enumerate(f):
            g[k] += (d*j+gamma+d-1-k)*a
            g[k+1] += a
        out.append(trim(g))
    return out


def operator_apply(f, c):
    out = [0]*(len(f)+1)
    for j, x in enumerate(f):
        out[j] += (c-j)*x
        out[j+1] += x
    return trim(out)


def polynomial_of_operator(poly, c, seed):
    out = [0]
    v = list(seed)
    for x in poly:
        out = add(out, scale(v, x))
        v = operator_apply(v, c)
    return out


def check_mixed(n, d):
    c = 3*n+d
    M = (d+3)*n+1
    C = M+d
    kappa = Fraction(d, (d+3)*((d+3)*n+d+1)*((d+3)*n+d+2))
    m = Fraction((d-2)*(2*M+d*d+1), 2) - (M+d+2)
    Rs = [[1]]
    for i in range(1, d):
        Rs.append(mul(Rs[-1], [d*(n+i-1), 1]))
    # Each column is a pair of polynomials in the operator X, the
    # constant and linear coefficients of the outer scalar u.
    columns = []
    for j in range(d-1):
        L = Rs[1]
        for r in range(j):
            L = mul(L, [c-r, -1])
        k = d-1-j
        vk = 1 if k == 2 else C+2 if k == 1 else 0
        corr = [-((-1)**(d+1))*kappa*vk]
        columns.append([L, corr])
    for j in range(d-1):
        columns[j] = [scale(x, (-1)**j) for x in columns[j]]
        # Eliminate in descending polynomial degree; previous columns
        # have the exact R_i polynomial as their u^0 component.
        for i in reversed(range(j)):
            rem = add(columns[j][0], scale(Rs[j+1], -1))
            q = rem[i+1] if i+1 < len(rem) else 0
            columns[j] = [add(x, scale(y, -q))
                          for x, y in zip(columns[j], columns[i])]
        assert columns[j][0] == Rs[j+1]
        expected = -kappa if j == d-3 else -kappa*m if j == d-2 else 0
        assert columns[j][1] == [expected]
    s = (d-9)*n+Fraction((d-4)*(d+3), 2)
    assert m-d*sum(n+r for r in range(d-1)) == s
    # The polynomial operator representation retains the actual seed1.
    Bs = actual_Bs(n, d)
    Pi = [1]
    for ell in range(n):
        Pi = mul(Pi, [d*ell, 1])
    for i in range(1, d):
        assert polynomial_of_operator(mul(Pi, Rs[i]), c, [1]) == Bs[n+i]


def sharpness():
    prod = [1]
    for j in range(1, 9):
        prod = mul(prod, [-j, 1])
    prod[0] -= factorial(8)
    assert prod[0] == 0
    Q = prod[1:]
    assert vp(Q[0], 2) == 4
    assert vp(Q[1], 2) == 2
    assert vp(Q[2], 2) >= 2
    assert vp(Q[6], 2) >= 1
    assert vp(Q[3], 2) == vp(Q[7], 2) == 0
    value = operator_apply(polynomial_of_operator(Q, 8, [1]), 8)
    assert value == [-factorial(8)] + [0]*7 + [1]
    return {"n": 1, "c": 8, "u": factorial(8), "Q": Q,
            "coefficient_valuations": [vp(x, 2) for x in Q]}


def series_mul(a, b, degree):
    return (mul(a, b)+[0]*(degree+1))[:degree+1]


def series_div(a, b, degree):
    out = [Fraction(0)]*(degree+1)
    for k in range(degree+1):
        out[k] = (Fraction(a[k] if k < len(a) else 0)
                  - sum(b[j]*out[k-j] for j in range(1, min(k, len(b)-1)+1)))/b[0]
    return out


def ah_inverse(degree):
    e = [Fraction(1)]
    for m in range(1, degree+1):
        powers = []
        k = 1
        while k <= m:
            powers.append(k)
            k *= 2
        e.append(-sum(e[m-k] for k in powers)/m)
    assert all(x == 0 or vp(x, 2) >= 0 for x in e)
    return e


def leaves(c, degree=2, initial_seed=None):
    initial_degree = 8*degree+7
    seed = [Fraction(1)]+[Fraction(0)]*initial_degree if initial_seed is None else (
        list(map(Fraction, initial_seed))+[Fraction(0)]*(initial_degree+1))[:initial_degree+1]
    states = [{"c": c, "F": {1: Fraction(1)}, "h": seed,
               "j": 0, "e": 0, "z_scale": Fraction(1)}]
    scale_controls = []
    for level in range(3):
        next_states = []
        for state in states:
            parent_degree = len(state["h"])-1
            child_degree = (parent_degree-1)//2
            b = Fraction(state["F"].get(2, 0)-1)
            assert vp(b, 2) == 0
            gauged = series_mul(ah_inverse(parent_degree), state["h"], parent_degree)
            newF = {1: Fraction(1)}
            power = 2
            while power <= child_degree:
                newF[power] = (state["F"].get(2*power, 0)-1)*2**(power-1)/b**power
                assert newF[power] == 0 or vp(newF[power], 2) >= 1
                power *= 2
            scale_controls.append([level, str(b)])
            for parity in [0, 1]:
                e = (state["c"]-parity) % 2
                h = [gauged[2*k+parity]*(2/b)**k for k in range(child_degree+1)]
                next_states.append({
                    "c": (state["c"]-parity-e)//2, "F": newF, "h": h,
                    "j": state["j"]+2**level*parity,
                    "e": state["e"]+2**level*e,
                    "z_scale": b*state["z_scale"]**2/2,
                })
        states = next_states
    for state in states:
        assert state["z_scale"] == Fraction(-27, 128)
        assert 8*state["c"]+state["e"] == c-state["j"]
        assert state["e"] == (c-state["j"]) % 8
        if initial_seed is None:
            assert vp(state["h"][0], 2) == 0
            assert all(x == 0 or vp(x, 2) >= 1 for x in state["h"][1:])
    return sorted(states, key=lambda x: x["j"])


def leaf_apply(h, a, F, degree):
    out = [Fraction((a-8*k)*h[k]) for k in range(degree+1)]
    for power, coeff in F.items():
        for k in range(power, degree+1):
            out[k] += 8*coeff*h[k-power]
    return out


def leaf_poly(poly, h, a, F, degree):
    v = h[:]
    out = [Fraction(0)]*(degree+1)
    for coef in poly:
        out = [x+coef*y for x, y in zip(out, v)]
        v = leaf_apply(v, a, F, degree)
    return out


def leaf_controls(n, degree=2):
    c = 3*n+8
    ls = [s for s in leaves(c, degree) if s["e"]]
    nodes = [c-s["j"] for s in ls]
    Q0 = [1]
    for a in nodes:
        Q0 = mul(Q0, [-a, 1])
    Gs = []
    for i in range(7):
        g = [1]
        for j, a in enumerate(nodes):
            if i != j:
                g = mul(g, [-a, 1])
        Gs.append(g)
    U = [[[Fraction(0)]*7 for _ in range(7)] for _ in range(degree+1)]
    forcing = [[Fraction(0)]*7 for _ in range(degree+1)]
    for i, s in enumerate(ls):
        a = nodes[i]
        H = [x/s["h"][0] for x in s["h"]]
        # Normalize by the constant factor after each multiplication.
        # It is a unit after removing its exact 2^v(e) content.
        for ell in range(n):
            H = [x/(a+8*ell) for x in leaf_apply(H, a+8*ell, s["F"], degree)]
        assert H[0] == 1
        assert all(x == 0 or vp(x, 2) >= 1 for x in H[1:])
        wi = ev(Gs[i], a)
        for j, g in enumerate(Gs):
            row = series_div(leaf_poly(g, H, a, s["F"], degree), H, degree)
            for k in range(degree+1):
                U[k][i][j] = row[k]*(-27)**k/wi
            assert U[0][i][j] == int(i == j)
            assert all(U[k][i][j] == 0 or vp(U[k][i][j], 2) >= 1
                       for k in range(1, degree+1))
        row = series_div(leaf_poly(Q0, H, a, s["F"], degree), H, degree)
        for k in range(degree+1):
            forcing[k][i] = row[k]*(-27)**k/wi
        assert forcing[0][i] == 0
        assert all(forcing[k][i] == 0 or vp(forcing[k][i], 2) >= 3
                   for k in range(1, degree+1))
    bs = [[Fraction(0)]*7 for _ in range(degree+1)]
    for k in range(1, degree+1):
        for i in range(7):
            bs[k][i] = -forcing[k][i]-sum(
                U[r][i][j]*bs[k-r][j] for r in range(1, k+1) for j in range(7))
            assert bs[k][i] == 0 or vp(bs[k][i], 2) >= 3
    Qs = [Q0]
    for k in range(1, degree+1):
        poly = [0]
        for i in range(7):
            poly = add(poly, scale(Gs[i], bs[k][i]))
        Qs.append(poly+[0]*(8-len(poly)))
    S0 = Q0[6]+22-n
    assert S0 == -22*n-6-(c % 8)
    record = {"n": n, "leaf_indices": [s["j"] for s in ls],
              "degree_in_v_checked": degree, "S0": S0,
              "S0_valuation": vp(S0, 2)}
    if n % 128 == 69 or n % 256 == 162:
        assert vp(ev(Q0, -8*n), 2) == 4
        assert all(vp(ev(g, -8*n), 2) >= 2 for g in Gs)
        for q in Qs[1:]:
            assert ev(q, -8*n) == 0 or vp(ev(q, -8*n), 2) >= 5
            assert q[6] == 0 or vp(q[6], 2) >= 3
        E = vp((11*n+9)*(11*n+10), 2)
        L = E+1-vp(S0, 2)
        record.update({"E": E, "forced_L": L})
        if L >= 8:
            kappa = Fraction(8, 11*(11*n+9)*(11*n+10))
            G = []
            for k in range(degree+1):
                upper = 0 if k == 0 else 128*kappa*(S0 if k == 1 else Qs[k-1][6])
                G.append(Fraction(ev(Qs[k], -8*n)+upper)*2**((L-7)*k)/16)
            assert vp(G[0], 2) == vp(G[1], 2) == 0
            assert all(x == 0 or vp(x, 2) >= 1 for x in G[2:])
            record["normalized_mixed_series_valuations"] = [vp(x, 2) for x in G]
    return record


def check_full_operator_split():
    # Independent exact comparison: split the original polynomial output
    # and compare with the composed leaf operators on the split seed1.
    controls = 0
    Q = [2, -3, 1]
    for n in range(1, 5):
        c = 3*n+8
        Pi = [1]
        for ell in range(n):
            Pi = mul(Pi, [8*ell, 1])
        poly = mul(Pi, Q)
        original = polynomial_of_operator(poly, c, [1])
        split_output = leaves(c, 2, initial_seed=original)
        split_seed = leaves(c, 2)
        for lhs, seed in zip(split_output, split_seed):
            rhs = leaf_poly(poly, seed["h"], c-seed["j"], seed["F"], 2)
            assert lhs["h"] == rhs
            controls += 1
    return controls


def main():
    mixed_controls = 0
    for d in [8]:
        for n in range(1, 13):
            check_mixed(n, d)
            mixed_controls += 1
    rank_controls = []
    for d, p in [(8, 2)]:
        for n in range(1, 41):
            v = [1]+[0]*(d-1)
            powers = [v]
            for _ in range(n+d):
                v = reduced_A(v, d, p, 3*n+d)
                powers.append(v)
            rank = rank_mod(powers[n+1:n+d], p)
            expected = d-min(n+1, d//p)
            assert rank == expected
            if n <= 4:
                rank_controls.append({"d": d, "p": p, "n": n, "rank": rank})
    radius_controls = 0
    for n in range(1, 101):
        for d, p, b in [(8, 2, 7)]:
            shape = 3*n+1
            rising = 1
            for k in range(1, n+1):
                for j in range(d*(k-1), d*k):
                    rising *= shape+j
                val = vp(comb(n, k)*rising, p)
                assert val >= b*k
                if d == 6 or n % 2 == 0 or n % 64 == 5:
                    assert val > b*k
                radius_controls += 1
    # Congruence equivalences, including the parity conditions.
    odd = [n for n in range(256) if n % 2 and (11*n+9) % 128 == 0]
    even = [n for n in range(256) if n % 2 == 0 and (11*n+10) % 256 == 0]
    assert odd == [69, 197] and even == [162]
    terminals = sorted((n+1) % 256 for n in odd+even)
    assert terminals == [70, 163, 198]
    # Derivative binomial coefficients used in the integral-ideal proof.
    assert [comb(j, 3) for j in range(3, 8)] == [1, 4, 10, 20, 35]
    assert [comb(j, 2) for j in range(2, 8)] == [1, 3, 6, 10, 15, 21]
    split_controls = check_full_operator_split()
    full_leaf_records = [leaf_controls(n) for n in [1, 2, 3, 4, 5, 6, 7, 8, 69, 197, 162, 930, 1954]]
    report = {
        "status": "PASS", "scope": "exact finite algebra controls; universal proof is in companion note",
        "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
        "mixed_triangularization_controls": mixed_controls,
        "actual_native_rank_controls": 40,
        "sample_actual_native_ranks": rank_controls,
        "actual_factorial_coefficient_radius_controls": radius_controls,
        "sharp_stage_control": sharpness(),
        "residue_stage_candidate_lower_indices_mod_256": sorted(odd+even),
        "residue_stage_candidate_terminal_indices_mod_256": terminals,
        "exact_three_stage_operator_controls": split_controls,
        "actual_full_leaf_interpolation_controls": full_leaf_records,
        "universal_closure_proved_by_finite_checks": False,
        "full_gap_eight_proof_location": "STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md, Sections 1-4",
    }
    path = Path(__file__).with_name("stepanov_lower_three_eighth_gap_closure_verification.json")
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    main()
