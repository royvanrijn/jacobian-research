#!/usr/bin/env python3
"""Generate an exact or finite-field Singular replay for the Case-1 rank-drop argument.

The input is the pinned exact replay archive.  This script independently
reconstructs the quotient-first degree-five equations, writes the 2x4 matrix
M=(a_i;b_i) with F_i=a_i*N+b_i*h.  The default exact replay
certifies the minimal adjacent-minor identity
    1 in (m_12,m_23,m_34,F_1),
where m_ij=a_i*b_j-a_j*b_i.  Multiplication by h and
    h*m_ij=a_i*F_j-a_j*F_i
then gives a direct certificate h in (F_1,...,F_4).  A separate exact
special-fibre Bezout identity gives N modulo h, so the hard residual ideal
is exactly (h,N).
"""
from __future__ import annotations

import argparse
import pickle
import re
from fractions import Fraction as Q
from pathlib import Path

import sympy as sp

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "plane-jc/external/zenodo-21479814/bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/release_bundle/exact_replay"
MINPOLY = "w^5-w^4+3*w^3+3*w^2+26"
VARIABLES = ("h", "u1", "u2")
DEFAULT_PAIRS = ((1, 2), (2, 3), (3, 4))
ALL_PAIRS = tuple((i, j) for i in range(1, 5) for j in range(i + 1, 5))


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def mul(a, b):
    if not a or not b:
        return []
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def mod_monic(a, modulus):
    a = list(a)
    modulus = trim(modulus)
    degree = len(modulus) - 1
    assert modulus[-1] == 1
    while len(a) > degree:
        exponent = len(a) - 1
        lead = a[exponent]
        if lead:
            for j in range(degree):
                a[exponent - degree + j] -= lead * modulus[j]
        a.pop()
    return trim(a)


FIELD_MODULUS = [Q(26), Q(0), Q(3), Q(3), Q(-1), Q(1)]


class L:
    __slots__ = ("a",)

    def __init__(self, value=0):
        if isinstance(value, L):
            self.a = value.a
        elif isinstance(value, (int, Q)):
            self.a = tuple(trim([Q(value)]))
        else:
            self.a = tuple(mod_monic([Q(x) for x in value], FIELD_MODULUS))

    def __add__(self, other):
        return L(add(self.a, L(other).a))

    __radd__ = __add__

    def __neg__(self):
        return L([-x for x in self.a])

    def __sub__(self, other):
        return self + (-L(other))

    def __mul__(self, other):
        return L(mod_monic(mul(self.a, L(other).a), FIELD_MODULUS))

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inv() ** (-exponent)
        result, base = L(1), self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def inv(self):
        if not self:
            raise ZeroDivisionError
        columns = []
        for j in range(5):
            col = list((self * L([0] * j + [1])).a) + [Q(0)] * 5
            columns.append(col[:5])
        matrix = sp.Matrix(
            [[sp.Rational(columns[j][i].numerator, columns[j][i].denominator) for j in range(5)] for i in range(5)]
        )
        solution = matrix.LUsolve(sp.Matrix([1, 0, 0, 0, 0]))
        return L([Q(int(x.p), int(x.q)) for x in solution])

    def __truediv__(self, other):
        return self * L(other).inv()

    def __bool__(self):
        return any(self.a)

    def __eq__(self, other):
        return self.a == L(other).a


def parse_degree35_eliminant():
    line = next(line for line in (ROOT / "firstblock_Q_exact.out").read_text().splitlines() if line.startswith("L[1]="))
    raw = {}
    for term in line.split("=", 1)[1].replace("-", "+-").split("+"):
        if not term:
            continue
        if "*a7^" in term:
            coefficient, exponent = term.split("*a7^")
            raw[int(exponent)] = Q(int(coefficient))
        else:
            raw[0] = Q(int(term))
    lead = raw[35]
    polynomial = {e: c / lead for e, c in raw.items()}
    assert set(polynomial) == {0, 7, 14, 21, 28, 35}
    return polynomial


H = parse_degree35_eliminant()
PHI = L([
    Q(-9725570295901, 12623962),
    Q(-1170753213563, 971074),
    Q(-387111042229, 12623962),
    Q(1578225240619, 12623962),
    Q(-469713794365, 6311981),
])


def decode_k(serialized):
    return {int(e): Q(int(n), int(d)) for e, (n, d) in serialized.items() if n}


def multiply_by_u_and_reduce(coefficient, exponent):
    polynomial = {e + exponent: value for e, value in coefficient.items()}
    while polynomial and max(polynomial) >= 35:
        e = max(polynomial)
        value = polynomial.pop(e)
        shift = e - 35
        for j in (0, 7, 14, 21, 28):
            target = shift + j
            polynomial[target] = polynomial.get(target, Q(0)) - value * H[j]
            if not polynomial[target]:
                polynomial.pop(target, None)
    return polynomial


def descend_coefficient(serialized, exponent):
    polynomial = multiply_by_u_and_reduce(decode_k(serialized), exponent)
    assert all(e % 7 == 0 for e in polynomial)
    result = L(0)
    for j in range(max((e // 7 for e in polynomial), default=0), -1, -1):
        result = result * PHI + L(polynomial.get(7 * j, Q(0)))
    return result


def transform_system(system):
    weights = (5, 6, 5, 6)
    characters = (6, 1, 0, 3, 2, 1, 0)
    transformed = []
    for equation, character in zip(system, characters):
        out = {}
        for monomial, coefficient in equation.items():
            exponent = sum(w * e for w, e in zip(weights, monomial)) - character
            value = descend_coefficient(coefficient, exponent)
            if value:
                out[tuple(monomial)] = value
        transformed.append(out)
    return transformed


def padd(a, b):
    out = dict(a)
    for monomial, coefficient in b.items():
        value = out.get(monomial, L(0)) + coefficient
        if value:
            out[monomial] = value
        else:
            out.pop(monomial, None)
    return out


def pmul(a, b):
    out = {}
    for m, c in a.items():
        for n, d in b.items():
            monomial = tuple(x + y for x, y in zip(m, n))
            out[monomial] = out.get(monomial, L(0)) + c * d
    return {m: c for m, c in out.items() if c}


def reconstruct_matrix():
    transformed = transform_system(pickle.loads((ROOT / "case1_branch1_after_w.pkl").read_bytes()))
    q0 = transformed[0]
    u3_coefficient = {m[:-1]: c for m, c in q0.items() if m[-1] == 1}
    rest = {m[:-1]: c for m, c in q0.items() if m[-1] == 0}
    assert set(u3_coefficient) == {(1, 0, 0)}
    bh = u3_coefficient[(1, 0, 0)]
    n_poly = {m: -c / bh for m, c in rest.items()}
    h_poly = {(1, 0, 0): L(1)}
    a_values, b_values, f_values = [], [], []
    for equation in transformed[1:]:
        a = {m[:-1]: c for m, c in equation.items() if m[-1] == 1}
        b = {m[:-1]: c for m, c in equation.items() if m[-1] == 0}
        f = padd(pmul(h_poly, b), pmul(a, n_poly))
        common_h = min(m[0] for m in f)
        if common_h:
            f = {(m[0] - common_h, m[1], m[2]): c for m, c in f.items()}
            b = {(m[0] - common_h, m[1], m[2]): c for m, c in b.items()}
            a = {(m[0] - common_h, m[1], m[2]): c for m, c in a.items()}
        assert common_h == 0
        a_values.append(a)
        b_values.append(b)
        f_values.append(f)

    stored = pickle.loads((ROOT / "hne0_polred.pkl").read_bytes())
    for index, (derived, serialized) in enumerate(zip(f_values, stored), start=1):
        expected = {
            tuple(m): L([Q(int(n), int(d)) if j in c else Q(0) for j in range(5) for n, d in [c.get(j, (0, 1))]])
            for m, c in serialized.items()
        }
        assert set(derived) == set(expected), f"support mismatch in residual {index}"
        monomial = next(iter(derived))
        ratio = derived[monomial] / expected[monomial]
        assert all(derived[m] == expected[m] * ratio for m in derived), f"coefficient mismatch in residual {index}"

    return n_poly, a_values[:4], b_values[:4], f_values[:4]


def special_fibre_bezout(a_values):
    """Return constants c1,c2 with c1*a1+c2*a2 = 1 modulo h."""
    allowed = {(0, 0, 0), (0, 0, 1)}
    for index in (0, 1):
        support = {m for m in a_values[index] if m[0] == 0}
        assert support <= allowed and support == allowed
    beta1 = a_values[0][(0, 0, 0)]
    alpha1 = a_values[0][(0, 0, 1)]
    beta2 = a_values[1][(0, 0, 0)]
    alpha2 = a_values[1][(0, 0, 1)]
    delta = alpha1 * beta2 - alpha2 * beta1
    assert delta
    c1 = -alpha2 / delta
    c2 = alpha1 / delta
    assert c1 * alpha1 + c2 * alpha2 == 0
    assert c1 * beta1 + c2 * beta2 == 1
    return c1, c2, delta


def rational_string(value):
    return str(value.numerator) if value.denominator == 1 else f"({value.numerator}/{value.denominator})"


def field_string(value, prime=None):
    terms = []
    for exponent, rational in enumerate(value.a):
        if prime is None:
            coefficient = rational
            if not coefficient:
                continue
            text = rational_string(coefficient)
        else:
            coefficient = (rational.numerator % prime) * pow(rational.denominator % prime, -1, prime) % prime
            if not coefficient:
                continue
            text = str(coefficient)
        variable = "" if exponent == 0 else ("w" if exponent == 1 else f"w^{exponent}")
        terms.append(text if not variable else (variable if text == "1" else f"{text}*{variable}"))
    return "+".join(terms) or "0"


def polynomial_string(poly, prime=None):
    terms = []
    for monomial, coefficient in sorted(poly.items(), key=lambda item: (sum(item[0]), item[0]), reverse=True):
        c = field_string(coefficient, prime)
        if c == "0":
            continue
        variable = "*".join(v if e == 1 else f"{v}^{e}" for v, e in zip(VARIABLES, monomial) if e)
        terms.append(c if not variable else (variable if c == "1" else f"({c})*{variable}"))
    return "+".join(terms) or "0"


def parse_pairs(text):
    if text == "all":
        return ALL_PAIRS
    pairs = []
    for token in text.split(","):
        token = token.strip().lower().removeprefix("m")
        if len(token) != 2 or not token.isdigit():
            raise argparse.ArgumentTypeError(f"invalid minor {token!r}; use e.g. 12,23,34")
        i, j = map(int, token)
        if not (1 <= i < j <= 4):
            raise argparse.ArgumentTypeError(f"invalid minor m{i}{j}")
        pairs.append((i, j))
    if not pairs:
        raise argparse.ArgumentTypeError("at least one minor is required")
    return tuple(pairs)


def write_singular(path, prime=None, lift=False, pairs=DEFAULT_PAIRS, write_certificate=None):
    n_poly, a_values, b_values, _ = reconstruct_matrix()
    sf_c1, sf_c2, _ = special_fibre_bezout(a_values)
    field = "(0,w)" if prime is None else f"({prime},w)"
    lines = [f"ring R={field},(h,u1,u2),dp;", f"minpoly={MINPOLY};", "option(redSB);"]
    lines.append(f"poly N={polynomial_string(n_poly, prime)};")
    for i in range(4):
        lines.append(f"poly a{i+1}={polynomial_string(a_values[i], prime)};")
        lines.append(f"poly b{i+1}={polynomial_string(b_values[i], prime)};")
        lines.append(f"poly f{i+1}=h*b{i+1}+a{i+1}*N;")
    for i, j in ALL_PAIRS:
        lines.append(f"poly m{i}{j}=a{i}*b{j}-a{j}*b{i};")
    selected = [f"m{i}{j}" for i, j in pairs]
    lines.append("ideal K=" + ",".join(selected) + ",f1;")
    lines += [
        f"number sf_c1={field_string(sf_c1, prime)};",
        f"number sf_c2={field_string(sf_c2, prime)};",
        "poly sfcheck=subst(sf_c1*a1+sf_c2*a2,h,0);",
        'if(sfcheck==1){print("CASE1_SPECIAL_FIBRE_BEZOUT_PASS");}else{print("CASE1_SPECIAL_FIBRE_BEZOUT_FAIL");}',
        "poly nmodh=subst(sf_c1*f1+sf_c2*f2-N,h,0);",
        'if(nmodh==0){print("CASE1_N_MOD_H_PASS");}else{print("CASE1_N_MOD_H_FAIL");}',
    ]
    if prime is None and not lift:
        lines += [
            'LIB "nfmodstd.lib";',
            "ideal G=nfmodStd(K);",
            'if(size(G)==1 && G[1]==1){print("RANKDROP_EXACT_UNIT_PASS");}else{print("RANKDROP_EXACT_UNIT_FAIL"); G;}',
        ]
    else:
        f1_index = len(pairs) + 1
        lines += [
            "matrix T; ideal G=liftstd(K,T);",
            "int jj=0; int k; for(k=1;k<=size(G);k++){if(G[k]==1||G[k]==-1){jj=k;}}",
            'if(jj==0){print("RANKDROP_LIFT_FAIL"); quit;}',
            "if(G[jj]==-1){for(k=1;k<=nrows(T);k++){T[k,jj]=-T[k,jj];}}",
            f"poly check1=0; for(k=1;k<={f1_index};k++){{check1=check1+K[k]*T[k,jj];}}",
            'if(check1==1){print("RANKDROP_UNIT_IDENTITY_PASS");}else{print("RANKDROP_UNIT_IDDITAY_FAIL");}',
            f"poly t1=h*T[{f1_index},jj]; poly t2=0; poly t3=0; poly t4=0;",
        ]
        for index, (i, j) in enumerate(pairs, start=1):
            lines.append(f"t{i}=t{i}-T[{index},jj]*a{j}; t{j}=t{j}+T[{index},jj]*a{i};")
        lines += [
            "poly checkh=t1*f1+t2*f2+t3*f3+t4*f4;",
            'if(checkh==h){print("CASE1_H_DETERMINANTAL_PASS");}else{print("CASE1_H_DETERMINANTAL_FAIL");}',
            'print("T1_TERMS="+string(size(t1))); print("T2_TERMS="+string(size(t2))); print("T3_TERMS="+string(size(t3))); print("T4_TERMS="+string(size(t4)));',
        ]
        if write_certificate is not None:
            destination = str(write_certificate).replace('\\', '/').replace('"', '\\"')
            lines.append(
                f'write("{destination}","t1="+string(t1)+";\\nt2="+string(t2)+";\\nt3="+string(t3)+";\\nt4="+string(t4)+";\\ne");'
            )
    lines.append("quit;")
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="generated Singular replay")
    parser.add_argument("--prime", type=int, help="reduce coefficients modulo this good prime")
    parser.add_argument("--lift", action="store_true", help="lift an explicit unit and h identity")
    parser.add_argument(
        "--pairs",
        type=parse_pairs,
        default=DEFAULT_PAIRS,
        help="comma-separated minors, default 12,23,34; use all for all six",
    )
    parser.add_argument("--write-certificate", type=Path, help="write the four final h multipliers")
    args = parser.parse_args()
    write_singular(args.output, args.prime, args.lift, args.pairs, args.write_certificate)
    print("CASE1_RECONSTRUCTION_PASS")
    print("CASE1_SPECIAL_FIBRE_PYTHON_PASS")
    print(args.output, args.output.stat().st_size)


if __name__ == "__main__":
    main()
