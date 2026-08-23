#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from __future__ import print_function

from sage.all import GF, FunctionField, PolynomialRing


F = GF(73, impl="modn")
KT = FunctionField(F, "T")
T = KT.gen()
R = PolynomialRing(KT, "V")
V = R.gen()

# Old pinned 2I7+3I2+4I1 model.
A = R([
    23,17,63,2,58,33,47,16,6
])
B = R([
    43,47,57,50,8,54,20,14,45,61,64,0,33
])

# Certified specialized horizontal P3.
PX = R([47,69,7,4,51])
PY0 = R([22,2,48,23,8,47,35])
assert PY0**2 == PX**3+A*PX+B

# Exact local rows already certified by the resolved/nodal preflight.
i7_roots = (F(6), F(53))
i2_data = {
    +1: {
        F(5): F(62),
        F(43): F(0),
        F(56): F(7),
    },
    -1: {
        F(5): F(11),
        F(43): F(0),
        F(56): F(66),
    },
}

TARGET = (15,74,600)


def branch_for(m, PY):
    return (
        m**4
        - 6*PX*m**2
        - 8*PY*m
        - 3*PX**2
        - 4*A
    )


def odd_part(poly):
    factorization = poly.factor()
    odd = R.one()
    profile = []
    for factor, exponent in factorization:
        exponent = int(exponent)
        profile.append(
            (str(factor), int(factor.degree()), exponent)
        )
        if exponent % 2:
            odd *= factor
    return odd, tuple(profile), factorization


def kodaira_data(ord_a,ord_b,ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        if n == 1:
            return 0,0,1,1,"I1"
        return n-1,n*(n-1),n,n,f"I{n}"
    if ord_delta == 2:
        return 0,0,1,2,"II"
    if ord_delta == 3:
        return 1,2,2,3,"III"
    if ord_delta == 4:
        return 2,6,3,4,"IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4,24,4,6,"I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta-6)
        rank = n+4
        return rank,2*rank*(rank-1),4,n+6,f"I{n}*"
    if ord_delta == 8:
        return 6,72,3,8,"IV*"
    if ord_delta == 9:
        return 7,126,2,9,"III*"
    if ord_delta == 10:
        return 8,240,1,10,"II*"
    raise ArithmeticError((ord_a,ord_b,ord_delta))


def classify_quartic(quartic,twist):
    coeffs = list(quartic.list())+[KT(0)]*5
    e,d,c,b,a = coeffs[:5]

    Iinv = 12*a*e-3*b*d+c**2
    Jinv = (
        72*a*c*e+9*b*c*d
        -27*a*d**2-27*b**2*e-2*c**3
    )

    jacA = twist**2*(-27*Iinv)
    jacB = twist**3*(-27*Jinv)
    delta = twist**6*(4*Iinv**3-Jinv**2)

    root_rank=root_count=euler=0
    root_det=1
    finite=[]
    factors=set()

    for value in (jacA,jacB,delta):
        for polynomial in (
            value.numerator(),value.denominator()
        ):
            factors.update(
                factor for factor,_ in polynomial.factor()
            )

    for factor in sorted(factors,key=str):
        raw = (
            int(jacA.valuation(factor)),
            int(jacB.valuation(factor)),
            int(delta.valuation(factor)),
        )
        scale = min(
            raw[0]//4,raw[1]//6,raw[2]//12
        )
        orders = (
            raw[0]-4*scale,
            raw[1]-6*scale,
            raw[2]-12*scale,
        )
        if orders[2] == 0:
            continue

        rank,count,determinant,local_euler,kind = (
            kodaira_data(*orders)
        )
        degree = int(factor.degree())
        root_rank += degree*rank
        root_count += degree*count
        root_det *= determinant**degree
        euler += degree*local_euler
        finite.append(
            (str(factor),degree,orders,kind)
        )

    raw_inf = tuple(
        int(
            value.denominator().degree()
            - value.numerator().degree()
        )
        for value in (jacA,jacB,delta)
    )
    scale = min(
        raw_inf[0]//4,
        raw_inf[1]//6,
        raw_inf[2]//12,
    )
    inf_orders = (
        raw_inf[0]-4*scale,
        raw_inf[1]-6*scale,
        raw_inf[2]-12*scale,
    )
    inf_kind = "smooth"
    if inf_orders[2] > 0:
        rank,count,determinant,local_euler,inf_kind = (
            kodaira_data(*inf_orders)
        )
        root_rank += rank
        root_count += count
        root_det *= determinant
        euler += local_euler

    return {
        "root_data": (
            int(root_rank),int(root_count),int(root_det)
        ),
        "MW": int(18-root_rank),
        "euler": int(euler),
        "finite": tuple(finite),
        "infinity_orders": tuple(inf_orders),
        "infinity": inf_kind,
        "A": str(jacA),
        "B": str(jacB),
    }


def interpolation_linear(x0,y0,x1,y1):
    assert x0 != x1
    slope = (y1-y0)/(x1-x0)
    intercept = y0-slope*x0
    return R([KT(intercept),KT(slope)])


hits = []
quartic_count = 0
tested = 0

for sign in (+1,-1):
    PY = PY0 if sign == 1 else -PY0

    for i7 in i7_roots:
        for i2,c2 in sorted(i2_data[sign].items(),key=lambda item:int(item[0])):
            assert i7 != i2

            per_block_quartics = 0
            per_block_targets = 0

            # The genuine connected-A6 quotient is one projective line in
            # the two-row componentwise span.  All nondegenerate such lines
            # are exactly evaluation rows (1,i7,i7^2,c7), c7 in F_73.
            #
            # With the I2 row (1,i2,i2^2,c2), a convenient exact kernel is:
            #   f0=(V-i7)(V-i2)
            #   f1=m+a(V),
            # where a(i7)=-c7 and a(i2)=-c2.
            # Hence T=f1/f0 and
            #   m=T*(V-i7)(V-i2)-a(V).
            for c7 in F:
                tested += 1

                p = (V-KT(i7))*(V-KT(i2))
                a = interpolation_linear(
                    i7,-c7,
                    i2,-c2,
                )
                m = T*p-a

                branch = branch_for(m,PY)
                assert branch.denominator() == 1
                square_class = branch.numerator()

                odd,profile,factorization = odd_part(
                    square_class
                )
                odd_degree = int(odd.degree())
                if odd_degree != 4:
                    continue

                quartic_count += 1
                per_block_quartics += 1

                quartic = odd.monic()
                twist = KT(factorization.unit())
                classification = classify_quartic(
                    quartic,twist
                )
                assert classification["euler"] == 24

                if (
                    classification["root_data"] == TARGET
                    and classification["MW"] == 3
                ):
                    per_block_targets += 1
                    hit = {
                        "sign":sign,
                        "i7":int(i7),
                        "i2":int(i2),
                        "c7":int(c7),
                        "c2":int(c2),
                        "a":str(a),
                        "quartic":str(quartic),
                        "twist":str(twist),
                        "classification":classification,
                        "profile":profile,
                    }
                    hits.append(hit)

                    print(
                        f"Q807774LINEHIT|sign={sign}|"
                        f"I7={int(i7)}|I2={int(i2)}|"
                        f"c7={int(c7)}|c2={int(c2)}|"
                        f"a={a}|quartic={quartic}|twist={twist}|"
                        f"finite={classification['finite']}|"
                        f"infinity={classification['infinity']}|"
                        f"root_data={classification['root_data']}|"
                        f"MW={classification['MW']}|"
                        "status=PASS_TRUE_A6_QUOTIENT_LINE",
                        flush=True,
                    )

            print(
                f"Q807774LINEBLOCK|sign={sign}|"
                f"I7={int(i7)}|I2={int(i2)}|c2={int(c2)}|"
                f"quartic_residues={per_block_quartics}|"
                f"target_residues={per_block_targets}|"
                "status=PASS_73_LINE_SCAN",
                flush=True,
            )


print(
    f"Q807774LINE|tested={tested}|"
    f"quartic_count={quartic_count}|"
    f"target_hits={len(hits)}|"
    f"hit_keys={tuple((h['sign'],h['i7'],h['i2'],h['c7']) for h in hits)}|"
    "status=PASS_CONNECTED_A6_LINE_SCAN",
    flush=True,
)

assert hits, "No projective A6 quotient line reached (15,74,600)/MW3"

# Pin the simplest target marking.  If automorphisms create multiple hits,
# preserve all hit keys in the previous line and choose deterministically.
chosen = min(
    hits,
    key=lambda h: (
        h["sign"] != 1,
        h["i7"],
        h["i2"],
        h["c7"],
    ),
)
classification = chosen["classification"]

print(
    f"Q807774LINESELECT|sign={chosen['sign']}|"
    f"I7={chosen['i7']}|I2={chosen['i2']}|"
    f"c7={chosen['c7']}|c2={chosen['c2']}|"
    f"a={chosen['a']}|"
    f"quartic={chosen['quartic']}|"
    f"twist={chosen['twist']}|"
    f"root_data={classification['root_data']}|"
    f"MW={classification['MW']}|"
    f"A={classification['A']}|"
    f"B={classification['B']}|"
    "status=PASS_PINNED_Q6_7774_CM24_EQUATION",
    flush=True,
)

print(
    "Q807774LINE|next=derive_transition_rule_for_selected_A6_line_"
    "then_q4_1938|"
    "status=PASS_READY_FOR_Q4_1938",
    flush=True,
)
"""


def main():
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    print(f"sage={sage}", flush=True)
    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_q6_7774_connected_a6_line_scan.sage"
        script.write_text(SAGE_CODE)
        subprocess.run([sage, str(script)], check=True)


if __name__ == "__main__":
    main()
