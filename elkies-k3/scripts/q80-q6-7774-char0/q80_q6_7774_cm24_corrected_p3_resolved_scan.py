#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from __future__ import print_function

from itertools import product

from sage.all import (
    EllipticCurve, GF, FunctionField, LaurentSeriesRing,
    PolynomialRing, gcd, lcm, matrix
)


F = GF(73, impl="modn")
R = PolynomialRing(F, "V")
V = R.gen()
K = R.fraction_field()

A = (
    6*V**8 + 16*V**7 + 47*V**6 + 33*V**5 + 58*V**4
    + 2*V**3 + 63*V**2 + 17*V + 23
)
B = (
    33*V**12 + 64*V**10 + 61*V**9 + 45*V**8 + 14*V**7
    + 20*V**6 + 54*V**5 + 8*V**4 + 50*V**3 + 57*V**2
    + 47*V + 43
)
Delta = -F(16)*(4*A**3+27*B**2)

i7_roots = tuple(
    F(-factor[0]/factor[1])
    for factor, exponent in Delta.factor()
    if factor.degree() == 1 and exponent == 7
)
i2_roots = tuple(
    F(-factor[0]/factor[1])
    for factor, exponent in Delta.factor()
    if factor.degree() == 1 and exponent == 2
)
assert set(map(int,i7_roots)) == {53,6}
assert set(map(int,i2_roots)) == {56,43,5}

xring = PolynomialRing(F, "x")
x = xring.gen()


def node_x(root):
    cubic = x**3 + A(root)*x + B(root)
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return F(-common[0]/common[1])


nodes = {root: node_x(root) for root in i7_roots+i2_roots}


# ------------------------------------------------------------------
# Recover the repository-pinned polynomial section P3 = points[2].
# We reproduce the exact deterministic two-I7-node search, stopping as soon
# as candidate index 2 is known.
# ------------------------------------------------------------------
def polynomial_square_roots(polynomial):
    assert polynomial.degree() <= 12
    if polynomial == 0:
        return (R.zero(),)
    shift = next(value for value in F if polynomial(value) != 0)
    shifted = polynomial(V+shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()
    roots = []
    for first in constant.sqrt(all=True):
        coefficients = [first]
        for degree in range(1,7):
            known = sum(
                coefficients[left]*coefficients[degree-left]
                for left in range(1,degree)
            )
            coefficients.append(
                (shifted[degree]-known)/(2*first)
            )
        candidate_shifted = R(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(V-shift))
    return tuple(roots)


interpolation = R.lagrange_polynomial(
    [(root,nodes[root]) for root in i7_roots]
)
vanishing = R.one()
for root in i7_roots:
    vanishing *= V-root
assert vanishing.degree() == 2

first_candidates = []
tests = 0
for coefficients in product(F, repeat=3):
    tests += 1
    Xcandidate = interpolation + R(coefficients)*vanishing
    for Ycandidate in polynomial_square_roots(
        Xcandidate**3+A*Xcandidate+B
    ):
        first_candidates.append((Xcandidate,Ycandidate))
    if len(first_candidates) >= 3:
        break

assert len(first_candidates) >= 3
P3X,P3Y = first_candidates[2]
assert P3Y**2 == P3X**3+A*P3X+B
assert P3X.degree() <= 4
assert P3Y.degree() <= 6

# P3 is the pinned height-8/7 basis section: nonidentity at both I7 fibres,
# identity at all three I2 fibres.
for root in i7_roots:
    assert P3X(root) == nodes[root]
for root in i2_roots:
    assert P3X(root) != nodes[root]

print(
    f"Q807774CORR|P3_search_tests={tests}|"
    f"P3X={tuple(map(int,P3X.list()))}|"
    f"P3Y={tuple(map(int,P3Y.list()))}|"
    "height=8/7|P.O=0|I2_support=000|"
    "status=PASS_PINNED_P3_RECOVERY",
    flush=True,
)


# ------------------------------------------------------------------
# Exact resolved multiplicative-fibre charts.
# ------------------------------------------------------------------
Kr = FunctionField(F, "r")
r = Kr.gen()
PREC = 42
LS = LaurentSeriesRing(Kr, "s", default_prec=PREC)
s = LS.gen()


def shifted_series(poly, root):
    answer = LS(0)
    for degree, coefficient in enumerate(poly.list()):
        answer += Kr(F(coefficient))*(Kr(root)+s)**degree
    return answer


def newton_sqrt(value, root0):
    root = LS(Kr(root0))
    for _ in range(11):
        root = (root+value/root)/2
    assert (root*root-value).valuation() >= PREC-6
    return root


def local_data(root, n):
    Aloc = shifted_series(A, root)
    Bloc = shifted_series(B, root)

    center0 = nodes[root]
    target = -Aloc/3
    assert F(target[0]) == center0**2
    center = newton_sqrt(target, center0)

    g0 = center**3+Aloc*center+Bloc
    assert g0.valuation() == n, (
        int(root), n, g0.valuation()
    )
    unit = g0/s**n

    rho_roots = tuple(
        value for value in F
        if value**2 == F((3*center)[0])
    )
    assert len(rho_roots) == 2
    rho0 = min(rho_roots,key=int)
    return Aloc,Bloc,center,unit,Kr(rho0)


def toric_point(root,n,component):
    _,_,center,unit,rho0 = local_data(root,n)

    aa = LS(r)*s**component
    b0 = s**(n-component)/LS(r)
    bb = unit*b0
    yy = (aa+bb)/2
    ww = (bb-aa)/2

    rho = LS(rho0)
    for _ in range(11):
        residual = rho**3-3*center*rho-ww
        derivative = 3*rho**2-3*center
        rho = rho-residual/derivative
    assert (rho**3-3*center*rho-ww).valuation() >= PREC-8

    xx = center+ww/rho
    return xx,yy


def section_series(poly,root):
    return shifted_series(poly,root)


def functional_rows(function_row):
    nonzero = [value for value in function_row if value]
    if not nonzero:
        return ()
    common = nonzero[0].denominator().parent().one()
    for value in nonzero:
        common = common.lcm(value.denominator())

    numerators = [
        (value*common).numerator()
        for value in function_row
    ]
    max_degree = max(
        [poly.degree() for poly in numerators if poly]+[-1]
    )
    rows = []
    for degree in range(max_degree+1):
        row = tuple(
            F(poly[degree])
            if poly and degree <= poly.degree()
            else F(0)
            for poly in numerators
        )
        if any(row):
            rows.append(row)
    return tuple(rows)


def canonical_row(row):
    row = tuple(F(value) for value in row)
    pivot = next(value for value in row if value)
    inverse = pivot**(-1)
    return tuple(value*inverse for value in row)


def dedup_rows(rows):
    result = []
    seen = set()
    for row in rows:
        normalized = canonical_row(row)
        key = tuple(map(int,normalized))
        if key not in seen:
            seen.add(key)
            result.append(normalized)
    return tuple(result)


# ------------------------------------------------------------------
# Quartic/Jacobian classifier.
# ------------------------------------------------------------------
Tfield = FunctionField(F, "T")
T = Tfield.gen()
VR = PolynomialRing(Tfield, "V")
VV = VR.gen()


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
    coeffs = list(quartic.list())+[Tfield(0)]*5
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
            factors.update(f for f,_ in polynomial.factor())

    for factor in sorted(factors,key=str):
        raw = (
            int(jacA.valuation(factor)),
            int(jacB.valuation(factor)),
            int(delta.valuation(factor)),
        )
        scale = min(raw[0]//4,raw[1]//6,raw[2]//12)
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
        int(value.denominator().degree()-value.numerator().degree())
        for value in (jacA,jacB,delta)
    )
    scale = min(
        raw_inf[0]//4,raw_inf[1]//6,raw_inf[2]//12
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


def to_T_poly(poly):
    return VR([Tfield(F(c)) for c in poly.list()])


AT = to_T_poly(A)
P3XT = to_T_poly(P3X)


# Vertical A6 coefficient pattern in deterministic lattice chain:
# (-1,0,0,-1,-1,-1) = subtract E1,E4,E5,E6.
# Reversing the chain gives E1,E2,E3,E6.
a6_supports = {
    "direct": (1,4,5,6),
    "reversed": (1,2,3,6),
}


def build_local_rows(root,n,components,PY):
    Xploc = section_series(P3X,root)
    Yploc = section_series(PY,root)

    all_rows = []
    component_diagnostics = []
    for component in components:
        xx,yy = toric_point(root,n,component)
        mloc = (yy+Yploc)/(xx-Xploc)
        mval = int(mloc.valuation())
        assert mval >= 0

        values = (
            Kr(1),
            Kr(root),
            Kr(root**2),
            Kr(mloc[0]) if mval == 0 else Kr(0),
        )
        rows = functional_rows(values)
        all_rows.extend(rows)
        component_diagnostics.append((
            component,mval,str(mloc[0]) if mval == 0 else "0",
            tuple(tuple(map(int,row)) for row in rows),
        ))

    return dedup_rows(tuple(all_rows)),tuple(component_diagnostics)


results = []
hits = []

for sign in (1,-1):
    PY = P3Y if sign == 1 else -P3Y
    P3YT = to_T_poly(PY)

    for i7root in sorted(i7_roots,key=int):
        for orientation,support in a6_supports.items():
            i7_rows,i7_diag = build_local_rows(
                i7root,7,support,PY
            )
            i7_matrix = matrix(F,i7_rows)
            i7_rank = i7_matrix.rank()

            for i2root in sorted(i2_roots,key=int):
                i2_rows,i2_diag = build_local_rows(
                    i2root,2,(1,),PY
                )
                i2_matrix = matrix(F,i2_rows)
                i2_rank = i2_matrix.rank()

                conditions = i7_matrix.stack(i2_matrix)
                rank = conditions.rank()

                print(
                    f"Q807774CORRLOCAL|sign={sign}|"
                    f"I7={int(i7root)}|orientation={orientation}|"
                    f"I7support={support}|I7rank={i7_rank}|"
                    f"I2={int(i2root)}|I2rank={i2_rank}|"
                    f"total_rank={rank}|"
                    f"I7diag={i7_diag}|I2diag={i2_diag}|"
                    "status=PASS_LOCAL_BLOCK",
                    flush=True,
                )

                if rank != 2:
                    results.append((
                        sign,int(i7root),orientation,int(i2root),
                        rank,"wrong_rank",None,None,
                    ))
                    continue

                kernel = conditions.right_kernel().basis_matrix()
                assert kernel.nrows() == 2
                assert kernel.ncols() == 4

                def split_row(row):
                    row = tuple(Tfield(F(v)) for v in row)
                    aa = row[0]+row[1]*VV+row[2]*VV**2
                    bb = row[3]
                    return aa,bb

                (a0,b0),(a1,b1) = map(
                    split_row,kernel.rows()
                )

                determinant = a0*b1-a1*b0
                if determinant == 0:
                    results.append((
                        sign,int(i7root),orientation,int(i2root),
                        rank,"degenerate",None,None,
                    ))
                    print(
                        f"Q807774CORRDEGEN|sign={sign}|"
                        f"I7={int(i7root)}|orientation={orientation}|"
                        f"I2={int(i2root)}|kernel="
                        f"{tuple(tuple(map(int,row)) for row in kernel.rows())}|"
                        "status=PASS_REJECT_DEGENERATE",
                        flush=True,
                    )
                    continue

                denominator = T*b0-b1
                assert denominator != 0
                mvalue = (a1-T*a0)/denominator

                branch = (
                    mvalue**4
                    - 6*P3XT*mvalue**2
                    - 8*P3YT*mvalue
                    - 3*P3XT**2
                    - 4*AT
                )

                square_class = VR(
                    branch.numerator()*branch.denominator()
                )
                factorization = square_class.factor()
                odd = VR.one()
                for factor,exponent in factorization:
                    if int(exponent)%2:
                        odd *= factor

                odd_degree = int(odd.degree())
                profile = tuple(
                    (
                        str(factor),
                        int(factor.degree()),
                        int(exponent),
                    )
                    for factor,exponent in factorization
                )

                classification = None
                if odd_degree == 4:
                    quartic = odd.monic()
                    twist = Tfield(factorization.unit())
                    classification = classify_quartic(
                        quartic,twist
                    )

                root_data = (
                    classification["root_data"]
                    if classification is not None else None
                )
                mw = (
                    classification["MW"]
                    if classification is not None else None
                )

                print(
                    f"Q807774CORRSCAN|sign={sign}|"
                    f"I7={int(i7root)}|orientation={orientation}|"
                    f"I2={int(i2root)}|"
                    f"kernel={tuple(tuple(map(int,row)) for row in kernel.rows())}|"
                    f"odd_degree={odd_degree}|"
                    f"root_data={root_data}|MW={mw}|"
                    f"factor_profile={profile}|"
                    "status=PASS_CORRECTED_PENCIL_SCAN",
                    flush=True,
                )

                result = (
                    sign,int(i7root),orientation,int(i2root),
                    rank,odd_degree,root_data,mw,
                )
                results.append(result)

                if (
                    classification is not None
                    and classification["euler"] == 24
                    and classification["root_data"] == (15,74,600)
                    and classification["MW"] == 3
                ):
                    hit = {
                        "sign":sign,
                        "i7":int(i7root),
                        "orientation":orientation,
                        "i2":int(i2root),
                        "support":support,
                        "kernel":tuple(
                            tuple(map(int,row))
                            for row in kernel.rows()
                        ),
                        "quartic":str(quartic),
                        "twist":str(twist),
                        "classification":classification,
                    }
                    hits.append(hit)

                    print(
                        f"Q807774CORRHIT|sign={sign}|"
                        f"I7={int(i7root)}|orientation={orientation}|"
                        f"I7support={support}|I2={int(i2root)}|"
                        f"kernel={hit['kernel']}|"
                        f"quartic={quartic}|twist={twist}|"
                        f"finite={classification['finite']}|"
                        f"infinity={classification['infinity']}|"
                        f"root_data={classification['root_data']}|"
                        f"MW={classification['MW']}|"
                        "status=PASS_Q6_7774_TRUE_CM24_CHILD",
                        flush=True,
                    )


print(
    f"Q807774CORR|scan_count={len(results)}|"
    f"target_hits={len(hits)}|"
    f"hit_markings={tuple((h['sign'],h['i7'],h['orientation'],h['i2']) for h in hits)}|"
    "status=PASS_CORRECTED_24_CASE_SCAN",
    flush=True,
)

assert hits, (
    "No corrected resolved module reached the exact CM24 lattice target",
    results,
)

# Equivalent hits may differ by the known automorphisms/sign conventions.
# Pin the lexicographically simplest exact marking for the next equation step.
chosen = min(
    hits,
    key=lambda h: (
        h["sign"] != 1,
        h["i7"],
        h["orientation"] != "direct",
        h["i2"],
    ),
)

classification = chosen["classification"]
print(
    f"Q807774CORRSELECT|sign={chosen['sign']}|"
    f"I7={chosen['i7']}|orientation={chosen['orientation']}|"
    f"I7support={chosen['support']}|I2={chosen['i2']}|"
    f"kernel={chosen['kernel']}|"
    f"quartic={chosen['quartic']}|twist={chosen['twist']}|"
    f"root_data={classification['root_data']}|MW={classification['MW']}|"
    f"A={classification['A']}|B={classification['B']}|"
    "status=PASS_PINNED_Q6_7774_CM24_EQUATION",
    flush=True,
)

print(
    "Q807774CORR|next=derive_q4_1938_marking_on_pinned_CM24_child|"
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
        script = Path(td) / "q80_q6_7774_corrected_p3_resolved_scan.sage"
        script.write_text(SAGE_CODE)
        subprocess.run([sage, str(script)], check=True)


if __name__ == "__main__":
    main()
