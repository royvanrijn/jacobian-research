#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from __future__ import print_function

from sage.all import (
    GF, FunctionField, LaurentSeriesRing, PolynomialRing,
    gcd, matrix
)


F = GF(73, impl="modn")
PREC = 40

# ------------------------------------------------------------------
# Old pinned 2I7+3I2+4I1 CM24 model, now proved isomorphic to the
# orbit1222 CM24 child.
# ------------------------------------------------------------------
K = FunctionField(F, "T")
T = K.gen()
VR = PolynomialRing(K, "V")
V = VR.gen()
VF = VR.fraction_field()

A = VR([
    23,17,63,2,58,33,47,16,6
])
B = VR([
    43,47,57,50,8,54,20,14,45,61,64,0,33
])

i7_roots = (F(53), F(6))


# Four representatives, one from each +/- pair in the exact maximal-height
# q6_7774 specialized MW census.
candidates = (
    {
        "name": "pair1",
        "mw": (-5,4,0),
        "labels": ((0,2),(1,0,1)),
        "X_num": (37,48,50,12,68,33,0,44,72,13,41),
        "X_den": (1,28,9,8,9,41,1),
        "Y_num": (70,10,51,43,29,24,51,3,62,47,21,45,42,26,27,19),
        "Y_den": (72,31,21,3,55,6,19,69,25,1),
    },
    {
        "name": "pair2",
        "mw": (-5,4,1),
        "labels": ((5,0),(1,0,1)),
        "X_num": (69,40,25,32,64,49,4,62,2,67,4),
        "X_den": (12,38,59,36,11,13,1),
        "Y_num": (16,27,26,60,27,0,13,18,42,15,1,0,10,28,12,11),
        "Y_den": (66,58,46,37,29,1,65,16,56,1),
    },
    {
        "name": "pair3",
        "mw": (-1,-4,-1),
        "labels": ((0,2),(1,0,1)),
        "X_num": (40,57,24,38,67,58,65,36,39,69,53),
        "X_den": (3,14,14,58,64,65,1),
        "Y_num": (21,21,7,66,10,61,10,29,69,45,28,4,70,14,9,4),
        "Y_den": (10,70,30,14,6,26,27,47,61,1),
    },
    {
        "name": "pair4",
        "mw": (-1,-4,0),
        "labels": ((5,0),(1,0,1)),
        "X_num": (22,0,68,61,4,28,51,41,47,21,45),
        "X_den": (71,49,12,17,65,7,1),
        "Y_num": (65,4,59,53,21,12,70,20,55,67,28,19,0,40,5,18),
        "Y_den": (49,6,15,33,51,54,67,52,47,1),
    },
)


def kpoly(coeffs):
    return VR([K(F(c)) for c in coeffs])


def fpoly(coeffs):
    R = PolynomialRing(F, "v")
    return R([F(c) for c in coeffs])


# ------------------------------------------------------------------
# Local multiplicative-fibre resolution helper:
# for an I_n germ, split to a*b0=t^n and use
# E_i: ord(a)=i, ord(b0)=n-i.
# ------------------------------------------------------------------
Kr = FunctionField(F, "r")
r = Kr.gen()
LS = LaurentSeriesRing(Kr, "s", default_prec=PREC)
s = LS.gen()


def newton_sqrt(value, root0):
    parent = value.parent()
    root = parent(root0)
    for _ in range(10):
        root = (root + value/root)/parent(2)
    assert (root*root-value).valuation() >= PREC-5
    return root


def polynomial_series(parent, variable, coeffs):
    return parent(sum(
        parent.base_ring()(c)*variable**i
        for i,c in enumerate(coeffs)
    ))


def shifted_polynomial_series(poly, root):
    # poly(root+s)
    answer = LS(0)
    for degree, coefficient in enumerate(poly.list()):
        answer += Kr(F(coefficient))*(Kr(root)+s)**degree
    return answer


def local_critical_data(root):
    Aloc = shifted_polynomial_series(
        PolynomialRing(F,"v")(A.list()), root
    )
    Bloc = shifted_polynomial_series(
        PolynomialRing(F,"v")(B.list()), root
    )

    # Find the node x-coordinate directly in the special cubic.
    xring = PolynomialRing(F, "x")
    x = xring.gen()
    cubic = (
        x**3
        + F(Aloc[0])*x
        + F(Bloc[0])
    )
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    center0 = F(-common[0]/common[1])

    target = -Aloc/LS(3)
    assert center0**2 == F(target[0])
    center = newton_sqrt(target, Kr(center0))

    g0 = center**3 + Aloc*center + Bloc
    assert g0.valuation() == 7

    rho0_candidates = tuple(
        value for value in F
        if value*value == F((3*center)[0])
    )
    assert len(rho0_candidates) == 2
    rho0 = min(rho0_candidates, key=int)

    return Aloc, Bloc, center, g0/s**7, Kr(rho0)


def toric_point_i7(root, component):
    Aloc, Bloc, center, unit, rho0 = local_critical_data(root)

    aa = LS(r)*s**component
    b0 = s**(7-component)/LS(r)
    bb = unit*b0
    yy = (aa+bb)/2
    ww = (bb-aa)/2

    rho = LS(rho0)
    for _ in range(10):
        value = rho**3 - 3*center*rho - ww
        derivative = 3*rho**2 - 3*center
        rho = rho-value/derivative
    assert (rho**3-3*center*rho-ww).valuation() >= PREC-7

    xx = center + ww/rho
    return xx, yy


def section_local_series(row, root):
    # Substitute V=root+s into rational section coordinates.
    def value(num, den):
        numerator = LS(0)
        denominator = LS(0)
        for degree,c in enumerate(num):
            numerator += Kr(F(c))*(Kr(root)+s)**degree
        for degree,c in enumerate(den):
            denominator += Kr(F(c))*(Kr(root)+s)**degree
        return numerator/denominator

    return (
        value(row["X_num"], row["X_den"]),
        value(row["Y_num"], row["Y_den"]),
    )


def functional_rows(function_row):
    nonzero = [value for value in function_row if value]
    common = nonzero[0].denominator().parent().one()
    for value in nonzero:
        common = common.lcm(value.denominator())

    numerators = [
        (value*common).numerator()
        for value in function_row
    ]
    max_degree = max(
        [poly.degree() for poly in numerators if poly] + [-1]
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


def dedup_rows(rows):
    result = []
    seen = set()
    for row in rows:
        key = tuple(map(int,row))
        if key not in seen:
            seen.add(key)
            result.append(tuple(F(v) for v in row))
    return tuple(result)


# ------------------------------------------------------------------
# Kodaira classifier.
# ------------------------------------------------------------------
def kodaira_data(ord_a,ord_b,ord_delta):
    if ord_a == 0 or ord_b == 0:
        n=int(ord_delta)
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
        n=int(ord_delta-6)
        rank=n+4
        return rank,2*rank*(rank-1),4,n+6,f"I{n}*"
    if ord_delta == 8:
        return 6,72,3,8,"IV*"
    if ord_delta == 9:
        return 7,126,2,9,"III*"
    if ord_delta == 10:
        return 8,240,1,10,"II*"
    raise ArithmeticError((ord_a,ord_b,ord_delta))


def classify_quartic(quartic,twist):
    coeffs=list(quartic.list())+[K(0)]*5
    e,d,c,b,a=coeffs[:5]
    Iinv=12*a*e-3*b*d+c**2
    Jinv=72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3

    jacA=twist**2*(-27*Iinv)
    jacB=twist**3*(-27*Jinv)
    delta=twist**6*(4*Iinv**3-Jinv**2)

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
        raw=(
            int(jacA.valuation(factor)),
            int(jacB.valuation(factor)),
            int(delta.valuation(factor)),
        )
        scale=min(
            raw[0]//4,raw[1]//6,raw[2]//12
        )
        orders=(
            raw[0]-4*scale,
            raw[1]-6*scale,
            raw[2]-12*scale,
        )
        if orders[2] == 0:
            continue
        rank,count,determinant,local_euler,kind = (
            kodaira_data(*orders)
        )
        degree=int(factor.degree())
        root_rank += degree*rank
        root_count += degree*count
        root_det *= determinant**degree
        euler += degree*local_euler
        finite.append(
            (str(factor),degree,orders,kind)
        )

    raw_inf=tuple(
        int(value.denominator().degree()-value.numerator().degree())
        for value in (jacA,jacB,delta)
    )
    scale=min(
        raw_inf[0]//4,raw_inf[1]//6,raw_inf[2]//12
    )
    inf_orders=(
        raw_inf[0]-4*scale,
        raw_inf[1]-6*scale,
        raw_inf[2]-12*scale,
    )
    inf_kind="smooth"
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
            int(root_rank),
            int(root_count),
            int(root_det),
        ),
        "MW": int(18-root_rank),
        "euler": int(euler),
        "finite": tuple(finite),
        "infinity_orders": tuple(inf_orders),
        "infinity": inf_kind,
        "A": str(jacA),
        "B": str(jacB),
    }


# ------------------------------------------------------------------
# Process the four +/- pair representatives.
# ------------------------------------------------------------------
summary = []

for row in candidates:
    Xp = VF(kpoly(row["X_num"])) / VF(kpoly(row["X_den"]))
    Yp = VF(kpoly(row["Y_num"])) / VF(kpoly(row["Y_den"]))
    assert Yp**2 == Xp**3 + VF(A)*Xp + VF(B)

    # Denominator divisor h^2,h^3, with deg(h)=3=P.O.
    xden = PolynomialRing(F,"v")(row["X_den"])
    yden = PolynomialRing(F,"v")(row["Y_den"])
    h = xden.sqrt()
    assert h**2 == xden
    assert h.degree() == 3
    assert h.monic() == h
    assert yden == h**3
    assert gcd(h,h.derivative()) == 1

    Nx = PolynomialRing(F,"v")(row["X_num"])
    Ny = PolynomialRing(F,"v")(row["Y_num"])

    # 9D raw ambient:
    # Acoef/h^2 with deg A<=6; Bcoef/h*m with deg B<=1.
    ambient = tuple(
        [("A",i) for i in range(7)]
        + [("B",i) for i in range(2)]
    )
    assert len(ambient) == 9

    modulus = h**2
    smooth_columns = []
    v = h.parent().gen()

    # m=(y+Yp)/(x-Xp), hence p=-Yp/Xp at P.O collisions:
    # A*Nx - B*Ny == 0 mod h^2.
    for kind,exponent in ambient:
        if kind == "A":
            residue = (v**exponent*Nx).mod(modulus)
        else:
            residue = (-v**exponent*Ny).mod(modulus)
        smooth_columns.append(residue)

    smooth_matrix = matrix(
        F,
        6,
        9,
        lambda rr,cc: F(smooth_columns[cc][rr]),
    )
    assert smooth_matrix.rank() == 6
    assert smooth_matrix.right_kernel().dimension() == 3

    print(
        f"Q807774RR|candidate={row['name']}|mw={row['mw']}|"
        f"h={h}|smooth_rank=6|smooth_kernel=3|"
        "status=PASS_SMOOTH_PO3_MODULE",
        flush=True,
    )

    # Identify nonzero I7 and endpoint correction.
    labels_i7 = row["labels"][0]
    if labels_i7[0]:
        fibre_index = 0
        root = i7_roots[0]
        component = labels_i7[0]
    else:
        fibre_index = 1
        root = i7_roots[1]
        component = labels_i7[1]

    assert component in (2,5)
    endpoint = 1 if component == 2 else 6

    Xploc,Yploc = section_local_series(row,root)
    xloc,yloc = toric_point_i7(root,endpoint)
    mloc = (yloc+Yploc)/(xloc-Xploc)
    assert mloc.valuation() >= 0

    # Restriction of the 9 ambient generators to the endpoint component.
    Hloc = polynomial_series(
        LS,s,tuple(h.list())
    )
    # h(root+s), not h(s)
    Hloc = LS(0)
    for degree,c in enumerate(h.list()):
        Hloc += Kr(F(c))*(Kr(root)+s)**degree
    assert Hloc.valuation() == 0

    endpoint_functions = []
    for kind,exponent in ambient:
        basepower = (Kr(root)+s)**exponent
        if kind == "A":
            value = basepower/Hloc**2
        else:
            value = basepower/Hloc*mloc
        assert value.valuation() >= 0
        endpoint_functions.append(
            Kr(value[0])
            if value.valuation() == 0
            else Kr(0)
        )

    endpoint_rows = dedup_rows(
        functional_rows(tuple(endpoint_functions))
    )
    endpoint_matrix = matrix(F,endpoint_rows)
    assert endpoint_matrix.rank() == 1

    conditions = smooth_matrix.stack(endpoint_matrix)
    rank = conditions.rank()
    kernel = conditions.right_kernel().basis_matrix()

    assert rank == 7
    assert kernel.nrows() == 2 and kernel.ncols() == 9
    assert (
        conditions*kernel.transpose()
        == matrix(F,conditions.nrows(),2)
    )

    print(
        f"Q807774RR|candidate={row['name']}|"
        f"I7={fibre_index+1}|root={int(root)}|component={component}|"
        f"endpoint=Theta{endpoint}|"
        f"endpoint_row={tuple(tuple(map(int,rw)) for rw in endpoint_matrix.rows())}|"
        f"rank=7|kernel={tuple(tuple(map(int,rw)) for rw in kernel.rows())}|"
        "status=PASS_RESOLVED_9D_TO_2D_KERNEL",
        flush=True,
    )

    def kernel_function(krow):
        coeffs = tuple(map(K,krow))
        Apoly = sum(coeffs[i]*V**i for i in range(7))
        Bpoly = coeffs[7] + coeffs[8]*V
        return (
            VF(Apoly / VF(kpoly(tuple(h.list())))**2),
            VF(Bpoly / VF(kpoly(tuple(h.list())))),
        )

    (a0,b0),(a1,b1) = map(kernel_function,kernel.rows())
    determinant = a0*b1-a1*b0
    assert determinant != 0

    denominator = K(T)*b0-b1
    assert denominator != 0
    mvalue = (a1-K(T)*a0)/denominator

    cover = (
        mvalue**4
        - 6*Xp*mvalue**2
        - 8*Yp*mvalue
        - 3*Xp**2
        - 4*VF(A)
    )

    square_class = VR(
        cover.numerator()*cover.denominator()
    )
    factorization = square_class.factor()
    odd = VR.one()
    for factor,exponent in factorization:
        if int(exponent)%2:
            odd *= factor

    profile = tuple(
        (str(factor),int(factor.degree()),int(exponent))
        for factor,exponent in factorization
    )
    odd_degree = int(odd.degree())

    print(
        f"Q807774RR|candidate={row['name']}|"
        f"factor_profile={profile}|odd_degree={odd_degree}|"
        "status=PASS_KERNEL_COVER",
        flush=True,
    )

    classification = None
    if odd_degree == 4:
        quartic = odd.monic()
        twist = K(factorization.unit())
        classification = classify_quartic(quartic,twist)

        print(
            f"Q807774RRCHILD|candidate={row['name']}|mw={row['mw']}|"
            f"quartic={quartic}|twist={twist}|"
            f"root_data={classification['root_data']}|"
            f"MW={classification['MW']}|"
            f"finite={classification['finite']}|"
            f"infinity={classification['infinity']}|"
            f"euler={classification['euler']}|"
            "status=PASS_CHILD_CLASSIFICATION",
            flush=True,
        )

        assert classification["euler"] == 24

    summary.append((
        row["name"],
        row["mw"],
        labels_i7,
        int(root),
        component,
        endpoint,
        odd_degree,
        classification["root_data"] if classification else None,
        classification["MW"] if classification else None,
    ))

print(
    f"Q807774RR|summary={tuple(summary)}|"
    f"quartic_count={sum(1 for row in summary if row[6] == 4)}|"
    f"distinct_root_data={tuple(sorted(set(row[7] for row in summary if row[7] is not None)))}|"
    "status=PASS_FOUR_PAIR_RESOLVED_DEGREE_TWO_SCAN",
    flush=True,
)
"""


def main():
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    print(f"sage={sage}", flush=True)
    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_q6_7774_resolved_degree_two_scan.sage"
        script.write_text(SAGE_CODE)
        subprocess.run([sage, str(script)], check=True)


if __name__ == "__main__":
    main()
