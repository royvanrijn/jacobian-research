#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from itertools import product

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, matrix, vector,
    ceil, sqrt
)


F = GF(73)
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
E = EllipticCurve(K, [0,0,0,K(A),K(B)])

i7_roots = (F(-20), F(-67))
i2_roots = (F(-17), F(-30), F(-68))

xring = PolynomialRing(F, "x")
x = xring.gen()


def node_x(root):
    cubic = x**3 + A(root)*x + B(root)
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return -common[0]/common[1]


nodes = {root: node_x(root) for root in i7_roots+i2_roots}


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


def constrained_x_family(selected_roots):
    interpolation = R.lagrange_polynomial(
        [(root,nodes[root]) for root in selected_roots]
    )
    vanishing = R.one()
    for root in selected_roots:
        vanishing *= V-root
    return interpolation, vanishing


# Reproduce the old verifier's deterministic 30 polynomial sections exactly.
interpolation, vanishing = constrained_x_family(i7_roots)
assert vanishing.degree() == 2

two_node_candidates = []
for coefficients in product(F, repeat=3):
    X = interpolation + R(coefficients)*vanishing
    for Y in polynomial_square_roots(X**3+A*X+B):
        assert Y**2 == X**3+A*X+B
        two_node_candidates.append((X,Y))

assert len(two_node_candidates) == 30
points = tuple(
    E(K(X),K(Y))
    for X,Y in two_node_candidates
)
assert len(set(points)) == 30

# Pinned basis from search_q80_third_child_polynomial_sections_gf73.sage.
P1 = points[20]
P2 = points[0]
P3 = points[2]

basis_labels = (
    ((5,2),(1,0,1)),
    ((1,3),(1,1,0)),
    ((5,5),(0,0,0)),
)

G = matrix(QQ, [
    [QQ(1)/7, QQ(1)/14, 0],
    [QQ(1)/14, QQ(3)/7, -QQ(1)/7],
    [0, -QQ(1)/7, QQ(8)/7],
])
assert G.det() == QQ(3)/49

generic_height = QQ(219)/28
generic_correction = QQ(61)/28
assert generic_height == 10-generic_correction
assert generic_correction == QQ(10)/7 + QQ(3)/4


def labels_for_coeffs(a,b,c):
    coeffs = (a,b,c)
    i7 = tuple(
        sum(
            coeffs[k]*basis_labels[k][0][j]
            for k in range(3)
        ) % 7
        for j in range(2)
    )
    i2 = tuple(
        sum(
            coeffs[k]*basis_labels[k][1][j]
            for k in range(3)
        ) % 2
        for j in range(3)
    )
    return i7,i2


def correction(labels):
    i7,i2 = labels
    return (
        sum(QQ(j*(7-j))/7 for j in i7)
        + QQ(sum(i2))/2
    )


def section_pole_from_height(height, labels):
    value = (height + correction(labels) - 4)/2
    assert value in QQ
    return value


# Exact coordinate bounds from v_i^2 <= h * (G^-1)_ii.
Ginv = G.inverse()
bounds = tuple(
    int(ceil(sqrt(generic_height*Ginv[i,i])))
    for i in range(3)
)
assert bounds == (8,5,3)

# Tighten each by exact square test.
tight_bounds = []
for i,bound in enumerate(bounds):
    while bound > 0 and QQ(bound**2) > generic_height*Ginv[i,i]:
        bound -= 1
    tight_bounds.append(bound)
tight_bounds = tuple(tight_bounds)
assert tight_bounds == (7,4,2)

all_po3 = []
for a in range(-tight_bounds[0], tight_bounds[0]+1):
    for b in range(-tight_bounds[1], tight_bounds[1]+1):
        for c in range(-tight_bounds[2], tight_bounds[2]+1):
            if (a,b,c) == (0,0,0):
                continue
            coeff = vector(QQ,(a,b,c))
            height = QQ(coeff*G*coeff)
            if height > generic_height:
                continue
            labels = labels_for_coeffs(a,b,c)
            pole = section_pole_from_height(height,labels)
            if pole == 3:
                all_po3.append(
                    (a,b,c,height,labels,correction(labels))
                )

assert len(all_po3) == 126
best_height = max(row[3] for row in all_po3)
best = [row for row in all_po3 if row[3] == best_height]

assert best_height == QQ(53)/7
assert generic_height-best_height == QQ(1)/4
assert len(best) == 8

for row in best:
    _,_,_,height,labels,corr = row
    i7,i2 = labels
    assert sorted(i7) in ([0,2],[0,5])
    assert i2 == (1,0,1)
    assert corr == QQ(17)/7
    assert height == 10-corr

print(
    f"Q807774MW|bounds={tight_bounds}|"
    f"po3_under_generic_height={len(all_po3)}|"
    f"best_height={best_height}|height_drop={generic_height-best_height}|"
    f"best_count={len(best)}|"
    "status=PASS_EXACT_SPECIALIZED_MW_CENSUS",
    flush=True,
)


def rational_coefficients(value):
    return (
        tuple(map(int,value.numerator().list())),
        tuple(map(int,value.denominator().list())),
    )


def section_pole(point):
    assert not point.is_zero()
    xx = point[0]
    denominator_degree = xx.denominator().degree()
    numerator_degree = xx.numerator().degree()
    assert denominator_degree % 2 == 0
    infinity_excess = max(
        0, numerator_degree-denominator_degree-4
    )
    assert infinity_excess % 2 == 0
    return denominator_degree//2 + infinity_excess//2


# Direct node support sanity on I2 fibres.
def hits_node(point, root):
    if point.is_zero():
        return False
    xx,yy = point[0],point[1]
    if xx.denominator()(root) == 0:
        return False
    return xx(root) == nodes[root] and yy(root) == 0


actual = []
for index,row in enumerate(best,1):
    a,b,c,height,labels,corr = row
    point = a*P1+b*P2+c*P3
    assert not point.is_zero()
    assert section_pole(point) == 3

    observed_i2 = tuple(
        1 if hits_node(point,root) else 0
        for root in i2_roots
    )
    assert observed_i2 == labels[1]

    actual.append((row,point))

    print(
        f"Q807774MWCAND|index={index}|coeffs={(a,b,c)}|"
        f"height={height}|labels={labels}|corr={corr}|P.O=3|"
        f"X={rational_coefficients(point[0])}|"
        f"Y={rational_coefficients(point[1])}|"
        "status=PASS_CANDIDATE_SECTION",
        flush=True,
    )

# Pair candidates by negation in MW coordinates.
unpaired = {tuple(row[:3]) for row,_ in actual}
pairs = []
while unpaired:
    coeffs = min(unpaired)
    neg = tuple(-x for x in coeffs)
    assert neg in unpaired
    unpaired.remove(coeffs)
    unpaired.remove(neg)
    pairs.append((coeffs,neg))

assert len(pairs) == 4
print(
    f"Q807774MW|pm_pairs={tuple(pairs)}|"
    "status=PASS_FOUR_PM_PAIRS",
    flush=True,
)

# Local vertical possibilities at the nonzero I7.
# For label 2, the minimal nef square-correcting vector is -Theta1;
# for label 5 it is the reflected -Theta6.
for row,_ in actual:
    a,b,c,height,labels,corr = row
    i7,_ = labels
    nonzero_index = 0 if i7[0] else 1
    component = i7[nonzero_index]
    assert component in (2,5)
    endpoint = 1 if component == 2 else 6
    print(
        f"Q807774MWVERT|coeffs={(a,b,c)}|"
        f"nonzero_I7={nonzero_index+1}|component={component}|"
        f"minimal_nef_vertical=-Theta{endpoint}|L1=1|twist=0|"
        "status=PASS_MINIMAL_SPECIAL_VERTICAL_OPTION",
        flush=True,
    )

print(
    "Q807774MW|next=test_four_pm_pairs_with_resolved_degree_two_module|"
    "status=PASS_Q6_7774_SPECIALIZED_CANDIDATE_SET",
    flush=True,
)
"""


def main():
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    print(f"sage={sage}", flush=True)
    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_q6_7774_cm24_mw_candidates.sage"
        script.write_text(SAGE_CODE)
        subprocess.run([sage, str(script)], check=True)


if __name__ == "__main__":
    main()
