"""Exact, bounded-cost interfaces for the written all-degree PF126 proof."""
from fractions import Fraction
import sympy as S


def gauss_frame(x, y, t):
    g = S.Poly(S.gcd(S.diff(x, t), S.diff(y, t)), t).monic().as_expr()
    gx, gy = S.cancel(S.diff(x, t)/g), S.cancel(S.diff(y, t)/g)
    a, b, unit = S.gcdex(gx, gy, t)
    if unit != 1:
        raise ValueError("not a polynomial Gauss frame")
    return g, (gx, gy), (-b, a)


def framing_data(multiplicities, unordered_contacts):
    if any(m < 1 for m in multiplicities) or any(k < 1 for k in unordered_contacts):
        raise ValueError("positive local multiplicities required")
    defect = sum(m-1 for m in multiplicities)
    delta = sum(m*(m-1)//2 for m in multiplicities)+sum(unordered_contacts)
    return dict(defect=defect, delta=delta, linking=2*delta+defect,
                preferred_critical_shear=(defect-(2*delta+defect)) % 2)


def binary_rank(columns):
    pivots = {}
    for column in columns:
        while column:
            k = column.bit_length()-1
            if k in pivots:
                column ^= pivots[k]
            else:
                pivots[k] = column
                break
    return len(pivots)


def residual_data(a, columns):
    """Columns label both preimages of each odd-contact double value."""
    if a < 2 or len(columns) % 2 or any(not 0 < u < 2**(a-1) for u in columns):
        raise ValueError("even number of nonzero hyperplane vectors required")
    rank = binary_rank(columns)
    infinity = 0
    for u in columns:
        infinity ^= u
    eps, w, b = int(infinity != 0), 2**rank, 2**a
    genus = 1-w+Fraction(w*(len(columns)+eps), 4)
    if genus.denominator != 1 or genus < 0:
        raise ValueError("Riemann-Hurwitz forbids these labels")
    return dict(dimension=a, odd_contact_values=len(columns)//2, rank=rank,
                component_degree=w, component_count=(b//2)//w,
                infinity_translation=infinity, infinity_inertia=2**eps,
                genus=int(genus), infinity_points=w//(2**eps),
                residual_longitude_cycles={2: b//4} if eps else {1: b//2},
                critical_tori=b//4, total_tori=b//2 if eps else 3*b//4,
                all_shears_zero=True)
