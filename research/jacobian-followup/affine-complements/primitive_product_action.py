"""PF120 exact product-action and table arithmetic; no classification algorithm."""
from itertools import product
from math import comb,gcd


def wreath_permutation(b,base,top):
    s=len(top)
    if b<3 or sorted(top)!=list(range(s)) or len(base)!=s:raise ValueError('Cartesian product action required')
    if any(sorted(g)!=list(range(b)) for g in base):raise ValueError('permutations required')
    points=list(product(range(b),repeat=s));index={x:i for i,x in enumerate(points)}
    return tuple(index[tuple(base[j][x[top[j]]] for j in range(s))] for x in points)


def classical_table(A,B,C,S):
    """Table2 arithmetic with row-local parameters; the proof/checker states their power ranges."""
    half=S.Rational(1,2)
    return {
      'L_field2_points':(2*A-1,A),
      'L_field3_points':(half*(3*A-1),A-1),
      'U4_field3':(S.Integer(112),S.Integer(72)),
      'Sp_field2_points':(2*A-1,A),
      'Sp_plus_forms':(A*(2*A+1),A*(A+1)),
      'Sp_minus_forms':(A*(2*A-1),A*(A-1)),
      'odd_O_field3_P1':(half*(9*B**2-1),B*(3*B-1)),
      'odd_O_field3_Nminus':(half*3*B*(3*B-1),3*B**2-2*B-1),
      'Oplus_field2_singular':((2*A-1)*(A+1),A*(A+1)),
      'Ominus_field2_singular':((2*A+1)*(A-1),A*(A-1)),
      'Ominus_field3_P1':(half*(3*C+1)*(C-1),C*(C-1)),
      'Oplus_field2_nonsingular':(A*(2*A-1),A*(A-1)),
      'Ominus_field2_nonsingular':(A*(2*A+1),A*(A+1)),
      'Oplus_field3_N1':(half*C*(3*C-1),C*(C-1)),
      'Ominus_field3_N1':(half*C*(3*C+1),C**2-1)}


def binomial_composite_certificate(m,k):
    if not 2<=k<=m-2:raise ValueError('interior binomial coefficient required')
    value=comb(m,k);factor=m//gcd(m,k)
    assert 1<factor<value and value%factor==0
    return dict(m=m,k=k,degree=value,factor=factor,cofactor=value//factor)
