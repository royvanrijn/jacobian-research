"""Factor-free Gauss coordinates followed by minimisation at explicitly listed primes.

No unrestricted hyperellminimalmodel call or complete discriminant factorization.
Exact polynomial identities certify the map; local minimality is not needed for
point/rank correctness. Different finite coordinate boxes require new exposure.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import pari
from search_observability import transform,multiply
from research_runtime.projective_box_change import classify
CAS=Path(__file__).resolve().parent
base=SourceFileLoader('bounded_prime_base_mapping',str(CAS/'factor_free_pari_mapping.sage')).load_module()
PRIMES=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]

def refine(original,primes):
    primes=list(primes)
    if primes!=PRIMES:raise ValueError('fixed25-prime refinement only')
    coeff=tuple(map(F,original['discriminant_quartic']))
    if any(c.denominator!=1 for c in coeff):raise ArithmeticError('integral input quartic required')
    polynomial='+'.join(f'({c})*x^{i}' for i,c in enumerate(coeff))
    vector='['+','.join(map(str,primes))+']'
    result=pari('my(m1,m2,C1,C2);C1=hyperellminimalmodel(['+polynomial+',0],&m1,'+vector+');C2=hyperellred(C1,&m2);[C2,m1,m2]')
    matrices=[tuple(F(str(m[i,j])) for i in range(2) for j in range(2)) for m in (result[1][1],result[2][1])]
    second=multiply(*matrices);combined=multiply(tuple(map(F,original['matrix'])),second)
    P=[F(str(result[0][0].polcoef(i))) for i in range(5)];Q=[F(str(result[0][1].polcoef(i))) for i in range(3)]
    disc=[4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)]
    raw=tuple(map(F,original['raw_coefficients']));transformed=transform(raw,combined);k=next(i for i,v in enumerate(disc) if v);ratio=transformed[k]/disc[k]
    if ratio<=0 or tuple(transformed)!=tuple(ratio*v for v in disc) or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator:raise ArithmeticError('exact restricted-prime map identity failed')
    return dict(centre=original['centre'],raw_coefficients=original['raw_coefficients'],denominator_clearing=original['denominator_clearing'],first_matrix=original['matrix'],second_matrix=list(map(str,second)),matrix=list(map(str,combined)),reduced_P=list(map(str,P)),reduced_Q=list(map(str,Q)),discriminant_quartic=list(map(str,disc)),square_ratio=str(ratio),coordinate_policy=dict(kind='raw',matrix=list(map(str,combined))),preparation_policy='integral-gauss-hyperellred-then-restricted25-prime-minimisation-v1',minimization_primes=primes,unrestricted_hyperellminimalmodel_called=False,comparison_to_factor_free=classify(original['matrix'],combined,125000),before_coefficient_bits=max(abs(int(c)).bit_length() for c in coeff),after_coefficient_bits=max(abs(int(c)).bit_length() for c in disc))

def mapping(model,points,centre):
    return refine(base.mapping(model,points,centre),PRIMES)
