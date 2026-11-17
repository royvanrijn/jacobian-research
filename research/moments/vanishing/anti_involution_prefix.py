"""A directed search chart for nonsingular first12-moment outsiders.

The linear part is z+w, so any solution has non-null H1 and is outside
Rodrigues. This is a candidate-search chart, not an exhaustion chart.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path

import sympy as S

az, aw, ez, ew, bz, bw = S.symbols('az aw ez ew bz bw')
PARAMETERS = (az, aw, ez, ew, bz, bw)


def odd_factorial(k):
    return math.prod(range(k, 0, -2))


@lru_cache(maxsize=None)
def height_powers(maximum):
    """Powers of 4+4*b*t+e*(5*t^2-1), as (b,e,t) exponent dictionaries."""
    factors = {(0, 0, 0): 4, (1, 0, 1): 4, (0, 1, 2): 5, (0, 1, 0): -1}
    powers = [{(0, 0, 0): 1}]
    for degree in range(maximum):
        nxt = defaultdict(int)
        for left, coefficient in powers[-1].items():
            for right, factor in factors.items():
                key = tuple(i+j for i, j in zip(left, right))
                nxt[key] += coefficient*factor
        powers.append({k: v for k, v in nxt.items() if v})
    return tuple(powers)


def exact_moment(order):
    if order % 2:
        return {}
    heights = height_powers(order)
    result = defaultdict(Fraction)
    for i in range(order+1):
        for j in range(order-i+1):
            k, ell = order//2-2*i+j, order//2+i-2*j
            if k < 0 or ell < 0:
                continue
            N = order//2+i+j
            multinomial = math.factorial(order)//(math.factorial(i)*math.factorial(j)*math.factorial(k)*math.factorial(ell))
            for (left_b, left_e, left_t), left in heights[k].items():
                for (right_b, right_e, right_t), right in heights[ell].items():
                    height = left_t+right_t
                    if height % 2:
                        continue
                    moment = Fraction(math.factorial(N)*odd_factorial(height-1), odd_factorial(2*N+height+1))
                    key = (i, j, left_e, right_e, left_b, right_b)
                    result[key] += multinomial*left*right*moment/Fraction(4**(k+ell))
    return {k: v for k, v in result.items() if v}


def primitive_moment(order):
    rational = exact_moment(order)
    denominator = math.lcm(*(value.denominator for value in rational.values()))
    integer = {key: int(value*denominator) for key, value in rational.items()}
    content = math.gcd(*(abs(value) for value in integer.values()))
    return {key: value//content for key, value in integer.items()}, Fraction(content, denominator)


def expression(terms, prime=0):
    pieces = []
    for exponents, raw in sorted(terms.items(), reverse=True):
        coefficient = raw % prime if prime else raw
        if coefficient == 0:
            continue
        factors = [str(variable) if exponent == 1 else str(variable)+'^'+str(exponent)
                   for variable, exponent in zip(PARAMETERS, exponents) if exponent]
        pieces.append(str(coefficient)+(' * '+'*'.join(factors) if factors else ''))
    return '+'.join(pieces).replace('+-', '-') or '0'


def phase():
    from conformal_euler_cubic import z, w, t
    return z+w+t*(bz*z+bw*w)+az*z**3+aw*w**3+ez*z*(t*t-z*w/2)+ew*w*(t*t-z*w/2)


def export(root, prime=41):
    root = Path(root)
    polynomials = {}
    scales = {}
    for n in range(2, 13, 2):
        terms, scale = primitive_moment(n)
        polynomials[n] = terms
        scales[n] = str(scale)
        print('ANTI MOMENT', n, 'TERMS', len(terms), flush=True)
    records = {str(n): [[list(exponents), coefficient] for exponents, coefficient in terms.items()]
               for n, terms in polynomials.items()}
    (root/'anti_involution_moments12.json').write_text(json.dumps({
        'variables': list(map(str, PARAMETERS)), 'moment_scales': scales,
        'primitive_moments': records,
        'phase': str(phase()),
        'cyclic_determinant': '(bz-bw)*(bz+bw)^2/8',
        'scope': 'Candidate-search chart only; no global exhaustion or outsider asserted.'}, separators=(',', ':'))+'\n')
    text = 'ring r='+str(prime)+',('+','.join(map(str, PARAMETERS))+'),dp;\noption(redSB);\n'
    text += '\n'.join('poly m'+str(n)+'='+expression(terms, prime)+';' for n, terms in polynomials.items())
    text += '\nideal I='+','.join('m'+str(n) for n in polynomials)+';\n'
    text += 'print("PREFIX12_BEGIN");\nideal G=std(I);\nprint("PREFIX12_DIM");print(dim(G));\n'
    text += 'print("PREFIX12_LENGTH");print(vdim(G));\nwrite("'+str(root/'anti_involution_groebner.sing')+'",G);\n'
    text += 'print("PREFIX12_END");\nquit;\n'
    path = root/('anti_involution_prefix12_p'+str(prime)+'.sing')
    path.write_text(text)
    return path
