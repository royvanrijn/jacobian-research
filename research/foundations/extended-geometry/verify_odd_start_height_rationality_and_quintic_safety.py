#!/usr/bin/env python3
"""Exact support for odd-start rationality, zero-center and quintic safety.

The literal integer Bezout and Z[a] Sylvester identities prove their
stated algebraic eliminations. Raw endpoint comparisons are finite
controls; general rational descent and all-marker safety need the
separate written arguments.
"""
if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")
from fractions import Fraction as Q
from itertools import permutations
from math import gcd,lcm,comb,factorial
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
PROOF_NAME="ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md"


def unit_degree_box_packet():

    def trim(a):
        a = list(map(Q, a))
        while len(a) > 1 and (not a[-1]):
            a.pop()
        return a

    def add(a, b):
        return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b)))])

    def scale(a, c):
        return trim([c * x for x in a])

    def mul(a, b):
        out = [Q(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return trim(out)

    def divide(a, b):
        a = trim(a)
        b = trim(b)
        q = [Q(0)] * max(1, len(a) - len(b) + 1)
        while a != [0] and len(a) >= len(b):
            i = len(a) - len(b)
            v = a[-1] / b[-1]
            q[i] += v
            a = add(a, [Q(0)] * i + scale(b, -v))
        return (trim(q), a)

    def bezout(a, b):
        r0, r1 = (trim(a), trim(b))
        s0, s1 = ([Q(1)], [Q(0)])
        t0, t1 = ([Q(0)], [Q(1)])
        while r1 != [0]:
            q, r = divide(r0, r1)
            r0, r1 = (r1, r)
            s0, s1 = (s1, add(s0, scale(mul(q, s1), -1)))
            t0, t1 = (t1, add(t0, scale(mul(q, t1), -1)))
        assert len(r0) == 1 and r0[0]
        s = scale(s0, 1 / r0[0])
        t = scale(t0, 1 / r0[0])
        assert add(mul(a, s), mul(b, t)) == [Q(1)]
        denominator = lcm(*(v.denominator for v in s + t))
        si = [int(v * denominator) for v in s]
        ti = [int(v * denominator) for v in t]
        common = gcd(denominator, gcd(*si, *ti))
        denominator //= common
        si = [v // common for v in si]
        ti = [v // common for v in ti]
        assert add(mul(a, si), mul(b, ti)) == [Q(denominator)]
        quotient, remainder = divide(b, a)
        return {'first': list(a), 'second': list(b), 'Bezout_first_multiplier': si, 'Bezout_second_multiplier': ti, 'nonzero_constant': denominator, 'second_remainder_mod_first': [str(x) for x in remainder]}
    unit2 = [1, 22, 55]
    unit4 = [comb(4, 4 - j) * comb(21, j) for j in range(5)]
    assert unit4 == [1, 84, 1260, 5320, 5985]
    packet = bezout(unit2, unit4)
    assert packet['nonzero_constant'] == 7301194
    return {'constant_H_second_coefficient': 55, 'linear_H_Bezout': packet}

def zero_center_packet():

    def trim(a):
        a = list(a)
        while len(a) > 1 and (not a[-1]):
            a.pop()
        return a

    def add(a, b):
        return trim([(a[j] if j < len(a) else 0) + (b[j] if j < len(b) else 0) for j in range(max(len(a), len(b)))])

    def scale(a, c):
        return trim([c * x for x in a])

    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return trim(out)

    def product(items):
        out = [1]
        for p in items:
            out = mul(out, p)
        return out

    def build_packet():
        factors2 = [[1, 2], [3, 2], [5, 2]]
        factors4 = [[1 + 2 * j, 4] for j in range(5)]
        quadratic = [product(factors2[1:]), scale(product([factors2[0], factors2[2]]), 2), product(factors2[:2])]
        quartic = [scale(product((f for k, f in enumerate(factors4) if k != j)), c) for j, c in enumerate((1, 4, 6, 4, 1))]
        zero = [0]
        matrix = []
        for shift in range(4):
            matrix.append([zero] * shift + quadratic + [zero] * (3 - shift))
        for shift in range(2):
            matrix.append([zero] * shift + quartic + [zero] * (1 - shift))
        assert len(matrix) == 6 and all((len(row) == 6 for row in matrix))
        determinant = [0]
        nonzero_terms = 0
        for permutation in permutations(range(6)):
            inversions = sum((permutation[i] > permutation[j] for i in range(6) for j in range(i + 1, 6)))
            term = product((matrix[i][permutation[i]] for i in range(6)))
            if term != [0]:
                nonzero_terms += 1
            determinant = add(determinant, scale(term, (-1) ** inversions))
        P = [10125, 115200, 518064, 1195296, 1565824, 1209600, 544768, 132096, 13312]
        assert all((c > 0 for c in P))
        assert determinant == scale(P, 4096)
        return {'coefficient_convention': 'All arrays ascend in a; quadratic/quartic coefficient lists descend in z.', 'D2_factors': factors2, 'D4_factors': factors4, 'A2_coefficients_in_z': quadratic, 'A4_coefficients_in_z': quartic, 'Sylvester_matrix': matrix, 'determinant_coefficients': determinant, 'positive_polynomial_P': P, 'constant_factor': 4096, 'nonzero_determinant_permutation_terms': nonzero_terms, 'identity': 'det Sylvester(D2*P2,D4*P4)=4096P(a) in Z[a].', 'quantifier_scope': 'Polynomial coefficient identity, not finite sampling in a. Positive coefficients make P(a) nonzero for every positive integer a; characteristic zero preserves nonvanishing.', 'limits': 'Applies to zero-center phases t^a(C+t^2D)/x with positive odd integera and C(0)=D(0)=0. It does not classify arbitrary polynomials of heighta+2.'}
    return build_packet()


def endpoint_controls():
    def multiply(A,B):
        out={}
        for (x,t),c in A.items():
            for (y,s),d in B.items():out[x+y,t+s]=out.get((x+y,t+s),Q(0))+c*d
        return {key:value for key,value in out.items()if value}
    controls=[]
    for a in(1,3,5):
        for h in(0,1,2):
            for b in map(Q,(0,1,-1,Q(1,2))):
                base={(0,a):2*b,(1,a):b*b,(-1,a):Q(1),(-1,a+2):Q(-1)}
                phase=multiply(base,{(0,0):Q(1),(1,0):Q(h)})
                power={(0,0):Q(1)}
                for m in range(1,7):
                    power=multiply(power,phase)
                    if m%2:continue
                    raw=sum(c/Q(t+1)for(x,t),c in power.items()if x==0 and t%2==0)
                    kappa=Q(2**m*factorial(m),1)
                    for j in range(m+1):kappa/=a*m+2*j+1
                    coefficient=sum(Q(comb(m,i)*comb((a+2)*m+1,m-i))*h**i*b**(m-i)for i in range(m+1))
                    assert raw==kappa*coefficient
                    controls.append({'a':a,'H_linear_coefficient':h,'J_linear_coefficient':str(b),
                        'm':m,'direct_angular_moment':str(raw),'kappa':str(kappa),
                        'endpoint_coefficient':str(coefficient)})
    assert 1+22+55==78
    wrong=1+2*9+comb(9,2)
    assert wrong==55 and wrong!=78
    return {'cases':controls,'case_count':len(controls),
            'wrong_height_exponent_rejected':True,
            'scope':'Independent direct Laurent-x/height integral versus endpoint coefficient at odd starts1,3,5. Finite controls only.'}

def contraction_weight_controls():
    def contract(alpha,beta):
        if any(a>b for a,b in zip(alpha,beta)):return None
        coefficient=1
        for a,b in zip(alpha,beta):coefficient*=factorial(b)//factorial(b-a)
        return coefficient,tuple(b-a for a,b in zip(alpha,beta))
    # Nonzero negative input weight can survive: E2(z1)=z1.
    assert contract((0,0),(1,0))==(1,(1,0))
    rows=[]
    # Fixed marker xi1*z2^2 times powers of xi2*z1: output degree1.
    for m in range(1,5):
        alpha=(1,m);beta=(m,2);out=contract(alpha,beta)
        if out is not None:
            coefficient,exponents=out
            assert sum(exponents)==1 and -exponents[0]==1-m
            assert -1<=-exponents[0]<=0
        if m>=3:assert out is None
        rows.append({'m':m,'input_weight':1-m,'output':out})
    assert any(row['input_weight']<0 and row['output']is not None for row in rows)
    # Every positive input weight vanishes because no z-output has it.
    for a in range(4):
        for b in range(4):
            if a>b:assert contract((a,0),(b,3))is None
    return {'E2_z1_is_z1':True,'negative_branch_fixed_degree_controls':rows,
        'positive_weights_rejected':True,
        'scope':'Regression against the false assertion that all nonzero original torus weights vanish. Eventual safety uses the written output-degree argument.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'odd_start_height_rationality_and_quintic_safety_verification.json')
    args=parser.parse_args()
    result={'passed':True,'unit_center_degree_box':unit_degree_box_packet(),
        'all_degree_zero_center':zero_center_packet(),'endpoint_controls':endpoint_controls(),
        'contraction_weight_controls':contraction_weight_controls(),
        'source_sha256':{Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            PROOF_NAME:hashlib.sha256((HERE/PROOF_NAME).read_bytes()).hexdigest()},
        'scope':'Exact literal Bezout and universal Z[a] Sylvester coefficient identities, with finite raw endpoint controls. Written arguments supply all-degree rationality, purity equivalence, safety and infinity restriction.',
        'external_review':False,'formal_verification':False}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_endpoint_cases':result['endpoint_controls']['case_count'],
        'Sylvester_nonzero_terms':result['all_degree_zero_center']['nonzero_determinant_permutation_terms'],
        'unit_Bezout_constant':result['unit_center_degree_box']['linear_H_Bezout']['nonzero_constant']}))

if __name__=='__main__':main()
