#!/usr/bin/env python3
"""Exact controls for the joint 2,3,5,17 elliptic local sieve."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import platform


A2, A4 = -10656, 27091730
Q0 = 2709173
P = (F(373248,49),F(59459616,343))
S = (F(92179201,6400),F(553221611649,512000))
T = (F(0),F(0))


def add(left,right):
    if left is None:
        return right
    if right is None:
        return left
    x,y = left
    u,v = right
    if x == u:
        if y+v == 0:
            return None
        slope = (3*x*x+2*A2*x+A4)/(2*y)
    else:
        slope = (v-y)/(u-x)
    z = slope*slope-A2-x-u
    return z,-y+slope*(x-z)


def on_curve(point):
    if point is None:
        return True
    x,y = point
    return y*y == x*(x*x+A2*x+A4)


def valuation(value,prime):
    value = F(value)
    assert value
    a,b = value.numerator,value.denominator
    result = 0
    while a % prime == 0:
        a //= prime
        result += 1
    while b % prime == 0:
        b //= prime
        result -= 1
    return result


def residue(value,prime):
    value = F(value)
    return value.numerator*pow(value.denominator,-1,prime) % prime


def parameter(point):
    return -(point[0]-3552)/point[1]


def add_mod(left,right,prime):
    if left is None:
        return right
    if right is None:
        return left
    x,y = left
    u,v = right
    if x == u:
        if (y+v) % prime == 0:
            return None
        slope = (3*x*x+2*A2*x+A4)*pow(2*y,-1,prime) % prime
    else:
        slope = (v-y)*pow(u-x,-1,prime) % prime
    z = (slope*slope-A2-x-u) % prime
    return z,(-y+slope*(x-z)) % prime


def run():
    checks = []

    def check(name,assertion):
        if not assertion:
            raise AssertionError(name)
        checks.append(name)

    check('exact_known_generators_on_curve',all(on_curve(point) for point in (P,S,T)))
    R = add(P,T)
    check('exact_R_on_curve',on_curve(R))
    check('T_has_order_two',add(T,T) is None)
    check('T_addition_reciprocal_x_for_P_and_S',
          all(add(point,T)[0] == A4/point[0] for point in (P,S)))
    check('five_adic_coefficient_valuations',valuation(A4,5) == 1 and valuation(A2,5) == 0)

    affine = [(x,y) for x,y in product(range(5),repeat=2)
              if (y*y-x*(x*x+A2*x+A4)) % 5 == 0]
    singular = [(x,y) for x,y in affine
                if (2*y) % 5 == 0 and (3*x*x+2*A2*x+A4) % 5 == 0]
    smooth = [None]+[point for point in affine if point not in singular]
    check('complete_reduction_singular_point',singular == [(0,0)])
    check('complete_nonsingular_reduction_group',set(smooth) == {None,(1,0),(2,2),(2,3)})
    pbar = tuple(residue(value,5) for value in P)
    check('P_reduces_to_generator',pbar == (2,2))
    chain,current = [None],pbar
    while current is not None:
        assert current not in chain
        chain.append(current)
        current = add_mod(current,pbar,5)
    check('complete_order_four_chain',chain == [None,(2,2),(1,0),(2,3)])
    check('nonsingular_group_closure',all(add_mod(a,b,5) in smooth for a,b in product(smooth,repeat=2)))

    check('S_has_negative_x_valuation_minus_two',valuation(S[0],5) == -2)
    check('S_formal_five_depth_one',valuation(parameter(S),5) == 1)
    check('S_formal_five_normalized_unit_one',residue(parameter(S)/5,5) == 1)
    check('R_and_S_formal_two_depth_four',all(valuation(parameter(point),2) == 4 for point in (R,S)))
    # The E4 addition estimate makes t/16 modulo four additive.
    units = [residue(parameter(point)/16,4) for point in (R,S)]
    check('normalized_two_adic_coefficients_mod_four',units == [3,3])

    # Enumerate the exact two-primary quotient needed by the new
    # C4 factor. Residual formal root masses retain the higher levels.
    before = {a:F(0) for a in (1,2,5*Q0,10*Q0)}
    after = before.copy()
    formal_five = before.copy()
    rows = []
    for m,n,c in product(range(4),range(4),range(2)):
        u = (units[0]*m+units[1]*n) % 4
        a = (5*Q0 if m % 2 else 1)*(10*Q0 if c else 1)
        if a == 50*Q0*Q0:
            a = 2
        two = F(1) if c else (F(4,7) if u == 0 else F(0))
        three = {(0,0):F(1,13),(1,0):F(0),(0,1):F(1,3),(1,1):F(1)}[(c,n % 2)]
        seventeen = F(205,307) if (m+c) % 2 == 0 else F(1)
        weight = two*three*seventeen/32
        in_formal_five = m == 0 and c == 0
        five = F(1,31) if in_formal_five else F(1)
        before[a] += weight
        after[a] += weight*five
        if in_formal_five:
            formal_five[a] += weight
        rows.append({'m_mod4':m,'n_mod4':n,'c':c,'u_mod4':u,
                     'class':a,'formal_five':in_formal_five,
                     'before_weight':str(weight),'after_weight':str(weight*five)})
    check('thirty_two_joint_two_primary_classes',len(rows) == 32)
    check('old_four_field_coefficients_recovered',before == {
        1:F(205,111748),2:F(205,2456),5*Q0:F(1,84),10*Q0:F(1,8)})
    check('only_class_one_has_formal_five_mass',formal_five == {
        1:before[1]/2,2:F(0),5*Q0:F(0),10*Q0:F(0)})
    check('formal_five_geometric_mass',F(4,125)/(1-F(1,125)) == F(1,31))
    check('class_one_retention_factor',F(1,2)+F(1,2)*F(1,31) == F(16,31))
    check('new_exact_four_field_coefficients',after == {
        1:F(820,866047),2:F(205,2456),5*Q0:F(1,84),10*Q0:F(1,8)})
    total = sum(after.values())
    check('new_exact_total',total == F(2300089,10392564))
    check('strict_improvement_over_three_prime_sieve',
          sum(before.values())-total == F(3075,3464188) > 0)

    # Direct finite-module check of the Haar-coordinate change:
    # for either parity of m, each formal coordinate of the same
    # parity constraint has equal numbers of m=0 or2 (or1 or3) lifts.
    for ar,ass in product((1,3),repeat=2):
        for u in range(4):
            for parity in range(2):
                lifts = [(m,n) for m,n in product(range(4),repeat=2)
                         if m % 2 == parity and (ar*m+ass*n) % 4 == u]
                assert len(lifts) == 2
                assert {m for m,n in lifts} == {parity,parity+2}
    check('unit_coefficient_coordinate_change_full_finite_domain',True)

    base = Path(__file__).resolve().parent
    proof = base/'ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md'
    sources = [
        'ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md',
        'ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md',
        'ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md',
        'ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md',
        'sources/milne-elliptic-curves-2006.pdf',
    ]
    check('retained_inputs_exist',all((base/name).is_file() for name in sources))
    return {
        'status':'pass','checks_passed':len(checks),'exact_controls':checks,
        'proof_path':proof.name,'proof_sha256':sha256(proof.read_bytes()).hexdigest(),
        'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'command':'.venv/bin/python research/arithmetic/verify_elliptic_exception_joint_2_3_5_17.py --write-receipt',
        'python_version':platform.python_version(),
        'source_hashes':{name:sha256((base/name).read_bytes()).hexdigest() for name in sources},
        'finite_reduction_five':{'all_affine_points':affine,'singular_points':singular,
                                 'nonsingular_group':smooth,'generator_chain':chain},
        'finite_joint_domain':{'m_modulus':4,'n_modulus':4,'c_modulus':2,
                               'angular_unit_pairs':[[a,b] for a,b in product((1,3),repeat=2)],
                               'rows':rows},
        'old_coefficients':{str(a):str(value) for a,value in before.items()},
        'new_coefficients':{str(a):str(value) for a,value in after.items()},
        'new_total':str(total),
        'unrestricted_local_and_joint_proofs_not_finite_controls':True,
        'formal_verification':False,'external_review':False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.write_receipt:
        path = Path(__file__).with_name('elliptic_exception_joint_2_3_5_17_checks.json')
        path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(f"PASS: {receipt['checks_passed']} exact controls; coefficient {receipt['new_total']}")
