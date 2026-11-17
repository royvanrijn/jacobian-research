#!/usr/bin/env python3
"""Exact arithmetic controls for the joint 2,3,5,11,17 elliptic sieve."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import platform
import runpy

HERE = Path(__file__).resolve().parent
SOURCE = 'ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md'


def run():
    old = runpy.run_path(str(HERE/'verify_elliptic_exception_joint_2_3_5_17.py'))
    rational = runpy.run_path(str(HERE/'verify_elliptic_exception_one_prime_positivity.py'))
    add, add_mod = old['add'], old['add_mod']
    parameter, valuation, residue = old['parameter'], old['valuation'], old['residue']
    P, S, T = old['P'], old['S'], old['T']
    R = add(P, T)
    q0, a2, a4 = old['Q0'], old['A2'], old['A4']
    checks = []

    def check(name, assertion):
        if not assertion:
            raise AssertionError(name)
        checks.append(name)

    check('all_exact_points_on_curve', all(old['on_curve'](V) for V in (P, R, S, T)))
    disc = 16*a4*a4*(a2*a2-4*a4)
    check('eleven_is_good_and_inert_mod_three', disc % 11 != 0 and 11 % 3 == 2)
    points = [None]+[(x,y) for x,y in product(range(11), repeat=2)
                    if (y*y-x*(x*x+a2*x+a4)) % 11 == 0]
    check('complete_reduction_has_ten_points', len(points) == 10)
    rb, sb = [tuple(residue(v,11) for v in V) for V in (R,S)]
    check('R_and_S_same_reduction_generator', rb == sb == (9,5))
    chain, current = [], None
    for k in range(5):
        chain.append(current)
        current = add_mod(current, rb, 11)
    check('exact_order_five_chain', chain == [None,(9,5),(4,9),(4,2),(9,6)] and current is None)
    tb = (0,0)
    check('torsion_translation_fills_group', set(chain) | {add_mod(V,tb,11) for V in chain} == set(points))
    check('complete_group_closure', all(add_mod(U,V,11) in points for U,V in product(points,repeat=2)))

    def chi(x):
        return 1 if pow(x % 11,5,11) == 1 else -1

    def character(V):
        return 1 if V is None else chi(a4 if V == tb else V[0])

    check('complete_local_character_kernel', {V for V in points if character(V) == 1} == set(chain))
    check('complete_character_homomorphism', all(character(add_mod(U,V,11)) == character(U)*character(V)
                                              for U,V in product(points,repeat=2)))
    class_signs = {1:chi(1),2:chi(2),5*q0:chi(5*q0),10*q0:chi(10*q0)}
    check('four_global_class_signs', class_signs == {1:1,2:-1,5*q0:1,10*q0:-1})
    check('zero_branch_squareclass_nonsquare', chi(51730) == -1)
    check('no_finite_zero_branch', all((128*V[0]-25865) % 11 != 0 for V in points if V is not None))

    theta = []
    for name,V,want in [('R',R,2),('S',S,1)]:
        U = rational['mul'](V,8,a2)
        check('exact_eightfold_on_curve_'+name, old['on_curve'](U))
        check('five_projection_depth_one_'+name, valuation(parameter(U),5) == 1)
        theta.append(residue(parameter(U)/40,5))
        check('five_projection_first_digit_'+name, theta[-1] == want)
    check('first_digit_determinant_one', (theta[0]-theta[1]) % 5 == 1)
    image = {((theta[0]*m+theta[1]*n) % 5,(m+n) % 5) for m,n in product(range(5),repeat=2)}
    check('complete_joint_first_digit_image', image == set(product(range(5),repeat=2)))
    U = rational['mul'](P,10,a2)
    check('exact_tenfold_on_curve', old['on_curve'](U))
    check('tenfold_P_reduces_to_O', U[0].denominator % 11 == 0)
    check('eleven_formal_depth_one', valuation(parameter(U),11) == 1)
    check('eleven_first_digit_ten', residue(parameter(U),121) == 110)

    soluble = F(10,11**3)/(1-F(1,11**3))
    check('exact_formal_annulus_mass', soluble == F(1,133))
    retain = F(4,5)+soluble/5
    check('class_dependent_retention_factor', retain == F(533,665))
    before = {1:F(820,866047),2:F(205,2456),5*q0:F(1,84),10*q0:F(1,8)}
    after = {a:c*(retain if class_signs[a] == 1 else 1) for a,c in before.items()}
    check('all_four_exact_new_coefficients', after == {
        1:F(6724,8860327),2:F(205,2456),5*q0:F(533,55860),10*q0:F(1,8)})
    check('new_total', sum(after.values()) == F(116302241,531619620))
    check('strict_improvement', sum(after.values()) < sum(before.values()) == F(2300089,10392564))
    inputs = [SOURCE,'verify_elliptic_exception_joint_2_3_5_17.py',
              'verify_elliptic_exception_one_prime_positivity.py',
              'ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md',
              'ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md',
              'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
              'ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md',
              'sources/milne-elliptic-curves-2006.pdf',Path(__file__).name]
    return dict(status='pass',checks_passed=len(checks),exact_controls=checks,
                command='.venv/bin/python research/arithmetic/verify_elliptic_exception_joint_2_3_5_11_17.py --write-receipt',
                input_sha256={s:sha256((HERE/s).read_bytes()).hexdigest() for s in inputs},
                reduction_points=points,order_five_chain=chain,
                first_digit_matrix=[theta,[1,1]],class_signs=class_signs,
                coefficients={str(a):str(c) for a,c in after.items()},
                new_total=str(sum(after.values())),python_version=platform.python_version(),
                full_joint_image_and_counting_proofs_not_finite_controls=True,
                formal_verification=False,external_review=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.write_receipt:
        (HERE/'elliptic_exception_joint_2_3_5_11_17_checks.json').write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(f"PASS: {receipt['checks_passed']} exact controls; coefficient {receipt['new_total']}")
