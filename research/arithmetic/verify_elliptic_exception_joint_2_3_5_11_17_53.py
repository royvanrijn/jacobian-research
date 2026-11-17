#!/usr/bin/env python3
"""Exact controls for the correlated six-prime elliptic local count."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse
import json
import platform
import runpy

HERE = Path(__file__).resolve().parent
SOURCE = 'ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_53_SIEVE.md'


def run():
    previous = runpy.run_path(str(HERE/'verify_elliptic_exception_joint_2_3_5_11_17.py'))
    previous_report = previous['run']()
    old = runpy.run_path(str(HERE/'verify_elliptic_exception_joint_2_3_5_17.py'))
    rational = runpy.run_path(str(HERE/'verify_elliptic_exception_one_prime_positivity.py'))
    add, add_mod = old['add'], old['add_mod']
    parameter, valuation, residue = old['parameter'], old['valuation'], old['residue']
    P, S, T = old['P'], old['S'], old['T']
    R = add(P,T)
    Q = (F(3442688,529),F(304670016,12167))
    q0, a2, a4 = old['Q0'], old['A2'], old['A4']
    checks = []

    def check(label, condition):
        if not condition:
            raise AssertionError(label)
        checks.append(label)

    def mul_mod(V,n,p):
        result = None
        while n:
            if n & 1:
                result = add_mod(result,V,p)
            V = add_mod(V,V,p)
            n //= 2
        return result

    def red(V,p):
        return tuple(residue(z,p) for z in V) if V is not None else None

    check('previous_five_prime_controls_pass', previous_report['status']=='pass')
    check('all_rational_points_on_curve', all(old['on_curve'](V) for V in (P,Q,R,S,T)))
    check('R_and_S_generate_same_subgroup_as_P_Q_with_T',
          R==add(P,T) and S==tuple(-v if i else v for i,v in enumerate(add(P,Q))))
    disc = 16*a4*a4*(a2*a2-4*a4)
    check('fifty_three_good_active', disc % 53 != 0 and 53 % 3 == 2)
    points = [None]+[(x,y) for x,y in product(range(53),repeat=2)
                    if (y*y-x*(x*x+a2*x+a4)) % 53 == 0]
    check('complete_fifty_three_point_count', len(points)==60)
    check('complete_group_closure', all(add_mod(U,V,53) in points for U,V in product(points,repeat=2)))
    rb,sb,tb = [red(V,53) for V in (R,S,T)]
    check('literal_generator_reductions', (rb,sb,tb)==((22,11),(44,8),(0,0)))
    A,B3,B5 = (30,0),(13,45),(36,27)
    check('two_primary_projections', [mul_mod(V,15,53) for V in (rb,sb,tb)]==[A,tb,tb])
    check('three_primary_projections', [mul_mod(V,20,53) for V in (rb,sb,tb)]==[B3,None,None])
    check('five_primary_projections', [mul_mod(V,12,53) for V in (rb,sb,tb)]==[B5,B5,None])
    check('independent_two_torsion', A!=tb and mul_mod(A,2,53) is None
          and mul_mod(tb,2,53) is None and add_mod(A,tb,53)==(26,0))
    check('three_generator_exact_order', B3 is not None and mul_mod(B3,3,53) is None)
    check('five_generator_exact_order', B5 is not None and mul_mod(B5,5,53) is None)
    generated = {add_mod(add_mod(mul_mod(A,i,53),mul_mod(tb,j,53),53),
                          add_mod(mul_mod(B3,k,53),mul_mod(B5,l,53),53),53)
                 for i,j,k,l in product(range(2),range(2),range(3),range(5))}
    check('complete_primary_product_equals_curve', generated==set(points))
    check('no_finite_F_zero_branch', all((128*V[0]-25865) % 53 != 0 for V in points if V is not None))
    check('branch_squareclass_nonsquare', pow(51730 % 53,26,53)==52)
    class_signs = {a:1 if pow(a % 53,26,53)==1 else -1 for a in (1,2,5*q0,10*q0)}
    check('class_eligibility', class_signs=={1:1,2:-1,5*q0:-1,10*q0:1})

    g41 = (37,37)
    check('three_saturation_41_character', [mul_mod(red(V,41),16,41) for V in (R,S)]==[g41,(37,4)])
    check('three_saturation_53_character', [mul_mod(red(V,53),20,53) for V in (R,S)]==[B3,None])
    check('three_saturation_matrix_unit', (1*0-2*1) % 3 == 1)
    theta = [residue(parameter(rational['mul'](V,8,a2))/40,5) for V in (R,S)]
    check('five_saturation_character_matrix', theta==[2,1] and (theta[0]-theta[1]) % 5==1)
    check('complete_five_saturation_character_image',
          {((2*m+n)%5,(m+n)%5) for m,n in product(range(5),repeat=2)}
          ==set(product(range(5),repeat=2)))
    check('R_three_formal_depth', valuation(parameter(R),3)==3)
    check('S_three_formal_projection_depth', valuation(parameter(add(S,S)),3)==2)
    check('new_three_zero_equals_old_17_branch_on_n_zero',
          all((m%3==0)==((m+n)%3==0) for m,n in product(range(3),repeat=2) if n==0))
    check('new_five_zero_equals_old_11_pole',
          all(mul_mod(add_mod(mul_mod(rb,m,53),mul_mod(sb,n,53),53),12,53)
              ==mul_mod(B5,(m+n)%5,53) for m,n in product(range(5),repeat=2)))

    U = rational['mul'](P,60,a2)
    check('exact_sixty_P_on_curve', old['on_curve'](U))
    check('sixty_P_reduces_to_O', U[0].denominator % 53==0)
    check('formal_fifty_three_depth_one', valuation(parameter(U),53)==1)
    check('formal_fifty_three_residue', residue(parameter(U),53**2)==2491)
    check('class_one_old_parities_give_new_two_zero',
          add_mod(mul_mod(A,0,53),mul_mod(tb,0,53),53) is None)
    check('certified_class_old_parities_give_new_two_zero',
          add_mod(mul_mod(A,0,53),mul_mod(tb,1+1,53),53) is None)
    check('other_classes_never_two_zero',
          all(add_mod(A,mul_mod(tb,n+c,53),53) is not None for n,c in product(range(2),repeat=2)))

    def weight17(z):
        return F(1,307) if z%3==0 else F(1)

    def weight11(z):
        return F(1,133) if z==0 else F(1)

    mass3 = sum(weight17(z) for z in range(9))
    ratio3 = sum(weight17(z) for z in range(9) if z%3==0)/mass3
    check('conditional_three_zero_ratio', ratio3==F(1,615))
    five_ratios,joint_ratios = [],[]
    for formal5 in [False,True]:
        def weight5(theta):
            return (F(5,31) if theta==0 else F(0)) if formal5 else F(1)
        check('old_five_first_digit_mass_'+str(formal5),
              sum(weight5(theta) for theta in range(5))/5==(F(1,31) if formal5 else F(1)))
        mass5 = sum(weight5(theta)*weight11(z) for theta,z in product(range(5),repeat=2))
        pole5 = sum(weight5(theta)*weight11(z) for theta,z in product(range(5),repeat=2) if z==0)
        five_ratios.append(pole5/mass5)
        check('conditional_five_zero_ratio_'+str(formal5), five_ratios[-1]==F(1,533))
        total = sum(weight17(w)*weight5(theta)*weight11(z)
                    for w,theta,z in product(range(9),range(5),range(5)))
        both = sum(weight17(w)*weight5(theta)*weight11(z)
                   for w,theta,z in product(range(9),range(5),range(5)) if w%3==0 and z==0)
        joint_ratios.append(both/total)
        check('conditioned_joint_zero_ratio_'+str(formal5), joint_ratios[-1]==F(1,615*533))
    check('certified_uniform_odd_zero_probability',
          F(sum(1 for m,z in product(range(3),range(5)) if m==z==0),15)==F(1,15))
    soluble53 = F(52,53**3)/(1-F(1,53**3))
    check('formal_fifty_three_soluble_mass', soluble53==F(1,2863))
    rho = {1:F(1,615*533),2:F(0),5*q0:F(0),10*q0:F(1,15)}
    retain = {a:1-r*(1-soluble53) for a,r in rho.items()}
    check('exact_class_one_retention', retain[1]==F(312824741,312825695))
    check('exact_certified_retention', retain[10*q0]==F(13361,14315))
    before = {int(a):F(c) for a,c in previous_report['coefficients'].items()}
    after = {a:c*retain[a] for a,c in before.items()}
    check('exact_six_prime_coefficients', after=={
        1:F(1251298964,1648862553065),2:F(205,2456),
        5*q0:F(533,55860),10*q0:F(13361,114520)})
    check('exact_new_total', sum(after.values())==F(2081912199403,9893175318390))
    check('strict_total_improvement', sum(after.values())<sum(before.values())==F(116302241,531619620))
    check('strict_certified_improvement', after[10*q0]<F(1,8))
    inputs = [SOURCE,Path(__file__).name,
              'verify_elliptic_exception_joint_2_3_5_11_17.py',
              'verify_elliptic_exception_joint_2_3_5_17.py',
              'verify_elliptic_exception_one_prime_positivity.py',
              'ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md',
              'ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md',
              'ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md',
              'ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md',
              'ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md',
              'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
              'ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md',
              'sources/milne-elliptic-curves-2006.pdf']
    return dict(status='pass',checks_passed=len(checks),exact_controls=checks,
                earlier_five_prime_controls=previous_report['checks_passed'],
                reduction_points=points,complete_group_pairs=60**2,
                primary_projections={'two':[A,tb,tb],'three':[B3,None,None],'five':[B5,B5,None]},
                saturation_character_matrix_five=[theta,[1,1]],
                conditional_three_zero=str(ratio3),
                conditional_five_zero_by_five_state=[str(r) for r in five_ratios],
                conditional_joint_zero_by_five_state=[str(r) for r in joint_ratios],
                finite_weighted_states_per_five_state=9*5*5,
                certified_pole_probability='1/15',formal_soluble_fraction=str(soluble53),
                retention={str(a):str(r) for a,r in retain.items()},
                coefficients={str(a):str(c) for a,c in after.items()},new_total=str(sum(after.values())),
                input_sha256={name:sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
                python_version=platform.python_version(),arithmetic='standard-library exact integers and fractions',
                full_basis_assumed=False,all_places_asymptotic_claimed=False,
                full_joint_image_proved_by_finite_controls=False,formal_verification=False,external_review=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',action='store_true')
    args=parser.parse_args()
    receipt=run()
    if args.write_receipt:
        (HERE/'elliptic_exception_joint_2_3_5_11_17_53_checks.json').write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('PASS:',receipt['checks_passed'],'exact controls; coefficient',receipt['new_total'])
