#!/usr/bin/env python3
"""Exact finite controls for FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md.

The finite cases here are those explicitly reduced to residue computations,
plus bounded controls of the uniform Newton-hull inequalities. They do not
prove an unrestricted integer-target nonexistence claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []

    def require(name, condition, **details):
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True, **details})

    def equal(name, left, right=0):
        difference = sp.factor(left-right)
        if difference != 0:
            raise AssertionError(f"{name}: {difference}")
        checks.append({"name": name, "passed": True})

    W, pi, b, c = sp.symbols("W pi b c")
    beta = -sp.Rational(5, 4)/pi-sp.Rational(9, 32)*b/pi**2
    P = W**5-sp.Rational(3, 2)*W**4+sp.Rational(3, 2)*W**3/pi**2+beta*W**2+sp.Rational(9, 16)*W/pi-sp.Rational(9, 32)*c
    Hmod9 = W**3-sp.Rational(3, 2)*W**2+sp.Rational(3, 2)*W/pi**2-sp.Rational(5, 4)/pi
    equal("mod_nine_factor_difference", P-W**2*Hmod9,
          -sp.Rational(9, 32)*b*W**2/pi**2+sp.Rational(9, 16)*W/pi-sp.Rational(9, 32)*c)
    shifted = sp.Poly(sp.expand(Hmod9.subs(W, W-pi)), W)
    equal("shifted_linear_coefficient", shifted.nth(1), 3*(pi**2+pi+1/(2*pi**2)))
    equal("shifted_constant", shifted.nth(0), -pi**3-sp.Rational(3, 2)*pi**2-sp.Rational(11, 4)/pi)
    aa, bb, cc = sp.symbols("aa bb cc")
    generic_discriminant = sp.discriminant(W**3+aa*W**2+bb*W+cc, W)
    equal("cubic_discriminant_formula", generic_discriminant,
          aa**2*bb**2-4*bb**3-4*aa**3*cc-27*cc**2+18*aa*bb*cc)
    rem = sp.Poly(sp.expand((generic_discriminant.subs({aa:3*aa,bb:3*bb,cc:3*cc})+108*bb**3)/81), aa, bb, cc)
    require("other_discriminant_terms_divisible_eighty_one", all(coef.q == 1 for coef in rem.coeffs()))

    def fraction_mod(value, modulus):
        numerator, denominator = sp.fraction(sp.cancel(value))
        return (int(numerator)*pow(int(denominator), -1, modulus)) % modulus

    unit_classes = [1, 2, 4, 5, 7, 8]
    constant_table = [fraction_mod(shifted.nth(0).subs(pi,p), 9) for p in unit_classes]
    require("all_unit_constant_classes", constant_table == [6,6,3,3,0,0],
            unit_classes=unit_classes, constant_residues=constant_table)
    require("shifted_linear_normalized_units",
            all(fraction_mod((shifted.nth(1)/3).subs(pi,p),3) == p%3 for p in unit_classes))
    require("shifted_cubic_discriminant_normalized_units",
            all((-4*(p%3)**3)%3 == (-p)%3 for p in unit_classes))

    def convolution(first, second):
        out = [0]*(len(first)+len(second)-1)
        for i, x in enumerate(first):
            for j, y in enumerate(second):
                out[i+j] += x*y
        return out

    def hensel_factors_mod(p, target_b, target_c, exponent=9):
        modulus = 3**exponent
        target = [(-9*target_c*pow(32,-1,modulus))%modulus,
                  (9*pow(16*p,-1,modulus))%modulus,
                  (-5*pow(4*p,-1,modulus)-9*target_b*pow(32*p*p,-1,modulus))%modulus,
                  (3*pow(2*p*p,-1,modulus))%modulus,
                  (-3*pow(2,-1,modulus))%modulus, 1]
        q = [0,0,1]
        h = [p%3,0,0,1]
        for kk in range(1,exponent):
            power = 3**kk
            prod = convolution(q,h)
            residual = [((target[ii]-prod[ii])//power)%3 for ii in range(5)]
            dq0 = residual[0]*pow(p,-1,3)%3
            dq1 = residual[1]*pow(p,-1,3)%3
            dh0 = residual[2]
            dh1 = (residual[3]-dq0)%3
            dh2 = (residual[4]-dq1)%3
            q[0] += power*dq0
            q[1] += power*dq1
            h[0] += power*dh0
            h[1] += power*dh1
            h[2] += power*dh2
        if any((x-y)%modulus for x,y in zip(convolution(q,h),target)):
            raise AssertionError("Hensel coefficient lift")
        return q,h,modulus

    def valuation(number, prime):
        if number == 0:
            return None
        number = abs(number)
        exponent = 0
        while number%prime == 0:
            exponent += 1
            number //= prime
        return exponent

    def cubic_disc_low(coeff):
        d, cc, bb, aa = coeff
        return cc**2*bb**2-4*cc**3*aa-4*d*bb**3-27*d**2*aa**2+18*d*cc*bb*aa

    lifted_cases = 0
    for p in unit_classes:
        for target_b in range(3):
            for target_c in range(27):
                q,h,modulus = hensel_factors_mod(p,target_b,target_c)
                require_q = q[0]%9 == 0 and q[1]%9 == 0 and valuation(q[1],3) == 2
                if not require_q:
                    raise AssertionError((p,target_b,target_c,"Q coefficients"))
                dh = cubic_disc_low(h)%modulus
                if valuation(dh,3) != 3 or (dh//27)%3 != (-p)%3:
                    raise AssertionError((p,target_b,target_c,"H discriminant",dh))
                dq = (q[1]**2-4*q[0])%modulus
                product = dq*dh%modulus
                vc = valuation(target_c,3)
                if vc is None or vc >= 3:
                    expected_v = 7
                elif vc == 0:
                    expected_v = 5
                elif vc == 1:
                    expected_v = 6
                    if (product//3**6)%3 != (target_c//3)%3:
                        raise AssertionError("normalized c valuation one unit")
                elif (target_c//9*p)%3 != 1:
                    expected_v = 7
                else:
                    expected_v = None
                if expected_v is not None and valuation(product,3) != expected_v:
                    raise AssertionError((p,target_b,target_c,"full discriminant valuation",product,expected_v))
                lifted_cases += 1
    require("finite_Hensel_discriminant_controls", lifted_cases == 486, residue_cases=lifted_cases,
            modulus=3**9)

    # A separate hull implementation checks bounded instances of the
    # unrestricted supporting-line proof of ILO3, including a missing S^2 term.
    def lower_hull(points):
        hull=[]
        for point in points:
            while len(hull)>=2:
                first,second=hull[-2:]
                cross=(second[0]-first[0])*(point[1]-second[1])-(second[1]-first[1])*(point[0]-second[0])
                if cross<=0:
                    hull.pop()
                else:
                    break
            hull.append(point)
        return hull

    hull_cases=0
    for e in range(1,13):
        for t in range(13):
            ks=[1] if e==1 else list(range(min(e,2),20))+[None]
            for kval in ks:
                points=[(0,2+t),(1,2),(3,1+e),(4,1+4*e),(5,5*e)]
                if kval is not None:
                    points.append((2,kval))
                hull=lower_hull(sorted(points))
                edges=list(zip(hull,hull[1:]))
                if ((3,1+e),(5,5*e)) not in edges:
                    raise AssertionError((e,t,kval,"right edge",hull))
                if not any(y[0]-x[0]==1 for x,y in edges):
                    raise AssertionError((e,t,kval,"adjacent edge",hull))
                hull_cases+=1
    require("divisible_pi_Newton_hull_controls", hull_cases>2500, hull_cases=hull_cases)
    require("supporting_line_index_two", all(min(e,2)>sp.Rational(3,2)-e for e in range(1,40)))
    require("supporting_line_index_one", all(2>2-3*e for e in range(1,40)))
    require("supporting_line_index_zero", all(2>sp.Rational(5,2)-5*e for e in range(1,40)))

    # Two-adic Newton factors and simple-root residue controls from ILO4/5.
    from fractions import Fraction
    two_hull_cases = 0
    for e in range(4):
        for h in list(range(13))+[None]:
            points=[(1,1),(2,0),(3,4+e),(4,4+4*e),(5,5+5*e)]
            if h is not None:
                points.append((0,h))
            edges=list(zip(lower_hull(sorted(points)),lower_hull(sorted(points))[1:]))
            if ((2,0),(5,5+5*e)) not in edges:
                raise AssertionError((e,h,'odd-b right edge',edges))
            if e in (0,1,3) and Fraction(5+5*e,3).denominator != 3:
                raise AssertionError('odd-b reduced cubic slope')
            two_hull_cases += 1
    for e,kval in [(0,1),(1,2),(2,1)]:
        for h in list(range(13))+[None]:
            points=[(1,1),(2,kval),(3,4+e),(4,4+4*e),(5,5+5*e)]
            if h is not None:
                points.append((0,h))
            hull=lower_hull(sorted(points))
            if ((2,kval),(5,5+5*e)) not in list(zip(hull,hull[1:])):
                raise AssertionError((e,kval,h,'even-b cubic edge',hull))
            if Fraction(5+5*e-kval,3).denominator != 3:
                raise AssertionError('even-b reduced cubic slope')
            two_hull_cases += 1
    require('two_adic_tame_cubic_Newton_controls',two_hull_cases==98,hull_cases=two_hull_cases,
            c_valuations='0 through12 and zero constant',pi_valuations='0 through3 and the three even-b cases')

    def two_coefficients(p, target_b, target_c):
        return [-9*target_c,18,-(40*p+9*target_b),48*p,-48*p**4,32*p**5]

    def scaled_reduction(p,target_b,target_c,r,power):
        coeff=[Fraction(a*2**power,2**(r*i))
               for i,a in enumerate(two_coefficients(p,target_b,target_c))]
        if any(a.denominator!=1 for a in coeff):
            raise AssertionError(('nonintegral scaled coefficient',p,target_b,target_c,r,power))
        return [int(a)%2 for a in coeff]

    automatic_cases=0
    odd_units=[-3,-1,1,3]
    odd_b=[-5,-3,-1,1,3,5]
    for e in [2,4,5,6,7,8,9]:
        for unit in odd_units:
            for target_b in odd_b:
                for target_c in range(-3,4):
                    r,power,expected=(5,10,[0,0,1,0,0,1]) if e==2 else (e+4,2*(e+4),[0,0,1,1,0,0])
                    actual=scaled_reduction(2**e*unit,target_b,target_c,r,power)
                    if actual!=expected:
                        raise AssertionError(('odd-b simple residual',e,unit,target_b,target_c,actual))
                    derivative=sum(i*actual[i] for i in range(1,6))%2
                    if sum(actual)%2 or derivative!=1:
                        raise AssertionError('simple residue root one')
                    automatic_cases+=1
    for unit in odd_units:
        for half_b in odd_b:
            for target_c in range(-3,4):
                actual=scaled_reduction(2*unit,2*half_b,target_c,3,5)
                if actual!=[0,0,1,0,0,1]:
                    raise AssertionError(('even-b simple residual',unit,half_b,target_c,actual))
                automatic_cases+=1
    for unit in odd_units:
        for quarter_b in range(-4,5):
            for half_c in range(-3,4):
                actual=scaled_reduction(unit,4*quarter_b,2*half_c,1,0)
                if actual!=[0,1,quarter_b%2,0,1,1]:
                    raise AssertionError(('simple root zero',unit,quarter_b,half_c,actual))
                automatic_cases+=1
    require('two_adic_automatic_simple_root_residues',automatic_cases==1596,residue_cases=automatic_cases,
            pi_unit_controls=odd_units,odd_b_unit_controls=odd_b,c_controls='-3 through3 or twice this range')

    Z=sp.symbols('Z')
    p5=sp.Poly(Z**5+Z**4+Z**2+Z+1,Z,modulus=2)
    require('two_adic_residual_quintic_irreducible',p5.is_irreducible)
    require('two_adic_residual_quintic_no_linear',p5.eval(0)!=0 and p5.eval(1)!=0)
    require('two_adic_residual_quintic_no_quadratic',sp.rem(p5,sp.Poly(Z**2+Z+1,Z,modulus=2)).as_expr()==1)
    irreducible_cases=0
    for unit in odd_units:
        for quarter_b in odd_b:
            for target_c in odd_b:
                if scaled_reduction(unit,4*quarter_b,target_c,1,0)!=[1,1,1,0,1,1]:
                    raise AssertionError(('irreducible residual',unit,quarter_b,target_c))
                irreducible_cases+=1
    require('two_adic_irreducible_quintic_residue_controls',irreducible_cases==144,residue_cases=irreducible_cases)
    equal('two_adic_tame_cubic_discriminant',sp.discriminant(Z**3-2,Z),-108)
    require('two_adic_resolvent_nonsplit_residue',sp.Poly(Z**2+Z+1,Z,modulus=2).is_irreducible)
    require('two_adic_unit_cube_simple_root',(1**3-1)%2==0 and (3*1**2)%2==1)

    # The proper S3 decomposition-group step in Corollary ILO3a.
    import itertools
    permutations = list(itertools.permutations(range(3)))
    identity = (0,1,2)
    compose = lambda first,second: tuple(first[second[i]] for i in range(3))
    group_controls = 0
    for mask in range(1<<6):
        subgroup = [permutations[i] for i in range(6) if mask & (1<<i)]
        if identity not in subgroup or len(subgroup) == 6:
            continue
        if not all(compose(g,h) in subgroup for g in subgroup for h in subgroup):
            continue
        def act(g, point):
            if point < 2:
                odd = sum(g[i]>g[j] for i in range(3) for j in range(i+1,3)) % 2
                return point ^ odd
            return 2+g[point-2]
        remaining=set(range(5))
        sizes=[]
        while remaining:
            point=min(remaining)
            orbit={act(g,point) for g in subgroup}
            sizes.append(len(orbit))
            remaining-=orbit
        if 2 in sizes and not (len(subgroup)==2 and sorted(sizes)==[1,2,2]):
            raise AssertionError((subgroup,sizes))
        group_controls+=1
    require("proper_S3_quadratic_orbit_controls", group_controls==5, proper_subgroups=group_controls)
    require("divisible_pi_theta_coefficient_valuations",
            all(e+sp.Rational(1,2)>0 and min(2*e-sp.Rational(3,2),e+sp.Rational(1,2))>0
                and 3*e>0 and 5*e-sp.Rational(1,2)>0 for e in range(1,50)))
    require("theta_unit_roots_simple_mod_three",
            all((t**5-t**3)%3==0 and (5*t**4-3*t**2)%3!=0 for t in (1,-1)))
    require("ramified_field_unit_minus_two_square", (-2)%3==1 and (2*1)%3!=0)
    require("unit_pi_resolvent_classes", 7%3==1 and (-8)%3==1)
    # A compatible low-prime example remains outside the Hasse locus.
    from math import isqrt
    escape = -32*W**5-48*W**4-48*W**3+4*W**2+18*W-108
    escape_derivative = sp.diff(escape,W)
    escape_disc = int(sp.discriminant(escape,W))
    escape_v2 = valuation(escape_disc,2)
    escape_v3 = valuation(escape_disc,3)
    hensel2 = (valuation(int(escape.subs(W,2)),2),valuation(int(escape_derivative.subs(W,2)),2))
    hensel3 = (valuation(int(escape.subs(W,5)),3),valuation(int(escape_derivative.subs(W,5)),3))
    require('explicit_integral_low_prime_escape_example',
            escape_disc == 2**24*3**6*50260753 and escape_v2==24 and escape_v3==6
            and (escape_disc//2**escape_v2)%8==1 and (escape_disc//3**escape_v3)%3==1
            and isqrt(50260753)**2 != 50260753 and hensel2==(3,1) and hensel3==(4,1)
            and all(int(escape.subs(W,a))%11 for a in range(11)),
            target=[1,4,12],discriminant='2^24*3^6*50260753',
            two_adic_Hensel_point_and_valuations=[2,*hensel2],
            three_adic_Hensel_point_and_valuations=[5,*hensel3],
            no_residue_roots_modulo=11)
    source=Path(__file__).with_name("FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md")
    report={
        "status":"passed",
        "scope":"Exact residue/discriminant controls and bounded Newton-hull controls; not a proof of unrestricted integral Hasse nonexistence or an independent whole-proof replay.",
        "python":platform.python_version(),"sympy":sp.__version__,
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "check_count":len(checks),"checks":checks,
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"check_count":len(checks),"output":str(args.output) if args.output else None}))


if __name__=="__main__":
    main()
