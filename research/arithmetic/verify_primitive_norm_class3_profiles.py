#!/usr/bin/env python3
"""Exact finite profiles for the discriminant -23 absorber argument.

Standard library only. These checks certify finite prime witnesses,
residue densities and assignment profiles. The class-group interpretation
and the analytic counting theorem require the separate written proof.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import json

PRIMES=(3,13,29,31)
Q1_WITNESSES=((0,1),(2,1),(-2,3),(-4,1))
ABSORBER_PRODUCT=prod(PRIMES)
QR={x*x%23 for x in range(1,23)}
assert ABSORBER_PRODUCT==35061
assert len(QR)==11 and 2 in QR and 22 not in QR


def chi(n):
    return 0 if n%23==0 else (1 if n%23 in QR else -1)


def q0_representations(n):
    """Complete signed search using (2x+y)^2+23y^2=4n."""
    result=[]
    y_bound=isqrt(4*n//23)
    z_bound=isqrt(4*n)
    for y in range(-y_bound,y_bound+1):
        for z in range(-z_bound,z_bound+1):
            if (z-y)%2==0 and z*z+23*y*y==4*n:
                x=(z-y)//2
                assert x*x+x*y+6*y*y==n
                result.append([x,y])
    return result


def prime_certificates():
    result=[]
    for p,(x,y) in zip(PRIMES,Q1_WITNESSES):
        assert all(p%d for d in range(2,isqrt(p)+1))
        assert pow(p,11,23)==1 and chi(p)==1
        assert gcd(x,y)==1 and 2*x*x+x*y+3*y*y==p
        assert not q0_representations(p)
        result.append(dict(
            prime=p,nonprincipal_form=[2,1,3],
            nonprincipal_form_witness=[x,y],principal_form=[1,1,6],
            principal_form_representations=[],
            complete_y_search_bound=isqrt(4*p//23),
            complete_z_search_bound=isqrt(4*p)))
    assert q0_representations(39)
    assert q0_representations(351)
    assert not any(gcd(x,y)==1 for x,y in q0_representations(351))
    assert 351==39*9 and gcd(9,ABSORBER_PRODUCT)>1
    return dict(
        certificates=result,
        exact_exponent_negative_example=dict(
            assigned_product=39,invalid_residual=9,value=351,
            residual_coprime_guard_rejected=True,
            primitive_principal_form_representations=[],
            imprimitive_principal_form_representations=q0_representations(351)))


def expected_gamma(t,target):
    eligible=target%2==(0 if t==1 else 1) and target%23!=0
    return F(46,121)*((22-2*chi(target))//4) if eligible else F(0)


def gamma_checks():
    cases=[]
    negatives={}
    for t in (1,2):
        period=46*t
        subgroup=[a for a in range(period)
                  if gcd(a,period)==1 and chi(a)==1]
        assert len(subgroup)==11*t
        counts=[0]*period
        for a in subgroup:
            for b in subgroup:
                counts[(a+t*b)%period]+=1
        rows=[]
        for target,count in enumerate(counts):
            actual=F(period*count,len(subgroup)**2)
            expected=expected_gamma(t,target)
            assert actual==expected,(t,target,actual,expected)
            rows.append([target,count,actual.numerator,actual.denominator])
        negatives[f"t{t}_nonzero_ramified_target_rejected"]=(
            expected_gamma(t,0 if t==1 else 23)==0)
        if t==2:
            actual=F(period*counts[1],len(subgroup)**2)
            assert actual>0 and actual/2!=expected_gamma(t,1)
            negatives["missing_t2_lift_factor_detected"]=True
        cases.append(dict(t=t,period=period,subgroup_size=len(subgroup),
                          columns=["target","pair_count","Gamma_numerator","Gamma_denominator"],
                          rows=rows))
    assert sum(len(case["rows"]) for case in cases)==138
    assert all(negatives.values())
    return dict(checks=138,cases=cases,negative_controls=negatives)


def expected_profile_mass(p,a,b,target):
    if (a,b)==(0,0):
        return F(p-1 if target%p==0 else p-2,p)
    if (a,b) in ((1,0),(0,1)):
        return F(2*(p-1),p*p) if target%p else F(0)
    return (F(0) if target%p else
            F(4*(p-1 if target==0 else p-2),p*p))


def profile_mass_checks():
    checks=0
    tables=[]
    for p in PRIMES:
        modulus=p*p
        residues={0:[n for n in range(modulus) if n%p],
                  1:[n for n in range(modulus) if n%p==0 and n]}
        for t in (1,2):
            for a,b in product((0,1),repeat=2):
                counts=[0]*modulus
                for x in residues[a]:
                    for y in residues[b]:
                        counts[(x+t*y)%modulus]+=2**(a+b)
                examples={}
                for target,count in enumerate(counts):
                    actual=F(count,modulus)
                    expected=expected_profile_mass(p,a,b,target)
                    assert actual==expected,(p,t,a,b,target,actual,expected)
                    checks+=1
                    if target in (0,1,p):
                        examples[str(target)]=str(actual)
                # The receipt records every numerator: the common
                # denominator is modulus, rather than only sample masses.
                tables.append(dict(
                    p=p,t=t,profile=[a,b],modulus=modulus,
                    weighted_line_counts=counts,
                    target_representatives=examples))
    assert checks==15840
    assert expected_profile_mass(3,1,1,3)==F(4,9)
    assert F(1,9)!=expected_profile_mass(3,1,1,3)
    assert expected_profile_mass(3,1,0,3)==0
    return dict(
        checks=checks,tables=tables,
        negative_controls=dict(
            missing_exact_one_weights_detected=True,
            one_sided_assignment_at_dividing_target_rejected=True))


def valid_template(profiles,dividing):
    if set(profiles)!=set(PRIMES):
        return False
    if any(tuple(pair) not in ((0,0),(1,0),(0,1),(1,1))
           for pair in profiles.values()):
        return False
    if any(sum(pair[side] for pair in profiles.values())<2
           for side in (0,1)):
        return False
    for p,pair in profiles.items():
        if tuple(pair)==(1,1) and p not in dividing:
            return False
        if tuple(pair) in ((1,0),(0,1)) and p in dividing:
            return False
    return True


def template_checks():
    templates=[]
    normalized_local_checks=0
    for pattern in product((False,True),repeat=4):
        dividing=[p for p,yes in zip(PRIMES,pattern) if yes]
        nondividing=[p for p,yes in zip(PRIMES,pattern) if not yes]
        profiles={p:[0,0] for p in PRIMES}
        if len(dividing)>=2:
            for p in dividing[:2]:profiles[p]=[1,1]
        elif dividing:
            profiles[dividing[0]]=[1,1]
            profiles[nondividing[0]]=[1,0]
            profiles[nondividing[1]]=[0,1]
        else:
            for p in nondividing[:2]:profiles[p]=[1,0]
            for p in nondividing[2:]:profiles[p]=[0,1]
        assert valid_template(profiles,dividing)
        assert all(sum(pair[side] for pair in profiles.values())==2
                   for side in (0,1))
        first=[p for p in PRIMES if profiles[p][0]]
        second=[p for p in PRIMES if profiles[p][1]]
        a,b=prod(first),prod(second)
        common=gcd(a,b)
        u,v=a//common,b//common
        assert gcd(u,v)==1 and max(a,b)<=899
        # Exact exponent one requires residuals coprime to the whole
        # fixed pool; a selected coefficient itself is squarefree.
        for p in first:assert a%p==0 and (a//p)%p
        for p in second:assert b%p==0 and (b//p)%p
        local_tables={}
        for p in PRIMES:
            for t in (1,2):
                counts=[]
                for target in range(p):
                    count=sum((u*x+t*v*y-target)%p==0
                              for x in range(1,p) for y in range(1,p))
                    if u%p==0 or v%p==0:
                        expected=0 if target==0 else p-1
                        assert p not in dividing
                    else:
                        expected=p-1 if target==0 else p-2
                    assert count==expected
                    normalized_local_checks+=1
                    counts.append(count)
                local_tables[f"{p}:{t}"]=counts
        templates.append(dict(
            target_dividing_primes=dividing,profiles=profiles,
            first_forced_primes=first,second_forced_primes=second,
            first_coefficient=a,second_coefficient=b,common_coefficient=common,
            normalized_coefficients=[u,v],
            required_residual_coprimality_modulus=ABSORBER_PRODUCT,
            normalized_unit_pair_counts=local_tables))
    assert len(templates)==16
    missing={p:[0,0] for p in PRIMES}
    assert not valid_template(missing,[])
    incompatible={3:[1,1],13:[1,1],29:[0,0],31:[0,0]}
    assert not valid_template(incompatible,[])
    assert not valid_template({3:[3,1],13:[1,1],29:[0,0],31:[0,0]},[3,13])
    return dict(
        count=len(templates),templates=templates,
        normalized_local_checks=normalized_local_checks,
        negative_controls=dict(
            missing_active_primes_rejected=True,
            common_assignment_at_unit_target_rejected=True,
            exponent_three_relaxation_rejected=True,
            residual_pool_divisibility_rejected=True))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",type=Path,
        default=Path(__file__).with_name("primitive_norm_class3_profile_checks.json"))
    args=parser.parse_args()
    result=dict(
        scope="Complete finite checks of four prime witnesses, 138 Gamma "
              "cases, 15840 local profile identities and 16 templates. "
              "The class-group criterion, analytic estimates and unweighting "
              "remain separate written proof inputs.",
        absorber_primes=list(PRIMES),absorber_product=ABSORBER_PRODUCT,
        prime_certificates=prime_certificates(),
        Gamma=gamma_checks(),local_profiles=profile_mass_checks(),
        assignment_templates=template_checks(),
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(
        Gamma_checks=result["Gamma"]["checks"],
        local_profile_checks=result["local_profiles"]["checks"],
        templates=result["assignment_templates"]["count"],
        normalized_local_checks=result["assignment_templates"]["normalized_local_checks"],
        output=str(args.output)),indent=2))


if __name__=="__main__":
    main()
