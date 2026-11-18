#!/usr/bin/env python3
"""Exact finite genus, profile and ramification checks for discriminant -56.

Standard library only. The written class-group and analytic proofs remain
separate inputs; these finite checks do not establish an eventual cutoff.
"""


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import json

L=56
PRIMES=(3,5,13,19)
ABSORBER_PRODUCT=prod(PRIMES)
Q1_WITNESSES=((1,0),(0,1),(-2,1),(-1,2))
QR7={1,2,4}


def psi(n):
    return 0 if n%2==0 else (1 if n%8 in (1,7) else -1)


def chi7(n):
    return 0 if n%7==0 else (1 if n%7 in QR7 else -1)


def chi(n):
    return psi(n)*chi7(n)


def subgroup(period):
    return [a for a in range(period) if gcd(a,period)==1 and chi(a)==1]


H=subgroup(L)
A={sign:[a for a in H if psi(a)==sign] for sign in (1,-1)}
assert H==[1,3,5,9,13,15,19,23,25,27,39,45]
assert A[1]==[1,9,15,23,25,39]
assert A[-1]==[3,5,13,19,27,45]
assert ABSORBER_PRODUCT==3705


def q0_representations(n):
    result=[]
    for y in range(-isqrt(n//14),isqrt(n//14)+1):
        rem=n-14*y*y
        x=isqrt(rem)
        if x*x==rem:
            result.append([x,y])
            if x:result.append([-x,y])
    assert all(x*x+14*y*y==n for x,y in result)
    return result


def class_and_prime_certificates():
    prime_rows=[]
    for p,(x,y) in zip(PRIMES,Q1_WITNESSES):
        assert all(p%d for d in range(2,isqrt(p)+1))
        assert chi(p)==1 and psi(p)==-1
        assert gcd(x,y)==1 and 3*x*x+2*x*y+5*y*y==p
        assert not q0_representations(p)
        prime_rows.append(dict(
            prime=p,chi_minus56=1,genus_chi8=-1,
            odd_class_form=[3,2,5],witness=[x,y],
            principal_form_representations=[],
            complete_y_search_bound=isqrt(p//14)))
    assert not q0_representations(2) and not q0_representations(7)
    assert [0,1] in q0_representations(14)
    assert [2,1] in q0_representations(18)
    assert (2+1)%3==0 and (2-1)%3!=0
    support={sum(signs)%4 for signs in product((-1,1),repeat=2)}
    assert support=={0,2}
    compensation=[]
    for prefix in (1,2,7,14):
        prefix_class=2*((prefix%2==0)+(prefix%7==0))%4
        for residual_class in (0,2):
            choices=sorted(s for s in support
                           if (s+prefix_class+residual_class)%4==0)
            assert choices
            compensation.append(dict(
                norm_prefix=prefix,prefix_class=prefix_class,
                residual_class=residual_class,absorber_choices=choices))
        for wrong_genus_class in (1,3):
            assert not any((s+prefix_class+wrong_genus_class)%4==0
                           for s in support)
    ordinary=q0_representations(450)
    assert ordinary and not any(gcd(x,y)==1 for x,y in ordinary)
    assert 450==2*(3*5)**2 and psi(225)==1
    assert gcd(15,ABSORBER_PRODUCT)>1
    return dict(
        scope="Finite integer and abstract C4 checks. Identifying the actual "
              "ideal class group requires the separate Minkowski/ideal proof.",
        prime_certificates=prime_rows,abstract_absorber_support=sorted(support),
        ramified_prefix_compensation=compensation,
        ideal_relation_integer_checks=dict(
            norm_sqrt_minus14=14,norm_2_plus_sqrt_minus14=18,
            theta_mod3_selected_root=1,other_root=-1),
        negative_controls=dict(
            wrong_genus_cannot_be_absorbed=True,
            no_global_norm2_or_norm7_generator=True,
            restoration_guard_rejects_450=True),
        restoration_guard_example=dict(
            value=450,ramified_prefix=2,assigned_product=15,
            invalid_residual=15,core_genus=1,
            imprimitive_representations=ordinary,
            primitive_representations=[]))


def finite_density(t,u,v,alpha,beta):
    period=t*L
    units=subgroup(period)
    counts=[0]*period
    for a in units:
        if psi(a)!=alpha:continue
        for b in units:
            if psi(b)==beta:counts[(u*a+t*v*b)%period]+=1
    return [F(period*n,len(units)**2) for n in counts],counts


def expected_plus_q(t,target):
    if target%7==0:return 0
    square_pairs=(6-2*chi7(target))//4
    if t==1:
        dyadic=2 if target%8==0 else (1 if target%8 in (2,6) else 0)
    else:
        dyadic=1 if target%2 else 0
    return dyadic*square_pairs


def density_and_transport_checks():
    cases=[]
    base_cache={}
    checks=coefficient_checks=common_checks=0
    for t in (1,2):
        period=t*L
        units=subgroup(period)
        assert len(units)==12*t
        base={}
        for alpha,beta in product((1,-1),repeat=2):
            density,counts=finite_density(t,1,1,alpha,beta)
            base[alpha,beta]=density
            cases.append(dict(
                t=t,alpha=alpha,beta=beta,period=period,
                subgroup_size=len(units),pair_counts=counts,
                Gamma_numerators=[value.numerator for value in density],
                Gamma_denominators=[value.denominator for value in density]))
            checks+=period
        base_cache[t]=base
        for (alpha,beta),density in base.items():
            for u in units:
                for v in units:
                    transformed,_=finite_density(
                        t,u,v,alpha*psi(u),beta*psi(v))
                    assert transformed==density
                    coefficient_checks+=1
            for common in units:
                transformed=base[alpha*psi(common),beta*psi(common)]
                inverse=pow(common,-1,period)
                for target in range(period):
                    assert density[target]==transformed[target*inverse%period]
                    common_checks+=1
        for target,value in enumerate(base[1,1]):
            q=sum((target-t*b)%L in A[1] for b in A[1])
            assert value==F(7,18)*q
            assert q==expected_plus_q(t,target)
    assert checks==672 and coefficient_checks==2880 and common_checks==13440
    original=base_cache[1][1,1][2]
    reduced=2*pow(3,-1,56)%56
    assert reduced==38
    assert original!=base_cache[1][1,1][reduced]
    assert original==base_cache[1][-1,-1][reduced]
    actual=base_cache[2][1,1][1]
    assert actual>0 and actual/2!=F(7,18)*expected_plus_q(2,1)
    return dict(
        Gamma_checks=checks,coefficient_transports=coefficient_checks,
        common_target_transports=common_checks,H56=H,A_plus=A[1],A_minus=A[-1],
        cases=cases,
        negative_controls=dict(
            unchanged_genus_common_division_rejected=True,
            missing_t2_lift_factor_detected=True),
        sign_transport_witness=dict(
            original_target=2,dividing_residue=3,reduced_target_mod56=38,
            original_Gamma=str(original),
            wrong_unchanged_sign_Gamma=str(base_cache[1][1,1][reduced]),
            corrected_Gamma=str(base_cache[1][-1,-1][reduced]))),base_cache


def local_crt_and_target_selection(base_cache):
    modulus=784
    cores=[n for n in range(modulus)
           if gcd(n,14)==1 and chi(n)==1 and psi(n)==1]
    assert len(cores)==84
    expected_values={prefix*n%modulus for prefix in (1,2,7,14)
                     for n in cores}
    actual_values={(x*x+14*y*y)%modulus
                   for x in range(modulus) for y in range(modulus)
                   if gcd(gcd(x,y),14)==1}
    assert actual_values==expected_values and len(actual_values)==144
    sums={(a+b)%modulus for a in actual_values for b in actual_values}
    targets={n for n in range(modulus) if n%49}
    assert sums==targets and len(targets)==768
    rows=[]
    for target in sorted(targets):
        seven=1 if target%7==0 else 0
        reduced=target//(7**seven)
        assert reduced%7
        if reduced%2:
            t,common_two,analytic=2,0,reduced
        elif reduced%8==4:
            t,common_two,analytic=1,1,reduced//2
            assert analytic%8 in (2,6)
        else:
            t,common_two,analytic=1,0,reduced
            assert reduced%8 in (0,2,6)
        period=t*L
        value=base_cache[t][1,1][analytic%period]
        assert value>=F(7,18)
        units=subgroup(period)
        a,b=next((a,b) for a in units if psi(a)==1
                 for b in A[1] if (a+t*b-analytic)%period==0)
        first_prefix=7**seven*2**common_two
        second_prefix=t*first_prefix
        assert first_prefix in (1,2,7,14) and second_prefix in (1,2,7,14)
        rows.append(dict(
            target_mod784=target,common_ramified_seven=seven,
            common_ramified_two=common_two,t=t,
            analytic_target_mod_tL=analytic%period,
            core_residue_witness=[a,b],
            norm_prefixes=[first_prefix,second_prefix],
            Gamma_numerator=value.numerator,Gamma_denominator=value.denominator))
    assert len(rows)==768
    assert base_cache[1][1,1][4]==0
    assert base_cache[1][1,1][2]>0
    return dict(
        modulus=modulus,individual_local_residue_count=len(actual_values),
        exact_individual_residues=sorted(actual_values),
        complete_additive_admissibility="49 does not divide N",
        admissible_target_classes=len(targets),target_cases=rows,
        negative_controls=dict(
            omitted_common_two_stratum_detected=True,
            nonprincipal_ramification_requires_reorientation=True))


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
    assert checks==4512
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
    if any(sum(pair[side] for pair in profiles.values())!=2
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
        assert psi(a)==psi(b)==1
        common=gcd(a,b)
        u,v=a//common,b//common
        assert gcd(u,v)==1 and max(a,b)<=247
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
    incompatible={3:[1,1],5:[1,1],13:[0,0],19:[0,0]}
    assert not valid_template(incompatible,[])
    three_active={3:[1,1],5:[1,1],13:[1,1],19:[0,0]}
    assert not valid_template(three_active,[3,5,13])
    assert not valid_template({3:[2,1],5:[1,1],13:[0,0],19:[0,0]},[3,5])
    return dict(
        count=len(templates),templates=templates,
        normalized_local_checks=normalized_local_checks,
        negative_controls=dict(
            missing_active_primes_rejected=True,
            three_odd_class_primes_wrong_genus_rejected=True,
            common_assignment_at_unit_target_rejected=True,
            nonbinary_exponent_profile_rejected=True,
            residual_pool_divisibility_rejected=True))



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",type=Path,
        default=Path(__file__).with_name("primitive_norm_class4_profile_checks.json"))
    args=parser.parse_args()
    density,base_cache=density_and_transport_checks()
    result=dict(
        scope="Complete finite C4/genus, CRT and exact-profile checks for "
              "discriminant -56. The actual ideal-class identification, "
              "infinite analytic estimates and unweighting require the written proof.",
        absorber_primes=list(PRIMES),absorber_product=ABSORBER_PRODUCT,
        class_and_prime_certificates=class_and_prime_certificates(),
        finite_genus_density=density,
        local_CRT_and_target_selection=local_crt_and_target_selection(base_cache),
        local_profiles=profile_mass_checks(),
        assignment_templates=template_checks(),
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(
        Gamma_checks=density["Gamma_checks"],
        coefficient_transports=density["coefficient_transports"],
        common_target_transports=density["common_target_transports"],
        admissible_target_classes=result["local_CRT_and_target_selection"]["admissible_target_classes"],
        local_profile_checks=result["local_profiles"]["checks"],
        templates=result["assignment_templates"]["count"],
        normalized_local_checks=result["assignment_templates"]["normalized_local_checks"],
        output=str(args.output)),indent=2))


if __name__=="__main__":
    main()
