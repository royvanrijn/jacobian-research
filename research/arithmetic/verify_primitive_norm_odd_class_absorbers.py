#!/usr/bin/env python3
"""Finite odd-class absorbers and listed-prime densities; no analytic claim."""
import argparse
from collections import defaultdict
from hashlib import sha256
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, lcm
from pathlib import Path
import json


class Group:
    def __init__(self, moduli):
        self.moduli=tuple(moduli)
        self.zero=tuple(0 for _ in moduli)
        self.elements=set(product(*(range(n) for n in moduli)))
    def add(self,a,b):
        return tuple((x+y)%n for x,y,n in zip(a,b,self.moduli))
    def neg(self,a):
        return tuple(-x%n for x,n in zip(a,self.moduli))
    def multiple(self,k,a):
        return tuple(k*x%n for x,n in zip(a,self.moduli))
    def order(self,a):
        return lcm(*(n//gcd(x,n) for x,n in zip(a,self.moduli)))


def signed_counts(group, variables):
    counts={group.zero:1}
    for variable in variables:
        new=defaultdict(int)
        for value,count in counts.items():
            for orientation in (variable,group.neg(variable)):
                new[group.add(value,orientation)]+=count
        counts=dict(new)
    assert sum(counts.values())==2**len(variables)
    return counts


def block_audit(group,generator):
    h=group.order(generator)
    assert h%2==1
    m=h-1
    expected={group.multiple(j,generator) for j in range(h)}
    base=signed_counts(group,[generator]*m)
    assert set(base)==expected
    checks=0
    for labels in product((-1,1),repeat=m):
        variables=[group.multiple(sign,generator) for sign in labels]
        assert signed_counts(group,variables)==base
        checks+=1
    return dict(order=h,assigned_count=m,initial_label_patterns_checked=checks,
                orientation_counts=[base[group.multiple(j,generator)] for j in range(h)])


def validate_pools(pools, orders):
    if len(pools)!=len(orders) or any(h<1 or h%2==0 for h in orders):
        raise ValueError('pool orders must be positive odd and match the pools')
    labels=[label for pool in pools for label in pool]
    if len(labels)!=len(set(labels)):
        raise ValueError('overlapping prime pool labels')
    if any(len(pool)!=2*(h-1) for pool,h in zip(pools,orders)):
        raise ValueError('wrong pool size')


def assign(pool,m,dividing):
    div=[q for q in pool if q in dividing]
    others=[q for q in pool if q not in dividing]
    if len(div)>=m:
        first=second=set(div[:m])
    else:
        k=len(div)
        first=set(div+others[:m-k])
        second=set(div+others[m-k:2*(m-k)])
    return first,second


def validate_profile(pools,orders,dividing,first,second,residual_coprime=True):
    validate_pools(pools,orders)
    all_labels={q for pool in pools for q in pool}
    if not residual_coprime:
        raise ValueError('residual factors can change assigned exponents')
    if (first|second)-all_labels:
        raise ValueError('prime outside reserved pool')
    for pool,h in zip(pools,orders):
        if len(first&set(pool))!=h-1 or len(second&set(pool))!=h-1:
            raise ValueError('missing absorber variables')
    if (first&second)-dividing:
        raise ValueError('common assignment at target unit')
    if ((first^second)&dividing):
        raise ValueError('exclusive assignment at dividing target')


def case_audit(name,moduli,generators):
    group=Group(moduli)
    generators=[tuple(g) for g in generators]
    blocks=[block_audit(group,g) for g in generators]
    orders=[group.order(g) for g in generators]
    pools=[]
    label_classes={}
    offset=0
    for g,h in zip(generators,orders):
        pool=list(range(offset,offset+2*(h-1)))
        offset+=len(pool)
        pools.append(pool)
        for q in pool:
            label_classes[q]=group.multiple(-1 if q%2 else 1,g)
    validate_pools(pools,orders)
    assert set(signed_counts(group,[g for g,h in zip(generators,orders) for _ in range(h-1)]))==group.elements
    digest=sha256()
    compensation_pairs=0
    for mask in range(1<<offset):
        dividing={q for q in range(offset) if mask>>q&1}
        first,second=set(),set()
        for pool,h in zip(pools,orders):
            a,b=assign(pool,h-1,dividing)
            first|=a
            second|=b
        validate_profile(pools,orders,dividing,first,second)
        left=signed_counts(group,[label_classes[q] for q in sorted(first)])
        right=signed_counts(group,[label_classes[q] for q in sorted(second)])
        assert set(left)==set(right)==group.elements
        # Shared rational-prime variables have independent ideal orientations
        # in the two values; challenge every pair of residual classes.
        for c,d in product(group.elements,repeat=2):
            assert left.get(group.neg(c),0)>0 and right.get(group.neg(d),0)>0
            compensation_pairs+=1
        first_mask=sum(1<<q for q in first)
        second_mask=sum(1<<q for q in second)
        row=[mask,first_mask,second_mask,first_mask&second_mask]
        digest.update((json.dumps(row,separators=(',',':'))+'\n').encode())
    return dict(name=name,moduli=list(moduli),generators=generators,
                group_size=len(group.elements),generator_blocks=blocks,
                globally_disjoint_pool_size=offset,
                complete_target_divisibility_patterns=1<<offset,
                residual_class_pairs_compensated=compensation_pairs,
                all_two_sided_signed_supports_complete=True,
                assignment_rows_sha256=digest.hexdigest())


def expect_rejection(action):
    try:
        action()
    except ValueError:
        return True
    raise AssertionError('invalid profile accepted')


def negative_audits():
    c2=Group((2,))
    assert set(signed_counts(c2,[(1,)]))!={c2.zero,(1,)}
    c5=Group((5,))
    assert len(signed_counts(c5,[(1,)]*3))==4<5
    pools=[[0,1,2,3]]
    orders=[3]
    assert expect_rejection(lambda:validate_pools([[0,1]],[2]))
    assert expect_rejection(lambda:validate_pools(pools,[]))
    assert expect_rejection(lambda:validate_pools([[0,1,2,3],[3,4,5,6]],[3,3]))
    assert expect_rejection(lambda:validate_profile(pools,orders,set(),{0,1},{0,1}))
    assert expect_rejection(lambda:validate_profile(pools,orders,{0},{0,1},{2,3}))
    assert expect_rejection(lambda:validate_profile(pools,orders,set(),{0},{2,3}))
    assert expect_rejection(lambda:validate_profile(pools,orders,set(),{0,1},{2,3},False))
    c3x3=Group((3,3))
    assert set(signed_counts(c3x3,[(1,0)]*2)) != c3x3.elements
    c3=Group((3,))
    assert c3.zero not in signed_counts(c3,[(0,),(1,)])
    assert sum((a+b)%2==1 for a in (1,) for b in (1,))==0
    return [
        'even order loses invertibility of two and is rejected',
        'mismatched prime pools and group orders rejected',
        'nongenerating classes leave a residual quotient obstruction',
        'too few same-class exact-one variables cannot span cyclic block',
        'overlapping prime-pool labels rejected',
        'common assignment at unit target rejected',
        'exclusive assignment at dividing target rejected',
        'missing absorber variable rejected',
        'missing residual coprimality guard rejected',
        'class3 exponent changed from one to three destroys absorption',
        'dyadic unit-unit local obstruction retained',
    ]


def local_unit_checks():
    checked=0
    for p in (3,5,7,11):
        for u,v in ((1,1),(0,1),(1,0)):
            for target in range(p):
                count=sum((u*a+v*b)%p==target for a in range(1,p) for b in range(1,p))
                expected=(p-1 if target==0 else p-2) if u*v else (0 if target==0 else p-1)
                assert count==expected
                checked+=1
    return checked



def prime_character(p,n):
    if n%p==0:
        return 0
    value=pow(n,(p-1)//2,p)
    assert value in (1,p-1)
    return 1 if value==1 else -1


def validate_prime_gamma_mode(p,t):
    if p<7 or p%4!=3 or any(p%d==0 for d in range(2,isqrt(p)+1)):
        raise ValueError('requires a prime p>=7 with p=3 mod4')
    if t not in (1,2) or (t==2 and prime_character(p,2)!=1):
        raise ValueError('t2 requires split two, chi_p(2)=1')


def raw_prime_line_counts(p,t):
    """Finite residues, also usable to challenge an inadmissible t2 mode."""
    period=2*p*t
    subgroup=[a for a in range(period)
              if gcd(a,period)==1 and prime_character(p,a)==1]
    counts=[0]*period
    for a in subgroup:
        for b in subgroup:
            counts[(a+t*b)%period]+=1
    return subgroup,counts


def prime_gamma_audit():
    """OC.10 is tested only at these eight primes; no extrapolation."""
    primes=(7,11,19,23,31,43,59,67)
    records=[]
    qr_cases=gamma_cases=coefficient_permutation_checks=0
    invalid_t2=[]
    for p in primes:
        validate_prime_gamma_mode(p,1)
        qr={a*a%p for a in range(1,p)}
        assert len(qr)==(p-1)//2 and prime_character(p,-1)==-1
        assert all((a in qr)==(prime_character(p,a)==1) for a in range(p))
        pair_counts=[0]*p
        for a in qr:
            for b in qr:
                pair_counts[(a+b)%p]+=1
        for target,count in enumerate(pair_counts):
            expected=0 if target==0 else (p-1-2*prime_character(p,target))//4
            assert count==expected
            if target:
                assert count>=1
            qr_cases+=1
        modes=[1,2] if prime_character(p,2)==1 else [1]
        cases=[]
        for t in modes:
            validate_prime_gamma_mode(p,t)
            period=2*p*t
            subgroup,counts=raw_prime_line_counts(p,t)
            assert len(subgroup)==t*(p-1)//2
            for u in subgroup:
                assert {u*a%period for a in subgroup}==set(subgroup)
                coefficient_permutation_checks+=1
            masses=[]
            for target,count in enumerate(counts):
                actual=Fraction(period*count,len(subgroup)**2)
                eligible=target%2==(0 if t==1 else 1) and target%p!=0
                expected=(Fraction(8*p,(p-1)**2)*pair_counts[target%p]
                          if eligible else Fraction(0))
                assert actual==expected,(p,t,target,actual,expected)
                masses.append([actual.numerator,actual.denominator])
                gamma_cases+=1
            if t==2:
                assert masses[1][0]>0
                actual=Fraction(*masses[1])
                assert actual/2 != Fraction(8*p,(p-1)**2)*pair_counts[1]
            cases.append(dict(t=t,period=period,subgroup_size=len(subgroup),
                              pair_counts=counts,Gamma_fractions=masses))
        if prime_character(p,2)==-1:
            assert expect_rejection(lambda:validate_prime_gamma_mode(p,2))
            subgroup,counts=raw_prime_line_counts(p,2)
            actual=Fraction(4*p*counts[p],len(subgroup)**2)
            assert actual==Fraction(4*p,p-1)>0
            invalid_t2.append(dict(p=p,target=p,mode_rejected=True,
                                  inadmissible_mode_ramified_mass=str(actual)))
        records.append(dict(p=p,chi2=prime_character(p,2),
                            square_pair_counts=pair_counts,cases=cases))
    # These challenge the prime threshold, Legendre sign, parity and t2 lifts.
    assert expect_rejection(lambda:validate_prime_gamma_mode(3,1))
    qr3={1}
    assert sum((a+b)%3==1 for a in qr3 for b in qr3)==0
    assert records[0]['square_pair_counts'][1] != (7-1+2*prime_character(7,1))//4
    assert raw_prime_line_counts(7,1)[1][1]==0 and records[0]['square_pair_counts'][1]>0
    assert len(invalid_t2)==5
    assert qr_cases==260 and gamma_cases==764 and coefficient_permutation_checks==184
    return dict(
        scope='Exact finite QR-pair, subgroup and Gamma checks only at the eight listed primes. The unrestricted elementary identity remains a written proof input.',
        primes=list(primes),square_pair_target_cases=qr_cases,
        Gamma_target_cases=gamma_cases,
        coefficient_subgroup_permutation_checks=coefficient_permutation_checks,
        cases=records,invalid_t2_when_two_is_inert=invalid_t2,
        negative_controls=['p3 is rejected and has a zero unit-target square-pair count',
                           'wrong Legendre sign changes the count',
                           'dropping the parity filter admits false targets',
                           'dropping the t2 lift factor changes Gamma',
                           't2 is rejected when chi2=-1; its ramified mass would be nonzero'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_odd_class_absorber_checks.json'))
    args=parser.parse_args()
    cases=[case_audit(name,moduli,generators) for name,moduli,generators in (
        ('trivial',(1,),[]),
        ('C3',(3,),[(1,)]),
        ('C5',(5,),[(1,)]),
        ('C7',(7,),[(1,)]),
        ('C3xC3',(3,3),[(1,0),(0,1)]),
        ('C3xC3_redundant',(3,3),[(1,0),(0,1),(1,1)]),
    )]
    composite=block_audit(Group((9,)),(1,))
    result=dict(
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Finite abstract group, assignment and local residue audits, plus QR/Gamma checks at the listed eight primes only. No class-prime existence theorem, analytic additive theorem or numerical cutoff is certified.',
        cases=cases,
        composite_odd_order_C9_block=composite,
        local_unit_residue_checks=local_unit_checks(),
        negative_controls=negative_audits(),
        prime_Gamma_audit=prime_gamma_audit(),
        total_complete_assignment_profiles=sum(c['complete_target_divisibility_patterns'] for c in cases),
        total_residual_class_pairs_compensated=sum(c['residual_class_pairs_compensated'] for c in cases))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(
        group_cases=len(cases),
        complete_assignment_profiles=result['total_complete_assignment_profiles'],
        residual_class_pairs_compensated=result['total_residual_class_pairs_compensated'],
        local_unit_residue_checks=result['local_unit_residue_checks'],
        group_negative_controls=len(result['negative_controls']),
        square_pair_target_cases=result['prime_Gamma_audit']['square_pair_target_cases'],
        Gamma_target_cases=result['prime_Gamma_audit']['Gamma_target_cases'],
        Gamma_mode_rejections=len(result['prime_Gamma_audit']['invalid_t2_when_two_is_inert']),
        output=str(args.output)),indent=2))


if __name__=='__main__':
    main()
