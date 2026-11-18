#!/usr/bin/env python3
"""Exact primitive discriminant -7 norm sumsets: complete finite proof checks.

Standard library only. Every invocation runs the full fixed proof domains.
The source-backed spectral and arithmetic proofs remain separate inputs;
these computations do not establish those infinite theorems by themselves.
"""
import argparse
from array import array
from bisect import bisect_right
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd,isqrt,prod
from pathlib import Path
import json
import sys

B=1<<24
X0=1<<18
R=1<<19
SCALE=1<<64
QR7={1,2,4}

def prime_flags(limit):
    flags=bytearray(b'\1')*(limit+1);flags[:2]=b'\0\0'
    for p in range(2,isqrt(limit)+1):
        if flags[p]:flags[p*p:limit+1:p]=b'\0'*((limit-p*p)//p+1)
    return flags


def iroot(n,k):
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while lo<hi:
        mid=(lo+hi+1)//2
        if mid**k<=n:lo=mid
        else:hi=mid-1
    return lo


def sigma_power(p,a):return sum(p**j for j in range(a+1))


def m_frontier(split,upper):
    data=[]
    def dfs(index,cap,n,value):
        data.append((n,value))
        if index>=len(split):return
        p=split[index]
        for a in range(1,cap+1):
            n*=p
            if n>upper:return
            dfs(index+1,a,n,value*(a+1))
    prefix=1;count=0
    for p in split:
        prefix*=p;count+=1
        if prefix>upper:break
    assert prefix>upper and 2**upper.bit_length()>upper
    dfs(0,upper.bit_length(),1,1)
    data.sort();out=[];best=0
    for n,value in data:
        if value>best:out.append((n,value));best=value
    return out,len(data),count


def rational_root_upper(value,k=6,scale=10**6):
    integer=(value.numerator*scale**k+value.denominator-1)//value.denominator
    root=iroot(integer,k)
    if root**k<integer:root+=1
    assert F(root,scale)**k>=value
    return F(root,scale)


def primes(limit):
    flags=prime_flags(limit)
    return flags,[p for p in range(2,limit+1) if flags[p]]



def coordinate_flags(bound):
    flags=bytearray(bound+1)
    for y in range(isqrt(4*bound//7)+1):
        root=isqrt(4*bound-7*y*y)
        for x in range(-(y//2),(root-y)//2+1):
            if gcd(x,y)==1:
                flags[x*x+x*y+2*y*y]=1
    flags[0]=0
    return flags


def semigroup_flags(bound):
    primes=prime_flags(bound)
    support=bytearray(b'\1')*(bound+1)
    for p in range(2,bound+1):
        if primes[p] and p%7 in (3,5,6):
            support[p:bound+1:p]=b'\0'*(bound//p)
    support[0:bound+1:49]=b'\0'*(bound//49+1)
    units=bytearray(support)
    units[0:bound+1:7]=b'\0'*(bound//7+1)
    restored=bytearray(units)
    restored[7:bound+1:7]=units[1:bound//7+1]
    assert restored==support
    return support,units


def trial_member(n,unit=False):
    if n<1:
        return False
    exponent=0
    while n%7==0:
        n//=7
        exponent+=1
    if exponent>(0 if unit else 1):
        return False
    p=2
    while p*p<=n:
        if n%p==0:
            if p%7 not in QR7:
                return False
            while n%p==0:
                n//=p
        p=3 if p==2 else p+2
    return n==1 or n%7 in QR7


def pack(flags):
    result=bytearray((len(flags)+7)//8)
    for n,yes in enumerate(flags):
        if yes:
            result[n>>3]|=1<<(n&7)
    return result


def positions(value):
    result=[]
    while value:
        low=value&-value
        result.append(low.bit_length()-1)
        value^=low
    return result


def require_equal(actual,expected):
    if actual!=expected:
        raise AssertionError('exact flag/mask mismatch')


def rejected(action):
    try:
        action()
    except AssertionError:
        return True
    raise AssertionError('corruption accepted')


def sumscan(flags,exclude_modulus,short_bound):
    bound=len(flags)-1
    packed=pack(flags)
    bits=int.from_bytes(packed,'little')
    anchors=[n for n in range(1,short_bound+1) if flags[n]]
    covered=0
    for a in anchors:
        covered|=bits<<a
    covered&=(1<<(bound+1))-1
    local=bytearray(b'\xff')*len(packed)
    local[0]&=~3
    for n in range(exclude_modulus,bound+1,exclude_modulus):
        local[n>>3]&=~(1<<(n&7))
    local[-1]&=(1<<((bound%8)+1))-1
    mask=int.from_bytes(local,'little')
    assert mask.bit_count()==bound-1-bound//exclude_modulus
    assert not covered&~mask
    residual=positions(mask&~covered)
    exceptions=[];witnesses=[]
    for n in residual:
        a=next((a for a in range(1,n//2+1) if flags[a] and flags[n-a]),None)
        if a is None:
            exceptions.append(n)
        else:
            witnesses.append([n,a,n-a])
            covered|=1<<n
    expected=mask
    for n in exceptions:
        expected&=~(1<<n)
    require_equal(covered,expected)
    output=covered.to_bytes(len(packed),'little')
    negatives=[]
    assert 2 not in exceptions
    assert rejected(lambda:require_equal(covered&~(1<<2),expected))
    negatives.append('deleted represented target2')
    first_exception=exceptions[0]
    assert rejected(lambda:require_equal(covered|(1<<first_exception),expected))
    negatives.append(f'injected exception{first_exception}')
    assert rejected(lambda:require_equal(covered|(1<<exclude_modulus),expected))
    negatives.append(f'injected locally forbidden target{exclude_modulus}')
    record=dict(excluded_target_modulus=exclude_modulus,short_bound=short_bound,short_count=len(anchors),
                value_count=bits.bit_count(),locally_admissible_target_count=mask.bit_count(),
                locally_excluded_target_count=bound//exclude_modulus,
                short_sweep_uncovered_local_targets=residual,residual_witnesses=witnesses,
                residual_witness_count=len(witnesses),all_residuals_resolved=True,
                exceptions=exceptions,exception_count=len(exceptions),largest_exception=max(exceptions,default=None),
                exhaustive_nonrepresentation_checks=True,represented_target_count=covered.bit_count(),
                forbidden_targets_covered=0,exact_target_mask_equal=True,negative_controls=negatives,
                value_flags_sha256=sha256(flags).hexdigest(),local_mask_sha256=sha256(local).hexdigest(),
                coverage_sha256=sha256(output).hexdigest())
    return record,output



def finite_coverage():
    short=5000
    direct=coordinate_flags(B)
    support,units=semigroup_flags(B)
    require_equal(direct,support)
    for n in range(1,10001):
        assert bool(direct[n])==trial_member(n)
        assert bool(units[n])==trial_member(n,unit=True)
    square=bytearray(1001)
    ycap=isqrt(4000//7)+1
    xcap=isqrt(1000)+ycap+1
    for x in range(-xcap,xcap+1):
        for y in range(-ycap,ycap+1):
            n=x*x+x*y+2*y*y
            if 1<=n<=1000 and gcd(x,y)==1:square[n]=1
    require_equal(square,direct[:1001])
    assert direct[4] and direct[7] and not direct[9] and not direct[49]
    p_record,p_mask=sumscan(direct,49,short)
    c_record,c_mask=sumscan(units,7,short)
    assert p_record['exceptions']==[7,77,203]
    assert c_record['exceptions']==[11,29]
    assert p_record['residual_witness_count']==c_record['residual_witness_count']==0
    assert p_record['short_count']==1050 and c_record['short_count']==903
    assert p_record['value_count']==2490663 and c_record['value_count']==2161776
    assert p_record['locally_admissible_target_count']==16434823
    assert c_record['locally_admissible_target_count']==14380470
    assert p_record['represented_target_count']==16434820
    assert c_record['represented_target_count']==14380468
    for n in range(1,B//7+1):
        nn=7*n
        assert ((p_mask[nn>>3]>>(nn&7))&1)==((c_mask[n>>3]>>(n&7))&1)
    for a,b in ((8,113),(2,317),(32,809)):
        assert units[a] and units[b]
    assert [(a+b) for a,b in ((8,113),(2,317),(32,809))]==[121,319,841]
    assert direct[4] and direct[7] and direct[22]
    for n,a,b in ((11,4,7),(29,7,22)):
        assert a+b==n and direct[a] and direct[b]
    bad=bytearray(direct);bad[9]=1
    assert rejected(lambda:require_equal(bad,support))
    bad=bytearray(direct);bad[4]=0
    assert rejected(lambda:require_equal(bad,support))
    bad=bytearray(direct);bad[49]=1
    assert rejected(lambda:require_equal(bad,support))
    return dict(bound=B,short_bound=short,
                complete_coordinate_vs_character_factor_flags_bound=B,
                independent_trial_factor_bound=10000,independent_square_box_bound=1000,
                primitive_P=p_record,unramified_C=c_record,
                ramified_value_identity='P intersect7Z=7C',
                ramified_sum_identity='(P+P) intersect7Z=7(C+C)',
                ramified_sum_mask_quotient_bound=B//7,
                exceptional_descent_product_witnesses=[[121,8,113],[319,2,317],[841,32,809]],
                unit_exception_P_witnesses=[[11,4,7],[29,7,22]],
                norm_corruption_controls=['inert norm9 inserted','split square4 deleted','ramified square49 inserted'])


def uint32_hash(values):
    encoded=array('I',values)
    assert encoded.itemsize==4
    if sys.byteorder!='little':encoded.byteswap()
    return sha256(encoded.tobytes()).hexdigest()


def validate_ordinary_counts(actual,character_sums):
    assert actual[0]==1
    for n in range(1,len(actual)):
        assert actual[n]==2*character_sums[n],n


def point_tables():
    ordinary=array('I',[0])*(X0+1)
    for y in range(isqrt(4*X0//7)+1):
        h=isqrt(4*X0-7*y*y)
        for x in range(-((h+y)//2),(h-y)//2+1):
            ordinary[x*x+x*y+2*y*y]+=1 if y==0 else 2
    character_sums=array('i',[0])*(X0+1)
    for d in range(1,X0+1):
        value=0 if d%7==0 else (1 if d%7 in QR7 else -1)
        if value:
            for n in range(d,X0+1,d):character_sums[n]+=value
    validate_ordinary_counts(ordinary,character_sums)
    primitive4=sum(x*x+x*y+2*y*y==4 and gcd(x,y)==1
                   for x in range(-4,5) for y in range(-4,5))
    assert ordinary[4]==6 and primitive4==4
    wrong=array('I',ordinary);wrong[4]=primitive4
    assert rejected(lambda:validate_ordinary_counts(wrong,character_sums))
    F1=array('I',[0])*(X0+1)
    current={1:0,2:0,4:0}
    cumulative=0
    for n in range(1,X0+1):
        cumulative+=ordinary[n]
        assert cumulative<=3*n
        if n%7 in current:current[n%7]+=ordinary[n]
        F1[n]=max(current.values())
    assert list(F1[1:8])==[2,4,4,6,6,6,6]
    assert all(F1[n]<=2*n for n in range(1,X0+1))
    limit=1024
    sigma=[0]*(limit+1)
    for d in range(1,limit+1):
        for n in range(d,limit+1,d):sigma[n]+=d
    theta=[sum(ordinary[j]*ordinary[n-j] for j in range(n+1)) for n in range(limit+1)]
    for n in range(1,limit+1):
        assert theta[n]==4*(sigma[n]-7*(sigma[n//7] if n%7==0 else 0))
    vectors=[]
    for x in range(-6,7):
        for y in range(-6,7):
            n=x*x+x*y+2*y*y
            if 1<=n<=18 and n%7:vectors.append((x,y,n))
    bad18=sum(a[2]+b[2]==18 and
              ((a[0]%3==a[1]%3==0) or (b[0]%3==b[1]%3==0))
              for a in vectors for b in vectors)
    bad3=sum(a[2]+b[2]==3 and
             ((a[0]%3==a[1]%3==0) or (b[0]%3==b[1]%3==0))
             for a in vectors for b in vectors)
    assert bad18==4 and bad3==0
    ratio=F(bad18,theta[18])
    assert ratio==F(1,39) and ratio<=F(1,9)
    def validate_cost(cost):assert ratio<=cost
    assert rejected(lambda:validate_cost(F(0)))
    record=dict(point_bound=X0,ordinary_coefficient_regression_bound=X0,
                ordinary_counts_sha256=uint32_hash(ordinary),F1_sha256=uint32_hash(F1),
                hash_encoding='unsigned32 little endian',F1_first_seven=list(F1[1:8]),
                maximal_order_norm_shell_regression_bound=limit,
                norm_shell_formula='4*(sigma(n)-7*sigma(n/7))',
                norm_shell_coefficients_through32=theta[:33],
                negative_controls=dict(primitive_point_array_rejected=True,
                    ordinary_rQ4=6,primitive_rQ4=primitive4,
                    forgotten_target_square_cost_rejected=True,
                    inert_prime=3,target=18,actual_bad_good7_count=bad18,
                    actual_spoiled_ratio=str(ratio),safe_target_square_cost='1/9',
                    valuation_one_target3_bad_count=bad3))
    return F1,record


def make_e_frontiers(small):
    caps={p:(4 if p==3 else 3 if p==5 else 1) for p in small}
    states={key:[] for key in product((0,1,2),repeat=2)}
    for exponents in product(*(range(caps[p]+1) for p in small)):
        n=prod(p**a for p,a in zip(small,exponents))
        value=prod(F((a+1)**6*p**(5*a),sigma_power(p,a)**6)
                   for p,a in zip(small,exponents))
        states[min(exponents[0],2),min(exponents[1],2)].append((n,value))
    assert sum(len(rows) for rows in states.values())==5120
    def frontier(rows):
        out=[];best=F(0)
        for n,value in sorted(rows):
            if value>best:out.append((n,value));best=value
        return out
    global_front=frontier([row for rows in states.values() for row in rows])
    result={key:frontier(rows) for key,rows in states.items()}
    return global_front,result,caps


def arithmetic_context():
    flags,prime_list=primes(R)
    assert [p for p in prime_list if p<=1000]==[p for p in range(2,1001)
        if all(p%d for d in range(2,isqrt(p)+1))]
    inert=[p for p in prime_list if p%7 in (3,5,6)]
    split=[p for p in prime_list if p%7 in QR7]
    small=[p for p in inert if p<64]
    assert small==[3,5,13,17,19,31,41,47,59,61]
    assert F(10,9)**6<2 and F(6,5)**6<3 and F(5,4)**6<5
    assert F(4,3)**6<11 and F(3,2)**6<13
    assert all(a+1<=2**a for a in range(1,101))
    CM6=prod(max(F((a+1)**6,p**a) for a in range(9)) for p in split if p<64)
    CE6=prod(max(F((a+1)**6*p**(5*a),sigma_power(p,a)**6) for a in range(5)) for p in small)
    assert CM6<8**6 and CE6<F(91,20)**6
    all_rows=[]
    for p in prime_list:
        if p>=64:break
        if p==7:continue
        cap={2:8,3:4,5:3,11:2}.get(p,1)
        values=[F((a+1)**6*p**(5*a),sigma_power(p,a)**6) for a in range(cap+1)]
        maximum=max(values)
        all_rows.append([p,cap,values.index(maximum),maximum.numerator,maximum.denominator])
    CEall6=prod(F(row[3],row[4]) for row in all_rows)
    assert CEall6<15**6
    assert rejected(lambda:require_equal(CEall6<14**6,True))
    Eglobal,Ecat,caps=make_e_frontiers(small)
    require_equal(Eglobal[-1][1],CE6)
    assert rejected(lambda:require_equal(Eglobal[-1][1],2*CE6))
    Mfront,Mstates,prefix_count=m_frontier(split,1<<54)
    assert Mstates==2556 and prod(split[:prefix_count])>1<<54
    Mkeys=[n for n,v in Mfront];Ekeys=[n for n,v in Eglobal]
    running=1
    for n in range(1,10001):
        rem=n;value=1;inert_only=True;sig=1;tau=1;p=2
        while p*p<=rem:
            a=0
            while rem%p==0:rem//=p;a+=1
            if a:
                tau*=a+1;sig*=sigma_power(p,a)
                if p%7 in QR7:value*=a+1
                if p%7 not in (3,5,6):inert_only=False
            p=3 if p==2 else p+2
        if rem>1:
            tau*=2;sig*=rem+1
            if rem%7 in QR7:value*=2
            if rem%7 not in (3,5,6):inert_only=False
        running=max(running,value)
        assert running==Mfront[bisect_right(Mkeys,n)-1][1]
        if inert_only:
            assert F(tau**6*n**5,sig**6)<=Eglobal[bisect_right(Ekeys,n)-1][1]
    pre2=[0]
    for p in inert:pre2.append(pre2[-1]+(SCALE+p*p-1)//(p*p))
    starts=[]
    for residue in (3,5,13):
        offset=(residue-R)%14
        starts.append(R+(offset if offset else 14))
    far=sum((F(1,a*a)+F(1,14*a) for a in starts),F(0))
    Bmu=F(sum((2*SCALE+p*(p-1)-1)//(p*(p-1)) for p in inert),SCALE)+2*F(R+1,R)*far
    assert Bmu<F(47,100)
    assert max(F(sum(1 for a in (3,5,13) if a<=r))-F(3*r,14) for r in range(14))==F(13,14)
    record=dict(prime_limit=R,prime_flags_sha256=sha256(flags).hexdigest(),
                fixed_point_scale=SCALE,inert_prime_count=len(inert),B_upper=str(Bmu),
                reciprocal_square_progression_tail=str(far),progression_tail_starts=starts,
                M_sixth_power=str(CM6),core_E_sixth_power=str(CE6),
                all_7_unit_E_sixth_power=str(CEall6),all_prime_E_maxima_columns=['p','cap','max_exponent','numerator','denominator'],
                all_prime_E_maxima=all_rows,M_bound_constant=8,core_E_bound_constant='91/20',all_E_bound_constant=15,
                core_small_prime_caps=caps,M_frontier_upper=1<<54,M_frontier_states=Mstates,
                split_prefix_completeness_length=prefix_count,split_prefix_product=prod(split[:prefix_count]),
                M_frontier=Mfront,E_global_frontier=[[n,str(v)] for n,v in Eglobal],
                E_category_frontiers={str(k):[[n,str(v)] for n,v in rows] for k,rows in Ecat.items()},
                independent_frontier_regression_bound=10000,
                negative_controls=['corrupted Euler product rejected','all-target E bound14 rejected'])
    return dict(inert=inert,split=split,pre2=pre2,far=far,Bmu=Bmu,
                Mfront=Mfront,Mkeys=Mkeys,Eglobal=Eglobal,Ecat=Ecat,record=record)



def validate_intervals(intervals,start,end):
    assert intervals and intervals[0][0]==start and intervals[-1][1]==end-1
    expected=start
    for lower,upper in intervals:
        assert lower==expected and lower<=upper
        expected=upper+1
    assert expected==end


def fine_profile_key(row):
    return (row['chi'],*row['v3_v5'],row['boundary_tau'] or 0)


def validate_fine_grid(rows,intervals):
    expected={(sign,a,b,0) for sign in (-1,1) for a,b in product((0,1,2),repeat=2)}
    expected.add((1,-1,-1,4))
    groups={pair:[] for pair in intervals}
    assert len(rows)==len(intervals)*19
    for row in rows:
        key=(row['L'],row['U'])
        assert key in groups
        groups[key].append(row)
        assert row['passed'] and F(row['margin'])>F(1,64)
    for group in groups.values():
        assert len(group)==19 and {fine_profile_key(row) for row in group}==expected


def fine_certificates(context,F1):
    inert=context['inert'];pre2=context['pre2'];Bmu=context['Bmu']
    Mfront=context['Mfront'];Mkeys=context['Mkeys'];Efront=context['Ecat']
    Ekeys={key:[n for n,v in rows] for key,rows in Efront.items()}
    cache=[None,None]
    def qr_suffix(upper):
        if cache[0]==upper:return cache[1]
        assert isqrt(upper)<R
        endpoint=bisect_right(inert,isqrt(upper))
        minimum=max(2,bisect_right(inert,isqrt(upper//X0)))
        suffix=[0]*(endpoint+1)
        for i in range(endpoint-1,minimum-1,-1):
            index=upper//inert[i]**2
            assert 0<=index<=X0
            suffix[i]=suffix[i+1]+F1[index]
        cache[0]=upper;cache[1]=(endpoint,minimum,suffix)
        return cache[1]

    def best_case(lower,upper,sign,category=None,boundary_tau=None):
        M=Mfront[bisect_right(Mkeys,upper)-1][1]
        if boundary_tau is not None:
            E=F(boundary_tau,isqrt(lower));D=F(1);present=0;correction=F(0)
            cat=(-1,-1)
        else:
            cat=category
            front=Efront[cat]
            sixth=front[bisect_right(Ekeys[cat],upper)-1][1]
            E=rational_root_upper(sixth)/iroot(lower,3)
            D=F(1);present=0;correction=F(0)
            for p,a in zip((3,5),cat):
                if a:
                    present+=1
                    correction+=(F(1,2) if sign==1 else 1)*F(2,p*(p-1))
                    if a==2:correction-=F(1,p*p)
                    D*=F(p**a,sigma_power(p,a))
        gamma=F(1,2) if sign==1 else F(1)
        baseline=gamma*(1-Bmu)+correction
        endpoint=bisect_right(inert,isqrt(upper))
        MD=M*D
        minimum=2
        if sign==1:
            endpoint,minimum,qr_tail=qr_suffix(upper)
            step=F(7,4)*E
            def saving(i):return MD*F(F1[upper//inert[i]**2],lower)
        else:
            step=E
            def saving(i):return 3*MD*F((SCALE+inert[i]**2-1)//inert[i]**2,SCALE)
        # The savings are nonincreasing with p, so this locates the
        # exact minimizing prime step among p<=sqrt(upper).
        lo,hi=minimum,endpoint
        while lo<hi:
            mid=(lo+hi)//2
            if saving(mid)>step:lo=mid+1
            else:hi=mid
        count=max(minimum,lo)
        assert present<=count<=endpoint and isqrt(upper)<R
        T=F(pre2[endpoint]-pre2[count],SCALE)
        if sign==1:
            large=MD*F(qr_tail[count],lower)
            small_error=E*(F(1,2)+F(7,4)*(count-present))
        else:
            large=3*MD*T
            small_error=E*(count-present)
        left=small_error+large
        return dict(L=lower,U=upper,chi=sign,v3_v5=list(cat),boundary_tau=boundary_tau,
                    M=M,E=str(E),sigma_discount=str(D),P=inert[count-1],pi=count,
                    prime_tail_count=endpoint-count,
                    small_error=str(small_error),large_error=str(large),left=str(left),baseline=str(baseline),
                    margin=str(baseline-left),passed=left<baseline)

    cuts={1<<24,1<<37}
    for k in range(24,37):
        step=1<<(k-4)
        cuts.update((1<<k)+j*step for j in range(17))
    cuts.update(n for n,v in Mfront if 1<<24<n<1<<37)
    for front in Efront.values():cuts.update(n for n,v in front if 1<<24<n<1<<37)
    cuts=sorted(cuts)
    intervals=[(lower,end-1) for lower,end in zip(cuts,cuts[1:])]
    assert len(intervals)==257
    validate_intervals(intervals,1<<24,1<<37)
    rows=[]
    for lower,upper in intervals:
        for sign in (-1,1):
            for category in product((0,1,2),repeat=2):
                rows.append(best_case(lower,upper,sign,category=category))
        rows.append(best_case(lower,upper,1,boundary_tau=4))
    assert len(rows)==4883
    validate_fine_grid(rows,intervals)
    wrong_intervals=list(intervals)
    wrong_intervals[1]=(wrong_intervals[1][0]+1,wrong_intervals[1][1])
    assert rejected(lambda:validate_intervals(wrong_intervals,1<<24,1<<37))
    wrong_rows=list(rows)
    boundary=next(i for i,row in enumerate(rows) if row['boundary_tau']==4)
    wrong_rows[boundary]=dict(rows[0])
    assert rejected(lambda:validate_fine_grid(wrong_rows,intervals))
    profile_keys=sorted({fine_profile_key(row) for row in rows})
    profiles=[]
    for key in profile_keys:
        example=next(row for row in rows if fine_profile_key(row)==key)
        profiles.append(dict(key=list(key),baseline=example['baseline'],sigma_discount=example['sigma_discount']))
    profile_ids={key:i for i,key in enumerate(profile_keys)}
    interval_ids={pair:i for i,pair in enumerate(intervals)}
    compact=[]
    for row in rows:
        e=F(row['E']);left=F(row['left'])
        compact.append([interval_ids[row['L'],row['U']],profile_ids[fine_profile_key(row)],
                        e.numerator,e.denominator,row['P'],row['pi'],row['prime_tail_count'],
                        left.numerator,left.denominator])
    worst=min(rows,key=lambda row:F(row['margin']))
    return dict(domain=[1<<24,(1<<37)-1],interval_count=len(intervals),case_count=len(rows),
                intervals=intervals,profiles=profiles,
                case_columns=['interval_id','profile_id','E_numerator','E_denominator','P_last_inert',
                              'pi_inert','large_prime_count','left_numerator','left_denominator'],
                cases=compact,uniform_margin='1/64',worst_case=worst,
                every_point_tail_argument_checked=True,every_sqrt_upper_below_prime_limit=True,
                negative_controls=['interval gap rejected','boundary omitted at unchanged row count rejected'])


def higher_bin_certificates(context):
    inert=context['inert'];pre2=context['pre2'];far=context['far']
    Bmu=context['Bmu'];Mfront=context['Mfront'];Mkeys=context['Mkeys']
    Efront=context['Eglobal'];Ekeys=[n for n,v in Efront]
    baseline=(1-Bmu)/2
    rows=[];worst=None
    for k in range(37,54):
        lower=1<<k;upper=(1<<(k+1))-1
        cube=iroot(lower,3)
        assert 4*cube<=isqrt(lower) # explicit tau<=4 / split-prime boundary domination
        C=rational_root_upper(Efront[bisect_right(Ekeys,upper)-1][1])
        assert C>=1
        E=C/cube
        M=Mfront[bisect_right(Mkeys,upper)-1][1]
        lo,hi=0,len(inert)
        while lo<hi:
            mid=(lo+hi)//2
            saved=3*M*F((SCALE+inert[mid]**2-1)//inert[mid]**2,SCALE)
            if saved>F(7,4)*E:lo=mid+1
            else:hi=mid
        count=lo
        assert 0<count<len(inert) and inert[count-1]<R
        tail=F(pre2[-1]-pre2[count],SCALE)+far
        left=E*(F(1,2)+F(7,4)*count)+3*M*tail
        margin=baseline-left
        assert margin>F(1,64)
        rows.append([k,M,cube,E.numerator,E.denominator,inert[count-1],count,
                     tail.numerator,tail.denominator,left.numerator,left.denominator])
        if worst is None or margin<worst[0]:worst=(margin,k,left)
    assert len(rows)==17 and [row[0] for row in rows]==list(range(37,54))
    return dict(domain=[1<<37,(1<<54)-1],bin_count=17,uniform_margin='1/64',
                columns=['k','M','floor_cuberoot_lower','E_numerator','E_denominator','P','pi_inert',
                         'tail_numerator','tail_denominator','left_numerator','left_denominator'],
                bins=rows,baseline=str(baseline),worst_bin=worst[1],worst_margin=str(worst[0]),
                tau_le4_boundary_dominance_explicitly_checked=True)


def coarse_tail_and_counts(context):
    assert context['Bmu']<F(47,100)
    assert F(91,20)>=4 # all tau<=4 targets satisfy the same coarse E bound
    assert 45**12<2**66
    leading=F(21,32)*F(91,20)+F(144,49)
    second=F(17,8)*F(91,20)+F(1152,49)
    coarse_left=leading/F(45,2)+second/(1<<18)
    coarse_margin=F(53,200)-coarse_left
    assert coarse_margin>F(1,1024)
    assert coarse_margin==F(9574319,6165626880)
    assert 45**2<2**11
    all_left=(F(3*15,8)+F(36,7))/45+(F(25*15,8)+72)/(1<<22)
    all_margin=F(53,200)-all_left
    assert all_left==F(281051653,1174405120)
    assert all_margin==F(150828519,5872025600)>F(1,40)
    assert 7*(1<<66)<1<<69 and 7**2<4**3
    assert 4*8**2==256 and 40*256//4==2560
    assert 2560*4<20000
    return dict(inert_core_and_tau_le4_onset=1<<54,
                coarse_error_leading=str(leading),coarse_error_second=str(second),
                coarse_error_bound=str(coarse_left),coarse_margin=str(coarse_margin),
                coarse_margin_exceeds='1/1024',
                all_7_unit_target_onset=1<<66,all_target_error_bound=str(all_left),
                all_target_normalized_margin=str(all_margin),all_target_margin_exceeds='1/40',
                unit_value_pair_lower_bound='n^(2/3)/2560',
                primitive_P_global_count_onset=1<<69,
                primitive_P_global_ordered_value_pair_lower_bound='N^(2/3)/20000',
                far_tail_has_no_upper_endpoint=True,prime_boundary_far_tail_covered=True)



def small_joint_group_checks():
    def projective(vector,p):
        inverse=pow(next(x for x in vector if x%p),-1,p)
        return tuple(x*inverse%p for x in vector)
    def r_mul(a,b):
        x,y=a;z,t=b
        return (x*z-2*y*t,x*t+y*z+y*t)
    def r_conj(a):return (a[0]+a[1],-a[1])
    def q_product(a,b,p):
        alpha,beta=a[:2],a[2:];gamma,delta=b[:2],b[2:]
        x=r_mul(alpha,gamma);y=r_mul(beta,r_conj(delta))
        z=r_mul(alpha,delta);t=r_mul(beta,r_conj(gamma))
        return projective((x[0]-y[0],x[1]-y[1],z[0]+t[0],z[1]+t[1]),p)
    def norm(a):
        x,y,z,t=a
        return x*x+x*y+2*y*y+z*z+z*t+2*t*t
    def character(value,p):
        result=pow(value%p,(p-1)//2,p)
        return 0 if result==0 else (1 if result==1 else -1)
    def f49_product(a,b):return ((a[0]*b[0]-a[1]*b[1])%7,(a[0]*b[1]+a[1]*b[0])%7)
    ram=[];value=(1,0)
    for k in range(8):
        ram.append(projective(value,7))
        value=f49_product(value,(1,2))
    assert len(set(ram))==8 and projective(value,7)==(1,0) and ram[4]==(0,1)
    assert all(character(a*a+b*b,7)==(-1)**k for k,(a,b) in enumerate(ram))
    gamma_fraction={sign:F(sum(k not in (0,4) and (-1)**k==sign for k in range(8)),
                           sum((-1)**k==sign for k in range(8))) for sign in (-1,1)}
    assert gamma_fraction=={-1:F(1),1:F(1,2)}
    delta=[1,0,0,0]
    ram_means=[F(sum(delta[k]*character7 for k,character7 in enumerate(values)),4)
               for values in ([1,1,1,1],[1,-1,1,-1])]
    assert ram_means==[F(1,4),F(1,4)]
    ram_perp=F(1,4)-sum(x*x for x in ram_means)
    assert ram_perp==F(1,8) and ram_perp/F(1,4)==F(1,2)
    rows=[]
    for p in (3,5):
        G=sorted({projective(a,p) for a in product(range(p),repeat=4) if norm(a)%p})
        K={a for a in G if a[:2]==(0,0) or a[2:]==(0,0)}
        j=(0,0,1,0);one=(1,0,0,0)
        assert len(G)==p*(p*p-1) and len(K)==2*(p+1)
        assert j in K and q_product(j,j,p)==one
        assert all(q_product(a,b,p) in K for a in K for b in K)
        mu=F(len(K),len(G))
        assert mu==F(2,p*(p-1))
        assert sum(character(norm(a),p)==1 for a in K)==len(K)//2
        # Quotient by the diagonal global unit (ramified exponent4,j).
        def gamma_action(pair):
            k,a=pair
            return ((k+4)%8,q_product(j,a,p))
        quotient={min((k,a),gamma_action((k,a))) for k in range(8) for a in G}
        assert len(quotient)==4*len(G)
        for pair in quotient:assert gamma_action(gamma_action(pair))==pair
        def indicators(pair):
            k,a=pair;inside=k in (0,4) and a in K
            outside=k not in (0,4) and a in K
            return int(inside),int(outside)
        assert all(indicators(pair)==indicators(gamma_action(pair)) for pair in quotient)
        chars=[(0,0),(1,0),(0,1),(1,1)]
        def char(pair,key):
            k,a=pair
            return ((-1)**k)**key[0]*character(norm(a),p)**key[1]
        for a in chars:
            for b in chars:
                assert F(sum(char(pair,a)*char(pair,b) for pair in quotient),len(quotient))==int(a==b)
        projections=[]
        for side in (0,1):
            projections.append([F(sum(indicators(pair)[side]*char(pair,key) for pair in quotient),len(quotient)) for key in chars])
        a=mu/4
        assert projections[0]==[a,a,0,0] and projections[1]==[3*a,-a,0,0]
        u_perp=a-sum(x*x for x in projections[0])
        v_perp=3*a-sum(x*x for x in projections[1])
        assert u_perp==a-2*a*a and v_perp==3*a-10*a*a
        C2=u_perp*v_perp/(a*a)
        assert C2==(1-mu/2)*(3-5*mu/2)<3<F(7,4)**2
        for s7,sq in product((-1,1),repeat=2):
            eigenvalues=[1,s7,sq,s7*sq]
            main=sum(x*y*z for x,y,z in zip(projections[0],projections[1],eigenvalues))/a
            assert main==mu*F(3-s7,4)
        assert rejected(lambda:require_equal(mu/2,mu))
        assert rejected(lambda:require_equal(a-a*a,u_perp))
        rows.append(dict(q=p,projective_group_size=len(G),normalizer_size=len(K),
                         joint_quotient_size=len(quotient),mu=str(mu),
                         u_projections=[str(x) for x in projections[0]],
                         v_projections=[str(x) for x in projections[1]],
                         u_complement_norm_square=str(u_perp),v_complement_norm_square=str(v_perp),
                         Cq_square=str(C2),Cq_upper='7/4'))
    return dict(scope='Finite groups, character projections and normalization only; the spectral transfer remains a written input.',
                ramified_cyclic_order=8,ramified_projective_generator=[1,2],
                ramified_global_unit_exponent=4,ramified_delta_complement_norm_square='1/8',
                ramified_return_error_coefficient='1/2',cases=rows,
                negative_controls=['missing normalizer half rejected','omitted ramified norm character rejected'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus7_checks.json'))
    args=parser.parse_args()
    coverage=finite_coverage()
    point_array,points=point_tables()
    context=arithmetic_context()
    fine=fine_certificates(context,point_array)
    higher=higher_bin_certificates(context)
    tail=coarse_tail_and_counts(context)
    joint=small_joint_group_checks()
    assert coverage['bound']==1<<24
    assert fine['domain']==[1<<24,(1<<37)-1]
    assert higher['domain']==[1<<37,(1<<54)-1]
    assert tail['inert_core_and_tau_le4_onset']==1<<54
    result=dict(checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Complete fixed finite checks supporting the exact -7 classification and explicit count onset. The source-backed class arithmetic, joint spectral bounds and descent proof are separate inputs; this is not formal verification or external review.',
                proof_domains_fixed=True,
                identities=dict(C_plus_C='all n>=2 with7 not dividing n, except11,29',
                                P_plus_P='all N>=2 with49 not dividing N, except7,77,203'),
                finite_coverage=coverage,ordinary_point_tables=points,
                arithmetic_envelopes=context['record'],fine_interval_certificates=fine,
                higher_interval_certificates=higher,unrestricted_tail_and_counts=tail,
                finite_joint_normalization=joint)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(source_sha256=result['checker_sha256'],finite_bound=B,
                         C_exceptions=coverage['unramified_C']['exceptions'],
                         P_exceptions=coverage['primitive_P']['exceptions'],
                         point_bound=X0,fine_intervals=fine['interval_count'],fine_cases=fine['case_count'],
                         higher_bins=higher['bin_count'],far_tail_onset=tail['inert_core_and_tau_le4_onset'],
                         all_target_count_onset=tail['primitive_P_global_count_onset'],
                         output=str(args.output),status='pass'),indent=2))


if __name__=='__main__':
    main()
