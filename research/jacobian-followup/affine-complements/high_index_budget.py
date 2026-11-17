"""PF112 necessary actual-Keller filters and exact seventeen-sheet controls."""
import math
import sympy as S
import full_peripheral_linking as F


def universal_index_cap(d):
    if not isinstance(d,int) or isinstance(d,bool) or d<5:
        raise ValueError('the new all-degree cap is stated for d>=5')
    return d-3


def high_index_reserve(d,e,hidden_generator_count=0):
    if any(not isinstance(x,int) or isinstance(x,bool) for x in (d,e,hidden_generator_count)) or d<2 or e<2 or hidden_generator_count<0:
        raise ValueError('positive field degree/index and nonnegative hidden generator count required')
    active=2*e>d
    minimum_others=max(1,hidden_generator_count) if active else hidden_generator_count
    reserve=1 if active else 0
    return dict(strict_high_index=active,residue_one_forced=active,
                positive_own_excess_forced=active,
                distinct_nonproper_image_forced=active,
                minimum_other_primes=minimum_others,minimum_excess_reserve=reserve,
                necessary_budget_passes=(e+minimum_others+reserve<=d-1),
                necessary_only=True,realization_claimed=False)


def next_strip_profile(d):
    if not isinstance(d,int) or isinstance(d,bool) or d<7:
        raise ValueError('e=d-3 is strictly high only for d>=7')
    return dict(field_degree=d,ramified_index=d-3,ramified_residue_degree=1,
                omitted_prime_count=2,other_index=1,other_residue_degree=1,
                own_excess=1,other_excess=0,exceptional_local_length=d-2,
                ramified_image_normalization_injective=True,
                ramified_nonimmersive_parameter_count=1,
                other_image_distinct=True,other_image_normalization_immersive=True,
                other_image_has_collision=True,source_smoothness_claimed=False,
                necessary_only=True,realization_claimed=False)


def same_image_mass_filter(d,primes,*,whole_nonproper_set_irreducible=False):
    """Known same-image prime sublists suffice; unknown primes are not removed."""
    primes=list(primes)
    if not isinstance(d,int) or isinstance(d,bool) or d<2 or not primes or any(
            not isinstance(e,int) or isinstance(e,bool) or not isinstance(f,int) or isinstance(f,bool) or e<1 or f<1
            for e,f in primes):
        raise ValueError('field degree and nonempty positive integral same-image prime data required')
    mass=sum(e*f for e,f in primes);b=sum(f==1 for e,f in primes)
    high=2*mass>d
    return dict(known_mass=mass,known_residue_one_prime_count=b,strict_high_mass=high,
                own_excess_units_forced=b if high else 0,
                outside_image_prime_forced=high,
                necessary_budget_passes=(mass+b+1<=d-1 if high else mass<=d-1) and
                    (not whole_nonproper_set_irreducible or 2*mass<=d),
                irreducible_total_mass_cap=d//2,
                irreducible_residue_normalized_excess_lower_bound=d-1-d//2,
                necessary_only=True,realization_claimed=False)


def gf16_mul(a,b):
    out=0
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a&16:a^=19
    return out


def gf16_power(a,n):
    out=1
    while n:
        if n&1:out=gf16_mul(out,a)
        a=gf16_mul(a,a);n>>=1
    return out


def gf16_inverse(a):
    if not a:raise ValueError('nonzero field element required')
    return gf16_power(a,14)


def projective_action(a,b,c,d):
    if gf16_mul(a,d)^gf16_mul(b,c)!=1:
        raise ValueError('determinant one over F16 required')
    out=[]
    for x in list(range(16))+[None]:
        num,den=(a,c) if x is None else (gf16_mul(a,x)^b,gf16_mul(c,x)^d)
        out.append(16 if not den else gf16_mul(num,gf16_inverse(den)))
    return tuple(out)


def seventeen_generators():
    t=2^gf16_inverse(2)^1
    return [projective_action(1,t,0,1),projective_action(0,1,1,1)]


def valid_cables(cables):
    return all(isinstance(r,int) and not isinstance(r,bool) and isinstance(k,int) and not isinstance(k,bool) and
               r>1 and k>1 and r%15==1 and k%15==2 and math.gcd(r,k)==1 for r,k in cables)


def seventeen_cover(cables=(),reverse=False):
    cables=list(cables)
    if not valid_cables(cables):raise ValueError('coprime ordinary r1/k2 mod15 cable sector required')
    rho=seventeen_generators();base=F.peripheral_cover(2,3,[],rho)
    mu=base['meridian_monodromy'];lam=base['longitude_monodromy']
    root=tuple(lam[mu[mu[i]]] for i in range(17))
    return F.peripheral_cover(2,3,cables,rho+[root]*len(cables),reverse)


def seventeen_transfer(cables):
    cables=list(cables)
    if not valid_cables(cables):raise ValueError('proved cable sector required')
    R=math.prod(r for r,k in cables);base=S.Rational(6,17)*S.ones(3)+S.diag(0,-6,-6)
    matrix=R*R*base+S.diag(S.Rational(6,15)*(1-R*R),0,0)
    q=17//math.gcd(17,R);hidden=51//q
    return dict(winding_product=R,matrix=[[str(x) for x in row] for row in matrix.tolist()],
                closed_H1=dict(free=0,torsion=[51]),core_orders=[q]*3,
                exterior_H1=dict(free=3,torsion=[hidden]),
                hidden_quotient=dict(free=0,torsion=[hidden]),
                exact_row_defects=[1]*3,all_degree_index_cap=14,
                known_ramified_index=15,actual_Keller_realization_excluded=True,
                arbitrary_uncolored_representation_excluded=False)
