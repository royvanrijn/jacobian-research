"""PF113 necessary marked-Keller local capacity and full-branch filters."""
import math
import sympy as S


def local_partition_gate(d,r,blocks):
    blocks=list(blocks)
    if any(not isinstance(n,int) or isinstance(n,bool) for n in [d,r]+blocks) or d<2 or r<1 or r>=d or not blocks or any(n<1 for n in blocks) or sum(blocks)!=d:
        raise ValueError('d>1, 1<=r<d and a positive integral partition of d required')
    cap=d-r
    return dict(local_cluster_cap=cap,normalization_fibre_count_lower_bound=math.ceil(d/cap),
                necessary_capacity_passes=max(blocks)<=cap,
                local_action_could_be_transitive_under_Keller_hypotheses=False,
                source_restriction_singleton_fibre_excluded=False,
                source_surjectivity_claimed=False,realization_claimed=False)


def positive_weighted_branch_gate(poly,weights,*,full_branch_union_certified):
    if not isinstance(poly,S.Poly) or len(poly.gens)!=2 or poly.total_degree()<1:
        raise ValueError('nonconstant polynomial in two named variables required')
    if len(weights)!=2 or any(not isinstance(n,int) or isinstance(n,bool) or n<=0 for n in weights):
        raise ValueError('strictly positive integer TARGET weights required')
    degrees={sum(a*b for a,b in zip(monom,weights)) for monom,coefficient in poly.terms() if coefficient!=0}
    homogeneous=len(degrees)==1 and next(iter(degrees))>0
    return dict(strictly_positive_target_weights=True,weighted_homogeneous=homogeneous,
                weighted_degrees=sorted(degrees),full_branch_union_certified=full_branch_union_certified is True,
                actual_Keller_cover_excluded=homogeneous and full_branch_union_certified is True,
                arbitrary_finite_cover_excluded=False,
                source_action_lift_assumed=False,whole_marked_Keller_plane_required=True)


def universal_index_cap(d):
    if not isinstance(d,int) or isinstance(d,bool) or d<7:
        raise ValueError('the strengthened universal index cap is stated for d>=7')
    return d-4


def next_strip_profile(d):
    if not isinstance(d,int) or isinstance(d,bool) or d<9:
        raise ValueError('e=d-4 is strictly high only for d>=9')
    return dict(field_degree=d,omitted_prime_count=2,ramified_prime_count=2,unramified_prime_count=0,
                first_index=d-4,second_index=2,first_residue=1,second_residue=1,
                first_excess=1,second_excess=0,first_special_local_length=d-3,
                distinct_ramified_images=True,first_image_injective=True,
                first_nonimmersive_parameter_count=1,second_image_immersive=True,
                second_image_has_collision=True,source_smoothness_assumed=False,
                necessary_only=True,realization_claimed=False)


def high_index_reserve(d,e,hidden_generator_count=0):
    if any(not isinstance(n,int) or isinstance(n,bool) for n in (d,e,hidden_generator_count)) or d<2 or e<2 or hidden_generator_count<0:
        raise ValueError('positive degree/index and nonnegative unseen generator count required')
    active=2*e>d;others=max(1,hidden_generator_count) if active else hidden_generator_count
    return dict(strict_high_index=active,minimum_other_primes=others,
                own_excess_reserve=1 if active else 0,
                additional_ramification_cost=1 if active else 0,
                distinct_ramified_image_forced=active,
                necessary_budget_passes=e+others+(2 if active else 0)<=d-1,
                realization_claimed=False)


def image_mass_gate(d,primes,*,irreducible_full_branch=False):
    primes=list(primes)
    if not isinstance(d,int) or isinstance(d,bool) or d<2 or not primes or any(
            not isinstance(e,int) or isinstance(e,bool) or not isinstance(f,int) or isinstance(f,bool) or e<1 or f<1
            for e,f in primes):
        raise ValueError('nonempty positive integral same-image data required')
    w=sum(e*f for e,f in primes);b=sum(f==1 for e,f in primes);high=2*w>d
    allowed=(w<=d-1)
    if high:allowed=not any(e==1 for e,f in primes) and w+b+2<=d-1 and not irreducible_full_branch
    return dict(known_same_image_mass=w,known_residue_one_count=b,strict_high_mass=high,
                minimum_outside_ramified_image_cost=2 if high else 0,
                necessary_budget_passes=allowed,irreducible_full_branch_mass_cap=d//2,
                no_total_mass_bound_outside_branch_claimed=True,
                necessary_only=True,realization_claimed=False)


def kummer_control(m):
    if not isinstance(m,int) or isinstance(m,bool) or m<2:
        raise ValueError('finite Kummer degree m>=2 required')
    z,t=S.symbols('z t')
    return dict(z=z,t=t,P=z**m,Q=t,Jacobian=m*z**(m-1),degree=m,
                generic_index=m,singleton_normalization_fibre_length=m,
                etale_open='C* times C',whole_plane_Keller_map=False)


def cubic_control():
    s,t=S.symbols('s t');u,v=S.symbols('u v')
    return dict(s=s,t=t,u=u,v=v,P=s,Q=t**3+s*t,Jacobian=s+3*t*t,
                branch=4*u**3+27*v**2,branch_parameter=[-3*t*t,-2*t**3],
                degree=3,generic_index=2,special_local_length=3,
                etale_open='C* times C',whole_plane_Keller_map=False)


def singular_quotient_control():
    A,B,C=S.symbols('A B C')
    return dict(A=A,B=B,C=C,relation=A*C-B*B,
                P=C*C+A*B,Q=C**3+S.Rational(3,2)*B**3+S.Rational(3,8)*A**3,
                degree=12,generic_index=9,source_singularity='A1',
                etale_open='C* times C',whole_plane_Keller_map=False)
