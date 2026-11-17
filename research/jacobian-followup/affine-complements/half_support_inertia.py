"""PF118 exact support bounds; existence/exclusion follows from the written theorem."""


def projective_bound(q,n):
    if not isinstance(q,int) or not isinstance(n,int) or q<2 or n<2:
        raise ValueError('integer q>=2,n>=2; group application requires q a prime power')
    degree=(q**n-1)//(q-1);hyperplane=(q**(n-1)-1)//(q-1)
    fixed=hyperplane+(q!=2)
    return dict(q=q,dimension=n,degree=degree,fixed_upper=int(fixed),
                support_lower=degree-fixed,half_support_excluded=2*(degree-fixed)>degree,
                prime_power_hypothesis_not_verified_here=True)


def affine_bound(q,n):
    if not isinstance(q,int) or not isinstance(n,int) or q<2 or n<1:
        raise ValueError('integer q>=2,n>=1; group application requires q a prime power')
    degree=q**n;fixed=q**(n-1)
    return dict(q=q,dimension=n,degree=degree,fixed_upper=fixed,
                support_lower=degree-fixed,half_support_excluded=2*(degree-fixed)>degree,
                prime_power_hypothesis_not_verified_here=True)


def support_obstruction(degree,minimum_moved):
    if degree<2 or not 1<=minimum_moved<=degree:raise ValueError('faithful nontrivial action required')
    return dict(degree=degree,minimum_moved=minimum_moved,
                excluded_for_plane_Keller=minimum_moved>degree//2,
                full_irreducible_half_support_necessary=2*minimum_moved==degree,
                higher_dimensional_exclusion_claimed=False)


def binary_affine_profile(n):
    if not isinstance(n,int) or n<2:raise ValueError('binary affine half-support case needs n>=2')
    d=2**n
    return dict(degree=d,full_branch_components=1,branch_inertia_support=d//2,
                inertia_cycles={'2':d//4,'1':d//2},ramified_prime_count=d//4,
                all_ramified_indices=2,all_residue_degrees=1,normalization_smooth=True,
                unramified_primes_plus_normalized_excess=d//2-1,
                branch_normalization_max_fibre=2,branch_double_values_exist=True,
                source_fibre_at_double_values=0,unramified_prime_over_branch_allowed=False,
                actual_plane_realization_constructed=False,binary_affine_family_excluded=False)


def cycles(p):
    todo=set(range(len(p)));out=[]
    while todo:
        a=min(todo);c=[];j=a
        while j in todo:
            todo.remove(j);c.append(j);j=p[j]
        out.append(c)
    return out
