"""Exact PF119 necessary balances; none constructs a plane Keller map."""


def half_support_cases(r):
    if not isinstance(r,int) or r<1:raise ValueError('positive integer ramified-prime count required')
    return [dict(degree=4*r,ramified_primes=r,nonimmersion_defect=M,
                 ordinary_cusps=M,total_excess=r*M,unramified_primes=2*r-1-r*M,
                 source_points_at_double_values=0,
                 actual_plane_realization=False,binary_affine_family_excluded=False)
            for M in [0,1]]


def source_fibre(degree,ramified_images,unramified_branch_counts):
    """ramified_images contains (number_of_primes, list_of_branch_multiplicities)."""
    if degree<1:raise ValueError('positive degree required')
    lengths=[]
    for r,ms in ramified_images:
        if r<1 or any(m<1 for m in ms):raise ValueError('positive ramified counts and multiplicities required')
        lengths.extend([m+1 for m in ms for _ in range(r)])
    if any(b<0 for b in unramified_branch_counts):raise ValueError('nonnegative branch counts required')
    return dict(ramified_local_lengths=lengths,unramified_omitted_points=sum(unramified_branch_counts),
                source_points=degree-sum(lengths)-sum(unramified_branch_counts))
