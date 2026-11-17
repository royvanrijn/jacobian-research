"""PF110 exact source controls; no function certifies a finite Keller target."""
import sympy as S


def circle_bundle_indices(p):
    if not isinstance(p,int) or isinstance(p,bool) or p<2:
        raise ValueError('oriented negative circle-bundle order p>=2 required')
    return [k+1 for k in range(1,p-1) if (k-p)%2==0]


def three_hole_factorization(c):
    c=S.sympify(c)
    if not c.is_Integer or c not in (0,1):
        raise ValueError('positive multiplicities require c=0 or1')
    return dict(single=c,pair=1-c,triple=c,
                length=3*c+3*(1-c)+c,
                b2_if_simply_connected=c)


def source_control(p):
    p=S.sympify(p)
    if not p.is_Integer or p<3:
        raise ValueError('the globally regular displayed primitive requires p>=3')
    return dict(p=int(p),Hirzebruch_index=int(p-2),section_square=int(p),
                section_minimal_section_intersection=1,source_smooth=True,
                source_affine=True,whole_plane_open=True,omitted_prime_count=1,
                omitted_prime='A1',canonical_coefficient=int(p-2),
                numerical_index=int(p-1),boundary_Chern_residue=int(2%p),
                intersection_form=[[-int(p)]],oriented_end='negative-Euler circle bundle L(p,1)',
                globally_exact_volume=True,finite_target_supplied=False,
                actual_Keller_realization_claimed=False)


def torsion_planarity_gate(*,integral_boundary_Chern_zero,fill_Chern_nonzero):
    """Apply the imported OSS Cor1.5, retaining its integral-zero hypothesis."""
    return dict(planarity_excluded=integral_boundary_Chern_zero and fill_Chern_nonzero,
                rational_boundary_vanishing_alone_suffices=False,
                graph_base_genus_zero_implies_contact_planarity=False)
