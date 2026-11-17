"""Finite exact controls for the graph-necessity written theorem."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    here = Path(__file__).resolve().parent
    # Source: unit u,constant genus g,h,and irreducible genus a,b.
    # Two global generators make the a,b plane irreducible/semisimple,
    # while each individual local matrix is nonsemisimple.
    plus = sp.eye(5)
    minus = sp.eye(5)
    plus[3, 4] = 1
    minus[4, 3] = 1
    unit = sp.eye(5)[:, 0]
    kernel = sp.eye(5)[:, 1]-unit
    q = sp.Matrix([[1, 1, 0, 0, 0], [0, 0, 1, 0, 0],
                   [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    # Projection E -> K fixes K and kills the prescribed unit.
    retraction_K = kernel*sp.Matrix([[0, 1, 0, 0, 0]])
    assert retraction_K*kernel == kernel
    assert retraction_K*unit == sp.zeros(5, 1)
    assert retraction_K**2 == retraction_K
    for a in [plus, minus]:
        assert a*kernel == kernel
        assert a*retraction_K == retraction_K*a
    inclusion = sp.Matrix([[1, 0, 0, 0], [0, 0, 0, 0],
                           [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    assert q*inclusion == sp.eye(4)
    assert inclusion*q == sp.eye(5)-retraction_K
    assert q*unit == sp.eye(4)[:, 0]
    for a in [plus, minus]:
        quotient = q*a*inclusion
        assert q*a == quotient*q
        assert quotient[0, :] == sp.Matrix([[1, 0, 0, 0]])
        assert quotient-sp.eye(4) != sp.zeros(4)
        assert (quotient-sp.eye(4))**2 == sp.zeros(4)
    # The two elementary transvections generate the full 2x2 matrix
    # algebra on genus, despite neither local transvection being SS.
    e12 = sp.Matrix([[0, 1], [0, 0]])
    e21 = e12.T
    matrices = [e12, e21, e12*e21, e21*e12]
    assert sp.Matrix.hstack(*(a.reshape(4, 1) for a in matrices)).rank() == 4

    # Negative control with admissible LOCAL relative-weight data.
    # Original W: u4,a5,b5,t6; relative M: u4,a4,b6,t6.
    N = sp.zeros(4)
    N[1, 2] = N[1, 3] = 1
    original_W = [4, 5, 5, 6]
    relative_M = [4, 4, 6, 6]
    assert N**2 == sp.zeros(4)
    assert all(relative_M[i] <= relative_M[j]-2
               for i in range(4) for j in range(4) if N[i, j] != 0)
    assert any(original_W[i] > original_W[j]-2
               for i in range(4) for j in range(4) if N[i, j] != 0)
    # Source unit has a local retraction, quotient unit does not.
    source_ell = sp.Matrix([[1, 0, 0, 0]])
    assert source_ell*N == sp.zeros(1, 4)
    K = sp.Matrix([[-1, 0], [1, 0], [0, 1], [0, -1]])
    quotient_map = sp.Matrix([[1, 1, 0, 0], [0, 0, 1, 1]])
    quotient_N = sp.Matrix([[0, 1], [0, 0]])
    assert N*K == sp.zeros(4, 2)
    assert quotient_map*K == sp.zeros(2)
    assert quotient_map*N == quotient_N*quotient_map
    assert quotient_map*sp.eye(4)[:, 0] == sp.eye(2)[:, 0]
    assert quotient_N*sp.eye(2)[:, 1] == sp.eye(2)[:, 0]
    # Limit F3 span(b,t) maps under N into F2=all, as type(-1,-1).
    F3 = sp.eye(4)[:, 2:4]
    assert (N*F3).rank() == 1
    # Cyclic induction's exact coinvariant map at ramification e.
    ramified = []
    A = sp.Matrix([[1, 1], [0, 1]])
    c = sp.Matrix([0, 1])
    for e in (2, 3):
        S = sp.zeros(2*e)
        for j in range(e-1):
            S[2*(j+1):2*(j+2), 2*j:2*j+2] = sp.eye(2)
        S[:2, -2:] = A
        Phi = sp.Matrix.hstack(*(sp.eye(2) for _ in range(e)))
        induced_c = sp.zeros(2*e, 1)
        induced_c[:2, :] = c
        assert Phi*induced_c == c
        assert (Phi*(S-sp.eye(2*e))).columnspace() == (A-sp.eye(2)).columnspace()
        assert (2*e)-(S-sp.eye(2*e)).rank() == 2-(A-sp.eye(2)).rank()
        ramified.append({'ramification': e, 'intrinsic_class_factor': 1,
                          'coinvariant_dimension': 1})
    # A forest has no top graph cycles even with arbitrary vertex genus.
    incidence = sp.Matrix([[-1, 0, 0], [1, -1, 0],
                            [0, 1, -1], [0, 0, 1]])
    assert 3-incidence.rank() == 0
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'unit_preserving_kernel_projection_and_global_complement_exact': True,
        'global_semisimple_genus_with_NONSEMISIMPLE_local_matrices_exact': True,
        'unit_retraction_passes_actual_style_nonHodge_quotient_control_exact': True,
        'relative_limit_M_not_original_W_control_exact': True,
        'local_split_source_unit_can_fail_in_quotient_without_global_SS_control_exact': True,
        'negative_control_only_local_not_actual_global_phase': True,
        'ramified_coinvariant_factor_one_controls': ramified,
        'forest_top_graph_rank_zero_control_exact': True,
        'universal_theorem': 'Under eventual purity GrW4 H2(generic fibre)=0 implies ALL original marked safety,with arbitrary boundary genus/compact middle/nonsemisimple local genus;EVERY eventual-pure failure requires nonzero Tate graph',
        'general_SS_quotient_gate': 'For V=E/L GLOBAL semisimple,eventual-pure safety iff intrinsic reduced xi_boundary,0=0;all H retain all branches and Bad_boundary root criterion',
        'universal_geometric_and_Borel_proofs_machine_verified': False,
        'all_exceptional_phases_claimed_top_graph_free': False,
        'all_boundary_genus_claimed_zero': False,
        'global_SS_replaced_by_local_SS_or_finite_monodromy': False,
        'reduced_xi_identified_with_unreduced_xi_without_proof': False,
        'arbitrary_actual_kernel_claimed_Hodge': False,
        'source_local_split_claimed_sufficient_without_global_SS_quotient': False,
        'hitting_bad_value_claimed_producer_of_eventual_purity': False,
        'primitive_purity_reflection_claimed': False,
        'coefficient_ring_topological_promotion_claimed': False,
        'new_degree_support_or_prefix_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (here/'haar_graph_necessity_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
