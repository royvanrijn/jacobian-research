"""Exact finite controls for full intrinsic LOCAL unit reflection."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    # Hom(P_nontriv,A) invariants evaluate to the nontrivial socle.
    P1 = sp.Matrix([[1, 1], [0, 1]])
    P2 = P1.T
    As = [sp.diag(1, *([1]*2), *([1]*2)) for _ in range(2)]
    for A, P in zip(As, [P1, P2]):
        A[1:3, 1:3] = P
    xs = sp.symbols('x0:10')
    X = sp.Matrix(5, 2, xs)
    equations = list(As[0]*X-X*P1)+list(As[1]*X-X*P2)
    solutions = sp.linsolve(equations, xs)
    row = next(iter(solutions))
    parameters = set().union(*(x.free_symbols for x in row))
    assert len(parameters) == 1
    parameter = next(iter(parameters))
    evaluation = sp.Matrix(5, 2, row).subs(parameter, 1)
    assert evaluation.columnspace() == [sp.eye(5)[:, 1], sp.eye(5)[:, 2]]
    # Pure genus removal preserves the local unit image test.
    removal = []
    for lam in (0, 1):
        for kappa in (1, 2):
            N = sp.zeros(4)
            N[1, 2] = 1
            N[0, 3], N[1, 3] = lam, kappa
            unit = sp.eye(4)[:, 0]
            C = sp.eye(4)[:, 1:3]
            assert (N*C).columnspace() == [sp.eye(4)[:, 1]]
            assert N**2 == sp.zeros(4)
            in_image = N.rank() == N.row_join(unit).rank()
            Nbar = sp.Matrix([[0, lam], [0, 0]])
            quotient_image = Nbar.rank() == Nbar.row_join(sp.eye(2)[:, 0]).rank()
            assert in_image == quotient_image == (lam != 0)
            removal.append({'lambda': lam, 'kappa': kappa,
                            'unit_in_source_and_genus_quotient_image': in_image})
    # Remaining constant kernel meets M4 only in the unit, which it avoids.
    # Basis u4,g5,h5,t6. A constant non-Hodge K can identify g with u,
    # but admissible local N has image in L, not in g.
    q = sp.Matrix([[1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    K = sp.Matrix([-1, 1, 0, 0])
    assert q*K == sp.zeros(3, 1)
    assert sp.Matrix.hstack(K, sp.eye(4)[:, 0]).rank() == 2
    for lam in (0, 1):
        N = sp.zeros(4)
        N[0, 3] = lam
        Nbar = sp.zeros(3)
        Nbar[0, 2] = lam
        assert q*N == Nbar*q
        assert (N.rank() == N.row_join(sp.eye(4)[:, 0]).rank()) == bool(lam)
        assert (Nbar.rank() == Nbar.row_join(sp.eye(3)[:, 0]).rank()) == bool(lam)
    # Local-only countercontrol cannot support a global nontrivial genus
    # quotient: its image of the removed genus would contain the unit.
    local_q = sp.Matrix([[1, 1, 0, 0], [0, 0, 1, 1]])
    local_C = sp.eye(4)[:, 1:3]
    assert (local_q*local_C).rank() == 2
    assert local_q*local_C[:, 0] == sp.eye(2)[:, 0]
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'Hom_invariant_evaluation_extracts_nontrivial_socle_exact': True,
        'pure_genus_M4_equals_N_image_controls_exact': True,
        'pure_genus_removal_unit_reflection_controls': removal,
        'constant_kernel_weight4_intersection_controls_exact': True,
        'nonHodge_kernel_still_cannot_create_local_unit_image_controls_exact': True,
        'local_only_counter_has_unit_in_nontrivial_genus_image_exact': True,
        'universal_result': 'ANY GLOBAL unit-retaining quotient of canonical boundary E reflects LOCAL unit splitting at EVERY puncture,before purity;under eventual purity full safety iff reduced intrinsic xi_boundary,0=0',
        'global_semisimplicity_of_V_assumed': False,
        'global_source_unit_splitting_assumed': False,
        'actual_kernel_claimed_Hodge': False,
        'arbitrary_local_quotient_reflection_claimed': False,
        'full_unreduced_xi_substituted_for_boundary_xi': False,
        'local_splitting_used_as_blanket_global_splitting_without_purity': False,
        'top_Tate_invariant_IMAGE_identified_with_entire_top_Tate': False,
        'exceptional_pure_locus_xi_vanishing_proved': False,
        'pure_composite_producer_implemented': False,
        'universal_Hodge_fixed_part_and_Borel_proofs_machine_verified': False,
        'coefficient_ring_topological_promotion_claimed': False,
        'new_degree_support_or_prefix_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (p/'haar_local_source_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
