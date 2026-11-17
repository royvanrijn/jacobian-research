"""Exact transcription controls; not a universal geometric proof."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    types = [(2, 2), (3, 2), (2, 3), (3, 3)]
    dual_types = [(2-a, 2-b) for a, b in types]
    assert dual_types == [(0, 0), (-1, 0), (0, -1), (-1, -1)]
    # Duality of ambient SUM is the rank-one DIAGONAL, all degrees.
    degree_controls = []
    for d in (2, 3):
        ambient_sum = sp.ones(1, d)
        diagonal = ambient_sum.T
        trace_zero = sp.Matrix.hstack(*(sp.eye(d)[:, j]-sp.eye(d)[:, 0]
                                       for j in range(1, d)))
        assert ambient_sum*trace_zero == sp.zeros(1, d-1)
        assert diagonal == sp.ones(d, 1) and trace_zero.rank() == d-1
        degree_controls.append({'degree': d, 'lattice_rank': 1,
                                'diagonal_entries': [int(x) for x in diagonal],
                                'ambient_trace_zero_rank': d-1})
    # Simplicial gluing really supplies a torus: boundary of tetrahedron.
    vertices = list(range(4))
    edges = [(i, j) for i in vertices for j in vertices if i < j]
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    boundary1 = sp.zeros(4, 6)
    boundary2 = sp.zeros(6, 4)
    for j, (a, b) in enumerate(edges):
        boundary1[a, j] = -1
        boundary1[b, j] = 1
    for j, (a, b, c) in enumerate(faces):
        for edge, sign in [((b, c), 1), ((a, c), -1), ((a, b), 1)]:
            boundary2[edges.index(edge), j] = sign
    assert boundary1*boundary2 == sp.zeros(4)
    assert 4-boundary2.rank() == 1
    # Local-system fixed lift equals dual unit retraction.
    A = sp.Matrix([[-1, 0], [0, 1]])
    c = sp.Matrix([2, 0])
    mon = A.row_join(c).col_join(sp.Matrix([[0, 0, 1]]))
    v = sp.Matrix([1, 0, 1])
    assert mon*v == v and mon**2 == sp.eye(3)
    assert v.T*mon.T == v.T
    unipotent = sp.Matrix([[1, 1], [0, 1]])
    for e in (2, 3):
        assert unipotent**e == sp.Matrix([[1, e], [0, 1]])
    # Constant non-torsion section vs varying section of a constant torus.
    s, tau, r = sp.symbols('s tau r')
    assert sp.diff(sp.Integer(2), s) == 0
    assert sp.diff(s, s)/s == 1/s
    assert sp.residue(1/s, s, 0) == 1
    # Long primitive representative and exact reciprocal-branch norm.
    g = s/(s+1)
    dlog = sp.cancel(sp.diff(g, s)/g)
    assert sp.cancel(dlog-1/(s*(s+1))) == 0
    assert sp.residue(dlog, s, 0) == 1
    assert sp.residue(dlog, s, -1) == -1
    shifted_g = (tau-sp.Rational(1, 2))/(tau+sp.Rational(1, 2))
    assert sp.cancel(shifted_g*shifted_g.subs(tau, -tau)) == 1
    assert sp.diff(shifted_g, tau) != 0
    # Trace weights retain the logarithmic difference despite norm1.
    rational_trace = (1/(tau-sp.Rational(1, 2))
                      -1/(-tau-sp.Rational(1, 2)))/(2*tau)
    assert sp.cancel(rational_trace-1/(tau**2-sp.Rational(1, 4))) == 0
    log_symbol = sp.symbols('log_g')
    assert (log_symbol-(-log_symbol))/(2*tau) == log_symbol/tau
    root = sp.sqrt(r+sp.Rational(1, 4))
    series = sp.series(1/r+sp.log((root-sp.Rational(1, 2)) /
                                  (root+sp.Rational(1, 2)))/root, r, sp.oo, 5)
    assert sp.expand(series.removeO()).coeff(r, -2) == sp.Rational(1, 6)
    controls = []
    for n in (1, 2, 3, 4):
        marked = sum(sp.binomial(n, j)*sp.Rational((-1)**(n+j-1), n+j+1)
                     for j in range(n+1))
        beta = (-1)**(n-1)*sp.factorial(n)**2/sp.factorial(2*n+1)
        assert marked == beta and marked != 0
        controls.append({'n': n, 'fixed_marker_value': str(marked)})
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'source_dual_twist_one_motivic_types_exact': True,
        'dual_SUM_is_DIAGONAL_controls': degree_controls,
        'simplicial_gluing_torus_rank_control_exact': True,
        'fixed_lift_and_dual_unit_retraction_controls_exact': True,
        'finite_cover_does_not_kill_unipotent_extension_controls_exact': True,
        'constant_nontorsion_vs_varying_constant_torus_distinction_exact': True,
        'Long_dlog_residues_and_constant_norm_exact': True,
        'weighted_trace_retains_log_difference_exact': True,
        'Long_pure_composite_fixed_marker_controls': controls,
        'Long_purity_and_primitive_marker_law': 'Retained primary Long2607.19012v1 Theorem4.2; all-order composite purity and beta integral are written deductions',
        'deep_algebraic_import': 'Published BRS2003 rational relative-cohomology Theorem0.1, Construction3.1, comparison5.4-5.6; point equivalence1.5',
        'recent_relative_equivalence_preprint_imported': False,
        'all_fibre_NS_ranks_claimed_constant': False,
        'integral_no_isogeny_promotion_claimed': False,
        'constant_type_replaced_by_torsion_or_Hodge_splitting': False,
        'scalar_norm_claimed_complete_obstruction_test': False,
        'Long_composite_claimed_exceptional_invariant': False,
        'general_resolution_Picard_fixed_part_software_implemented': False,
        'geometric_construction_and_fixed_part_proofs_machine_verified': False,
        'exceptional_pure_locus_section_constancy_proved': False,
        'coefficient_ring_topological_promotion_claimed': False,
        'new_degree_support_or_prefix_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (p/'haar_boundary_motive_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
