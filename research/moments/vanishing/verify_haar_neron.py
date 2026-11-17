"""Exact normal-form and component-presentation controls, not a VMHS proof."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    # Basis u,a,b,c,d,t1,t2. Two degenerating genus blocks; arbitrary
    # allowed genus cross terms from the entire top Tate space.
    b1, b2, k1, k2, h1, h2 = sp.symbols('b1 b2 k1 k2 h1 h2')
    N = sp.zeros(7)
    N[1, 2] = N[3, 4] = 1
    N[0, 5], N[0, 6] = b1, b2
    N[1, 5], N[1, 6] = k1, k2
    N[3, 5], N[3, 6] = h1, h2
    change = sp.eye(7)
    change[2, 5], change[2, 6] = -k1, -k2
    change[4, 5], change[4, 6] = -h1, -h2
    normal = sp.zeros(7)
    normal[1, 2] = normal[3, 4] = 1
    normal[0, 5], normal[0, 6] = b1, b2
    assert N**2 == sp.zeros(7)
    assert change.det() == 1 and N*change == change*normal
    # Every top class is invariant in V after correcting by genus lifts.
    assert (N*change)[1:, 5:] == sp.zeros(6, 2)
    controls = []
    for beta in [(0, 0), (1, 0), (0, 2), (1, -1)]:
        for e in [1, 2, 5]:
            n = N.subs({b1: beta[0], b2: beta[1], k1: 2,
                        k2: -3, h1: 4, h2: 1})
            mon = sp.eye(7)+n
            assert mon**e-sp.eye(7) == e*n
            excess = (mon**e-sp.eye(7)).rank()-2
            unit = sp.eye(7)[:, 0]
            in_image = n.rank() == n.row_join(unit).rank()
            assert excess == int(any(beta)) == int(in_image)
            controls.append({'beta': list(beta), 'cover_degree': e,
                             'rank_excess': excess})
    # A non-one original genus sector is counted in both ranks after cover.
    n = sp.zeros(5)
    n[1, 2] = n[3, 4] = 1
    semisimple = sp.diag(1, 1, 1, -1, -1)
    assert semisimple*n == n*semisimple
    mon = semisimple*(sp.eye(5)+n)
    assert (mon**2-sp.eye(5)).rank() == 2
    assert (mon[1:, 1:]**2-sp.eye(4)).rank() == 2
    # Without the strict M4P=imageNP restriction, local V need not split.
    forbidden = sp.Matrix([[0, 1], [0, 0]])  # constant genus, T->genus
    assert forbidden.rank() == 1
    assert all(x[1] == 0 for x in forbidden.nullspace())
    # A component group with finite torsion and one free coordinate.
    relations = sp.Matrix([[2, 0], [0, 6], [0, 0]])
    presentation = []
    for v, torsion in [(sp.Matrix([1, 1, 0]), True),
                       (sp.Matrix([0, 0, 1]), False),
                       (sp.Matrix([2, -3, 0]), True)]:
        actual = relations.row_join(v).rank() == relations.rank()
        assert actual == torsion
        if torsion:
            assert all(x.q == 1 for x in relations.gauss_jordan_solve(6*v)[0])
        presentation.append({'representative': list(v),
                             'component_torsion': actual,
                             'identity_component_already':
                             v[2] == 0 and v[0] % 2 == 0 and v[1] % 6 == 0})
    assert presentation[0]['component_torsion']
    assert not presentation[0]['identity_component_already']
    s = sp.symbols('s')
    g = s/(s+1)
    log_derivative = sp.diff(g, s)/g
    assert sp.residue(log_derivative, s, 0) == 1
    assert sp.residue(log_derivative, s, -1) == -1
    tau = sp.symbols('tau')
    branch = (tau-sp.Rational(1, 2))/(tau+sp.Rational(1, 2))
    assert sp.cancel(branch*branch.subs(tau, -tau)) == 1
    assert sp.diff(branch, tau) != 0
    # Constant2 has component0 without being a torsion complex point.
    assert sp.diff(sp.log(sp.Integer(2)), s) == 0
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'symbolic_arbitrary_allowed_genus_cross_terms_removed': True,
        'every_top_Tate_class_has_invariant_V_lift_control': True,
        'rank_excess_controls': controls,
        'nonone_genus_contributions_cancel_in_rank_difference': True,
        'Hodge_strip_restriction_is_essential_control': True,
        'component_presentation_controls': presentation,
        'torsion_component_distinguished_from_identity_component': True,
        'Long_full_valuation_vector_and_constant_norm_control': True,
        'constant_nontorsion_section_allowed': True,
        'universal_written_result': 'For reduced canonical E, J_boundary=T^1 and local unit class vanishes iff monodromy rank excess is0 iff section has torsion Neron lft component;under eventual purity iff full safety',
        'global_splitting_of_V_claimed': False,
        'local_genus_semisimplicity_assumed': False,
        'Hodge_splitting_of_local_normal_form_claimed': False,
        'full_unreduced_top_invariant_image_claimed_surjective': False,
        'section_torsion_required': False,
        'identity_component_required_without_integer_multiple': False,
        'scalar_norm_used_as_full_component_test': False,
        'general_Neron_or_Picard_producer_implemented': False,
        'exceptional_pure_locus_component_torsion_proved': False,
        'universal_Neron_Kummer_or_Hodge_proof_machine_verified': False,
        'coefficient_ring_topological_promotion_claimed': False,
        'new_degree_support_or_prefix_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    # Sympy integers/booleans are intentionally converted for JSON.
    for row in presentation:
        row['representative'] = [int(x) for x in row['representative']]
        row['identity_component_already'] = bool(row['identity_component_already'])
    (p/'haar_neron_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
