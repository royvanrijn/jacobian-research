"""Exact controls for cohomological pullback and a uniform original marker."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    d, k = sp.symbols('d k')
    target = sp.Matrix([[0, 1], [0, 0]])
    source = sp.zeros(4)  # u,a,b,t; genus b->a, toric t->d*u+k*a
    source[1, 2], source[0, 3], source[1, 3] = 1, d, k
    alpha = sp.Matrix([[d, 0], [0, 0], [0, -k], [0, 1]])
    assert source*alpha == alpha*target
    assert source**2 == sp.zeros(4)
    assert alpha[:, 0] == d*sp.eye(4)[:, 0]
    for degree in [1, 2, 3]:
        n = source.subs({d: degree, k: 3})
        assert n.rank()-1 == 1
        assert n.rank() == n.row_join(sp.eye(4)[:, 0]).rank()
    # d=0 must be excluded: this map kills the target unit and obstruction.
    zero = sp.zeros(2)
    zero_alpha = sp.Matrix([[0, 0], [0, 1]])
    assert zero*zero_alpha == zero_alpha*target
    assert zero_alpha[:, 0] == sp.zeros(2, 1)
    # A nontrivial global simple block has no map into a trivial lowest line.
    r1, r2 = sp.symbols('r1 r2')
    row = sp.Matrix([[r1, r2]])
    generators = [sp.Matrix([[1, 1], [0, 1]]),
                  sp.Matrix([[1, 0], [1, 1]])]
    assert sp.solve(list(row*generators[0]-row)+list(row*generators[1]-row),
                    [r1, r2]) == {r1: 0, r2: 0}
    s, t = sp.symbols('s t')
    cs = sp.symbols('c0:7')
    R = sum(cs[i]*s**i for i in range(len(cs)))
    interval = sp.integrate(R.subs(s, -t), (t, 0, 1))
    pure = cs[0]
    marked_Q = sum(cs[i]*sp.Rational((-1)**(i-1), i+1)
                   for i in range(1, len(cs)))
    assert sp.expand(pure-marked_Q-interval) == 0
    # No H-specific marker: the same functional evaluates every composition.
    examples = []
    for H in [s, s+2, s*(s+1), s*(s+1)*(s+sp.Rational(1, 2))]:
        values = []
        for m in range(1, 7):
            hm = sp.Poly(H**m, s)
            direct = sp.integrate(H.subs(s, -t)**m, (t, 0, 1))
            coefficient = sum(c*sp.Rational((-1)**mon[0], mon[0]+1)
                              for mon, c in hm.terms())
            assert direct == coefficient
            values.append(str(direct))
        examples.append({'H': str(H), 'marked_values_orders1_to6': values})
    odd_control = (s*(s+1)*(s+sp.Rational(1, 2))).subs(s, -t)
    assert sp.expand(odd_control.subs(t, 1-t)+odd_control) == 0
    for m in [1, 3, 5]:assert sp.integrate(odd_control**m, (t, 0, 1)) == 0
    for m in [2, 4, 6]:assert sp.integrate(odd_control**m, (t, 0, 1)) > 0
    degrees = [{'Gamma': group, 'n': n, 'm': m, 'order': order,
                'integral_H3_degree': n*m//order}
               for group, n, m, order in [('2O', 8, 12, 48),
                                          ('2I', 12, 20, 120)]]
    assert all(x['integral_H3_degree'] == 2 for x in degrees)
    # Distinguish normalized J=2W from the raw nm factor.
    for x in degrees:assert sp.Rational(x['n']*x['m'], x['order']) == 2
    g = s/(s+1)
    dg = sp.diff(g, s)/g
    assert sp.residue(dg, s, 0) == 1
    assert sp.residue(dg, s, -1) == -1
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'symbolic_nonzero_unit_map_and_local_obstruction_transport': True,
        'zero_H3_degree_negative_control': True,
        'nontrivial_global_block_has_no_unit_map_control': True,
        'symbolic_full_polynomial_uniform_marker_identity': True,
        'same_marker_composition_controls': examples,
        'odd_orders_may_vanish_even_orders_survive_control': True,
        'integral_H3_degrees': degrees,
        'normalized_and_raw_volume_factors_distinguished': True,
        'actual_Long_local_classes_at0_and_minus1_pinned': True,
        'universal_written_result': 'NonzeroH3 regular pullback preserves nonzero reduced local component. For every pullback g of LongF one fixed original C=(J/d)(1-Q)pullback satisfies h(C R(g))=integral_0^1 R(-t)dt for allR;MugerTuset implies infinitelymany markednonzeroorders for every nonzeroH(g)',
        'proper_or_etale_map_assumed': False,
        'generic_field_degree_identified_with_H3_degree': False,
        'reverse_local_obstruction_implication_claimed': False,
        'marker_depends_on_H_degree_or_order': False,
        'all_positive_marked_orders_claimed_nonzero_for_every_H': False,
        'Jacobian_weight_dropped_from_pure_moments': False,
        'invariance_alone_claimed_to_force_component_torsion': False,
        'eventual_pure_Bezout_phase_or_composite_constructed': False,
        'all_Bezout_phases_or_composites_excluded_from_purity': False,
        'sufficient_pullback_family_claimed_to_cover_full_quotient': False,
        'arbitrary_coefficient_ring_promotion_claimed': False,
        'universal_Hodge_or_interval_proof_machine_verified': False,
        'new_degree_support_or_parameter_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (p/'haar_pullback_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
