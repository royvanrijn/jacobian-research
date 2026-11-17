"""Exact Rodrigues, density, marked-tail and local-log controls."""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    t, s, w = sp.symbols('t s w')
    def C(m):
        return sum((-1)**j*sp.binomial(m, j)*sp.binomial(m+2*j, m)
                   *t**(2*j) for j in range(m+1))
    def B(m, marker_power):
        return sum((-1)**j*sp.binomial(m, j)
                   *sp.binomial(m+2*j, m-marker_power)*t**(2*j)
                   for j in range(m+1))
    for m in [0, 1, 2, 3, 4]:
        assert sp.expand(C(m)-sp.diff(t**m*(1-t*t)**m, t, m)
                         /sp.factorial(m)) == 0
    pairing = sp.Matrix([[sp.integrate(t**k*C(m), (t, 0, 1))
                          for m in range(4)] for k in range(4)])
    assert pairing.is_lower
    assert pairing.det() != 0
    for m in range(4):
        assert pairing[m, m] == (-1)**m*sp.integrate(
            t**m*(1-t*t)**m, (t, 0, 1))
    examples = []
    for q in [sp.Integer(1), 2*t, 3*t**2, 2*(1-t),
              (1+sp.I*t)/(1+sp.I/2)]:
        q = sp.expand(q)
        assert sp.integrate(q, (t, 0, 1)) == 1
        degree = sp.degree(q, t)
        first = min(i for (i,), coeff in sp.Poly(q, t).terms() if coeff != 0)
        coeff = sp.expand(q).coeff(t, first)
        m = int(degree)+2
        assert sp.integrate(q*C(m), (t, 0, 1)) == 0
        assert sp.expand(B(m, first+1)-t**(-first-1)
                         *sp.diff(t**m*(1-t*t)**m, t, m-first-1)
                         /sp.factorial(m-first-1)) == 0
        tail = sp.integrate(q*B(m, first+1), (t, 0, 1))
        expected = coeff*sp.integrate(t**first*(1-t*t)**m, (t, 0, 1))
        assert sp.simplify(tail-expected) == 0 and tail != 0
        initial = sp.Matrix([sp.integrate(q*C(j), (t, 0, 1))
                             for j in range(int(degree)+1)])
        small = pairing[:int(degree)+1, :int(degree)+1]
        reconstructed = small.T.inv()*initial
        assert sp.expand(sum(reconstructed[j]*t**j
                             for j in range(len(reconstructed)))-q) == 0
        examples.append({'density': str(q), 'degree': int(degree),
                         'first_nonzero_index': int(first),
                         'marker_power': int(first+1),
                         'checked_tail_order': m, 'nonzero_marked_value': str(tail)})
    # Polynomial-density Cauchy decomposition, retaining both parities.
    ce0, ce1, co0, co1 = sp.symbols('ce0 ce1 co0 co1')
    qe = ce0+ce1*s
    qo = co0+co1*s
    q = qe.subs(s, t*t)+t*qo.subs(s, t*t)
    correction = sp.cancel((qe-qe.subs(s, t*t))/(s-t*t)
                           +t*(qo-qo.subs(s, t*t))/(s-t*t))
    assert sp.cancel(q/(s-t*t)-qe/(s-t*t)-qo*t/(s-t*t)+correction) == 0
    dL1 = 1/(w+1)-1/(w-1)
    dL2 = 2/w-2*w/(w*w-1)
    assert sp.residue(dL1, w, 1) == sp.residue(dL2, w, 1) == -1
    a, b = qe.subs(s, w*w)/(2*w), qo.subs(s, w*w)/2
    c = sp.cancel(b/a)
    assert sp.cancel(c.subs(w, -w)+c) == 0
    assert sp.diff(c, w) != 0
    # An endpoint zero in q does NOT remove the differential-module log.
    endpoint_q = 2*(1-t)
    endpoint_a, endpoint_b = 1/w, -1
    assert endpoint_q.subs(t, 1) == 0
    assert sp.cancel((endpoint_a+endpoint_b).subs(w, 1)) == 0
    assert sp.diff(endpoint_b/endpoint_a, w) == -1
    # Base uniform-even profile: any meromorphic unit retraction would
    # require a primitive of this form, whose residue is nonzero.
    forcing = -1/(2*sp.sqrt(s)*(s-1))
    assert sp.residue(forcing, s, 1) == -sp.Rational(1, 2)
    # Proper compact fibre weight3 has no unit N-image; the punctured
    # sphere-height fibre supplies the weight4 Kummer negative control.
    genus_N = sp.zeros(3)
    genus_N[1, 2] = 1
    unit = sp.eye(3)[:, 0]
    assert genus_N.rank() < genus_N.row_join(unit).rank()
    tate_N = sp.Matrix([[0, 1], [0, 0]])
    assert tate_N.rank() == tate_N.row_join(sp.eye(2)[:, 0]).rank()
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'Rodrigues_identity_and_triangular_density_pairing_controls': True,
        'density_pairing_determinant': str(pairing.det()),
        'complex_density_reconstruction_and_eventual_marker_controls': examples,
        'even_odd_Cauchy_decomposition_exact': True,
        'both_log_derivatives_have_nonzero_residue_at1': True,
        'parity_forces_mixed_log_coefficient_ratio_nonconstant': True,
        'density_endpoint_zero_does_not_remove_unit_log_control': True,
        'uniform_even_profile_meromorphic_primitive_obstruction': True,
        'proper_genus_weight3_and_punctured_Tate_weight4_distinguished': True,
        'universal_written_result': 'Proper ORIGINAL-volume quadric model gives locally split actual UNIT Cauchy module at everypuncture;oddNglobal. EventualHopfpurity iff polynomial-density even-heightlaw;allnonzero suchdensitylaws haveunitlogobstructionat1. Thus BOTHexceptional torus-equivariantcomplexHopf phases areNOTeventuallypure inANYdegree',
        'height_or_density_assumed_real_or_positive': False,
        'local_genus_semisimplicity_assumed': False,
        'actual_period_kernel_claimed_Hodge': False,
        'all_original_marker_forms_assumed_to_extend_to_Y': False,
        'arbitrary_local_quotient_reflection_claimed': False,
        'proper_original_volume_model_assumed_for_every_full_phase': False,
        'purity_of_polynomial_composite_reflected_to_primitive': False,
        'mixed_weight_Bezout_sections_claimed_excluded': False,
        'all_polynomial_composites_of_Hopf_phase_claimed_excluded': False,
        'all_Bezout_phases_or_composites_claimed_excluded': False,
        'full_exceptional_pure_locus_component_vanishing_proved': False,
        'arbitrary_coefficient_ring_promotion_claimed': False,
        'universal_Hodge_or_MZ_proof_machine_verified': False,
        'new_degree_support_or_parameter_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (p/'haar_proper_unit_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
