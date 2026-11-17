"""Exact controls; universal module/geometry arguments are in the proof note."""
from pathlib import Path
import hashlib
import json
import sympy as S


def fixed(gens):
    n = gens[0].rows
    return S.Matrix.vstack(*[g - S.eye(n) for g in gens]).nullspace()


def span_equal(a, b):
    return a.rank() == b.rank() == a.row_join(b).rank()


def quotient(gens, q, section):
    assert q * section == S.eye(q.rows)
    out = [q * g * section for g in gens]
    assert all(q * g == h * q for g, h in zip(gens, out))
    return out


def unit_splits(gens, u):
    return any((v.T * u)[0] != 0 for v in fixed([g.T for g in gens]))


def cyclic(a, e):
    n = a.rows
    g = S.zeros(e * n)
    g[:n, (e - 1) * n:] = a
    for j in range(1, e):
        g[j * n:(j + 1) * n, (j - 1) * n:j * n] = S.eye(n)
    return g


def run():
    # Unit, two fixed-genus directions, two varying-genus directions, top.
    ug = S.eye(6); ug[3, 4] = 1
    vg = S.eye(6); vg[4, 3] = 1
    trans = S.eye(6); trans[1, 5] = 1
    finite = S.diag(1, 1, 1, 1, 1, -1)
    gens = [ug, vg, trans, finite]
    u = S.eye(6)[:, 0]
    assert unit_splits(gens, u)
    w5 = [g[:5, :5] for g in gens]
    f_basis = S.Matrix.hstack(*fixed(w5))
    assert f_basis.cols == 3
    assert span_equal(f_basis, S.eye(5)[:, :3])
    q = S.eye(6)[[0, 1, 2, 5], :]
    red = quotient(gens, q, q.T)
    assert len(fixed([g[:3, :3] for g in red])) == 3
    r = (red[2] - S.eye(4))[:, 3:4]
    assert r == S.eye(4)[:, 1:2]
    assert (red[2] - S.eye(4))**2 == S.zeros(4)
    assert unit_splits(red, q * u)

    # An ACTUAL relation can mix the unit with a radical fixed-genus direction.
    # This is an abstract matrix model, not an exceptional phase realization.
    bad_q = S.Matrix([[1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    bad_sec = S.Matrix([[1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 1]])
    bad = quotient(red, bad_q, bad_sec)
    assert bad_q * q * u != S.zeros(3, 1)
    assert not unit_splits(bad, bad_q * q * u)
    k_bad = S.Matrix([1, -1, 0, 0])
    assert bad_q * k_bad == S.zeros(3, 1)
    assert r.row_join(k_bad).row_join(q * u).rank() == r.row_join(k_bad).rank()
    good_q = S.eye(4)[[0, 1, 3], :]
    good = quotient(red, good_q, good_q.T)
    assert unit_splits(good, good_q * q * u)
    k_good = S.eye(4)[:, 2:3]
    assert r.row_join(k_good).row_join(q * u).rank() > r.row_join(k_good).rank()

    # No fixed genus: varying-genus radical directions cannot become a unit.
    a = S.eye(4); a[1, 2] = 1
    b = S.eye(4); b[2, 1] = 1
    c = S.eye(4); c[1, 3] = 1
    no_fixed = [a, b, c]
    assert fixed([g[1:3, 1:3] for g in no_fixed]) == []
    q0 = S.eye(4)[[0, 3], :]
    zero_trace = quotient(no_fixed, q0, q0.T)
    assert all(g == S.eye(2) for g in zero_trace)

    # Recomposition SUM identifies the ENTIRE fixed bottom, not just its unit.
    ram = []
    primitive = red[2]
    for e in [2, 3]:
        merge = S.zeros(3 + e, 4 * e)
        sec = S.zeros(4 * e, 3 + e)
        for k in range(3):
            sec[k, k] = 1
        for j in range(e):
            for k in range(3):
                merge[k, 4 * j + k] = 1
            merge[3 + j, 4 * j + 3] = 1
            sec[4 * j + 3, 3 + j] = 1
        assert merge * sec == S.eye(3 + e)
        ind = cyclic(primitive, e)
        comp = merge * ind * sec
        assert merge * ind == comp * merge
        delta = comp**e - S.eye(3 + e)
        assert delta.rank() == 1 and delta**2 == S.zeros(3 + e)
        assert span_equal(delta, S.eye(3 + e)[:, 1:2])
        ram.append({'degree': e, 'fixed_bottom_rank': 3,
                    'translation_range_rank': 1, 'radical_square_zero': True})

    # Lagrange trace identity in the generic finite function-field algebra.
    t, x, z, aa, bb = S.symbols('t x z aa bb')
    trace_controls = []
    for d, h in [(2, t**2 + aa * t + bb), (3, t**3 + aa * t + bb)]:
        poly = S.Poly(h - z, t)
        companion = S.zeros(d)
        for j in range(d - 1):
            companion[j + 1, j] = 1
        for j in range(d):
            companion[j, d - 1] = -poly.nth(j)
        deriv = S.diff(h, t)
        dh = sum((coef * companion**power[0] for power, coef in S.Poly(deriv, t).terms()), S.zeros(d))
        for j in range(d):
            lhs = (companion**j * dh.inv() * (companion - x * S.eye(d)).inv()).trace()
            rhs = x**j / (z - h.subs(t, x))
            assert S.cancel(lhs - rhs) == 0
        trace_controls.append({'degree': d, 'all_d_trace_basis_identities_exact': True,
                               'rational_trace_kernel_rank': d - 1})

    # Genuine Long fibre isomorphism; no expansion to a deeper moment prefix.
    s, uu, vv = S.symbols('s uu vv')
    a0 = (uu * vv + s) / (s + 1)
    yy = a0 - 1
    zz = (s - a0) / a0**2
    xx, ww = uu / a0, vv / a0
    assert S.cancel(xx * ww - yy * zz - 1) == 0
    assert S.cancel((1 + yy) * (1 + zz + yy * zz) - s) == 0
    cq = 1 / s - S.log(1 + 1 / s)
    diff = S.diff(cq - 1 / s, s)
    assert S.simplify(diff - 1 / (s * (s + 1))) == 0
    assert S.residue(diff, s, 0) == 1 and S.residue(diff, s, -1) == -1
    return {
        'completed': True, 'sympy_version': S.__version__,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'split_source_varying_genus_removal_unit_Hom_controls_exact': True,
        'constant_genus_actual_relation_can_create_residual_obstruction_exact': True,
        'zero_fixed_genus_selected_model_closes_residual_exact': True,
        'translation_span_and_radical_square_zero_controls_exact': True,
        'fixed_bottom_SUM_recomposition_controls': ram,
        'actual_function_trace_controls': trace_controls,
        'genuine_Long_fibre_isomorphism_identities_exact': True,
        'genuine_Long_Cauchy_nonrationality_residues_exact': True,
        'Long_all_order_purity_and_marked_law_imported_from_primary': True,
        'universal_constant_genus_and_actual_semisimplicity_reflection_machine_verified': False,
        'global_split_source_claimed_general_safety_without_actual_relation_test': False,
        'constant_genus_rank_claimed_zero_universally': False,
        'reduced_actual_quotient_claimed_scalar_function_submodule': False,
        'algebraic_1_motive_or_relation_producer_implemented': False,
        'Long_phase_claimed_right_exceptional_invariant': False,
        'abstract_models_claimed_actual_exceptional_phases': False,
        'pure_composite_existence_claimed': False,
        'all_exceptional_intrinsic_and_residual_obstructions_proved_zero': False,
        'all_coefficient_ring_topological_promotion_claimed': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'new_degree_support_or_prefix_campaign': False,
    }


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('haar_constant_genus_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
