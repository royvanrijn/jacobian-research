"""Selected exact controls for the written universal boundary-reduction proof."""
from pathlib import Path
import hashlib
import json
import sympy as S


def fixed(gens):
    n = gens[0].rows
    return S.Matrix.vstack(*[g - S.eye(n) for g in gens]).nullspace()


def span_equal(a, b):
    return a.rank() == b.rank() == a.row_join(b).rank()


def matrix_basis(mats):
    flat = S.Matrix.hstack(*[x.reshape(x.rows**2, 1) for x in mats])
    return [mats[i] for i in flat.rref()[1]]


def algebra(gens):
    out = matrix_basis([S.eye(gens[0].rows)] + gens)
    while True:
        more = matrix_basis(out + [a * g for a in out for g in gens])
        if len(more) == len(out):
            return more
        out = more


def cyclic(a, e):
    n = a.rows
    g = S.zeros(e * n)
    g[:n, (e - 1) * n:] = a
    for j in range(1, e):
        g[j * n:(j + 1) * n, (j - 1) * n:j * n] = S.eye(n)
    return g


def ambient_sum(n, e):
    out = S.zeros(1 + e * (n - 1), n * e)
    sec = S.zeros(n * e, out.rows)
    sec[0, 0] = 1
    for j in range(e):
        out[0, j * n] = 1
        for k in range(1, n):
            out[1 + j * (n - 1) + k - 1, j * n + k] = 1
            sec[j * n + k, 1 + j * (n - 1) + k - 1] = 1
    assert out * sec == S.eye(out.rows)
    return out, sec


def run():
    # Pure compact-middle Sym2 has no invariant under the ENTIRE group.
    u = S.Matrix([[1, 2, 1], [0, 1, 1], [0, 0, 1]])
    v = S.Matrix([[1, 0, 0], [1, 1, 0], [1, 2, 1]])
    assert fixed([u, v]) == []
    assert fixed([u])  # A local invariant is not a global fixed class.
    gu, gv = S.diag(1, u, 1), S.diag(1, v, 1)
    w4 = S.eye(5)[:, :4]
    c = S.eye(5)[:, 1:4]
    generated = S.Matrix.hstack((gu - S.eye(5)) * w4,
                                (gv - S.eye(5)) * w4)
    assert span_equal(generated, c)
    invariant_w4 = fixed([gu[:4, :4], gv[:4, :4]])
    assert len(invariant_w4) == 1 and invariant_w4[0] == S.eye(4)[:, 0]
    projection = S.eye(5)[[0, 4], :]
    section = S.eye(5)[:, [0, 4]]
    unit = S.eye(5)[:, 0]
    controls = []
    for name, target in [('unit_boundary_nonsplit', 0),
                         ('compact_extension_removed', 1)]:
        trans = S.eye(5)
        trans[target, 4] = 1
        gens = [gu, gv, trans]
        red = [projection * g * section for g in gens]
        assert all(projection * g == h * projection for g, h in zip(gens, red))
        full_rows = fixed([g.T for g in gens])
        red_rows = fixed([g.T for g in red])
        a = S.Matrix.hstack(*full_rows)
        b = S.Matrix.hstack(*[projection.T * row for row in red_rows])
        assert span_equal(a, b)
        split = any((row.T * unit)[0] != 0 for row in full_rows)
        assert split == (name == 'compact_extension_removed')
        # Both full sources are nonsemisimple, but only one obstructs its unit.
        mono = algebra(gens)
        tr = S.Matrix([[(x * y).trace() for y in mono] for x in mono])
        assert tr.rank() < len(mono)
        controls.append({'model': name, 'unit_splits_before_and_after': split,
                         'full_monodromy_algebra_dimension': len(mono),
                         'full_trace_rank': tr.rank()})

    # The local geometric-strip model loses its quadratic blocks on quotient.
    n = S.zeros(7)
    n[1, 2] = n[2, 3] = n[4, 5] = n[0, 6] = 1
    r = S.eye(7)[[0, 4, 5, 6], :]
    sec = r.T
    nr = r * n * sec
    assert r * n == nr * r
    assert n**3 == S.zeros(7) and (n**2).rank() == 1
    assert nr**2 == S.zeros(4) and nr.rank() == 2
    # Actual unit stays nonzero because the removed module has no trivial factor.
    assert r * S.eye(7)[:, 0] != S.zeros(4, 1)

    # Canonical compact removal commutes with finite-cover ambient SUM.
    ram = []
    prim = S.eye(5)
    prim[0, 4] = 1
    for e in [2, 3]:
        sum_full, sf = ambient_sum(5, e)
        sum_red, sr = ambient_sum(2, e)
        induced_reduction = S.diag(*([projection] * e))
        reduce_composite = S.zeros(1 + e, 1 + 4 * e)
        reduce_composite[0, 0] = 1
        for j in range(e):
            reduce_composite[1 + j, 1 + 4 * j + 3] = 1
        assert reduce_composite * sum_full == sum_red * induced_reduction
        ind_full = cyclic(prim, e)
        ind_red = cyclic(projection * prim * section, e)
        comp_full = sum_full * ind_full * sf
        comp_red = sum_red * ind_red * sr
        assert reduce_composite * comp_full == comp_red * reduce_composite
        c_ind_gens = [cyclic(u, e), S.diag(*([v] * e))]
        assert fixed(c_ind_gens) == []
        ram.append({'degree': e, 'full_composite_rank': 1 + 4 * e,
                    'removed_compact_rank': 3 * e,
                    'reduced_rank': 1 + e,
                    'induced_compact_global_invariants': 0})

    # Genus-weight zero: a nonzero radical image is exactly the ambient line.
    t = S.eye(3); t[0, 1] = 1
    j = t - S.eye(3)
    assert j.rank() == 1 and j**2 == S.zeros(3)
    assert span_equal(j, S.eye(3)[:, :1])
    quotient = S.eye(3)[[0, 1], :]
    assert quotient * t == S.Matrix([[1, 1], [0, 1]]) * quotient
    assert all(row[0] == 0 for row in fixed([t.T]))

    # A local square-zero bound is not a global algebra square-zero argument.
    a = S.zeros(3); a[0, 1] = 1
    b = S.zeros(3); b[1, 2] = 1
    assert a**2 == b**2 == S.zeros(3) and a * b != S.zeros(3)
    return {
        'completed': True, 'sympy_version': S.__version__,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'global_fixed_part_and_canonical_nontrivial_compact_controls_exact': True,
        'unit_Hom_spaces_preserved_by_compact_quotient_exact': True,
        'full_nonsemisimple_source_can_have_safe_reduced_unit_gate_exact': True,
        'unit_controls': controls,
        'local_log_square_zero_after_boundary_reduction_exact': True,
        'finite_cover_reduction_and_ambient_sum_commute_exact': True,
        'finite_cover_controls': ram,
        'genus_weight_zero_nonsemisimple_radical_image_is_ambient_exact': True,
        'local_square_zero_not_used_as_global_radical_square_zero': True,
        'universal_lowest_weight_Leray_and_boundary_reduction_proofs_machine_verified': False,
        'nonproper_Leray_MHS_theorem_independently_reproved': False,
        'full_Cauchy_quadratic_logs_claimed_absent': False,
        'reduced_actual_quotient_claimed_scalar_function_submodule': False,
        'general_algebraic_1_motive_or_curve_presentation_implemented': False,
        'all_exceptional_reduced_unit_extensions_proved_split': False,
        'pure_composites_claimed_constructed_from_nonsplit_source': False,
        'abstract_controls_claimed_actual_exceptional_phases': False,
        'all_coefficient_ring_topological_promotion_claimed': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'new_degree_support_or_prefix_campaign': False,
    }


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('haar_boundary_reduction_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
