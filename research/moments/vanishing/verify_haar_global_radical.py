"""Exact linear/period boundary controls; not universal geometric verification."""
from pathlib import Path
import hashlib
import itertools
import json
import sympy as S


def basis(mats):
    if not mats:
        return []
    piv = S.Matrix.hstack(*[x.reshape(x.rows**2, 1) for x in mats]).rref()[1]
    return [mats[i] for i in piv]


def algebra(gens):
    out = basis([S.eye(gens[0].rows)] + gens)
    while True:
        more = basis(out + [a * g for a in out for g in gens])
        if len(more) == len(out):
            return more
        out = more


def orbit(gens, g):
    out = S.Matrix.hstack(g)
    while True:
        more = S.Matrix.hstack(out, *[a * out for a in gens])
        piv = more.rref()[1]
        more = more[:, list(piv)]
        if more.cols == out.cols:
            return more
        out = more


def fixed(gens):
    n = gens[0].rows
    return S.Matrix.vstack(*[a - S.eye(n) for a in gens]).nullspace()


def lift_exists(gens, w, q):
    mats = [(a - S.eye(a.rows)) * w for a in gens]
    lhs = S.Matrix.vstack(*mats, q * w)
    rhs = S.zeros(lhs.rows, 1)
    rhs[-1] = 1
    return lhs.rank() == lhs.row_join(rhs).rank()


def run():
    # Semisimple GLOBAL Sym2(SL2) representation with local size3 blocks.
    u = S.Matrix([[1, 2, 1], [0, 1, 1], [0, 0, 1]])
    v = S.Matrix([[1, 0, 0], [1, 1, 0], [1, 2, 1]])
    assert (u - S.eye(3))**3 == S.zeros(3)
    assert (u - S.eye(3))**2 != S.zeros(3)
    gens = [S.diag(1, u), S.diag(1, v)]
    a = algebra(gens)
    assert len(a) == 10  # C direct_sum Mat3(C).
    trace = S.Matrix([[(x * y).trace() for y in a] for x in a])
    assert trace.rank() == 10
    gamma = S.Matrix([1, 1, 0, 0])
    w = orbit(gens, gamma)
    q = S.Matrix([[1, 0, 0, 0]])
    assert w.cols == 4 and lift_exists(gens, w, q)
    assert gens[0] * gamma != gamma or gens[1] * gamma != gamma
    # No purity/phase realization is claimed for this specified gamma.

    # A full split homology extension can lack a fixed lift IN its orbit.
    nh = S.zeros(4)
    nh[2, 1] = nh[3, 2] = 1
    t = S.eye(4) + nh + nh**2 / 2
    g = S.Matrix([1, 1, 0, 0])
    wg = orbit([t], g)
    assert wg.cols == 3
    assert lift_exists([t], S.eye(4), q)
    assert not lift_exists([t], wg, q)
    assert nh**2 * g != S.zeros(4, 1)

    # Separately fixed lifts need not be a simultaneous fixed lift.
    ta = S.Matrix([[1, 0], [1, -1]])
    tb = S.Matrix([[1, 0], [-1, -1]])
    q2 = S.Matrix([[1, 0]])
    assert lift_exists([ta], S.eye(2), q2)
    assert lift_exists([tb], S.eye(2), q2)
    assert not lift_exists([ta, tb], S.eye(2), q2)
    assert ta * tb * tb * ta == S.eye(2)  # infinity can be identity.

    # Finite-cover cohomology quotient identifies all ambient copies by SUM.
    ramification = []
    for e in [2, 3]:
        for name, prim in [('semisimple', S.diag(1, -1)),
                           ('nonsplit', S.Matrix([[1, 1], [0, 1]]))]:
            ind = S.zeros(2 * e)
            ind[:2, 2 * (e - 1):] = prim
            for k in range(1, e):
                ind[2 * k:2 * k + 2, 2 * (k - 1):2 * k] = S.eye(2)
            r = S.zeros(e + 1, 2 * e)
            sec = S.zeros(2 * e, e + 1)
            sec[0, 0] = 1
            for k in range(e):
                r[0, 2 * k] = 1
                r[1 + k, 2 * k + 1] = 1
                sec[2 * k + 1, 1 + k] = 1
            comp = r * ind * sec
            assert r * sec == S.eye(e + 1)
            assert r * ind == comp * r
            assert len(r.nullspace()) == e - 1
            for k in range(e):
                branch = r[:, 2 * k:2 * k + 2]
                assert branch.rank() == 2
                assert comp**e * branch == branch * prim
            if name == 'semisimple':
                assert comp**(2 * e) == S.eye(e + 1)
            else:
                delta = comp**e - S.eye(e + 1)
                assert delta.rank() == 1 and delta**2 == S.zeros(e + 1)
            ramification.append({'degree': e, 'primitive_type': name,
                                  'quotient_rank': e + 1,
                                  'ambient_kernel_rank': e - 1})

    # Weight4/5/6 upper-triangular algebra: radical cube zero, unit in J2M.
    e01 = S.zeros(3); e01[0, 1] = 1
    e12 = S.zeros(3); e12[1, 2] = 1
    e02 = e01 * e12
    upper = algebra([S.eye(3) + e01, S.eye(3) + e12,
                     S.diag(1, -1, -1), S.diag(1, 1, -1)])
    assert len(upper) == 6
    j = [e01, e12, e02]
    assert all(x * y * z == S.zeros(3) for x, y, z in itertools.product(j, repeat=3))
    assert e01 * e12 != S.zeros(3)
    unit = S.Matrix([1, 0, 0])
    jm = S.Matrix.hstack(*j)
    assert jm.row_join(unit).rank() == jm.rank()
    assert e01 * e12 * S.Matrix([0, 0, 1]) == unit
    invariant_rows = fixed([x.T for x in [S.eye(3) + e01, S.eye(3) + e12]])
    assert all((row.T * unit)[0] == 0 for row in invariant_rows)
    # The one-step boundary extension also places the unit in JM, with J2=0.
    assert e02**2 == S.zeros(3) and e02 * S.Matrix([0, 0, 1]) == unit

    # Principal parts commute with the Weyl action modulo polynomials.
    s = S.symbols('s')
    def pp0(f):
        ex = S.series(f, s, 0, 1).removeO().expand()
        return S.Add(*[x for x in S.Add.make_args(ex) if x.as_powers_dict().get(s, 0) < 0])
    r = s**2 + 3 / s + 2 / s**2 + 1 / (s - 2)
    assert S.simplify(pp0(S.diff(r, s)) - S.diff(pp0(r), s)) == 0
    assert S.simplify(pp0(s * r) - s * pp0(r)).has(s) is False
    # The earlier entire nonsplit control remains nonsplit as a Cauchy system.
    log = S.log(s / (s - 1))
    assert S.simplify((1 - s) * S.diff(log, s) - 1 / s) == 0
    return {
        'completed': True, 'sympy_version': S.__version__,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'global_semisimplicity_can_coexist_with_local_size_three_blocks': True,
        'full_fixed_lift_outside_actual_orbit_is_insufficient_exact': True,
        'separate_local_fixed_lifts_need_not_be_simultaneous_exact': True,
        'finite_cover_ambient_sum_quotient_and_primitive_injections_exact': True,
        'finite_cover_controls': ramification,
        'boundary_weight_lowering_radical_cube_zero_exact': True,
        'unit_can_lie_in_first_or_second_radical_layer_exact': True,
        'invariant_covectors_kill_radical_unit_exact': True,
        'principal_part_Weyl_retraction_controls_exact': True,
        'entire_countercontrol_Cauchy_log_retained_exact': True,
        'global_semisimple_algebra_dimension': 10,
        'global_semisimple_trace_rank': 10,
        'full_split_dimension': 4, 'restricted_orbit_dimension': 3,
        'universal_orbit_retraction_and_geometric_pure_middle_closure_machine_verified': False,
        'all_exceptional_global_unit_extensions_proved_split': False,
        'all_exceptional_generic_H2_proved_pure_weight_two': False,
        'global_semisimplicity_identified_with_local_matrix_semisimplicity': False,
        'full_extension_fixed_lift_claimed_sufficient_without_orbit_condition': False,
        'ordinary_generic_or_exponential_localization_substituted_for_Cauchy_functions': False,
        'abstract_controls_claimed_actual_exceptional_phases': False,
        'all_coefficient_ring_topological_promotion_claimed': False,
        'general_monodromy_or_orbit_producer_implemented': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'new_degree_support_or_prefix_campaign': False,
    }


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('haar_global_radical_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
