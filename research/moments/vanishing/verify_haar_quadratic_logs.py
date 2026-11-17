"""Exact transcription controls; universal geometric proofs are in the note."""
from pathlib import Path
import hashlib
import json
import sympy as S


def span_equal(a, b):
    return a.rank() == b.rank() == a.row_join(b).rank()


def cyclic(a, e):
    n = a.rows
    out = S.zeros(n * e)
    out[:n, (e - 1) * n:] = a
    for j in range(1, e):
        out[j * n:(j + 1) * n, (j - 1) * n:j * n] = S.eye(n)
    return out


def run():
    # A relative-limit model, NOT an exceptional phase realization.
    # compact J3; ambient L; genus J2; boundary Tate mapping into L.
    n = S.zeros(7)
    n[0, 1] = n[1, 2] = n[4, 5] = n[3, 6] = 1
    w = [4, 4, 4, 4, 5, 5, 6]
    m = [2, 4, 6, 4, 4, 6, 6]
    p = [1, 2, 3, 2, 2, 3, 3]
    for i in range(7):
        for j in range(7):
            if n[i, j]:
                assert w[i] <= w[j]
                assert m[i] == m[j] - 2 and p[i] == p[j] - 1
    assert n**3 == S.zeros(7) and (n**2).rank() == 1
    m2 = S.eye(7)[:, :1]
    w4 = S.eye(7)[:, :4]
    assert span_equal(n**2, m2)
    assert span_equal(n**2, n**2 * w4)
    top = S.eye(7)[:, [j for j in range(7) if m[j] == 6]]
    assert top.cols == 3 and (n**2 * top).rank() == 1
    ambient = S.eye(7)[:, 3:4]
    assert (n**2).row_join(ambient).rank() == 2
    keep = [0, 1, 2, 4, 5, 6]
    nk = n.extract(keep, keep)
    assert (nk**2).rank() == (n**2).rank()
    assert (nk[:3, :3]**2).rank() == (n**2).rank()

    # Split intrinsic extension, but the specified lift has quadratic logs.
    split = n.copy()
    split[3, 6] = 0
    nh = -split.T
    g = S.eye(7)[:, 3] + S.eye(7)[:, 0]
    assert nh * S.eye(7)[:, 3] == S.zeros(7, 1)
    assert nh * g != S.zeros(7, 1) and nh**2 * g != S.zeros(7, 1)
    t = S.eye(7) + nh + nh**2 / 2
    assert (t - S.eye(7))**3 == S.zeros(7)
    assert nh == t - S.eye(7) - (t - S.eye(7))**2 / 2
    assert nh**2 == (t - S.eye(7))**2
    for j in [-3, 0, 1, 2, 5]:
        assert t**j * g == g + j * nh * g + S.Rational(j * j, 2) * nh**2 * g
    z, ell = S.symbols('z ell')
    h = S.eye(7)
    h[0, 1] = z
    assert h.det() == 1
    period = h * (g + ell * nh * g + ell**2 * nh**2 * g / 2)
    assert period.diff(ell, 3) == S.zeros(7, 1)
    assert period.diff(ell, 2) == h * nh**2 * g
    assert period.diff(ell).subs(ell, 0) == h * nh * g
    assert period[3] == 1  # One row misses all motion; no phase purity claim.

    # Quadratic-log vanishing alone leaves genuine first-order motion.
    g1 = S.eye(7)[:, 3] + S.eye(7)[:, 4]
    assert nh * g1 != S.zeros(7, 1) and nh**2 * g1 == S.zeros(7, 1)
    # A length4 chain needs four distinct p levels and cannot fit1..3.
    assert all(not all(1 <= start - j <= 3 for j in range(4))
               for start in range(-3, 8))

    # Cyclic ramification keeps the bound for ALL eigenvalue sectors.
    j3 = S.zeros(3)
    j3[0, 1] = j3[1, 2] = 1
    a1 = S.eye(3) + j3 + j3**2 / 2
    a = S.diag(a1, S.Matrix([[-1, -1], [0, -1]]))
    ram = []
    for e in [2, 3]:
        ap = cyclic(a, e)
        cover = ap**(2 * e)
        d = cover - S.eye(ap.rows)
        assert d**3 == S.zeros(ap.rows)
        log = d - d**2 / 2
        assert log**3 == S.zeros(ap.rows)
        assert (log**2).rank() == e
        ram.append({'root_multiplicity': e, 'all_sector_log_cube_zero': True,
                    'size_three_blocks_after_full_cover': e})
    return {
        'completed': True, 'sympy_version': S.__version__,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'Hodge_strip_and_relative_weight_matrix_control_exact': True,
        'image_N_squared_equals_M2_and_compact_middle_image_exact': True,
        'ambient_line_disjoint_from_image_N_squared_exact': True,
        'full_top_limit_grade_to_bottom_need_not_be_injective_exact': True,
        'split_extension_specified_lift_can_have_quadratic_logs': True,
        'quadratic_log_zero_does_not_force_first_log_zero': True,
        'original_one_block_log_and_all_integer_continuations_exact': True,
        'finite_meromorphic_row_basis_detects_both_log_coefficients_exact': True,
        'one_row_invisible_motion_not_promoted_to_actual_phase': True,
        'ramification_controls': ram,
        'universal_geometric_limit_and_marker_basis_proofs_machine_verified': False,
        'meromorphic_log_coefficients_claimed_algebraic': False,
        'N_squared_claimed_isomorphism_on_full_top_limit_grade': False,
        'general_original_basis_or_canonical_frame_producer_implemented': False,
        'all_exceptional_actual_logarithms_proved_zero': False,
        'arbitrary_coefficient_ring_topological_promotion_claimed': False,
        'abstract_models_claimed_actual_exceptional_phases': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'new_degree_support_or_prefix_campaign': False,
    }


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('haar_quadratic_logs_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
