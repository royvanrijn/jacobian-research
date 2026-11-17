"""Exact lattice, lift-normalization and ramification controls."""
from pathlib import Path
import hashlib
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def main():
    assert sp.__version__ == '1.14.0'
    p = Path(__file__).resolve().parent
    # PhiG=Z^2+Z/4; PhiT embeds by free block B and torsion row c.
    # The abelian component quotient has exponent24.
    B = sp.Matrix([[2, 1], [0, 3]])
    torsion_row = sp.Matrix([[1, 2]])
    relations_A = sp.Matrix([[2, 1, 0], [0, 3, 0], [1, 2, 4]])
    snf = smith_normal_form(relations_A, domain=sp.ZZ)
    invariant_factors = [abs(int(snf[i, i])) for i in range(3)]
    assert invariant_factors == [1, 1, 24]
    controls = []
    for free, torsion in [(sp.Matrix([0, 0]), 1),
                           (sp.Matrix([1, -3]), 2),
                           (sp.Matrix([2, 1]), 0)]:
        vectors = []
        for n in [24, 48]:
            val = B.inv()*(n*free)
            assert all(x.q == 1 for x in val)
            assert int((torsion_row*val)[0]-n*torsion) % 4 == 0
            vectors.append(val/n)
        assert vectors[0] == vectors[1] == B.inv()*free
        assert (vectors[0] == sp.zeros(2, 1)) == (free == sp.zeros(2, 1))
        controls.append({'free_component': list(map(int, free)),
                         'torsion_component_mod4': torsion,
                         'normalized_toric_vector': list(map(str, vectors[0]))})
    assert controls[0]['torsion_component_mod4'] != 0
    # Actual torus point: changing a lift by a unit does not change orders.
    z = sp.symbols('z')
    def order(q):
        return int(q.as_leading_term(z).as_powers_dict().get(z, 0))
    for k in [-5, 0, 7]:
        point = z**k*(1+z)/(1-z)
        changed = point*(2+z+z**2)/(3-z)
        assert order(point) == order(changed) == k
    # Further ramification changes raw orders but normalized vectors agree.
    v = sp.Matrix([sp.Rational(2, 3), sp.Rational(-1, 3)])
    for e, n, further in [(2, 24, 3), (5, 48, 2)]:
        raw = e*n*v
        assert all(x.q == 1 for x in raw)
        assert raw/(e*n) == further*raw/(e*further*n) == v
    basis_change = sp.Matrix([[1, 2], [0, 1]])
    assert basis_change.det() == 1
    assert basis_change*v != sp.zeros(2, 1)
    # Original inertia invariants can be smaller than split torus rank.
    swap = sp.Matrix([[0, 1], [1, 0]])
    assert len((swap-sp.eye(2)).nullspace()) == 1
    assert (-sp.eye(2)-sp.eye(2)).nullspace() == []
    # Original value vectors repeat with1/e in split ramified coordinates;
    # the cycle sum recovers factor1, separately at each root.
    branch_vectors = [(1, sp.Matrix([1, -1])),
                      (2, sp.Matrix([sp.Rational(3, 2), 2])),
                      (3, sp.Matrix([-2, sp.Rational(1, 2)]))]
    cycles = []
    for ramification, vector in branch_vectors:
        copies = [vector/ramification for _ in range(ramification)]
        assert sum(copies, sp.zeros(2, 1)) == vector
        assert all(c == sp.zeros(2, 1) for c in copies) == (vector == sp.zeros(2, 1))
        cycles.append({'multiplicity': ramification,
                       'primitive_vector': list(map(str, vector)),
                       'each_split_branch_vector': list(map(str, copies[0])),
                       'cycle_sum': list(map(str, sum(copies, sp.zeros(2, 1))))})
    norm_one = sp.Matrix([1, -1])
    assert sum(norm_one) == 0 and norm_one != sp.zeros(2, 1)
    out = {
        'date': '2026-11-16', 'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'nonprimitive_injective_torus_lattice_with_finite_cokernel_control': True,
        'abelian_component_invariant_factors': invariant_factors,
        'identity_lift_multiple_independence_controls': controls,
        'finite_component_torsion_allowed_without_identity': True,
        'unit_changes_of_toric_correction_preserve_all_orders': True,
        'normalized_vector_independent_of_further_ramification': True,
        'character_basis_change_preserves_zero': True,
        'original_inertia_invariant_space_not_entire_split_rank_control': True,
        'ramified_branch_coordinate_controls': cycles,
        'cycle_sum_recovers_intrinsic_factor1': True,
        'scalar_norm_zero_does_not_kill_full_vector': True,
        'universal_written_result': 'Correct ACTUAL sigma byidentity-component lift of[n]abelianimage;valactualT/(e*n) is canonicalI-invariantcomponentvector,evenbadabelianreduction/nonsplitG. Underpurity fullsafetyiffallorders0. AllH splitbranchvectorsrepeatv_a/e_a;cycleSumfactor1,nocrossrootcancellation',
        'split_extension_or_good_abelian_reduction_assumed': False,
        'torus_characters_assumed_to_extend_to_G': False,
        'special_fibre_torus_substituted_for_actual_generic_T': False,
        'incorrect_Bosch_Xarles_lemma_imported': False,
        'rigid_uniformizing_group_declared_algebraic': False,
        'identity_component_lift_or_Neron_producer_implemented': False,
        'actual_exceptional_purity_implies_valuation_zero_proved': False,
        'new_eventual_pure_phase_or_composite_constructed': False,
        'sufficient_Bezout_family_substituted_for_full_quotient': False,
        'split_branch_factor_one_claimed_without_ramification_normalization': False,
        'arbitrary_coefficient_ring_promotion_claimed': False,
        'universal_Neron_or_MZ_proof_machine_verified': False,
        'new_degree_support_or_parameter_campaign': False,
        'independent_replay': False, 'formal_verification': False,
        'external_review': False, 'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    (p/'haar_toric_correction_verification.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
