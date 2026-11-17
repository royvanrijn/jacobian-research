"""Exact finite controls for the internal intrinsic-source safety proof.

This does not machine-prove fixed part, admissibility, or Borel splitting.
Run with the repository's SymPy 1.14.0 environment; writes only its receipt.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp


def main():
    here = Path(__file__).resolve().parent
    assert sp.__version__ == '1.14.0'
    # Semisimple fixed/dual pairing and canonical evaluation projection.
    a = sp.diag(1, 1, -1, -1)
    b = sp.diag(1, 1, 2, sp.Rational(1, 2))
    fixed = sp.Matrix.vstack(a-sp.eye(4), b-sp.eye(4)).nullspace()
    dual = sp.Matrix.vstack(a.inv().T-sp.eye(4), b.inv().T-sp.eye(4)).nullspace()
    F, Fdual = sp.Matrix.hstack(*fixed), sp.Matrix.hstack(*dual)
    pairing = Fdual.T*F
    projection = F*pairing.inv()*Fdual.T
    assert pairing.det() != 0 and projection**2 == projection
    assert projection*F == F and projection*a == a*projection
    C = sp.Matrix.hstack(a-sp.eye(4), b-sp.eye(4)).columnspace()
    assert len(C) == 2 and all(projection*x == sp.zeros(4, 1) for x in C)

    # In the residual VMHS, relative M=W permits only weight drop >=2.
    # Order: unit4, two constant genus5, two Tate6.
    weights = [4, 5, 5, 6, 6]
    variables = sp.symbols('a0:2 b0:2 c0:2')
    local_N = sp.zeros(5)
    for j in range(2):
        local_N[0, 3+j] = variables[j]
        local_N[1, 3+j] = variables[2+j]
        local_N[2, 3+j] = variables[4+j]
    forbidden = [local_N[i, j] for i in range(5) for j in range(5)
                 if weights[i] > weights[j]-2 and local_N[i, j] != 0]
    assert set(forbidden) == set(variables[2:])
    allowed = local_N.subs(dict.fromkeys(forbidden, 0))
    p1, p2 = sp.symbols('p1 p2')
    retraction = sp.Matrix([[1, p1, p2, 0, 0]])
    equations = list(retraction*allowed)
    solution = sp.solve(equations, variables[:2], dict=True)
    assert solution == [{variables[0]: 0, variables[1]: 0}]
    assert allowed.subs(solution[0]) == sp.zeros(5)

    # Earlier obstructing local translation violates that weight-drop rule.
    bad_N = sp.zeros(5)
    bad_N[1, 3] = 1
    unit = sp.eye(5)[:, 0]
    radical = sp.eye(5)[:, 1]
    relation = radical-unit
    assert sp.Matrix.hstack(radical, relation).rank() == 2
    assert sp.Matrix.hstack(radical, relation, unit).rank() == 2
    assert bad_N[1, 3] != 0 and weights[1] > weights[3]-2

    # Finite LOCAL monodromy does not imply finite GLOBAL monodromy.
    # Four involutions have product I, but their global algebra is nonsplit
    # after a unit-retaining abstract quotient. No actual pure phase asserted.
    s1 = sp.diag(1, 1, 1, -1)
    s2 = s1.copy()
    s2[1, 3] = 2
    source = [s1, s2, s2, s1]
    q = sp.Matrix([[1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    n1 = sp.diag(1, 1, -1)
    n2 = n1.copy()
    n2[0, 2] = 2
    actual_abstract = [n1, n2, n2, n1]
    assert all(x*x == sp.eye(4) for x in source)
    assert all(x*x == sp.eye(3) for x in actual_abstract)
    assert source[0]*source[1]*source[2]*source[3] == sp.eye(4)
    assert all(q*x == y*q for x, y in zip(source, actual_abstract))
    assert q*sp.Matrix([-1, 1, 0, 0]) == sp.zeros(3, 1)
    assert q*sp.eye(4)[:, 0] == sp.eye(3)[:, 0]
    assert (source[0]*source[1]-sp.eye(4))**2 == sp.zeros(4)
    assert source[0]*source[1] != sp.eye(4)
    ell1 = sp.Matrix([[1, 0, 0]])
    ell2 = sp.Matrix([[1, 0, 1]])
    assert ell1*n1 == ell1 and ell2*n2 == ell2
    assert (n1*n2-sp.eye(3))[:, 2] == sp.Matrix([2, 0, 0])
    assert ell1*sp.eye(3)[:, 0] == sp.ones(1)
    assert ell2*sp.eye(3)[:, 0] == sp.ones(1)
    x, y = sp.symbols('x y')
    ell = sp.Matrix([[1, x, y]])
    assert sp.solve(list(ell*n1-ell)+list(ell*n2-ell), [x, y]) == []

    # Principal part intertwines s and d/ds after polynomial quotient.
    s = sp.symbols('s')
    coeffs = sp.symbols('c0:9')
    laurent = sum(c*s**k for c, k in zip(coeffs, range(-4, 5)))

    def pp(expr):
        return sp.Add(*(term for term in sp.expand(expr).as_ordered_terms()
                        if term.as_powers_dict().get(s, 0) < 0))

    assert pp(sp.diff(laurent, s)) == sp.diff(pp(laurent), s)
    assert pp(s*laurent) == pp(s*pp(laurent))
    assert pp(sum(coeffs[i]*s**i for i in range(9))) == 0
    for k in range(1, 5):
        assert pp(s**(-k)) == s**(-k)

    output = {
        'date': '2026-11-16',
        'completed': True,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sympy_version': sp.__version__,
        'semisimple_fixed_and_dual_pairing_projection_exact': True,
        'relative_weight_drop_and_split_unit_force_local_N_zero_exact': True,
        'earlier_obstructing_local_translation_is_not_admissible_exact': True,
        'finite_local_involutions_can_have_infinite_global_translation_exact': True,
        'abstract_actual_kernel_retains_unit_but_destroys_global_splitting_exact': True,
        'local_unit_covectors_survive_that_abstract_quotient_exact': True,
        'principal_part_Weyl_quotient_and_polynomial_kernel_controls_exact': True,
        'mixed_fixed_part_import': 'Steenbrink--Zucker1985 Proposition4.19, hypotheses3.13; original scanned pages visually read; André1992 Section5 context',
        'new_universal_argument': 'Admissible residual source has finite LOCAL monodromy; compose actual Cauchy quotient with a local meromorphic functional, then PP0, to split actual Borel unit under eventual purity',
        'actual_kernel_needs_Hodge_compatibility': False,
        'actual_quotient_embedded_as_scalar_functions': False,
        'finite_local_monodromy_used_to_infer_finite_global_monodromy': False,
        'global_source_unit_splitting_proved_for_every_exceptional_pure_phase': False,
        'constant_genus_rank_or_global_R_claimed_zero': False,
        'arbitrary_abstract_K_claimed_harmless_without_eventual_purity': False,
        'all_coefficient_ring_topological_promotion_claimed': False,
        'universal_fixed_part_admissibility_Borel_proofs_machine_verified': False,
        'new_degree_support_or_moment_prefix_campaign': False,
        'independent_replay': False,
        'formal_verification': False,
        'external_review': False,
        'external_review_next_task': False,
        'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    }
    target = here/'haar_intrinsic_safety_verification.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
