"""Two exact first-image matrices and complete residual singular-value matching."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as s
from image_norm import families, norm_matrix, matching_rows

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    U, V, a = s.symbols('U V a')
    definitions, models = [], {}
    for name, parameter, p, q, t in families():
        matrix = norm_matrix(p, q, t, U, V)
        pp = s.Poly(p, t)
        check(name+' chart has only its declared coefficient parameter',
              p.free_symbols | q.free_symbols == {parameter, t})
        companion = s.zeros(6)
        for j in range(5):
            companion[j+1, j] = 1
        for i in range(6):
            companion[i, 5] = (U if i == 0 else 0)-pp.nth(i)
        check(name+' fixed monic normalization degrees', pp.degree() == 6 and pp.LC() == 1
              and s.Poly(q, t).degree() == 10 and s.Poly(q, t).LC() == 1)
        check(name+' exact first column by independent monic polynomial division',
              s.expand(sum(matrix[i, 0]*t**i for i in range(6))-s.rem(q-V, p-U, t)) == 0)
        check(name+' whole multiplication matrix commutes with t',
              all(s.expand(x) == 0 for x in matrix*companion-companion*matrix))
        check(name+' target V enters as minus identity', matrix.diff(V) == -s.eye(6))
        model = (parameter, p, q, t, matrix)
        models[name] = model
        for parameter_value in ([0, 1] if name == 'L' else [1, 2]):
            point_p = s.Poly(p.subs(parameter, parameter_value), t)
            point_q = s.Poly(q.subs(parameter, parameter_value), t)
            point_matrix = matrix.subs(parameter, parameter_value)
            for target_U, target_V in [(0, 1), (1, 0), (-2, 3)]:
                specialized = point_matrix.subs({U: target_U, V: target_V})
                f, g = point_p.as_expr()-target_U, point_q.as_expr()-target_V
                check(f'{name}{parameter_value} norm/resultant {target_U},{target_V}',
                      specialized.det() == s.resultant(f, g, t))
                check(f'{name}{parameter_value} corank/gcd {target_U},{target_V}',
                      6-specialized.rank() == s.Poly(s.gcd(f, g), t).degree())
            target_U = point_p.eval(0)
            target_V = point_q.eval(0)
            gcd = s.Poly(s.gcd(point_p.as_expr()-target_U, point_q.as_expr()-target_V), t)
            check(f'{name}{parameter_value} immersed origin has one normalization point', gcd.degree() == 1)
            check(f'{name}{parameter_value} origin image matrix rank5',
                  point_matrix.subs({U: target_U, V: target_V}).rank() == 5)
        domain = ('3*d**3+2 !=0' if name == 'L' else 'x*(3*x**3+2)*(15*x**3-4) !=0')
        definitions.append(dict(family=name, parameter=str(parameter), immersed_domain=domain,
            normalized_target_coordinates=['U', 'V'], p=str(p), q=str(q),
            multiplication_basis=['1', 't', 't^2', 't^3', 't^4', 't^5'],
            multiplication_matrix_entries=[[str(matrix[i, j]) for j in range(6)] for i in range(6)],
            image_equation='det(multiplication_matrix)=0; reduced irreducible, monic degree6 in V',
            projective_image_degree=10,
            first_omitted_prime=dict(name='R', e=2, f=1, ramifies=True, omitted=True,
                                     branch_count='6-rank(matrix_at_target)',
                                     ambient_length_per_normalization_point=2,
                                     total_omitted_length='2*(6-rank(matrix_at_target))'),
            actual_global_family_parameter=None, actual_original_target_coordinate_changes=None))

    d, p, q, t, matrix = models['L']
    # A rational double-image value: H3(v,d)=0 at v=1,d=-8/9.
    d_value = -s.Rational(8, 9)
    pair = t*t-t+s.Rational(10, 27)
    node_U = s.rem(p.subs(d, d_value), pair, t)
    node_V = s.rem(q.subs(d, d_value), pair, t)
    check('rational collision target coordinates are constant', not node_U.has(t) and not node_V.has(t))
    check('rational collision has two distinct normalization parameters', s.discriminant(pair, t) != 0)
    check('rational collision matrix rank4', matrix.subs({d: d_value, U: node_U, V: node_V}).rank() == 4)
    check('rational collision normalization gcd exactly quadratic',
          s.Poly(s.gcd(p.subs(d, d_value)-node_U, q.subs(d, d_value)-node_V), t).monic().as_expr() == pair)

    # All three triple values at once, over a squarefree algebra; no root sampling.
    triple_base = 81*a**3+135*a**2+63*a+7
    monic_base = s.Poly(triple_base, a).monic().as_expr()
    triple_matrix = matrix.subs({d: 0, U: a*a+a, V: 0}).applyfunc(lambda x: s.rem(s.expand(x), monic_base, a))
    kernel = s.zeros(6, 3)
    for j in range(3):
        kernel[j, j] = a+1
        kernel[j+3, j] = 1
    check('triple base squarefree', s.gcd(triple_base, s.diff(triple_base, a)) == 1)
    check('three triple target values remain distinct', s.Poly(s.gcd(triple_base, triple_base.subs(a, -1-a)), a).degree() == 0)
    check('three independent kernel vectors at every triple value', kernel[3:6, :] == s.eye(3)
          and all(s.rem(s.expand(x), monic_base, a) == 0 for x in triple_matrix*kernel))
    witness = None
    for rows in combinations(range(6), 3):
        for cols in combinations(range(6), 3):
            minor = s.rem(s.expand(triple_matrix.extract(rows, cols).det()), monic_base, a)
            if s.Poly(s.gcd(minor, triple_base), a).degree() == 0:
                witness = dict(rows=list(rows), columns=list(cols), minor_mod_base=str(minor))
                break
        if witness:
            break
    check('uniform rank at least3 over all triple roots', witness is not None)
    triple_targets = s.Poly(s.resultant(triple_base, U-a*a-a, a), U).clear_denoms()[1].primitive()[1]
    check('triple target eliminant has three distinct roots', triple_targets.degree() == 3
          and s.Poly(s.gcd(triple_targets.as_expr(), s.diff(triple_targets.as_expr(), U)), U).degree() == 0)

    matching = []
    for c in range(7):
        for group in ['A_distinct', 'A_H=C', 'A_H=B', 'B', 'C']:
            completions = matching_rows(c, group)
            if c >= 4:
                check(f'{group} forbids {c} C branches', completions == [])
            matching.append(dict(C_branches=c, C_matrix_rank=6-c,
                                 grouping=group, necessary_completions=completions))
    for c in [2, 3]:
        check(f'rowC second image avoids {c}-branch C values',
              all(r['b'] == 0 for r in matching_rows(c, 'C')))
        check(f'rowB cusp avoids {c}-branch C values',
              all(not r['cusp'] for r in matching_rows(c, 'B')))
        for group in ['A_H=C', 'A_H=B']:
            check(f'{group} second image avoids {c}-branch C values',
                  all(r['b'] == 0 for r in matching_rows(c, group)))
    check('rowC each smooth C/B intersection is a singleton',
          [r['retained_length'] for r in matching_rows(1, 'C') if r['b'] > 0] == [1])
    check('rowB cusp on smooth C is a singleton',
          [r['retained_length'] for r in matching_rows(1, 'B') if r['cusp']] == [1])
    check('H=C impossible at any triple C value', matching_rows(3, 'A_H=C') == [])
    check('all other triple-value completions contain no other omitted branch',
          all(r['b'] == 0 and r['h'] == 0 and not r['cusp'] and r['retained_length'] == 0
              for group in ['A_distinct', 'A_H=B', 'B', 'C'] for r in matching_rows(3, group)))
    strata_path = ROOT/'jacobian/jc2-r7-collisions/strata.json'
    strata = json.loads(strata_path.read_text())
    factors = {k: s.sympify(v, locals={'d': d}) for k, v in strata['exceptional_factors'].items()}
    excluded = s.expand(d*factors['U']*factors['M'])
    check('H=C excludes13 distinct polynomial-family parameters', s.degree(excluded, d) == 13
          and s.Poly(s.gcd(excluded, s.diff(excluded, d)), d).degree() == 0)
    check('all13 sharing-excluded parameters are immersed', s.Poly(s.gcd(excluded, factors['C']), d).degree() == 0)
    for row in strata['singularities']:
        if row['cusps']:
            continue
        check(row['stratum']+' complete immersed delta12 catalogue',
              row['nodes']+2*row['tacnodes']+3*row['triple_points'] == 12)
    # Immersion is essential for translating cokernel length to branch count.
    cusp_matrix = s.Matrix([[0, 0], [0, 0]])
    check('nonimmersed cusp countercontrol has corank2 but one branch',
          s.Poly(s.gcd(t*t, t**3), t).degree() == 2 and 2-cusp_matrix.rank() == 2)

    ledger = dict(as_of='2026-11-14', claim='PF2OMIT11',
        hypotheses='PF2OMIT6 residual immersed r7 degree6; first image normalized by PF2R7FAM1/2',
        first_image_families=definitions,
        branch_count_formula='6-rank of explicit six-by-six multiplication matrix',
        first_image_singular_locus_set_theoretic='rank<=4 on each immersed chart',
        forbidden_for_Keller_first_image='rank<=2 at any target value',
        singular_value_matching=matching,
        polynomial_family_immersed_strata=[r for r in strata['singularities'] if r['cusps'] == 0],
        H_equals_C_excluded_polynomial_parameter_factor=str(excluded),
        H_equals_C_excluded_polynomial_parameter_count=13,
        triple_values_at_d0=dict(target_V=0, target_U_polynomial=str(triple_targets.as_expr()),
                                 normalization_cube_parameter_polynomial=str(triple_base),
                                 uniform_rank3_minor=witness,
                                 R_omitted_local_lengths=[2, 2, 2], retained_length=0,
                                 every_other_omitted_image_avoids=True),
        rational_double_value_control=dict(parameter=str(d_value), U=str(node_U), V=str(node_V),
                                           normalized_parameter_equation=str(pair), C_matrix_rank=4),
        actual_second_ramified_image_equation=None, actual_optional_unramified_image_equation=None,
        actual_global_inverse=None, actual_pole_orders=None,
        Laurent_family_equisingularity_classification_complete=False,
        global_Keller_realization_supplied=False, full_goal_complete=False)
    target = HERE/'image_norm_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/n for n in ['IMAGE_NORM_LEDGER.md', 'image_norm.py', 'RESIDUE_LEDGER.md',
              'SIMPLE_RAMIFICATION_AND_RADIAL.md']]
    inputs += [ROOT/'jacobian/jc2-r7-image-family'/n for n in ['family.py', 'families.json', 'PROOF.md']]
    inputs += [ROOT/'jacobian/jc2-r7-collisions'/n for n in ['PROOF.md', 'strata.json']]
    inputs += [ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target), runtime_seconds=round(time.monotonic()-started, 6),
        python=platform.python_version(), sympy=s.__version__,
        proof_boundary='Complete symbolic multiplication matrices, uniform triple-rank certificate and matching controls checked. Universal norm/rank, imported family completeness/immersion/collisions and global finite-flat interpretation remain written proofs. No global source cover is constructed.',
        formal_verification=False, external_review=False, full_goal_complete=False)
    (HERE/'image_norm_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))
    print(json.dumps(ledger['triple_values_at_d0']))


if __name__ == '__main__':
    main()
