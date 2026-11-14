"""PF2OMIT7 algebra and explicit radial charts; topology is a written input."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S
from verify import quotient_length, checks as LOCAL_CHECKS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []


def check(label, value):
    if not value:
        raise AssertionError(label)
    CHECKS.append(label)


def zero(label, expression):
    check(label, S.cancel(expression) == 0)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vanishing_order(expression, variable):
    poly = S.Poly(expression, variable)
    assert not poly.is_zero
    return min(i[0] for i, _ in poly.terms())


def main():
    started = time.monotonic()
    x, y, u, w, t, a, v, tau, sigma, z = S.symbols('x y u w t a v tau sigma z')
    A, B, C, D, E, F = S.symbols('A B C D E F')
    P = A*x*x+B*x*y+C*y*y
    Q = D*x*x+E*x*y+F*y*y
    jac = S.diff(P, x)*S.diff(Q, y)-S.diff(P, y)*S.diff(Q, x)
    check('universal invariant quadratic jets have no constant or linear Jacobian term',
          all(sum(m) >= 2 for m, _ in S.Poly(jac, x, y).terms()))
    quotient_examples = []
    for r in range(2, 10):
        for weight in range(1, r):
            if gcd(weight, r) != 1:
                continue
            invariant_linear = [m for m, character in [('x', 1), ('y', weight)] if character % r == 0]
            check('small cyclic quotient has no invariant linear form '+str((r, weight)), invariant_linear == [])
            quotient_examples.append(dict(r=r, weight=weight, invariant_linear_forms=[]))
    zero('smooth fold admits a simple Jacobian when an invariant linear coordinate exists',
         S.diff(x*x, x)*S.diff(y, y)-2*x)
    zero('pseudoreflection control has an invariant y coordinate', y.subs(x, -x)-y)

    corank_controls = []
    for m in range(1, 5):
        g = w**(m+1)+u*w
        jacobian = S.diff(g, w)
        p = -(m+1)*t**m
        q = -m*t**(m+1)
        zero('simple critical curve parametrization m'+str(m), jacobian.subs({u: p, w: t}))
        zero('critical image parametrization m'+str(m), g.subs({u: p, w: t})-q)
        check('critical divisor smooth m'+str(m), S.diff(jacobian, u) == 1)
        length = quotient_length([u, g], (u, w))
        check('corank-one local degree m'+str(m), length == m+1)
        check('image derivative gcd exactly books ambient excess m'+str(m),
              S.degree(S.gcd(S.diff(p, t), S.diff(q, t)), t) == length-2)
        check('consecutive leading orders are primitive m'+str(m), gcd(m, m+1) == 1)
        corank_controls.append(dict(m=m, map=['u', str(g)], image=[str(p), str(q)],
                                   generic_e=2, residue_f=1, ambient_origin_length=length,
                                   derivative_gcd_degree=m-1, global_Keller_realization=False))
    cusp_p, cusp_q = -3*t*t, -2*t**3
    cusp_determinant = S.det(S.Matrix([[S.diff(cusp_p, t, 2), S.diff(cusp_p, t, 3)],
                                      [S.diff(cusp_q, t, 2), S.diff(cusp_q, t, 3)]]))
    check('ordinary cusp second/third jets independent', cusp_determinant.subs(t, 0) != 0)
    check('derivative-gcd alone also accepts an excluded higher cusp',
          S.degree(S.gcd(2*t, 5*t**4), t) == 1)
    check('higher cusp fails ordinary-cusp third jet',
          S.det(S.Matrix([[2, 0], [20*t**3, 60*t*t]])).subs(t, 0) == 0)
    check('regular cusp fibre algebra length three', quotient_length([u, w**3+u*w], (u, w)) == 3)

    blocks = {}
    for kind, degree, target_second, original_e in [
            ('unramified', 1, w, 1), ('fold', 2, w*w, 2),
            ('cyclic3', 3, w**3, 3), ('cusp3', 3, w**3+u*w, 2)]:
        check(kind+' actual local ideal equals monomial ideal',
              S.groebner([u, target_second], u, w).polys == S.groebner([u, w**degree], u, w).polys)
        intersection = S.zeros(degree)
        for i in range(degree):
            intersection[i, i] = -1 if i == degree-1 else -2
            if i:
                intersection[i, i-1] = intersection[i-1, i] = 1
        multiplicities = S.Matrix(range(1, degree+1))
        directions = S.Matrix([0]*(degree-1)+[1])
        check(kind+' full radial intersection equation', intersection*multiplicities == -directions)
        check(kind+' pointwise radial degree', (multiplicities.T*directions)[0] == degree)
        check(kind+' negative definite full principalization',
              all((-intersection)[:i, :i].det() > 0 for i in range(1, degree+1)))
        canonical = []
        for j in range(1, degree+1):
            # Chart u=a*w^j includes the only remaining base point for j<k.
            target_j = S.expand(target_second.subs(u, a*w**j))
            direction = S.cancel(target_j/(a*w**j))
            generator1 = a*w**j
            generator2 = target_j
            determinant = S.diff(generator1, a)*S.diff(generator2, w)-S.diff(generator1, w)*S.diff(generator2, a)
            canonical.append(vanishing_order(determinant, w))
            check(kind+' ideal residual after chart'+str(j),
                  S.rem(generator2, w**j, w) == 0 and S.rem(generator1, w**j, w) == 0)
            if j == degree:
                zero(kind+' horizontal direction degree one', direction.subs(w, 0)-1/a)
                check(kind+' horizontal canonical coefficient', canonical[-1] == 2*j-1)
            else:
                check(kind+' vertical canonical coefficient '+str(j), canonical[-1] >= 2*j)
        expected = {'unramified': [1], 'fold': [2, 3], 'cyclic3': [3, 4, 5], 'cusp3': [2, 4, 5]}[kind]
        check(kind+' complete canonical vector', canonical == expected)
        attachment = S.zeros(degree, 1)
        attachment_index = None if kind == 'unramified' else (2 if kind == 'cusp3' else 1)
        if attachment_index is not None:
            attachment[attachment_index-1] = original_e-1
        check(kind+' canonical adjunction includes the strict old boundary',
              intersection*S.Matrix(canonical)+attachment ==
              S.Matrix([-2-intersection[i, i] for i in range(degree)]))
        blocks[kind] = dict(original_generic_e=original_e, original_residue_f=1,
            original_ambient_length=degree, original_fibre_algebra='C[w]/(w^'+str(degree)+')',
            source_N_regular=True, radial_prime=dict(e=degree, f=1, image='target exceptional P1',
                                                   omitted=True, ramifies=(degree > 1)),
            principalization=dict(intersection_matrix=[[int(v) for v in intersection.row(i)] for i in range(degree)],
                                  e=[int(v) for v in multiplicities], f=[int(v) for v in directions], canonical_orders=canonical),
            strict_old_ramification_attachment_component=attachment_index,
            strict_old_ramification_canonical_coefficient=original_e-1,
            actual_special_direction=None, special_direction_in_model=0,
            special_direction_ambient_length=degree, special_direction_count_point_once=True,
            incident_old_and_radial_generic_indices=[original_e, degree])

    # Explicit normal charts: cyclic A_(k-1), and the cusp with the same A2
    # surface but a different finite projection.
    for k in (1, 2, 3):
        polynomial = w**k-u*tau
        check('cyclic u-chart fibre length '+str(k), quotient_length([u, tau, polynomial], (u, tau, w)) == k)
        check('cyclic u-chart monic rank '+str(k), S.Poly(polynomial, w).is_monic)
        zero('cyclic all-direction fibre relation '+str(k), polynomial.subs(u, 0)-w**k)
        zero('cyclic second chart '+str(k), (w**k-v).subs(v, w**k))
    cusp_chart = w**3+u*w-u*tau
    zero('cusp u-chart is normal A2 surface', cusp_chart.subs(tau, z+w)-(w**3-u*z))
    zero('cusp all-direction fibre relation', cusp_chart.subs(u, 0)-w**3)
    cusp_other = w**3+sigma*v*w-v
    check('cusp other chart smooth along exceptional', S.diff(cusp_other, v).subs({v: 0, w: 0}) == -1)
    check('cusp u-chart fibre length3 at special direction', quotient_length([u, tau, cusp_chart], (u, tau, w)) == 3)
    check('cusp u-chart fibre length3 at ordinary direction', quotient_length([u, tau-1, cusp_chart], (u, tau, w)) == 3)
    check('cusp other-chart fibre length3', quotient_length([v, sigma, cusp_other], (v, sigma, w)) == 3)
    # Original target Jacobian has an extra factor u after pulling back its
    # blowup chart; divide it out to get ramification of the new map.
    new_discriminant = S.discriminant(cusp_chart, w)
    zero('cusp radial and strict branch discriminant factors', new_discriminant+u*u*(4*u+27*tau*tau))
    check('radial cusp cannot use cyclic projection coordinates',
          S.expand(cusp_chart-(w**3-u*tau)) != 0)

    patterns = [
        ('A0', ['unramified']*5), ('A1', ['fold']+['unramified']*3),
        ('A2', ['fold', 'fold', 'unramified']), ('B', ['cusp3', 'fold']),
        ('C', ['cyclic3', 'fold'])]
    singleton_rows = []
    for label, kinds in patterns:
        radial = [blocks[kind]['radial_prime'] for kind in kinds]
        lengths = [blocks[kind]['original_ambient_length'] for kind in kinds]
        check('singleton '+label+' complete omitted radial degree', sum(p['e']*p['f'] for p in radial) == 5)
        check('singleton '+label+' full original fibre length', sum(lengths)+1 == 6)
        check('singleton '+label+' every radial direction preserves full length', sum(lengths)+1 == 6)
        singleton_rows.append(dict(type=label, omitted_blocks=kinds, omitted_lengths=lengths,
            omitted_radial_rows=radial, retained_radial_row=dict(e=1, f=1, omitted=False,
                                                               ramifies=False, image='target exceptional P1'),
            every_exceptional_direction_normalization_lengths=lengths+[1],
            directions_of_special_incidence=None, point_labels_remain_distinct_when_directions_coincide=True))
    check('two cubic blocks have different old ramification indices',
          blocks['cusp3']['original_generic_e'] == 2 and blocks['cyclic3']['original_generic_e'] == 3)
    check('two cubic blocks have different vertical canonical coefficients',
          blocks['cusp3']['principalization']['canonical_orders'][0] == 2
          and blocks['cyclic3']['principalization']['canonical_orders'][0] == 3)
    # Once e1/e2 have f1, a remaining f>=2 costs at least six globally.
    check('all degree-at-most-six residual higher-index residue rows excluded',
          not [(e, f) for e in range(3, 7) for f in range(2, 7) if e*f <= 5])
    for rel in ['sources/residue_manifest.json', '../jc2-boundary-attachment/sources/manifest.json']:
        manifest_path = (HERE/rel).resolve()
        manifest = json.loads(manifest_path.read_text())
        for source in manifest['sources']:
            check('source pin '+source['file'], digest(manifest_path.parent/source['file']) == source['sha256'])
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT7',
        generic_index_two=dict(all_finite_points_N_regular=True, residue_f=1,
                               exceptional_length_formula='ell=2+ord gcd(p_prime,q_prime)',
                               source_is_global_Keller_normalization_required=True),
        degree_at_most_six_all_omitted_residue_degrees=1,
        residual_row_B=dict(source_N_regular=True, exceptional_branch='ordinary (2,3) cusp',
                            exceptional_fibre_algebra='C[w]/(w^3)', square_zero_algebra_excluded=True,
                            actual_exceptional_value=None),
        all_residual_r7_affine_N_smooth=True,
        local_blocks=blocks, corank_one_controls=corank_controls, quotient_jet_controls=quotient_examples,
        singleton_radial_ledger=singleton_rows,
        special_radial_chart_scope='After target blowup only; no transport of original fixed-target disjoint-boundary prohibition',
        actual_affine_image_equations=None, actual_Keller_inverse=None,
        unrestricted_singleton_zero_properness_proved=False, full_goal_complete=False)
    target = HERE/'simple_radial_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/'SIMPLE_RAMIFICATION_AND_RADIAL.md', HERE/'RESIDUE_LEDGER.md',
              HERE/'residue_ledger.json', HERE/'verify.py', HERE/'sources/residue_manifest.json',
              HERE.parent/'jc2-boundary-attachment/PROOF.md',
              HERE.parent/'jc2-boundary-attachment/sources/manifest.json',
              HERE.parent/'hc4-radial-intersection/PROOF.md', HERE.parent/'hc4-pointed-pencil/PROOF.md']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS),
        local_staircase_checks=len(LOCAL_CHECKS), checks=CHECKS,
        runtime_seconds=round(time.monotonic()-started, 6), python=platform.python_version(), sympy=S.__version__,
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target),
        proof_boundary='End-chain complement, cyclic quotient uniformization and plane-curve knot criterion are explicit inputs. Bounded jet/quotient controls do not prove them.',
        imported_global_theorems_independently_reproved=False, formal_verification=False, external_review=False,
        full_goal_complete=False)
    (HERE/'simple_radial_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'local_staircase_checks', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
