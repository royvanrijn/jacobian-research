"""Exact local controls and complete residual integer ledgers for PF2OMIT6."""
from pathlib import Path
from itertools import product
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


def order(expression, u, v):
    poly = S.Poly(expression, u, v)
    assert not poly.is_zero
    return min(sum(m) for m, _ in poly.terms())


def main():
    started = time.monotonic()
    s, t, u, v, w, h, T = S.symbols('s t u v w h T')
    controls = []
    for e in (1, 2, 3):
        for f in (1, 2, 3, 4):
            length = quotient_length([s**e, t**f], (s, t))
            check('product control local degree '+str((e, f)), length == e*f)
            controls.append(dict(map=['s^'+str(e), 't^'+str(f)], boundary='s=0', e=e, f=f,
                                 residue_local_degree=f, ambient_origin_length=length,
                                 ambient_excess=length-e, residue_cost=e*(f-1), residual_eta=0,
                                 global_Keller_boundary=False))
    # A residue cubic with two distinct simple critical points, rather than
    # one total finite ramification point.
    g = t**3-3*t
    zero('first residue collision factorization', g-g.subs(t, 1)-(t-1)**2*(t+2))
    zero('second residue collision factorization', g-g.subs(t, -1)-(t+1)**2*(t-2))
    check('two residue critical local factors cost four plus two',
          quotient_length([s**2, (t-1)**2], (s, t))+
          quotient_length([s**2, t+2], (s, t)) == 6)
    check('two finite residue ramification costs add correctly', 2*(2-1)+2*(2-1) == 2*(3-1))
    cusp_v = w**3+u*w
    critical = u+3*w*w
    zero('smooth cubic critical curve', S.diff(critical, u)-1)
    zero('cubic cusp critical parametrization', cusp_v.subs({u: -3*t*t, w: t})+2*t**3)
    check('ordinary cusp ambient length three', quotient_length([u, cusp_v], (u, w)) == 3)
    check('generic simple ramified factor length two', quotient_length([u+3, (w-1)**2], (u, w)) == 2)
    check('abstract cusp double cover ambient length two',
          quotient_length([u, v, h*h-u**3+v*v], (u, v, h)) == 2)
    check('cusp normalization nonimmersive at zero',
          S.diff(t*t, t).subs(t, 0) == S.diff(t**3, t).subs(t, 0) == 0)
    check('curvilinear length-three algebra', quotient_length([w**3], (w,)) == 3)
    check('noncurvilinear length-three algebra', quotient_length([s*s, s*t, t*t], (s, t)) == 3)

    old = json.loads((HERE/'infinity_ledger.json').read_text())
    triples = [(row['other_ramified_e'], row['other_ramified_excess'], row['unramified_omitted_primes'])
               for row in old['residual_immersed_r7_ledger']]
    check('prior residual table retained exactly', triples == [(2, 0, 1), (2, 1, 0), (3, 0, 0)])
    rows = []
    for label, (e, excess, count) in zip('ABC', triples):
        possible_f = [f for f in range(1, 7) if e*(f-1) <= excess]
        check('row'+label+' forces residue degree one', possible_f == [1])
        singleton = []
        for c, b, hh, delta in product(range(4), range(4), range(7), range(2)):
            if count == 0 and hh != 0:
                continue
            if label != 'B' and delta != 0:
                continue
            if label == 'B' and delta == 1 and b == 0:
                continue
            if 6-2*c-e*b-count*hh-delta == 1:
                singleton.append(dict(c=c, b=b, h=hh if count else None, at_exceptional_value=bool(delta)))
        if label == 'A':
            expected = sorted([(0, 0, 5), (0, 1, 3), (1, 0, 3), (0, 2, 1), (1, 1, 1), (2, 0, 1)])
            check('rowA complete singleton triples', sorted((a['c'], a['b'], a['h']) for a in singleton) == expected)
            check('rowA shared unramified image permits only1,1,1',
                  [(a['c'], a['b'], a['h']) for a in singleton if a['h'] == a['c']] == [(1, 1, 1)]
                  and [(a['c'], a['b'], a['h']) for a in singleton if a['h'] == a['b']] == [(1, 1, 1)])
        elif label == 'B':
            check('rowB complete singleton pairs are exceptional',
                  [(a['c'], a['b']) for a in singleton] == [(0, 2), (1, 1)]
                  and all(a['at_exceptional_value'] for a in singleton))
        else:
            check('rowC complete singleton pair', [(a['c'], a['b']) for a in singleton] == [(1, 1)])
        # Enumerate precisely the degree-bounded local length rows. These
        # are numeric gates; no listed integer row is asserted realizable.
        target_rows = []
        for c, b, hh, delta in product(range(4), range(4), range(7), range(2)):
            if (count == 0 and hh) or (label != 'B' and delta) or (delta and b == 0):
                continue
            source = 6-2*c-e*b-count*hh-delta
            if source < 0:
                continue
            omitted = [2]*c+[e]*b+[1]*(count*hh)
            if delta:
                omitted[c] += 1
            check('local length balance '+str((label, c, b, hh, delta)), source+sum(omitted) == 6)
            if label == 'B':
                check('rowB parity '+str((c, b, delta)), source % 2 == delta)
                if delta:
                    check('exceptional odd cardinality '+str((c, b)), source in (1, 3))
            target_rows.append(dict(c=c, b=b, h=hh if count else None,
                                    at_exceptional_value=bool(delta), omitted_local_lengths=omitted,
                                    source_length=source))
        normalized = [1, e-1]
        odd = ['C']+(['B'] if (e-1) % 2 else [])
        check('row'+label+' normalized trace discriminant weights', normalized == ([1, 2] if label == 'C' else [1, 1]))
        rows.append(dict(row=label,
            ramified_primes=[dict(name='R', e=2, f=1, image='C', image_equation=None, excess=0, immersive=True),
                            dict(name='S', e=e, f=1, image='B', image_equation=None, excess=excess,
                                 immersive_everywhere=(excess == 0), nonimmersion_parameter_count=excess,
                                 generic_ambient_length=e,
                                 exceptional_points=[] if not excess else [dict(point=None, target=None, ambient_length=3,
                                    Artin_algebra_type=None, regular_source_forces='C[s]/(s^3)')])],
            omitted_unramified_primes=[] if not count else [dict(name='D', e=1, f=1, image='H', image_equation=None,
                                                                 immersive=True, all_ambient_lengths=1,
                                                                 may_share_image_with=['C', 'B'])],
            ramified_images_distinct=True, normalized_discriminant_weights_C_B=normalized,
            raw_discriminant_odd_part=odd, singleton_branch_patterns=singleton,
            all_degree_bounded_numeric_target_rows=target_rows, numeric_rows_are_realizations=False,
            affine_N_smooth_everywhere=(excess == 0)))
    check('index-three ramified image invisible in odd discriminant part', rows[2]['raw_discriminant_odd_part'] == ['C'])

    jet_controls = []
    for pole in (1, 2, 3):
        roots = [(u+j*v)**pole for j in range(1, 6)]
        inverse = S.Poly((T-u-7*v)*S.prod(1-z*T for z in roots), T)
        orders = [order(inverse.nth(i), u, v) for i in range(2, 7)]
        disc_order = 2*sum(order(a-b, u, v) for i, a in enumerate(roots) for b in roots[i+1:])
        check('equal-pole singleton coefficient staircase '+str(pole), orders == [j*pole for j in range(1, 6)])
        check('equal-pole singleton raw discriminant bound sharp '+str(pole), disc_order == 20*pole)
        jet_controls.append(dict(pole_order=pole, coefficient_orders=orders, discriminant_order=disc_order,
                                 polynomial_index_order=disc_order//2, formal_split_only=True,
                                 global_Keller_realization=False))
    manifest = json.loads((HERE/'sources/residue_manifest.json').read_text())
    for row in manifest['sources']:
        check('primary source hash '+row['file'], digest(HERE/'sources'/row['file']) == row['sha256'])
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT6',
        hypotheses='degree6 plane Keller map with an immersed e2,f1 birational degree(6,10) gap21 prime; additional ramification allowed',
        residue_formula='E_D=e_D*(f_D-1)+sum eta_q; eta_q=ell_q-e_D*q_D(q)>=0',
        all_omitted_residue_degrees_known=True, rows=rows,
        residue_controls=controls, cubic_cusp_control_global_Keller=False,
        abstract_double_cusp_control_global_Keller_boundary=False,
        off_branch_singleton=dict(row='A', ramified_rows=[[2, 1], [2, 1]], unramified_rows=[[1, 1]],
                                 distinct_target_components=3, omitted_lengths=[1]*5,
                                 unramified_image_branches_at_zero=5, pole_order=None,
                                 exact_leading_coefficient_order='5*p', raw_discriminant_min_order='20*p',
                                 polynomial_index_min_order='10*p', weighted_controls=jet_controls),
        actual_ramified_image_B_equation=None, actual_unramified_image_H_equation=None,
        actual_exceptional_value=None, actual_Keller_inverse=None,
        unrestricted_singleton_zero_properness_proved=False, full_goal_complete=False)
    target = HERE/'residue_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/'RESIDUE_LEDGER.md', HERE/'GLOBAL_LEDGER.md', HERE/'R7_INFINITY_COMPLETION.md',
              HERE/'infinity_ledger.json', HERE/'DISCRIMINANT_LEDGER.md', HERE/'verify.py',
              HERE/'sources/residue_manifest.json',
              ROOT/'foundations/plane-jc/OREVKOV_RESIDUE_DEGREE_BUDGET.md']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), local_staircase_checks=len(LOCAL_CHECKS),
        checks=CHECKS, runtime_seconds=round(time.monotonic()-started, 6),
        python=platform.python_version(), sympy=S.__version__, checker_sha256=digest(Path(__file__)),
        input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs}, ledger_sha256=digest(target),
        visual_primary_source_review='Printed pages7--9 rendered with Poppler and inspected; Lemma4.2, Lemma5.2 and contraction scope verified',
        imported_theorems_independently_reproved=False, formal_verification=False, external_review=False,
        proof_scope='Flat local degree and Riemann-Hurwitz are written arguments; Orevkov5.2 is an explicit global import. Finite controls do not prove these statements.',
        full_goal_complete=False)
    (HERE/'residue_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'local_staircase_checks', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
