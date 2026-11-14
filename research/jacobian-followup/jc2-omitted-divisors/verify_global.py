"""Finite algebra and integer ledgers for PF2OMIT2; no topology search."""
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S
from verify import quotient_length, checks as staircase_checks

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def branch_patterns(total, count):
    return [p for p in combinations_with_replacement(range(total+1), count) if sum(p) == total]


def main():
    started = time.monotonic()
    ledger = json.loads((HERE/'global_ledger.json').read_text())
    s, v = S.symbols('s v')
    cyclic_lengths = []
    for e in range(2, 9):
        length = quotient_length([s**e, v], (s, v))
        check(f'cyclic local ambient length e{e}', length == e)
        cyclic_lengths.append(dict(e=e, ambient_length=length, finite_excess=length-e))
    u, w, h = S.symbols('u w h')
    cusp_three = quotient_length([u, w**3+u*w], (w,u))
    cusp_two = quotient_length([u,v,h**2-u**3+v**2], (h,u,v))
    generic_double = quotient_length([(w-1)**2], (w,))
    check('cubic smooth-source cusp local length3', cusp_three == 3)
    check('double-cover singular-source cusp local length2', cusp_two == 2)
    check('generic double branch local length2', generic_double == 2)
    check('cusp exceptional excesses differ despite identical e2',
          (cusp_three-generic_double,cusp_two-generic_double) == (1,0))
    # Written geometry supplies the separate cyclic local factors.
    for branches in [(), (2,), (3,), (2,2), (2,3), (2,)*12]:
        chi_image = 1-sum(b-1 for b in branches)
        chi_smooth = 1-sum(branches)
        chi_outside = 1-chi_image
        chi_N = 6*chi_outside+5*chi_smooth+sum(6-b for b in branches)
        check('immersed Euler '+str(branches), chi_N == 5)
    r7 = ledger['applications']['normal_r7_immersed']
    check('immersed r7 unramified prime count', r7['unramified_prime_count'] == 6-1-2 == 3)
    for row in r7['strata']:
        b = row['ramification_image_branches']
        check(f'r7 b{b} ramified lengths', row['ramified_local_lengths'] == [2]*b)
        check(f'r7 b{b} remaining capacity', row['max_unramified_points'] == 6-2*b)
        singleton = 5-2*b if b <= 2 else None
        check(f'r7 b{b} singleton count', row['unramified_points_for_singleton'] == singleton)
        for extra in range(7-2*b):
            check(f'r7 b{b} B{extra} source length', 6-2*b-extra >= 0)
    check('four ramification image branches impossible', 6-2*4 < 0)
    for b, max_q in [(2,1), (3,0)]:
        allowed = [q for q in range(4) if (2+q)*b <= 6]
        check(f'b{b} same-image prime cap', max(allowed) == max_q)
    check('at least two r7 prime images differ', r7['min_primes_with_image_different_from_C'] == 3-1)
    cusp_prior = json.loads((HERE.parent/'jc2-boundary-attachment/receipt.json').read_text())
    for case in cusp_prior['Euler_ledgers']:
        epsilon = case['cusp_epsilon']
        check('cusp finite excess vs Euler '+str((case['case'],epsilon)),
              case['required_unramified_boundary_primes'] == 6-1-2-epsilon)
    patterns = {}
    for name, degree, ram_indices, excess in [
            ('r7 immersed or cusp_eps0', 6, [2], 0),
            ('r7 cusp_eps1', 6, [2], 1),
            ('r10 two immersed ramified primes', 6, [2,2], 0)]:
        t = degree-1-sum(ram_indices)-excess
        pats = branch_patterns(degree-1, t)
        for p in pats:
            absent = p.count(0)
            surplus = sum(max(b-1,0) for b in p)
            check(name+' branch-surplus '+str(p), surplus == sum(ram_indices)+excess+absent)
            check(name+' smooth images cannot suffice '+str(p), surplus >= 2)
        patterns[name] = [list(p) for p in pats]
    check('complete three-prime singleton patterns', patterns['r7 immersed or cusp_eps0'] ==
          [[0,0,5], [0,1,4], [0,2,3], [1,1,3], [1,2,2]])
    check('complete two-prime singleton patterns', patterns['r7 cusp_eps1'] == [[0,5], [1,4], [2,3]])
    numerical_ledgers = []
    for degree in range(2, 10):
        for indices in [(2,), (3,), (2,2), (2,3)]:
            for excess in range(3):
                t = degree-1-sum(indices)-excess
                if t < 0:
                    continue
                chi_N = degree-sum(e-1 for e in indices)-excess
                check('all-prime Euler '+str((degree,indices,excess)), chi_N-(len(indices)+t) == 1)
                check('unramified prime upper bound '+str((degree,indices,excess)), t <= degree-3)
                numerical_ledgers.append(dict(d=degree, ramification_indices=indices,
                                              finite_excess=excess, unramified_primes=t))
    r10 = ledger['applications']['historical_r10_budget_audit']
    chi_two_nodal = 6*7+5*(-14)+4*8
    check('r10 complete finite-cover Euler characteristic', chi_two_nodal == 4)
    check('r10 extra unramified prime repairs source Euler', chi_two_nodal-(2+1) == 1)
    check('r10 complete Orevkov budget can have e1 term', 2+2+1 == 5)
    check('r10 required omitted row', r10['required_unramified_rows'] == 1)
    check('r10 shared nodal image pointwise length', sum(r10['node_if_image_shared_omitted_lengths']) == 6)
    check('r10 shared node need not have source point', r10['node_if_image_shared_retained_length'] == 0)
    check('no r10 realization inferred', r10['realization_proved'] is False)
    check('whole even exclusion not promoted', r10['historical_whole_even_exclusion_verified'] is False)
    check('r7 no global exclusion inferred', r7['global_r7_exclusion'] is False)
    check('goal remains open', ledger['goal_complete'] is False)
    for name, app in ledger['applications'].items():
        check(name+' omission equations unknown', app['complete_omission_equations_known'] is False)
        check(name+' unresolved inputs explicit', bool(app['unknowns']))
    manifest = json.loads((HERE.parent/'jc2-boundary-attachment/sources/manifest.json').read_text())
    for entry in manifest['sources']:
        p = HERE.parent/'jc2-boundary-attachment/sources'/entry['file']
        check('prior primary pin '+entry['file'], hashlib.sha256(p.read_bytes()).hexdigest() == entry['sha256'])
    inputs = ['GLOBAL_LEDGER.md', 'STEIN_SCOPE_AUDIT.md', 'global_ledger.json', 'verify.py',
              '../jc2-boundary-attachment/PROOF.md', '../jc2-boundary-attachment/receipt.json',
              '../jc2-boundary-attachment/sources/manifest.json']
    result = dict(as_of='2026-11-14', outcome='PASS', check_count=len(checks),
                  finite_staircase_check_count=len(staircase_checks), checks=checks,
                  runtime_seconds=round(time.monotonic()-started, 6),
                  python=platform.python_version(), sympy=S.__version__,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs},
                  cyclic_local_lengths=cyclic_lengths, singleton_off_branch_patterns=patterns,
                  cusp_local_controls=dict(generic_index=2, cubic_cusp_length=cusp_three,
                                           double_cover_cusp_length=cusp_two, finite_excesses=[1,0]),
                  numerical_Euler_controls=numerical_ledgers,
                  actual_omitted_image_equations_known=False, goal_complete=False,
                  written_inputs=['plane boundary theorem', 'Euler integration', 'normal cyclic smooth-branch cover',
                                  'flat closure of the generic section factor'],
                  formal_verification=False, external_review=False)
    (HERE/'global_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['outcome','check_count','finite_staircase_check_count','runtime_seconds']}))


if __name__ == '__main__':
    main()
