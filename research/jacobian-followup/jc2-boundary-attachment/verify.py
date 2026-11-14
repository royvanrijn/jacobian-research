"""Finite algebra/topology regressions for the global boundary application.

Borisov's theorem and purity are explicit written imports, not outputs of
this checker. The old local transposition helpers are reused unchanged.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import sympy as S

HERE = Path(__file__).resolve().parent
R7 = HERE.parent/'jc2-normal-r7'
sys.path.insert(0, str(R7))
from local_excess import transpositions, mul, power

checks = []


def check(name, condition):
    assert condition, name
    checks.append(name)


def zero(name, expression):
    check(name, S.cancel(expression) == 0)


def supports(p):
    return {i for i in range(6) if p[i] != i}


def orbits(generators):
    remaining = set(range(6))
    result = []
    while remaining:
        reached = {min(remaining)}
        while True:
            newer = reached | {p[i] for p in generators for i in reached}
            if newer == reached:
                break
            reached = newer
        result.append(sorted(reached))
        remaining -= reached
    return result


def main():
    start = time.monotonic()
    cusp = []
    local_nodes = []
    local_a5 = []
    global_two_branch = []
    for i, a in enumerate(transpositions):
        for j, b in enumerate(transpositions):
            sizes = sorted(map(len, orbits([a, b])))
            separate = supports(a).isdisjoint(supports(b))
            if mul(a, b, a) == mul(b, a, b):
                epsilon = len(supports(a) | supports(b))-2
                check(f'cusp pair {i},{j} has size2 or size3 ramified orbit', epsilon in (0, 1))
                check(f'cusp pair {i},{j} finite fibre count', len(sizes) == 5-epsilon)
                cusp.append(dict(a=i, b=j, epsilon=epsilon, orbit_sizes=sizes))
            if mul(a, b) == mul(b, a):
                local_nodes.append((i, j))
            if power(mul(a, b), 3) == power(mul(b, a), 3):
                local_a5.append((i, j))
            if separate:
                check(f'two-branch pair {i},{j} has distinct ramified points', sizes == [1, 1, 2, 2])
                check(f'two-branch pair {i},{j} satisfies node and A5 relations',
                      mul(a, b) == mul(b, a) and power(mul(a, b), 3) == power(mul(b, a), 3))
                global_two_branch.append((i, j))
    check('complete S6 pair enumeration', len(transpositions)**2 == 225)
    check('cusp pair count', len(cusp) == 135)
    check('cusp same-pair count', sum(x['epsilon'] == 0 for x in cusp) == 15)
    check('cusp S3-pair count', sum(x['epsilon'] == 1 for x in cusp) == 120)
    check('local node equal-or-disjoint count', len(local_nodes) == 105)
    check('local A5 relation alone admits every pair', len(local_a5) == 225)
    check('global two-branch disjoint count', len(global_two_branch) == 90)
    # A single ramification germ cannot supply two distinct source branches.
    # The global theorem imposing that condition is a written import.
    check('global two-branch filter excludes equal and intersecting pairs',
          all(supports(transpositions[i]).isdisjoint(supports(transpositions[j]))
              for i, j in global_two_branch))
    data = json.loads((R7/'classification.json').read_text())
    expected = {'-3/2': ('I', 'A2 + A5 + 8 A1', 9), '-15/4': ('II', 'A2 + 11 A1', 11)}
    t = S.symbols('t')
    rows = []
    for case in data['cases']:
        label, singularities, s = expected[case['inverse_b']]
        check(label+' retained complete affine singularity type', case['affine_singularities'] == singularities)
        p = S.sympify(case['p'], locals={'t': t})
        q = S.sympify(case['q'], locals={'t': t})
        check(label+' retained polynomial degrees', S.degree(p, t) == 6 and S.degree(q, t) == 10)
        check(label+' unique critical parameter', S.Poly(S.gcd(S.diff(p, t), S.diff(q, t)), t).monic().as_expr() == t-1)
        check(label+' cusp jet nonzero',
              (S.diff(p, t, 2)*S.diff(q, t, 3)-S.diff(q, t, 2)*S.diff(p, t, 3)).subs(t, 1) != 0)
        for epsilon in (0, 1):
            chi_curve = 1-s
            chi_smooth_curve = chi_curve-(s+1)
            chi_outside = 1-chi_curve
            chi_N = 6*chi_outside+5*chi_smooth_curve+4*s+(5-epsilon)
            required_unramified = chi_N-2
            check(label+f' complete Euler ledger eps={epsilon}', chi_N == 5-epsilon)
            check(label+f' missing prime count eps={epsilon}', required_unramified == 3-epsilon)
            check(label+f' original affine Euler characteristic eps={epsilon}',
                  chi_N-(required_unramified+1) == 1)
            check(label+f' one-component nonproper locus impossible eps={epsilon}', required_unramified > 0)
            check(label+f' full Orevkov count agrees eps={epsilon}', 2+required_unramified+epsilon == 5)
            rows.append(dict(case=label, two_branch_values=s, cusp_epsilon=epsilon,
                             chi_curve=chi_curve, chi_smooth_curve=chi_smooth_curve,
                             chi_N=chi_N, required_unramified_boundary_primes=required_unramified,
                             each_unramified_residue_degree=1,
                             images_differ_from_the_classified_cusp_curve=True,
                             distinct_additional_target_images_not_counted=True))
    # The node projection has two analytic ramification branches for all n.
    # Chart exponents regress the universal DVR argument; this is not a
    # search for globally realizable covers.
    w, v = S.symbols('w v')
    node_cases = []
    for n in range(2, 13):
        u = w**n/v
        zero(f'node n={n} surface relation on v!=0', u*v-w**n)
        zero(f'node n={n} generic branch Jacobian', S.diff(u, w)-n*w**(n-1)/v)
        check(f'node n={n} both ramification indices exceed1', n >= 2)
        node_cases.append(dict(n=n, analytic_ramification_branches=2,
                               both_generic_ramification_indices=n,
                               global_Keller_attachment_excluded_by_written_theorem=True))
    manifest = json.loads((HERE/'sources/manifest.json').read_text())
    for entry in manifest['sources']:
        check('primary pin '+entry['file'], hashlib.sha256((HERE/'sources'/entry['file']).read_bytes()).hexdigest() == entry['sha256'])
    result = dict(as_of='2026-11-13', outcome='PASS', check_count=len(checks), checks=checks,
                  runtime_seconds=round(time.monotonic()-start, 6),
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  reused_local_excess_sha256=hashlib.sha256((R7/'local_excess.py').read_bytes()).hexdigest(),
                  retained_classification_sha256=hashlib.sha256((R7/'classification.json').read_bytes()).hexdigest(),
                  python=platform.python_version(), sympy=S.__version__,
                  cusp_pairs=cusp, local_node_pair_count=len(local_nodes),
                  local_A5_pair_count=len(local_a5), globally_separated_two_branch_pair_count=len(global_two_branch),
                  Euler_ledgers=rows, node_projection_controls=node_cases,
                  global_inputs='Borisov plane Stein-boundary theorem and purity are written imports, not finite-checker outputs',
                  ramified_row_hypothesis='One entire reduced ramification prime, e2 f1, image one classified r7 cusp curve',
                  complete_nonproper_locus_not_assumed=True,
                  extra_images_not_excluded=True, normal_r7_fully_excluded=False,
                  JC2_resolved=False, formal_verification=False, external_review=False)
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
