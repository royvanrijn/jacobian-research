"""Exact finite-algebra and accounting controls; global imports stay written.

No network, elimination campaign or missing historical-artifact reconstruction.
"""
from collections import deque
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def zero(name, expression):
    check(name, S.cancel(expression) == 0)


def quotient_length(equations, variables):
    """Count standard monomials of a zero-dimensional rational quotient.

    A nonzero finite quotient must have a pure leading power in every
    variable. BFS then visits its finite staircase, with a hard guard
    against accidentally running on a positive-dimensional input.
    """
    G = S.groebner(equations, *variables, order='lex', domain=S.QQ)
    leads = [tuple(p.LM(order=G.order).exponents) for p in G.polys]
    if any(not any(e) for e in leads):
        return 0
    check('zero-dimensional leading ideal', G.is_zero_dimensional)
    initial = (0,) * len(variables)
    todo, seen = deque([initial]), {initial}
    standard = set()
    while todo:
        mon = todo.popleft()
        if any(all(x >= y for x, y in zip(mon, lm)) for lm in leads):
            continue
        standard.add(mon)
        check('bounded finite staircase', len(standard) <= 1000)
        for i in range(len(variables)):
            newer = tuple(e + (j == i) for j, e in enumerate(mon))
            if newer not in seen:
                seen.add(newer)
                todo.append(newer)
    return len(standard)


def marked_case(label, n, P, t, a):
    m = n-2
    u, v = S.symbols('u v')
    degree = S.degree(P, t)
    Q = S.Poly(S.expand(u**n * P.subs(t, 1/u)), u).as_expr()
    infinity_k = int(n-degree)
    check(label+' coefficient T=1', S.expand(P).coeff(t, 1) == 1)
    finite_length = quotient_length([P, a**m-S.diff(P, t)], (a, t))
    check(label+' finite chart algebra length', finite_length == degree*m)
    infinity_length = 0
    if infinity_k:
        # Localize at u=0 by replacing Q by u^k. The complementary factor
        # is a unit there; retaining Q' in the marking relation is exact.
        infinity_length = quotient_length([u**infinity_k, v**m+S.diff(Q, u)], (v, u))
    check(label+' infinity local algebra length', infinity_length == infinity_k*m)
    factors = S.factor_list(P, t)[1]
    finite_simple = sum(S.degree(f, t) for f, k in factors if k == 1)
    finite_multiple = sum(S.degree(f, t)*k for f, k in factors if k >= 2)
    zero_root = int(P.subs(t, 0) == 0)
    simple = int(finite_simple) + int(infinity_k == 1)
    multiple = int(finite_multiple) + (infinity_k if infinity_k >= 2 else 0)
    retained = m*simple-(m-1)*zero_root
    omitted = m*multiple+(m-1)*zero_root
    check(label+' complete scheme length', finite_length+infinity_length == n*m)
    check(label+' omission identity', retained+omitted == n*m)
    return dict(label=label, n=n, m=m, polynomial=str(P),
                finite_chart_length=finite_length, infinity_multiplicity=int(infinity_k),
                infinity_local_length=infinity_length, retained_length=retained,
                omitted_length=omitted, simple_projective_roots=simple)


def completeness_gaps(scope):
    """Necessary data gates; an empty result would still require a proof.

    In particular this function never promotes a numerical balance to a
    completeness certificate for the normal algebra or reconstruction open.
    """
    gaps = []
    for field in ('complete_omission_list', 'exceptional_lengths_complete'):
        if scope[field] is not True:
            gaps.append(field)
    for row in scope['rows']:
        for field in ('target_image', 'ramification_index', 'residue_degree',
                      'reconstruction_pole', 'exceptional_ambient_lengths'):
            if row[field] is None:
                gaps.append(row['prime']+': '+field)
    gaps.extend(scope.get('unknowns', []))
    return gaps


def main():
    started = time.monotonic()
    ledger = json.loads((HERE/'ledger.json').read_text())
    for scope in [ledger['controls']['fonollosa'], *ledger['applications'].values()]:
        for row in scope['rows']:
            check('independent row fields '+row['prime'], set(ledger['row_fields']) <= row.keys())
            if row['ramification_index'] == 1:
                check('e1 does not ramify '+row['prime'], row['ramifies'] is False)
    completion_audit = {}
    for name, scope in ledger['applications'].items():
        check(name+' unresolved omission list', scope['complete_omission_list'] is False)
        check(name+' unresolved exceptional algebras', scope['exceptional_lengths_complete'] is False)
        check(name+' actual unknowns retained', bool(scope['unknowns']))
        check(name+' unknown target images stay null', any(r['target_image'] is None for r in scope['rows']))
        completion_audit[name] = completeness_gaps(scope)
        check(name+' completeness rejected on missing evidence', bool(completion_audit[name]))
    check('full research goal remains open', ledger['full_research_goal_complete'] is False)

    t, a, u, v, b = S.symbols('t a u v b')
    cases = []
    for n in (4, 5, 7):
        cases.extend([
            marked_case('simple generic n'+str(n), n, t**n+t-1, t, a),
            marked_case('zero root n'+str(n), n, t**n+t, t, a),
            marked_case('double finite root n'+str(n), n, -S.expand((t-1)**2*(t**(n-2)+1))/2, t, a),
            marked_case('simple infinity n'+str(n), n, t**(n-1)+t-1, t, a),
            marked_case('double infinity n'+str(n), n, t**(n-2)+t-1, t, a),
            marked_case('maximal infinity plus zero n'+str(n), n, t, t, a),
            marked_case('triple finite plus infinity n'+str(n), n, S.expand((t-1)**3)/3, t, a),
            marked_case('double finite plus zero n'+str(n), n, t*(t-1)**2, t, a),
        ])
        m = n-2
        check(f'n{n} discriminant DVR degree sum', 2*m+(n-2)*m == n*m)
        check(f'n{n} zero-hyperplane DVR degree sum', (m-1)+1+(n-1)*m == n*m)
        Q0 = u**(n-1)-b*u**n
        r = -v**m-S.diff(Q0, u)
        q = u*v**m+u*S.diff(Q0, u)-Q0
        zero(f'n{n} infinity-chart Jacobian',
             S.det(S.Matrix([r, q]).jacobian([u, v]))-m*v**(2*m-1))
        zero(f'n{n} Kummer gluing derivative',
             S.diff(u**n*(t**n+t-b).subs(t, 1/u), u)
             -n*u**(n-1)*(t**n+t-b).subs(t, 1/u)
             +u**(n-2)*S.diff(t**n+t-b, t).subs(t, 1/u))
        # Ecrit is t!=0 in this finite chart; the zero root is always simple.
        check(f'n{n} zero markings reduced', S.gcd(a**m-1, m*a**(m-1)) == 1)
    for n in (4, 5, 7):
        m = n-2
        for k in range(2, n+1):
            # A k-fold finite root at t=1 with coefficient T=1, so the
            # local equation below is the exact whole root algebra.
            P = S.expand((t-1)**k) / (k*(-1)**(k-1))
            length = quotient_length([P, a**m-S.diff(P, t)], (a, t))
            check(f'n{n} k{k} single multiple-root length', length == m*k)
            check(f'n{n} k{k} only reduced marking zero', S.diff(P, t).subs(t, 1) == 0)

    for case in cases:
        if case['label'].startswith('maximal infinity plus zero'):
            n, m = case['n'], case['m']
            check(f'n{n} marked-root singleton zero', case['retained_length'] == 1)
            check(f'n{n} singleton omitted length includes unramified',
                  case['omitted_length'] == m*(n-1)+(m-1) == n*m-1)
    check('n4 zero-target length ledger', ledger['controls']['fonollosa']['singleton_zero_target']['n4_lengths'] == [1, 6, 1])
    x, y, z3, z4, T = S.symbols('x y z3 z4 T')
    p = 1+x*y
    Phi4 = ((1-y*T)**3*(p*T-x)+(p*T-x)**2
            *(4*y*(y*T+S.Rational(1, 2))+z3*(p*T+x/2)+z4*T**2))
    ordinary_path = {x: b, y: 0, z3: 0, z4: 0}
    escaping_path = {x: -b, y: 2/b, z3: 12/b**2, z4: -8/b**3}
    zero('literal polynomial retained path Phi4=T-b', Phi4.subs(ordinary_path)-(T-b))
    zero('literal polynomial escaping path Phi4=T-b', Phi4.subs(escaping_path)-(T-b))
    zero('escaping source path marking -1', (1/p).subs(escaping_path)+1)
    zero('escaping source path root b', (x/p).subs(escaping_path)-b)
    zero('escaping source path has pole in y', (b*y).subs(escaping_path)-2)

    # First-order unit tests for omitted and retained zero-root charts:
    # t=0, a=zeta has invertible (t,a) Jacobian even if deleted.
    zeta = S.symbols('zeta', nonzero=True)
    for m in (2, 3, 5):
        J = S.Matrix([[1, 0], [-2, m*zeta**(m-1)]])
        zero(f'm{m} unramified zero-branch determinant', J.det()-m*zeta**(m-1))
    # At a=1 the sum is a unit and removes the apparent reconstruction pole.
    for m in (2, 3, 5):
        zero(f'm{m} reconstruction cancellation',
             (a-1)*sum(a**i for i in range(m))-(a**m-1))
        check(f'm{m} retained branch unit', sum(a**i for i in range(m)).subs(a, 1) == m)

    # An immersed nodal polynomial normalization has TWO length-one points
    # above its node. This models the image of one etale omitted prime;
    # it does not claim a global Keller realization of the image.
    p, q = t**2-1, t*(t**2-1)
    check('immersed nodal normalization', S.gcd(S.diff(p, t), S.diff(q, t)) == 1)
    node_gcd = S.gcd(p, q)
    check('one prime may contribute two points', S.degree(node_gcd, t) == 2)
    check('node normalization fibre reduced', S.gcd(node_gcd, S.diff(node_gcd, t)) == 1)
    check('node local inverse parameter proves birationality', S.cancel(q/p) == t)

    r7 = ledger['applications']['normal_r7_nonimmersed']
    samples = []
    for row in r7['strata']:
        L = sum(row['ramified_omitted_lengths'])
        capacity = 6-L
        check(row['name']+' base length', row['source_length_before_unramified_omission'] == capacity)
        check(row['name']+' branch capacity', row['max_unramified_branches'] == capacity)
        check(row['name']+' singleton branch count', row['unramified_branches_for_singleton'] == capacity-1)
        for branch_sum in range(capacity+1):
            source = 6-L-branch_sum
            check(row['name']+f' fibre identity B={branch_sum}', source+L+branch_sum == 6 and source >= 0)
        check(row['name']+' overfull omitted fibre rejected', 6-L-(capacity+1) < 0)
        samples.append(dict(stratum=row['name'], ramified_omitted_length=L,
                            max_unramified_branch_sum=capacity, singleton_branch_sum=capacity-1))
    # Replay both prior Euler calculations without rerunning their classification.
    prior = json.loads((HERE.parent/'jc2-boundary-attachment/receipt.json').read_text())
    for row in prior['Euler_ledgers']:
        eps, s = row['cusp_epsilon'], row['two_branch_values']
        chi_N = 6*s-10*s+4*s+5-eps
        check('r7 prior Euler '+str((s, eps)), chi_N-(1+3-eps) == 1)
        check('r7 extra prime count '+str((s, eps)), row['required_unramified_boundary_primes'] == 3-eps)
        check('r7 cusp singleton numerical compatibility '+str((s, eps)), 6-(2+eps)-(3-eps) == 1)

    # Finite-cover local model: discriminant pullback has a ramified and
    # unramified factor. The ramification ideal cannot choose the open.
    w, pi = S.symbols('w pi')
    z = w**3+pi**7*w
    zero('normal-r7 local cubic factors', -4*pi**21-27*z**2+(3*w**2+pi**7)**2*(3*w**2+4*pi**7))
    # At pi=-3,w=54 the companion factor vanishes but Jacobian is nonzero.
    check('local cubic companion generically etale witness',
          (3*w**2+4*pi**7).subs({pi: -3, w: 54}) == 0
          and S.diff(z, w).subs({pi: -3, w: 54}) != 0)

    manifest = json.loads((HERE/'sources/manifest.json').read_text())
    for entry in manifest['sources']:
        path = HERE/'sources'/entry['file']
        check('primary pin '+entry['file'], path.is_file() and
              hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'])
    inputs = ['PROOF.md', 'ledger.json', 'sources/manifest.json']
    inputs += ['../jc2-boundary-attachment/PROOF.md', '../jc2-boundary-attachment/receipt.json',
               '../hc4-relative-row-orbit/PROOF.md']
    result = dict(as_of='2026-11-14', outcome='PASS', check_count=len(checks),
                  checks=checks, runtime_seconds=round(time.monotonic()-started, 6),
                  python=platform.python_version(), sympy=S.__version__,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs},
                  marked_root_algebra_cases=cases, r7_pointwise_ledgers=samples,
                  completion_audit=completion_audit,
                  global_theorems='Written imports; not outputs of finite algebra checks',
                  reconstruction_criterion='Written normality and inverse-identity proof',
                  normal_r7_global_omission_list_complete=False,
                  singleton_zero_omission_list_complete=False,
                  full_research_goal_complete=False, formal_verification=False, external_review=False)
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
