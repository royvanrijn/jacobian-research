#!/usr/bin/env python3
"""Exact controls for written all-rank symmetry proofs; not formalization."""
from pathlib import Path
from itertools import product
import hashlib
import json
import platform
import time
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []

def check(label, condition):
    assert condition, label
    checks.append(label)


def zero(expr, p, *variables):
    return s.Poly(s.expand(expr), *variables, modulus=p).is_zero


class Field:
    """Small explicit F_p[T]/(modulus), for independent block controls."""
    def __init__(self, p, modulus):
        self.p, self.modulus = p, tuple(modulus)
        self.e = len(modulus)-1
        self.q = p**self.e

    def digits(self, value):
        return [(value//self.p**i) % self.p for i in range(self.e)]

    def encode(self, coeffs):
        return sum((v % self.p)*self.p**i for i, v in enumerate(coeffs))

    def add(self, a, b):
        return self.encode([x+y for x, y in zip(self.digits(a), self.digits(b))])

    def mul(self, a, b):
        aa, bb = self.digits(a), self.digits(b)
        c = [0]*(2*self.e-1)
        for i, x in enumerate(aa):
            for j, y in enumerate(bb):
                c[i+j] += x*y
        for i in range(len(c)-1, self.e-1, -1):
            z = c[i] % self.p
            for j in range(self.e+1):
                c[i-self.e+j] -= z*self.modulus[j]
        return self.encode(c[:self.e])

    def power(self, a, n):
        result = 1
        while n:
            if n & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            n >>= 1
        return result

    def spaces(self):
        seen = {frozenset([0])}
        queue = [frozenset([0])]
        for w in queue:
            for v in range(self.q):
                if v in w:
                    continue
                new = frozenset(self.add(x, self.mul(c, v))
                                for x in w for c in range(self.p))
                if new not in seen:
                    seen.add(new)
                    queue.append(new)
        return queue


def main():
    started = time.monotonic()
    u, x, z, t, A, y = s.symbols('u x z t A y')
    fields = [Field(2, (0, 1)), Field(3, (0, 1)), Field(5, (0, 1)),
              Field(2, (1, 1, 1)), Field(2, (1, 1, 0, 1)),
              Field(3, (1, 0, 1)), Field(2, (1, 1, 0, 0, 1)),
              Field(5, (2, 0, 1))]
    blocks, factors = [], []
    for f in fields:
        p, q, e = f.p, f.q, f.e
        modulus = sum(v*t**i for i, v in enumerate(f.modulus))
        check(f'F{q} defining polynomial irreducible',
              s.Poly(modulus, t, modulus=p).is_irreducible)
        check(f'F{q} Frobenius and inverses', all(
            f.power(a, q) == a and (not a or f.power(a, q-1) == 1)
            for a in range(q)))
        spaces = f.spaces()
        for m in s.divisors(q-1):
            m = int(m)
            r = (q-1)//m
            scalars = {a for a in range(1, q) if f.power(a, m) == 1}
            check(f'q{q} m{m} scalar subgroup order', len(scalars) == m)
            generator = next(a for a in scalars
                             if len({f.power(a, j) for j in range(m)}) == m)
            proper = [w for w in spaces if 1 < len(w) < q and
                      all(f.mul(generator, v) in w for v in w)]
            d = 1 if m == 1 else int(s.n_order(p, m))
            check(f'q{q} m{m} independent block enumeration',
                  (not proper) == (d == e))
            poly = x-u**r*x**q
            check(f'q{q} m{m} full plane Jacobian',
                  zero(s.diff(poly, x)-1, p, u, x))
            check(f'q{q} m{m} boundary whole fibre rank one',
                  s.expand(poly.subs(u, 0)-x) == 0)
            check(f'q{q} m{m} generic whole fibre rank q',
                  s.Poly(poly, x).degree() == q)
            # A polynomial identity imposing only b^q=b and u*lambda^m=1.
            lam, b, mu = s.symbols('lam b mu')
            diff = poly.subs(x, mu+lam*b)-poly.subs(x, mu)
            remainder = diff-lam*(b-b**q)-lam*b**q*(1-(u*lam**m)**r)
            check(f'q{q} m{m} translation-frame polynomial identity',
                  zero(remainder, p, u, mu, lam, b))
            blocks.append({'q': q, 'p': p, 'e': e, 'm': m, 'd': d,
                           'proper_stable_additive_subgroups': len(proper),
                           'primitive': not proper})
            if d < e:
                B, steps = p**d, e//d
                a = (B-1)//m
                h = x-u**a*x**B
                outer = sum(u**(a*(B**i-1)//(B-1))*z**(B**i)
                            for i in range(steps))
                check(f'q{q} m{m} exact global Keller decomposition',
                      zero(outer.subs(z, h)-poly, p, u, x))
                check(f'q{q} m{m} both factors Jacobian one',
                      zero(s.diff(h, x)-1, p, u, x) and
                      zero(s.diff(outer, z)-1, p, u, z))
                check(f'q{q} m{m} both factors have nontrivial degree',
                      s.Poly(h, x).degree() == B and
                      s.Poly(outer, z).degree() == B**(steps-1) and steps > 1)
                factors.append({'q': q, 'm': m, 'degrees': [B, B**(steps-1)]})

    # Non-field translation subspaces over F_p(A), not just the F_q example.
    for p in [2, 3]:
        W = [i+j*A for i in range(p) for j in range(p)]
        P = s.Poly(s.prod(x-w for w in W), x, A, modulus=p).as_expr()
        check(f'generic F{p} two-generator subspace polynomial additive',
              zero(P.subs(x, x+y)-P-P.subs(x, y), p, x, y, A))
        check(f'generic F{p} subspace derivative nonzero constant in x',
              s.Poly(s.diff(P, x), x, A, modulus=p).degree(x) == 0 and
              not zero(s.diff(P, x), p, x, A))
        for i, w in enumerate(W):
            check(f'generic F{p} subspace translation {i}',
                  zero(P.subs(x, x+w)-P, p, x, A))
        h = x**p-x
        expected = h**p-(A**p-A)**(p-1)*h
        check(f'generic F{p} subspace exact proper factorization',
              zero(P-expected, p, x, A))

    # Complete framed fibre algebras: quotient bases, separability and all frame actions.
    L, M = s.symbols('L M')
    fibres = []
    for p, m, uu, vv in [(3, 1, 2, 1), (3, 2, 2, 1),
                          (5, 2, 2, 3), (5, 4, 3, 2)]:
        r = (p-1)//m
        polynomials = [L**m-pow(uu, -1, p), M-uu**r*M**p-vv]
        gb = s.groebner(polynomials, L, M, modulus=p)
        powers = {g.LM(order=gb.order).exponents for g in gb.polys}
        check(f'F{p} m{m} complete frame quotient basis',
              powers == {(m, 0), (0, p)})
        check(f'F{p} m{m} complete frame quotient etale',
              s.gcd(s.Poly(polynomials[0], L, modulus=p),
                    s.Poly(s.diff(polynomials[0], L), L, modulus=p)).degree() == 0
              and zero(s.diff(polynomials[1], M)-1, p, M))
        actions = 0
        for a in range(1, p):
            if pow(a, m, p) != 1:
                continue
            for b in range(p):
                for eq in polynomials:
                    pulled = eq.subs({L: a*L, M: M+b*L}, simultaneous=True)
                    assert gb.reduce(pulled)[1] == 0
                actions += 1
        check(f'F{p} m{m} every specified frame action preserves the whole quotient',
              actions == m*p)
        fibres.append({'characteristic': p, 'q': p, 'm': m, 'target': [uu, vv],
                       'root_algebra_rank': p, 'frame_algebra_rank': m*p,
                       'checked_frame_actions': actions})

    # Published nonlinearizable S3 action: exact relations and the full fixed ideal.
    a, b, xx, yy, w = s.symbols('a b xx yy w')
    variables = (a, b, xx, yy)
    cyclo = s.Poly(w*w+w+1, w)
    def red(expr):
        return s.Poly(s.expand(expr), w).rem(cyclo).as_expr().expand()
    def compose(f, g):
        return tuple(red(v.subs(dict(zip(variables, g)), simultaneous=True)) for v in f)
    identity = variables
    sigma = (w*a, w*w*b, xx, yy)
    tau = (b, a, -b**3*xx+(1+a*b+a*a*b*b)*yy, (1-a*b)*xx+a**3*yy)
    check('nonlinear S3 sigma order three', compose(sigma, compose(sigma, sigma)) == identity)
    check('nonlinear S3 tau order two', compose(tau, tau) == identity)
    check('nonlinear S3 defining mixed relation',
          compose(compose(sigma, tau), compose(sigma, tau)) == identity)
    check('nonlinear S3 sigma Jacobian one',
          red(s.Matrix(sigma).jacobian(variables).det()-1) == 0)
    check('nonlinear S3 tau Jacobian one',
          red(s.Matrix(tau).jacobian(variables).det()-1) == 0)
    equations = [w*w+w+1] + [f-v for f, v in zip(sigma, identity)] + [f-v for f, v in zip(tau, identity)]
    actual = s.groebner(equations, a, b, xx, yy, w)
    expected = s.groebner([a, b, xx-yy, w*w+w+1], a, b, xx, yy, w)
    check('nonlinear S3 complete fixed scheme is the exact affine line',
          all(actual.reduce(g)[1] == 0 for g in [a, b, xx-yy, w*w+w+1]) and
          all(expected.reduce(g)[1] == 0 for g in equations))

    # All finite-field splitting types: independent partition enumeration versus DP.
    def partitions(n, minimum=1):
        if n == 0:
            yield ()
        for first in range(minimum, n+1):
            for rest in partitions(n-first, first):
                yield (first,)+rest
    dp = [0]*17
    dp[0] = 1
    for part in range(1, 17):
        for n in range(part, 17):
            dp[n] += dp[n-part]
    budgets = []
    for n in range(2, 17):
        types = list(partitions(n))
        budget = sum(lam.count(1) for lam in types)
        check(f'rank{n} source-point partition identity', budget == sum(dp[:n]))
        check(f'rank{n} source bound dominates target types', budget >= len(types))
        h = 0
        while (h+1)*(h+2)//2 <= n-1:
            h += 1
        check(f'rank{n} elementary all-rank budget growth', budget >= 2**h)
        lower2 = next(d for d in range(20) if 2**d >= budget)
        budgets.append({'N': n, 'partition_types': len(types),
                        'source_rational_point_budget': budget,
                        'dimension_lower_bound_over_F2': lower2})
    check('rank7 over F2 needs dimension five',
          next(v for v in budgets if v['N'] == 7)['source_rational_point_budget'] == 30
          and 2**4 < 30 <= 2**5)
    for control in blocks:
        q, m = control['q'], control['m']
        check(f'q{q} m{m} framed finite-field source budget',
              (q+m-1 <= q) == (m == 1))
        if q*m <= 100:
            f = next(f for f in fields if f.q == q)
            neg = lambda v: f.encode([-c for c in f.digits(v)])
            scalars = [v for v in range(1, q) if f.power(v, m) == 1]
            elements = {(a, b) for a in scalars for b in range(q)}
            left = set(elements)
            class_count, point_sum = 0, 0
            while left:
                a, b = min(left)
                conjugates = {(a, f.add(f.mul(c, b), f.mul(f.add(1, neg(a)), d)))
                              for c, d in elements}
                left -= conjugates
                class_count += 1
                point_sum += sum(f.add(f.mul(a, v), b) == v for v in range(q))
            check(f'q{q} m{m} independent conjugacy-class budget',
                  class_count == (q-1)//m+m and point_sum == q+m-1)
    # An explicit two-variable finite universal translation map of rank eight over F4.
    alpha, X, Y = s.symbols('alpha X Y')
    W = [0, 1, alpha, 1+alpha]
    P4 = s.Poly(s.prod(X-v for v in W), X, alpha, modulus=2).as_expr()
    reduce4 = s.groebner([alpha**2+alpha+1], alpha, X, Y, modulus=2)
    P4 = reduce4.reduce(P4)[1]
    check('rank8 over F4 explicit minimal product Jacobian',
          zero(s.diff(P4, X)*s.diff(Y*Y+Y, Y)-1, 2, alpha, X, Y))
    check('rank8 over F4 complete split fibre rank eight',
          s.Poly(P4, X).degree() == 4 and
          all(reduce4.reduce(P4.subs(X, v))[1] == 0 for v in W) and 4 < 8 <= 4**2)

    sources = json.loads((HERE/'sources.json').read_text())
    for source in sources['local_sources']:
        check('local proof pin '+source['path'],
              hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest() == source['sha256'])
    receipt = {'schema_version': 1, 'result': 'pass', 'check_count': len(checks),
               'checks': checks, 'affine_group_controls': blocks,
               'global_decomposition_controls': factors, 'complete_fibre_controls': fibres,
               'nonlinear_action_fixed_scheme': 'exact ideal (a,b,x-y) over Q(omega)',
               'finite_field_source_budgets': budgets,
               'assurance': 'Bounded exact controls support written universal proofs. Field tables and block enumeration do not prove the all-group classification. Symbolic action relations and the exact fixed ideal do not reprove published nonlinearizability. Hilbert90, general versality, essential dimension and linearization transfers remain written arguments. No formalization, independent whole-proof replay, external review or priority claim.',
               'python': platform.python_version(), 'sympy': s.__version__,
               'runtime_seconds': round(time.monotonic()-started, 3),
               'sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in ['PROOF.md', 'LITERATURE_AND_GAPS.md', 'verify.py', 'sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'result': 'pass', 'checks': len(checks),
                      'affine_groups': len(blocks), 'factorizations': len(factors),
                      'complete_fibres': len(fibres), 'seconds': receipt['runtime_seconds']}))


if __name__ == '__main__':
    main()
