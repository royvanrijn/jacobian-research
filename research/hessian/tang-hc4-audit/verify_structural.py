"""Exact algebra supporting the written structural review and applications.

This checker does not certify the imported line-family classification or
the written case exhaustion. Run with a 60-second external worker cap.
"""
from pathlib import Path
import hashlib
import itertools
import json
import platform
import sys
import time
import flint
import sympy as S

HERE = Path(__file__).resolve().parent
checks = []


def check(name, value):
    assert value, name
    checks.append(name)


def zero(name, expression):
    check(name, S.cancel(expression) == 0)


def determinant(matrix, ctx):
    total = ctx.constant(0)
    for permutation in itertools.permutations(range(len(matrix))):
        term = ctx.constant((-1) ** sum(
            permutation[i] > permutation[j]
            for i in range(len(matrix)) for j in range(i + 1, len(matrix))))
        for i, j in enumerate(permutation):
            term *= matrix[i][j]
        total += term
    return total


def rank_one_update():
    pairs = list(itertools.combinations_with_replacement(range(4), 2))
    names = [f's{i}{j}' for i, j in pairs] + ['v0', 'v1', 'v2', 'v3', 't']
    ctx = flint.fmpq_mpoly_ctx.get(names, 'lex')
    g = ctx.gens()
    coefficients = dict(zip(pairs, g[:10]))
    matrix = [[coefficients[min(i, j), max(i, j)] for j in range(4)] for i in range(4)]
    v, t = g[10:14], g[14]
    updated = [[matrix[i][j] + t*v[i]*v[j] for j in range(4)] for i in range(4)]
    contraction = ctx.constant(0)
    for i in range(4):
        for j in range(4):
            minor = [[matrix[k][ell] for ell in range(4) if ell != i]
                     for k in range(4) if k != j]
            contraction += (-1)**(i+j)*v[i]*v[j]*determinant(minor, ctx)
    check('universal symmetric 4x4 rank-one determinant update',
          determinant(updated, ctx) == determinant(matrix, ctx) + t*contraction)


def line_frames():
    lam, a, d, alpha, gamma, rho, sigma = S.symbols('lam a d alpha gamma rho sigma')
    frame = S.Matrix([[a/(1-lam), gamma/lam, 0],
                      [alpha/(1-lam), d/lam, 0],
                      [rho/(1-lam), sigma/lam, 0]])
    zero('finite-endpoint trace numerator', lam*(1-lam)*S.trace(frame) - (a*lam+d*(1-lam)))
    c2 = sum(frame.extract([i, j], [i, j]).det() for i, j in itertools.combinations(range(3), 2))
    zero('finite-endpoint second characteristic coefficient',
         c2.subs({a:0, d:0}) + alpha*gamma/(lam*(1-lam)))
    s, t = S.symbols('s t')
    b, q = S.Function('b')(s, t), S.Function('q')(t)
    aa = -S.diff(b, t)
    cs = q*S.diff(aa, s) + S.diff(q, t)*S.diff(b, s)
    ct = q*S.diff(aa, t)
    zero('infinite-endpoint compatibility leaves q_second*b_s',
         S.diff(cs, t)-S.diff(ct, s)-S.diff(q,t,2)*S.diff(b,s))
    covector = S.Matrix([t*S.diff(q,t)-q, -S.diff(q,t), 1])
    ws = S.Matrix([S.diff(aa,s), t*S.diff(aa,s)+S.diff(b,s), cs])
    wt = S.Matrix([S.diff(aa,t), aa+t*S.diff(aa,t)+S.diff(b,t), ct])
    zero('infinite-endpoint covector annihilates W_s', (covector.T*ws)[0])
    zero('infinite-endpoint covector annihilates W_t', (covector.T*wt)[0])


def weights():
    m, p, q, beta, ell, D, eta, M = S.symbols('m p q beta ell D eta M')
    E = ell+m*beta
    zero('one-phase determinant weight', 4*E-2*(1+beta+E-p+E-q)-2*(p+q-1-beta))
    zero('three-phase slope', 4*m-2*3-(4*m-6))
    zero('two-phase plus partner slope', 4*m-2*(2+m)-(2*m-4))
    zero('binary-common fibre-raise slope', 4-2*2)
    zero('three-affine-partner slope', 4-2*3+2)
    zero('negative-active path determinant weight',
         4*D-2*(-eta+1+1+eta+1+2*eta)-(4*D-6-4*eta))
    zero('negative-active first-wall lower bound', (4*D-6-4*eta).subs(eta,D-2)-2)
    zero('lower-weight pure-power determinant weight',
         4*M*beta-2*(1+beta+M*beta-p+M*beta-q)-2*(p+q-1-beta))
    zero('ordered first two crossings',
         (p+1)/M-p/(M-1)-(M-p-1)/(M*(M-1)))
    zero('rank-one initial endpoint', 4*D-2*(1+3*(2*D-1)/3))
    zero('rank-two second-wall lower bound',
         4*D-2*(2+beta+D-1)-2*(D-beta-1))


def meng_yang():
    x,y,p,q,r,L,M,N = S.symbols('x y p q r L M N')
    variables = (x,y,p,q,r)
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y**2*u*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x**2*y)*r
    psi = S.Poly(L*A**2+M*A+N*B, *variables)
    degree_by_r = {}
    for exponent in psi.monoms():
        degree_by_r[exponent[4]] = max(degree_by_r.get(exponent[4], -1), sum(exponent[:4]))
    check('complete Meng-Yang support base degrees by r power', degree_by_r == {0:14, 1:10, 2:6})
    d = S.Symbol('d')
    zero('middle support degree is arithmetic mean', d+10-(14+2*d+6)/2)
    check('general degree-seven polynomial pullback bound', psi.total_degree() == 14 and 101//psi.total_degree() == 7)
    check('degree-47 potential bound is 100', max(14,2*47+6) == 100)
    check('degree-48 endpoint potential bound', max(14,2*48+6) == 102)
    r2 = S.Poly(psi.as_expr(),r).coeff_monomial(r**2)
    zero('unique degree-48 highest part is a square',r2-L*x**6)
    check('degree-48 next weighted layer is strictly lower', degree_by_r[1]+48 < 102 and degree_by_r[0]<102)
    check('degree-49 first uncovered degree after proper-power refinement', max(14,2*49+6)-2 == 102 > 100)
    v2 = psi.as_expr().subs({L:1,M:13,N:2})
    expected = S.Matrix([0,0,-S.Rational(1,2),0,0])
    grad = S.Matrix([S.diff(v2,v) for v in variables])
    points = [(1,-S.Rational(3,2),0,0,0),(-1,S.Rational(3,2),0,0,0)]
    for i, point in enumerate(points):
        check(f'marked v2 gradient at point {i}', grad.subs(dict(zip(variables,point))) == expected)
    jets = S.symbols('Rx Ry Rp Rq')
    graph_differential = S.eye(5,4)
    graph_differential[4,:] = S.Matrix([jets])
    check('collision transfer independent of all graph first derivatives',
          graph_differential.T*expected == expected[:4,0])
    return {'base_degree_by_r_power': degree_by_r,
            'all_graph_degrees_potential_bound': 'max(14,2*d+6)',
            'plane_degree_import': 100, 'potential_degree_consequence': 101,
            'proper_power_top_potential_degree_consequence':102,
            'marked_collision_graph_degree_excluded': 48,
            'unmarked_graphs_excluded': False,
            'general_pullback_degree_excluded_with_collision_jet_match': 7}


def main():
    start = time.monotonic()
    receipt = json.loads((HERE/'source-receipt.json').read_text())
    check('Tang source version hash', hashlib.sha256((HERE/'source-v1.pdf').read_bytes()).hexdigest() == receipt['sha256'])
    for row in json.loads((HERE/'structural-sources/sources.json').read_text()):
        check('import source hash: '+row['file'], row['outcome']=='downloaded' and
              hashlib.sha256((HERE/'structural-sources'/row['file']).read_bytes()).hexdigest()==row['sha256'])
    rank_one_update()
    line_frames()
    weights()
    applications = meng_yang()
    print(json.dumps(dict(outcome='PASS',checks=checks,check_count=len(checks),
        applications=applications,elapsed_seconds=round(time.monotonic()-start,3),
        runtime={'python':platform.python_version(),'sympy':S.__version__,'flint':flint.__version__},
        geometric_classification='imported De Poi Theorem0.1; written hypothesis review',
        full_structural_coverage='written in STRUCTURAL_REVIEW.md, not certified by this finite checker',
        formal_verification=False,external_review=False),indent=2))


if __name__ == '__main__':
    main()
