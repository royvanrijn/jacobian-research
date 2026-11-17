#!/usr/bin/env python3
"""Frozen exact controls for the written torus/Whitehead corollary.

Requires Sympy 1.14.0 for symbolic identities. Finite field arithmetic and
matrix-algebra spans below use ordinary integer operations. These checks
are not a formal proof or independent verification of the source theorem.
"""
import hashlib
import json
from math import gcd
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent


def ident(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def mm(a, b, r):
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) % r
                       for col in zip(*b)) for row in a)


def scale(a, c, r):
    return tuple(tuple(x*c % r for x in row) for row in a)


def inverse(a, r):
    n = len(a)
    rows = [list(x)+list(y) for x, y in zip(a, ident(n))]
    for j in range(n):
        k = next(k for k in range(j, n) if rows[k][j])
        rows[j], rows[k] = rows[k], rows[j]
        rows[j] = [x*pow(rows[j][j], -1, r) % r for x in rows[j]]
        for k in range(n):
            if k != j:
                c = rows[k][j]
                rows[k] = [(x-c*y) % r for x, y in zip(rows[k], rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def mpow(a, exponent, r):
    if exponent < 0:
        return mpow(inverse(a, r), -exponent, r)
    out = ident(len(a))
    while exponent:
        if exponent & 1:
            out = mm(out, a, r)
        a = mm(a, a, r)
        exponent //= 2
    return out


def add_basis(basis, vector, r):
    row = list(vector)
    for j, b in sorted(basis.items()):
        c = row[j]
        row = [(x-c*y) % r for x, y in zip(row, b)]
    j = next((j for j, x in enumerate(row) if x), None)
    if j is None:
        return False
    c = pow(row[j], -1, r)
    basis[j] = tuple(x*c % r for x in row)
    return True


def algebra_dimension(gens, r):
    n = len(gens[0])
    queue, basis = [ident(n)], {}
    add_basis(basis, sum(queue[0], ()), r)
    for a in queue:
        for g in gens:
            b = mm(a, g, r)
            if add_basis(basis, sum(b, ()), r):
                queue.append(b)
    assert len(queue) <= n*n
    return len(basis)


def companion(coefficients, r):
    n = len(coefficients)-1
    return tuple(tuple((int(i == j+1) if j < n-1 else -coefficients[i]) % r
                       for j in range(n)) for i in range(n))


def poly_roots(roots, r):
    coeff = [1]
    for root in roots:
        new = [0]*(len(coeff)+1)
        for j, c in enumerate(coeff):
            new[j] = (new[j]-root*c) % r
            new[j+1] = (new[j+1]+c) % r
        coeff = new
    return coeff


def torus_controls(plan):
    results = []
    for p, q in plan['torus_pairs']:
        assert 2 <= p < q and gcd(p, q) == 1
        r = next(r for r in range(4*p*q+2, 5001)
                 if (r-1) % q == 0 and sp.isprime(r))
        g = int(sp.primitive_root(r))
        zeta = pow(g, (r-1)//q, r)
        delta = [pow(zeta, j, r) for j in range(p)]
        c = (-1)**(p+1)
        for x in delta:
            c = c*x % r
        a = pow(q, -1, p)
        b = (1-a*q)//p
        u, v = pow(-a, -1, p), pow(b, -1, q)
        f = companion([-1]+[0]*(p-1)+[1], r)
        samples = []
        for s in range(1, r):
            mu = c*pow(s, p, r) % r
            if mu == 1 or any(pow(s*x, p, r) == 1 for x in delta):
                continue
            cc = companion(poly_roots([s*x % r for x in delta], r), r)
            dd = scale(cc, pow(s, -1, r), r)
            x = scale(mpow(f, u, r), pow(s, q, r), r)
            y = scale(mpow(dd, v, r), pow(s, p, r), r)
            center = scale(ident(p), pow(s, p*q, r), r)
            assert mpow(x, p, r) == mpow(y, q, r) == center
            meridian = mm(mpow(x, a, r), mpow(y, b, r), r)
            assert meridian == mm(inverse(f, r), cc, r)
            basis = {}
            for i, row in enumerate(meridian):
                add_basis(basis, [(z-int(i == j)) % r for j, z in enumerate(row)], r)
            assert len(basis) == 1
            longitude = mm(center, mpow(meridian, -p*q, r), r)
            lam = pow(s, p*q, r)*pow(mu, -p*q, r) % r
            diff = tuple(tuple((z-int(i == j)) % r for j, z in enumerate(row))
                         for i, row in enumerate(meridian))
            assert mm(longitude, diff, r) == scale(diff, lam, r)
            assert lam*pow(mu, (p-1)*q, r) % r == (-1)**(p+1) % r
            dim = algebra_dimension([f, cc], r)
            assert dim == p*p
            samples.append({'s':s, 'mu':mu, 'lambda':lam, 'algebra_dimension':dim})
            if len(samples) == plan['limits']['samples_per_torus_pair']:
                break
        assert len(samples) == 3
        results.append({'p':p, 'q':q, 'prime':r, 'samples':samples})
    return results


def symbolic_torus(plan):
    s, tau = sp.symbols('s tau')
    relation = tau*tau-tau-1
    def red(e):
        return sp.cancel(sp.rem(sp.expand(e), relation, tau))
    def matred(a):
        return a.applyfunc(red)
    def power(a, n):
        if n < 0:
            a, n = matred(a.inv()), -n
        out = sp.eye(3)
        for _ in range(n):
            out = matred(out*a)
        return out
    f = sp.Matrix([[0,0,1],[1,0,0],[0,1,0]])
    cc = sp.Matrix([[0,0,s**3],[1,0,(tau-1)*s**2],[0,1,(1-tau)*s]])
    dd = cc/s
    assert power(dd, 5) == sp.eye(3)
    x, y = s**5*f, s**3*power(dd, 3)
    assert power(x, 3) == power(y, 5) == s**15*sp.eye(3)
    m = matred(power(x, -1)*power(y, 2))
    expected = sp.Matrix([[1,0,(tau-1)*s**2],[0,1,(1-tau)*s],[0,0,s**3]])
    assert matred(m-expected) == sp.zeros(3)
    longitude = s**15*power(m, -15)
    assert matred(longitude*(m-sp.eye(3))-s**-30*(m-sp.eye(3))) == sp.zeros(3)
    assert sp.diff(s**3, s).subs(s, -1) == 3
    assert sp.diff(s**-30, s).subs(s, -1) == 30
    from quadratic_integer_lifts import IDENTITY, matmul
    def specialize(a):
        out = []
        for e in a:
            z = sp.Poly(red(e.subs(s, -1)), tau)
            out.append((int(z.nth(0)), int(z.nth(1))))
        return tuple(out)
    xx, yy, mm0 = map(specialize, [x,y,m])
    assert matmul(mm0, mm0) == IDENTITY
    found, queue = {IDENTITY}, [IDENTITY]
    for a in queue:
        for b in [xx, yy]:
            z = matmul(a,b)
            if z not in found:
                found.add(z)
                queue.append(z)
                assert len(found) <= plan['limits']['H3_group_order_cap']
    assert len(found) == 120
    return {'coefficient_ring':'Q(tau)(s), tau^2=tau+1', 'identities':'PASS',
            'mu':'s^3', 'lambda':'s^-30', 'eta':'s^-60',
            'specialization_s':-1, 'specialization_group_order':len(found),
            'd_mu_at_minus1':3, 'd_lambda_at_minus1':30}


def symbolic_whitehead():
    m, l, d = sp.symbols('mu lambda d')
    eta, a = l*l, m*d
    r = eta*(1-m)*a*d+(m*eta**2-1)*(a-1)*(d-1)
    x = sp.Matrix([[a,d-1],[a-1,d]])
    ll = sp.diag(l,1/l)
    y = ll.inv()*x*ll
    def zero_mod(e):
        num = sp.fraction(sp.cancel(e))[0]
        return sp.cancel(sp.rem(num, r, d)) == 0
    comm = y*x.inv()*y.inv()*x-sp.diag(m,1/m)
    assert all(zero_mod(e) for e in comm)
    long = x.inv()*ll*x*ll.inv()*x.inv()*ll.inv()*x*ll
    lam = (eta*(a-1)+d-1)/(a-1+eta*(d-1))
    assert all(zero_mod(e) for e in long*sp.ones(2,1)-lam*sp.ones(2,1))
    assert zero_mod(x.det()-(a+d-1))
    assert (x-sp.eye(2)).det() == 0
    return {'boundary_commutator':True, 'longitude_word_eigenvalue':True,
            'rank_one_meridian':True, 'method':'rational polynomial remainder modulo R, with a=mu*d'}


def finite_whitehead(plan):
    results, tried = [], 0
    for r in [5,7,11,13]:
        count = 0
        for mu in range(2,r):
            for lam in range(1,r):
                eta = lam*lam % r
                for d in range(r):
                    tried += 1
                    a = mu*d % r
                    rr = eta*(1-mu)*a*d+(mu*eta*eta-1)*(a-1)*(d-1)
                    bad = (eta-1)*(mu*eta+1)*(mu*mu*eta*eta-mu*eta*eta+4*mu*eta-mu+1)
                    if rr % r or bad % r == 0:
                        continue
                    new_mu = (a+d-1) % r
                    assert new_mu not in [0,1] and (a-1)*(d-1) % r
                    x = ((a,(d-1)%r),((a-1)%r,d))
                    ll = ((lam,0),(0,pow(lam,-1,r)))
                    li = inverse(ll,r)
                    y = mm(mm(li,x,r),ll,r)
                    comm = mm(mm(mm(y,inverse(x,r),r),inverse(y,r),r),x,r)
                    assert comm == ((mu,0),(0,pow(mu,-1,r)))
                    word = [inverse(x,r),ll,x,li,inverse(x,r),li,x,ll]
                    longitude = ident(2)
                    for z in word:
                        longitude = mm(longitude,z,r)
                    den = a-1+eta*(d-1)
                    assert den % r
                    new_lam = (eta*(a-1)+d-1)*pow(den,-1,r) % r
                    assert all(sum(row) % r == new_lam for row in longitude)
                    assert new_lam != 0
                    count += 1
        results.append({'prime':r, 'valid_blocks':count})
    assert tried == plan['limits']['whitehead_triples']
    assert sum(row['valid_blocks'] for row in results) > 0
    return {'triples_examined':tried, 'valid':results}


def main():
    plan = json.loads((HERE/'torus_whitehead_plan.json').read_text())
    assert plan['frozen_before_execution'] and sp.__version__ == '1.14.0'
    result = {'status':'PASS', 'sympy':sp.__version__,
              'torus_finite_controls':torus_controls(plan),
              'symbolic_T35':symbolic_torus(plan),
              'symbolic_whitehead':symbolic_whitehead(),
              'finite_whitehead':finite_whitehead(plan),
              'assurance':'Exact controls for a written working proof; no independent full-source audit, formal verification, novelty or full MRC claim'}
    inputs = ['verify_torus_whitehead.py','torus_whitehead_plan.json','quadratic_integer_lifts.py',
              'TORUS_WHITEHEAD_PERIPHERAL_CURVES.tex','torus_whitehead_source_audit.json',
              'sources/cornwell-knot-groups-v5.pdf','sources/cornwell-kch-v3.pdf','sources/whitehead-v1.pdf']
    result['input_sha256'] = {p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    (HERE/'torus_whitehead_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'torus_finite_controls' and k != 'input_sha256'},indent=2))


if __name__ == '__main__':
    main()
