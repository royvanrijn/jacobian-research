"""Bounded exact audit and constructive follow-up to arXiv:2610.24889v1."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parent
X = s.symbols('x1:19')


def homogeneous(f, degree, variables=X):
    return s.Add(*(co * s.prod(v**e for v, e in zip(variables, mon))
                   for mon, co in s.Poly(s.expand(f), *variables).terms()
                   if sum(mon) == degree))


def cubic12(scale11=1):
    a, b, c = X[:3]
    v = 1 + 2*a*b
    base = [a-3*a*a*b-a**3*c,
            3*v*v*a*c+(9*v*v-6*v-2)*b,
            v**3*c+(6*v*v+2*v)*b*b]
    f = base + [X[3]+a*a*c+3*a*b, X[4]+b**3, X[5]+a*a,
                X[6]+a*c, X[7]+b*b, X[8]+a*X[4],
                X[9]+a*X[7], X[10]+scale11*a*X[8],
                X[11]+a*b+a*X[9]]
    k = [f[0]+f[5]*f[6],
         f[1]+12*f[0]*f[7]+12*f[3]*f[9]-12*f[6]*f[11],
         f[2]+8*f[0]*f[4]-12*f[3]*f[7]+8*f[3]*f[8]
         +12*f[6]*f[9]-s.Rational(8,scale11)*f[6]*f[10]] + f[3:]
    p = [0,1,-8,0,-1,0,0,-1,0,0,0,0]
    q = [-1,1,4,-1,-1,-1,4,-1,-1,-1,-scale11,0]
    return base, list(map(s.expand, k)), [p,q]


def homogeneous18(k):
    # This sign convention is the author's: tail = w+C, top = x+tQ-t^2w.
    t = X[17]
    q = [homogeneous(f,2) for f in k]
    c = [homogeneous(f,3) for f in k]
    assert all(f == 0 for f in c[5:])
    h = [t*f-(t*t*X[12+i] if i < 5 else 0) for i,f in enumerate(q)]
    return list(map(s.expand,h+c[:5]+[s.Integer(0)])), c[:5]


def coefficient_rank(polynomials, variables=X):
    polys = [s.Poly(f,*variables) for f in polynomials]
    mons = sorted(set(m for p in polys for m,_ in p.terms()))
    return s.Matrix([[p.coeff_monomial(m) for m in mons] for p in polys]).rank()


def evaluate(polynomials, point, variables=X):
    return [s.expand(f.subs(dict(zip(variables,point)), simultaneous=True)) for f in polynomials]


def grouped18(h, scale11=1, mask=511):
    y = s.symbols('y1:6')
    z = [2*X[8]-3*X[7],3*X[9]-s.Rational(2,scale11)*X[10],
         3*X[1]-X[6],X[1]+X[9],3*X[1]+6*X[9]-s.Rational(4,scale11)*X[10]]
    g = list(h)
    replacements = [
        (2,4*y[0]*X[3]*X[17]), (2,4*y[1]*X[6]*X[17]),
        (11,y[3]*X[0]*X[17]), (12,-y[2]*X[0]**2),
        (13,12*y[2]*y[3]*X[0]), (14,2*y[4]*X[0]*X[2]),
        (14,4*y[0]*y[2]*X[0]), (1,12*y[3]*X[3]*X[17]),
        (10,scale11*s.Rational(1,2)*y[0]*X[0]*X[17])]
    for bit,(i,term) in enumerate(replacements):
        if mask & (1<<bit):
            g[i] = term+s.expand(g[i]-term.subs(dict(zip(y,z)), simultaneous=True))
    assert all(s.expand(a-b.subs(dict(zip(y,z)),simultaneous=True)) == 0 for a,b in zip(h,g))
    return g,y,z


def main():
    base,k,points = cubic12()
    assert s.factor(s.Matrix(base).jacobian(X[:3]).det()) == 1
    assert evaluate(k,points[0]) == evaluate(k,points[1]) == [0,1]+[0]*10
    degrees = [s.Poly(f,*X).total_degree() for f in k]
    assert degrees == [3]*5+[2]*7
    h,c = homogeneous18(k)
    lifted = [p+[-v for v in evaluate(c,p)]+[1] for p in points]
    gh = [a+b for a,b in zip(X,h)]
    assert evaluate(gh,lifted[0]) == evaluate(gh,lifted[1])
    assert all(homogeneous(f,3)==f for f in h)
    result = {
        'source_sha256':hashlib.sha256((ROOT/'source-v1.tar').read_bytes()).hexdigest(),
        'base_jacobian':1,'component_degrees':degrees,
        'cubic_output_rank':coefficient_rank(c),
        'homogeneous_dimension':18,'homogeneous_output_rank':coefficient_rank(h),
        'collision12':points, 'collision18':lifted,
        'collision_image18':evaluate(gh,lifted[0]),
    }
    (ROOT/'audit.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps(result,indent=2,default=str))


if __name__ == '__main__':
    main()
