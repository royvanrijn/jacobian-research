"""Certified upper/lower parameter-height bounds for a retained pointed chart.

Sage produces exact real-algebraic extrema and finite residue trees. A separate
Fraction/SymPy checker verifies the packet. No point enumeration or minimisation
is performed here. Local bounds are safe upper bounds, not asserted optimal.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from sage.all import AA, GF, PolynomialRing, QQ, RealIntervalField, ZZ, gcd, lcm

R = PolynomialRing(QQ, 't')
t = R.gen()
BALL = RealIntervalField(192)


def coefficients(f, degree=4):
    return [f[i] for i in range(degree+1)]


def reverse(f):
    return R(coefficients(f)[::-1])


def substitute(values, matrix):
    a,b,c,d = map(QQ, matrix)
    if a*d == b*c:
        raise ArithmeticError('singular projective change')
    return sum(QQ(v)*(a*t+b)**i*(c*t+d)**(4-i) for i,v in enumerate(values))


def integral_pair(n, d):
    scale = lcm(v.denominator() for v in coefficients(n)+coefficients(d))
    content = gcd(ZZ(v*scale) for v in coefficients(n)+coefficients(d))
    return R(n*scale/content), R(d*scale/content)


def bezout(f, g):
    h,u,v = f.xgcd(g)
    if h != 1:
        raise ArithmeticError('covering map has a base point')
    k = lcm(c.denominator() for c in list(u)+list(v))
    return {'constant': str(k), 'left': list(map(str,u*k)),
            'right': list(map(str,v*k))}


def fixed_nonsquare(q, p):
    """Exclude a whole Z_p ball using its constant valuation and unit."""
    values = [ZZ(c) for c in q.list()]
    v = min(c.valuation(p) for c in values if c)
    digits = 3 if p == 2 else 1
    fixed = q[0] and ZZ(q[0]).valuation(p) == v and all(
        not c or c.valuation(p) >= v+digits for c in values[1:])
    if not fixed:
        return False
    unit = ZZ(q[0])//p**v
    return bool(v % 2 or (unit % 8 != 1 if p == 2 else pow(int(unit),int((p-1)//2),int(p)) != 1))


def local_tree(n, d, q, p, cap, node_limit):
    field = GF(p)
    polys = PolynomialRing(field, 'z')
    count = 0

    def walk(f, g, q, accumulated):
        nonlocal count
        count += 1
        if count > node_limit:
            # The original Bezout bound remains valid on an unfinished ball.
            return {'kind': 'budget-fallback', 'upper': int(cap)}
        if fixed_nonsquare(q,p):
            return {'kind': 'nonsquare', 'upper': -1}
        e = min(ZZ(c).valuation(p) for c in list(f)+list(g) if c)
        at = accumulated+e
        if at >= cap:
            return {'kind': 'bezout-cap', 'upper': int(cap)}
        f,g = f/p**e,g/p**e
        common = polys(f).gcd(polys(g))
        roots = sorted(map(int,common.roots(multiplicities=False))) if common.degree() > 0 else []
        children = [{'residue':r, 'node':walk(f(r+p*t),g(r+p*t),q(r+p*t),at)} for r in roots]
        return {'kind':'split', 'content_valuation':int(e), 'roots':roots,
                'children':children, 'upper':max([int(at)]+[c['node']['upper'] for c in children])}

    affine = walk(n,d,q,0)
    infinity = walk(reverse(n)(p*t),reverse(d)(p*t),reverse(q)(p*t),0)
    upper = max(affine['upper'],infinity['upper'])
    if upper < 0:
        raise ArithmeticError('pointed curve unexpectedly locally insoluble')
    return {'prime':str(p), 'bezout_cap':int(cap), 'upper':upper,
            'nodes':count, 'affine':affine, 'infinity':infinity}


def real_bounds(n,d,q):
    values = []
    candidate_count = 0
    for f,g,square in ((n,d,q),(reverse(n),reverse(d),reverse(q))):
        candidates = [AA(-1),AA(1)]
        # On a compact chart extrema occur at endpoints, branch points,
        # crossings of |f| and |g|, or a stationary active polynomial.
        for polynomial in (f.derivative(),g.derivative(),f-g,f+g,square):
            if polynomial:
                candidates += polynomial.roots(AA,multiplicities=False)
        for x in candidates:
            if -1 <= x <= 1 and square(x) >= 0:
                values.append(max(abs(f(x)),abs(g(x))))
                candidate_count += 1
    low,high = min(values),max(values)
    if low <= 0:
        raise ArithmeticError('nonpositive real map bound')
    # Round outward, with a small extra margin for efficient independent
    # rational root isolation. These are bounds, not floating estimates.
    margin = QQ(1)/2**32
    low = BALL(low).lower().exact_rational()*(1-margin)
    high = BALL(high).upper().exact_rational()*(1+margin)
    return {'lower':str(low), 'upper':str(high), 'candidate_count':candidate_count,
            'domain':'max(|U|,|V|)=1 and discriminant_quartic(U,V)>=0'}


def interval(value):
    return [str(value.lower().exact_rational()),str(value.upper().exact_rational())]


def build(payload, node_limit=4096):
    started = time.process_time()
    A,B = map(QQ,payload['curve'][3:])
    if list(map(QQ,payload['curve'][:3])) != [0,0,0] or 4*A**3+27*B**2 == 0:
        raise ArithmeticError('nonsingular short Weierstrass model required')
    mapping = payload['mapping']
    raw = list(map(QQ,mapping['raw_coefficients']))
    a,b = -raw[2]/6,-raw[1]/8
    if raw != [-3*a*a-4*A,-8*b,-6*a,0,1] or b*b != a**3+A*a+B:
        raise ArithmeticError('pointed input does not agree with the equation')
    nq = [a**3+4*B,4*a*b,6*a*a+4*A,4*b,a]
    n = substitute(nq,mapping['matrix'])
    d = substitute(raw,mapping['matrix'])
    q = R(list(map(QQ,mapping['discriminant_quartic'])))
    ratio = QQ(mapping['square_ratio'])
    if ratio <= 0 or not ratio.is_square() or d != ratio*q:
        raise ArithmeticError('exact discriminant transport failed')
    if q.gcd(q.derivative()).degree() > 0:
        raise ArithmeticError('singular quartic')
    n,d = integral_pair(n,d)
    packets = [bezout(n,d),bezout(reverse(n),reverse(d))]
    uniform = lcm(ZZ(p['constant']) for p in packets)
    real_start = time.process_time()
    real = real_bounds(n,d,q)
    real_cpu = time.process_time()-real_start
    remaining = ZZ(uniform)
    locals_ = []
    local_start = time.process_time()
    for p in map(ZZ,payload['certified_primes']):
        cap = remaining.valuation(p)
        if cap:
            remaining //= p**cap
            locals_.append(local_tree(n,d,q,p,cap,node_limit))
    divisor = remaining
    for row in locals_:
        divisor *= ZZ(row['prime'])**row['upper']
    lower = -BALL(QQ(real['upper'])).log()/4
    upper = (BALL(divisor).log()-BALL(QQ(real['lower'])).log())/4
    H = ZZ(payload['search_height'])
    radius = H**4*QQ(real['lower'])/divisor
    result = {'schema':'pointed-chart-height-bounds.v1',
        'status':'PRODUCED_PENDING_INDEPENDENT_REPLAY',
        'input':payload, 'anchor':list(map(str,(a,b))),
        'numerator':list(map(str,coefficients(n))),
        'denominator':list(map(str,coefficients(d))),
        'bezout':packets, 'uniform_bezout_divisor':str(uniform),
        'real':real, 'local_trees':locals_, 'unfactored_residual':str(remaining),
        'finite_gcd_divisor':str(divisor), 'B1_interval':interval(lower),
        'B2_interval':interval(upper),
        'guaranteed_multiplicative_x_height_radius':str(radius),
        'display':{'B1':float(lower.center()), 'B2':float(upper.center()),
            'log_guaranteed_x_height_radius':float(BALL(radius).log().center()),
            'bezout_only_B2':float(((BALL(uniform).log()-BALL(QQ(real['lower'])).log())/4).center())},
        'timing':{'real_cpu_seconds':real_cpu,'local_cpu_seconds':time.process_time()-local_start,
                  'total_cpu_seconds':time.process_time()-started},
        'boundary':'Uniform height bounds on this marked chart, including infinity. '
            'Finite upper bounds may include locally unattainable cases. No claim of optimal '
            'constants, global minimality, new points, search completeness under timeout, '
            'independence, or CPU improvement. Unfactored residuals retain a safe Bezout bound.'}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve the previous bound packet')
    data = json.loads(args.input.read_text())
    result = build(data)
    result['input_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['display']),flush=True)
