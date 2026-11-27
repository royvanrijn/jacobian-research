#!/usr/bin/env python3
"""Full-domain interval certificate for the finite-amplitude three-cell bound.

See THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md, FA1--FA7. The proof is global
in the partition amplitude, on 11/25 <= rho <= 221/500. This program is an
internal directed-arithmetic certificate, not independent or formal review.
All tables are rebuilt; no exploratory cache is an input.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

import numpy as np
from flint import arb, ctx

SOURCE_DIR = Path(__file__).resolve().parent / 'sources/cut_bounds_6b720fe57bfe'
sys.path.insert(0, str(SOURCE_DIR))
import check_positive_spectral_extension as s

I = s.Interval
s.LAMBDAS = (Q(0), Q(1,16), Q(1,4), Q(1), Q(4), Q(16), Q(64), Q(256))
s.COEFFICIENTS = tuple(Q(n, 1000000) for n in
    (-19606173,19029196,1042060,528898,3055,7376,-4335,-76))
ctx.prec = 100
pi = arb.pi()
j1 = (9-pi*arb(3).sqrt())/(4*pi*pi)
cw = (6-pi)/(4*pi*pi.sqrt())
cs, ls, M = s.trial_constants()
MOM = []
for n, row in enumerate(s.positive.CENTER):
    p = arb(n+1)/2
    MOM.append(sum((c*p.gamma()/(1+l)**p for c,l in zip(cs,ls)),arb(0))
               -M*(1+arb(n)/2).gamma())
U = ((1/arb(2).sqrt(),1/arb(6).sqrt()),
     (-1/arb(2).sqrt(),1/arb(6).sqrt()),(arb(0),-2/arb(6).sqrt()))


class Jet:

    def __init__(self, v, d=0, dd=0):
        self.v = I.coerce(v)
        self.d = I.coerce(d)
        self.dd = I.coerce(dd)

    @staticmethod
    def co(x):
        return x if isinstance(x, Jet) else Jet(x)

    def __add__(self, o):
        o = self.co(o)
        return Jet(self.v + o.v, self.d + o.d, self.dd + o.dd)
    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d, -self.dd)

    def __sub__(self, o):
        return self + -self.co(o)

    def __rsub__(self, o):
        return self.co(o) + -self

    def __mul__(self, o):
        o = self.co(o)
        return Jet(self.v * o.v, self.d * o.v + self.v * o.d, self.dd * o.v + 2 * self.d * o.d + self.v * o.dd)
    __rmul__ = __mul__

    def inverse(self):
        inv = I(1) / self.v
        inv2 = inv.square()
        return Jet(inv, -self.d * inv2, 2 * self.d.square() * inv2 * inv - self.dd * inv2)

    def __truediv__(self, o):
        return self * self.co(o).inverse()

    def sq(self):
        return self * self

def table_jet(tab, key, q):
    value = tab.query(key, q.v)
    mid = (q.v.lo + q.v.hi) / 2
    ix = np.clip(np.searchsorted(tab.edges, mid, side='right') - 1, 0, len(tab.edges) - 2)
    left = tab.edges[ix]
    right = tab.edges[ix + 1]
    vl, vh = tab.values[key]
    slope = (I(vl[ix + 1], vh[ix + 1]) - I(vl[ix], vh[ix])) / I(s.base.sub_lo(right, left), s.base.sub_hi(right, left))
    last = np.maximum(ix + 1, np.searchsorted(tab.edges, q.v.hi, side='left'))
    curv = tab.curvatures[key][last]
    distance = np.maximum((q.v - left).magnitude(), (q.v - right).magnitude())
    error = s.base.mul_pos_hi(curv, distance)
    der = I(s.base.sub_lo(slope.lo, error), s.base.add_hi(slope.hi, error))
    second = I(-curv, curv)
    return Jet(value, der * q.d, second * q.d.square() + der * q.dd)

def exp_minus(et, z):
    z.v.lo = np.maximum(0, z.v.lo)
    value = I(*et.query('E', z.v.lo, z.v.hi))
    return Jet(value, -value * z.d, value * (z.d.square() - z.dd))

class Mesh:

    def __init__(self, N):
        edges = np.linspace(0, 7, N + 1)
        mids = (edges[:-1] + edges[1:]) / 2
        self.width = I(s.base.sub_lo(edges[1:], edges[:-1]), s.base.sub_hi(edges[1:], edges[:-1]))
        self.x = Jet(I(edges[:-1], edges[1:]), 1)
        self.mid = Jet(I(mids), 1)
        cs, ls, M = s.trial_constants()

        def hjet(lo, hi):
            vv = []
            dd = []
            ss = []
            for l, u in zip(lo, hi):
                x = s.base.ball(float(l), float(u))
                v = sum((c * (-a * x * x).exp() for c, a in zip(cs, ls)), s.arb(0)) - M * x
                d = sum((-2 * a * c * x * (-a * x * x).exp() for c, a in zip(cs, ls)), s.arb(0)) - M
                d2 = sum((c * (4 * a * a * x * x - 2 * a) * (-a * x * x).exp() for c, a in zip(cs, ls)), s.arb(0))
                vv.append(s.base.bounds(v))
                dd.append(s.base.bounds(d))
                ss.append(s.base.bounds(d2))
            return Jet(I(*np.array(vv).T), I(*np.array(dd).T), I(*np.array(ss).T))
        self.h = hjet(edges[:-1], edges[1:])
        self.hm = hjet(mids, mids)

def sqrtjet(z):
    val = I.arb(s.base.ball(float(z.v.lo), float(z.v.hi)).sqrt())
    return Jet(val, z.d / (2 * val), z.dd / (2 * val) - z.d.square() / (4 * val * val * val))

def power(z, n):
    out = Jet(1)
    for _ in range(n):
        out = out * z
    return out

def center(rho):
    ans = Jet(0)
    for moment, row in zip(MOM, s.positive.CENTER):
        c = Jet(0)
        for k, v in enumerate(row):
            c += Jet(I.arb(s.arb(v))) * power(rho, k)
        ans += c * Jet(I.arb(moment))
    return ans

def components(rho, x, h, table, et, refU=None, delta=None):
    d = 1 - rho.sq()
    root = sqrtjet(d)
    alpha = rho / root
    kappa = rho.sq() / d
    q = alpha * x
    q.v.lo = np.maximum(0, q.v.lo)
    jb = alpha * (table_jet(table, 'HPLUS', q) + Jet(I.arb(j1 / 2)) * table_jet(table, 'S', q))
    jb.v.lo = np.maximum(jb.v.lo, (I.arb(cw) * alpha.v).lo)
    g = Jet(0)
    for c, l in zip(cs, ls):
        w = 1 + Jet(I.arb(l)) * d
        v = q / sqrtjet(w)
        v.v.lo = np.maximum(0, v.v.lo)
        term = Jet(I.arb(pi.sqrt())) / (root * sqrtjet(w)) * table_jet(table, 'L', v) - Jet(I.arb(4 * j1)) / (root * w) * q * table_jet(table, 'R', v)
        g += Jet(I.arb(c)) * term
    g -= Jet(I.arb(M)) * (table_jet(table, 'P', q) - Jet(I.arb(2 * j1)) * q * table_jet(table, 'S', q))
    em = exp_minus(et, kappa * x.sq())
    gauss = exp_minus(et, x.sq()) * em
    U = g - center(rho) * em.inverse() - 2 * jb * h
    if refU is not None:
        err = s.base.mul_pos_hi(U.d.magnitude(), delta)
        lo = s.base.sub_lo(refU.lo, err)
        hi = s.base.add_hi(refU.hi, err)
        U.v = I(np.maximum(U.v.lo, lo), np.minimum(U.v.hi, hi))
        assert np.all(U.v.lo <= U.v.hi)
    return (jb, gauss, U)

def integrands(jb, gauss, U, h):
    return [2 * gauss * jb * h.sq(), 2 * gauss * h * U, 2 * gauss * U.sq() / jb]

def tails(rlo, rhi):
    rho = s.base.ball(rlo, rhi)
    d = 1 - rho * rho
    alpha = rho / d.sqrt()
    k = rho * rho / d
    mass = sum((abs(c) for c in cs), s.arb(0))
    ft = mass / 7 + abs(M)
    jm = s.arb(1) / 2 + j1 / (2 * s.arb(1).exp())
    km = 2 + 2 * j1 / s.arb(1).exp()
    c0 = s.centered_trial_integral(rlo, rhi)
    U = km * (mass * pi.sqrt() / d.sqrt() + abs(M)) + abs(c0) + 2 * jm * ft
    gt = (-s.arb(49)).exp()
    return (s.base.bounds(jm * ft * ft * gt)[1], s.base.bounds(ft * U * gt)[1], s.base.bounds(U * U * (-(1 - k) * 49).exp() / (cw * alpha * (1 - k) * 7))[1])

def derivative_bounds(lo, hi, mesh, table, et):
    mid = (lo + hi) / 2
    delta = max(float(s.base.sub_hi(mid, lo)), float(s.base.sub_hi(hi, mid)))
    x = Jet(mesh.x.v)
    h = Jet(mesh.h.v)
    ux = components(Jet(I(mid)), mesh.x, mesh.h, table, et)[2]
    um = components(Jet(I(mid)), mesh.mid, mesh.hm, table, et)[2].v
    errx = s.base.mul_pos_hi(ux.d.magnitude(), s.base.div_pos_hi(mesh.width.hi, 2))
    ref = I(np.maximum(ux.v.lo, s.base.sub_lo(um.lo, errx)), np.minimum(ux.v.hi, s.base.add_hi(um.hi, errx)))
    assert np.all(ref.lo <= ref.hi)
    comp = components(Jet(I(lo, hi), 1), x, h, table, et, ref, delta)
    ints = integrands(*comp, h)
    return [float(s.integrate(I(0, f.d.magnitude()), mesh.width).hi) for f in ints]

def point_integrals(r, mesh, table, et):
    rho = Jet(I(r))
    mid = integrands(*components(rho, mesh.mid, mesh.hm, table, et), mesh.hm)
    whole = integrands(*components(rho, mesh.x, mesh.h, table, et), mesh.h)
    out = []
    for f, w in zip(mid, whole):
        err = s.base.div_pos_hi(s.base.mul_pos_hi(w.dd.magnitude(), mesh.width.square().hi), 24)
        out.append(s.integrate(I(s.base.sub_lo(f.v.lo, err), s.base.add_hi(f.v.hi, err)), mesh.width))
    return out

def cover(lo, hi, mesh, table, et, panels=100):
    L = derivative_bounds(lo, hi, mesh, table, et)
    tail = tails(lo, hi)
    print('LIPSCHITZ', L, 'TAILS', tail, flush=True)
    edges = np.linspace(lo, hi, panels + 1)
    result = []
    for l, h in zip(edges[:-1], edges[1:]):
        mid = (l + h) / 2
        delta = max(float(s.base.sub_hi(mid, l)), float(s.base.sub_hi(h, mid)))
        v = point_integrals(mid, mesh, table, et)
        v = [w + I(-float(s.base.mul_pos_hi(c, delta)), float(s.base.mul_pos_hi(c, delta))) for w, c in zip(v, L)]
        D = v[0] + I(0, tail[0])
        II = v[1] + I(-tail[1], tail[1])
        B = v[2] + I(0, tail[2])
        a = I(2) + II / D
        r2 = float(s.base.div_pos_hi(B.hi, D.lo))
        result.append(dict(rho=[float(l), float(h)], D=[float(D.lo), float(D.hi)], a=[float(a.lo), float(a.hi)], r2=r2))
    return dict(derivatives=L, tail=tail, panels=panels, minD=min((v['D'][0] for v in result)), maxD=max((v['D'][1] for v in result)), mina=min((v['a'][0] for v in result)), maxa=max((v['a'][1] for v in result)), maxr2=max((v['r2'] for v in result)), slabs=result)

def ba(l, h):
    lo = arb(l.numerator) / l.denominator
    hi = arb(h.numerator) / h.denominator
    return lo.union(hi)

def upper(x):
    return math.nextafter(float(x.upper()), math.inf)

def lower(x):
    return math.nextafter(float(x.lower()), -math.inf)

def sinc(x):
    z = x * x
    p = arb((-1) ** 7) / math.factorial(15)
    for k in reversed(range(7)):
        p = p * z + arb((-1) ** k) / math.factorial(2 * k + 1)
    return p + (pi / 6) ** 16 / math.factorial(17) * arb(0, 1)

def gamma_box(box, c=Q(1, 200)):
    (pl, ph), (ql, qh), (rl, rh), (sl, sh) = box
    p3l = max(Q(1, 6), 1 - ph - qh)
    p3h = min(Q(1, 3), 1 - pl - ql)
    q3l = max(Q(1, 6), 1 - rh - sh)
    q3h = min(Q(1, 2), 1 - rl - sl)
    if p3l > p3h or q3l > q3h or ph < ql or (qh < p3l):
        return None
    p = (ba(pl, ph), ba(ql, qh), ba(p3l, p3h))
    q = (ba(rl, rh), ba(sl, sh), ba(q3l, q3h))
    a = [(pi * x).sin() for x in p]
    b = [(pi * x).sin() for x in q]
    cb = [(pi * (x + y) / 2).cos() * sinc(pi * (x - y) / 2) for x, y in zip(p, q)]

    def one(p, ratios):
        eta = [(ratios[i] - ratios[(i + 1) % 3]) * (ratios[i] - ratios[(i + 2) % 3]) / ratios[i] for i in range(3)]
        M = max([0.0] + [upper(v) for v in eta])
        den = upper(2 * (sum(ratios, arb(0)) + arb(M)))
        G = [[cb[i] * cb[j] * (arb(1) if i == j else -(pi * p[3 - i - j]).cos()) for j in range(3)] for i in range(3)]
        return [[sum((U[i][v] * G[i][j] * U[j][w] for i in range(3) for j in range(3)), arb(0)) / den for w in range(2)] for v in range(2)]
    A = one(p, [y / x for x, y in zip(a, b)])
    B = one(q, [x / y for x, y in zip(a, b)])
    M = [[(A[i][j] + B[i][j]) / 2 for j in range(2)] for i in range(2)]
    cc = arb(c.numerator) / c.denominator
    f = M[0][0] - cc
    g = M[1][1] - cc
    det = f * g - M[0][1] * M[1][0]
    return (lower(f), lower(g), lower(det))

def gamma_cover(limit=100000):
    start = time.monotonic()
    root = ((Q(1, 3), Q(1, 2)), (Q(1, 4), Q(5, 12)), (Q(1, 6), Q(1, 2)), (Q(1, 6), Q(1, 2)))
    pending = [root]
    accepted = skipped = attempts = 0
    maxdepth = 0
    mindet = 1.0
    minlead = 1.0
    while pending:
        box = pending.pop()
        ans = gamma_box(box)
        attempts += 1
        if ans is None:
            skipped += 1
        elif min(ans) > 0:
            accepted += 1
            mindet = min(mindet, ans[2])
            minlead = min(minlead, ans[0], ans[1])
        else:
            if attempts >= limit:
                return dict(status='incomplete', attempts=attempts, pending=len(pending), accepted=accepted, skipped=skipped, box=[[str(a), str(b)] for a, b in box], margins=ans, elapsed=time.monotonic() - start)
            i = max(range(4), key=lambda j: box[j][1] - box[j][0])
            lo, hi = box[i]
            mid = (lo + hi) / 2
            x = list(box)
            y = list(box)
            x[i] = (lo, mid)
            y[i] = (mid, hi)
            pending.extend([tuple(y), tuple(x)])
        if attempts % 10000 == 0:
            print(json.dumps(dict(progress=attempts, pending=len(pending), accepted=accepted, elapsed=time.monotonic() - start)), flush=True)
    return dict(status='pass', attempts=attempts, accepted=accepted, skipped=skipped, minimum_leading_margin=minlead, minimum_determinant_margin=mindet, elapsed=time.monotonic() - start)


def interval_covers(records, target=(Q(11, 25), Q(221, 500))):
    """Exact rational audit of the retained binary64 interval endpoints."""
    current = target[0]
    for item in sorted(records, key=lambda z: z['rho'][0]):
        lo, hi = map(Q.from_float, item['rho'])
        if lo > current or lo >= hi:
            return False
        current = max(current, hi)
    return current >= target[1]


def amplitude_constants():
    """Arb checks of the unrestricted analytic assembly in FA3--FA6."""
    rat = lambda n, d=1: arb(n) / d
    rho = s.base.ball(math.nextafter(.44, -math.inf),
                      math.nextafter(.442, math.inf))
    d = 1-rho*rho
    alpha = rho/d.sqrt()
    kappa = rho*rho/d
    x0 = rat(1, 10)
    q = alpha*x0
    S = lambda z: pi.sqrt()*(1+2*z*z)*(z*z).exp()*z.erfc()-2*z
    weighted_j = j1*alpha*(-(1+kappa)*x0*x0).exp()*(S(q/2)+S(q)/2)
    positive_slope = sum((c*l for c, l in zip(cs, ls) if c > 0), arb(0))
    pole = (sum(cs, arb(0))-positive_slope*x0*x0)/rat(271, 1000).sqrt()
    constant = abs(M)/rat(268, 1000).sqrt()
    A0, c, C, j = rat(27, 1000), rat(47, 25), rat(14, 25), rat(107, 1000)
    gap, e, R2, P = rat(7, 5), rat(31, 10000), rat(1, 50000), rat(97, 20)
    H2 = ((R2.sqrt()+(R2+gap*e).sqrt())/gap)**2
    B = rat(1, 3)+A0*C/arb(6).sqrt()
    X = A0*c/(arb(6).sqrt()*B)
    Phi = j*c**3*A0/(9*arb(6).sqrt()*B)
    zeta = rat(1, 6)-rat(2, 3).sqrt()*A0*P
    anchor_margin = zeta*zeta/20-(2*P*A0)**2/400
    rest_margin = rat(3, 5)**2-(rat(217, 50)-rat(399, 200)**2)
    determinant = rat(13, 20000)*rat(89, 135)-R2
    assert weighted_j > j and pole > c and constant < C
    assert X < x0 and Phi > H2 and zeta > 0 and anchor_margin > 0
    assert j1 > rat(9, 100) and rest_margin > 0 and determinant > 0
    values = dict(weighted_J_local=weighted_j, pole_lower=pole,
                  pole_constant_upper=constant, feasibility_H_squared=H2,
                  feasibility_Phi=Phi, feasibility_margin=Phi-H2,
                  integration_endpoint=X, outside_cube_margin=anchor_margin,
                  orthogonal_block_margin=rest_margin,
                  final_determinant=determinant, j1=j1)
    return {name: list(s.base.bounds(value)) for name, value in values.items()}


def angular_band(exp_table, panels=2048):
    """Positive interval Riemann sum over the entire dyadic annular square."""
    rho = I(math.nextafter(.44, -math.inf), math.nextafter(.442, math.inf))
    d = 1-rho.square()
    edges = np.linspace(5/16, 2., panels+1)
    xx = I(edges[:-1], edges[1:])
    widths = I(s.base.sub_lo(edges[1:], edges[:-1]),
               s.base.sub_hi(edges[1:], edges[:-1]))
    values = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = s.base.ball(float(lo), float(hi))
        h = sum((c*(-l*x*x).exp() for c, l in zip(cs, ls)), arb(0))-M*x
        values.append(s.base.bounds(h))
    hh = I(*np.array(values).T)
    normalized = (hh/xx)/I.arb((arb(268)/1000).sqrt())
    phi_upper = float(np.max(normalized.magnitude()))
    lower_sum, upper_sum = 0., 0.
    for start in range(0, panels, 64):
        stop = start+64
        x = I(xx.lo[start:stop, None], xx.hi[start:stop, None])
        y = I(xx.lo[None, :], xx.hi[None, :])
        h = I(hh.lo[start:stop, None], hh.hi[start:stop, None])
        k = I(hh.lo[None, :], hh.hi[None, :])
        z = (x.square()+y.square()+2*rho*x*y)/d
        exponential = I(*exp_table.query('E', z.lo, z.hi))
        integrand = (8*rho/(d*d))*exponential*(y*h-x*k).square()
        area = I(widths.lo[start:stop, None], widths.hi[start:stop, None])*I(
            widths.lo[None, :], widths.hi[None, :])
        value = integrand*area
        lo, _ = s.base.sum_enclosure(value.lo)
        _, hi = s.base.sum_enclosure(value.hi)
        lower_sum = float(s.base.add_lo(lower_sum, lo))
        upper_sum = float(s.base.add_hi(upper_sum, hi))
    norm = I.arb(arb(271)/1000)
    normalized_lower = float(s.base.div_pos_lo(lower_sum, norm.hi))
    assert Q.from_float(phi_upper) < Q(97, 20)
    assert Q.from_float(normalized_lower) > Q(3, 2)
    return dict(band=[5/16, 2.], panels=panels,
                raw_integral=[lower_sum, upper_sum],
                normalized_dirichlet_lower=normalized_lower,
                normalized_trial_absolute_upper=phi_upper)


def numerical_controls(table):
    """Supplementary controls; the unrestricted error argument is in FA7."""
    checks = s.finite_controls(table)
    # Midpoint error is width^3 sup|F''|/24. For x^2 on [0,1] it is exact.
    assert Q(1, 4)+Q(2, 24) == Q(1, 3)
    assert Q(1, 4)+Q(2, 25) < Q(1, 3)
    assert gamma_box(((Q(1, 3), Q(1, 3)),)*4, Q(1)) is not None
    assert min(gamma_box(((Q(1, 3), Q(1, 3)),)*4, Q(1))) < 0
    checks.update(midpoint_remainder_control=True,
                  smaller_midpoint_remainder_rejected=True,
                  excessive_angular_coercivity_rejected=True)
    return checks


def audit_result(spectral, hs, gamma, band):
    assert gamma['status'] == 'pass'
    assert gamma['attempts'] == 2*(gamma['accepted']+gamma['skipped'])-1
    assert interval_covers(hs)
    children = []
    for parent in spectral:
        assert interval_covers(parent['slabs'], tuple(map(Q.from_float, parent['rho'])))
        for child in parent['slabs']:
            assert parent['rho'][0] <= child['rho'][0] < child['rho'][1] <= parent['rho'][1]
            lo, hi = map(Q.from_float, child['D'])
            assert Q(268, 1000) <= lo <= hi <= Q(271, 1000)
            lo, hi = map(Q.from_float, child['a'])
            assert Q(399, 200) <= lo <= hi <= Q(20031, 10000)
            assert Q.from_float(child['r2']) <= Q(1, 50000)
        children.extend(parent['slabs'])
    assert interval_covers(children)
    assert not interval_covers(children[1:])
    assert not interval_covers(hs[1:])
    assert all(Q.from_float(item['hs_squared']) <= Q(217, 50) for item in hs)
    return dict(exact_rational_endpoint_cover=True, child_parent_containment=True,
                all_spectral_caps=True, all_hilbert_schmidt_caps=True,
                incomplete_spectral_cover_rejected=True,
                incomplete_hilbert_schmidt_cover_rejected=True,
                gamma_full_bisection_cover=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    start = time.monotonic()
    print('Rebuilding validated one-variable tables', flush=True)
    table = s.InterpolatedTable(4., 8192)
    et = s.base.build_exp_table(64., 131072, 2.)
    controls = numerical_controls(table)
    constants = amplitude_constants()
    gamma = gamma_cover()
    print('Angular cube', json.dumps(gamma), flush=True)
    band = angular_band(et)
    print('Angular annulus', json.dumps(band), flush=True)
    mesh = Mesh(16384)
    spectral = []
    for lo, hi in zip(np.linspace(.44, .442, 5)[:-1], np.linspace(.44, .442, 5)[1:]):
        lo, hi = math.nextafter(float(lo), -math.inf), math.nextafter(float(hi), math.inf)
        parent = cover(lo, hi, mesh, table, et, 25)
        parent['rho'] = [lo, hi]
        spectral.append(parent)
        print('Spectral parent', json.dumps({k: v for k, v in parent.items() if k != 'slabs'}), flush=True)
    print('Rebuilding Hilbert--Schmidt tables', flush=True)
    qt = s.positive.build_positive_q_table(4., 1024, 64, 2.)
    at = s.positive.build_positive_a_table(64., 4096, 64, 2.)
    upstream_constants = s.positive.build_constants()
    edges = np.linspace(.44, .442, 21)
    hs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        lo, hi = math.nextafter(float(lo), -math.inf), math.nextafter(float(hi), math.inf)
        value = s.positive.certify_slab(lo, hi, 7., 6400, 48, at, qt, et, upstream_constants)
        item = dict(rho=[lo, hi], hs_squared=value.total, finite=value.finite,
                    tail=value.tail, kernel_tail=value.kernel_tail,
                    additive_tail=value.additive_tail)
        hs.append(item)
        print('Hilbert--Schmidt parent', json.dumps(item), flush=True)
    controls.update(audit_result(spectral, hs, gamma, band))
    summary = dict(minimum_trial_norm=min(v['minD'] for v in spectral),
                   maximum_trial_norm=max(v['maxD'] for v in spectral),
                   minimum_rayleigh=min(v['mina'] for v in spectral),
                   maximum_rayleigh=max(v['maxa'] for v in spectral),
                   maximum_raw_residual_squared=max(v['maxr2'] for v in spectral),
                   maximum_hilbert_schmidt_squared=max(v['hs_squared'] for v in hs))
    input_files = ('plur_cert.py', 'plur_cert_positive.py', 'check_positive_spectral_extension.py')
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result = dict(status='pass', exact_correlation_range=['11/25', '221/500'],
                  coefficients=[str(c) for c in s.COEFFICIENTS],
                  exponents=[str(l) for l in s.LAMBDAS],
                  analytic_trial_mean='0', summary=summary,
                  spectral_cover=spectral, hilbert_schmidt_cover=hs,
                  angular_cube=gamma, angular_band=band,
                  analytic_assembly=constants, controls=controls,
                  input_sources=[dict(path=str(SOURCE_DIR.name+'/'+name),
                                      sha256=digest(SOURCE_DIR/name)) for name in input_files],
                  checker_sha256=digest(Path(__file__)),
                  environment=dict(python=sys.version.split()[0],
                                   python_flint=__import__('flint').__version__, numpy=np.__version__),
                  elapsed_seconds=time.monotonic()-start)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status='pass', summary=summary, elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
