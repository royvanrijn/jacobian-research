"""Directed spectral refinement of the positive radial certificate.

The retained upstream modules are imported without modification. This program
does not assert a mathematical result unless every covering slab passes. The
trial is an auxiliary Hilbert-space vector, not a partition mass profile.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import sys
import time

import numpy as np
from flint import arb, acb, ctx

import plur_cert as base
import plur_cert_positive as positive


LAMBDAS = (Fraction(0), Fraction(1, 4), Fraction(1), Fraction(4), Fraction(16))
COEFFICIENTS = (Fraction(151, 250), Fraction(-1), Fraction(28, 125),
                Fraction(-33, 1250), Fraction(121, 2500))


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = np.asarray(lo if hi is None else hi, dtype=float)

    @classmethod
    def arb(cls, value):
        return cls(*base.bounds(value))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Interval(base.add_lo(self.lo, other.lo), base.add_hi(self.hi, other.hi))

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Interval(*base.mul_interval(self.lo, self.hi, other.lo, other.hi))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if np.any(other.lo <= 0):
            raise ArithmeticError('positive denominator required')
        return self * Interval(base.down(1 / other.hi), base.up(1 / other.lo))

    def square(self):
        hi = base.abs_square_upper(self.lo, self.hi)
        lo = np.where((self.lo <= 0) & (self.hi >= 0), 0,
                      np.minimum(base.down(self.lo * self.lo), base.down(self.hi * self.hi)))
        return Interval(np.maximum(0, lo), hi)

    def magnitude(self):
        return np.maximum(abs(self.lo), abs(self.hi))


def rational(value):
    return arb(value.numerator) / value.denominator


def trial_constants():
    pi = arb.pi()
    cs = [rational(c) for c in COEFFICIENTS]
    ls = [rational(lam) for lam in LAMBDAS]
    mean = pi.sqrt() * sum((c / (1 + lam).sqrt() for c, lam in zip(cs, ls)), arb(0))
    return cs, ls, mean


def point_values(q):
    """Validated complex integration of entire functions; q is a point ball."""
    pi = arb.pi()
    sqrtpi = pi.sqrt()
    end = arb(2) / 3
    opts = dict(abs_tol=arb('1e-20'), rel_tol=arb('1e-20'), eval_limit=10000)

    def erfcx_minus(z):
        return (z * z).exp() * (-z).erfc()

    def lint(v, analytic):
        z = q * (pi * v).cos()
        return (2 + 3 * v) * erfcx_minus(z) / 2

    def pint(v, analytic):
        z = q * (pi * v).cos()
        return (2 + 3 * v) * (1 + sqrtpi * z * erfcx_minus(z)) / 2

    def hint(v, analytic):
        z = -q * (pi * v).cos()
        s = sqrtpi * (1 + 2 * z * z) * (z * z).exp() * z.erfc() - 2 * z
        return pi * v * (1 - 3 * v / 2) * (pi * v).sin() * s / 2

    out = {}
    for name, fn in [('L', lint), ('P', pint), ('HPLUS', hint)]:
        value = acb.integral(fn, 0, end, **opts)
        if not value.is_finite() or not value.imag.contains(0):
            raise ArithmeticError('failed point integral ' + name)
        out[name] = value.real
    f = (q * q).exp() * q.erfc()
    out['S'] = sqrtpi * (1 + 2 * q * q) * f - 2 * q
    out['R'] = 1 - sqrtpi * q * f
    return out


def curvature_bounds(q):
    """Increasing explicit upper bounds for absolute second derivatives."""
    pi = arb.pi()
    root = pi.sqrt()
    exponential = (q * q).exp()
    e1 = 4 * q * exponential + 2 / root
    e2 = (2 + 4 * q * q) * 2 * exponential + 4 * q / root
    s2 = root * (6 + 24 * q * q + 8 * q**4) * 2 * exponential + 20 * q + 8 * q**3
    return {'L': e2, 'P': root * (2 * e1 + q * e2),
            'HPLUS': pi * s2 / 27, 'S': s2,
            'R': root * (2 * e1 + q * e2)}


class InterpolatedTable:
    """Point enclosures plus Lagrange remainder, including extrapolation."""
    def __init__(self, maximum, bins):
        self.edges = base.graded_edges(maximum, bins, 2.0)
        self.values = {key: (np.empty(bins + 1), np.empty(bins + 1))
                       for key in ('L', 'P', 'HPLUS', 'S', 'R')}
        self.curvatures = {key: np.empty(bins + 1) for key in self.values}
        for i, edge in enumerate(self.edges):
            q = arb(float(edge))
            values = point_values(q)
            bounds = curvature_bounds(q)
            for key in self.values:
                self.values[key][0][i], self.values[key][1][i] = base.bounds(values[key])
                self.curvatures[key][i] = base.bounds(bounds[key])[1]
        self.maximum = float(self.edges[-1])

    def query(self, key, q):
        if np.any(q.lo < 0) or np.any(q.hi > self.maximum):
            raise ArithmeticError('query outside interpolation table')
        midpoint = (q.lo + q.hi) / 2
        index = np.clip(np.searchsorted(self.edges, midpoint, side='right') - 1,
                        0, len(self.edges) - 2)
        left = self.edges[index]
        right = self.edges[index + 1]
        width = Interval(base.sub_lo(right, left), base.sub_hi(right, left))
        vals_lo, vals_hi = self.values[key]
        value_left = Interval(vals_lo[index], vals_hi[index])
        value_right = Interval(vals_lo[index + 1], vals_hi[index + 1])
        linear = value_left + ((q - left) / width) * (value_right - value_left)
        # The interpolation remainder holds on the convex hull of q and
        # the two nodes even if q crosses an interpolation cell boundary.
        last = np.maximum(index + 1, np.searchsorted(self.edges, q.hi, side='left'))
        second = self.curvatures[key][last]
        error = base.mul_pos_hi(0.5, base.mul_pos_hi(second,
                   base.mul_pos_hi((q - left).magnitude(), (q - right).magnitude())))
        return Interval(np.maximum(0, base.sub_lo(linear.lo, error)),
                        base.add_hi(linear.hi, error))


def centered_trial_integral(rho_lo, rho_hi):
    cs, ls, mean = trial_constants()
    rho = base.ball(rho_lo, rho_hi)
    answer = arb(0)
    for n, row in enumerate(positive.CENTER):
        coefficient = sum((arb(c) * rho**i for i, c in enumerate(row)), arb(0))
        power = arb(n + 1) / 2
        first = sum((c * power.gamma() / (1 + lam)**power
                     for c, lam in zip(cs, ls)), arb(0))
        second = mean * (1 + arb(n) / 2).gamma()
        answer += coefficient * (first - second)
    return answer


def integrate(values, widths):
    product = values * widths
    lo, _ = base.sum_enclosure(product.lo)
    _, hi = base.sum_enclosure(product.hi)
    return Interval(lo, hi)


def temple_bound(hs_squared, rayleigh_lo, rayleigh_hi, residual_squared):
    if not rayleigh_lo > 0:
        raise ArithmeticError('leading eigenvalue not isolated')
    rest_square = base.sub_hi(hs_squared, base.mul_pos_lo(rayleigh_lo, rayleigh_lo))
    rest = base.sqrt_pos_hi(max(0, float(rest_square)))
    if not rayleigh_lo > rest:
        raise ArithmeticError('leading eigenvalue not isolated')
    upper = base.add_hi(rayleigh_hi, base.div_pos_hi(residual_squared,
                        base.sub_lo(rayleigh_lo, rest)))
    return rest, float(upper)


def finite_controls(table):
    # A self-adjoint two-dimensional positive control with an exactly
    # normalized rational vector, and a rejected absent spectral gap.
    a = arb(19) / 10
    b = arb(1) / 5
    w = arb(99)**2 / 101**2
    rayleigh = w * a + (1 - w) * b
    raw_residual = w * (a - 2)**2 + (1 - w) * (b - 2)**2
    h2 = a*a + b*b
    rest, upper = temple_bound(base.bounds(h2)[1], *base.bounds(rayleigh),
                               base.bounds(raw_residual)[1])
    assert upper >= base.bounds(a)[1] and upper < 2
    try:
        temple_bound(4, 0.1, 0.2, 0.01)
    except ArithmeticError:
        pass
    else:
        raise AssertionError('absent isolation must be rejected')
    points = (0.0, float(table.edges[len(table.edges)//3]),
              float((table.edges[len(table.edges)//2] + table.edges[len(table.edges)//2+1])/2))
    for q in points:
        direct = point_values(arb(q))
        for key in direct:
            enclosed = table.query(key, Interval(q))
            lo, hi = base.bounds(direct[key])
            assert float(enclosed.lo) <= lo and float(enclosed.hi) >= hi, (key, q)
    try:
        table.query('L', Interval(-0.01))
    except ArithmeticError:
        pass
    else:
        raise AssertionError('negative table argument must be rejected')
    return dict(rational_matrix_spectral_upper=upper,
                absent_isolation_rejected=True, interpolation_control_points=len(points),
                negative_table_argument_rejected=True)


def verify_cover(results, lower, upper):
    position = Fraction(str(lower))
    for item in sorted(results, key=lambda r: r['rho_lo']):
        if Fraction(item['rho_lo']) > position or not item['spectral_upper'] < 2:
            return False
        position = max(position, Fraction(item['rho_hi']))
    return position >= Fraction(str(upper))


class TrialMesh:
    def __init__(self, cutoff, panels):
        edges = np.linspace(0, cutoff, panels + 1)
        edges[0] = 0.0
        edges[-1] = cutoff
        self.x = Interval(edges[:-1], edges[1:])
        self.width = Interval(base.sub_lo(edges[1:], edges[:-1]),
                              base.sub_hi(edges[1:], edges[:-1]))
        cs, ls, mean = trial_constants()
        hlo, hhi = np.empty(panels), np.empty(panels)
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            x = base.ball(float(lo), float(hi))
            h = sum((c * (-lam * x * x).exp() for c, lam in zip(cs, ls)), arb(0)) - mean * x
            hlo[i], hhi[i] = base.bounds(h)
        self.h = Interval(hlo, hhi)


def spectral_slab(rho_lo, rho_hi, cutoff, mesh, table, exp_table, constants, hs_squared):
    rho = base.ball(rho_lo, rho_hi)
    d = 1 - rho * rho
    rootd = d.sqrt()
    alpha = rho / rootd
    kappa = rho * rho / d
    pi = arb.pi()
    j1 = (9 - pi * arb(3).sqrt()) / (4 * pi * pi)
    cw = (6 - pi) / (4 * pi * pi.sqrt())
    cs, ls, mean = trial_constants()
    q = Interval.arb(alpha) * mesh.x
    q.lo = np.maximum(0, q.lo)  # q is analytically nonnegative.
    jbar = Interval.arb(alpha) * (table.query('HPLUS', q) +
                                  Interval.arb(j1 / 2) * table.query('S', q))
    jbar.lo = np.maximum(jbar.lo, base.bounds(cw * alpha)[0])
    g = Interval(0)
    for coefficient, lam in zip(cs, ls):
        weight = 1 + lam * d
        v = q / Interval.arb(weight.sqrt())
        v.lo = np.maximum(0, v.lo)
        term = Interval.arb(pi.sqrt() / (rootd * weight.sqrt())) * table.query('L', v)
        term -= Interval.arb(4 * j1 / (rootd * weight)) * q * table.query('R', v)
        g += Interval.arb(coefficient) * term
    g -= Interval.arb(mean) * (table.query('P', q) -
                               Interval.arb(2 * j1) * q * table.query('S', q))
    c0 = centered_trial_integral(rho_lo, rho_hi)
    square_x = mesh.x.square()
    kz = Interval.arb(kappa) * square_x
    kz.lo = np.maximum(0, kz.lo)
    emin = Interval(*exp_table.query('E', kz.lo, kz.hi))
    eplus = Interval(1) / emin
    gauss = Interval(*exp_table.query('E', square_x.lo, square_x.hi)) * emin
    residual = g - Interval.arb(c0) * eplus - 2 * jbar * mesh.h
    denom = integrate(2 * gauss * jbar * mesh.h.square(), mesh.width)
    correction = integrate(2 * gauss * mesh.h * residual, mesh.width)
    defect = integrate(2 * gauss * residual.square() / jbar, mesh.width)

    # Absolute row bound from |tau_+-4theta|<=2+2j1/e and
    # sup_r P_rho(t^(-1/2))<=sqrt(pi/d).
    coefficient_mass = sum((abs(c) for c in cs), arb(0))
    f_tail = coefficient_mass / cutoff + abs(mean)
    jmax = arb(1) / 2 + j1 / (2 * arb(1).exp())
    kmax = 2 + 2 * j1 / arb(1).exp()
    row_bound = kmax * (coefficient_mass * pi.sqrt() / rootd + abs(mean))
    residual_bound = row_bound + abs(c0) + 2 * jmax * f_tail
    gaussian_tail = (-arb(cutoff)**2).exp()
    d_tail = base.bounds(jmax * f_tail**2 * gaussian_tail)[1]
    i_tail = base.bounds(f_tail * residual_bound * gaussian_tail)[1]
    rate = 1 - kappa
    if not rate > 0:
        raise ArithmeticError('residual tail rate is not positive')
    r_tail = base.bounds(residual_bound**2 * (-(rate * cutoff**2)).exp()
                         / (cw * alpha * rate * cutoff))[1]
    denom += Interval(0, d_tail)
    correction += Interval(-i_tail, i_tail)
    defect += Interval(0, r_tail)
    if not denom.lo > 0:
        raise ArithmeticError('trial norm lower bound not positive')
    rayleigh = Interval(2) + correction / denom
    residual_squared = base.div_pos_hi(defect.hi, denom.lo)
    rest, upper = temple_bound(hs_squared, float(rayleigh.lo), float(rayleigh.hi),
                               float(residual_squared))
    return dict(rho_lo=rho_lo, rho_hi=rho_hi, hs_squared=float(hs_squared),
                trial_norm=[float(denom.lo), float(denom.hi)],
                rayleigh=[float(rayleigh.lo), float(rayleigh.hi)],
                residual_squared_upper=float(residual_squared), rest_spectrum_upper=rest,
                spectral_upper=float(upper), residual_tail_upper=r_tail)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rho-min', type=float, default=0.42)
    parser.add_argument('--rho-max', type=float, default=0.44)
    parser.add_argument('--rho-panels', type=int, default=20)
    parser.add_argument('--cutoff', type=float, default=7.0)
    parser.add_argument('--trial-panels', type=int, default=262144)
    parser.add_argument('--interpolation-bins', type=int, default=4096)
    parser.add_argument('--hs-panels', type=int, default=800)
    parser.add_argument('--a-bins', type=int, default=4096)
    parser.add_argument('--q-bins', type=int, default=1024)
    parser.add_argument('--exp-bins', type=int, default=131072)
    parser.add_argument('--v-panels', type=int, default=64)
    parser.add_argument('--min-rho-width', type=float, default=0.0000001)
    parser.add_argument('--max-attempts', type=int, default=20000)
    parser.add_argument('--output')
    args = parser.parse_args()
    assert 0 < args.rho_min < args.rho_max < 0.5
    ctx.prec = 100
    start = time.perf_counter()
    rhomax = arb(math.nextafter(args.rho_max, math.inf))
    d = 1 - rhomax * rhomax
    # Integer headroom also encloses the subsequent binary64 interval
    # operations, rather than only the exact real parameter expressions.
    qmax = float(math.ceil(base.bounds(rhomax * args.cutoff / d.sqrt())[1]))
    amax = float(math.ceil(base.bounds(2 * rhomax * args.cutoff**2 / d)[1]) + 1)
    zmax = float(math.ceil(base.bounds(max(arb(args.cutoff**2), rhomax**2 / d * 2 * args.cutoff**2))[1]) + 1)
    print('building validated interpolation points', flush=True)
    table = InterpolatedTable(qmax, args.interpolation_bins)
    controls = finite_controls(table)
    print('building upstream HS tables', flush=True)
    qt = positive.build_positive_q_table(qmax, args.q_bins, args.v_panels, 2.0)
    at = positive.build_positive_a_table(amax, args.a_bins, args.v_panels, 2.0)
    et = base.build_exp_table(zmax, args.exp_bins, 2.0)
    constants = positive.build_constants()
    mesh = TrialMesh(args.cutoff, args.trial_panels)
    edges = np.linspace(args.rho_min, args.rho_max, args.rho_panels + 1)
    edges[0] = args.rho_min
    edges[-1] = args.rho_max
    results = []
    parent_covers = []
    attempts = 0
    stopped = False
    for parent_lo, parent_hi in zip(edges[:-1], edges[1:]):
        lo = math.nextafter(float(parent_lo), -math.inf)
        hi = math.nextafter(float(parent_hi), math.inf)
        hs = positive.certify_slab(lo, hi, args.cutoff, args.hs_panels, 48,
                                  at, qt, et, constants)
        parent_covers.append(dict(rho_lo=lo, rho_hi=hi, hs_squared=hs.total,
                                  finite=hs.finite, tail=hs.tail,
                                  kernel_tail=hs.kernel_tail, additive_tail=hs.additive_tail))
        pending = [(float(parent_lo), float(parent_hi))]
        while pending:
            child_lo, child_hi = pending.pop()
            attempts += 1
            expanded_lo = math.nextafter(child_lo, -math.inf)
            expanded_hi = math.nextafter(child_hi, math.inf)
            try:
                result = spectral_slab(expanded_lo, expanded_hi, args.cutoff, mesh,
                                       table, et, constants, hs.total)
            except ArithmeticError as exc:
                if str(exc) != 'leading eigenvalue not isolated':
                    raise
                result = None
            if result is not None and result['spectral_upper'] < 2:
                result['hs_parent_index'] = len(parent_covers) - 1
                print(json.dumps(result), flush=True)
                results.append(result)
            elif child_hi-child_lo > args.min_rho_width and attempts < args.max_attempts:
                midpoint = (child_lo+child_hi)/2
                pending.append((midpoint, child_hi))
                pending.append((child_lo, midpoint))
            else:
                print('unaccepted slab '+repr((child_lo, child_hi, result)), flush=True)
                stopped = True
                break
        if stopped:
            break
    passed = not stopped and verify_cover(results, args.rho_min, args.rho_max)
    # A deliberately incomplete cover cannot be promoted by successful
    # individual slabs.
    if passed and len(results) > 1:
        assert not verify_cover(results[1:], args.rho_min, args.rho_max)
    controls['incomplete_cover_rejected'] = bool(passed and len(results) > 1)
    output = dict(status='pass' if passed else 'fail', parameters=vars(args),
                  coefficients=[str(c) for c in COEFFICIENTS],
                  exponents=[str(lam) for lam in LAMBDAS],
                  trial_mean='0', slabs=results, hilbert_schmidt_cover=parent_covers,
                  attempts=attempts, controls=controls,
                  environment=dict(python=sys.version.split()[0],
                                   python_flint=__import__('flint').__version__,
                                   numpy=np.__version__),
                  maximum_spectral_upper=max(x['spectral_upper'] for x in results),
                  elapsed_seconds=time.perf_counter() - start)
    if args.output:
        with open(args.output, 'w') as file:
            json.dump(output, file, indent=2)
            file.write('\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'slabs'}), flush=True)
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
