#!/usr/bin/env python3
"""Directed cube and finite-amplitude controls for SHARP_TWO_MODE_ENERGY.

The unrestricted sharp one-set energy formula has a written proof. This
checker certifies the specified four-dimensional mass cube, its two Gram
coercivity constants, and the explicit rho=3/5 radial-family comparison.
It makes no all-correlation or all-profile optimality claim.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time
from flint import arb, ctx
from check_three_cell_coupled_radial import certificate as scalar_certificate, enclose

ctx.prec = 100
pi = arb.pi()
rt = arb(3).sqrt()
def ball(l, h):
    return (arb(l.numerator) / l.denominator).union(arb(h.numerator) / h.denominator)

def lo(x):
    return math.nextafter(float(x.lower()), -math.inf)

def gram(p, m):
    a = [(m * pi * x).sin() ** 2 for x in p]
    return [a[0] + a[1] - a[2] / 2, rt * (a[0] - a[1]) / 2, 3 * a[2] / 2]

def det(A):
    return A[0] * A[2] - A[1] * A[1]

def inverse(A):
    d = det(A)
    if not d > 0:
        raise ArithmeticError()
    return [A[2] / d, -A[1] / d, A[0] / d]

def trans(A, B):
    da, db = (det(A), det(B))
    if not (da > 0 and db > 0):
        raise ArithmeticError()
    h = (da * db).sqrt()
    F = A[0] * B[0] + 2 * A[1] * B[1] + A[2] * B[2] + 2 * h
    if not F > 0:
        raise ArithmeticError()
    F = F.sqrt()
    Ai = inverse(A)
    Bi = inverse(B)
    return [[(b + h * a) / F for a, b in zip(Ai, B)], [(a + h * b) / F for a, b in zip(A, Bi)]]

def lambdas(X):
    return [X[0] + rt * X[1], X[0] - rt * X[1], (3 * X[2] - X[0]) / 2]

def check(box, ratio=Q(1, 4), coercivities=(Q(1, 10), Q(3, 40))):
    ps = []
    for v in [box[:2], box[2:]]:
        (a, b), (c, d) = v
        l = max(Q(17, 60), 1 - b - d)
        h = min(Q(23, 60), 1 - a - c)
        if l > h:
            return None
        ps.append([ball(a, b), ball(c, d), ball(l, h)])
    p, q = ps
    try:
        X = trans(gram(p, 1), gram(q, 1))
        Y = trans(gram(p, 2), gram(q, 2))
    except ArithmeticError:
        return (-1.0, -1.0, -1.0)
    r = arb(ratio.numerator) / ratio.denominator
    first = second = margin = 100.0
    for z, x, y in zip([p, q], X, Y):
        lx, ly = (lambdas(x), lambdas(y))
        for i in range(3):
            first = min(first, lo(lx[i]))
            second = min(second, lo(ly[i]))
            margin = min(margin, lo(lx[i] - r * (pi * z[i]).sin() ** 2 * ly[i]))
    coerc = []
    for folded, constant in zip((False, True), coercivities):
        pp, qq = ([1-2*z for z in p], [1-2*z for z in q]) if folded else (p, q)
        w = [(pi*a).sin()*(pi*b).sin()/(pi*pi) for a, b in zip(pp, qq)]
        beta = pi*pi/2-pi**4/(600 if folded else 2400)
        v = [beta*w[(i+1)%3]*w[(i+2)%3]/sum(w, arb(0)) for i in range(3)]
        c = arb(constant.numerator)/constant.denominator
        a = (v[0]+v[1])/2-c
        d = (v[0]+v[1]+4*v[2])/6-c
        b = (v[0]-v[1])/(2*rt)
        coerc.extend((lo(a), lo(d), lo(a*d-b*b)))
    return (first, second, margin, *coerc)

def run(limit=1000000):
    start = time.monotonic()
    root = ((Q(17, 60), Q(23, 60)),) * 4
    todo = [root]
    attempts = accepted = skipped = 0
    minima = [100.0] * 9
    depth = 0
    while todo:
        box = todo.pop()
        ans = check(box)
        attempts += 1
        if ans is None:
            skipped += 1
        elif min(ans) > 0:
            accepted += 1
            minima = [min(a, b) for a, b in zip(minima, ans)]
        else:
            if attempts >= limit:
                raise RuntimeError(('incomplete', attempts, len(todo), ans))
            i = max(range(4), key=lambda i: box[i][1] - box[i][0])
            l, h = box[i]
            m = (l + h) / 2
            L = list(box)
            R = list(box)
            L[i] = (l, m)
            R[i] = (m, h)
            todo.extend([tuple(R), tuple(L)])
        if attempts % 50000 == 0:
            print('progress', attempts, len(todo), accepted, time.monotonic() - start, flush=True)
    return dict(status='pass', attempts=attempts, accepted=accepted, infeasible=skipped, minimum_first_coefficient=minima[0], minimum_second_coefficient=minima[1], minimum_absorption_margin=minima[2], minimum_coercivity_margins=minima[3:], elapsed_seconds=time.monotonic() - start)


def family_certificate():
    """The whole stated amplitude interval uses the analytic cube theorem."""
    value = scalar_certificate(panels=8192, precision=120, rho_value=Q(3, 5))
    ctx.prec = 120
    rho = arb(3)/5
    d = 1-rho*rho
    rate = 1/(1-rho)
    coefficients = (-arb(21)/34, arb(1), arb(2))
    exponents = (0, 1, 16)
    mass = 1/rate**2
    first = sum((c/(rate+l)**2 for c, l in zip(coefficients, exponents)), arb(0))
    second = sum((c*b/(rate+l+m)**2
                  for c, l in zip(coefficients, exponents)
                  for b, m in zip(coefficients, exponents)), arb(0))
    second_loss = (arb(19)/40)*(2*rho*rho/d**3)*(mass*second-first*first)
    interval = lambda z: arb(z[0]).union(arb(z[1]))
    scalar = interval(value['numerator_minus_twice_denominator'])
    first_loss = (arb(8)/5)*interval(value['closure_correction_enclosure'])
    upper = scalar-first_loss-second_loss
    cube_deviation = (arb(1)/40)*(arb(2)/3).sqrt()*(arb(81)/34)
    assert scalar > 0 and second_loss > 0
    assert upper < -arb(1)/500 and cube_deviation < arb(1)/20
    return dict(rho='3/5', exponents=list(exponents),
                coefficients=['-21/34', '1', '2'], exact_radial_mean='0',
                label_direction='every centered unit vector in R^3',
                amplitude_interval=['-1/40', '1/40'],
                scalar_quadratic_coefficient=enclose(scalar),
                first_mode_loss_lower=enclose(first_loss),
                second_polynomial_loss=enclose(second_loss),
                corrected_quadratic_upper=enclose(upper),
                maximum_coordinate_deviation=enclose(cube_deviation),
                guaranteed_gap='epsilon^2/500', scalar_replay=value)


def controls():
    point = ((Q(1, 3), Q(1, 3)),)*4
    assert min(check(point)) > 0
    assert min(check(point, ratio=Q(2))) < 0
    assert min(check(point, coercivities=(Q(1, 4), Q(3, 40)))) < 0
    try:
        run(limit=1)
    except RuntimeError:
        pass
    else:
        raise AssertionError('an incomplete mass cover must fail')
    return dict(balanced_point=True, excessive_second_weight_rejected=True,
                excessive_coercivity_rejected=True, incomplete_cover_rejected=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    start = time.monotonic()
    cube = run()
    assert cube['attempts'] == 2*(cube['accepted']+cube['infeasible'])-1
    tests = controls()
    print('Complete mass cube', json.dumps(cube), flush=True)
    family = family_certificate()
    print('Finite-amplitude family', json.dumps({k: v for k, v in family.items()
                                               if k != 'scalar_replay'}), flush=True)
    source = Path(__file__).with_name('check_three_cell_coupled_radial.py')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    result = dict(status='pass', exact_mass_cube=['17/60', '23/60'],
                  probability_constraints='p3=1-p1-p2; q3=1-q1-q2; all coordinates in cube',
                  maximum_second_to_first_weight='1/4',
                  gram_coercivities=['1/10', '3/40'], cube_arb_precision_bits=100,
                  scalar_arb_precision_bits=120, cube=cube, family=family,
                  controls=tests, checker_sha256=digest(Path(__file__)),
                  input_source=dict(path=source.name, sha256=digest(source)),
                  environment=dict(python=sys.version.split()[0],
                                   python_flint=__import__('flint').__version__),
                  elapsed_seconds=time.monotonic()-start)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status='pass', elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
