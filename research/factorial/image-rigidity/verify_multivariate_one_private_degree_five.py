"""Exact scalar estimates for the private-degree-five phase comparison."""

from fractions import Fraction as Q
from math import comb, factorial


def round_up(x, scale=10**12):
    return Q(-(-x * scale // 1), scale)


def exp_partial(x):
    assert 0 <= x < 32
    return sum(x**k / factorial(k) for k in range(33))


def exp_lower(x):
    return -round_up(-exp_partial(Q(x)))


def exp_upper(x):
    x = Q(x)
    return round_up(exp_partial(x) + x**33 / (factorial(33) * (1 - x / 34)))


def log_upper(x):
    def series(y):
        t = (y - 1) / (y + 1)
        return 2 * (sum(t**(2*k+1) / (2*k+1) for k in range(12))
                    + t**25 / (25 * (1 - t*t)))

    x = Q(x)
    power = 0
    while x > 2:
        x /= 2
        power += 1
    return round_up(power * series(Q(2)) + series(x))


def q_upper(n):
    loss = 5 * (1 + log_upper(n) / 2) / n
    return round_up(exp_upper(loss) / exp_lower(Q(5)))


def geometric_moments(ratio, start):
    mass = ratio**start / (1 - ratio)
    offset = ratio / (1 - ratio)
    return (mass, mass * (start + offset),
            mass * ((start + offset)**2 + offset*(1 + offset)))


def factorial_mean_lower(m):
    scale = 10000
    target = factorial(m) * scale**m
    lo, hi = scale, m * scale
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**m <= target:
            lo = mid
        else:
            hi = mid
    return Q(lo, scale)


def private_envelope(n, radius, start, positive_start):
    moments = [Q(0)] * 3
    for m in (n, n+1):
        mean = factorial_mean_lower(m)
        terms = []
        for k in range(start, m//2 + 1):
            index = m-k if positive_start else k
            terms.append(comb(m, k) * (factorial(index) / mean**index)**5 * radius**k)
        # Round outward before adding terms with different denominators.
        for j in range(3):
            bound = sum((start+i)**j * value for i, value in enumerate(terms))
            moments[j] += round_up(bound, 10**10)
    return moments


def uniform_light(n, radius, ratio):
    head = radius * exp_upper(5*Q(1)) * Q(9, 4)**4 / n**4
    half = Q(n, 2)
    exponent = 1 + 4 * (1 + log_upper(half)/2) / half
    assert radius * exp_upper(exponent) / 16 < ratio
    tails = geometric_moments(ratio, 5)
    return tuple(sum(k**j * head**k for k in range(1, 5)) + tails[j]
                 for j in range(3))


def check_direct(n, cutoff, first, intervals, tail, light):
    radius = 1 + Q(1, n)
    x = radius * first
    change = x / n
    individual = list(tail)
    adjacent = list(tail[:2])
    for k in range(2, cutoff+1):
        eta = n - k - Q(1, 2)
        residual_bounds, derivative_bounds = [], []
        for upper_first, upper_curvature in intervals:
            f = upper_curvature*k*(k-1) / (2*eta)
            exp_f = exp_upper(f)
            term = (radius*upper_first)**k / factorial(k)
            residual_bounds.append(term * (exp_f - 1))
            derivative_bounds.append(term * (Q(k, n)*(exp_f-1) + exp_f*f/eta))
        for j in range(3):
            individual[j] += k**j * max(residual_bounds)
            if j < 2:
                adjacent[j] += k**j * max(derivative_bounds)

    for j in range(3):
        first_tail = x**(cutoff+1) * (cutoff+1)**j / factorial(cutoff+1)
        ratio = x / (cutoff+2) * Q(cutoff+2, cutoff+1)**j
        assert ratio < 1
        individual[j] += first_tail / (1-ratio)
        if j < 2:
            adjacent[j] += first_tail / (1-ratio)

    e_x = exp_upper(x)
    j0, j1, j2 = individual
    a0, a1 = adjacent
    e0 = e_x*j0
    e1 = e_x*(j1+x*j0)
    e2 = e_x*(j2+2*x*j1+(x+x*x)*j0)
    diff0 = e_x*(a0+change*j0)
    diff1 = e_x*(a1+x*a0+change*(j1+(1+x)*j0))
    l0, l1, l2 = light
    assert max(e0, l0) < 1
    h, s = 1-e0, 1-l0
    mass = x + e0/h + l0/s
    first_euler = x + e1/h + l1/s
    second_euler = x + e2/h + (e1/h)**2 + l2/s + (l1/s)**2
    log_difference = change + diff0/h + l0/s
    euler_difference = change + diff1/h + e1*diff0/h**2 + l1/s + l1*l0/s**2
    assert mass < 1-Q(1, 2*n), n
    assert first_euler < n, n
    radius_drift = (Q(1, n)+log_difference) / (1-first_euler/n)
    assert euler_difference + second_euler*radius_drift/n < Q(1, 2), n


def main():
    assert (5*(1+log_upper(14)/2) + Q(28, 13) + Q(41, 6)) < 28
    for n in range(14, 100):
        radius = 1 + Q(1, n)
        cap = (n+1) * q_upper(n)
        if n < 70:
            assert cap < Q(3, 5)
        first = min(Q(3, 5), cap)
        intervals = [(min(first, 4*cap/lo), hi) for lo, hi in ((4, 6), (6, 9), (9, 14))]
        tail = private_envelope(n, radius, 5, positive_start=True)
        light = private_envelope(n, radius, 1, positive_start=False)
        check_direct(n, 4, first, intervals, tail, light)

    n = 100
    radius = Q(101, 100)
    ratio = Q(9, 25)
    assert radius * Q(3, 40) * exp_upper(1 + 14*Q(7, 183)) < ratio
    assert radius * exp_upper(4 * (1 + log_upper(50)/2) / 50) / 8 < ratio
    check_direct(n, 8, Q(3, 5), [(Q(3, 5), 14)], geometric_moments(ratio, 9),
                 uniform_light(n, radius, Q(11, 50)))

    q = Q(43, 5000)
    assert Q(71, 70)*q_upper(70) < q
    eta = Q(777, 1000)
    assert 1-25*q-Q(1, 140) > eta
    assert 25*q+Q(1, 70) < Q(1, 2)
    h_times_b = 4*q*Q(28, 27)
    h = h_times_b / Q(3, 5)
    extra = exp_upper(1/(70*eta))
    u, beta_times_b = Q(283, 1000), Q(43, 1000)
    assert 25*q*extra/eta < u
    assert extra/(4*eta) < Q(1, 3)
    assert h_times_b*exp_upper(h+Q(1, 8)) < beta_times_b
    p = Q(7, 20)
    assert (Q(256, 81)*beta_times_b/(1-u)
            * exp_upper(25*beta_times_b*(Q(1, 2)+u/(6*(1-u))))) < p
    assert 1+8/(exp_lower(Q(1))*eta) < 5
    assert h_times_b*exp_upper(h) < Q(19, 500)

    for radius, cap in ((Q(6, 5), Q(7, 5)), (Q(21, 10), Q(37, 10))):
        norm = radius + Q(19, 500)*radius**2/2 + p**2*radius**3/(3*(1-p*radius))
        assert norm < cap

    assert exp_upper(1 + 25*4*q/(2*eta))/10 < Q(19, 40)
    assert exp_upper(4 * (1 + log_upper(35)/2) / 35)/8 < Q(19, 40)
    assert exp_upper(Q(11, 20)) < Q(7, 4)
    ratio = Q(57, 100)
    assert ratio*exp_upper(Q(7, 50)) < 1
    residual = (1/exp_lower(Q(519, 100))
                + ratio**15/(1-ratio)*exp_upper(Q(21, 25)))
    assert residual/(1-residual) < Q(1, 125)

    light_mass, _, _ = uniform_light(70, Q(6, 5), Q(29, 100))
    assert light_mass/(1-light_mass) < Q(1, 250)
    assert Q(7, 5)*q + (Q(1, 125)+Q(1, 250))/70 < Q(1, 80)
    adjacent = Q(23, 200)
    assert 5*q*Q(6, 5)/(1-p*Q(6, 5)) + 2*(Q(1, 125)+Q(1, 250)) < adjacent
    assert exp_upper(Q(13, 100)) < Q(8, 7)
    assert 1/(exp_lower(Q(1))*Q(13, 100)) < 3
    assert 4/(exp_lower(Q(1))*Q(13, 100))**2 < 33
    assert exp_upper(Q(1, 80)) < Q(21, 20)
    drift = (Q(1, 20)+adjacent)/(1-Q(3, 80))
    assert 3*adjacent + Q(33, 80)*drift < Q(1, 2)
    print("Private degree five coefficient, first-factor and phase bounds passed.")


if __name__ == "__main__":
    main()
