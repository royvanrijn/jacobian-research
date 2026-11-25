#!/usr/bin/env python3
"""Exact paired-ray separation for terminal-twelve retained profiles."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from flint import acb_poly, arb, ctx

from verify_multivariate_terminal_twelve_shared_light import (
    ORDERS, coefficient_model, exact_cap, independent_exponents,
    nominal_exponent, scaled_native,
)
from verify_multivariate_terminal_nine_retained_rows import power_bounds as powq
from verify_multivariate_terminal_nine_shared_matched_upper import ANNULUS, replay_python
from verify_multivariate_terminal_six_intermediate_cayley import ONE, interval
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_six_private_angles import multiply_interval, require
from verify_multivariate_terminal_six_retained_rows import ball
from verify_multivariate_terminal_eleven_retained_profiles import fixed_powers, verify_grouped


CAP = Q(1, 8)
SHARED_CUTOFF = Q(8)
PRIVATE_CUTOFF = Q(6)
TAIL_EXPONENT = Q(19, 22)
ATTENUATION = 194
RAY_SCALE = 10**8
RAY_POSITIONS = (
    Q(1, 2), Q(2, 5), Q(3, 5), Q(3, 10), Q(7, 10), Q(1, 5),
    Q(4, 5), Q(1, 10), Q(9, 10), Q(1, 16), Q(15, 16),
)
SECTOR_ORDERS = (12, 11) * 5 + (12,)


def cores(private):
    cutoff = PRIVATE_CUTOFF if private else SHARED_CUTOFF
    eligible = []
    gap_bound = int(cutoff * (1 + 1 / (11 * CAP)))
    for gap in range(1, gap_bound + 1):
        base_bound = int((gap / CAP - 1) / 11)
        for base in range(1, base_bound + 1):
            strength = Q(gap * gap, base + gap)
            if Q(gap, 11 * base + 1) > CAP and strength <= cutoff:
                eligible.extend((base, gap, sign, strength) for sign in (-1, 1))

    # Compare both orientations in the same lexicographic order.
    if not private:
        eligible.sort()

    def collect(index, rows, budget):
        if index == len(eligible):
            if private or rows and rows <= tuple(sorted((a, d, -sign) for a, d, sign in rows)):
                yield rows, budget
            return
        base, gap, sign, strength = eligible[index]
        for copies in range(1 + int(budget / strength)):
            yield from collect(index + 1, rows + ((base, gap, sign),) * copies,
                               budget - copies * strength)

    return sorted(collect(0, (), cutoff))


@lru_cache(maxsize=32)
def core_coefficients(rows, private):
    return {order: tuple(powq(value, Q(1, order)) for value in fixed_powers(order, rows, private))
            for order in ORDERS}


def anchor_angles(rows, private, mean):
    upper_roots = []
    for order in ORDERS:
        factors = fixed_powers(order, rows, private)
        coefficients = [arb(1)] + [
            comb(order, k) * ball(factors[k])**ball(Q(1, order))
            * ball(mean)**ball(nominal_exponent(order, k))
            for k in range(1, order)
        ] + [arb(1)]
        try:
            roots = acb_poly(coefficients).roots(tol=arb(2)**-58, maxprec=512)
        except (ArithmeticError, ValueError):
            return None
        upper = [z.real / abs(z) for z in roots if z.imag > 0]
        if len(roots) != order or len(upper) != order // 2 or sum(z.imag < 0 for z in roots) != order // 2:
            return None
        upper_roots.extend((cosine, order) for cosine in upper)
    upper_roots.sort(key=lambda item: float(item[0].mid()), reverse=True)
    if (tuple(order for _, order in upper_roots) != SECTOR_ORDERS
        or any(not left > right for (left, _), (right, _) in zip(upper_roots, upper_roots[1:]))):
        return None
    return tuple(cosine for cosine, _ in upper_roots)


def bounds(core, means, cap):
    tail = {}
    for order in ORDERS:
        if cap:
            interior = [(powq(means[0], hi)[0], powq(means[1], lo)[1])
                        for lo, hi in independent_exponents(order, cap)]
            tail[order] = [(Q(1), Q(1)), *interior, (Q(1), Q(1))]
        else:
            tail[order] = [(Q(1), Q(1))] * (order + 1)
    result = {order: [(comb(order, k) * lo, comb(order, k) * hi) for k, (lo, hi) in enumerate(
                  multiply_interval(c, t) for c, t in zip(core[order], tail[order], strict=True))]
              for order in ORDERS}
    for k in range(1, 12):
        left = tail[11][k][1] * powq(means[1], Q(k * k, ATTENUATION))[1]
        right = tail[11][k - 1][1] * powq(means[1], Q((12-k)**2, ATTENUATION))[1]
        lower, upper = result[12][k]
        result[12][k] = lower, min(upper, comb(12, k) * core[12][k][1] * min(left, right))
    require(all(lo <= hi for order in ORDERS for lo, hi in result[order]), "coefficient bounds")
    return result


def cap_at(lower, declared=CAP):
    cap = min(declared, (lower**-2 - 1) / 10)
    if cap == 0:
        return cap
    # Round outward once to reuse the coefficient models across core lists.
    scaled = cap * 10000
    rounded = min(declared, Q(-(-scaled.numerator // scaled.denominator), 10000))
    return exact_cap(rounded, lower)


def verify_rays(core, coefficients, means, cap, targets, grouped, replay, *,
                model_builder=coefficient_model):
    factors = {order: tuple(interval(*pair) for pair in core[order]) for order in ORDERS}
    boxes = {order: ((ONE, tuple(interval(*pair) for pair in coefficients[order])),) for order in ORDERS}
    count = 0
    for left, right in zip(targets, targets[1:]):
        left_midpoint, right_midpoint = float(left.mid()), float(right.mid())
        proposed = set()
        for position in RAY_POSITIONS:
            c = Q(round(RAY_SCALE * (left_midpoint * float(position) +
                                     right_midpoint * float(1-position))), RAY_SCALE)
            if c in proposed or not left > ball(c) > right:
                continue
            proposed.add(c)
            cells = []
            for order in ORDERS:
                answer = replay(order, (c, c), boxes[order], Q(0))
                if answer is None and grouped and cap:
                    answer = verify_grouped(factors[order], order, means, cap, c, replay, model_builder)
                if answer is None:
                    break
                cells.append(answer)
            if len(cells) == 2:
                count += sum(cells)
                break
        else:
            return None
    return count


def enclosing_annulus(coefficients):
    inner, outer = ANNULUS
    for order in ORDERS:
        box = coefficients[order]
        caps = tuple(hi for _, hi in box)
        dominance = (sum(caps[k]*inner**k for k in range(1, order+1)) < 1
                     and sum(caps[k]*outer**k for k in range(order)) < outer**order)
        ratios = all(a > inner*d and b < outer*c for (a,b), (c,d) in zip(box, box[1:]))
        if not (dominance or ratios):
            return False
    return True


def verify_interval(rows, private, outer, declared, replay, *, model_builder=coefficient_model,
                    grouped_depth=5):
    core = core_coefficients(rows, private)
    global_annulus = enclosing_annulus(bounds(core, outer, cap_at(outer[0], declared)))
    pending = [(outer, 0)]
    leaves = cells = 0
    while pending:
        means, depth = pending.pop()
        cap = cap_at(means[0], declared)
        midpoint = sum(means) / 2
        targets = anchor_angles(rows, private, midpoint)
        count = None
        if targets is not None:
            coefficients = bounds(core, means, cap)
            if global_annulus or enclosing_annulus(coefficients):
                count = verify_rays(core, coefficients, means, cap, targets, depth >= grouped_depth, replay,
                                    model_builder=model_builder)
        if count is not None:
            leaves += 1
            cells += count
            continue
        require(depth < 14, ("unresolved core", rows, private, means, cap, targets))
        pending.extend((((means[0], midpoint), depth + 1), ((midpoint, means[1]), depth + 1)))
    return leaves, cells


def verify_parameters():
    ctx.prec = 192
    _, log_outer_upper = log_bounds(Q(9, 4))
    require(276 * TAIL_EXPONENT * log_outer_upper < ATTENUATION, "factorial tail attenuation")
    require(log_bounds(Q(217,96))[0] / log_outer_upper > 1 + Q(1, ATTENUATION),
            "Gaussian tail attenuation")
    require(Q(1,110) > Q(12,11*ATTENUATION), "Gaussian endpoint exponent")


def verify(private, selected, replay):
    verify_parameters()
    items = cores(private)
    chosen = list(enumerate(items)) if selected is None else [(i, items[i]) for i in selected]
    leaf_total = cell_total = 0
    for index, (rows, budget) in chosen:
        lower = powq(Q(4, 9), TAIL_EXPONENT * budget)[0]
        leaves, cells = verify_interval(rows, private, (lower, Q(1)), CAP, replay)
        leaf_total += leaves
        cell_total += cells
        print(index, rows, leaves, cells, flush=True)
    print(f"Verified {len(chosen)} of {len(items)} lists, {leaf_total} outer cells and "
          f"{cell_total} radial cells.", flush=True)


def main():
    args = [x for x in argv[1:] if x != "--native"]
    require(args and args[0] in ("private", "shared"), "expected profile family")
    private = args[0] == "private"
    selected = tuple(map(int, args[1:])) or None
    if "--native" not in argv:
        verify(private, selected, replay_python)
        return
    with TemporaryDirectory(prefix="sfc-twelve-retained-") as directory:
        library = Path(directory) / "bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(private, selected, scaled_native(CDLL(str(library)), compress=True))


if __name__ == "__main__":
    main()
