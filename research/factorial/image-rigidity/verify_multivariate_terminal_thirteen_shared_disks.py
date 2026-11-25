#!/usr/bin/env python3
"""Reproduce the grouped Rouché disk certificate for terminal thirteen."""
from fractions import Fraction as Q
from functools import lru_cache
from sys import argv

from flint import acb_poly, arb, ctx

from verify_multivariate_terminal_thirteen_retained_disks import choose_centers
from verify_multivariate_terminal_thirteen_shared_light import (
    GRID, ORDERS, coefficient_model, exact_cap, ordinary_ratio_upper,
)

# Integers p and q denote the closed interval [p/100,q/100].
BANDS = (
    (65, 66, Q(1, 40), 249),
    (69, 70, Q(1, 50), 128),
    (70, 71, Q(1, 50), 128),
    (72, 73, Q(1, 50), 155),
    (73, 74, Q(1, 50), 256),
    (74, 75, Q(1, 50), 248),
    (75, 76, Q(1, 50), 304),
)


@lru_cache(maxsize=4096)
def outward(lo, hi):
    return arb(lo + hi, hi - lo) / (2 * GRID)


def midpoint(bounds):
    return arb(bounds.lo + bounds.hi) / (2 * GRID)


def disk_model(model):
    groups = tuple(
        (outward(weight.lo, weight.hi),
         acb_poly([outward(entry.lo, entry.hi) for entry in vector]))
        for weight, vector in model if weight.lo or weight.hi
    )
    reference = tuple(
        sum(midpoint(weight) * midpoint(vector[k]) for weight, vector in model)
        for k in range(len(model[0][1]))
    )
    return reference, groups


def radius_at(center, groups, order, radii):
    shift = acb_poly([center, 1])
    translated = sum((weight * polynomial(shift) for weight, polynomial in groups),
                     acb_poly())
    upper = tuple(abs(translated[k]).upper() for k in range(order + 1))
    derivative_lower = abs(translated[1]).lower()
    for radius in radii:
        if derivative_lower * radius > upper[0] + sum(
            radius**k * upper[k] for k in range(2, order + 1)
        ):
            return radius
    return None


def check_pair(means, signed, cap, radii):
    models = {}
    for order in ORDERS:
        ceiling = ordinary_ratio_upper(cap, means[0]) if order == ORDERS[1] else None
        raw = coefficient_model(
            order, cap, means, signed, True, precise=True, monotone=True,
            coupled=True, ratio_ceiling=ceiling, loss_transport=True,
            center_ratio=True, center_skew=True,
        )
        if raw is None:
            # At least one bound makes this common parameter cell empty.
            return True
        models[order] = raw
    if not all(models.values()):
        return False

    models = {order: disk_model(raw) for order, raw in models.items()}
    proposed = choose_centers({order: data[0] for order, data in models.items()})
    if proposed is None:
        return False
    upper, (odd, negative) = proposed

    disks = []
    for order, center in upper:
        radius = radius_at(center, models[order][1], order, radii)
        if radius is None or not center.imag > radius:
            return False
        disks.append((center, radius))
    for (a, r), (b, s) in zip(disks, disks[1:]):
        if not (a.conjugate() * b).imag > abs(a) * s + abs(b) * r + r * s:
            return False

    radius = radius_at(negative, models[odd][1], odd, radii)
    return radius is not None and negative.real < -radius


def check_band(lo, hi, cap, radii):
    cap = exact_cap(cap, lo)
    pending = [((lo, hi), (Q(0), Q(1)), 0, 0)]
    accepted = 0
    while pending:
        means, signed, mean_depth, skew_depth = pending.pop()
        if check_pair(means, signed, cap, radii):
            accepted += 1
            continue
        if mean_depth + skew_depth >= 20:
            raise ArithmeticError(("unresolved rectangle", means, signed, cap))
        if mean_depth < skew_depth + 2:
            mid = sum(means) / 2
            pending.extend((
                ((means[0], mid), signed, mean_depth + 1, skew_depth),
                ((mid, means[1]), signed, mean_depth + 1, skew_depth),
            ))
        else:
            mid = sum(signed) / 2
            pending.extend((
                (means, (signed[0], mid), mean_depth, skew_depth + 1),
                (means, (mid, signed[1]), mean_depth, skew_depth + 1),
            ))
    return accepted


def main():
    requested = {int(value) for value in argv[1:]}
    if requested - {band[0] for band in BANDS}:
        raise ValueError("Select bands by their integer lower endpoint.")
    ctx.prec = 112
    radii = tuple(arb(k) / 2**p for p in range(20, 4, -1) for k in (4, 5, 6, 7))
    total = 0
    for begin, end, cap, expected in BANDS:
        if requested and begin not in requested:
            continue
        checked = check_band(Q(begin, 100), Q(end, 100), cap, radii)
        if checked != expected:
            raise ArithmeticError(("unexpected disk count", begin, checked, expected))
        print(begin, end, cap, checked, flush=True)
        total += checked
    print(f"Exact accepted rectangles: {total}.", flush=True)


if __name__ == "__main__":
    main()
