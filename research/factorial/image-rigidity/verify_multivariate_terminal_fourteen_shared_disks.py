#!/usr/bin/env python3
"""Exact grouped Rouché tests for fully shared terminal-fourteen profiles."""
from fractions import Fraction as Q
from sys import argv

from flint import acb, acb_poly, arb, ctx

from verify_multivariate_terminal_fourteen_retained_profiles import configure_shared_model
from verify_multivariate_terminal_thirteen_shared_disks import disk_model, radius_at

BANDS = {64: 128, 65: 128, 66: 128, 67: 135, 68: 139, 69: 134, 70: 131, 71: 144, 72: 224, 73: 293, 74: 336, 75: 386}


def propose_centers(references):
    upper = []
    real = []
    for order, reference in references.items():
        try:
            roots = acb_poly(list(reference)).roots(tol=arb(2)**-45, maxprec=384)
        except (ArithmeticError, ValueError):
            return None
        positive = [z for z in roots if z.imag > 0]
        negative = [z for z in roots if z.imag < 0]
        if len(positive) != order//2 or len(negative) != order//2:
            return None
        upper.extend((order, acb(z.real.mid(), z.imag.mid())) for z in positive)
        real.extend((order, acb(z.real.mid())) for z in roots if z.imag.contains(0))
    if len(real) != 1 or real[0][0] % 2 != 1:
        return None
    upper.sort(key=lambda entry: float((entry[1].real/abs(entry[1])).mid()), reverse=True)
    return upper, real[0]


def check_rectangle(model, means, skew, cap, radii):
    models = {}
    for order in model.ORDERS:
        ceiling = model.ordinary_ratio_upper(cap, means[0]) if order == model.UPPER else None
        raw = model.coefficient_model(
            order, cap, means, skew, True, precise=True, monotone=True,
            coupled=True, ratio_ceiling=ceiling, loss_transport=True,
            center_ratio=True, center_skew=True,
        )
        if raw is None:
            return True
        models[order] = raw
    if not all(models.values()):
        return False

    models = {order: disk_model(raw) for order, raw in models.items()}
    centers = propose_centers({order: data[0] for order, data in models.items()})
    if centers is None:
        return False
    upper, (odd, negative) = centers

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


def check_band(model, lo, hi, cap, radii):
    cap = model.exact_cap(cap, lo)
    pending = [((lo, hi), (Q(0), Q(1)), 0, 0)]
    accepted = 0
    while pending:
        means, signed, mean_depth, skew_depth = pending.pop()
        if check_rectangle(model, means, signed, cap, radii):
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
    if requested - BANDS.keys():
        raise ValueError("Select bands by their integer lower endpoint.")
    ctx.prec = 112
    model = configure_shared_model()
    radii = tuple(arb(k) / 2**p for p in range(20, 4, -1) for k in (4, 5, 6, 7))
    total = 0
    for lo, expected in BANDS.items():
        if requested and lo not in requested:
            continue
        checked = check_band(model, Q(lo, 100), Q(lo+1, 100), Q(1, 60), radii)
        if checked != expected:
            raise ArithmeticError(("unexpected disk count", lo, checked, expected))
        print(lo, lo+1, checked, flush=True)
        total += checked
    print(f"Exact accepted rectangles: {total}.", flush=True)


if __name__ == "__main__":
    main()
