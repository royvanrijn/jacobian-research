"""Enclose the GH6 and GH10 criterion on a rational correlation cover."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from flint import arb, ctx


ORIGINAL_SHA256 = "827b923ff66c2e976527b6489db79d3a1008db87d9ed8e4832f0595671e9f4cb"


def load_original():
    script = Path(__file__).resolve().with_name("check_three_cell_global_scalar_through_half.py")
    content = script.read_bytes()
    assert sha256(content).hexdigest() == ORIGINAL_SHA256
    namespace = {"__file__": str(script), "__name__": "_gh10_original"}
    exec(compile(content, str(script), "exec"), namespace)
    return namespace


def enclose(lower, upper, original, angular, panels, row_order, radius):
    rho = (arb(lower.numerator) / lower.denominator).union(arb(upper.numerator) / upper.denominator)
    finite, row, d, kappa, rate = original["finite_bound"](rho, angular, panels, row_order, radius)
    first = len(angular[0])
    small_rate = 4 * rate
    small_ratio = small_rate / (first + 1)
    eta = 1 + 7 * kappa / 4
    large_ratio = rate / eta
    assert 0 < small_ratio < 1 and 0 < large_ratio < 1
    small_tail = 16 * small_rate**first / (row[0] * factorial(first) * (1 - small_ratio))
    large_tail = 288 * large_ratio**first / (7 * d * eta * (1 - large_ratio))
    return finite, small_tail + large_tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=Fraction, nargs="+")
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--denominator", type=int)
    parser.add_argument("--cap", type=Fraction, default=Fraction(2))
    parser.add_argument("--panels", type=int, default=16384)
    parser.add_argument("--modes", type=int, default=80)
    parser.add_argument("--rows", type=int, default=16)
    parser.add_argument("--radius", type=int, default=12)
    args = parser.parse_args()
    ctx.prec = 100
    original = load_original()
    angular = original["angular_moments"](args.modes)
    upper_float = original["upper_float"]

    settings = dict(checker_sha256=ORIGINAL_SHA256, precision_bits=ctx.prec,
                    panels=args.panels, modes=args.modes, rows=args.rows, radius=args.radius)
    if args.points:
        for point in args.points:
            finite, tail = enclose(point, point, original, angular, args.panels, args.rows, args.radius)
            print(json.dumps(dict(settings, rho=str(point), finite_upper=upper_float(finite),
                                  tail_upper=upper_float(tail), total_upper=upper_float(finite + tail))), flush=True)
        return

    assert args.start < args.stop
    assert Fraction(1, 2) <= Fraction(args.start, args.denominator)
    maximum_upper = 0.0
    maximum_tail = 0.0
    for index in range(args.start, args.stop):
        lower = Fraction(index, args.denominator)
        upper = Fraction(index + 1, args.denominator)
        finite, tail = enclose(lower, upper, original, angular, args.panels, args.rows, args.radius)
        total = finite + tail
        assert total < arb(args.cap.numerator) / args.cap.denominator, (str(lower), str(upper), total)
        maximum_upper = max(maximum_upper, upper_float(total))
        maximum_tail = max(maximum_tail, upper_float(tail))

    print(json.dumps(dict(settings, status="pass", cap=str(args.cap),
                          rho_bounds=[str(Fraction(i, args.denominator)) for i in (args.start, args.stop)],
                          correlation_slabs=args.stop - args.start,
                          slab_width=str(Fraction(1, args.denominator)), tail_upper=maximum_tail,
                          total_upper=maximum_upper), indent=2))


if __name__ == "__main__":
    main()
