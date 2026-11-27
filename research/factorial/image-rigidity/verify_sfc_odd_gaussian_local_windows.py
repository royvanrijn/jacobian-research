#!/usr/bin/env python3
"""Exact modular checks for the quadratic Gaussian moment gate."""

import argparse
import subprocess
from math import isqrt


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def scaled_terms(order, modulus):
    """Generate the homogeneous terms of order**order times the moment."""
    factorial = [1]
    double_factorial = [1]
    powers = [1]
    for k in range(1, order + 1):
        factorial.append(factorial[-1] * k % modulus)
        double_factorial.append(double_factorial[-1] * (2 * k - 1) % modulus)
        powers.append(powers[-1] * order % modulus)
    inv = [pow(f, -1, modulus) for f in factorial]

    for a in range(order + 1):
        for b in range((order - a) // 2 + 1):
            c = order - a - 2 * b
            coef = (
                factorial[order]
                * inv[a]
                * inv[2 * b]
                * inv[c]
                * double_factorial[order - a - b]
                * powers[a + b]
                % modulus
            )
            yield a, b, c, coef


def chart_form(order, prime, chart):
    terms = []
    for a, b, c, coef in scaled_terms(order, prime):
        if not coef or (chart == "infinity" and c) or (chart == "constant_zero" and a):
            continue
        if chart == "quadratic":
            i, j = a, b
        elif chart == "infinity":
            i, j = a, 0
        elif chart == "constant":
            i, j = b, c
        else:
            i, j = 0, c
        terms.append(f"{coef}*a^{i}*u^{j}")
    return "+".join(terms) or "0"


def hensel_row(order, prime, point):
    modulus = prime * prime
    b_square, quadratic = point
    b_powers = [pow(b_square, k, modulus) for k in range(order + 1)]
    c_powers = [pow(quadratic, k, modulus) for k in range(order + 1)]
    value = derivative_b = derivative_c = 0
    for _, b, c, coefficient in scaled_terms(order, modulus):
        value += coefficient * b_powers[b] * c_powers[c]
        if b:
            derivative_b += b * coefficient * b_powers[b - 1] * c_powers[c]
        if c:
            derivative_c += c * coefficient * b_powers[b] * c_powers[c - 1]
    if value % prime:
        raise RuntimeError("The supplied point is not a residue zero.")
    return derivative_b % prime, derivative_c % prime, (value % modulus) // prime


def check_hensel_rows(orders, prime, point):
    (a, b, c), (d, e, f), (g, h, i) = rows = [
        hensel_row(order, prime, point) for order in orders
    ]
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if determinant % prime == 0:
        raise RuntimeError(f"No Hensel obstruction from the rows {rows}.")
    print(f"Hensel rows at {point}: {rows}, determinant {determinant % prime}.")


def verify(start, prime, intermediate, executable, residue_point=None):
    orders = (start, start + 1, start + 2)
    if intermediate:
        n = start + 1
        if not n + 1 < prime < 2 * n - 1:
            raise ValueError("The intermediate prime lies outside the local theorem.")
        charts = ("constant", "constant_zero")
        term_order = "dp"
    else:
        if prime <= 2 * max(orders):
            raise ValueError("Use a prime greater than twice the largest order.")
        charts = ("quadratic", "infinity")
        term_order = "wp(1,2)"

    blocks = [f"ring R={prime},(a,u),{term_order};", "option(redSB);"]
    for i, chart in enumerate(charts):
        forms = ",".join(chart_form(j, prime, chart) for j in orders)
        blocks.append(f"ideal I{i}={forms};")
        blocks.append(f"ideal H{i}=std(I{i});")
    first_chart = "H0[1]==1" if residue_point is None else "vdim(H0)==1"
    blocks.extend([f"{first_chart} && H1[1]==1;", "quit;"])
    result = subprocess.run(
        [executable, "-q"],
        input="\n".join(blocks),
        text=True,
        capture_output=True,
        check=True,
    )
    if result.stdout.strip() != "1":
        raise RuntimeError(
            f"Window {orders} failed at {prime}: {result.stdout}\n{result.stderr}"
        )
    if residue_point is not None:
        check_hensel_rows(orders, prime, residue_point)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--last", type=int, default=127)
    parser.add_argument("--prime", type=int, default=65537)
    parser.add_argument("--intermediate", action="store_true")
    parser.add_argument("--residue-point", nargs=2, type=int, metavar=("B_SQUARED", "C"))
    parser.add_argument("--singular", default="Singular")
    args = parser.parse_args()
    if args.first < 1 or args.first % 2 == 0 or args.last < args.first:
        parser.error("The start range must begin at a positive odd order.")
    if not is_prime(args.prime):
        parser.error("The requested modulus must be prime.")

    if args.residue_point is not None and (not args.intermediate or args.first != args.last):
        parser.error("A residue point requires an intermediate prime and one start.")

    for start in range(args.first, args.last + 1, 2):
        verify(start, args.prime, args.intermediate, args.singular, args.residue_point)
        print(f"Verified odd start {start} at prime {args.prime}.", flush=True)


if __name__ == "__main__":
    main()
