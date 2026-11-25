"""Generate exact prefix equations for the two normalized even quintic charts."""

from argparse import ArgumentParser
from functools import cache
from math import comb, factorial, gcd, lcm
from pathlib import Path

from flint import fmpq, fmpq_mpoly_ctx, fmpq_poly


@cache
def trace(radial, height):
    """Return Lambda((xy)**radial * t**(2*height))."""
    return fmpq(
        2**radial * factorial(radial) * factorial(2 * height)
        * factorial(radial + height),
        factorial(height) * factorial(2 * radial + 2 * height + 1),
    )


def chart(d_zero, *, d_free=False, q_value=None):
    names = ("a4", "a3", "a2", "a1", "b2", "b1", "b0", "c0")
    if q_value is None:
        names += ("q",)
    names += ("a",)
    if not d_zero:
        names += ("b",)
    if d_free:
        names += ("D",)
    ring = fmpq_mpoly_ctx.get(names, ordering="lex")
    v = dict(zip(names, ring.gens()))
    d = v["D"] if d_free else (0 if d_zero else 1)
    b = 1 if d_zero else v["b"]
    q = v["q"] if q_value is None else ring.constant(q_value)
    # Angular winding, half t-degree, and coefficient. The core is centered.
    terms = (
        (1, 1, q), (1, 2, 1),
        (-5, 0, d), (-4, 0, v["a4"]),
        (-3, 0, v["a3"]), (-3, 1, b),
        (-2, 0, v["a2"]), (-2, 1, v["b2"]),
        (-1, 0, v["a1"]), (-1, 1, v["b1"]), (-1, 2, v["a"]),
        (0, 0, -v["b0"] / 3 - v["c0"] / 5),
        (0, 1, v["b0"]), (0, 2, v["c0"]),
    )
    return names, ring, terms[:2] + tuple(term for term in terms[2:] if term[2] != 0)


def generate_moments(ring, terms, max_order):
    factorials = [factorial(i) for i in range(max_order + 1)]
    coefficient_powers = [
        [coefficient**i for i in range(max_order + 1)]
        for _, _, coefficient in terms
    ]
    counts = [0] * len(terms)

    for order in range(1, max_order + 1):
        moment = ring.constant(0)

        def visit(index, remaining, radial, height, denominator):
            # A term of winding -j consumes j+1 slots including the j
            # compensating positive-winding factors. At a leaf choose how
            # many of those positive factors have height 1 or 2.
            if index == len(terms):
                for q_count in range(radial + 1):
                    counts[0] = q_count
                    counts[1] = radial - q_count
                    monomial = ring.constant(1)
                    for j, count in enumerate(counts):
                        if count:
                            monomial *= coefficient_powers[j][count]
                    multiplicity = factorials[order] // (
                        denominator * factorials[q_count] * factorials[radial - q_count]
                    )
                    moment.iadd(
                        monomial * trace(radial, height + 2 * radial - q_count)
                        * multiplicity
                    )
                return
            winding, half_t_degree, _ = terms[index]
            cost = 1 - winding
            possible = (
                range(remaining // cost + 1) if index < len(terms) - 1 else (remaining,)
            )
            for count in possible:
                counts[index] = count
                visit(
                    index + 1, remaining - count * cost, radial - count * winding,
                    height + count * half_t_degree, denominator * factorials[count],
                )

        visit(2, order, 0, 0, 1)
        yield moment


def literal_moments(ring, terms, max_order):
    """Small independent control by literal Laurent multiplication."""
    zero = ring.constant(0)
    laurent = {}
    for winding, height, coefficient in terms:
        y_degree = max(0, -winding)
        for i in range(y_degree + 1):
            key = (winding, height + i)
            scale = fmpq((-1)**i * comb(y_degree, i), 2**y_degree)
            laurent[key] = laurent.get(key, zero) + scale * coefficient
    power = {(0, 0): ring.constant(1)}
    for _ in range(max_order):
        product = {}
        for (winding, height), coefficient in power.items():
            for (shift, extra), factor in laurent.items():
                key = (winding + shift, height + extra)
                product[key] = product.get(key, zero) + coefficient * factor
        power = product
        yield sum(
            (coefficient / (2 * height + 1)
             for (winding, height), coefficient in power.items() if winding == 0),
            zero,
        )


def triangular_pivots():
    q = fmpq_poly([0, 1])
    def integral(j, a):
        return sum(
            comb(j, k) * trace(j, j + a + k) * q**(j - k) for k in range(j + 1)
        )
    def pivot(j, a):
        return (j + 1) * integral(j, a)

    for j in range(1, 4):
        assert pivot(j, 0).gcd(pivot(j, 1)) == 1
    p4 = pivot(4, 0)
    assert p4 == fmpq(8, 21879) * fmpq_poly([
        fmpq(1287, 15295), fmpq(1716, 3059), fmpq(198, 133), fmpq(36, 19), 1,
    ])
    return [pivot(j, 0) for j in range(1, 5)]


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Output directory for the equations")
    parser.add_argument("--max-order", type=int, default=12)
    parser.add_argument("--d-zero", action="store_true", help="Use the normalized B=C=1 chart")
    args = parser.parse_args()
    if args.max_order < 1:
        parser.error("--max-order must be positive")
    args.directory.mkdir(parents=True, exist_ok=True)
    names, ring, terms = chart(args.d_zero)
    print("Variables:", ",".join(names))
    print("Pure a_j pivots:", *(str(pivot) for pivot in triangular_pivots()))
    controls = list(literal_moments(ring, terms, min(args.max_order, 4)))
    for order, moment in enumerate(generate_moments(ring, terms, args.max_order), 1):
        if order <= len(controls):
            assert moment == controls[order - 1]
        if order == 1:
            assert not moment
        coefficients = list(moment.to_dict().values())
        scale = lcm(*(int(c.denominator) for c in coefficients))
        content = gcd(*(int(c * scale) for c in coefficients))
        primitive = moment * (fmpq(scale, content) if content else 1)
        (args.directory / f"M{order:02}.txt").write_text(f"{moment}\n", encoding="ascii")
        (args.directory / f"P{order:02}.txt").write_text(f"{primitive}\n", encoding="ascii")
        print(order, len(moment), scale, content, flush=True)


if __name__ == "__main__":
    main()
