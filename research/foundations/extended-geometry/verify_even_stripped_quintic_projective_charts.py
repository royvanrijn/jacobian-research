"""Replay two weighted-projective exclusions in the even B=1 chart."""

from argparse import ArgumentParser
from math import isqrt
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory

from flint import fmpz_mod_mpoly_ctx, nmod_poly

from compute_even_stripped_quintic_moments import chart, generate_moments, literal_moments


QUARTIC = (1287, 8580, 22770, 28980, 15295)
CASES = (
    (65521, 0, 13, {"b2": 8, "b1": 10, "b0": 6, "c0": 9, "a": 10, "b": 6}),
    (101, 69, 14, {"a4": 5, "b2": 8, "b1": 8, "b0": 5, "c0": 8, "a": 9, "b": 6}),
)


def equations(prime, q, max_order):
    names, rational_ring, terms = chart(False, d_free=True, q_value=q)
    rational = list(generate_moments(rational_ring, terms, max_order))
    assert rational[:4] == list(literal_moments(rational_ring, terms, 4))
    assert not rational[0]
    ring = fmpz_mod_mpoly_ctx.get(names, prime, ordering="lex")
    moments = []
    for moment in rational:
        coefficients = {}
        for powers, coefficient in moment.to_dict().items():
            assert int(coefficient.denominator) % prime != 0
            coefficients[powers] = (
                int(coefficient.numerator) * pow(int(coefficient.denominator), -1, prime)
                % prime
            )
        moments.append(ring.from_dict(coefficients))

    eliminated = ("a1", "a2", "a3") if q else ("a1", "a2", "a3", "a4")
    for index, name in enumerate(eliminated, 1):
        variable = names.index(name)
        current = moments[index]
        pivot = current.derivative(variable)
        assert pivot and pivot.is_constant()
        # Each eliminated variable has positive weight. The first equation
        # containing it is affine with a nonzero pivot in this residue field.
        substitute = list(ring.gens())
        substitute[variable] = -current.subs({name: 0}) / pivot.coefficient(0)
        assert not current.compose(*substitute)
        for order in range(index + 1, len(moments)):
            moments[order] = moments[order].compose(*substitute)

    remaining = tuple(name for name in names if name not in eliminated)
    positions = tuple(names.index(name) for name in remaining)
    weights = ((5,) if q else ()) + (3, 2, 1, 1, 2, 4, 6)
    if q:
        assert not moments[4].derivative(names.index("a4"))
    reduced = []
    start = len(eliminated) + 1
    for order, moment in enumerate(moments[start:], start + 1):
        coefficients = {}
        for exponent, coefficient in moment.to_dict().items():
            powers = tuple(exponent[i] for i in positions)
            assert sum(exponent) == sum(powers)
            assert sum(power * weight for power, weight in zip(powers, weights)) == order
            coefficients[powers] = coefficient
        reduced.append(coefficients)
    return remaining[:-1], reduced


def check_special_fiber(prime, q, max_order, expected_powers, solver, threads):
    assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    names, reduced = equations(prime, q, max_order)
    ring = fmpz_mod_mpoly_ctx.get(names, prime, ordering="lex")
    with TemporaryDirectory() as directory:
        for d_value in (0, 1):
            polys = []
            for moment in reduced:
                data = {
                    exponent[:-1]: coefficient
                    for exponent, coefficient in moment.items()
                    if d_value or not exponent[-1]
                }
                polys.append(str(ring.from_dict(data)))
            source = Path(directory) / f"d{d_value}.ms"
            output = Path(directory) / f"d{d_value}.out"
            source.write_text(
                ",".join(names) + f"\n{prime}\n" + ",\n".join(polys) + "\n",
                encoding="ascii",
            )
            subprocess.run(
                [solver, "-f", str(source), "-o", str(output), "-t", str(threads),
                 "-l", "2", "-g", "0" if d_value else "1"],
                check=True,
            )
            result = output.read_text(encoding="ascii").strip()
            if d_value:
                assert result == "[-1]:", "The D=1 special-fiber chart is nonempty"
            else:
                # The exact Gröbner leading ideal contains all six powers.
                # A weighted homogeneous ideal with this leading ideal has
                # only the origin as a geometric affine zero.
                powers = {
                    name: int(power)
                    for name, power in re.findall(
                        r"^\[?([a-z][a-z0-9]*)\^([0-9]+)[,\]]", result, re.MULTILINE
                    )
                }
                assert powers == expected_powers, powers
    print(f"Orders 1–{max_order}, q={q} mod {prime}. D=1 unit ideal. D=0 pure powers {powers}.")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--msolve", default="msolve")
    parser.add_argument("--threads", type=int, default=3)
    args = parser.parse_args()
    x = nmod_poly([0, 1], 31)
    ell31 = nmod_poly(list(QUARTIC), 31)
    assert ell31.degree() == 4 and ell31.gcd(x**(31**2) - x) == 1
    ell101 = nmod_poly(list(QUARTIC), 101)
    assert ell101(69) == 0 and ell101.derivative()(69) != 0
    for case in CASES:
        check_special_fiber(*case, args.msolve, args.threads)

