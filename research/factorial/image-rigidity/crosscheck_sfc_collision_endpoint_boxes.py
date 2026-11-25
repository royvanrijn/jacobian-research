"""Cross-check the modular collision sieve by direct projective Macaulay rank."""

import argparse
import math
import pathlib
import random
import subprocess

from flint import nmod_mat


def compositions(degree):
    for i in range(degree + 1):
        for j in range(degree - i + 1):
            yield (i, j, degree - i - j)


def macaulay_rank(forms, prime):
    degree = sum(forms) - 2
    index = {power: i for i, power in enumerate(compositions(degree))}
    rows = []
    for order, coefficients in forms.items():
        for shift in compositions(degree - order):
            row = [0] * len(index)
            for exponent, coefficient in coefficients.items():
                combined = tuple(a + b for a, b in zip(exponent, shift))
                row[index[combined]] = coefficient
            rows.append(row)
    return nmod_mat(rows, prime).rank(), len(index)


def parse(line):
    size, *columns = map(int, line.split())
    assert len(columns) == 3 * size
    return tuple(tuple(columns[i:i + 3]) for i in range(0, 3 * size, 3))


def factorial_forms(support, prime):
    result = {}
    for order in (2, 3, 4):
        form = {}
        for exponents in compositions(order):
            multinomial = math.factorial(order) // math.prod(map(math.factorial, exponents))
            moment = math.prod(
                math.factorial(sum(a * b for a, b in zip(exponents, row)))
                for row in support
            )
            form[exponents] = multinomial * moment % prime
        result[order] = form
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("first")
    parser.add_argument("second")
    parser.add_argument("survivors", nargs="+")
    args = parser.parse_args()
    endpoints = [
        ((5, 2, 3), (2, 5, 3), (5, 4, 4), (4, 5, 4), (5, 5, 4)),
        ((7, 7, 6), (5, 4, 4), (4, 5, 4), (5, 5, 4)),
    ]
    random_generator = random.Random(261125)
    controls = endpoints + [
        ((3, 3, 4), (1, 1, 7)),
        ((2, 2, 2), (1, 3, 1)),
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ]
    controls.extend(
        tuple(tuple(random_generator.randrange(10) for _ in range(3))
              for _ in range(random_generator.randrange(1, 8)))
        for _ in range(240)
    )
    survivor_lines = [line for path in args.survivors
                      for line in pathlib.Path(path).read_text().splitlines()]
    survivors = [parse(line) for line in survivor_lines]
    cases = controls + survivors
    text = "\n".join(" ".join(map(str, [len(case), *sum(case, ())])) for case in cases) + "\n"
    for prime, binary in [(65521, args.first), (65519, args.second)]:
        run = subprocess.run([binary, "filter"], check=True, input=text,
                             capture_output=True, text=True)
        flagged = {parse(line) for line in run.stdout.splitlines()}
        ranks = {case: macaulay_rank(factorial_forms(case, prime), prime) for case in cases}
        expected = {case for case, (rank, target) in ranks.items() if rank < target}
        assert expected == flagged, (prime, flagged ^ expected)
        assert all(case not in flagged for case in endpoints)
        assert all(case in flagged for case in controls[2:4])
        assert (all(case in flagged for case in survivors) if prime == 65521
                else all(case not in flagged for case in survivors))
        print(f"prime={prime}: checked={len(cases)}, rank deficient={len(flagged)}")


if __name__ == "__main__":
    main()

