#!/usr/bin/env python3
"""Verify exclusion certificates for de Smit's three curated abc catalogs.

Run with Python 3.11+ and the directory containing good.html, merit.html,
and unbeaten.html.
Unresolved rows get their reported quality displayed, without treating it
as a certified upper bound when a factor is hidden as BIG.
"""

from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
import re
from sys import set_int_max_str_digits


CATALOGS = {
    "good.html": (
        "https://numberdb.org/files/T62/abc-quality.html",
        "71c1e8057c32f46ead771392c7fd73c0bd10eeec91c33e5f57a809f6dc68b619",
        241,
    ),
    "merit.html": (
        "https://web.archive.org/web/20251211064253id_/https://pub.math.leidenuniv.nl/~smitbde/abc/?set=3",
        "9dceeafbe448482d3ac58d7fcd365cfc72fc0a9c33dcfb9aa7c94a7142d2d1ea",
        202,
    ),
    "unbeaten.html": (
        "https://web.archive.org/web/20251211063930id_/https://pub.math.leidenuniv.nl/~smitbde/abc/?set=4",
        "9e0a973a6f4c5ae47cce3d645ac72da802196000302805462f15c43bf9d8a806",
        160,
    ),
}
TRIAL_LIMIT = 1_000_000
EXTRA_WITNESSES = [2_259_353, 6_971_303, 11_837_159, 28_542_797, 500_120_767]
LARGE = 2**63
FACTOR = re.compile(r"(\d+)(?:<sup>(\d+)</sup>)?")
HIDDEN_FACTOR = re.compile(r"BIG(?:<sup>\d+</sup>)?")


def primes_through(bound):
    sieve = bytearray(b"\1") * (bound + 1)
    sieve[:2] = b"\0\0"
    for prime in range(2, isqrt(bound) + 1):
        if sieve[prime]:
            sieve[prime * prime :: prime] = b"\0" * ((bound - prime * prime) // prime + 1)
    return [prime for prime in range(2, bound + 1) if sieve[prime]]


def certified_prime(value, primes):
    if not 2 <= value <= TRIAL_LIMIT**2:
        return False
    for prime in primes:
        if prime * prime > value:
            return True
        if value % prime == 0:
            return False
    raise ValueError(f"Prime sieve too short for {value}")


def parse_leg(cell):
    cell = cell.replace("&middot;", " ").replace("&#x200B;", " ").replace("BIG <sup>", "BIG<sup>")
    factors = []
    hidden = False
    for token in cell.split():
        if HIDDEN_FACTOR.fullmatch(token):
            hidden = True
            continue
        match = FACTOR.fullmatch(token)
        if match is None:
            raise ValueError(f"Unrecognized factor: {token}")
        prime = int(match[1])
        if prime != 1:
            factors.append((prime, int(match[2] or "1")))
    return factors, hidden


def read_rows(path):
    for line in path.read_text(encoding="latin-1").splitlines():
        if not line.startswith("<tr>") or '<td class="abcnum">' not in line:
            continue
        number, printed_quality = line.split("<td>")[1:3]
        legs = [parse_leg(cell) for cell in line.split('<td class="abcnum">')[-3:]]
        values = [prod(p**e for p, e in factors) for factors, _ in legs]
        hidden_legs = [i for i, (_, hidden) in enumerate(legs) if hidden]
        if len(hidden_legs) > 1:
            raise ValueError(f"More than one hidden leg in {path.name} row {number}")
        if hidden_legs:
            index = hidden_legs[0]
            lower_part = values[index]
            a, b, c = values
            values[index] = (c - b, c - a, a + b)[index]
            if values[index] % lower_part:
                raise ValueError(f"Invalid reconstructed leg in {path.name} row {number}")
        a, b, c = values
        if not (a > 0 and b > 0 and a + b == c and gcd(a, b) == 1):
            raise ValueError(f"Invalid triple in {path.name} row {number}")
        yield int(number), printed_quality, values, legs, hidden_legs


def exclusion(values, legs, hidden_legs, primes, prime_set):
    known_prime_divisors = set()
    for value, (factors, _) in zip(values, legs):
        for prime, exponent in factors:
            if value % prime:
                raise ValueError(f"Reported factor {prime} does not divide its leg")
            if prime in prime_set:
                known_prime_divisors.add(prime)
            if exponent < 4 and value % prime**4 and (
                prime in prime_set or certified_prime(prime, primes)
            ):
                return "visible_valuation", prime

    known_radical = prod(known_prime_divisors)
    if known_radical**4 > prod(values):
        return "radical_lower_bound", None

    for index in hidden_legs:
        value = values[index]
        for prime in EXTRA_WITNESSES:
            if value % prime == 0 and value % prime**4:
                return "additional_valuation", prime
        for prime in primes:
            if prime > value:
                break
            if value % prime == 0 and value % prime**4:
                return "trial_valuation", prime
    return "unresolved", None


def check_catalogs(directory):
    primes = primes_through(TRIAL_LIMIT)
    prime_set = set(primes)
    if not all(certified_prime(prime, primes) for prime in EXTRA_WITNESSES):
        raise ValueError("Invalid additional witness")
    totals = set()
    exact = set()
    unresolved = []
    for name, (_, expected_hash, expected_rows) in CATALOGS.items():
        path = directory / name
        actual_hash = sha256(path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            raise ValueError(f"{name}: unexpected SHA256 {actual_hash}")
        counts = Counter()
        keys = set()
        for row, quality, values, legs, hidden_legs in read_rows(path):
            a, b, c = values
            key = (min(a, b), max(a, b), c)
            if key in keys:
                raise ValueError(f"{name}: duplicate row {row}")
            keys.add(key)
            counts["rows"] += 1
            counts["c_at_least_2^63"] += c >= LARGE
            counts["reconstructed"] += bool(hidden_legs)
            result, witness = exclusion(values, legs, hidden_legs, primes, prime_set)
            counts[result] += 1
            if result == "unresolved":
                unresolved.append((key, name, row, quality))
            else:
                exact.add(key)
                if result == "additional_valuation":
                    print(f"{name} row {row}: verified extra prime witness {witness}")
        if counts["rows"] != expected_rows:
            raise ValueError(f"{name}: expected {expected_rows} rows")
        totals.update(keys)
        print(f"{name}: {dict(counts)}")

    for key, name, row, quality in unresolved:
        if key not in exact:
            print(f"{name} row {row}: no certificate (reported q={quality})")
    print(
        f"Unique triples: {len(totals)}. "
        f"Large triples: {sum(c >= LARGE for _, _, c in totals)}. "
        f"Excluded exactly: {len(exact)}. "
        f"Excluded large triples exactly: {sum(c >= LARGE for _, _, c in exact)}. "
        f"Unresolved: {len(totals - exact)}."
    )


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path.cwd())
    args = parser.parse_args()
    set_int_max_str_digits(0)
    check_catalogs(args.directory)


if __name__ == "__main__":
    main()

