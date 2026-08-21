#!/usr/bin/env sage
"""Reconstruct the slope-8/87 q80 surface series from split-prime jets.

Each input is produced by ``extend_q80_cm24_branch_modp.sage`` from the same
characteristic-zero marked system.  The four surface coordinates are rational,
so their coefficients can be reconstructed independently of the chosen split
CM embedding.  Every reconstructed fraction is reduced back at every input
prime before the exact series artifact is written.
"""

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import CRT_list, QQ, ZZ, prod


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input-glob",
    action="append",
    required=True,
    help="glob for q80 split-prime formal-jet JSON files; repeat as needed",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-order28-series.json",
)
arguments = parser.parse_args()

paths = sorted(
    {
        Path(filename)
        for pattern in arguments.input_glob
        for filename in glob.glob(pattern)
    }
)
if not paths:
    raise ValueError("no input jets matched")

records = []
for path in paths:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "q80-cm24-split-prime-formal-branch-v1":
        raise ValueError(f"unexpected schema in {path}")
    if payload.get("kind") != "normalized_formal_jet":
        raise ValueError(f"expected a formal jet in {path}")
    if payload.get("characteristic_zero_slope") != "8/87":
        raise ValueError(f"wrong branch in {path}")
    records.append((path, payload))

orders = {payload["order"] for _, payload in records}
active_names_set = {tuple(payload["active_variables"]) for _, payload in records}
if len(orders) != 1 or len(active_names_set) != 1:
    raise ValueError("input jets do not have a common order/marking")
order = orders.pop()
active_names = active_names_set.pop()
surface_names = ("D", "P", "Q", "E")
surface_columns = tuple(active_names.index(name) for name in surface_names)
primes = tuple(ZZ(payload["prime"]) for _, payload in records)
if len(set(primes)) != len(primes):
    raise ValueError("duplicate primes")
modulus = prod(primes)


def reconstruct(values):
    residue = ZZ(CRT_list(list(map(ZZ, values)), list(primes)))
    value = QQ(residue.rational_reconstruction(modulus))
    for prime, expected in zip(primes, values):
        reduced = ZZ(value.numerator())*ZZ(value.denominator()).inverse_mod(prime)
        if reduced % prime != ZZ(expected) % prime:
            raise AssertionError((prime, expected, value))
    return value


series = {name: [] for name in surface_names}
for degree in range(order):
    for name, column in zip(surface_names, surface_columns):
        values = [payload["coefficients"][degree][column] for _, payload in records]
        try:
            series[name].append(reconstruct(values))
        except ArithmeticError as error:
            raise ArithmeticError(
                f"rational reconstruction failed for {name}[{degree}] "
                f"with a {modulus.nbits()}-bit modulus"
            ) from error

expected_centers = {
    "D": QQ(-1)/2,
    "P": QQ(9)/4,
    "Q": -QQ(9)/4,
    "E": -QQ(27)/32,
}
assert {name: values[0] for name, values in series.items()} == expected_centers
assert series["P"][1] == 1 and not any(series["P"][2:])

input_metadata = []
for path, payload in records:
    input_metadata.append(
        {
            "path": str(path),
            "prime": payload["prime"],
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    )

output = Path(arguments.output)
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(
    json.dumps(
        {
            "schema": "q80-cm24-qq-surface-series-v1",
            "status": "exact_crt_reconstruction_verified_at_all_input_primes",
            "scope": "exact_formal_series_not_global_algebraization",
            "slope": "8/87",
            "normalization": "P=9/4+h",
            "order": order,
            "crt_modulus_bits": modulus.nbits(),
            "inputs": input_metadata,
            "series": {
                name: [str(value) for value in values]
                for name, values in series.items()
            },
        },
        indent=2,
        sort_keys=True,
    )+"\n"
)

largest_bits = max(
    max(abs(value.numerator()).nbits(), value.denominator().nbits())
    for values in series.values()
    for value in values
)
print(
    "Q80SURFACECRT|slope=8/87|"
    f"primes={len(primes)}|modulus_bits={modulus.nbits()}|order={order}|"
    f"largest_numden_bits={largest_bits}|output={output}|"
    "status=PASS_EXACT_RATIONAL_RECONSTRUCTION",
    flush=True,
)
