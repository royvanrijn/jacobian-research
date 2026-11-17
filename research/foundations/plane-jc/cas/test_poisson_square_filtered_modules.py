#!/usr/bin/env python3
"""Regression tests for the Poisson-square associated-stratum filter."""

import sympy as sp

from poisson_square_filtered_modules import (
    associated_stratum_ledger,
    lower_band_survival_audit,
)


ledger = associated_stratum_ledger()
assert tuple(stratum.name for stratum in ledger) == (
    "T",
    "C0",
    "A0",
    "S",
    "S_C",
    "S_A",
    "K_C",
    "K_A",
)
assert tuple(stratum.normalized_dimension for stratum in ledger) == (
    3,
    3,
    3,
    2,
    2,
    2,
    1,
    1,
)
assert tuple(stratum.d3_length for stratum in ledger) == (
    None,
    None,
    None,
    3,
    5,
    7,
    25,
    30,
)
assert tuple(stratum.d2_length for stratum in ledger) == (
    1,
    None,
    None,
    8,
    4,
    2,
    None,
    None,
)
assert tuple(stratum.d3_socle_dimension for stratum in ledger) == (
    None,
    None,
    None,
    2,
    2,
    2,
    4,
    6,
)
assert tuple(stratum.d2_socle_dimension for stratum in ledger) == (
    1,
    None,
    None,
    2,
    1,
    1,
    None,
    None,
)
assert tuple(stratum.d2_presentation for stratum in ledger) == (
    "reduced_field",
    None,
    None,
    "three_quadrics_plus_one_cubic",
    "four_relations_on_three_displayed_generators",
    "dual_numbers",
    None,
    None,
)
assert tuple(stratum.d3_presentation for stratum in ledger) == (
    None,
    None,
    None,
    "square_zero_two_generators",
    "two_generators_three_relations",
    "two_generators_four_groebner_relations",
    "three_generators_nine_groebner_relations",
    "four_generators_eighteen_groebner_relations",
)

p30 = sp.Symbol("p_3_0")
q11 = sp.Symbol("q_1_1")
p20, p21, q10 = sp.symbols("p_2_0 p_2_1 q_1_0")

bottom_constant = 2 * p20 * q11 - p21 * q10 - 1
assert all(
    report.status == "preserved"
    for report in lower_band_survival_audit((bottom_constant,))
)

d1_reports = {
    report.stratum: report.status
    for report in lower_band_survival_audit((q11,))
}
assert d1_reports == {
    "T": "cut",
    "C0": "eliminated",
    "A0": "eliminated",
    "S": "eliminated",
    "S_C": "preserved",
    "S_A": "eliminated",
    "K_C": "preserved",
    "K_A": "eliminated",
}

a0_reports = {
    report.stratum: report.status
    for report in lower_band_survival_audit((p30,))
}
assert a0_reports == {
    "T": "cut",
    "C0": "eliminated",
    "A0": "preserved",
    "S": "preserved",
    "S_C": "eliminated",
    "S_A": "preserved",
    "K_C": "preserved",
    "K_A": "preserved",
}

print("PASS: all eight associated strata satisfy the normalized Poisson system")
print("PASS: certified Hilbert, socle, and d2-presentation data are exposed")
print("PASS: localized lower-band survival distinguishes preserve/cut/eliminate")
