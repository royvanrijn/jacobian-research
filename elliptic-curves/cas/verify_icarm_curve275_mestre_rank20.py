#!/usr/bin/env python3
"""Replay ICARM #275 as a generalized Fermigier rank-20 specialization.

The public model is recovered from ``(u,v)=(-3,-1/2)``.  In the primitive
integer-root chart its roots are ``(0,113,550,753,868,1058)`` and its
parameter is ``3413/11`` (the public ``10239/176`` is the native parameter).
The 54 pinned rational quartic abscissas below came from the declared
height-2,000,000, denominator-13,000 ``ratpoints`` box.

This verifier reconstructs the ordinates, injects Fermigier's thirteenth
generic section, maps everything to the short Jacobian, and gives an exact
mod-3 finite-reduction rank-20 certificate.  The selected basis consists of
twelve generic directions and eight exceptional directions.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from math import isqrt

from elliptic_candidate_record import (
    WeierstrassChange,
    change_weierstrass_model,
)
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from search_six_root_low_conductor_centers import pool_with_sources


Q = Fraction
ROOTS = tuple(map(Q, (0, 113, 550, 753, 868, 1058)))
PARAMETER = Q(3413, 11)
NATIVE_PARAMETER = Q(10239, 176)
CONSTRUCTION = SixRootMestreConstruction(ROOTS)
PUBLIC_MODEL = (
    Q(1), Q(0), Q(1), Q(-2034488389107661074627844285),
    Q(35847670110541831966937994064437784692732),
)
SHORT_TO_PUBLIC = (Q(18, 121), Q(27, 14641), Q(9, 121), Q(2916, 1771561))
CONDUCTOR = 42943483208607336815574462222443765285847682460232958112909965535488718

RATPOINTS_ABSCISSAS = tuple(map(Q, (
    "-94898/11", "-3413/11", "-2170/11", "1743/11", "2637/11",
    "3413/11", "4656/11", "4870/11", "6135/11", "8225/11",
    "9463/11", "11696/11", "12961/11", "15051/11", "79501/73",
    "92446/307", "245933/341", "136343/349", "-344923/451",
    "182760/451", "316171/451", "394366/451", "400168/451",
    "445865/451", "-250232/583", "-149269/583", "221537/605",
    "454737/605", "-267109/781", "379655/803", "335689/935",
    "823249/935", "1471257/979", "410241/1067", "366081/1177",
    "1273321/1177", "518019/1199", "575695/1529", "897911/1661",
    "1254167/1661", "870187/2179", "1954913/2273", "-899645/2497",
    "42437/2497", "1361808/2651", "1415224/3377", "1647627/3641",
    "208925/3839", "-309043/3905", "1750765/4961", "-1404297/5533",
    "1675479/5599", "1354257/9691", "1850215/9977",
)))

EXPECTED_PIVOTS = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 25, 28, 30,
)
EXPECTED_EXCEPTIONAL_X = tuple(map(Q, (
    "-94898/11", "79501/73", "92446/307", "245933/341",
    "136343/349", "-250232/583", "454737/605", "379655/803",
)))


def quartic_value(x_value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(
        CONSTRUCTION.primitive_quartic_coefficients(PARAMETER)
    ):
        answer = answer*x_value+coefficient
    return answer


def rational_square_root(value: Fraction) -> Fraction:
    value = Q(value)
    if value < 0:
        raise AssertionError("a pinned quartic value became negative")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator*numerator != value.numerator or denominator*denominator != value.denominator:
        raise AssertionError("a pinned quartic value is not a rational square")
    return Q(numerator, denominator)


def main() -> None:
    assert CONSTRUCTION.quartic_condition == 0
    assert Q(16, 3)*NATIVE_PARAMETER == PARAMETER
    short_model = CONSTRUCTION.primitive_jacobian_coefficients(PARAMETER)
    change = WeierstrassChange.from_values(SHORT_TO_PUBLIC)
    assert change_weierstrass_model(short_model, change) == PUBLIC_MODEL

    searched = tuple(
        (x_value, rational_square_root(quartic_value(x_value)))
        for x_value in RATPOINTS_ABSCISSAS
    )
    # Fermigier's A+B*t formula becomes
    # x=(-12430-61*T)/41 in this canonical chart.
    extra_x = Q(-12430, 41)-Q(61, 41)*PARAMETER
    assert extra_x == Q(-344923, 451)
    extra = (extra_x, rational_square_root(quartic_value(extra_x)))

    points, sources, visible_count, prescribed_count = pool_with_sources(
        CONSTRUCTION, PARAMETER, searched, (extra,)
    )
    assert (len(points), visible_count, prescribed_count) == (54, 12, 13)
    certificate = mod3_independence_certificate(
        short_model, points, prime_bound=499
    )
    assert certificate["combined_exact_rank_over_F3"] == 20
    assert tuple(certificate["independent_subset_indices_one_based"]) == EXPECTED_PIVOTS

    exceptional_x = tuple(
        Q(sources[index-1]["quartic_point"]["x"])
        for index in EXPECTED_PIVOTS
        if sources[index-1]["source"].startswith("searched-")
    )
    assert exceptional_x == EXPECTED_EXCEPTIONAL_X
    assert sum(
        sources[index-1]["source"].startswith(("visible-", "forced-"))
        for index in EXPECTED_PIVOTS
    ) == 12

    getcontext().prec = 40
    score = Decimal(CONDUCTOR).ln()
    print(
        "ICARM275MESTRE|u=-3|v=-1/2|native_T=10239/176|"
        "roots=0,113,550,753,868,1058|canonical_T=3413/11|"
        f"pool={len(points)}|generic=12|exceptional=8|rank_lower=20|"
        f"certificate_primes={','.join(map(str, certificate['certificate_primes']))}|"
        f"ln_conductor={score}|status=PASS_EXACT_RANK20_DECOMPOSITION",
        flush=True,
    )


if __name__ == "__main__":
    main()
