"""Exact finite-normalization exports for MBP pipeline regressions.

These are family-specific adapters.  The extractor in
``jcsearch.minimal_boundary`` does not import this module and therefore
cannot use the fixture or mechanism names while deciding predicates.
"""

from __future__ import annotations

from dataclasses import replace
from math import gcd

import sympy as sp

from jcsearch.boundary import (
    cancellation_boundary_cover_profile,
    weighted_boundary_cover_profile,
)
from jcsearch.minimal_boundary import (
    ChartRecord,
    CollisionRecord,
    ConormalRecord,
    CriticalCurveRecord,
    FiniteNormalizationRecord,
    LedgerPrime,
    LinkRecord,
    PrimeRecord,
    mutate_record,
)


P, S, Q, w, q, t = sp.symbols("P S Q w q t")


def _prime_from_profile(profile, *, critical: bool) -> PrimeRecord:
    return PrimeRecord(
        label=profile.label,
        boundary=True,
        critical_fitting=critical,
        quotient_image=(
            "critical-core" if critical else f"aux:{profile.target_divisor}"
        ),
        target_image=profile.target_divisor,
        ramification_index=profile.ramification_index,
        residue_degree=profile.residue_degree,
        different_exponent=profile.different_exponent,
        color="ramified" if profile.different_exponent else "unramified",
        completed_incidence=profile.inertia_cycle,
    )


def _link(orientation: int) -> LinkRecord:
    values = (0, 1, 2)
    return LinkRecord(
        label="rank-one-boundary-link",
        normal_source=True,
        normal_target=True,
        source_ufd=True,
        target_ufd=True,
        scalar_source_units=True,
        scalar_target_units=True,
        source_parameter_prime=True,
        target_parameter_prime=True,
        unit_lattice_rank=1,
        target_parameter_support_in_source=("a",),
        source_parameter_support_in_target=("d",),
        orientation_exponent=orientation,
        forward_generator_valuations=values if orientation == 1 else (),
        reverse_generator_valuations=values if orientation == -1 else (),
    )


def _line_curve() -> CriticalCurveRecord:
    return CriticalCurveRecord(
        geometrically_integral=True,
        smooth=True,
        genus=0,
        parameter=t,
        punctures=("infinity",),
        unit_generators=(),
    )


def _torus_curve() -> CriticalCurveRecord:
    return CriticalCurveRecord(
        geometrically_integral=True,
        smooth=True,
        genus=0,
        parameter=t,
        punctures=("zero", "infinity"),
        unit_generators=(t,),
    )


def _line_conormal() -> ConormalRecord:
    return ConormalRecord(
        height_one_content_vectors=((1,),),
        collisions=(CollisionRecord("generic-collision", 1, ((1,),)),),
        residue_coefficient=t,
    )


def _torus_conormal() -> ConormalRecord:
    return ConormalRecord(
        height_one_content_vectors=((1,),),
        collisions=(CollisionRecord("generic-collision", 1, ((1,),)),),
        residue_coefficient=t,
    )


def weighted_record(degree: int) -> FiniteNormalizationRecord:
    """Export the weighted tangent core of degree ``degree``."""

    if degree < 3:
        raise ValueError("weighted regression degree must be at least three")
    H = sp.expand(w**degree * (1 - w))
    h = sp.diff(H, w)
    profile = weighted_boundary_cover_profile(H, w)
    primes = tuple(
        _prime_from_profile(prime, critical=prime.label == "E_Delta")
        for prime in profile.primes
    )
    chart = ChartRecord(
        source_variables=(w, q),
        target_expressions=(q, sp.expand(w * q - H)),
        controlled_divisor=sp.expand(q - h),
        controlled_exponent=1,
        reciprocal_chart=False,
        fitting_support_rows=((1,),),
    )
    return FiniteNormalizationRecord(
        name=f"weighted-degree-{degree}",
        primes=primes,
        links=(_link(1),),
        ledger=(
            LedgerPrime("gamma", False, True, 2, 1, 3),
            LedgerPrime("E_zero_cluster", False, True, 3, 0, 3),
        ),
        unrecorded_graph_primes=(),
        critical_curve=_line_curve(),
        conormal=_line_conormal(),
        residue_parameter=t,
        residue_generators=(t,),
        chart=chart,
    )


def cancellation_record(m: int, r: int) -> FiniteNormalizationRecord:
    """Export one exact cancellation normalization profile."""

    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    profile = cancellation_boundary_cover_profile(m, r)
    primes = tuple(
        _prime_from_profile(prime, critical=prime.label == "E_Delta")
        for prime in profile.primes
    )
    D = sp.expand(1 - S * (Q - P * S) ** m)
    u = sp.Dummy("u")
    R = sp.expand(sp.integrate(D.subs(S, u) ** r, (u, 0, S)))
    auxiliary = tuple(
        LedgerPrime(prime.label, False, True, 0, 0, 0)
        for prime in profile.primes
        if prime.label != "E_Delta"
    )
    return FiniteNormalizationRecord(
        name=f"cancellation-m{m}-r{r}",
        primes=primes,
        links=(_link(-1),),
        ledger=(
            LedgerPrime("E_Delta", False, True, -r, r, 0),
        )
        + auxiliary,
        unrecorded_graph_primes=(),
        critical_curve=_torus_curve(),
        conormal=_torus_conormal(),
        residue_parameter=t,
        residue_generators=(t, 1 / t),
        chart=ChartRecord(
            source_variables=(P, S, Q),
            target_expressions=(P, Q, R),
            controlled_divisor=D,
            controlled_exponent=r,
            reciprocal_chart=True,
            fitting_support_rows=((1, -1),),
        ),
    )


def quadratic_gauge_record(degree: int) -> FiniteNormalizationRecord:
    """Export a root-engineered quadratic-incidence profile."""

    if degree < 3:
        raise ValueError("quadratic-gauge degree must be at least three")
    U = S + P * S**3
    if degree > 3:
        U += P**degree * S**degree
    beta = sp.cancel((sp.diff(U, S) - 1 - P * S**2) / S)
    B = sp.expand(Q + beta)
    C = sp.expand(2 * U - B * S**2)
    D = sp.expand(1 - S * Q + P * S**2)
    primes = [
        PrimeRecord(
            label="E_Delta",
            boundary=True,
            critical_fitting=True,
            quotient_image="critical-core",
            target_image="Z_Delta",
            ramification_index=2,
            residue_degree=1,
            different_exponent=1,
            color="ramified",
            completed_incidence=(2,) + (1,) * (degree - 2),
        )
    ]
    if degree >= 4:
        second_vertex_degree = degree - 3
        prime_count = gcd(second_vertex_degree, 2)
        ramification_index = second_vertex_degree // prime_count
        for index in range(prime_count):
            primes.append(
                PrimeRecord(
                    label=f"E_P_{index}",
                    boundary=True,
                    critical_fitting=False,
                    quotient_image="aux:Z_0",
                    target_image="Z_0",
                    ramification_index=ramification_index,
                    residue_degree=1,
                    different_exponent=ramification_index - 1,
                    color=(
                        "ramified"
                        if ramification_index > 1
                        else "unramified"
                    ),
                    completed_incidence=(ramification_index,)
                    + (1,) * (degree - ramification_index),
                )
            )
    auxiliary_ledger = tuple(
        LedgerPrime(prime.label, False, True, 0, 0, 0)
        for prime in primes
        if prime.label != "E_Delta"
    )
    return FiniteNormalizationRecord(
        name=f"quadratic-gauge-degree-{degree}",
        primes=tuple(primes),
        links=(_link(-1),),
        ledger=(LedgerPrime("E_Delta", False, True, -1, 1, 0),)
        + auxiliary_ledger,
        unrecorded_graph_primes=(),
        critical_curve=_torus_curve(),
        conormal=_torus_conormal(),
        residue_parameter=t,
        residue_generators=(t, 1 / t),
        chart=ChartRecord(
            source_variables=(P, S, Q),
            target_expressions=(P, B, C),
            controlled_divisor=D,
            controlled_exponent=1,
            reciprocal_chart=True,
            fitting_support_rows=((1, 0), (0, 1)),
        ),
    )


def countermodels() -> tuple[FiniteNormalizationRecord, ...]:
    """Return single-defect perturbations and a spectator-prime model."""

    weighted = weighted_record(4)
    assert weighted.chart is not None
    w_var, q_var = weighted.chart.source_variables
    perturbed_chart = replace(
        weighted.chart,
        target_expressions=(
            q_var,
            weighted.chart.target_expressions[1] + w_var**2 * q_var**2,
        ),
    )
    chart_perturbation = mutate_record(
        weighted,
        name="countermodel-chart-perturbation",
        chart=perturbed_chart,
    )

    cancellation = cancellation_record(2, 1)
    conormal_perturbation = mutate_record(
        cancellation,
        name="countermodel-imprimitive-conormal",
        conormal=replace(
            cancellation.conormal,
            height_one_content_vectors=((2,),),
            residue_coefficient=t**2,
        ),
    )
    contracted = mutate_record(
        cancellation,
        name="countermodel-contracted-residue",
        residue_generators=(sp.Integer(1),),
    )
    nonsaturated = mutate_record(
        cancellation,
        name="countermodel-nonsaturated-link",
        links=(replace(cancellation.links[0], orientation_exponent=-2),),
    )

    spectator_prime = replace(
        cancellation.primes[0],
        label="E_spectator",
        automorphism_orbit_size=1,
    )
    spectator = mutate_record(
        cancellation,
        name="countermodel-spectator-prime",
        primes=cancellation.primes + (spectator_prime,),
        ledger=cancellation.ledger
        + (LedgerPrime("E_spectator", False, False, 0, 0, 0),),
    )
    return (
        chart_perturbation,
        conormal_perturbation,
        contracted,
        nonsaturated,
        spectator,
    )
