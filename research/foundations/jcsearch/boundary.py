"""Boundary saturation data for weighted inverse primitives."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class BoundarySaturationProfile:
    """Multiplicity data controlling the pulled-back discriminant at C=0."""

    degree: int
    zero_multiplicity: int
    one_multiplicity: int
    extra_degree: int
    extra_repetition_excess: int
    saturation_exponent: int
    zero_leading_coefficient: sp.Expr

    def boundary_trace(self, A, B, c=sp.Integer(1)):
        """Return the boundary trace of Q up to a nonzero scalar."""
        if self.zero_multiplicity == 2:
            conic = B**2 - 4 * sp.sympify(c) * self.zero_leading_coefficient * A
            return sp.expand(B**self.extra_repetition_excess * conic)
        return B**self.saturation_exponent


@dataclass(frozen=True)
class BoundaryPrimeProfile:
    """Generic formal data of one canonical upstairs boundary prime."""

    label: str
    target_divisor: str
    ramification_index: int
    residue_degree: int
    different_exponent: int
    local_parameter_equation: str
    inertia_cycle: tuple[int, ...]
    arithmetic_factor: sp.Expr


@dataclass(frozen=True)
class UpstairsBoundaryCoverProfile:
    """Computable generic divisorial layer of the full boundary invariant."""

    degree: int
    target_divisors: tuple[str, ...]
    primes: tuple[BoundaryPrimeProfile, ...]
    degree_sums: tuple[tuple[str, int, int, int], ...]


@dataclass(frozen=True)
class BoundaryIntersectionBranchProfile:
    """One branch in a completed critical-normalization intersection."""

    label: str
    divisor_multiplicity: int
    contact_order: int
    length_contribution: int


@dataclass(frozen=True)
class CancellationBoundaryIntersectionProfile:
    """Exact generic higher-stratum contact for cancellation boundaries.

    At the generic point of ``P=Q=0`` the target completion is ``K[[P,Q]]``
    with ``K=k(R)``.  The regular critical-normalization plane
    ``K[[Y,Q]]`` pulls ``P`` back to ``(Q-Y)Y^m``.  The two displayed branch
    contributions retain information lost by the single target nilpotency
    index.
    """

    m: int
    r: int
    target_completion: str
    normalization_completion: str
    pulled_back_boundary_equation: sp.Expr
    branches: tuple[BoundaryIntersectionBranchProfile, ...]
    target_intersection_polynomial: sp.Expr
    nilpotency_index: int


@dataclass(frozen=True)
class BoundaryPrimeIntersection:
    """One scheme-theoretic edge in an upstairs boundary diagram."""

    left: str
    right: str
    support: str
    contact_length: int
    completed_intersection_ring: str


@dataclass(frozen=True)
class CancellationFullBoundaryDiagram:
    """Certified geometric intersection diagram of all cancellation primes.

    The labels describe geometric primes after strict henselization.  The
    ``E_A`` branch is retained because it belongs to the same finite
    ``K``-cluster as the boundary primes, although it lies in the affine
    source.  It makes the arithmetic length calculation transparent.
    """

    m: int
    r: int
    branch_polynomial: sp.Expr
    critical_coefficient: sp.Expr
    branch_discriminant: sp.Expr
    contact_resultant: sp.Expr
    contact_certified: bool
    central_stratum: str
    boundary_primes: tuple[str, ...]
    cluster_primes: tuple[str, ...]
    intersections: tuple[BoundaryPrimeIntersection, ...]
    higher_intersection_ring: str
    critical_cluster_total_length: int
    critical_boundary_total_length: int


def cancellation_branch_polynomial(m: int, r: int, variable=None):
    """Return the exact ``K_{m,r}`` polynomial for the ``P=0`` cluster."""
    m, r = int(m), int(r)
    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    w = variable if variable is not None else sp.Symbol("w")
    v = sp.Dummy("v")
    numerator = sp.integrate(v**r * (1 - v) ** (m * r), (v, 0, w))
    return sp.factor(sp.cancel(numerator / w ** (r + 1)))


def cancellation_critical_coefficient(m: int, r: int, variable=None):
    """Return the coefficient controlling critical/K-cluster contact.

    With ``p(w)=w(1-w)^m``, this is

    ``w^(-r-1) integral_0^w (p(w)-p(v))^r dv``.

    Its image in ``k[w]/(K_{m,r})`` is a unit exactly when the resultant in
    :func:`cancellation_full_boundary_diagram` is nonzero.
    """
    m, r = int(m), int(r)
    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    w = variable if variable is not None else sp.Symbol("w")
    v = sp.Dummy("v")
    p_w = w * (1 - w) ** m
    p_v = v * (1 - v) ** m
    numerator = sp.integrate((p_w - p_v) ** r, (v, 0, w))
    return sp.factor(sp.cancel(numerator / w ** (r + 1)))


def cancellation_full_boundary_diagram(m: int, r: int):
    """Extract every geometric prime intersection for a cancellation map.

    The calculation is exact over ``QQ``.  It certifies the uniform contact
    assertion precisely when ``Res(K_{m,r}, L_{m,r})`` is nonzero.  In that
    case the critical prime restricts to ``Q^m`` times a unit on the finite
    etale K-cluster, so each geometric cluster branch has length ``m``.
    Distinct cluster branches force ``Q=0`` and therefore meet in the reduced
    central stratum.  The cubic ``(m,r)=(1,1)`` has no ``P=0`` boundary prime,
    but its affine K-branch is still included in the length certificate.
    """
    m, r = int(m), int(r)
    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    w = sp.Symbol("w")
    K = cancellation_branch_polynomial(m, r, w)
    L = cancellation_critical_coefficient(m, r, w)
    K_poly = sp.Poly(K, w, domain=sp.QQ)
    L_poly = sp.Poly(L, w, domain=sp.QQ)
    discriminant = sp.factor(sp.discriminant(K_poly.as_expr(), w))
    resultant = sp.factor(sp.resultant(K_poly.as_expr(), L_poly.as_expr(), w))
    certified = discriminant != 0 and resultant != 0

    boundary_cluster = tuple(f"E_P_{index}" for index in range(m * r - 1))
    cluster = ("E_A",) + boundary_cluster
    boundary = ("E_Delta",) + boundary_cluster
    intersections: list[BoundaryPrimeIntersection] = []
    if certified:
        for prime in boundary_cluster:
            intersections.append(
                BoundaryPrimeIntersection(
                    left="E_Delta",
                    right=prime,
                    support="S={P=Q=0}",
                    contact_length=m,
                    completed_intersection_ring=f"k(R)[[Q]]/(Q^{m})",
                )
            )
        for left_index, left in enumerate(boundary_cluster):
            for right in boundary_cluster[left_index + 1 :]:
                intersections.append(
                    BoundaryPrimeIntersection(
                        left=left,
                        right=right,
                        support="S={P=Q=0}",
                        contact_length=1,
                        completed_intersection_ring="k(R)",
                    )
                )

    return CancellationFullBoundaryDiagram(
        m=m,
        r=r,
        branch_polynomial=K,
        critical_coefficient=L,
        branch_discriminant=discriminant,
        contact_resultant=resultant,
        contact_certified=certified,
        central_stratum="S={P=Q=0}, with R free",
        boundary_primes=boundary,
        cluster_primes=cluster,
        intersections=tuple(intersections),
        higher_intersection_ring="k(R)",
        critical_cluster_total_length=m * m * r,
        critical_boundary_total_length=m * (m * r - 1),
    )


def weighted_boundary_cover_profile(primitive, variable):
    """Extract generic boundary primes, differents, and inertia for a seed.

    Factor degrees are arithmetic residue degrees.  Over an algebraic closure
    each irreducible factor splits into geometric primes of residue degree one.
    The displayed completed-DVR equations are understood up to multiplication
    by a unit; characteristic zero makes every extension tame, so the
    different exponent is ``e-1``.
    """
    saturation = boundary_saturation_profile(primitive, variable)
    polynomial = sp.Poly(sp.sympify(primitive), variable)
    reduced = polynomial
    for _ in range(saturation.zero_multiplicity):
        reduced = reduced.exquo(sp.Poly(variable, variable))
    reduced = reduced.exquo(sp.Poly(variable - 1, variable))

    primes = [
        BoundaryPrimeProfile(
            label="E_Delta",
            target_divisor="Z_Delta",
            ramification_index=2,
            residue_degree=1,
            different_exponent=1,
            local_parameter_equation="delta = unit*q^2",
            inertia_cycle=(2,) + (1,) * (saturation.degree - 2),
            arithmetic_factor=sp.Integer(1),
        )
    ]
    for index, (factor, multiplicity) in enumerate(reduced.factor_list()[1]):
        factor_expression = sp.factor(factor.as_expr())
        c = -sp.diff(polynomial.as_expr(), variable).subs(variable, 1)
        escape_numerator = sp.expand(c * variable + sp.diff(polynomial.as_expr(), variable))
        if sp.resultant(factor_expression, escape_numerator, variable) == 0:
            raise ValueError(
                "an extra-root factor fails the boundary-clean condition"
            )
        primes.append(
            BoundaryPrimeProfile(
                label=f"E_extra_{index}",
                target_divisor="Z_0",
                ramification_index=int(multiplicity),
                residue_degree=int(factor.degree()),
                different_exponent=int(multiplicity) - 1,
                local_parameter_equation=(
                    f"C = unit*q^{int(multiplicity)}"
                ),
                inertia_cycle=(int(multiplicity),)
                + (1,) * (saturation.degree - int(multiplicity)),
                arithmetic_factor=factor_expression,
            )
        )
    if saturation.zero_multiplicity >= 3:
        exponent = saturation.zero_multiplicity - 1
        primes.append(
            BoundaryPrimeProfile(
                label="E_zero_cluster",
                target_divisor="Z_0",
                ramification_index=exponent,
                residue_degree=1,
                different_exponent=exponent - 1,
                local_parameter_equation=f"C = unit*q^{exponent}",
                inertia_cycle=(exponent,)
                + (1,) * (saturation.degree - exponent),
                arithmetic_factor=variable,
            )
        )

    delta_boundary_degree = 2
    zero_boundary_degree = sum(
        prime.ramification_index * prime.residue_degree
        for prime in primes
        if prime.target_divisor == "Z_0"
    )
    return UpstairsBoundaryCoverProfile(
        degree=saturation.degree,
        target_divisors=("Z_Delta", "Z_0"),
        primes=tuple(primes),
        degree_sums=(
            (
                "Z_Delta",
                delta_boundary_degree,
                saturation.degree - delta_boundary_degree,
                saturation.degree,
            ),
            (
                "Z_0",
                zero_boundary_degree,
                saturation.degree - zero_boundary_degree,
                saturation.degree,
            ),
        ),
    )


def cancellation_boundary_cover_profile(m: int, r: int):
    """Return the geometric generic boundary profile of cancellation type.

    This is the divisorial layer over an algebraic closure.  Arithmetic
    grouping of the unramified ``P=0`` primes is separately determined by the
    factorization of the cancellation parameter polynomial.
    """
    m, r = int(m), int(r)
    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    degree = r * (m + 1) + 1
    critical_e = r + 1
    primes = [
        BoundaryPrimeProfile(
            label="E_Delta",
            target_divisor="Z_Delta",
            ramification_index=critical_e,
            residue_degree=1,
            different_exponent=r,
            local_parameter_equation=f"delta = unit*q^{critical_e}",
            inertia_cycle=(critical_e,) + (1,) * (degree - critical_e),
            arithmetic_factor=sp.Integer(1),
        )
    ]
    target_divisors = ["Z_Delta"]
    degree_sums = [
        ("Z_Delta", critical_e, degree - critical_e, degree)
    ]
    if (m, r) != (1, 1):
        target_divisors.append("Z_0")
        for index in range(m * r - 1):
            primes.append(
                BoundaryPrimeProfile(
                    label=f"E_P_{index}",
                    target_divisor="Z_0",
                    ramification_index=1,
                    residue_degree=1,
                    different_exponent=0,
                    local_parameter_equation="P = unit*q",
                    inertia_cycle=(1,) * degree,
                    arithmetic_factor=sp.Integer(1),
                )
            )
        degree_sums.append(
            ("Z_0", m * r - 1, r + 2, degree)
        )
    return UpstairsBoundaryCoverProfile(
        degree=degree,
        target_divisors=tuple(target_divisors),
        primes=tuple(primes),
        degree_sums=tuple(degree_sums),
    )


def cancellation_boundary_intersection_profile(m: int, r: int):
    """Return the exact completed generic contact of ``Z_Delta`` and ``Z_0``.

    This packages the all-parameter completed-local calculation, including
    its two normalization branches.  It is the computed critical-value
    higher stratum; it does not pretend to supply still-unknown intersections
    among every prime of the full source normalization cover.
    """
    m, r = int(m), int(r)
    if m < 1 or r < 1:
        raise ValueError("cancellation parameters must be positive")
    Y, Q, R, C = sp.symbols("Y Q R C")
    zero_contact = m * r
    diagonal_contact = m * r
    zero_contribution = m * zero_contact
    diagonal_contribution = diagonal_contact
    total = zero_contribution + diagonal_contribution
    return CancellationBoundaryIntersectionProfile(
        m=m,
        r=r,
        target_completion="k(R)[[P,Q]]",
        normalization_completion="k(R)[[Y,Q]]",
        pulled_back_boundary_equation=sp.expand((Q - Y) * Y**m),
        branches=(
            BoundaryIntersectionBranchProfile(
                label="Y=0",
                divisor_multiplicity=m,
                contact_order=zero_contact,
                length_contribution=zero_contribution,
            ),
            BoundaryIntersectionBranchProfile(
                label="Y=Q",
                divisor_multiplicity=1,
                contact_order=diagonal_contact,
                length_contribution=diagonal_contribution,
            ),
        ),
        target_intersection_polynomial=sp.expand(
            Q**total * ((r + 1) * R * Q**m - C)
        ),
        nilpotency_index=total,
    )


def boundary_saturation_profile(primitive, variable) -> BoundarySaturationProfile:
    """Compute the exact C-saturation exponent from the roots of ``H``.

    This applies to admissible primitives: the zero at 0 has multiplicity at
    least two and the distinguished zero at 1 is simple.  Extra factors may be
    simple or repeated and need not split over the coefficient field.
    """
    polynomial = sp.Poly(sp.sympify(primitive), variable)
    if polynomial.domain.characteristic() != 0:
        raise ValueError("boundary saturation requires characteristic zero")

    reduced = polynomial
    zero_multiplicity = 0
    zero_factor = sp.Poly(variable, variable)
    while reduced.eval(0) == 0:
        reduced = reduced.exquo(zero_factor)
        zero_multiplicity += 1

    one_multiplicity = 0
    one_factor = sp.Poly(variable - 1, variable)
    while reduced.eval(1) == 0:
        reduced = reduced.exquo(one_factor)
        one_multiplicity += 1

    if zero_multiplicity < 2:
        raise ValueError("an admissible primitive must have a double zero at 0")
    if one_multiplicity != 1:
        raise ValueError("the distinguished zero at 1 must be simple")

    extra_repetition_excess = sum(
        factor.degree() * (multiplicity - 1)
        for factor, multiplicity in reduced.sqf_list()[1]
    )
    saturation_exponent = zero_multiplicity + extra_repetition_excess
    zero_leading_coefficient = polynomial.coeff_monomial(variable**zero_multiplicity)

    return BoundarySaturationProfile(
        degree=polynomial.degree(),
        zero_multiplicity=zero_multiplicity,
        one_multiplicity=one_multiplicity,
        extra_degree=reduced.degree(),
        extra_repetition_excess=extra_repetition_excess,
        saturation_exponent=saturation_exponent,
        zero_leading_coefficient=zero_leading_coefficient,
    )
