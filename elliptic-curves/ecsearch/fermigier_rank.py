"""Rank-evaluation helpers for Fermigier family specializations."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt, lcm
from pathlib import Path
import re
import shutil
import subprocess
from typing import Sequence

from .fermigier import (
    FermigierQuartic,
    evaluate_polynomial,
    fermigier_canonical_coefficients,
    fermigier_quartic,
    quartic_point_to_canonical_point,
    thirteen_visible_points,
)
from .rank_certification import (
    AffinePoint,
    IndependenceCertificate,
    build_independence_certificate,
    negate_rational_point,
    subtract_rational_points,
)


@dataclass(frozen=True)
class FermigierRankSpecialization:
    adapter_parameter: Fraction
    quartic_model: FermigierQuartic
    quartic_points: tuple[AffinePoint, ...]
    canonical_model: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
    canonical_points: tuple[AffinePoint, ...]
    section_differences: tuple[AffinePoint, ...]


def specialize_fermigier_rank_sections(
    adapter_parameter: Fraction | int,
) -> FermigierRankSpecialization:
    """Specialize the thirteen quartic points and twelve Jacobian differences."""

    adapter_parameter = Fraction(adapter_parameter)
    quartic_model = fermigier_quartic(2 * adapter_parameter)
    quartic_points = thirteen_visible_points(quartic_model)
    canonical_model = fermigier_canonical_coefficients(adapter_parameter)
    canonical_points = tuple(
        quartic_point_to_canonical_point(quartic_model, point)
        for point in quartic_points
    )
    origin = canonical_points[0]
    differences: list[AffinePoint] = []
    for point in canonical_points[1:]:
        difference = subtract_rational_points(canonical_model, point, origin)
        if difference is None:
            raise ArithmeticError("two specialized quartic points had the same image")
        differences.append(difference)
    return FermigierRankSpecialization(
        adapter_parameter=adapter_parameter,
        quartic_model=quartic_model,
        quartic_points=quartic_points,
        canonical_model=canonical_model,
        canonical_points=canonical_points,
        section_differences=tuple(differences),
    )


def certify_fermigier_rank_sections(
    specialization: FermigierRankSpecialization,
    *,
    relation_primes: Sequence[int] = (5, 3, 7),
    maximum_reduction_prime: int = 2000,
) -> IndependenceCertificate:
    """Try small relation primes until the twelve differences are certified.

    The covariant images are divisible by two, so a mod-2 certificate cannot
    see their independence; mod 5 is the useful deterministic first attempt.
    """

    failures: list[str] = []
    for relation_prime in relation_primes:
        try:
            return build_independence_certificate(
                specialization.canonical_model,
                specialization.section_differences,
                relation_prime=relation_prime,
                maximum_reduction_prime=maximum_reduction_prime,
            )
        except ArithmeticError as error:
            failures.append(f"ell={relation_prime}: {error}")
    raise ArithmeticError("; ".join(failures))


def search_quartic_points_with_gp(
    model: FermigierQuartic,
    height_bound: int,
) -> tuple[AffinePoint, ...]:
    """Run PARI's bounded ``hyperellratpoints`` search on a quartic model.

    The bound is a search limit, not a completeness or rank certificate.
    PARI returns both ordinate signs when both are present.
    """

    if not isinstance(height_bound, int) or height_bound <= 0:
        raise ValueError("height_bound must be a positive integer")
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP executable 'gp' is required for point search")
    terms = [
        f"({coefficient.numerator}/{coefficient.denominator})*x^{degree}"
        for degree, coefficient in enumerate(model.quartic)
    ]
    program = (
        "x='x;r="
        + "+".join(terms)
        + f";V=hyperellratpoints(r,{height_bound});"
        + 'for(i=1,#V,print(V[i][1],"\\t",V[i][2]));\n'
    )
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise RuntimeError(completed.stdout + completed.stderr)
    points: list[AffinePoint] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise RuntimeError(f"unexpected PARI point-search output: {line!r}")
        points.append((Fraction(fields[0]), Fraction(fields[1])))
    return tuple(points)


def _rational_square_root(value: Fraction) -> Fraction:
    if value < 0:
        raise ArithmeticError("a negative rational is not a rational square")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if (
        numerator * numerator != value.numerator
        or denominator * denominator != value.denominator
    ):
        raise ArithmeticError("the rational value is not a square")
    return Fraction(numerator, denominator)


def parse_ratpoints_output(
    model: FermigierQuartic, output: str
) -> tuple[AffinePoint, ...]:
    """Parse quiet abscissa-only `ratpoints` output and check both ordinates."""

    pattern = re.compile(r"\((-?\d+) : (\d+)\)")
    points: list[AffinePoint] = []
    for line in output.splitlines():
        match = pattern.fullmatch(line.strip())
        if match is None:
            if line.strip():
                raise RuntimeError(f"unexpected ratpoints output: {line!r}")
            continue
        x_coordinate = Fraction(int(match[1]), int(match[2]))
        y_coordinate = _rational_square_root(
            evaluate_polynomial(model.quartic, x_coordinate)
        )
        points.append((x_coordinate, y_coordinate))
        if y_coordinate:
            points.append((x_coordinate, -y_coordinate))
    return tuple(points)


def search_quartic_points_with_ratpoints(
    model: FermigierQuartic,
    height_bound: int,
    *,
    denominator_bound: int | None = None,
    executable: str = "ratpoints",
) -> tuple[AffinePoint, ...]:
    """Run an optional `ratpoints` bounded search and reconstruct ordinates.

    `ratpoints` is not a repository dependency.  The caller must install it
    and, when using a nonstandard shared-library location, configure the
    process environment.  Its output supplies abscissas; ordinates are checked
    and reconstructed with exact rational arithmetic here.
    """

    if not isinstance(height_bound, int) or height_bound <= 0:
        raise ValueError("height_bound must be a positive integer")
    if denominator_bound is not None and (
        not isinstance(denominator_bound, int) or denominator_bound <= 0
    ):
        raise ValueError("denominator_bound must be a positive integer")
    resolved = shutil.which(executable)
    if resolved is None:
        candidate = Path(executable)
        if not candidate.is_file():
            raise RuntimeError(f"ratpoints executable {executable!r} was not found")
        resolved = str(candidate)

    common_denominator = 1
    for coefficient in model.quartic:
        common_denominator = lcm(common_denominator, coefficient.denominator)
    integral_coefficients = [
        int(coefficient * common_denominator**2)
        for coefficient in model.quartic
    ]
    command = [
        resolved,
        " ".join(map(str, integral_coefficients)),
        str(height_bound),
    ]
    if denominator_bound is not None:
        command.extend(("-du", str(denominator_bound)))
    command.extend(("-q", "-y"))
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=True,
    )
    return parse_ratpoints_output(model, completed.stdout)


def section_and_point_cloud_differences(
    specialization: FermigierRankSpecialization,
    quartic_points: Sequence[AffinePoint],
) -> tuple[AffinePoint, ...]:
    """Put searched quartic points after the twelve baseline differences.

    Exact duplicates and inverses are removed because they cannot increase
    the generated subgroup rank.  A ramification point where the covariant
    formula has zero denominator is skipped.
    """

    curve = specialization.canonical_model
    origin = specialization.canonical_points[0]
    differences = list(specialization.section_differences)
    seen = set(differences)
    seen.update(negate_rational_point(curve, point) for point in differences)
    for quartic_point in quartic_points:
        try:
            canonical_point = quartic_point_to_canonical_point(
                specialization.quartic_model, quartic_point
            )
        except ArithmeticError:
            continue
        difference = subtract_rational_points(curve, canonical_point, origin)
        if difference is None or difference in seen:
            continue
        differences.append(difference)
        seen.add(difference)
        seen.add(negate_rational_point(curve, difference))
    return tuple(differences)


def write_json_exclusively(path: Path, text: str) -> None:
    """Write a generated evaluator result without silently refreshing it."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as output:
        output.write(text)
