#!/usr/bin/env python3
"""Exact direct-reduction backend for pointed half-lattice quartics.

The historical half-lattice backend first called ``hyperellminimalmodel`` and
then ``hyperellred`` under one wall timeout.  On the large ICARM/MW16
prospective quartics the first operation consumes the entire budget before a
bounded point search starts.  This experimental backend instead:

1. clears denominators exactly;
2. removes the whole coefficient content only when it is provably a square;
3. calls PARI ``hyperellred`` directly; and
4. runs the same bounded ``hyperellratpoints`` search on that reduced model.

Every returned point is transported through PARI's exact inverse map, through
the recorded ordinate scaling, and through the exact pointed-quartic map to
the original short elliptic curve.  A bounded miss remains only a bounded
miss.  No minimal-model, Selmer, saturation, or rank-upper-bound claim is
made.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from math import gcd, isqrt, lcm
from pathlib import Path
import subprocess
import time
from typing import Any, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
base = SourceFileLoader("half_lattice_direct_reduction_base", str(BASE_SOURCE)).load_module()


BACKEND_NAME = "exact_square_content_then_direct_hyperellred_v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _timeout_text(value) -> str:
    if value is None:
        return ""
    return value.decode(errors="replace") if isinstance(value, bytes) else str(value)


def _common_record(
    *,
    mask: int,
    representative: Sequence[int],
    base_point,
    cover,
    denominator: int,
    integral_coefficients: Sequence[int],
    content: int,
    content_root: int | None,
    ordinate_scale: int,
    normalized_coefficients: Sequence[int],
    height_bound: int,
    timeout_seconds: float,
    wall_seconds: float,
) -> dict[str, Any]:
    return {
        "mask": mask,
        "hex": f"0x{mask:05x}",
        "height_bound": height_bound,
        "timeout_seconds": timeout_seconds,
        "wall_seconds": wall_seconds,
        "representative": list(map(int, representative)),
        "base_point": base.point_record(base_point),
        "raw_quartic_coefficients_ascending": [
            base.rational_to_string(value) for value in cover.coefficients
        ],
        "raw_rational_coefficient_maximum_bits": max(
            base.bit_height(value) for value in cover.coefficients
        ),
        "denominator_clearing_factor_bits": denominator.bit_length(),
        "integral_model_maximum_coefficient_bits": max(
            abs(value).bit_length() for value in integral_coefficients
        ),
        "coefficient_content_bits": content.bit_length(),
        "coefficient_content_is_square": content_root is not None,
        "removed_square_content_bits": (
            content.bit_length() if content_root is not None else 0
        ),
        "ordinate_scale": str(ordinate_scale),
        "normalized_integral_model_maximum_coefficient_bits": max(
            abs(value).bit_length() for value in normalized_coefficients
        ),
        "backend": BACKEND_NAME,
        "hyperellminimalmodel_called": False,
    }


def run_quartic_search(
    *,
    mask: int,
    representative: Sequence[int],
    short_model,
    generic_points,
    height_bound: int,
    timeout_seconds: float,
    stack_bytes: int,
):
    """Run one exact direct-reduction pointed-quartic search."""

    base_point = base.exact_linear_combination(
        Fraction(short_model[3]), generic_points, representative
    )
    if base_point is None:
        raise ArithmeticError("a nonzero selected class produced the point at infinity")
    cover = base.alternate_cover(short_model, base_point)

    denominator = 1
    for coefficient in cover.coefficients:
        denominator = lcm(denominator, Fraction(coefficient).denominator)
    integral_coefficients = tuple(
        Fraction(coefficient) * denominator * denominator
        for coefficient in cover.coefficients
    )
    if any(value.denominator != 1 for value in integral_coefficients):
        raise ArithmeticError("quartic denominator clearing failed")
    integral_coefficients = tuple(int(value) for value in integral_coefficients)

    content = 0
    for coefficient in integral_coefficients:
        content = gcd(content, abs(coefficient))
    if content <= 0:
        raise ArithmeticError("integral quartic has zero coefficient content")
    root = isqrt(content)
    content_root = root if root * root == content else None
    if content_root is not None:
        if denominator % content_root:
            raise ArithmeticError("square-content ordinate scale is nonintegral")
        normalized_coefficients = tuple(
            coefficient // content for coefficient in integral_coefficients
        )
        ordinate_scale = denominator // content_root
    else:
        normalized_coefficients = integral_coefficients
        ordinate_scale = denominator

    polynomial = base.gp_polynomial(
        tuple(Fraction(value) for value in normalized_coefficients)
    )
    x_base, y_base = base_point
    program = f"""
print("PHASE|reduction");
C0=[{polynomial},0];
gettime(); C2=hyperellred(C0,&m2); redms=gettime();
print("REDMS|",redms);
print("REDP|",Vec(C2[1]));
print("REDQ|",Vec(C2[2]));
print("REDDISC|",hyperelldisc(C2));
print("PHASE|search");
gettime(); R=hyperellratpoints(C2,{height_bound}); searchms=gettime();
print("SEARCHMS|",searchms);
print("SEARCHCOUNT|",#R);
for(i=1,#R,p=R[i];z=m2[2][2,1]*p[1]+m2[2][2,2];if(z==0,print("TRIVIAL_INFINITY|stage=red"),p0=[(m2[2][1,1]*p[1]+m2[2][1,2])/z,(m2[1]*p[2]+subst(m2[3],x,p[1]))/z^2];if(!hyperellisoncurve(C0,p0),error("inverse model map left C0"));ex=(p0[1]^2-{base.gp_rational(x_base)}+p0[2]/{ordinate_scale})/2;ey=p0[1]*(ex-{base.gp_rational(x_base)})-{base.gp_rational(y_base)};print("POINT|",p0[1],"|",p0[2]/{ordinate_scale},"|",ex,"|",ey)));
quit
"""

    started = time.monotonic()
    try:
        completed = subprocess.run(
            ["gp", "-q", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
        wall_seconds = time.monotonic() - started
    except subprocess.TimeoutExpired as error:
        wall_seconds = time.monotonic() - started
        stdout = _timeout_text(error.stdout)
        phase = "search" if "PHASE|search" in stdout else "reduction"
        record = _common_record(
            mask=mask,
            representative=representative,
            base_point=base_point,
            cover=cover,
            denominator=denominator,
            integral_coefficients=integral_coefficients,
            content=content,
            content_root=content_root,
            ordinate_scale=ordinate_scale,
            normalized_coefficients=normalized_coefficients,
            height_bound=height_bound,
            timeout_seconds=timeout_seconds,
            wall_seconds=wall_seconds,
        )
        record.update(
            {
                "status": f"bounded_{phase}_timeout",
                "timeout_phase": phase,
            }
        )
        return base.QuarticSearchResult(record, ())

    common = _common_record(
        mask=mask,
        representative=representative,
        base_point=base_point,
        cover=cover,
        denominator=denominator,
        integral_coefficients=integral_coefficients,
        content=content,
        content_root=content_root,
        ordinate_scale=ordinate_scale,
        normalized_coefficients=normalized_coefficients,
        height_bound=height_bound,
        timeout_seconds=timeout_seconds,
        wall_seconds=wall_seconds,
    )
    if completed.returncode != 0 or "***" in completed.stderr:
        common.update(
            {
                "status": "pari_failure",
                "error": completed.stderr.strip()[-2000:],
            }
        )
        return base.QuarticSearchResult(common, ())

    markers: dict[str, str] = {}
    curve_points = []
    raw_points = []
    infinity_count = 0
    for line in completed.stdout.splitlines():
        if line.startswith("POINT|"):
            unused, raw_x, raw_y, curve_x, curve_y = line.split("|", 4)
            raw_point = (Fraction(raw_x), Fraction(raw_y))
            curve_point = (Fraction(curve_x), Fraction(curve_y))
            if raw_point[1] ** 2 != cover.value(raw_point[0]):
                raise ArithmeticError("mapped PARI point left the raw quartic")
            if cover.cover_point_to_curve(raw_point) != curve_point:
                raise ArithmeticError("PARI/Python quartic maps disagree")
            if not base.point_on_short_curve(short_model, curve_point):
                raise ArithmeticError("mapped quartic point left E")
            raw_points.append(raw_point)
            curve_points.append(curve_point)
        elif line.startswith("TRIVIAL_INFINITY|"):
            infinity_count += 1
        elif "|" in line:
            key, value = line.split("|", 1)
            markers[key] = value.strip()
    required = ("REDMS", "REDP", "REDQ", "REDDISC", "SEARCHMS", "SEARCHCOUNT")
    for marker in required:
        if marker not in markers:
            raise ArithmeticError(f"PARI omitted marker {marker} for mask {mask:#x}")

    red_p = list(reversed([int(value) for value in base.parse_gp_vector(markers["REDP"])]))
    red_q = list(reversed([int(value) for value in base.parse_gp_vector(markers["REDQ"])]))
    local_profile = [
        base.modular_square_density(red_p, red_q, prime)
        for prime in base.SQUARE_SIEVE_PRIMES
    ]
    common.update(
        {
            "status": "bounded_search_complete",
            "reduction_milliseconds": int(markers["REDMS"]),
            "reduced_model": {
                "P_coefficients_ascending": red_p,
                "Q_coefficients_ascending": red_q,
                "maximum_coefficient_bits": max(
                    abs(value).bit_length() for value in red_p + red_q
                ),
                "discriminant": markers["REDDISC"],
                "discriminant_bits": abs(int(markers["REDDISC"])).bit_length(),
            },
            "local_stage": {
                "solubility_filter": "not_applicable_pointed_model_birational_to_E",
                "reason": "the pointed quartic has rational points at infinity",
                "affine_modular_square_sieve_profile": local_profile,
                "joint_independent_density_product": str(
                    np.prod(
                        [
                            row["affine_x_survivors"] / row["prime"]
                            for row in local_profile
                        ]
                    )
                ),
            },
            "search_milliseconds": int(markers["SEARCHMS"]),
            "signed_affine_points_reported": int(markers["SEARCHCOUNT"]),
            "trivial_points_mapping_to_infinity": infinity_count,
            "finite_raw_points": [base.point_record(point) for point in raw_points],
            "finite_curve_points": [base.point_record(point) for point in curve_points],
        }
    )
    return base.QuarticSearchResult(common, tuple(curve_points))


def provenance() -> dict[str, str]:
    return {
        "backend": BACKEND_NAME,
        "backend_source": str(Path(__file__).resolve().relative_to(ROOT)),
        "backend_sha256": digest(Path(__file__).resolve()),
        "base_engine_source": str(BASE_SOURCE.relative_to(ROOT)),
        "base_engine_sha256": digest(BASE_SOURCE),
    }
