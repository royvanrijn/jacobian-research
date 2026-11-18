"""Replay exact local controls for quadratic radial coefficient roots."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from verify_gaussian_quadratic_radial_logarithmic_cover_finite_images import (
    source_controls,
    wick,
)


if not __debug__:
    raise RuntimeError("Exact assertions must be enabled")


def phase(leading_coefficient, linear_coefficient, slope, y, r):
    return -(linear_coefficient - r) ** 2 / (4 * r * leading_coefficient) - (
        y - slope / r
    ) ** 2 / 2


def scaled_faces():
    y, u = sp.symbols("y u", nonzero=True)
    r = sp.Symbol("r", positive=True)
    b, a = sp.symbols("b a", nonzero=True)
    expected_top_functions = [
        (4, 1, sp.Rational(1, 2), -1 / (4 * a * u**2) - b**2 / 2),
        (4, 2, sp.Rational(3, 4), -1 / (4 * a * u**4) - b**2 / 2),
        (6, 2, sp.Rational(1, 2), -(u**2 - 1) ** 2 / (4 * a * u**6) - b**2 / 2),
    ]
    for m, n, y_scale, expected in expected_top_functions:
        scaled = r**2 * phase(a * y**m, y**n, b, y, r).subs(y, r**y_scale * u)
        assert sp.cancel(sp.limit(scaled, r, 0, dir="+") - expected) == 0

    # For H=Y the h-coordinate is exact, centered at the moving root Y=r.
    moving = r**2 * phase(a * y**5, y, b, y, r).subs(y, r + r**2 * u)
    assert sp.cancel(sp.limit(moving, r, 0, dir="+") + u**2 / (4 * a) + b**2 / 2) == 0

    z = sp.Symbol("z", nonzero=True)
    resonant = -(z - 1) ** 2 / (4 * a * z**3) - b**2 / 2
    assert sp.cancel(sp.diff(resonant, z) + (z - 1) * (3 - z) / (4 * a * z**4)) == 0
    degenerate = resonant.subs(a, -2 / (27 * b**2))
    assert sp.cancel(degenerate + b**2 * (z - 3) ** 2 * (4 * z - 3) / (8 * z**3)) == 0
    assert sp.diff(degenerate, z, 2).subs(z, 3) == -b**2 / 12
    return {
        "nonresonant_and_resonant_faces": len(expected_top_functions),
        "moving_root_coordinate": True,
        "genuine_double_resonant_root": 3,
        "double_root_second_derivative": "-b^2/12",
    }


def volume_powers():
    m = sp.Symbol("m", positive=True)
    n = sp.Symbol("n", positive=True)
    epsilon = sp.Symbol("epsilon", positive=True)
    volume_order = lambda scale: sp.Rational(1, 2) + (1 - m / 2) * scale

    annulus = volume_order((1 + epsilon) / m)
    assert sp.cancel(annulus - (1 / m - epsilon * (sp.Rational(1, 2) - 1 / m))) == 0
    inverse_h_weight = sp.Rational(1, 2) - m / (2 * n) - (1 - 1 / n)
    h_radius = (1 + m / n - epsilon) / 2
    assert sp.cancel(inverse_h_weight + h_radius - (1 / n - epsilon / 2)) == 0
    assert sp.cancel(volume_order(1 / m) - 1 / m) == 0
    assert sp.cancel(volume_order(1 / (m - 2 * n)) - (1 - n) / (m - 2 * n)) == 0
    assert sp.cancel(volume_order(3 / m) - (3 / m - 1)) == 0
    assert sp.cancel(volume_order(1 / n).subs(m, 3 * n) - (1 / n - 1)) == 0
    assert sp.cancel(inverse_h_weight + (m / n - 1) / 2 - (1 / n - 1)) == 0

    # Substituting z=(-C/T)^(1/m) in sqrt(r) z^-m/2 dz.
    z_exponent = 1 - m / 2
    c_exponent = z_exponent / m
    t_exponent = -1 - z_exponent / m
    assert sp.cancel(t_exponent + sp.Rational(1, 2) + 1 / m) == 0
    assert sp.cancel(sp.Rational(1, 2) - c_exponent - (1 - 1 / m)) == 0
    assert sp.cancel(sp.Rational(1, 2) + c_exponent - 1 / m) == 0
    return {
        "zero_slope_truncated_orders": ["1/m-epsilon*(1/2-1/m)", "1/n-epsilon/2"],
        "nonzero_slope_resonant_order_after_shrink": "1/n-epsilon/2",
        "descending_phase_power": "-1/2-1/m",
        "descending_tail_scalar_orders": ["1-1/m", "1/m"],
    }


def original_wick_checks():
    s, y = sp.symbols("s y")
    checks = 0
    fixtures = [
        (-2 * y**3 / 27, y, y),
        (y**5, y, y),
        ((y - 1) ** 4 * (y + 1) ** 3, (y - 1) ** 2, sp.Integer(0)),
        ((y - 1) ** 6, sp.Integer(0), y + 3),
    ]
    for a, h, boundary in fixtures:
        p = a * s**2 + h * s + boundary
        powers = [sp.expand(p**degree) for degree in range(4)]
        derivative = 2 * a * s + h
        for marker in (sp.Integer(1), y, y**2, s, s * y):
            for degree, power in enumerate(powers):
                # Multiply the Taylor coefficient of the Gamma endpoint identity by degree!.
                amplitude = (marker - sp.diff(marker, s)) * power
                if degree:
                    amplitude -= degree * derivative * marker * powers[degree - 1]
                expected = wick(marker.subs(s, 0) * boundary**degree)
                assert wick(amplitude - expected) == 0
                checks += 1
    return {"phases": len(fixtures), "orders": [0, 3], "original_endpoint_controls": checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[2]
    files = [
        here / name
        for name in (
            Path(__file__).name,
            "GAUSSIAN_QUADRATIC_RADIAL_UNRESTRICTED_FINITE_IMAGES.md",
            "GAUSSIAN_QUADRATIC_RADIAL_LOGARITHMIC_COVER_FINITE_IMAGES.md",
            "verify_gaussian_quadratic_radial_logarithmic_cover_finite_images.py",
            "GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md",
            "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
            "GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md",
            "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
            "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
            "GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md",
            "sources/sabbah_stokes_0912.2762v5.pdf",
        )
    ]
    primary = root / "results/gaussian-boundary"
    manifest = primary / "sources/laurent_primary_manifest.json"
    files += [manifest, primary / "analytic_core.tex"]
    for source in json.loads(manifest.read_text(encoding="utf-8"))["sources"]:
        path = primary / source["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"]
        files.append(path)
    result = {
        "status": "passed",
        "sympy": sp.__version__,
        "source": source_controls(),
        "scaled_faces": scaled_faces(),
        "volume_orders": volume_powers(),
        "original": original_wick_checks(),
        "sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "scope": "Exact finite controls of local identities. The proof supplies the uniform chain estimates and the Stokes filtration argument.",
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    output = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
