"""Bounded exact original-Wick controls for invariant-core boundary images.

The unrestricted endpoint and finite-module statements are written proofs.
No finite packet is used to infer finiteness or vacuum primitivity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O so exact assertions execute")

S, S2, Y, t = sp.symbols("S S2 Y t")
Z, W, Z2, W2 = sp.symbols("Z W Z2 W2")


def zero(value):
    assert sp.expand(value) == 0, value


def ordinary(power):
    return 0 if power % 2 else factorial(power)//(2**(power//2)*factorial(power//2))


def wick(polynomial, two_pairs=False):
    polynomial = sp.expand(polynomial.subs({S: Z*W, S2: Z2*W2}))
    variables = (Z, W, Z2, W2, Y) if two_pairs else (Z, W, Y)
    answer = sp.Integer(0)
    for powers, coefficient in sp.Poly(polynomial, *variables).terms():
        if powers[0] != powers[1]:
            continue
        factor = factorial(powers[0])
        if two_pairs:
            if powers[2] != powers[3]:
                continue
            factor *= factorial(powers[2])
        answer += coefficient*factor*ordinary(powers[-1])
    return sp.expand(answer)


def boundary_operator(phase, amplitude, variable=S):
    return sp.expand(amplitude-sp.diff(amplitude, variable)
                     -t*sp.diff(phase, variable)*amplitude)


def raw_time_amplitude(phase, amplitude, order, two_pairs=False):
    result = sp.Integer(0)
    for (power,), coefficient in sp.Poly(sp.expand(amplitude), t).terms():
        if power <= order:
            result += sp.Rational(factorial(order), factorial(order-power))*wick(
                coefficient*phase**(order-power), two_pairs)
    return sp.expand(result)


def spatial_degree(polynomial):
    value = sp.Poly(sp.expand(polynomial.subs(S, Z*W)), Z, W, Y, t)
    return max((a+b+c for (a,b,c,power), coefficient in value.terms()), default=0)


def checks():
    phases = [S**2*Y**3+S*Y**2+S**3+Y,
              S**2*(Y**3+Y**2)+S*(Y+1)+2*Y**2-3*Y+1,
              S*(Y**2+Y)+Y**3+Y**2-2*Y,
              3+2*Y]
    markers = [sp.Integer(1), S, Y, S*Y, S**2, Y**2, S**2*Y**2]
    controls = 0
    for phase in phases:
        boundary = phase.subs(S, 0)
        for marker in markers:
            amplitude = boundary_operator(phase, marker)
            degree_bound = spatial_degree(marker)+max(0, spatial_degree(phase)-2)
            assert spatial_degree(amplitude) <= degree_bound
            assert sp.Poly(amplitude, t).degree() <= 1
            for n in range(5):
                zero(raw_time_amplitude(phase, amplitude, n)
                     -wick(marker.subs(S, 0)*boundary**n))
                controls += 1
        amplitude = S**2*Y+t*S*Y**2+t**2*Y
        left = sp.diff(boundary_operator(phase, amplitude), t)+phase*boundary_operator(phase, amplitude)
        right = boundary_operator(phase, sp.diff(amplitude, t)+phase*amplitude)
        zero(left-right)
    remainder = Z*(1+Y**2)+Z**2*W*Y
    phase = phases[0]
    strict_controls = 0
    for marker in markers[:5]:
        amplitude = boundary_operator(phase, marker)
        for n in range(5):
            zero(raw_time_amplitude(phase+remainder, amplitude, n)
                 -wick(marker.subs(S, 0)*phase.subs(S, 0)**n))
            strict_controls += 1
    # The zero boundary need not imply a finite image. Its exact periods
    # have unbounded pole order; bounded coefficients check the formula.
    for marker_degree in range(5):
        for n in range(7):
            expected = factorial(n+marker_degree)
            zero(wick(S**(n+marker_degree))-expected)
            zero(raw_time_amplitude(S, 1-t, n)-(1 if n == 0 else 0))
    # An independent affine boundary has its primitive generator exactly.
    affine = 2*Y+3
    for n in range(7):
        actual = wick(affine**n)/factorial(n)
        expected = sp.expand(sp.series(sp.exp(3*t+2*t**2), t, 0, 7).removeO()).coeff(t, n)
        zero(actual-expected)
    # Multiple radial endpoints: all cross derivatives and time factors
    # are contracted in the literal original two-Witt-pair source.
    phase2 = S*S2*Y**3+S**2*Y+S2*Y**2+Y
    multi_controls = 0
    for marker in (sp.Integer(1), Y, Y**2, S+S2*Y):
        first = boundary_operator(phase2, boundary_operator(phase2, marker, S2), S)
        second = boundary_operator(phase2, boundary_operator(phase2, marker, S), S2)
        zero(first-second)
        assert sp.Poly(first, t).degree() <= 2
        for n in range(4):
            expected = wick(marker.subs({S:0,S2:0})*Y**n, True)
            zero(raw_time_amplitude(phase2, first, n, True)-expected)
            multi_controls += 1
    return dict(status="passed", coefficient_field="QQ",
                one_radial_phases=len(phases), original_marker_count=len(markers),
                one_radial_time_orders=[0,4], exact_boundary_controls=controls,
                strict_remainder_controls=strict_controls,
                parameter_commutators=len(phases),
                two_radial_controls=multi_controls, two_radial_time_orders=[0,3],
                scope="Finite exact original-Wick and operator controls; unrestricted identities and finite-image consequences are the written proofs.",
                independent_replay=False, formal_verification=False, external_review=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    imports = ["GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md",
               "GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md",
               "GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md",
               "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
               "GAUSSIAN_W_AFFINE_ARBITRARY_FIBRE_FINITE_IMAGES.md",
               "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md"]
    pins = {name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
            for name in imports}
    pins[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result.update(python=platform.python_version(), sympy=sp.__version__, sha256=pins)
    output = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
