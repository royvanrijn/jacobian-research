#!/usr/bin/env python3
"""Exact finite controls for the written conformal-Euler classification."""
from pathlib import Path
import hashlib
import json
import math

from conformal_euler_cubic import (
    S, z, w, t, q, W, V, Z, Y, harmonic_basis, harmonic_parts,
    laplace, fields, conformal_certificate, sphere_average,
    gaussian_average, shear, hopf_matrix, biform_to_phase,
)


def main():
    root = Path(__file__).resolve().parent
    a, b, c, u, v = S.symbols('a b c u v')
    basis = harmonic_basis()
    for degree, level in basis.items():
        assert len(level) == 2*degree+1
        for charge, f in level.items():
            assert laplace(f) == 0
            assert S.expand(fields(f)[1]-charge*f) == 0
    assert fields(q) == (2*q, 0, 0, 0)

    matrix, source_basis, target_basis = hopf_matrix()
    determinant = matrix.det()
    assert matrix.shape == (16, 16) and determinant != 0
    coefficients = S.symbols('p0:16')
    generic = sum(coefficient*monomial for coefficient, monomial in zip(coefficients, target_basis))
    phase = biform_to_phase(generic)
    parts = harmonic_parts(phase)
    s = W*Z+V*Y
    hopf = {z: W*Y, w: 2*V*Z, t: W*Z-V*Y}
    lifted = sum(s**(3-j)*f.subs(hopf, simultaneous=True) for j, f in parts.items())
    assert S.expand(lifted-generic) == 0
    assert all(laplace(f) == 0 for f in parts.values())

    patterns = []
    for k1 in range(-1, 2):
        for k2 in range(-2, 3):
            for k3 in range(-3, 4):
                if k2-k1 > 0 and k3-k2 == k2-k1:
                    patterns.append((k1, k2, k3))
    assert patterns == [(-1, 0, 1), (-1, 1, 3), (0, 1, 2), (1, 2, 3)]
    forbidden1 = c*w+b*(t*t-z*w)+a*z*(t*t-z*w/2)
    forbidden2 = c*t+b*z*t+a*z*z*t
    assert sphere_average(forbidden1**2) == b*b/5
    assert sphere_average(forbidden2**2) == c*c/3

    opposite_charge_checks = []
    for lower, upper in [(1, 2), (2, 3)]:
        for k in range(1, lower+1):
            for l in range(-upper, 0):
                divisor = math.gcd(k, -l)
                r, ss = -l//divisor, k//divisor
                product = S.expand(basis[lower][k]**r*basis[upper][l]**ss)
                average_square = sphere_average(product*product)
                assert average_square > 0
                moment = sphere_average((basis[lower][k]+basis[upper][l])**(2*(r+ss)))
                assert moment == S.binomial(2*(r+ss), 2*r)*average_square
                opposite_charge_checks.append([lower, upper, k, l, str(average_square)])

    shear_family = a*z**3+b*z*t+c*w
    invariant = c*c*(2*a*c+b*b)/2
    assert sphere_average(shear_family**4) == 4*c*c*(2*a*c+b*b)/35
    for n in range(1, 4):
        expected = S.factorial(4*n)*invariant**n/S.factorial(n)
        assert S.expand(gaussian_average(shear_family**(4*n))-expected) == 0
        assert S.expand(sphere_average(shear_family**(4*n))-expected/S.factorial2(6*n+1)) == 0

    assert S.expand(shear(q, u)-q) == 0
    assert S.expand(shear(shear(z+w+t, u), v)-shear(z+w+t, u+v)) == 0
    assert S.expand(shear(w, u)-(w+u*z*t-u*u*z**3/2)) == 0
    derivation_checks = 0
    for degree in range(7):
        for i in range(degree+1):
            for j in range(degree-i+1):
                monomial = z**i*w**j*t**(degree-i-j)
                derivative = z*(t*S.diff(monomial, w)-z*S.diff(monomial, t))
                assert sphere_average(derivative) == 0
                derivation_checks += 1

    f0 = shear(w, S.Integer(1))
    certificate = conformal_certificate(f0)
    assert certificate['recognized']
    assert (certificate['alpha'], certificate['beta'], certificate['gamma'], certificate['delta']) == (-S.Rational(1, 2), 0, 0, S.Rational(3, 2))
    rotated = f0.subs({z: z, t: t+2*z, w: w-2*t-2*z}, simultaneous=True)
    assert conformal_certificate(rotated)['recognized']
    unsupported_safe = z+z*t+z**3
    assert not conformal_certificate(unsupported_safe)['recognized']
    for m in range(1, 9):
        assert sphere_average(f0**m) == sphere_average(unsupported_safe**m) == 0
    nonpure = w+z*t+z**3
    assert conformal_certificate(nonpure)['recognized']
    assert sphere_average(nonpure**4) != 0

    sharp = {}
    for degree in range(1, 5):
        actual = sphere_average(w**degree*f0**(3*degree))
        expected = (-S.Rational(1, 2))**degree*S.factorial(3*degree)/S.factorial2(6*degree+1)
        assert actual == expected and actual != 0
        assert sphere_average(w**degree*f0**(3*degree+1)) == 0
        sharp[str(degree)] = str(actual)

    rodrigues = (W**3*Z**2*Y+W**3*Y**3-W**2*V*Z**3
                 -2*W*V**2*Z**2*Y-W*V**2*Y**3-V**3*Z*Y**2)
    assert S.expand(biform_to_phase(rodrigues)-(z**3+z*t-w/2)) == 0
    assert S.expand(rodrigues-(W*Y)**3-s*(W*Y)*(W*Z-V*Y)+s*s*V*Z) == 0
    assert conformal_certificate(biform_to_phase(rodrigues))['recognized']

    report = {
        'claim': 'SIC2B33EULER20261116', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'adapter_sha256': hashlib.sha256((root/'conformal_euler_cubic.py').read_bytes()).hexdigest(),
        'sympy_version': S.__version__,
        'full_biform_adapter_determinant': str(determinant),
        'three_summand_patterns_up_to_charge_reversal': patterns,
        'opposite_charge_two_summand_checks': opposite_charge_checks,
        'sphere_replica_derivation_checks': derivation_checks,
        'sharp_spherical_marked_values_at_m_equal_3b': sharp,
        'unsupported_safe_control': 'z+z*t+z**3',
        'boundary': 'The supplied conformal-Euler architecture is classified by the written '
                    'all-order proof. No assertion that every semistable pure phase admits '
                    'a grading; the full (3,3) Image problem remains open. No external, '
                    'independent specialist or formal verification.',
    }
    (root/'conformal_euler_cubic_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
