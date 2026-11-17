"""Eight exact transcription controls; the all-order proof is in the note."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as s

x, y, z, w = s.symbols('x y z w')
a, b, c, d = s.symbols('a b c d')
entries = (x, y, z, w)
coords = (a, b, c, d)
delta = x*w-y*z
quaternion = {x:a+s.I*b, y:c+s.I*d, z:-c+s.I*d, w:a-s.I*b}


def haar(poly):
    result = 0
    for (ax, ay, az, aw), coefficient in s.Poly(s.expand(poly), *entries).terms():
        if ax == aw and ay == az:
            result += coefficient*s.Rational((-1)**ay*factorial(ax)*factorial(ay), factorial(ax+ay+1))
    return s.simplify(result)


def gaussian(poly):
    result = 0
    for powers, coefficient in s.Poly(s.expand(poly.subs(quaternion)), *coords).terms():
        if all(k % 2 == 0 for k in powers):
            moment = s.Integer(1)
            for k in powers:
                moment *= factorial(k)//(2**(k//2)*factorial(k//2))
            result += coefficient*moment
    return s.simplify(result)


def lift(poly, degree):
    return s.expand(sum(coef*x**i*y**j*z**k*w**l*delta**((degree-i-j-k-l)//2)
                        for (i,j,k,l), coef in s.Poly(poly, *entries).terms()))


def main():
    assert s.expand(delta.subs(quaternion)-(a*a+b*b+c*c+d*d)) == 0
    cases = [
        (1, 1, 0, 0, 3),
        (x*w, 1, 2, 0, 2),
        (y*z, x*w, 2, 2, 3),
        (x*x, 1, 2, 0, 3),
        (1+x*w, 1+y*z, 2, 2, 2),
        (1+x*w+x*x*w*w, 1+y*z, 4, 2, 2),
        (x*w-y*z, x*w+y*z, 2, 2, 4),
        ((x+y)*(z+w), x*w, 2, 2, 3),
    ]
    results = []
    for phase, marker, D, E, m in cases:
        phase, marker = s.sympify(phase), s.sympify(marker)
        pt, qt = lift(phase, D), lift(marker, E)
        radial = 2**((D*m+E)//2)*factorial((D*m+E)//2+1)
        left = gaussian(qt*pt**m)
        right = radial*haar(marker*phase**m)
        assert s.expand(left-right) == 0
        results.append(dict(phase=str(phase), marker=str(marker), D=D, E=E, m=m,
                            gaussian=str(left), haar=str(haar(marker*phase**m)), radial_factor=radial))
    here = Path(__file__).resolve().parent
    receipt = dict(passed=True, exact_identities=len(results), controls=results,
                   proof_scope='Finite even-sector transcription controls only; all-order radial and exceptional Reynolds proofs are written in the canonical note.',
                   user_reported_controls='Eight controls reported; original inputs not supplied. These are separate reproducible controls.',
                   checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   independent_replay=False, formal_verification=False, external_review=False)
    (here/'gaussian_radial_transfer_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('PASS: 8 exact radial identities; determinant/quaternion normalization verified.')


if __name__ == '__main__':
    main()
