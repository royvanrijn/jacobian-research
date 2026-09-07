#!/usr/bin/env sage-python
"""Distinguish the two MW14 triangle j-maps by their branch divisors.

At1013, interpolate the binary discriminant in the target j-coordinate.
Its degree is at most46. A nonaffine PGL2 control preserves the normalized
polynomial. Unequal branch divisors exclude equivalence of the two j-maps.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
from sage.all import GF, PolynomialRing

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
INPUTS = [ART/(name+'-mod1013-v1.json') for name in
          ['elkies-k3-curve302-triangle-mw14','elkies-k3-curve302-triangle-6-8-mw14']]
OUT = ART/'elkies-k3-curve302-triangle-branch-separation-v1.json'


def binary_discriminant(f, degree=24):
    if f.degree() == degree:
        return f.discriminant()
    if f.degree() == degree-1:
        return f.leading_coefficient()**2*f.discriminant()
    return f.base_ring()(0)


def branch_polynomial(n,d):
    R = n.parent(); J = R.gen()
    values = [R.base_ring()(i) for i in range(50)]
    if d.degree() == 24:
        leading_drop = n[24]/d[24]
        if leading_drop not in values:
            values[-1] = leading_drop
    samples = [(a,binary_discriminant(n-a*d)) for a in values]
    samples.sort(key=lambda row:int(row[0]))
    h = R.lagrange_polynomial(samples[:47])
    assert h and h.degree() <= 46
    assert all(h(a) == value for a,value in samples)
    return h.monic(), [[int(a),int(b)] for a,b in samples]


def build():
    field = GF(1013); R = PolynomialRing(field,'s'); s = R.gen()
    records = []; polynomials = []
    for path in INPUTS:
        data = json.loads(path.read_text())
        n = R(data['j_numerator_coefficients_low_to_high'])
        d = R(data['j_denominator_coefficients_low_to_high'])
        assert max(n.degree(),d.degree()) == 24 and n.gcd(d) == 1
        assert all(n(a) == field(j)*d(a) and d(a) for a,j in data['samples'])
        h,samples = branch_polynomial(n,d)
        residual, rem = h.quo_rem(s**16*(s-field(1728))**12)
        assert not rem and h.degree() == 43 and residual.degree() == 15
        polynomials.append(h)
        records.append({'input':str(path.relative_to(ROOT)), 'binary_discriminant_samples':samples,
                        'normalized_branch_polynomial':list(map(int,h.list())),
                        'residual_branch_polynomial':list(map(int,residual.list())),
                        'degree':int(h.degree())})
    # A nonaffine base change with determinant1. Its polynomial binary
    # discriminant must be identical, including the exceptional leading-term drop.
    old = json.loads(INPUTS[0].read_text())
    n = R(old['j_numerator_coefficients_low_to_high']); d = R(old['j_denominator_coefficients_low_to_high'])
    def change(f):
        return sum(f[i]*(2*s+3)**i*(s+2)**(24-i) for i in range(25))
    changed_n,changed_d = change(n),change(d)
    assert changed_n.gcd(changed_d) == 1
    control, control_samples = branch_polynomial(changed_n,changed_d)
    assert control == polynomials[0]
    assert polynomials[0] != polynomials[1]
    return {'schema':'curve302.triangle-branch-separation.v1','status':'TWO_TRIANGLE_J_MAPS_HAVE_DIFFERENT_BRANCH_DIVISORS',
            'prime':1013,'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__)]+INPUTS},
            'records':records,'interpolation_samples':47,'additional_checks':3,
            'PGL2_control':{'matrix':[[2,3],[1,2]],'determinant':1,'normalized_branch_polynomial_preserved':True,
                             'binary_discriminant_samples':control_samples},
            'conclusion':'With the certified degree-preserving reductions of both characteristic-zero j-maps, their unequal normalized branch discriminants exclude PGL2 equivalence over Qbar, hence also Q. These are different elliptic fibrations on the same K3, not different surface isomorphism classes or302 parents.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args = parser.parse_args()
    signal.alarm(120)
    result = build()
    if args.check:
        assert result == json.loads(OUT.read_text())
    else:
        assert not OUT.exists(), 'Preserve certificate'
        OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)
