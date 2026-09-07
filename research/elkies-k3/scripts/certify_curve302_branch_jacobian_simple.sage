#!/usr/bin/env sage-python
"""Small Rabin certificate for Q-simplicity of the branch-carrier Jacobian.

The Frobenius polynomial comes from the separately replayed finite-field
point counts. Irreducibility modulo3 at the good prime47 proves Q-simplicity;
no absolute simplicity or Mordell-Weil rank bound is asserted.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ, gcd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'artifacts/generated-results/elkies-k3-curve302-branch-carrier-probe-v1.json'
OUT = ROOT/'artifacts/generated-results/elkies-k3-curve302-branch-jacobian-simple-v1.json'


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def compute():
    d = json.loads(SOURCE.read_text())
    for path, expected in d['input_sha256'].items():
        assert digest(ROOT/path) == expected
    rows = d['jacobian_arithmetic']['records']
    assert [r['prime'] for r in rows] == [47, 53]
    zr = PolynomialRing(ZZ, 'T')
    orders = []
    for row in rows:
        p = row['prime']
        polynomial = zr(row['frobenius_low_to_high'])
        a1, a2, a3 = [ZZ(p**n+1-row['extension_point_counts_n1_n2_n3'][n-1]) for n in [1, 2, 3]]
        c1, c2 = -a1, (a1*a1-a2)//2
        assert (a1*a1-a2) % 2 == 0
        assert (a3+c1*a2+c2*a1) % 3 == 0
        c3 = -(a3+c1*a2+c2*a1)//3
        assert polynomial == zr([p**3, p**2*c1, p*c2, c3, c2, c1, 1])
        assert polynomial(1) == ZZ(row['jacobian_order'])
        orders.append(polynomial(1))
    assert gcd(orders) == 64
    ring = PolynomialRing(GF(3), 'T')
    T = ring.gen()
    f = ring(rows[0]['frobenius_low_to_high'])
    assert f == T**6+T**5+T+2
    quotient = ring.quotient(f, 't')
    t = quotient.gen()
    assert t**(3**6) == t
    witnesses = []
    for degree in [2, 3]:
        remainder = (t**(3**degree)-t).lift()
        common, u, v = f.xgcd(remainder)
        assert common == 1 and u*f+v*remainder == 1
        witnesses.append({'exponent_degree': degree,
                          'remainder_low_to_high': list(map(int, remainder.list())),
                          'bezout_for_f_low_to_high': list(map(int, u.list())),
                          'bezout_for_remainder_low_to_high': list(map(int, v.list()))})
    return {'schema': 'curve302.branch-jacobian-simple.v1',
            'status': 'Q_SIMPLE_JACOBIAN_WITH_RATIONAL_TORSION_Z2_POWER6',
            'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), SOURCE]},
            'good_reduction_prime': 47, 'irreducibility_prime': 3,
            'mod3_polynomial_low_to_high': list(map(int, f.list())),
            'T_power729_minus_T_remainder': [], 'rabin_witnesses': witnesses,
            'jacobian_orders': list(map(str, orders)), 'order_gcd': 64,
            'deduction': 'Rabin proves the degree-six Frobenius polynomial irreducible mod3 and hence overQ. A proper positive-dimensional Q-abelian subvariety would factor that good-reduction Frobenius polynomial, a contradiction. Thus the genus-three Jacobian is Q-simple. Its64 rational two-torsion points exhaust rational torsion by the two good-reduction orders.',
            'conditional_rank_zero_consequence': 'If rank J(Q)=0, every rational point on the genus-three curve is Weierstrass: [P-W] then has order dividing2, so2P belongs to the unique hyperelliptic degree-two linear system.',
            'boundary': 'No rank-zero proof, global rational-point determination, absolute simplicity proof or parent construction. Point counts are supplied by the probe and its direct enumeration replay.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    text = json.dumps(compute(), indent=2, sort_keys=True)+'\n'
    if args.check:
        assert OUT.read_text() == text
    else:
        assert not OUT.exists(), 'Use --check; preserve the certificate'
        OUT.write_text(text)
    print('PASS_RABIN_Q_SIMPLE_AND_TORSION', flush=True)
