#!/usr/bin/env sage-python
"""One uniform inherited RR point; algebra only, with a 25-second limit.

Only generic RR coefficients and the already certified universal equation
are used. No exceptional point, control outcome, or search coordinate enters.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, version

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_rr_inherited_point_v1.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    net_path = ART / 'det1092_first_centre_rr_net_v1.json'
    universal_path = ART / 'det1092_universal_rr_descent_preflight_v1.json'
    net = json.loads(net_path.read_text())
    data = json.loads(universal_path.read_text())['universal_genus2']
    R = PolynomialRing(QQ, 'T')
    A = [R(row) for row in net['A']]
    B = [R(row) for row in net['B']]
    h = R(data['h'])
    assert h.degree() == 3
    a, ra = A[2].quo_rem(h)
    b, rb = B[2].quo_rem(h)
    assert not ra and not rb and a.degree() == 0 and b.degree() == 1
    alpha, beta0, beta1 = a[0], b[0], b[1]
    UV = PolynomialRing(QQ, ['u', 'v'])
    u, v = UV.gens()
    t0, z0 = -(beta0 + alpha*u), beta1 + alpha*v

    def hom_eval(poly, degree):
        assert poly.degree() <= degree
        return sum(UV(c)*t0**i*z0**(degree-i)
                   for i, c in enumerate(poly.list()))

    H = hom_eval(h, 3)
    # f1 has degree <= 6; u+vT becomes (uZ+vT)/Z.
    assert A[1].degree() <= 5 and B[1].degree() <= 6
    F1 = hom_eval(B[1], 6) + (u*z0+v*t0)*hom_eval(A[1], 5)
    G, remainder = F1.quo_rem(H**2)
    assert not remainder
    W = H*G**2
    Q = sum(QQ(row['coefficient'])*u**row['u']*v**row['v'] *
            t0**row['T']*z0**(6-row['T']) for row in data['sparse_q'])
    scale = QQ(data['scale'])
    assert W**2 == scale*Q
    assert z0(0, 0) != 0 and W(0, 0) != 0
    fixed_T = t0(0, 0)/z0(0, 0)
    fixed_s = W(0, 0)/z0(0, 0)**3
    fixed_q = R([sum(QQ(row['coefficient']) for row in data['sparse_q']
                     if row['T'] == i and row['u'] == row['v'] == 0)
                 for i in range(7)])
    assert fixed_q.degree() == 6
    assert fixed_q.gcd(fixed_q.derivative()) == 1
    assert fixed_s**2 == scale*fixed_q(fixed_T)
    exceptional = [-beta0/alpha, -beta1/alpha]
    assert B[2]+(exceptional[0]+exceptional[1]*R.gen())*A[2] == 0

    def sparse(poly):
        return [{'u': int(i), 'v': int(j), 'coefficient': str(c)}
                for (i, j), c in sorted(UV(poly).dict().items())]

    return {
        'classification': 'verified application and new deduction',
        'status': 'PASS_UNIFORM_INHERITED_RATIONAL_RR_POINT',
        'equation': 's^2=scale*q(T;u,v), from the immutable universal certificate',
        'point_formula': {
            'weights': [1, 1, 3],
            'coordinates': '[T:Z:s]=[t0:z0:+/-W]',
            'alpha': str(alpha), 'beta0': str(beta0), 'beta1': str(beta1),
            't0': sparse(t0), 'z0': sparse(z0),
            'H': sparse(H), 'G': sparse(G), 'W': sparse(W),
            'identity': 'F1(t0,z0)=H^2*G; W=H*G^2; W^2=scale*q_hom(t0,z0;u,v)',
            'unique_base_parameter': list(map(str, exceptional)),
            'base_parameter_member': 'reducible O+P_w, excluded from the smooth integral locus',
            'fixed_0_0_point': [str(fixed_T), str(fixed_s)],
            'degrees_H_G_W': [int(H.total_degree()), int(G.total_degree()), int(W.total_degree())]
        },
        'divisor_argument': {
            'D': '2O+5F+phi(w)=O+P_w', 'O_squared': -2,
            'O_dot_F': 1, 'O_dot_phi_w': 0, 'D_dot_O': 1,
            'conclusion': 'Every smooth geometrically integral Q-rational RR member intersects O in a rational point.'
        },
        'arithmetic_consequences': {
            'rationally_soluble_on_smooth_integral_locus': True,
            'everywhere_locally_soluble_on_smooth_integral_locus': True,
            'curve_fake_2_Selmer_set_nonempty': True,
            'full_Jacobian_2_Selmer_group': 'NOT_COMPUTED',
            'Cassels_Tate_pairing': 'NOT_COMPUTED',
            'incidence_discriminator': 'Raw RR-curve solubility or curve-Selmer-set nonemptiness cannot distinguish 302 from any smooth rational control member.',
            'new_MW_point': False,
            'scope': 'Inherited point maps to O or P_w at its own parameter, not to a new point at an arbitrarily marked elliptic fibre.'
        },
        'limits': {'wall_seconds': 25, 'uniform_point_identities': 1,
                   'point_searches': 0, 'Selmer_runs': 0, 'class_group_runs': 0,
                   'control_sweeps': 0, 'pilot_changes': 0},
        'software': version(),
        'inputs': {str(p.relative_to(ROOT)): sha(p)
                   for p in [net_path, universal_path, Path(__file__)]}
    }


if __name__ == '__main__':
    signal.alarm(25)
    result = build()
    payload = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if OUT.exists():
        assert OUT.read_text() == payload
    else:
        OUT.write_text(payload)
    print(result['status'], result['point_formula']['degrees_H_G_W'], flush=True)
