#!/usr/bin/env sage-python
"""Independent identities and ramification proof; no producer import.

Manual tangent/chord arithmetic replaces the producer's elliptic group law.
The full original-coordinate composition is verified by generic identities;
finite-field rational functions certify that the chosen curve is not in the
exceptional locus of the birational map. Run each multiple under timeout25s.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from sage.all import QQ, GF, PolynomialRing, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_pencil_multiples_v2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(n, write=False):
    started = time.monotonic()
    protocol = json.loads((OUT/'protocol.json').read_text())
    for path, expected in protocol['inputs'].items():
        assert sha(ROOT/path) == expected
    assert protocol['multiples'] == [2, 3]
    data = json.loads((ART/'det1092_norm8_seed_cover_v2/generic.json').read_text())
    parent = json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
    result = json.loads((OUT/('multiple%d.json'%n)).read_text())
    assert result['inputs'] == protocol['inputs']
    R = PolynomialRing(QQ, 'z'); z = R.gen(); K = R.fraction_field()
    def val(record):
        return K(R(record['numerator']))/R(record['denominator'])
    old = EllipticCurve(K, [val(row) for row in parent['a_invariants']])
    a, b = -old.c4()/48, -old.c6()/864
    h, nx, ny, shift = [R(data[name]) for name in ['pole_h', 'nx', 'ny', 'shift']]
    assert ny*ny == nx**3+a*nx*h**4+b*h**6
    # Reconstruct the generic two-variable pencil identity before specializing.
    L = PolynomialRing(R, 'v'); v = L.gen()
    M = L(h*h)*v-L(shift)
    left = sum((L([QQ(c)*z**i for c in row])
                for i,row in enumerate(data['quartic_t_coefficients_in_z'])), L.zero())
    assert h**6*left == M**4-6*nx*M*M-8*ny*M-3*nx**2-4*L(R(a))*h**4
    # The inverse map is universal algebra, checked without large compositions.
    A = PolynomialRing(QQ, ['c', 'd', 'm', 'h', 'a'])
    c, d, mm, hh, alpha = A.gens()
    V = PolynomialRing(A.fraction_field(), 'W'); ww = V.gen()
    XX = (hh*ww-c+mm*mm)/2
    YY = mm*(XX-c)-d
    beta = d*d-c**3-alpha*c
    relation = hh*hh*ww*ww-(mm**4-6*c*mm*mm-8*d*mm-3*c*c-4*alpha)
    assert (YY*YY-XX**3-alpha*XX-beta)%relation == 0
    # Manual group law on the pointed quartic's monic cubic.
    S = PolynomialRing(K, 't'); t = S.gen()
    f = S([R(row) for row in data['quartic_t_coefficients_in_z']])
    t0, s = [val(data['sections'][0][key]) for key in ['t_of_z','W_of_z']]
    q0,q1,q2,q3,q4 = f(t+t0).list()
    assert q0 == s*s
    aa4 = q1*q3-4*s*s*q4
    aa6 = s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4
    bx = q1*q1/(4*s*s)-q2
    by = -(q1*bx+2*s*s*q3)/(2*s)
    assert by*by == bx**3+q2*bx*bx+aa4*bx+aa6
    slope = (3*bx*bx+2*q2*bx+aa4)/(2*by)
    xx = slope*slope-q2-2*bx
    yy = -by+slope*(bx-xx)
    if n==3:
        slope = (yy-by)/(xx-bx)
        newx = slope*slope-q2-bx-xx
        yy = -by+slope*(bx-newx); xx = newx
    assert yy*yy == xx**3+q2*xx*xx+aa4*xx+aa6
    T, W = val(result['t_of_z']), val(result['W_of_z'])
    u = T-t0
    assert u == (2*s*yy+q1*xx+2*s*s*q3)/(xx*xx-4*s*s*q4)
    assert W == xx*u*u/(2*s)-s-q1*u/(2*s)
    assert W*W == f(T)
    print('n', n, 'exact generic maps replayed', time.monotonic()-started, flush=True)
    N,D = T.numerator(), T.denominator()
    deg = max(N.degree(),D.degree())
    assert deg == result['parameter_map_degree'] and deg>1
    delta = old.discriminant(); assert delta.denominator().degree()==0
    delta = delta.numerator()
    assert delta.degree()==24 and delta.gcd(delta.derivative()).degree()==0
    assert max(old.j_invariant().numerator().degree(),old.j_invariant().denominator().degree())==24
    ram = N.derivative()*D-N*D.derivative()
    support = R(result['ramification_support'])
    assert support == (ram//ram.gcd(ram.derivative())).monic()
    assert support.degree()==2*deg-2
    cert = result['ramification_certificate']; p = cert['prime']
    Rp = PolynomialRing(GF(p), 'z'); kp = Rp.fraction_field(); zp = Rp.gen()
    nn,dd,rr,ds = map(Rp,[N,D,support,delta])
    assert all(x.degree()==y.degree() for x,y in [(nn,N),(dd,D),(rr,support),(ds,delta)])
    assert nn.gcd(dd)==1 and rr.gcd(rr.derivative())==1
    # Horner homogeneous substitution differs from the producer's power sum.
    bad = Rp.zero()
    for i in range(ds.degree(),-1,-1):
        bad = bad*nn+ds[i]*dd**(ds.degree()-i)
    bad *= dd
    g, v1, v2 = rr.xgcd(bad)
    assert v1*rr+v2*bad == 1 and g==1
    assert list(map(int,rr.list())) == cert['ramification_support_mod_p']
    assert cert['bad_support_gcd_mod_p']==[1]
    assert cert['smooth_ramification_support_lower_bound']==int(support.degree())
    # Certify nonzero denominators and recovery of the parameter on the
    # selected rational curve. A nonzero reduction implies nonzero over Q.
    tp = kp(nn)/dd
    wp = kp(Rp(W.numerator()))/Rp(W.denominator())
    hp = Rp(h)(tp)
    cp,dp = Rp(nx)(tp)/(hp*hp),Rp(ny)(tp)/(hp**3)
    sp = Rp(shift)(tp); mp = hp*zp-sp/hp
    xp = (hp*wp-cp+mp*mp)/2
    yp = mp*(xp-cp)-dp
    assert hp and xp-cp
    assert (hp*(yp+dp)/(xp-cp)+sp)/(hp*hp)==zp
    report = {
        'classification': 'independent verified application and new construction',
        'status': 'PASS_INDEPENDENT_RANK18_RATIONAL_BASE_CHANGE',
        'n': n, 'degree': int(deg), 'ground_field': 'Q',
        'smooth_ramification_points': int(support.degree()), 'certificate_prime': p,
        'generic_rank_lower_bound': 18, 'new_independent_directions': 1,
        'checks': ['generic quartic identity', 'universal elliptic map identity',
                   'manual tangent/chord multiple', 'exact stored quartic point',
                   'original discriminant degree24 squarefree and j degree24',
                   'full-degree squarefree ramification', 'modular Bezout gcd1',
                   'nonexceptional image and birational parameter recovery'],
        'seconds': time.monotonic()-started,
        'checker_sha256': sha(Path(__file__)),
        'construction_sha256': sha(OUT/('multiple%d.json'%n)),
        'boundary': 'Uses the certified original full generic MW17. Does not prove generic rank exactly18, a specialized rank, rank19 on a common cover, or any302 incidence.',
    }
    if write:
        with (OUT/('replay%d.json'%n)).open('x') as stream:
            json.dump(report,stream,indent=2,sort_keys=True); stream.write('\n')
    print(json.dumps(report,sort_keys=True),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--multiple',type=int,choices=[2,3],required=True)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); verify(args.multiple,args.write)
