#!/usr/bin/env sage-python
"""Exact generic plane-cubic conversion and completion of ten triangle inverses.

No fibre normalization or interpolation is used here. Over Q(s), compute
one explicit birational plane cubic and its Jacobian; settle the remaining
302 inverse by a primitive polynomial obstruction. Over Fp(s), independently
reconstruct all sixteen saved j-maps. Build180 seconds, replay120 seconds.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, vector, identity_matrix, gcd, lcm, prime_range, EllipticCurve
from sage.env import SAGE_VERSION
from sage.schemes.toric.weierstrass import WeierstrassForm

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
PROTOCOL = ART/'elkies-k3-curve302-shortword-triangles-protocol-v1.json'
SCREEN = ART/'elkies-k3-curve302-shortword-triangles-v1.json'
PUBLIC = ROOT/'elliptic-curves/cas/icarm_curve302.py'
OUT = ART/'elkies-k3-curve302-triangle-generic-conversion-v1.json'


def digest(p): return sha256(p.read_bytes()).hexdigest()
def enc_s(a): return {'n':list(map(str,a.numerator().list())), 'd':list(map(str,a.denominator().list()))}
def enc_u(a): return {'n':[enc_s(c) for c in a.numerator().list()], 'd':[enc_s(c) for c in a.denominator().list()]}
def dec_s(K,a):
    R = K.ring(); return K(R(a['n'])/R(a['d']))
def dec_u(Ku,a):
    R = Ku.ring(); Ks = R.base_ring()
    return Ku(R([dec_s(Ks,c) for c in a['n']])/R([dec_s(Ks,c) for c in a['d']]))


def convert(field,model,provided=None):
    Ks = PolynomialRing(field,'s').fraction_field(); s = Ks.gen()
    R = PolynomialRing(Ks,'u'); u = R.gen(); Ku = R.fraction_field(); X = PolynomialRing(Ku,'x')
    coefficients = []
    for row in model['raw_cubic_coefficients']:
        coefficients.append(sum(Ku(R([field(QQ(c)) for c in a['n']]))/R([field(QQ(c)) for c in a['d']])*s**j for j,a in enumerate(row)))
    f = X(coefficients).monic()
    assert all(c.denominator() == 1 for c in f.list())
    f = f.change_ring(R); A,B,C = f[2],f[1],f[0]
    assert [c.degree() for c in f.list()] == [12,8,4,0]
    D = R(f.discriminant()); q = D.gcd(D.derivative()).monic()
    assert q.degree() == 9
    reduced = D//q**2
    assert D == q**2*reduced and reduced.degree() == 6
    assert reduced.gcd(reduced.derivative()) == 1 and reduced.gcd(q) == 1
    if provided is None:
        r = ((9*C-A*B)*(2*(A*A-3*B)).inverse_mod(q))%q
        a = (A+r)%q; b = (-r*(A+2*r))%q
        P = matrix(R,[[q,0,0],[0,q,0],[b,a,1]])
        shifts = [0,4,8]; steps = 0
        while True:
            info = []
            for row in P.rows():
                deg = max(c.degree()+shifts[i] for i,c in enumerate(row) if c)
                pos = max(i for i,c in enumerate(row) if c and c.degree()+shifts[i] == deg)
                info.append((deg,pos))
            pairs = [(i,j) for i in range(3) for j in range(i) if info[i][1] == info[j][1]]
            if not pairs: break
            i,j = pairs[0]
            if info[i][0]<info[j][0]: i,j = j,i
            pos = info[i][1]; power = info[i][0]-info[j][0]
            P[i] -= P[i,pos].leading_coefficient()/P[j,pos].leading_coefficient()*u**power*P[j]
            steps += 1; assert steps <= 100
        assert P.det() == q**2
        assert sorted(d-q.degree() for d,pos in info) == [0,1,2]
        chosen = next(i for i,(d,pos) in enumerate(info) if d-q.degree() == 1)
        g = [Ku(a/q) for a in P[chosen]]
    else:
        g = [dec_u(Ku,a) for a in provided]
    MX = matrix(Ku,[[0,0,-C],[1,0,-B],[0,1,-A]])
    I = identity_matrix(Ku,3)
    MG = g[0]*I+g[1]*MX+g[2]*MX**2
    plane = MG.charpoly('v')
    assert sum(c*MG**i for i,c in enumerate(plane.list())) == matrix(Ku,3)
    assert all(c.denominator() == 1 for c in plane.list())
    plane = plane.change_ring(R)
    assert [a.degree() for a in plane.list()] == [3,2,1,0]
    # Invert the generator change: 1,v,v^2 span the same cubic algebra.
    one = vector(Ku,[1,0,0])
    T = matrix(Ku,[one,MG*one,MG**2*one]).transpose()
    assert T.det()
    inverse = T.solve_right(vector(Ku,[0,1,0]))
    assert inverse[0]*I+inverse[1]*MG+inverse[2]*MG**2 == MX
    plane_D = R(plane.discriminant())
    assert plane_D.degree() == 6 and plane_D.gcd(plane_D.derivative()) == 1
    assert Ku(T.det()**2*D) == plane_D
    H = PolynomialRing(Ks,('U','V','W')); U,V,W = H.gens()
    cubic = sum(H(a)*U**k*V**i*W**(3-i-k) for i,c in enumerate(plane.list()) for k,a in enumerate(c.list()) if a)
    ja,jb = WeierstrassForm(cubic)
    ja,jb = Ks(ja),Ks(jb)
    assert 4*ja**3+27*jb**2
    j = Ks(1728*4*ja**3/(4*ja**3+27*jb**2))
    assert (j.numerator().degree(),j.denominator().degree()) == (24,21)
    return {'field':Ks,'j':j,'forward':g,'inverse':list(inverse),'plane':plane,
            'jacobian_A':ja,'jacobian_B':jb,'index_polynomial':q}


def packet(data):
    return {'forward_plane_v_coefficients_in_old_x':list(map(enc_u,data['forward'])),
            'inverse_old_x_coefficients_in_plane_v':list(map(enc_u,data['inverse'])),
            'plane_cubic_coefficients_v_then_u':[[enc_s(a) for a in c.list()] for c in data['plane'].list()],
            'jacobian_A':enc_s(data['jacobian_A']),'jacobian_B':enc_s(data['jacobian_B']),
            'j_map':enc_s(data['j']),
            'normalization_index_polynomial':[enc_s(a) for a in data['index_polynomial'].list()]}


def modular_checks(protocol,screen):
    records = []; source_paths = []
    for row in screen['records']:
        model = next(m for m in protocol['models'] if m['index'] == row['index'])
        if 'control_source' in row:
            path = ROOT/row['control_source']; source_paths.append(path)
            assert digest(path) == row['control_sha256']
            old = json.loads(path.read_text())
            results = [{'prime':old['prime'],'samples':old['samples'],
                        'j_numerator':old['j_numerator_coefficients_low_to_high'],
                        'j_denominator':old['j_denominator_coefficients_low_to_high']}]
        else: results = row['results']
        for result in results:
            field = GF(result['prime']); data = convert(field,model)
            j = data['j']; R = j.parent().ring()
            assert j == R(result['j_numerator'])/R(result['j_denominator'])
            assert all(j(a) == value for a,value in result['samples'])
            records.append({'index':row['index'],'prime':result['prime'],
                            'entire_j_map_equal':True,'sample_values_equal':len(result['samples'])})
            print('GENERIC_MODULAR_REPLAY',row['index'],result['prime'],flush=True)
    assert len(records) == 16
    return records,source_paths


def inverse(j):
    public = runpy.run_path(str(PUBLIC))
    target = EllipticCurve(QQ,list(map(QQ,public['GENERAL_WEIERSTRASS_COEFFICIENTS']))).j_invariant()
    f = j.numerator()*target.denominator()-j.denominator()*target.numerator()
    f *= lcm([a.denominator() for a in f]); f /= gcd([ZZ(a) for a in f])
    if f.leading_coefficient()<0: f = -f
    assert f.degree() == 24
    for p in prime_range(5,2000):
        reduced = f.change_ring(GF(p))
        if reduced.degree() == 24 and all(reduced(a) for a in GF(p)):
            return {'primitive_comparison':list(map(str,f.list())),'prime':int(p),
                    'comparison_mod_prime':list(map(int,reduced.list())),
                    'finite_root_count':0,'infinity_root':False}
    raise AssertionError('No exact witness found: keep inverse UNKNOWN')


def build(check=False):
    protocol = json.loads(PROTOCOL.read_text()); screen = json.loads(SCREEN.read_text())
    assert screen['protocol_sha256'] == digest(PROTOCOL)
    assert screen['excluded'] == 9 and screen['unknown'] == 1
    pending = [r['index'] for r in screen['records'] if r['status'] == 'UNKNOWN']
    assert pending == [17]
    model = next(m for m in protocol['models'] if m['index'] == 17)
    old = json.loads(OUT.read_text()) if check else None
    data = convert(QQ,model,old['QQ_conversion']['forward_plane_v_coefficients_in_old_x'] if check else None)
    print('QQ_BIRATIONAL_CONVERSION_VERIFIED',flush=True)
    qq = packet(data)
    if check: assert qq == old['QQ_conversion']
    witness = inverse(data['j'])
    checks,paths = modular_checks(protocol,screen)
    result = {'schema':'curve302.triangle-generic-conversion.v1',
              'status':'ALL_TEN_FROZEN_MW14_TRIANGLES_EXCLUDE302',
              'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),PROTOCOL,SCREEN,PUBLIC]+paths},
              'software':{'SageMath':SAGE_VERSION},'QQ_index':17,'QQ_conversion':qq,
              'QQ_inverse':witness,'independent_generic_modular_checks':checks,
              'excluded_indices':[m['index'] for m in protocol['models']],
              'boundary':'All ten pencils in the frozen289-word dictionary, not all degree-three fibrations or all302 parents. The QQ Jacobian and birational plane cubic are explicit for row17. No new302 fibre or full Weierstrass-coordinate MW basis is produced. Generic algebraic conversion verifies the earlier sixteen j-maps without their fibre-normalization algorithm.'}
    if check: assert result == old
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args = parser.parse_args()
    signal.alarm(120 if args.check else 180)
    result = build(args.check)
    if not args.check:
        assert not OUT.exists(), 'Preserve certificate'
        OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['QQ_inverse']['prime'],flush=True)
