#!/usr/bin/env sage-python
"""Exact generic model and section transport between08234 and published R17."""
import argparse
from pathlib import Path
import sys
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, pari, RealField
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ATLAS = ROOT/'artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json'
MODEL = ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
SECTIONS = ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'

def setup():
    R = PolynomialRing(QQ, 't'); t = R.gen(); K = R.fraction_field()
    a = next(r for r in cert.read(ATLAS)['families'] if r['family'] == '08234')
    b = cert.read(MODEL)
    A = R(list(map(QQ,b['A_coefficients_low_to_high'])))
    B = R(list(map(QQ,b['B_coefficients_low_to_high'])))
    s = -26*t-50
    aa = R(list(map(QQ,a['A_coefficients_low_to_high'])))(s)
    bb = R(list(map(QQ,a['B_coefficients_low_to_high'])))(s)
    ra, rb = K(aa)/A, K(bb)/B
    if ra not in QQ or rb not in QQ:
        raise ArithmeticError('model scale is not constant')
    u2 = QQ(rb/ra)
    if not u2.is_square():
        raise ArithmeticError('model twist is not a rational isomorphism')
    u = u2.sqrt()
    if aa != u**4*A or bb != u**6*B:
        raise ArithmeticError('generic model identity failed')
    E = EllipticCurve(K, [0,0,0,A,B])
    old = []
    for row in cert.read(SECTIONS)['sections']:
        x = R(list(map(QQ,row['x_coefficients_low_to_high'])))
        if 'y_coefficients_low_to_high' in row:
            y = R(list(map(QQ,row['y_coefficients_low_to_high'])))
        else:
            chord = row['chord']; q = old[chord['reference_basis_index']]
            y = q[1]+R(list(map(QQ,chord['slope_coefficients_low_to_high'])))*(x-q[0])
        old.append(E(x,y))
    def rational(row):
        return K(R(list(map(QQ,row['numerator_coefficients_low_to_high'])))(s))/R(list(map(QQ,row['denominator_coefficients_low_to_high'])))(s)
    new = [E(rational(row['X'])/u**2,rational(row['Y'])/u**3) for row in a['sections']]
    if len(old) != 17 or len(new) != 17:
        raise ArithmeticError('section count differs')
    return E, old, new, u, A, B

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ATLAS,MODEL,SECTIONS)}

def verify(E, old, new, words):
    if len(words) != 17 or any(len(row) != 17 or any(type(x) is not int for x in row) for row in words):
        raise ArithmeticError('invalid integral section matrix')
    for i, word in enumerate(words):
        point = sum((int(k)*p for k,p in zip(word,old) if k), E(0))
        if point != new[i]:
            raise ArithmeticError('generic section word failed')
    det = matrix(ZZ, words).det()
    if abs(det) != 1:
        raise ArithmeticError('same integral generic subgroup not established')
    return int(det)

def build(output):
    if output.exists():
        raise FileExistsError('preserve transport certificate')
    E, old, new, u, A, B = setup()
    words = []
    for p in new:
        found = None
        for j,q in enumerate(old):
            for sign in (1,-1):
                if p == sign*q:
                    found = [sign if k == j else 0 for k in range(17)]
        words.append(found)
    unmatched = [i for i,word in enumerate(words) if word is None]
    if unmatched:
        E1 = EllipticCurve(QQ,[0,0,0,A(1),B(1)])
        points = [E1(p[0](1),p[1](1)) for p in old+new]
        pari.default('realprecision',90)
        raw = pari(E1).ellheightmatrix([list(p.xy()) for p in points])
        RF = RealField(256)
        H = matrix(RF,17,17,lambda i,j:RF(str(raw[i,j])))
        for i in unmatched:
            v = H.inverse()*matrix(RF,17,1,lambda j,k:RF(str(raw[j,17+i])))
            word = [int(v[j,0].round()) for j in range(17)]
            if max(abs(v[j,0]-word[j]) for j in range(17)) > RF('1e-30'):
                raise ArithmeticError('integral word proposal unresolved')
            words[i] = word
    det = verify(E,old,new,words)
    checkpoint(output, {'schema':'elliptic-curves.compact-published-r17-generic-transport.v1',
        'sources':sources(),'status':'PASS_EXACT_GENERIC_TRANSPORT',
        'published_to_compact_parameter':{'constant':-50,'linear':-26},
        'compact_to_published_x_scale':str(1/u**2),'compact_to_published_y_scale':str(1/u**3),
        'compact_model_scale_u':str(u),'compact_sections_in_published_basis':words,
        'matrix_determinant':det,'numeric_word_proposals':len(unmatched),
        'claim_boundary':'Exact Q(t) coefficient identities and all17 group words with unimodular matrix prove the same integral generic subgroup under this base/model isomorphism. No extra rank direction from these duplicated presentations.'})
    print('EXACT GENERIC R17 TRANSPORT', 'u',u,'det',det,'nontrivial words',len(unmatched),flush=True)

def check(path):
    data = cert.read(path)
    if data['sources'] != sources():
        raise ArithmeticError('transport sources changed')
    E,old,new,u,_,_ = setup()
    if str(u) != data['compact_model_scale_u'] or verify(E,old,new,data['compact_sections_in_published_basis']) != data['matrix_determinant']:
        raise ArithmeticError('generic transport replay failed')
    print('REPLAYED GENERIC R17 TRANSPORT WITHOUT NUMERICAL HEIGHTS',flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path);p.add_argument('--check',type=Path)
    a = p.parse_args(); check(a.check.resolve()) if a.check else build(a.output.resolve())
