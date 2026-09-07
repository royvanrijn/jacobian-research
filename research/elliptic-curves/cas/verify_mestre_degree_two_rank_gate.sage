#!/usr/bin/env sage-python
"""Independent witness replay: no producer import or short-vector search.

Check every parity class, three integral norm-minus-two roots orthogonal
to each admitted fibre, and their independence modulo the fibre.
"""
import hashlib
import json
from pathlib import Path
import signal
import numpy as np
from sage.all import ZZ, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def verify():
    p = ART/'mestre_degree_two_rank_gate_v1.json'
    d = json.loads(p.read_text())
    source_path = ART/'mestre_rational_ns_gram_v2.json'
    source = json.loads(source_path.read_text())
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == d['input_sha256'][str(source_path.relative_to(ROOT))]
    first = source['rows'][0]
    assert [r['outer_u'] for r in source['rows']] == ['11','13','17','19','23','29']
    assert all(r['rational_NS_Gram'] == first['rational_NS_Gram'] for r in source['rows'])
    G = matrix(ZZ, first['rational_NS_Gram'])
    e = vector(ZZ,d['old_fibre']); h = vector(ZZ,d['hyperbolic_partner'])
    V = matrix(ZZ,d['frame_basis']); L = matrix(ZZ,d['frame_gram'])
    change = matrix(ZZ,[e,h]+V.rows())
    assert abs(change.det()) == 1 and e*G*e == h*G*h == 0 and e*G*h == 1
    assert all(c == 0 for c in V*G*e) and all(c == 0 for c in V*G*h)
    assert -V*G*V.transpose() == L and L.is_positive_definite() and L.det() == 468
    assert e == vector(ZZ,[1]+[0]*17)
    assert h == vector(ZZ,[1,1]+[0]*16)
    wp = ROOT/d['witness_file']
    assert hashlib.sha256(wp.read_bytes()).hexdigest() == d['witness_sha256']
    with np.load(wp,allow_pickle=False) as data:
        masks = data['masks'].astype(np.int64)
        numerator = data['root_numerators'].astype(np.int64)
        columns = data['nonzero_minor_columns'].astype(np.int64)
    assert numerator.shape == (32509,3,16) and columns.shape == (32509,3)
    assert np.all((columns>=0)&(columns<16))
    gram = np.array(L.rows(),dtype=np.int64)
    oldgram = np.array(G.rows(),dtype=np.int64)
    basis = np.array(V.rows(),dtype=np.int64)
    ee = np.array(list(e),dtype=np.int64); hh = np.array(list(h),dtype=np.int64)
    # Explicit bounds keep all following integer matrix arithmetic below 2^63.
    assert max(np.max(np.abs(a)) for a in [gram,oldgram,basis,numerator]) < 1000
    words = ((np.arange(65536,dtype=np.int64)[:,None]>>np.arange(16))&1)
    pair = words@gram; norm = np.sum(pair*words,axis=1)
    divisibility = np.gcd.reduce(np.column_stack((np.full(65536,2),norm//4,pair)),axis=1)
    expected = np.flatnonzero((norm%4 == 0)&(divisibility == 1))
    assert np.array_equal(masks,expected) and len(masks) == 32509
    w = words[masks]; k = norm[masks]//4
    fibres = k[:,None]*ee+2*hh+w@basis
    assert np.max(np.abs(fibres)) < 1000000
    assert np.all(np.sum((fibres@oldgram)*fibres,axis=1) == 0)
    assert np.all(np.gcd.reduce(fibres@oldgram,axis=1) == 1)
    assert np.all(fibres@oldgram@ee == 2)
    yy = numerator.reshape((-1,16))
    assert np.all(np.sum((yy@gram)*yy,axis=1) == 8)
    parity = numerator%2
    old = np.all(parity == 0,axis=2)
    new = np.all(parity == w[:,None,:],axis=2)
    assert np.all(old|new) and not np.any(old&new)
    b = new.astype(np.int64)
    wx = np.sum(numerator*pair[masks,None,:],axis=2)
    assert np.all((wx+2*b*k[:,None])%4 == 0)
    a = (wx+2*b*k[:,None])//4
    x = (numerator+b[:,:,None]*w[:,None,:])//2
    roots = a[:,:,None]*ee+b[:,:,None]*hh+x@basis
    assert max(np.max(np.abs(roots)),np.max(np.abs(fibres))) < 1000000
    flat = roots.reshape((-1,18))
    assert np.all(np.sum((flat@oldgram)*flat,axis=1) == -2)
    assert np.all(np.sum((roots@oldgram)*fibres[:,None,:],axis=2) == 0)
    # The quotient coordinate x-(b/2)w is numerator/2. A nonzero
    # 3x3 minor proves the three roots independent modulo the fibre.
    small = np.take_along_axis(numerator,columns[:,None,:],axis=2)
    v = small
    det = (v[:,0,0]*(v[:,1,1]*v[:,2,2]-v[:,1,2]*v[:,2,1])
           -v[:,0,1]*(v[:,1,0]*v[:,2,2]-v[:,1,2]*v[:,2,0])
           +v[:,0,2]*(v[:,1,0]*v[:,2,1]-v[:,1,1]*v[:,2,0]))
    assert np.all(det%101 != 0)
    assert d['root_rank_lower_bound'] == 3 and d['arithmetic_MW_rank_upper_bound'] == 18-2-3
    print('PASS_ALL_PARITY_CLASSES',len(words),flush=True)
    print('PASS_DIVISIBILITY_ONE_CLASSES',len(masks),flush=True)
    print('PASS_INDEPENDENT_VERTICAL_ROOTS',len(masks)*3,flush=True)
    print('PASS_MW_RANK_AT_MOST13_FOR_EVERY_OLD_DEGREE_TWO_Q_JACOBIAN_FIBRATION',flush=True)


if __name__ == '__main__':
    signal.alarm(120)
    verify()
