#!/usr/bin/env sage
"""Bounded positive-witness solvers for the first labelled tangent conic.

--worker norm uses PARI relative S-unit norm equations. --worker lift uses
Hensel lifting in O_K/(p^e), followed by projective Minkowski reconstruction.
Neither a miss nor an incomplete backend is a solubility decision.
Run workers through research_runtime.supervisor; --run does this automatically.
"""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import sys
import time

from sage.all import (GF, QQ, ZZ, NumberField, PolynomialRing, RealField,
                      diagonal_matrix, matrix, pari, vector)

import run_fixed_field_tangent_conics as gate
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = gate.ROOT
WORK = ROOT/'artifacts/local/fixed-field-conic-norm-reconstruction'
SUMMARY = ROOT/'artifacts/generated-results/elliptic-curves/fixed_field_conic_solver_comparison_v1.json'
EVIDENCE = SUMMARY.with_name('fixed_field_conic_solver_comparison_evidence_v1.json.gz')
LONG_SUMMARY = ROOT/'artifacts/generated-results/elliptic-curves/fixed_field_conic_long_search_v1.json'
LONG_EVIDENCE = LONG_SUMMARY.with_name('fixed_field_conic_long_search_evidence_v1.json.gz')
LONG_WORK = {
    'norm_2700': ROOT/'artifacts/local/fixed-field-conic-norm-2700s-v1',
    'deep_2700': ROOT/'artifacts/local/fixed-field-conic-deep-2700s-v1',
    'exhaustive_p37': ROOT/'artifacts/local/fixed-field-conic-exhaustive-p37-v1',
}
PRIMES = [2, 3, 5, 7, 13, 17, 19, 31, 79, 1049, 71889448247,
          40200713707633, 491007790268548705232623905732119,
          550485889697549509378456967440666703691905279903,
          65908956571, 22407310087431607913]


def enc(x):
    return list(map(str, x.list()))


def context(pari_stack_bytes=256_000_000):
    summary = json.loads(gate.SUMMARY.read_text())
    assert gate.models.sha(gate.EVIDENCE) == summary['evidence_sha256']
    e = json.loads(gzip.decompress(gate.EVIDENCE.read_bytes()))
    # Replay the selected algebraic input directly. The historical whole-gate
    # checker pins a point-realization source since moved to a runtime snapshot;
    # do not weaken or rewrite that historical hash check for this experiment.
    assert e['models_evidence_sha256'] == gate.models.sha(gate.models.EVIDENCE)
    for p in PRIMES:
        pari.addprimes(p)
    pari.allocatemem(pari_stack_bytes, silent=True)
    K = NumberField(PolynomialRing(QQ, 't')(e['field_polynomial']), 'theta')
    dec = lambda z: K(list(map(QQ, z)))
    row = next(z for z in e['conics'] if z['mask'] == 1047173)
    C = matrix(K, 3, [dec(z) for z in row['conic_matrix']])
    models = json.loads(gzip.decompress(gate.models.EVIDENCE.read_bytes()))['covers']
    model = next(z for z in models if z['cover']['mask'] == 1047173
                 and not z['cover']['translated_by_universal_point'])
    H = [matrix(K, 4, list(map(QQ,z))) for z in model['quadric_model']['quadric_matrices']]
    pencil = dec(row['lambda'])*H[0]+H[1]
    assert pencil.rank() == 3 and pencil.det() == 0
    assert C == pencil.matrix_from_rows_and_columns(row['indices'],row['indices'])
    assert C.det() and K.defining_polynomial().change_ring(GF(e['irreducibility_prime'])).is_irreducible()
    source = json.loads(gate.models.base.SOURCE.read_text())
    anchor_B, anchor_A = map(QQ,source['anchor']['base_polynomial_ascending'][:2])
    anchor_theta = 9*K.gen()+6
    assert anchor_theta**3+anchor_A*anchor_theta+anchor_B == 0
    initial = e['norm_conic_start']
    A, B, d = dec(initial['a']), dec(initial['b']), QQ(initial['ordinate_scale'])
    assert A == -C[1,1]/2 and C[0,0] == 2 and C[0,1] == C[0,2] == 0
    assert B == (-C[2,2]/2-C[1,2]**2/(4*A))*d*d
    T0 = matrix(K, [[0, 0, A], [1, C[1, 2]*d/(2*A), 0], [0, d, 0]])
    cp = next(z for z in e['local_runs'] if z['seed'] == 4)['checkpoint']
    a, b = dec(cp['a']), dec(cp['b'])
    T = T0*matrix(K, 3, [dec(z) for z in cp['map_to_initial_norm_conic']])
    M = T.transpose()*C*T
    assert T.det() and M == M[0, 0]*diagonal_matrix(K, [1, -a, -b])
    return K, a, b, C, T


def witness(K, a, b, p, C=None, T=None):
    p = vector(K, p)
    assert any(p), 'zero vector is not a projective point'
    assert p[0]**2-a*p[1]**2-b*p[2]**2 == 0, 'nonzero exact residual'
    p /= next(z for z in p if z)
    result = {'normalized_point': [enc(z) for z in p], 'exact_residual': ['0']*3}
    if C is not None:
        q = T*p
        assert any(q) and q*C*q == 0
        q /= next(z for z in q if z)
        result['tangent_conic_point'] = [enc(z) for z in q]
    return result


def norm_worker(K, a, b, emit):
    # Variable x has higher PARI priority than y. Reuse the prepared maximal
    # order, rather than initializing from an unlabelled absolute sextic.
    nf = K.pari_nf('y')
    ay = pari(a.polynomial().change_variable_name('y')).Mod(nf.nf_get_pol())
    by = pari(b.polynomial().change_variable_name('y')).Mod(nf.nf_get_pol())
    rel = pari('x')**2-ay
    emit({'phase': 'relative_norm_initialization', 'relative_polynomial': str(rel),
          'base_maximal_basis': list(map(str, nf.nf_get_zk()))})
    data = pari.rnfisnorminit(nf, rel, 1)
    emit({'phase': 'relative_norm_solve'})
    z, remainder = pari.rnfisnorm(data, by)
    emit({'phase': 'relative_norm_returned', 'remainder': str(remainder)})
    if remainder != 1:
        return None
    coefficients = z.lift()
    decode = lambda c: K(list(map(QQ, c.lift().Vec().python()[::-1])))
    u, v = (decode(coefficients.polcoef(i)) for i in range(2))
    return vector(K, [u, v, 1])


def reduced_basis(K):
    RR = RealField(512)
    embeddings = K.embeddings(RR)
    basis = list(K.integral_basis())
    A = matrix(ZZ, [[ZZ((f(z)*2**256).round()) for f in embeddings] for z in basis])
    _, U = A.LLL(transformation=True)
    basis = [sum((U[i,j]*basis[j] for j in range(3)), K(0)) for i in range(3)]
    assert abs(U.det()) == 1
    return basis, embeddings


def local_points(K, a, b, p, seed):
    F = GF(p**3, name='z', modulus=K.defining_polynomial().change_ring(GF(p)))
    reduce = lambda z: sum((F(c)*F.gen()**i for i,c in enumerate(z.list())), F(0))
    fa, fb = reduce(a), reduce(b)
    assert fa and fb
    rng = random.Random(seed)
    for _ in range(256):
        v = sum((F(rng.randrange(p))*F.gen()**i for i in range(3)), F(0))
        u2 = fa*v*v+fb
        if u2 and u2.is_square():
            u = u2.sqrt()
            if rng.randrange(2):
                u = -u
            yield K(u.polynomial().change_ring(ZZ)), K(v.polynomial().change_ring(ZZ))


def lift_root(K, a, b, p, u, v, exponent):
    # Inert p, prime to the equation index and denominators: arithmetic in
    # (Z/p^e)[theta] is the full unramified local ring at this precision.
    m = ZZ(p)
    target = ZZ(p)**exponent
    reduce = lambda z, n: K([ZZ(c.numerator())*ZZ(c.denominator()).inverse_mod(n) % n
                            for c in z.list()])
    while m < target:
        n = min(m*m, target)
        u = reduce(u-(u*u-a*v*v-b)/(2*u), n)
        m = n
    assert all(QQ(c).valuation(p) >= exponent for c in u*u-a*v*v-b)
    return u, v


def reconstruct(K, a, b, basis, embeddings, p, exponent, u, v):
    m = ZZ(p)**exponent
    # Row lattice consists exactly of triples congruent modulo p^e to
    # an O_K multiple of (u,v,1). This reconstructs all three coordinates;
    # it never substitutes a lifted residue for an exact K-point.
    B = matrix(QQ, [z.list() for z in basis])
    def residue(z):
        c = vector(QQ, z.list())*B.inverse()
        return vector(ZZ, [ZZ(x.numerator())*ZZ(x.denominator()).inverse_mod(m) % m for x in c])
    rows = []
    for j in range(6):
        row = [ZZ(0)]*9
        row[j] = m
        rows.append(row)
    for j,z in enumerate(basis):
        rows.append(list(residue(u*z))+list(residue(v*z))+[ZZ(i==j) for i in range(3)])
    L = matrix(ZZ, rows)
    assert abs(L.det()) == m**6
    RR = embeddings[0].codomain()
    metric = matrix(RR, 9, 9)
    for coordinate, coefficient in enumerate([K(1), a, b]):
        for i,z in enumerate(basis):
            for j,f in enumerate(embeddings):
                metric[3*coordinate+i, 3*coordinate+j] = f(z)*abs(f(coefficient)).sqrt()
    approx = L.change_ring(RR)*metric
    maxentry = max(abs(z) for z in approx.list())
    scale = RR(2)**384/maxentry
    integral = matrix(ZZ, [[ZZ((x*scale).round()) for x in row] for row in approx.rows()])
    _, change = integral.LLL(transformation=True)
    assert abs(change.det()) == 1
    reduced = change*L
    points = [vector(K, [sum((row[3*j+i]*basis[i] for i in range(3)),K(0))
                         for j in range(3)]) for row in reduced.rows()]
    # Nine rows and the two signs of each pair: 81 exact candidates per cell.
    proposals = points+[points[i]+s*points[j] for i in range(9) for j in range(i) for s in (-1,1)]
    residuals = []
    for point in proposals:
        residual = point[0]**2-a*point[1]**2-b*point[2]**2
        assert all(QQ(c).valuation(p) >= exponent for c in residual)
        residuals.append(enc(residual))
        if any(point) and not residual:
            return point, {'candidate_count': len(residuals), 'residuals': residuals,
                           'reconstruction_basis': [[enc(z) for z in point] for point in points]}
    return None, {'candidate_count': len(residuals), 'residuals': residuals,
                  'reconstruction_basis': [[enc(z) for z in point] for point in points]}


def lift_worker(K, a, b, emit, wall_seconds, control=False):
    from sage.rings.number_field.number_field import NumberField_generic
    def forbidden(*args, **kwargs):
        raise RuntimeError('BNF forbidden in local reconstruction')
    NumberField_generic.pari_bnf = forbidden
    basis, embeddings = reduced_basis(K)
    emit({'phase':'basis', 'basis':[enc(z) for z in basis]})
    start = time.monotonic()
    cells = 0
    for p in [11, 23, 29, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 83, 89, 97, 101, 103, 107, 109]:
        if not K.defining_polynomial().change_ring(GF(p)).is_irreducible():
            continue
        if any(QQ(c).valuation(p)<0 for z in [a,b,*basis] for c in z.list()):
            continue
        if K.discriminant()%p == 0 or a.norm().valuation(p) or b.norm().valuation(p):
            continue
        for seed in range(1,9):
            if control:
                # A retained exact control residue tests lifting/reconstruction
                # at target scale; it is explicitly a seeded recovery control.
                u0, v0 = K(2), K(1)
            else:
                u0,v0 = next(local_points(K,a,b,p,seed))
            for exponent in [4, 8, 16, 32, 64]:
                if time.monotonic()-start > wall_seconds:
                    return None
                u,v = lift_root(K,a,b,p,u0,v0,exponent)
                point, record = reconstruct(K,a,b,basis,embeddings,p,exponent,u,v)
                cells += 1
                emit({'phase':'reconstruction_cell', 'cell':cells, 'prime':p,
                      'seed':seed, 'exponent':exponent, 'local_point':[enc(u),enc(v),enc(K(1))],
                      **record})
                if point is not None:
                    return point
            if control:
                break
    return None


def deep_lift_worker(K, a, b, emit, wall_seconds, starts):
    """A broader residue/reconstruction sweep with compact checkpoints.

    Each local residue is independent and merely proposes a finite list of
    global vectors.  Every vector is still checked in K before being omitted
    from the checkpoint.  The hash retains an integrity handle for the full
    transient list without claiming that this samples a complete local disk.
    """
    from sage.rings.number_field.number_field import NumberField_generic
    def forbidden(*args, **kwargs):
        raise RuntimeError('BNF forbidden in deep local reconstruction')
    NumberField_generic.pari_bnf = forbidden
    basis, embeddings = reduced_basis(K)
    emit({'phase':'basis', 'basis':[enc(z) for z in basis], 'profile':'deep'})
    start = time.monotonic()
    cells = 0
    transcript = sha256()
    inert_primes = []
    for p in [11, 23, 29, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 83, 89, 97, 101, 103, 107, 109]:
        if not K.defining_polynomial().change_ring(GF(p)).is_irreducible():
            continue
        if any(QQ(c).valuation(p)<0 for z in [a,b,*basis] for c in z.list()):
            continue
        if K.discriminant()%p == 0 or a.norm().valuation(p) or b.norm().valuation(p):
            continue
        inert_primes.append(p)
    emit({'phase':'eligible_primes','primes':inert_primes,'starts_per_prime':starts,
          'exponents':[16,32,64,128]})
    for p in inert_primes:
        for local_seed in range(1,starts+1):
            if time.monotonic()-start > wall_seconds:
                return None
            u0,v0 = next(local_points(K,a,b,p,local_seed))
            for exponent in [16,32,64,128]:
                if time.monotonic()-start > wall_seconds:
                    return None
                u,v = lift_root(K,a,b,p,u0,v0,exponent)
                point, raw = reconstruct(K,a,b,basis,embeddings,p,exponent,u,v)
                cells += 1
                # `raw` includes all exact residuals and the reconstruction
                # basis. Store only its SHA-256 to make the 45-minute run
                # checkpointable without exploding retained evidence size.
                packed=json.dumps({'local_point':[enc(u),enc(v),enc(K(1))],
                                   'reconstruction_basis':raw['reconstruction_basis'],
                                   'residuals':raw['residuals']},sort_keys=True,separators=(',',':'))
                cell_hash=sha256(packed.encode()).hexdigest()
                transcript.update((str(p)+'/'+str(local_seed)+'/'+str(exponent)+'/'+cell_hash+'\n').encode())
                emit({'phase':'deep_reconstruction_cell','cell':cells,'prime':p,
                      'local_seed':local_seed,'exponent':exponent,
                      'candidate_count':raw['candidate_count'],'candidate_sha256':cell_hash,
                      'transcript_sha256':transcript.hexdigest()})
                if point is not None:
                    return point
    return None


def exhaustive_affine_worker(K, a, b, emit, wall_seconds):
    """Enumerate the W=1 residue chart at the smallest eligible inert prime."""
    from sage.rings.number_field.number_field import NumberField_generic
    def forbidden(*args, **kwargs):
        raise RuntimeError('BNF forbidden in exhaustive local reconstruction')
    NumberField_generic.pari_bnf = forbidden
    basis, embeddings = reduced_basis(K)
    p = ZZ(37)
    assert K.defining_polynomial().change_ring(GF(p)).is_irreducible()
    assert K.discriminant()%p and not a.norm().valuation(p) and not b.norm().valuation(p)
    F = GF(p**3, name='z', modulus=K.defining_polynomial().change_ring(GF(p)))
    reduce = lambda z: sum((F(c)*F.gen()**i for i,c in enumerate(z.list())), F(0))
    fa,fb=reduce(a),reduce(b)
    assert fa and fb
    emit({'phase':'basis','basis':[enc(z) for z in basis],'profile':'exhaustive-affine'})
    emit({'phase':'finite_chart','prime':int(p),'chart':'W=1','field_order':int(F.order()),
          'exponent':128})
    start=time.monotonic();cells=0;finite_points=0;transcript=sha256()
    for v0 in F:
        u2=fa*v0*v0+fb
        if not u2.is_square():
            continue
        roots=[u2.sqrt()]
        if roots[0]:
            roots.append(-roots[0])
        for u0 in roots:
            finite_points += 1
            if time.monotonic()-start > wall_seconds:
                emit({'phase':'finite_chart_progress','finite_points':finite_points,'cells':cells,
                      'complete':False,'transcript_sha256':transcript.hexdigest()})
                return None
            u=K(u0.polynomial().change_ring(ZZ));v=K(v0.polynomial().change_ring(ZZ))
            u,v=lift_root(K,a,b,p,u,v,128)
            point,raw=reconstruct(K,a,b,basis,embeddings,p,128,u,v)
            cells += 1
            packed=json.dumps({'local_point':[enc(u),enc(v),enc(K(1))],
                               'reconstruction_basis':raw['reconstruction_basis'],
                               'residuals':raw['residuals']},sort_keys=True,separators=(',',':'))
            cell_hash=sha256(packed.encode()).hexdigest()
            transcript.update((str(cells)+'/'+cell_hash+'\n').encode())
            emit({'phase':'exhaustive_reconstruction_cell','cell':cells,'prime':int(p),
                  'exponent':128,'candidate_count':raw['candidate_count'],
                  'candidate_sha256':cell_hash,'transcript_sha256':transcript.hexdigest()})
            if point is not None:
                return point
    emit({'phase':'finite_chart_progress','finite_points':finite_points,'cells':cells,
          'complete':True,'transcript_sha256':transcript.hexdigest()})
    return None


def worker(args):
    K, a, b, C, T = context(args.pari_stack_bytes)
    if args.control:
        if args.worker == 'norm':
            K=NumberField(PolynomialRing(QQ,'t')([-1,-1,0,1]),'theta')
            a,b=K(2),K(7)
        else:
            b = K(4)-a
        C = T = None
    path = args.workdir/f'{args.worker}{"-control" if args.control else ""}.json'
    record = {'target_mask':None if args.control else 1047173, 'control':args.control, 'method':args.worker,
              'field_polynomial':list(map(str,K.defining_polynomial().list())),
              'a':enc(a), 'b':enc(b), 'events':[], 'point':None,
              'status':'RUNNING', 'source_evidence_sha256':gate.models.sha(gate.EVIDENCE)}
    def emit(event):
        record['events'].append(event)
        is_cell=event['phase'] in ('reconstruction_cell','deep_reconstruction_cell',
                                   'exhaustive_reconstruction_cell')
        if not is_cell or len(record['events']) % args.checkpoint_stride == 0:
            checkpoint(path,record)
        if not is_cell or event.get('cell',0) % args.checkpoint_stride == 0:
            print(json.dumps({k:v for k,v in event.items() if k not in
                             ('residuals','reconstruction_basis','local_point','basis','base_maximal_basis')}),flush=True)
    checkpoint(path,record)
    try:
        point = (norm_worker(K,a,b,emit) if args.worker=='norm' else
                 exhaustive_affine_worker(K,a,b,emit,args.seconds-5)
                 if args.profile=='exhaustive-affine' else
                 deep_lift_worker(K,a,b,emit,args.seconds-5,args.local_starts)
                 if args.profile=='deep' else
                 lift_worker(K,a,b,emit,args.seconds-5,args.control))
        if point is not None:
            record['point'] = witness(K,a,b,point,C,T)
            record['status'] = 'EXACT_POINT'
        else:
            record['status'] = 'NO_WITNESS_WITHIN_PROTOCOL'
    except Exception as exc:
        record['status'] = 'INCOMPLETE'
        record['failure'] = repr(exc)
        raise
    finally:
        checkpoint(path,record)


def replay(evidence):
    """Check retained cells and witnesses without local search, LLL, or BNF."""
    K,a,b,C,T=context()
    assert evidence['source_evidence_sha256'] == gate.models.sha(gate.EVIDENCE)
    stats={'target_conic_points':0,'target_reconstruction_cells':0,'exact_target_candidates':0,
           'seeded_reconstruction_control_points':0,'small_field_norm_control_points':0,
           'genuine_higher_covers':0,'higher_cover_searches':0,
           'point_or_sha':{str(m):'UNKNOWN' for m in gate.models.MASKS},'curve_rank_lower_bound':1}
    assert set(evidence['runs']) == {'lift','norm','lift-control','norm-control'}
    for label,run_record in evidence['runs'].items():
        r=run_record['result']; supervisor=run_record['supervisor']
        assert r['source_evidence_sha256']==evidence['source_evidence_sha256']
        assert r['method'] == label.split('-')[0]
        assert r['control'] == label.endswith('-control')
        assert r['target_mask'] == (None if r['control'] else 1047173)
        field = (NumberField(PolynomialRing(QQ,'t')([-1,-1,0,1]),'theta')
                 if label=='norm-control' else K)
        assert r['field_polynomial'] == list(map(str,field.defining_polynomial().list()))
        dec=lambda z: field(list(map(QQ,z)))
        aa,bb = dec(r['a']),dec(r['b'])
        if label=='norm-control':
            assert (aa,bb)==(field(2),field(7))
        else:
            assert aa==a and bb==(4-a if r['control'] else b)
        if r['point']:
            assert r['status']=='EXACT_POINT' and supervisor['outcome']=='completed'
            pt=vector(field,[dec(z) for z in r['point']['normalized_point']])
            assert r['point']==witness(field,aa,bb,pt,None if r['control'] else C,None if r['control'] else T)
            stats[('small_field_norm_control_points' if label=='norm-control' else
                   'seeded_reconstruction_control_points' if r['control'] else 'target_conic_points')]+=1
        else:
            assert r['status']!='EXACT_POINT'
            if supervisor['outcome']=='completed':
                assert r['status']=='NO_WITNESS_WITHIN_PROTOCOL'
            else:
                assert supervisor['outcome'] in {'strict_wall_timeout','strict_rss_limit','backend_failure'}
        if r['method']!='lift':
            continue
        basis_events=[z for z in r['events'] if z['phase']=='basis']
        assert len(basis_events)==1
        basis=[dec(z) for z in basis_events[0]['basis']]
        B=matrix(QQ,[z.list() for z in basis])
        U=B*matrix(QQ,[z.list() for z in field.integral_basis()]).inverse()
        assert all(z in ZZ for z in U.list()) and abs(U.det())==1
        seen=set()
        for cell in (z for z in r['events'] if z['phase']=='reconstruction_cell'):
            p,exponent=ZZ(cell['prime']),cell['exponent']; modulus=p**exponent
            assert p.is_prime() and exponent in [4,8,16,32,64] and cell['seed'] in range(1,9)
            assert field.defining_polynomial().change_ring(GF(p)).is_irreducible()
            assert field.discriminant()%p and not aa.norm().valuation(p) and not bb.norm().valuation(p)
            assert all(c.valuation(p)>=0 for z in [aa,bb,*basis] for c in z)
            key=(int(p),cell['seed'],exponent)
            assert key not in seen;seen.add(key)
            local=vector(field,[dec(z) for z in cell['local_point']])
            assert local[2]==1
            residual=local[0]**2-aa*local[1]**2-bb
            assert all(c.valuation(p)>=exponent for c in residual)
            rows=[vector(field,[dec(z) for z in row]) for row in cell['reconstruction_basis']]
            assert len(rows)==9 and all(len(row)==3 for row in rows)
            lattice=matrix(QQ,[sum((list(vector(QQ,z.list())*B.inverse()) for z in row),[]) for row in rows])
            assert all(z in ZZ for z in lattice.list()) and abs(lattice.det())==modulus**6
            for row in rows:
                assert all(c.valuation(p)>=exponent for j in (0,1) for c in row[j]-local[j]*row[2])
            proposals=rows+[rows[i]+s*rows[j] for i in range(9) for j in range(i) for s in (-1,1)]
            count=cell['candidate_count']
            assert 1<=count<=81
            residuals=[q[0]**2-aa*q[1]**2-bb*q[2]**2 for q in proposals[:count]]
            assert all(any(q) for q in proposals[:count])
            assert all(all(c.valuation(p)>=exponent for c in r0) for r0 in residuals)
            if count<81:
                assert not residuals[-1] and r['point']
            elif not r['point']:
                assert all(residuals)
            if not r['control']:
                stats['target_reconstruction_cells']+=1
                stats['exact_target_candidates']+=count
        if not r['control'] and supervisor['outcome']=='completed' and not r['point']:
            # These are finite recorded candidate lists, not projective-height
            # boxes or a complete enumeration of points in a local residue disk.
            expected={(p,s,e) for p in [37,41,53,73,83,89,97,101]
                      for s in range(1,9) for e in [4,8,16,32,64]}
            assert seen==expected
    return stats


def collect(workdir):
    evidence={'source_evidence_sha256':gate.models.sha(gate.EVIDENCE),'runs':{}}
    for label in ('lift','norm','lift-control','norm-control'):
        result=json.loads((workdir/f'{label}.json').read_text())
        for event in result['events']:
            event.pop('residuals',None)  # Independently recomputed from the nine rows.
        supervisor=json.loads((workdir/f'{label}.supervisor.json').read_text())
        evidence['runs'][label]={'result':result,'supervisor':supervisor,
                                 'log':(workdir/f'{label}.log').read_text()}
    arithmetic=replay(evidence)
    assert arithmetic['target_conic_points']==0, 'a positive target witness requires a new construction record'
    EVIDENCE.write_bytes(gzip.compress((json.dumps(evidence,sort_keys=True)+'\n').encode(),mtime=0))
    from sage.version import version
    summary={'schema':'elliptic-curves.fixed-field-conic-solver-comparison.v1',
             'status':'NO_TARGET_CONIC_POINT_NO_HIGHER_COVER', 'arithmetic':arithmetic,
             'checker_sha256':gate.models.sha(__file__),'evidence_sha256':gate.models.sha(EVIDENCE),
             'software':{'sage':version,'pari':str(pari.version())},
             'declared_budget':{'new_solver_seconds_ceiling':900,'per_worker_rss_bytes':2*1024**3,
                                'norm_seconds':180,'lift_seconds':480},
             'retained_runs':{label:{k:row['supervisor'][k] for k in
                                    ('outcome','wall_seconds','peak_observed_rss_bytes')}
                              for label,row in evidence['runs'].items()},
             'claim_boundary':[
                 'The target conic is known soluble; these methods returned no coordinates.',
                 'The reconstruction control is seeded with a known local point and does not calibrate target discovery probability.',
                 'The norm positive control uses a small cubic field, not the target arithmetic.',
                 'Local residues and nonzero reconstruction residuals are not global obstructions.',
                 'No higher-cover construction or cover point search was reached.']}
    checkpoint(SUMMARY,summary)
    print(json.dumps(arithmetic))


def _transcript_step(digest, *parts):
    digest.update(('/'.join(map(str,parts))+'\n').encode())


def long_replay(evidence):
    """Audit structure, dispatch bindings and compact exact-check transcripts.

    This does not rerun the 55-minute discovery calculation. Each retained
    transcript row is the SHA-256 of an exact residual list and reconstruction
    basis generated by ``reconstruct`` during the worker run; a target point
    would instead require a full explicit witness and immediate substitution.
    """
    K,a,b,C,T=context()
    assert evidence['source_evidence_sha256']==gate.models.sha(gate.EVIDENCE)
    assert set(evidence['runs'])==set(LONG_WORK)
    stats={'relative_norm_seconds':0,'deep_reconstruction_cells':0,
           'deep_exact_candidates':0,'exhaustive_w1_residues':0,
           'exhaustive_exact_candidates':0,'target_conic_points':0,
           'genuine_higher_covers':0,'higher_cover_searches':0,
           'point_or_sha':{str(m):'UNKNOWN' for m in gate.models.MASKS},
           'curve_rank_lower_bound':1}
    expected_field=list(map(str,K.defining_polynomial().list()))
    for label,row in evidence['runs'].items():
        protocol,row_result,supervisor=row['protocol'],row['result'],row['supervisor']
        assert protocol['source_evidence_sha256']==evidence['source_evidence_sha256']
        assert row_result['source_evidence_sha256']==evidence['source_evidence_sha256']
        assert row_result['target_mask']==1047173 and not row_result['control']
        assert row_result['field_polynomial']==expected_field
        assert row_result['point'] is None and row_result['status']!='EXACT_POINT'
        if label=='norm_2700':
            assert row_result['method']=='norm' and supervisor['outcome']=='strict_wall_timeout'
            assert protocol['wall_seconds']==2700 and protocol['rss_bytes']==4*1024**3
            assert supervisor['wall_seconds']>=2700 and supervisor['peak_observed_rss_bytes']<4*1024**3
            assert len(row_result['events'])==1 and row_result['events'][0]['phase']=='relative_norm_initialization'
            stats['relative_norm_seconds']=round(supervisor['wall_seconds'],6)
            continue
        assert row_result['method']=='lift' and supervisor['outcome']=='completed'
        assert protocol['rss_bytes']==4*1024**3 and supervisor['peak_observed_rss_bytes']<4*1024**3
        events=row_result['events']
        basis_events=[x for x in events if x['phase']=='basis']
        assert len(basis_events)==1
        basis=[K(list(map(QQ,z))) for z in basis_events[0]['basis']]
        U=matrix(QQ,[z.list() for z in basis])*matrix(QQ,[z.list() for z in K.integral_basis()]).inverse()
        assert all(z in ZZ for z in U.list()) and abs(U.det())==1
        if label=='deep_2700':
            assert protocol['wall_seconds']==2700
            prime_event=[x for x in events if x['phase']=='eligible_primes']
            assert len(prime_event)==1 and prime_event[0]['primes']==[37,41,53,73,83,89,97,101]
            assert prime_event[0]['starts_per_prime']==4096 and prime_event[0]['exponents']==[16,32,64,128]
            cells=[x for x in events if x['phase']=='deep_reconstruction_cell']
            assert len(cells)==8*4096*4
            digest=sha256()
            for n,cell in enumerate(cells,1):
                group=(n-1)//(4096*4); offset=(n-1)%(4096*4)
                p=[37,41,53,73,83,89,97,101][group]
                assert cell['cell']==n and cell['prime']==p
                assert cell['local_seed']==offset//4+1 and cell['exponent']==[16,32,64,128][offset%4]
                assert cell['candidate_count']==81 and len(cell['candidate_sha256'])==64
                _transcript_step(digest,p,cell['local_seed'],cell['exponent'],cell['candidate_sha256'])
                assert cell['transcript_sha256']==digest.hexdigest()
            stats['deep_reconstruction_cells']=len(cells)
            stats['deep_exact_candidates']=81*len(cells)
            continue
        assert label=='exhaustive_p37' and protocol['wall_seconds']==900
        chart=[x for x in events if x['phase']=='finite_chart']
        assert len(chart)==1 and chart[0]=={'phase':'finite_chart','prime':37,'chart':'W=1',
                                             'field_order':37**3,'exponent':128}
        cells=[x for x in events if x['phase']=='exhaustive_reconstruction_cell']
        progress=[x for x in events if x['phase']=='finite_chart_progress']
        assert len(progress)==1 and not progress[0]['complete']
        digest=sha256()
        for n,cell in enumerate(cells,1):
            assert cell['cell']==n and cell['prime']==37 and cell['exponent']==128
            assert cell['candidate_count']==81 and len(cell['candidate_sha256'])==64
            _transcript_step(digest,n,cell['candidate_sha256'])
            assert cell['transcript_sha256']==digest.hexdigest()
        assert progress[0]['cells']==len(cells) and progress[0]['finite_points']==len(cells)+1
        assert progress[0]['transcript_sha256']==digest.hexdigest()
        stats['exhaustive_w1_residues']=len(cells)
        stats['exhaustive_exact_candidates']=81*len(cells)
    return stats


def long_collect():
    evidence={'source_evidence_sha256':gate.models.sha(gate.EVIDENCE),'runs':{}}
    for label,workdir in LONG_WORK.items():
        worker='norm' if label=='norm_2700' else 'lift'
        result=json.loads((workdir/(worker+'.json')).read_text())
        supervisor=json.loads((workdir/(worker+'.supervisor.json')).read_text())
        protocol=json.loads((workdir/(worker+'.protocol.json')).read_text())
        evidence['runs'][label]={'result':result,'supervisor':supervisor,'protocol':protocol,
                                 'log':(workdir/(worker+'.log')).read_text()}
    arithmetic=long_replay(evidence)
    LONG_EVIDENCE.write_bytes(gzip.compress((json.dumps(evidence,sort_keys=True)+'\n').encode(),mtime=0))
    from sage.version import version
    summary={'schema':'elliptic-curves.fixed-field-conic-long-search.v1',
             'status':'NO_TARGET_CONIC_POINT_NO_HIGHER_COVER','arithmetic':arithmetic,
             'checker_sha256':gate.models.sha(__file__),'evidence_sha256':gate.models.sha(LONG_EVIDENCE),
             'software':{'sage':version,'pari':str(pari.version())},
             'declared_budget':{'relative_norm_seconds':2700,'deep_reconstruction_seconds':2700,
                                'exhaustive_affine_seconds':900,'per_worker_rss_bytes':4*1024**3,
                                'relative_norm_pari_stack_bytes':1_000_000_000,
                                'reconstruction_pari_stack_bytes':512_000_000},
             'claim_boundary':[
                 'The target conic is known soluble; none of these finite procedures returned coordinates.',
                 'The deep reconstruction transcripts record exact checks performed during discovery; their compact hashes are integrity records, not a new global obstruction.',
                 'The prime-37 enumeration covers only the W=1 local residue chart and stopped before completing that chart.',
                 'The relative norm initialization timeout is an algorithmic resource result, not a norm or solubility decision.',
                 'No target conic point, higher-cover map, cover search, or global obstruction results.']}
    checkpoint(LONG_SUMMARY,summary)
    print(json.dumps(arithmetic))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,default=WORK)
    parser.add_argument('--worker',choices=['norm','lift'])
    parser.add_argument('--control',action='store_true')
    parser.add_argument('--seconds',type=int,default=180)
    parser.add_argument('--rss-bytes',type=int,default=2*1024**3,
                        help='per-worker RSS limit; recorded by --run')
    parser.add_argument('--pari-stack-bytes',type=int,default=256_000_000,
                        help='PARI stack allocation for a discovery worker')
    parser.add_argument('--profile',choices=['standard','deep','exhaustive-affine'],default='standard')
    parser.add_argument('--local-starts',type=int,default=8,
                        help='deterministic local residues per eligible prime in deep mode')
    parser.add_argument('--checkpoint-stride',type=int,default=1,
                        help='checkpoint every N cell records')
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--collect',action='store_true')
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--long-collect',action='store_true')
    parser.add_argument('--long-verify',action='store_true')
    args=parser.parse_args()
    args.workdir.mkdir(parents=True,exist_ok=True)
    if args.long_collect:
        long_collect()
    elif args.long_verify:
        summary=json.loads(LONG_SUMMARY.read_text())
        assert summary['checker_sha256']==gate.models.sha(__file__)
        assert summary['evidence_sha256']==gate.models.sha(LONG_EVIDENCE)
        assert long_replay(json.loads(gzip.decompress(LONG_EVIDENCE.read_bytes())))==summary['arithmetic']
        print('PASS_LONG_SEARCH_TRANSCRIPTS; NO_TARGET_CONIC_POINT; NO_HIGHER_COVER')
    elif args.collect:
        collect(args.workdir)
    elif args.verify:
        s=json.loads(SUMMARY.read_text())
        assert s['checker_sha256']==gate.models.sha(__file__) and s['evidence_sha256']==gate.models.sha(EVIDENCE)
        assert replay(json.loads(gzip.decompress(EVIDENCE.read_bytes())))==s['arithmetic']
        print('PASS_EXACT_RECONSTRUCTION_REPLAY; NO_TARGET_CONIC_POINT; NO_HIGHER_COVER')
    elif args.run:
        assert args.worker and 5 < args.seconds <= 7200
        assert args.rss_bytes >= 256_000_000 and args.pari_stack_bytes >= 64_000_000
        assert args.local_starts > 0 and args.checkpoint_stride > 0
        label=args.worker+('-control' if args.control else '')
        command=[sys.executable,str(Path(__file__).resolve()),'--worker',args.worker,
                 '--seconds',str(args.seconds),'--workdir',str(args.workdir.resolve())]
        command += ['--pari-stack-bytes',str(args.pari_stack_bytes),'--profile',args.profile,
                    '--local-starts',str(args.local_starts),'--checkpoint-stride',str(args.checkpoint_stride)]
        if args.control:
            command.append('--control')
        checkpoint(args.workdir/f'{label}.protocol.json',
                   {'command':command,'checker_sha256':gate.models.sha(__file__),
                    'source_evidence_sha256':gate.models.sha(gate.EVIDENCE),
                    'wall_seconds':args.seconds,'rss_bytes':args.rss_bytes,
                    'pari_stack_bytes':args.pari_stack_bytes})
        result=run(command,limits=Limits(wall_seconds=args.seconds,rss_bytes=args.rss_bytes,
                                          pari_stack_bytes=args.pari_stack_bytes),
                   log_path=args.workdir/f'{label}.log',cwd=ROOT,
                   checkpoint_path=args.workdir/f'{label}.supervisor.json')
        print(json.dumps(result))
    else:
        assert args.worker
        worker(args)


if __name__=='__main__':
    main()
