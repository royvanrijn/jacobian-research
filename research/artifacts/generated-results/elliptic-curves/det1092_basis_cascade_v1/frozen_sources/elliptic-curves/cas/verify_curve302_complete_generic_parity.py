#!/usr/bin/env python3
"""Exact complete parity/norm and quartic identity audit; no CVP optimality claim."""
import gzip,heapq,math
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
from half_lattice_pointed_sieve import linear_combination
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'artifacts/local/elliptic-curves/curve302-complete-generic-parity-v1'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_complete_generic_parity_geometry_v1.json'

def main():
    assert not OUT.exists()
    p=cert.read(D/'protocol.json');data=cert.read(D/'maps.json');ledger=cert.read(D/'ledger.json')
    assert ledger['status']=='PASS' and data['status']=='COMPLETE_DECLARED_MAPS'
    for name,h in p['sources'].items():assert cert.hashed(ROOT/name)==h
    seedpath=ROOT/'artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2/curve302-generic17/seed.json'
    assert cert.hashed(seedpath)==p['seed_sha256']
    seed=cert.read(seedpath);g=data['rounded_gram'];u=data['change_of_basis'];heap=[];mask=0
    oldpath=seedpath.parent/'maps.json';old=cert.read(oldpath)
    assert old['rounded_gram']==g
    oldsample={c['parity']:c for c in old['sample']};matched=0
    for b in data['blocks']:
        path=ROOT/b['path'];assert cert.hashed(path)==b['sha256']
        entries=__import__('json').loads(gzip.decompress(path.read_bytes()));assert len(entries)==b['rows']
        for c in entries:
            mask+=1;assert c['parity']==mask
            w=c['representative'];v=c['reduced_representative'];assert len(w)==len(v)==17
            assert all(w[j]%2==((mask>>j)&1) for j in range(17))
            assert w==[sum(v[i]*u[i][j] for i in range(17)) for j in range(17)]
            nonzero=[(i,x) for i,x in enumerate(w) if x]
            assert sum(x*g[i][j]*y for i,x in nonzero for j,y in nonzero)==c['metric_norm']
            if mask in oldsample:
                assert c==oldsample[mask];matched+=1
            key=(c['metric_norm'],-mask)
            if len(heap)<49:heapq.heappush(heap,(key,c))
            elif key>heap[0][0]:heapq.heapreplace(heap,(key,c))
    assert mask==131071 and matched==2048
    centres=[c for key,c in sorted(heap,reverse=True)];assert centres==data['centres']
    model=tuple(map(F,seed['curve']));points=[tuple(map(F,P)) for P in seed['points']]
    assert len(data['rows'])==49
    for c,row in zip(centres,data['rows']):
        assert c==row['centre'];x,y=linear_combination(model,points,c['representative'])
        raw=[-3*x*x-4*model[3],-8*y,-6*x,F(0),F(1)]
        assert raw==list(map(F,row['raw_coefficients']))
        a,b,c1,d=map(F,row['matrix']);assert a*d-b*c1
        result=[F(0)]*5
        for i in range(5):
            for j in range(i+1):
                for k in range(5-i):
                    result[j+k]+=raw[i]*math.comb(i,j)*a**j*b**(i-j)*math.comb(4-i,k)*c1**k*d**(4-i-k)
        ratio=F(row['square_ratio']);assert ratio>0
        assert math.isqrt(ratio.numerator)**2==ratio.numerator and math.isqrt(ratio.denominator)**2==ratio.denominator
        disc=list(map(F,row['discriminant_quartic']))
        P=list(map(F,row['reduced_P']));Q=list(map(F,row['reduced_Q']))
        assert disc==[4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0<=i-j<3) for i in range(5)]
        assert result==[ratio*v for v in disc]
    reads=cert.read(D/'data-access.json')
    assert all((ROOT/n).is_relative_to(D) or ROOT/n==seedpath for n in reads)
    overlap=len({c['parity'] for c in centres}&{c['parity'] for c in old['centres']})
    files=[D/n for n in ['protocol.json','maps.json','metric.json','ledger.json','data-access.json']]+[oldpath]
    checkpoint(OUT,dict(schema='elliptic-curves.curve302-complete-generic-parity-geometry.v1',status='PASS',
        checked_classes=mask,checked_maps=49,old_sample_reproduced=matched,old_selected_overlap=overlap,
        old_norm_range=[min(c['metric_norm'] for c in old['centres']),max(c['metric_norm'] for c in old['centres'])],
        new_norm_range=[min(c['metric_norm'] for c in centres),max(c['metric_norm'] for c in centres)],
        geometry_seconds=ledger['supervision']['wall_seconds'],
        inputs={str(q.relative_to(ROOT)):cert.hashed(q) for q in files},checker_sha256=cert.hashed(Path(__file__)),
        scope='All nonzero parities of the original17-section subgroup, exact rounded norms, transport and all49 factor-free quartic identities verified. The original2048 candidates reproduce exactly;48 of the top49 complete-enumeration choices were absent from that sample. Complete subgroup parity coverage does not cover every unknown ambient class, prove CVP optimality or produce a new rational point. No point worker or active48-fibre rule changed.'))
    print('PASS131071 parity/norm checks,49 exact maps,48 new choices; no point search',flush=True)

if __name__=='__main__':main()
