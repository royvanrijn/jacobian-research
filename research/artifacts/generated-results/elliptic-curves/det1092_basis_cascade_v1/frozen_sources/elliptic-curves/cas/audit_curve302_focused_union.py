#!/usr/bin/env python3
"""Post-search full31 union certificate and retrospective public-point visibility.

Neither the union nor the oracle coordinates feed the frozen point campaign.
Finite rank31 of this union is not an exact-rank or rational-span upper bound.
"""
import argparse
from dataclasses import asdict
from pathlib import Path
import curve302_focused_point_exposure_v2 as batch
import audit_recorded_point_mod2_rank_v3 as finite
from memory_rank_certificate import checked_rank
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint
from mod2_reduction_independence import _primes_up_to
from search_observability import point_visibility

ROOT,CAS,ART,D=batch.ROOT,batch.CAS,batch.ART,batch.BATCH
OUT=ART/'curve302_focused_union_mod2_v1.json'
VIS=ART/'curve302_focused_visibility_v1.json'

def build():
    assert not OUT.exists() and not VIS.exists()
    assert batch.cert.read(D/'run-ledger.json')['status']=='PASS'
    seedpath=D/'curve302-full31/seed.json';seed=batch.cert.read(seedpath)
    model=tuple(map(batch.cert.F,seed['curve']));pts=[tuple(map(batch.cert.F,P)) for P in seed['points']]
    assert len(pts)==31 and not any(model[:3])
    seen={(x,abs(y)) for x,y in pts};paths=[seedpath,D/'protocol.json'];runs=[]
    for ident in ['curve302-generic17','curve302-full31']:
        src=D/ident/'result.json';r=batch.cert.read(src);paths.append(src)
        assert len(r['charts'])==49 and all(c['search']['status']=='bounded_search_complete' for c in r['charts'])
        assert r['curve']==seed['curve'];runs.append((ident,r))
        # Include original generic section values as well as every raw output.
        candidates=[tuple(map(batch.cert.F,P)) for P in r['initial_state']['state']['reductions']['points']]
        candidates += [(batch.cert.F(P['x']),batch.cert.F(P['y'])) for c in r['charts'] for P in c['search']['finite_curve_points']]
        for P in candidates:
            key=P[0],abs(P[1])
            if key not in seen:seen.add(key);pts.append(P)
    assert all(batch.cert.is_on_weierstrass_curve(model,P) for P in pts)
    cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];processed=[]
    disc=-16*(4*model[3]**3+27*model[4]**2)
    for p in _primes_up_to(997):
        if p==2 or disc.numerator%p==0:continue
        sig=finite.signature(cache,model,pts,p);before=len(pivots)
        for row in sig.rows:finite.insert(pivots,row)
        processed.append(p)
        if len(pivots)>before:signatures.append(asdict(sig))
    assert len(processed)>100
    indices=sorted(pivots);points=[pts[i] for i in indices]
    proof=checked_rank(model,points,[s['prime'] for s in signatures],seed['rank_certificate']['no_rational_2_torsion_prime'])
    inputs={str(p.relative_to(ROOT)):batch.cert.hashed(p) for p in paths}
    checkpoint(OUT,dict(schema='elliptic-curves.focused302-recorded-union.v1',status='COMPLETE_DECLARED_FINITE_AUDIT',
        sources=finite.sources(),producer_sha256=batch.cert.hashed(Path(__file__)),inputs=inputs,
        curve=seed['curve'],family='det1092',parameter='0',points=[list(map(str,P)) for P in pts],
        signatures=signatures,rank_certificate=proof,rank_lower_bound=len(indices),
        independent_column_indices=indices,independent_points=[list(map(str,P)) for P in points],
        processed_primes=processed,original_rank_lower_bound=31,
        claim_boundary='All supplied31 public points,17 generic images and every point returned by the two completed302 searches, with exact sign deduplication. Finite-certified lower bound only; rank31 does not prove every point in the rational span or an upper bound. No new search.'))
    visibility=[]
    for ident,r in runs:
        rows=[];counts={};minimum=[None]*31;inside=set()
        for ci,c in enumerate(r['charts']):
            record=dict(c['search'],completed_denominator=125000)
            for j,P in enumerate(seed['points']):
                x,y=map(batch.cert.F,P)
                for sign in [1,-1]:
                    v=point_visibility(record,(x,sign*y));counts[v['status']]=counts.get(v['status'],0)+1
                    h=v['minimum_affine_height']
                    if h is not None:minimum[j]=h if minimum[j] is None else min(h,minimum[j])
                    if v.get('in_requested_box'):inside.add(j)
                    rows.append(dict(chart=ci,public_index=j,sign=sign,**v))
        checkpoint(D/ident/'public-visibility-after-search.json',dict(rows=rows))
        visibility.append(dict(id=ident,observations=len(rows),status_counts=counts,
            minimum_affine_heights=minimum,public_indices_inside_any_box=sorted(inside),
            details_sha256=batch.cert.hashed(D/ident/'public-visibility-after-search.json')))
    checkpoint(VIS,dict(status='PASS',inputs=inputs,producer_sha256=batch.cert.hashed(Path(__file__)),rows=visibility,
        scope='Retrospective only after all four runs and proofs. Exact original31 signed representatives in frozen charts; these are not all translates or all representatives of their quotient directions. Completed GP boxes retain their original trust boundary. No success threshold or subsequent coordinates chosen from this oracle.'))
    print('UNION',len(pts),'points; rank >=',len(indices),'visibility',[(r['id'],r['status_counts']) for r in visibility],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.parse_args();build()
