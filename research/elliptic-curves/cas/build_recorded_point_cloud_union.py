#!/usr/bin/env python3
"""Exact finite independence on a hash-pinned union of existing point artifacts."""
import argparse
from pathlib import Path
from dataclasses import asdict
import audit_recorded_point_mod2_rank_v3 as finite
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint
from mod2_reduction_independence import _primes_up_to
ROOT=Path(__file__).resolve().parents[2]
def main(inp,out):
    assert not out.exists();p=cert.read(inp);points=[];seen=set();model=None
    for source in p['inputs']:
        path=ROOT/source['path'];assert cert.hashed(path)==source['sha256'];d=cert.read(path)
        m=tuple(map(cert.F,d['curve']))
        if model is None:model=m;torsion=d['rank_certificate']['no_rational_2_torsion_prime']
        assert m==model and not any(model[:3])
        for raw in d['points']:
            P=tuple(map(cert.F,raw));assert cert.is_on_weierstrass_curve(model,P)
            key=P[0],abs(P[1])
            if key not in seen:seen.add(key);points.append(P)
    cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];primes=[]
    disc=-16*(4*model[3]**3+27*model[4]**2);assert disc.denominator==1
    for q in _primes_up_to(p['prime_bound']):
        if q==2 or disc.numerator%q==0:continue
        sig=finite.signature(cache,model,points,q);before=len(pivots)
        for row in sig.rows:finite.insert(pivots,row)
        primes.append(q)
        if len(pivots)>before:signatures.append(asdict(sig))
    assert primes and signatures
    indices=sorted(pivots);basis=[points[i] for i in indices]
    proof=checked_rank(model,basis,[s['prime'] for s in signatures],torsion)
    checkpoint(out,dict(schema='elliptic-curves.recorded-cloud-union.v1',status='COMPLETE_DECLARED_FINITE_AUDIT',
        sources=finite.sources(),producer_sha256=cert.hashed(Path(__file__)),input_path=str(inp.relative_to(ROOT)),
        input_sha256=cert.hashed(inp),inputs=p['inputs'],curve=list(map(str,model)),family=p['family'],parameter=p['parameter'],
        points=[list(map(str,P)) for P in points],signatures=signatures,processed_primes=primes,
        independent_column_indices=indices,independent_points=[list(map(str,P)) for P in basis],
        rank_certificate=proof,rank_lower_bound=len(basis),
        claim_boundary='Exact union of the specified point lists with sign deduplication and certified lower bound only. No new search and no rank or rational-span upper bound.'))
    print('UNION',len(points),'points; rank >=',len(basis),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.input.resolve(),a.output.resolve())
