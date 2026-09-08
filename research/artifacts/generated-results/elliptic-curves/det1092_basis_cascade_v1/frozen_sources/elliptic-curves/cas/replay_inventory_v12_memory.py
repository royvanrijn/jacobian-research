#!/usr/bin/env python3
"""Exact V12 inventory/CSV replay with an explicit ephemeral finite cache."""
import argparse,csv,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import export_new_high_rank_curve_index_v12 as original
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'

def replay(path,output):
    if output.exists():raise FileExistsError('preserve memory inventory replay')
    data=cert.read(path)
    if data['sources']!=original.sources() or data['previous_index_sha256']!=cert.hashed(original.OLD):raise ArithmeticError('original inventory sources changed')
    expected,previous=original.expected_rows()
    if data['curves']!=expected or data['catalogue']!=cert.read(original.ADDITIONAL)['catalogue']:raise ArithmeticError('inventory extraction or stable IDs differ')
    inputs={}
    for name,h in data['source_certificate_hashes'].items():
        if cert.hashed(ART/name)!=h:raise ArithmeticError('source certificate changed')
        inputs[name]=cert.read(ART/name)['curves']
    seen=set();ranks=[]
    for r in data['curves']:
        source=inputs[r['source_certificate']][r['source_curve_index']]
        if any(r[k]!=source[k] for k in ('parameter','curve','points','rank_certificate')):raise ArithmeticError('source point proof differs')
        model=tuple(map(cert.F,r['curve']));points=[tuple(map(cert.F,P)) for P in r['points']];old=r['rank_certificate'];actual=checked_rank(model,points,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.dumps(actual,sort_keys=True)!=json.dumps(old,sort_keys=True) or r['rank_lower_bound']!=len(points):raise ArithmeticError('exact finite certificate differs')
        inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
        if str(j)!=r['j_invariant'] or j in seen:raise ArithmeticError('j identity or inventory distinctness differs')
        if any(cert.isomorphic(model,q['ainvs']) for q in data['catalogue']['equations']):raise ArithmeticError('catalogue match in new inventory')
        seen.add(j);ranks.append({'id':r['id'],'rank_lower_bound':len(points)});print('MEMORY REPLAY V12',r['id'],'rank >=',len(points),flush=True)
    csv_path=path.with_suffix('.csv');rows=list(csv.reader(csv_path.open(newline='')));header=['id','rank_lower_bound','family','parameter','a1','a2','a3','a4','a6','source_certificate','source_curve_index'];expected_csv=[header]+[[str(v) for v in [r['id'],r['rank_lower_bound'],r['family'],r['parameter'],*r['curve'],r['source_certificate'],r['source_curve_index']]] for r in data['curves']]
    if rows!=expected_csv:raise ArithmeticError('equation CSV differs from proved inventory')
    sources={str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(original.__file__).resolve(),ROOT/'elliptic-curves/cas/memory_rank_certificate.py',path,csv_path)}
    checkpoint(output,{'schema':'elliptic-curves.inventory-v12-memory-replay.v1','status':'PASS','sources':sources,'curves_checked':len(seen),'rank_lower_bounds':ranks,'csv_checked':True,'claim_boundary':'Same original V12 extraction, source bindings, independent finite point certificates, catalogue exclusion and distinctness checks, with explicit MemoryFactStore for ephemeral exact finite computations. Earlier frozen exporters and their outcomes remain unchanged.'});print('EXACT MEMORY INVENTORY PASS',len(seen),'AND CSV',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=ART/'new_high_rank_curve_index_v12.json');p.add_argument('--output',type=Path,required=True);a=p.parse_args();replay(a.input.resolve(),a.output.resolve())
