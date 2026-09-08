#!/usr/bin/env sage-python
"""Refine every fixed map and witness exposure beyond declared same-centre history."""
import argparse,sys
from pathlib import Path
from collections import Counter
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import bounded_prime_point_portfolio as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.projective_box_change import classify,transition,apply,primitive
from pointed_quartic_search import PointedQuarticSearch,point_record
from half_lattice_pointed_sieve import linear_combination
import pari_pointed_backend as backend
mapper=SourceFileLoader('bounded_prime_mapper',str(CAS/'bounded_prime_pari_mapping_v2.sage')).load_module()
def main(index):
 batch.configure(index);p=batch.protocol();r=batch.ROW;old=cert.read(batch.ROOT/r['base_maps']);seed=cert.read(batch.SEED);out=batch.D/'maps.json';assert not out.exists();assert len(old['rows'])==49 and [c['parity'] for c in old['sample']]==batch.masks(p)
 prior=[(cert.read(batch.ROOT/h['maps_path']),h['height']) for h in r['prior_same_centre_boxes']]
 assert all(m['centres']==old['centres'] for m,h in prior)
 model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points']);data={k:v for k,v in old.items() if k!='rows'};data.update(status='RUNNING_MAPS',protocol_hash=digest(p),rows=[]);checkpoint(out,data);mapper.pari.allocatemem(256000000,silent=True)
 for i,base in enumerate(old['rows']):
  m=mapper.refine(base,p['minimization_primes']);history=[(base,125000)]+[(v['rows'][i],h) for v,h in prior];comparisons=[dict(height=h,transition=list(transition(oldmap['matrix'],m['matrix']))) for oldmap,h in history]
  candidates=[(1,0),(0,1),(1,1),(-1,1)]
  for oldmap,h in history:
   c=classify(oldmap['matrix'],m['matrix'],p['height'])
   if c['new_coordinate_witness']:candidates.append(tuple(c['new_coordinate_witness']))
  witness=next((primitive(q) for q in candidates if max(map(abs,primitive(q)))<=p['height'] and all(max(map(abs,apply(c['transition'],primitive(q))))>c['height'] for c in comparisons)),None)
  m['history_comparison']=dict(status='PROVED_NEW_VS_DECLARED_HISTORY' if witness else 'NO_NEW_HISTORY_WITNESS',comparisons=comparisons,witness=list(witness) if witness else None)
  C=linear_combination(model,points,m['centre']['representative']);search=PointedQuarticSearch(curve=seed['curve'],subgroup=[],centre={'point':point_record(C)},coordinate_policy=m['coordinate_policy']);backend.validate_map(search,m)
  data['rows'].append(m);checkpoint(out,data)
 data['history_counts']=dict(Counter(m['history_comparison']['status'] for m in data['rows']));data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('FROZEN49',data['history_counts'],flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,required=True);a=parser.parse_args();main(a.index)
