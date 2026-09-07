#!/usr/bin/env python3
"""Retrospective exact visibility of the17 own sections before/after refinement."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/bounded-prime-point-portfolio-v1';ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'det1092_refined_visibility_v1.json'
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__),ROOT/'elliptic-curves/cas/search_observability.py']}
def freeze():
 p=cert.read(D/'protocol.json');paths=[D/'protocol.json'];rows=[]
 for r in p['rows']:
  if r['initial_rank']!=17:continue
  old=(ROOT/r['base_maps']).with_name('result.json');new=D/r['id']/'result.json';seed=D/r['id']/'seed.json';paths += [old,new,seed];rows.append(dict(id=r['id'],seed=str(seed.relative_to(ROOT)),before=str(old.relative_to(ROOT)),after=str(new.relative_to(ROOT))))
 assert len(rows)==2 and not (D/'visibility-protocol.json').exists();checkpoint(D/'visibility-protocol.json',dict(rows=rows,sources=sources(),inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},seconds=120,rss_bytes=2147483648,maximum_workers=1,maximum_observations=6664,boundary='Retrospective own-section visibility only. No public point or new search; no entire-coset coverage or exceptional sensitivity claim. All experimental maps and point outcomes already frozen.'))
def expected():
 p=cert.read(D/'visibility-protocol.json');assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());rows=[];total=0
 for r in p['rows']:
  seed=cert.read(ROOT/r['seed']);assert len(seed['points'])==17
  for phase in ['before','after']:
   data=cert.read(ROOT/r[phase]);assert len(data['charts'])==49;counts=Counter();minima=[None]*17;hashes=[]
   for i,c in enumerate(data['charts']):
    record={**c['search'],'completed_denominator':125000};assert record['status']=='bounded_search_complete' and record['height_bound']==125000
    for k,P in enumerate(seed['points']):
     x,y=map(cert.F,P)
     for sign in (-1,1):
      v=point_visibility(record,(x,sign*y));counts[v['status']]+=1;entry=dict(chart=i,section=k,sign=sign,visibility=v);hashes.append(digest(entry));total+=1
      H=v.get('minimum_affine_height')
      if H is not None and (minima[k] is None or H<minima[k]['height']):minima[k]=dict(height=H,chart=i,sign=sign,coordinate=v['coordinate'])
   assert counts['VISIBLE_NOT_RECORDED']==0
   rows.append(dict(id=r['id'],phase=phase,counts=dict(counts),observations=len(hashes),observation_digest=digest(hashes),section_minima=minima,smallest_known_section_coordinate_height=min(x['height'] for x in minima if x)))
 assert total==6664
 return dict(schema='elliptic-curves.det1092-refined-visibility.v1',status='PASS',rows=rows,observations=total,protocol_sha256=cert.hashed(D/'visibility-protocol.json'),sources=sources(),inputs=p['inputs'],boundary=p['boundary'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','run','check']);a=p.parse_args()
 if a.mode=='freeze':freeze()
 else:
  r=expected()
  if a.mode=='check':assert cert.read(OUT)==r
  else:assert not OUT.exists();checkpoint(OUT,r)
  print('PASS6664 exact own-section visibility checks',[(a['id'],a['phase'],a['counts'],len(str(a['smallest_known_section_coordinate_height']))) for a in r['rows']],flush=True)
