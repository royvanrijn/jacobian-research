#!/usr/bin/env python3
"""Frozen training-prime comparison on seven already retained rank26 curves."""
import json,argparse
from pathlib import Path
from fractions import Fraction as Q
import score_retained_native19 as scalar
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=scalar.ROOT;ART=scalar.ART;D=ROOT/'artifacts/local/elliptic-curves/retained26-gap-scores-v1';OUT=ART/'retained26_gap_scores_v1.json';GAP=ART/'retained26_source_gap_v1.json';INDEX=ART/'new_high_rank_curve_index_v22.json'
def freeze():
 assert not (D/'protocol.json').exists();gap=cert.read(GAP);inventory=cert.read(INDEX);rows=[];paths=[Path(__file__).resolve(),Path(scalar.__file__),GAP,INDEX,ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']
 for ident in gap['eligible_ids']:
  r=next(r for r in inventory['curves'] if r['id']==ident);assert r['rank_lower_bound']==26 and not r['current_catalogue_matches'];model=[Q(c) for c in r['curve']];assert not any(model[:3])
  model=[Q(0),Q(0),Q(0),12**4*model[3],12**6*model[4]];assert all(c.denominator==1 for c in model)
  rows.append({'id':ident,'family':r['family'],'parameter':r['parameter'],'model':list(map(str,model)),'model_coefficient_bits':max(abs(c.numerator).bit_length() for c in model)})
 checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained26-gap-scores.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'primes':scalar.PRIMES,'direct_check_primes':scalar.CHECKS,'gp_sha256':cert.hashed(scalar.GP),'seconds_per_curve':30,'rss_bytes':536870912,'workers':1,
  'score':'Reuse exact good-prime local scaling and rounded sum(2-a_p)log(p)/(p+1-a_p) from5 through32749. Fixed scale12 makes short coefficients integral and is invertible at every scoring prime. No validation primes are computed.',
  'selection':'Sort by decreasing score, then coefficient bits and id. Take the first two rows from different fibration labels. All seven already have certified26 and the same source-centre rank17; no new parameter, public target or outcome-based refill.',
  'following_exposure':'At most49 factor-free own26 charts per selected curve, using the calibrated2048-mask policy and unchanged125000/10-second boxes. Freeze both map sets before any points. No automatic following wave.'})
 print('FROZEN7 retained26 equations and fixed two-fibration selection',flush=True)
def run(check=False):
 p=cert.read(D/'protocol.json');assert all(cert.hashed(ROOT/n)==h for n,h in p['sources'].items()) and cert.hashed(scalar.GP)==p['gp_sha256'];scalar.D=D
 rows=[]
 for row in p['rows']:
  rows.append(scalar.evaluate(row,not check));print(row['id'],rows[-1]['score_units'],flush=True)
 ordered=sorted(rows,key=lambda r:(-r['score_units'],r['model_coefficient_bits'],r['id']));selected=[];families=set()
 for r in ordered:
  if r['family'] not in families:selected.append(r['id']);families.add(r['family'])
  if len(selected)==2:break
 assert len(selected)==2
 out={'schema':'elliptic-curves.retained26-gap-scores.v1','status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'ordering':[r['id'] for r in ordered],'selected_ids':selected,'direct_checks':sum(len(r['direct_checks']) for r in rows),'total_worker_seconds':sum(r['wall_seconds'] for r in rows),'scope':p['selection']+' This is a scheduling comparison on already selected rank26 curves, not an unbiased estimate of score efficacy or rank yield.'}
 if check:assert out==cert.read(OUT)
 else:
  assert not OUT.exists();checkpoint(OUT,out)
 print('PASS selected',selected,'from7 retained equations;',out['total_worker_seconds'],'seconds',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','run','check']);args=p.parse_args();freeze() if args.mode=='freeze' else run(args.mode=='check')
