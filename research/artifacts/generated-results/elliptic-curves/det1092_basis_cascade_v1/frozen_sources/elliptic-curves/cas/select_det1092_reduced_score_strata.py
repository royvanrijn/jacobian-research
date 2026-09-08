#!/usr/bin/env python3
"""Small fixed equation-only panel; training scores and height-matched strata."""
import argparse,sys,math
from pathlib import Path
from fractions import Fraction as Q
import certify_compact_r17_candidates as cert
import score_retained_native19 as scalar
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/det1092-reduced-score-strata-v1'
SOURCE=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json';CHART=ART/'det1092_reduced_parameter_chart_v1/chart-search.json'
def freeze():
 assert not (D/'protocol.json').exists()
 parent=cert.read(SOURCE);chart=cert.read(CHART)['selected']['state'];aa,bb,cc,dd=map(Q,chart['parameter_matrix'])
 def value(v,t):
  def ev(co):
   ans=Q(0)
   for c in reversed(co):ans=ans*t+Q(c)
   return ans
  return ev(v['numerator'])/ev(v['denominator'])
 rows=[];excluded=[];seen=set()
 parameters=sorted({Q(m,n) for n in range(1,9) for m in range(-8,9)},key=lambda t:(max(abs(t.numerator),t.denominator),t.denominator,t.numerator))
 for i,t in enumerate(parameters):
  original=(aa*t+bb)/(cc*t+dd)
  if t==1 or original in [0,1,Q(1009,101)]:excluded.append(dict(parameter=str(t),reason='previous fixed fibre or construction anchor'));continue
  a=[value(v,t) for v in parent['a_invariants']];assert not any(a[:3]);A,B=a[-2:];disc=4*A**3+27*B**2
  if not disc:excluded.append(dict(parameter=str(t),reason='singular'));continue
  j=6912*A**3/disc
  if j in seen:excluded.append(dict(parameter=str(t),reason='internal j duplicate, conservative computational deduplication'));continue
  seen.add(j);jb=abs(j.numerator).bit_length()
  # Exact denominator clearing, without factorization or minimal-model assumptions.
  q=t.denominator;F=A*q**8;G=B*q**12;assert F.denominator==G.denominator==1
  rows.append(dict(id=f'fibre-{i:03}',family='det1092-reduced',parameter=str(t),original_parameter=str(original),model=['0','0','0',str(F),str(G)],model_coefficient_bits=max(abs(F.numerator).bit_length(),abs(G.numerator).bit_length()),j_numerator_bits=jb,j_denominator_bits=j.denominator.bit_length(),height_bin=jb//64))
 assert len(parameters)<=100
 paths=[SOURCE,CHART,Path(__file__).resolve(),Path(scalar.__file__),ART/'blind_factor_free_28_control_v1.json',ART/'det1092_reduced_chart_point_pilot_v2.json']
 checkpoint(D/'protocol.json',dict(schema='elliptic-curves.det1092-reduced-score-strata.v1',sources={str(x.relative_to(ROOT)):cert.hashed(x) for x in paths},rows=rows,parameter_count=len(parameters),excluded=excluded,parameter_height=8,primes=scalar.PRIMES,direct_check_primes=scalar.CHECKS,gp_sha256=cert.hashed(scalar.GP),seconds_per_curve=30,rss_bytes=536870912,maximum_workers=1,gate='The independently certified reduced chart reaches low invariant-height fibres, and the full original27-only49-chart factor-free regression recovered28. The first two fixed new fibres completed98 boxes but were not selected for score. Compare six source-only fibres at comparable arithmetic height and identical point exposure, with three already retained rank26 source gaps in the same following portfolio.',selection='Complete rational |m|,n<=8 panel, at most100 addresses. Omit prior fixed fibres and original construction anchor; internally deduplicate j without public catalogue input. Select the most populated64-bit j-numerator bin among256<=bits<640, ties to smaller bin. Within it order by descending training score then model bits and id. Strong: first two. Moderate: floor(N/3),floor(N/3)+1. Lower: the two smallest SHA256([det1092-lower-fixed-v1,id]) among the bottom third. Require N>=9 and six distinct choices; otherwise stop without enlargement. No section ranks, points, validation scores or catalogue matches select parameters.',score='Exact good-prime short-model scaling before reduction at every training prime5..32749. Sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)). Reuse calibrated scorer with independent direct character-sum checks. Validation primes65537..131071 are neither computed nor read.',following_exposure='Exactly49 factor-free own-subgroup charts per selected new fibre and each of the three remaining retained26 gaps; all maps before any points. Height125000, ten seconds per chart, no rank stop, parameter refill, or automatic next wave. Specialized span is certified before exposure; generic17 is never assumed to survive specialization.'))
 print('FROZEN',len(parameters),'fixed addresses;',len(rows),'equations',flush=True)
def expected(create):
 p=cert.read(D/'protocol.json');assert all(cert.hashed(ROOT/n)==h for n,h in p['sources'].items()) and cert.hashed(scalar.GP)==p['gp_sha256'];scalar.D=D
 rows=[]
 for r in p['rows']:
  rows.append(scalar.evaluate(r,create))
  if create:checkpoint(D/'progress.json',dict(status='RUNNING',rows=rows))
  print('SCORED',r['id'],flush=True)
 bins={k:[r for r in rows if r['height_bin']==k] for k in range(4,10)};chosen_bin=min(bins,key=lambda k:(-len(bins[k]),k));band=sorted(bins[chosen_bin],key=lambda r:(-r['score_units'],r['model_coefficient_bits'],r['id']));N=len(band);assert N>=9
 selected=[dict(band[i],stratum='strong') for i in [0,1]]+[dict(band[i],stratum='moderate') for i in [N//3,N//3+1]]+[dict(r,stratum='lower_fixed') for r in sorted(band[2*N//3:],key=lambda r:digest(['det1092-lower-fixed-v1',r['id']]))[:2]]
 assert len({r['id'] for r in selected})==6
 return dict(status='PASS',protocol_sha256=cert.hashed(D/'protocol.json'),rows=rows,bin_counts={str(k):len(v) for k,v in bins.items()},selected_bin=chosen_bin,within_bin_order=[r['id'] for r in band],selected=selected,total_score_seconds=sum(r['wall_seconds'] for r in rows),direct_check_count=sum(len(r['direct_checks']) for r in rows),boundary=p['selection']+' Scores schedule computation; no rank estimate, efficacy conclusion, or novelty claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['freeze','run','check']);a=ap.parse_args()
 if a.mode=='freeze':freeze()
 else:
  r=expected(a.mode=='run')
  if a.mode=='check':assert r==cert.read(D/'result.json')
  else:assert not (D/'result.json').exists();checkpoint(D/'result.json',r)
  print('PASS strata',[(s['parameter'],s['stratum'],s['j_numerator_bits']) for s in r['selected']],flush=True)
