#!/usr/bin/env python3
"""Bounded recorded-exposure audit for the eighteen retained rank26 curves."""
import argparse,json,re,subprocess
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'retained26_source_gap_v1.json'
EXPECTED=['new-20260906-73','new-20260906-63','new-20260906-42','new-20260906-74','new-20260906-49','new-20260906-50','new-20260906-75']
def build():
 assert not OUT.exists();index=ART/'new_high_rank_curve_index_v22.json';rows=[r for r in cert.read(index)['curves'] if r['rank_lower_bound']==26];assert len(rows)==18
 files=sorted(p for depth in [1,2,3] for p in LOCAL.glob('*/'*depth+'result.json') if not any(s in str(p.relative_to(LOCAL)).lower() for s in ['portable','workspace','standalone']))
 pattern='"parameter"\\s*:\\s*"(?:'+'|'.join(re.escape(p) for p in sorted({r['parameter'] for r in rows}))+')"'
 result=subprocess.run(['rg','-l',pattern,*map(str,files)],capture_output=True,text=True);assert result.returncode in [0,1] and not result.stderr
 matched={r['id']:[] for r in rows};paths=[Path(__file__).resolve(),index]
 for name in result.stdout.splitlines():
  path=Path(name);data=cert.read(path)
  for row in rows:
   if data.get('family')==row['family'] and str(data.get('parameter'))==row['parameter']:
    paths.append(path);dims=sorted({len(c.get('centre',{}).get('representative',[])) for c in data.get('charts',[])})
    matched[row['id']].append({'path':str(path.relative_to(ROOT)),'recorded_status':data.get('status'),'recorded_initial_dimension':data.get('initial_dimension'),'centre_word_dimensions':dims,'chart_records':len(data.get('charts',[]))})
 records=[];eligible=[]
 for row in rows:
  source=ART/row['source_certificate'];q=cert.read(source)['curves'][row['source_curve_index']];paths.append(source)
  assert q['curve']==row['curve'] and q['rank_lower_bound']==26 and not row['current_catalogue_matches']
  later=any(any(dim>=26 for dim in r['centre_word_dimensions']) for r in matched[row['id']]);entry={k:row[k] for k in ['id','family','parameter','source_certificate','source_curve_index']};entry.update(records=matched[row['id']],recorded_own26_exposure=later)
  if not later:
   assert len(matched[row['id']])==1 and len(q['generic_points'])==17 and q['points'][:17]==q['generic_points'] and len(q['points'])==26
   witness=q['discovery_witness'];path=ROOT/witness['path'];data=cert.read(path);assert cert.hashed(path)==witness['sha256'] and data['status']=='COMPLETE_DECLARED_POINT_ATTEMPT'
   for c in data['charts']:
    rep=c['centre']['representative'];word=c['search']['input']['centre']['coefficients'];assert len(rep)==17 and word[:17]==rep and not any(word[17:]) and c['search']['status']=='bounded_search_complete'
   proof=q['rank_certificate'];assert json.loads(json.dumps(checked_rank(tuple(map(cert.F,q['curve'])),[tuple(map(cert.F,p)) for p in q['points']],[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])))==proof
   entry.update(generic_prefix_rank=17,unused_independent_directions_in_source_centres=9,completed_source_boxes=len(data['charts']),eligible_for_fixed_comparison=True);eligible.append(row['id'])
  else:entry['eligible_for_fixed_comparison']=False
  records.append(entry)
 assert eligible==EXPECTED
 output={'schema':'elliptic-curves.retained26-source-gap.v1','status':'PASS','inventory_rank26_count':18,'eligible_ids':eligible,'rows':records,'layout_search':{'result_files_examined':len(files),'matched_parameter_files':len(result.stdout.splitlines()),'relative_files':[str(p.relative_to(LOCAL)) for p in files],'depths':[1,2,3],'excluded_path_terms':['portable','workspace','standalone'],'matching':'exact stored family and parameter'},
 'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'scope':'Seven complete retained discovery runs use only17 generic directions despite26-point final proofs. The bounded live-artifact layout audit locates later own26 exposure for eleven other retained curves and excludes them. This is not a repository-wide absence theorem across renamed parameters, isomorphic models, deeper directories or archives. Recorded RUNNING labels are historical data, not claims of currently live processes or completed exposure. No duplicate search or automatic population expansion is authorized by this audit.'}
 OUT.write_text(json.dumps(output,indent=2)+'\n');print('PASS seven source gaps; eleven with later own26 records; searched',len(files),'result paths',flush=True)
def check():
 d=cert.read(OUT);assert d['status']=='PASS' and d['eligible_ids']==EXPECTED
 for n,h in d['sources'].items():assert cert.hashed(ROOT/n)==h
 # Replay the affirmative source-gap witnesses; do not infer global absence
 # from the historical layout snapshot or rescan changing directories.
 for row in d['rows']:
  if not row['eligible_for_fixed_comparison']:continue
  q=cert.read(ART/row['source_certificate'])['curves'][row['source_curve_index']];raw=cert.read(ROOT/q['discovery_witness']['path'])
  assert len(q['points'])==26 and len(q['generic_points'])==17
  for c in raw['charts']:assert len(c['centre']['representative'])==17 and not any(c['search']['input']['centre']['coefficients'][17:]) and c['search']['status']=='bounded_search_complete'
 print('PASS retained26 source-gap bindings and seven exact centre-prefix audits')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else build()
