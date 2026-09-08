#!/usr/bin/env python3
"""Post-terminal comparison and presentation ledger, never worker input."""
import collections,gzip,hashlib,json,shutil,time,statistics
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1';PKG=ART/'det1092_basis_cascade_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(d):return hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def common_prefix(a,b):return next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),min(len(a),len(b)))
roster=read(PKG/'roster.json');commit=roster['starting_commit']
def write(p,d):
 with p.open('x') as f:json.dump({'starting_commit':commit,**d},f,sort_keys=True,indent=2);f.write('\n')
assert all((D/r['arm_id']/'supervisor-result.json').exists() for r in roster['arms'])
# This historical comparison is strictly after every arm is terminal.
historical={};historical_points=set();comparison_sources={}
for name in ['adaptive-visibility-cascade-v1','adaptive-visibility-cascade-v2','adaptive-visibility-cascade-v3']:
 base=ROOT/'artifacts/local/elliptic-curves'/name
 for path in sorted(base.rglob('chart-*.json')):
  raw=read(path)
  if not isinstance(raw,dict) or 'mapping' not in raw or 'search' not in raw:continue
  if not all(k in raw['mapping'] for k in ['reduced_P','reduced_Q']):continue
  key=digest({k:raw['mapping'][k] for k in ['reduced_P','reduced_Q']})
  historical.setdefault(key,[]).append(str(path.relative_to(ROOT)))
  comparison_sources[str(path.relative_to(ROOT))]=sha(path)
  for p in raw['search'].get('finite_curve_points',[]):historical_points.add((tuple(raw['search']['short_model']),p['x'],str(abs(F(p['y'])))))
new_equations={};summary=[];bundles={};raw_bindings={};new_points={}
for row in roster['arms']:
 aid=row['arm_id'];arm=D/aid;supervisor=read(arm/'supervisor-result.json');terminal=read(arm/'run/terminal.json') if (arm/'run/terminal.json').exists() else {};replay=read(PKG/f'replay-{aid}.json') if (PKG/f'replay-{aid}.json').exists() else None
 if supervisor['exit_code']==0:assert replay and replay['status']=='PASS_EXACT_REPLAY'
 coefficient_bits=[];maximum_certified_rank=17;stages=[];search_status=collections.Counter();charts_total=0;search_wall=0;point_records=0;path=[];bundle={'starting_commit':commit,'arm_id':aid,'fixture':read(arm/'fixture.json'),'protocol':read(arm/'protocol.json'),'terminal':terminal,'supervisor':supervisor,'resources':read(arm/'resources.json') if (arm/'resources.json').exists() else None,'stages':[]}
 for wd in sorted((arm/'run').glob('epoch-*')):
  charts=[read(p) for p in sorted(wd.glob('chart-*.json'))];stage=read(wd/'stage.json') if (wd/'stage.json').exists() else None
  selection=read(wd/'selection.json') if (wd/'selection.json').exists() else None
  packed={'epoch':int(wd.name.split('-')[1]),'stage':stage,'charts':charts,'selected_centres':selection['centres'] if selection else None,'metric_transport':read(wd/'metric-transport.json') if (wd/'metric-transport.json').exists() else None}
  if stage:
   audit=read(wd/stage['audit']);maximum_certified_rank=max(maximum_certified_rank,audit['rank_lower_bound']);packed['audit']=audit;packed['modl']=read(wd/'modl.json') if (wd/'modl.json').exists() else None;stages.append(stage)
   if packed['modl']:maximum_certified_rank=max([maximum_certified_rank]+[v['finite_column_rank'] for v in packed['modl']['audits']])
   packed['final_state']=read(wd/f"cloud-{stage['charts']-1:03d}.json")['final_state']
  bundle['stages'].append(packed)
  for c in charts:
   coefficient_bits.append(c['search']['maximum_coefficient_bits']);charts_total+=1;search_status[c['search']['status']]+=1;search_wall+=c['search']['wall_seconds'];point_records+=len(c['search']['finite_curve_points']);path.append(c['centre']['point'])
   mapping={k:c['mapping'][k] for k in ['reduced_P','reduced_Q']};key=digest(mapping)
   ledger=new_equations.setdefault(key,{'equation_id':key,'equation':'v^2+Q(u)*v=P(u)','coefficients_ascending':mapping,'classification':'EXACT_Q_BIRATIONAL_PRESENTATION_OF_KNOWN_CURVE302','new_exact_equation_in_compared_chart_catalogue':key not in historical,'historical_matches':historical.get(key,[])[:4],'occurrences':[],'map_witness':{k:c['mapping'][k] for k in ['raw_coefficients','matrix','square_ratio','discriminant_quartic']},'centre':c['centre']['point']})
   ledger['occurrences'].append({'arm':aid,'epoch':packed['epoch'],'chart':c['index']})
   for point in c['search']['finite_curve_points']:
    keyp=(tuple(c['search']['short_model']),point['x'],str(abs(F(point['y']))))
    if keyp not in historical_points:
     canonical={'x':point['x'],'y':str(abs(F(point['y'])))}
     rec=new_points.setdefault(digest(canonical),{'point_modulo_sign':canonical,'classification':'ABSENT_MODULO_SIGN_FROM_COMPARED_SAME_CURVE_CHART_TRANSCRIPTS; SUBGROUP_NOVELTY_NOT_INFERRED','occurrences':[]});rec['occurrences'].append({'arm':aid,'epoch':packed['epoch'],'chart':c['index']})
 bundles[aid]=bundle
 data=json.dumps(bundle,sort_keys=True,separators=(',',':')).encode();target=PKG/f'evidence-{aid}.json.gz'
 with target.open('xb') as f:f.write(gzip.compress(data,mtime=0))
 for name in ['fixture.json','protocol.json','resources.json','supervisor-result.json','failure.json','launch.json']:
  if (arm/name).exists():
   q=PKG/'arms'/aid/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(arm/name,q)
 summary.append({'arm_id':aid,'status':supervisor['classification'],'final_rank_lower_bound':terminal.get('final_rank_lower_bound'),'maximum_retained_cloud_certified_rank':maximum_certified_rank,'above_known31_review_required':maximum_certified_rank>31,'rank_path':[17]+[s['after'] for s in stages],'charts_per_stage':[s['charts'] for s in stages],'charts':charts_total,'selected_charts':sum(len(e['selected_centres'] or []) for e in bundle['stages']),'selected_but_unsearched_charts':sum(len(e['selected_centres'] or [])-len(e['charts']) for e in bundle['stages']),'full_cosets_scored_in_completed_epochs':sum(s['full_cosets_scored'] for s in stages),'search_status_counts':dict(search_status),'point_records':point_records,'chart_maximum_coefficient_bits_median':statistics.median(coefficient_bits) if coefficient_bits else None,'chart_maximum_coefficient_bits_range':[min(coefficient_bits),max(coefficient_bits)] if coefficient_bits else None,'search_wall_seconds':search_wall,'resources':bundle['resources'],'replay_status':replay['status'] if replay else 'UNREPLAYED','evidence_sha256':sha(target),'chart_centre_path_sha256':digest(path)})
 for p in arm.rglob('*'):
  if p.is_file():raw_bindings[str(p.relative_to(ROOT))]={'sha256':sha(p),'bytes':p.stat().st_size}
comparisons=[];base=bundles['original']
for aid,bundle in bundles.items():
 result={'arm_id':aid,'epochs':[],'first_selected_chart_divergence':None,'first_tested_chart_divergence':None,'first_admitted_point_divergence':None}
 for i,(a,b) in enumerate(zip(base['stages'],bundle['stages'])):
  s0=[r['point'] for r in a['selected_centres'] or []];s1=[r['point'] for r in b['selected_centres'] or []];t0=[c['centre']['point'] for c in a['charts']];t1=[c['centre']['point'] for c in b['charts']]
  x=set(map(tuple,s0));y=set(map(tuple,s1));same_extra=a.get('audit',{}).get('independent_points',[])[17:]==b.get('audit',{}).get('independent_points',[])[17:]
  result['epochs'].append({'epoch':i,'ordered_selected_centres_equal':s0==s1,'selected_common_prefix_length':common_prefix(s0,s1),'selected_set_intersection':len(x&y),'selected_set_union':len(x|y),'ordered_tested_centres_equal':t0==t1,'tested_common_prefix_length':common_prefix(t0,t1),'admitted_extension_points_equal':same_extra})
  for key,same in [('first_selected_chart_divergence',s0==s1),('first_tested_chart_divergence',t0==t1),('first_admitted_point_divergence',same_extra)]:
   if not same and result[key] is None:result[key]=i
 comparisons.append(result)
write(PKG/'presentation_ledger.json',{'status':'POST_RUN_EXACT_PRESENTATION_CLASSIFICATION','equations':list(new_equations.values()),'distinct_presentations':len(new_equations),'new_equation_presentations_in_compared_catalogue':sum(r['new_exact_equation_in_compared_chart_catalogue'] for r in new_equations.values()),'new_elliptic_curve_isomorphism_classes':0,'explanation':'Every tested equation is connected by exact nonsingular PGL2/square scaling to a pointed chord quartic of302, birational to302. No parameter sweep or new Jacobian equation producer was run. Catalogue absence establishes equation-text novelty only.','historical_comparison_source_bindings':comparison_sources})
write(PKG/'point_novelty_ledger.json',{'points_absent_modulo_sign':list(new_points.values()),'count':len(new_points),'boundary':'Absence modulo sign from V1/V2/V3 chart transcripts on the same short curve only, not a new Mordell-Weil direction or global novelty.'})
write(PKG/'raw_archive.json',{'files':raw_bindings,'total_bytes':sum(x['bytes'] for x in raw_bindings.values())})
write(PKG/'summary.json',{'status':'TERMINAL_POST_RUN_EVALUATION','arms':summary,'comparisons_to_original':comparisons,'presentation_ledger_sha256':sha(PKG/'presentation_ledger.json'),'point_novelty_ledger_sha256':sha(PKG/'point_novelty_ledger.json'),'raw_archive_sha256':sha(PKG/'raw_archive.json'),'package_script_sha256':sha(Path(__file__)),'claim_boundary':'Fixed curve302 visibility comparison of four unimodular presentations of the same generic17 subgroup. All displayed ranks are certified lower bounds. No new generic rank, new elliptic curve, upper rank bound, or prospective family-transfer theorem.'})
print(json.dumps(summary,indent=2))
