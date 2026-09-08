#!/usr/bin/env python3
"""Bind exact coefficient bounds, skew-box retention and terminal point witnesses."""
import argparse,json
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_r17_parameter_box_skew as bound
import scan_skew_r17_boxes as population
import extend_skew_r17_scores as extension
import skew8_r17_pari_batch as batch
import audit_higher24_visibility_cost as costs
from research_runtime.store import checkpoint
ROOT=batch.ROOT;ART=batch.ART;OUT=ART/'skew_r17_experiment_v1.json'

def expected():
    paths={Path(__file__).resolve(),Path(costs.__file__).resolve(),Path(costs.minimal.__file__).resolve()}
    def read(p):paths.add(p);return cert.read(p)
    def terminal(p):
        s=read(p)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required terminal replay failed or censored')
    b=read(bound.OUT)
    if b!=bound.expected():raise ArithmeticError('exact rectangle bound differs')
    pop=read(population.D/'result.json');scores=read(extension.D/'result.json');ep=read(extension.D/'protocol.json');p=read(batch.D/'protocol.json');v=read(batch.D/'verification-ledger.json');o=read(batch.D/'odd-cloud-audit/ledger.json');proof=read(ART/'skew8_r17_results_v1.json')
    for folder,label in [(population.D,'replay'),(extension.D,'benchmark-check'),(extension.D,'replay'),(batch.D/'post-batch','certify'),(batch.D/'post-batch','proof-replay'),(batch.D/'post-batch','geometry')]:terminal(folder/(label+'.supervisor.json'))
    if pop['status']!='COMPLETE_FROZEN_SKEW_POPULATION' or scores['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or len(pop['rows'])!=2048 or len(scores['rows'])!=2048 or v['status']!='PASS' or o['status']!='PASS':raise ArithmeticError('complete finite population and point checks required')
    total=sum(r['summary'][3] for r in pop['shards']);small=sum(r['summary'][3] for r in pop['small_checks'])
    if (total,small,pop['old_S1_retained_overlap'],ep['cached_rows'],ep['new_rows'])!=(40733526,2634,654,656,1392):raise ArithmeticError('fixed population/cache counts differ')
    vk={r['id']:r for r in v['rows']};ok={r['id']:r for r in o['rows']};qk={r['id']:r for r in proof['curves']};rows=[]
    if len(qk)!=8 or set(vk)!=set(qk) or set(ok)!=set(qk):raise ArithmeticError('fixed eight-curve roster differs')
    for row in p['rows']:
        name=row['id'];d=batch.D/name;raw=read(d/'result.json');q=qk[name];odd=ok[name];cloud=read(ROOT/vk[name]['cloud_certificate']);ell=read(ROOT/odd['output']);rank=q['rank_lower_bound'];attempted=len(raw['charts']);declared=len(p['generic_masks'][row['family']]);completed=sum(c['search']['status']=='bounded_search_complete' for c in raw['charts'])
        if abs(cert.F(row['parameter']).numerator)<=32768 or (q['family'],q['parameter'])!=(row['family'],row['parameter']):raise ArithmeticError('new rectangle strip differs')
        if cert.hashed(d/'result.json')!=vk[name]['input_sha256'] or cert.hashed(ROOT/vk[name]['cloud_certificate'])!=vk[name]['cloud_sha256'] or cert.hashed(ROOT/odd['output'])!=odd['output_sha256']:raise ArithmeticError('bound cloud inputs changed')
        if cloud['rank_lower_bound']!=rank or odd['mod2_lower_bound']!=rank or {str(a['modulus']):a['finite_column_rank'] for a in ell['audits']}!=odd['odd_lower_bounds']:raise ArithmeticError('complete-cloud bounds differ')
        if any(c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('uniform point budget differs')
        if raw['status']=='COMPLETE_DECLARED_POINT_ATTEMPT':
            if attempted!=declared:raise ArithmeticError('short complete attempt')
        elif raw['status']=='TARGET_REACHED_PENDING_REPLAY':
            if raw['rank_lower_bound']<28:raise ArithmeticError('target stop without bound')
        else:raise ArithmeticError('nonterminal point status')
        rows.append({'id':name,'family':row['family'],'parameter':row['parameter'],'rank_lower_bound':rank,'attempted_charts':attempted,'declared_charts':declared,'completed_boxes':completed,'retained_points':len(cloud['points']),'odd_modulus_lower_bounds':odd['odd_lower_bounds'],'catalogue_matches':q['icarm_matches'],'previous_matches':q['previous_matches'],'model_size':costs.size(q['curve'])})
    if sum(r['declared_charts'] for r in rows)!=368:raise ArithmeticError('fixed368 chart policy differs')
    earlier=read(ART/'higher24_r17_results_v1.json');oldcosts=[{'id':r['id'],'family':r['family'],'parameter':r['parameter'],**costs.size(r['curve'])} for r in earlier['curves'] if r['family'] in ('08234','08f72')]
    if len(oldcosts)!=8:raise ArithmeticError('fixed two-family old cost comparator differs')
    cost_comparison={'old_square_rows':oldcosts,'old_square_statistics':costs.statistics(oldcosts,'normalized_coefficient_bits'),'skew_statistics':costs.statistics([r['model_size'] for r in rows],'normalized_coefficient_bits'),'boundary':'Post-selection descriptive comparison of the eight earlier and eight skew candidates in the same two families. Normalization is at the existing search scale; global minimality is asserted only where the exact bounded invariant proof closes it. Different selected populations confound causality.'}
    return {'schema':'elliptic-curves.skew-r17-experiment.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'observed_model_size_comparison':cost_comparison,'coefficient_bound_families':[{k:r[k] for k in ('family','selected_k','selected_numerator_bound','selected_denominator_bound','weighted_bound_improvement')} for r in b['rows'] if r['selected_k']],'primitive_addresses':total,'small_exhaustive_scores':small,'retained_addresses':2048,'outside_old_square':sum(abs(r['numerator'])>32768 for r in pop['rows']),'old_S1_retained_overlap':654,'cached_trace_rows':656,'new_trace_rows':1392,'fresh_extension_traces':1392*5978,'reused_extension_traces':656*5978,'rows':rows,'rank_lower_bound_counts':dict(sorted(Counter(r['rank_lower_bound'] for r in rows).items())),'declared_point_boxes':368,'attempted_point_boxes':sum(r['attempted_charts'] for r in rows),'completed_point_boxes':sum(r['completed_boxes'] for r in rows),'within_batch_isomorphic_pairs':proof['within_batch_isomorphic_pairs'],'catalogue_snapshot_curves':proof['catalogue']['curve_count'],'previous_address_equations':len(proof['previous_equations']),'claim_boundary':'A finite equal-area rectangle optimization for conservative coefficient bounds, followed by unchanged S1 retention and point exposure only on addresses outside the earlier square. These are parameter-population and lower-bound measurements, not minimal discriminant, conductor, rank-density, causal model-size or selector-superiority claims. Numerical search geometry is not an optimality proof. New inventory entries require separate point-proof and distinctness promotion. No exact rank, rank upper bound, saturation, point absence or universal novelty.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('skew experiment aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve skew experiment aggregate')
        checkpoint(OUT,d)
    print('SKEW EXPERIMENT PASS',d['rank_lower_bound_counts'],d['completed_point_boxes'],'completed point boxes',flush=True)
