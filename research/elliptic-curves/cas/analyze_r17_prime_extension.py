#!/usr/bin/env python3
"""Descriptive held-out-prime and frozen prior-exposure diagnostics."""
from pathlib import Path
import statistics
import certify_compact_r17_candidates as cert
from compare_bounded_prime_selectors import correlation,ranks
import extend_retained_r17_prime_scores as extension
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'r17_retained_prime_extension_diagnostics_v1.json';OLD=ART/'compact_r17_wide_results_v1.json'

def main():
    if OUT.exists():raise FileExistsError('preserve prime diagnostics')
    protocol=extension.protocol();data=cert.read(extension.D/'result.json');gate=cert.read(extension.D/'replay.supervisor.json')
    if data['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or gate['outcome']!='completed' or gate['returncode']!=0:raise ArithmeticError('trace replay incomplete')
    old=cert.read(OLD);rows=[];by={(r['family'],r['parameter']):r['rank_lower_bound'] for r in old['curves']}
    for f,choices in data['selection'].items():
        current=[r for r in data['rows'] if r['family']==f];remaining=[r for r in current if r['retained_index']>=4];index={r['retained_index']:r for r in current};vr=dict(zip([r['retained_index'] for r in remaining],ranks([r['validation_units'] for r in remaining])));arms={}
        for name,ids in choices.items():
            arms[name]={'indices':ids,'parameters':[index[i]['parameter'] for i in ids],'mean_validation_score':statistics.mean(index[i]['validation_units']/10**12 for i in ids),'mean_validation_percentile':statistics.mean(vr[i]/(len(remaining)-1) for i in ids)}
        measured=[r for r in current if r['retained_index']<4]
        if len(measured)!=4 or any((f,r['parameter']) not in by for r in measured):raise ArithmeticError('prior fixed exposure roster differs')
        lower=[by[(f,r['parameter'])] for r in measured]
        rows.append({'family':f,'arms':arms,'validation_score_difference_extended_minus_original':arms['extended_top_two']['mean_validation_score']-arms['original_next_two']['mean_validation_score'],'prior_four_lower_bounds':lower,'prior_four_spearman_with_original_score':correlation([r['score_units'] for r in measured],lower),'prior_four_spearman_with_extended_score':correlation([r['combined_selection_units'] for r in measured],lower),'remaining_pool_selection_vs_validation_spearman':correlation([r['combined_selection_units'] for r in remaining],[r['validation_units'] for r in remaining])})
    checkpoint(OUT,{'schema':'elliptic-curves.r17-prime-extension-diagnostics.v1','source_sha256':cert.hashed(Path(__file__).resolve()),'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),OLD,extension.D/'protocol.json',extension.D/'result.json',extension.D/'replay.supervisor.json')},'rows':rows,'candidate_count':len(data['rows']),'additional_prime_count':len(protocol['prime_roster']),'trace_lookups':len(data['rows'])*len(protocol['prime_roster']),'families_with_higher_extended_arm_validation':sum(r['validation_score_difference_extended_minus_original']>0 for r in rows),'claim_boundary':'Descriptive finite diagnostics after selection and paired point roster were frozen. Disjoint validation primes are not actual ranks. The24 previous lower bounds came from a different point policy and a truncated score pool; four per family is not a general rank-classifier validation. No policy promotion, p-value or selection completeness claim.'})
    print('EXTENDED PRIME DIAGNOSTICS',len(rows),'families; validation wins',sum(r['validation_score_difference_extended_minus_original']>0 for r in rows),'of6',flush=True)
if __name__=='__main__':main()
