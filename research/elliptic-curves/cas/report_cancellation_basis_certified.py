#!/usr/bin/env sage -python
"""Post-protocol report adapter; the frozen comparison and all workers are unchanged."""
import json
from cancellation_basis_certified import OUT, RAW, read, sha, guard, write

def report():
    from cancellation_cloud_programme import compare
    from retain_cancellation_basis_certified_v2 import bindings
    audit=bindings(); plan,inputs=guard(); cases={c['id']:c for c in inputs}; rows=[]
    for receipt in read(OUT/'supervision.json')['records']:
        r=read(RAW/'arms'/receipt['case']/receipt['arm']/'result.json')
        rows.append({**{k:r[k] for k in ('case','family','stratum','arm','success','new_directions','later_cloud_directions','first_cloud_gain','unknowns','status')},
            'validation_block':cases[r['case']]['validation_block'],'calls':len(r['calls']),'epochs':len(r['epochs']),
            'cpu_seconds':receipt['charged_cpu_seconds'],'components':r['components']})
    def comparison(selected,later):
        renamed=[{**r,'arm':'factor_free' if r['arm']=='fixed_bank' else 'adaptive_cloud',
            'new_directions':r['later_cloud_directions'] if later else r['new_directions']} for r in selected]
        totals,outcome=compare(renamed)
        return {'totals':{('fixed_bank' if k=='factor_free' else 'basis_refresh'):v for k,v in totals.items()},'comparison':outcome}
    all_gain=comparison(rows,False); later=comparison(rows,True)
    blocks={str(b):{'all':comparison([r for r in rows if r['validation_block']==b],False),
                   'later':comparison([r for r in rows if r['validation_block']==b],True)} for b in (0,1)}
    unknowns=sum(r['unknowns'] for r in rows); a,b=all_gain['totals']['basis_refresh'],all_gain['totals']['fixed_bank']
    performance=bool(later['comparison']['gate_passed'] and a['directions']>=b['directions'] and a['target_completions']>=b['target_completions'] and not unknowns)
    repeat=all(x['later']['comparison']['uncapped_direction_rate_ratio'] is not None and x['later']['comparison']['uncapped_direction_rate_ratio']>1 for x in blocks.values())
    mechanism=read(OUT/'mechanism.json'); representatives=read(OUT/'representative-audit.json')
    causal=bool(mechanism['status']=='PASS' and representatives['status']=='PASS_FINITE_REPRESENTATIVE_SENSITIVITY' and
        all(representatives['blocks'][str(i)]['both']>0 for i in (0,1)))
    shared=read(OUT/'preflight.json')['cpu_seconds']; diagnostic=sum(read(OUT/n)['full_child_cpu_seconds'] for n in ('mechanism-process.json','representative-process.json'))
    development=read(OUT/'integration-prerequisite.json')['summary']
    dev_cpu=development['full_child_cpu_seconds']+development['negative_check_component_cpu_seconds']
    la,lb=later['totals']['basis_refresh'],later['totals']['fixed_bank']
    def ratio(extra):return (la['directions']/(la['cpu_seconds']+extra))/(lb['directions']/lb['cpu_seconds']) if lb['directions'] else None
    summary={'status':'COMPLETE_FIXED_TWO_BLOCK_CERTIFIED_BANK_CONTROL','rows':rows,'all':all_gain,'later':later,'blocks':blocks,
        'aggregate_performance_gate':performance,'repeatability_rate_gate':repeat,'actual_mechanism_gate':causal,
        'promotion_gate':performance and repeat and causal,'preparation_unknowns':unknowns,
        'shared_preflight_component_cpu_seconds':shared,'postmortem_full_child_cpu_seconds':diagnostic,
        'transcript_development_recorded_cpu_seconds':dev_cpu,
        'conservative_later_rate_ratio':ratio(shared),'conservative_ratio_including_postmortem':ratio(shared+diagnostic),
        'conservative_ratio_including_recorded_integration':ratio(shared+diagnostic+dev_cpu),
        'exposure':{k:audit[k] for k in ('totals','paired_initial_exposure','point_call_statuses')},
        'protocol_sha256':sha(OUT/'protocol.json'),'fresh_fibres_run':0,'rank32':'UNKNOWN','boundary':plan['boundary']}
    write(OUT/'summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('all','later','aggregate_performance_gate','repeatability_rate_gate','actual_mechanism_gate','promotion_gate','preparation_unknowns')},indent=2),flush=True)


if __name__=='__main__':
    report()
