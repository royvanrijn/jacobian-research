#!/usr/bin/env python3
"""Exact generic/quotient parity rosters, chart maps and complete-cloud provenance."""
import certify_compact_r17_candidates as cert
import followup_scaled13_rank25 as follow
from replay_retention24_geometry import geometry,cloud_check,tuples
from research_runtime.store import digest
IDS=('11952-300',)
def main():
    count=0
    for identifier in IDS:
        follow.configure(identifier);p=follow.protocol();D=follow.D;data=cert.read(D/'result.json');seed=cert.read(follow.SEED);maps=cert.read(D/'maps.json');initial=tuples(seed['final_state']['state']['reductions']['points'])
        census=cert.read(follow.CENSUS);pool=sorted(census['records'][1:],key=lambda r:(-cert.F(r['norm']),r['mask']))[:301]
        if p['generic_pool']!=[{'mask':r['mask'],'computed_generic_norm':r['norm']} for r in pool]:raise ArithmeticError('generic pool selection differs')
        words=sorted(range(1,256),key=lambda w:(w.bit_count(),w));expected={(r['mask'],words[i%len(words)]) for i,r in enumerate(pool)}
        if len(maps['centres'])!=301 or {(c['generic_mask'],c['quotient_word']) for c in maps['centres']}!=expected:raise ArithmeticError('quotient parity roster differs')
        for c in maps['centres']:
            if c['parity']!=c['generic_mask']|(c['quotient_word']<<17) or len(c['representative'])!=25 or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(25)):raise ArithmeticError('exact centre parity differs')
        if len(initial)!=25 or tuples(data['initial_state']['state']['reductions']['points'])!=initial or data['protocol_hash']!=digest(p) or maps['protocol_hash']!=digest(p):raise ArithmeticError('initial subgroup/protocol differs')
        geometry(data,maps,initial,{**p,'maps_path':D/'maps.json'})
        cloud_check(data,seed['charts']+data['charts'],follow.ART/('scaled13_rank25_'+identifier.replace('-','_')+'_all_retained_mod2_v1.json'),D/'all-retained-point-cloud-only.json')
        count+=len(data['charts'])
        if data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=301:raise ArithmeticError('declared complete attempt is short')
    print('EXACT RETENTION RANK25 FOLLOWUP GEOMETRY AND CLOUD PROVENANCE',count,'charts',flush=True)
if __name__=='__main__':main()
