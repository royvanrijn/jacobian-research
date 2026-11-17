#!/usr/bin/env python3
"""Frozen exhaustive GL3(F2) core quotients and binary lift tests."""
import hashlib
import json
from pathlib import Path
import time
from gl3f2 import group_model,all_subgroups,tuple_orbits,generated,steinberg,cocycles,root_cocycles,matmul
from joint_a5 import core_extension
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'gl3f2_plan.json').read_text());assert plan['frozen_before_execution']
    matrices,table,inverses,orders,autos,dual,inner,core,stdpair=group_model()
    catalogue_path=HERE/'gl3f2_catalogue.json';stamp=time.monotonic()
    if catalogue_path.exists():
        catalogue=json.loads(catalogue_path.read_text())
        assert catalogue['model_sha256']==hashlib.sha256((HERE/'gl3f2.py').read_bytes()).hexdigest()
        triples=[tuple(v) for v in catalogue['generating_triple_orbits']]
        capacity=catalogue['pair_capacity']
    else:
        subgroups=all_subgroups(table);assert len(subgroups)==179
        counts={}
        for arity in [2,3]:
            values={}
            for h in subgroups:values[h]=len(h)**arity-sum(v for k,v in values.items() if k<h)
            counts[arity]=values[subgroups[-1]]
        pairs=tuple_orbits(2,table,orders,autos);triples=tuple_orbits(3,table,orders,autos)
        assert counts[2]==336*len(pairs) and counts[3]==336*len(triples)
        assert len(triples)<=20000
        capacity=len(pairs)
        catalogue={'status':'COMPLETE','matrices':matrices,'subgroups':[sorted(h) for h in subgroups],
                   'std_pair':stdpair,'automorphisms':autos,'dual_automorphism':dual,
                   'pair_capacity':capacity,'generating_tuple_counts':counts,'generating_pair_orbits':pairs,
                   'generating_triple_orbits':triples,'model_sha256':hashlib.sha256((HERE/'gl3f2.py').read_bytes()).hexdigest()}
        catalogue_path.write_text(json.dumps(catalogue,indent=2)+'\n')
        assert time.monotonic()-stamp<120
        print(json.dumps({'catalogue':'COMPLETE','subgroups':len(subgroups),'pair_capacity':capacity,'triple_orbits':len(triples),'seconds':time.monotonic()-stamp}),flush=True)
    eight,building=steinberg(matrices,inverses)
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases'][:16]
    old=json.loads((HERE/'joint_a5_results.json').read_text())['results'][:16]
    assert [r['name'] for r in panel]==plan['targets']
    paths=['gl3f2_plan.json','gl3f2.py','audit_gl3f2.py','gl3f2_catalogue.json','joint_a5_panel.json',
           'joint_a5_results.json','joint_a5_verification.json','joint_a5.py','parabolic_lifts.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    output=HERE/'gl3f2_results.json'
    report=json.loads(output.read_text()) if output.exists() else {'status':'RUNNING','results':[],
        'input_sha256':hashes,'Steinberg_building':building,'Steinberg_matrices':eight}
    assert report['input_sha256']==hashes
    done={r['name'] for r in report['results']};start=time.monotonic()
    for case,source in zip(panel,old):
        if case['name'] in done:continue
        assert time.monotonic()-start<180;stamp=time.monotonic()
        crossings=oriented_diagram(case['pd']);seeds=source['seeds'];trace=source['seed_trace']
        assert source['normalized_arc']==seeds[0]
        covered=set(seeds)
        for o,a,b in trace:assert o in covered and a in covered and b not in covered;covered.add(b)
        assert covered==set(range(13))
        words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in crossings]+[[seeds[0]+1]]
        knotwords=[[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)] for o,a,b,sign in crossings]
        hits=[]
        for orbit,values in enumerate(triples):
            labels=core_extension(13,[(o,a,b) for o,a,b,_ in crossings],seeds,trace,values,core)
            if labels is None:continue
            assert len(hits)<2000
            tests=[]
            for name,model,indices in [('natural',matrices,labels),('dual',matrices,[dual[i] for i in labels]),('Steinberg',eight,labels)]:
                tests.append({'module':name,**cocycles([model[i] for i in indices],[model[inverses[i]] for i in indices],words)})
            implementers=[a for a in autos if all(a[i]==inverses[i] for i in labels)]
            assert len(implementers)<=1
            deck=None;roots=[]
            if implementers:
                a=implementers[0];inner_values=[g for g in range(168) if all(table[table[g][i]][inverses[g]]==a[i] for i in range(168))]
                if inner_values:
                    assert len(inner_values)==1;deck=inner_values[0];assert orders[deck]==2
                    meridians=[table[i][deck] for i in labels]
                    assert all(orders[i]==2 for i in meridians)
                    for label,indices in [('natural',meridians),('dual',[dual[i] for i in meridians])]:
                        roots.append({'module':label,'meridian_images':indices,
                                      **root_cocycles([matrices[i] for i in indices],knotwords)})
            hits.append({'source_group':source['source_group'],'orbit_index':orbit,'seed_images':values,'arc_images':labels,
                         'ordinary_cocycles':tests,'deck_type':'inner' if deck is not None else ('outer' if implementers else 'none'),
                         'deck_element':deck,'root_lifts':roots})
        result={'name':case['name'],'source_group':source['source_group'],'seeds':seeds,'trace':trace,
                'GL3_quotient_kernels':len(hits),'maps':hits,'pair_capacity':capacity,
                'sharp_by_kernel_count':len(hits)>capacity,'sharp_by_ordinary_cocycles':any(t['sharp_cover_rank3'] for h in hits for t in h['ordinary_cocycles']),
                'sharp_by_root_lifts':any(t['sharp_four'] for h in hits for t in h['root_lifts']),
                'seconds':time.monotonic()-stamp}
        report['results'].append(result);output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'name':result['name'],'kernels':len(hits),'decks':[h['deck_type'] for h in hits],
            'cocycles':[[t['Z1_dimension'] for t in h['ordinary_cocycles']] for h in hits],
            'root_dimensions':[[t['dimension'] for t in h['root_lifts']] for h in hits],
            'sharp':any(result[k] for k in ['sharp_by_kernel_count','sharp_by_ordinary_cocycles','sharp_by_root_lifts']),
            'seconds':result['seconds']}),flush=True)
        assert result['seconds']<20
    assert len(report['results'])==16;report['status']='COMPLETE'
    output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
