#!/usr/bin/env python3
"""Full pinned uncertainty panel for the singular affine architecture."""
import hashlib
import json
from pathlib import Path
import time
from gl3f2 import group_model,generated,cocycles,root_cocycles,matmul,matvec,rank,kernel,xor
from joint_a5 import core_extension
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def torus_control(matrices,table,inverses,orders):
    x,y=next((x,y) for x in range(168) for y in range(168)
             if orders[x]==4 and orders[y]==7
             and orders[table[inverses[x]][table[y][y]]]==2
             and len(generated([x,y],table))==168)
    xx,yy=matrices[x],matrices[y];mx=matrices[inverses[x]]
    m=matmul(mx,matmul(yy,yy));diff=tuple(row^(1<<j) for j,row in enumerate(m))
    root=next(sum(((row>>k)&1)<<j for j,row in enumerate(diff)) for k in range(3) if any(row>>k&1 for row in diff))
    def homogeneous(a,v):return tuple(row|(((v>>j)&1)<<3) for j,row in enumerate(a))+(8,)
    def power(a,n):
        out=tuple(1<<j for j in range(len(a)))
        for _ in range(n):out=matmul(out,a)
        return out
    meridian_vectors=[]
    for j in range(6):
        a,b=homogeneous(xx,(1<<j) if j<3 else 0),homogeneous(yy,(1<<(j-3)) if j>=3 else 0)
        meridian=matmul(power(a,3),power(b,2))
        meridian_vectors.append(sum(((r>>3)&1)<<i for i,r in enumerate(meridian[:3])))
        assert power(a,4)==power(b,7)==(1,2,4,8)
    functionals=[f for f in range(1,8) if (f&root).bit_count()%2==0]
    rows=[sum(((f&v).bit_count()%2)<<j for j,v in enumerate(meridian_vectors)) for f in functionals]
    basis=kernel(rows,6);assert len(basis)==4
    principal=[]
    for v in [1,2,4]:principal.append((matvec(xx,v)^v)|((matvec(yy,v)^v)<<3))
    assert rank(principal,6)==3
    chosen=next(v for v in basis if rank(principal+[v],6)==4)
    a,b=homogeneous(xx,chosen&7),homogeneous(yy,chosen>>3)
    mu=matmul(power(a,3),power(b,2))
    assert rank([r^(1<<j) for j,r in enumerate(mu)],4)==1
    assert rank([r^(1<<j) for g in [a,b] for j,r in enumerate(g)],4)==4
    found={(1,2,4,8)};queue=[(1,2,4,8)]
    for g in queue:
        for h in [a,b]:
            z=matmul(g,h)
            if z not in found:found.add(z);queue.append(z);assert len(found)<=1344
    assert len(found)==1344
    invariant_forms=[]
    for bits in range(1,8):
        f=( ((bits&1)<<1)|(((bits>>1)&1)<<2), (bits&1)|(((bits>>2)&1)<<2), ((bits>>1)&1)|(((bits>>2)&1)<<1) )
        from gl3f2 import transpose
        if all(matmul(matmul(transpose(g),f),g)==f for g in [xx,yy]):invariant_forms.append(f)
    assert not invariant_forms
    return {'knot':'T(4,7)','X_index':x,'Y_index':y,'meridian_index':table[inverses[x]][table[y][y]],
            'admissibility_rows':rows,'kernel_basis':basis,'principal_basis':principal,'chosen_cocycle':chosen,
            'homogeneous_X':a,'homogeneous_Y':b,'homogeneous_meridian':mu,
            'affine_group_order':len(found),'linear_invariant_nonzero_alternating_forms':invariant_forms,
            'meridian_difference_rank':1,'global_fixed_codimension':4,'admissible_cocycle_dimension':4}


def main():
    plan=json.loads((HERE/'gl3f2_full_panel_plan.json').read_text());assert plan['frozen_before_execution']
    matrices,table,inverses,orders,autos,dual,_,core,_=group_model()
    catalogue=json.loads((HERE/'gl3f2_catalogue.json').read_text())
    triples=catalogue['generating_triple_orbits'];capacity=catalogue['pair_capacity']
    first=json.loads((HERE/'gl3f2_results.json').read_text());assert first['status']=='COMPLETE'
    eight=[tuple(a) for a in first['Steinberg_matrices']]
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases']
    source=json.loads((HERE/'joint_a5_results.json').read_text())['results']
    assert len(panel)==len(source)==293
    assert [r['name'] for r in panel]==json.loads((HERE/'sources/current13/data/still_open.json').read_text())
    paths=['gl3f2_full_panel_plan.json','audit_gl3f2_full_panel.py','gl3f2.py','gl3f2_catalogue.json','gl3f2_results.json',
           'joint_a5_panel.json','joint_a5_results.json','joint_a5.py','parabolic_lifts.py',
           'sources/current13/search/glsearch.py','sources/current13/search/w3.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    output=HERE/'gl3f2_full_panel_results.json'
    report=json.loads(output.read_text()) if output.exists() else {'status':'RUNNING','results':list(first['results']),
               'input_sha256':hashes,'T47_control':torus_control(matrices,table,inverses,orders)}
    assert report['input_sha256']==hashes
    done={r['name'] for r in report['results']};start=time.monotonic()
    for case,old in zip(panel,source):
        if case['name'] in done:continue
        assert time.monotonic()-start<180;stamp=time.monotonic()
        crossings=oriented_diagram(case['pd']);seeds=old['seeds'];trace=old['seed_trace']
        assert len(seeds)==4 and seeds[0]==old['normalized_arc']
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
            assert len(implementers)<=1;deck=None;roots=[]
            if implementers:
                action=implementers[0]
                candidates=[g for g in range(168) if all(table[table[g][i]][inverses[g]]==action[i] for i in range(168))]
                if candidates:
                    assert len(candidates)==1;deck=candidates[0];assert orders[deck]==2
                    meridians=[table[i][deck] for i in labels];assert all(orders[i]==2 for i in meridians)
                    for name,indices in [('natural',meridians),('dual',[dual[i] for i in meridians])]:
                        roots.append({'module':name,'meridian_images':indices,**root_cocycles([matrices[i] for i in indices],knotwords)})
            hits.append({'source_group':old['source_group'],'orbit_index':orbit,'seed_images':values,'arc_images':labels,
                         'ordinary_cocycles':tests,'deck_type':'inner' if deck is not None else ('outer' if implementers else 'none'),
                         'deck_element':deck,'root_lifts':roots})
        result={'name':case['name'],'source_group':old['source_group'],'seeds':seeds,'trace':trace,
                'GL3_quotient_kernels':len(hits),'maps':hits,'pair_capacity':capacity,'sharp_by_kernel_count':len(hits)>capacity,
                'sharp_by_ordinary_cocycles':any(t['sharp_cover_rank3'] for h in hits for t in h['ordinary_cocycles']),
                'sharp_by_root_lifts':any(t['sharp_four'] for h in hits for t in h['root_lifts']),
                'seconds':time.monotonic()-stamp}
        report['results'].append(result);output.write_text(json.dumps(report,indent=2)+'\n')
        if len(report['results'])%25==0 or result['sharp_by_root_lifts'] or result['sharp_by_ordinary_cocycles']:
            print(json.dumps({'completed':len(report['results']),'last':case['name'],'last_kernels':len(hits),
                             'sharp_roots':result['sharp_by_root_lifts'],'sharp_ordinary':result['sharp_by_ordinary_cocycles']}),flush=True)
        assert result['seconds']<20
    assert len(report['results'])==293
    report['status']='COMPLETE';output.write_text(json.dumps(report,indent=2)+'\n')
    from collections import Counter
    print(json.dumps({'status':'COMPLETE','kernel_histogram':dict(Counter(r['GL3_quotient_kernels'] for r in report['results'])),
          'deck_types':dict(Counter(m['deck_type'] for r in report['results'] for m in r['maps'])),
          'root_dimensions':dict(Counter(t['dimension'] for r in report['results'] for m in r['maps'] for t in m['root_lifts'])),
          'sharp_names':[r['name'] for r in report['results'] if any(r[k] for k in ['sharp_by_kernel_count','sharp_by_ordinary_cocycles','sharp_by_root_lifts'])],
          'T47_control_group_order':report['T47_control']['affine_group_order']},indent=2))


if __name__=='__main__':main()
