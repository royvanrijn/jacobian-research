#!/usr/bin/env python3
"""Adaptive equation-only relation completion with a positive class stopping gate."""
import argparse
import hashlib
from math import gcd,prod
from pathlib import Path
import runpy
import subprocess
import sys
import time
import retrospective as r
import reference_class_targeted_relations as fixed
import reference_partial_class_extraction as extract
import seeded_reference_class as seed

PROTOCOL=Path(__file__).with_name('REFERENCE_ADAPTIVE_CLASS_COMPLETION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_reference_adaptive_class_completion_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_reference_adaptive_class_completion_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-adaptive-class-completion-v1'

def prepare():
    inp=r.read(fixed.INPUT);previous=r.read(fixed.OUTPUT)
    r.write_new(INPUT,{'schema':'rank-jump.reference-adaptive-class-completion-inputs.v1',
        'form':{k:inp[k] for k in ['cubic_ascending','field_discriminant','S_finite','fixed_a','w_power_basis',
                                  'sl2_matrix','binary_cubic_descending','slope_scale']},
        'previous_fixed_chunks':previous['chunks'],
        'bindings':seed.bindings([Path(__file__),PROTOCOL,fixed.INPUT,fixed.OUTPUT,seed.POOL,
                                extract.ACCEPTED,seed.REFERENCE,fixed.LEGACY,Path(extract.__file__)])})

def worker():
    from sage.all import QQ,ZZ,PolynomialRing,pari,prime_range
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    policy=r.read(PROTOCOL);inp=r.read(INPUT);form=inp['form'];bounds=policy['bounds']
    pari.allocatemem(64000000,bounds['pari_stack_bytes'],silent=True)
    for name,sha in inp['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    legacy=runpy.run_path(str(fixed.LEGACY));state=legacy['prior_state'](6)
    history=legacy['history'](6)
    R=PolynomialRing(QQ,'z');f=R(form['cubic_ascending']);nf=pari.nfinit([pari(f),form['S_finite']])
    assert str(nf.disc())==form['field_discriminant']
    cols=state.base['columns'];lookup={(c['p'],c['hnf']):i for i,c in enumerate(cols)}
    elements=list(r.read(extract.ACCEPTED)['elements_GP'])
    for item in inp['previous_fixed_chunks']:
        path=r.ROOT/item['path'];assert r.digest(path.read_bytes())==item['sha256'];chunk=r.read(path)
        target=chunk['target'];history[target['column']]=max(history.get(target['column'],0),8192)
        for row in chunk['relations']:
            state.add([[lookup[(q['p'],q['hnf'])],q['valuation']] for q in row['ideal_factorization']])
            alpha=pari.Mod(pari(R(row['alpha_ascending'])),pari(f));elements.append(str(pari.nfalgtobasis(nf,alpha)))
    # Include the six equation-only BNF relations in selection as well as extraction.
    pool_keys={tuple(row['alpha_ascending']) for row in r.read(seed.POOL)['relations']}
    for row in r.read(extract.ACCEPTED)['elements_GP']:
        alpha=pari.nfbasistoalg(nf,pari(row));fac=pari.idealfactor(nf,alpha)
        state.add([[lookup[(int(fac[j,0][0]),str(pari.idealhnf(nf,fac[j,0])))],int(fac[j,1])]
                   for j in range(fac.nrows())])
    legacy['selection'].__globals__['history']=lambda wave:dict(history)
    initial=state.report();r.write_new(WORK/'initial-matrix.json',initial)
    w=pari.Mod(pari(R(form['w_power_basis'])),pari(f));a=int(form['fixed_a']);M=form['sl2_matrix']
    c0,c1,c2,c3=map(int,form['binary_cubic_descending']);scale=int(form['slope_scale'])
    primorial=prod(map(int,prime_range(bounds['smooth_bound']+1)))
    completed=[];skips=[];monitors=[];total=0;hits=0;started=time.monotonic()
    def monitor():
        number=len(monitors);path=WORK/('monitor_%03d_elements.json'%number)
        r.write_new(path,{'basis_GP':list(map(str,nf.nf_get_zk())),'elements_GP':list(elements)})
        extract.ACCEPTED=path
        result=extract.compute();result['bindings'].update(seed.bindings([Path(__file__),PROTOCOL,INPUT]))
        out=WORK/('monitor_%03d_result.json'%number);r.write_new(out,result)
        monitors.append({'path':str(out.relative_to(r.ROOT)),'sha256':r.digest(out.read_bytes()),
            'extra':result['additional_independence_witnesses']})
        print('MONITOR',number,'ELEMENTS',result['total_distinct_principal_elements'],'STRICT_RANK',
            result['strict_character_rank'],'EXTRA',result['additional_independence_witnesses'],flush=True)
        return bool(result['additional_independence_witnesses'])
    for level in bounds['strip_levels']:
        selected=legacy['selection'](state,6,level,bounds['targets_per_level_max'])
        r.write_new(WORK/('level_%d_selection.json'%level),selected)
        print('LEVEL',level,'TARGETS',len(selected['targets']),'FORMAL_DIM',state.dimension(),flush=True)
        for ordinal,target in enumerate(selected['targets']):
            if target['column'] in state.inside:
                skips.append({'level':level,'ordinal':ordinal,'column':target['column'],'reason':'now a pivot'});continue
            v1,v2=target['lattice_basis'];pairs=[];values=[];seen=set();begin=time.monotonic()
            for v in range(history.get(target['column'],0)+1,level+1):
                for slope in target['root_slopes_scaled']:
                    center=slope*v//scale
                    for u in range(center-1,center+2):
                        if (u,v) in seen or gcd(u,v)!=1:continue
                        seen.add((u,v));m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                        if gcd(m,n)!=1:continue
                        value=((c0*m+c1*n)*m+c2*n*n)*m+c3*n*n*n
                        assert value and value%target['p']==0
                        pairs.append((u,v,m,n));values.append(value)
            digest=hashlib.sha256();relations=[]
            if values:
                for pair,value,residue in zip(pairs,values,fixed.residues(list(map(abs,values)),primorial)):
                    remainder=fixed.strip(abs(value),gcd(abs(value),residue))
                    digest.update(('%s,%s,%s,%s\n'%(pair[0],pair[1],value,remainder)).encode())
                    if remainder!=1:continue
                    u,v,m,n=pair;alpha=a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
                    assert pari.nfeltnorm(nf,alpha)==a*a*value
                    fac=pari.idealfactor(nf,alpha);factors=[]
                    assert pari.idealhnf(nf,pari.idealfactorback(nf,fac))==pari.idealhnf(nf,alpha)
                    for j in range(fac.nrows()):
                        P,e=fac[j,0],int(fac[j,1]);p=int(P[0]);assert e>0 and p<=bounds['smooth_bound']
                        factors.append([lookup[(p,str(pari.idealhnf(nf,P)))],e])
                    assert any(i==target['column'] for i,e in factors)
                    state.add(factors);elements.append(str(pari.nfalgtobasis(nf,alpha)))
                    relations.append({'u':u,'v':v,'m':m,'n':n,'norm':str(a*a*value),
                        'alpha_ascending':[str(pari.lift(alpha).polcoef(i)) for i in range(3)],'legacy_ideal_factorization':factors})
            path=WORK/('level_%d_target_%04d.json'%(level,ordinal))
            record={'target':target,'level':level,'previous_vmax':history.get(target['column'],0),
                'candidate_count':len(values),'population_digest_sha256':digest.hexdigest(),'relations':relations,
                'ending_formal_dimension':state.dimension(),'elapsed_seconds':time.monotonic()-begin}
            r.write_new(path,record);history[target['column']]=level
            completed.append({'path':str(path.relative_to(r.ROOT)),'sha256':r.digest(path.read_bytes())})
            total+=len(values);hits+=len(relations)
            print('TARGET',len(completed),'DIM',state.dimension(),'RELATIONS',hits,'SECONDS',round(time.monotonic()-started,1),flush=True)
            if len(completed)%bounds['monitor_interval']==0 and monitor():
                r.write_new(WORK/'worker.json',{'status':'CANDIDATE_FOUND','completed':completed,'skips':skips,'monitors':monitors});return
        if monitor():
            r.write_new(WORK/'worker.json',{'status':'CANDIDATE_FOUND','completed':completed,'skips':skips,'monitors':monitors});return
    r.write_new(WORK/'worker.json',{'status':'BOUNDS_COMPLETED','completed':completed,'skips':skips,'monitors':monitors})

def capture():
    WORK.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            reason=None if p.returncode==0 else 'worker failure'
        except subprocess.TimeoutExpired:reason='bounded timeout'
    terminal=r.read(WORK/'worker.json') if (WORK/'worker.json').exists() else {'status':'PARTIAL','reason':reason}
    terminal.update(schema='rank-jump.reference-adaptive-class-completion.v1',elapsed_seconds=time.monotonic()-started,
        bindings=seed.bindings([Path(__file__),PROTOCOL,INPUT,WORK/'worker.log']),
        retained_target_chunks=len(list(WORK.glob('level_*_target_*.json'))),
        retained_monitor_results=len(list(WORK.glob('monitor_*_result.json'))))
    r.write_new(OUTPUT,terminal);print(terminal['status'],terminal['retained_target_chunks'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture','worker']);a=p.parse_args();globals()[a.mode]()
