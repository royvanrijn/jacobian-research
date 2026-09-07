#!/usr/bin/env python3
"""Construct virtual units using small representatives of unresolved ideal classes."""
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
import reference_adaptive_class_completion as adaptive
import reference_partial_class_extraction as extract
import seeded_reference_class as seed

PROTOCOL=Path(__file__).with_name('REFERENCE_SMALL_REPRESENTATIVE_CLASS_PROTOCOL.json')
ANCHORS=r.OUT/'rank_jump_generic_only_class_anchors_v1.json'
INPUT=r.OUT/'rank_jump_reference_small_representative_class_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_reference_small_representative_class_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-small-representative-class-v1'

def indices(mask):
    while mask:
        bit=mask&-mask;yield bit.bit_length()-1;mask^=bit

class Matrix:
    def __init__(self,base,anchors,bound):
        self.base=base;self.anchors=sum(1<<i for i in anchors)
        self.canonical={max(row['parity_columns']):sum(1<<i for i in row['parity_columns']) for row in base['canonical_rational_relations']}
        assert not set(anchors)&self.canonical.keys()
        self.inside_mask=sum(1<<i for i,c in enumerate(base['columns']) if c['p']<=bound)
        self.free=[i for i,c in enumerate(base['columns']) if c['p']<=bound and i not in self.canonical]
        self.rows=[];self.inside={};self.outside={}

    def add(self,factors):
        row=sum(1<<i for i,e in factors if e%2)
        for i,e in factors:
            if e%2 and i in self.canonical:row^=self.canonical[i]
        self.rows.append(row);inside=row&self.inside_mask;outside=row&~self.inside_mask
        while outside:
            p=outside.bit_length()-1
            if p not in self.outside:self.outside[p]=(outside,inside);return
            o,j=self.outside[p];outside^=o;inside^=j
        while inside:
            free=inside&~self.anchors
            assert free,'principal relation contradicts independently certified anchors'
            p=free.bit_length()-1
            if p not in self.inside:self.inside[p]=inside;return
            inside^=self.inside[p]

    def projection(self,column):
        assert column not in self.canonical
        row=(1<<column)&self.inside_mask;outside=(1<<column)&~self.inside_mask
        while outside:
            p=outside.bit_length()-1
            if p not in self.outside:return 0
            o,j=self.outside[p];outside^=o;row^=j
        row&=~self.anchors;result=0
        while row:
            p=row.bit_length()-1
            if p in self.inside:row^=self.inside[p]
            else:result^=1<<p;row^=1<<p
            row&=~self.anchors
        return result

    def report(self):
        total=r.rank(self.rows);outside=r.rank([x&~self.inside_mask for x in self.rows])
        assert total-outside==len(self.inside) and outside==len(self.outside)
        return {'formal_quotient_dimension':len(self.free)-len(self.inside),'protected_generic_anchor_dimension':self.anchors.bit_count(),
                'supported_rank':len(self.inside),'outside_rank':outside,'all_rank':total}

def prepare():
    terminal=r.read(adaptive.OUTPUT)
    assert terminal['status'] in ['PARTIAL','BOUNDS_COMPLETED'],'do not continue a successful or live experiment'
    chunks=[{'path':str(p.relative_to(r.ROOT)),'sha256':r.digest(p.read_bytes())} for p in sorted(adaptive.WORK.glob('level_*_target_*.json'))]
    r.write_new(INPUT,{'schema':'rank-jump.reference-small-representative-class-inputs.v1',
        'form':r.read(adaptive.INPUT)['form'],'previous_fixed_chunks':r.read(fixed.OUTPUT)['chunks'],
        'previous_adaptive_chunks':chunks,'bindings':seed.bindings([Path(__file__),PROTOCOL,ANCHORS,adaptive.OUTPUT,
            adaptive.INPUT,fixed.OUTPUT,seed.POOL,seed.REFERENCE,extract.ACCEPTED,Path(extract.__file__),fixed.LEGACY])})

def worker():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,pari,prime_range
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    inp=r.read(INPUT);bounds=r.read(PROTOCOL)['bounds'];form=inp['form']
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,bounds['pari_stack_bytes'],silent=True)
    legacy=runpy.run_path(str(fixed.LEGACY));old=legacy['prior_state'](6);history=legacy['history'](6)
    R=PolynomialRing(QQ,'z');f=R(form['cubic_ascending']);nf=pari.nfinit([pari(f),form['S_finite']])
    cols=old.base['columns'];lookup={(c['p'],c['hnf']):i for i,c in enumerate(cols)}
    reference=r.read(seed.REFERENCE);controls=r.read(ANCHORS)['cases'][0]
    betas=[R(c['beta_ascending']) for c in reference['generic_classes']]
    masks=list(map(int,controls['generic_coefficient_masks']));anchor_cols=[];anchor_labels=[];anchor_records=[]
    for i,c in enumerate(cols):
        p=c['p']
        if p<101 or p>bounds['anchor_prime_bound'] or i in old.canonical or c['f']!=1 or c['e']!=1 or f.discriminant()%p==0:continue
        rr=[int(x) for x in f.change_ring(GF(p)).roots(multiplicities=False)
            if str(pari.idealhnf(nf,p,pari.Mod(pari(R.gen()-int(x)),pari(f))))==c['hnf']]
        assert len(rr)==1;root=rr[0];vals=[int(b.change_ring(GF(p))(root)) for b in betas]
        if not all(vals):continue
        raw=r.pack([int(pow(v,(p-1)//2,p)==p-1) for v in vals]);label=r.pack([(raw&m).bit_count()%2 for m in masks])
        if r.rank(anchor_labels+[label])==len(anchor_labels):continue
        anchor_cols.append(i);anchor_labels.append(label);anchor_records.append({'column':i,'p':p,'root':root,'hnf':c['hnf'],'label':str(label)})
        if len(anchor_cols)==len(masks):break
    assert len(anchor_cols)==len(masks)
    state=Matrix(old.base,anchor_cols,bounds['inside_bound'])
    for row in old.rows:state.add([[i,1] for i in indices(row)])
    elements=list(r.read(extract.ACCEPTED)['elements_GP'])
    for key in ['previous_fixed_chunks','previous_adaptive_chunks']:
        for item in inp[key]:
            path=r.ROOT/item['path'];assert r.digest(path.read_bytes())==item['sha256'];chunk=r.read(path)
            col=chunk['target']['column'];history[col]=max(history.get(col,0),chunk.get('level',8192))
            for row in chunk['relations']:
                fac=row.get('legacy_ideal_factorization')
                if fac is None:fac=[[lookup[(q['p'],q['hnf'])],q['valuation']] for q in row['ideal_factorization']]
                state.add(fac);alpha=pari.Mod(pari(R(row['alpha_ascending'])),pari(f));elements.append(str(pari.nfalgtobasis(nf,alpha)))
    for encoded in r.read(extract.ACCEPTED)['elements_GP']:
        fac=pari.idealfactor(nf,pari.nfbasistoalg(nf,pari(encoded)))
        state.add([[lookup[(int(fac[j,0][0]),str(pari.idealhnf(nf,fac[j,0])))],int(fac[j,1])] for j in range(fac.nrows())])
    for j in range(16):
        state.add([[lookup[(v['p'],v['hnf'])],1] for v in controls['bad_prime_valuations'] if v['valuations'][j]%2])
    r.write_new(WORK/'initial-matrix.json',{'matrix':state.report(),'generic_anchors':anchor_records})
    w=pari.Mod(pari(R(form['w_power_basis'])),pari(f));a=int(form['fixed_a']);M=form['sl2_matrix']
    c0,c1,c2,c3=map(int,form['binary_cubic_descending']);scale=int(form['slope_scale'])
    normroots=R([c3,c2,c1,c0]).roots(AA,multiplicities=False);primorial=prod(map(int,prime_range(bounds['smooth_bound']+1)))
    olddata=legacy['old']['setup']()[0];completed=[];monitors=[];started=time.monotonic()
    def monitor():
        k=len(monitors);path=WORK/('monitor_%03d_elements.json'%k)
        r.write_new(path,{'basis_GP':list(map(str,nf.nf_get_zk())),'elements_GP':list(elements)})
        extract.ACCEPTED=path;result=extract.compute();result['bindings'].update(seed.bindings([Path(__file__),PROTOCOL,INPUT]))
        out=WORK/('monitor_%03d_result.json'%k);r.write_new(out,result)
        monitors.append({'path':str(out.relative_to(r.ROOT)),'sha256':r.digest(out.read_bytes()),'extra':result['additional_independence_witnesses']})
        print('MONITOR',k,'ELEMENTS',result['total_distinct_principal_elements'],'STRICT_RANK',result['strict_character_rank'],'EXTRA',result['additional_independence_witnesses'],flush=True)
        return bool(result['additional_independence_witnesses'])
    def finish(status):r.write_new(WORK/'worker.json',{'status':status,'completed':completed,'monitors':monitors,'matrix':state.report()})
    if monitor():finish('CANDIDATE_FOUND');return
    for level in bounds['strip_levels']:
        eligible=[]
        for i,c in enumerate(cols):
            p=c['p']
            if i in state.canonical or i in anchor_cols or p>bounds['smooth_bound'] or c['f']!=1 or c['e']!=1 or history.get(i,0)>=level:continue
            if f.discriminant()%p==0 or a%p==0:continue
            projection=state.projection(i)
            if projection:eligible.append((i,projection))
        spans={};selected=[]
        for i,v in eligible:
            reduced=r.reduce(v,spans)
            if reduced:spans[reduced.bit_length()-1]=reduced;selected.append(i)
            if len(selected)>=bounds['targets_per_level_max']:break
        selected_set=set(selected)
        selected += [i for i,v in eligible if i not in selected_set][:bounds['targets_per_level_max']-len(selected)]
        selected.sort();r.write_new(WORK/('level_%d_selection.json'%level),{'columns':selected,'eligible':len(eligible),'span_rank':len(spans),'matrix':state.report()})
        print('LEVEL',level,'TARGETS',len(selected),'DIM',state.report()['formal_quotient_dimension'],flush=True)
        for ordinal,column in enumerate(selected):
            if not state.projection(column):continue
            c=cols[column];p=c['p'];root=next(int(x) for x in f.change_ring(GF(p)).roots(multiplicities=False)
                if str(pari.idealhnf(nf,p,pari.Mod(pari(R.gen()-int(x)),pari(f))))==c['hnf'])
            wbar=int(R(form['w_power_basis']).change_ring(GF(p))(root));bm=(a*M[0][0]+wbar*M[1][0])%p;bn=(a*M[0][1]+wbar*M[1][1])%p
            slope=-bn*pow(bm,-1,p)%p;lattice=legacy['old']['reduce_lattice'](p,slope,olddata['reduced_hessian']);v1,v2=lattice
            slopes=[int((((x*v2[1]-v2[0])/(v1[0]-x*v1[1]))*scale).floor()) for x in normroots]
            target={'column':column,'p':p,'hnf':c['hnf'],'root':root,'lattice_basis':lattice,'root_slopes_scaled':slopes,'projection_hex':hex(state.projection(column))}
            pairs=[];values=[];seen=set()
            for v in range(history.get(column,0)+1,level+1):
                for s in slopes:
                    center=s*v//scale
                    for u in range(center-1,center+2):
                        if (u,v) in seen or gcd(u,v)!=1:continue
                        seen.add((u,v));m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                        if gcd(m,n)!=1:continue
                        value=((c0*m+c1*n)*m+c2*n*n)*m+c3*n*n*n
                        assert value and value%p==0;pairs.append((u,v,m,n));values.append(value)
            digest=hashlib.sha256();relations=[]
            for pair,value,residue in zip(pairs,values,fixed.residues(list(map(abs,values)),primorial) if values else []):
                remainder=fixed.strip(abs(value),gcd(abs(value),residue));digest.update(('%s,%s,%s,%s\n'%(pair[0],pair[1],value,remainder)).encode())
                if remainder!=1:continue
                u,v,m,n=pair;alpha=a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
                assert pari.nfeltnorm(nf,alpha)==a*a*value;fac=pari.idealfactor(nf,alpha)
                assert pari.idealhnf(nf,pari.idealfactorback(nf,fac))==pari.idealhnf(nf,alpha)
                factors=[]
                for j in range(fac.nrows()):
                    P,e=fac[j,0],int(fac[j,1]);q=int(P[0]);assert e>0 and q<=bounds['smooth_bound']
                    factors.append([lookup[(q,str(pari.idealhnf(nf,P)))],e])
                assert any(i==column for i,e in factors);state.add(factors);elements.append(str(pari.nfalgtobasis(nf,alpha)))
                relations.append({'m':m,'n':n,'u':u,'v':v,'norm':str(a*a*value),'alpha_ascending':[str(pari.lift(alpha).polcoef(i)) for i in range(3)],'legacy_ideal_factorization':factors})
            path=WORK/('level_%d_target_%04d.json'%(level,ordinal))
            r.write_new(path,{'target':target,'level':level,'previous_vmax':history.get(column,0),'candidate_count':len(values),
                'population_digest_sha256':digest.hexdigest(),'relations':relations,'ending_formal_dimension':len(state.free)-len(state.inside)})
            history[column]=level;completed.append({'path':str(path.relative_to(r.ROOT)),'sha256':r.digest(path.read_bytes())})
            print('TARGET',len(completed),'P',p,'DIM',len(state.free)-len(state.inside),'RELATIONS',len(relations),'SECONDS',round(time.monotonic()-started,1),flush=True)
            if len(completed)%bounds['monitor_interval']==0 and monitor():finish('CANDIDATE_FOUND');return
        if monitor():finish('CANDIDATE_FOUND');return
    finish('BOUNDS_COMPLETED')

def capture():
    WORK.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            reason=None if p.returncode==0 else 'worker failure'
        except subprocess.TimeoutExpired:reason='bounded timeout'
    terminal=r.read(WORK/'worker.json') if (WORK/'worker.json').exists() else {'status':'PARTIAL','reason':reason}
    terminal.update(schema='rank-jump.reference-small-representative-class.v1',elapsed_seconds=time.monotonic()-started,
        bindings=seed.bindings([Path(__file__),PROTOCOL,INPUT,WORK/'worker.log']),
        retained_target_chunks=len(list(WORK.glob('level_*_target_*.json'))),retained_monitor_results=len(list(WORK.glob('monitor_*_result.json'))))
    r.write_new(OUTPUT,terminal);print(terminal['status'],terminal['retained_target_chunks'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture','worker']);a=p.parse_args();globals()[a.mode]()
