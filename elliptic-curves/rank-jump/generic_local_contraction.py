#!/usr/bin/env python3
"""Masked generic-only certificate of maximal relative Selmer contraction."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('GENERIC_LOCAL_CONTRACTION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_generic_local_contraction_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_local_contraction_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-generic-local-contraction-v1'
LOCAL_SOURCE=r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,LOCAL_SOURCE,INPUT)}


def export():
    # This is the only mode allowed to read the previous retrospective inputs.
    # The worker below accepts an explicit whitelist and reads only the masked file.
    import production_minus_twist_completion as old
    import remaining_bad_primes as rem
    cases=[]
    for row in old.merged():
        i=row['case_index'];source=rem.bad.cases()[i]
        kept=r.read(rem.bad.INPUT)['cases'][i]
        selected=kept['selected_input_indices'][:kept['generic_dimension']]
        assert all(j<len(source['generic_points']) for j in selected)
        model,points=r.short(source['model'],source['generic_points'])
        scale=r.F(kept['elliptic_scaling_d'])
        points=[[str(r.F(points[j][0])*scale**2),str(r.F(points[j][1])*scale**3)] for j in selected]
        cases.append({'case_index':i,'cubic_ascending':kept['integral_cubic_ascending'],
            'discriminant_factors':r.read(rem.INPUT)['cases'][i]['factor']['factors'],
            'generic_points':points,'local':[{'place':x['place'],
                 'twist_x_witnesses':[w['short_x'] for w in x['witnesses']]} for x in row['local']]})
    r.write_new(INPUT,{'schema':'rank-jump.generic-local-contraction-inputs.v1',
        'mask':'Generic points and local twist witnesses only; retrospective selection, no exceptional point input.',
        'cases':cases})


def coordinates(v,columns):
    piv={}
    for i,b in enumerate(columns):
        mask=1<<i
        while b:
            k=b.bit_length()-1
            if k not in piv:piv[k]=(b,mask);break
            b^=piv[k][0];mask^=piv[k][1]
    mask=0
    for k in sorted(piv,reverse=True):
        if v>>k&1:v^=piv[k][0];mask^=piv[k][1]
    if v:raise ValueError('outside column span')
    return mask


def xor_selected(mask,columns):
    value=0
    for i,x in enumerate(columns):
        if mask>>i&1:value^=x
    return value


def linear_certificate(generic,original,twist):
    tb=r.basis(twist);residual=[r.reduce(x,tb) for x in generic]
    d=r.rank(original+twist)-r.rank(twist);e=r.rank(residual)
    complete=e==d;corrections=[]
    if complete:
        for x in original:
            mask=coordinates(r.reduce(x,tb),residual)
            assert r.reduce(x^xor_selected(mask,generic),tb)==0
            corrections.append(mask)
    return {'local_change_dimension':d,'generic_quotient_dimension':e,
            'generic_surjectivity_certified':complete,'original_basis_generic_correction_masks':corrections}


def compute(index):
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari,GF,VectorSpace
    from sage.version import version
    sys.path.insert(0,str(LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    row=next(x for x in r.read(INPUT)['cases'] if x['case_index']==index)
    assert set(row)=={'case_index','cubic_ascending','discriminant_factors','generic_points','local'}
    assert all(set(x)=={'place','twist_x_witnesses'} for x in row['local'])
    R=PolynomialRing(QQ,'z');f=R(list(map(QQ,row['cubic_ascending'])))
    assert f[3]==1 and f[2]==0
    primes=[x['place'] for x in row['local'] if x['place']!='infinity']
    product=ZZ(1)
    for p,e in row['discriminant_factors']:
        assert p in primes and ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert product==abs(16*f.discriminant())
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f))
    points=[tuple(map(QQ,P)) for P in row['generic_points']]
    assert all(y*y==f(x) for x,y in points)
    betas=[pari(x)-theta for x,y in points];roots=f.roots(AA,multiplicities=False)
    generic=[0]*len(points);original=[];twist=[];offset=0;locals=[]
    for local in row['local']:
        p=local['place'];witnesses=[]
        if p=='infinity':
            sigs=[r.pack(int(x<a) for a in roots) for x,y in points]
            dim=int(len(roots)==3);width=len(roots);tb=[3] if dim else []
            assert not local['twist_x_witnesses']
        else:
            chars=LocalSquareclasses(nf,p);dim=chars.point_kummer_dimension
            sigs=[r.pack(chars.signature(beta)) for beta in betas];width=len(chars.signature(betas[0]))
            if p%4==1:
                assert not local['twist_x_witnesses'];tb=list(r.basis(sigs).values())
            else:
                tb=[]
                for raw in local['twist_x_witnesses']:
                    x=QQ(raw);value=x**3+f[1]*x-f[0];assert value!=0
                    v=int(value.valuation(p));unit=value/QQ(p)**v
                    assert v%2==0 and pari.issquare(pari(unit)+pari(f'O({p}^8)'))==1
                    beta=pari(x)+theta;sig=r.pack(chars.signature(beta));tb.append(sig)
                    witnesses.append({'x':raw,'valuation':v,'signature':sig})
        ob=list(r.basis(sigs).values())
        if len(ob)!=dim or r.rank(tb)!=dim:
            return {'bindings':bindings(),'case_index':index,'status':'UNKNOWN','incomplete_place':p}
        for i,v in enumerate(sigs):generic[i]|=v<<offset
        original.extend(v<<offset for v in ob);twist.extend(v<<offset for v in tb)
        locals.append({'place':p,'width':width,'dimension':dim,'generic_signatures':sigs,
                       'original_basis':ob,'twist_basis':tb,'twist_witnesses':witnesses})
        offset+=width
    result=linear_certificate(generic,original,twist)
    # Independent matrix implementation, also checking the explicit correction map.
    V=VectorSpace(GF(2),offset)
    vec=lambda x:V([(x>>i)&1 for i in range(offset)])
    span=lambda xs:V.subspace([vec(x) for x in xs])
    O,T,G=span(original),span(twist),span(generic);C=O.intersection(T)
    assert G.is_subspace(O)
    assert O.dimension()-C.dimension()==result['local_change_dimension']
    assert G.dimension()-G.intersection(T).dimension()==result['generic_quotient_dimension']
    for x,mask in zip(original,result['original_basis_generic_correction_masks']):
        assert vec(x^xor_selected(mask,generic)) in C
    return {'bindings':bindings(),'case_index':index,'status':'PASS' if result['generic_surjectivity_certified'] else 'GATE_NOT_MET',
        'software':{'sage':version,'pari':str(pari.version())},'local':locals,
        'generic_point_count':len(points),'joint_generic_dimension':int(G.dimension()),
        'local_product_dimension':int(O.dimension()),'joint_generic_signatures':generic,
        'original_product_basis':original,'twist_product_basis':twist,**result,
        'boundary':'Proves an exact relative Selmer contraction only when the gate passes. No exceptional points, derivative class, global class group or CT values were read.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in r.read(PROTOCOL)['cases']:
        if check:
            assert compute(i)==next(x for x in r.read(OUTPUT)['rows'] if x['case_index']==i)
            print('PASS independent masked replay',i,flush=True);continue
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--index',str(i),
                        '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:failure='30-second timeout'
                if failure and not path.exists():r.write_new(path,{'bindings':bindings(),'case_index':i,'status':'UNKNOWN','reason':failure})
        record=r.read(path);assert record['bindings']==bindings();rows.append(record)
        print('checkpoint',i,record['status'],record.get('generic_quotient_dimension'),flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.generic-local-contraction.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check'])
    p.add_argument('--index',type=int);p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index))
    else:capture(args.mode=='check')
