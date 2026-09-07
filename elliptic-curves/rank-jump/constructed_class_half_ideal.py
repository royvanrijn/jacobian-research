#!/usr/bin/env python3
"""Ordinary half-ideal images of the two equation-only constructed strict classes."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('CONSTRUCTED_CLASS_HALF_IDEAL_PROTOCOL.json')
REFERENCE=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
CONSTRUCTION=r.OUT/'rank_jump_reference_additional_strict_verification_v1.json'
EXTRACTION=r.OUT/'rank_jump_reference_large_support_extraction_v1.json'
ANCHORS=r.OUT/'rank_jump_generic_only_class_anchors_v1.json'
OUTPUT=r.OUT/'rank_jump_constructed_class_half_ideal_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-constructed-class-half-ideal-v1'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    ref=r.read(REFERENCE);c=r.read(CONSTRUCTION);ext=r.read(EXTRACTION);anc=r.read(ANCHORS)['cases'][0]
    assert c['status']=='PASS' and c['additional_independent_strict_classes']==2
    for source in [c,ext,r.read(ANCHORS)]:
        for path,sha in source['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,r.read(PROTOCOL)['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    gammas=[pari.Mod(pari(R(x['beta_ascending'])),pari(f)) for x in ref['generic_classes']]
    bad=[(p,str(pari.idealhnf(nf,P)),P) for p in ref['S_finite'] for P in pari.idealprimedec(nf,p)]
    generic_gcd=[]
    for g,row in zip(gammas,ref['generic_classes']):
        b=ZZ(row['norm']).sqrt();assert b*b==pari.nfeltnorm(nf,g)
        generic_gcd.append(pari.idealadd(nf,b,g))
    columns=c['columns'];atoms={a['index']:a for a in c['atoms']};factors={('generic',i):a for i,a in enumerate(gammas)}
    atom_valuations={}
    for i,a in atoms.items():
        factors[('projected_atom',i)]=ZZ(a['norm'])*pari.Mod(pari(R(a['alpha_ascending'])),pari(f))
        values=dict(a['ideal_factorization']);normvals={}
        for j,e in values.items():
            col=columns[j];p=col['p'];normvals[p]=normvals.get(p,0)+col['f']*e
        atom_valuations[i]={j:values.get(j,0)+col['e']*normvals[col['p']] for j,col in enumerate(columns) if col['p'] in normvals}
    coords={key:list(map(ZZ,pari.nfalgtobasis(nf,value))) for key,value in factors.items()}
    def generic_labels(mask):return [{'kind':'generic','index':i} for i in range(16) if mask>>i&1]
    class_labels=[generic_labels(int(m)) for m in ext['generic_strict_masks']]+[z['factor_labels'] for z in c['class_representatives']]
    character_labels=[generic_labels(int(m)) for m in anc['generic_coefficient_masks']]+[z['factor_labels'] for z in c['class_representatives']]
    assert len(class_labels)==8 and len(character_labels)==12
    ideals=[];cache={};S=set(ref['S_finite'])
    directions=[[0,0,0]]+[[0,a,b] for a in range(-4,5) for b in range(-4,5) if a or b]
    for number,labels in enumerate(class_labels):
        generic_ids=[x['index'] for x in labels if x['kind']=='generic'];exponents={}
        for x in labels:
            if x['kind']=='projected_atom':
                for j,e in atom_valuations[x['index']].items():exponents[j]=exponents.get(j,0)+e
        ideal=pari.idealhnf(nf,1);steps=[]
        def multiply(I,description):
            nonlocal ideal
            product=pari.idealmul(nf,ideal,I);J,a=pari.idealred(nf,[product,1])
            assert pari.idealmul(nf,J,a)==product
            steps.append({'input':description,'reduced_hnf':str(J),'principal_multiplier_GP':str(a)})
            ideal=J
        for i in generic_ids:multiply(generic_gcd[i],{'generic_gcd':i})
        for j,e in sorted(exponents.items()):
            col=columns[j]
            if col['p'] in S:continue
            assert e%2==0
            key=(col['p'],col['hnf'])
            if key not in cache:
                P=next(P for P in pari.idealprimedec(nf,col['p']) if str(pari.idealhnf(nf,P))==col['hnf'])
                assert int(P[2])==col['e'] and int(P[3])==col['f'];cache[key]=P
            multiply(pari.idealpow(nf,cache[key],e//2),{'column':j,'exponent':e//2})
        corrections=[]
        for p,h,P in bad:
            atomval=sum(e for j,e in exponents.items() if columns[j]['p']==p and columns[j]['hnf']==h)
            genval=sum(int(pari.idealval(nf,gammas[i],P)) for i in generic_ids)
            gcdval=sum(int(pari.idealval(nf,generic_gcd[i],P)) for i in generic_ids)
            assert (atomval+genval)%2==0;correction=(atomval+genval)//2-gcdval
            corrections.append({'p':p,'hnf':h,'atom_valuation':atomval,'generic_valuation':genval,'generic_gcd_valuation':gcdval,'exponent':correction})
            if correction:multiply(pari.idealpow(nf,P,correction),{'bad_prime':p,'hnf':h,'exponent':correction})
        trials=[];success=None
        for direction in directions:
            H,a=pari.idealred(nf,[ideal,1],direction);assert pari.idealmul(nf,H,a)==ideal
            N=ZZ(pari.idealnorm(nf,H));cyclic=bool(N>0 and N%2 and H[0,0]==N and H[1,1]==H[2,2]==1)
            trial={'direction':direction,'hnf':str(H),'norm':str(N),'cyclic_odd':cyclic}
            if cyclic:
                residues={key:ZZ(v[0]-H[0,1]*v[1]-H[0,2]*v[2])%N for key,v in coords.items()}
                nonunit=next((key for key,value in residues.items() if value.gcd(N)!=1),None)
                if nonunit is None:
                    symbols={key:int(pari.kronecker(value,N)) for key,value in residues.items()}
                    assert all(v in (-1,1) for v in symbols.values());bits=[]
                    for character in character_labels:bits.append(sum(int(symbols[(x['kind'],x['index'])]==-1) for x in character)%2)
                    trial.update(artin_bits=bits,principal_multiplier_GP=str(a));success=trial;trials.append(trial);break
                trial['nonunit_factor']=list(nonunit)
            trials.append(trial)
        record={'column':number,'role':'generic_strict' if number<6 else 'constructed_extra',
            'factor_labels':labels,'half_ideal_reduction_steps':steps,'bad_corrections':corrections,
            'final_reduced_half_ideal_hnf':str(ideal),'artin_trials':trials,'status':'PASS' if success else 'UNKNOWN',
            'artin_bits':success['artin_bits'] if success else None}
        ideals.append(record);r.write_new(WORK/('column_%02d.json'%number),record)
        print('COLUMN',number,record['status'],record['artin_bits'],flush=True)
    full=all(x['status']=='PASS' for x in ideals)
    grank=r.rank([r.pack(x['artin_bits']) for x in ideals[:6]]) if all(x['status']=='PASS' for x in ideals[:6]) else None
    rank=r.rank([r.pack(x['artin_bits']) for x in ideals]) if full else None
    return {'schema':'rank-jump.constructed-class-half-ideal.v1','status':'PASS' if full else 'PARTIAL',
        'character_labels':character_labels,'columns':ideals,'generic_half_ideal_artin_rank':grank,
        'combined_half_ideal_artin_rank':rank,'detected_relative_half_ideal_rank':rank-grank if full else None,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,REFERENCE,CONSTRUCTION,EXTRACTION,ANCHORS]},
        'boundary':'Ordinary ideal-class two-torsion images tested by certified unramified characters. A positive relative rank proves new half-ideal information; zero does not prove a unit representation. This is not CT and does not decide rational versus Sha.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()
