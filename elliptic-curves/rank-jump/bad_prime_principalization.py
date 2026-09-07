#!/usr/bin/env python3
"""Reduced-ideal principal-square probe on the frozen 103b2 pair."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base
from prime_square_class_constructor import directions

PROTOCOL=Path(__file__).with_name('BAD_PRIME_PRINCIPALIZATION_PROTOCOL.json')
PRIOR=r.OUT/'rank_jump_matched103b2_class_boundary_v1.json'
INPUT=r.OUT/'rank_jump_bad_prime_principalization_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_bad_prime_principalization_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-bad-prime-principalization-v1'


def bind(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    from sage.all import QQ,PolynomialRing
    cases=[];R=PolynomialRing(QQ,'z');z=R.gen()
    for old in r.read(PRIOR)['rows']:
        token=old['token'];red=old['reduction'];f=R(red['reduced_cubic_ascending'])
        image=R(red['original_root_in_reduced']);original,pts,scale=base.model_data(token)
        assert R(original.list())(image)%f==0
        controls=[list(map(str,((QQ(x)-image)%f).list())) for x,y in pts[:3]]
        cases.append({'token':token,'cubic':list(map(str,f.list())),
            'basis':red['transported_maximal_order_basis'],'field_discriminant':red['field_discriminant'],
            'S_finite':old['boundary']['S_finite'],'generic_controls':controls})
    assert [c['token'] for c in cases]==r.read(PROTOCOL)['cases']
    r.write_new(INPUT,{'schema':'rank-jump.bad-prime-principalization-inputs.v1','cases':cases,
        'bindings':bind([Path(__file__),PROTOCOL,PRIOR,base.INPUT,Path(base.__file__)])})


def worker(token):
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    pari.allocatemem(64000000,268435456,silent=True)
    data=next(c for c in r.read(INPUT)['cases'] if c['token']==token)
    R=PolynomialRing(QQ,'z');f=R(data['cubic']);theta=pari.Mod('z',pari(f))
    nf=pari.nfinit([pari(f),[pari(R(b)) for b in data['basis']]])
    assert str(nf.disc())==data['field_discriminant']
    zk=list(nf.nf_get_zk());T=matrix(ZZ,3,3,lambda i,j:ZZ(pari.nfelttrace(nf,zk[i]*zk[j])))
    assert T.is_positive_definite();ds=directions();assert len(ds)==49
    enc=lambda a:[str(pari.lift(a).polcoef(i)) for i in range(3)]
    mat=lambda A:[[str(A[i,j]) for j in range(3)] for i in range(3)]
    targets=[]
    for p in data['S_finite']:
        for j,P in enumerate(pari.idealprimedec(nf,p)):
            I=pari.idealpow(nf,P,2)
            targets.append((f'{p}:{j}','candidate',I,{'p':p,'prime_index':j,'e':int(P[2]),'f':int(P[3]),'prime_hnf':mat(pari.idealhnf(nf,P))}))
    for j,cs in enumerate(data['generic_controls']):
        alpha=pari(R(cs))(theta);I=pari.idealhnf(nf,alpha)
        targets.append((f'generic:{j}','control',I,{'generic_index':j}))
    rows=[]
    for label,kind,I,meta in targets:
        J,a=pari.idealred(nf,[I,1]);multiplier=pari.nfbasistoalg(nf,a)
        assert pari.idealmul(nf,J,multiplier)==I
        H=matrix(QQ,3,3,lambda i,j:QQ(J[i,j]));gram=H.transpose()*T*H
        U=matrix(ZZ,pari.qflllgram(pari(gram)));assert abs(U.det())==1
        N=QQ(pari.idealnorm(nf,I));JN=QQ(pari.idealnorm(nf,J));assert N.is_square()
        norms=[];hits=[]
        for v in ds:
            coords=H*U*vector(ZZ,v)
            a=pari.nfbasistoalg(nf,pari.Col(list(coords)));norm=QQ(pari.nfeltnorm(nf,a));norms.append(str(norm))
            if abs(norm)!=JN:continue
            alpha=a*multiplier;sign=1 if pari.nfeltnorm(nf,alpha)>0 else -1;alpha*=sign
            assert pari.idealhnf(nf,alpha)==I and QQ(pari.nfeltnorm(nf,alpha))==N
            hits.append({'direction':list(v),'sign':sign,'alpha':enc(alpha),'norm':str(N)})
        rows.append({'id':label,'kind':kind,**meta,'target_hnf':mat(I),'target_norm':str(N),
            'reduced_hnf':mat(J),'reduced_norm':str(JN),'multiplier':enc(multiplier),
            'lll_matrix':mat(U),'tested_norms':norms,'hits':hits})
        print(token,label,len(hits),flush=True)
    return {'token':token,'status':'PASS','basis':[enc(b) for b in zk],'trace_gram':mat(T),
        'directions':list(map(list,ds)),'rows':rows,
        'candidate_generators':sum(len(z['hits']) for z in rows if z['kind']=='candidate'),
        'control_generators':sum(len(z['hits']) for z in rows if z['kind']=='control'),
        'strict_class_status':'NOT_YET_TESTED: exact generator output is not a class-creation assertion'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(PROTOCOL)['cases']:
        out=WORK/f'{token}.json';log=WORK/f'{token}.log'
        if not out.exists():
            with log.open('x') as stream:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=stream,stderr=stream,timeout=60)
                    error=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:error='60-second timeout'
            if error:r.write_new(out,{'token':token,'status':'UNKNOWN','reason':error})
        row=r.read(out);rows.append(row);print(token,row['status'],row.get('candidate_generators'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.bad-prime-principalization.v1','cases':rows,
        'bindings':bind([Path(__file__),PROTOCOL,INPUT,Path(__file__).with_name('prime_square_class_constructor.py')]),
        'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['export','capture','worker']);parser.add_argument('--token');args=parser.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'{args.token}.json',worker(args.token))
    else:globals()[args.mode]()
