#!/usr/bin/env python3
"""Small exact local repairs and factorization-free zero-kernel certificates."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base

PROTOCOL=Path(__file__).with_name('FRESH_GOVERNING_COMPLETION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_governing_completion_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-governing-completion-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,base.OUTPUT,base.INPUT,Path(base.__file__),Path(r.__file__),base.LOCAL)}


def subset_local(token):
    from sage.all import QQ,ZZ,AA,pari,EllipticCurve,GF,matrix
    import local_collision as lc
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    f,pts,_=base.model_data(token);E=EllipticCurve([0,0,0,f[1],f[0]])
    primes=[2]+[p for p in r.primes(r.read(PROTOCOL)['limits']['trial_prime_bound']) if f.discriminant()%p==0 and E.local_data(p,proof=True).conductor_valuation()>0]
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f))
    betas=[pari(x)-theta for x,y in pts];columns=[];locals=[]
    for p in primes:
        chars=LocalSquareclasses(nf,p);sigs=[list(chars.signature(beta)) for beta in betas]
        for j in range(len(sigs[0])):columns.append(r.pack(s[j] for s in sigs))
        locals.append({'place':p,'signatures':sigs,'point_kummer_dimension':chars.point_kummer_dimension})
    roots=f.roots(AA,multiplicities=False);signs=[[int(x<a) for a in roots] for x,y in pts]
    for j in range(len(roots)):columns.append(r.pack(s[j] for s in signs))
    locals.append({'place':'infinity','signatures':signs,'point_kummer_dimension':int(len(roots)==3)})
    masks=lc.orthogonal(columns,len(pts))
    assert matrix(GF(2),[[v>>i&1 for i in range(len(pts))] for v in columns]).right_kernel().dimension()==len(masks)
    # Generic mod2 independence cannot be inferred solely from a nonzero kernel.
    model=['0','0','0',str(f[1]),str(f[0])];blocks=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['trial_prime_bound']):
        roots_p=r.roots_at(model[3],model[4],p)
        if roots_p:blocks.append((p,roots_p))
        signatures=[r.point_signature(model,list(map(str,P)),blocks) for P in pts]
        if r.rank(signatures)==len(pts):break
    assert r.rank(signatures)==len(pts)
    return {'status':'PASS_ZERO_KERNEL' if not masks else 'PARTIAL_KERNEL_UPPER_BOUND',
      'tested_bad_primes':primes,'local':locals,'generic_dimension':len(pts),'tested_kernel_masks':masks,
      'strict_generic_dimension':0 if not masks else None,'strict_generic_dimension_upper_bound':len(masks),
      'minus_twist_CT_rank':0 if not masks else None,'independence_blocks':blocks,'independence_signatures':signatures,
      'galois':r.galois(model),'boundary':'Only zero tested kernel certifies the full strict kernel without the remaining bad primes.'}


def repair(token,old):
    from sage.all import QQ,ZZ,pari,matrix
    f,pts,_=base.model_data(token);prior=next(x for x in r.read(base.OUTPUT)['rows'] if x['token']==token)
    primes=[p for p,e in prior['factor']['factors']];nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f))
    betas=[pari.Mod(pari(f.parent()(list(map(QQ,x['beta_ascending'])))),pari(f)) for x in old['class_records']]
    A=[list(row) for row in old['Artin_matrix']];records=[]
    for j,c in enumerate(old['artin_columns']):
        if not c['cyclic']:return {'status':'UNKNOWN','reason':'noncyclic residual ring'}
        H=pari(matrix(QQ,c['coprime_ideal_hnf']));N=ZZ(c['norm'])
        for i,entry in enumerate(c['evaluations']):
            if A[i][j] is not None:continue
            gcd=ZZ(entry['nonunit_gcd'])
            if gcd.nbits()>r.read(PROTOCOL)['limits']['nonunit_gcd_maximum_bits']:
                return {'status':'UNKNOWN','reason':'nonunit gcd exceeds frozen small-factor cap'}
            support=[p for p,e in gcd.factor(proof=True)];cofactor=N;bit=0;local=[]
            for p in support:
                e=int(N.valuation(p));cofactor//=p**e
                contributions=[]
                for P in pari.idealprimedec(nf,p):
                    exponent=int(pari.idealval(nf,H,P))
                    if not exponent:continue
                    assert int(P[3])==1 # cyclic rational residue ring
                    valuation=int(pari.idealval(nf,betas[i],P));assert valuation%2==0
                    frob=int(pari.nfislocalpower(nf,P,betas[i],2)==0);bit^=(exponent%2)*frob
                    # Independent rational-point character product at this labelled root.
                    assert f.discriminant()%p
                    tc=pari.nfalgtobasis(nf,theta);root=int((tc[0]-H[0,1]*tc[1]-H[0,2]*tc[2])%p)
                    assert pari.idealval(nf,theta-root,P)>0
                    mask=old['class_records'][i]['generic_mask'];point_bit=0
                    for k,point in enumerate(pts):
                        if mask>>k&1:point_bit^=r.point_signature(['0','0','0',str(f[1]),str(f[0])],list(map(str,point)),[(int(p),[root])])
                    assert point_bit==frob
                    contributions.append({'exponent':exponent,'valuation':valuation,'frobenius_bit':frob,'independent_generic_character_bit':point_bit,'root':root})
                local.append({'prime':int(p),'norm_exponent':e,'contributions':contributions})
            coords=pari.nfalgtobasis(nf,betas[i]);residue=ZZ(coords[0]-H[0,1]*coords[1]-H[0,2]*coords[2])%cofactor
            assert residue.gcd(cofactor)==1
            jacobi=int(pari.kronecker(residue,cofactor));assert jacobi in (-1,1)
            bit^=int(jacobi==-1);A[i][j]=bit
            records.append({'row':i,'column':j,'local':local,'cofactor':str(cofactor),'residue':str(residue),'jacobi':jacobi,'artin_bit':bit})
    assert all(x is not None for row in A for x in row)
    k=len(A);M=[[A[i][j]^A[j][i] for j in range(k)] for i in range(k)]
    arank=r.rank(list(map(r.pack,A)));mrank=r.rank(list(map(r.pack,M)))
    assert int(matrix(__import__('sage.all',fromlist=['GF']).GF(2),M).rank())==mrank
    return {'status':'PASS_REPAIRED','repairs':records,'Artin_matrix':A,'minus_twist_CT_matrix':M,
             'strict_generic_dimension':k,'Artin_rank':arank,'minus_twist_CT_rank':mrank}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for old in r.read(base.OUTPUT)['rows']:
        token=old['token'];path=WORK/f'{token}.json'
        if old['local']['status']=='PASS':continue
        if not path.exists():
            with (WORK/f'{token}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['worker_seconds'])
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='45-second timeout'
            if reason:r.write_new(path,{'bindings':bindings(),'token':token,'status':'UNKNOWN','reason':reason})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(token,row['status'],row.get('strict_generic_dimension'),row.get('minus_twist_CT_rank'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-governing-completion.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        old=next(x['local'] for x in r.read(base.OUTPUT)['rows'] if x['token']==args.token)
        result=repair(args.token,old) if old['status']=='PARTIAL' else subset_local(args.token)
        r.write_new(WORK/f'{args.token}.json',{'bindings':bindings(),'token':args.token,**result})
