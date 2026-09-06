#!/usr/bin/env python3
"""Independent equation-defined norm projections on the completed masked panel."""
import argparse
from itertools import product
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base
import matched103b2_class_boundary as original
import fresh_retained_factors as supplement

PROTOCOL=Path(__file__).with_name('FRESH_NORM_PROJECTION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_norm_projection_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-norm-projection-v1'


def cases():
    a=[x['token'] for x in r.read(base.OUTPUT)['rows'] if x['factor']['status']=='PASS']
    b=[x['token'] for x in r.read(supplement.OUTPUT)['rows'] if x['boundary']['status']=='PASS']
    return sorted(set(a+b))


def bindings():
    paths=(Path(__file__),PROTOCOL,base.INPUT,base.OUTPUT,supplement.OUTPUT,Path(base.__file__),
           Path(original.__file__),Path(supplement.__file__),base.LOCAL,Path(r.__file__))
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def setup(token):
    from sage.all import ZZ,pari
    row=next(x for x in r.read(base.OUTPUT)['rows'] if x['token']==token)
    if row['factor']['status']!='PASS':row=next(x for x in r.read(supplement.OUTPUT)['rows'] if x['token']==token)
    factor=row['factor'];assert factor['status']=='PASS';f,pts,_=base.model_data(token);n=ZZ(1)
    for p,e in factor['factors']:assert ZZ(p).is_prime(proof=True);n*=ZZ(p)**e
    assert n==abs(16*f.discriminant());primes=[p for p,e in factor['factors']]
    return row,f,pts,primes,pari.nfinit([pari(f),primes])


def worker(token):
    from sage.all import QQ,ZZ,AA,GF,matrix,pari
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old,f,pts,primes,nf=setup(token);S=old['local']['S_finite'];m=len(pts);k=old['local']['strict_generic_dimension']
    theta=pari.Mod('z',pari(f));basis=list(nf.nf_get_zk());assert basis[0]==1
    enc=lambda b:[str(pari.lift(b).polcoef(i)) for i in range(3)]
    coords=[v for v in product((-1,0,1),repeat=3) if any(v) and next(x for x in v if x)!=-1 and v!=(1,0,0)]
    assert len(coords)==12
    candidates=[];betas=[pari(x)-theta for x,y in pts]
    for v in coords:
        alpha=sum(v[i]*basis[i] for i in range(3));alpha=pari.Mod(alpha,pari(f));N=ZZ(pari.nfeltnorm(nf,alpha));assert N
        beta=pari(N)*alpha;assert pari.nfeltnorm(nf,beta)==N**4
        candidates.append({'basis_coordinates':list(v),'alpha_ascending':enc(alpha),'norm_alpha':str(N),'beta_ascending':enc(beta)})
        betas.append(beta)
    constraints=[];locals=[]
    for p in S:
        chars=LocalSquareclasses(nf,p);sigs=[list(chars.signature(b)) for b in betas]
        assert sigs[:m]==next(x['signatures'] for x in old['local']['local'] if x['place']==p)
        constraints.extend(list(zip(*sigs)));locals.append({'place':p,'signatures':sigs})
    roots=f.roots(AA,multiplicities=False)
    signs=[[int(f.parent()(list(map(QQ,enc(b))))(a)<0) for a in roots] for b in betas]
    constraints.extend(list(zip(*signs)));locals.append({'place':'infinity','signatures':signs})
    initial_kernel=matrix(GF(2),constraints).right_kernel();outside=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['probe_prime_bound']):
        if p in S or not any(int(c['norm_alpha'])%p==0 for c in candidates):continue
        for j,P in enumerate(pari.idealprimedec(nf,p)):
            vals=[int(pari.idealval(nf,b,P)) for b in betas];assert all(v%2==0 for v in vals[:m])
            constraints.append([v%2 for v in vals]);outside.append({'prime':p,'prime_index':j,'valuations':vals})
    kernel=matrix(GF(2),constraints).right_kernel();assert int(matrix(GF(2),[x[:m] for x in constraints]).right_kernel().dimension())==k
    selected=[];projected=[]
    for v in kernel.basis():
        projection=r.pack(v[m:])
        if r.rank(projected+[projection])==len(projected):continue
        projected.append(projection);b=pari.Mod(1,pari(f));mask=r.pack(v)
        for bit,beta in zip(v,betas):
            if bit:b*=beta
        selected.append({'coefficient_mask':mask,'beta_ascending':enc(b),'is_global_square':bool(pari.nfeltissquare(nf,b))})
    assert len(selected)==int(kernel.dimension())-k
    # Good-prime field characters on the generic subgroup and retained candidate products.
    all_coeff=[[QQ(pari.lift(b).polcoef(i)) for i in range(3)] for b in betas[:m]]
    all_coeff += [list(map(QQ,b['beta_ascending'])) for b in selected]
    signatures=[0]*len(all_coeff);places=[];offset=0
    for p in r.primes(r.read(PROTOCOL)['limits']['probe_prime_bound']):
        if p==2 or int(f.discriminant())%p==0:continue
        roots_p=r.roots_at(str(f[1]),str(f[0]),p)
        if not roots_p:continue
        if any(int(c.denominator())%p==0 for coeff in all_coeff for c in coeff):continue
        residues=[[(sum(int(c.numerator())*pow(int(c.denominator()),-1,p)*pow(x,i,p) for i,c in enumerate(coeff)))%p for x in roots_p] for coeff in all_coeff]
        if any(v==0 for row in residues for v in row):continue
        for j,x in enumerate(roots_p):
            for i,row in enumerate(residues):signatures[i]|=int(pow(row[j],(p-1)//2,p)==p-1)<<offset
            places.append([p,x]);offset+=1
        if r.rank(signatures[:m])==m and r.rank(signatures)==m+sum(not b['is_global_square'] for b in selected):break
    assert r.rank(signatures[:m])==m
    return {'status':'PASS','token':token,'maximal_order_basis_ascending':[enc(pari.Mod(b,pari(f))) for b in basis],
        'candidates':candidates,'generic_dimension':m,'generic_strict_dimension':k,'local':locals,'outside_S_probes':outside,
        'local_only_kernel_dimension':int(initial_kernel.dimension()),'final_kernel_dimension':int(kernel.dimension()),
        'final_kernel_masks':[r.pack(v) for v in kernel.basis()],
        'candidate_coefficient_quotient_dimension':len(selected),'survivors':selected,
        'character_places':places,'character_signatures':signatures,
        'certified_global_squareclass_increment_beyond_G':r.rank(signatures)-m,
        'complete_additional_strict_classes':'UNKNOWN' if any(not b['is_global_square'] for b in selected) else 0,
        'boundary':'Finite signatures certify only global squareclass independence. Any nonsquare survivor still needs full outside-S support certification and rational solubility.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in cases():
        path=WORK/f'{token}.json'
        if not path.exists():
            with (WORK/f'{token}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['worker_seconds_per_case'])
                    error=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:error='bounded timeout'
            if error:r.write_new(path,{'bindings':bindings(),'token':token,'status':'UNKNOWN','reason':error})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(token,row['status'],row.get('candidate_coefficient_quotient_dimension'),row.get('certified_global_squareclass_increment_beyond_G'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-norm-projection.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        r.write_new(WORK/f'{args.token}.json',{'bindings':bindings(),**worker(args.token)})
