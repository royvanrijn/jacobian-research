#!/usr/bin/env python3
"""Independent exact local-square and Sage matrix replay, no local enumeration."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import production_minus_twist as first
import production_minus_twist_completion as completed
import scalar_cup

OUTPUT=r.OUT/'rank_jump_production_minus_twist_verification_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-production-minus-twist-verification-v1-retry2'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),first.INPUT,first.OUTPUT,completed.INPUT,completed.OUTPUT,scalar_cup.OUTPUT)}


def verify(index):
    from sage.all import QQ, GF, AA, PolynomialRing, VectorSpace, pari
    from sage.version import version
    sys.path.insert(0,str(first.rem.bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    assert r.read(first.INPUT)['bindings']==first.bindings()
    case=next(x for x in completed.merged() if x['case_index']==index)
    assert case['complete']
    old=r.read(first.rem.bad.INPUT)['cases'][index]; source=first.rem.bad.cases()[index]
    _,allpoints=r.short(source['model'],source['generic_points']+source['points'])
    points=[allpoints[i] for i in old['selected_input_indices']];m=old['generic_dimension']
    R=PolynomialRing(QQ,'z');f=R(list(map(QQ,case['integral_cubic_ascending'])))
    primes=[x['place'] for x in case['local'] if x['place']!='infinity']
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f));scale=QQ(old['elliptic_scaling_d'])
    betas=[pari(QQ(P[0])*scale**2)-theta for P in points]
    derivative=-pari(f.discriminant())*pari(f.derivative())(theta);betas.append(derivative)
    roots=f.roots(AA,multiplicities=False);local_records=[]
    joint=[0]*len(betas);original=[];twist=[];offset=0
    for local in case['local']:
        p=local['place'];basis=local['twist_basis'];dim=local['point_dimension'];width=local['width']
        V=VectorSpace(GF(2),width)
        vec=lambda x:V([(x>>i)&1 for i in range(width)])
        def span(xs):return V.subspace([vec(x) for x in xs])
        assert span(basis).dimension()==dim
        if p=='infinity':
            signatures=[r.pack(int(QQ(P[0])*scale**2<a) for a in roots) for P in points]
            signatures.append(r.pack(int(-f.discriminant()*f.derivative()(a)<0) for a in roots))
            assert basis==([3] if len(roots)==3 else [])
            verified_witnesses=0;independent_nonzero_subsets=0
        else:
            chars=LocalSquareclasses(nf,p)
            signatures=[r.pack(chars.signature(beta)) for beta in betas]
            assert dim==len(chars.primes)-1+(p==2)
            representatives=[]
            if p%4==1:
                assert int(pari.issquare(pari(-1)+pari(f'O({p}^8)')))==1
                selected=[]
                for i,sig in enumerate(signatures[:m]):
                    if vec(sig) not in span([signatures[j] for j in selected]):selected.append(i)
                for sig in basis:
                    # A basis vector is a combination, not necessarily an input column.
                    mask=first.lc.coordinates(sig,[signatures[j] for j in selected])
                    beta=pari.Mod(1,pari(f))
                    for i,j in enumerate(selected):
                        if mask>>i&1:beta*=betas[j]
                    representatives.append(beta)
                verified_witnesses=0
            else:
                assert len(local['witnesses'])==dim
                for witness in local['witnesses']:
                    x=QQ(witness['short_x']);value=x**3+f[1]*x-f[0]
                    assert value==QQ(witness['cubic_value']) and value!=0
                    # PARI p-adic square testing, separate from the producer's
                    # rational valuation/unit-residue criterion.
                    valuation=int(value.valuation(p));unit=value/QQ(p)**valuation
                    assert valuation%2==0 and int(pari.issquare(pari(unit)+pari(f'O({p}^8)')))==1
                    assert first.local_square(value,p)==witness['square_witness']
                    beta=pari(x)+theta
                    assert r.pack(chars.signature(beta))==witness['signature']
                    representatives.append(beta)
                assert [r.pack(chars.signature(x)) for x in representatives]==basis
                verified_witnesses=len(representatives)
            # nfislocalpower is a different API from LocalSquareclasses' residue logs.
            # Every nonempty subset must be nonsquare in at least one completion.
            for mask in range(1,1<<dim):
                value=pari.Mod(1,pari(f))
                for j,beta in enumerate(representatives):
                    if mask>>j&1:value*=beta
                den=pari.denominator(pari.nfalgtobasis(nf,value));value*=den**2
                assert any(int(pari.nfislocalpower(nf,P,value,2))==0 for P in chars.primes)
            independent_nonzero_subsets=(1<<dim)-1
        assert signatures==local['global_signatures']
        assert span(signatures[:m])==span(local['original_basis'])
        assert span(local['original_basis']).intersection(span(basis)).dimension()==local['intersection_dimension']
        for j,sig in enumerate(signatures):joint[j]|=sig<<offset
        original.extend(v<<offset for v in local['original_basis']);twist.extend(v<<offset for v in basis)
        offset+=width
        local_records.append({'place':p,'local_square_witnesses':verified_witnesses,
            'independent_nonzero_subsets':independent_nonzero_subsets,'status':'PASS'})
    assert joint==case['joint_global_signatures'] and original==case['original_product_basis'] and twist==case['twist_product_basis']
    V=VectorSpace(GF(2),offset)
    def span(xs):return V.subspace([V([(x>>i)&1 for i in range(offset)]) for x in xs])
    Lambda=span(joint);L0=span(original);Lm=span(twist)
    A=Lambda.intersection(L0);B=Lambda.intersection(Lm);C=A.intersection(B)
    assert Lambda.dimension()==case['relaxed_boundary_dimension']==L0.dimension()==Lm.dimension()
    assert A==span(joint[:-1]) and B.is_subspace(A)
    generic_derivative=span(joint[:m]+[joint[-1]])
    expected=next(x for x in r.read(completed.OUTPUT)['rows'] if x['case_index']==index)
    assert B.dimension()==expected['twist_boundary_dimension'] and C.dimension()==expected['common_boundary_dimension']
    k=case['known_strict_dimension']
    cup=next(x for x in r.read(scalar_cup.OUTPUT)['production_cases'] if x['case_index']==index)
    from sage.all import matrix
    M=matrix(GF(2),cup['scalar_cup_matrix']);ct_rank=int(M.rank())
    assert M==M.transpose() and not any(M[i,i] for i in range(M.nrows()))
    assert ct_rank==cup['detected_scalar_cup_rank']
    assert cup['strict_dimension']==k
    drop=int(A.dimension()-B.dimension())
    local_cost=L0.dimension()-L0.intersection(Lm).dimension()
    assert drop==local_cost
    return {'bindings':bindings(),'case_index':index,'id':case['id'],'status':'PASS',
        'software':{'sage':version,'pari':str(pari.version())},'local':local_records,
        'full_twist_Selmer_is_subspace_of_original':True,'Selmer_dimension_drop':drop,
        'full_twist_Selmer_dimension':f'{k+B.dimension()} + epsilon',
        'twist_Sha_2_dimension_lower_bound':ct_rank,
        'full_twist_rank_upper_bound':f'{k+B.dimension()-ct_rank} + epsilon',
        'rank_drop_lower_bound':f'{drop+ct_rank} - dim Sha(original)[2]',
        'generic_and_derivative_boundary_dimension':int(generic_derivative.dimension()),
        'full_relaxed_boundary_dimension':int(Lambda.dimension()),
        'boundary_certificate_needs_exceptional_points':generic_derivative!=Lambda,
        'boundary':'epsilon is the same unknown strict excess for each original/twist pair; separate across cases. No numerical twist rank or original exact rank is claimed.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in [0,4,5]:
        path=WORK/f'case-{i}.json'
        # Check mode independently reruns in memory, without touching artifacts.
        if check:
            expected=next(x for x in r.read(OUTPUT)['rows'] if x['case_index']==i)
            assert verify(i)==expected;print('PASS independent replay',i,flush=True);continue
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker',
                    '--index',str(i),'--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=60)
                if proc.returncode:raise RuntimeError(f'independent verification failed; see {path.with_suffix(".log")}')
        row=r.read(path);assert row['bindings']==bindings();rows.append(row);print('PASS',i,flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.production-minus-twist-verification.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);p.add_argument('--index',type=int)
    p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='worker':r.write_new(args.destination,verify(args.index))
    else:capture(args.mode=='check')
