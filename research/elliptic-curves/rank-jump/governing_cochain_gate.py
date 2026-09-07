#!/usr/bin/env python3
"""Finite S3 hypotheses and a complete two-class governing-cochain model."""
import argparse
from collections import Counter
from itertools import product
from pathlib import Path
import retrospective as r
from fibre_discrimination import hash_file
from quadratic_norm_blocks import cubic_certificate

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_quadratic_norm_block_inputs_v1.json'
PRIOR=r.OUT/'rank_jump_quadratic_norm_blocks_v1.json'
VERIFIED=r.OUT/'rank_jump_quadratic_norm_block_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_governing_cochain_gate_v1.json'


def mm(a,b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%2 for i in range(2) for j in range(2))


def act(g,x):
    a,b=x&1,x>>1
    return ((g[0]*a+g[1]*b)%2)|(((g[2]*a+g[3]*b)%2)<<1)


def weil(x,y):return ((x&1)*(y>>1)+(x>>1)*(y&1))%2


def compute():
    inp=r.read(INPUT);prior=r.read(PRIOR);verified=r.read(VERIFIED)
    assert prior['status']==verified['status']=='PASS'
    for doc in (prior,verified):
        for path,sha in doc['bindings'].items():assert hash_file(r.ROOT/path)==sha
    for path,sha in inp['source_hashes'].items():assert hash_file(r.ROOT/path)==sha
    mats=list(product(range(2),repeat=4));G=[g for g in mats if (g[0]*g[3]-g[1]*g[2])%2]
    I=(1,0,0,1);zero=(0,0,0,0);identity=G.index(I)
    gt=[[G.index(mm(g,h)) for h in G] for g in G]
    central=[g for g in mats if all(mm(g,h)==mm(h,g) for h in G)]
    assert central==[zero,I]
    assert all({act(g,v) for g in G}=={1,2,3} for v in (1,2,3))
    fixed_free=[i for i,g in enumerate(G) if all(act(g,v)!=v for v in (1,2,3))]
    assert len(fixed_free)==2 and all(gt[gt[i][i]][i]==identity for i in fixed_free)
    cocycles=[]
    for f in product(range(4),repeat=6):
        if all(f[gt[i][j]]==f[i]^act(G[i],f[j]) for i in range(6) for j in range(6)):cocycles.append(f)
    coboundaries=[tuple(act(g,v)^v for g in G) for v in range(4)]
    assert set(cocycles)==set(coboundaries) and len(cocycles)==4
    # Gamma=(V x V) semidirect S3, and omega=a cup b.
    base=list(product(range(4),range(4),range(6)));idx={x:i for i,x in enumerate(base)};N=len(base)
    bt=[];omega=[]
    for a,b,g in base:
        bt.append([idx[(a^act(G[g],aa),b^act(G[g],bb),gt[g][h])] for aa,bb,h in base])
        omega.append([weil(a,act(G[g],bb)) for aa,bb,h in base])
    for i in range(N):
        for j in range(N):
            ij=bt[i][j]
            for k in range(N):
                jk=bt[j][k]
                assert bt[ij][k]==bt[i][jk]
                assert omega[i][j]^omega[ij][k]==omega[j][k]^omega[i][jk]
    def mult(x,y):
        i,z=x//2,x%2;j,w=y//2,y%2
        return 2*bt[i][j]+(z^w^omega[i][j])
    unit=2*idx[(0,0,identity)];H=range(2*N)
    inverses=[]
    for x in H:
        candidates=[y for y in H if mult(x,y)==unit and mult(y,x)==unit];assert len(candidates)==1
        inverses.append(candidates[0])
    def psi(x):
        a,b,g=base[x//2];assert g in fixed_free
        p=next(v for v in range(4) if act(G[g],v)^v==a)
        return weil(p,b)^(x%2)
    fixed=[x for x in H if base[x//2][2] in fixed_free]
    for x in fixed:
        assert psi(x^1)!=psi(x)
        for h in H:assert psi(mult(mult(h,x),inverses[h]))==psi(x)
    counts=dict(sorted(Counter(str(psi(x)) for x in fixed).items()));assert counts=={'0':32,'1':32}
    rows=[]
    for old in inp['rows']:
        cert=cubic_certificate(old['cubic_ascending'],503);assert cert['group']=='S3'
        k=len(old['CT_matrix']);assert k in (6,8,10)
        target=[[int(i!=j and i//2==j//2) for j in range(k)] for i in range(k)]
        assert r.rank([r.pack(row) for row in target])==k
        rows.append({'id':old['id'],'cubic_certificate':cert,'strict_rational_block_dimension':k,
            'proposed_non_degenerate_CT_form':target,'proposed_CT_rank':k,
            'existential_twist_conclusion':'Morgan Theorem5.8 permits this form on the retained block while preserving all local Kummer images and the full Selmer group. No twist character or governing number field is computed.',
            'explicit_twist':'UNKNOWN','original_full_rank':'UNKNOWN'})
    return {'schema':'rank-jump.governing-cochain-gate.v1','status':'PASS',
        'S3_hypotheses':{'matrices':[list(g) for g in G],'endomorphism_centralizer':[list(g) for g in central],
            'simple_module':True,'fixed_point_free_element_indices':fixed_free,
            'one_cocycles':[list(f) for f in cocycles],'one_coboundaries':[list(f) for f in coboundaries],'H1_dimension':0},
        'two_class_model':{'semidirect_order':N,'central_extension_order':2*N,'two_cocycle_checks':N**3,
            'fixed_point_free_lifts':len(fixed),'governing_value_counts':counts,'conjugacy_checks':len(fixed)*len(H),
            'central_toggle_changes_value':True,'omega':'<a,g*b_next>','psi':'<(g-1)^(-1)*a,b>+z'},
        'rows':rows,
        'external_theorem':{'url':'https://arxiv.org/pdf/2309.02374v2','sections':['2 setup','Lemma2.8','Proposition3.3','Section5.2 Theorem5.8'],
            'paper_local_sha256':'ec4d757f37167d4c42a43580ddd3c208e6386eb6e16d15cc45fa02d7128835e1','scope':'Section2 resets A to an arbitrary principally polarized abelian variety; the Section1 Kummer-surface dimension restriction is not used by Section5.2.'},
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,PRIOR,VERIFIED,Path(__file__),HERE/'quadratic_norm_blocks.py',HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
        'boundary':'Finite group verification plus a written application of an external theorem. This is neither an explicit arithmetic governing extension nor an executed twist experiment. CT vanishing means second-descent lifting, not rational points. Existing original-family generic sections need not remain rational after twisting.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: H1(S3,V)=0; End=F2; 884736 cocycle checks; 12288 conjugacy checks; 3 production applicability gates')
