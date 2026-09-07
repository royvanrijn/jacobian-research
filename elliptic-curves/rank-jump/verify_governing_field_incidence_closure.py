#!/usr/bin/env python3
"""Independent quotient-module proof and replay of the 48 local exclusions."""
import argparse
from itertools import product
from pathlib import Path
import retrospective as r
import governing_field_incidence_closure as run

OUTPUT=r.OUT/'rank_jump_governing_field_incidence_closure_verification_v1.json'


def compute():
    inp=r.read(run.INPUT);out=r.read(run.OUTPUT)
    for obj in (inp,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    for path in (run.panel_verifier.OUTPUT,run.SUPPLEMENT_REPLAY):assert r.read(path)['status']=='PASS'
    V=list(product((0,1),repeat=2));zero=(0,0);I=(1,0,0,1)
    G=[g for g in product((0,1),repeat=4) if (g[0]*g[3]+g[1]*g[2])%2]
    add=lambda a,b:tuple(x^y for x,y in zip(a,b))
    act=lambda g,v:tuple(sum(g[2*i+j]*v[j] for j in range(2))%2 for i in range(2))
    e=lambda a,b:(a[0]*b[1]+a[1]*b[0])%2
    def mul(x,y):
        a,b,z=x;c,d,w=y
        return add(a,c),add(b,d),z^w^e(a,d)
    P=list(product(V,V,(0,1)));unit=(zero,zero,0);center=(zero,zero,1)
    inverse={x:next(y for y in P if mul(x,y)==unit and mul(y,x)==unit) for x in P}
    comm={mul(mul(mul(x,y),inverse[x]),inverse[y]) for x in P for y in P}
    assert comm=={unit,center} and {mul(x,x) for x in P}=={unit,center}
    homs=[]
    for entries in product((0,1),repeat=8):
        phi=lambda a,b:tuple(sum(entries[4*i+j]*v for j,v in enumerate(a+b))%2 for i in range(2))
        if all(phi(act(g,a),act(g,b))==act(g,phi(a,b)) for g in G for a in V for b in V):homs.append(entries)
    assert len(homs)==4
    mat=lambda g,h:tuple(sum(g[2*i+k]*h[2*k+j] for k in range(2))%2 for i in range(2) for j in range(2))
    multiplication=[[G.index(mat(g,h)) for h in G] for g in G]
    cocycles=[]
    for vals in product(V,repeat=6):
        if all(vals[multiplication[i][j]]==add(vals[i],act(G[i],vals[j])) for i in range(6) for j in range(6)):cocycles.append(vals)
    coboundaries=[tuple(add(act(g,v),v) for g in G) for v in V]
    assert set(cocycles)==set(coboundaries) and len(cocycles)==4
    # Direct semidirect-product extension of each equivariant hom gives all four H1 classes.
    assert out['finite_group']['H1_dimension']==2 and out['finite_group']['Z1_dimension']==4
    base={x['token']:x for x in r.read(run.panel.OUTPUT)['rows']}
    patch={x['token']:x for x in r.read(run.completion.OUTPUT)['rows']}
    supp={x['token']:x for x in r.read(run.supplement.OUTPUT)['rows']}
    count=0
    for row,result in zip(inp['rows'],out['rows']):
        token=row['token'];assert token==result['token']
        source=[]
        for local in (base[token]['local'],patch.get(token,{}),supp.get(token,{}).get('local',{})):
            source.extend(local.get('local',[]))
        assert [w['word'] for w in row['local_witnesses']]==[1,2,3]
        for w in row['local_witnesses']:
            assert any(b['place']==w['place'] and b['signatures'][:2]==w['pair_signatures'] for b in source)
            a,b=w['pair_signatures'];sig=[((a[i] if w['word']&1 else 0)+(b[i] if w['word']&2 else 0))%2 for i in range(len(a))]
            assert any(sig) and r.pack(sig)==w['nonzero_signature'];count+=1
        assert 6*4**2==row['joint_class_field_degree']<row['governing_field_degree']<6*4**3
        assert result['pair_governing_strict_class_dimension']==result['pair_governing_new_class_dimension']==0
    assert count==48
    return {'schema':'rank-jump.governing-field-incidence-closure-verification.v1','status':'PASS',
        'kernel_order':32,'kernel_frattini_and_commutator_order':2,
        'equivariant_homomorphisms_to_V':len(homs),'S3_one_cocycles':len(cocycles),'S3_H1_dimension':0,
        'governing_H1_dimension':2,'local_exclusion_witnesses_verified':48,'fibres':16,
        'third_class_degree_exclusion_verified':True,
        'method':'Independent tuple kernel, all 256 linear maps, all 4096 S3 cocycles, retained local-character provenance and XOR replay.',
        'bindings':run.bindings([Path(__file__),run.INPUT,run.OUTPUT,run.panel_verifier.OUTPUT,run.SUPPLEMENT_REPLAY,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent H1 calculation and all 48 strict-class exclusions')
