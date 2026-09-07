#!/usr/bin/env python3
"""Independent norm-one identities, low-pivot parity and finite characters."""
import argparse
from pathlib import Path
import retrospective as r
import retained_radical_closure as run

OUTPUT=r.OUT/'rank_jump_retained_radical_closure_verification_v1.json'


def low_pivot_rank(rows):
    piv={}
    for v in rows:
        while v:
            p=v&-v
            if p not in piv:piv[p]=v;break
            v^=piv[p]
    return len(piv)


def compute():
    from sage.all import QQ,GF,PolynomialRing,matrix
    out=r.read(run.OUTPUT);data=r.read(run.pool.INPUT);inp=r.read(run.root.INPUT);affine=r.read(run.affine.OUTPUT)
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    blocks={}
    for i,col in enumerate(data['columns']):blocks.setdefault(col['p'],[]).append(i)
    S=set(data['S_finite']);rows=[]
    for rec in data['relations']:
        vals=dict(rec['ideal_factorization']);mask=0
        for p,ids in blocks.items():
            if p in S:continue
            if not any(i in vals for i in ids):continue
            n=sum(data['columns'][i]['f']*vals.get(i,0) for i in ids)
            for i in ids:
                # Valuation of a^3/Norm(a), independently of the pi(a) worker formula.
                v=3*vals.get(i,0)-data['columns'][i]['e']*n
                if v%2:mask^=1<<i
        rows.append(mask)
    j=out['replaced_dictionary_index'];assert rows[j]==0
    assert [data['relations'][j][k] for k in ('m','n')]==out['parent_address']
    selected=rows[:j]+rows[j+1:]
    assert low_pivot_rank(selected)==4133
    target=sum((v%2)<<i for i,v in enumerate(affine['projection_valuations_by_column']) if v is not None)
    assert low_pivot_rank(selected+[target])==4134
    f=list(map(QQ,inp['cubic_ascending']));assert f[-1]==1
    C=matrix(QQ,[[0,0,-f[0]],[1,0,-f[1]],[0,1,-f[2]]]);I=matrix.identity(QQ,3)
    def elt(cs):return sum((QQ(c)*C**i for i,c in enumerate(cs)),matrix(QQ,3,3))
    def D(M):return M**3/M.det()
    a=elt(inp['alpha_ascending']);w=elt(inp['root_ascending']);gammas=[elt(c) for c in inp['generic_classes_ascending']]
    u=elt(out['norm_one_replaced_generator']);h=elt(out['norm_one_root']);vs=[elt(c) for c in out['norm_one_generic_generators']]
    assert u==D(a) and h==D(w) and vs==list(map(D,gammas))
    assert u.det()==h.det()==1 and all(v.det()==1 for v in vs)
    mask=inp['generic_product_mask'];rhs=u
    for i,v in enumerate(vs):
        if (mask>>i)&1:rhs*=v
    assert rhs==h*h
    # The old basis in the replacement basis: u_j=h^2*product(v_i)^-1.
    inclusion=matrix.identity(QQ,17);inclusion[0,0]=2
    for i in range(16):inclusion[0,i+1]=-int((mask>>i)&1)
    assert inclusion.det()==2 and inclusion.change_ring(GF(2)).rank()==16
    sigs=[0]*16
    for b,block in enumerate(out['generic_independence']['blocks']):
        p=block['p'];R=PolynomialRing(GF(p),'z');pol=R(f)
        roots=sorted(map(int,pol.roots(multiplicities=False)))
        assert roots==block['roots'] and len(roots)==3 and pol.discriminant()!=0
        for i,cs in enumerate(inp['generic_classes_ascending']):
            vals=[R(list(map(QQ,cs)))(x) for x in roots];assert all(vals)
            bits=[int(not v.is_square()) for v in vals];assert sum(bits)%2==0
            value=r.pack(bits);assert value==block['signatures'][i];sigs[i]|=value<<(3*b)
    assert sigs==out['generic_independence']['signatures'] and low_pivot_rank(sigs)==16
    return {'schema':'rank-jump.retained-radical-closure-verification.v1','status':'PASS',
        'dictionary_rows_rebuilt_from_norm_one_valuations':len(rows),
        'independent_low_pivot_ranks':[4133,4134],
        'generic_finite_character_blocks':len(out['generic_independence']['blocks']),
        'generic_character_rank':16,'norm_one_generator_identities':18,
        'exact_norm_one_root_relation':True,'basis_inclusion_determinant':2,
        'replacement_squareclass_rank':4134+16,
        'saturation_scope':'The accompanying torsion-free group proof promotes the exact index and squareclass independence to full 2-saturation and exhausts the squareclass image of all in-field radical operations.',
        'bindings':run.bindings([Path(__file__),run.OUTPUT,run.pool.INPUT,run.root.INPUT,run.affine.OUTPUT,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent parity ranks, norm-one relation and index-two inclusion')
