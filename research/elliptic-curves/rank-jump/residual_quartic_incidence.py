#!/usr/bin/env python3
"""Exact generic theta pairing and its residual quartic Galois scheme."""
import argparse
from collections import Counter
from itertools import permutations
from pathlib import Path
import retrospective as r
import branch_divisibility_capacity as source

PROTOCOL=Path(__file__).with_name('RESIDUAL_QUARTIC_INCIDENCE_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_residual_quartic_incidence_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_residual_quartic_incidence_v1.json'


def bindings(paths):
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    old=r.read(source.INPUT)
    r.write_new(INPUT,{'schema':'rank-jump.residual-quartic-incidence-inputs.v1',
        **{k:old[k] for k in ('A','B','sections')},
        'bindings':bindings([Path(__file__),PROTOCOL,source.INPUT]),
        'boundary':'Only generic polynomials. No specialized points, parameters, covers or ranks.'})


def q4(x):
    a,b,c,d=[(x>>i)&1 for i in range(4)]
    return (a*b+c+c*d+d)%2


def action(cols,x):
    ans=0
    for j,col in enumerate(cols):
        if (x>>j)&1:ans^=col
    return ans


def residual_group():
    e=lambda a,b:q4(a^b)^q4(a)^q4(b)
    base=[1,2,4,8]
    orth=[]
    for cols in permutations(range(1,16),4):
        if any(q4(v)!=q4(b) for v,b in zip(cols,base)):continue
        if any(e(cols[i],cols[j])!=e(base[i],base[j]) for i in range(4) for j in range(i)):continue
        assert r.rank(cols)==4
        assert all(q4(action(cols,x))==q4(x) for x in range(16))
        orth.append(cols)
    iso=[x for x in range(1,16) if q4(x)==0];assert len(iso)==5 and len(orth)==120
    radical=1;four=[x for x in iso if x!=radical]
    stab=[c for c in orth if c[0]==radical];assert len(stab)==24
    perm_to_cols={tuple(four.index(action(c,x)) for x in four):c for c in stab}
    assert len(perm_to_cols)==24
    perms=sorted(perm_to_cols);indices={p:i for i,p in enumerate(perms)};identity=indices[(0,1,2,3)]
    mult=[[indices[tuple(p[q[i]] for i in range(4))] for q in perms] for p in perms]
    def cycles(p):
        seen=set();out=[]
        for i in range(4):
            if i in seen:continue
            j=i;n=0
            while j not in seen:seen.add(j);n+=1;j=p[j]
            out.append(n)
        return tuple(sorted(out))
    rows=[]
    for p in perms:
        c=perm_to_cols[p];fixed=[x for x in range(16) if action(c,x)==x]
        rows.append({'permutation':list(p),'columns':list(c),'cycles':list(cycles(p)),
            'fixed_dimension':r.rank(fixed),'trace_mod_two':sum((c[i]>>i)&1 for i in range(4))%2,
            'rank_minus_identity':r.rank([c[i]^base[i] for i in range(4)])})
    def close(gens):
        ans={identity};pending=[identity]
        while pending:
            a=pending.pop()
            for g in gens:
                b=mult[a][g]
                if b not in ans:ans.add(b);pending.append(b)
        return frozenset(ans)
    trivial=frozenset([identity]);groups={trivial};pending=[trivial]
    while pending:
        h=pending.pop()
        for g in range(24):
            if g in h:continue
            k=close(list(h)+[g])
            if k not in groups:groups.add(k);pending.append(k)
    assert len(groups)==30
    admissible=[]
    for h in sorted(groups,key=lambda s:(len(s),sorted(s))):
        if not any(rows[i]['trace_mod_two']==1 for i in h):continue
        if not any(rows[i]['rank_minus_identity']==1 for i in h):continue
        fixed=[x for x in range(16) if all(action(perm_to_cols[perms[i]],x)==x for i in h)]
        rational=[x for x in four if x in fixed]
        assert len(h) in (6,24) and len(rational)==int(len(h)==6)
        assert r.rank(fixed)==1+len(rational)
        admissible.append({'elements':sorted(h),'group_order':len(h),'fixed_dimension':r.rank(fixed),
            'fixed_quartic_letters':rational,'global_pool_dimension':16+r.rank(fixed)})
    assert Counter(x['group_order'] for x in admissible)=={6:4,24:1}
    return {'quadratic_form':'a*b+c+c*d+d over F2','isotropic_nonzero':iso,
        'chosen_radical':radical,'quartic_letters':four,'orthogonal_group_order':len(orth),
        'radical_stabilizer_order':len(stab),'subgroups_checked':len(groups),
        'elements':rows,'admissible_galois_images':admissible}


def compute():
    from sage.all import QQ,GF,PolynomialRing,matrix
    data=r.read(INPUT);R=PolynomialRing(QQ,'t');t=R.gen();A=R(data['A']);B=R(data['B'])
    D=-4*A**3-27*B**2
    assert [A.degree(),B.degree(),D.degree()]==[8,12,24] and D.is_squarefree() and D.gcd(A)==1
    s=[(R(z['x']),R(z['y'])) for z in data['sections']];assert len(s)==17
    def coeff(p,n):return [p[i] for i in range(n)]
    singles=[];pairs=[];E=matrix(GF(2),17)
    for i,(x,y) in enumerate(s):
        assert y*y==x**3+A*x+B and x.degree()==4 and y.degree()==6
        assert y.is_squarefree() and y.gcd(3*x*x+A)==1
        # Multiply L(D_i+F) by x-x_i: degree-five polynomials plus x*degree-one polynomials, vanishing at D_i.
        cols=[coeff(t**j%y,6) for j in range(6)]+[coeff(t**j*x%y,6) for j in range(2)]
        dim=8-matrix(QQ,cols).rank();assert dim==2
        singles.append({'index':i,'theta_translate_h0':dim,'theta_value':dim%2})
        for j in range(i):
            xx,yy=s[j];assert y.gcd(yy)==1 and (x-xx).degree()==4
            # L(D_i+D_j-F) multiplied by (x-x_i)(x-x_j) lies in L(7F-D_i-D_j).
            cols=[]
            for k in range(12):
                a=t**k if k<8 else R(0);b=t**(k-8) if k>=8 else R(0)
                cols.append(coeff((a+b*x)%y,6)+coeff((a+b*xx)%yy,6))
            dim=12-matrix(QQ,cols).rank();E[i,j]=E[j,i]=dim%2
            pairs.append({'i':i,'j':j,'theta_translate_h0':dim,'weil_pairing_bit':dim%2})
    assert E.rank()==16
    radical=E.right_kernel().basis()[0];mask=r.pack(radical);assert mask==27687
    rows=[r.pack(row) for row in E.rows()]
    def q(v):return sum(((v&((1<<i)-1)&rows[i]).bit_count()&1) for i in range(17) if (v>>i)&1)%2
    assert q(mask)==0
    # Omit the first marked vector; the radical has nonzero first coefficient.
    H=E.matrix_from_rows_and_columns(range(1,17),range(1,17));assert H.rank()==16
    zeros=sum(q(v<<1)==0 for v in range(1<<16));assert zeros==32640
    return {'schema':'rank-jump.residual-quartic-incidence.v1','status':'PASS',
        'theta':{'divisor':'3F_infinity','h0':4,'arf':0},'singletons':singles,'pairs':pairs,
        'weil_pairing_rows':rows,'generic_pairing_rank':16,'generic_radical_mask':mask,
        'generic_radical_theta_value':0,'nondegenerate_generic_indices':list(range(1,17)),
        'nondegenerate_generic_zero_count':zeros,'nondegenerate_generic_arf':1,
        'residual_dimension':4,'residual_arf':1,'residual_group':residual_group(),
        'degree_four_etale_scheme':'{v in H_perp: v != 0, v != radical, q_theta(v)=0}',
        'coefficient_quartic':'UNKNOWN','rational_point_on_degree_four_scheme':'UNKNOWN',
        'global_pool_dimension':'17 (S4) or 18 (S3); actual alternative UNKNOWN',
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT,Path(r.__file__)]),
        'boundary':'Exact fixed-root-curve incidence reduction; no rational h, twist point, specialized class-group basis or candidate score constructed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);a=p.parse_args()
    if a.mode=='export':export()
    else:
        result=compute()
        if a.mode=='build':r.write_new(OUTPUT,result)
        else:assert result==r.read(OUTPUT)
        print('PASS 136 theta pairings; residual degree-four scheme; Galois image S3 or S4, alternative UNKNOWN')
