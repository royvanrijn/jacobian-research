#!/usr/bin/env sage-python
"""Bounded effective three-section pencils on the explicit u11 parent."""
import hashlib,json
from itertools import combinations,product
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-degree-three-pencils-v1'
OUT=D/'roster.json'
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
    heights=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    G=matrix(QQ,ns['rational_NS_Gram']);H=matrix(QQ,heights['seed_height_gram'])
    std=matrix.identity(QQ,18).rows();F,O=std[:2];C=G[2:7,2:7]
    phi=[]
    for i,p in enumerate(ns['section_profiles']):
        v=std[7+i]-O-(2+p['zero_section_intersection'])*F
        v[2:7]-=C.inverse()*vector(QQ,G[2:7,7+i].column(0));phi.append(v)
    words=[]
    for size in range(1,4):
        for indices in combinations(range(11),size):
            for signs in product((-1,1),repeat=size):
                w=vector(ZZ,11)
                for i,s in zip(indices,signs):w[i]=s
                words.append(w)
    assert len(words)==1562
    sections=[]
    for w in words:
        comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
        correction=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4
        oo=(w*H*w-4+correction)/2
        assert oo in ZZ and oo>=0
        v=O+(2+oo)*F+sum((w[i]*phi[i] for i in range(11)),vector(QQ,18))
        cross=vector(QQ,comp[:2]+[int(comp[2]==k) for k in (1,2,3)])
        v[2:7]+=C.inverse()*cross
        assert v*G*v==-2 and v*G*O==oo and v*G*F==1 and all(z in ZZ for z in v)
        sections.append({'word':list(map(int,w)),'class':list(map(int,v)),'O_intersection':int(oo),'components':list(map(int,comp))})
    curves=[O]+std[2:7]+[F-std[2],F-std[3],F-sum(std[4:7])]+[vector(QQ,r['class']) for r in sections]
    V=matrix(QQ,curves);VG=V*G
    eligible=[i for i,r in enumerate(sections) if r['O_intersection']==1]
    rows=[]
    for i,j in combinations(eligible,2):
        P,Q=[vector(QQ,sections[k]['class']) for k in (i,j)]
        if P*G*Q!=1:continue
        fibre=O+P+Q;assert fibre*G*fibre==0 and fibre*G*F==3
        cross=VG*fibre;assert min(cross)>=0
        ones=[k for k,v in enumerate(cross) if v==1]
        if not ones:continue
        zero=[k for k,v in enumerate(cross) if v==0]
        vertical_rank=V.matrix_from_rows(zero).rank();upper=17-vertical_rank
        # Intersection one with an actual rational curve proves primitivity
        # and a degree-one multisection. O,P,Q are effective, so nefness holds:
        # each component has zero intersection and every other curve >=0.
        rows.append({'section_indices':[i,j],'words':[sections[k]['word'] for k in (i,j)],
                     'fibre_class':list(map(int,fibre)), 'section_witness_curve_index':ones[0],
                     'vertical_curve_indices':zero,'vertical_class_rank':int(vertical_rank),
                     'generic_Q_MW_upper_bound_from_visible_curves':int(upper)})
    rows.sort(key=lambda r:(-r['generic_Q_MW_upper_bound_from_visible_curves'],sum(abs(z) for w in r['words'] for z in w),r['words']))
    result={'schema':'elliptic-curves.mestre-degree-three-roster.v1','status':'PASS','outer_u':'11',
            'signed_word_count':len(words),'O_intersection_one_sections':len(eligible),
            'eligible_pencils':len(rows),'sections':sections,'visible_curve_classes':[[int(z) for z in v] for v in curves],
            'rows':rows,'selected_row':0 if rows else None,'sources':protocol['sources'],
            'scope':'Fixed signed support-at-most-three section words, effective D=O+P+Q with all three pairwise intersections1 and a visible rational degree-one multisection. Exact bounded class roster only; no all-degree-three completeness, pairwise inequivalence, new parent or point search. Visible upper bounds are not attained-rank claims.'}
    with OUT.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS signed words',len(words),'O.P1',len(eligible),'pencils',len(rows),'upper histogram',
          {u:sum(r['generic_Q_MW_upper_bound_from_visible_curves']==u for r in rows) for u in sorted(set(r['generic_Q_MW_upper_bound_from_visible_curves'] for r in rows))},flush=True)
    if rows:print('SELECTED',rows[0],flush=True)
if __name__=='__main__':main()
