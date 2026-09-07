#!/usr/bin/env sage-python
"""Use exact NS classes to select bisection translates before equation work."""
import hashlib,json
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-bisection-pencil-admission-v1'
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
    data=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    G=matrix(QQ,ns['rational_NS_Gram']);inv=G.inverse();H=matrix(QQ,data['seed_height_gram']);std=matrix.identity(QQ,18).rows();F,O=std[:2];C=G[2:7,2:7]
    phi=[]
    for i,p in enumerate(ns['section_profiles']):
        v=std[7+i]-O-(2+p['zero_section_intersection'])*F
        v[2:7]-=C.inverse()*vector(QQ,G[2:7,7+i].column(0));phi.append(v)
    def section(w):
        if not any(w):return O,[0,0,0]
        comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
        corr=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4;oo=(w*H*w-4+corr)/2
        v=O+(2+oo)*F+sum((w[i]*phi[i] for i in range(11)),vector(QQ,18))
        cross=vector(QQ,comp[:2]+[int(comp[2]==k) for k in (1,2,3)]);v[2:7]+=C.inverse()*cross
        assert v*G*v==-2 and all(c in ZZ for c in v)
        return v,comp
    translations=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-translated-bisections-v1/result.json').read_text())
    options=[(i,j,a,b,c) for i,j in product(range(3),repeat=2) for a,b,c in product(range(3),repeat=3) if a+b+c<=2]
    bases=[]
    for source_index in (27,35):
        rows=[r for r in translations['rows'] if r['source_curve_index']==source_index]
        cross=[2,rows[0]['old_O_intersection'],*[rows[1+2*i]['old_O_intersection'] for i in range(11)]]
        possibilities=[]
        for comp in options:
            v=inv*vector(QQ,cross[:2]+list(comp)+cross[2:])
            if all(c in ZZ for c in v) and v*G*v==-2:possibilities.append((v,list(comp)))
        assert len(possibilities)==1
        v,comp=possibilities[0]
        for r in rows:
            S,_=section(-vector(ZZ,r['translation_word']));assert v*G*S==r['old_O_intersection']
        bases.append({'source_index':source_index,'class':list(map(int,v)),'components':list(map(int,comp))})
    domain=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-all-section-triangles-v1/height-domain.json').read_text())
    words=[s*vector(ZZ,w) for w in domain['section_words_up_to_sign'] for s in (-1,1)]
    sec=[section(w) for w in words]
    oldpool=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-degree-three-pencils-v1/roster.json').read_text())
    visible=[vector(QQ,v) for v in oldpool['visible_curve_classes']]
    V=matrix(QQ,visible);VG=V*G;rows=[]
    for base in bases:
        R=vector(QQ,base['class']);com=base['components'];full=[[2-com[0],com[0]],[2-com[1],com[1]],[2-sum(com[2:]),*com[2:]]]
        for w,(S,sc) in zip(words,sec):
            if R*G*S!=2:continue
            shifted=[[a[(i+c)%len(a)] for i in range(len(a))] for a,c in zip(full,sc)]
            cmp=[shifted[0][1],shifted[1][1],*shifted[2][1:]]
            cross=[R*G*section(w+vector(ZZ,[int(i==j) for i in range(11)]))[0] for j in range(11)]
            v=inv*vector(QQ,[2,2,*cmp,*cross]);assert all(c in ZZ for c in v) and v*G*v==-2
            f=O+v;assert f*G*f==0 and f*G*F==3
            ints=VG*f;assert min(ints)>=0
            zero=[i for i,c in enumerate(ints) if c==0];one=[i for i,c in enumerate(ints) if c==1]
            vertical=matrix(QQ,[visible[i] for i in zero]+[v]);rank=vertical.rank()
            rows.append({'source_index':base['source_index'],'subtract_section_word':list(map(int,w)),
                         'bisection_class':list(map(int,v)),'fibre_class':list(map(int,f)),
                         'visible_vertical_indices':zero,'visible_section_indices':one,'visible_MW_upper_bound':int(17-rank)})
    rows.sort(key=lambda r:(-r['visible_MW_upper_bound'],sum(abs(c) for c in r['subtract_section_word']),r['source_index'],r['subtract_section_word']))
    out={'schema':'mestre-bisection-pencil-admission.v1','status':'PASS','source_bisections':bases,'rows':rows,'selected':rows[0] if rows else None,'sources':protocol['sources'],
         'scope':'Exact classes of two rational bisections determined uniquely by old basis intersections, component nonnegativity and square-minus-two. A retained height<=6 section domain supplies translates with O.R2, each an effective nef old-degree3 pencil. Visible section curves certify Jacobian pencils when present. Upper bounds are not ranks. No equation compiler, point search, new parent or complete classification of all degree3 pencils.'}
    with (D/'result.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS exact bisection classes',bases,'pencils',len(rows),'upper histogram',{u:sum(r['visible_MW_upper_bound']==u for r in rows) for u in sorted(set(r['visible_MW_upper_bound'] for r in rows))},flush=True)
    print('SELECTED',out['selected'],flush=True)
if __name__=='__main__':main()
