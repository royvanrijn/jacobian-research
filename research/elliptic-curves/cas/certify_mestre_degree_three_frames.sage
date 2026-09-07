#!/usr/bin/env sage-python
"""Complete fibre roots for the three admitted old-degree3 pencils."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-degree-three-frames-v1'
def main():
    p=json.loads((D/'protocol.json').read_text())
    for name,h in p['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
    h=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    admission=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-bisection-pencil-admission-v1/result.json').read_text())
    oldpool=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-degree-three-pencils-v1/roster.json').read_text())
    glue=ns['eligible_index_two_geometric_glues'][0];G=matrix(ZZ,glue['geometric_Gram']);action=matrix(ZZ,glue['Galois_action'])
    change=matrix.identity(QQ,19);change[glue['replace_basis_index']]=vector(QQ,glue['fixed_numerator']+[1])/2
    togeo=change.inverse();assert abs(G.det())==468
    rows=[]
    selected=[r for r in admission['rows'] if r['visible_MW_upper_bound']==14];assert len(selected)==3
    for row in selected:
        f=vector(QQ,row['fibre_class']+[0])*togeo
        S=vector(QQ,oldpool['visible_curve_classes'][row['visible_section_indices'][0]]+[0])*togeo
        assert f*G*f==0 and S*G*S==-2 and f*G*S==1 and f*action==f and S*action==S
        kernel=matrix(ZZ,[f*G,S*G]).right_kernel_matrix();H=-kernel*G*kernel.transpose();assert H.is_positive_definite() and H.det()==468
        U=matrix(ZZ,pari(H).qflllgram());assert abs(U.det())==1;reduced=U.transpose()*H*U
        enumeration=pari(reduced).qfminim(2,100000,2);short=matrix(ZZ,enumeration[2]).transpose();assert 2*short.nrows()==int(enumeration[0])
        roots=short*U.transpose()*kernel;root_rank=roots.rank();fixed_rank=root_rank-(roots*(action-matrix.identity(ZZ,19))).rank()
        mwQ=16-fixed_rank;mwbar=17-root_rank
        r={**row,'geometric_fibre_class':list(map(int,f)),'geometric_zero_section_class':list(map(int,S)),
           'frame_basis':[list(map(int,r)) for r in kernel.rows()],'frame_gram':[list(map(int,r)) for r in H.rows()],
           'reduction_change':[list(map(int,r)) for r in U.rows()],'reduced_frame_gram':[list(map(int,r)) for r in reduced.rows()],
           'root_vectors_up_to_sign':[list(map(int,r)) for r in short.rows()],
           'geometric_root_count':int(enumeration[0]),'geometric_root_rank':int(root_rank),'rational_root_rank':int(fixed_rank),
           'generic_Q_MW_rank':int(mwQ),'generic_geometric_MW_rank':int(mwbar)}
        rows.append(r);print('FRAME',len(rows),'roots',int(enumeration[0]),'root ranks',fixed_rank,root_rank,'MW',mwQ,mwbar,flush=True)
    # Verify the selected actual curve, including intersections, beyond its class.
    first=rows[0];base=next(r for r in json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-rational-bisections-v1/result.json').read_text())['rows'] if r['index']==first['source_index'])
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();Z=PolynomialRing(QQ,'z');z=Z.gen();L=Z.fraction_field()
    t=L(base['old_T']);A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients']);E=EllipticCurve(L,[A(t),B(t)])
    seed=[E([K(v)(t) for v in coords]) for coords in ns['section_points']]
    point=E([L(base['old_x']),L(base['old_y'])])-sum((a*P for a,P in zip(first['subtract_section_word'],seed)),E(0))
    x,y=point.xy();deck=L(base['deck']);assert x(deck)!=x or y(deck)!=y
    def oi(x):
        dx=x.denominator();dt=t.denominator()
        n=dx.degree()-dx.gcd(dt**4).degree()+max(0,x.numerator().degree()-dx.degree()-4*max(0,t.numerator().degree()-dt.degree()))
        assert n%2==0;return n//2
    assert oi(x)==2
    oldG=matrix(ZZ,ns['rational_NS_Gram']);v=vector(ZZ,first['bisection_class'])
    for i,P in enumerate(seed):assert oi((point-P)[0])==(v*oldG)[7+i]
    result={'schema':'mestre-degree-three-frames.v1','status':'ROOT_ENUMERATION_PENDING_INDEPENDENT_REPLAY','rows':rows,
        'geometric_NS_gram':[list(map(int,r)) for r in G.rows()],'Galois_action':[list(map(int,r)) for r in action.rows()],
        'selected_bisection':{'old_T':str(t),'old_x':str(x),'old_y':str(y),'deck':str(deck),'O_intersection':2},
        'sources':p['sources'],'scope':'Complete roots in the frame of three effective nef rational old-degree3 pencils, with rational sections. Generic ranks depend on independent enumeration completeness replay. The first translated bisection is verified as an exact rational curve and against all eleven old section intersections. Explicit new Weierstrass equations and section bases are not constructed; no new point search, parent or fibre rank claim.'}
    with (D/'result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS selected actual bisection O.R2 and all eleven basis intersections',flush=True)
if __name__=='__main__':main()
