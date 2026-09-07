#!/usr/bin/env sage-python
"""Complete the fixed200 retained pencils with visible rank ceilings12 or13."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-retained-pencil-frames-v1'
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    bundle=json.loads((ART/'mestre_degree_three_proof_bundle_v1.json').read_text());ns=bundle['ns'];h=bundle['heights']
    admission=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-bisection-pencil-admission-v1/result.json').read_text())
    selected=[r for r in admission['rows'] if r['visible_MW_upper_bound'] in (12,13)];assert len(selected)==200
    glue=ns['eligible_index_two_geometric_glues'][0];G=matrix(ZZ,glue['geometric_Gram']);action=matrix(ZZ,glue['Galois_action'])
    change=matrix.identity(QQ,19);change[glue['replace_basis_index']]=vector(QQ,glue['fixed_numerator']+[1])/2;togeo=change.inverse()
    rows=[]
    for i,row in enumerate(selected):
        f=vector(QQ,row['fibre_class']+[0])*togeo
        S=vector(QQ,bundle['visible_curve_classes'][row['visible_section_indices'][0]]+[0])*togeo
        assert f*G*f==0 and S*G*S==-2 and f*G*S==1 and f*action==f and S*action==S
        kernel=matrix(ZZ,[f*G,S*G]).right_kernel_matrix();H=-kernel*G*kernel.transpose();assert H.is_positive_definite() and H.det()==468
        U=matrix(ZZ,pari(H).qflllgram());assert abs(U.det())==1;reduced=U.transpose()*H*U
        enumeration=pari(reduced).qfminim(2,100000,2);short=matrix(ZZ,enumeration[2]).transpose();assert 2*short.nrows()==int(enumeration[0])
        roots=short*U.transpose()*kernel;rg=roots.rank();rq=rg-(roots*(action-matrix.identity(ZZ,19))).rank()
        r={**row,'retained_index':i,'geometric_fibre_class':list(map(int,f)),'geometric_zero_section_class':list(map(int,S)),
           'frame_basis':[list(map(int,r)) for r in kernel.rows()],'frame_gram':[list(map(int,r)) for r in H.rows()],
           'reduction_change':[list(map(int,r)) for r in U.rows()],'reduced_frame_gram':[list(map(int,r)) for r in reduced.rows()],
           'root_vectors_up_to_sign':[list(map(int,r)) for r in short.rows()],
           'geometric_root_count':int(enumeration[0]),'geometric_root_rank':int(rg),'rational_root_rank':int(rq),
           'generic_Q_MW_rank':int(16-rq),'generic_geometric_MW_rank':int(17-rg)}
        write(D/('frame'+str(i)+'.json'),r);rows.append(r)
        if i%10==0:print('FRAME',i,'MW',16-rq,17-rg,flush=True)
    rows.sort(key=lambda r:(-r['generic_Q_MW_rank'],sum(abs(c) for c in r['subtract_section_word']),r['source_index'],r['subtract_section_word']))
    first=rows[0];base=next(r for r in bundle['base_bisections'] if r['index']==first['source_index'])
    R=PolynomialRing(QQ,'T');K=R.fraction_field();Z=PolynomialRing(QQ,'z');L=Z.fraction_field();t=L(base['old_T'])
    A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients']);E=EllipticCurve(L,[A(t),B(t)])
    seed=[E([K(v)(t) for v in coords]) for coords in ns['section_points']]
    point=E([L(base['old_x']),L(base['old_y'])])-sum((a*P for a,P in zip(first['subtract_section_word'],seed)),E(0))
    x,y=point.xy();deck=L(base['deck']);assert x(deck)!=x or y(deck)!=y
    result={'schema':'mestre-retained-pencil-frames.v1','status':'ROOT_ENUMERATION_PENDING_INDEPENDENT_REPLAY','rows':rows,
            'geometric_NS_gram':[list(map(int,r)) for r in G.rows()],'Galois_action':[list(map(int,r)) for r in action.rows()],
            'selected_bisection':{'old_T':str(t),'old_x':str(x),'old_y':str(y),'deck':str(deck),'O_intersection':2},'sources':protocol['sources'],
            'scope':'Exactly200 previously retained pencils whose visible ceilings12or13 exceed the current MW11 baseline. Three prior ceiling14 cases are not rerun. No new divisor roster, parameter scan or point search. Exact generic ranks await independent enumeration replay. Best true generic rank selected before any equation work; ties use L1 word length, source index and lex word.'}
    write(D/'result.json',result)
    print('PASS200 histogram',{k:sum(r['generic_Q_MW_rank']==k for r in rows) for k in sorted(set(r['generic_Q_MW_rank'] for r in rows))},flush=True)
    print('BEST',first['generic_Q_MW_rank'],'source',first['source_index'],'word',first['subtract_section_word'],flush=True)
if __name__=='__main__':main()
