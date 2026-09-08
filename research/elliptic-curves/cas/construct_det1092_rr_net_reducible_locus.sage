#!/usr/bin/env sage-python
"""One-coset exact enumeration classifying reducible members of the fixed net.

Limit200000 exact LDL nodes; one coset only; no rational point search.
"""
import hashlib,json,runpy,sys
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
from visibility_lattice_v2 import ExactParity
OUT=ART/'det1092_rr_net_reducible_locus_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';loader=CAS/'load_curve302_recovered_parent.sage'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text());G=matrix(ZZ,parent['generic_height_gram']);w=-vector(ZZ,net['trace_word'])
    U=G.LLL_gram();assert abs(U.det())==1
    Gr=U.transpose()*G*U;wr=vector(ZZ,U.inverse()*w)
    result=ExactParity(Gr.rows()).solve(wr,wr,node_limit=200000)
    assert result['norm']==10
    minima=sorted({tuple(U*vector(ZZ,v)) for v in result['minima']})
    pairs=[];seen=set()
    for v in minima:
        x=vector(ZZ,[(a+b)//2 for a,b in zip(w,v)]);y=w-x
        key=tuple(sorted([tuple(x),tuple(y)]))
        if key in seen:continue
        seen.add(key);assert x*G*x+y*G*y==10
        pairs.append({'words':[list(map(int,z)) for z in key],
                      'heights':[int(vector(ZZ,z)*G*vector(ZZ,z)) for z in key]})
    assert len(pairs)*2==len(minima)
    # Explicit canonical pair O+C, lifted as a vertical chord through -C.
    E,basis,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    T=-sum((n*P for n,P in zip(w,basis)),E(0));R=E.base_ring().ring();t=R.gen()
    A=list(map(R,net['A']));B=list(map(R,net['B']))
    V=[-R(T[0].numerator()),R(T[0].denominator()),R(0)]
    def flat(row):return vector(QQ,[f[j] for f,bound in zip(row,[10,6,4]) for j in range(bound+1)])
    coeff=matrix(QQ,[flat(A),flat([t*f for f in A]),flat(B)]).transpose().solve_right(flat(V))
    assert all(V[i]==coeff[0]*A[i]+coeff[1]*t*A[i]+coeff[2]*B[i] for i in range(3))
    assert coeff[2] and V[0]+V[1]*T[0]==0
    return {'classification':'verified application and new deduction','status':'PASS_COMPLETE_FIXED_NET_REDUCIBLE_LOCUS',
        'centre_word':list(map(int,w)),'LLL_columns':[list(map(int,row)) for row in U.rows()],
        'exact_enumeration':{'minimum':10,'nodes':result['nodes'],'node_limit':200000,'cosets':1,'short_representatives':[list(map(int,v)) for v in minima]},
        'section_pairs':pairs,'section_pair_count':len(pairs),
        'canonical_O_plus_C_member':{'line_coefficients':[list(map(str,f.list())) for f in V],
            'coefficients_in_A_tA_B':list(map(str,coeff))},
        'vertical_locus':'Projective line span(A,tA): C_min plus one fibre, including infinity.',
        'classification_argument':'A horizontal split P_x+P_(w-x) with k vertical fibres has norm(2x-w)=10-4k. The exact coset minimum10 forces k=0; minimum representatives enumerate all pairs. An irreducible bisection with k vertical fibres has arithmetic genus2-2k, forcing k<=1; k=1 is the unique C_min.',
        'conclusion':'Every reducible member is C_min+F_tau or one of the listed pairs of known generic sections. None supplies an independent specialized point on a smooth fibre.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,np,loader,CAS/'visibility_lattice_v2.py',Path(__file__)]},
        'limits':{'parity_cosets':1,'exact_nodes':200000,'point_searches':0,'parameter_sweeps':0}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'pairs',d['section_pair_count'],'nodes',d['exact_enumeration']['nodes'],flush=True)
