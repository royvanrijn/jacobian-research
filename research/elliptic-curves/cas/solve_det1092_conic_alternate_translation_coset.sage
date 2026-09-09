#!/usr/bin/env sage-python
"""One full-MW12 parity coset, not a new source/subgroup-state census.

Construct the integral Mordell-Weil image of the19 NS generators by HNF,
then solve only the norm<=5 parity class forced by the fixed47755 conic.
Each action25s, exact CVP<=100000 nodes, no maps/points/target addresses.
"""
import argparse,hashlib,json,signal,time
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,identity_matrix,block_diagonal_matrix,lcm
from visibility_lattice_v2 import ExactParity
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_conic_alternate_translation_coset_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';FRAME=ART/'det1092_genus1_picard_image_v1/frame.json'
CONIC=ART/'det1092_rational_bisection_index_v1/orbit-47755.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(A):return [list(map(str,r)) for r in A.rows()]
def save(name,row):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(row,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    save('protocol.json',dict(classification='generic-only global minimum in one fixed translation orbit',
        rule='Use the complete integral NS image in the alternate MW12 lattice and the single parity class L=phi(C47755-F). Solve that coset with the original L as feasible radius. No selection from its target incidence.',
        limits=dict(seconds_per_action=25,parity_cosets=1,scaled_radius=10,node_cap=100000,point_searches=0,target_inputs=0,parameter_evaluations=0),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,FRAME,CONIC,Path(__file__),
            Path(__file__).with_name('visibility_lattice_v2.py'),
            ART/'det1092_norm8_singular_members_v1/independent-replay.json',
            ART/'det1092_short_alternate_translations_v1/independent-replay.json']}))
def protocol():
    p=read(OUT/'protocol.json')
    for name,h in p['inputs'].items():assert sha(ROOT/name)==h
    return p
def geometry():
    parent=read(PARENT);frame=read(FRAME);G=matrix(QQ,parent['generic_height_gram'])
    NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G);I=identity_matrix(QQ,19)
    D=vector(QQ,frame['D']);O=vector(QQ,frame['new_zero_section']);F=I.column(1)
    dot=lambda a,b:a*NS*b
    roots=[vector(QQ,r['NS_coordinates']) for r in frame['vertical_old_sections'] if dot(vector(QQ,r['NS_coordinates']),O)==0]
    V=matrix(QQ,[D,O,*roots]).transpose();perp=I-V*(V.transpose()*NS*V).inverse()*V.transpose()*NS
    C=vector(QQ,[2,4]+read(CONIC)['word']);L=perp*(C-F)
    return G,NS,D,O,F,roots,perp,C,L
def frame():
    p=protocol();G,NS,D,O,F,roots,perp,C,L=geometry()
    den=lcm([v.denominator() for v in perp.list()]);integers=matrix(ZZ,den*perp.transpose())
    H,U=integers.hermite_form(transformation=True)
    assert U*integers==H and abs(U.det())==1 and H.rank()==12
    assert all(not any(H.row(i)) for i in range(12,19))
    basis=matrix(QQ,H[:12]).transpose()/den;gram=-basis.transpose()*NS*basis
    assert gram.det()==QQ(273)/8 and all(2*v in ZZ for v in gram.list())
    lword=vector(ZZ,basis.solve_right(L));assert basis*lword==L and lword*gram*lword==5
    scaled=matrix(ZZ,2*gram);lll=scaled.LLL_gram();assert abs(lll.det())==1
    save('frame.json',dict(status='EXACT_FULL_INTEGRAL_MW12_IMAGE_AND_ONE_PARITY',
        projection=rows(perp),denominator=str(den),HNF=rows(H),HNF_unimodular=rows(U),
        basis_NS=rows(basis),gram=rows(gram),scaled_gram=rows(scaled),LLL_columns=rows(lll),
        L_NS=list(map(str,L)),L_word=list(map(int,lword)),L_height='5',
        C_root_intersections=[str(C*NS*r) for r in roots],F_root_intersections=[str(F*NS*r) for r in roots],
        protocol_sha256=sha(OUT/'protocol.json'),
        argument='NS generated integrally by its19 displayed classes; normalized restriction to the elliptic generic fibre is surjective onto MW, with kernel the trivial lattice. Orthogonal projection kills this kernel and has no torsion kernel here. HNF retains the exact image lattice, not its arbitrary rational saturation.'))
    print('full MW12 determinant',gram.det(),'L height',lword*gram*lword,'projection denominator',den,flush=True)
def solve():
    p=protocol();data=read(OUT/'frame.json');G,NS,D,O,F,roots,perp,C,L=geometry()
    basis=matrix(QQ,data['basis_NS']);gram=matrix(QQ,data['gram']);scaled=matrix(ZZ,data['scaled_gram'])
    U=matrix(ZZ,data['LLL_columns']);lword=vector(ZZ,data['L_word']);reduced=U.transpose()*scaled*U
    seed=vector(ZZ,U.inverse()*lword);assert seed*reduced*seed==10
    answer=ExactParity(reduced.rows()).solve(seed,seed,node_limit=p['limits']['node_cap'])
    trivial=matrix(QQ,[D,O,*roots]).transpose();pp=(trivial.transpose()*NS).right_kernel().basis_matrix().transpose()
    allbasis=trivial.augment(pp);results=[]
    for small in answer['minima']:
        r=U*vector(ZZ,small);q=vector(ZZ,(r-lword)/2);phi=basis*q;height=q*gram*q
        possibilities=[]
        for bits in product([0,1],repeat=5):
            S=O+(height/2+QQ(sum(bits))/4)*D-sum((b*v/2 for b,v in zip(bits,roots)),vector(QQ,19))+phi
            if all(v in ZZ for v in S):possibilities.append((S,bits))
        assert len(possibilities)==1
        S,bits=possibilities[0]
        assert S*NS*S==-2 and S*NS*D==1 and [S*NS*v for v in roots]==list(bits)
        if height:assert S*NS*O>=0
        images=[D,S,*[v if not bit else D-v for bit,v in zip(bits,roots)],
            *[v-(v*NS*phi)*D for v in pp.columns()]]
        A=matrix(QQ,images).transpose()*allbasis.inverse()
        assert all(v in ZZ for v in A.list()) and A.transpose()*NS*A==NS
        image=A*C;degree=image*NS*F
        expected=(QQ(answer['norm'])/2-1)/2
        assert degree==expected and degree>=1
        results.append(dict(q_word=list(map(int,q)),q_height=str(height),q_section_NS=list(map(str,S)),component_bits=list(bits),
            image_NS=list(map(str,image)),original_degree=int(degree),parity_minimum_word=list(map(int,r))))
    save('minimum.json',dict(status='EXACT_GLOBAL_ORIGINAL_DEGREE_MINIMUM_IN_FULL_ALTERNATE_TRANSLATION_ORBIT',
        scaled_minimum=answer['norm'],minimum_height=str(QQ(answer['norm'])/2),
        original_degree=results[0]['original_degree'],multiplicity=len(results),nodes=answer['nodes'],minimizers=results,
        frame_sha256=sha(OUT/'frame.json'),protocol_sha256=sha(OUT/'protocol.json'),
        boundary='One fixed conic and all generic alternate translations. This is neither an automorphism-group census nor rational incidence at302.'))
    print('minimum height',QQ(answer['norm'])/2,'original degree',results[0]['original_degree'],'multiplicity',len(results),'nodes',answer['nodes'],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','frame','solve']);args=parser.parse_args()
    signal.alarm(25);start=time.monotonic()
    if args.action=='freeze':freeze()
    elif args.action=='frame':frame()
    else:solve()
    print('seconds',round(time.monotonic()-start,3),flush=True)
