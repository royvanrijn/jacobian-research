#!/usr/bin/env sage-python
"""Calibrated exact quadratic-shell recognition and K3 height-two gate.

Replay finite numerical-height selections with exact point integrality,
linear kernels and every integer form in two declared adaptive boxes.
The height-two obstruction uses the separately proved2-saturation of D.
One worker,120 seconds; no elliptic parent is inferred from a height form.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import runpy
import signal
from sage.all import QQ,ZZ,GF,RealField,PolynomialRing,EllipticCurve,matrix,vector,pari

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results'
INPUT=ART/'elliptic-curves/curve302_integral_shell_inputs_v1.json'
SOURCE=ART/'elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
HEIGHTS=ART/'elliptic-curves/record_height_lattices_28_29_273_302_v1.json'
PARITY=ART/'elliptic-curves/curve302_displayed_group_parity_domain_v1.json'
SPACES=ART/'elliptic-curves/curve302_mw14_16_packet_protocol_v1.json'
CORE=ART/'elliptic-curves/record_rank17_core_candidates_v1.json'
PUBLIC=ROOT/'elliptic-curves/cas/icarm_curve302.py'
OUT=ART/'elliptic-curves/curve302_integral_shell_lift_v1.json'

def digest(p):return sha256(p.read_bytes()).hexdigest()
def rows(m):return [list(map(str,r)) for r in m.rows()]
def normalize(v):return tuple(v) if next(c for c in v if c)>0 else tuple(-c for c in v)
def pairs(n):return [(i,j) for i in range(n) for j in range(i,n)]
def gram(v,n):
    G=matrix(QQ,n)
    for a,(i,j) in zip(v,pairs(n)):G[i,j]=G[j,i]=a
    return G

def select(E,basis,H):
    minimum=RealField(192)(pari(H).qfminim(None,1,2)[1])
    data=pari(H).qfminim(QQ(7)/5*minimum,50000,2)
    V=matrix(ZZ,data[2]).transpose();assert int(data[0])==2*V.nrows()
    chosen=sorted(V.rows(),key=lambda v:(v*H*v,tuple(v)))[:4000]
    integral=[]
    for v in chosen:
        P=sum((a*Q for a,Q in zip(v,basis) if a),E(0))
        if P[0].denominator()==1:integral.append(list(map(int,v)))
    return len(chosen),integral

def kernel(words,n):
    M=matrix(QQ,[[v[i]*v[j]*(1 if i==j else 2) for i,j in pairs(n)]+[-1] for v in words])
    return M.right_kernel().basis_matrix()

def control(saved):
    src=json.loads(SOURCE.read_text());t0=ZZ(10)**40+7;R=PolynomialRing(QQ,'u')
    A,B=[R(src['weierstrass_model'][k+'_coefficients_low_to_high'])(t0) for k in ['A','B']]
    E0=EllipticCurve(QQ,[A,B]);E=E0.integral_model();iso=E0.isomorphism_to(E)
    def value(r):return R(r['numerator_coefficients_low_to_high'])(t0)/R(r['denominator_coefficients_low_to_high'])(t0)
    points=[iso(E0(value(r['X']),value(r['Y']))) for r in src['sections']['records']]
    H=E.height_pairing_matrix(points,precision=192);U=matrix(ZZ,pari(H).qflllgram())
    basis=[sum((a*P for a,P in zip(v,points) if a),E(0)) for v in U.columns()]
    count,words=select(E,basis,U.transpose()*H*U)
    assert count==saved['selected_count']==1313
    assert {normalize(v) for v in words}=={normalize(v) for v in saved['integral_words']}
    K=kernel(words,17);assert K.nrows()==1 and K[0,-1]
    G=gram(4*K[0]/K[0,-1],17)
    # The generic height matrix is used only after selection and recovery.
    truth=U.transpose()*matrix(QQ,src['sections']['height_gram'])*U
    assert G==truth and G==matrix(QQ,saved['recovered_gram']) and G.det()==948
    return dict(selected_rays=count,integral_rays=len(words),kernel_dimension=1,
                recovered_gram=rows(G),matches_withheld_generic_gram=True)

def geometric_gate():
    # All root types of rank<=4, with maximal local diagonal corrections.
    types=[('A1',1,QQ(1)/2),('A2',2,QQ(2)/3),('A3',3,QQ(1)),('A4',4,QQ(6)/5),('D4',4,QQ(1))]
    configs=[]
    for multiplicities in product(range(5),repeat=5):
        rank=sum(m*r for m,(_,r,c) in zip(multiplicities,types))
        if rank>4:continue
        correction=sum(m*c for m,(_,r,c) in zip(multiplicities,types))
        configs.append(dict(types=[name for m,(name,r,c) in zip(multiplicities,types) for _ in range(m)],rank=rank,correction=str(correction)))
    extreme=[r for r in configs if QQ(r['correction'])>=2]
    assert len(extreme)==1 and extreme[0]['types']==['A1']*4
    assert all(QQ(r['correction'])<2 for r in configs if r['rank']<=3)
    # In the forced4A1 case, an even height meets0 or4 nonidentity components.
    profiles=[list(v) for v in product(range(2),repeat=4) if sum(v)%4==0]
    assert profiles==[[0,0,0,0],[1,1,1,1]]
    # Three universal radical vectors supported on the four root rows.
    null=matrix(GF(2),3,20)
    for i in range(3):null[i,2+i]=1;null[i,5]=1
    assert null.rank()==3
    # Linearity in each section component bit proves the identity for all bits.
    for label in range(-1,14):
        N=matrix(GF(2),20);N[0,1]=N[1,0]=1
        for j in range(14):
            N[0,6+j]=N[6+j,0]=1
            for i in range(4):N[2+i,6+j]=N[6+j,2+i]=int(j==label)
        assert N*null.transpose()==0 and N.rank()%2==0
    return dict(root_correction_profiles=configs,even_height_component_profiles=profiles,
                universal_radical_vectors=rows(null),radical_rank=3,
                even_Gram_forced_2_length_at_least=4,primitive_rank20_K3_NS_2_length_at_most=2,
                conclusion='Rank>=15 forbids height2 outright. Rank14 plus an even section lattice2-saturated in geometric MW forbids height2 by the discriminant-length contradiction. The302 primitive image and2-saturation of D enforce that saturation for any proposed lift.')

def build():
    data=json.loads(INPUT.read_text());calibration=control(data['control'])
    pub=runpy.run_path(str(PUBLIC));E=EllipticCurve(QQ,list(map(QQ,pub['GENERAL_WEIERSTRASS_COEFFICIENTS'])))
    public=[E(QQ(x),QQ(y)) for x,y in pub['POINTS']]
    h=next(r for r in json.loads(HEIGHTS.read_text())['curves'] if r['label']=='curve302')
    ambient=matrix(RealField(192),h['height_gram'])
    parity=json.loads(PARITY.read_text());assert parity['status']=='PASS_DISPLAYED_GROUP_2_SATURATED'
    assert matrix(GF(2),parity['binary_matrix_rows']).rank()==31 and parity['rational_torsion_order']==1
    expected={r['candidate_id']:matrix(ZZ,r['basis_rows']) for r in json.loads(SPACES.read_text())['candidate_inputs']}
    c=next(r for r in json.loads(CORE.read_text())['curves'] if r['label']=='curve302')
    expected['core17']=matrix(ZZ,c['saturated_basis_columns_in_public_point_coordinates']).transpose()
    targets=[];kernels={}
    for saved in [dict(data['target_core'],candidate_id='core17',rank=17)]+data['target_spaces']:
        label=saved['candidate_id'];n=saved['rank'];B=matrix(ZZ,saved['public_embedding'])
        if label=='core17':B=B.transpose()
        assert B.row_module()==expected[label].row_module()
        assert list(B.smith_form()[0].diagonal())==[1]*n
        basis=[sum((a*P for a,P in zip(v,public) if a),E(0)) for v in B.rows()]
        count,words=select(E,basis,B*ambient*B.transpose())
        assert count==saved['selected_count']
        assert {normalize(v) for v in words}=={normalize(v) for v in saved['integral_words']}
        K=kernel(words,n)
        assert K==matrix(QQ,len(saved['kernel_basis']),n*(n+1)//2+1,sum(saved['kernel_basis'],[]))
        kernels[label]=K
        targets.append(dict(candidate_id=label,rank=n,selected_rays=count,integral_rays=len(words),
                            point_word_rank=matrix(ZZ,words).rank(),kernel_dimension=K.nrows(),
                            admits_nonzero_constant=any(K.column(-1)),primitive_public_embedding=True))
        print('EXACT_SHELL',label,len(words),K.nrows(),bool(any(K.column(-1))),flush=True)
    boxes=[]
    for box in data['integer_boxes']:
        label=box['candidate_id'];saved=next(r for r in data['target_spaces'] if r['candidate_id']==label)
        n=saved['rank'];K=kernels[label];assert K[0,-1]==1 and all(v==0 for v in K.column(-1)[1:])
        positive=[];total=0
        for coefficients in product(*box['box']):
            total+=1;v=4*K[0]+sum((a*b for a,b in zip(coefficients,K.rows()[1:])),K[0]*0)
            G=gram(v,n)
            if not G.is_positive_definite():continue
            assert all(a in ZZ for a in G.list()) and all(a%2==0 for a in G.diagonal())
            G=matrix(ZZ,G)
            assert all(vector(ZZ,w)*G*vector(ZZ,w)==4 for w in saved['integral_words'])
            found=pari(G).qfminim(2,1,2);assert int(found[0])>0
            w=vector(ZZ,matrix(ZZ,found[2]).column(0));assert w*G*w==2
            positive.append(dict(coefficients=list(coefficients),gram=rows(G),determinant=str(G.det()),height_two_word=list(map(int,w))))
        expected_count=142 if label==32 else 2
        assert len(positive)==expected_count
        boxes.append(dict(candidate_id=label,form_count=total,positive_even_forms=positive,
                          excluded_as_K3_height_lifts=len(positive)))
        print('INTEGER_BOX',label,total,len(positive),'ALL_K3_LIFTS_EXCLUDED',flush=True)
    return dict(schema='curve302.integral-shell-lift.v1',status='CONTROL_GRAM_RECOVERED_144_FINITE_K3_HEIGHT_LIFTS_EXCLUDED',
                input_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),INPUT,SOURCE,HEIGHTS,PARITY,SPACES,CORE,PUBLIC]},
                control=calibration,targets=targets,integer_boxes=boxes,height_two_gate=geometric_gate(),
                boundary='One native-family calibration and six bounded target point selections. The integer boxes are adaptive and finite; other forms, omitted points, different subspaces and non-K3 elliptic surfaces remain open. Integrality does not imply generic zero-section disjointness. No302 parent or generic basis is recovered.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args()
    signal.alarm(120);result=build()
    # Explicit conversion handles Sage integers in compact certificate fields.
    serialized=json.dumps(result,sort_keys=True,indent=2,default=lambda x:int(x))+'\n'
    if args.build:
        assert not OUT.exists();OUT.write_text(serialized)
    assert json.loads(serialized)==json.loads(OUT.read_text())
    print(result['status'],flush=True)
