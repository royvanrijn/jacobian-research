#!/usr/bin/env sage-python
"""Independent HNF witnesses, universal degree identity and affine CVP.

One parity coset only. The producer enumerates parity vectors; this check
enumerates integer translations in a closed affine ellipsoid. <=25s and
100000 nodes. No producer import, point search or exceptional input.
"""
import hashlib,json,signal
from fractions import Fraction as Rat
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,identity_matrix,block_diagonal_matrix
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_conic_alternate_translation_coset_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rat(x):
    x=QQ(x);return Rat(int(x.numerator()),int(x.denominator()))
def save(name,row):
    path=OUT/name
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(row,f,indent=2,sort_keys=True);f.write('\n')
protocol=read(OUT/'protocol.json');data=read(OUT/'frame.json');answer=read(OUT/'minimum.json')
for name,h in protocol['inputs'].items():assert sha(ROOT/name)==h
assert data['protocol_sha256']==answer['protocol_sha256']==sha(OUT/'protocol.json')
assert answer['frame_sha256']==sha(OUT/'frame.json')
save('replay-protocol.json',dict(checker_sha256=sha(Path(__file__)),seconds=25,node_cap=100000,
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'frame.json',OUT/'minimum.json',
        ART/'det1092_two_fibration_action_v1/independent-replay.json',
        ART/'det1092_bifibration_conic_sources_v1/independent-replay.json']}))
for path in [ART/'det1092_two_fibration_action_v1/independent-replay.json']:
    for name,h in read(path)['inputs'].items():assert sha(ROOT/name)==h
parent=read(ART/'curve302_recovered_mw17_parent_v1.json');frame=read(ART/'det1092_genus1_picard_image_v1/frame.json')
G=matrix(QQ,parent['generic_height_gram']);NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G)
I=identity_matrix(QQ,19);F=I.column(1);D=vector(QQ,frame['D']);O=vector(QQ,frame['new_zero_section'])
roots=[vector(QQ,r['NS_coordinates']) for r in frame['vertical_old_sections'] if vector(QQ,r['NS_coordinates'])*NS*O==0]
trivial=matrix(QQ,[D,O,*roots]);perp=matrix(QQ,data['projection'])
assert perp*perp==perp and perp.transpose()*NS==NS*perp and perp.rank()==12
assert perp*trivial.transpose()==0 and (trivial*NS)*perp==0
den=ZZ(data['denominator']);assert den>0
generators=matrix(ZZ,den*perp.transpose());H=matrix(ZZ,data['HNF']);U=matrix(ZZ,data['HNF_unimodular'])
assert abs(U.det())==1 and U*generators==H
assert H[:12].rank()==12 and H[12:]==0
basis=matrix(QQ,data['basis_NS']);assert basis==matrix(QQ,H[:12]).transpose()/den
gram=matrix(QQ,data['gram']);assert gram==-basis.transpose()*NS*basis and gram.det()==QQ(273)/8
scaled=matrix(ZZ,data['scaled_gram']);assert scaled==2*gram
w=vector(QQ,D[2:]);a=vector(QQ,read(ART/'det1092_rational_bisection_index_v1/orbit-47755.json')['word'])
C=vector(QQ,[2,4,*a]);L=perp*(C-F);ell=vector(ZZ,data['L_word'])
assert basis*ell==L==vector(QQ,data['L_NS']) and ell*gram*ell==5
assert [C*NS*r for r in roots]==[2,1,1,1,0] and all(F*NS*r==1 for r in roots)
assert C*NS*F==2 and C*NS*D==F*NS*D==2
# Balanced original fibre: root components do not cause a periodic term.
phiF=perp*F
assert F==2*O+5*D-sum((r/2 for r in roots),vector(QQ,19))+phiF
# Universal identity in twelve indeterminates. The written proof identifies
# this vector with the actual pullback under arbitrary generic translation.
P=PolynomialRing(QQ,12,'q');q=vector(P,P.gens());qphi=basis.change_ring(P)*q
norm=q*gram.change_ring(P)*q;pairF=-(vector(P,phiF)*NS.change_ring(P)*qphi)
pullback=vector(P,F)+(-pairF+norm)*vector(P,D)-2*qphi
assert pullback*NS.change_ring(P)*pullback==0
assert vector(P,D)*NS.change_ring(P)*pullback==2
assert all(vector(P,r)*NS.change_ring(P)*pullback==1 for r in roots)
degree=vector(P,C)*NS.change_ring(P)*pullback
parity=2*q+vector(P,ell)
assert degree==(parity*gram.change_ring(P)*parity-1)/2
# Independent closed affine enumeration q=V*z with radius10/4 in 2*Gram.
V=matrix(ZZ,data['LLL_columns']);assert abs(V.det())==1
reduced=V.transpose()*scaled*V;center=[rat(v/2) for v in V.inverse()*ell];n=12
mu=[[Rat(int(i==j)) for j in range(n)] for i in range(n)];diag=[]
for j in range(n):
    d=rat(reduced[j,j])-sum(mu[j][k]**2*diag[k] for k in range(j));assert d>0;diag.append(d)
    for i in range(j+1,n):mu[i][j]=(rat(reduced[i,j])-sum(mu[i][k]*mu[j][k]*diag[k] for k in range(j)))/d
assert all(rat(reduced[i,j])==sum(mu[i][k]*mu[j][k]*diag[k] for k in range(n)) for i in range(n) for j in range(n))
z=[0]*n;found=[];nodes=0
def visit(j,used):
    global nodes
    nodes+=1;assert nodes<=100000,'UNKNOWN_INDEPENDENT_NODE_CAP'
    if j<0:
        word=V*vector(ZZ,z);r=2*word+ell;nn=r*scaled*r
        assert rat(nn)/4==used and nn<=10
        found.append((int(nn),tuple(map(int,word))));return
    shift=center[j]+sum((mu[k][j]*(Rat(z[k])+center[k]) for k in range(j+1,n)),Rat(0))
    allowance=(Rat(5,2)-used)/diag[j]
    if allowance<0:return
    b,a0=shift.denominator,shift.numerator
    radius=isqrt((allowance*b*b).numerator//(allowance*b*b).denominator)
    for value in range(-((radius+a0)//b),(radius-a0)//b+1):
        z[j]=value;visit(j-1,used+diag[j]*(Rat(value)+shift)**2)
visit(11,Rat(0))
assert len(found)==2 and sorted(found)==[(10,tuple(map(int,-ell))),(10,tuple([0]*12))] or False
expected={(10,tuple(map(int,-ell))),(10,tuple([0]*12))}
assert set(found)==expected
assert answer['scaled_minimum']==10 and answer['minimum_height']=='5' and answer['original_degree']==2
assert answer['multiplicity']==2 and {tuple(r['q_word']) for r in answer['minimizers']}=={q for _,q in found}
# The two curves are exactly the same opposite pair as in the original
# translation minimum: q_original=w-a. This checks the intersection of
# the two minimum-degree orbits without constructing more rational maps.
originalq=w-a
companion=vector(QQ,[2,4+a*G*originalq+originalq*G*originalq,*(a+2*originalq)])
assert companion*NS*companion==-2 and companion*NS*F==companion*NS*D==2
images=[]
for row in answer['minimizers']:
    q=vector(QQ,row['q_word']);phi=basis*q;height=q*gram*q;S=vector(QQ,row['q_section_NS']);bits=row['component_bits']
    assert height==QQ(row['q_height']) and set(bits)<={0,1}
    assert S==O+(height/2+QQ(sum(bits))/4)*D-sum((b*r/2 for b,r in zip(bits,roots)),vector(QQ,19))+phi
    assert all(v in ZZ for v in S) and S*NS*S==-2 and S*NS*D==1
    pp=(trivial*NS).right_kernel().basis_matrix().transpose();bb=trivial.transpose().augment(pp)
    ii=[D,S,*[r if not b else D-r for b,r in zip(bits,roots)],*[v-(v*NS*phi)*D for v in pp.columns()]]
    action=matrix(QQ,ii).transpose()*bb.inverse()
    assert all(v in ZZ for v in action.list()) and action.transpose()*NS*action==NS
    image=action*C;assert image==vector(QQ,row['image_NS']) and image in [C,companion]
    images.append(tuple(image))
assert set(images)=={tuple(C),tuple(companion)}
save('independent-replay.json',dict(status='PASS_FULL_MW12_CVP_AND_TWO_SIDED_CONIC_MINIMUM',
    classification='new generic-only global translation-orbit obstruction',MW12_determinant='273/8',
    L_height='5',universal_degree_identity='degree(t_Q C)= (height(2Q+L)-1)/2',
    minimum_parity_height='5',minimum_original_degree=2,multiplicity=2,
    independent_affine_nodes=nodes,minimizing_Q_words=[list(q) for _,q in sorted(found)],
    minimum_curves=[list(map(str,C)),list(map(str,companion))],
    no_original_section_in_entire_alternate_translation_orbit=True,
    degree_two_plateau_under_both_translation_groups_has_two_curves=True,
    boundary='Any alternating path to a section must leave bidegree<=2. No exclusion of longer/higher-degree paths or any302 seed incidence.',
    replay_protocol_sha256=sha(OUT/'replay-protocol.json')))
print('PASS_FULL_MW12_CVP_AND_TWO_SIDED_CONIC_MINIMUM','nodes',nodes,'two minimum curves, degree2',flush=True)
