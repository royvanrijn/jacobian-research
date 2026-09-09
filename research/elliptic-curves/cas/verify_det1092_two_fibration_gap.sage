#!/usr/bin/env sage-python
"""Independent Picard eigenspace and exact short-shell proof.

No PARI, floating point, producer imports, target parameter, or point search.
Exact rational ellipsoid enumeration is capped at 200000 nodes and 25s outside.
"""
import hashlib,json
from fractions import Fraction as Rat
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,identity_matrix,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_two_fibration_action_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FRAME=ART/'det1092_genus1_picard_image_v1/frame.json';MAPS=ART/'det1092_translated_sections_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
    path=OUT/name
    if path.exists():assert read(path)==d
    else:
        with path.open('x') as s:json.dump(d,s,indent=2,sort_keys=True);s.write('\n')
for name in ['protocol.json','degree-minima-protocol-v2.json','nontrivial-gap-protocol.json']:
    for path,expected in read(OUT/name)['inputs'].items():assert sha(ROOT/path)==expected
G=matrix(QQ,read(PARENT)['generic_height_gram']);ns=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G)
assert G.det()==1092 and all(G[i,i] in 2*ZZ for i in range(17))
frame=read(FRAME);action=read(OUT/'action.json');gap=read(OUT/'nontrivial-gap.json')
D=vector(QQ,frame['D']);origin=vector(QQ,frame['new_zero_section'])
F=vector(QQ,[0,1,*([0]*17)]);O=vector(QQ,[1,0,*([0]*17)])
roots=[vector(QQ,r['NS_coordinates']) for r in frame['vertical_old_sections']]
for row,S in zip(frame['vertical_old_sections'],roots):
    q=vector(QQ,row['word']);assert S==vector(QQ,[1,q*G*q/2,*q])
    assert S*ns*S==-2 and S*ns*D==0
pairing=matrix(QQ,[[a*ns*b for b in roots] for a in roots]);pairs=[]
for i in range(10):
    js=[j for j in range(10) if j!=i and pairing[i,j]!=0]
    assert len(js)==1 and pairing[i,js[0]]==2 and roots[i]+roots[js[0]]==D
    if i<js[0]:pairs.append([i,js[0]])
assert pairs==action['component_pairs']
nonidentity=[]
for i,j in pairs:
    vals=[roots[i]*ns*origin,roots[j]*ns*origin];assert sorted(vals)==[0,1]
    nonidentity.append(roots[[i,j][vals.index(0)]])
trivial=matrix(QQ,[D,origin,*nonidentity]);perp=(trivial*ns).right_kernel().basis_matrix()
assert trivial.nrows()==7 and perp.nrows()==12 and trivial.rank()==7
full=trivial.stack(perp).transpose();assert full.det()!=0
inv=full*block_diagonal_matrix(identity_matrix(QQ,7),-identity_matrix(QQ,12))*full.inverse()
assert inv==matrix(QQ,action['alternate_inversion']) and inv*inv==1 and inv.transpose()*ns*inv==ns
def original_translate(S,q):
    a,b=S[:2];v=vector(QQ,S[2:])
    return vector(QQ,[a,b+v*G*q+a*q*G*q/2,*(v+a*q)])
w=vector(QQ,[0]*14+[1,-1,0]);old_inverse=block_diagonal_matrix(identity_matrix(QQ,2),-identity_matrix(QQ,17))
deck=matrix(QQ,[original_translate(old_inverse*e,w) for e in identity_matrix(QQ,19).columns()]).transpose()
assert deck*deck==1 and deck*F==F and deck*D==D
tau=matrix(QQ,action['alternate_B_translation']);assert tau==deck*inv and tau.transpose()*ns*tau==ns
for k in range(6):
    mp=read(MAPS/('map-%02d.json'%k));rp=read(MAPS/('independent-%02d.json'%k))
    assert rp['status']=='PASS_GENERIC_RANK18_AND_EXACT302_EXCLUSION'
    for path,expected in rp['inputs'].items():assert sha(ROOT/path)==expected
    i=mp['case']['section'];sign=mp['case']['translation'];q=vector(QQ,[int(j==i) for j in range(17)])
    S=original_translate(O,q);image=(tau**sign)*S
    assert list(image)==action['six_degree_regressions'][k]['NS'] and image*ns*F==mp['degree']
C=vector(QQ,action['source']['NS']);assert C*ns*F==13
closed=[[-1,1,[6,-19,-6,-6,0,-12,0,0,6,-6,12,-6,0,12,9,3,0]],
        [1,63,[32,-19,-32,-32,0,-64,0,0,32,-32,64,-32,0,64,48,16,0]]]
for sign,d0,rv in closed:
    K=(tau**(-sign))*F;a=C[0];c=K[0];v=vector(QQ,C[2:]);r=vector(QQ,K[2:])
    assert a*c==247 and C*ns*K==d0 and a*r-c*v==vector(QQ,rv)

# Independently derive an exact LDL decomposition, using only saved unimodular
# columns as an efficiency choice. Every branch bound is integer/rational exact.
U=matrix(ZZ,read(OUT/'degree-minima-v2-plus.json')['LLL_columns']);assert abs(U.det())==1
reduced=U.transpose()*G*U;n=17
A=[[Rat(int(x)) for x in row] for row in reduced.rows()]
L=[[Rat(int(i==j)) for j in range(n)] for i in range(n)];diag=[]
for j in range(n):
    pivot=A[j][j]-sum(L[j][k]**2*diag[k] for k in range(j));assert pivot>0;diag.append(pivot)
    for i in range(j+1,n):L[i][j]=(A[i][j]-sum(L[i][k]*L[j][k]*diag[k] for k in range(j)))/pivot
assert all(A[i][j]==sum(L[i][k]*diag[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
word=[0]*n;found=[];nodes=0
def enumerate_shell(j,used):
    global nodes
    nodes+=1
    if nodes>200000:raise RuntimeError('INDEPENDENT_NODE_LIMIT_NO_CERTIFICATE')
    if j<0:
        v=U*vector(ZZ,word);norm=v*G*v
        assert norm==QQ(used.numerator)/used.denominator
        if norm:assert norm==4;found.append(tuple(map(int,v)))
        return
    shift=sum((L[k][j]*word[k] for k in range(j+1,n)),Rat(0))
    allowance=(Rat(4)-used)/diag[j]
    if allowance<0:return
    b=shift.denominator;a=shift.numerator
    radius=isqrt((allowance*b*b).numerator//(allowance*b*b).denominator)
    lower=-((radius+a)//b);upper=(radius-a)//b
    for x in range(lower,upper+1):
        word[j]=x;enumerate_shell(j-1,used+diag[j]*(Rat(x)+shift)**2)
enumerate_shell(16,Rat(0))
saved=read(OUT/'norm-four-shell.json');assert sorted(found)==[tuple(q) for q in saved['words']]
assert len(found)==len(set(found))==saved['count']==2436
proofs=[]
for branch,(sign,d0,rlist) in zip(gap['branches'],closed):
    r=vector(QQ,rlist);rnorm=r*G*r
    products=[r*G*vector(QQ,q) for q in found];maximum=max(products)
    degree=d0+494-maximum;winners=[list(q) for q,value in zip(found,products) if value==maximum]
    assert degree==branch['minimum_nonzero_translation_degree']
    assert sorted(winners)==sorted(branch['minimizers']) and len(winners)==branch['number_of_minimizers']
    # For n>=6, f(n)=d0+247*n/2-sqrt(rnorm*n) increases.
    assert 247**2*6>rnorm
    test=d0+741-degree;assert test>0 and test*test>6*rnorm
    proofs.append({'sign':sign,'r_norm':int(rnorm),'maximum_short_shell_pairing':int(maximum),
        'minimum_nonzero_translation_degree':int(degree),'multiplicity':len(winners),
        'norm_ge_six_exclusion_squared_margin':int(test*test-6*rnorm)})
paths=[PARENT,FRAME,OUT/'action.json',OUT/'norm-four-shell.json',OUT/'nontrivial-gap.json',Path(__file__)]
save('independent-replay.json',{'status':'PASS_EXACT_GLOBAL_TWO_FIBRATION_DEGREE_GAP',
    'classification':'new verified generic-only geometric obstruction',
    'component_pairs':pairs,'independent_shell_nodes':nodes,'short_vectors':len(found),
    'branches':proofs,'ellipsoid_radius_only_63_replay_unnecessary':True,
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
    'boundary':'Global over all original generic translations, for this single starting curve and final +/-B. No general automorphism-group exclusion, first-seed incidence, or new specialized rank claimed.'})
print('PASS_EXACT_GLOBAL_TWO_FIBRATION_DEGREE_GAP','shell nodes',nodes,'short vectors',len(found),'minima',[p['minimum_nonzero_translation_degree'] for p in proofs],flush=True)
