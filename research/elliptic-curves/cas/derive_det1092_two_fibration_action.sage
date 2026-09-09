#!/usr/bin/env sage-python
"""Generic Picard actions; one frozen 70-word, lattice-only diagnostic.

No specialized parameter, exceptional point, or point search. External 25s cap.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,identity_matrix,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_two_fibration_action_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FRAME=ART/'det1092_genus1_picard_image_v1/frame.json'
MAPS=ART/'det1092_translated_sections_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
    OUT.mkdir(exist_ok=True)
    p=OUT/name
    if p.exists():assert read(p)==d
    else:
        with p.open('x') as s:json.dump(d,s,indent=2,sort_keys=True);s.write('\n')
def ints(v):return list(map(int,v))
sources=[PARENT,FRAME,*[MAPS/('map-%02d.json'%i) for i in range(6)],Path(__file__)]
save('protocol.json',{'classification':'generic-only Picard action and bounded word calculation',
    'rule':'Reconstruct alternate inversion from the five vertical component pairs. Verify all six previous maps. Select the unique smallest-degree previous map, then examine final translation by each sign of B after original translation by zero or either sign of a displayed generic basis section. No outcome-based refill.',
    'limits':{'seconds':25,'lattice_words':70,'new_rational_maps':0,'parameter_evaluations':0,'point_searches':0,'exceptional_inputs':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in sources}})
p=read(PARENT);f=read(FRAME);G=matrix(QQ,p['generic_height_gram'])
NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G)
I=identity_matrix(QQ,19);F=I.column(1);O=I.column(0)
D=vector(QQ,f['D']);Z0=vector(QQ,f['new_zero_section']);w=vector(QQ,[0]*14+[1,-1,0])
def dot(a,b):return a*NS*b
def sec(v):return vector(QQ,[1,(v*G*v)/2,*v])
vertical=[vector(QQ,r['NS_coordinates']) for r in f['vertical_old_sections']]
assert len(vertical)==10 and all(dot(a,a)==-2 and dot(a,D)==0 for a in vertical)
pairs=[]
for i,a in enumerate(vertical):
    for j in range(i+1,len(vertical)):
        b=vertical[j]
        if dot(a,b):
            assert dot(a,b)==2 and a+b==D
            pairs.append([i,j])
assert len(pairs)==5 and sorted(sum(pairs,[]))==list(range(10))
nonidentity=[]
for i,j in pairs:
    vals=[dot(Z0,vertical[i]),dot(Z0,vertical[j])]
    assert sorted(vals)==[0,1]
    nonidentity.append(vertical[[i,j][vals.index(0)]])
T=matrix(QQ,[D,Z0,*nonidentity]).transpose();TT=T.transpose()*NS*T
assert T.rank()==7 and TT.det()!=0
projection=T*TT.inverse()*T.transpose()*NS
invalt=2*projection-I
assert all(a in ZZ for a in invalt.list()) and invalt*invalt==I and invalt.transpose()*NS*invalt==NS
assert invalt*D==D and invalt*Z0==Z0 and all(invalt*a==a for a in vertical)
invold=block_diagonal_matrix(identity_matrix(QQ,2),-identity_matrix(QQ,17))
def translation(v):
    ans=[]
    for e in I.columns():
        a,b=e[:2];x=vector(QQ,e[2:])
        ans.append(vector(QQ,[a,b+x*G*v+a*(v*G*v)/2,*(x+a*v)]))
    A=matrix(QQ,ans).transpose()
    assert A.transpose()*NS*A==NS and A*F==F and A*O==sec(v)
    return A
deck=translation(w)*invold
assert deck*deck==I and deck*D==D and deck*F==F
tau=deck*invalt
assert tau*Z0==sec(-vector(QQ,[int(j==15) for j in range(17)]))
regressions=[]
for k in range(6):
    mp=read(MAPS/('map-%02d.json'%k));i=mp['case']['section'];sign=mp['case']['translation']
    S=sec(vector(QQ,[int(j==i) for j in range(17)]));C=(tau**sign)*S
    degree=dot(C,F);assert degree==mp['degree'] and dot(C,C)==-2
    regressions.append({'case':k,'degree':int(degree),'NS':ints(C)})
best=min(regressions,key=lambda r:r['degree']);assert best['case']==1 and best['degree']==13
C=vector(QQ,best['NS']);rows=[]
words=[(None,0,vector(QQ,17))]+[(i,s,vector(QQ,[s*int(j==i) for j in range(17)])) for i in range(17) for s in [-1,1]]
for sign in [-1,1]:
    targetF=(tau**(-sign))*F
    a,b=C[:2];v=vector(QQ,C[2:]);c,e=targetF[:2];r=vector(QQ,targetF[2:])
    for index,s,q in words:
        out=(tau**sign)*translation(q)*C
        degree=dot(out,F)
        formula=dot(C,targetF)+(c*v-a*r)*G*q+a*c*(q*G*q)/2
        assert degree==formula and degree>=0 and dot(out,out)==-2
        rows.append({'original_section':index,'original_sign':s,'alternate_sign':sign,'degree':int(degree),'NS':ints(out)})
save('action.json',{'status':'EXACT_GENERIC_PICARD_ACTION_PENDING_INDEPENDENT_REPLAY',
    'component_pairs':pairs,'trivial_Gram':[[str(c) for c in row] for row in TT.rows()],
    'alternate_inversion':[[int(c) for c in row] for row in invalt.rows()],
    'alternate_B_translation':[[int(c) for c in row] for row in tau.rows()],
    'six_degree_regressions':regressions,'source':best,'word_ledger':rows,
    'degree_histogram':{str(d):sum(r['degree']==d for r in rows) for d in sorted(set(r['degree'] for r in rows))},
    'boundary':'Lattice action and geometric degree, not rational splitting at any fixed fibre. No new generic curve equations are asserted.'})
print('component pairs',pairs,'trivial determinant',TT.det(),'six degree regressions',[r['degree'] for r in regressions])
print('70-word degrees',sorted(set(r['degree'] for r in rows)))
for r in rows:
    if 1<r['degree']<=15:print({k:v for k,v in r.items() if k!='NS'})
