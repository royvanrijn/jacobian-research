#!/usr/bin/env sage-python
"""Eight fixed generic parity-kernel basis classes, one exact LLL each.

No enumeration of their256 combinations, unit group, class group, or point
search. An ideal reducing to O proves a unit relation; other output does
not prove nonprincipality. Run under timeout25s.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,NumberField,pari,matrix,lcm,GF
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_generic_virtual_units_v1'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
PROJ=ART/'det1092_seed_half_ideal_v1/generic-projection.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,obj):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==obj
    else:
        with path.open('x') as stream:json.dump(obj,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(A):return [[str(x) for x in row] for row in A.rows()]
def coeff(a):return [str(a[i]) for i in range(3)]
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'generic-only bounded reduction, not full principal-ideal testing',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,PROJ,Path(__file__)]},
    'limits':{'seconds':25,'classes':8,'LLL_dimension':3,'vectors_per_class':1,
              'combination_enumeration':0,'class_groups':0,'unit_groups':0,'point_searches':0},
    'selection':'Canonical F2 kernel basis from only the17 generic classes; one unweighted trace-LLL of each inverse square-root ideal.'})
d=json.loads(ARITH.read_text());pr=json.loads(PROJ.read_text())
R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);K=NumberField(f,'theta')
nf=pari.nfinit([pari(f),d['S_finite']]);assert ZZ(nf.disc())==ZZ(d['field_discriminant'])
zk=[K(R(str(a))) for a in nf.nf_get_zk()]
def pa(a):return pari.Mod(pari(R(a.list())),pari(f))
polys=[K(R(x['beta_ascending'])) for x in d['generic_classes']]
places=[(p,i,P) for p in d['S_finite'] for i,P in enumerate(pari.idealprimedec(nf,p))]
valuations=[[int(pari.idealval(nf,pa(a),P))%2 for p,i,P in places] for a in polys]
words=matrix(GF(2),valuations).left_kernel().basis_matrix()
assert rows(words)==rows(matrix(GF(2),pr['generic_kernel_basis'])) and words.nrows()==8
halves=[pari.idealhnf(nf,pa(a),pari(ZZ(a.norm()).sqrt())) for a in polys]
cases=[]
for index,word in enumerate(words.rows()):
    alpha=K.one();I=pari.idealhnf(nf,1)
    for bit,a,half in zip(word,polys,halves):
        if bit:alpha*=a;I=pari.idealmul(nf,I,half)
    J=I
    for p,i,P in places:
        v=int(pari.idealval(nf,pa(alpha),P));w=int(pari.idealval(nf,I,P));assert v%2==0
        if v//2!=w:J=pari.idealmul(nf,J,pari.idealpow(nf,P,v//2-w))
    assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,pa(alpha))
    inv=matrix(QQ,pari.idealinv(nf,J).sage())
    basis=[sum((inv[i,j]*zk[i] for i in range(3)),K.zero()) for j in range(3)]
    G=matrix(QQ,3,3,lambda i,j:(basis[i]*basis[j]).trace())
    scale=lcm([q.denominator() for q in G.list()])
    U=matrix(ZZ,pari.qflllgram(pari(matrix(ZZ,scale*G))).sage());assert abs(U.det())==1
    gamma=sum((U[i,0]*basis[i] for i in range(3)),K.zero());beta=alpha*gamma**2
    Jred=pari.idealmul(nf,J,pa(gamma));assert pari.idealpow(nf,Jred,2)==pari.idealhnf(nf,pa(beta))
    norm=ZZ(pari.idealnorm(nf,Jred));assert norm>0
    row={'index':index,'word':list(map(int,word)),'alpha':coeff(alpha),'gamma':coeff(gamma),
         'beta':coeff(beta),'ideal':rows(matrix(QQ,Jred.sage())),
         'norm':str(norm),'unit_found':norm==1}
    save('basis-%02d.json'%index,row);cases.append(row)
    print('basis',index,'unit',norm==1,'norm digits',len(str(norm)),flush=True)
save('summary.json',{'status':'PASS_EIGHT_FIXED_GENERIC_VIRTUAL_UNIT_REDUCTIONS',
    'classification':'verified application','checker_sha256':sha(Path(__file__)),
    'unit_indices':[r['index'] for r in cases if r['unit_found']],
    'equal_reduced_ideal_pairs':[[i,j] for i in range(8) for j in range(i) if cases[i]['ideal']==cases[j]['ideal']],
    'boundary':'No nonprincipal-ideal assertion follows from a nonunit reduced ideal; no256-class or class-group computation was attempted.'})
