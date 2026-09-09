#!/usr/bin/env sage-python
"""Exactly two Artin-compatible relative ideals, one exact LLL each <=25s."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,pari,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_seed_artin_v2';ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
SEED=ART/'det1092_seed_half_ideal_v1/virtual-unit.json'
GEN=[ART/('det1092_generic_virtual_units_v1/basis-%02d.json'%i) for i in range(8)]
COLS=[OUT/('column-%02d.json'%i) for i in range(9)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(A):return [[str(x) for x in row] for row in A.rows()]
def coeff(a):return [str(a[i]) for i in range(3)]
save('relative-protocol.json',{'classification':'retrospective exact two-candidate reduction',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,SEED,*GEN,*COLS,Path(__file__)]},
    'limits':{'seconds':25,'relative_ideal_candidates':2,'LLL_vectors_per_class':1,
              'point_searches':0,'unit_groups':0,'class_groups':0,'factorizations':0},
    'selection':'The exact affine F2 solution set to all seven frozen Artin equations; no further ideal search.'})
d=json.loads(ARITH.read_text());seed=json.loads(SEED.read_text());gen=[json.loads(p.read_text()) for p in GEN]
cols=[json.loads(p.read_text()) for p in COLS];M=matrix(GF(2),[c['artin_bits'] for c in cols])
A=M[:8,:];v=M.row(8);word=A.transpose().solve_right(v);kernel=A.left_kernel().basis_matrix()
assert kernel.nrows()==1 and A.rank()==7 and word*A==v
words=[word,word+kernel.row(0)]
R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);K=NumberField(f,'theta')
nf=pari.nfinit([pari(f),d['S_finite']]);zk=[K(R(str(v))) for v in nf.nf_get_zk()]
def pa(a):return pari.Mod(pari(R(a.list())),pari(f))
results=[]
for index,w in enumerate(words):
    J=pari(matrix(QQ,seed['reduced_square_root_ideal']));alpha=K(R(seed['beta']))
    for bit,g in zip(w,gen):
        if bit:J=pari.idealdiv(nf,J,pari(matrix(QQ,g['ideal'])));alpha/=K(R(g['beta']))
    assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,pa(alpha))
    inv=matrix(QQ,pari.idealinv(nf,J).sage())
    basis=[sum((inv[i,j]*zk[i] for i in range(3)),K.zero()) for j in range(3)]
    G=matrix(QQ,3,3,lambda i,j:(basis[i]*basis[j]).trace())
    scale=lcm([x.denominator() for x in G.list()])
    U=matrix(ZZ,pari.qflllgram(pari(matrix(ZZ,scale*G))).sage());assert abs(U.det())==1
    gamma=sum((U[i,0]*basis[i] for i in range(3)),K.zero());beta=alpha*gamma**2
    reduced=pari.idealmul(nf,J,pa(gamma));N=ZZ(pari.idealnorm(nf,reduced))
    assert N>0 and pari.idealpow(nf,reduced,2)==pari.idealhnf(nf,pa(beta))
    record={'index':index,'generic_virtual_basis_word':list(map(int,w)),
            'unreduced_ideal':rows(matrix(QQ,J.sage())),'alpha':coeff(alpha),
            'gamma':coeff(gamma),'beta':coeff(beta),'ideal':rows(matrix(QQ,reduced.sage())),
            'norm':str(N),'unit_found':N==1,
            'boundary':'A nonunit reduction is not a proof of nonprincipality.'}
    save('relative-%02d.json'%index,record);results.append(record)
    print('relative',index,'word',list(map(int,w)),'unit',N==1,'norm digits',len(str(N)),flush=True)
save('relative-summary.json',{'status':'PASS_TWO_ARTIN_COMPATIBLE_IDEAL_REDUCTIONS',
    'classification':'retrospective verified application','generic_kernel':list(map(int,kernel.row(0))),
    'compatible_words':[list(map(int,w)) for w in words],
    'unit_found':any(r['unit_found'] for r in results),
    'criterion':'The seed ideal is in the inherited ideal-class image iff at least one of these two relative ideals is principal. Neither direction follows merely from its reduced norm.',
    'checker_sha256':sha(Path(__file__))})
