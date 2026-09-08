#!/usr/bin/env sage-python
"""Recover a quadratic section via the saved finite quotient algebra."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,QuadraticField,lcm
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/kihara-quotient-section-extraction-v1'
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def main():
    p=json.loads((D/'protocol.json').read_text())
    for name,h in p['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    old=ROOT/'artifacts/local/elliptic-curves/kihara-quotient-minimal-sections-v1'
    eq=json.loads((old/'equations.json').read_text());saved=json.loads((old/'groebner.json').read_text())
    R=PolynomialRing(QQ,names=eq['unknown_order'],order='degrevlex');generators=R.gens();gb=[R(f) for f in saved['basis']]
    leading=[tuple(f.lm().exponents()[0]) for f in gb]
    def standard(v):return not any(all(a>=b for a,b in zip(v,w)) for w in leading)
    exponents={(0,)*5};pending=[(0,)*5]
    while pending:
        v=pending.pop()
        for i in range(5):
            w=tuple(v[j]+int(i==j) for j in range(5))
            if w not in exponents and standard(w):exponents.add(w);pending.append(w)
        assert len(exponents)<=256
    exponents=sorted(exponents,key=lambda v:(sum(v),v));indices={v:i for i,v in enumerate(exponents)}
    monomials=[R.monomial(*v) for v in exponents];n=len(monomials);print('QUOTIENT DIMENSION',n,flush=True)
    assert n==56
    theta=generators[2]+generators[4] # c+a separates opposite ordinate signs.
    M=matrix(QQ,n,n)
    for j,monomial in enumerate(monomials):
        reduced=(theta*monomial).reduce(gb)
        for exponent,coefficient in reduced.dict().items():M[indices[tuple(exponent)],j]=coefficient
    den=lcm(c.denominator() for c in M.list());Z=matrix(ZZ,den*M)
    print('CHARPOLY START max entry bits',max(abs(c).nbits() for c in Z.list()),flush=True)
    cp=Z.charpoly();S=PolynomialRing(QQ,'v');v=S.gen();poly=S(cp(den*v)/den**n)
    fac=poly.factor();write(D/'characteristic.json',{'dimension':n,'theta':'c+a','polynomial':list(map(str,poly.list())),'factors':[{'coefficients':list(map(str,f.list())),'multiplicity':int(m)} for f,m in fac]})
    print('FACTORS',[(f.degree(),m) for f,m in fac],flush=True)
    chosen=next(f for f,m in fac if f.degree()==2 and m==1 and QQ(-f.discriminant()/3).is_square())
    K=QuadraticField(-3,'r');r=K.gen();alpha=(-K(chosen[1])+QQ(-chosen.discriminant()/3).sqrt()*r)/(2*chosen[2])
    assert chosen(alpha)==0
    kernel=(matrix(K,M)-alpha*matrix.identity(K,n)).left_kernel();assert kernel.dimension()==1
    functional=kernel.basis()[0];functional/=functional[0]
    solution=[functional[indices[tuple(int(i==j) for i in range(5))]] for j in range(5)]
    assert all(R(f)(*solution)==0 for f in eq['equations'])
    T=PolynomialRing(K,'s');s=T.gen();ee,dd,cc,bb,aa=solution
    X=12*s*s+aa*s+bb;Y=cc*s*s+dd*s+ee;A=T(eq['A']);B=T(eq['B']);assert Y*Y==X**3+A*X+B
    assert any(c not in QQ for c in solution)
    write(D/'result.json',{'schema':'kihara-quotient-extra-section.v1','status':'EXACT_SECTION_FOUND','field_polynomial':'r^2+3','path_parameter':'3/2','A':list(map(str,A.list())),'B':list(map(str,B.list())),'X':list(map(str,X.list())),'Y':list(map(str,Y.list())),'coordinate_scale':eq['coordinate_scale'],'base_shift':eq['base_shift'],'solution_unknown_order':eq['unknown_order'],'solution':list(map(str,solution)),'sources':p['sources'],'scope':'One exactly substituted non-rational height3/2 section on the normalized first Kihara quotient over Q(sqrt(-3))(s). Independence, exact field interpretation and generic parent consequences require separate replay; no new rational direction or point search.'})
    print('PASS exact non-rational section over Q(sqrt(-3))',flush=True)
if __name__=='__main__':main()
