#!/usr/bin/env sage-python
"""Extract the rational component of the saved anti-section equations."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing
ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'artifacts/local/elliptic-curves/kihara-universal-anti-extraction-v1'
OLD=ROOT/'artifacts/local/elliptic-curves/kihara-universal-anti-v1'
def write(p,o):
    with p.open('x') as f:json.dump(o,f,indent=2);f.write('\n')
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for p,h in protocol['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    V=PolynomialRing(QQ,'v');v=V.gen();F=V.fraction_field()
    P=PolynomialRing(F,names=('k','l'),order='lex');k,l=P.gens()
    gb=[P(s) for s in json.loads((OLD/'groebner.json').read_text())['basis']]
    L=PolynomialRing(F,'l');f=L(gb[1](0,L.gen()));factors=list(f.factor())
    linear=[g for g,e in factors if g.degree()==1];assert len(linear)==1
    lv=-linear[0][0]/linear[0][1];kv=-F(gb[0](0,lv));assert all(g(kv,lv)==0 for g in gb)
    U=PolynomialRing(P,'s');eq=json.loads((OLD/'equations.json').read_text())
    assert all(P(s)(kv,lv)==0 for s in eq['equations'])
    R=PolynomialRing(F,'T');T=R.gen()
    X=R([F(c(kv,lv)) for c in U(eq['X']).list()]);Y=R([F(c(kv,lv)) for c in U(eq['Y_div_sqrt_minus3']).list()])
    model=json.loads((OLD/'model.json').read_text());a=R(list(map(F,model['normalized_quotient_A'])));b=R(list(map(F,model['normalized_quotient_B'])))
    assert -3*Y**2==X**3+a*X+b
    lam=F(model['scale']);shift=F(model['base_shift'])
    rawX=R(lam**2*X(T*T-shift));rawY=R(lam**3*Y(T*T-shift))
    A=R(list(map(F,model['raw_A'])));B=R(list(map(F,model['raw_B'])))
    assert -3*rawY**2==rawX**3+A*rawX+B
    quad=next(g for g,e in factors if g.degree()==2)
    result={'schema':'kihara-universal-anti-extraction.v1','status':'EXACT_IDENTITY_PENDING_REPLAY','sources':protocol['sources'],
        'anti_X':list(map(str,rawX.list())),'anti_Y_div_sqrt_minus3':list(map(str,rawY.list())),
        'normalized_quotient_X':list(map(str,X.list())),'normalized_quotient_Y_div_sqrt_minus3':list(map(str,Y.list())),
        'k':str(kv),'l':str(lv),'factor_degrees':[int(g.degree()) for g,e in factors],
        'other_component_discriminant_factorization':str(quad.discriminant().factor()),
        'scope':'One explicit section over Q(v)(sqrt(-3))(T) on the full parent-ratio line. Neither a new rational direction nor a claim of full geometric MW or constant field at every specialization. Exceptional zero/pole parameters require separate checking.'}
    write(D/'result.json',result)
    print('PASS universal quadratic section; other component discriminant',result['other_component_discriminant_factorization'],flush=True)
if __name__=='__main__':main()
