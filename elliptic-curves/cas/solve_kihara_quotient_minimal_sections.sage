#!/usr/bin/env sage-python
"""Scholten's height3/2 polynomial-section equations on one fixed quotient."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/kihara-quotient-minimal-sections-v1'
def write(p,obj):
    with p.open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    data=json.loads((ROOT/'artifacts/local/elliptic-curves/kihara-quotient-characters-v1/result.json').read_text())
    row=next(r for r in data['rows'] if r['path_parameter']=='3/2')
    R=PolynomialRing(QQ,'s');s=R.gen();A=R(row['quotient_A']);B=R(row['quotient_B'])
    node=-3*B[6]/(2*A[4]);scale=(3*node).sqrt()/6
    aa=A/scale**4;bb=B/scale**6;assert aa[4]==-432 and bb[6]==3456
    shift=-aa[3]/(4*aa[4]);aa=R(aa(s+shift));bb=R(bb(s+shift));assert aa[3]==0 and bb[5]==0
    P=PolynomialRing(QQ,names=('e','d','c','b','a'),order='degrevlex');e,d,c,b,a=P.gens()
    U=PolynomialRing(P,'s');ss=U.gen();X=12*ss**2+a*ss+b;Y=c*ss**2+d*ss+e
    polynomial=Y*Y-X**3-U(aa)*X-U(bb);assert polynomial.degree()<=4
    equations=list(polynomial);I=P.ideal(equations)
    write(D/'equations.json',{'schema':'kihara-quotient-minimal-sections-equations.v1','path_parameter':'3/2','coordinate_scale':str(scale),'base_shift':str(shift),'A':list(map(str,aa.list())),'B':list(map(str,bb.list())),'unknown_order':['e','d','c','b','a'],'equations':list(map(str,equations)),'section_x':'12*s^2+a*s+b','section_y':'c*s^2+d*s+e','sources':protocol['sources']})
    print('START five height3/2 section equations',flush=True)
    gb=I.groebner_basis();print('GROEBNER',len(gb),'dimension',I.dimension(),flush=True)
    write(D/'groebner.json',{'order':'degrevlex','basis':list(map(str,gb)),'dimension':int(I.dimension())})
    assert I.dimension()==0
    lex=P.change_ring(order='lex');J=I.transformed_basis('fglm',other_ring=lex)
    univariate=[f for f in J.groebner_basis() if len(f.variables())==1]
    factors=[]
    for f in univariate:
        z=f.variables()[0];S=PolynomialRing(QQ,str(z));g=S([f.monomial_coefficient(z**i) for i in range(f.degree(z)+1)])
        factors.append({'variable':str(z),'polynomial':list(map(str,g.list())),'factors':[{'coefficients':list(map(str,h.list())),'multiplicity':int(k)} for h,k in g.factor()]})
    write(D/'result.json',{'schema':'kihara-quotient-minimal-sections.v1','status':'ELIMINATION_COMPLETE','lex_basis':list(map(str,J.groebner_basis())),'univariate_factors':factors,'sources':protocol['sources'],'scope':'Exact polynomial section ideal on one fixed quotient. Factor degrees or candidate coefficients alone do not identify a section or its field; full substitution and independence remain required.'})
    print('PASS univariate factors',[(r['variable'],[(len(f['coefficients'])-1,f['multiplicity']) for f in r['factors']]) for r in factors],flush=True)
if __name__=='__main__':main()
