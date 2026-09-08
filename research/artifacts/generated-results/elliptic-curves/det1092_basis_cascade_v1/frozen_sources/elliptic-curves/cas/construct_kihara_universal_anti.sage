#!/usr/bin/env sage-python
"""Seek the primitive quadratic anti-section over the full parent-ratio line."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,prod,matrix,vector
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/kihara-universal-anti-v1'
def write(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    V=PolynomialRing(QQ,'v');v=V.gen();F=V.fraction_field();R=PolynomialRing(F,'T');T=R.gen();S=PolynomialRing(R,'x');x=S.gen()
    roots=[F(0),(2*v*v+v+2)**2,2*(v+1)**2*(2*v*v+v+1),4*v*v-v+4,
           v*(2*v-1)*(2*v*v+4*v+5),4*v**4+8*v**3+9*v*v-2*v+2]
    product=prod((x-a-T)*(x-a+T) for a in roots);square=x**6
    for j in range(5,-1,-1):square+=(product[6+j]-(square*square)[6+j])/2*x**j
    remainder=square*square-product;assert remainder.degree()==4
    q=[]
    for c in remainder.list():
        a,b=c.quo_rem(T*T);assert not b;q.append(a)
    e,d,c,b,a=q;I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    A=-27*I;B=-27*J;assert A(-T)==A and B(-T)==B
    aa=R([A[2*i] for i in range(5)]);bb=R([B[2*i] for i in range(7)])
    lam=F(a[2]).sqrt();assert lam in F and lam*lam==a[2]
    aa=aa/lam**4;bb=bb/lam**6;assert aa[4]==-432 and bb[6]==3456
    shift=-aa[3]/(4*aa[4]);aa=R(aa(T+shift));bb=R(bb(T+shift));assert aa[3]==0 and bb[5]==0
    raw={'schema':'kihara-universal-anti-model.v1','parent_parameter':'v=p/q, q=1','roots':list(map(str,roots)),
         'quartic_coefficients':[list(map(str,c.list())) for c in q],'raw_A':list(map(str,A.list())),'raw_B':list(map(str,B.list())),
         'scale':str(lam),'base_shift':str(shift),'normalized_quotient_A':list(map(str,aa.list())),'normalized_quotient_B':list(map(str,bb.list())),'sources':protocol['sources']}
    write(D/'model.json',raw);print('PASS full parent-ratio model; leading scale',lam,flush=True)
    P=PolynomialRing(F,names=('k','l'),order='lex');k,l=P.gens();U=PolynomialRing(P,'s');s=U.gen()
    X=-36*s*s+k*s+l;Y=96*s**3-6*k*s*s
    residual=-3*Y*Y-X**3-U(aa)*X-U(bb)
    Y+=residual[4]/576*s
    residual=-3*Y*Y-X**3-U(aa)*X-U(bb);Y+=residual[3]/576
    residual=-3*Y*Y-X**3-U(aa)*X-U(bb);assert residual.degree()<=2
    equations=list(residual);write(D/'equations.json',{'unknowns':['k','l'],'equations':list(map(str,equations)),'X':str(X),'Y_div_sqrt_minus3':str(Y)})
    ideal=P.ideal(equations);gb=ideal.groebner_basis();write(D/'groebner.json',{'basis':list(map(str,gb))});print('ANTI GROEBNER',len(gb),[g.total_degree() for g in gb],flush=True)
    linear=[g for g in gb if g.total_degree()==1];assert len(linear)>=2
    M=matrix(F,[[g.monomial_coefficient(k),g.monomial_coefficient(l)] for g in linear]);rhs=vector(F,[-g.constant_coefficient() for g in linear]);solution=M.solve_right(rhs)
    xpoly=R([F(c(*solution)) for c in X.list()]);ypoly=R([F(c(*solution)) for c in Y.list()])
    assert -3*ypoly**2==xpoly**3+aa*xpoly+bb and xpoly[2]==-36 and ypoly[3]==96
    rawX=R(lam**2*xpoly(T*T-shift));rawY=R(lam**3*ypoly(T*T-shift))
    assert -3*rawY**2==rawX**3+A*rawX+B
    write(D/'result.json',{'schema':'kihara-universal-anti.v1','status':'EXACT_IDENTITY_PENDING_REPLAY','parent_parameter':'v=p/q, q=1',
          'anti_X':list(map(str,rawX.list())),'anti_Y_div_sqrt_minus3':list(map(str,rawY.list())),
          'normalized_quotient_X':list(map(str,xpoly.list())),'normalized_quotient_Y_div_sqrt_minus3':list(map(str,ypoly.list())),
          'solution':list(map(str,solution)),'scale':str(lam),'base_shift':str(shift),'sources':protocol['sources'],
          'scope':'Exact anti-section identity over Q(v)(sqrt(-3))(T) on the full Kihara parent-ratio line, not only the rank14 path. Requires separate independent reconstruction, nondegeneracy and height checks. No new rational direction or fibre search; exceptional zero/pole parent parameters are not covered implicitly.'})
    print('PASS universal anti-section identity over full parent-ratio line',flush=True)
if __name__=='__main__':main()
