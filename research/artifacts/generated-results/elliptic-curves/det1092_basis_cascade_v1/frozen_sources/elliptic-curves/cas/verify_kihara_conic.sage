#!/usr/bin/env sage-python
"""Independent conic rank and saturation proof by deck symmetry.
Uses the previously independently proved full rank12 parent as an explicit
input theorem. No constructor imports or repeated full height calculation.
"""
import json,argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,lcm,gcd
def finite_rank(model,points,primes,ell):
    E=EllipticCurve(QQ,model);P=[E([QQ(x),QQ(y)]) for x,y in points];projective=[]
    for point in P:
        den=lcm([c.denominator() for c in point]);v=[ZZ(c*den) for c in point];g=gcd(v);projective.append([c//g for c in v])
    rows=[];records=[]
    for prime in primes:
        F=GF(prime);e=EllipticCurve(F,[F(c) for c in model])
        if not e.discriminant():raise ArithmeticError('good finite specialization required')
        key=lambda point:tuple(int(c) for c in point)
        elements=e.points();multiples={key(ell*P):ell*P for P in elements};mask={key(P):0 for P in multiples.values()};reps=[e(0)]
        while len(mask)<len(elements):
            P=next(P for P in elements if key(P) not in mask);old=list(reps);size=len(old)
            for digit in range(1,ell):
                for i,R in enumerate(old):
                    rep=R+digit*P;reps.append(rep)
                    for T in multiples.values():
                        k=key(rep+T)
                        if k in mask:raise ArithmeticError('quotient cosets overlap')
                        mask[k]=i+size*digit
        dimension=ZZ(len(reps)).valuation(ell)
        if ell**dimension!=len(reps) or dimension>2:raise ArithmeticError('elliptic quotient dimension differs')
        reduced=[e([F(c) for c in P]) for P in projective]
        for j in range(dimension):rows.append([(mask[key(P)]//ell**j)%ell for P in reduced])
        records.append({'prime':prime,'group_order':len(elements),'ell_multiple_subgroup_order':len(multiples),'quotient_dimension':int(dimension)})
    rank=int(matrix(GF(ell),rows).rank())
    if rank!=13:raise ArithmeticError('selected generic seed is not injective in finite ell quotients')
    return {'modulus':ell,'rank':rank,'groups':records}

def main(path):
 d=json.loads(path.read_text());c=d['construction'];control=d['control'];parent=d['parent'];rank=d['parent_rank'];assert d['parent_replay']['status']=='PASS' and rank['generic_Q_MW_rank']==12 and rank['remaining_index']==1
 R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();P=PolynomialRing(K,'x');x=P.gen();quartic=P([K(R(a)) for a in parent['quartic_coefficients']]);aa=R(quartic[4]);constant=QQ(c['conic_constant']);lam=QQ(c['leading_square_root'])
 assert aa==lam*lam*(T*T+constant) and constant!=0
 S=PolynomialRing(QQ,'s');s=S.gen();L=S.fraction_field();t=L(c['base_T']);w=L(c['quartic_infinity_ordinate_leading_coefficient']);scale=L(c['integral_function_field_scale'])
 assert t==(s*s-constant)/(2*s) and w==lam*(s*s+constant)/(2*s) and w*w==aa(t)
 A=S(c['A']);B=S(c['B']);assert A==scale**4*R(parent['raw_A'])(t) and B==scale**6*R(parent['raw_B'])(t)
 delta=-16*(4*A**3+27*B**2);assert (A.degree(),B.degree(),delta.degree())==(16,24,44) and delta.gcd(A)==1
 assert delta.valuation(s)==4;rest=delta//s**4;assert rest.degree()==40 and rest.gcd(rest.derivative())==1
 # chi4, fortyI1 and twoI4. Hence any nonzero geometric section has height>=6.
 assert c['chi']==4 and c['finite_fibres']=='40I1+I4 at s0' and c['infinity_fibre']=='I4'
 E=EllipticCurve(L,[A,B]);points=[E([L(X),L(Y)]) for X,Y in c['sections']];old=points[:12];Q=points[12]
 for p,(X,Y) in zip(old,rank['basis']):assert p==E([scale**2*K(X)(t),scale**3*K(Y)(t)])
 X0,Y0=[K(v) for v in parent['quartic_points'][0]];shift=quartic(x+X0);ee,dd,cc,bb,a=shift.list()
 assert Q==E([scale**2*(18*Y0(t)*w+3*cc(t)),scale**3*27*(Y0(t)*bb(t)+dd(t)*w)])
 phis=[L(s),L(-s),L(constant/s),L(-constant/s)]
 def act(point,phi):
  ratio=scale/scale(phi);assert A==ratio**4*A(phi) and B==ratio**6*B(phi)
  return E([ratio**2*point[0](phi),ratio**3*point[1](phi)])
 assert points[12:]==[act(Q,phi) for phi in phis]
 deck=phis[3];assert t(deck)==t and w(deck)==-w and all(act(p,deck)==p for p in old)
 assert act(Q,deck)==old[0]-Q
 anti=2*Q-old[0];assert anti!=E(0) and act(anti,deck)==-anti
 four=4*anti;n,den=four[0].numerator(),four[0].denominator();twice=max(den.degree(),n.degree()-8);assert twice>=0 and twice%2==0
 antiheight=QQ(8+twice)/16;assert antiheight==10
 # Rank13: apply1-deck to a relation. Primitive anti direction: if nP=anti
 # for n>=2 then h(P)<=10/4<6, impossible. This plus the full parent group
 # proves saturation of the13-dimensional span, even if ambient rank is higher.
 H=2*matrix(QQ,rank['basis_height_gram']);G=matrix(QQ,13);G[:12,:12]=H
 for i in range(12):G[i,12]=G[12,i]=H[i,0]/2
 G[12,12]=(H[0,0]+antiheight)/4
 assert G==matrix(QQ,c['basis_height_gram']) and G.is_positive_definite() and G.det()==1935360
 for rel in c['exact_relations']:
  assert rel['multiplier']*points[rel['target']]==sum((ZZ(a)*p for a,p in zip(rel['coefficients'],points[:13])),E(0))
 print('PASS explicit conic, all infinity images, primitive rank13 subgroup and anti-height10',flush=True)
 dc=ZZ(control['conic_squarefree_part']);rs=QQ(control['conic_scale']);v=ZZ(control['balanced_integer_parameter']);s0=QQ(control['conic_parameter_s'])
 assert constant==dc*rs*rs and dc.is_squarefree() and rs>0 and (v-1)**2<dc<v*v and s0==rs*v
 assert QQ(control['old_parent_T'])==t(s0)
 Ef=EllipticCurve(QQ,control['curve']);u=QQ(control['point_transport_from_polynomial_model'])
 assert Ef.a4()==u**4*A(s0) and Ef.a6()==u**6*B(s0) and Ef.discriminant()!=0
 specialization=[Ef([u*u*p[0](s0),u**3*p[1](s0)]) for p in points[:13]]
 assert [[str(z) for z in p.xy()] for p in specialization]==control['points']
 proof=control['rank_certificate'];finite_rank(control['curve'],control['points'],[r['prime'] for r in proof['signatures']],2)
 prime=proof['no_rational_2_torsion_prime'];F=GF(prime);Rp=PolynomialRing(F,'z');z=Rp.gen();assert (z**3+F(Ef.a4())*z+F(Ef.a6())).is_irreducible()
 assert max(abs(ZZ(a)).nbits() for a in Ef.a_invariants())==control['coefficient_bits']==998
 print('PASS balanced v774161 control,998-bit model and independent rational rank13',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);a=p.parse_args();main(a.input)
