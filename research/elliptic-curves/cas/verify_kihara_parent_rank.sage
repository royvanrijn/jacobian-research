#!/usr/bin/env sage-python
"""Standalone independent Kihara Picard, heights and full rational basis replay.
No repository imports. Finite groups and surface fibres are recounted in Sage;
heights use multiples meeting identity components, not producer local profiles.
"""
import json,argparse
from pathlib import Path
from sage.all import QQ,ZZ,AA,GF,PolynomialRing,EllipticCurve,matrix,vector,lcm,gcd,cyclotomic_polynomial,prime_range
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
    if rank!=12:raise ArithmeticError('selected generic seed is not injective in finite ell quotients')
    return {'modulus':ell,'rank':rank,'groups':records}

def height(P):
 if P.is_zero():return QQ(0)
 Q=4*P;assert not Q.is_zero();x=Q[0];n,d=x.numerator(),x.denominator();twice=max(d.degree(),n.degree()-4);assert twice>=0 and twice%2==0
 return QQ(4+twice)/16

def recount(A,B,p,k,stored):
 if k==1:F=GF(p)
 else:
  d=stored['quadratic_nonsquare'];R=PolynomialRing(GF(p),'z');assert not GF(p)(d).is_square();F=GF(p*p,'z',modulus=R.gen()**2-d)
 R=PolynomialRing(F,'T');a,b=R(A),R(B);delta=-16*(4*a**3+27*b**2);q=p**k
 assert delta.degree()==20 and delta.gcd(delta.derivative())==1 and delta.gcd(a)==1
 counts=[]
 for i in range(q+1):
  if i==q:av,bv=a[8],b[12]
  else:
   t=F(i) if k==1 else F(i%p)+F(i//p)*F.gen();av,bv=a(t),b(t)
  if 4*av**3+27*bv**2:n=EllipticCurve(F,[av,bv]).cardinality(algorithm='pari')
  else:
   # Singular nodal cubic: q if split, q+2 if nonsplit (including infinity).
   node=-3*bv/(2*av);assert node**3+av*node+bv==0 and 3*node*node+av==0
   n=q if (3*node).is_square() else q+2
  counts.append(int(n))
 assert counts==stored['weierstrass_fibre_counts']
 assert stored['resolution_corrections']==[{'base_index':q,'correction':3*q}]
 assert sum(counts)+3*q==stored['surface_point_count']
 return sum(counts)+3*q

def main(path):
 d=json.loads(path.read_text());h=d['height'];r=d['rank'];parent=d['parent'];seed=d['seed'];counts=d['counts']
 R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(h['raw_A']);B=R(h['raw_B']);delta=-16*(4*A**3+27*B**2)
 assert (A.degree(),B.degree(),delta.degree())==(8,12,20) and delta.gcd(delta.derivative())==1 and delta.gcd(A)==1
 assert A(-T)==A and B(-T)==B;node=-3*B[12]/(2*A[8]);assert (3*node).is_square()
 E=EllipticCurve(K,[A,B]);P=[E([K(x),K(y)]) for x,y in h['sections']];assert len(P)==24
 assert P[12:]==[E([x(-T),y(-T)]) for x,y in [point.xy() for point in P[:12]]]
 H=matrix(QQ,12);hs=[height(p) for p in P[:12]]
 for i in range(12):
  H[i,i]=hs[i]
  for j in range(i):H[i,j]=H[j,i]=(height(P[i]+P[j])-hs[i]-hs[j])/2
 assert H==matrix(QQ,h['basis_height_gram']) and H.is_positive_definite() and H.det()==6804
 for rel in h['exact_relations']:
  assert rel['multiplier']*P[rel['target']]==sum((ZZ(c)*p for c,p in zip(rel['coefficients'],P[:12])),E(0))
 M=matrix(QQ,h['involution_action']);assert M*M==matrix.identity(QQ,12) and M.transpose()*H*M==H
 for i in range(12):assert sum((M[j,i]*P[j] for j in range(12)),E(0))==P[12+i]
 assert 12-(M-matrix.identity(QQ,12)).rank()==6 and 12-(M+matrix.identity(QQ,12)).rank()==6
 a=R(h['quotient_A']);b=R(h['quotient_B']);dd=-16*(4*a**3+27*b**2)
 assert A==a(T*T) and B==b(T*T) and (a.degree(),b.degree(),dd.degree())==(4,6,10)
 assert dd.gcd(dd.derivative())==1 and dd(0)!=0
 # Rational elliptic surface10I1+I2: geometric MW=10-2-1=7.
 assert h['quotient_geometric_MW_rank']==7 and h['geometric_NS_lower_bound']==18
 print('PASS78 independent generic heights and exact6+6 symmetry; geometric NS>=18',flush=True)
 new=[E([K(x),K(y)]) for x,y in r['basis']];assert new[1:]==P[1:12] and 6*new[0]==sum(P[:11],E(0))
 C=matrix(QQ,r['basis_change']);assert C.det()==QQ(1)/6
 for i in range(12):
  den=lcm([c.denominator() for c in C.row(i)]);assert den*new[i]==sum((ZZ(den*c)*p for c,p in zip(C.row(i),P[:12])),E(0))
 Hnew=C*H*C.transpose();assert Hnew==matrix(QQ,r['basis_height_gram']) and Hnew.det()==189
 # Specialization is independently transported to the already recorded fibre.
 oldE=EllipticCurve(QQ,seed['old_model']);oldP=[oldE([QQ(x),QQ(y)]) for x,y in seed['old_points']]
 s=seed['seed'];Ef=EllipticCurve(QQ,s['curve']);original=[Ef([QQ(x),QQ(y)]) for x,y in s['original_points']];full=[Ef([QQ(x),QQ(y)]) for x,y in s['points']]
 t=QQ(parent['control_T']);Et=EllipticCurve(QQ,[A(t),B(t)]);special=[Et([x(t),y(t)]) for x,y in [p.xy() for p in P[:12]]]
 iso1=next(iso for iso in Et.isomorphisms(oldE) if all(iso(p)==q-oldP[0] for p,q in zip(special,oldP[1:13])))
 iso2=next(iso for iso in oldE.isomorphisms(Ef) if all(iso(p)==q for p,q in zip(oldP,original)))
 specialized=[iso2(iso1(Et([c(t) for c in p.xy()]))) for p in new]
 matrix14=matrix(ZZ,r['specialized_basis_in_full14_coefficients'])
 for i,p in enumerate(specialized):assert p==sum((c*q for c,q in zip(matrix14.row(i),full)),Ef(0))
 for ell in [2,3]:
  primes=[v['prime'] for v in (seed['finite']['signatures'] if ell==2 else next(a for a in d['odd']['audits'] if a['modulus']==3)['signatures'])]
  finite_rank(s['curve'],[[str(c) for c in p.xy()] for p in specialized],primes,ell)
  print('PASS complete finite quotient rank12 modulo',ell,flush=True)
 def good(p):
  try:
   Rp=PolynomialRing(GF(p),'z');ap,bp=Rp(A),Rp(B);dp=-16*(4*ap**3+27*bp**2)
  except (ValueError,ZeroDivisionError):return False
  if (ap.degree(),bp.degree(),dp.degree())!=(8,12,20) or dp.gcd(dp.derivative())!=1 or dp.gcd(ap)!=1 or not dp(0):return False
  node=-3*bp[12]/(2*ap[8]);return not (3*node*node+ap[8]) and not (node**3+ap[8]*node+bp[12]) and (3*node).is_square()
 selected=[int(p) for p in prime_range(31,252) if good(p)][:2];assert selected==counts['primes']==[53,83]
 for row in counts['rows']:
  p=row['prime'];F=GF(p);Rp=PolynomialRing(F,'s');aq,bq=Rp(a),Rp(b);dq=-16*(4*aq**3+27*bq**2)
  assert dq.degree()==10 and dq.gcd(dq.derivative())==1
  quotient=[]
  for i in range(p+1):
   av,bv=(aq[4],bq[6]) if i==p else (aq(F(i)),bq(F(i)))
   if 4*av**3+27*bv**2:n=EllipticCurve(F,[av,bv]).cardinality()
   else:
    node=-3*bv/(2*av);n=p if (3*node).is_square() else p+2
   quotient.append(int(n))
  assert quotient==row['quotient']['weierstrass_fibre_counts']
  total=sum(quotient)+p;assert total==row['quotient']['surface_point_count']
  sign=(total-1-p*p-9*p)/p;assert sign==row['quotient']['remaining_geometric_quotient_divisor_sign']==-1
  n1,n2=[recount(A,B,p,k,c) for k,c in enumerate(row['counts'],1)]
  a1=ZZ(n1-1-p*p-16*p);a2=ZZ(n2-1-p**4-18*p*p);b2=QQ(a1*a1-a2)/2;assert b2.denominator()==1 and b2!=0
  z=R.gen();f=z**4-a1*z**3+b2*z*z-p*p*a1*z+p**4
  candidate=row['frobenius']['candidates'];assert len(candidate)==1 and list(map(str,f.list()))==candidate[0]['residual_polynomial']
  g=z*z-QQ(a1)/p*z+QQ(b2)/(p*p)-2;roots=g.roots(AA);assert sum(m for x,m in roots)==2 and all(abs(x)<=2 for x,m in roots)
  normalized=R(f(p*z)/p**4)
  for n in [1,2,3,4,5,6,8,10,12]:assert normalized.gcd(R(cyclotomic_polynomial(n))).degree()==0
  print('PASS p',p,'complete Fp/Fp2 counts; rational/geometric Picard upper17/18',flush=True)
 assert r['rational_NS_rank']==17 and r['geometric_NS_rank']==18 and r['generic_Q_MW_rank']==12 and r['generic_Qbar_MW_rank']==13
 det=4*Hnew.det();assert det==756
 possible=[n for n in range(1,28) if ZZ(det)%(n*n)==0];assert possible==r['possible_indices_before_finite_checks']==[1,2,3,6]
 assert r['remaining_index']==r['geometric_torsion_order']==1 and r['all_Q_Jacobian_fibrations_generic_MW_upper_bound']==15
 print('PASS full rational NSdet756; MW12/13; every Q-Jacobian fibration MW<=15',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);args=p.parse_args();main(args.input)
