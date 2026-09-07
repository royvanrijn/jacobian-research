#!/usr/bin/env sage-python
"""Independent identities and complete semistable fibre accounting, nine pencils."""
import json,hashlib,argparse
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
from itertools import combinations,product
def compute(bundle):
 ns=bundle['ns'];results=[]
 verify_roster(bundle)
 for index,row in enumerate(bundle['presentations']):
  R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(row['old_A']);B=R(row['old_B']);E=EllipticCurve(K,[A,B]);P=E([K(c) for c in row['pole_point']])
  points=[E([K(c) for c in coords]) for coords in ns['section_points']]
  if 'selection' in row:assert P==sum((int(w)*point for w,point in zip(row['selection']['word'],points)),E(0))
  else:assert P==points[0]
  q=R(row['q']);a=R(row['a']);b=R(row['b']);c=R(row['c']);assert q.degree()==2 and P[0]==a/q**2 and P[1]==b/q**3
  assert b*b==a**3+A*a*q**4+B*q**6 and (a*c-b)%(q*q)==0
  Z=PolynomialRing(QQ,'z');z=Z.gen();L=Z.fraction_field();TR=PolynomialRing(L,'T');tt=TR.gen();F=TR.fraction_field()
  quartic=sum(TR(Z(coeff))*tt**i for i,coeff in enumerate(row['quartic_T_coefficients']));qp=TR(q);ap=TR(a);bp=TR(b);cp=TR(c)
  AA=TR(A);BB=TR(B);xp=ap/qp**2;yp=bp/qp**3;m=qp*z-cp/qp
  Wring=PolynomialRing(F,'W');W=Wring.gen();xx=(qp*W-xp+m*m)/2;yy=m*(xx-xp)-yp
  assert (yy*yy-xx**3-AA*xx-BB)%(W*W-quartic)==0
  assert yy+yp==m*(xx-xp) and (qp*m+cp)/qp**2==z
  for t,y in row['quartic_sections']:assert quartic(L(t))==L(y)**2
  e,d,cc,bb,aa=[Z(c) for c in row['quartic_T_coefficients']];I=12*aa*e-3*bb*d+cc*cc;J=72*aa*cc*e+9*bb*cc*d-27*aa*d*d-27*bb*bb*e-2*cc**3
  newA=Z(row['Jacobian_A']);newB=Z(row['Jacobian_B']);assert newA==-27*I and newB==-27*J
  delta=-16*(4*newA**3+27*newB**2);assert newA.degree()==8 and newB.degree()==12 and delta.degree()==22 and delta.gcd(newA)==1
  reconstructed=Z(QQ(row['discriminant_unit']));fibres=[];geometric=0;arithmetic=0;total=2
  factors=[Z(f['coefficients']) for f in row['discriminant_factors']]
  assert all(f.gcd(g)==1 for f,g in combinations(factors,2))
  for item in row['discriminant_factors']:
   f=Z(item['coefficients']);n=item['multiplicity'];reconstructed*=f**n;total+=f.degree()*n
   if n==1:
    assert f.gcd(f.derivative())==1
    continue
   assert f.degree()==1;r=-f[0]/f[1];node=-3*newB(r)/(2*newA(r));assert 3*node**2+newA(r)==0 and node**3+newA(r)*node+newB(r)==0
   split=bool((3*node).is_square());g=n-1;rat=g if split else n//2
   geometric+=g;arithmetic+=rat;fibres.append({'base':str(r),'type':'I'+str(n),'split':split,'geometric_component_rank':g,'rational_component_rank':rat,'tangent_squareclass_representative':str(3*node)})
  assert reconstructed==delta and total==24
  node=-3*newB[12]/(2*newA[8]);assert 3*node**2+newA[8]==0 and node**3+newA[8]*node+newB[12]==0
  geometric+=1;arithmetic+=1;fibres.append({'base':'infinity','type':'I2','split':bool((3*node).is_square()),'geometric_component_rank':1,'rational_component_rank':1})
  rankQ=18-2-arithmetic;rankgeom=19-2-geometric
  assert rankQ<=11
  results.append({'id':('baseline' if 'pencil_index' not in row else 'pencil'+str(row['pencil_index'])),'reducible_fibres':fibres,'geometric_generic_MW_rank':rankgeom,'rational_generic_MW_rank':rankQ,'rational_point_search_admission':False})
  expected=bundle['geometry']['rows'][index].copy();expected.pop('source');expected.pop('source_sha256');assert results[-1]==expected
  print(results[-1]['id'],'generic Q rank',rankQ,'geometric',rankgeom,flush=True)
 return results

def verify_roster(bundle):
 ns=bundle['ns'];G=matrix(QQ,ns['rational_NS_Gram']);C=G[2:7,2:7]
 H=-(G[7:,7:]-G[7:,:7]*G[:7,:7].inverse()*G[:7,7:])
 std=matrix.identity(QQ,18).rows();F,O=std[:2];phi=[]
 for i,p in enumerate(ns['section_profiles']):
  v=std[7+i]-O-(2+p['zero_section_intersection'])*F;v[2:7]-=C.inverse()*vector(QQ,G[2:7,7+i].column(0));phi.append(v)
 visible=[std[i] for i in range(1,18)]+[F-std[2],F-std[3],F-sum(std[4:7])]
 words=[];found=[]
 for support in (1,2,3):
  for indices in combinations(range(11),support):
   for signs in product((-1,1),repeat=support-1):
    w=vector(ZZ,11);w[indices[0]]=1
    for i,c in zip(indices[1:],signs):w[i]=c
    words.append(w)
    comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
    corr=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4;oo=(w*H*w-4+corr)/2
    if oo!=2:continue
    v=O+4*F+sum((w[i]*phi[i] for i in range(11)),vector(QQ,18));v[2:7]+=C.inverse()*vector(QQ,comp[:2]+[int(comp[2]==k) for k in (1,2,3)])
    assert v*G*v==-2 and v*G*O==2 and all(c.denominator()==1 for c in v)
    D=O+v;vertical=[p for p in visible+[v] if D*G*p==0];bound=17-matrix(QQ,vertical).rank()
    found.append((int(bound),[int(c) for c in w],[int(c) for c in comp]))
 rows=[{'upper_bound_from_visible_curves':a,'word':b,'components':c} for a,b,c in sorted(found,reverse=True)]
 assert len(words)==781 and len(rows)==151 and bundle['roster']=={'word_count':781,'rows':rows}
 selected=[r for r in rows if r['upper_bound_from_visible_curves']>11][:8]
 assert selected==bundle['admission']['rows'] and len(selected)==8
 assert [p['selection'] for p in bundle['presentations'][1:]]==selected
 print('PASS all781 words,151 eligible divisors and fixed eight-pencil admission',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path('mestre_chord_pencils_v1.json'));args=p.parse_args()
 result=compute(json.loads(args.input.read_text()));assert len(result)==9
 print('PASS9 INDEPENDENT BIRATIONAL MAPS AND EXACT GENERIC RANKS',flush=True)
