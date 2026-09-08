#!/usr/bin/env sage-python
"""Independent identities and complete semistable fibre accounting, nine pencils."""
import json,hashlib,argparse
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/mestre-chord-pencils-v1'
OUT=ART/'mestre_chord_pencil_geometry_v1.json'
def compute():
 ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
 paths=[ART/'mestre_u11_visible_two_neighbor_v1.json']+[D/('pencil'+str(i)+'.json') for i in range(8)];results=[]
 for path in paths:
  row=json.loads(path.read_text());R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(row['old_A']);B=R(row['old_B']);E=EllipticCurve(K,[A,B]);P=E([K(c) for c in row['pole_point']])
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
  results.append({'source':str(path.relative_to(ROOT)),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'id':('baseline' if 'pencil_index' not in row else 'pencil'+str(row['pencil_index'])),'reducible_fibres':fibres,'geometric_generic_MW_rank':rankgeom,'rational_generic_MW_rank':rankQ,'rational_point_search_admission':False})
  print(results[-1]['id'],'generic Q rank',rankQ,'geometric',rankgeom,flush=True)
 return {'schema':'elliptic-curves.mestre-chord-pencil-geometry.v1','status':'PASS','rows':results,'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),ART/'mestre_rational_ns_gram_v2.json',ART/'mestre_468_portable_replay_v2.json']},'scope':'Exact birational chord inverse, rational quartic sections, Jacobian invariants and complete minimal semistable fibre decomposition. Picard ranks18/19 transport through the birational K3 identification and the rational quartic section; Shioda-Tate yields exact rational and geometric generic ranks. Nine presentations, not claimed pairwise inequivalent. No rational point search or high-rank addition.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();r=compute()
 if args.check:assert r==json.loads(OUT.read_text())
 else:
  assert not OUT.exists();OUT.write_text(json.dumps(r,indent=2)+'\n')
 print('PASS9 COMPLETE PENCIL GEOMETRIES',flush=True)
