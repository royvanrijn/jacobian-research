#!/usr/bin/env sage-python
"""Fixed conic base change on the first new Kihara parent; all infinity images."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/kihara-split-infinity-v1';OUT=ART/'kihara_split_infinity_v1.json'
def height(P):
 if P.is_zero():return QQ(0)
 Q=4*P;assert not Q.is_zero();x=Q[0];n,d=x.numerator(),x.denominator();twice=max(d.degree(),n.degree()-8);assert twice>=0 and twice%2==0
 return QQ(8+twice)/16

def main():
 protocol=json.loads((D/'protocol.json').read_text())
 for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
 parent=json.loads((ART/'kihara_five_parent_distinctness_v1.json').read_text())['rows'][1]
 original=json.loads((ART/'kihara_first_parent_rank_v1.json').read_text());assert parent['path_parameter']=='3/2'
 R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();P=PolynomialRing(K,'x');x=P.gen();quartic=P([K(R(f)) for f in parent['quartic_coefficients']]);a=R(quartic[4]);assert a.degree()==2 and a[1]==0 and a[2].is_square()
 lam=a[2].sqrt();c=a[0]/a[2];assert c and not (-c).is_square()
 S=PolynomialRing(QQ,'s');s=S.gen();L=S.fraction_field();t=(s*s-c)/(2*s);w=lam*(s*s+c)/(2*s);scale=L((2*s)**2)
 A=S(scale**4*R(parent['raw_A'])(t));B=S(scale**6*R(parent['raw_B'])(t));delta=-16*(4*A**3+27*B**2)
 assert (A.degree(),B.degree(),delta.degree())==(16,24,44)
 assert delta.valuation(s)==4 and (delta//s**4).gcd((delta//s**4).derivative())==1
 assert delta.gcd(A)==1
 E=EllipticCurve(L,[A,B]);old=[E([scale**2*K(X)(t),scale**3*K(Y)(t)]) for X,Y in original['basis']]
 X0,Y0=[K(v) for v in parent['quartic_points'][0]];shift=quartic(x+X0);ee,dd,cc,bb,aa=shift.list()
 Qraw=[18*Y0(t)*w+3*cc(t),27*(Y0(t)*bb(t)+dd(t)*w)]
 Q=E([scale**2*Qraw[0],scale**3*Qraw[1]])
 # Automorphisms over u=T^2: s, -s, c/s, -c/s. Transport weighted short models exactly.
 maps=[L(s),L(-s),L(c/s),L(-c/s)];infinity=[]
 for phi in maps:
  ratio=scale/scale(phi);assert A==ratio**4*A(phi) and B==ratio**6*B(phi)
  infinity.append(E([ratio**2*Q[0](phi),ratio**3*Q[1](phi)]))
 cloud=old+infinity;H=matrix(QQ,len(cloud));hs=[height(p) for p in cloud]
 for i in range(len(cloud)):
  H[i,i]=hs[i]
  for j in range(i):
   hp,hm=height(cloud[i]+cloud[j]),height(cloud[i]-cloud[j]);assert hp+hm==2*(hs[i]+hs[j]);H[i,j]=H[j,i]=(hp-hs[i]-hs[j])/2
  print('HEIGHT',i+1,'of',len(cloud),flush=True)
 assert H[:12,:12]==2*matrix(QQ,original['basis_height_gram'])
 indices=list(H.pivots());G=H.matrix_from_rows_and_columns(indices,indices);assert G.is_positive_definite();basis=[cloud[i] for i in indices];relations=[]
 for j,Q in enumerate(cloud):
  v=G.solve_right(vector(QQ,[H[i,j] for i in indices]));den=lcm([z.denominator() for z in v]);word=[ZZ(den*z) for z in v]
  assert den*Q==sum((a*p for a,p in zip(word,basis)),E(0));relations.append({'target':j,'multiplier':int(den),'coefficients':list(map(int,word))})
 result={'schema':'elliptic-curves.kihara-split-infinity.v1','status':'PASS','parent_parameter':'3/2','conic_constant':str(c),'leading_square_root':str(lam),'base_T':str(t),'quartic_infinity_ordinate_leading_coefficient':str(w),'integral_function_field_scale':str(scale),
  'A':list(map(str,A.list())),'B':list(map(str,B.list())),'chi':4,'finite_fibres':'40I1+I4 at s0','infinity_fibre':'I4','curve_coefficient_bits':max(abs(z.numerator()).nbits() for z in list(A)+list(B)),
  'sections':[[str(v) for v in p.xy()] for p in cloud],'height_gram':[[str(c) for c in r] for r in H.rows()],'basis_indices':indices,'basis_height_gram':[[str(c) for c in r] for r in G.rows()],
  'basis_height_determinant':str(G.det()),'exact_relations':relations,'supplied_generic_rank':len(indices),'new_independent_generic_directions':len(indices)-12,
  'full_generic_rank':'UNKNOWN','full_saturation':'UNKNOWN','sources':protocol['sources'],
  'scope':'Explicit rational conic base change of one proved new K3 parent, with all four infinity images under its V4 action retained before rank calculations. Exact generic section identities and height/group-law proofs only. The chi4 surface is not another K3 fibration or a new Q-isomorphic parent; no fibre rank gain, point search or parameter population is claimed.'}
 assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n');print('PASS GENERIC RANK',len(indices),'GAIN',len(indices)-12,'CHI4','HEIGHT DET',G.det(),flush=True)
if __name__=='__main__':main()
