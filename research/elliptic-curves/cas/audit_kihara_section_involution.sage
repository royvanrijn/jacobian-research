#!/usr/bin/env sage-python
"""Fixed first Kihara parent: exact heights and its base-involution closure."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/kihara-section-involution-v1'
SOURCE=ART/'kihara_five_parent_distinctness_v1.json';OUT=ART/'kihara_section_involution_v1.json'
def inf(f,w):
 if not f:return QQ(0)
 n,d=f.numerator(),f.denominator();v=d.degree()-n.degree()+w
 return None if v<0 else QQ(0) if v>0 else n.leading_coefficient()/d.leading_coefficient()
def main():
 protocol=json.loads((D/'protocol.json').read_text())
 for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
 row=json.loads(SOURCE.read_text())['rows'][1];assert row['path_parameter']=='3/2'
 R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(row['raw_A']);B=R(row['raw_B']);delta=-16*(4*A**3+27*B**2)
 assert (A.degree(),B.degree(),delta.degree())==(8,12,20) and delta.gcd(delta.derivative())==1 and delta.gcd(A)==1
 assert A(-T)==A and B(-T)==B
 node=-3*B[12]/(2*A[8]);tangent=QQ(3*node).sqrt();assert tangent in QQ and tangent and A[6]*node+B[10]==0
 E=EllipticCurve(K,[A,B]);P=[E([K(x),K(y)]) for x,y in row['generic_sections']]
 def profile(P):
  if P.is_zero():return {'height':'0','zero':True,'component':0}
  x,y=P.xy();n,d=x.numerator(),x.denominator();twice=max(d.degree(),n.degree()-4);assert twice>=0 and twice%2==0
  component=0
  if inf(x,4)==node and inf(y,6)==0:
   f=x-node*T**4;vx=10**9 if not f else 4+f.denominator().degree()-f.numerator().degree();vy=10**9 if not y else 6+y.denominator().degree()-y.numerator().degree()
   if min(vx,vy)>=2:component=2
   else:
    slope=inf(y/f,2);assert slope in [tangent,-tangent];component=1 if slope==tangent else 3
  h=4+twice-QQ(component*(4-component))/4;assert h>0
  return {'height':str(h),'zero':False,'component':component,'O_intersection':int(twice//2)}
 # Every reflected supplied section is retained before any rank calculation.
 Q=P+[E([c(-T) for c in point.xy()]) for point in P];profiles=[profile(point) for point in Q];G=matrix(QQ,24,24)
 for i in range(24):
  G[i,i]=QQ(profiles[i]['height']);doubled=profile(2*Q[i]);assert QQ(doubled['height'])==4*G[i,i] and doubled['component']==2*profiles[i]['component']%4
  for j in range(i):
   plus,minus=profile(Q[i]+Q[j]),profile(Q[i]-Q[j]);assert QQ(plus['height'])+QQ(minus['height'])==2*(G[i,i]+G[j,j])
   assert plus['component']==(profiles[i]['component']+profiles[j]['component'])%4
   G[i,j]=G[j,i]=(QQ(plus['height'])-G[i,i]-G[j,j])/2
  print('HEIGHT ROW',i+1,'of24',flush=True)
 indices=list(G.pivots());H=G.matrix_from_rows_and_columns(indices,indices);assert H.is_positive_definite();basis=[Q[i] for i in indices];relations=[]
 for j,point in enumerate(Q):
  v=H.solve_right(vector(QQ,[G[i,j] for i in indices]));den=lcm([x.denominator() for x in v]);word=[ZZ(den*x) for x in v]
  assert den*point==sum((a*P for a,P in zip(word,basis)),E(0));relations.append({'target':j,'multiplier':int(den),'coefficients':list(map(int,word))})
 action=matrix(QQ,len(indices));columns=[]
 for j,i in enumerate(indices):
  reflected=i+12 if i<12 else i-12;r=relations[reflected];v=vector(QQ,r['coefficients'])/r['multiplier'];action.set_column(j,v)
 assert action*action==matrix.identity(QQ,len(indices)) and action.transpose()*H*action==H
 plus=len(indices)-(action-matrix.identity(QQ,len(indices))).rank();minus=len(indices)-plus
 a=R([A[2*i] for i in range(5)]);b=R([B[2*i] for i in range(7)]);dd=-16*(4*a**3+27*b**2)
 assert (a.degree(),b.degree(),dd.degree())==(4,6,10) and dd.gcd(dd.derivative())==1 and dd(0)
 # Rational elliptic quotient: ten I1 and split I2, hence geometric MW7.
 assert plus<=7
 result={'schema':'elliptic-curves.kihara-section-involution.v1','status':'PASS','parent_parameter':'3/2','input_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
  'raw_A':row['raw_A'],'raw_B':row['raw_B'],'sections':[[str(c) for c in p.xy()] for p in Q], 'section_profiles':profiles,'height_gram':[[str(c) for c in r] for r in G.rows()],
  'basis_indices':indices,'basis_height_gram':[[str(c) for c in r] for r in H.rows()],'basis_height_determinant':str(H.det()),'exact_relations':relations,
  'involution_action':[[str(c) for c in r] for r in action.rows()],'invariant_rank':int(plus),'anti_invariant_rank':int(minus),
  'supplied_rational_MW_rank':len(indices),'rational_NS_lower_bound':len(indices)+5,'known_rational_NS_absolute_determinant':str(4*H.det()),
  'quotient_A':list(map(str,a.list())),'quotient_B':list(map(str,b.list())),'quotient_geometric_MW_rank':7,'geometric_NS_lower_bound':int(12+minus),
  'full_rational_NS_rank':'UNKNOWN','full_geometric_NS_rank':'UNKNOWN','scope':'Fixed first fresh Kihara K3 parent and closure of all twelve retained rational sections under T->-T. Exact local heights, parallelogram identities and function-field group relations; quotient geometry gives a geometric lower bound. No new fibre or point search, no full saturation or NS identification, and no generic-rank-only yield claim.'}
 assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n');print('PASS MW',len(indices),'NS det',4*H.det(),'involution',plus,minus,'geom NS >=',12+minus,flush=True)
if __name__=='__main__':main()
