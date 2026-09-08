#!/usr/bin/env sage-python
"""Independent exact maps, prior-box inclusion and invariant coefficient bounds."""
import json,hashlib
from pathlib import Path
from sage.all import QQ,ZZ,EllipticCurve,PolynomialRing,matrix,lcm,gcd,ceil
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/bounded-prime-point-portfolio-v1';OUT=ROOT/'artifacts/generated-results/elliptic-curves/bounded_prime_history_and_height_v1.json'
def read(p):return json.loads(p.read_bytes())
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(v):
 v=list(map(QQ,v));d=lcm([x.denominator() for x in v]);n=[ZZ(x*d) for x in v];g=gcd(n);n=[x//g for x in n]
 if next(x for x in n if x)<0:n=[-x for x in n]
 return list(map(int,n))
def bits(q):
 q=QQ(q);return max(abs(q.numerator()).nbits(),q.denominator().nbits())
def main():
 assert not OUT.exists();p=read(D/'protocol.json');R=PolynomialRing(QQ,'z');z=R.gen();rows=[];inputs={}
 paths=[D/'protocol.json',D/'ledger.json'];assert read(D/'ledger.json')['status']=='PASS'
 for row in p['rows']:
  folder=D/row['id'];seed=read(folder/'seed.json');maps=read(folder/'maps.json');old=read(ROOT/row['base_maps']);paths += [folder/'seed.json',folder/'maps.json',ROOT/row['base_maps']]
  prior=[(old,125000)]+[(read(ROOT/h['maps_path']),h['height']) for h in row['prior_same_centre_boxes']]
  E=EllipticCurve(QQ,seed['curve']);P=[E([QQ(x),QQ(y)]) for x,y in seed['points']];j=E.j_invariant();nb=int(abs(j.numerator()).nbits());db=int(j.denominator().nbits());short_lb=max(1,int(ceil(QQ(nb-13)/3)),int(ceil(QQ(db-5)/3)));quartic_lb=max(1,int(ceil(QQ(nb-25)/6)),int(ceil(QQ(db-16)/6)))
  new=covered=0;quartic_bits=[]
  for i,m in enumerate(maps['rows']):
   C=sum((int(c)*Q for c,Q in zip(m['centre']['representative'],P)),E(0));assert C
   x,y=C.xy();raw=[-3*x*x-4*E.a4(),-8*y,-6*x,0,1];assert list(map(str,raw))==m['raw_coefficients']
   a,b,c,d=map(QQ,m['matrix']);phi=matrix(QQ,2,[a,b,c,d]);assert phi.det()!=0
   F=R(m['discriminant_quartic']);actual=sum(raw[k]*(a*z+b)**k*(c*z+d)**(4-k) for k in range(5));ratio=QQ(m['square_ratio']);assert ratio>0 and ratio.is_square() and actual==ratio*F
   vals=[F[k] for k in range(5)];assert all(v.denominator()==1 for v in vals);aa,bb,cc,dd,ee=vals
   I=12*aa*ee-3*bb*dd+cc**2;J=72*aa*cc*ee+9*bb*cc*dd-27*aa*dd**2-27*bb**2*ee-2*cc**3
   assert 4*I**3-J**2!=0 and 6912*I**3/(4*I**3-J**2)==j
   height=max(int(abs(v.numerator()).nbits()) for v in vals);assert height>=quartic_lb;quartic_bits.append(height)
   inclusion=False;transitions=[]
   for (prior_maps,H),record in zip(prior,m['history_comparison']['comparisons']):
    oldmap=prior_maps['rows'][i];assert oldmap['centre']==m['centre'] and list(map(QQ,oldmap['raw_coefficients']))==raw
    T=matrix(QQ,2,list(map(QQ,oldmap['matrix']))).inverse()*phi;flat=primitive(T.list());assert flat==record['transition'] and H==record['height'];T=matrix(ZZ,2,flat);transitions.append((T,H))
    # Primitive normalization can only decrease the resulting integer height.
    inclusion |= max(sum(abs(T[k,l]) for l in range(2)) for k in range(2))*p['height']<=H
   witness=m['history_comparison']['witness']
   if witness is not None:
    q=primitive(witness);assert q==witness and max(map(abs,q))<=p['height']
    for T,H in transitions:assert max(map(abs,primitive((T*matrix(ZZ,2,1,q)).list())))>H
    assert not inclusion;new+=1
   if inclusion:covered+=1
  rows.append(dict(id=row['id'],j_numerator_bits=nb,j_denominator_bits=db,integral_short_coefficient_bits_lower_bound=short_lb,integral_quartic_coefficient_bits_lower_bound=quartic_lb,observed_refined_quartic_bits_min=min(quartic_bits),observed_refined_quartic_bits_max=max(quartic_bits),proved_new_history_boxes=new,proved_covered_by_prior_boxes=covered))
  print(row['id'],'NEW',new,'COVERED',covered,'QUARTIC lower/observed',quartic_lb,min(quartic_bits),max(quartic_bits),flush=True)
 result=dict(schema='elliptic-curves.bounded-prime-history-height.v1',status='PASS',rows=rows,inputs={str(f.relative_to(ROOT)):hashed(f) for f in paths},source_sha256=hashed(Path(__file__)),inequalities=dict(short='For integral |A|,|B|<2^k, j=6912 A^3/(4A^3+27B^2) has numerator bits<=3k+13 and denominator bits<=3k+5 after cancellation.',quartic='For integral binary-quartic coefficients of absolute value<2^k, |I|<16*2^(2k), |J|<137*2^(3k). From j=6912 I^3/(4I^3-J^2), numerator bits<=6k+25 and denominator bits<=6k+16. Cancellation can only decrease both.'),boundary='Universal coefficient lower bounds for integral short models and integral binary quartics with this j, not rank bounds or bounds on point-search height. Exact new coordinate witnesses need not lift to curve points. Prior coverage is limited to the explicitly bound same-centre boxes.')
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
