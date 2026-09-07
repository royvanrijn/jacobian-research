#!/usr/bin/env sage-python
"""Two fixed good primes on the first retained Kihara K3 parent."""
import argparse,sys,json,hashlib
from pathlib import Path
import numpy as np
from sage.all import QQ,ZZ,AA,GF,PolynomialRing,EllipticCurve,prime_range,cyclotomic_polynomial
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/kihara-picard-count-v1'
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint
SOURCE=ART/'kihara_section_involution_v1.json'
def good_model(A,B,p):
 try:
  R=PolynomialRing(GF(p),'T');a,b=R(A),R(B)
 except (ValueError,ZeroDivisionError):return None
 delta=-16*(4*a**3+27*b**2)
 if (a.degree(),b.degree(),delta.degree())!=(8,12,20) or delta.gcd(delta.derivative())!=1 or delta.gcd(a)!=1 or not delta(0):return None
 node=-3*b[12]/(2*a[8])
 if 3*node**2+a[8] or node**3+a[8]*node+b[12] or not (3*node).is_square():return None
 return a,b,delta,node
def legendre_table(p):
    chars=np.full(p,-1,dtype=np.int64);chars[0]=0
    for a in range(1,p):chars[a*a%p]=1
    return chars

def count(A,B,p,extension,progress):
    aa,bb,delta,_=good_model(A,B,p);chars=legendre_table(p)
    d=next(n for n in range(2,p) if chars[n]==-1)
    if extension==1:
        F=GF(p);q=p;x0=np.arange(p,dtype=np.int64);x1=None;z=None
    else:
        R=PolynomialRing(GF(p),'z');F=GF(p*p,'z',modulus=R.gen()**2-d);z=F.gen();q=p*p
        index=np.arange(q,dtype=np.int64);x0=index%p;x1=index//p
        cube0=(x0*x0*x0+3*d*x0*x1*x1)%p;cube1=(3*x0*x0*x1+d*x1*x1*x1)%p
    R=PolynomialRing(F,'T');a=R(aa);b=R(bb);disc=-16*(4*a**3+27*b**2)
    def pair_eval(coeffs,t0,t1):
        r0=r1=0
        for c in reversed(coeffs):r0,r1=(r0*t0+d*r1*t1+int(c))%p,(r0*t1+r1*t0)%p
        return r0,r1
    counts=[];repairs=[];smooth=0;count_total=0
    for index in range(q+1):
        if index==q:
            t=None;a0,a1=int(aa[8]),0;b0,b1=int(bb[12]),0;av,bv=F(a0),F(b0)
        else:
            t0=index%p;t1=0 if extension==1 else index//p
            t=F(t0) if extension==1 else F(t0)+F(t1)*z
            a0,a1=pair_eval(aa.list(),t0,t1);b0,b1=pair_eval(bb.list(),t0,t1)
            av,bv=a(t),b(t)
            if extension==1:
                if av!=F(a0) or bv!=F(b0) or a1 or b1:raise ArithmeticError('prime-field coefficient evaluation differs')
            elif av!=F(a0)+F(a1)*z or bv!=F(b0)+F(b1)*z:raise ArithmeticError('quadratic-field coefficient evaluation differs')
        if extension==1:
            values=(x0*x0*x0+a0*x0+b0)%p;naive=q+1+int(chars[values].sum())
        else:
            f0=(cube0+a0*x0+d*a1*x1+b0)%p;f1=(cube1+a0*x1+a1*x0+b1)%p
            norm=(f0*f0-d*f1*f1)%p;naive=q+1+int(chars[norm].sum())
        if 4*av**3+27*bv**2:
            independent=int(EllipticCurve(F,[av,bv]).cardinality(algorithm='pari'))
            if independent!=naive:raise ArithmeticError('independent PARI cardinality differs from norm-character count')
            smooth+=1;correction=0
        elif index==q:correction=3*q
        else:
            multiplicity=0;temp=disc;T=R.gen()
            while temp(t)==0:temp=temp//(T-t);multiplicity+=1
            if multiplicity not in (1,2):raise ArithmeticError('semistable reduction type changed')
            correction=(multiplicity-1)*q
        if correction:repairs.append({'base_index':index,'correction':correction})
        counts.append(naive);count_total+=naive+correction
        if (index+1)%256==0:progress(index+1,q+1)
    if repairs!=[{'base_index':q,'correction':3*q}]:raise ArithmeticError('only split I4 repair expected')
    return {'prime':p,'extension_degree':extension,'field_order':q,'quadratic_nonsquare':None if extension==1 else d,
      'base_order':'index=a+p*b represents a+b*z; infinity has index q',
      'weierstrass_fibre_counts':counts,'resolution_corrections':repairs,'surface_point_count':count_total,
      'smooth_fibres_independently_counted':smooth,'all_smooth_counts_verified_by':'PARI via Sage elliptic cardinality'}

def quotient_count(A,B,p):
 F=GF(p);R=PolynomialRing(F,'s');a=R([A[2*i] for i in range(5)]);b=R([B[2*i] for i in range(7)]);d=-16*(4*a**3+27*b**2)
 assert (a.degree(),b.degree(),d.degree())==(4,6,10) and d.gcd(d.derivative())==1
 chars=legendre_table(p);values=np.arange(p,dtype=np.int64);counts=[]
 for i in range(p+1):
  av,bv=(a[4],b[6]) if i==p else (a(F(i)),b(F(i)))
  n=p+1+int(chars[(values**3+int(av)*values+int(bv))%p].sum())
  if 4*av**3+27*bv**2:assert n==EllipticCurve(F,[av,bv]).cardinality()
  counts.append(n)
 total=sum(counts)+p;sign=(total-1-p*p-9*p)/p;assert sign in [-1,1]
 return {'weierstrass_fibre_counts':counts,'infinity_correction':p,'surface_point_count':total,'remaining_geometric_quotient_divisor_sign':int(sign)}
def frobenius(p,n1,n2,eps):
 R=PolynomialRing(QQ,'z');z=R.gen();a=ZZ(n1-1-p*p-(17+eps)*p);s2=ZZ(n2-1-p**4-18*p*p);b=QQ(a*a-s2)/2;assert b.denominator()==1
 candidates=[]
 for sign in [1,-1]:
  if sign==-1 and b:continue
  f=z**4-a*z**3+b*z*z-sign*p*p*a*z+sign*p**4
  if sign==1:
   roots=(z*z-QQ(a)/p*z+QQ(b)/(p*p)-2).roots(AA)
   if sum(m for x,m in roots)!=2 or any(abs(x)>2 for x,m in roots):continue
  elif abs(a)>2*p:continue
  remaining=R(f(p*z)/p**4);cycles=[]
  for n in [1,2,3,4,5,6,8,10,12]:
   phi=R(cyclotomic_polynomial(n));m=0
   while remaining%phi==0:remaining//=phi;m+=1
   if m:cycles.append({'order':n,'multiplicity':m,'degree':int(phi.degree())})
  candidates.append({'determinant_sign':sign,'residual_polynomial':list(map(str,f.list())),'cyclotomic_factors':cycles,
    'geometric_NS_upper_bound':18+sum(c['degree']*c['multiplicity'] for c in cycles),
    'rational_NS_upper_bound':17+int(eps==1)+sum(c['multiplicity'] for c in cycles if c['order']==1)})
 assert candidates
 return {'known_rational_NS_rank':17,'extra_quotient_eigenvalue_sign':eps,'residual_trace1':int(a),'residual_trace2':int(s2),'candidates':candidates,'geometric_NS_upper_bound':max(c['geometric_NS_upper_bound'] for c in candidates),'rational_NS_upper_bound':max(c['rational_NS_upper_bound'] for c in candidates)}
def main():
 p=json.loads((D/'protocol.json').read_text())
 for name,h in p['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
 source=json.loads(SOURCE.read_text());R=PolynomialRing(QQ,'T');A=R(source['raw_A']);B=R(source['raw_B']);primes=[int(q) for q in prime_range(31,252) if good_model(A,B,q) is not None][:2];assert len(primes)==2
 out=D/'counts.json';assert not out.exists();result={'status':'RUNNING','prime_selection':'first two good primes in31..251','primes':primes,'rows':[]};checkpoint(out,result)
 for prime in primes:
  row={'prime':prime,'quotient':quotient_count(A,B,prime),'counts':[]};result['rows'].append(row);checkpoint(out,result)
  for k in [1,2]:
   row['counts'].append(count(A,B,prime,k,lambda i,n:print(prime,k,i,'of',n,flush=True)));checkpoint(out,result)
  row['frobenius']=frobenius(prime,*[c['surface_point_count'] for c in row['counts']],row['quotient']['remaining_geometric_quotient_divisor_sign']);checkpoint(out,result)
  print('PRIME',prime,'FROBENIUS',row['frobenius'],flush=True)
 result['status']='PASS';checkpoint(out,result)
 print('PASS two fixed-prime Kihara Picard witnesses',flush=True)
if __name__=='__main__':main()
