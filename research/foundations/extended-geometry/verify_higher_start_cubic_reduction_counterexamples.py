#!/usr/bin/env python3
"""Exact coefficient, derivative, valuation and moment controls; standard library."""
import argparse,json,sys
from fractions import Fraction as Q
from hashlib import sha256
from math import comb,factorial
from pathlib import Path
if not __debug__:raise SystemExit("Run without -O: assertions are required.")
HERE=Path(__file__).resolve().parent

def trim(p):
 p=list(p)
 while p and not p[-1]:p.pop()
 return p
def add(p,q,s=1):
 r=[Q(0)]*max(len(p),len(q))
 for i,c in enumerate(p):r[i]+=c
 for i,c in enumerate(q):r[i]+=s*c
 return trim(r)
def scale(p,s):return trim([c*s for c in p])
def mul(p,q,limit=None):
 if not p or not q:return []
 n=len(p)+len(q)-1
 if limit is not None:n=min(n,limit+1)
 r=[Q(0)]*n
 for i,c in enumerate(p):
  for j,d in enumerate(q):
   if i+j<n:r[i+j]+=c*d
 return trim(r)
def power(p,n,limit=None):
 assert isinstance(n,int) and n>=0
 r=[Q(1)]
 while n:
  if n&1:r=mul(r,p,limit)
  p=mul(p,p,limit);n//=2
 return r
def derivative(p):return trim([i*p[i]for i in range(1,len(p))])
def compose(p,w):
 out=[]
 for c in reversed(p):out=add(mul(out,w),[c])
 return out
def divide(p,q):
 r=list(p);out=[Q(0)]*max(0,len(p)-len(q)+1)
 while r and len(r)>=len(q):
  k=len(r)-len(q);c=r[-1]/q[-1];out[k]=c
  r=add(r,[Q(0)]*k+scale(q,c),-1)
 return trim(out),trim(r)
def valuation(p,q):
 n=0
 while p:
  out,rem=divide(p,q)
  if rem:return n
  n+=1;p=out
 raise AssertionError("Valuation of zero is not a finite control")
def gcd(p,q):
 while q:_,r=divide(p,q);p,q=q,r
 return scale(p,1/p[-1])if p else []
def rational_derivative(N,D):return add(mul(derivative(N),D),mul(N,derivative(D)),-1),power(D,2)
def kappa(a,m):
 k=Q(2**m*factorial(m))
 for j in range(m+1):k/=a*m+2*j+1
 return k
def binom_integer(n,k):
 r=1
 for i in range(k):r=r*(n-i)//(i+1)
 return r
def raw(C,D,a,m,G):
 if a*m%2:return Q(0)
 out=Q(0)
 for i in range(m+1):
  P=mul(mul(G,power(C,m-i,m),m),power(D,i,m),m)
  out+=Q(comb(m,i),a*m+2*i+1)*(P[m]if len(P)>m else 0)
 return out
def lagrange_phi_one(m,j):
 if not j or m<j:return Q(0)
 n=m-j
 inv=[Q(binom_integer(-2*m,k))*Q(2,3)**k for k in range(n+1)]
 c=mul(power([1,1],3*m,n),inv,n)
 return Q(j,m)*(c[n]if len(c)>n else 0)
def mu(n,lam,r):
 out=lam**(-1)*(-n*n*lam*lam)**r
 for j in range(r):out*=Q(1,2*n)+j;out*=Q(1)-Q(1,2*n)+j;out/=Q(1,2)+j;out/=j+1
 return out
def digest(p):return sha256(p.read_bytes()).hexdigest()

def first_family(a,maxm):
 assert a>=3 and a%2==1
 U=[1,1];V=[1,Q(2,3)];X=mul(U,V)
 C=mul(power(V,a-2),power(U,a+3))
 D=scale(mul(power(V,a),power(U,a+5)),-1)
 assert not add(D,mul(C,power(X,2)))
 N=mul([0,1],power(V,2));Den=power(U,3)
 dN,dD=rational_derivative(N,Den)
 assert mul(power(X,a-1),dD)==mul(C,dN)
 assert mul(mul([0,1],power(X,a)),Den)==mul(C,N)
 assert mul(C,V)==mul(power(X,a-1),power(U,4))
 assert not divide(C,power(X,a-1))[1]==[]
 vals={}
 for name,L in [('U',U),('V',V)]:
  vals[name]={'P':valuation(X,L),'H':valuation(C,L),'D':valuation(D,L),
              'u_order':valuation(N,L)-valuation(Den,L),
              'omega_order':(a-1)*valuation(X,L)-valuation(C,L)}
 assert vals['V']=={'P':1,'H':a-2,'D':a,'u_order':2,'omega_order':1}
 assert vals['U']=={'P':1,'H':a+3,'D':a+5,'u_order':-3,'omega_order':-4}
 count=0;nonzero=[]
 for m in range(1,maxm+1):
  for j in range(4):
   actual=raw(C,D,a,m,[0]*j+[1])
   expected=Q(0)if m%2 else kappa(a,m)*lagrange_phi_one(m,j)
   assert actual==expected,(a,m,j,actual,expected)
   count+=1
   if j==1 and actual:nonzero.append(m)
 assert len(C)-1==2*a+1 and len(D)-1==2*a+5
 assert add(scale(Den,Q(-4,9)),N)==scale([4,3],Q(-1,9))
 return {'a':a,'degree':3*a+6,'C_degree':len(C)-1,'D_degree':len(D)-1,
         'geometric_valuations':vals,'raw_controls':count,'x_nonzero_tested_orders':nonzero,
         'C_coefficients':list(map(str,C)),'D_coefficients':list(map(str,D))}

def hyperbolic_pullback(a,n,maxm):
 assert a>=3 and a%2==1 and n>=1
 lam=Q(27,4);A=[1,1];B=[1,3];Qp=[1,Q(3,4)]
 w=mul([0,1],power(A,2));wp=derivative(w);Z0=add([1],scale(w,lam))
 assert wp==mul(A,B)and Z0==mul(power(B,2),Qp)
 z=[0,1];Zbar=[1,lam]
 Sbar=[Q(comb(2*n,j+1))*lam**j for j in range(2*n)]
 S=compose(Sbar,w)
 Xbar=scale(mul(Zbar,Sbar),Q(1,2*n))
 Cbar=scale(mul(power(Zbar,n+a),power(Sbar,a-1)),Q(1,2*n)**(a-1))
 Dbar=scale(mul(power(Zbar,n+a+2),power(Sbar,a+1)),-Q(1,2*n)**(a+1))
 X=scale(mul(mul(mul(A,B),Qp),S),Q(1,2*n))
 C=scale(mul(mul(mul(power(A,a-2),power(B,2*n+a)),power(Qp,n+a)),power(S,a-1)),Q(1,2*n)**(a-1))
 D=scale(mul(mul(mul(power(A,a),power(B,2*n+a+2)),power(Qp,n+a+2)),power(S,a+1)),-Q(1,2*n)**(a+1))
 assert not add(D,mul(C,power(X,2)))
 # Full coefficient reconstruction of the actual admissible source pullback.
 assert mul(C,power(wp,a))==mul(compose(Cbar,w),power(A,2*(a-1)))
 assert mul(D,power(wp,a+2))==mul(compose(Dbar,w),power(A,2*(a+1)))
 assert mul(mul(X,[0,1]),wp)==mul(w,compose(Xbar,w))
 # u=ubar(w) and its NONCONSTANT Jacobian T=2/(1+Z0^(2n)).
 N=scale(add(power(Z0,2*n),[1],-1),Q(1,2*n)/lam);Den=power(Z0,n)
 assert mul(mul([0,1],power(X,a)),Den)==mul(C,N)
 dN,dD=rational_derivative(N,Den);Tden=add([1],power(Z0,2*n))
 assert mul(mul(power(X,a-1),dD),Tden)==scale(mul(C,dN),2)
 assert len(Tden)>1
 assert mul(C,A)==mul(power(X,a-1),mul(power(B,2*n+1),power(Qp,n+1)))
 assert divide(C,power(X,a-1))[1]
 # Every root of S is genuinely simple and disjoint from A,B,Q.
 assert gcd(S,derivative(S))==[1]
 assert gcd(S,mul(mul(A,B),Qp))==[1]
 vals={}
 for name,L in [('A',A),('B',B),('Q',Qp),('S',S)]:
  vals[name]={'P':valuation(X,L),'H':valuation(C,L),'D':valuation(D,L),
              'u_order':valuation(N,L)-valuation(Den,L),
              'omega_order':(a-1)*valuation(X,L)-valuation(C,L)}
 assert vals['A']=={'P':1,'H':a-2,'D':a,'u_order':2,'omega_order':1}
 assert vals['B']['H']==2*n+a and vals['Q']['H']==n+a and vals['S']['H']==a-1
 count=0;marked=[]
 for m in range(1,maxm+1):
  pure=raw(C,D,a,m,[1]);assert pure==0
  actual=raw(C,D,a,m,w);target=raw(Cbar,Dbar,a,m,[0,1])
  expected=Q(0)if m%2 else kappa(a,m)*mu(n,lam,m//2)
  assert actual==target==expected,(a,n,m,actual,target,expected)
  count+=3
  if actual:marked.append({'m':m,'value':str(actual)})
 c=3*(2*a-1)*n+1;e=3*(2*a+3)*n+1;degree=e+a+1
 assert len(C)-1==c and len(D)-1==e and len(X)-1==6*n
 return {'a':a,'n':n,'degree':degree,'C_degree':c,'D_degree':e,'X_degree':len(X)-1,
         'nonconstant_T':True,'geometric_valuations':vals,'raw_controls':count,
         'nonzero_even_marker_w_controls':marked,
         'C_coefficients':list(map(str,C)),'D_coefficients':list(map(str,D))}

# Sparse original polynomial operations, independent of angular integration.
def sa(P,Q):
 R=dict(P)
 for e,c in Q.items():R[e]=R.get(e,0)+c
 return {e:c for e,c in R.items()if c}
def ss(P,c):return {e:v*c for e,v in P.items()if v*c}
def sm(P,Q):
 R={}
 for e,c in P.items():
  for f,d in Q.items():
   k=tuple(a+b for a,b in zip(e,f));R[k]=R.get(k,0)+c*d
 return {e:c for e,c in R.items()if c}
def sp(P,n):
 R={(0,0,0,0):Q(1)}
 while n:
  if n&1:R=sm(R,P)
  P=sm(P,P);n//=2
 return R
def image_scalar(P,j=0):
 return sum(c*factorial(e[0]+j)*factorial(e[1])for e,c in P.items()if e[0]+j==e[2]and e[1]==e[3]+j)
def original_control():
 x1,x2,z1,z2=[{tuple(int(k==i)for k in range(4)):Q(1)}for i in range(4)]
 R=sa(sm(x1,z1),sm(x2,z2));Z=sm(x1,z2);W=ss(sm(x2,z1),2);T=sa(sm(x1,z1),ss(sm(x2,z2),-1))
 U=sa(R,Z);V=sa(R,ss(Z,Q(2,3)));P=sm(U,V);C=sm(V,sp(U,6))
 sq=sa(sp(P,2),ss(sp(R,4),-1))
 quotient={}
 for e,c in sq.items():
  # Symbolic division by Z=xi1*z2 is literal termwise in these coordinates.
  assert e[0]>0 and e[3]>0
  f=(e[0]-1,e[1],e[2],e[3]-1);quotient[f]=c
 assert sm(Z,quotient)==sq
 B=sa(ss(sm(W,sp(R,4)),2),ss(sm(sp(T,2),quotient),-1));F=sm(sm(sp(T,3),C),B);d=15
 assert all(e[0]+e[1]==e[2]+e[3]==d for e in F)
 F2=sm(F,F)
 assert image_scalar(F)==image_scalar(F2)==0
 marked=image_scalar(F2,1);expected=factorial(2*d+2)*kappa(3,2)*Q(5,3)
 assert marked==expected!=0
 assert image_scalar(F2,2)==factorial(2*d+3)*kappa(3,2)!=0
 assert marked!=factorial(2*d+1)*kappa(3,2)*Q(5,3)
 radial=sum(c*((-1)**e[3])for e,c in F.items())
 assert radial!=0
 return {'a':3,'degree':d,'monomials':len(F),'F_squared_monomials':len(F2),
         'radial_point_value':str(radial),'Z_F_squared_normalized':str(Q(marked,factorial(2*d+2))),
         'wrong_factorial_rejected':True}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--repo-root',type=Path,default=Path.cwd())
 parser.add_argument('--output',type=Path,default=HERE/'higher_start_cubic_reduction_counterexamples_verification.json')
 args=parser.parse_args()
 ip=HERE/'higher_start_cubic_reduction_counterexamples_inputs.json'
 lp=HERE/'higher_start_cubic_reduction_counterexamples_source_ledger.json'
 data=json.loads(ip.read_text(encoding='utf-8-sig'));ledger=json.loads(lp.read_text(encoding='utf-8-sig'))
 pins=[]
 for src in ledger['sources']:
  p=args.repo_root/src['path'];assert digest(p)==src['sha256'];pins.append({'path':src['path'],'sha256':digest(p),'role':src['role']})
 first=[first_family(a,data['maximum_moment'])for a in data['phi_one_starting_heights']]
 second=[hyperbolic_pullback(a,n,data['maximum_moment'])for a,n in data['hyperbolic_parameters']]
 assert any(r['degree']==59 for r in second)
 original=original_control()
 # A genuine geometric valuation contradiction to universal divisibility.
 assert all(r['geometric_valuations']['V']['H']<(r['a']-1)*r['geometric_valuations']['V']['P']for r in first)
 names=['HIGHER_START_CUBIC_REDUCTION_COUNTEREXAMPLES.md',Path(__file__).name,ip.name,lp.name]
 receipt={'status':'PASS','scope':'Counterexamples over every characteristic-zero field; all odd a>=3. Finite controls supplement the all-order endpoint/pullback proofs, not a moment cutoff.',
 'phi_one_family':first,'nonconstant_phi_pullbacks':second,'original_control':original,
 'raw_control_count':sum(r['raw_controls']for r in first+second),'source_pins':pins,
 'packet_sha256':{name:digest(HERE/name)for name in names},
 'runtime':{'python':sys.version.split()[0],'standard_library_only':True},
 'assurance':'Internal exact algebra; no formal verification, external review or prime-degree priority claim.'}
 args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':'PASS','raw_controls':receipt['raw_control_count'],'prime59':True,'output':str(args.output)},sort_keys=True))
if __name__=='__main__':main()

