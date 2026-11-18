#!/usr/bin/env python3
"""Complete degree26 reciprocal census over F3 and mod9 obstruction.

Exact enumeration covers both unit lambda classes and every normalized
monic numerator satisfying the linear zero-branch/reciprocity conditions.
The note supplies the reduction and characteristic-zero consequence.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import argparse,hashlib,json,platform
import sympy as sp

P=3
def trim(a):
 a=[x%P for x in a]
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def add(a,b,sign=1):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=sign*x
 return trim(c)
def rem(a,m):
 a=trim(a[:]);inv=pow(m[-1],-1,P)
 while len(a)>=len(m) and a!=[0]:
  j=len(a)-len(m);q=a[-1]*inv%P
  for k,x in enumerate(m):a[j+k]=(a[j+k]-q*x)%P
  a=trim(a)
 return a
def mul(a,b,m):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return rem(c,m)
def gcd(a,b):
 while b!=[0]:a,b=b,rem(a,b)
 return [x*pow(a[-1],-1,P)%P for x in a]
def substitute(a,b,m):
 out=[0]
 for x in reversed(a):out=add(mul(out,b,m),[x])
 return rem(out,m)
def affine(rows,nvars):
 rows=[[x%P for x in row]for row in rows];pivots=[]
 for j in range(nvars):
  k=next((k for k in range(len(pivots),len(rows))if rows[k][j]),None)
  if k is None:continue
  i=len(pivots);rows[i],rows[k]=rows[k],rows[i]
  scale=pow(rows[i][j],-1,P);rows[i]=[(scale*x)%P for x in rows[i]]
  for k in range(len(rows)):
   if k==i:continue
   scale=rows[k][j]
   rows[k]=[(a-scale*b)%P for a,b in zip(rows[k],rows[i])]
  pivots.append(j)
 assert all(any(row[:nvars])or not row[-1]for row in rows)
 free=[i for i in range(nvars)if i not in pivots]
 return rows,pivots,free
def rank(rows):
 if not rows:return 0
 _,pivots,_=affine([row+[0]for row in rows],len(rows[0]))
 return len(pivots)
def enumerate_models(lam):
 eq=[[1]*13+[2],list(range(13))+[(-13)%3]]
 for j in range(6):
  c=pow(lam,6-j,3)
  eq.append([((2*k-1)*((k>=j+1)-c*(k>=13-j)))%3 for k in range(13)]+[(c-1)%3])
 rows,pivots,free=affine(eq,13)
 counts=Counter();models=[]
 for values in product(range(3),repeat=len(free)):
  N=[0]*14;N[13]=1
  for i,value in zip(free,values):N[i]=value
  for row,j in zip(rows,pivots):N[j]=(row[-1]-sum(row[i]*N[i]for i in free))%3
  S=[0]*13;partial=0
  for j in reversed(range(1,14)):
   partial=(partial+(2*j-1)*N[j])%3;S[j-1]=partial
  counts["linear_candidates"]+=1
  if sum(S)%3==0 or gcd(S,[(-lam)%3,0,1])!=[1]:continue
  derivative=trim([(i*S[i])%3 for i in range(1,13)])
  if gcd(S,derivative)!=[1]or gcd(N,S)!=[1]:continue
  counts["nondegenerate"]+=1
  invu=[(-S[i])%3 for i in range(1,13)]
  M=mul(mul(rem(N,S),rem(N,S),S),invu,S)
  if substitute(M,[(lam*x)%3 for x in invu],S)!=M:continue
  counts["invariant"]+=1
  powers=[];current=[1]
  for _ in range(7):
   powers.append(current+[0]*(12-len(current)));current=mul(current,M,S)
  if rank(powers)!=6:continue
  counts["models"]+=1;models.append((N[:],S[:],M[:]))
 return {"counts":dict(counts),"free_indices":free,"linear_rank":len(pivots)},models

def constraints(values,mod=9):
 N=values[:13]+[1];rho=values[13]%mod
 S=[0]*13;partial=0
 for j in reversed(range(1,14)):
  partial=(partial+(2*j-1)*N[j])%mod;S[j-1]=partial
 out=[sum(N)%mod,sum(j*c for j,c in enumerate(N))%mod]
 out.extend((S[j]-pow(rho,12-2*j,mod)*S[12-j])%mod for j in range(6))
 def product_raw(a,b):
  c=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%mod
  return c
 reflected=[N[13-i]*pow(rho,2*(13-i),mod)%mod for i in range(14)]
 first=[0]*24+[rho*rho*x%mod for x in product_raw(N,N)]
 second=product_raw(reflected,reflected)
 E=[(first[i]-(second[i] if i<len(second)else 0))%mod for i in range(len(first))]
 inverse=pow(S[-1],-1,mod)
 for j in reversed(range(12,len(E))):
  quotient=E[j]*inverse%mod
  for k,c in enumerate(S):E[j-12+k]=(E[j-12+k]-quotient*c)%mod
 return out+E[:12]

def nonlift_certificate(N):
 point=N[:13]+[1]
 residual=constraints(point)
 assert all(x%3==0 for x in residual)
 columns=[]
 for j in range(14):
  shifted=point[:];shifted[j]+=3
  derivative=[((a-b)%9)//3 for a,b in zip(constraints(shifted),residual)]
  columns.append(derivative)
 rhs=[(-x//3)%3 for x in residual]
 assert rank(columns)==13 and rank(columns+[rhs])==14
 equations=[column+[0]for column in columns]+[rhs+[1]]
 rows,pivots,free=affine(equations,20)
 witness=[0]*20
 for row,j in zip(rows,pivots):witness[j]=row[-1]
 assert all(sum(a*b for a,b in zip(witness,column))%3==0 for column in columns)
 assert sum(a*b for a,b in zip(witness,rhs))%3==1
 return {"rank":13,"augmented_rank":14,"rhs":rhs,"left_kernel_witness":witness}

def polynomial_certificates(model):
 N,S,M=model
 u,y,z,s=sp.symbols("u y z s")
 Nu=sum(c*u**i for i,c in enumerate(N));Su=sum(c*u**i for i,c in enumerate(S))
 assert sp.Poly(2*u*sp.diff(Nu,u)-Nu-(u-1)*Su,u,modulus=3).is_zero
 assert sp.Poly(Nu,u,modulus=3).eval(1)==0
 assert sp.Poly(sp.diff(Nu,u),u,modulus=3).eval(1)==0
 invu=[(-S[i])%3 for i in range(1,13)]
 reciprocalN=substitute(N,invu,S)
 plus=gcd(S,add(rem(N,S),mul([0,1],reciprocalN,S)))
 minus=gcd(S,add(rem(N,S),mul([0,1],reciprocalN,S),-1))
 assert (len(plus)-1,len(minus)-1)in((8,4),(6,6))
 columns=[]
 for i in range(12):
  column=mul(M,[0]*i+[1],S);columns.append(column+[0]*(12-len(column)))
 matrix=sp.Matrix(12,12,lambda i,j:columns[j][i])
 chi=sp.Poly(matrix.charpoly(y).as_expr(),y,modulus=3)
 Q=sp.Poly(y**6,y,modulus=3)
 for k in range(1,7):
  coefficient=(int(chi.nth(12-k))-int((Q*Q).nth(12-k)))*2%3
  Q+=sp.Poly(coefficient*y**(6-k),y,modulus=3)
 assert Q.degree()==6 and Q.eval(0)!=0
 assert sp.gcd(Q,Q.diff()).degree()==0 and Q*Q==chi
 H=Nu.subs(u,z*z)-1
 D=sp.Poly(H-s*z+1,z,domain=sp.GF(3).poly_ring(s)).discriminant()
 assert sp.Poly(sp.sympify(str(D))-s*s*Q.as_expr().subs(y,s*s)**2,s,modulus=3).is_zero
 return {"N_coefficients_ascending":N,"S_coefficients_ascending":S,
         "M_coefficients_ascending":M,"mixed_sign_degrees":[len(plus)-1,len(minus)-1],
         "Q6":str(Q.as_expr()),"full_discriminant":str(D),
         "nonlift":nonlift_certificate(N)}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path)
 args=parser.parse_args()
 assert sp.__version__=="1.14.0"
 census={};models={}
 for lam in(1,2):
  census[str(lam)],models[lam]=enumerate_models(lam)
 assert census["1"]["counts"]=={"linear_candidates":2187,"nondegenerate":728,"invariant":2,"models":2}
 assert census["2"]["counts"]=={"linear_candidates":243}
 assert len(models[1])==2 and not models[2]
 certificates=[polynomial_certificates(model)for model in models[1]]
 expected=[1,1,1,1,0,0,1,1,0,0,1,0,1,1]
 assert models[1][0][0]==expected
 assert certificates[0]["Q6"]=="y**6 - y**5 + y**4 - y**3 - y**2 - 1"
 report={
  "schema_version":1,"kind":"complete_finite_reciprocal_degree26_census",
  "python":platform.python_version(),"sympy":sp.__version__,
  "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  "arithmetic":"Complete F3 enumeration; Z/9 polynomial constraints; exact finite-field discriminants",
  "census":census,"models":certificates,
  "scope":"Finite F3 feasibility and no mod9 lift with good paired reciprocal reduction. The written proof excludes normalized integral reciprocal models with good paired reduction at3, without excluding all Q or Q3 models.",
  "independent_replay":False,"formal_verification":False,
  "external_review":False,"result":"PASS"
 }
 output=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if args.write_receipt:args.write_receipt.write_text(output)
 print(output,end="")
if __name__=="__main__":main()

