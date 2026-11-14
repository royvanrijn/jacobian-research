"""Exact cofactor coverage and two-engine finite-free chart unit replay.
No symbolic resultant, Groebner search or factorization is called.
"""
from pathlib import Path
import hashlib,json,math,time,sympy as S
from flint import fmpz_poly,fmpq_poly,fmpz_mat,nmod_poly
from build_remaining_charts import Algebra
HERE=Path(__file__).resolve().parent;START=time.monotonic();P=1000000007;checks=0
if not __debug__:raise RuntimeError('Run without -O')
def yes(c):
 global checks
 assert c;checks+=1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(stage,**kw):print(json.dumps({'stage':stage,'seconds':time.monotonic()-START,**kw}),flush=True)
def bareiss(M):
 W=[r[:] for r in M];prev=1;sign=1
 for k in range(len(W)-1):
  j=next((i for i in range(k,len(W)) if W[i][k]),None)
  if j is None:return 0
  if j!=k:W[k],W[j]=W[j],W[k];sign=-sign
  v=W[k][k]
  for i in range(k+1,len(W)):
   z=W[i][k]
   for j in range(k+1,len(W)):
    u,r=divmod(W[i][j]*v-z*W[k][j],prev);assert not r;W[i][j]=u
   W[i][k]=0
  prev=v
 return sign*W[-1][-1]
def sylv(f,g):
 m,n=len(f)-1,len(g)-1;M=[[0]*(m+n) for _ in range(m+n)]
 for i in range(n):
  for j,c in enumerate(reversed(f)):M[i][i+j]=c
 for i in range(m):
  for j,c in enumerate(reversed(g)):M[n+i][i+j]=c
 return M
class Plain:
 def __init__(self,h,f):
  self.h=list(h);v=pow(self.h[-1],-1,P);self.h=[c*v%P for c in self.h];self.n=len(f)-1;self.f=[tuple(c) for c in f];self.zero=tuple([()]*self.n);self.one=tuple([(1,)]+[()]*(self.n-1))
 def br(self,v):
  v=[c%P for c in v]
  while v and not v[-1]:v.pop()
  while len(v)>=len(self.h):
   k=len(v)-len(self.h);c=v[-1]
   for j,d in enumerate(self.h):v[k+j]=(v[k+j]-c*d)%P
   while v and not v[-1]:v.pop()
  return tuple(v)
 def ba(self,u,v,sg=1):return self.br([(u[i] if i<len(u) else 0)+sg*(v[i] if i<len(v) else 0) for i in range(max(len(u),len(v)))])
 def bm(self,u,v):
  if not u or not v:return ()
  out=[0]*(len(u)+len(v)-1)
  for i,c in enumerate(u):
   for j,d in enumerate(v):out[i+j]+=c*d
  return self.br(out)
 def add(self,u,v):return tuple(self.ba(c,d) for c,d in zip(u,v))
 def sub(self,u,v):return tuple(self.ba(c,d,-1) for c,d in zip(u,v))
 def mul(self,u,v):
  out=[()]*(2*self.n-1)
  for i,c in enumerate(u):
   for j,d in enumerate(v):out[i+j]=self.ba(out[i+j],self.bm(c,d))
  for k in range(len(out)-1,self.n-1,-1):
   for j in range(self.n):out[k-self.n+j]=self.ba(out[k-self.n+j],self.bm(out[k],self.f[j]),-1)
  return tuple(out[:self.n])
 def const(self,c):return tuple([self.br([c])]+[()]*(self.n-1))
 def decode(self,u):
  yes(len(u)==self.n and all(len(c)<len(self.h) and all(0<=v<P for v in c) and (not c or c[-1]) for c in u));return tuple(tuple(c) for c in u)
 def evaluate(self,expr,variables,values):
  terms=S.Poly(expr,*variables).terms();powers=[[self.one] for _ in values]
  for j,value in enumerate(values):
   for _ in range(max(ex[j] for ex,c in terms)):powers[j].append(self.mul(powers[j][-1],value))
  out=self.zero
  for ex,c in terms:
   term=self.const(int(c))
   for j,k in enumerate(ex):term=self.mul(term,powers[j][k])
   out=self.add(out,term)
  return out
 def tinverse(self,q,derivative,inverse):
  out=[self.zero]*(len(derivative)+len(inverse)-1)
  for i,c in enumerate(derivative):
   for j,d in enumerate(inverse):out[i+j]=self.add(out[i+j],self.mul(c,d))
  for k in range(len(out)-1,len(q)-2,-1):
   for j,c in enumerate(q[:-1]):out[k-len(q)+1+j]=self.sub(out[k-len(q)+1+j],self.mul(out[k],c))
   out[k]=self.zero
  yes(out[0]==self.one and all(v==self.zero for v in out[1:]))
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());src=json.loads((HERE/'quadruple-input.json').read_text());cp=json.loads((HERE/'coordinate-poles-candidate.json').read_text());fac=json.loads((HERE/'coordinate-pole-resultant-candidate.json').read_text());charts=json.loads((HERE/'remaining-charts-candidate.json').read_text());s,a,b=S.symbols('s a b');L=S.Poly(S.sympify(old['linear_equations'][0]),b);PP,QQ=[S.Poly(v,a,s) for v in [L.nth(1),L.nth(0)]];G=S.Poly(S.sympify(old['equations'][3]['polynomial']),a,s);f2=S.sympify(src['equations'][0]['polynomial'])
yes(cp['source_sha256']==sha(HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json'));yes(charts['prime']==P and all(P%d for d in range(2,math.isqrt(P)+1)))
U,H,num,den,rnum,rden,common=[fmpz_poly(list(map(int,cp[k]))) for k in ['resultant_coefficients','projection_coefficients','a_numerator','a_denominator','raw_numerator','raw_denominator','common_factor']];yes(rnum==num*common and rden==den*common and H.degree()==87)
product=fmpz_poly([int(fac['factor_content'])])
for row in fac['factors']:product*=fmpz_poly(list(map(int,row['coefficients'])))**row['multiplicity']
yes(product==U and U==fmpz_poly(list(map(int,fac['coefficients']))));yes(H==fmpz_poly(list(map(int,next(r['coefficients'] for r in fac['factors'] if len(r['coefficients'])==88)))))
# Full actual-degree8 Sylvester matrix and cofactor identities.
def coefs(F):return [fmpz_poly([int(F.coeff_monomial(a**j*s**k)) for k in range(F.degree(s)+1)]) for j in range(F.degree(a)+1)]
pp,qq=coefs(PP),coefs(QQ);yes([PP.degree(a),QQ.degree(a),PP.degree(s),QQ.degree(s)]==[3,5,14,16]);yes(cp['resultant_degree_bound']==118 and cp['cofactor_degree_bound']==102)
for sv in range(-59,60):
 M=sylv([int(c(sv)) for c in pp],[int(c(sv)) for c in qq]);yes(len(M)==8);value=int(U(sv));yes(int(fmpz_mat(M).det())==value);yes(bareiss(M)==value)
for sv in range(-51,52):
 M=sylv([int(c(sv)) for c in pp],[int(c(sv)) for c in qq])
 for col,cofactor in [(6,rnum),(7,rden)]:
  W=[[c for j,c in enumerate(r) if j!=col] for i,r in enumerate(M) if i!=5];value=int(cofactor(sv));sign=(-1)**(5+col);yes(sign*int(fmpz_mat(W).det())==value);yes(sign*bareiss(W)==value)
yes(U.degree()<=118 and rnum.degree()<=102 and rden.degree()<=102)
# Positive characteristic-zero identities for every retained base root.
Hq=fmpq_poly(list(map(int,H.coeffs())))
for row in [pp,qq]:
 degree=len(row)-1;value=sum((c*num**j*den**(degree-j) for j,c in enumerate(row)),fmpz_poly([]));yes(not (fmpq_poly(list(map(int,value.coeffs())))%Hq))
# Base units are replayed independently using coefficient-list arithmetic.
hp=nmod_poly([int(c)%P for c in H.coeffs()],P);yes(hp.degree()==87);baseplain=Plain(list(map(int,hp.coeffs())),[[],[1]])
for name,value in [('projection_derivative',H.derivative()),('raw_cofactor_denominator',rden),('reduced_denominator',den),('D',fmpz_poly([0,-2,0,0,3])),('ell_poles',fmpz_poly([-1,0,0,3])*fmpz_poly([1,0,0,3]))]:
 v=nmod_poly([int(c)%P for c in value.coeffs()],P);g,u,w=v.xgcd(hp);yes(g.degree()==0 and int(g[0])!=0);u=u*pow(int(g[0]),-1,P)%hp;yes(v*u%hp==1);yes(baseplain.bm(tuple(map(int,v.coeffs())),tuple(map(int,u.coeffs())))==(1,))
# Every leading-pole curve really has a-degree7 over the cubic base.
for htext,leading in [(3*s**3-1,S.Rational(-38440,9)),(3*s**3+1,S.Integer(-10728))]:
 remainder=S.Poly(S.rem(S.Poly(G.as_expr(),s,a),S.Poly(htext,s,a)).as_expr(),a);yes(remainder.degree()==7 and S.expand(remainder.LC()-leading)==0);yes(int(leading.p)%P!=0 and int(leading.q)%P!=0)
results=[]
for row in charts['charts']:
 name=row['name'];h=nmod_poly(row['base_modulus'],P);A=Algebra(h,[nmod_poly(c,P) for c in row['extension_modulus']]);B=Plain(row['base_modulus'],row['extension_modulus']);yes(A.f[-1]==1 and row['rank']==int(h.degree())*A.n)
 sv=A.base(nmod_poly([0,1],P));sp=B.decode(A.encode(sv));av=tuple(nmod_poly(c,P) for c in row['a_mod']);bv=tuple(nmod_poly(c,P) for c in row['b_mod']);ap=B.decode(row['a_mod']);bp=B.decode(row['b_mod'])
 if name=='coordinate-poles':
  yes(h==hp and A.n==2);nn=nmod_poly([int(c)%P for c in num.coeffs()],P);dd=nmod_poly([int(c)%P for c in den.coeffs()],P);yes(av==A.base(nn*A.binv(dd)%h));yes(bv==tuple([A.bz,A.bo]));fp=S.Poly(f2,b);co=[A.evaluate(fp.nth(j),[s,a],[sv,av])[0] for j in range(3)];iv=A.binv(co[-1]);yes(A.f==[c*iv%h for c in co])
 else:
  yes(name in ['leading-minus','leading-plus'] and A.n==7 and h==nmod_poly([-1 if name=='leading-minus' else 1,0,0,3],P));yes(av==tuple([A.bz,A.bo]+[A.bz]*5));co=[A.evaluate(S.Poly(G.as_expr(),a).nth(j),[s],[sv])[0] for j in range(9)];yes(not co[-1]);iv=A.binv(co[7]);yes(A.f==[c*iv%h for c in co[:8]]);pv=A.evaluate(L.nth(1),[s,a],[sv,av]);qv=A.evaluate(L.nth(0),[s,a],[sv,av]);pi=tuple(nmod_poly(c,P) for c in row['b_denominator_inverse_mod']);yes(A.mul(pv,pi)==A.one and A.add(A.mul(pv,bv),qv)==A.zero);pvplain=B.evaluate(L.nth(1),[s,a],[sp,ap]);yes(B.mul(pvplain,B.decode(row['b_denominator_inverse_mod']))==B.one);yes(B.add(B.mul(pvplain,bp),B.evaluate(L.nth(0),[s,a],[sp,ap]))==B.zero)
 # A is a successive monic finite-free presentation of actual degrees87*2 or3*7.
 defining_derivative=tuple(A.f[j+1]*(j+1) for j in range(A.n));fi=tuple(nmod_poly(c,P) for c in row['extension_derivative_inverse_mod']);yes(A.mul(defining_derivative,fi)==A.one);yes(B.mul(B.decode(A.encode(defining_derivative)),B.decode(row['extension_derivative_inverse_mod']))==B.one)
 quotient=[];qp=[]
 for text,encoded in zip(src['quotient_coefficients'],row['quotient_coefficients_mod']):
  expression=S.sympify(text['polynomial']);u=A.evaluate(expression,[s,a,b],[sv,av,bv]);yes(A.encode(u)==encoded);v=B.evaluate(expression,[s,a,b],[sp,ap,bp]);yes(v==B.decode(encoded));quotient.append(u);qp.append(v)
 yes(quotient[-1]==A.one and qp[-1]==B.one)
 inverse=[B.decode(c) for c in row['derivative_inverse_mod']];derivative=[B.mul(B.const(j+1),qp[j+1]) for j in range(8)];B.tinverse(qp,derivative,inverse)
 # Native replay of the same finite polynomial inverse.
 ni=[tuple(nmod_poly(c,P) for c in v) for v in row['derivative_inverse_mod']];nd=[A.mul(A.const(j+1),quotient[j+1]) for j in range(8)];out=[A.zero]*(len(ni)+7)
 for i,c in enumerate(nd):
  for j,d in enumerate(ni):out[i+j]=A.add(out[i+j],A.mul(c,d))
 for k in range(len(out)-1,7,-1):
  for j,c in enumerate(quotient[:-1]):out[k-8+j]=A.sub(out[k-8+j],A.mul(out[k],c))
  out[k]=A.zero
 yes(out[0]==A.one and all(v==A.zero for v in out[1:]));marked=B.zero
 for c in reversed(qp):marked=B.add(B.mul(marked,sp),c)
 yes(marked==B.decode(row['marked_value_mod']));yes(B.mul(marked,B.decode(row['marked_inverse_mod']))==B.one)
 results.append({'name':name,'rank':row['rank'],'residual_squarefree':True,'marked_root_exactly_quadruple':True});emit('chart_pass',name=name,rank=row['rank'])
receipt={'result':'PASS','scope':'Full D!=0 coordinate-pole chart:174 distinct marked targets, and both ell-leading pole charts:21 each. All216 targets have one quadruple root and eight other simple roots. Complete source coverage by exact cofactor adjugate recovery and generic quadratic-linear equivalence; squarefree and marked-root units by literal two-engine replay in actual finite-free presentations.','checks':checks,'charts':results,'integer_resultant_points':119,'integer_cofactor_points':103,'integer_determinant_engines':['FLINT','plain Python Bareiss'],'finite_arithmetic_engines':['FLINT nmod_poly nested quotient','separate Python coefficient lists'],'checker_sha256':sha(Path(__file__)),'native_helper_sha256':sha(HERE/'build_remaining_charts.py'),'hashes':{f:sha(HERE/f) for f in ['coordinate-poles-candidate.json','coordinate-pole-resultant-candidate.json','remaining-charts-candidate.json','quadruple-input.json']},'seconds':time.monotonic()-START}
(HERE/'remaining-charts-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
