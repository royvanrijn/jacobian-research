"""Exact internal replay; universal claims rely on the retained written proof."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,prod,inf
import hashlib,json,platform
import sympy as s
HERE=Path(__file__).resolve().parent
REPO=next(x for x in HERE.parents if (x/"research/MATH_STATUS.json").is_file())
RESEARCH=REPO/"research/factorial/image-rigidity"
PINS={
 "LOWER_TWO_PRIME_GAP_SFC.tex":"8453b5b0ccd220573ef106fa26790e1f691885970b8c2116f036ed893ab2c200",
 "DIAGONAL_PRIME_GAP_SFC.tex":"db034b5fb29ece9c3a0cb8cbce4c82f692b933e3679b27ee823db5ca5e7f98eb",
 "LOWER_FOUR_QUADRATIC_GAP_SFC.tex":"ea2c5885f6f03a06ca2aa12840eddb191c128c20848ca8cf83d472268b89aca4",
 "lower_exponent_initial_and_boundary.tex":"e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45",
}
for name,digest in PINS.items():
 assert hashlib.sha256((RESEARCH/name).read_bytes()).hexdigest()==digest
counts={}
# The full symbolic d3,h4 upper identity modulo the lower monic ODE.
z,n,M=s.symbols("z n M")
def opadd(P,Q):
 out=dict(P)
 for j,a in Q.items():out[j]=s.factor(out.get(j,0)+a)
 return {j:a for j,a in out.items() if a!=0}
def compose(P,Q):
 out={}
 for i,a in P.items():
  for j,b in Q.items():
   for k in range(i+1):
    d=i-k+j;out[d]=out.get(d,0)+s.binomial(i,k)*a*s.diff(b,z,k)
 return {j:s.factor(a) for j,a in out.items() if a!=0}
def risingE(r):
 out={0:s.Integer(1)}
 for j in range(r):out=compose({0:M+j,1:-z},out)
 return out
C=M+3;Ah=s.rf(C,4);V={}
for j in range(4):V=opadd(V,{k:s.rf(C+j+1,3-j)*a for k,a in risingE(j).items()})
upper=opadd(risingE(7),{k:z**3*a for k,a in risingE(4).items()})
target=opadd({k:Ah*a for k,a in risingE(3).items()},{0:z**3*Ah})
target=opadd(target,{k:-3*(n+1)*z**3*a for k,a in V.items()})
diff=opadd(upper,{k:-a for k,a in target.items()})
L={0:3*n/z,1:-(M*(M-1)*(M+1)/z**3+1),2:3*M*(M-1)/z**2,3:-3*(M-1)/z,4:s.Integer(1)}
while diff and max(diff)>=4:
 j=max(diff);correction=compose({j-4:diff[j]},L)
 diff=opadd(diff,{k:-a for k,a in correction.items()})
assert not diff
counts["symbolic_d3_h4_Ore_remainder_coefficients"]=4
# Exact augmented3x3 determinant equals the explicit cubic2x2 minor.
eta,a,Delta=s.symbols("eta a Delta")
bi,bj,ci,cj,di,dj=s.symbols("bi bj ci cj di dj")
Z=s.Matrix([[di-a*ci,-ci,bi],[dj-a*cj,-cj,bj],
            [-eta*a*Delta,-eta*Delta,1-eta]])
pair=((1-eta)*ci-eta*Delta*bi)*dj-((1-eta)*cj-eta*Delta*bj)*di
assert s.expand(Z.det()-pair)==0
counts["symbolic_fulljet_cubic_elimination_identity"]=1
# Polynomial arithmetic modulo p, with gamma left formal.
def trim(A):
 A=list(A)
 while len(A)>1 and A[-1]==0:A.pop()
 return A
def padd(A,B,p,sgn=1):
 out=[0]*max(len(A),len(B))
 for i,a in enumerate(A):out[i]=(out[i]+a)%p
 for i,b in enumerate(B):out[i]=(out[i]+sgn*b)%p
 return trim(out)
def pmul(A,B,p):
 out=[0]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):out[i+j]=(out[i+j]+a*b)%p
 return trim(out)
def pscale(A,c,p):return trim([c*a%p for a in A])
def prf(offset,slope,m,p):
 out=[1]
 for j in range(m):out=pmul(out,[(offset+j)%p,slope%p],p)
 return out
def deriv(A,p):return trim([j*A[j]%p for j in range(1,len(A))] or [0])
def hasse(A,r,p):return trim([comb(j,r)*A[j]%p for j in range(r,len(A))] or [0])
def peval(A,x,p):
 out=0
 for a in reversed(A):out=(out*x+a)%p
 return out
kernel_polys=0;sigma_polys=0
for p0 in (2,3,5,7,11,13):
 k={m:deriv(prf(0,1,m,p0),p0) for m in range(1,p0+1)}
 assert k[p0]==[p0-1]
 for lam in range(1,p0):
  total=[0]
  for j in range(p0-1):
   term=pmul(prf(lam+1,-1,j,p0),k[p0-1-j],p0)
   total=padd(total,pscale(term,((-1)**j)*lam,p0),p0)
  assert padd(total,[1],p0)==[0]
  kernel_polys+=1
 if p0%2:
  assert hasse(prf(0,1,p0+1,p0),2,p0)==[p0-1]
  sigma_polys+=1
assert hasse(prf(0,1,2,2),2,2)==[1]
assert deriv(deriv(prf(0,1,2,2),2),2)==[0]
counts["formal_gamma_kernel_polynomial_identities"]=kernel_polys
counts["formal_gamma_odd_prime_plus_one_Hasse_identities"]=sigma_polys
counts["binary_Hasse_not_ordinary_second_derivative_control"]=1

def matvec(A,v,p):return [sum(a*b for a,b in zip(row,v))%p for row in A]
def matmul(A,B,p):
 return [[sum(A[i][k]*B[k][j] for k in range(len(B)))%p for j in range(len(B[0]))] for i in range(len(A))]
def determinant(A,p):
 A=[row[:] for row in A];out=1;N=len(A)
 for j in range(N):
  pivot=next((i for i in range(j,N) if A[i][j]%p),None)
  if pivot is None:return 0
  if pivot!=j:A[pivot],A[j]=A[j],A[pivot];out=-out
  a=A[j][j]%p;out=out*a%p;ai=pow(a,-1,p)
  for i in range(j+1,N):
   b=A[i][j]*ai%p
   A[i]=[(x-b*y)%p for x,y in zip(A[i],A[j])]
 return out%p
def Z_reduced(p,gamma,n):
 A=[[0]*p for _ in range(p)]
 for j in range(p):
  A[j][j]=(gamma-1-j)%p
  if j<p-1:A[j+1][j]=1
 # Here ubar0, exactly the positive-valuation case used in the determinant.
 powers=[[1]+[0]*(p-1)]
 for j in range(n+p):powers.append(matvec(A,powers[-1],p))
 B=powers[n];L0=powers[n+1];Ls=[L0]
 for j in range(p-2):
  AL=matvec(A,Ls[-1],p)
  Ls.append([((gamma-1-j)*v-a)%p for v,a in zip(Ls[-1],AL)])
 columns=[[( -x)%p for x in Ls[p-1-m]] for m in range(1,p)]+[B]
 return [[columns[j][i] for j in range(p)] for i in range(p)]
augmented=0;kernel_states=0
for p0 in (2,3,5,7):
 for gamma in range(p0):
  kval=[peval(deriv(prf(0,1,m,p0),p0),gamma,p0) for m in range(1,p0+1)]
  for n0 in range(1,p0):
   Zm=Z_reduced(p0,gamma,n0)
   assert not any(matvec(Zm,kval,p0))
   ep=[0]*(p0-1)+[1]
   chosen=None
   for omit in range(p0):
    I=[i for i in range(p0) if i!=omit]
    c0=determinant([Zm[i] for i in I]+[ep],p0)
    if c0:chosen=(I,c0);break
   assert chosen
   kernel_states+=1
   I,c0=chosen
   for h in range(1,p0+2):
    v=[]
    for m in range(1,p0+1):
     v.append(prod((gamma+m+1+j)%p0 for j in range(h-1-m))%p0 if m<h else 0)
    sigma=peval(hasse(prf(0,1,h,p0),2,p0),gamma,p0)
    assert sum(a*b for a,b in zip(v,kval))%p0==sigma
    linear=determinant([Zm[i] for i in I]+[[(-a)%p0 for a in v]],p0)
    assert linear==c0*sigma%p0
    if sigma:assert linear
    augmented+=1
counts["complete_sample_prime_residue_positive_power_kernel_states"]=kernel_states
counts["augmented_affine_minor_coefficient_identities"]=augmented

def rf_int(x,m):return prod(range(x,x+m))
def Bfamily(n,p,gamma):
 B=[s.Integer(1)];out={0:B}
 for r in range(n+p-1):
  D=[s.Integer(0)]*(len(B)+1)
  for j,a in enumerate(B):
   D[j]+=(gamma+p*r+p-1-j)*a;D[j+1]+=a
  B=D
  if r+1>=n:out[r+1]=B
 return out
def pcyclic(B,row,p,U):
 return sum(a*U**(j//p) for j,a in enumerate(B) if j%p==row)
U,X=s.symbols("U X")
def affine_minus_check(coeff,p,c):
 assert c and coeff[0]%p==c and coeff[1]%p==(-c)%p
 assert all(a%p==0 for a in coeff[2:])
literal=[]
for p0,n0 in ((3,2),(3,3),(3,7),(3,11),(5,4),(5,6)):
 gamma=(p0+1)*n0+1;M0=gamma+p0*n0;C0=M0+p0
 A0=rf_int(C0,p0+1);q=s.Rational(A0,p0*(n0+1))
 assert q.q==1 and int(q)%p0==0
 fam=Bfamily(n0,p0,gamma)
 B=fam[n0];Ls=[fam[n0+1]]
 for j in range(p0-2):
  P=Ls[-1];D=[s.Integer(0)]*(len(P)+1)
  for k,a in enumerate(P):D[k]+=(k-j)*a;D[k+1]-=a
  Ls.append(D)
 Zcols=[]
 for m in range(1,p0):Zcols.append([-pcyclic(Ls[p0-1-m],row,p0,U) for row in range(p0)])
 Zcols.append([pcyclic(B,row,p0,U) for row in range(p0)])
 Zm=s.Matrix(p0,p0,lambda i,j:Zcols[j][i])
 Z0=[[int(Zm[i,j].subs(U,0))%p0 for j in range(p0)] for i in range(p0)]
 ep=[0]*(p0-1)+[1]
 for omit in range(p0):
  I=[i for i in range(p0) if i!=omit]
  c0=determinant([Z0[i] for i in I]+[ep],p0)
  if c0:break
 assert c0
 v=[rf_int(C0+m+1,p0-m) for m in range(1,p0+1)]
 rows=[[Zm[i,j].subs(U,q*X) for j in range(p0)] for i in I]
 rows.append([ep[j]-X*v[j] for j in range(p0)])
 F=s.Poly(s.Matrix(rows).det(method="domain-ge"),X,domain=s.ZZ)
 coeff=[int(F.nth(j)) for j in range(F.degree()+1)]
 affine_minus_check(coeff,p0,c0)
 literal.append({"p":p0,"n":n0,"rows":I,"minor_degree":int(F.degree()),
                 "full_coefficient_reduction":"c(1-X)"})
counts["full_literal_augmented_minors"]=len(literal)
# Deliberately invalidate the critical unit linear coefficient.
bad=coeff.copy();bad[1]+=c0
try:affine_minus_check(bad,p0,c0)
except AssertionError:pass
else:raise AssertionError("Corrupted unit linear coefficient accepted")
counts["corrupted_linear_coefficient_rejected"]=1

# Actual ramified coefficients in Q(sqrt(3)): q may lie in m but outside3O.
def vp(x):
 x=Q(x)
 if not x:return inf
 a,b=abs(x.numerator),x.denominator;v=0
 while a%3==0:a//=3;v+=1
 while b%3==0:b//=3;v-=1
 return v
def qa(a=0,b=0):return Q(a),Q(b)
def qadd(x,y):return x[0]+y[0],x[1]+y[1]
def qneg(x):return -x[0],-x[1]
def qmul(x,y):return x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def qscale(x,a):return x[0]*a,x[1]*a
def qvalue(x):return min(vp(x[0]),vp(x[1])+Q(1,2))
def qreduce(x):
 v=qvalue(x);assert v>=0
 if v>0:return 0
 assert vp(x[0])==0 and vp(x[1])+Q(1,2)>0
 return x[0].numerator*pow(x[0].denominator,-1,3)%3
def fqadd(A,B,sgn=1):
 C=[qa()] * max(len(A),len(B))
 for j,a in enumerate(A):C[j]=qadd(C[j],a)
 for j,b in enumerate(B):C[j]=qadd(C[j],qscale(b,sgn))
 return C
def fqmul(A,B):
 C=[qa()]*(len(A)+len(B)-1)
 for j,a in enumerate(A):
  for k,b in enumerate(B):C[j+k]=qadd(C[j+k],qmul(a,b))
 return C
def fqscale(A,c):return [qmul(a,c) for a in A]
def fqvar(A,q):
 power=qa(1);out=[]
 for a in A:out.append(qmul(a,power));power=qmul(power,q)
 return out
ramified=[]
for n0,beta0,coeff in ((2,10,3),(3,1,3),(3,1,9)):
 gamma=qa(4*n0+beta0,coeff)
 Delta=qadd(gamma,qa(3*n0+6))
 B=[qa(1)];fam={}
 for r in range(n0+2):
  D=[qa()]*(len(B)+1)
  for j,a in enumerate(B):
   D[j]=qadd(D[j],qmul(qadd(gamma,qa(3*r+2-j)),a))
   D[j+1]=qadd(D[j+1],a)
  B=D
  if r+1>=n0:fam[r+1]=B
 A=qa(1)
 for j in range(3*n0+3,3*n0+7):A=qmul(A,qadd(gamma,qa(j)))
 q=qscale(A,Q(1,3*(n0+1)))
 assert qvalue(q)==(Q(1,2) if coeff==3 else Q(3,2))
 for i,j in ((0,1),(0,2),(1,2)):
  ci=fam[n0+1][i];cj=fam[n0+1][j]
  di=fam[n0+2][i];dj=fam[n0+2][j]
  if qreduce(qadd(qmul(ci,dj),qneg(qmul(cj,di)))):break
 else:raise AssertionError("No unit cubic row pair")
 b=[fam[n0][k::3] or [qa()] for k in (i,j)]
 c=[fam[n0+1][k::3] or [qa()] for k in (i,j)]
 d=[fam[n0+2][k::3] or [qa()] for k in (i,j)]
 b,c,d=[[fqvar(P,q) for P in V] for V in (b,c,d)]
 H=[fqadd(fqadd(cj,[qa()]+cj,-1),[qa()]+fqscale(bj,Delta),-1) for cj,bj in zip(c,b)]
 F=fqadd(fqmul(H[0],d[1]),fqmul(H[1],d[0]),-1)
 reductions=[qreduce(a) for a in F]
 assert reductions[0] and reductions[1]==(-reductions[0])%3 and not any(reductions[2:])
 ramified.append({"beta":str(beta0)+"+"+str(coeff)+"sqrt(3)","n":n0,"rows":[i,j],
                  "q_valuation":str(qvalue(q)),
                  "all_minor_coefficients_integral":True,"minor_reduction":reductions,
                  "higher_coefficient_valuations":[str(qvalue(a)) for a in F[2:]]})
counts["actual_ramified_minor_controls"]=len(ramified)
# The uniform lower-four corollary uses these exact characteristic-zero values.
g=s.symbols("g")
N4=s.prod(g+j for j in range(4))
sigma4=s.expand(N4).diff(g,2)/2
assert s.expand(sigma4-(6*g*g+18*g+11))==0
assert [sigma4.subs(g,a) for a in (0,-1,-2,-3)]==[11,-1,-1,11]
counts["exact_lower_four_Hasse_root_values"]=4
out={
 "passed":True,
 "scope":"Exact internal operator/kernel/Hasse/minor controls; universal claims rely on written proofs",
 "pointwise_domain":"primep,1<=h<=p+1,positive real algebraic gamma integral at one p-place; sigma_h nonzero or v(kappa_h)>=0",
 "prime_plus_one_family":"allprimep, alln, beta positive real algebraic integral at one p-place",
 "lower_four_prime_family":"oddprimep!=11 with the same algebraic integral-at-one-place Gamma domain; ordinaryp2 imported",
 "ordinary_support_four_seven":True,
 "all_positive_real_Gamma_shapes_claimed":False,
 "unrestricted_binomial_SFC_claimed":False,
 "controls":counts,
 "literal_augmented_minor_controls":literal,
 "ramified_controls":ramified,
 "source_sha256":{
  "PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex":hashlib.sha256((HERE/"PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex").read_bytes()).hexdigest(),
  Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  **PINS},
 "python":platform.python_version(),"sympy":s.__version__,
 "command":"python research/factorial/image-rigidity/prime_plus_one_exponent_gap_sfc.py",
 "independent_replay":False,"formal_verification":False,"external_review":False}
(HERE/"prime_plus_one_exponent_gap_sfc_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
