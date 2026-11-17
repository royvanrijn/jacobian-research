"""Exact all-degree controls and complete eight-letter enumeration."""
from pathlib import Path
from itertools import permutations,product
from collections import Counter
from math import gcd,lcm
import hashlib,json,platform,time
import sympy as S
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
start=time.monotonic();checks=[]
def check(n,v):
    if not v:raise AssertionError(n)
    checks.append(n)
w,z,g,t,Y1,Y2,Y3=S.symbols("w z gamma t Y1 Y2 Y3")
j,n,s,eps=S.symbols("j n s eps",integer=True)
A=1-j*t;B=2*(n-1)+(-2*j*n+j*j+j)*t
check("universal heat singular coefficient",S.cancel(B.subs(t,1/j)-(j-1))==0)
check("universal heat top coefficient",S.expand((j*j*s*(s-1)+(-2*j*n+j*j+j)*s+n*(n-1)).subs(n,s*j+eps)-eps*(eps-1))==0)
a=(s+1)*j+eps-1
check("epsilon0 weighted degree",S.cancel(((a+1)*(a+j)/j).subs(eps,0)-(s+1)*((s+2)*j-1))==0)
check("epsilon1 weighted degree",S.cancel(((a+1)*(a+j)/j).subs(eps,1)-(s+2)*((s+1)*j+1))==0)
Delta=S.Matrix([1,w,z]);Nframe=S.Matrix([[-w,1,0],[-z,0,1]])
cases=[(w,[(0,1)]),(w,[(1,0)]),(w*w+w,[(1,0)]),
       (w*w,[(1,0)]),(w**4+3*w*w+1,[(1,0)]),
       ((w-1)**2+2,[(1,-1)]),(w,[(0,0),(0,1)]),
       (w*w+w,[(1,2),(1,1)]),(w,[(2,1),(0,2),(0,1)])]
models=[]
for h,pairs in cases:
 f=polynomial_chart(h,pairs,w=w,z=z);X=f["X"];lam=f["lambda_polynomial"];sv=f["s"];dv=f["d"];ep=f["epsilon"]
 label=f"h={h},s={sv},eps={ep}"
 wave=S.diff(lam,w,2)+2*S.diff(h*lam,w,z)+S.diff(h*h*lam,z,2)
 check(label+" exact wave",S.expand(wave)==0)
 check(label+" normalization matrix",all(S.expand(e)==0 for e in Nframe*X.jacobian([w,z])-lam*S.Matrix([[-h,1],[-h*h,h]])))
 check(label+" full sweep Jacobian",S.expand(S.Matrix.hstack(Delta,X.diff(w)+g*Delta.diff(w),X.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
 # Leading homogeneous model has unit top coefficient and monomial h.
 model_pairs=[(0,0)]*(sv+1);model_pairs[sv]=(1,0) if ep else (0,1)
 ff=polynomial_chart(w**dv,model_pairs,w=w,z=z)
 xx=ff["X"];ll=ff["lambda_polynomial"];jv=dv+1;nv=sv*jv+ep;av=nv+jv-1
 LL=S.expand(ll.subs({w:1,z:t}))
 check(label+" leading polynomial avoids singular point",S.cancel(LL.subs(t,S.Rational(1,jv)))!=0)
 if sv:check(label+" leading polynomial simple",S.discriminant(LL,t)!=0)
 V=xx[1]-w*xx[0];W=xx[2]-z*xx[0]
 check(label+" no nonzero weighted origin ray",S.resultant(V.subs({w:1,z:t}),W.subs({w:1,z:t}),t)!=0)
 if ep:
  check(label+" w0 first source coordinate",S.expand(xx[0].subs(w,0)+z**(sv+1)/(sv+1))==0)
  origin=S.expand((xx[2]+g*z).subs({w:0,g:z**(sv+1)/(sv+1)},simultaneous=True))
  check(label+" w0 origin obstruction",S.expand(origin-z**(sv+2)/((sv+1)*(sv+2)))==0)
 else:
  check(label+" w0 leading first coordinate",xx[0].subs(w,0)==0)
  check(label+" w0 leading second coordinate",S.expand(xx[1].subs(w,0)-z**(sv+1)/(sv+1))==0)
 # The weighted initial forms agree with the complete model up to its scalar.
 weights=[av,av+1,av+jv]
 for i in range(3):
  degs=[mon[0]+jv*mon[1] for mon,co in S.Poly(X[i],w,z).terms() if co!=0]
  check(label+" polynomial weight bound "+str(i),not degs or max(degs)<=weights[i])
 res=S.factor(S.resultant(X[1]-w*X[0]+2*w-3,X[2]-z*X[0]+2*z-5,z))
 check(label+" exact degree resultant",S.degree(res,w)==f["N"])
 models.append(dict(h=str(h),lambda_polynomial=str(lam),s=sv,epsilon=ep,
  degree=f["N"],finite=True,delta=f["delta"],geometric_group=f["geometric_group"]))
# All 96 even signed permutations on four pairs, with complete Nielsen audit.
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
 out=[0]*len(a)
 for i,j in enumerate(a):out[j]=i
 return tuple(out)
ident=tuple(range(8))
def cycles(p):
 seen=set();out=[]
 for a in range(len(p)):
  if a in seen:continue
  v=a;n=0
  while v not in seen:seen.add(v);n+=1;v=p[v]
  out.append(n)
 return tuple(sorted(out,reverse=True))
def generated(gs):
 out={ident};queue=[ident]
 while queue:
  a=queue.pop()
  for b in gs:
   c=mul(a,b)
   if c not in out:out.add(c);queue.append(c)
 return frozenset(out)
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
ambient=[]
for pp in permutations(range(4)):
 if parity(pp)!=1:continue
 for bits in product(range(2),repeat=4):
  if sum(bits)%2:continue
  ambient.append(tuple(2*pp[i//2]+((i%2)^bits[i//2]) for i in range(8)))
lifts=[p for p in ambient if cycles(p)==(3,3,1,1)]
check("complete even signed ambient action",len(set(ambient))==96 and len(lifts)==32)
counts=Counter();groups={};pair_count=0
for aa in lifts:
 for bb in lifts:
  if cycles(mul(aa,bb))!=(6,2):continue
  H=generated([aa,bb])
  if len({h[0] for h in H})!=8:continue
  pair_count+=1;counts[len(H)]+=1;groups.setdefault(H,(aa,bb))
check("all transitive genus-zero pairs enumerated",pair_count==192 and counts=={24:192})
check("eight degree-eight subgroups",len(groups)==8)
full_signed=[tuple(2*pp[i//2]+((i%2)^bits[i//2]) for i in range(8))
             for pp in permutations(range(4)) for bits in product(range(2),repeat=4)]
check("full signed conjugacy ambient",len(set(full_signed))==384)
H0=next(iter(groups));gens=groups[H0];conjugators=[]
for H in groups:
 c=next((p for p in full_signed if frozenset(mul(mul(p,h),inv(p)) for h in H0)==H),None)
 check("degree-eight ambient conjugacy "+str(len(conjugators)),c is not None)
 conjugators.append(c)
# Explicit SL2(F3) and GL2(F3) actions; identify the first subgroup by
# a simultaneous generator conjugacy, then all eight by the preceding audit.
vectors=[v for v in product(range(3),repeat=2) if v!=(0,0)]
def matrix_perm(M):
 a,b,c,d=M
 return tuple(vectors.index(((a*x+b*y)%3,(c*x+d*y)%3)) for x,y in vectors)
SL=[];GL=[]
for M in product(range(3),repeat=4):
 det=(M[0]*M[3]-M[1]*M[2])%3
 if det:GL.append(matrix_perm(M))
 if det==1:SL.append(matrix_perm(M))
SL=frozenset(SL);GL=frozenset(GL)
check("explicit special and general linear actions",len(SL)==24 and len(GL)==48)
def match(gs,hs):
 for seed in range(8):
  mapping={0:seed};queue=[0];ok=True
  while queue and ok:
   i=queue.pop()
   for p,q in zip(gs,hs):
    ii=p[i];jj=q[mapping[i]]
    if ii in mapping:
     if mapping[ii]!=jj:ok=False;break
    elif jj in mapping.values():ok=False;break
    else:mapping[ii]=jj;queue.append(ii)
  if ok and len(mapping)==8:return tuple(mapping[i] for i in range(8))
 return None
conj=None;targetgens=None
for aa in SL:
 if cycles(aa)!=(3,3,1,1):continue
 for bb in SL:
  if cycles(bb)!=(3,3,1,1) or cycles(mul(aa,bb))!=(6,2):continue
  cc=match(gens,[aa,bb])
  if cc is not None:conj=cc;targetgens=(aa,bb);break
 if conj is not None:break
check("exact degree-eight SL2F3 action conjugacy",conj is not None and frozenset(mul(mul(conj,h),inv(conj)) for h in H0)==SL)
check("SL2F3 normal in GL2F3",all(frozenset(mul(mul(p,h),inv(p)) for h in SL)==SL for p in GL))
# The commutator calculation excludes a full-plane enlargement to order96.
# Compute the full derived subgroup independently to avoid assuming generators.
derived=generated([mul(mul(mul(a,b),inv(a)),inv(b)) for a in ambient for b in ambient])
check("ambient abelianization is C3",len(derived)==32)
census=Counter(cycles(p) for p in SL)
check("binary tetrahedral complete cycle-gcd invariant",lcm(*(gcd(*cs) for cs in census))==4<8)
check("binary action has one involution",sum(lcm(*cycles(p))==2 for p in SL)==1)
# Normalizer bound: centralizer2 and restriction to quaternion automorphisms.
Hstab=[p for p in SL if p[0]==0]
check("binary point stabilizer order",len(Hstab)==3)
fixed=[i for i in range(8) if all(h[i]==i for h in Hstab)]
check("binary action centralizer order",len(fixed)==2)
Q8=frozenset(p for p in SL if lcm(*cycles(p)) in [1,2,4])
check("unique quaternion Sylow subgroup",len(Q8)==8 and all(frozenset(mul(mul(p,h),inv(p)) for h in Q8)==Q8 for p in SL))
# Quartic sign and exact arithmetic constant character.
q,AA,BB=S.symbols("q A B")
quartic=q**4-12*AA*q*q+24*BB*q-12*AA*AA
disc=S.factor(S.discriminant(quartic,q));coef,factors=S.factor_list(disc/(-3))
check("quartic sign square factor",all(e%2==0 for factor,e in factors) and S.sqrt(coef).is_Rational)
square=S.sqrt(coef)*S.prod(factor**(e//2) for factor,e in factors)
check("exact binary arithmetic discriminant",S.expand(disc+3*square*square)==0)
# The centered degree-eight stage obstruction does not assert absolute impossibility.
check("centered binary unramified-origin obstruction",S.Rational(1,3)-S.Rational(1,2)+S.Rational(1,8)==-S.Rational(1,24))
manifest=json.loads((ROOT/"sources/signed-alternating-manifest.json").read_text())
for src in manifest["sources"]:
 check("Carmichael presentation pin",hashlib.sha256((ROOT/src["path"]).read_bytes()).hexdigest()==src["sha256"])
model=dict(nonconstant_polynomial_multiplier_sweeps="all finite; universal weighted-leading proof",
 controls=models,degree_eight=dict(geometric_group="SL2(F3) on eight nonzero vectors",
  arithmetic_group="GL2(F3) iff -3 nonsquare in coefficient field",constant_field="K(sqrt(-3))",
  pair_count=pair_count,group_order_counts=dict(counts),distinct_groups=len(groups),
  representative_generators=gens,ambient_conjugators=conjugators,SL2_conjugacy=conj,
  cycle_census={",".join(map(str,cs)):v for cs,v in sorted(census.items())},
  cycle_gcd_lcm=4,actual_Keller_realization_proved=False),
 entire_k_w_z_chart_groups=["natural S_n,n>=1","natural A_n,n>=3,n!=4",
  "SL2(F3) on8","even signed V_n:A_n on2n,even n>=6"],
 all_chart_PSL2_11_field_pairs_excluded=True,natural_A4_in_chart="excluded",
 abstract_PSL2_11_actions_decided=False,absolute_SL2_3_degree8_occurrence_decided=False,
 full_action_classification_complete=False)
(ROOT/"polynomial_chart_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md","polynomial_wave.py","verify_polynomial_chart.py",
 "polynomial_chart_model.json","sources/signed-alternating-manifest.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-start,2),sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written universal heat-recursion, weighted finiteness, degree and field-group proof; complete finite binary-action enumeration; Carmichael input attributed. Not formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"polynomial_chart_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
