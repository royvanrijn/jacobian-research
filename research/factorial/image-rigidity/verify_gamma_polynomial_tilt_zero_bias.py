#!/usr/bin/env python3
"""Exact/source controls of the positive Gamma-tilt zero-bias theorem and obstruction."""
import argparse,ast,copy,hashlib,itertools,json,sys
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
if not __debug__:raise SystemExit("Optimized Python is unsupported.")
HERE=Path(__file__).resolve().parent
PN="GAMMA_POLYNOMIAL_TILT_ZERO_BIAS_CLOSURE_AND_OBSTRUCTION.md"
CN="verify_gamma_polynomial_tilt_zero_bias.py"
IN="gamma_polynomial_tilt_inputs.json"
LN="gamma_polynomial_tilt_source_ledger.json"
RN="gamma_polynomial_tilt_verification.json"
EXPECTED={'barriers': {'all_order_resultant_proved': False, 'factorial_SFC_counterexample': False, 'law_class_positivity_implies_ordered_separation': False, 'no_fixed_finite_band_recurrence_in_THIS_Appell_basis': True, 'other_basis_or_multi_index_excluded_claim': False, 'pure_product_or_Gamma_power_counterexample': False}, 'closure': {'coefficientwise_rescaling': '(gamma)_(K+d) A N', 'finite_product_Gamma_closed_claim': False, 'law_identity_uses_full_test_functions': True, 'moment_determinacy_assumed': False, 'moment_matching_alone_identifies_law': False, 'nonnegative_coefficients': True, 'translated_zero_bias_law': 'X^d under tilt X^d N(X)'}, 'probability_domain': {'all_moments_exist': True, 'base_rate': 1, 'd': 'positive integer', 'gamma': 'real gamma>0', 'tilt': 'nonzero polynomial with nonnegative real coefficients'}, 'schema': 1, 'source_boundary': {'GR_primary': 'classical zero-bias definition/existence/uniqueness only', 'old_convolution_note': 'authenticated retained-example credit only; new witness verified directly', 'unrelated_multiple_orthogonality_or_branched_fraction_PDFs_imported': False}, 'verification': {'bounded_exact_controls_only': True, 'external_review': False, 'finite_controls_prove_universal_closure': False, 'formal_verification': False, 'independent_replay': False, 'worldwide_priority': False}, 'witness': {'base_rate': 1, 'degrees': [4, 5], 'gamma': '5014/361', 'mean': '264/19', 'mixture_nodes': [0, 1, 24, 25, 88, 89], 'new_scope': 'exact finite positive polynomial tilt realization of the retained collision', 'old_Gamma_convolution_collision_rediscovered_claim': False, 'old_collision_scale': 2, 'shared_factor': 'u^2-2u+76', 'shared_factor_discriminant': -300, 'tilt_degree': 89, 'variance': '14'}}
HARD={"proof":"28b7302c8b63accdb41872bcb76655881bcb912742355a04acac3b0feb425113",
 "zero_bias_primary":"e435d58b7cba5567e53a11c269d7cf305c5af0870b51ae6cdb6e79ed25f9b983",
 "retained_collision":"672dd8c896c0cda82227209cb11a05b313cbd710c41d3c4be4c98119c5929790"}
SPECS=[
 ("proof","here",PN,"WRITTEN_CLOSURE_AND_EXACT_FINITE_REALIZATION"),
 ("checker","here",CN,"BOUNDED_EXACT_CONTROLS_NOT_UNIVERSAL_QUANTIFIER_PROOF"),
 ("inputs","here",IN,"EXACT_TYPED_SCOPE"),
 ("zero_bias_primary","repo","research/factorial/image-rigidity/sources/goldstein_reinert_zero_bias_author2001.pdf","CLASSICAL_ZERO_BIAS_IDENTIFICATION_AND_ATTRIBUTION_ONLY"),
 ("retained_collision","repo","research/factorial/image-rigidity/GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex","AUTHENTICATED_RETAINED_COLLISION_CREDIT_ONLY")]
PRIMARY={"title":"Stein's Method and the Zero Bias Transformation with Application to Simple Random Sampling","authors":"Larry Goldstein and Gesine Reinert","published":"Ann. Appl. Probab.7(1997),935-952","doi":"10.1214/aoap/1043862419","retained_version":"Author manuscript dated28Nov2001","url":"https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/zer.pdf","locators":["Definition1.1 printedp2","Existence/uniqueness printedpp2-3","Lemma2.1 printedp4"],"scope":"Centered finite nonzero variance zero bias. Translate back to Y; the Gamma-law identity here is proved with full compactly-supported C1 tests, not from moments alone."}
FLAGS={"external_review":False,"formal_verification":False,"independent_replay":False,"worldwide_priority":False}
def need(v,m):
 if not v:raise ValueError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def exact(a,b):return json.dumps(a,sort_keys=True,separators=(",",":"))==json.dumps(b,sort_keys=True,separators=(",",":"))
def input_check(d):need(exact(d,EXPECTED),"Exact typed scope")
def ledger_check(d):
 need(set(d)=={"schema","sources","primary","assurance","logical_boundary"},"Ledger keys")
 need(type(d["schema"]) is int and d["schema"]==1,"Ledger schema")
 need(exact(d["primary"],PRIMARY),"Primary scope/version")
 need(exact(d["assurance"],FLAGS),"Assurance flags")
 need(exact(d["logical_boundary"],EXPECTED["source_boundary"]),"Logical source boundary")
 rows=d["sources"]
 need(type(rows) is list and len(rows)==len(SPECS),"Mandatory role count")
 need(len({r.get("role") for r in rows if isinstance(r,dict)})==len(SPECS),"Duplicate/missing role")
 by={r["role"]:r for r in rows}
 for role,root,path,b in SPECS:
  r=by.get(role)
  need(type(r) is dict and set(r)=={"role","root","path","sha256","boundary"},"Source row schema")
  need((r["root"],r["path"],r["boundary"])==(root,path,b),"Source identity/scope")
  h=r["sha256"];need(type(h) is str and len(h)==64 and all(c in "0123456789abcdef" for c in h),"SHA format")
 return rows
def sources(repo):
 input_check(load(HERE/IN));rows=ledger_check(load(HERE/LN));pins={}
 for r in rows:
  base=HERE if r["root"]=="here" else repo;p=(base/r["path"]).resolve()
  try:p.relative_to(base)
  except ValueError:raise ValueError("Source path escaped root")
  need(p.is_file(),"Missing source "+r["role"]);h=sha(p)
  need(h==r["sha256"],"Current source hash "+r["role"])
  if r["role"] in HARD:need(h==HARD[r["role"]],"Reviewed source mutation "+r["role"])
  pins[r["role"]]=h
 return rows,pins
def rf(x,n):
 a=F(1)
 for j in range(n):a*=x+j
 return a
def vmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def appell(m,mom,shift=F(0)):
 out=[F(0)]*(m+1)
 for k in range(m+1):
  for j in range(m-k+1):out[j]+=comb(m,k)*mom[k]*comb(m-k,j)*shift**(m-k-j)
 return out
# Formal sparse polynomials in independent commuting parameters; no CAS.
def cl(p):return {m:F(v) for m,v in p.items() if v}
def co(c,n):return {} if not c else {(0,)*n:F(c)}
def va(i,n):
 m=[0]*n;m[i]=1;return {tuple(m):F(1)}
def ad(*ps):
 o={}
 for p in ps:
  for m,v in p.items():o[m]=o.get(m,F(0))+v
 return cl(o)
def sc(p,c):return cl({m:v*c for m,v in p.items()})
def mu(a,b):
 o={}
 for m,v in a.items():
  for n,w in b.items():
   k=tuple(i+j for i,j in zip(m,n));o[k]=o.get(k,F(0))+v*w
 return cl(o)
def pw(p,n,nv):
 o=co(1,nv)
 for j in range(n):o=mu(o,p)
 return o
def rising(p,n,nv):
 o=co(1,nv)
 for j in range(n):o=mu(o,ad(p,co(j,nv)))
 return o
def sum_poly(ps,nv):
 o={}
 for p in ps:o=ad(o,p)
 return o
def poly_record(p):return [[list(m),str(v)] for m,v in sorted(p.items())]
def formal_closure():
 rows=[];ncoeff=0
 for K in range(4):
  for d in range(1,5):
   nv=K+2;g=va(0,nv);ww=[va(r+1,nv) for r in range(K+1)]
   aa=[mu(ww[r],rising(g,r,nv)) for r in range(K+1)]
   rr=[rising(ad(g,co(r,nv)),d,nv) for r in range(K+1)]
   A=sum_poly(aa,nv);B=sum_poly([mu(aa[r],rr[r]) for r in range(K+1)],nv)
   def part(vv,u):return sum_poly(vv[:max(0,min(K,u)+1)],nv)
   hats=[]
   for j in range(K+d):
    Aj=part(aa,j);Au=part(aa,j-d);Bu=part([mu(aa[r],rr[r]) for r in range(K+1)],j-d)
    numerator=ad(mu(B,Aj),sc(mu(A,Bu),-1))
    positive=mu(B,ad(Aj,sc(Au,-1)))
    for r in range(K+1):
     for s in range(K+1):
      if s<=j-d<r:positive=ad(positive,mu(mu(aa[r],aa[s]),ad(rr[r],sc(rr[s],-1))))
    need(numerator==positive,"Formal positive cross-sum identity Z5")
    hat=mu(rising(ad(g,co(j+1,nv)),K+d-j-1,nv),numerator)
    need(all(v>=0 for v in hat.values()),"Formal coefficientwise positivity")
    hats.append(hat);ncoeff+=len(hat)
   norm=rising(g,K+d,nv)
   for j in range(K+d+1):
    left=mu(ad(g,co(j,nv)),hats[j]) if j<len(hats) else {}
    if j:left=ad(left,sc(hats[j-1],-1))
    rhs={}
    if j<=K:rhs=ad(rhs,mu(B,ww[j]))
    if 0<=j-d<=K:rhs=ad(rhs,sc(mu(A,ww[j-d]),-1))
    need(left==mu(norm,rhs),"Formal rescaled Stein coefficient identity")
   need(hats[-1]==mu(mu(norm,A),ww[-1]),"Formal exact top coefficient")
   rows.append([K,d,K+d-1,len(hats),sum(len(p) for p in hats)])
 # Genuine d2 denominator has nonzero numerator at its pole, for every gamma>0.
 g=va(0,1)
 atpole=mu(mu(ad(sc(g,-1),co(-1,1)),sc(g,-1)),co(2,1))
 need(atpole==sc(mu(g,ad(g,co(1,1))),2),"d2 denominator noncancellation polynomial")
 return {"formal_K_d_domains":rows,"nonnegative_formal_coefficients":ncoeff,
  "d2_denominator_pole_numerator":"2*gamma*(gamma+1)>0",
  "scope":"Exact formal identities for gamma and every original tilt coefficient in the listed K,d domains; written proof supplies all K,d."}
def numerical_closure():
 count=0;remainders=0
 for gamma in [F(1,2),F(3,2),F(5,3)]:
  for d in [1,2,3]:
   for w in [[F(1)],[F(0),F(1)],[F(2),F(0),F(1)]]:
    K=len(w)-1;A=sum(w[r]*rf(gamma,r) for r in range(K+1))
    moms=[sum(w[r]*rf(gamma,r+d*j) for r in range(K+1))/A for j in range(7)]
    mean=moms[1];V=moms[2]-mean**2;need(V>0,"Positive variance")
    N=[];prev=F(0)
    for j in range(K+d):
     val=(prev+mean*(w[j] if j<=K else 0)-(w[j-d] if 0<=j-d<=K else 0))/(gamma+j)
     N.append(val);prev=val
    need(all(v>=0 for v in N) and N[-1]==w[-1],"Numerical kernel positivity/top")
    norm=sum(v*rf(gamma,j+d) for j,v in enumerate(N));need(norm==A*V/d,"Law normalization")
    stars=[]
    for j in range(5):
     star=sum(v*rf(gamma,i+d*(j+1)) for i,v in enumerate(N))/norm
     need(star==(moms[j+2]-mean*moms[j+1])/((j+1)*V),"Actual zero-bias moment/chain")
     stars.append(star)
    for m in range(2,7):
     pm=appell(m,moms);other=vmul([mean,F(1)],appell(m-1,moms))
     left=[pm[j]-other[j] for j in range(m+1)]
     right=[(m-1)*V*v for v in appell(m-2,stars)]+[F(0)]*2
     need(left==right,"Literal first remainder uses new law")
     remainders+=1
    count+=1
 # Missing X^d changes the normalized law, even in the simplest d2 kernel.
 gamma=F(1);N=[F(2),F(1)]
 right=sum(v*rf(gamma,i+4) for i,v in enumerate(N))/sum(v*rf(gamma,i+2) for i,v in enumerate(N))
 wrong=sum(v*rf(gamma,i+2) for i,v in enumerate(N))/sum(v*rf(gamma,i) for i,v in enumerate(N))
 need(right!=wrong,"Missing X^d tilt countercontrol")
 return {"normalized_closure_cases":count,"literal_first_remainder_identities":remainders}
def cycles(p):
 seen=set();n=0
 for i in range(len(p)):
  if i in seen:continue
  n+=1;j=i
  while j not in seen:seen.add(j);j=p[j]
 return n
def transitive(ps):
 reached={0}
 while True:
  bigger=reached|{p[i] for p in ps for i in reached}
  if bigger==reached:return len(reached)==len(ps[0])
  reached=bigger
def cumulant_controls():
 domains=[];trans=[];recurrence_coeffs=0
 for r in [1,2,3]:
  nv=r+1;bs=[va(i,nv) for i in range(r)];x=va(r,nv)
  moments=[co(1,nv)]
  for n in range(1,7):
   a=co(1,nv)
   for b in bs:a=mu(a,rising(b,n,nv))
   moments.append(a)
  kap=[{}]
  for n in range(1,7):
   a=moments[n]
   for j in range(1,n):a=ad(a,sc(mu(kap[j],moments[n-j]),-comb(n-1,j-1)))
   need(a and all(v>0 for v in a.values()),"Positive product-Gamma cumulants")
   kap.append(a);domains.append([r,n,len(a)])
  def P(n):return sum_poly([sc(mu(moments[k],pw(x,n-k,nv)),comb(n,k)) for k in range(n+1)],nv)
  for n in range(6):
   identity=ad(mu(x,P(n)),sc(P(n+1),-1))
   for j in range(n+1):identity=ad(identity,sc(mu(kap[j+1],P(n-j)),comb(n,j)))
   need(not identity,"Literal formal fixed-law Appell recurrence Z11")
   recurrence_coeffs+=n+2
  for n in range(1,5 if r<=2 else 4):
   perms=list(itertools.permutations(range(n)));poly={}
   for ps in itertools.product(perms,repeat=r):
    if transitive(ps):
     m=tuple([cycles(p) for p in ps]+[0]);poly[m]=poly.get(m,F(0))+1
   need(poly==kap[n],"Transitive tuple/exponential-formula cumulant")
   trans.append([r,n,int(sum(poly.values())),len(poly)])
 return {"positive_cumulant_domains":domains,"transitive_permutation_tuple_domains":trans,
  "literal_Appell_recurrence_coefficients":recurrence_coeffs,
  "scope":"Finite exact combinatorial/formal controls; the written labeled-component proof supplies all n and all numbers of Gamma factors."}
def exact_witness():
 gamma=F(5014,361);alpha=F(125,9);beta=F(1,3249)
 weights=[(0,F(23959174167747319,24009255465548716)),(1,F(787871138391581,404852698201373433)),(24,F(337387137375,4512464765741536)),(25,F(578110771745,9270172181795112)),(88,F(14673828944125,11828739703970562912)),(89,F(70917583625,48018510931097432))]
 need(gamma==alpha+beta and all(p>0 for n,p in weights) and sum(p for n,p in weights)==1,"Literal positive mixture")
 moments=[sum(p*rf(gamma+n,k) for n,p in weights) for k in range(6)]
 old=[sum(F(comb(k,j))*rf(alpha,j)*rf(beta,k-j)*19**(k-j) for j in range(k+1)) for k in range(6)]
 need(moments==old,"Rate-one doubled old collision moments")
 mean=moments[1];V=moments[2]-mean**2;need((mean,V)==(F(264,19),F(14)),"Rate-one mean/variance")
 q=[F(76),F(-2),F(1)];p4=appell(4,moments,-mean);p5=appell(5,moments,-mean)
 need(p4==vmul(q,[F(12),F(2),F(1)]),"Fourth factorization")
 need(p5==vmul([F(4),F(1)],vmul(q,q)),"Fifth factorization")
 need(q[1]**2-4*q[0]*q[2]==-300,"Nonreal factor discriminant")
 w=[F(0)]*90
 for n,p in weights:w[n]=p/rf(gamma,n)
 need(sum(w[j]*rf(gamma,j) for j in range(90))==1,"Literal degree89 tilt")
 N=[];previous=F(0)
 for j in range(90):
  now=(previous+mean*w[j]-(w[j-1] if j else 0))/(gamma+j);N.append(now);previous=now
 need(all(v>=0 for v in N) and N[-1]==w[-1],"Degree89 zero-bias kernel")
 for j in range(91):
  left=(gamma+j)*N[j] if j<90 else F(0)
  if j:left-=N[j-1]
  right=mean*w[j] if j<90 else F(0)
  if j:right-=w[j-1]
  need(left==right,"Every degree89 Stein coefficient")
 norm=sum(v*rf(gamma,j+1) for j,v in enumerate(N));need(norm==V,"Exact law normalization")
 star=[sum(v*rf(gamma,j+k+1) for j,v in enumerate(N))/norm for k in range(4)]
 for k in range(4):need(star[k]==(moments[k+2]-mean*moments[k+1])/((k+1)*V),"Literal zero-bias moment")
 b3=appell(3,star,-mean);need(b3==vmul(q,[F(38,7),F(1)]),"Cross-law remainder collision")
 need([p5[j]-(p4[j-1] if j else 0) for j in range(6)]==[4*V*v for v in b3]+[F(0)]*2,"Exact Euclidean first remainder")
 # Class distinctions and honest corrupt controls.
 need(3*F(32)**2!=2*F(14)*F(324),"Not a single scaled Gamma")
 wrong_rate_moments=[sum(p*rf(gamma+n,k)/2**k for n,p in weights) for k in range(6)]
 need(wrong_rate_moments!=moments and wrong_rate_moments[1]==mean/2,"Wrong rate2/doubling countercontrol")
 need(sum(p for n,p in weights[:-1])!=1,"Omitted positive mixture atom")
 need(sum(p+(F(1,1000) if n==0 else 0) for n,p in weights)!=1,"Corrupted atom")
 return {"base_shape":str(gamma),"base_rate":1,"tilt_degree":89,"nodes_and_weights":[[n,str(p)] for n,p in weights],
 "raw_moments_0_through_5":[str(x) for x in moments],"mean":str(mean),"variance":str(V),
 "shared_centered_factor":"u^2-2u+76","discriminant":-300,
 "centered_P4":[str(v) for v in p4],"centered_P5":[str(v) for v in p5],
 "Stein_coefficients":91,"zero_bias_moments":4,
 "first_remainder":"56*(u^2-2u+76)*(u+38/7)",
 "mathematical_negative_controls":["omitted positive atom breaks normalization","corrupted atom breaks normalization","rate2 law is not the displayed rate1 law","single scaled Gamma cumulant identity fails"],
 "retained_collision_credit":"Exactly twice the retained Gamma convolution, matched only through moment5; finite positive polynomial realization is the new scope."}
def scope_controls(rows):
 n=0
 def reject(fn):
  nonlocal n
  try:fn()
  except (ValueError,TypeError,KeyError,IndexError):n+=1;return
  raise AssertionError("Bad scope/schema accepted")
 for path,value in [
 (("probability_domain","base_rate"),2),
 (("closure","moment_matching_alone_identifies_law"),True),
 (("closure","moment_determinacy_assumed"),True),
 (("closure","law_identity_uses_full_test_functions"),False),
 (("closure","finite_product_Gamma_closed_claim"),True),
 (("witness","old_collision_scale"),1),
 (("witness","base_rate"),2),
 (("witness","variance"),"7/2"),
 (("witness","old_Gamma_convolution_collision_rediscovered_claim"),True),
 (("barriers","other_basis_or_multi_index_excluded_claim"),True),
 (("barriers","law_class_positivity_implies_ordered_separation"),True),
 (("barriers","pure_product_or_Gamma_power_counterexample"),True),
 (("barriers","all_order_resultant_proved"),True),
 (("barriers","factorial_SFC_counterexample"),True),
 (("verification","finite_controls_prove_universal_closure"),True),
 (("verification","formal_verification"),True),
 (("schema",),True)]:
  a=copy.deepcopy(EXPECTED);q=a
  for k in path[:-1]:q=q[k]
  q[path[-1]]=value;reject(lambda a=a:input_check(a))
 template={"schema":1,"sources":rows,"primary":PRIMARY,"assurance":FLAGS,"logical_boundary":EXPECTED["source_boundary"]}
 for role,root,path,b in SPECS:
  a=copy.deepcopy(template);a["sources"]=[r for r in a["sources"] if r["role"]!=role];reject(lambda a=a:ledger_check(a))
 for mutate in [lambda a:a.update(schema=True),lambda a:a.update(schema=2),lambda a:a.update(sources=[]),lambda a:a["sources"].append(a["sources"][0]),lambda a:a["assurance"].update(formal_verification=1),lambda a:a["primary"].update(scope="Positive zero bias proves Appell separation"),lambda a:a["logical_boundary"].update(unrelated_multiple_orthogonality_or_branched_fraction_PDFs_imported=0)]:
  a=copy.deepcopy(template);mutate(a);reject(lambda a=a:ledger_check(a))
 return n
def run(repo):
 rows,pins=sources(repo)
 out={"schema":1,"status":"PASS_SOURCE_GUARDED_GAMMA_POLYNOMIAL_TILT_CLOSURE_AND_OBSTRUCTION",
 "scope":EXPECTED,"proof_sha256":pins["proof"],"checker_sha256":pins["checker"],"inputs_sha256":pins["inputs"],
 "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),"pins":rows,"primary":PRIMARY,
 "formal_closure_controls":formal_closure(),"numerical_normalized_closure_cases":numerical_closure(),
 "cumulant_and_no_fixed_band_controls":cumulant_controls(),"exact_witness":exact_witness(),
 "source_scope_negative_controls":scope_controls(rows),"assurance":FLAGS,
 "verification_boundary":"Universal closure/positivity and all-order no-fixed-band theorem are WRITTEN proofs. Finite formal controls and the complete rational degree89 witness have distinct scopes. No product-Gamma/factorial SFC counterexample or separation theorem is inferred."}
 return out
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--repo-root",type=Path,required=True);p.add_argument("--output",type=Path,default=HERE/RN);a=p.parse_args()
 repo=a.repo_root.resolve();need(repo.is_dir(),"Repository root")
 out=run(repo);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":out["status"],"roles":len(out["pins"]),"source_controls":out["source_scope_negative_controls"],"receipt_sha256":sha(a.output)},sort_keys=True))
if __name__=="__main__":
 try:main()
 except (ValueError,TypeError,KeyError,AssertionError) as exc:raise SystemExit("FAILED: "+str(exc))

