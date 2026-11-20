#!/usr/bin/env python3
"""Literal cubic index17 collision and source reconstruction; no CAS required."""
if not __debug__:raise SystemExit("Refusing optimized Python: exact assertions required.")
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,copy,hashlib,json
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(p):return {m:c for m,c in p.items() if c}
def add(*ps):
 r={}
 for p in ps:
  for m,c in p.items():r[m]=r.get(m,F(0))+c
 return clean(r)
def mul(p,q):
 r={}
 for m,c in p.items():
  for n,d in q.items():
   k=tuple(sorted(m+n));r[k]=r.get(k,F(0))+c*d
 return clean(r)
def jacobian(terms):
 r={}
 for t in terms:
  for i,k in Counter(t["inputs"]).items():
   mon=list(t["inputs"]);mon.remove(i);mon=tuple(mon)
   p=r.setdefault((t["output"],i),{})
   p[mon]=p.get(mon,F(0))+F(t["coefficient"])*k
 return {ij:clean(p) for ij,p in r.items() if clean(p)}
def mmul(a,b):
 rows={}
 for (i,j),p in b.items():rows.setdefault(i,[]).append((j,p))
 r={}
 for (i,k),p in a.items():
  for j,q in rows.get(k,[]):r[(i,j)]=add(r.get((i,j),{}),mul(p,q))
 return {ij:p for ij,p in r.items() if p}
def ev(p,x):
 r=F(0)
 for m,c in p.items():
  for i in m:c*=x[i]
  r+=c
 return r
def validate(d):
 assert set(d)=={"schema_version","dimension","H_terms","collision_points","common_image",
                 "nilpotency_index","sixteenth_power_witness"}
 assert d["schema_version"]==1 and d["dimension"]==18 and d["nilpotency_index"]==17
 terms=d["H_terms"];assert len(terms)==45
 seen=set()
 for t in terms:
  assert set(t)=={"output","inputs","coefficient"}
  i=t["output"];m=t["inputs"]
  assert type(i) is int and 0<=i<18
  assert len(m)==3 and m==sorted(m) and all(type(k) is int and 0<=k<18 for k in m)
  assert (i,tuple(m)) not in seen;seen.add((i,tuple(m)))
  assert isinstance(t["coefficient"],str) and F(t["coefficient"])
 assert len(d["collision_points"])==2 and all(len(x)==18 for x in d["collision_points"])
 assert len(d["common_image"])==18 and d["collision_points"][0]!=d["collision_points"][1]
 w=d["sixteenth_power_witness"]
 assert w=={"row":3,"column":17,"inputs":[0]*10+[1]*4+[2]*3+[17]*15,"coefficient":"288"}
 return terms
def reconstruct(base):
 # The pinned12-input has scale11=1. Diagonal conjugation x11 -> x11/2
 # gives the scale11=2 map used in POWER89 and the source companion.
 assert len(base["coefficients"])==12
 diagonal=[F(1)]*12;diagonal[10]=F(1,2)
 K=[]
 for i,row in enumerate(base["coefficients"]):
  p={}
  for t in row:
   e=t["exponent"];c=F(t["coefficient"])/diagonal[i]
   assert len(e)==12 and all(type(n) is int and n>=0 for n in e)
   for a,n in zip(diagonal,e):c*=a**n
   mon=tuple(j for j,n in enumerate(e) for k in range(n))
   assert mon not in p and c;p[mon]=c
  K.append(p)
 Q=[{m:c for m,c in p.items() if len(m)==2} for p in K]
 C=[{m:c for m,c in p.items() if len(m)==3} for p in K]
 assert all(not p for p in C[5:])
 H=[]
 for i in range(12):
  p={tuple(sorted(m+(17,))):c for m,c in Q[i].items()}
  if i<5:p=add(p,{(12+i,17,17):F(-1)})
  H.append(p)
 H+=C[:5]+[{}]
 points=[[F(x)/a for x,a in zip(p,diagonal)] for p in base["points"]]
 images=[[x+ev(add(Q[i],C[i]),p) for i,x in enumerate(p)] for p in points]
 assert points[0]!=points[1] and images[0]==images[1]
 lifted=[p+[-ev(q,p) for q in C[:5]]+[F(1)] for p in points]
 return H,lifted
def numeric_power(j,n):
 out=[[F(i==k) for k in range(18)] for i in range(18)]
 for _ in range(n):
  out=[[sum((out[i][k]*j[k][l] for k in range(18)),F(0))
        for l in range(18)] for i in range(18)]
 return out
def jet_terms(seed,n,M):
 out={}
 for t in seed:
  a,b,c=t["inputs"];co=F(t["coefficient"])
  for grade in range(1,M+1):
   for i in range(grade):
    for j in range(grade-i):
     k=grade-1-i-j
     mon=tuple(sorted((a+n*i,b+n*j,c+n*k)))
     key=(t["output"]+n*grade,mon)
     out[key]=out.get(key,F(0))+co
 return [{"output":i,"inputs":list(mon),"coefficient":str(c)}
         for (i,mon),c in sorted(out.items()) if c]
def jet_operator(seedj,n,M):
 out={}
 for (row,col),p in seedj.items():
  for r in range(M+1):
   for k in range(r):
    q={}
    degree=r-k-1
    for mon,co in p.items():
     assert len(mon)==2
     a,b=mon
     for i in range(degree+1):
      m=tuple(sorted((a+n*i,b+n*(degree-i))))
      q[m]=q.get(m,F(0))+co
    if clean(q):out[(row+n*r,col+n*k)]=clean(q)
 return out
def series_mul(a,b,M):
 out=[{} for _ in range(M+1)]
 for i,p in enumerate(a):
  for k,q in enumerate(b[:M+1-i]):out[i+k]=add(out[i+k],mul(p,q))
 return out
def series_pow(a,k,M):
 out=[{():F(1)}]+[{} for _ in range(M)]
 for _ in range(k):out=series_mul(out,a,M)
 return out
def jet_controls(terms,j):
 from math import comb
 tests=[]
 for M in range(5):
  t=jet_terms(terms,18,M);jj=jacobian(t)
  assert jj==jet_operator(j,18,M)
  assert all(r//18>c//18 for r,c in jj)
  for a in t:assert sum(2*(i//18)+1 for i in a["inputs"])==2*(a["output"]//18)+1
  power=jj
  for _ in range(M):power=mmul(power,jj)
  assert not power
  tests.append({"M":M,"dimension":18*(M+1),"jet_cubic_terms":len(t),"full_jacobian_operator_entries":len(jj)})
 # Scalar x+x^3 has nonterminating inverse and exact jet inverse degrees2M+1.
 scalar=[]
 for M in range(1,5):
  V=[{(i,):F(1)} for i in range(M+1)]
  inverse=[{} for _ in range(M+1)]
  for k in range(M+1):
   coeff=F((-1)**k*comb(3*k,k),2*k+1)
   P=series_pow(V,2*k+1,M)
   for r in range(k,M+1):inverse[r]=add(inverse[r],{m:coeff*c for m,c in P[r-k].items()})
  cube=series_pow(inverse,3,M)
  image=[inverse[0]]+[add(inverse[r],cube[r-1]) for r in range(1,M+1)]
  assert image==V
  degree=max(len(mon) for p in inverse for mon in p)
  assert degree==2*M+1
  scalar.append({"M":M,"inverse_degree":degree})
 # The nonzero invertible seed(x,y)->(x,y+x^3) has jet inverse degree3
 # at every positive M, despite unbounded available grading weights.
 safe=[]
 for M in range(1,5):
  seed=[{"output":1,"inputs":[0,0,0],"coefficient":"1"}]
  jt=jet_terms(seed,2,M);n=2*(M+1)
  H=[{} for _ in range(n)]
  for t in jt:H[t["output"]][tuple(t["inputs"])]=F(t["coefficient"])
  V=[{(i,):F(1)} for i in range(n)]
  inverse=[add(V[i],{m:-c for m,c in H[i].items()}) for i in range(n)]
  assert all(inverse[i]==V[i] for i in range(0,n,2))
  assert max(len(mon) for p in inverse for mon in p)==3
  for i,p in enumerate(H):
   # Every nonlinear input has even index and is fixed by inverse.
   assert all(all(k%2==0 for k in mon) for mon in p)
   assert add(inverse[i],p)==V[i]
  safe.append({"M":M,"inverse_degree":3})
 # Gap-window endpoints used for the quantitative bound.
 for M in range(16,1017):
  K=(M+2)//3
  assert K>=1 and K<=3*K-2<=M
  assert 2*K+1<=2*M+1
 return {"literal_jet_operator_checks":tests,"scalar_noninvertible_jet_controls":scalar,
         "invertible_seed_bounded_jet_controls":safe,"gap_endpoint_checks":1001}

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument("--repo-root",type=Path,required=True)
 ap.add_argument("--input",type=Path,default=HERE/"cubic_index_seventeen_shuffle_failure_inputs.json")
 ap.add_argument("--ledger",type=Path,default=HERE/"cubic_index_seventeen_shuffle_failure_source_ledger.json")
 ap.add_argument("--output",type=Path,default=HERE/"cubic_index_seventeen_shuffle_failure_verification.json")
 a=ap.parse_args();repo=a.repo_root.resolve()
 led=json.loads(a.ledger.read_text());assert led["schema_version"]==1
 required={"proof","checker","input","retained_base_input","retained_companion_proof",
           "retained_shuffle_proof","primary_shuffle"}
 assert len(led["sources"])==len(required) and {x["role"] for x in led["sources"]}==required
 paths={}
 for x in led["sources"]:
  root=a.ledger.parent if x["scope"]=="packet" else repo
  p=(root/x["path"]).resolve()
  assert p.is_file() and sha(p)==x["sha256"],"Source pin mismatch: "+str(p)
  paths[x["role"]]=p
 assert paths["checker"]==Path(__file__).resolve() and paths["input"]==a.input.resolve()
 d=json.loads(a.input.read_text());terms=validate(d)
 H,points=reconstruct(json.loads(paths["retained_base_input"].read_text()))
 actual=[{} for _ in range(18)]
 for t in terms:actual[t["output"]][tuple(t["inputs"])]=F(t["coefficient"])
 assert H==actual and [[str(x) for x in p] for p in points]==d["collision_points"]
 target=[F(x) for x in d["common_image"]]
 assert [[x+ev(p,X) for x,p in zip(X,H)] for X in points]==[target,target]
 j=jacobian(terms);power=j;profiles=[]
 for n in range(1,18):
  if n>1:power=mmul(power,j)
  profiles.append({"power":n,"nonzero_entries":len(power),"nonzero_terms":sum(map(len,power.values()))})
  if n<17:assert power
  if n==16:
   w=d["sixteenth_power_witness"]
   assert power[(w["row"],w["column"])]=={tuple(w["inputs"]):F(w["coefficient"])}
 assert not power
 bad=copy.deepcopy(terms);bad[0]["coefficient"]=str(F(bad[0]["coefficient"])+1)
 jb=jacobian(bad);point=[F(1)]*18
 nj=[[ev(jb.get((i,k),{}),point) for k in range(18)] for i in range(18)]
 badpower=numeric_power(nj,17)
 assert any(c for row in badpower for c in row)
 negative=next((i,k,str(c)) for i,row in enumerate(badpower) for k,c in enumerate(row) if c)
 rejected=0
 for mode in ["duplicate","zero","degree","index"]:
  b=copy.deepcopy(d)
  if mode=="duplicate":b["H_terms"][-1]=copy.deepcopy(b["H_terms"][-2])
  elif mode=="zero":b["H_terms"][0]["coefficient"]="0"
  elif mode=="degree":b["H_terms"][0]["inputs"].append(0)
  else:b["H_terms"][0]["output"]=18
  try:validate(b)
  except (AssertionError,ValueError):rejected+=1
 assert rejected==4
 result={"schema_version":1,"status":"PASS","source_ledger_sha256":sha(a.ledger),
  "dimension":18,"literal_cubic_terms":45,"exact_nilpotency_index":17,
  "power_profiles":profiles,"sixteenth_power_witness":d["sixteenth_power_witness"],
  "all_source_coefficients_reconstructed":True,"distinct_rational_collision":True,
  "common_image":d["common_image"],"coefficient_perturbation_nonzero_J17_at_ones":negative,
  "malformed_inputs_rejected":rejected,
  "jet_controls":jet_controls(terms,j),
  "scope":"Direct QQ cubic identity/collision; infinite jet-growth and shuffle consequences use the written formal/graded arguments, not finite inverse-prefix sampling."}
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n");print(json.dumps(result))
if __name__=="__main__":main()
