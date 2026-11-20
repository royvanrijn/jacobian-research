#!/usr/bin/env python3
"""Guard the completed actual finite box and its all-terminal proof."""
if not __debug__:raise SystemExit("Unoptimized Python is required for exact verification.")
from pathlib import Path
from math import comb,factorial,isqrt
from hashlib import sha256
import argparse,copy,importlib.util,json,struct,sys
HERE=Path(__file__).resolve().parent
PRIMARY=1000003;OTHER=1000033
ROLES={"proof","checker","input","finite_checker","finite_receipt","analytic_proof","lower_two_proof"}
def digest(p):return sha256(p.read_bytes()).hexdigest()
def repo_root(given):
 if given is not None:return given.resolve()
 for p in [HERE,*HERE.parents]:
  if (p/"research/MATH_STATUS.json").is_file():return p
 raise AssertionError("Supply --repo-root for a relocated packet")
def pin_sources(ledger,root):
 assert ledger["schema_version"]==1
 rows=ledger["sources"];assert len(rows)==len(ROLES)and {r["role"]for r in rows}==ROLES
 out={}
 for row in rows:
  assert row["scope"]in{"packet","repo"}
  p=(HERE if row["scope"]=="packet"else root)/row["path"]
  assert p.is_file()and digest(p)==row["sha256"],"Mandatory source pin mismatch: "+str(p)
  out[row["role"]]=p.resolve()
 assert out["checker"]==Path(__file__).resolve()
 return out
def validate_input(v):
 assert v=={"schema_version":1,"support_gap_range":[9,63],"lower_exponent_minimum":3,
  "relative_gap_factor":3,"absolute_terminal_onset":4000,"linear_terminal_factor":128,
  "support_pairs":532,"terminal_cases":3694040,"maximum_terminal":10751,"maximum_factorial":903084,
  "primary_prime":PRIMARY,"alternative_prime":OTHER,
  "finite_source_sha256":"e33fe9e57513b9a0ff74c1e5a5334718e57179450d87cdae2964da27d34ebe32",
  "completed_receipt_sha256":"95bf53dba99f375ea908aedcaa80a36162a9824a0e193840150d81663f6d9036",
  "fallbacks":[[11,41,6142],[4,48,4491],[13,61,3147],[19,61,3843],[19,62,2254]]}
def verify_completed(r,inp,stream=True):
 assert r["status"]=="PASS"
 assert r["expected_support_pairs"]==r["completed_support_pairs"]==inp["support_pairs"]
 assert r["expected_cases"]==r["completed_cases"]==inp["terminal_cases"]
 assert r["max_terminal"]==inp["maximum_terminal"]and r["max_factorial_argument"]==inp["maximum_factorial"]
 assert r["source_sha256"]==inp["finite_source_sha256"]
 assert r["primary_prime"]==PRIMARY and r["alternative_prime"]==OTHER
 assert r["versions"]=={"python":"3.12.3","python_flint":"0.9.0","FLINT":"3.6.0"}
 assert r["actual_polynomials"]=="P_m(t)=sum binom(m,k)*(a*m+1)_(d*k)*t^k, consecutive m=N-1,N, each full degree retained"
 cases=[(a,d,max(4000,128*(a+d))-1)for d in range(9,64)for a in range(3,d//3+1)]
 assert len(cases)==532 and sum(last-1 for a,d,last in cases)==3694040
 assert len(r["rows"])==len(cases)
 observed=[];stream_cases=0
 for row,(a,d,last)in zip(r["rows"],cases):
  assert row["lower_exponent"]==a and row["gap"]==d
  assert row["first_terminal"]==2 and row["last_terminal"]==last and row["tested"]==last-1
  assert row["max_factorial"]==(a+d)*last<=903084<PRIMARY<OTHER
  assert row["primary_prime"]==PRIMARY and row["all_actual_full_degree_pairs_certified"]is True
  failures={}
  for f in row["fallbacks"]:
   N=f["terminal"]
   assert 2<=N<=last and N not in failures
   assert f=={"terminal":N,"primary_gcd_degree":1,"success_prime":OTHER,"full_Bezout":True}
   failures[N]=OTHER;observed.append([a,d,N])
  if stream:
   h=sha256()
   for N in range(2,last+1):h.update(struct.pack("<IIIII",a,d,N,failures.get(N,PRIMARY),0))
   assert h.hexdigest()==row["case_stream_sha256"]
   stream_cases+=last-1
 assert observed==inp["fallbacks"]
 assert max(last for a,d,last in cases)==10751
 assert max((a+d)*last for a,d,last in cases)==903084
 return len(cases),stream_cases
def mod_poly(m,a,d,p):
 # A fresh ratio recurrence, independent of the executed factorial-table builder.
 from flint import nmod_poly
 cs=[1];v=1
 for k in range(m):
  v=v*(m-k)*pow(k+1,-1,p)%p
  for j in range(a*m+d*k+1,a*m+d*(k+1)+1):v=v*j%p
  cs.append(v)
 out=nmod_poly(cs,p);assert out.degree()==m and int(out[0])==1 and int(out[m])!=0
 return out
def verify_fallbacks(inp):
 from flint import nmod_poly
 for p in [PRIMARY,OTHER]:assert all(p%j for j in range(2,isqrt(p)+1))
 records=[]
 for a,d,N in inp["fallbacks"]:
  lo=mod_poly(N-1,a,d,PRIMARY);hi=mod_poly(N,a,d,PRIMARY)
  assert hi.gcd(lo).degree()==1
  lo=mod_poly(N-1,a,d,OTHER);hi=mod_poly(N,a,d,OTHER)
  gg,s,t=hi.xgcd(lo);assert gg==nmod_poly([1],OTHER)and s*hi+t*lo==gg
  records.append({"a":a,"d":d,"N":N,"primary_gcd_degree":1,"alternate_Bezout_verified":True})
 return records
def controls(r,inp,finite_module):
 literal=0
 for a,d in [(3,9),(7,24),(21,63)]:
  for p in [PRIMARY,OTHER]:
   table=finite_module.tables(p,(a+d)*9)
   for n in [1,2,5,8]:
    expected=[comb(n,k)*(factorial(a*n+d*k)//factorial(a*n))%p for k in range(n+1)]
    assert list(map(int,finite_module.polynomial(n,a,d,p,table)))==expected
    assert list(map(int,mod_poly(n,a,d,p)))==expected;literal+=1
 negative=0
 bad=copy.deepcopy(r);bad["rows"].pop()
 try:verify_completed(bad,inp,False)
 except AssertionError:negative+=1
 else:raise AssertionError("Missing support not rejected")
 bad=copy.deepcopy(r);bad["completed_cases"]-=1
 try:verify_completed(bad,inp,False)
 except AssertionError:negative+=1
 else:raise AssertionError("Partial count not rejected")
 bad=copy.deepcopy(r);bad["source_sha256"]="0"*64
 try:verify_completed(bad,inp,False)
 except AssertionError:negative+=1
 else:raise AssertionError("Wrong actual source not rejected")
 bad=copy.deepcopy(r);next(x for x in bad["rows"]if x["fallbacks"])["fallbacks"][0]["success_prime"]=PRIMARY
 try:verify_completed(bad,inp,False)
 except AssertionError:negative+=1
 else:raise AssertionError("Wrong fallback prime not rejected")
 # Check the stream path on a corruption in its first support.
 bad=copy.deepcopy(r);bad["rows"][0]["case_stream_sha256"]="0"*64
 try:verify_completed(bad,inp,True)
 except AssertionError:negative+=1
 else:raise AssertionError("Changed terminal stream not rejected")
 from flint import nmod_poly
 assert nmod_poly([1,1],PRIMARY).gcd(nmod_poly([1,2,1],PRIMARY)).degree()==1
 negative+=1
 return literal,negative
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument("--repo-root",type=Path)
 ap.add_argument("--input",type=Path,default=HERE/"relative_gap_three_completion_inputs.json")
 ap.add_argument("--ledger",type=Path,default=HERE/"relative_gap_three_completion_source_ledger.json")
 ap.add_argument("--output",type=Path,default=HERE/"relative_gap_three_completion_verification.json")
 a=ap.parse_args();root=repo_root(a.repo_root)
 ledger=json.loads(a.ledger.read_text());paths=pin_sources(ledger,root)
 assert a.input.resolve()==paths["input"]
 inp=json.loads(a.input.read_text());validate_input(inp)
 assert digest(paths["finite_checker"])==inp["finite_source_sha256"]
 assert digest(paths["finite_receipt"])==inp["completed_receipt_sha256"]
 result=json.loads(paths["finite_receipt"].read_text());supports,cases=verify_completed(result,inp)
 spec=importlib.util.spec_from_file_location("reviewed_relative_gap_finite_source",paths["finite_checker"])
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 literal,negative=controls(result,inp,module);fallbacks=verify_fallbacks(inp)
 out={"schema_version":1,"status":"PASS","source_ledger_sha256":digest(a.ledger),
  "scope":"Pinned completed finite artifact, all532support/3694040terminal stream entries,24literal coefficient vectors and five freshly reconstructed alternate full Bezout identities; the full3694040GCD run is authenticated, not rerun here. Analytical and all-terminal conclusions use the separate written proofs.",
  "completed_support_pairs":supports,"completed_terminal_stream_entries":cases,"maximum_terminal":10751,
  "maximum_factorial_argument":903084,"literal_coefficient_vectors":literal,
  "fresh_fallback_Bezout_replays":fallbacks,"negative_controls":negative,
  "actual_finite_source_sha256":digest(paths["finite_checker"]),
  "actual_completed_receipt_sha256":digest(paths["finite_receipt"]),
  "no_completed_secondary_run_claim":True}
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","supports":supports,"terminal_stream_entries":cases,
  "fresh_full_Bezout_fallbacks":len(fallbacks),"literal_coefficient_vectors":literal,"negative_controls":negative}))
if __name__=="__main__":main()
