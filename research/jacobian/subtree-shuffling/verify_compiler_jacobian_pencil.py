#!/usr/bin/env python3
"""Exact retained companion resolvent/profile and all-affine obstruction controls."""
if not __debug__:raise SystemExit("Assertions required")
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,importlib.util,copy
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def validate_input(d):
 assert set(d)=={"schema_version","seed_input","arithmetic_module","mixed_dimension","base_dimension","cubic_rank","tau_index","resolvent_degree_profile","inverse_jacobian_degree","affine_centers","top_times_Q_witness","top_chain_columns"}
 assert d["schema_version"]==1 and [d[k] for k in ("mixed_dimension","base_dimension","cubic_rank","tau_index")]==[17,12,5,17]
 assert d["resolvent_degree_profile"]==[15,14,14,13] and d["inverse_jacobian_degree"]==15
 assert d["top_chain_columns"]==[7,16]
 assert d["affine_centers"]==[[0]*12,[1]*12,list(range(1,13))]
 assert d["top_times_Q_witness"]=={"row":3,"column":0,"inputs":[0]*10+[1]*4+[2]*3,"coefficient":"-288"}
 assert d["seed_input"]=="research/jacobian/subtree-shuffling/cubic_index_seventeen_shuffle_failure_inputs.json"
 assert d["arithmetic_module"]=="research/jacobian/subtree-shuffling/verify_cubic_index_seventeen_shuffle_failure.py"

def main():
 ap=argparse.ArgumentParser();ap.add_argument("--repo-root",type=Path,required=True)
 ap.add_argument("--input",type=Path,default=HERE/"compiler_jacobian_pencil_inputs.json")
 ap.add_argument("--ledger",type=Path,default=HERE/"compiler_jacobian_pencil_source_ledger.json")
 ap.add_argument("--output",type=Path,default=HERE/"compiler_jacobian_pencil_verification.json")
 a=ap.parse_args();repo=a.repo_root.resolve();base=repo/"research/jacobian/subtree-shuffling"
 ledger=json.loads(a.ledger.read_text());assert ledger["schema_version"]==1
 roles={"proof","checker","input","seed_input","arithmetic_module","companion_proof","index17_proof"}
 assert len(ledger["sources"])==len(roles) and {e["role"] for e in ledger["sources"]}==roles
 paths={}
 for e in ledger["sources"]:
  root=a.ledger.parent if e["scope"]=="packet" else repo
  path=(root/e["path"]).resolve()
  assert path.is_file() and sha(path)==e["sha256"],"Source pin mismatch: "+str(path)
  paths[e["role"]]=path
 assert paths["checker"]==Path(__file__).resolve() and paths["input"]==a.input.resolve()
 inp=json.loads(a.input.read_text());validate_input(inp)
 assert paths["seed_input"]==(repo/inp["seed_input"]).resolve()
 assert paths["arithmetic_module"]==(repo/inp["arithmetic_module"]).resolve()
 data=base/"cubic_index_seventeen_shuffle_failure_inputs.json"
 module=base/"verify_cubic_index_seventeen_shuffle_failure.py"
 assert hashlib.sha256(data.read_bytes()).hexdigest()=="f305ecd8db14df064ff706c9c74242b8d004024d0ef25747c6f2334a3f98e26b"
 assert hashlib.sha256(module.read_bytes()).hexdigest()=="4aed8de2304ea9d4df6c25208ecf8cbe3dad5332843e92d8947c84ecd9d11b49"
 spec=importlib.util.spec_from_file_location("pinned_polynomial_ops",module)
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 d=json.loads(data.read_text());J=m.jacobian(d["H_terms"])
 def term(*i,c=1):return {tuple(sorted(i)):F(c)}
 def scale(p,c):return m.clean({x:y*c for x,y in p.items()})
 def summat(*args):
  r={}
  for mat in args:
   for ij,p in mat.items():r[ij]=m.add(r.get(ij,{}),p)
  return {ij:p for ij,p in r.items() if p}
 def mscale(A,p):return {ij:m.mul(q,p) for ij,q in A.items() if m.mul(q,p)}
 def identity(n):return {(i,i):term() for i in range(n)}
 def tauone(p):
  r={}
  for mon,c in p.items():
   key=tuple(i for i in mon if i!=17);r[key]=r.get(key,F(0))+c
  return m.clean(r)
 A={ij:tauone(p) for ij,p in J.items() if max(ij)<17}
 powers=[identity(17)]
 for k in range(1,17):powers.append(m.mmul(powers[-1],A))
 assert powers[15] and not powers[16]
 assert len(powers[15])==10 and {j for i,j in powers[15]}==set(inp["top_chain_columns"])
 assert {j:p for (i,j),p in A.items() if i==7}=={1:term(1,c=2)}
 assert {j:p for (i,j),p in A.items() if i==16}=={1:term(1,1,c=3)}
 S=[{ij:scale(p,(-1)**k) for ij,p in P.items() if max(ij)<12} for k,P in enumerate(powers[:16])]
 Q={}
 C={}
 for t in d["H_terms"]:
  i=t["output"];mon=tuple(t["inputs"]);co=F(t["coefficient"])
  if i<12 and mon.count(17)==1:Q[(i,0)]=m.add(Q.get((i,0),{}),tauone({mon:co}))
  if 12<=i<17:C[i-12]=m.add(C.get(i-12,{}),{mon:co})
 Cterms=[{"output":i,"inputs":list(mon),"coefficient":str(c)} for i,p in C.items() for mon,c in p.items()]
 DC=m.jacobian(Cterms);Iota={(i,i):term() for i in range(5)}
 DCfull={(i,j):p for (i,j),p in DC.items()}
 SI=[m.mmul(s,Iota) for s in S]
 DCS=[m.mmul(DC,s) for s in S]
 DCSI=[m.mmul(z,Iota) for z in DCS]
 def deg(series):return max(k for k,q in enumerate(series) if q)
 assert [deg(S),deg(SI),deg(DCS),deg(DCSI)]==[15,14,14,13]
 lam=18
 R={}
 for k in range(16):
  factor=term(*([lam]*k))
  for (i,j),p in S[k].items():R[(i,j)]=m.add(R.get((i,j),{}),m.mul(p,factor))
  for (i,j),p in SI[k].items():R[(i,j+12)]=m.add(R.get((i,j+12),{}),m.mul(p,term(*([lam]*(k+1)))))
  for (i,j),p in DCS[k].items():R[(i+12,j)]=m.add(R.get((i+12,j),{}),scale(m.mul(p,term(*([lam]*(k+1)))),-1))
  for (i,j),p in DCSI[k].items():R[(i+12,j+12)]=m.add(R.get((i+12,j+12),{}),scale(m.mul(p,term(*([lam]*(k+2)))),-1))
 R=summat(R,{(i+12,i+12):term() for i in range(5)})
 L=summat(identity(17),mscale(A,term(lam)))
 assert m.mmul(L,R)==m.mmul(R,L)==identity(17)
 wrong_resolvent={ij:(scale(p,-1) if ij[0]<12 and ij[1]>=12 else p) for ij,p in R.items()}
 assert m.mmul(L,wrong_resolvent)!=identity(17)
 for (i,j),p in J.items():
  if max(i,j)<17:
   weight=1 if i<12 and j<12 else 2 if i<12 else 0
   assert p==m.mul(A[(i,j)],term(*([17]*weight)))
 b={i:p for (i,j),p in J.items() if j==17}
 expected={i:p for (i,z),p in Q.items()}
 for i in range(5):expected[i]=m.add(expected.get(i,{}),term(17,12+i,c=-2))
 assert {i:p for i,p in expected.items() if p}==b
 full=J
 for k in range(2,18):
  full=m.mmul(full,J)
  if k==16:
   assert len(full)==5 and sum(map(len,full.values()))==10 and all(j==17 for i,j in full)
 assert not full
 Nd=S[15]
 assert Nd and not m.mmul(Nd,DCfull) and not m.mmul(Nd,Iota)
 NdQ=m.mmul(Nd,Q)
 witness=(0,)*10+(1,)*4+(2,)*3
 assert NdQ[(3,0)]=={witness:F(-288)}
 assert {j for i,j in Nd}=={7}
 # The highest chain cannot be cancelled by ANY constanttau^3 target shear.
 assert max(len(mon) for (i,j),p in powers[15].items() if i==3 for mon in p)<=15
 assert len(witness)==17
 # Exact affine-normalization identities at several centers, plus NdDC=0
 # supplies the coefficient identity for ALL centers in the written proof.
 JK=summat(identity(12),{ij:p for ij,p in m.jacobian(
  [{"output":i,"inputs":list(mon),"coefficient":str(c)} for (i,z),p in Q.items() for mon,c in p.items()]).items()},DCfull)
 def inv(M):
  n=len(M);X=[list(map(F,row))+[F(i==j) for j in range(n)] for i,row in enumerate(M)]
  for i in range(n):
   k=next(k for k in range(i,n) if X[k][i]);X[i],X[k]=X[k],X[i]
   t=X[i][i];X[i]=[v/t for v in X[i]]
   for k in range(n):
    if k!=i:
     t=X[k][i]
     if t:X[k]=[v-t*w for v,w in zip(X[k],X[i])]
  return [r[n:] for r in X]
 wrong_order_controls=0
 for center in inp["affine_centers"]:
  E=[[m.ev(JK.get((i,j),{}),center) for j in range(12)] for i in range(12)]
  EI=inv(E)
  Epoly={(i,j):term(c=c) for i,row in enumerate(E) for j,c in enumerate(row) if c}
  Lpoly={(i,j):term(c=c) for i,row in enumerate(EI) for j,c in enumerate(row) if c}
  avec={(i,0):term(c=c) for i,c in enumerate(center) if c}
  Qshift=summat(Q,m.mmul(DCfull,avec))
  Qa=m.mmul(Lpoly,Qshift);Nda=m.mmul(Nd,Epoly);Ia=m.mmul(Lpoly,Iota)
  assert m.mmul(Nda,Qa)==NdQ and m.mmul(Nda,Ia)==m.mmul(Nd,Iota)
  if m.mmul(m.mmul(Epoly,Nd),Qa)!=NdQ:wrong_order_controls+=1
 assert wrong_order_controls>=1
 # Universal two-shear product splitting, including its complete collision lift.
 Fi,f,g,U,V=[term(i) for i in range(5)]
 after=m.add(Fi,scale(m.mul(m.add(U,scale(f,-1)),m.add(V,scale(g,-1))),-1))
 formula=m.add(Fi,scale(m.mul(f,g),-1),scale(m.mul(U,V),-1),m.mul(U,g),m.mul(V,f))
 assert after==formula
 lifted={}
 for mon,c in after.items():
  key=tuple(sorted(1 if i==3 else 2 if i==4 else i for i in mon))
  lifted[key]=lifted.get(key,F(0))+c
 assert m.clean(lifted)==Fi
 rejected=0
 for name in ["degree","column","rank","selector"]:
  bad=copy.deepcopy(inp)
  if name=="degree":bad["inverse_jacobian_degree"]=14
  elif name=="column":bad["top_chain_columns"]=[7,15]
  elif name=="rank":bad["cubic_rank"]=6
  else:bad["seed_input"]="missing.json"
  try:validate_input(bad)
  except (AssertionError,ValueError):rejected+=1
 assert rejected==4
 result={"status":"PASS","source_ledger_sha256":sha(a.ledger),"mixed_dimension":17,"full_dimension":18,
  "inverse_jacobian_degree":15,"resolvent_block_series_degrees":[15,14,14,13],
  "mixed_companion_exact_index":16,"full_exact_index":17,
  "J16_only_tau_column":True,"J16_entries":5,"J16_monomials":10,
  "top_inverse_jacobian_columns":[7],"top_times_cubic_inclusion_zero":True,
  "top_times_DC_zero":True,"top_times_Q_row3":{"coefficient":"-288","inputs":list(witness)},
  "all_constant_target_tau3_shears_obstructed_by_degree":True,
  "affine_center_controls":3,"both_resolvent_inverse_orders_verified":True,
  "A15_entries":10,"A15_support_columns":inp["top_chain_columns"],"complete_carrier_rows_verified":2,
  "wrong_affine_matrix_order_controls":wrong_order_controls,"wrong_resolvent_sign_rejected":True,
  "universal_product_splitting_and_lift_verified":True,"malformed_inputs_rejected":rejected}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+"\n")
 print(json.dumps(result))
if __name__=="__main__":main()
