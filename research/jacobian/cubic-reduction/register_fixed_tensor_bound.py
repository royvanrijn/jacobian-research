#!/usr/bin/env python3
"""Register the fixed-tensor derivative bound without changing global minima."""
from pathlib import Path
import copy,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE="jacobian/cubic-reduction/"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())
def save(p,d):(ROOT/p).write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
v=load(BASE+"fixed-tensor-bound-validation.json")
assert v["status"]=="pass" and v["rank"]==34
assert v["source_sha256"]==sha(HERE/"verify_fixed_tensor_bound.py")
assert v["primary_diagnostic_sha256"]==sha(HERE/"fixed-tensor-bound-primary.json")
ledger=load("MATH_STATUS.json");known={e["id"]:e for e in ledger["entries"]}
identifier="KCUB6";assert identifier not in known or known[identifier]["canonical_source"]==BASE+"FIXED_TENSOR_BOUND.md"
entry=copy.deepcopy(known["KCUB5"])
entry.update(id=identifier,title="Fixed padded companion requires at least34 shared cubes",
 scope="For the particular18-output homogeneous cubic tensor P=(tau Q-tau^2 iota W,C1,...,C5,tau^3) used by KCUB5, every simultaneous cube representation P=sum b_a ell_a^3 over any characteristic-zero extension requires at least34 forms. Its324 partial derivatives span a34-dimensional quadratic space (77 nonzero derivatives,51 occurring quadratic monomials). Every such derivative lies in the span of the selected squared forms, proving the lower bound. Independent standard-library Fraction elimination reconstructs the tensor from the scale11=1 base coefficient artifact by diagonal conjugation; its full coefficient hash agrees with a separate SymPy reconstruction and FLINT rank computation. The saved89-form identity supplies the upper bound, so the fixed-tensor count lies in[34,89]. Invertible linear source or output changes preserve this bound. No claim about the minimum dimension of arbitrary cubic-linear Keller maps, different companions or padding choices, ordinary cubic dimensions or three-variable sparsity is made. Warm-started178-form and weighted858-form searches remain bounded misses at89; both cube spans have rank121. A primary-only first Koszul diagnostic also gives34 and is not an independently checked stronger result.",
 canonical_source=BASE+"FIXED_TENSOR_BOUND.md",dependencies=["BCRDB1","KCUB5"],
 checker=BASE+"verify_fixed_tensor_bound.py",artifact_hash="sha256:"+sha(HERE/"verify_fixed_tensor_bound.py"),
 software_lock=[BASE+n for n in ["verify.py","de-bondt12-witness.json","compute_fixed_tensor_bound.py",
 "fixed-tensor-bound-primary.json","fixed-tensor-bound-validation.json","requirements.txt","graph-requirements.txt"]]+["jacobian/de-bondt-2610.24889/verify.py"],
 supersedes_notes=[],supersedes=[],replaced_by=[],closes_problems=[],narrows_problems=[],consumers=[],invalidates_assumptions=[])
if identifier in known:known[identifier].clear();known[identifier].update(entry)
else:ledger["entries"].append(entry);known[identifier]=entry
save("MATH_STATUS.json",ledger)
lessons=load("knowledge/lessons.json")
for r in lessons["lessons"]:
 if r["id"]=="METHOD-KELLER-SHARED-CUBE-BASIS-EXCHANGE-20261114":
  addition=" A warm start from89 made16980 exchanges in60 seconds in178 forms without improvement. Weighted sums with coefficients -2 through2 enlarged the pool to858 forms but preserved rank121;1966 exchanges in90 seconds also stayed at89. KCUB6 independently proves a34-form derivative-space lower bound for the fixed padded tensor only. The first modular Koszul diagnostic gives the same34, not a stronger result."
  if addition not in r["boundary"]:r["boundary"]+=addition
  r["reopen_when"]="Mix factors across the existing term groups, exploit a joint quadratic decomposition, or change the companion beyond its existing linear orbit. Do not repeat either terminal warm-start cap. Retain the verified89 seed and require full matrix, projection, determinant and collision verification before a smaller cubic-linear claim."
  r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"FIXED_TENSOR_BOUND.md"]))
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]))
  r["implementation"]=list(dict.fromkeys(r["implementation"]+[BASE+n for n in ["power_basis_from89.py","power_basis_weighted.py","compute_fixed_tensor_bound.py","verify_fixed_tensor_bound.py","fixed_tensor_koszul.py"]]))
  addition2=" The36-term base gives91 cubes in a50-restart termwise test, but its companion is linearly conjugate to the old tensor and inherits89 cubes. Exact primary bridge checks verify every conjugacy coefficient and the unchanged transported factor product A. Lower monomial count alone is not shared-cube dimension progress."
  if addition2 not in r["boundary"]:r["boundary"]+=addition2
  r["implementation"]=list(dict.fromkeys(r["implementation"]+[BASE+"power36_terms.py",BASE+"power36_bridge.py"]))
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/lessons.json",lessons)
work=load("knowledge/work_items.json")
for r in work["actions"]:
 if r["id"]=="WORK-GLDC":
  addition=" KCUB6 gives the fixed padded tensor shared-cube interval[34,89]. Warm-start178 and weighted858 pools both terminated at89 with cube-span rank121; do not repeat those caps. Next test mixed factor groups or a joint quadratic decomposition, or change the companion; the34 bound does not constrain all Keller maps."
  if addition not in r["next_step"]:r["next_step"]+=addition
  r["claims"]=list(dict.fromkeys(r["claims"]+[identifier]));r["sources"]=list(dict.fromkeys(r["sources"]+[BASE+"FIXED_TENSOR_BOUND.md"]))
  addition2=" The36-term base's companion is in the old tensor's linear orbit; the primary bridge transports the same89-form product. A50-restart termwise91-form output is dominated. Seek changes beyond that linear orbit, or use conjugacy only to enlarge the finite search pool."
  if addition2 not in r["next_step"]:r["next_step"]+=addition2
  r.setdefault("reviewed_claims",{})[identifier]=claim_fingerprint(entry)
save("knowledge/work_items.json",work)
print("Registered KCUB6: fixed tensor shared-cube interval[34,89]; all global minima remain open.")
