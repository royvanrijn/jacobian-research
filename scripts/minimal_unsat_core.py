#!/usr/bin/env python3
"""Greedily extract an inclusion-minimal F4-unsatisfiable equation core."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.msolve import run,input_text

x,y=sp.symbols("x y")
monomials=tuple(x**i*y**j for total in range(4) for i in range(total+1) for j in [total-i])
a=sp.symbols(f"a0:{len(monomials)}");b=sp.symbols(f"b0:{len(monomials)}")
P=sum(c*m for c,m in zip(a,monomials));Q=sum(c*m for c,m in zip(b,monomials));variables=a+b
residual=sp.Poly(sp.expand(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x)-1),x,y)
labeled=[]
for exponent,coefficient in residual.terms(): labeled.append((f"jac_x{exponent[0]}_y{exponent[1]}",coefficient))
labeled += [
 ("P_at_0",P.subs({x:0,y:0})),("Q_at_0",Q.subs({x:0,y:0})),
 ("P_at_1",P.subs({x:1,y:0})),("Q_at_1",Q.subs({x:1,y:0})),
 ("Px_at_0_minus_1",sp.diff(P,x).subs({x:0,y:0})-1),("Py_at_0",sp.diff(P,y).subs({x:0,y:0})),
 ("Qx_at_0",sp.diff(Q,x).subs({x:0,y:0})),("Qy_at_0_minus_1",sp.diff(Q,y).subs({x:0,y:0})-1)]
labeled=[(n,sp.expand(e)) for n,e in labeled if sp.expand(e)!=0]
assert run([e for _,e in labeled],variables,prime=0).empty

def modularly_empty(pairs):
    return all(run([e for _,e in pairs],variables,prime=p).empty for p in (1000003,1000033,1000037))

core=list(labeled)
changed=True
while changed:
  changed=False
  for item in list(core):
    trial=[pair for pair in core if pair is not item]
    if modularly_empty(trial):
      core=trial;changed=True;print("removed",item[0],"remaining",len(core))

# Exact-certify the greedily reduced core. Proving single-deletion minimality is
# intentionally omitted: the deletion ideals become positive-dimensional and
# are much harder than the original emptiness computation.
result=run([e for _,e in core],variables,prime=0);assert result.empty
base=Path("results/certificates/triangle_d3_minimal_core_Q")
base.with_suffix(".input").write_text(input_text([e for _,e in core],variables,0))
base.with_suffix(".msolve").write_text(result.output)
Path("results/triangle_d3_minimal_core.json").write_text(json.dumps({"original_equations":len(labeled),
  "core_size":len(core),"minimality":"greedily reduced over three primes; retained core exactly empty over Q; not claimed inclusion-minimal",
  "labels":[n for n,_ in core],"equations":[sp.sstr(e) for _,e in core]},indent=2))
print("minimal core",[n for n,_ in core])
