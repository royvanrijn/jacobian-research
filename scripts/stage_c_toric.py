#!/usr/bin/env python3
"""First bounded two-divisor toric chart search with exact F4 certificates."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.toric import catalog,cone_support,x,y,u,v
from jcsearch.msolve import run,input_text
from jcsearch.laurent import laurent_dict

def build(chart,support):
    aa=sp.symbols(f"a0:{len(support)}");bb=sp.symbols(f"b0:{len(support)}")
    A=sum(c*u**i*v**j for c,(i,j) in zip(aa,support));B=sum(c*u**i*v**j for c,(i,j) in zip(bb,support))
    P,Q=sp.expand(chart.pullback(A)),sp.expand(chart.pullback(B))
    assert sp.denom(P)==sp.denom(Q)==1
    residual=sp.expand(sp.Matrix([A,B]).jacobian((u,v)).det()-chart.reciprocal_outer_jacobian)
    equations=list(laurent_dict(residual,(u,v)).values())
    equations += [P.subs({x:0,y:0}),Q.subs({x:0,y:0}),P.subs({x:1,y:0}),Q.subs({x:1,y:0}),
      sp.diff(P,x).subs({x:0,y:0})-1,sp.diff(P,y).subs({x:0,y:0}),
      sp.diff(Q,x).subs({x:0,y:0}),sp.diff(Q,y).subs({x:0,y:0})-1]
    return A,B,P,Q,aa+bb,[sp.expand(e) for e in equations if sp.expand(e)!=0]

records=[]
for index,chart in enumerate(catalog()):
  for limit in (10,15):
    support=cone_support(chart,limit=limit);A,B,P,Q,variables,equations=build(chart,support)
    record={"chart":chart.expressions,"signature":chart.signature,"support":support,"terms_per_coordinate":len(support),"runs":[]}
    for prime in (1000003,1000033,1000037,0):
      result=run(equations,variables,prime=prime)
      record["runs"].append({"prime":prime,"empty":result.empty,"positive_dimensional":result.positive_dimensional,"returncode":result.returncode})
      if prime==0:
        base=Path(f"results/certificates/stage_c_chart{index}_n{limit}_Q")
        base.with_suffix(".input").write_text(input_text(equations,variables,0));base.with_suffix(".msolve").write_text(result.output)
    records.append(record);print(record)
Path("results/stage_c_toric.json").write_text(json.dumps(records,default=sp.sstr,indent=2))
print("wrote results/stage_c_toric.json")
