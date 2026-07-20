#!/usr/bin/env python3
"""Exact search in the non-toric translated two-divisor chart."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.translated import TranslatedTwoDivisorChart,cancellation_system,kernel_expressions,support_families,x,y,u,v
from jcsearch.laurent import laurent_dict
from jcsearch.msolve import run,input_text

chart=TranslatedTwoDivisorChart();records=[]
assert sp.factor(chart.jacobian+1/(y+1/x))==0
assert chart.pullback(u*v)==x and chart.pullback(u-1/(u*v))==y
assert sp.factor(sp.Matrix([u*v,u-1/(u*v)]).jacobian((u,v)).det()+u)==0

for name,support in support_families().items():
    cancellation=cancellation_system(support)
    basis=kernel_expressions(cancellation);dimension=len(basis)
    aa=sp.symbols(f"a0:{dimension}");bb=sp.symbols(f"b0:{dimension}")
    A=sp.expand(sum(c*f for c,f in zip(aa,basis)));B=sp.expand(sum(c*f for c,f in zip(bb,basis)))
    P,Q=chart.pullback(A),chart.pullback(B)
    assert sp.denom(P)==sp.denom(Q)==1
    residual=sp.expand(sp.Matrix([A,B]).jacobian((u,v)).det()+u)
    equations=list(laurent_dict(residual,(u,v)).values())
    equations += [P.subs({x:0,y:0}),Q.subs({x:0,y:0}),P.subs({x:1,y:0}),Q.subs({x:1,y:0}),
      sp.diff(P,x).subs({x:0,y:0})-1,sp.diff(P,y).subs({x:0,y:0}),
      sp.diff(Q,x).subs({x:0,y:0}),sp.diff(Q,y).subs({x:0,y:0})-1]
    equations=[sp.expand(e) for e in equations if sp.expand(e)!=0];variables=aa+bb
    record={"family":name,"support_size":len(support),"support":support,"kernel_dimension":dimension,
      "pole_rank":cancellation["certificate"].rank,"denominator_orders":cancellation["denominator_orders"],
      "earliest_remainder_monomial":cancellation["earliest_remainder_monomial"],"equations":len(equations),"runs":[]}
    for prime in (1000003,1000033,1000037,0):
      result=run(equations,variables,prime=prime)
      record["runs"].append({"prime":prime,"empty":result.empty,"positive_dimensional":result.positive_dimensional,"returncode":result.returncode})
      if prime==0:
        base=Path(f"results/certificates/translated_{name}_Q")
        base.with_suffix(".input").write_text(input_text(equations,variables,0));base.with_suffix(".msolve").write_text(result.output)
    records.append(record);print(record)
Path("results/stage_c_translated.json").write_text(json.dumps(records,default=sp.sstr,indent=2))
print("wrote results/stage_c_translated.json")
