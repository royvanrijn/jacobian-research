#!/usr/bin/env python3
"""Exhaust named 10--30 coefficient sparse families with F4/F4SAT."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.triangular import TriangularChart,p_catalog,x,y
from jcsearch.msolve import run,input_text

def triangle(d): return tuple(x**i*y**j for total in range(d+1) for i in range(total+1) for j in [total-i])

def families():
    pool=triangle(6)
    required=(sp.Integer(1),x,y)
    def select(key,n=15):
      chosen=list(required)
      for monomial in sorted(pool,key=key):
        if monomial not in chosen: chosen.append(monomial)
        if len(chosen)==n: break
      return tuple(chosen)
    return {
      "triangle_d3":triangle(3),                 # 10 terms/coordinate
      "triangle_d4":triangle(4),                 # 15 terms/coordinate
      "x_heavy_15":select(lambda m:(sp.degree(m,y)*3+sp.degree(m,x),sp.degree(m,y))),
      "y_heavy_15":select(lambda m:(sp.degree(m,x)*3+sp.degree(m,y),sp.degree(m,x))),
      "balanced_15":select(lambda m:(abs(sp.degree(m,x)-sp.degree(m,y)),sp.total_degree(m))),
    }

def build(monomials):
    a=sp.symbols(f"a0:{len(monomials)}");b=sp.symbols(f"b0:{len(monomials)}")
    P=sum(c*m for c,m in zip(a,monomials));Q=sum(c*m for c,m in zip(b,monomials))
    equations=sp.Poly(sp.expand(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x)-1),x,y).coeffs()
    equations += [P.subs({x:0,y:0}),Q.subs({x:0,y:0}),P.subs({x:1,y:0}),Q.subs({x:1,y:0}),
      sp.diff(P,x).subs({x:0,y:0})-1,sp.diff(P,y).subs({x:0,y:0}),
      sp.diff(Q,x).subs({x:0,y:0}),sp.diff(Q,y).subs({x:0,y:0})-1]
    equations=[sp.expand(e) for e in equations if sp.expand(e)!=0]
    top=max(map(sp.total_degree,monomials))
    leading=[c for c,m in zip(a+b,monomials+monomials) if sp.total_degree(m)==top]
    saturation=sp.prod(leading)
    return P,Q,a+b,equations,saturation

records=[]
for name,monomials in families().items():
    P,Q,variables,equations,saturation=build(monomials)
    chart_realizations={}
    for p in p_catalog(4):
      chart=TriangularChart(p)
      terms=set()
      for m in monomials:
        terms.update(sp.Add.make_args(sp.expand(chart.to_chart(m))))
      chart_realizations[sp.sstr(p)]=len(terms)
    record={"family":name,"terms_per_coordinate":len(monomials),"coefficient_variables":len(variables),
            "equations":len(equations),"chart_term_counts":chart_realizations,"runs":[]}
    for prime in (1000003,1000033,1000037):
      result=run(equations,variables,prime=prime,saturate=saturation)
      record["runs"].append({"prime":prime,"saturated":True,"empty":result.empty,
        "positive_dimensional":result.positive_dimensional,"returncode":result.returncode})
    exact=run(equations,variables,prime=0)
    record["runs"].append({"prime":0,"saturated":False,"empty":exact.empty,
      "positive_dimensional":exact.positive_dimensional,"returncode":exact.returncode})
    cert=Path(f"results/certificates/{name}_Q.msolve");cert.parent.mkdir(parents=True,exist_ok=True);cert.write_text(exact.output)
    cert.with_suffix(".input").write_text(input_text(equations,variables,0))
    records.append(record);print(record)
path=Path("results/stage_a_support_families.json");path.write_text(json.dumps(records,indent=2));print("wrote",path)
