#!/usr/bin/env python3
"""Stage A: complete coefficient solve for named sparse polynomializable supports."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.msolve import run,available,input_text

x,y=sp.symbols("x y")

def dense_polynomial(degree,prefix):
    monomials=[x**i*y**j for total in range(degree+1) for i in range(total+1) for j in [total-i]]
    coefficients=sp.symbols(f"{prefix}0:{len(monomials)}")
    return sum(c*m for c,m in zip(coefficients,monomials)),coefficients

def system(degree):
    P,a=dense_polynomial(degree,"a");Q,b=dense_polynomial(degree,"b")
    equations=sp.Poly(sp.expand(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x)-1),x,y).coeffs()
    equations += [P.subs({x:0,y:0}),Q.subs({x:0,y:0}),P.subs({x:1,y:0}),Q.subs({x:1,y:0}),
      sp.diff(P,x).subs({x:0,y:0})-1,sp.diff(P,y).subs({x:0,y:0}),
      sp.diff(Q,x).subs({x:0,y:0}),sp.diff(Q,y).subs({x:0,y:0})-1]
    return P,Q,a+b,[sp.expand(e) for e in equations if sp.expand(e)!=0]

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--max-degree",type=int,default=3);args=parser.parse_args()
    if not available(): raise SystemExit("msolve not installed")
    records=[]
    for degree in range(1,args.max_degree+1):
      P,Q,variables,equations=system(degree)
      for prime in (1000003,1000033,1000037,0):
        result=run(equations,variables,prime=prime)
        records.append({"degree":degree,"prime":prime,"equations":len(equations),"variables":len(variables),
          "empty":result.empty,"positive_dimensional":result.positive_dimensional,
          "contains_one":result.contains_one,"returncode":result.returncode})
        print(records[-1])
        if prime==0:
          name="triangle_d3" if degree==3 else f"stage_a_degree_{degree}"
          cert=Path(f"results/certificates/{name}_Q.msolve")
          cert.parent.mkdir(parents=True,exist_ok=True);cert.write_text(result.output)
          cert.with_suffix(".input").write_text(input_text(equations,variables,0))
    path=Path("results/stage_a_sparse.json");path.write_text(json.dumps(records,indent=2));print("wrote",path)

if __name__=="__main__":main()
