#!/usr/bin/env python3
"""Exhaust pole cancellation for deg p<=4 in the requested Laurent box."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.triangular import TriangularChart,p_catalog,rectangular_support,sparse_weighted_support,polynomialization_system
from jcsearch.laurent import integer_matrix,rank_mod_prime

parser=argparse.ArgumentParser()
parser.add_argument("--deg-p",type=int,default=4);parser.add_argument("--imax",type=int,default=6)
parser.add_argument("--jmin",type=int,default=-8);parser.add_argument("--jmax",type=int,default=8)
parser.add_argument("--sparse",action="store_true");parser.add_argument("--weight",type=int,default=1)
parser.add_argument("--output")
args=parser.parse_args()
support=(sparse_weighted_support(args.imax,args.jmin,args.jmax,args.weight) if args.sparse
         else rectangular_support(args.imax,args.jmin,args.jmax))
primes=(1000003,1000033,1000037)
records=[]
for p in p_catalog(args.deg_p):
    chart=TriangularChart(p); system=polynomialization_system(chart,support)
    matrix,rhs,scale=integer_matrix(system["matrix"],system["rhs"])
    record={"p":sp.sstr(p),"signature":chart.signature,"support_size":len(support),
      "equations":len(system["equations"]),"rank_Q":system["certificate"].rank,
      "kernel_dimension":len(support)-system["certificate"].rank,
      "earliest_forbidden":system["earliest_forbidden"],
      "ranks_mod_p":{str(q):rank_mod_prime(matrix,q) for q in primes}}
    record["bad_primes"]=[q for q in primes if record["ranks_mod_p"][str(q)]!=record["rank_Q"]]
    records.append(record); print(record)
name=args.output or f"results/pole_kernels_p{args.deg_p}_i{args.imax}_j{args.jmin}_{args.jmax}_{'sparse' if args.sparse else 'box'}.json"
path=Path(name);path.parent.mkdir(exist_ok=True)
path.write_text(json.dumps({"scope":{"deg_p":args.deg_p,"i":[0,args.imax],"j":[args.jmin,args.jmax],"sparse":args.sparse,"weight":args.weight},"records":records},indent=2))
print("wrote",path)
