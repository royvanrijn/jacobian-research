#!/usr/bin/env python3
"""Stratify pole-matrix ranks as every lower coefficient of p varies mod q."""
import argparse,itertools,json,math,random,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from jcsearch.triangular import rectangular_support

def convolution(a,b,q):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
      for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%q
    return out

def rank_mod(rows,ncols,q):
    rows=[dict((j,v%q) for j,v in row.items() if v%q) for row in rows]
    rank=0
    for col in range(ncols):
      pivot=next((i for i in range(rank,len(rows)) if rows[i].get(col,0)),None)
      if pivot is None:continue
      rows[rank],rows[pivot]=rows[pivot],rows[rank]
      inv=pow(rows[rank][col],-1,q);rows[rank]={j:(v*inv)%q for j,v in rows[rank].items()}
      for i in range(len(rows)):
        factor=rows[i].get(col,0)
        if i!=rank and factor:
          for j,v in rows[rank].items():
            rows[i][j]=(rows[i].get(j,0)-factor*v)%q
            if not rows[i][j]:rows[i].pop(j,None)
      rank+=1
      if rank==len(rows):break
    return rank

def cancellation_rank(coefficients,q,support):
    imax=max(i for i,_ in support);powers=[[1]]
    for _ in range(imax):powers.append(convolution(powers[-1],coefficients,q))
    rowdata={}
    for col,(i,j) in enumerate(support):
      for k in range(i+1):
        factor=math.comb(i,k)%q
        for ell,value in enumerate(powers[i-k]):
          exponent=ell+j
          if exponent>0 and value:
            rowdata.setdefault((k,exponent),{})[col]=(factor*value)%q
    return rank_mod(list(rowdata.values()),len(support),q)

parser=argparse.ArgumentParser();parser.add_argument("--random-large",type=int,default=200);args=parser.parse_args()
support=rectangular_support(6,-8,8);records=[]
for q in (3,5):
  for degree in range(1,5):
    ranks=Counter()
    for lower in itertools.product(range(q),repeat=degree):
      ranks[cancellation_rank(list(lower)+[1],q,support)]+=1
    records.append({"field":q,"degree":degree,"samples":q**degree,"rank_distribution":dict(ranks)})
    print(records[-1])
rng=random.Random(20260720)
for q in (1000003,1000033,1000037):
  for degree in range(1,5):
    ranks=Counter()
    for _ in range(args.random_large):
      lower=[rng.randrange(q) for _ in range(degree)]
      ranks[cancellation_rank(lower+[1],q,support)]+=1
    records.append({"field":q,"degree":degree,"samples":args.random_large,"rank_distribution":dict(ranks)})
    print(records[-1])
Path("results/parameter_rank_scan.json").write_text(json.dumps(records,indent=2))

