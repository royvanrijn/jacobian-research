#!/usr/bin/env python3
"""Bounded modular Koszul-flattening diagnostic of the fixed derivative system."""
from pathlib import Path
from itertools import combinations
from math import ceil
import hashlib,json,time
from flint import nmod_mat
from power_search import cubic12,homogeneous18,X,sp
HERE=Path(__file__).resolve().parent
v=HERE/"fixed-tensor-bound-validation.json";d=json.loads(v.read_text());assert d["status"]=="pass"
_,K,_=cubic12(2);H,_=homogeneous18(K);H[-1]=X[-1]**3
labels=d["independent_derivative_rows"];qs=[sp.Poly(sp.diff(H[i],X[j]),*X) for i,j in labels]
pairs=list(combinations(range(18),2));pair_id={p:i for i,p in enumerate(pairs)}
M=[[0]*(18*len(qs)) for _ in range(len(pairs)*18)]
for j,q in enumerate(qs):
    entries=[]
    for e,c in q.terms():
        if not c:continue
        ids=[k for k,power in enumerate(e) for _ in range(power)];assert len(ids)==2
        i,k=ids
        entries.append((i,k,int(c)*(2 if i==k else 1)))
        if i!=k:entries.append((k,i,int(c)))
    for a in range(18):
        for i,k,c in entries:
            if a==i:continue
            pair=tuple(sorted((a,i)));sign=1 if a<i else -1
            M[pair_id[pair]*18+k][j*18+a]+=sign*c
prime=1000003;t=time.monotonic();rank=nmod_mat(M,prime).rank()
receipt=dict(status="modular-lower-bound-diagnostic",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 input_sha256=hashlib.sha256(v.read_bytes()).hexdigest(),prime=prime,rank_mod_prime=rank,
 rows=len(M),columns=len(M[0]),per_cube_rank=17,implied_cube_lower_bound=ceil(rank/17),seconds=round(time.monotonic()-t,3),
 scope="Primary-only diagnostic. A single shared cube has Koszul rank17, so modular rank divided by17 bounds the fixed tensor cube count. No independent certificate or global Keller bound.")
(HERE/"fixed-tensor-koszul.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
