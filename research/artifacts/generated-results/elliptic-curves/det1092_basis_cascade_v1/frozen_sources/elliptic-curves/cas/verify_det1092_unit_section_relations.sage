#!/usr/bin/env sage-python
"""At most256 signed finite-kernel words; exact rational group identities only."""
import json,itertools,hashlib
from pathlib import Path
from sage.all import QQ,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/det1092-unit-section-relations-v1'
p=json.loads((D/'protocol.json').read_bytes());assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==p['script_sha256']
def read(n):
 path=ROOT/n;assert hashlib.sha256(path.read_bytes()).hexdigest()==p['inputs'][n];return json.loads(path.read_bytes())
seed=read(p['seed']);diag=read(p['diagnosis']);E=EllipticCurve(QQ,seed['curve']);points=[E(P) for P in seed['points']];relations=[];attempts=0
for base in diag['kernel_words']:
 ids=[i for i,v in enumerate(base) if v];assert len(ids)<=8
 found=[]
 for signs in itertools.product([-1,1],repeat=len(ids)-1):
  word=[0]*17
  for i,v in zip(ids,[1,*signs]):word[i]=v
  attempts+=1;assert attempts<=256
  if sum((word[i]*points[i] for i in ids),E(0)).is_zero():found.append(word)
 assert len(found)==1;relations.append(found[0])
assert len(relations)==2 and matrix(QQ,relations).rank()==2
out=dict(status='PASS',curve=seed['curve'],points=seed['points'],relations=relations,attempted_signed_words=attempts,independent_indices=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,16],section_span_rank=15,boundary='Two independent exact integral relations bound the specialized17-section span above by15. Independent finite lower bound is supplied by the separate seed and cloud certificates. No upper bound on the whole fibre rank.')
for j in [14,15]:assert sum(abs(r[j]) for r in relations)==1
assert not (D/'result.json').exists();(D/'result.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',attempts,'signed words; relations',relations,flush=True)
