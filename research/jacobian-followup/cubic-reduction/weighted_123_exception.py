#!/usr/bin/env python3
"""Retain and test the sole nonzero root of the lifted boundary denominator."""
from pathlib import Path
import hashlib,json
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=HERE/'weighted-sparse7/1_-2_-3/slices/both/boundary/input.json';d=json.loads(source.read_text());P=sp.Symbol('c3');W=sp.symbols(' '.join(d['variables']),seq=True);eqs=[sp.expand(sum(sp.Rational(c)*sp.prod(v**a for v,a in zip((P,)+W,e)) for e,c in terms).subs(P,-sp.Rational(810,2021))) for terms in d['equations']];eqs=[f for f in eqs if f];out=source.parent.parent/'exception';out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,W)),specialization={'c3':'-810/2021'},equations=[[[list(e),str(c)] for e,c in sp.Poly(f,*W).terms() if c] for f in eqs],scope='Exact rational nonzero denominator root. A separate unit replay is required; generic unit alone does not exclude it.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(out)
