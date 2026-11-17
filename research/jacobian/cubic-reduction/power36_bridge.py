#!/usr/bin/env python3
"""Exact primary bridge: the36-term companion lies in the old tensor's GL orbit."""
from pathlib import Path
import hashlib,json
from power_search import cubic12,homogeneous18,sp,X
from power36_terms import make_terms
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
_,K,_=cubic12(2);old,_=homogeneous18(K);old[-1]=X[-1]**3
_,new=make_terms()
V=sp.eye(18);V[8,8]=sp.Rational(1,2);V[8,7]=sp.Rational(3,2);V[10,9]=3
z=sp.Matrix(X);transformed=V.inv()*sp.Matrix(old).subs(dict(zip(X,V*z)),simultaneous=True)
assert all(sp.expand(a-b)==0 for a,b in zip(transformed,new))
w=HERE/"power89/witness.json";d=json.loads(w.read_text())
L=sp.Matrix(d["forms"]);B=sp.Matrix(d["B"]);L[d["padding_index"],:]=sp.zeros(1,18)
assert (L*V)*(V.inv()*B)==sp.Matrix(d["A"])
receipt=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),
 witness89_sha256=sha(w),base36_sha256=sha(HERE/"cubic12-integer-witness.json"),
 source36_sha256=sha(HERE/"power36_terms.py"),
 source_matrix_nonidentity_entries=[dict(row=i,column=j,value=str(V[i,j])) for i in range(18) for j in range(18) if V[i,j]!=int(i==j)],
 determinant_source_matrix=str(V.det()),all18_tensor_conjugacy_identities=True,
 transformed_factor_product_same_A=True,scope="Primary exact bridge only. Linear conjugacy preserves the simultaneous cube minimum. The raw transported factors have the same product A; no new smaller Keller witness is claimed.")
(HERE/"power36-bridge.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
