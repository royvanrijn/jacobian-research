#!/usr/bin/env python3
"""Primary FLINT/SymPy derivative rank of the fixed padded companion tensor."""
from pathlib import Path
import hashlib,json
from power_search import cubic12,homogeneous18,X,sp,fmpq_mat,DB
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    _,K,_=cubic12(2);H,_=homogeneous18(K);H[-1]=X[-1]**3
    tensor=[{e:c for e,c in sp.Poly(h,*X).terms() if c} for h in H]
    tensor_hash=hashlib.sha256(json.dumps([[[list(e),str(c)] for e,c in sorted(p.items())] for p in tensor],separators=(",",":")).encode()).hexdigest()
    derivatives=[{e:c for e,c in sp.Poly(sp.diff(h,x),*X).terms() if c} for h in H for x in X]
    monomials=sorted(set().union(*(set(d) for d in derivatives)))
    M=fmpq_mat([[str(d.get(e,0)) for e in monomials] for d in derivatives]);rank=M.rank()
    receipt=dict(status="primary-exact-pass",source_sha256=sha(Path(__file__)),
       base_source_sha256=sha(DB/"verify.py"),tensor_sha256=tensor_hash,rank=rank,derivatives=len(derivatives),
       nonzero_derivatives=sum(bool(d) for d in derivatives),ambient_quadratic_monomials=len(monomials),
       scope="Derivative-space rank for the fixed padded18-output tensor; independent replay is separate.")
    (HERE/"fixed-tensor-bound-primary.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()
