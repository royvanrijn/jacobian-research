#!/usr/bin/env python3
"""Reproduce the rejected formal-root cache collision and verify the repaired domains."""
from pathlib import Path
import hashlib,json
import sympy as sp
HERE=Path(__file__).resolve().parent
def monic_modulus(K):
    coefficients=K.mod.to_list()
    return [sp.Rational(c/coefficients[0]) for c in coefficients]
def main():
    w,t,a=sp.symbols("w t a")
    bad1=sp.QQ.algebraic_field(sp.AlgebraicNumber((sp.Poly(w-2,w),a)))
    bad2=sp.QQ.algebraic_field(sp.AlgebraicNumber((sp.Poly(w-14,w),a)))
    assert bad1==bad2 and monic_modulus(bad2)!=[1,-14]
    fixed=[]
    for value in (2,14):
        alpha=sp.Dummy("alpha")
        K=sp.QQ.algebraic_field(sp.AlgebraicNumber((sp.Poly(w-value,w),alpha)))
        assert monic_modulus(K)==[1,-value]
        ring=K.poly_ring(t).ring
        assert monic_modulus(ring.domain)==[1,-value]
        theta=K.unit*K.one
        assert theta-K.convert(value)==K.zero
        fixed.append(K)
    assert fixed[0]!=fixed[1]
    receipt=dict(status="pass",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 sympy=sp.__version__,moduli=["w-2","w-14"],
                 rejected_shared_formal_root_alias=True,
                 distinct_generators_and_moduli_verified=True,
                 boundary="This diagnoses the discarded exploratory specialization implementation; predecessor sparse proofs used no such number-field construction.")
    (HERE/"field-regression.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

