#!/usr/bin/env python3
"""Exact filters and a small enumerator for rational C^2 -> C^2 pairs.

This is search infrastructure, not a claim of a 2D counterexample. Candidate
pairs can be added to `candidate_stream`; every survivor is checked exactly.
"""
from dataclasses import dataclass
from typing import Iterable
import sympy as sp

x, y, u, v, s, tau = sp.symbols("x y u v s tau")

@dataclass(frozen=True)
class Pair:
    name: str
    A: tuple[sp.Expr, sp.Expr]
    B: tuple[sp.Expr, sp.Expr]

def jac(coords, variables):
    return sp.factor(sp.Matrix(coords).jacobian(variables).det())

def compose(pair: Pair):
    sub = {u: pair.A[0], v: pair.A[1]}
    return tuple(sp.cancel(f.subs(sub, simultaneous=True)) for f in pair.B)

def is_polynomial(expr, variables=(x, y)):
    num, den = sp.fraction(sp.cancel(expr))
    return sp.expand(den) in (1, -1) and sp.Poly(num/den, *variables) is not None

def analyze(pair: Pair):
    C = compose(pair)
    ja = jac(pair.A, (x, y))
    jb = jac(pair.B, (u, v))
    chain = sp.factor(ja * jb.subs({u: pair.A[0], v: pair.A[1]}, simultaneous=True))
    polynomial = all(is_polynomial(f) for f in C)
    constant_keller = polynomial and chain.is_number and chain != 0
    return {"name": pair.name, "A": pair.A, "B": pair.B, "composition": C,
            "jac_A": ja, "jac_B_after_A": sp.factor(jb.subs({u: pair.A[0], v: pair.A[1]}, simultaneous=True)),
            "chain_jacobian": chain, "polynomial": polynomial,
            "constant_keller": constant_keller}

def branch_limit(composition, branch: tuple[sp.Expr, sp.Expr]):
    """Limit along a branch (x(s,tau),y(s,tau)) as s -> 0."""
    return tuple(sp.simplify(sp.limit(f.subs({x: branch[0], y: branch[1]}, simultaneous=True), s, 0))
                 for f in composition)

def branches_identified(composition, branch1, branch2):
    l1, l2 = branch_limit(composition, branch1), branch_limit(composition, branch2)
    distinct = any(sp.simplify(a-b) != 0 for a, b in zip(branch1, branch2))
    finite = all(not value.has(sp.oo, -sp.oo, sp.zoo, sp.nan) for value in l1+l2)
    same = all(sp.simplify(a-b) == 0 for a, b in zip(l1, l2))
    return distinct and finite and same, l1, l2

def candidate_stream(max_m: int = 3) -> Iterable[Pair]:
    # Positive controls: poles and reciprocal Jacobians cancel, but the result
    # is a triangular automorphism, hence never a counterexample.
    for m in range(1, max_m + 1):
        for h in (0, u, u**2, u**3-u):
            yield Pair(f"triangular_m{m}_h{sp.sstr(h)}",
                       (x, y/x**m), (u, u**m*v + h))

    # A deliberately failing candidate exercises rejection of nonconstant J.
    yield Pair("reject_nonconstant", (x**2, y/x), (u, u*v))

def main():
    survivors = []
    for pair in candidate_stream():
        result = analyze(pair)
        status = "SURVIVES" if result["constant_keller"] else "reject"
        print(f"{status:8} {pair.name}: comp={result['composition']}, J={result['chain_jacobian']}")
        if result["constant_keller"]:
            survivors.append(result)
    assert survivors, "positive controls failed"
    assert all(r["composition"][0] == x for r in survivors)
    print(f"\n{len(survivors)} cancellation controls survive; all are triangular automorphisms.")
    print("No 2D counterexample is claimed. Extend candidate_stream with sparse ansatz solutions.")

if __name__ == "__main__":
    main()

