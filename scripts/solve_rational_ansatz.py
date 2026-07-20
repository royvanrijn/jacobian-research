#!/usr/bin/env python3
"""Solve a small Laurent ansatz for reciprocal-Jacobian rational pairs.

The default ansatz is intentionally small enough for exact symbolic solving:

  A=(x, (a0+a1*x+a2*y+a3*x*y)/x)
  B=(u, (b0+b1*u+b2*v+b3*u*v)/u)

We impose (i) exact polynomiality of B o A and (ii) constant nonzero chain
Jacobian k. Solutions are then classified; invertible triangular solutions are
expected controls, not counterexamples. Enlarge supports/pole orders to search.
"""
import argparse
import sympy as sp

x, y, u, v = sp.symbols("x y u v")
a0, a1, a2, a3, b0, b1, b2, b3, k = sp.symbols(
    "a0 a1 a2 a3 b0 b1 b2 b3 k")

def numerator_coefficients(expr):
    num, den = sp.fraction(sp.cancel(expr))
    return sp.Poly(sp.expand(num), x, y).coeffs(), sp.factor(den)

def equations():
    A = (x, (a0+a1*x+a2*y+a3*x*y)/x)
    B = (u, (b0+b1*u+b2*v+b3*u*v)/u)
    sub = {u: A[0], v: A[1]}
    C = tuple(sp.cancel(f.subs(sub, simultaneous=True)) for f in B)

    eqs = []
    # Polynomiality: after putting over a monomial x^d, coefficients of terms
    # below degree d in x must vanish.
    for f in C:
        num, den = sp.fraction(f)
        den_poly = sp.Poly(den, x, y)
        if len(den_poly.terms()) != 1 or den_poly.degree(y) != 0:
            raise ValueError(f"ansatz produced unsupported denominator {den}")
        pole_order = den_poly.degree(x)
        poly = sp.Poly(sp.expand(num), x, y)
        for (ix, _iy), coeff in poly.terms():
            if ix < pole_order:
                eqs.append(coeff)

    ja = sp.Matrix(A).jacobian((x, y)).det()
    jb = sp.Matrix(B).jacobian((u, v)).det().subs(sub, simultaneous=True)
    chain = sp.cancel(ja*jb-k)
    coeffs, _ = numerator_coefficients(chain)
    eqs.extend(coeffs)
    return A, B, C, list(dict.fromkeys(map(sp.factor, eqs)))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--show-equations", action="store_true")
    args = parser.parse_args()
    A, B, C, eqs = equations()
    print("A =", A)
    print("B =", B)
    print("B o A =", C)
    print(f"generated {len(eqs)} coefficient equations")
    if args.show_equations:
        for e in eqs:
            print("  ", e, "= 0")

    # Normalize the nonzero constant Jacobian to one. Solve across both maps;
    # otherwise generic A-parameters hide conditional nonzero families.
    unknowns = (a0, a1, a2, a3, b0, b1, b2, b3, k)
    raw = sp.solve(eqs + [k-1], unknowns, dict=True, simplify=False)
    solutions = []
    for solution in raw:
        if solution not in solutions:
            solutions.append(solution)
    print(f"found {len(solutions)} symbolic solution families with k=1")
    for i, solution in enumerate(solutions, 1):
        print(f"[{i}]", solution)
    print("This tiny ansatz yields only triangular/birational families after")
    print("specialization. Enlarge supports, recheck exactly, and then test")
    print("branch identification and generic function-field degree separately.")

if __name__ == "__main__":
    main()
