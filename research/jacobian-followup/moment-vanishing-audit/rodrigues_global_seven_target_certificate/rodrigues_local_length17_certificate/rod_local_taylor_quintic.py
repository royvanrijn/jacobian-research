"""Bounded extension of the exact transverse IFT calculation to orders 5/6."""
import time
started = time.monotonic()
from rod_local_taylor_quartic import *

def prepare(expressions):
    return {k: sp.Poly(sp.expand(f), a, b).as_dict()
            for k, f in expressions.items()}

def jet(n, prepared):
    return sp.expand(sum(
        sp.binomial(n, k)*sum(H(n-k, coeff)*a**i*b**j
                             for (i, j), coeff in terms.items())
        for k, terms in prepared.items() if n >= k))

def solve_correction(prepared):
    v = sp.Matrix([jet(n, prepared) for n in range(1, 11)])
    correction = sp.expand(sum(c*q for c, q in zip(-Jinv*v, directions)))
    assert all(sp.expand(n*parameter_H(n-1, correction)+v[n-1]) == 0
               for n in range(1, 11))
    return correction

parts4 = prepare({2: 2*T1*T3+T2*T2,
                  3: 3*T1*T1*T2,
                  4: T1**4})
T4 = solve_correction(parts4)
parts5 = prepare({2: 2*T1*T4+2*T2*T3,
                  3: 3*T1*T1*T3+3*T1*T2*T2,
                  4: 4*T1**3*T2,
                  5: T1**5})
for n in (11, 12, 13, 14):
    print(f"R{n} quintic:", sp.factor(jet(n, parts5)), flush=True)
print("Elapsed through quintics:", time.monotonic()-started, flush=True)

T5 = solve_correction(parts5)
parts6 = prepare({2: 2*T1*T5+2*T2*T4+T3*T3,
                  3: 3*T1*T1*T4+6*T1*T2*T3+T2**3,
                  4: 4*T1**3*T3+6*T1*T1*T2*T2,
                  5: 5*T1**4*T2,
                  6: T1**6})
for n in (11, 12, 13, 14):
    print(f"R{n} sextic:", sp.factor(jet(n, parts6)), flush=True)
monomials = [(i, d-i) for d in range(7) for i in range(d+1)]
residuals = {n: sp.Poly(jet(n, parts4)+jet(n, parts5)+jet(n, parts6), a, b)
             for n in (11, 12, 13, 14)}
def rows_for(indices):
    rows = []
    for n in indices:
        for degree in range(3):
            for i in range(degree+1):
                product = residuals[n] * sp.Poly(a**i*b**(degree-i), a, b)
                rows.append([product.coeff_monomial(a**j*b**k)
                             for j, k in monomials])
    return rows
rank12 = sp.Matrix(rows_for((11, 12))).rank()
rank14 = sp.Matrix(rows_for((11, 12, 13, 14))).rank()
print("Dimensions modulo (a,b)^7:",
      {"A12": len(monomials)-rank12, "A14": len(monomials)-rank14}, flush=True)
alpha = residuals[13].coeff_monomial(a**3*b**2) / residuals[11].coeff_monomial(a**3*b)
rem13 = residuals[13]-sp.Poly(alpha*b, a, b)*residuals[11]
print("R13 after subtracting its quintic b*R11 multiple, degree6 terms:",
      [(ij, sp.factor(c)) for ij, c in rem13.as_dict().items() if sum(ij) == 6], flush=True)
print("Elapsed through sextics:", time.monotonic()-started, flush=True)
