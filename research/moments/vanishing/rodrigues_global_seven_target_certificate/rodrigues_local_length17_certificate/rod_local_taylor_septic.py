"""Exact seventh-order extension; no weighted truncation of the Taylor series."""
import time
started7 = time.monotonic()
import rod_local_taylor_quintic as previous
from rod_local_taylor_quintic import *
from math import factorial, comb, prod

@lru_cache(None)
def odd_df(n):
    return prod(range(n, 0, -2))

@lru_cache(None)
def fast_monomial(n, iz, it, iw):
    # Angular balance reduces the multinomial expansion of F0**n to one sum.
    twice_N = iw+n-iz
    if twice_N < 0 or twice_N % 2:
        return sp.S.Zero
    N = twice_N//2
    if (it+N) % 2:
        return sp.S.Zero
    numerator = 0
    for j in range(max(0, N-n), min(n, N//2)+1):
        v = N-2*j
        numerator += ((-1)**j * 2**(n-j) * comb(n, j) * comb(n-j, v)
                      * factorial(iz+N+j) * odd_df(it+N-2*j-1))
    return sp.Rational(numerator, 2**n * odd_df(2*iz+it+3*N+1))

@lru_cache(None)
def fast_H(n, marker):
    return sum(coeff*fast_monomial(n, *ijk)
               for ijk, coeff in poly_dict(marker).items())

for n in (0, 1, 3, 7, 11, 13):
    for marker in (Qa, Qb, Qa*Qb, Qa**2, z**2*w**3*t**4):
        assert fast_H(n, marker) == H(n, marker)

# Existing helpers close over their defining module; replace only their moment
# evaluator with the checked, equivalent one-dimensional multinomial sum.
previous.H = fast_H

print("Starting seventh-order implicit correction", flush=True)
T6 = solve_correction(parts6)
parts7 = prepare({
    2: 2*T1*T6+2*T2*T5+2*T3*T4,
    3: 3*T1*T1*T5+6*T1*T2*T4+3*T1*T3*T3+3*T2*T2*T3,
    4: 4*T1**3*T4+12*T1*T1*T2*T3+4*T1*T2**3,
    5: 5*T1**4*T3+10*T1**3*T2*T2,
    6: 6*T1**5*T2,
    7: T1**7,
})
for n in (11, 12, 13, 14):
    print(f"R{n} septic:", sp.factor(jet(n, parts7)), flush=True)

monomials7 = [(i, d-i) for d in range(8) for i in range(d+1)]
residuals7 = {n: sp.Poly(sum(jet(n, pp) for pp in (parts4, parts5, parts6, parts7)), a, b)
              for n in (11, 12, 13, 14)}

def rows7(indices):
    rows = []
    for n in indices:
        for degree in range(4):
            for i in range(degree+1):
                product7 = residuals7[n]*sp.Poly(a**i*b**(degree-i), a, b)
                rows.append([product7.coeff_monomial(a**j*b**k)
                             for j, k in monomials7])
    return rows

for label, indices in (("A12", (11, 12)), ("A13", (11, 12, 13)),
                       ("A14", (11, 12, 13, 14))):
    rank7 = sp.Matrix(rows7(indices)).rank()
    print(f"{label} dimension modulo (a,b)^8:", len(monomials7)-rank7, flush=True)
print("Elapsed including earlier orders:", time.monotonic()-started7, flush=True)
