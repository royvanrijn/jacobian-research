"""All residual moments through 28, ordinary jets modulo (a,b)^8.

Exact full local thickness additionally uses the separately proved bounds:
(a,b)^8 is contained in the first-13 ideal; moments above 28 vanish there.
The numerical rank calculation below makes neither inference on its own.
"""
import time
start_full = time.monotonic()
from rod_local_taylor_septic import *

parts2 = prepare({2: T1*T1})
parts3 = prepare({2: 2*T1*T2, 3: T1**3})
all_parts = (parts2, parts3, parts4, parts5, parts6, parts7)
echelon = {}

def add_row(row):
    row = list(row)
    for pivot in sorted(echelon):
        c = row[pivot]
        if c:
            row = [v-c*w for v, w in zip(row, echelon[pivot])]
    pivot = next((i for i, v in enumerate(row) if v), None)
    if pivot is not None:
        c = row[pivot]
        echelon[pivot] = [v/c for v in row]

def add_ideal_generator(poly):
    if poly.is_zero:
        return
    order = min(sum(ij) for ij, c in poly.as_dict().items() if c)
    for degree in range(8-order):
        for i in range(degree+1):
            product_poly = poly*sp.Poly(a**i*b**(degree-i), a, b)
            add_row([product_poly.coeff_monomial(a**j*b**k)
                     for j, k in monomials7])

for n in (11, 12, 13, 14):
    add_ideal_generator(residuals7[n])
assert len(monomials7)-len(echelon) == 17

full_residuals = dict(residuals7)
for n in range(15, 29):
    residual = sp.Poly(sum(jet(n, pp) for pp in all_parts), a, b)
    full_residuals[n] = residual
    old_dimension = len(monomials7)-len(echelon)
    add_ideal_generator(residual)
    dimension = len(monomials7)-len(echelon)
    order = min((sum(ij) for ij, c in residual.as_dict().items() if c), default=None)
    print(f"A{n} dimension modulo (a,b)^8: {dimension}; residual order: {order}", flush=True)
    if dimension != old_dimension:
        print(f"New residual R{n}:", sp.factor(residual.as_expr()), flush=True)

print("Final prefix-28 local jet dimension:", len(monomials7)-len(echelon), flush=True)
print("Elapsed including prior orders:", time.monotonic()-start_full, flush=True)
