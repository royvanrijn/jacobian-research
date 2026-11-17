"""Rebuild exact local jets and independently verify the finite certificate.

This checks rational arithmetic and finite linear algebra.  The accompanying
mathematical proof must separately justify the formal slice, the filtered tail
bound, and passage from these truncated rings to completed moment ideals.
"""
from pathlib import Path
import json
import time
import sympy as sp

started = time.monotonic()
import rod_local_full_prefix28 as calc

a, b, z, t, w = calc.a, calc.b, calc.z, calc.t, calc.w
destination = Path(__file__).resolve().parent

# Independent coordinates for the twelve slice directions and four orbit
# directions, in the sixteen-element harmonic basis used by the jet scripts.
columns = []
for index in (0, 1, 2, 3, 6, 7, 8, 13, 14, 15):
    column = sp.zeros(16, 1)
    column[index] = 1
    columns.append(column)
for coefficients in (
    {8: 1, 12: -sp.Rational(1, 10), 2: sp.Rational(17, 70)},
    {7: 1, 11: sp.Rational(2, 5), 1: -sp.Rational(1, 5)},
    {3: 1, 5: 1, 9: -sp.Rational(1, 2)},
    {3: -1, 5: 1, 9: -sp.Rational(3, 2)},
    {4: 1, 2: -1},
    {6: -1, 10: sp.Rational(3, 2)},
):
    column = sp.zeros(16, 1)
    for index, coefficient in coefficients.items():
        column[index] = coefficient
    columns.append(column)
transverse_matrix = sp.Matrix.hstack(*columns)
assert transverse_matrix.det() == -sp.Rational(3, 50)
assert transverse_matrix.rank() == 16
assert calc.J0.det() == sp.Rational(33554432, 789441476538345505875)
orbit_polynomials = [sum(column[i]*calc.basis[i] for i in range(16))
                     for column in columns[-4:]]
assert sp.expand(orbit_polynomials[0]-calc.F0) == 0
assert sp.expand(orbit_polynomials[1]-(z*sp.diff(calc.F0,z)-w*sp.diff(calc.F0,w))) == 0
assert sp.expand(orbit_polynomials[2]-(z*sp.diff(calc.F0,t)-t*sp.diff(calc.F0,w))) == 0
assert sp.expand(orbit_polynomials[3]-(w*sp.diff(calc.F0,t)-t*sp.diff(calc.F0,z))) == 0
print("Transverse determinant -3/50; rank 16; first-ten Jacobian certified.", flush=True)

residuals = calc.full_residuals
weight = lambda exponents: (3*exponents[0]+2*exponents[1], exponents[0])

def initial(poly):
    terms = poly.as_dict()
    least = min(weight(exponents)[0] for exponents in terms)
    return sp.expand(sum(c*a**i*b**j for (i,j),c in terms.items()
                         if weight((i,j))[0] == least))

initials = {n: initial(residuals[n]) for n in (11,12,13,14)}
assert all(min(weight(e)[0] for e in residuals[n].as_dict()) == n
           for n in (11,12,13,14))
A = residuals[11].coeff_monomial(a**3*b)
C = residuals[12].coeff_monomial(a**4)
f = sp.expand(initials[11]/A)
g = sp.expand(initials[12]/C)
alpha = sp.Rational(6488776,18982155)
beta = sp.Rational(3109877708,1048504325)
gamma = sp.Rational(4973565856,29245471875)
assert f == a**3*b+alpha*a*b**4
assert g == a**4+beta*a**2*b**3+gamma*b**6

h13 = sp.expand(initials[13]-sp.Poly(initials[13],a,b).coeff_monomial(a**3*b**2)*b*f)
c13 = sp.Rational(2538438656,2305418019694453125)
assert h13 == c13*a*b**5
first_h14 = sp.expand(initials[14]-sp.Poly(initials[14],a,b).coeff_monomial(a**4*b)*b*g)
relation = sp.expand(b*g-a*f)
h14 = sp.expand(first_h14-sp.Poly(first_h14,a,b).coeff_monomial(a**2*b**4)/(beta-alpha)*relation)
c14 = -sp.Rational(809733013864448,95603262363699979798828125)
assert h14 == c14*b**7

expected_leads = {
    12: {(4,0),(3,1),(2,4),(1,7),(0,10)},
    13: {(4,0),(3,1),(2,4),(1,5),(0,8)},
    14: {(4,0),(3,1),(2,4),(1,5),(0,7)},
}
weighted_footprints = {}
for prefix, generators in ((12,[f,g]),(13,[f,g,h13]),(14,[f,g,h13,h14])):
    groebner = sp.groebner(generators,a,b,order=weight,domain=sp.QQ)
    leads = {poly.LM(order=weight).exponents for poly in groebner.polys}
    assert leads == expected_leads[prefix]
    footprint = [(i,j) for i in range(4) for j in range(10)
                 if not any(i >= u and j >= v for u,v in leads)]
    assert len(footprint) == {12:22,13:18,14:17}[prefix]
    weighted_footprints[prefix] = footprint
    print(f"Weighted initial ideal {prefix}: leads={sorted(leads)}; length={len(footprint)}", flush=True)
print(f"Reduced R13 coefficient: {c13}; reduced R14 coefficient: {c14}", flush=True)

# Build a second set of exact rational matrices directly from coefficient
# dictionaries, independently of the incremental echelon implementation.
monomials = [(i,d-i) for d in range(8) for i in range(d+1)]
assert len(monomials) == 36

def make_rows(indices):
    rows, labels = [], []
    for n in indices:
        terms = residuals[n].as_dict()
        if residuals[n].is_zero:
            continue
        order = min(i+j for (i,j),c in terms.items() if c)
        for degree in range(8-order):
            for u in range(degree+1):
                v = degree-u
                rows.append([terms.get((i-u,j-v),sp.S.Zero) for i,j in monomials])
                labels.append([n,u,v])
    return rows, labels

rows14, labels14 = make_rows(range(11,15))
tail_rows, tail_labels = make_rows(range(15,29))
assert len(rows14) == 32 and len(tail_rows) == 18
M14 = sp.Matrix(rows14)
M28 = sp.Matrix(rows14+tail_rows)
assert M14.shape == (32,36) and M28.shape == (50,36)
reduced, pivots = M14.rref()
assert len(pivots) == M14.rank() == M28.rank() == 19
assert all(residuals[n].is_zero for n in range(22,29))
for prefix, expected in ((12,20),(13,18),(14,17)):
    rows, _ = make_rows(range(11,prefix+1))
    assert 36-sp.Matrix(rows).rank() == expected
print("Exact ranks: 32x36 and augmented 50x36 both rank 19.", flush=True)

# Multiplication in the finite ring Q[a,b]/(R11,R12,R13,R14,(a,b)^8).
free = [index for index in range(36) if index not in pivots]
assert len(free) == 17
quotient_monomials = [monomials[index] for index in free]
reductions = {}
for j,index in enumerate(free):
    vector = sp.zeros(17,1)
    vector[j] = 1
    reductions[monomials[index]] = vector
for row,index in enumerate(pivots):
    reductions[monomials[index]] = sp.Matrix([-reduced[row,j] for j in free])

degree7_monomials = [exponent for exponent in monomials if sum(exponent) == 7]
weight13_monomials = [exponent for exponent in monomials
                      if weight(exponent)[0] >= 13]
assert all(reductions[exponent] == sp.zeros(17,1)
           for exponent in degree7_monomials)
assert all(reductions[exponent] == sp.zeros(17,1)
           for exponent in weight13_monomials)
assert reductions[(0,6)] != sp.zeros(17,1)
print(f"Direct row membership: all {len(degree7_monomials)} degree-seven monomials "
      f"and all {len(weight13_monomials)} truncated monomials of weight >=13; "
      "b^6 is not in the row span.",flush=True)

def multiply_by(step):
    output = []
    for i,j in quotient_monomials:
        exponent = (i+step[0],j+step[1])
        output.append(reductions.get(exponent,sp.zeros(17,1)))
    return sp.Matrix.hstack(*output)

Ma, Mb = multiply_by((1,0)), multiply_by((0,1))
assert Ma*Mb == Mb*Ma
powers_a = [sp.eye(17)]
powers_b = [sp.eye(17)]
for degree in range(1,8):
    powers_a.append(powers_a[-1]*Ma)
    powers_b.append(powers_b[-1]*Mb)
ideal_dimensions = []
for degree in range(8):
    images = [powers_a[i]*powers_b[degree-i] for i in range(degree+1)]
    ideal_dimensions.append(sp.Matrix.hstack(*images).rank())
assert ideal_dimensions == [17,16,14,11,7,4,1,0]
assert powers_b[6] != sp.zeros(17) and all(
    powers_a[i]*powers_b[7-i] == sp.zeros(17) for i in range(8))
hilbert = [ideal_dimensions[i]-ideal_dimensions[i+1] for i in range(7)]
assert hilbert == [1,2,3,4,3,3,1]
print(f"Finite quotient Hilbert function: {hilbert}; maximal-ideal nilpotence index 7.", flush=True)

def encode_scalar(value):
    value = sp.sympify(value)
    assert value.is_Rational and not value.has(sp.Float)
    return str(value)

def encode_matrix(matrix):
    return [[encode_scalar(v) for v in row] for row in matrix.tolist()]

data = {
    "field":"Q",
    "truncation":"(a,b)^8",
    "monomials":monomials,
    "residuals":{str(n):[[list(exponents),encode_scalar(coefficient)]
                            for exponents,coefficient in sorted(poly.as_dict().items()) if coefficient]
                 for n,poly in sorted(residuals.items())},
    "first14_matrix":{"row_labels":labels14,"entries":encode_matrix(M14),"rank":19},
    "first28_matrix":{"row_labels":labels14+tail_labels,"entries":encode_matrix(M28),"rank":19},
    "weighted_leading_monomials":{str(n):sorted(leads) for n,leads in expected_leads.items()},
    "weighted_footprints":weighted_footprints,
    "checked_degree7_membership":degree7_monomials,
    "checked_weight13_membership":weight13_monomials,
    "reduced_initial_coefficients":{"R13_ab5":str(c13),"R14_b7":str(c14)},
    "quotient_basis":quotient_monomials,
    "multiplication_a":encode_matrix(Ma),
    "multiplication_b":encode_matrix(Mb),
    "maximal_ideal_power_dimensions":ideal_dimensions,
    "hilbert_function":hilbert,
    "nilpotence_index":7,
}
output = destination/"certificate_data.json"
output.write_text(json.dumps(data,separators=(",",":"))+"\n")
print(f"Exported exact certificate: {output.name} ({output.stat().st_size} bytes)", flush=True)
print("ALL ASSERTIONS PASSED; full-prefix finite quotient length 17.", flush=True)
print("Elapsed seconds:",time.monotonic()-started,flush=True)
