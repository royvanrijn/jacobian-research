"""Direct exact contraction checks for Gaussian/sphere moment interpolation.

The all-degree and surjectivity claims use the accompanying written proof.
This checker independently reconstructs finite Jacobians, without inverse shear.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import factorial, prod
from pathlib import Path
import json
import platform
import sys

if not __debug__:
    raise RuntimeError("Exact certificates require assertions enabled")

@lru_cache(maxsize=None)
def odd_factorial(n):
    return prod(range(1, n + 1, 2))

def monomial_mean(a, b, c, functional):
    if a != b or c % 2:
        return Fraction(0)
    value = Fraction(factorial(a) * odd_factorial(c - 1))
    if functional == "sphere":
        value /= odd_factorial(a + b + c + 1)
    return value

def direct_mixed(moment_power, marker, functional, shear_order=1):
    """Literal multinomial contraction, keeping exactly the balanced terms."""
    a, b, c = marker
    r = shear_order
    numerator = moment_power+b-a
    if numerator < 0 or numerator % (r+1):
        return Fraction(0)
    h = numerator//(r+1)
    if (c+h) % 2:
        return Fraction(0)
    value = Fraction(0)
    for k in range(max(0, h-moment_power), min(moment_power, h//2)+1):
        j = h-2*k
        i = moment_power-h+k
        coefficient = Fraction(factorial(moment_power),
                               factorial(i)*factorial(j)*factorial(k))
        coefficient *= Fraction(-1, 2)**k
        value += coefficient * monomial_mean(a+r*j+(2*r+1)*k, b+i, c+j,
                                              functional)
    return value

def direct_mixed_full(moment_power, marker, functional, shear_order=1):
    """Separate unfiltered expansion checks the balance-sum implementation."""
    a, b, c = marker
    value = Fraction(0)
    # P_r^s=(W+Z^rY-Z^(2r+1)/2)^s, by the multinomial theorem.
    for j in range(moment_power + 1):
        for k in range(moment_power - j + 1):
            i = moment_power - j - k
            coefficient = Fraction(factorial(moment_power),
                                   factorial(i)*factorial(j)*factorial(k))
            coefficient *= Fraction(-1, 2)**k
            value += coefficient * monomial_mean(a+shear_order*j+(2*shear_order+1)*k, b+i, c+j,
                                                  functional)
    return value

def markers(degree, shear_order=1):
    r = shear_order
    size = (2*r+1)*(degree-r+1)+r
    available = {}
    # Prefer maximal W and then maximal Y exponent; for r=1 these
    # are exactly the original W^a, ZW^a, YW^a markers.
    for k in range(degree, -1, -1):
        for j in range(degree-k, -1, -1):
            for i in range(degree-j-k+1):
                weight = i+(r+1)*j+(2*r+1)*k
                available.setdefault(weight, (i, k, j))
    assert all(n in available for n in range(size))
    assert size not in available
    return [available[n] for n in range(size)]

def projected_markers(degree, shear_order):
    r = shear_order
    assert r >= 3 and degree >= 2*r+1
    basis = markers(degree, r)
    # Tuples retain ORIGINAL spatial exponents (Z,W,Y).
    basis.append((0, degree-r+2, 1))
    basis.extend((j-1, degree-r+1, 1) for j in range(1, r))
    basis.append((0, degree-r+3, 0))
    basis.extend((j-1, degree-r+2, 0) for j in range(1, r))
    assert len(basis) == len(set(basis))
    assert all(sum(marker) <= degree for marker in basis)
    return basis

def verify_matrix(matrix, degree, functional, shear_order=1, projected=False):
    r = shear_order
    basis = projected_markers(degree, r) if projected else markers(degree, r)
    base_size = (2*r+1)*(degree-r+1)+r
    size = len(basis)
    assert len(matrix) == size and all(len(row) == size for row in matrix)
    expected = []
    for n, (_, k, _) in enumerate(basis):
        if projected and n == base_size:
            coefficient = r*k*(2*k-1)*Fraction(-1, 2)**(k-1)
        elif projected and n == base_size+r:
            coefficient = 2*r*k*(k-1)*Fraction(-1, 2)**(k-1)
        else:
            coefficient = Fraction(-1, 2)**k
        value = factorial(n+1) * coefficient
        if functional == "sphere":
            special = projected and n in (base_size, base_size+r)
            value /= odd_factorial(2*n+3 if special else 2*n+1)
        expected.append(value)
    assert all(matrix[i][j] == 0 for i in range(size) for j in range(i))
    assert [matrix[i][i] for i in range(size)] == expected
    determinant = prod(matrix[i][i] for i in range(size))
    exponent = sum(k for _, k, _ in basis)
    if shear_order == 1:
        assert exponent == degree*(3*degree-1)//2
    formula = prod(factorial(j) for j in range(1, size+1))
    formula *= Fraction(-1, 2)**exponent
    if projected:
        b1, b2 = degree-r+2, degree-r+3
        formula *= 4 * (r*b1*(2*b1-1)) * (2*r*b2*(b2-1))
    if functional == "sphere":
        formula /= prod(odd_factorial(2*n+1) for n in range(size))
        if projected:
            formula /= (2*base_size+3)*(2*(base_size+r)+3)
    assert determinant == formula != 0
    return determinant

def run():
    records = []
    direct_comparisons = 0
    for functional in ("gaussian", "sphere"):
        assert direct_mixed(0, (0, 0, 0), functional) == 1
        assert direct_mixed(1, (0, 0, 0), functional) == 0
        assert monomial_mean(1, 1, 0, functional) == (
            1 if functional == "gaussian" else Fraction(1, 3))
        for r in range(1, 5):
            for s in range(9):
                for marker in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),
                               (0,2,0),(0,1,1),(1,1,1),(0,3,1)]:
                    assert direct_mixed(s, marker, functional, r) == direct_mixed_full(
                        s, marker, functional, r)
                    direct_comparisons += 1
        cases = [(1,d,False) for d in range(3,7)]
        cases += [(2,5,False),(2,6,False),(3,7,False)]
        cases += [(3,7,True),(3,8,True),(4,9,True)]
        for shear_order, degree, projected in cases:
            basis = (projected_markers(degree, shear_order) if projected
                     else markers(degree, shear_order))
            base_size = (2*shear_order+1)*(degree-shear_order+1)+shear_order
            assert len(set(basis)) == base_size+(2*shear_order if projected else 0)
            assert all(sum(marker) <= degree for marker in basis)
            matrix = [[j*direct_mixed(j-1, marker, functional, shear_order)
                       for marker in basis] for j in range(1, len(basis)+1)]
            determinant = verify_matrix(matrix, degree, functional, shear_order, projected)
            corrupted = [row[:] for row in matrix]
            corrupted[-1][-1] += 1
            try:
                verify_matrix(corrupted, degree, functional, shear_order, projected)
            except AssertionError:
                pass
            else:
                raise AssertionError("Corrupted Jacobian was accepted")
            # These finite controls support the implementation, not purity at all orders.
            assert all(direct_mixed(j, (0, 0, 0), functional, shear_order) == 0
                       for j in range(1, len(basis)+1))
            records.append({
                "functional": functional, "spatial_degree": degree,
                "shear_order": shear_order,
                "projected_jacobian": projected,
                "first_moment_count": len(basis),
                "original_markers_Z_W_Y": [list(marker) for marker in basis],
                "diagonal": [str(matrix[i][i]) for i in range(len(basis))],
                "determinant": str(determinant),
                "matrix_sha256": sha256(json.dumps(
                    [[str(value) for value in row] for row in matrix],
                    separators=(",", ":")).encode()).hexdigest(),
                "corrupted_matrix_rejected": True})
    return {
        "status": "passed", "arithmetic": "fractions.Fraction over Q",
        "python": platform.python_version(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "balanced_vs_full_multinomial_comparisons": direct_comparisons,
        "records": records,
        "scope": "Direct finite Jacobians at r=1,d=3..6 and (r,d)=(2,5),(2,6),(3,7) for each functional. "
                 "Projected finite Jacobians at (r,d)=(3,7),(3,8),(4,9), with 44,51,66 moments, "
                 "including exact subleading sphere normalizations. "
                 "General-degree nonvanishing, all-order purity and global "
                 "surjectivity are written proofs, not finite-check conclusions.",
        "independent_replay": False, "formal_verification": False,
        "external_review": False}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name(
                            "gaussian_moment_interpolation_verification.json"))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "matrices": len(result["records"]),
                      "output": str(args.output)}))

if __name__ == "__main__":
    main()
