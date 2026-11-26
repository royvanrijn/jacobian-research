"""Exact Sylvester-interpolation checker for two increasing unit-gap rows."""

from argparse import ArgumentParser
from collections import defaultdict
from hashlib import sha256
from math import comb, factorial


def coefficients(n, S, P):
    terms = [1]
    rising_product = 1
    for k in range(1, n+1):
        rising_product *= k*k + n*k*S + n*n*P
        terms.append(comb(n,k)*rising_product)
    return terms


def bareiss(matrix):
    n = len(matrix)
    sign = 1
    previous = 1
    for i in range(n-1):
        if matrix[i][i] == 0:
            swap = next((j for j in range(i+1,n) if matrix[j][i]), None)
            if swap is None:
                return 0
            matrix[i], matrix[swap] = matrix[swap], matrix[i]
            sign = -sign
        pivot = matrix[i][i]
        for j in range(i+1,n):
            for k in range(i+1,n):
                numerator = pivot*matrix[j][k] - matrix[j][i]*matrix[i][k]
                matrix[j][k], remainder = divmod(numerator, previous)
                assert remainder == 0
            matrix[j][i] = 0
        previous = pivot
    return sign*matrix[-1][-1]


def resultant_value(n,S,P):
    first = coefficients(n,S,P)
    second = coefficients(n+1,S,P)
    matrix = []
    for shift in range(n+1):
        matrix.append([0]*shift + first + [0]*(n-shift))
    for shift in range(n):
        matrix.append([0]*shift + second + [0]*(n-shift-1))
    return bareiss(matrix)


def differences(values):
    values = list(values)
    while values:
        yield values[0]
        values = [right-left for left,right in zip(values,values[1:])]


def resultant_polynomial(n):
    degree = n*(n+1)
    # Interpolate first in S, then in P on the triangular grid.
    rows = [
        list(differences(resultant_value(n, S, P) for S in range(degree - P + 1)))
        for P in range(degree + 1)
    ]
    factorial_d = factorial(degree)
    # Represent each binomial polynomial with the common denominator degree!.
    fallings = [[1]]
    scaled = [[factorial_d]]
    for k in range(1,degree+1):
        old = fallings[-1]
        current = [0]*(k+1)
        for i,c in enumerate(old):
            current[i] -= (k-1)*c
            current[i+1] += c
        fallings.append(current)
        scaled.append([c*(factorial_d//factorial(k)) for c in current])
    out = defaultdict(int)
    for i in range(degree+1):
        for j,c in enumerate(differences(rows[P][i] for P in range(degree-i+1))):
            if c:
                for x, u in enumerate(scaled[i]):
                    for y, v in enumerate(scaled[j]):
                        out[x,y] += c*u*v
    result = {}
    for pos,c in out.items():
        quo, rem = divmod(c, factorial_d**2)
        assert rem == 0
        if quo:
            result[pos] = quo
    assert result[0,0] > 0
    return result


def divide_unit_factor(poly):
    rest = defaultdict(int, poly)
    out = {}
    degree = max(sum(k) for k in poly)
    for i in range(degree,0,-1):
        for j in range(degree-i+1):
            c = rest.pop((i,j),0)
            if c:
                out[i-1,j] = c
                rest[i-1,j] -= c
                rest[i-1,j+1] -= c
    assert not any(rest.values())
    return out


def ordinary_coordinates(poly):
    out = defaultdict(int)
    for (i,j),c in poly.items():
        for l in range(i+1):
            out[l+j,i-l+j] += c*comb(i,l)
    return dict(out)


def shifted_coordinates(poly):
    shifted_S = defaultdict(int)
    for (i,j),c in poly.items():
        for l in range(i+1):
            shifted_S[l,j] += c*comb(i,l)*2**(i-l)
    out = defaultdict(int)
    for (i,j),c in shifted_S.items():
        for l in range(j+1):
            for m in range(j-l+1):
                out[i+m,l] += c*comb(j,l)*comb(j-l,m)
    return dict(out)


def value(poly, S, P):
    return sum(c*S**i*P**j for (i,j),c in poly.items())


def verify_shift_witnesses():
    witnesses = (
        (1, 16, 27, "e494da3b000cfbeaa14c22c67e21673a3d0cf2ae9961e0fb37f326d4cc487a0e"),
        (2, 21, 220, "84831d03f827e5a69a3d7a869d06cbcc1ff5ff59e3c10fbc6fd97ef8e9e5831d"),
    )
    for shift, n, shifted_product, expected_digest in witnesses:
        shifted_sum = 0
        original_sum = shifted_sum + 2*shift
        original_product = shifted_product + shift*shift
        resultant = resultant_value(n, original_sum, original_product)
        quotient, remainder = divmod(resultant, 1 + original_sum + original_product)
        assert remainder == 0
        assert resultant < 0 and quotient < 0
        digest = sha256(str(resultant).encode("ascii")).hexdigest()
        assert digest == expected_digest
        print("negative shifted resultant", shift, n, shifted_sum, shifted_product, digest)


def main():
    if not __debug__:
        raise RuntimeError("Run without Python optimization to enable the exact checks.")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--shift-witnesses-only", action="store_true",
        help="check only the two negative shifted point evaluations",
    )
    args = parser.parse_args()
    verify_shift_witnesses()
    if args.shift_witnesses_only:
        return

    resultant_witness = -int(
        "950844528495324832888575815277770346318989613491315656077122140236"
        "3311847372186386362265625000"
    )
    quotient_witness = -int(
        "587555943571407548820183814836214577759238428301808919279143064688"
        "86463756364419595179404296875"
    )
    for n in range(1,8):
        resultant = resultant_polynomial(n)
        quotient = divide_unit_factor(resultant)
        degree = n*(n+1)
        assert len(resultant) == (degree+1)*(degree+2)//2
        assert len(quotient) == degree*(degree+1)//2
        assert resultant[0,degree] == (2*n+1)**degree
        assert value(resultant, -1, 2) == resultant_value(n, -1, 2)
        assert value(resultant, 2, -1) == resultant_value(n, 2, -1)
        neg_resultant = sorted(k for k,c in resultant.items() if c<0)
        neg_quotient = sorted(k for k,c in quotient.items() if c<0)
        expected = [(0,j) for j in range(36,40)] if n == 7 else []
        assert neg_resultant == expected
        assert neg_quotient == [(i,j-1) for i,j in expected]
        for polynomial, size in ((resultant, degree), (quotient, degree-1)):
            ordinary = ordinary_coordinates(polynomial)
            shifted = shifted_coordinates(polynomial)
            assert len(ordinary) == (size+1)**2
            assert len(shifted) == (size+1)*(size+2)//2
            assert min(ordinary.values()) > 0
            assert min(shifted.values()) > 0
        print(n, len(resultant), len(quotient), neg_resultant, neg_quotient)
        if n == 7:
            assert resultant[0,39] == resultant_witness
            assert quotient[0,38] == quotient_witness


if __name__ == "__main__":
    main()
