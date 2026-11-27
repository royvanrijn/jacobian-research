#!/usr/bin/env python3
"""Certify the balanced tetrahedral standard block; see TETRAHEDRAL_STANDARD_GAP.md."""

from functools import cache
from math import comb, factorial
import json
import platform

import flint
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat


DEGREE = 32
COMPARISON_RHO = fmpq(19, 20)
STANDARD_GAP = fmpq(1, 1000)
# (left endpoint, right endpoint, strict upper bound, method)
SLABS = (
    (fmpq(0), fmpq(1, 2), fmpq(19, 20), "Hilbert–Schmidt"),
    (fmpq(1, 2), fmpq(3, 5), fmpq(113, 100), "Temple"),
    (fmpq(3, 5), fmpq(16, 25), fmpq(31, 25), "Temple"),
    (fmpq(16, 25), fmpq(33, 50), fmpq(13, 10), "Temple"),
    (fmpq(33, 50), fmpq(67, 100), fmpq(133, 100), "Temple"),
    (fmpq(67, 100), fmpq(169, 250), fmpq(1345, 1000), "Temple"),
    (fmpq(169, 250), fmpq(17, 25), fmpq(1356, 1000), "Temple"),
)


@cache
def hermite(n, variance):
    """Integer coefficients, low degree first, for variance-scaled Hermite."""
    if n == 0:
        return (1,)
    previous = hermite(n - 1, variance)
    if n == 1:
        return (0, 1)
    coefficients = [0, *previous]
    for i, c in enumerate(hermite(n - 2, variance)):
        coefficients[i] -= (n - 1) * variance * c
    return tuple(coefficients)


def accumulate(target, u, t, coefficient):
    if coefficient:
        target[u, t] = target.get((u, t), 0) + coefficient


@cache
def shifted_hermites(j, normal):
    """Even-in-t part is taken after multiplying by H_i(t)."""
    polynomial = {}
    for d, c in enumerate(hermite(j, 2)):
        for e, a in enumerate(hermite(normal, 2)):
            for v in range(d + 1):
                for w in range(e + 1):
                    accumulate(polynomial, d + e - v - w, v + w,
                               c * a * comb(d, v) * comb(e, w) * (-1)**w)
    return polynomial


def coefficient_matrices():
    # Nonconstant even-in-t monomials, in ascending total degree.
    basis = [(d - j, j) for d in range(1, DEGREE + 1)
             for j in range(0, d + 1, 2)]
    index = {key: k for k, key in enumerate(basis)}
    by_degree = [fmpq_mat(len(basis), len(basis)) for _ in range(DEGREE + 1)]
    for m in range(1, DEGREE + 1):
        matrix = by_degree[m]
        for normal in range(0, m + 1, 2):
            for j in range(0, m - normal + 1, 2):
                i = m - normal - j
                at_zero = hermite(normal, 2)[0]
                factor = fmpq(at_zero, 2**(1 + j + normal) *
                              factorial(i) * factorial(j) * factorial(normal))
                left = {(u, t): c * a * 2**t
                        for u, c in enumerate(hermite(i, 1))
                        for t, a in enumerate(hermite(j, 2))
                        if c and a and (u, t) != (0, 0)}
                right = {key: at_zero * coefficient for key, coefficient in left.items()}
                for (u, t), c in shifted_hermites(j, normal).items():
                    for d, a in enumerate(hermite(i, 1)):
                        if (t + d) % 2 == 0 and (u, t + d) != (0, 0):
                            accumulate(right, u, t + d, c * a)
                for x, a in left.items():
                    for y, c in right.items():
                        matrix[index[x], index[y]] += factor * a * c
        if matrix != matrix.transpose():
            raise ArithmeticError(f"Hermite coefficient {m} is asymmetric")
    return basis, by_degree


def moments(basis):
    root2 = arb(2).sqrt()
    alpha = root2.atan()
    cosine, sine = 1 / arb(3).sqrt(), (arb(2) / 3).sqrt()

    @cache
    def angular(i, j):
        if i >= 2:
            return ((i - 1) * angular(i - 2, j) + cosine**(i - 1) * sine**(j + 1)) / (i + j)
        if i == 1:
            return sine**(j + 1) / (j + 1)
        if j == 0:
            return alpha
        return ((j - 1) * angular(i, j - 2) - cosine * sine**(j - 1)) / j

    @cache
    def moment(i, j):
        return root2**i * (1 + arb(i + j) / 2).gamma() * angular(i, j) / alpha

    gram = arb_mat([[moment(i + u, j + t) - moment(i, j) * moment(u, t)
                     for u, t in basis] for i, j in basis])
    return alpha, gram


def require(proposition, message):
    if not proposition:
        raise ArithmeticError(message)


def certify():
    ctx.prec = 256
    basis, by_degree = coefficient_matrices()
    alpha, gram = moments(basis)

    # Exact low-degree controls from the first two centered kernel coefficients.
    u, u2, t2 = (basis.index(key) for key in [(1, 0), (2, 0), (0, 2)])
    require(by_degree[1][u, u] == fmpq(1, 2), "Incorrect degree-one control")
    require((by_degree[2][u2, u2], by_degree[2][u2, t2], by_degree[2][t2, t2])
            == (fmpq(1, 4), fmpq(1, 4), fmpq(5, 4)), "Incorrect degree-two control")

    rows = []
    b = arb(COMPARISON_RHO)
    require(0 < b and b < 1, "Invalid comparison correlation")
    sigma = (1 - b*b).sqrt()
    tail_prefactor = arb.pi() / (2 * alpha * b * sigma)
    hs_prefactor = arb.pi() / (alpha * b * sigma**3)
    last = fmpq(0)
    for left, right, target, method in SLABS:
        require(last == left and left < right < COMPARISON_RHO, "Invalid slab cover")
        last = right
        rho = arb(right)
        matrix = fmpq_mat(len(basis), len(basis))
        for coefficient in reversed(by_degree[1:]):
            matrix = right * matrix + coefficient
        action = arb_mat(matrix) * gram
        hs2 = sum((action[i, j] * action[j, i]
                   for i in range(len(basis)) for j in range(len(basis))), arb(0))
        tail = tail_prefactor * (rho / b)**DEGREE
        hs = hs2.sqrt()
        hs_tail = hs_prefactor * (rho / b)**DEGREE
        row = {"left": str(left), "right": str(right), "method": method,
               "truncated_hilbert_schmidt": hs.str(12), "operator_tail_bound": tail.str(12),
               "full_hilbert_schmidt_lower_bound": hs.str(12),
               "full_hilbert_schmidt_upper_bound": (hs + hs_tail).str(12)}
        upper = hs
        if method == "Temple":
            trial = arb_mat(len(basis), 1)
            trial[u, 0] = 1
            for _ in range(14):
                trial = action * trial
            dual = trial.transpose() * gram
            norm2 = (dual * trial)[0, 0]
            require(norm2 > 0, f"Vanishing trial vector at rho={right}")
            applied = action * trial
            rayleigh = (dual * applied)[0, 0] / norm2
            residual2 = ((applied.transpose() * gram * applied)[0, 0] /
                         norm2 - rayleigh**2)
            complement = (hs2 - rayleigh**2).sqrt()
            require(rayleigh > complement, f"No spectral isolation at rho={right}")
            require(residual2 >= 0, f"Residual enclosure too wide at rho={right}")
            upper = rayleigh + residual2 / (rayleigh - complement)
            row.update(truncated_rayleigh=rayleigh.str(12), residual_squared=residual2.str(12),
                       full_rayleigh_lower_bound=rayleigh.str(12),
                       full_rayleigh_upper_bound=(rayleigh + tail).str(12),
                       complementary_spectral_bound=complement.str(12))

        limit = 1 / (1 - arb(left)**2).sqrt()
        require(upper + tail < arb(target), f"Upper bound failed at rho={right}: {upper + tail}")
        require(arb(target + STANDARD_GAP) < limit, f"Coercivity threshold failed at rho={right}")
        row.update(operator_upper_bound=(upper + tail).str(12), certified_upper=str(target))
        rows.append(row)

    return {"command": ".venv/bin/python research/plurality-stability/check_tetrahedral_standard_gap.py",
            "python": platform.python_version(), "python_flint": flint.__version__,
            "precision_bits": ctx.prec, "degree": DEGREE,
            "comparison_correlation": str(COMPARISON_RHO),
            "standard_gap_in_units_of_ell": str(STANDARD_GAP), "slabs": rows}


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
