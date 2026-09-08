"""Target-free, dimension-independent half-lattice scheduling utilities.

Version3 retains explicit PARI precision and uses a scale-aware floating consistency guard. No fixture loads.
Floating CVP and heights schedule searches; exact norm checks do not certify
CVP optimality or rank. All rank admissions need separate exact certificates.
"""
from typing import Sequence
from fractions import Fraction
from decimal import Decimal, getcontext
from fpylll import Enumeration, GSO, IntegerMatrix
from research_runtime.floating_norm_check import checked_distance_error
from sage.all import EllipticCurve, QQ, pari
Point = tuple[Fraction, Fraction]

class CosetOracle:
    """Floating fplll CVP decisions with exact rounded-form recomputation."""

    def __init__(self, gram: Sequence[Sequence[int]], degree: int = 2) -> None:
        self.gram = tuple(tuple(int(value) for value in row) for row in gram)
        self.dimension = len(self.gram)
        if not self.dimension or any(len(row) != self.dimension for row in self.gram):
            raise ValueError("the CVP Gram matrix must be nonempty and square")
        self.degree = int(degree)
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(self.gram),
            gram=True,
            float_type="dd",
            update=True,
        )
        self.mu = tuple(
            tuple(self.gso.get_mu(i, j) if i > j else 0.0 for j in range(self.dimension))
            for i in range(self.dimension)
        )
        self.distance_bound = (
            (degree - 1) ** 2
            * sum(abs(value) for row in self.gram for value in row)
            / (degree * degree)
            + 1.0
        )

    def solve(self, residue: Sequence[int]) -> tuple[int, tuple[int, ...], float]:
        residue = tuple(int(value) for value in residue)
        if len(residue) != self.dimension or any(value not in (0, 1) for value in residue):
            raise ValueError("a parity residue has the wrong shape")
        target = [
            -(
                residue[i]
                + sum(residue[j] * self.mu[j][i] for j in range(i + 1, self.dimension))
            )
            / self.degree
            for i in range(self.dimension)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0,
            self.dimension,
            self.distance_bound,
            0,
            target=target,
        )
        if not solutions:
            raise ArithmeticError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(value)) for value in coordinates)
        if any(abs(value - integer) > 1.0e-7 for value, integer in zip(coordinates, closest)):
            raise ArithmeticError("CVP enumeration returned nonintegral coordinates")
        representative = tuple(
            residue[index] + self.degree * closest[index]
            for index in range(self.dimension)
        )
        norm = sum(
            representative[i] * self.gram[i][j] * representative[j]
            for i in range(self.dimension)
            for j in range(self.dimension)
        )
        # Exact representative, parity and rounded norm remain authoritative.
        # This only guards numerical consistency; it does not certify CVP optimality.
        error = checked_distance_error(norm, reported_distance, self.degree)
        return norm, representative, error


def canonical_height_gram(model, basis: Sequence[Point]):
    getcontext().prec = 110
    pari.default("realprecision", 110)
    curve = pari(EllipticCurve(QQ, list(model)))
    raw = curve.ellheightmatrix([list(point) for point in basis], precision=384)
    dimension = len(basis)
    gram = tuple(
        tuple(Decimal(str(raw[i, j])) for j in range(dimension))
        for i in range(dimension)
    )
    maximum_asymmetry = max(
        abs(gram[i][j] - gram[j][i])
        for i in range(dimension)
        for j in range(dimension)
    )
    if maximum_asymmetry > Decimal("1e-90"):
        raise ArithmeticError(f"canonical-height Gram is asymmetric by {maximum_asymmetry}")
    return gram, maximum_asymmetry


def quadratic_decimal(gram, vector: Sequence[int]) -> Decimal:
    return sum(
        Decimal(vector[i]) * gram[i][j] * Decimal(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def rounded_gram(gram, scale: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
        for row in gram
    )


