"""Exact incremental binary linear algebra with retained provenance."""

from dataclasses import dataclass


def pack(row):
    bits = tuple(row)
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError("nonbinary matrix entry")
    return sum(int(bit) << i for i, bit in enumerate(bits))


def unpack(mask, width):
    if type(mask) is not int or mask < 0 or mask >> width:
        raise ValueError("mask outside the declared binary space")
    return tuple((mask >> i) & 1 for i in range(width))


@dataclass(frozen=True)
class BinaryBasis:
    width: int
    pivots: tuple[tuple[int, int, int], ...] = ()
    column_count: int = 0

    @property
    def rank(self):
        return len(self.pivots)

    def reduce(self, value):
        unpack(value, self.width)
        combination = 0
        for pivot, vector, provenance in self.pivots:
            if (value >> pivot) & 1:
                value ^= vector
                combination ^= provenance
        return value, combination

    def append(self, value):
        residual, combination = self.reduce(value)
        provenance = combination ^ (1 << self.column_count)
        pivots = self.pivots
        if residual:
            pivots = tuple(sorted((*pivots, (residual.bit_length()-1, residual, provenance)), reverse=True))
        return BinaryBasis(self.width, pivots, self.column_count+1), (None if residual else provenance)


def kernel_masks(rows, *, width=None):
    rows = tuple(tuple(row) for row in rows)
    if width is None:
        width = len(rows[0]) if rows else 0
    if any(len(row) != width for row in rows):
        raise ValueError("inconsistent matrix width")
    basis = BinaryBasis(width)
    dependencies = []
    for row in rows:
        basis, dependency = basis.append(pack(row))
        if dependency is not None:
            dependencies.append(dependency)
    return tuple(dependencies)


def combine(mask, rows):
    unpack(mask, len(rows))
    value = 0
    for i, row in enumerate(rows):
        if (mask >> i) & 1:
            value ^= row
    return value


def quotient_rows(subspace, rows, *, width=None):
    """Coordinates in a deterministically extended quotient basis, in O(n r²).

    The ordering agrees with first-independent-row extension, without trying
    any of its 2^r combinations. The supplied subspace must be independent.
    """
    subspace, rows = tuple(map(tuple, subspace)), tuple(map(tuple, rows))
    if width is None:
        width = len((subspace + rows)[0]) if subspace or rows else 0
    basis = BinaryBasis(width)
    coordinates = []
    for i, row in enumerate(subspace + rows):
        if len(row) != width:
            raise ValueError("inconsistent quotient-map width")
        residual, combination = basis.reduce(pack(row))
        if residual:
            combination = 1 << basis.rank
            basis, _ = basis.append(pack(row))
        elif i < len(subspace):
            raise ValueError("dependent local image generators")
        if i >= len(subspace):
            coordinates.append(combination >> len(subspace))
    dimension = basis.rank - len(subspace)
    return [list(unpack(c, dimension)) for c in coordinates], dimension
