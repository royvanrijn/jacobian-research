"""Narrow PARI/GP backend for canonical heights and bounded SVP enumeration."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re
import subprocess
from typing import Sequence

from .elliptic import EllipticCurve, Point, point_complexity
from .integer import canonical_unoriented, content, rational_rank


def gp_rational(value: int | str | Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_point(point: Point) -> str:
    if point is None:
        return "[0]"
    return "[" + ",".join(gp_rational(value) for value in point) + "]"


def gp_vector(values: Sequence[int | str | Fraction]) -> str:
    return "[" + ",".join(gp_rational(value) for value in values) + "]"


def gp_matrix(rows: Sequence[Sequence[int]]) -> str:
    literal = "[" + ";".join(",".join(str(int(value)) for value in row) for row in rows) + "]"
    # A one-row GP literal is a vector, not a matrix.  Mat is harmless for
    # genuine matrix literals and removes that shape ambiguity.
    return f"Mat({literal})"


def run_gp(program: str, *, timeout: float = 300.0) -> list[str]:
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=True,
    )
    stderr = re.sub(r"\x1b\[[0-9;]*m", "", completed.stderr)
    errors = [
        line
        for line in stderr.splitlines()
        if line.strip()
        and "Warning: new" not in line
        and "increasing stack size" not in line
    ]
    if errors:
        raise RuntimeError("PARI/GP: " + "\n".join(errors))
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def primitive_column_closure(
    rows: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    """Saturate a full-column-rank matrix inside its ambient ``Z^r``.

    The returned matrix has the same rational column span.  PARI's double
    integer-kernel construction is exact and does not invoke Mordell--Weil
    saturation assumptions.
    """

    if not rows or not rows[0]:
        raise ValueError("closure input must be a nonempty matrix")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("closure matrix rows have inconsistent widths")
    if width == len(rows):
        # A full-rank rational subspace is the whole ambient space, whose
        # primitive closure has the standard basis.  PARI represents the
        # intervening zero-column kernel as a 0x0 matrix, so handle it here.
        return tuple(
            tuple(1 if row == column else 0 for column in range(width))
            for row in range(width)
        )
    program = f"""
B={gp_matrix(rows)};K=matkerint(B~);S=matkerint(K~);
print("SIZE|",matsize(S)[1],"|",matsize(S)[2]);print("BEGIN");
for(i=1,matsize(S)[1],for(j=1,matsize(S)[2],if(j>1,print1("|"));print1(S[i,j]));print());
print("END");
"""
    lines = run_gp(program)
    size = next(line.split("|")[1:] for line in lines if line.startswith("SIZE|"))
    height, returned_width = map(int, size)
    answer = tuple(
        tuple(int(value) for value in line.split("|"))
        for line in _parse_block(lines, "BEGIN", "END")
    )
    if height != len(rows) or returned_width != width or len(answer) != height:
        raise ArithmeticError(
            "primitive closure changed matrix rank or shape: "
            f"input={len(rows)}x{width}, output={height}x{returned_width}, rows={len(answer)}"
        )
    oriented = [list(row) for row in answer]
    for column in range(width):
        first = next((oriented[row][column] for row in range(height) if oriented[row][column]), 0)
        if first < 0:
            for row in range(height):
                oriented[row][column] *= -1
    return tuple(tuple(row) for row in oriented)


def row_embedding_smith_invariant_factors(
    rows: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """Return nonzero Smith factors of a full-row-rank embedding matrix.

    A ``k x n`` matrix records the images of a basis of ``Z^k`` in
    ``Z^n``.  Its image is primitive exactly when every returned factor is
    one.  PARI's integer Smith form is exact; no Mordell--Weil saturation
    assumption enters this calculation.
    """

    if not rows or not rows[0]:
        raise ValueError("embedding matrix must be nonempty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("embedding matrix rows have inconsistent widths")
    if rational_rank(rows) != len(rows):
        raise ValueError("embedding matrix is not full row rank")
    lines = run_gp(f"print(matsnf({gp_matrix(rows)}));")
    if len(lines) != 1 or not lines[0].startswith("["):
        raise ArithmeticError("PARI did not return Smith invariants")
    factors = tuple(
        abs(int(value.strip()))
        for value in lines[0].strip("[]").split(",")
        if int(value.strip()) != 0
    )
    if len(factors) != len(rows):
        raise ArithmeticError("embedding matrix Smith rank changed")
    return factors


def row_embedding_saturation_indices(
    matrices: Sequence[Sequence[Sequence[int]]],
    *,
    batch_size: int = 128,
    timeout: float = 300.0,
) -> tuple[int, ...]:
    """Batch the saturation index of each full-row-rank embedding matrix.

    For a ``k x n`` row matrix, the product of its nonzero Smith factors is
    the index of the generated row lattice in its primitive closure.  One GP
    process handles the declared batch so large candidate ledgers do not pay
    subprocess startup once per matrix.
    """

    matrices = tuple(
        tuple(tuple(map(int, row)) for row in matrix) for matrix in matrices
    )
    if not matrices:
        return ()
    row_counts = []
    for matrix in matrices:
        if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
            raise ValueError("batch embedding matrices must be nonempty and rectangular")
        row_counts.append(len(matrix))
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    parsed = []
    for start in range(0, len(matrices), batch_size):
        batch = matrices[start : start + batch_size]
        literals = ",".join(gp_matrix(matrix) for matrix in batch)
        program = f"""
V=[{literals}];
for(i=1,#V,s=matsnf(V[i]);z=1;r=0;for(j=1,#s,if(s[j],z*=abs(s[j]);r++));print(i,"|",r,"|",z));
"""
        lines = run_gp(program, timeout=timeout)
        batch_parsed = {}
        for line in lines:
            index, rank, value = line.split("|")
            batch_parsed[int(index)] = (int(rank), int(value))
        if set(batch_parsed) != set(range(1, len(batch) + 1)):
            raise ArithmeticError("PARI batch Smith output indices changed")
        parsed.extend(batch_parsed[index] for index in range(1, len(batch) + 1))
    if any(rank != row_counts[index] for index, (rank, _value) in enumerate(parsed)):
        raise ValueError("batch embedding matrix is not full row rank")
    return tuple(value for _rank, value in parsed)


def exact_rational_ranks(
    matrices: Sequence[Sequence[Sequence[int]]],
    *,
    batch_size: int = 128,
    timeout: float = 300.0,
) -> tuple[int, ...]:
    """Batch exact rational ranks of rectangular integer matrices in PARI."""

    matrices = tuple(
        tuple(tuple(map(int, row)) for row in matrix) for matrix in matrices
    )
    if not matrices:
        return ()
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    for matrix in matrices:
        if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
            raise ValueError("rank matrices must be nonempty and rectangular")
    answer = []
    for start in range(0, len(matrices), batch_size):
        batch = matrices[start : start + batch_size]
        literals = ",".join(gp_matrix(matrix) for matrix in batch)
        lines = run_gp(
            f'V=[{literals}];for(i=1,#V,print(i,"|",matrank(V[i])));',
            timeout=timeout,
        )
        parsed = {}
        for line in lines:
            index, rank = line.split("|")
            parsed[int(index)] = int(rank)
        if set(parsed) != set(range(1, len(batch) + 1)):
            raise ArithmeticError("PARI batch rank output indices changed")
        answer.extend(parsed[index] for index in range(1, len(batch) + 1))
    return tuple(answer)


def row_embedding_is_primitive(rows: Sequence[Sequence[int]]) -> bool:
    """Certify that a rectangular row embedding has primitive image."""

    return all(value == 1 for value in row_embedding_smith_invariant_factors(rows))


def _parse_block(lines: Sequence[str], begin: str, end: str) -> list[str]:
    start = lines.index(begin) + 1
    stop = lines.index(end, start)
    return list(lines[start:stop])


def height_gram(
    curve: EllipticCurve,
    points: Sequence[Point],
    *,
    digits: int = 100,
    timeout: float = 300.0,
) -> list[list[str]]:
    if not points or any(point is None for point in points):
        raise ValueError("height Gram input must contain finite points")
    if any(not curve.is_on_curve(point) for point in points):
        raise ValueError("height Gram input contains an off-curve point")
    program = f"""
default(realprecision,{int(digits)});
E=ellinit({gp_vector(curve.coefficients)});P=[{','.join(map(gp_point, points))}];
H=ellheightmatrix(E,P);print("BEGIN");
for(i=1,#P,for(j=1,#P,if(j>1,print1("|"));print1(H[i,j]));print());
print("END");
"""
    lines = run_gp(program, timeout=timeout)
    return [line.split("|") for line in _parse_block(lines, "BEGIN", "END")]


@dataclass(frozen=True)
class ExactEmbedding:
    """Columns expressing target points in a displayed independent basis."""

    columns: tuple[tuple[int, ...], ...]
    max_abs_coordinate: int
    nonzero_coordinates: int
    numerical_residual_max: str

    def rows(self) -> tuple[tuple[int, ...], ...]:
        if not self.columns:
            return ()
        return tuple(tuple(column[row] for column in self.columns) for row in range(len(self.columns[0])))

    def smith_invariant_factors(self) -> tuple[int, ...]:
        """Return the nonzero Smith factors of the column embedding."""

        rows = self.rows()
        if not rows:
            return ()
        return row_embedding_smith_invariant_factors(rows)

    def is_primitive(self) -> bool:
        return all(value == 1 for value in self.smith_invariant_factors())


def recover_exact_embedding(
    curve: EllipticCurve,
    basis: Sequence[Point],
    targets: Sequence[Point],
    *,
    digits: int = 120,
    timeout: float = 300.0,
) -> ExactEmbedding:
    """Recover an integral embedding by a height-dual solve and exact replay.

    This is a certificate constructor once both point lists are supplied.  It
    is not a blind latent-subgroup selector and must not be scored as one.
    """

    if not basis or not targets or any(point is None for point in (*basis, *targets)):
        raise ValueError("embedding inputs must be nonempty finite point lists")
    if any(not curve.is_on_curve(point) for point in (*basis, *targets)):
        raise ValueError("embedding input contains an off-curve point")
    point_text = ",".join(map(gp_point, basis))
    target_text = ",".join(map(gp_point, targets))
    n_basis = len(basis)
    n_targets = len(targets)
    program = f"""
default(realprecision,{int(digits)});
E=ellinit({gp_vector(curve.coefficients)});P=[{point_text}];Q=[{target_text}];
A=concat(P,Q);H=ellheightmatrix(E,A);HP=H[1..#P,1..#P];
X=HP^-1*H[1..#P,#P+1..#A];Z=matrix(#P,#Q,i,j,round(X[i,j]));
ok=1;for(j=1,#Q,T=[0];for(i=1,#P,if(Z[i,j],T=elladd(E,T,ellmul(E,P[i],Z[i,j]))));if(ellsub(E,T,Q[j])!=[0],ok=0));
print("EXACT|",ok);
print("RESIDUAL|",vecmax(abs(X-Z)));
print("BEGIN");for(j=1,#Q,for(i=1,#P,if(i>1,print1("|"));print1(Z[i,j]));print());print("END");
"""
    lines = run_gp(program, timeout=timeout)
    residual = next(line.split("|", 1)[1] for line in lines if line.startswith("RESIDUAL|"))
    exact = next(line.split("|", 1)[1] for line in lines if line.startswith("EXACT|"))
    if exact != "1":
        raise ArithmeticError("height-dual coordinates failed PARI exact group-law replay")
    columns = tuple(
        tuple(int(value) for value in line.split("|"))
        for line in _parse_block(lines, "BEGIN", "END")
    )
    if len(columns) != n_targets or any(len(column) != n_basis for column in columns):
        raise ArithmeticError("PARI returned an embedding matrix of the wrong shape")
    # Retain an implementation-independent replay when its rational arithmetic
    # is small enough for a routine calibration run.
    work = sum(abs(value) for column in columns for value in column)
    if work <= 500:
        for column, target in zip(columns, targets):
            if curve.linear_combination(basis, column) != target:
                raise ArithmeticError("height-dual coordinates failed Python exact group-law replay")
    return ExactEmbedding(
        columns=columns,
        max_abs_coordinate=max(abs(value) for column in columns for value in column),
        nonzero_coordinates=sum(value != 0 for column in columns for value in column),
        numerical_residual_max=residual,
    )


@dataclass(frozen=True)
class ShortVectorRecord:
    coordinates: tuple[int, ...]
    canonical_height: str
    point: Point
    arithmetic: dict[str, int | bool]

    def to_record(self) -> dict[str, object]:
        return {
            "coordinates": list(self.coordinates),
            "canonical_height": self.canonical_height,
            "point": None
            if self.point is None
            else [str(self.point[0]), str(self.point[1])],
            "arithmetic": self.arithmetic,
        }


def enumerate_short_vectors(
    curve: EllipticCurve,
    basis: Sequence[Point],
    *,
    height_bound: float,
    digits: int = 100,
    maximum_lines: int = 100_000,
    materialize_points: bool = True,
    timeout: float = 600.0,
) -> tuple[ShortVectorRecord, ...]:
    """Enumerate a complete primitive unoriented shell through a height bound."""

    if not basis or any(point is None for point in basis):
        raise ValueError("short-vector basis must contain finite points")
    if any(not curve.is_on_curve(point) for point in basis):
        raise ValueError("short-vector basis contains an off-curve point")
    materialize = """
T=[0];for(i=1,#P,if(V[i,j],T=elladd(E,T,ellmul(E,P[i],V[i,j]))));
print1("|",T[1],"|",T[2]);
""" if materialize_points else ""
    program = f"""
default(parisizemax,4000000000);default(parisize,500000000);
default(realprecision,{int(digits)});
E=ellinit({gp_vector(curve.coefficients)});P=[{','.join(map(gp_point, basis))}];
H=ellheightmatrix(E,P);U=qflllgram(H);R=U~*H*U;
Q=qfminim(R,{float(height_bound):.17g},{int(maximum_lines)},2);V=U*Q[3];
print("COUNT|",matsize(V)[2]);print("BEGIN");
for(j=1,matsize(V)[2],{{for(i=1,matsize(V)[1],if(i>1,print1(","));print1(V[i,j]));
print1("|",V[,j]~*H*V[,j]);{materialize}print();}});
print("END");
"""
    lines = run_gp(program, timeout=timeout)
    _raw_count = int(next(line.split("|", 1)[1] for line in lines if line.startswith("COUNT|")))
    records: dict[tuple[int, ...], ShortVectorRecord] = {}
    for line in _parse_block(lines, "BEGIN", "END"):
        fields = line.split("|")
        raw = tuple(map(int, fields[0].split(",")))
        if content(raw) != 1:
            continue
        coordinates = canonical_unoriented(raw)
        point: Point = None
        if materialize_points:
            point = Fraction(fields[2]), Fraction(fields[3])
            # Canonicalizing a negative PARI representative changes the point.
            if coordinates != raw:
                point = curve.negate(point)
        arithmetic = point_complexity(point) if materialize_points else {}
        record = ShortVectorRecord(coordinates, fields[1], point, arithmetic)
        if coordinates in records:
            raise ArithmeticError("PARI returned both orientations of a primitive vector")
        records[coordinates] = record
    return tuple(sorted(records.values(), key=lambda item: (float(item.canonical_height), item.coordinates)))
