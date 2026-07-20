"""Sparse Laurent bookkeeping and exact linear obstruction certificates."""
from __future__ import annotations
from dataclasses import dataclass
from functools import reduce
from math import gcd
import sympy as sp

@dataclass
class LinearCertificate:
    consistent: bool
    rank: int
    augmented_rank: int
    equations: int
    unknowns: int
    witness_row: list[str] | None

def term_exponents(term: sp.Expr, variables):
    powers = term.as_powers_dict()
    return tuple(int(powers.get(v, 0)) for v in variables)

def laurent_dict(expr: sp.Expr, variables):
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        exponents = term_exponents(term, variables)
        monomial = sp.prod(v**e for v, e in zip(variables, exponents))
        coefficient = sp.cancel(term/monomial)
        out[exponents] = sp.simplify(out.get(exponents, 0) + coefficient)
    return {e: c for e, c in out.items() if c != 0}

def forbidden_pole_equations(expr: sp.Expr, pole_variable, other_variables=()):
    variables = (pole_variable,) + tuple(other_variables)
    return [coefficient for exponents, coefficient in laurent_dict(expr, variables).items()
            if exponents[0] < 0]

def linear_certificate(equations, unknowns):
    equations = [sp.expand(e) for e in equations if sp.expand(e) != 0]
    if not equations:
        return LinearCertificate(True, 0, 0, 0, len(unknowns), None), sp.zeros(0, len(unknowns)), sp.zeros(0, 1)
    matrix, rhs = sp.linear_eq_to_matrix(equations, unknowns)
    rank = matrix.rank()
    augmented = matrix.row_join(rhs)
    augmented_rank = augmented.rank()
    witness = None
    if augmented_rank > rank:
        rref, _ = augmented.rref()
        for row in range(rref.rows):
            if all(rref[row, col] == 0 for col in range(matrix.cols)) and rref[row, -1] != 0:
                witness = [sp.sstr(rref[row, col]) for col in range(rref.cols)]
                break
    return LinearCertificate(rank == augmented_rank, rank, augmented_rank,
                             len(equations), len(unknowns), witness), matrix, rhs

def integer_matrix(matrix, rhs):
    entries = list(matrix) + list(rhs)
    denominators = [sp.denom(sp.cancel(e)) for e in entries]
    scale = sp.ilcm(*[int(d) for d in denominators]) if denominators else 1
    return matrix.applyfunc(lambda e: int(e*scale)), rhs.applyfunc(lambda e: int(e*scale)), scale

def rank_mod_prime(matrix: sp.Matrix, prime: int):
    rows = [[int(matrix[i, j]) % prime for j in range(matrix.cols)] for i in range(matrix.rows)]
    rank, col = 0, 0
    while rank < len(rows) and col < matrix.cols:
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][col], -1, prime)
        rows[rank] = [(v*inv) % prime for v in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                factor = rows[i][col]
                rows[i] = [(a-factor*b) % prime for a, b in zip(rows[i], rows[rank])]
        rank += 1
        col += 1
    return rank

