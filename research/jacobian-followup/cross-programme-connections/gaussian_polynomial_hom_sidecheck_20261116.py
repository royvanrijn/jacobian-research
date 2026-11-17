"""Exact Q[t]<D> polynomial Hom from a finite module presentation.

Coordinate annihilators MUST carry verified left relation witnesses.
This computes Hom of the SUPPLIED module. Gaussian safety additionally
requires exact source identification or a proved surjection covering all
amplitudes plus the one-dimensional time-zero quotient certificate.
"""
from __future__ import annotations
import json
import sys
from math import comb, factorial
from pathlib import Path
import sympy as sp

Z = sp.Symbol("z")


def parse_op(terms):
    result = {}
    for a, b, c in terms:
        if not isinstance(a, int) or not isinstance(b, int) or a < 0 or b < 0:
            raise ValueError("Weyl exponents must be nonnegative integers")
        v = sp.Rational(c)
        result[a, b] = result.get((a, b), sp.S.Zero) + v
    return {k: v for k, v in result.items() if v}


def op_terms(op):
    return [[a, b, str(c)] for (a, b), c in sorted(op.items()) if c]


def add(*ops):
    out = {}
    for op in ops:
        for mon, c in op.items():
            out[mon] = out.get(mon, sp.S.Zero) + c
    return {k: v for k, v in out.items() if v}


def scale(op, c):
    return {k: sp.Rational(c)*v for k, v in op.items() if c*v}


def falling(n, j):
    if j > n:
        return 0
    return factorial(n)//factorial(n-j)


def mul(left, right):
    """Exact PBW multiplication, with D*t=t*D+1."""
    out = {}
    for (a, b), c in left.items():
        for (u, v), e in right.items():
            for j in range(min(b, u)+1):
                mon = (a+u-j, b+v-j)
                out[mon] = out.get(mon, sp.S.Zero) + c*e*comb(b, j)*falling(u, j)
    return {k: v for k, v in out.items() if v}


def degree_certificate(op):
    if not op:
        raise ValueError("Zero is not a coordinate annihilator")
    shift = max(a-b for a, b in op)
    expr = sp.S.Zero
    for (a, b), c in op.items():
        if a-b == shift:
            f = sp.prod(Z-j for j in range(b))
            expr += c*f
    poly = sp.Poly(sp.expand(expr), Z, domain=sp.QQ)
    if poly.is_zero:
        raise AssertionError("Distinct falling-factorial degrees cannot cancel identically")
    roots = sorted(int(x) for x in sp.polys.polytools.ground_roots(poly)
                   if x.is_Integer and x >= 0)
    return dict(shift=shift, indicial_polynomial=str(poly.as_expr()),
                nonnegative_integer_roots=roots,
                degree_bound=max(roots, default=-1))


def apply_polynomial(op, coefficients):
    """Output coefficients of op applied to sum c[n]*t^n."""
    result = {}
    for n, c in enumerate(coefficients):
        if not c:
            continue
        for (a, b), q in op.items():
            if n >= b:
                k = n-b+a
                result[k] = result.get(k, sp.S.Zero) + c*q*falling(n, b)
    return {k: sp.cancel(v) for k, v in result.items() if v}


def serialize_polynomial(coefficients):
    t = sp.Symbol("t")
    return dict(coefficients=[str(x) for x in coefficients],
                expression=str(sum(c*t**n for n, c in enumerate(coefficients))))


def solve(data, max_unknowns=300):
    rank = data["generator_count"]
    if not isinstance(rank, int) or rank < 1:
        raise ValueError("generator_count must be a positive integer")
    rows = [[parse_op(op) for op in row] for row in data["relations"]]
    if any(len(row) != rank for row in rows):
        raise ValueError("Every relation must have generator_count entries")
    certs = data["coordinate_annihilators"]
    if len(certs) != rank:
        raise ValueError("One certified annihilator is required per generator")
    bounds = []
    for i, cert in enumerate(certs):
        op = parse_op(cert["operator"])
        witnesses = [parse_op(w) for w in cert["left_witnesses"]]
        if len(witnesses) != len(rows):
            raise ValueError("Each annihilator needs one LEFT witness per relation")
        for coordinate in range(rank):
            actual = add(*(mul(w, row[coordinate]) for w, row in zip(witnesses, rows)))
            expected = op if coordinate == i else {}
            if actual != expected:
                raise ValueError(f"False coordinate-annihilator witness at ({i},{coordinate})")
        bounds.append(degree_certificate(op))
    unknowns = [(i, n) for i, cert in enumerate(bounds)
                for n in range(cert["degree_bound"]+1)]
    base_result = dict(coordinate_degree_certificates=bounds,
                       verified_coordinate_annihilator_identities=True,
                       unknown_coefficient_count=len(unknowns),
                       gaussian_safety_claimed=False,
                       source_contract="Exact source, or all-amplitude surjection plus dim(N/tN)=1, must be proved separately")
    if max_unknowns is not None and len(unknowns) > max_unknowns:
        return dict(base_result, outcome="UNKNOWN", reason="Explicit coefficient budget exceeded")
    index = {key: j for j, key in enumerate(unknowns)}
    equations = []
    for row in rows:
        out = {}
        for i, op in enumerate(row):
            for n in range(bounds[i]["degree_bound"]+1):
                for (a, b), c in op.items():
                    if n >= b:
                        v = out.setdefault(n-b+a, {})
                        col = index[i, n]
                        v[col] = v.get(col, sp.S.Zero) + c*falling(n, b)
        for k in sorted(out):
            equation = [out[k].get(j, sp.S.Zero) for j in range(len(unknowns))]
            if any(equation):
                equations.append(equation)
    matrix = sp.Matrix(equations) if equations else sp.zeros(0, len(unknowns))
    kernel = matrix.nullspace() if unknowns else []
    anchor = [parse_op(op) for op in data.get("unit_anchor", [[[0, 0, "1"]]]+[[]]*(rank-1))]
    if len(anchor) != rank:
        raise ValueError("unit_anchor must have generator_count entries")
    basis, values = [], []
    for vector in kernel:
        polys = [[vector[index[i, n]] for n in range(bounds[i]["degree_bound"]+1)]
                 for i in range(rank)]
        # Replay every defining relation on the returned literal polynomials.
        for row in rows:
            total = {}
            for op, poly in zip(row, polys):
                for degree, value in apply_polynomial(op, poly).items():
                    total[degree] = total.get(degree, sp.S.Zero)+value
            if any(sp.cancel(c) for c in total.values()):
                raise AssertionError("Returned map does not satisfy the defining relations")
        value = sum(apply_polynomial(op, poly).get(0, sp.S.Zero)
                    for op, poly in zip(anchor, polys))
        values.append(value)
        basis.append([serialize_polynomial(poly) for poly in polys])
    normalized = None
    if len(kernel) == 1 and values[0] != 0:
        normalized = [
            serialize_polynomial([sp.Rational(x)/values[0] for x in poly["coefficients"]])
            for poly in basis[0]]
    return dict(base_result, outcome="COMPLETE", polynomial_hom_dimension=len(kernel),
                polynomial_hom_basis=basis,
                unit_values_at_zero=[str(v) for v in values],
                normalized_unique_map=normalized,
                defining_relations_replayed=True,
                linear_equation_count=matrix.rows)



def find_coordinate_annihilators(data, max_weyl_degree=6, max_columns=500):
    """Enumerate exact finite LEFT witnesses. Uncapped search terminates
    for holonomic D1 modules; a cap yields UNKNOWN, never a zero Hom claim.
    """
    rank = data["generator_count"]
    if not isinstance(rank, int) or rank < 1:
        raise ValueError("generator_count must be positive")
    rows = [[parse_op(op) for op in row] for row in data["relations"]]
    if any(len(row) != rank for row in rows):
        raise ValueError("Wrong relation width")
    result = []
    stages = []
    for coordinate in range(rank):
        degree = 0
        while max_weyl_degree is None or degree <= max_weyl_degree:
            mons = [(a, total-a) for total in range(degree+1) for a in range(total+1)]
            column_count = len(mons)*(len(rows)+1)
            if max_columns is not None and column_count > max_columns:
                return dict(outcome="UNKNOWN", reason="Witness column budget exceeded",
                            completed_coordinates=coordinate, searched_degree=degree)
            columns = []
            # Unknown annihilator A: coefficient of -A*e_coordinate.
            for mon in mons:
                columns.append({(coordinate, *mon): -sp.S.One})
            # Unknown W_j: coefficients of W_j times each supplied relation.
            for row in rows:
                for mon in mons:
                    image = {}
                    for i, op in enumerate(row):
                        for (a,b), value in mul({mon: sp.S.One}, op).items():
                            image[i,a,b] = value
                    columns.append(image)
            equations = sorted(set().union(*(column.keys() for column in columns)))
            entries = {(i,j): column[key]
                       for j,column in enumerate(columns)
                       for i,key in enumerate(equations) if key in column}
            matrix = sp.MutableSparseMatrix(len(equations), len(columns), entries)
            found = None
            for vector in matrix.nullspace():
                if any(vector[j] for j in range(len(mons))):
                    found = vector
                    break
            if found is not None:
                lead = next(found[j] for j in range(len(mons)) if found[j])
                coeffs = [x/lead for x in found]
                op = {mon: coeffs[j] for j,mon in enumerate(mons) if coeffs[j]}
                witnesses = []
                for j in range(len(rows)):
                    start = (j+1)*len(mons)
                    w = {mon: coeffs[start+k] for k,mon in enumerate(mons) if coeffs[start+k]}
                    witnesses.append(op_terms(w))
                result.append(dict(operator=op_terms(op),left_witnesses=witnesses))
                stages.append(degree)
                break
            degree += 1
        else:
            return dict(outcome="UNKNOWN", reason="Witness Weyl-degree budget exceeded",
                        completed_coordinates=coordinate, searched_degree=max_weyl_degree)
    return dict(outcome="COMPLETE", coordinate_annihilators=result,
                witness_search_degrees=stages)


def solve_presentation(data, max_unknowns=300, max_weyl_degree=6, max_witness_columns=500):
    if "coordinate_annihilators" in data:
        return solve(data, max_unknowns=max_unknowns)
    found = find_coordinate_annihilators(data, max_weyl_degree, max_witness_columns)
    if found["outcome"] != "COMPLETE":
        return dict(found, gaussian_safety_claimed=False)
    enriched = dict(data, coordinate_annihilators=found["coordinate_annihilators"])
    result = solve(enriched, max_unknowns=max_unknowns)
    result["generated_coordinate_annihilators"] = found["coordinate_annihilators"]
    result["witness_search_degrees"] = found["witness_search_degrees"]
    return result


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: gaussian_polynomial_hom_sidecheck_20261116.py input.json result.json")
    data = json.loads(Path(sys.argv[1]).read_text())
    result = solve_presentation(data, max_unknowns=data.get("max_unknowns", 300),
                                max_weyl_degree=data.get("max_weyl_degree", 6),
                                max_witness_columns=data.get("max_witness_columns", 500))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
