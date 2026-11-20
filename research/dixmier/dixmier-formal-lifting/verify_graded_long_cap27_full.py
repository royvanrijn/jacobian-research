#!/usr/bin/env python3
"""Replay the full-parity canonical Long degree-27 exclusion."""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import factorial, isqrt
import json
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import flint
import check_graded_long_cap23 as columns
import check_graded_long_cap23_full_jet as moyal
import check_graded_long_cap26_full_jet as parameters
import check_long_graph_curvature as polynomials


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "dixmier-audit" / "verify_long_pbw.py"
CAP = 27
PRIME = 1000003
SPATIAL = ("X", "Q", "P", "Z")


def rational(c):
    return Fraction(int(c.numerator), int(c.denominator))


def mod(c):
    return c.numerator * pow(c.denominator, -1, PRIME) % PRIME


def coefficient_basis(order):
    return sorted(
        (index, (a, b, c, d))
        for index, weight in enumerate(columns.WEIGHTS)
        for b in range(CAP + 1)
        for c in range(CAP + 1 - b)
        for d in range(CAP + 1 - b - c)
        if (a := b + 4*c + 2*d + 3*order - weight) >= 0
        and a + b + c + d <= CAP
    )


def leading_columns(matrix, rank):
    return [next(j for j in range(matrix.ncols()) if matrix[i, j])
            for i in range(rank)]


class LinearElimination:
    """Pivot on constants and check the entire characteristic-zero kernel."""

    def __init__(self, images, order):
        self.basis = coefficient_basis(order)
        self.columns = columns.linear_columns(polynomials, images, self.basis)
        rows = sorted(set().union(*self.columns))
        row_index = {row: i for i, row in enumerate(rows)}
        n, m = len(rows), len(self.basis)
        transpose = flint.nmod_mat(m, n, PRIME)
        for i, column in enumerate(self.columns):
            for row, c in column.items():
                transpose[i, row_index[row]] = mod(c)
        reduced, self.rank = transpose.rref()
        keep = leading_columns(reduced, self.rank)
        restricted = flint.nmod_mat(
            self.rank, m, [transpose[j, i] for i in keep for j in range(m)], PRIME
        )
        reduced, rank_check = restricted.rref()
        assert rank_check == self.rank
        self.pivots = leading_columns(reduced, self.rank)
        self.free = sorted(set(range(m)) - set(self.pivots))
        self.pivot_rows = [rows[i] for i in keep]

        square = flint.fmpq_mat(self.rank, self.rank)
        for i, row in enumerate(self.pivot_rows):
            for j, column in enumerate(self.pivots):
                c = self.columns[column].get(row)
                if c:
                    square[i, j] = flint.fmpq(c.numerator, c.denominator)
        self.inverse = square.inv()
        self.kernel = []
        for j, solution in zip(self.free, self.particulars(
                [{row: -c for row, c in self.columns[j].items()} for j in self.free])):
            solution[j] = Fraction(1)
            assert not columns.product(self.columns, solution)
            self.kernel.append(self.unpack(solution))

    def particulars(self, right_sides):
        matrix = flint.fmpq_mat(self.rank, len(right_sides))
        for j, side in enumerate(right_sides):
            for i, row in enumerate(self.pivot_rows):
                if c := side.get(row):
                    matrix[i, j] = flint.fmpq(c.numerator, c.denominator)
        solution = self.inverse * matrix
        result = []
        for j in range(len(right_sides)):
            vector = [Fraction(0)] * len(self.basis)
            for i, pivot in enumerate(self.pivots):
                vector[pivot] = rational(solution[i, j])
            result.append(vector)
        return result

    def unpack(self, vector):
        image = [{} for _ in SPATIAL]
        for (index, exponents), coefficient in zip(self.basis, vector):
            if coefficient:
                image[index][exponents] = coefficient
        return image

    def solve(self, rhs):
        """Return an affine family and all equations on its earlier parameters."""
        frame = {}
        residuals = {}
        for parameter, vector in zip(rhs, self.particulars(list(rhs.values()))):
            image = self.unpack(vector)
            if any(image):
                frame[parameter] = image
            defect = moyal.difference(rhs[parameter], columns.product(self.columns, vector))
            for row, coefficient in defect.items():
                residuals.setdefault(row, {})[parameter] = coefficient
        return frame, list(row_span(residuals.values()).values())


def minus_multiple(row, pivot_row, multiplier):
    for monomial, coefficient in pivot_row.items():
        value = row.get(monomial, 0) - multiplier * coefficient
        if value:
            row[monomial] = value
        else:
            row.pop(monomial, None)


def row_span(equations):
    """Reduce only over Q, without dividing by a parameter."""
    basis = {}
    for equation in equations:
        row = dict(equation)
        while row:
            pivot = max(row)
            if pivot in basis:
                minus_multiple(row, basis[pivot], row[pivot])
            else:
                c = row[pivot]
                basis[pivot] = {m: q/c for m, q in row.items()}
                break
    for pivot in sorted(basis):
        for other, row in basis.items():
            if other != pivot and (c := row.get(pivot)):
                minus_multiple(row, basis[pivot], c)
    return basis


def contained_in_span(equation, span):
    row = dict(equation)
    while row:
        pivot = max(row)
        if pivot not in span:
            return False
        minus_multiple(row, span[pivot], row[pivot])
    return True


def rhs_coefficient(names, frames, order):
    """Accelerated evaluation of the complete coefficient, in all six pairs."""
    ctx = flint.fmpq_mpoly_ctx.get((*SPATIAL, *names))
    images = [
        [ctx.from_dict({space + mon: flint.fmpq(q.numerator, q.denominator)
                        for mon, tup in frame.items() for space, q in tup[i].items()})
         for i in range(4)]
        for frame in frames
    ]

    @lru_cache(None)
    def derivative(k, i, powers):
        poly = images[k][i]
        for variable, count in enumerate(powers):
            for _ in range(count):
                poly = poly.derivative(variable)
        return poly

    result = {}
    for pair, (i, j) in enumerate(columns.PAIRS):
        poly = ctx.constant(0)
        for a, b, ell, _ in moyal.rhs_indices(order, len(frames)):
            for k in range(ell + 1):
                for l in range(ell - k + 1):
                    for m in range(ell - k - l + 1):
                        n = ell - k - l - m
                        left = derivative(a, i, (l, n, k, m))
                        right = derivative(b, j, (k, m, l, n))
                        if left and right:
                            denominator = (2**(ell-1)*factorial(k)*factorial(l)
                                           * factorial(m)*factorial(n))
                            poly -= flint.fmpq((-1)**(l+n), denominator)*left*right
        for exponents, coefficient in poly.to_dict().items():
            exponents = tuple(map(int, exponents))
            result.setdefault(exponents[4:], {})[(pair, exponents[:4])] = rational(coefficient)
    return result


def compare_moyal_at_point(names, frames, order, fast):
    """Compare with the independent Fraction derivative engine at the same input."""
    values = tuple(Fraction((-1)**i) for i, _ in enumerate(names))
    jets = moyal.evaluate_frames(polynomials, frames, values)
    direct_frames = [{(): [jet[j] for jet in jets]} for j in range(len(frames))]
    independent = moyal.rhs_coefficients(columns, polynomials, direct_frames, order).get((), {})
    combined = {}
    for exponents, row in fast.items():
        scale = (-1)**sum(exponents[1::2])
        for key, c in row.items():
            combined[key] = combined.get(key, 0) + scale*c
    assert independent == {k: v for k, v in combined.items() if v}, order


def term(names, **powers):
    return tuple(powers.get(name, 0) for name in names)


def expression(names, *items):
    return {term(names, **powers): Fraction(coefficient) for coefficient, powers in items}


def divisors(basis):
    """Normal forms using monic polynomials already certified as consequences."""
    @lru_cache(None)
    def reduce(monomial):
        for row in basis:
            lead = row[0][0]
            extra = tuple(a-b for a, b in zip(monomial, lead))
            if min(extra, default=0) < 0:
                continue
            output = {}
            for m, c in row[1:]:
                minus_multiple(output, reduce(moyal.plus(m, extra)), c)
            return output
        return {monomial: Fraction(1)}
    return reduce


def reduce_polynomial(poly, reduce):
    out = {}
    for m, c in poly.items():
        minus_multiple(out, reduce(m), -c)
    return out


def reduce_frames(frames, reduce):
    result = []
    for frame in frames:
        new = {}
        for mon, images in frame.items():
            for m, c in reduce(mon).items():
                moyal.add_tuple(polynomials, new, m,
                                moyal.tuple_scale(polynomials, images, c))
        result.append(new)
    return result


def reduce_rhs(rhs, reduce):
    result = {}
    for mon, row in rhs.items():
        for m, c in reduce(mon).items():
            minus_multiple(result.setdefault(m, {}), row, -c)
    return {m: row for m, row in result.items() if row}


def change_variables(names, new, frames, equations, replacements=None):
    replacements = replacements or {}
    monomials = {
        name: {term(new, **{name: 1}): Fraction(1)} if name in new else {}
        for name in names
    }
    monomials.update(replacements)
    images = [monomials[name] for name in names]
    moved = [parameters.substitute(polynomials, moyal, f, images, len(new)) for f in frames]
    equations = [
        {m: image[0][(0, 0, 0, 0)] for m, image in parameters.substitute(
            polynomials, moyal, {m: [{(0, 0, 0, 0): c}, {}, {}, {}]
                                  for m, c in p.items()}, images, len(new)).items()}
        for p in equations
    ]
    return moved, list(row_span(equations).values())


def groebner_with_lift(names, equations, directory, singular):
    """Check polynomial ideal membership by literal matrix multiplication."""
    def render(p):
        return "+".join(
            f"({c})" + "".join(f"*{n}^{e}" for n, e in zip(names, m) if e)
            for m, c in p.items()
        )

    script = directory / "elimination.sing"
    script.write_text(
        f'ring r=0,({",".join(names)}),dp;\n'
        'option(redSB);\n'
        'ideal I=' + ',\n'.join(render(p) for p in equations if p) + ';\n'
        'ideal G=slimgb(I);\n'
        'matrix M=lift(I,G);\n'
        'matrix Irow[1][size(I)]=I;\n'
        'matrix Grow[1][size(G)]=G;\n'
        'if (Irow*M!=Grow) {print("BAD LIFT");quit;}\n'
        'ideal reduction=reduce(I,G);\n'
        'if (size(reduction)!=0) {print("BAD REDUCTION");quit;}\n'
        'print("LIFT CHECKED");\npoly p; int j;\n'
        'for(j=1;j<=size(G);j++){p=G[j];'
        'while(p!=0){print(string(leadcoef(p))+" "+string(leadexp(p)));p=p-lead(p);}'
        'print("END");}\nquit;\n'
    )
    result = run([singular, "-q", str(script)], check=True, capture_output=True, text=True)
    entries = result.stdout.splitlines()
    assert entries[0] == "LIFT CHECKED" and "?" not in result.stdout, result.stdout[:500]
    basis = []
    row = []
    for line in entries[1:]:
        if line == "END":
            if row:
                normalizer = row[0][1]
                basis.append([(m, c/normalizer) for m, c in row])
            row = []
        else:
            coefficient, exponents = line.split()
            row.append((tuple(int(e) for e in exponents.split(",")), Fraction(coefficient)))
    assert not row
    return basis


def first_kernel_generators(images):
    add, multiply, scale = polynomials.add, polynomials.multiply, polynomials.scale

    def power(poly, degree):
        return columns.power(polynomials, poly, degree)

    x, q, p, z = [{term(SPATIAL, **{name: 1}): Fraction(1)} for name in SPATIAL]
    one = {term(SPATIAL): Fraction(1)}
    r, t_image, d_image, s_image = images
    t = multiply(x, q)
    rs = multiply(r, s_image)
    b = add(scale(multiply(power(x, 2), p), 3),
            scale(multiply(add(one, scale(t, -3)), z), 2))
    a0 = add(q, scale(multiply(x, power(q, 2)), 3),
             scale(multiply(x, b), Fraction(-1, 3)))
    y = add(multiply(r, a0), scale(t, -3))
    ty = multiply(t, y)
    ty2 = multiply(t, power(y, 2))
    t2 = power(t, 2)
    t2y = multiply(t2, y)
    h4 = add(power(y, 4), scale(ty2, 4), scale(ty, Fraction(-16, 3)), scale(t2, 10))
    shared = add(scale(power(y, 5), 9), scale(multiply(t, power(y, 3)), 45),
                 scale(ty2, -60))
    n5 = add(shared, scale(t2y, 45), scale(ty, 110), scale(t2, -720))
    h6 = add(multiply(power(r, 4), d_image), scale(n5, Fraction(1, 25)))
    h5 = add(shared, scale(t2y, 135), scale(ty, -10), scale(t2, -555))
    hamiltonians = {
        "_t": t, "_rs": rs, "a": y, "b": multiply(power(r, 2), t_image),
        "c": h4, "d": h6, "e": multiply(power(x, 2), z), "f": h5,
        "eta": power(rs, 2), "j": multiply(power(x, 4), p),
        "k": multiply(multiply(power(x, 3), q), z),
    }
    return {name: [polynomials.bidifferential(h, f, 1) for f in images]
            for name, h in hamiltonians.items()}


def solve_all(images, matrices, directory, singular):
    canonical = first_kernel_generators(images)
    first = matrices[0]

    def rank_in_free_coordinates(jets):
        matrix = flint.fmpq_mat(len(jets), len(first.free))
        for i, jet in enumerate(jets):
            vector = columns.vector(first.basis, jet)
            for j, k in enumerate(first.free):
                c = vector[k]
                matrix[i, j] = flint.fmpq(c.numerator, c.denominator)
        return matrix.rref()[1]

    assert rank_in_free_coordinates(list(canonical.values())) == 11
    for jet in first.kernel:
        if rank_in_free_coordinates([*canonical.values(), jet]) > len(canonical):
            canonical[f"x{len(canonical)-11}"] = jet
    assert rank_in_free_coordinates(list(canonical.values())) == len(first.free) == 14
    assert all(not columns.product(first.columns, columns.vector(first.basis, u))
               for u in canonical.values())

    # The two exact linear shears from the cap-26 proof supply these gauges.
    del canonical["_t"], canonical["_rs"]
    names = list(canonical)
    frames = [{term(names): images},
              {term(names, **{name: 1}): jet for name, jet in canonical.items()}]
    raw3 = rhs_coefficient(names, frames, 3)
    compare_moyal_at_point(names, frames, 3, raw3)
    _, equations3 = matrices[1].solve(raw3)
    span = row_span(equations3)
    for name in ("x0", "x1", "x2"):
        assert contained_in_span(expression(names, (1, {name: 2})), span)
    new = [n for n in names if n not in {"x0", "x1", "x2"}]
    frames, eq = change_variables(names, new, frames, equations3)
    span = row_span(eq)
    assert contained_in_span(expression(new, (4, {"j": 2}), (4, {"j": 1,"k": 1}),
                                         (1, {"k": 2})), span)
    conic = ((1, {"d": 2}), (-250, {"d": 1, "f": 1}), (35625, {"f": 2}))
    assert contained_in_span(expression(new, *conic), span)
    old = new
    names = "a b c d e f eta z u v gamma delta chi psi kappa sigma".split()
    frames, _ = change_variables(old, names, frames, [], {
        "j": expression(names, (3, {"z": 1})),
        "k": expression(names, (-6, {"z": 1})),
        "e": expression(names, (1, {"e": 1}), (2, {"z": 1})),
    })
    normal_conic = divisors([list(expression(names, *conic).items())])
    rhs3 = reduce_rhs(rhs_coefficient(names, frames, 3), normal_conic)
    frame2, defects = matrices[1].solve(rhs3)
    assert not defects

    # Use the retained positive even particular at the constant parameter.
    even_witness = json.loads((HERE / "graded_long_cap27_even_witness.json").read_text())
    even = even_witness["first_family"]
    particular = parameters.parse(columns, even["particular"])
    assert columns.product(matrices[1].columns, columns.vector(matrices[1].basis, particular)) == rhs3[term(names)]
    frame2[term(names)] = particular
    known_kernel = [parameters.parse(columns, rows) for rows in even["kernels"]]
    for name, jet in zip(names[8:], known_kernel, strict=True):
        assert not columns.product(matrices[1].columns, columns.vector(matrices[1].basis, jet))
        frame2[term(names, **{name: 1})] = jet
    witnesses = [columns.vector(matrices[1].basis, jet) for jet in known_kernel]
    identity = flint.fmpq_mat([
        [flint.fmpq(v[j].numerator, v[j].denominator) for j in matrices[1].free]
        for v in witnesses
    ])
    assert identity.det()
    frames.append(frame2)

    rhs4 = reduce_rhs(rhs_coefficient(names, frames, 4), normal_conic)
    frame3, equations4 = matrices[2].solve(rhs4)
    row_space = row_span(equations4)
    assert contained_in_span(expression(names, (1, {"d": 1, "f": 2}), (375, {"f": 3})), row_space)
    assert contained_in_span(expression(names, (1, {"e": 2, "z": 1})), row_space)
    # f!=0 would give d=-375*f and then 270000*f^2=0 by the conic.
    assert 375**2+250*375+35625 == 270000
    # Then d=f=0, and e^2*z=0 gives e*z=0 on field-valued points.
    new = [name for name in names if name not in {"d", "f"}] + ["omega", "tau"]
    frames, equations4 = change_variables(names, new, [*frames, frame3], equations4)
    names = new
    equations4.append(expression(names, (1, {"e": 1, "z": 1})))
    for name, jet in zip(("omega", "tau"), matrices[2].kernel, strict=True):
        frames[3][term(names, **{name: 1})] = jet

    basis4 = groebner_with_lift(names, equations4, directory, singular)
    normal4 = divisors(basis4)
    assert not reduce_polynomial(expression(names,
        (1, {"z": 1, "kappa": 1}), (-162, {"z": 3})), normal4)
    frames = reduce_frames(frames, normal4)
    rhs5 = reduce_rhs(rhs_coefficient(names, frames, 5), normal4)
    frame4, equations5 = matrices[3].solve(rhs5)
    row_space = row_span(equations5)

    def require5(*monomials):
        assert contained_in_span(reduce_polynomial(expression(names, *monomials), normal4), row_space)

    require5((1, {"sigma": 2}))
    require5((1, {"kappa": 2}), (-162, {"kappa": 1, "z": 2}))
    require5((-27, {"delta": 1, "kappa": 1}), (4374, {"delta": 1, "z": 2}),
             (11, {"kappa": 2}))
    require5((1, {"e": 3, "eta": 1}))
    # z*(kappa-162*z^2)=kappa*(kappa-162*z^2)=0, so their difference
    # forces (kappa-162*z^2)^2=0. The next equation then gives kappa=z=0.
    new = [name for name in names if name not in {"sigma", "kappa", "z"}] + ["free5_0"]
    frames, equations = change_variables(names, new, [*frames, frame4], equations4 + equations5)
    names = new
    frames[4][term(names, free5_0=1)], = matrices[3].kernel
    i, j = names.index("e"), names.index("eta")
    # e^3*eta=0 implies e*eta=0 over every field.
    equations = [{m: c for m, c in p.items() if not (m[i] and m[j])} for p in equations]
    equations.append(expression(names, (1, {"e": 1, "eta": 1})))
    gamma_relation = ((3, {"gamma": 1}), (10, {"chi": 1}), (9, {"psi": 1}),
                      (-576, {"c": 2}), (3888, {"c": 1, "e": 1}))
    assert contained_in_span(expression(names, *gamma_relation), row_span(equations))
    new = [n for n in names if n != "gamma"]
    gamma = expression(new, (-Fraction(10,3), {"chi": 1}), (-3, {"psi": 1}),
                       (192, {"c": 2}), (-1296, {"c": 1, "e": 1}))
    frames, equations = change_variables(names, new, frames, equations, {"gamma": gamma})
    names = new
    basis5 = groebner_with_lift(names, equations, directory, singular)
    normal5 = divisors(basis5)
    frames = reduce_frames(frames, normal5)
    rhs6 = reduce_rhs(rhs_coefficient(names, frames, 6), normal5)
    frame5, equations6 = matrices[4].solve(rhs6)
    assert not matrices[4].kernel

    basis6 = groebner_with_lift(names, equations + equations6, directory, singular)
    assert len(basis6) == 13
    normal6 = divisors(basis6)
    frames = reduce_frames([*frames, frame5], normal6)
    raw7 = rhs_coefficient(names, frames, 7)
    compare_moyal_at_point(names, frames, 7, raw7)
    rhs7 = reduce_rhs(raw7, normal6)
    _, equations7 = matrices[5].solve(rhs7)
    assert contained_in_span({term(names): Fraction(1)}, row_span(equations7))
    return {
        "h3_initial_residual_rank": len(equations3),
        "h4_frame_constraints": len(basis4),
        "h5_frame_constraints": len(basis5),
        "h6_frame_constraints": len(basis6),
        "h6_parameter_names": names,
        "h6_parameter_basis": [[[list(m), str(c)] for m, c in row] for row in basis6],
        "h7_RHS_parameter_monomials": len(rhs7),
        "h7_new_residual_rank": len(equations7),
        "h7_residual_span_contains_one": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    assert all(PRIME % j for j in range(2, isqrt(PRIME)+1))
    source = columns.load(SOURCE, "canonical_long_complete_cap27")
    images = [{e: Fraction(c) for e, c in p.items()} for p in source.long_images()]
    assert all(sum(e) <= CAP and sum(x*y for x, y in zip(e, (-1, 1, 4, 2))) == w
               for p, w in zip(images, columns.WEIGHTS) for e in p)
    generators = [{term(SPATIAL, **{name: 1}): Fraction(1)} for name in SPATIAL]
    for i, j in columns.PAIRS:
        assert polynomials.bidifferential(images[i], images[j], 1) == polynomials.bidifferential(
            generators[i], generators[j], 1
        )
    assert coefficient_basis(11) == []
    assert (1, (26, 0, 1, 0)) in coefficient_basis(8)
    assert (2, (26, 0, 0, 0)) in coefficient_basis(10)
    assert polynomials.bidifferential({(26, 0, 1, 0): Fraction(1)},
                                     {(26, 0, 0, 0): Fraction(1)}, 1)
    matrices = [LinearElimination(images, m) for m in range(1, 11)]
    counts = [len(L.basis) for L in matrices]
    ranks = [L.rank for L in matrices]
    assert sum(counts) == 2054
    assert ranks == [647, 480, 345, 235, 152, 91, 49, 22, 7, 1]
    print("Exact coefficient counts and ranks:", list(zip(counts, ranks)), flush=True)
    with TemporaryDirectory(prefix="long-cap27-") as work:
        evidence = solve_all(images, matrices, Path(work), args.singular)

    sources = [Path(__file__), SOURCE, *(HERE / name for name in (
        "GRADED_LONG_CAP27_FULL_OBSTRUCTION.md",
        "check_graded_long_cap23.py", "check_graded_long_cap23_full_jet.py",
        "check_graded_long_cap26_full_jet.py", "check_long_graph_curvature.py",
        "graded_long_cap27_even_witness.json",
    ))]
    report = {"schema_version": 1, "status": "passed", "cap": CAP,
              "parity": "unrestricted", "FLINT_python_version": flint.__version__,
              "raw_coefficient_count": sum(counts), "symbol_order_ceiling": 10,
              "commutator_order_ceiling": 19, "component_sizes": [
                  [sum(i == j for i, _ in L.basis) for j in range(4)] for L in matrices],
              "linear_ranks_Q": ranks, **evidence,
              "input_hashes": {p.name: sha256(p.read_bytes()).hexdigest() for p in sources}}
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(report, indent=2) + "\n")
    print("PASS: the exact residual span contains 1 at h^7 with all parity branches retained.")


if __name__ == "__main__":
    main()
