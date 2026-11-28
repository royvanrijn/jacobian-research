#!/usr/bin/env python3
"""Certify the ordinary-balanced six-cell standard sector on rational slabs."""

import argparse
from functools import cache
from hashlib import sha256
from pathlib import Path
from itertools import combinations_with_replacement, permutations
import json
from math import factorial
import platform
import sys

import flint
from flint import acb, arb, arb_mat, ctx, fmpq, fmpq_mat, fmpz_mat, fmpz_poly

import check_tetrahedral_standard_gap as standard_base
import check_tetrahedral_weighted_gap as weighted_base
from check_tetrahedral_standard_gap import hermite, require
from check_tetrahedral_weighted_gap import finite_matrix, planar_tail_coefficients, rational_upper


GAP = fmpq(1, 1000)
SLABS = (
    (fmpq(1, 3), 6), (fmpq(1, 2), 8), (fmpq(3, 5), 10),
    (fmpq(13, 20), 12), (fmpq(7, 10), 14),
)

# Raw coordinates against four mutually orthogonal integer vectors.
# RIGHT gives the coordinates after the adjacent transport (23).
VARIANCES = (12, 2, 6, 12, 2)
LEFT = ((1, 1, 1, 1), (-1, 1, 0, 0), (-1, -1, 2, 0), (-1, -1, -1, 3), (0, 0, 0, 0))
RIGHT = ((-2, 1, 1, 1), (0, 1, 0, 0), (0, -1, 2, 0), (0, -1, -1, 3), (1, 0, 0, 0))



def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first, *rest)


def coefficient_numerators(basis, maximum):
    """Interaction numerators over 6*12**m*m!, planar over 12**m*m!."""
    index = {exponents: i for i, exponents in enumerate(basis)}
    radix = maximum + 1
    packed_index = {
        i + radix*(j + radix*(k + radix*l)): index[tuple(sorted((i, j, k, l)))]
        for degree in range(1, maximum + 1) for i, j, k, l in compositions(degree, 4)
    }

    @cache
    def packed_hermite(linear, degree, variance):
        p = fmpz_poly()
        for order, a in enumerate(hermite(degree, variance)):
            if not a:
                continue
            for i, j, k, l in compositions(order, 4):
                b = linear[0]**i * linear[1]**j * linear[2]**k * linear[3]**l
                if b:
                    coefficient = factorial(order)//(factorial(i)*factorial(j)*factorial(k)*factorial(l))
                    p[i + radix*(j + radix*(k + radix*l))] = a*b*coefficient
        return p

    def project(p, size):
        values = [0]*size
        for position, a in enumerate(p.coeffs()):
            if position and a:
                values[packed_index[position]] += int(a)
        return values

    def sparse(matrix):
        return {(i, j): int(matrix[i, j]) for i in range(matrix.nrows())
                for j in range(matrix.ncols()) if matrix[i, j]}

    interaction = []
    planar = []
    for m in range(1, maximum + 1):
        size = sum(sum(exponents) <= m for exponents in basis)
        left_rows = []
        right_rows = []
        plane_rows = []
        scaled_plane_rows = []
        for degrees in compositions(m, 5):
            if degrees[-1] % 2:
                continue
            x, y = fmpz_poly([1]), fmpz_poly([1])
            divisor = 1
            for l, r, degree, variance in zip(LEFT, RIGHT, degrees, VARIANCES):
                x *= packed_hermite(l, degree, variance)
                y *= packed_hermite(r, degree, variance)
                divisor *= variance**degree * factorial(degree)
            fx = project(x, size)
            if not any(fx):
                continue
            fy = project(y, size)
            factor, remainder = divmod(12**m*factorial(m), divisor)
            require(remainder == 0, f"Nonintegral coefficient factor of degree {m}")
            left_rows.append(fx)
            right_rows.append([factor*(2*a + 4*b) for a, b in zip(fx, fy)])
            if degrees[-1] == 0:
                plane_rows.append(fx)
                scaled_plane_rows.append([factor*a for a in fx])
        kernel = fmpz_mat(left_rows).transpose() * fmpz_mat(right_rows)
        require(kernel == kernel.transpose(), f"Asymmetric Hermite coefficient {m}")
        interaction.append(sparse(kernel))
        planar.append(sparse(fmpz_mat(plane_rows).transpose() * fmpz_mat(scaled_plane_rows)))
    require(interaction[0] == {(0, 0): 48}, "Incorrect degree-one standard coefficient")
    require(planar[0] == {(0, 0): 16}, "Incorrect degree-one planar coefficient")
    return interaction, planar


def bound_weights(a, c, degree):
    sigma_a = (1 - arb(a*a)).sqrt()
    sigma_c = rational_upper((1 - arb(c*c)).sqrt())
    weights = []
    corrections = []
    for m, tail in enumerate(planar_tail_coefficients(c, degree), 1):
        critical = fmpq(m-1, m)
        upper = sigma_c*c**(m-1)
        if a*a >= critical:
            upper = max(upper, rational_upper(sigma_a * arb(a**(m-1))))
        elif c*c >= critical:
            peak = arb(critical).sqrt()**(m-1) / arb(m).sqrt()
            upper = max(upper, rational_upper(peak))
        weights.append(upper/(6*12**m*factorial(m)))
        corrections.append(tail/6**m)
    return weights, corrections


def covariance_key(matrix):
    return tuple(tuple(matrix[i, j] for j in range(matrix.ncols()))
                 for i in range(matrix.nrows()))


@cache
def orthant_probability(covariance):
    dimension = len(covariance)
    require(dimension <= 4, "The orthant formula requires dimension at most four")
    if dimension == 4:
        correlations = {
            (i, j): arb(covariance[i][j]) / arb(covariance[i][i]*covariance[j][j]).sqrt()
            for i in range(4) for j in range(4) if i != j
        }
        def derivative(t, analytic):
            total = acb(0)
            for i in range(4):
                for j in range(i):
                    k, l = [q for q in range(4) if q not in (i, j)]
                    r = t * correlations[i, j]
                    determinant = 1-r**2

                    def conditioned_entry(u, v):
                        a, b = (t*correlations[u, i], t*correlations[u, j])
                        c, d = (t*correlations[v, i], t*correlations[v, j])
                        prior = 1 if u == v else t*correlations[u, v]
                        return prior-(a*c+b*d-r*(a*d+b*c))/determinant

                    cond = (conditioned_entry(k,l)
                            / conditioned_entry(k,k).sqrt(analytic=analytic)
                            / conditioned_entry(l,l).sqrt(analytic=analytic))
                    if analytic and not abs(cond) < 1:
                        return acb(arb("nan"))
                    total += correlations[i, j] / determinant.sqrt(analytic=analytic) * (
                        arb(1)/4 + cond.asin()/(2*arb.pi()))
            return total/(2*arb.pi())

        accuracy = arb(2)**(-125)
        remainder = acb.integral(derivative, 0, 1, abs_tol=accuracy, rel_tol=accuracy)
        require(remainder.real.is_finite() and remainder.imag.contains(0), "Orthant integral failed")
        return arb(1)/16+remainder.real
    coefficients = (
        (arb(covariance[i][j]) / arb(covariance[i][i]*covariance[j][j]).sqrt()).asin()/arb.pi()
        for i in range(dimension) for j in range(i)
    )
    return (1+2*sum(coefficients, arb(0)))/2**dimension


@cache
def conditioned(covariance, j):
    remaining = [i for i in range(len(covariance)) if i != j]
    return tuple(tuple(covariance[k][l] - covariance[k][j]*covariance[j][l]/covariance[j][j]
                       for l in remaining) for k in remaining)


@cache
def face_density(covariance, j):
    return 1 / (2 * arb.pi() * arb(covariance[j][j])).sqrt()


@cache
def raw_moment(covariance, alpha):
    if not any(alpha):
        return orthant_probability(covariance)
    i = next(i for i, exponent in enumerate(alpha) if exponent)
    beta = tuple(x-(i == j) for j, x in enumerate(alpha))
    value = arb(0)
    for j, c in enumerate(covariance[i]):
        if not c:
            continue
        if beta[j]:
            reduced = tuple(x-(j == l) for l, x in enumerate(beta))
            term = beta[j] * raw_moment(covariance, reduced)
        else:
            boundary = raw_moment(conditioned(covariance, j), beta[:j] + beta[j+1:])
            term = face_density(covariance, j) * boundary
        value += arb(c)*term
    return value


@cache
def orbit(exponents):
    return tuple(sorted(set(permutations(exponents))))


def weighted_gram(basis, rho, allowance, covariance):
    remain = 1 - allowance
    require(remain > 0, "Tail and target gap exceed the local lower bound")
    precision = covariance.inv()
    base = covariance_key(covariance)
    p = orthant_probability(base)
    rate = rho*rho/(3*(1-rho*rho))
    if rho == 0:
        terms = [(base, 1/p)]
    else:
        rates = [fmpq(0), rate, fmpq(4, 3)*rate, 4*rate]
        terms = []
        for rr in combinations_with_replacement(rates, 4):
            tilted_precision = fmpq_mat([
                [precision[i, j] + (2*rr[i] if i == j else 0) for j in range(4)]
                for i in range(4)
            ])
            tilted = tilted_precision.inv()
            mass_ratio = arb(tilted.det()/covariance.det()).sqrt()
            weight = arb(6)*len(orbit(rr))*mass_ratio/(2*9**sum(bool(x) for x in rr))
            terms.append((covariance_key(tilted), weight))

    @cache
    def moment(alpha):
        alpha_orbit = orbit(alpha)
        return sum((w*sum((raw_moment(c, v) for v in alpha_orbit), arb(0))/len(alpha_orbit)
                    for c, w in terms), arb(0))

    mass = moment((0, 0, 0, 0))
    require(mass > 0, "Weighted facet mass must be positive")
    if rho == 0:
        require(mass.contains(1), "Unweighted facet measure is not a probability")
    moments = [moment(x) for x in basis]

    # Each basis function is the normalized permutation average of one monomial.
    size = len(basis)
    gram = arb_mat(size, size)
    for i, x in enumerate(basis):
        for j in range(i+1):
            right_orbit = orbit(basis[j])
            product = sum((moment(tuple(sorted(x[k]+y[k] for k in range(4))))
                           for y in right_orbit), arb(0))/len(right_orbit)
            value = (product - moments[i]*moments[j]/mass)/remain
            gram[i, j] = value
            gram[j, i] = value
    raw_moment.cache_clear()
    conditioned.cache_clear()
    orthant_probability.cache_clear()
    face_density.cache_clear()
    return gram


def certify(slabs, start=fmpq(0), ratio_target=fmpq(1)):
    ctx.prec = 512
    maximum = max(n for _, n in slabs)
    basis = sorted({tuple(sorted(x)) for d in range(1, maximum + 1) for x in compositions(d, 4)},
                   key=lambda x: (sum(x), x))
    covariance = fmpq_mat([[fmpq(1, 2)+(i == j) for j in range(4)] for i in range(4)])
    p = orthant_probability(covariance_key(covariance))
    require(arb(747)/5000 < p < arb(1)/6, "Incorrect six-cell facet mass")
    small_sigma = ((fmpq(2),),)
    require(raw_moment(small_sigma, (2,)).contains(1), "Incorrect half-line second moment")
    interaction, planar = coefficient_numerators(basis, maximum)

    left = start
    for right, degree in slabs:
        require(0 <= left <= right < 1, "Invalid correlation cover")
        if ratio_target > 1:
            require(left > 0, "The stationary-gap target requires positive correlation")
            factor = 1 - (ratio_target-1)*(1/left-1)
        else:
            factor = fmpq(1)
        small_basis = [x for x in basis if sum(x) <= degree]
        weights, corrections = bound_weights(left, right, degree)
        tail = arb(right)**degree/p
        gram = weighted_gram(small_basis, left, 1-arb(factor) + tail + arb(GAP), covariance)
        matrix = finite_matrix(interaction[:degree] + planar[:degree], len(small_basis),
                               weights + corrections)
        powered = matrix*gram
        for exponent in (2, 4, 8, 16, 32, 64):
            powered = powered*powered
            trace = powered.trace()
            if 0 < trace < 1:
                break
        require(0 < trace < 1, f"No trace certificate on [{left},{right}]: {trace.str(14)}")
        norm_bound = trace**(arb(1)/exponent)
        normalized_upper = tail + norm_bound*(arb(factor)-tail-arb(GAP))
        z_upper = min(fmpq(117, 100), 1 + fmpq(15, 8)*right**2/(1-right**2))
        ratio_lower = (1-arb(left)*normalized_upper)/(arb(1-left)*arb(z_upper))
        yield {
            "left": str(left), "right": str(right), "degree": degree, "size": len(small_basis),
            "normalized_factor": str(factor), "trace_exponent": exponent, "trace": trace.str(12),
            "norm_bound": norm_bound.str(12), "operator_tail_bound": tail.str(12),
            "normalized_interaction_upper": normalized_upper.str(12),
            "mass_mismatch_upper": str(z_upper), "stationary_ratio_lower": ratio_lower.str(12),
        }
        left = right


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slabs", type=int, default=len(SLABS), choices=range(1, len(SLABS)+1),
                        help="number of consecutive slabs to replay")
    parser.add_argument("--start", type=fmpq, default=fmpq(0))
    parser.add_argument("--stop", type=fmpq)
    parser.add_argument("--degree", type=int)
    parser.add_argument("--ratio-target", type=fmpq, default=fmpq(1))
    args = parser.parse_args()
    if args.stop is None:
        slabs = SLABS[:args.slabs]
    else:
        require(args.degree is not None and args.degree > 0, "Degree is required for a custom slab")
        slabs = ((args.stop, args.degree),)
    rows = []
    for row in certify(slabs, args.start, args.ratio_target):
        rows.append(row)
        print(f"Certified {row['left']} to {row['right']} (trace power {row['trace_exponent']})",
              flush=True, file=sys.stderr)
    sources = tuple(Path(m.__file__) for m in (standard_base, weighted_base))
    hashes = {p.name: sha256(p.read_bytes()).hexdigest()
              for p in (Path(__file__), *sources)}
    print(json.dumps({"status": "pass", "python": platform.python_version(),
                      "python_flint": flint.__version__, "precision_bits": ctx.prec,
                      "epsilon": str(GAP), "ratio_target": str(args.ratio_target),
                      "slabs": rows, "source_sha256": hashes}, indent=2))


if __name__ == "__main__":
    main()
