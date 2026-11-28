#!/usr/bin/env python3
"""Certify the six-cell standard gap by a counted complementary bound."""

import argparse
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import platform
import sys

import flint
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

import check_six_cell_weighted_gap as six
import check_tetrahedral_standard_gap as standard
import check_tetrahedral_weighted_gap as weighted


FACET_MASS_LOWER = fmpq(747, 5000)
RECIPROCAL_ROW_LOWER = 3 * FACET_MASS_LOWER
SHIFTS = tuple(fmpq(n, 100) for n in range(60, 100))
RADIAL = {(0, 0, 0, 2): 5, (0, 0, 1, 1): -3}
LINEAR = {(0, 0, 0, 1): 1}


def kernel_and_tail(left, right, degree):
    basis = sorted({tuple(sorted(x)) for d in range(1, degree + 1)
                    for x in six.compositions(d, 4)}, key=lambda x: (sum(x), x))
    covariance = fmpq_mat([[fmpq(1, 2) + (i == j) for j in range(4)] for i in range(4)])
    facet_mass = six.orthant_probability(six.covariance_key(covariance))
    standard.require(facet_mass > arb(FACET_MASS_LOWER), "Facet mass lower bound failed")
    tail = weighted.rational_upper(arb(right)**degree / facet_mass)
    standard.require(tail < 1, "The analytic tail alone exceeds the target")

    interaction, planar = six.coefficient_numerators(basis, degree)
    weights, corrections = six.bound_weights(left, right, degree)
    kernel = weighted.finite_matrix(interaction + planar, len(basis), weights + corrections)
    return basis, covariance, kernel, tail


def finite_forms(left, right, degree):
    """Return a signed kernel, a dominating feature Gram and the tail cap."""
    basis, covariance, kernel, tail = kernel_and_tail(left, right, degree)
    # For z=D^{-1} between a and 1, the secant line overbounds z/(1-tail*z).
    # Replace z by its existing upper envelope E and use the ordinary-balanced
    # covariance inequality G_(alpha E-beta) <= alpha G_E - beta G_mu.
    a = RECIPROCAL_ROW_LOWER
    denominator = (1 - tail)*(1 - tail*a)
    alpha = 1 / denominator
    beta = tail*a / denominator
    gram_e = six.weighted_gram(basis, left, 0, covariance)
    gram_mu = six.weighted_gram(basis, fmpq(0), 0, covariance)
    gram = gram_e*arb(alpha) - gram_mu*arb(beta)
    return basis, kernel, gram, tail


def count_shifts(operator):
    """Use even spectral moments, independently of any trial choice."""
    proved = {}
    power = operator
    for exponent in (2, 4, 8, 16):
        power = power*power
        trace = power.trace()
        for shift in SHIFTS:
            if shift not in proved and trace < 2*arb(shift)**exponent:
                proved[shift] = exponent, trace
    return proved


def rayleigh_and_residual(vector, image, gram):
    norm = (vector.transpose()*gram*vector)[0, 0]
    standard.require(norm > 0, "Trial norm is not strictly positive")
    rayleigh = (vector.transpose()*gram*image)[0, 0] / norm
    residual = image - vector*rayleigh
    squared = (residual.transpose()*gram*residual)[0, 0] / norm
    return rayleigh, squared


def feature_column(basis, terms):
    return arb_mat(len(basis), 1, [terms.get(exponent, 0) for exponent in basis])


def trial_upper(basis, operator, gram, shift):
    """Build radial and linear Krylov trials using only directed matrix products."""
    best = None
    for name, terms in (("radial", RADIAL), ("linear", LINEAR)):
        vector = feature_column(basis, terms)
        for steps in range(9):
            image = operator*vector
            rayleigh, residual_squared = rayleigh_and_residual(vector, image, gram)
            if rayleigh > arb(shift):
                enclosure = rayleigh + residual_squared/(rayleigh - arb(shift))
                upper = weighted.rational_upper(enclosure)
                if best is None or upper < best[0]:
                    best = upper, name, steps, rayleigh, residual_squared
            vector = image
    standard.require(best is not None, "No trial proves there is an eigenvalue above the shift")
    return best


def true_complement(left, right, degree):
    """Independently bound the full operator off one radial moment."""
    basis, covariance, kernel, tail = kernel_and_tail(left, right, degree)
    gram = six.weighted_gram(basis, left, 0, covariance)
    # Remove the representer of the ordinary radial moment in L²(E^{-1} dμ).
    # Its orthogonality condition is codimension one in L²(D dμ) as well.
    radial = feature_column(basis, RADIAL)
    column = gram*radial
    norm = (radial.transpose()*column)[0, 0]
    standard.require(norm > 0, "Radial norm is not strictly positive")
    projected_gram = gram - column*column.transpose()/norm
    power = kernel*projected_gram
    for _ in range(4):
        power = power*power
    trace = power.trace()
    standard.require(trace > 0, "Could not enclose the projected even trace")
    upper = weighted.rational_upper(arb(tail) + trace**(arb(1)/16))
    standard.require(upper < 1, "The count did not exclude a second dangerous mode")
    return {
        "degree": degree,
        "feature_count": len(basis),
        "trial": "radial moment complement",
        "analytic_tail_upper": str(tail),
        "trace_exponent": 16,
        "compressed_trace": trace.str(16),
        "true_second_eigenvalue_upper": str(upper),
        "centered_radial_envelope_moment_upper": str(weighted.rational_upper(4*norm/9)),
    }


def gaussian_density(x):
    return (-x*x/2).exp()/(2*arb.pi()).sqrt()


def gaussian_cdf(x):
    return (1 + (x/arb(2).sqrt()).erf())/2


def latent_cdf(w):
    return (1 + w.erf())/2 - (-w*w/2).exp()*gaussian_cdf(w)/arb(2).sqrt()


def maximum_density_product(w):
    """Bound the maximum of 4(a-w) phi(a) Phi(a)^3 over a >= w."""
    def log_derivative(a):
        return 1/(a-w) - a + 3*gaussian_density(a)/gaussian_cdf(a)

    low, high = w, max(w + 2, arb(3))
    standard.require(log_derivative(high) < 0, "Maximum bracket did not close")
    for _ in range(24):
        mid = (low + high)/2
        derivative = log_derivative(mid)
        if derivative > 0:
            low = mid
        elif derivative < 0:
            high = mid
        else:
            break  # the previous bracket still encloses the maximum
    closest_zero = low if low > 0 else high if high < 0 else arb(0)
    return 4*(high-w)*gaussian_density(closest_zero)*gaussian_cdf(high)**3


def radial_derivative_bound():
    # The conditional maximum decreases in the latent coordinate w.
    # The W law is that of (S+N)/sqrt(2), with S Gaussian and N Rayleigh.
    end, per_unit = 6, 64
    endpoints = [arb(k)/per_unit for k in range(-end*per_unit, end*per_unit+1)]
    masses = [latent_cdf(upper) - latent_cdf(lower)
              for lower, upper in zip(endpoints, endpoints[1:])]
    standard.require(all(mass > 0 for mass in masses), "A latent bin mass is not positive")
    body = sum((mass*maximum_density_product(w)
                for mass, w in zip(masses, endpoints)), arb(0))
    # For w >= 0 the integrand is at most 4 phi(0). For w <= 0 it
    # is at most 1+2|w|, and W is bounded below by S/sqrt(2).
    positive_tail = 4*gaussian_density(arb(0))*(1 - latent_cdf(arb(end)))
    negative_cutoff = -arb(end)*arb(2).sqrt()
    negative_tail = gaussian_cdf(negative_cutoff) + arb(2).sqrt()*gaussian_density(negative_cutoff)
    tail = positive_tail + negative_tail
    upper = weighted.rational_upper((body + tail)/arb(RECIPROCAL_ROW_LOWER))
    return {
        "latent_bin_count": len(masses),
        "latent_bin_width": str(fmpq(1, per_unit)),
        "latent_tail_contribution": tail.str(16),
        "radial_row_derivative_upper": str(upper),
    }


def envelope_mass_upper(rho):
    """Enclose the mass of E using the same four-dimensional tilts as its Gram."""
    covariance = fmpq_mat([[fmpq(1, 2) + (i == j) for j in range(4)] for i in range(4)])
    precision = covariance.inv()
    rate = rho*rho/(3*(1-rho*rho))
    rates = (fmpq(0), rate, fmpq(4, 3)*rate, 4*rate)
    mass = arb(0)
    for diagonal in combinations_with_replacement(rates, 4):
        tilted_precision = fmpq_mat([
            [precision[i, j] + (2*diagonal[i] if i == j else 0) for j in range(4)]
            for i in range(4)
        ])
        tilt = tilted_precision.inv()
        density_ratio = arb(tilt.det()/covariance.det()).sqrt()
        multiplicity = len(six.orbit(diagonal))
        weight = 3*multiplicity*density_ratio/9**sum(bool(rate_i) for rate_i in diagonal)
        mass += weight*six.orthant_probability(six.covariance_key(tilt))
    return weighted.rational_upper(mass)


def radial_certify(left, right, degree):
    """Use one true radial direction and an independently certified complement."""
    standard.require(0 < left <= right < 1, "Invalid positive correlation endpoints")
    standard.require(degree >= 2, "The radial trial requires degree two")
    ctx.prec = 512
    complement = true_complement(left, right, degree)
    complement_cap = fmpq(complement["true_second_eigenvalue_upper"])
    radial_moment = fmpq(complement["centered_radial_envelope_moment_upper"])
    mass = envelope_mass_upper(left)
    derivative = radial_derivative_bound()
    derivative_cap = fmpq(derivative["radial_row_derivative_upper"])
    error_squared = ((1-left*left)/left)**2 * derivative_cap**2*mass*radial_moment/256
    error = arb(error_squared).sqrt()
    angle = arb((1-RECIPROCAL_ROW_LOWER)/(1+RECIPROCAL_ROW_LOWER))

    # The complementary space is orthogonal in dμ, not in D dμ.
    # The angle term handles that cross product in the true block comparison.
    for n in range(1, 1000):
        candidate = fmpq(n, 1000)
        if candidate <= max(right, complement_cap):
            continue
        smallest = arb(candidate-right)
        if not smallest > error:
            continue
        residual_term = 2*arb(error_squared)/(smallest + (smallest**2-arb(error_squared)).sqrt())
        coupling = angle*(arb(candidate-left) + error).sqrt() + residual_term.sqrt()
        margin = arb(candidate-complement_cap) - coupling**2
        if margin > 0:
            break
    else:
        raise ArithmeticError("The one-mode block bound did not prove the gap")
    return {
        "left": str(left),
        "right": str(right),
        "independent_true_complement": complement,
        "radial_derivative": derivative,
        "row_envelope_mass_upper": str(mass),
        "squared_radial_residual_upper": str(error_squared),
        "schur_margin": margin.str(16),
        "normalized_interaction_upper": str(candidate),
        "standard_gap_factor_lower": str(1 - candidate),
        "stationary_ratio_lower": str((1-left*candidate)/((1-left)*fmpq(117, 100))),
    }


def certify(left, right, degree):
    standard.require(0 < left <= right < 1, "Invalid positive correlation endpoints")
    standard.require(degree >= 2, "The radial seed requires degree two")
    ctx.prec = 512
    basis, kernel, gram, tail = finite_forms(left, right, degree)
    print(f"Enclosed both facet Grams at degree {degree} ({len(basis)} monomials)",
          file=sys.stderr, flush=True)
    operator = kernel*gram
    counts = count_shifts(operator)
    standard.require(counts, "No independent count excludes a second eigenvalue")
    shift = min(counts)
    exponent, trace = counts[shift]
    bound, seed, steps, rayleigh, residual_squared = trial_upper(basis, operator, gram, shift)
    print(f"Complementary upper bound {arb(bound).str(12)} at shift {shift}",
          file=sys.stderr, flush=True)
    standard.require(bound < 1, "The counted upper bound did not prove the gap")
    interaction_upper = tail + (1 - tail)*bound
    return {
        "left": str(left),
        "right": str(right),
        "degree": degree,
        "feature_count": len(basis),
        "row_inverse_lower": str(RECIPROCAL_ROW_LOWER),
        "analytic_tail_upper": str(tail),
        "shift": str(shift),
        "count_exponent": exponent,
        "count_trace": trace.str(16),
        "count_threshold": str(2*shift**exponent),
        "trial_seed": seed,
        "krylov_steps": steps,
        "rayleigh": rayleigh.str(16),
        "squared_residual": residual_squared.str(16),
        "complementary_upper": str(bound),
        "normalized_interaction_upper": str(interaction_upper),
        "standard_gap_factor_lower": str(1 - interaction_upper),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=fmpq, default=fmpq(3, 4))
    parser.add_argument("--stop", type=fmpq, default=fmpq(3, 4))
    parser.add_argument("--degree", type=int)
    parser.add_argument("--radial", action="store_true",
                        help="use one true radial trial and an independent Hermite complement")
    args = parser.parse_args()
    if args.radial:
        result = radial_certify(args.start, args.stop, 10 if args.degree is None else args.degree)
    else:
        result = certify(args.start, args.stop, 13 if args.degree is None else args.degree)
    sources = (Path(__file__), *(Path(module.__file__) for module in (six, standard, weighted)))
    print(json.dumps({
        "status": "pass",
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "precision_bits": ctx.prec,
        "mode": "radial" if args.radial else "complementary",
        "certificate": result,
        "source_sha256": {path.name: sha256(path.read_bytes()).hexdigest() for path in sources},
    }, indent=2))


if __name__ == "__main__":
    main()
