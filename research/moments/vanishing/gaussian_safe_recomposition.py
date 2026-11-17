"""Exact original-marker recomposition maps at the safe cubic shear P*.

The general safe-recomposition classification is a written theorem. This
QQ utility computes its filtered Gaussian/sphere benchmark and finite jets
without expanding high powers of the recomposed spatial phase.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import factorial

import sympy as sp

from moment_jet_tail_certificates import Z, W, Y, T, PSTAR, linear_null_period, truncate

U = sp.Symbol("U")


def recomposition_data(polynomial):
    value = sp.Poly(sp.expand(polynomial), U, domain=sp.QQ)
    if value.is_zero or value.degree() < 1:
        raise ValueError("Use a nonconstant nonzero recomposition polynomial")
    if value.nth(0) != 0:
        raise ValueError("The polynomial-period backend requires f(0)=0")
    order = min(power[0] for power, coefficient in value.terms() if coefficient)
    return value, order


def inverse_shear(marker):
    return sp.expand(sp.sympify(marker).subs(
        {Z: Z, W: W-Z*Y-Z**3/2, Y: Y+Z**2}, simultaneous=True))


def recompose_polynomial_period(period, polynomial, source_time=T, target_time=T):
    """R_f(g): exact map from a polynomial P-period to an f(P)-period.

    Coefficients of g may include external parameters. A zero constant
    term of f is necessary for a polynomial output; nonzero constants
    would introduce an exponential factor instead.
    """
    f, order = recomposition_data(polynomial)
    g = sp.Poly(sp.expand(period), source_time)
    if g.is_zero:
        return sp.Integer(0)
    limit = g.degree()
    # Terms above the source-period degree can never affect a needed
    # coefficient. Truncate before every multiplication, even for high deg f.
    low_f = sp.Poly.from_dict({powers: coefficient for powers, coefficient in f.terms()
                              if powers[0] <= limit}, U, domain=sp.QQ)
    power = sp.Poly(1, U, domain=sp.QQ)
    answer = sp.Integer(0)
    for m in range(limit//order+1):
        if m:
            power = sp.Poly.from_dict({powers: coefficient
                                      for powers, coefficient in (power*low_f).terms()
                                      if powers[0] <= limit}, U, domain=sp.QQ)
        coefficient = sum(factorial(j)*g.nth(j)*power.nth(j)
                          for j in range(limit+1))
        answer += coefficient*target_time**m/factorial(m)
    return sp.expand(answer)


@lru_cache(maxsize=256)
def _original_period(marker, polynomial, functional):
    if functional not in ("gaussian", "sphere_3"):
        raise ValueError("Unknown normalized functional")
    sp.Poly(marker, Z, W, Y, domain=sp.QQ)
    g = linear_null_period(inverse_shear(marker), functional)
    return recompose_polynomial_period(g, polynomial)


def original_period(marker, polynomial, functional="gaussian"):
    """Phi_Q^{f(P*)} for an ORIGINAL fixed rational polynomial marker Q."""
    return _original_period(sp.expand(marker), sp.expand(polynomial), functional)


def marker_for_weight(weight):
    if type(weight) is not int or weight < 0:
        raise ValueError("The original shear weight must be a nonnegative integer")
    blocks, remainder = divmod(weight, 3)
    return W**blocks*(1 if remainder == 0 else Z if remainder == 1 else Y)


def filtered_basis(polynomial, degree, functional="gaussian"):
    """A triangular ORIGINAL marker basis of the exact degree<=D image."""
    if type(degree) is not int or degree < 0:
        raise ValueError("Original spatial degree must be nonnegative")
    _, order = recomposition_data(polynomial)
    answer = []
    for m in range(3*degree//order+1):
        weight = order*m
        marker = marker_for_weight(weight)
        period = original_period(marker, polynomial, functional)
        answer.append(dict(time_degree=m, marker=marker,
                           marker_degree=sp.Poly(marker, Z, W, Y).total_degree(),
                           period=period, leading_coefficient=sp.expand(period).coeff(T, m)))
    return answer


def marker_lift(target_period, polynomial, functional="gaussian"):
    """Lift a target polynomial with minimum possible original marker degree."""
    target = sp.Poly(sp.expand(target_period), T, domain=sp.QQ)
    if target.is_zero:
        return sp.Integer(0)
    _, order = recomposition_data(polynomial)
    degree = (order*target.degree()+2)//3
    basis = filtered_basis(polynomial, degree, functional)
    remainder, marker = target.as_expr(), sp.Integer(0)
    for entry in reversed(basis):
        m = entry["time_degree"]
        coefficient = sp.expand(remainder).coeff(T, m)/entry["leading_coefficient"]
        if coefficient:
            marker += coefficient*entry["marker"]
            remainder = sp.expand(remainder-coefficient*entry["period"])
    if remainder != 0:
        raise ArithmeticError("The triangular original-marker lift failed")
    return sp.expand(marker)


@lru_cache(maxsize=32)
def _moment_jet_items(perturbation, parameters, order, polynomial, functional):
    recomposition_data(polynomial)
    if type(order) is not int or order < 0:
        raise ValueError("Jet order must be nonnegative")
    if functional not in ("gaussian", "sphere_3"):
        raise ValueError("Unknown normalized functional")
    value = sp.Poly(perturbation, *parameters, Z, W, Y, domain=sp.QQ)
    if any(not sum(powers[:len(parameters)]) for powers, coefficient in value.terms() if coefficient):
        raise ValueError("Perturbations must vanish at the parameter origin")
    transformed = inverse_shear(perturbation)
    power, exponential = sp.Integer(1), sp.Integer(1)
    for j in range(1, order+1):
        power = truncate(sp.expand(power*transformed), parameters, order)
        if power == 0:
            break
        base_period = linear_null_period(power, functional)
        recomposed = recompose_polynomial_period(base_period, polynomial)
        exponential += T**j*recomposed/factorial(j)
    return tuple((m, sp.expand(coefficient*factorial(m)))
                 for (m,), coefficient in sp.Poly(sp.expand(exponential), T).terms()
                 if m > 0 and coefficient)


def moment_jets(perturbation, parameters, order, polynomial, functional="gaussian"):
    """Actual moment jets at f(P*)+H, with H of arbitrary bounded spatial degree."""
    parameters = tuple(parameters)
    reserved = {variable.name for variable in (Z, W, Y, T, U)}
    if (not parameters or any(not isinstance(value, sp.Symbol) for value in parameters)
            or len({value.name for value in parameters}) != len(parameters)
            or any(value.name in reserved for value in parameters)):
        raise ValueError("Use distinct coefficient symbols separate from Z,W,Y,t,U")
    return dict(_moment_jet_items(sp.expand(perturbation), parameters, order,
                                  sp.expand(polynomial), functional))


def jet_bound(polynomial, perturbation_degree, order):
    """Sharp uniform mixed-cutoff jet budget floor((3d/e+1)r)."""
    _, e = recomposition_data(polynomial)
    if (type(perturbation_degree) is not int or perturbation_degree < 0
            or type(order) is not int or order < 0):
        raise ValueError("Degrees and jet orders must be nonnegative integers")
    return ((3*perturbation_degree+e)*order)//e


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", nargs="+", required=True,
                        help="Rational coefficients of f, starting with its constant term")
    parser.add_argument("--degree", type=int, required=True, help="Original marker degree bound")
    parser.add_argument("--functional", choices=["gaussian", "sphere_3"], default="gaussian")
    args = parser.parse_args()
    polynomial = sum(sp.Rational(value)*U**i for i, value in enumerate(args.coefficients))
    f, order = recomposition_data(polynomial)
    basis = filtered_basis(polynomial, args.degree, args.functional)
    print(json.dumps(dict(polynomial=str(f.as_expr()), order_at_zero=order,
                          original_marker_degree=args.degree,
                          exact_max_time_degree=3*args.degree//order,
                          basis=[{key: str(value) if isinstance(value, sp.Basic) else value
                                  for key, value in entry.items()} for entry in basis]), indent=2))


if __name__ == "__main__":
    main()
