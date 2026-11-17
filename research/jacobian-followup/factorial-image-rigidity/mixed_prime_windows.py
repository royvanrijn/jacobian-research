"""Exact mixed-prime Newton data and square-divisor degree certificates.

The all-parameter lifting and indexed-family conclusions use
MIXED_PRIME_RESIDUAL_WINDOWS.tex. A failed residual or degree-intersection
test is UNKNOWN, not a counterexample to the Factorial conjectures.
"""
from math import factorial, comb
from fractions import Fraction
import sympy as s
import prime_power_windows as p


def _power_mod(base, exponent, modulus):
    result = s.Poly(1, *modulus.gens, domain=modulus.domain)
    base = base.rem(modulus)
    while exponent:
        if exponent & 1:
            result = (result*base).rem(modulus)
        base = (base*base).rem(modulus)
        exponent >>= 1
    return result


def _irreducible_verified(factor, prime):
    """Rabin's finite-field criterion, separate from the factorization call."""
    degree = factor.degree()
    if degree < 1:
        return False
    T = s.Poly(factor.gens[0], *factor.gens, domain=factor.domain).rem(factor)
    powers = [T]
    for j in range(degree):
        powers.append(_power_mod(powers[-1], prime, factor))
    return powers[-1] == T and all(s.gcd(factor, powers[degree//ell]-T).degree() == 0
                                   for ell in s.factorint(degree))


def square_divisor_degrees(rows, prime):
    """Necessary degrees, from residual multiplicities at every edge."""
    T = s.Symbol('T')
    possible = {0}
    branches = []
    for edge in rows:
        residual = s.Poly(sum(c*T**j for j, c in enumerate(edge['residual_coefficients'])), T, modulus=prime)
        unit, factors = residual.factor_list()
        reconstruction = s.Poly(unit, T, modulus=prime)
        for factor, multiplicity in factors:
            assert _irreducible_verified(factor, prime), 'Residual factor irreducibility failed'
            reconstruction *= factor**multiplicity
            step = edge['e']*factor.degree()
            capacity = multiplicity//2
            possible = {value+step*k for value in possible for k in range(capacity+1)}
            branches.append({'edge_left': edge['left'][0], 'e': edge['e'],
                             'residual_factor': [int(factor.nth(j)) % prime for j in range(factor.degree()+1)],
                             'residual_degree': factor.degree(), 'multiplicity': multiplicity,
                             'factor_irreducibility_verified': True,
                             'square_divisor_degree_step': step, 'max_copies_in_square_divisor': capacity})
        assert reconstruction == residual
    return {'possible_degrees': sorted(possible), 'branches': branches}


def mixed_window_certificate(terminal_index, degree, primes=(2, 3), shape=1, use_real_parity=True):
    """Combine p-adic budgets with the proved real-place restriction for moments."""
    if not isinstance(use_real_parity, bool):
        raise ValueError("use_real_parity must be a boolean")
    n, degree = p.integer(terminal_index, 2), p.integer(degree, 1)
    shape = s.sympify(shape)
    if shape.has(s.Float) or shape.is_Rational is not True or shape <= 0:
        raise ValueError('an exact positive rational Gamma shape is required')
    coefficients = [Fraction(comb(n, k)*s.rf(shape, degree*k)) for k in range(n+1)]
    primes = list(dict.fromkeys(p.integer(prime, 2) for prime in primes))
    if not primes or any(not s.isprime(prime) for prime in primes):
        raise ValueError('at least one prime is required')
    intersection = set(range(n//2+1))
    profiles = []
    for prime in primes:
        edges = p.newton_residuals(coefficients, prime)
        profile = square_divisor_degrees(edges, prime)
        intersection &= set(profile['possible_degrees'])
        profiles.append({'prime': prime, 'edges': edges, **profile})
    local_intersection = sorted(intersection)
    if use_real_parity:
        intersection = {value for value in intersection if value % 2 == 0}
    accepted = intersection == {0}
    return {'status': 'adjacent_window_certified' if accepted else 'unresolved_by_supplied_local_data',
            'terminal_index': n, 'moment_indices': [n-1, n], 'exponent': degree,
            'gamma_shape': str(shape),
            'phase': 'A+B*Y^exponent with arbitrary complex A,B, nonzero',
            'profiles': profiles, 'local_possible_square_divisor_degrees': local_intersection,
            'real_place_parity_used': use_real_parity,
            'possible_square_divisor_degrees': sorted(intersection),
            'real_root_structure': 'no real roots for even n; exactly one negative simple real root for odd n',
            'real_place_theorem': 'FIRREALROOTPARITY20261116',
            'all_starting_indices_claimed': False, 'counterexample_claimed': False,
            'ordinary_pair_profile': ({'complete_intersection_degrees': [n-1, n],
                'length': n*(n-1), 'highest_nonzero_degree': 2*n-3,
                'actual_coefficient_ideal_nilpotence_bound': 2*n-2,
                'universal_bound_is_sharp': True} if accepted else None),
            'theorem': 'FIRMULTIPRIMEWINDOW20261116' if accepted else None}


def lifting_certificate(n, cofactor, prime, insertion_power, shape=1):
    """Literal check of the all-a Newton-data lifting identity for supplied data."""
    n, cofactor = p.integer(n, 1), p.integer(cofactor, 1)
    prime, insertion_power, shape = p.integer(prime, 2), p.integer(insertion_power, 1), p.integer(shape, 1)
    if not s.isprime(prime):
        raise ValueError('prime required')
    inserted = prime**insertion_power
    if shape > inserted:
        raise ValueError('shape must be <=p^a for this lifting theorem')
    degree = cofactor*inserted
    shift = cofactor*(inserted-1)//(prime-1)
    base = [comb(n, k)*factorial(cofactor*k) for k in range(n+1)]
    lifted = [Fraction(comb(n, k)*int(s.rf(shape, degree*k)), prime**(shift*k)) for k in range(n+1)]
    base_rows, lifted_rows = p.newton_residuals(base, prime), p.newton_residuals(lifted, prime)
    assert len(base_rows) == len(lifted_rows)
    gamma = (-1)**shift % prime
    for source, target in zip(base_rows, lifted_rows):
        assert all(source[key] == target[key] for key in ('left', 'right', 'e', 'h'))
        expected = [value*pow(gamma, source['e']*j, prime) % prime
                    for j, value in enumerate(source['residual_coefficients'])]
        assert expected == target['residual_coefficients']
    return {'status': 'newton_data_lifting_identity_verified', 'n': n, 'cofactor': cofactor,
            'prime': prime, 'insertion_power': insertion_power, 'exponent': degree,
            'gamma_shape': shape, 'parameter_scale_denominator_exponent': shift,
            'base_edges': base_rows, 'lifted_edges': lifted_rows,
            'square_divisor_degree_profile_preserved': True,
            'theorem': 'FIRNEWTONLIFT20261116'}


def indexed_family_certificate(cofactor, index_multiplier, prime):
    """A finite kernel certifies ALL a,r>=1, not just a list of examples."""
    m, k, prime = p.integer(cofactor, 1), p.integer(index_multiplier, 1), p.integer(prime, 2)
    if not s.isprime(prime) or m*k >= prime:
        raise ValueError('prime p>m*k required')
    T = s.Symbol('T')
    residual = s.Poly(sum(comb(k, j)*(factorial(m*j) % prime)*T**j for j in range(k+1)), T, modulus=prime)
    gcd = s.gcd(residual, residual.diff())
    accepted = gcd.degree() == 0
    return {'status': 'unbounded_indexed_family_certified' if accepted else 'residual_kernel_not_squarefree_no_decision',
            'cofactor': m, 'index_multiplier': k, 'prime': prime,
            'exponents': 'd=m*p^a for EVERY a>=1',
            'gamma_shapes': 'EVERY integer1<=alpha<=p^a',
            'terminal_indices': 'N=k*p^r for EVERY r>=1',
            'moment_indices': ['k*p^r-1', 'k*p^r'], 'phase': 'A+B*Y^d, arbitrary complex coefficients',
            'ordinary_pair_profile': ({'length': 'N*(N-1)',
                'actual_coefficient_ideal_nilpotence_bound': '2*N-2',
                'universal_bound_is_sharp': True} if accepted else None),
            'kernel_residual_coefficients': [int(residual.nth(j)) % prime for j in range(k+1)],
            'kernel_gcd_coefficients': [int(gcd.nth(j)) % prime for j in range(gcd.degree()+1)],
            'all_starting_indices_claimed': False, 'counterexample_claimed': False,
            'theorem': 'FIRINDEXEDMIXEDSFC20261116' if accepted else None}
