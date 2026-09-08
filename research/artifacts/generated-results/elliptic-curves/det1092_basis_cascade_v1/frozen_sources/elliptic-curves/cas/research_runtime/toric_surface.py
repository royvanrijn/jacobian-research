"""Replay retained toric Frobenius output for a separated elliptic K3 twist.

The expensive exact CAS computation is retained as driver output, with its
source/build provenance. Replay derives the model, boundary quotient, Weil
condition and all rank factors. It does not rerun controlled reduction.
Supported geometry: squarefree K3 discriminant, smooth infinity and a disjoint
squarefree quadratic or quartic twist. Names of experiments play no role.
"""
from pathlib import Path
import sys

from .sage_surface import separated_twist_fibres, UnsupportedFibres
from .store import digest

TCR_COMMIT = '74cda9e8148cd8e9a3928fc15a558c9a70b67cc1'
PROOF_KIND = 'retained-toric-k3-separated-twist-v1'


def weil_circle(polynomial, prime):
    from sage.all import AA, PolynomialRing, QQ
    ring = PolynomialRing(QQ, 'z'); z = ring.gen(); n = polynomial.degree()
    normalized = ring(polynomial(prime*z)/prime**n)
    sign = normalized[0]
    if sign not in (-1, 1) or normalized != sign*z**n*normalized(1/z):
        raise ArithmeticError('Frobenius functional equation failed')
    core = normalized
    for endpoint in (-1, 1):
        while core.degree() and core(endpoint) == 0:
            core //= z-endpoint
    if core.degree() % 2 or core != z**core.degree()*core(1/z):
        raise ArithmeticError('nonreciprocal Weil core')
    m = core.degree()//2
    traces = PolynomialRing(QQ, 'w'); w = traces.gen(); chebyshev = [traces(2), w]
    for i in range(2, m+1):
        chebyshev.append(w*chebyshev[-1]-chebyshev[-2])
    trace = traces(core[m])+sum(core[m+i]*chebyshev[i] for i in range(1, m+1))
    if core != z**m*trace(z+1/z):
        raise ArithmeticError('Weil trace reconstruction failed')
    roots = trace.roots(AA)
    if sum(e for _, e in roots) != m or any(not -2 <= a <= 2 for a, _ in roots):
        raise ArithmeticError('Frobenius roots fail the Weil circle')


def replay_toric(surface, proof):
    from sage.all import GF, PolynomialRing, QQ, ZZ
    scripts = Path(__file__).resolve().parents[3]/'elkies-k3/scripts'
    sys.path.insert(0, str(scripts))
    from parse_toric_controlled_reduction_output import parse_readfile_output
    if proof.get('proof_kind') != PROOF_KIND or proof.get('surface_key') != surface.key:
        raise ValueError('toric witness belongs to a different surface or proof algorithm')
    provenance = proof['provenance']
    if provenance.get('source_commit') != TCR_COMMIT or provenance.get('nondegenerate_driver_completed') is not True:
        raise ValueError('unknown toric arithmetic backend or incomplete driver')
    for key in ('executable_sha256', 'raw_output_sha256'):
        value = provenance.get(key, '')
        if len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
            raise ValueError('missing toric build/output provenance')
    from hashlib import sha256
    if sha256(proof['raw_output'].encode()).hexdigest() != provenance['raw_output_sha256']:
        raise ArithmeticError('retained toric output changed')
    row = parse_readfile_output(proof['raw_output'])
    p = proof['prime']
    fibre = separated_twist_fibres(surface, p)
    if fibre['geometric_I1_count'] != 24 or fibre['geometric_I0star_count'] not in (2, 4):
        raise UnsupportedFibres('this toric comparison theorem covers quadratic/quartic K3 twists')
    d_degree = fibre['geometric_I0star_count']
    toric_degree = fibre['expected_L_degree']+2*d_degree
    pg = fibre['chi']-1
    if row['prime'] != p or row['hodge_numbers'] != [pg, toric_degree-2*pg, pg]:
        raise ArithmeticError('wrong toric prime or Hodge vector')
    ring = PolynomialRing(GF(p), 't')
    A, B, d = (ring(list(map(QQ, c))) for c in (surface.A, surface.B, surface.d))
    terms = {(0, 3, 0): -GF(p).one()}
    for polynomial, exponents in ((d, (0,2)), (-A, (1,0)), (-B, (0,0))):
        for i, coefficient in enumerate(polynomial):
            if coefficient:
                terms[(i, *exponents)] = coefficient
    monomials = sorted(terms)
    facets = [(0,0,1), (0,1,0), (-1,-4,-(6-d_degree//2)), (0,-2,-3), (1,0,0)]
    if (row['monomials'] != monomials or row['coefficients'] != [int(terms[m]) for m in monomials]
        or row['halfspace_A'] != facets or row['halfspace_b'] != [0,0,12,6,0]):
        raise ArithmeticError('toric output does not represent the exact surface')
    polynomials = PolynomialRing(ZZ, 'T'); T = polynomials.gen()
    full = polynomials(row['frobenius_coefficients'])
    if full.degree() != toric_degree or not full.is_monic():
        raise ArithmeticError('incomplete toric polynomial')
    numerator, denominator = polynomials(1), polynomials(1)
    for branch in fibre['branches']:
        e = branch['place_degree']; denominator *= T**e-p**e
        for degree in branch['residual_factor_degrees']:
            numerator *= T**(e*degree)-p**(e*degree)
    boundary, remainder = numerator.quo_rem(denominator)
    if remainder or boundary.degree() != 2*d_degree:
        raise ArithmeticError('invalid boundary permutation quotient')
    elliptic, remainder = full.quo_rem(boundary)
    if remainder or elliptic.degree() != fibre['expected_L_degree'] or not elliptic.is_monic():
        raise ArithmeticError('wrong elliptic Frobenius factor')
    if list(map(str, elliptic)) != list(map(str, proof['coefficients'])):
        raise ArithmeticError('supplied L-polynomial differs from toric boundary quotient')
    n = elliptic.degree()
    moments = [-elliptic[n-1], elliptic[n-1]**2-2*elliptic[n-2]]
    if list(map(str, moments)) != list(map(str, proof['moments'])):
        raise ArithmeticError('retained moment audit disagrees')
    weil_circle(elliptic, p)
    return True
