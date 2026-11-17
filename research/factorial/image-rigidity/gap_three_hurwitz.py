"""Exact generalized-Euclidean certificates for the coupled gap-three family.

All coefficients use Q[beta] or Q[n,beta]. This module makes no all-index
positivity assertion. The fixed-index bound belongs to the replay record.
"""
from functools import lru_cache
from itertools import permutations
import sympy as s
from gap_three_wronskian import centered_coefficients


def euclidean_rows(n, beta):
    """HKK step-three Euclidean rows, retaining the full polynomial identity.

    The input degree is N=2*n; f_i has degree N-i. Sparse dictionaries
    avoid introducing or cancelling the auxiliary z variable. A zero or
    wrong-degree leading coefficient raises an error rather than passing.
    """
    beta = s.sympify(beta)
    N = 2*n
    c = centered_coefficients(n, beta)
    rows = [{N-k: c[N-k] for k in range(r, N+1, 3)} for r in range(3)]
    for i in range(N-2):
        h0 = rows[i].get(N-i, s.Integer(0))
        h1 = rows[i+1].get(N-i-1, s.Integer(0))
        if h0 == 0 or h1 == 0:
            raise ValueError('degenerate Euclidean row')
        row = {k: s.cancel(rows[i].get(k, 0)-h0*rows[i+1].get(k-1, 0)/h1)
               for k in range(N-i-3, -1, -3)}
        if row.get(N-i-3, s.Integer(0)) == 0:
            raise ValueError('unexpected Euclidean degree drop')
        rows.append(row)
    return rows


def leading_coefficients(n, beta):
    return [s.cancel(row[2*n-i]) for i, row in enumerate(euclidean_rows(n, beta))]


def positive_fraction_certificate(value, variable):
    """Sufficient exact certificate of positivity on variable >= 0."""
    num, den = s.fraction(s.cancel(value))
    P, Q = s.Poly(num, variable), s.Poly(den, variable)
    positive = (P.TC() > 0 and Q.TC() > 0 and
                all(c >= 0 for c in P.all_coeffs()+Q.all_coeffs()))
    return {'positive_on_nonnegative_axis': bool(positive),
            'numerator_degree': int(P.degree()), 'denominator_degree': int(Q.degree()),
            'numerator_constant': str(P.TC()), 'denominator_constant': str(Q.TC()),
            'numerator_min_coefficient': str(min(P.all_coeffs())),
            'denominator_min_coefficient': str(min(Q.all_coeffs()))}


@lru_cache(maxsize=None)
def uniform_top_coefficients(order):
    """Return C_j(n,beta)=[z^(2n-j)]F_n, polynomially continued in n.

    A polynomial finite-difference identity and its initial value prove
    each formula for every integer n>=0. No interpolation is used.
    """
    n, beta, m = s.symbols('n beta m')
    C = [s.Integer(1)]
    for j in range(1, order+1):
        rhs = (4*m+beta+2)*C[j-1].subs({n:m, beta:beta+1}, simultaneous=True)
        if j >= 2:
            rhs += (4*m+beta+2)*(2*m+beta-1+j)*C[j-2].subs(
                {n:m, beta:beta+1}, simultaneous=True)
        C.append(s.factor(s.summation(s.expand(rhs), (m, 0, n-1))))
    return n, beta, C


def uniform_pivot_numerators():
    """Return unreduced, identity-friendly numerators for h3,...,h6."""
    n, beta, C = uniform_top_coefficients(8)
    D3 = s.expand(C[1]*C[3]-C[4])
    D4 = s.expand(C[2]*C[4]-C[1]*C[5])
    D5 = s.expand(C[5]*D3-C[2]*(C[1]*C[6]-C[7]))
    D6 = s.expand(D4*(C[1]*C[6]-C[7])-D3*(C[2]*C[7]-C[1]*C[8]))
    return n, beta, C, {3:D3, 4:D4, 5:D5, 6:D6}


def check_euclidean_rows(rows, input_coefficients, n):
    """Replay identities in a native exact fraction field, independently of generation."""
    N = 2*n
    if len(rows) != N+1:
        return False
    symbols = sorted(set().union(*(s.sympify(c).free_symbols for c in input_coefficients)), key=str)
    domain = s.QQ.frac_field(*symbols) if symbols else s.QQ
    exact = [{k:domain.from_sympy(s.sympify(v)) for k,v in row.items()} for row in rows]
    source = [domain.from_sympy(s.sympify(v)) for v in input_coefficients]
    zero = domain.zero
    for i in range(3):
        if any(exact[i].get(N-k,zero) != source[N-k] for k in range(i,N+1,3)):
            return False
    for i,row in enumerate(exact):
        if not row or max(row)!=N-i or row[N-i]==zero:
            return False
        if any(k<0 or (N-i-k)%3 for k in row):
            return False
    for i in range(N-2):
        c = exact[i][N-i]/exact[i+1][N-i-1]
        for k in set(exact[i]) | {k+1 for k in exact[i+1]} | set(exact[i+3]):
            if exact[i].get(k,zero)-c*exact[i+1].get(k-1,zero)-exact[i+3].get(k,zero) != zero:
                return False
    return True


def special_minor(index, C):
    """HKK special minor Delta_index from reversed coefficients C_j."""
    size = (index+1)//2
    first_row = 2 if index % 2 else 1
    def at(i,j):
        k = 3*(j+1)-(first_row+i)
        return C[k] if 0<=k<len(C) else s.Integer(0)
    result = 0
    for permutation in permutations(range(size)):
        inversions = sum(permutation[i]>permutation[j]
                         for i in range(size) for j in range(i+1,size))
        result += (-1)**inversions*s.prod(at(i,permutation[i]) for i in range(size))
    return s.expand(result)
