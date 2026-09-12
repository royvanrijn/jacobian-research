"""Deterministic equation-only skew and exact root-progression norm sieving."""
from math import gcd
import two_class_relation_core as core


def rectangle(coefficients, area_exponent=18):
    """Minimize a rigorous triangle bound, over fixed-area dyadic rectangles."""
    choices = []
    for a_bits in range(area_exponent+1):
        A, B = 2**a_bits, 2**(area_exponent-a_bits)
        bound = sum(abs(c)*A**i*B**(3-i) for i, c in enumerate(coefficients))
        choices.append((bound, A, B))
    bound, A, B = min(choices)
    return {'A': A, 'B': B, 'norm_triangle_upper_bound': str(bound),
            'lattice_cells': (2*A+1)*B, 'area_exponent': area_exponent}


def progressions(c, p, roots_a, roots_b, axis, fixed):
    """None means the entire line vanishes mod p; otherwise distinct residues."""
    if fixed % p:
        roots = roots_a if axis == 'a' else roots_b
        return sorted({fixed*r % p for r in roots})
    leading = c[3] if axis == 'a' else -c[0]
    return [0] if leading % p else None


def sieve_lines(c, box, roots):
    """Visit every cell; divide all supported powers, with no log-score cutoff.

    roots contains (p, roots of F(x,1), roots of F(1,x)), for F(a,b)
    equal to the norm numerator c3*a^3-c2*a^2*b+c1*a*b^2-c0*b^3.
    This is a small exact prototype, not optimized special-q / large-prime NFS.
    """
    A, B = box['A'], box['B']
    axis = 'a' if A >= B else 'b'
    fixed_values = range(1, B+1) if axis == 'a' else range(-A, A+1)
    lower, upper = (-A, A) if axis == 'a' else (1, B)
    for fixed in fixed_values:
        pairs = [(x, fixed) if axis == 'a' else (fixed, x) for x in range(lower, upper+1)]
        norms = [core.form_value(c, a, b) for a, b in pairs]
        if any(n == 0 for n in norms):
            raise ArithmeticError('unexpected zero norm in irreducible cubic')
        residuals = list(map(abs, norms))
        for p, ra, rb in roots:
            residues = progressions(c, p, ra, rb, axis, fixed)
            if residues is None:
                indices = range(len(norms))
            else:
                indices = (i for r in residues for i in range((r-lower) % p, len(norms), p))
            for i in indices:
                assert residuals[i] % p == 0
                while residuals[i] % p == 0:
                    residuals[i] //= p
        yield [(a, b, n, r) for (a, b), n, r in zip(pairs, norms, residuals) if gcd(a, b) == 1]
