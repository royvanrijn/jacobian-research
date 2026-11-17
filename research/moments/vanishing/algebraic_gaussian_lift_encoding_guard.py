"""Rational input guard for the retained branch-encoding constructor.

The sealed finite-transcription helper labels rank by polynomial degree.
That label is valid only for an irreducible polynomial over QQ(t). Use
this entry point whenever producing a new rank-labelled encoding model.
The normalization proof itself works for nonminimal input too; its rank
is then the degree of the selected minimal factor, not the input degree.
"""
import sympy as s
import verify_algebraic_gaussian_lifts as V


def encode_checked(F,a,k,v):
    F=s.expand(F)
    if F.free_symbols-{V.t,V.x}:
        raise ValueError('guard supports exact rational polynomial inputs only')
    polynomial=s.Poly(F,V.x,domain=s.QQ.frac_field(V.t))
    if polynomial.degree()<1 or not polynomial.is_irreducible:
        raise ValueError('rank-labelled encoding needs an irreducible input over QQ(t)')
    return V.encode(F,a,k,v)
