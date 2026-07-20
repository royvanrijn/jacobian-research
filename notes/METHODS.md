# Additional methods and the 2D rational-pair search

## Exact verification methods

- Expand the determinant over `Z[x,y,z]`; equality to `-2` is coefficientwise.
- Repeat modulo several good primes. This is not a characteristic-zero proof by
  itself unless combined with coefficient-height bounds, but it catches CAS and
  transcription errors cheaply.
- Verify the cubic identities instead of expanding the full inverse. They give a
  compact structural certificate.
- Use elimination ideals for generic fiber degree and projective homogenization
  for solutions at infinity.
- Sample rational targets away from the discriminant and reconstruct every root
  of the cubic, then substitute exactly.

## Rational factor search in dimension two

We seek rational maps `A: C^2 -->> C^2` and `B: C^2 -->> C^2` such that

1. `det DA * (det DB) o A` is a nonzero constant (reciprocal factors);
2. every denominator in `B o A` cancels exactly;
3. the composition extends polynomially across the pole divisor of `A` or `B`;
4. two distinct local/valuation branches at that divisor have the same finite
   image; and
5. ideally, the polynomial composition is not injective.

Condition 5 plus condition 1 would be an actual 2D counterexample, so no current
search should be described as routine or expected to succeed.

The script `search_rational_pairs.py` represents numerators and denominators as
SymPy expressions and applies exact filters. Its baseline pair

\[
A=(x,y/x^m),\qquad B=(u,u^m v+H(u))
\]

has reciprocal determinants and polynomial composition `(x,y+H(x))`. It is a
positive control for pole cancellation but is invertible and has no branch
identification. This prevents a search implementation from confusing mere
cancellation with a counterexample.

`solve_rational_ansatz.py` complements the filter: it clears denominators in a
small Laurent ansatz, generates coefficient equations for polynomial extension
and constant chain-rule determinant, and solves them symbolically. Its support
and pole orders are meant to be enlarged systematically. Any returned family
still needs the branch and function-field rejection tests below.

### Search ansätze

- **Laurent/weighted:** denominators powers of a line `l(x,y)`; bound numerator
  support by weights. This directly imitates the 3D mechanism.
- **Two divisors:** denominators `l1^m l2^n`, allowing distinct valuations to
  cancel only after composition.
- **Blow-up charts:** write maps in local coordinates `(s, s^k t)` and compare
  limits on exceptional divisors. Distinct `t` values encode branches.
- **Cubic primitive element:** impose a relation
  `c T^3-2 T^2+b T-2a=0`, seek rational reconstruction of `(x,y)` from
  `(a,b,T)`, and force the repeated-root locus to correspond to a pole rather
  than finite ramification.
- **Coefficient solving:** clear denominators, equate forbidden coefficients to
  zero, enforce a constant chain-rule determinant, then use Gröbner bases,
  modular solving, or homotopy continuation. Recheck survivors over `Q`.

### Necessary rejection tests

- Reject if the composition has an uncancelled denominator.
- Reject if its determinant is zero or nonconstant.
- Reject birational compositions (Keller's theorem makes them automorphisms).
- Reject Galois generic field extensions in the Keller setting.
- Test whether apparent distinct branches are merely two parameterizations of
  the same valuation.
- Projectivize: a branch collision at a finite critical point contradicts the
  nonzero determinant; the merger must occur through infinity/nonproperness.
