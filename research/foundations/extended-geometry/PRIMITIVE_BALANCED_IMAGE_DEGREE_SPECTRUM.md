# The primitive balanced-degree spectrum in two contraction pairs

Over every characteristic-zero field, the balanced degrees of radially
primitive, absolutely nonpower two-pair Image counterexamples are precisely
\(\{4,6,7,8,\ldots\}\). Every pure, radially primitive balanced quintic
is safe against every fixed polynomial marker, including when the phase
is a proper power.

## Exhaustion of quintics

For a radially primitive balanced phase of degree five, restriction to the
Hopf quadric \(t^2+2xy=1\) has least degree five. Purity and full-marker
safety transfer between this restriction and original contractions by the
[rank-one trace and marker formulas](RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md).
The leading degree-five symbol on the null cone gives a nonzero binary
decic. Its geometric root multiplicities exhaust the following cases.

- **All multiplicities at most four.** The [proper period obstruction](STABLE_LEADING_SPHERE_PHASE_PERIOD_OBSTRUCTION.md)
  gives infinitely many nonzero moments. It uses a sublevel-area estimate to
  extend the spherical volume to a proper model of the phase map.
- **A root of multiplicity exactly five.** The
  [leading coefficient calculation](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#4-proof-of-the-elementary-obstructions)
  makes the coefficient at the unique prime interval denominator a nonzero
  power of the central angular coefficient.
- **A root of multiplicity at least six.** The
  [unstable leading-form theorem](UNSTABLE_QUINTIC_SPHERE_PURITY_EXCLUSION.md)
  treats all lower coefficients. It proves that every pure phase in this case
  is safe. Its arithmetic input is the
  [rank-two CM ordinarity criterion](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md).
  The finite supporting-line census and the cyclic signatures are replayable
  by [an exact checker](verify_unstable_quintic_faces.py).

For a radially divisible balanced quintic, the restriction has smaller
least degree. In particular, multiplying the unsafe quartic by the radial
bilinear gives an unsafe balanced quintic of least sphere degree four.

## The other degrees

The [cubic Image theorem](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md) includes
the full degree-three coefficient space and arbitrary original markers.
Radial padding to degree three and the contraction normalization give the
same exclusion in lower degrees. Nonzero constants fail purity.

The [joint endpoint identity and lift](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md)
give integer examples in degrees \(4+2a+3b\), \(a,b\ge0\). They have a
common linear marker with a nonzero rational factorial value at every
positive order, a nonzero specialization on the radial zero set, and a
simple-root specialization which excludes absolute proper powers.

## Exact exceptional parameters and finite moment checks

The two \(T=A^2/(BC)\) charts (I) and (II) have exceptional polynomial \(1\)
by the [rank-two face theorem](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md#the-two-quadratic-face-families).
The even \(D=0,\ BC\ne0\) chart has the same exact polynomial. The full
even \(BD\ne0\) chart has empty geometric pure image, with reduced
elimination ideal \((1)\) in its two parameters. The latter two statements
follow from the
[even leading-form calculation](EVEN_STRIPPED_QUINTIC_PRYM_AND_EXACT_MOMENTS.md).
These are all-order geometric conclusions.

The [rational moment generator](compute_quintic_prime_face_moments.py)
returns exact sparse polynomials in Singular syntax. The first six moments
in (I) and the first four in (II) still have zero parameter elimination
ideal, as verified by triangular rational sections. The even computations
give finite-prefix exclusions in the \(q=0\) and quartic-pivot-zero charts
by exact special-fiber Gröbner bases and weighted projective properness.
The [even calculation](EVEN_STRIPPED_QUINTIC_PRYM_AND_EXACT_MOMENTS.md)
gives the equations, scope, and replay command. These computations and the
finite Newton census support the argument without a universal purity cutoff.

The [standalone manuscript](../../../unpublished/primitive-balanced-image-spectrum/main.tex)
assembles the three leading-form cases, the low-degree exclusions and the
common-marker integer family.
