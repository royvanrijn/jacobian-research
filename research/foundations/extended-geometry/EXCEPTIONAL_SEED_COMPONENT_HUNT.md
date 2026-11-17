# Low-degree exceptional-seed component hunt

This note stress-tests the exceptional-seed dimension formula against
scheme-theoretic elimination, rather than point samples.  The executable
audit is

```bash
.venv/bin/python scripts/hunt_exceptional_seed_components.py
```

It runs exact characteristic-zero eliminations in degrees four through seven.
For each full-contact partition it uses the quotient-coordinate incidence from
`jcsearch.discriminant_geometry`, eliminates the root coordinates with
Singular, takes the saturated projective closure in the normalized seed
coordinates, and compares its Hilbert data with the maximal `2/3` images.

## Results

No additional reduced component occurs through degree seven.  Every
full-contact partition image is contained in the union of the images indexed
by partitions using only twos and threes.  The maximal-component data are:

| degree | type | dimension | projective degree | Hilbert polynomial |
|---:|:---:|---:|---:|:---|
| 4 | `(2,2)` | 1 | 1 | `m+1` |
| 5 | `(3,2)` | 1 | 6 | `6m-9` |
| 6 | `(3,3)` | 1 | 6 | `6m-2` |
| 6 | `(2,2,2)` | 2 | 4 | `2m^2+2` |
| 7 | `(3,2,2)` | 2 | 24 | `12m^2-22m-111` |

The dimensions agree with `ell(lambda)-1`.  The degrees agree independently
with

\[
 \deg C_{a,b}=\binom{a+b}{a}
 \left(2^a3^b-(a+b+1)\right).
\]

The reduced union Hilbert polynomials are `m+1`, `6m-9`,
`2m^2+6m-24`, and `12m^2-22m-111` in degrees four through seven.
In degree six the lower-dimensional `(3,3)` component changes the linear and
constant terms of the union polynomial even though it does not change the
top-dimensional degree.  This is one reason a dimension or generic-point
sample alone is insufficient.

## The five proposed failure modes

### Nonreduced primitives

The maximal component ideals in degrees four, five, and six are radical and
have one minimal prime.  Thus no nonreduced primitive component appears in
those degrees.  The degree-seven component has the predicted irreducible
support and Hilbert data; a direct characteristic-zero `radical` call on its
large implicit ideal exceeds the practical memory budget of this audit, so
the executable deep-primary check stops at degree six by default.

There is nevertheless genuine nonreduced geometry in an intersection.  In
degree six, the projective scheme intersection

\[
 C_{2,0}\cap C_{0,2}
 =C_{(2,2,2)}\cap C_{(3,3)}
\]

has length 24 and reduced support of degree five.  After restricting to the
exact-degree admissible seed open, its reduced support has degree four but its
scheme length is eight.  Thus each admissible common-collision point carries
the expected double intersection structure.  This does not create another
irreducible component of the nonsurjective locus.

### Collisions at infinity and exact-degree loss

No maximal component is supported in the projective hyperplane at infinity:
intersecting each component with that hyperplane lowers dimension.

There is, however, a substantial closure-only degree-drop contribution in the
degree-six pairwise intersection.  Its total length 24 is entirely in the
affine homogenizing chart.  Saturating by the exact-degree factor reduces the
length to eight, while saturating only by the weighted forbidden divisor does
nothing.  The discarded length 16 is supported at the single reduced point

\[
 [h_3:h_4:h_5:z]=[-1:0:0:1],
 \qquad h_6=0.
\]

This is the lower-degree seed `H=W^2-W^3` viewed on the degree-six closure.  It
is best interpreted as support escaping to root infinity/degree drop.  It is
outside the exact-degree space `A_6`, so it does not contradict the claimed
dimension of `N_6`, but any scheme-theoretic compactification statement must
include its multiplicity.

### Hidden normalization charts

For every maximal type through degree seven, `Phi` is coprime to both the
collision discriminant and the normalized admissibility factors `D` and
`M''(1)-2D`.  Therefore removing those divisors removes no component of the
root hypersurface, and the raw and collision-retaining closures agree.  The
degree-six calculation above shows why the exact-degree factor must still be
remembered when two projective seed closures are intersected.

### Components appearing only after closure

The intersection of the projective ideals of all full-contact partition
images equals the intersection of only the maximal `2/3` ideals through
degree seven.  Every tested coarsened partition satisfies the reverse ideal
containment predicted by collision refinement.  Closure creates embedded or
nonreduced intersection structure, including the length-16 degree-drop point,
but no new minimal prime.

### Automorphism and symmetry components

Equal-multiplicity roots are represented by elementary symmetric coordinates.
Independent ordered-root eliminations through degree six give exactly the
same projective seed ideals as the quotient-coordinate eliminations.  Fixed
loci of the permutation action specialize to root-collision boundaries; they
do not produce a separate image component.  In characteristic zero the
ordered-to-quotient map is finite, so it cannot change dimension.

## Conclusion and remaining boundary

The scheme-level low-degree evidence supports

\[
 \dim N_N=\left\lfloor N/2\right\rfloor-1
\]

and the phase transition “generic quartic exceptional, generic degree
`N>=5` seed surjective.”  The search found no countercomponent.  Its positive
finding is instead a warning about strengthening the theorem: component
support is generated by the `2/3` atoms, but projective closures and component
intersections are not reduced.  Already in degree six, exact-degree loss adds
a length-16 boundary contribution and the admissible common-collision scheme
has length two at each reduced point.

Full degree-eight implicitization and primary decomposition remain outside
this script's current memory envelope.  The existing degree-eight
coefficient-comparison obstruction still rules out off-diagonal points, but
it is not a substitute for the Hilbert-polynomial comparison performed here
through degree seven.
