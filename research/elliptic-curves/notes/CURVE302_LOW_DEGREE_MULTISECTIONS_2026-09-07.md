# Curve 302: canonical low-degree multisection quotient

<!-- status-consumer: EC-CURVE302-PARENT-DEGREE2-MULTISECTION-QUOTIENT 2f3f052d8ad9cebe -->
<!-- status-consumer: EC-CURVE302-PARENT-CHEAPEST-LATTICE-BISECTION 77f2987cfea994a1 -->

The recovered determinant-1092 parent has a particularly useful feature for
this question: its full **geometric** Mordell--Weil lattice is known and is
rational.  The parent has 24 irreducible `I1` fibres, and

\[
 \operatorname{NS}(X)=U\oplus(-M),\qquad
 \operatorname{rank}M=17,\quad\det M=1092.
\]

Its fibre at `t=0` is literally curve 302, and the displayed rank-31 group
has `D/sp(M)=Z^14`.  This gives a canonical low-degree search space before an
equation ansatz is chosen.

## Quotient, not a height box

Write

\[
 C=dO+bF+\phi(w),\qquad
 C^2=-2d^2+2db-\langle w,w\rangle.
\]

For arithmetic genus `g`, integrality forces

\[
 b=d+\frac{\langle w,w\rangle+2g-2}{2d}.
\]

Translation by the section indexed by `x\in M` is

\[
 w\longmapsto w+dx,
 \qquad
 b\longmapsto b+\langle w,x\rangle+
 \frac d2\langle x,x\rangle.
\]

Thus the finite object at fixed degree is `M/dM`, not a coefficient or
height cutoff.  For a coset `c`, put

\[
 \mu_d(c)=\min_{w\in c}\langle w,w\rangle.
\]

Completing the square gives

\[
 \min_x C\mathbin{.}S_x=
 \frac{\mu_d(c)}{2d}-d+\frac{g-1}{d}.
\]

The all-section nonnegativity gate is therefore
`mu_d(c) >= 2d^2-2g+2`.  The genus/integrality condition remains
`mu_d(c)+2g-2 == 0 (mod 2d)`.

## Complete degree-two result

The reproducible census visits all `2^17=131072` translation orbits, using a
complete CVP search at double-double and MPFR precision and checking every
returned vector over the integral Gram form.  The complete minimum spectrum
is:

| minimum norm | translation orbits |
| ---: | ---: |
| 0 | 1 |
| 4 | 1,218 |
| 6 | 24,875 |
| 8 | 63,922 |
| 10 | 40,917 |
| 12 | 139 |

Consequently there are **40,917 rational-bisection translation orbits** and
**64,061 section-nonnegative genus-one-bisection lattice classes**.  All
40,917 rational classes have minimum norm 10.

For a rational degree-two survivor, Riemann--Roch gives an effective divisor
over `Q`.  The fibre degree is positive, so the opposite class cannot be
effective.  Rootlessness makes every vertical component a multiple of `F`.
The square `-2` then rules out a vertical summand; if there were two section
components, its intersection with either would be negative, contrary to the
all-section gate.  Hence each of these 40,917 classes is a geometrically
irreducible smooth rational bisection.  The genus-one count is deliberately
weaker: it is a complete lattice/divisor candidate count, not an assertion
that every class is nef or irreducible.

This already changes the scale of the question.  The earlier 178
conic/three-line twisted-cubic bisection equations were a particular
construction family, and all were nonsplit at `t=0`; they have not been
identified with, nor promoted to, a complete quotient census.

The certificate is
[`curve302_parent_degree2_multisection_lattice_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_parent_degree2_multisection_lattice_v1.json),
with one deterministic minimum representative per exported class in
[`curve302_parent_degree2_multisection_orbits_v1.tsv`](../../artifacts/generated-results/elliptic-curves/curve302_parent_degree2_multisection_orbits_v1.tsv).

## Equation-to-quotient gate

Lattice enumeration does not factor a fibre.  For each exported rational
orbit, the next certificate must:

1. construct the indicated bisection over `Q`;
2. verify generic irreducibility and the attached NS orbit;
3. factor its degree-two fibre over `t=0` exactly;
4. map each rational factor point into `D/sp(M)=Z^14`; and
5. prove any claimed new quotient direction independent of the specialized
   generic image.

There is a uniform construction target for the minimum rational shell.  If
`w^2=10`, then `P_{-w}.O=3` and

\[
 (2O+4F+\phi(w))+P_{-w}=3O+9F.
\]

So a resolved Riemann--Roch calculation in `H^0(3O+9F)` through the explicit
section `P_{-w}` should leave the required residual bisection.  This is the
right common equation interface: it works from the lattice representative and
the 17 explicit generic sections, rather than from an ad hoc configuration of
lines and conics.  A residual quadratic that does not split at zero is an
exact nonsplit result for that constructed orbit only; it says nothing about
unconstructed orbits.

### One blind equation-side control

The first deterministic representative under the rule “minimize `l1`, then
`linfinity`, then the orbit mask” is orbit `8044`, with

\[
 w=(0,-1,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0),\qquad w^2=10.
\]

This rule reads only the exported lattice quotient; it does not read `t=0`
or the exceptional curve-302 points. Exact interpolation through `P_{-w}`
has rank 19 in the 20-dimensional displayed linear system and produces its
primitive unique relation. Removing the known trace factor leaves a
generically irreducible quadratic. The residual quadratic on the fibre at
`t=0` is separable and irreducible over `Q`. Thus this bisection is an exact
nonsplit control, with no rational fibre point to attach to `D/sp(M)`.

The complete equation, elimination data, and exact zero-fibre factorization
are in
[`curve302_parent_cheapest_lattice_bisection_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_parent_cheapest_lattice_bisection_v1.json).
It demonstrates the common RR/elimination interface, but neither samples nor
discards any remaining orbit.

For orientation, the raw quotient sizes are `3^17=129140163`,
`4^17=17179869184`, and `5^17=762939453125`.  They are finite but not all
equally cheap.  The degree-three, four and five stages therefore need
checkpointed residue shards, a frozen ordering independent of the curve-302
exceptional points, and streamed equation construction; no finite prefix is
to be described as an exhaustive degree-`d` result.

## Reproduction

From `research/`:

```sh
sage -python elliptic-curves/cas/enumerate_curve302_parent_degree2_multisections.py

sage -python elliptic-curves/cas/enumerate_curve302_parent_degree2_multisections.py --check

sage -python elliptic-curves/cas/construct_curve302_parent_cheapest_lattice_bisection.sage --check

python3 -m unittest elliptic-curves.tests.test_curve302_parent_degree2_multisections
```

The result is a geometric supply theorem for degree two.  It does not yet
identify a split fibre, a missing quotient class, or a rank gain.
