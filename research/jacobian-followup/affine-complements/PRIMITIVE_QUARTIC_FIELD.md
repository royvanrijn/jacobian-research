# The quartic normalization field is primitive with symmetric monodromy

14 November 2026. PF2PRIM1 applies to both surviving quartic packets
without assuming a ruled completion or a smooth normalization.
It rules out every proper intermediate field, not just a polynomial
intermediate plane.

The geometric input is purity of the branch locus,
[Stacks Lemma58.21.4, tag0BMB](https://stacks.math.columbia.edu/tag/0BMB),
with a [pinned copy](sources/primitive-cover-manifest.json). Purity is
imported with its normal-source, regular-target, finite/quasi-finite and
codimension-one unramified hypotheses. We also use triviality of connected
finite etale covers of complex A2.

## A general field-primitivity statement

Let pi:N -> A2 be a finite dominant map from an integral normal complex
surface. Suppose its full ramification divisor is one prime R, with
generic ramification index two and residue degree one over its image C.
Put F=C(A2) and K=C(N).

Then there is no field F strictly contained in L strictly contained in K.

To prove this, normalize A2 in a proposed L, obtaining a normal finite
surface M. The normalization ring of M is finite over the target ring, which is
excellent, and is contained in the normalization ring of N. A finite set
of target-module generators of the latter is also a set of generators
over the former. Thus N -> M is finite.
Let h:M -> A2 and g:N -> M be the finite maps.

If h is unramified at every codimension-one point, purity makes it etale
everywhere. The source is normal, the target regular, the map finite,
and equal local dimensions hold for these finite surface covers.
The connected finite etale cover of A2 is trivial, so L=F.

Otherwise choose a ramification prime Q of h, with image C0. Every
prime of N above Q ramifies for the composite, since ramification indices
multiply. There can therefore be only the one prime R above Q, and C0=C.
Multiplicativity gives

\[
 2=e(R/Q)e(Q/C),\qquad
 1=f(R/Q)f(Q/C).
\]

Since e(Q/C)>1, we have e(Q/C)=2, e(R/Q)=1 and both residue degrees one.
The complete finite DVR degree identity at Q now gives

\[
 [K:L]=\sum_{R_i\mid Q}e(R_i/Q)f(R_i/Q)=1,
\]

because there is exactly one prime above Q. Hence L=K, contradicting
a proper intermediate field.

The completeness of that DVR sum is essential. It is taken in the finite
normalization N, not in the source open with deleted primes.

## Symmetric monodromy

The action of the Galois-closure group G on the d=[K:F] sheets is
primitive, by the intermediate-field result. A meridian at a general
point of C acts by a single transposition: the sole ramified prime has
index two and residue degree one; every other point over that general
branch value is unramified.

A primitive permutation group containing a transposition is S_d.
For completeness, form the graph whose edges are all G-conjugates of
that transposition. Its connected components form a G-invariant block
system. There is an edge, and primitivity forces the graph to be
connected. Edge transpositions of a connected graph generate S_d, so
G already contains S_d.

Thus the generic monodromy and Galois-closure group of either actual
quartic packet are S4. No clean-cusp or chosen matching hypothesis is
needed. This recovers the old local S4 conclusion in the clean
one-boundary case and extends the receiver restriction to both packets.

## What is and is not excluded

Any factorization of the quartic cover through a normal finite surface
with two generic degrees greater than one is excluded. In degree four,
this eliminates all degree-two intermediate fields, regardless of whether
their integral closures are polynomial rings.

The [prime-ramification polynomial-receiver theorem](PRIME_RAMIFICATION_FACTORIZATION.md)
also handles a different endpoint issue: a polynomial receiver with the
same function field as N would be a finite birational plane receiver and
would force N to be a plane, contradicting its canonical class.
Field primitivity alone does not handle that endpoint.

The sign/discriminant quadratic field in the S4 Galois closure is still
allowed. It is not an intermediate field of the non-Galois quartic field:
the point stabilizer S3 is not contained in A4. Therefore this result does
not remove the old Cox orientation/discriminant base-change programme.
Cubic resolvent fields likewise belong to the Galois closure, rather than
the quartic field itself.

Actual ruled compatibility, higher-pole map realizability and either
full quartic packet remain unresolved.

The [checker](verify_primitive_cover.py) verifies the complete degree-two
DVR profile list, enumerates all subgroups of S4, checks the primitive
transposition consequence and distinguishes the sign field. These are
finite controls. Normalization, purity, generic inertia and the
all-degree field argument above are written proofs, without formal
verification, independent whole-proof implementation or external review.

