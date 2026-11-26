# Conductor-eight ambient families and period obligations

For an actual plane Keller map in the monic source frame, the global
coordinate-order conductor has even degree. Degrees two, four and six are
excluded; the [hyperelliptic octic branch](HYPERELLIPTIC_OCTIC_CONDUCTORS_AND_FOUR_PERIOD_RIGIDITY.md)
is also excluded. This note fixes the remaining degree-eight research tasks.
It asserts no new exclusion.

The [Nakazawa–Oka Main Theorem 1.9, printed p. 670](sources/nakazawa-oka-genus-four-1997.pdf)
classifies **minimal smooth affine plane curves of compactification genus
four with one place at infinity**, in ambient coordinates. Their four
normal-form families are parameterized, not four fixed equations. Write
`h0=v` after the allowed linear change, with each `ξi≠0`. Coefficients of
each displayed `R` vary freely in the indicated linear space, subject to
the theorem's smoothness conditions.

| Ambient class | Full parameterized normal form at infinity | Period task on a singular rational specialization |
|---|---|---|
| N1, `k=1`, `P1=(5,2)` | `h=h0^5+ξ1 u²+R1`, with `R1` in the span of `u³,u⁴,u⁵,h0u²,h0u³,h0u⁴,h0²u²,h0²u³,h0³u,h0³u²` | Derive all collision partitions and the physical action vector for the entire ten-parameter space; test nonzero action and the Keller period derivative on every singular stratum. |
| N2, `k=1`, `P1=(9,7)` | `h=h0^9+ξ1 u⁷+R1`, with `R1` spanned by `h0^a u^b` where `a+b≤9`, `9b+7a>63`, `a<8` (13 coefficients) | Identify exactly which specializations have the already excluded order `K[u²,u g(u²)]`; compute actions for any remaining ambient embeddings rather than equating the class with that order. |
| N3, `k=2`, `P1=(3,1), P2=(2,7)` | `h1=h0³+ξ1u+R1`, `R1∈span{u²,u³,h0u,h0u²}`; `h=h1²+ξ2u⁴h0+R2`, `R2∈span{u⁵,u⁶,h0u⁵,h0²u⁴}` | Use the approximate-root frame to derive the matching endpoints and four actions, retaining all four `R2` parameters and every singular coefficient locus. |
| N4, `k=2`, `P1=(3,1), P2=(3,16)` | The same four-parameter `h1`; `h=h1³+ξ2u⁸h0+R2`, `R2∈span{u⁹,h1u⁶}` | Compute the two-parameter upper layer and the lower `h1` dependence together; test whether action ratios obey the required common Keller derivative. |

The formula for N2 is the monomial condition printed in Theorem 1.9;
the table must be used with the source's affine degree and minimality
hypotheses. The [Assi–García-Sánchez approximate-root and δ-sequence
source](sources/assi-garcia-sanchez-1407.0490v1.pdf), Propositions 2 and 4
and Section 6, supplies a way to enumerate possible one-place semigroups
and to track distinct embeddings. Its numerical-semigroup genus is not by
itself a proof that a particular Keller order belongs to one of N1–N4.

## Transfer gate before applying the list

An exact specialized Keller order is an immersed singular rational plane
curve with normalization `C[x]` and defect four. Nakazawa–Oka's hypothesis
instead concerns a **smooth** affine curve of geometric genus four.
The [sextic proof](GLOBAL_SEXTIC_CONDUCTOR_EXCLUSION_AND_PERIOD_POLE_REPAIR.md)
uses a good-at-infinity smoothing bridge for defect three. Reprove the
corresponding bridge at defect four: a smoothing of all finite singularities
must preserve the embedded one-place type at infinity, have genus four,
and land in a *minimal* ambient normal form after admissible polynomial
target changes. Then determine which singular parameter loci in N1–N4
admit normalization by one affine line and conductor degree eight. Until
this gate is checked, the four families are candidate tasks, not an
exhaustive classification of Keller octic orders.

## Physical action and descent tests

For each admitted singular specialization, compute the exact conductor
and matched endpoint scheme in the original source coordinate. On a
simple four-pair open set, for each matched pair `x_i^-,x_i^+` define

`I_i(y)=∫_{x_i^-}^{x_i^+} P(x,y) ∂_xQ(x,y) dx`.

The retained physical-action calculation gives `I_i'=-J(x_i^+-x_i^-)`.
Target affine and triangular automorphisms scale every action by the
same determinant. For N1, N3 and N4, derive the canonical action vector
as a function of *all* normal-form parameters, then test whether a
single nonzero rational scaling can satisfy these four derivative
equations. A vanishing action vector, colliding endpoints, or a
parameter denominator must be handled on its own stratum. For N2,
first compare the exact order with the hyperelliptic theorem before
using its existing obstruction.

Finally repeat the finite filtered-order descent of the sextic proof over
the parameter line. Establish regularity of the chosen target changes
across every pole, root collision and degree drop, using monicity and the
exact specialized conductor rather than a generic-family argument.
This last step is needed even if the generic period equations yield a
contradiction. The output of a successful class calculation is an exact
excluded family with its degeneration proof, not an unrestricted
conductor-eight conclusion until all admitted ambient classes are covered.
