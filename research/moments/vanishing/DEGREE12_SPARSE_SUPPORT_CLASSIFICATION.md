# All degree-twelve supports with at most four left weights

15 November 2026. Internal proof with exact rational polynomial identities.
No novelty claim, independent specialist replay, external review or formal
verification. The [checker](verify_degree12_sparse_supports.py),
[frozen certificates](degree12_sparse_certificates.json),
[replay receipt](degree12_sparse_verification.json), and
[source audit](degree12_sparse_source_audit.json) specify the evidence.

For **both** exceptional quotients SL_2/2O and SL_2/2I, this closes every
function in the actual V_12 having at most four nonzero left-weight
coefficients. Every pure-null member has weights of one strict sign and
is safe for all regular markers in the entire quotient algebra. The
least sufficient pure-only set prefix on this union of subspaces is11.
This is substantially larger than the single four-weight sector closed
in the [preceding theorem](DEGREE12_FOUR_WEIGHT_IDEALS_AND_QUARTIC_BARRIER.md).

Consequently a pure-null function outside the SL_2 nullcone, or an MZ
failure, lying in either actual V_12 must have **at least five weights in
every left SL_2 frame**. This does not classify the full V_12, the mixed
octahedral V_8+V_12, or either all-degree exceptional model. Those remain
OPEN. In particular, the ordinary-ideal quartic deficit from the preceding
theorem remains a deficit of ideals; no inequality of radicals follows.

## 1. Functions, weights, and the precise statement

Use g=[[x,y],[z,w]], det(g)=1, F=z partial_x+w partial_y, and

    f_O(X,Y)=X^2Y^2(X^4-Y^4)^2,
    f_I(X,Y)=XY(X^10+11X^5Y^5-Y^10),
    p_j=F^j f_Gamma(x,y)/j!,  j=0,...,12.

The checker re-verifies exact finite-group invariance, orders48/120,
and the one-dimensional fixed line in V_12. These are the complete
actual V_12 copies in the two quotient algebras. It verifies all ladder
weights 12-2j and the Weyl identity

    p_j(z,w,-x,-y)=(-1)^j p_(12-j)(x,y,z,w).             (1)

For P=sum a_j p_j let supp(P)={j:a_j!=0}. Assume |supp(P)|<=4. Then
the following conditions are equivalent:

1. h(P^m)=0 for every positive integer m.
2. h(P^m)=0 for 1<=m<=11.
3. supp(P) is contained in {0,...,5} or in {7,...,12}.

The zero function is included. Such P belongs to the SL_2 nullcone by
an explicit contracting diagonal subgroup. For every regular Q,

    h(QP^m)=0 whenever 2m>deg_matrix_entries(Q).         (2)

The coefficient1/2 in this uniform zero-offset bound is optimal on
each quotient. The theorem also applies after any common left SL_2
translate, since left translation preserves Haar integration, regular
markers and the matrix-entry degree filtration. Adding a constant to
this family adds its vanishing as the first-moment condition; it does
not change the sufficient set prefix11.

## 2. The easy supports and finite normalization

Every monomial in P^m has positive weight when all weights of P are
positive, and negative weight when all are negative. Haar invariance
therefore proves (3)=>(1) at every order. A marker of matrix-entry
degree d has weights in [-d,d]. The minimum absolute nonzero weight
here is2, proving (2) by the same argument.

If all weights are weakly of one sign and weight zero occurs, the only
weight-zero contribution to H_2=h(P^2) is a_6^2 h(p_6^2). Direct exact
integration gives

    h(p_6^2)=32/13 on2O, 275/13 on2I.                   (3)

Thus H_2=0 removes the zero-weight coefficient. This also covers
single supports and all supports without weights of both strict signs.

It remains to exclude each support i_0<...<i_(s-1) with i_0<6<i_(s-1)
and s=2,3,4, with **all** its coefficients nonzero. Scalar multiplication
and the left diagonal torus normalize the first and last coefficients
to1: the needed nonzero complex root exists because the two endpoint
weights differ. Write the normalized function

    P=p_(i_0)+sum_(r=1)^(s-2) t_r p_(i_r)+p_(i_(s-1)),
    every t_r !=0.                                      (4)

Scalar multiplication rescales H_m by its mth power, and the torus
preserves H_m. This normalization therefore preserves every zero-moment
condition. The Weyl identity (1) reflects supports by j->12-j, preserving
nonzero coefficients and permitting another normalization. The exact
coverage counts per model are:

| Support size | All crossing supports | Weyl representatives checked |
|---|---:|---:|
| 2 | 36 | 21 |
| 3 | 216 | 111 |
| 4 | 645 | 330 |

The checker constructs the sets of all supports and verifies equality
with the union of these representatives and their reflections. This
is coverage of every support, not a random sample of coefficients.

## 3. Exact moment polynomials and unit-ideal identities

For fixed support write w_r=6-i_r. If the interior exponents in a term
of P^m are b_r, its endpoint exponents are uniquely determined by

    a=(-w_last*m-sum_r (w_r-w_last)*b_r)/(w_first-w_last),
    d=m-a-sum_r b_r.                                    (5)

Only integral nonnegative a,d contribute. The coefficient of
product t_r^b_r is exactly

    m!/(a! d! product b_r!)
      *h(p_first^a product p_(i_r)^b_r p_last^d).        (6)

The imported sparse four-variable arithmetic multiplies these functions
and integrates every monomial using the retained exact SU_2 Haar rule.
No numerical integration, algebraic sampling or rounding occurs.

For every nonzero moment polynomial H_m in (4), the certificate stores

    H_m=c_m product t_r^v_(m,r) G_m(t), c_m in Q*,       (7)

where the removed monomial uses the componentwise minimum exponents,
and G_m is normalized to leading coefficient1. Since every t_r is
nonzero, H_m=0 is equivalent to G_m=0 on this normalized stratum.
For two weights there are no interior variables and the product below
is1.

For each representative the frozen certificate supplies explicit
rational polynomials A_m(t,z),B(t,z) satisfying

    sum_m A_m(t,z)G_m(t)
       +B(t,z)*(z product_r t_r-1)=1.                   (8)

These are polynomial identities, replayed by direct sparse
multiplication and coefficient equality. Taking z=(product t_r)^(-1)
proves that the indicated moment equations have **no** common zero
with all coefficients in (4) nonzero. The generator for saturation is
z-1 when there are no interior variables.

The constructor uses a tracked exact Buchberger procedure with rational
arithmetic and a product criterion. The default checker does not trust
that procedure: it recomputes (5)--(7) and directly checks the frozen
identity (8). Hence a mistake in the constructor cannot produce a
passing false certificate. The replay still shares the previously
audited Haar arithmetic and is not an independent specialist replay.

There are462 identities per model,924 in total. The largest orders
needed are11 for two weights,10 for three, and9 for four. Thus if all
first11 moments vanish, no crossing support of size<=4 is possible.
Together with (3), this proves (2)=>(3) and the whole equivalence.
This is a finite exact inconsistency proof, not an inference that
uncomputed moments vanish after a long zero prefix.

## 4. Sharp bounds and the frame-independent obstruction

The preceding theorem's P_*=p_0+p_11 already lies in this family.
Its weights12,-10 force its first ten pure moments to vanish, and

    h(P_*^11)=-43987763200000/23813095609334142593573 on2O,
               23788050537109375/807223579977428562494 on2I. (9)

The new checker replays the weight balances and both nonzero values.
Consequently the least sufficient pure-only **set** prefix on this
sparse union is exactly11. This does not assert ordinary moment-ideal
generation through11, or sufficiency of11 on the whole V_12.

For sharpness of (2), choose P=p_5 and Q=p_12. Their weights are2,-12,
and deg(Q)=12. Direct integration gives

    h(p_12*p_5^6)=-536540715417600/436092044389001 on2O,
                   17610926850000000/39644731308091 on2I. (10)

These values exclude every smaller uniform zero-offset coefficient
than1/2. Powers of P alone vanish at every order by their positive
weight. No claim about optimal affine bounds with an extra constant
offset is made.

Finally suppose a pure-null P in either V_12 can be left-translated
to at most four weights. The theorem makes the translate torus-null
and safe for every regular marker. Translating back preserves both
properties. Therefore every pure-null P outside the SL_2 nullcone,
and every MZ counterexample in this V_12, has at least five nonzero
weights in every left frame. This is a necessary support condition,
not a converse or a classification of the five-or-more-weight strata.

## 5. Reproduction and remaining frontier

From the repository root, run

    python3 research/moments/vanishing/verify_degree12_sparse_supports.py

It uses only the Python standard library. The explicit --freeze flag
reconstructs the polynomial identities; it is not needed for replay.
The receipt records hashes, coverage, all orders, exact sharp controls,
and the open full-space boundary.

The earlier special sector retains its stronger ordinary-ideal and
prefix6 conclusions. The present theorem extends the safe-support
classification to all supports of size<=4 and to their common left
translates. Neither result handles arbitrary sums with five or more
weights, nor a nonzero V_8 coefficient on the octahedral mixed target.
A next exact calculation must address such a new stratum or prove a
full radical/nullcone criterion; rerunning already excluded supports
with longer prefixes adds no mathematical information. A genuine MZ
failure still requires a pure-null function with a certified infinite
marked tail (equivalently the actual entire-period rank excess in the
[period criterion](ENTIRE_HAAR_PERIODS_AND_EXCEPTIONAL_DECISION.md)).
