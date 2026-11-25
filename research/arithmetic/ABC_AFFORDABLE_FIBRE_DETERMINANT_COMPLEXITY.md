# Determinants of distinct targets in an affordable fibre

A pair of actual targets in one affordable record gives a new
primitive additive triple by crossing the two residual summands.
Counting those triples forces a lower bound on the radical of the
typical normalized cross determinant in a large matching. The
argument does not require a common-neighbor star. It gives a saving
when the determinants are sufficiently small-radical and identifies
the remaining large-radical case.

## Fixed records and normalized determinants

Use the iterated balanced box conventions with \(n=X^L\).
Take the inner height limit at fixed \(L\) before \(L\to0\).
On a fixed affordable pair record write
\[
 U_i=F u_i,\qquad V_i=G v_i,\qquad
 T_i=\sigma F u_i+\tau G v_i=D w_i
\]
with common signs \(\sigma,\tau\), fixed positive \(F,G\), and
an optional fixed positive target factor \(D\). The complete
letters are positive and pairwise coprime. Assume that the
targets are distinct. Affordability supplies a fixed positive
numerical weight on at least one of \(F,G\). In the fibre all
original letter root caps and numerical height caps are bounded
in the usual units. A matching in what follows means that
the \(u_i\)'s are distinct and the \(v_i\)'s are distinct.
Record one representative per target and refine to dyadic
numerical boxes for \(u_i,v_i\) at subpower cost. Denote their
numerical height exponents by \(E_u,E_v\).

For \(i\ne j\) set
\[
 \Delta_{ij}=u_i v_j-u_j v_i,\qquad
 g_u=\gcd(u_i,u_j),\quad g_v=\gcd(v_i,v_j),\quad
 g_t=\gcd(T_i,T_j),\qquad
 d_{ij}=\frac{|\Delta_{ij}|}{g_u g_v}.
\]
The determinant is nonzero by primitivity of each pair and
distinctness of the targets. The identity
\[
     \gcd(u_i v_j,u_j v_i)=g_u g_v
\]
holds prime by prime. A prime present on both edges in the
same letter contributes its smaller valuation. A prime in
opposite letters puts both factors on the same side and
contributes no gcd. Thus \(d_{ij}\) is the positive third
term in a primitive additive relation between the two
cross-products divided by \(g_u g_v\). No global prime
colouring is needed. The normalized determinant is unchanged
if a further *fixed* divisor is transferred from each
residual summand into the corresponding fixed factor.

The target gcd divides even the normalized determinant:
\[
                         g_t\mid d_{ij}.                    \tag{DC1}
\]
Indeed \(u_iT_j-u_jT_i=\tau G\Delta_{ij}\), and
\(v_iT_j-v_jT_i=-\sigma F\Delta_{ij}\).
Since \(\gcd(F,G)=1\), the common target gcd divides
\(\Delta_{ij}\). It is coprime to \(g_u g_v\) by primitivity.

### A fixed-anchor bound

Fix \(i\). For prescribed \(\Delta_{ij}\), the integral
cofactor pairs lie on an affine line with primitive
integer difference direction \((u_i,v_i)\). Both
coordinates of the anchor are comparable to their upper
bounds in the dyadic box, so there are \(O(1)\) pairs
on this line. The three gcds are divisor-many choices
from the anchor values. For each such record the raw
determinant is a nonzero multiple of \(g_u g_v g_t\).

Assume \(D\) records the full parts on its prime support
and write \(\operatorname{rad}_{\not\mid D}\) for the
radical outside that support. Uniformly for a fixed anchor,
the number of neighbours satisfying
\[
 \operatorname{rad}_{\not\mid D}(d_{ij})\le n^\alpha,\qquad
 \operatorname{rad}(g_t/D)\ge n^\beta,\qquad
                   \log_X(g_u g_v g_t)\ge w
\]
is at most
\[
        X^{o_X(1)}
            \min\{n^{\alpha-\beta},\,X^{E_u+E_v-w}\}.
                                                               \tag{DC2}
\]
If \(\alpha<\beta\) there are no neighbours. A bound with
vanishingly thin height bins means the usual arbitrarily
small height allowances. For the first estimate fix the
gcds and strip the full part of \(d_{ij}\) supported on
\(g_t\) at subpower cost. By (DC1) its remaining radical
cap is \(n^{\alpha-\beta}\), since \((D,g_t/D)=1\).
Apply the one-variable radical count. Each resulting
\(d_{ij}\) determines \(\Delta_{ij}\) up to the sign.
The second estimate counts multiples in the raw determinant
height range. In particular, a fixed target divisor genuinely
reduces the available numerical height by \(\log_X D\)
when counting the quotient determinant. It also gives
\(d_{ij}/D\le X^{E_u+E_v-\log_X D+o(1)}\).
[Uniform radical and fixed-support counts](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

## Removing concentrated support from a large fibre

Fix losses \(\epsilon>0\) and thresholds \(\eta>0,C>3\).
After losing at most \(\epsilon\) in normalized exponent,
one can enlarge the fixed full factors \(F,G,D\) so that the
following holds on the resulting subfamily. For every fixed
integer \(J\le X^C\), all but an \(n^{-\delta}\) fraction
of edges have both
- total numerical weight less than \(\eta\) in their full
  \(J\)-supported parts of \((u_i,v_i,w_i)\),
- total normalized radical mass less than \(\eta\) in those
  same parts.

Here \(\delta>0\) depends only on \(\epsilon,\eta\), and the
supported-part costs are taken in the inner height limit. To
see this, if a pool violates the assertion, retain its large
exceptional class and fix the values of the three supported
parts at subpower cost. Include these full parts in \(F,G,D\).
First record the full parts on the supports of the initial
locks at subpower cost if necessary. No locked prime can
then reappear in a residual.
At each step the sum of the fixed numerical weights and
the normalized fixed radical masses increases by at least
\(\eta\), and that potential is bounded by \(4+o(1)\)
at the balanced budget. The process has bounded length.
Choose \(\delta\) small compared with \(\epsilon\) divided
by that length bound. The pools and the fixed values may
depend on \(X\) and on the current subfamily. Fixing each
pool here is an existential thinning, with no enumeration
of possible pools. The supported-part estimate is uniform
over all pools up to the indicated polynomial height.

For each anchor take \(J=u_i v_i w_i\). It follows that
all but a vanishing fraction of ordered pairs satisfy
\[
 \log_X(g_u g_v)<\eta,\qquad
                  \log_n\operatorname{rad}(g_u g_v)<\eta.
                                                               \tag{DC3}
\]
In fact one obtains the same combined upper bounds including
the *residual* target gcd \(g_t/D\). The assertion also
allows selection of any fixed number of edges with small
full overlaps across all residual letters, by taking \(C\)
large enough to include the product of their anchor pools.

## Positive height after division

The small-budget exceptional-set count uses a radical budget
relative to the height of its primitive triple. A fibre with
positive normalized exponent bounded below cannot have
\(E_u+E_v\to0\) after the support locks above. The assertion
holds uniformly in the iterated limit for fibres of exponent
at least any fixed \(r_0>0\), with a fixed lower margin for the
numerical weight of one initial \(F,G\).

Indeed otherwise pass to a sequence with \(E_u+E_v\to0\)
and rebin the residual target valuations to obtain their limit measure. If it has positive
nonatomic mass, choose a full target block of arbitrarily
small normalized root cost and fixed positive numerical
weight. For a fixed value of this block, the two cofactor
pairs in any two edges have their cross determinant
divisible by it, as in (DC1). The cross determinant has
vanishing numerical height. It is zero, so the primitive
cofactor pair is unique. This contradicts exponent at least
\(r_0\). If the residual target measure is purely atomic,
take a positive-weight full modulus from \(F\) or \(G\).
The original radical cap permits discarding its low orders
while retaining positive weight and fullness of order \(1/L\).
Group the residual target roots in shrinking windows around
finitely many atoms, freezing the remaining tail at arbitrarily
small normalized cost. The other free residual summand has
height \(X^{o(1)}\), so it belongs to the wider coordinate in
the amplified full-modulus count for the opposite ratio.
That count again gives arbitrarily small normalized exponent
and the ratio recovers the primitive triple.
[Amplified full-modulus count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

Thus choose \(\eta\) in (DC3) below half a uniform positive
lower bound for \(E_u+E_v\), retaining a fixed positive
fibre exponent in the thinning. On the resulting good pairs,
the two primitive cross-products after division by \(g_u g_v\)
still have height \(X^\rho\) for \(\rho\) bounded away from zero.
This follows from the dyadic cofactor caps and (DC3).
The count below may alternatively be used directly under
that positive normalized height hypothesis.

## Counting pairs with a small determinant radical

Let \(a_u,a_v\) denote the normalized radical caps of the
remaining \(u,v\), and put \(a=a_u+a_v\). Suppose the
normalized gcd radical exponent is at most \(s\), and
suppose \(\operatorname{rad}(d_{ij})\le n^\alpha\).
For pairs with the preceding positive lower height bound, the
pair count has upper normalized exponent
\[
                  \frac{2a+\alpha+\min(s,a)}{3}
                                                               \tag{DC4}
\]
up to arbitrarily small fixed losses in the iterated limit.
If one instead assumes
\(\operatorname{rad}_{\not\mid D}(d_{ij})\le n^\alpha\),
replace \(\alpha\) in (DC4) by
\(\alpha+\log_n\operatorname{rad}(D)\). This accounts for
the fixed target support in the third term rather than
dropping its radical from the relative budget.

For the proof put \(\tilde u_i=u_i/g_u\) and likewise
\(\tilde u_j,\tilde v_i,\tilde v_j\). Let \(q_u\) be the
larger of the two normalized radical exponents on the
reduced \(u\)'s, \(q_v\) the corresponding exponent on the
reduced \(v\)'s, and set \(q=q_u+q_v\le a\).
Group these exponents into bins of small fixed width.
The new primitive triple has absolute radical cap
\(n^{2q+\alpha}\) up to the bin losses. Group its absolute
height exponent into bins as well. Since the lower endpoints
of the height bins are bounded below by a positive constant,
the resulting *relative* budgets are \(O(L)\). The
unrestricted small-budget \(1/3\) bound counts the primitive
triples at exponent \((2q+\alpha)/3\) with arbitrarily small
losses: take the relative budget just above the absolute
radical exponent divided by the lower bound for the
height, then take the bin widths small. All budgets and
the finite number of bins are fixed in each inner limit.
[The limiting small-budget bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

Each triple specifies the reduced cross-products up to bounded
ordering choices. Splitting them into the individual reduced
cofactors costs divisor-many choices. For fixed split data,
the parts of \(g_u,g_v\) supported on their respective reduced
cofactor pairs are enumerable at subpower cost. The parts
outside those supports have radical caps \(n^{a_u-q_u}\)
and \(n^{a_v-q_v}\), because each outside prime is added to
both corresponding original factors. Their combined radical
cap is also \(n^s\). Thus the scales cost exponent at most
\(\min(s,a-q)\). The maximum of
\[
                      \frac{2q+\alpha}{3}+\min(s,a-q)
\]
over \(0\le q\le a\) occurs at \(q=a-\min(s,a)\).
This proves (DC4) and accounts for all lifts to the original
ordered distinct targets.

## Three orientations of the same quotient determinant

Write the three residual letters in a row as \(u_i,v_i,w_i\),
with \(T_i=Dw_i\) and \(D\) full on its prime support.
For two rows let \(g_w=\gcd(w_i,w_j)=g_t/D\), and set
\[
              H_{ij}=\frac{|\Delta_{ij}|}{Dg_ug_vg_w}.
\]
This is a positive integer by (DC1). Let \(a_U,a_V,a_T\)
be the residual radical caps after recording \(F,G,D\),
and \(b_U,b_V,b_T\) the logarithms of the radicals of those
fixed factors, all in base \(n\). For any two of the residual
letters one can form the primitive difference of their
cross-products after dividing by the two same-letter gcds.
The positive difference terms, using respectively \(UV,UT,VT\)
as the crossed products, are exactly
\[
                      Dg_w H_{ij},\qquad Gg_v H_{ij},
                                   \qquad Fg_u H_{ij}.
                                                               \tag{DC6}
\]
For example \(u_i w_j-u_j w_i=\tau G\Delta_{ij}/D\), and
the cross-product gcd for this orientation is \(g_u g_w\).
The same identity with \(v\) proves the third expression.

Assume a positive lower bound on the absolute height exponent
of the reduced products in the chosen orientation. On pairs
with combined residual gcd radical mass at most \(s\), the
count with \(\operatorname{rad}(H_{ij})\le n^\alpha\) is at most
\[
        n^{(2(a_A+a_B)+\alpha+b_C+2s)/3+o(1)}                \tag{DC7}
\]
when \(A,B\) are the crossed letters and \(C\) the omitted
letter. The \(o(1)\) allows arbitrarily small exponent
losses in the stated iterated limits. In fact the omitted
fixed factor costs \(b_C\) and its residual gcd costs at most
\(s\) in the radical of the corresponding term of (DC6).
The reconstruction cost of the two crossed gcds is bounded
by the proof of (DC4). Splitting the reduced crossed products
and fixing those two gcds recovers two complete letters,
hence the original triples including the omitted letter.

The positive-height condition for each orientation holds on
almost all pairs after support peeling in any fibre of
exponent bounded below with one fixed factor of positive
numerical weight. To check this, if the two residual numerical
heights tend to zero and the positive fixed factor is on the
omitted letter, the signed ratio of the two complete terms
modulo that factor directly determines their residual ratio.
The latter is reduced and has numerator-denominator product
\(X^{o(1)}\), so rational reconstruction applies. If the
positive fixed factor is on a crossed letter, apply the
diffuse or atomic argument of the preceding positive-height
section to the omitted letter. In the diffuse case a cheap
full block on the omitted letter fixes the primitive ratio of
the two small residuals. In the atomic case the crossed full
factor gives a modulus for the other two letters, and the
small residual of the other crossed letter joins the wide
factor in the amplified bounded-root count. In every case
this contradicts a positive fibre exponent. The gcd numerical
bound in (DC3), including the residual target gcd, then
preserves positive height in the reduced products.

Consequently a fibre of exponent \(r>0\), after arbitrarily
small thinning as above, has almost all pairs satisfying
\[
  \log_n\operatorname{rad}(H_{ij})\ge
       \max\!\left\{ r,\ \max_{C}
               \bigl(6r-2(a_A+a_B)-b_C\bigr)\right\}
                           -\text{allowance}.                \tag{DC8}
\]
The \(r\) bound already follows by fixing an anchor, recording
the three gcds as its divisors, and counting the positive
integer \(H_{ij}\) by its radical. The other bounds follow
from (DC7) and comparison to the number of ordered pairs.
In particular the lower bound concerns a quotient determinant
with the entire common factor of each of the three original
letters removed. It can use any of the three residual caps
when the pooled records shrink them unevenly.

## Consequence near a saturated affordable matching

Suppose the initial \(U,V\) full root costs are \(f,g\) and
\(t=f+g\), and a sparse matching in that record has exponent
\(r=1/3-t\) at the endpoint. The initial residual cap
satisfies \(a\le2/3-t\), and further support locks only
decrease \(a\). After arbitrarily small exponent loss and
support normalization, the pair count above shows that,
for every fixed positive margin, almost every pair has
\[
      \log_n\operatorname{rad}(d_{ij})
                   \ \ge\ \frac23-4t-\text{margin}            \tag{DC5}
\]
whenever the right side is positive. Indeed the remaining
set has at least \(n^{2r-o(1)}\) ordered pairs up to the
chosen thinning allowance, which may be made smaller
than the strict deficit in (DC4) using \(s=\eta\).
The same argument with \(a\) after peeling gives the
threshold \(6r-2a\), up to the thinning allowances.
The fixed-anchor bound independently forces an almost-all
lower threshold \(r+\log_n\operatorname{rad}(D)\) for
the full radical, or \(r\) outside the fixed target support.
The two-pair estimate improves the basic anchor threshold
\(r\) already for \(t<1/9\). At \(t=1/24\) its threshold
is \(1/2\), whereas the basic anchor threshold is \(7/24\).

The normalized determinant is invariant under the fixed
support expansions. Thus a uniform bound strictly below
the right side of (DC5) on all pairs in such an initial
matching excludes its saturation. One concrete way a
bound of this type can enter is through a fixed target
divisor close to the projective reconstruction scale:
if for each fixed \(L\)
\[
 E_u+E_v-\log_X D\le \zeta L+o_X(1),
\]
then the outside radical exponent is at most \(\zeta\);
the full exponent in (DC4) also charges the radical of
\(D\). A fixed positive numerical slack alone supplies no
small radical bound. Full support peeling controls shared
parts, but still allows the remaining determinants to
carry large radicals on new primes. No global upper
bound on that large-radical case is proved here.

