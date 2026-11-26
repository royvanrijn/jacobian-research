# Recovering common factors of splices from primitive outputs

Dividing a nonprimitive additive relation by its common gcd
does not require a separate charge for the gcd's numerical
height. The primes still present in the primitive relation
are known there. A prime that disappears from all three terms
can be recorded at a cost exactly paid by the reduction in
its three original radical contributions at the one-third
counting coefficient. This gives a joint radical bound for
actual near splices at every strict subterminal depth.

[Near splices of genuine rows (SRS1–SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md)
· [Splice gcds and near-splice radical bounds (HPS2–HPS4)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md)

## Primitive outputs and the exact radical discount

Fix \(L>0\) and put \(n=X^L\). Let \(P\) be a prescribed
positive integer of polynomial height, and write
\(\operatorname{rad}_{P^c}\) for the radical outside its
prime support. Consider a signed additive relation on
three nonzero integers of absolute value at most \(X^H\),
where \(H\) is fixed. After taking absolute values and
possibly permuting the terms, write it as
\(U+V=W\) with \(U,V,W>0\).
Put
\[
     G=\gcd(U,V,W),\qquad
     (u,v,w)=(U,V,W)/G,\qquad
     R=\operatorname{rad}_{P^c}(G).
\]
Divide \(R\) into two squarefree factors
\[
 R_{\rm vis}=\gcd(R,\operatorname{rad}(uvw)),
 \qquad R_{\rm miss}=R/R_{\rm vis}.
\]
Use normalized masses
\(z=\log_n R_{\rm vis}\), \(b=\log_n R_{\rm miss}\), and set
\[
   \mathcal R=\log_n\operatorname{rad}_{P^c}U
         +\log_n\operatorname{rad}_{P^c}V
         +\log_n\operatorname{rad}_{P^c}W .
\]
The exact identity is
\[
       \mathcal R-2z
             =\log_n\operatorname{rad}_{P^c}(uvw)+3b .
 \tag{ESL1}
\]
Indeed a prime outside \(P\) that divides two original terms
divides all three and belongs to \(R\).
A prime outside \(R\) has the same contribution before and after
division. Each prime in \(R_{\rm vis}\) loses two of its three
contributions, since the normalized terms are pairwise coprime.
Each prime in \(R_{\rm miss}\) loses all three.

Given the primitive output and the integer \(R_{\rm miss}\),
the visible part has divisor-many possible values, since it
divides the radical of \(uvw\). The whole common factor \(G\)
is then supported on the primes of \(P R_{\rm miss}R_{\rm vis}\).
The uniform fixed-support bound enumerates all \(G\le X^H\)
on this support at height-subpower cost. In particular this
counts arbitrary powers of both visible and missing primes.
No restriction on the powerful part of \(G\) is needed to
restore the output. A known prescribed divisor of \(G\)
may simply be required in this enumeration.
[Uniform fixed-support enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

**ESL2.** Suppose the normalized triples have numerical
height \(w\ge X^{h_0}\) for some fixed \(h_0>0\), and raw
radical at most \(n^{B_{\rm raw}}\) for a fixed
\(B_{\rm raw}>0\). For every fixed \(B,\eta>0\) and
all sufficiently small fixed \(L\), the number of
nonprimitive signed additive outputs with
\(\mathcal R-2z\le B\) is
\[
                         O(n^{B/3+\eta} X^\varepsilon).
 \tag{ESL2}
\]
The constants may depend on the fixed cap data, the polynomial
height bound for \(P\), \(L,\eta,\varepsilon\), uniformly in \(P\).
One can infer the raw cap from the outside cap if
\(\operatorname{rad}(P)\le n^{O(1)}\).
If each nonprimitive output represents at most
\(n^{e+o(1)}\) original configurations, add \(e\) to the
count exponent. A list of fixed pool choices is priced
separately.

Indeed (ESL1) gives \(b\le B/3\).
Bin \(b\) on a fixed mesh of width \(\tau\).
For a bin \([b_-,b_+]\) the normalized effective radical
cap outside \(P\) is at most \(n^{B-3b_-}\).
Apply (EB4) to the normalized primitive triples, pricing
the primes of \(P\) at zero. Their raw radical is bounded.
To pass to a fixed relative budget, split the height
range \([X^{h_0},X^H]\) into a fixed fine mesh of
logarithmic exponents and use dyadic shells within it.
The relative budgets range over a fixed bounded
multiple of \(L\) and can be rounded upwards on a finite
mesh at the fixed \(L\). This gives exponent
\((B-3b_-)/3\) in base \(n\), with an arbitrarily small
normalized allowance. The finite sign and ordering
choices are harmless. There are at most \(n^{b_+}\)
possible missing radicals. The visible part and all the
powers of \(G\) cost height-subpower by the preceding
argument. The total exponent is at most
\(B/3+b_+-b_-\) up to these allowances. Choose the
mesh and inner losses inside \(\eta\) and \(X^\varepsilon\).
[Effective one-third count on prescribed pools (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

### Two-letter costs after restoring the common factor

Write \(\mathcal R_T=\log_n\operatorname{rad}_{P^c}(T)\)
on each of the three physical terms, and
\(z_T=\log_n\gcd(R,\operatorname{rad}(T/G))\).
These three visible contributions sum to \(z\) by
primitivity. Define
\[
          c_T=\mathcal R_T-z+z_T.
\]
The primitive outside radical mass on this role is
exactly \(c_T-b\), since every prime of \(R\) that
does not survive there must be removed from
the original radical of \(T\).

Under the height and raw cap hypotheses of (ESL2),
a bound \(c_T+c_U\le B_2\) for any specified pair
of roles gives the nonprimitive output count
\[
                          O(n^{B_2/2+\eta}X^\varepsilon).
 \tag{ESL2a}
\]
One can use the minimum over the three choices.
The multiplicity factor of (ESL2) applies here as well.

For this assertion fix the full \(P\)-supported part
of each primitive term at height-subpower cost.
These parts can be of positive numerical height.
They enter the fixed coefficients, with the primitive
complete terms unchanged. On the remaining parts
use singleton anatomic layers with the outside raw
radical caps. The initial exceptional factors can
be fixed at an arbitrarily small normalized loss,
before taking the inner height limit. The arbitrary
coefficient estimate (MCB2) applies since the complete
height is at least \(X^{h_0}\), hence is at least
\(n^C\) for its prescribed threshold \(C\) at sufficiently
small fixed \(L\). For a bin \(b\in[b_-,b_+]\)
the sum of the two selected primitive caps is at most
\(B_2-2b_-\), so the count before restoring \(G\) costs
at most \(B_2/2-b_-\), up to the counting allowance.
Add the missing radical price \(b_+\) as before.
The bound \(b\le B_2/2\) makes the mesh finite.
Take the mesh and other errors inside \(\eta\).
This proves (ESL2a).
[Arbitrary-coefficient half bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

The coefficient \(3\) of \(b\) in (ESL1) pays precisely
for restoring a missing radical at the effective one-third
count. The two-letter bound pays the same record because
both charged letters lost that radical. Visible factors
reduce the compensated cost without a power-sized record.

## Application to near splices at all strict subterminal depths

Work with positive primitive rows \(A+B=C\) and the original
three radical caps. Fix a marked computed role \(W\),
a compared complete role \(Z\), and one of the splice
orientations with
\[
       W^*>0,\qquad \Delta=W^*-W_j=Z_i-Z_j\ne0 .
\]
The signed identity uses the compared role \(Z=A\) in a
sum splice of \(A_i+B_j\), or \(Z=C\) for the difference
splice \(C_i-B_j\) with \(W=A\). Other orientations are
handled by fixed changes of roles and signs. Here \(j\)
is the original anchor row and \(i\) the other row.
A known shallow divisor \(M\) divides
\(G=\gcd(W^*,W_j)=\gcd(W_j,\Delta)\).
The pool \(P\) contains the support of the fixed
full locks, hence those of \(M\).
Suppose
\[
 \log_n\operatorname{rad}(W_j)\le r_W+o(1),\qquad
 \log_n\operatorname{rad}_{P^c}(W_j)\le r_W-\mu_W+o(1).
\]
In the lock setting \(\mu_W\) can be taken to be
the normalized radical mass of the marked full lock.
Write
\[
    b_\Sigma=\log_n\operatorname{rad}_{P^c}(W^*)
                       +\log_n\operatorname{rad}_{P^c}(|\Delta|),
 \qquad
    z=\log_n\gcd\!\left(\operatorname{rad}_{P^c}G,\
       \operatorname{rad}\!\left(
                  \frac{W_j}{G}\frac{W^*}{G}\frac{|\Delta|}{G}\right)
                         \right).
\]

**ESL3.** For the splice pairs with \(b_\Sigma-2z\le B\),
positive lower numerical height after division by \(G\),
and bounded normalized raw radical of the normalized
triple, the count is at most
\[
        O\!\left(n^{(r_W-\mu_W+B)/3+u_W+u_Z+\eta}
                    X^\varepsilon\right).
 \tag{ESL3}
\]
Here \(u_W,u_Z\) are upper exponents for the maximum
complete degrees on \(W,Z\), and \(\eta>0\) is any
fixed counting allowance at sufficiently small fixed budgets.
The constants have the same uniformity as (ESL2).
The fixed pool's bounded raw mass suffices for the
raw output hypothesis on this subfamily.

Indeed restoration produces
the values \(W_j,W^*,|\Delta|\), with at most a bounded
choice of sign. Choose the anchor row among the
\(n^{u_W+o(1)}\) rows with that \(W_j\).
The identity for \(\Delta\) recovers the compared
label \(Z_i\) exactly from \(Z_j\).
Then there are at most \(n^{u_Z+o(1)}\)
possibilities for the second row. Apply (ESL2) with
total compensated cap \(r_W-\mu_W+B+o(1)\).
If \(i,j\) belong to the same complete \(W\)-fibre,
the second row is determined by \((W_j,Z_i)\).
Thus the second multiplicity factor can be dropped
for such a ripple.

In the fine-layer complete-matching setting of (HPS3),
let the combined depth be \(s<d=1/3\).
Assume the original compared marked value satisfies
\(W_j\ge X^{h_0}\) for a fixed \(h_0>0\), and
the outside primes of the marked role meet the
small-individual-prime hypothesis of (HPS3).
The fixed support pool has bounded normalized
radical mass, as for the union of the three locks.
Off a power-small fraction of the actual near splices,
that theorem gives
\[
      \frac{G/M}{\operatorname{rad}_{P^c}(G/M)}
                                  \le n^\zeta
\]
with \(\zeta>0\) arbitrarily small and fixed.
Since \(M\) is shallow and
\(\operatorname{rad}_{P^c}(G/M)\mid\operatorname{rad}(W_j)\),
we have \(G\le n^{\log_n M+r_W+\zeta+o(1)}\).
Hence \(W_j/G\ge X^{h_0/2}\) at sufficiently small fixed
budgets. This use of the gcd estimate is only for
normalized height. The possibly large squarefree gcd
does not cost its mass in (ESL3).

The full actual population has exponent \(2d-s\)
by (SRS2). Therefore outside a power-small fraction,
for every fixed \(\eta>0\) with positive threshold,
\[
              b_\Sigma-2z
                    >3(2d-s)-r_W+\mu_W-\eta .
 \tag{ESL4}
\]
At balanced caps and half depth \(s=d/2\) this gives
\(7/6+\mu_W-\eta\) for the compensated two-output
radical cost. It holds at every fixed \(s<d\), including
above half depth, with the corresponding threshold.
It does not require the first valuation layer of
the marked letter to have vanishing mass. For a more
general splice population of exponent \(t\) for which
the normalized height hypothesis has been proved,
the threshold supplied by (ESL3) is
\(3(t-u_W-u_Z)-r_W+\mu_W-\eta\).
Use \(u_W\) alone in place of \(u_W+u_Z\)
for a same-complete-target ripple. A size hypothesis
for the normalized outputs remains necessary in that
variant, where (HPS3) is not automatic.

There is a stronger two-letter condition from (ESL2a).
For the three terms \(W_j,W^*,|\Delta|\), write
\(z_W,z_*,z_\Delta\) for the masses of \(R_{\rm vis}\)
surviving in the respective primitive terms.
The smaller of the two compensated prices
\(c_W+c_*\) and \(c_W+c_\Delta\) is bounded above by their
average, giving the following count with conditional
degree exponent \(e\):
\[
   O\!\left(n^{(r_W-\mu_W)/2+B_2/4+e+\eta}
                 X^\varepsilon\right)
       \quad\text{if}\quad b_\Sigma-3z+z_W\le B_2.
 \tag{ESL4a}
\]
Here \(e=u_W+u_Z\), or \(e=u_W\) for a same-target
ripple. To check the price, note \(c_W=\mathcal R_W-z+z_W\)
and \(c_*+c_\Delta=b_\Sigma-z-z_W\).
Thus \((2c_W+c_*+c_\Delta)/4\) is at most
\((r_W-\mu_W)/2+B_2/4+o(1)\).
Bin the prices and select the cheaper pair in each bin
before applying (ESL2a), then restore the original rows
as in (ESL3).

On the complete matching this gives, off a power-small
exception, the stronger compensated joint floor
\[
     b_\Sigma-3z+z_W
          > 4(2d-s)-2(r_W-\mu_W)-\eta.
 \tag{ESL4b}
\]
For a general permitted population of exponent \(t\),
replace \(2d-s\) by \(t-e\).
At balanced half depth (ESL4b) gives
\(4/3+2\mu_W-\eta\) on the left, and it is a stronger
joint floor than (ESL4) at every \(s<d\) with \(r_W=d\).
Indeed the baseline difference is
\(d-s+\mu_W\), and \(3z-z_W\ge2z\).
The same two pair tests individually give
\[
 \begin{split}
 \log_n\operatorname{rad}_{P^c}(W^*)-2z+z_W+z_*
                       &>2(t-e)-r_W+\mu_W-\eta,\\
 \log_n\operatorname{rad}_{P^c}(|\Delta|)-2z+z_W+z_\Delta
                       &>2(t-e)-r_W+\mu_W-\eta .
 \end{split}
\]
All these bounds retain the primitive normalized height
condition. HPS3 supplies it on the stated actual splice
populations.

For a complete matching one can also combine the
second individual floor with the independent floor
for the positive splice from (SRS2). Let \(r_U,r_V\)
be its two original operand caps and \(\mu_U,\mu_V\)
the outside discounts from their recorded locks.
Both exceptional sets are power-small at each fixed
\(s<d\). On the remaining pairs, since
\(2z-z_W-z_\Delta=z+z_*\),
\[
       b_\Sigma-z-z_*
          >5(2d-s)-r_U-r_V-r_W+\mu_U+\mu_V+\mu_W-\eta.
\]
For balanced caps at half depth the right side is
\(3/2+\mu_U+\mu_V+\mu_W-\eta\).
This estimate combines the actual three-letter counts
for the original operand splice and for its normalized
comparison against a genuine target.

[Actual near-splice population (SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md#near-splices-in-both-residue-and-numerical-scales)
· [Subterminal gcd excess (HPS3)](ABC_HYBRID_CELL_PRIVATE_SUPPORT_AND_SPLICE_GCD.md#gcds-of-a-splice-with-its-compared-true-target)

## Physical differences of two original rows

For two distinct genuine primitive rows in the original
shell put
\[
 \Delta_T=T(y)-T(x),\qquad
 G=\gcd(|\Delta_A|,|\Delta_B|),\qquad
 u_T'=\Delta_T/G .
\]
Assume all three differences are nonzero. They satisfy
the original signed additive law before and after the
division. Apply the preceding definitions with
\[
  \mathcal R=\sum_{T=A,B,C}\log_n\operatorname{rad}_{P^c}
                                      (|\Delta_T|),
  \quad z=\log_n\gcd\!\left(\operatorname{rad}_{P^c}G,\
                    \operatorname{rad}(|u_A'u_B'u_C'|)\right).
\]
Thus (ESL2) counts the actual signed physical
displacements with cap \(\mathcal R-2z\le B\)
and a tall normalized direction, with no cap on
the powerful part or numerical size of \(G\).

**ESL5.** Suppose the original balanced chart has
\(N\ge n^{d-\kappa}\) rows with \(d=1/3\)
and sufficiently small cap tolerance \(\kappa\),
and on \(C\) each complete value occurs at most
\(n^{u_C+o(1)}\) times. Let \(P\) be fixed with bounded
normalized radical mass. For each fixed \(B,\eta>0\),
the number of ordered pairs with nonzero differences
and \(\mathcal R-2z\le B\) is bounded by
\[
 O\!\left(N^2n^{-\delta}+
        n^{\,B/3+(d+\kappa)/2+u_C+\eta} X^\varepsilon\right)
 \tag{ESL5}
\]
for a fixed positive \(\delta\) at sufficiently small
fixed budgets. The other fixed parameters determine
the permitted cap tolerance and budget threshold.
There is a strict saving on such a physical budget
when \(B<3/2-3u_C\) by a fixed margin and the
allowances are small enough.

The short normalized directions have a power-small pair
count by (SC1)–(SC2), keeping all three original caps.
Also discard pairs with
\(\max(C(x),C(y))/\gcd(C(x),C(y))\le X^\theta\)
for a sufficiently small fixed \(\theta>0\).
For fixed \(x\), enumerate the scalar gcd
\(a\mid C(x)\) at divisor cost and use the
integer feature \(C(y)/a\) in (CF1). This gives
a uniform fixed power loss after summing the anchors.
The raw normalized direction radical is bounded
because (ESL1) and the fixed pool bound give
a bound on its outside and inside masses.

For a fixed nonzero displacement, its scalar target
difference \(\Delta_C\) is fixed. Choose
\(a\mid|\Delta_C|\) at divisor cost and apply (MCB2)
to the two scalar quotients and the fixed shift
\(|\Delta_C|/a\). Their maximum has a positive
numerical height floor. The fixed shift has zero
dictionary cost. Both moving radical caps are at
most \(n^{d+\kappa}\). Anatomic factorization and
singleton boxes on the moving terms give
\(n^{(d+\kappa)/2+o(1)}\) possibilities for the complete
target pair. A choice of the first target leaves at
most \(n^{u_C+o(1)}\) first rows. The displacement
then determines the second row. Apply (ESL2) with
this multiplicity bound to prove (ESL5).
[Short primitive secants (SC1–SC2)](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#a-uniform-bound-for-queries-short-in-two-directions)
· [Small features on one original role (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)
· [Fixed-coefficient half bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

The direction cap in (ESL5) concerns the physical
differences with the stated visible-factor credit.
A bound only on the primitive direction radical
does not bound the exhausted mass \(b\) in (ESL1).
That distinction is relevant when comparing with
(NGT3) and with the one-coordinate primitive
direction exclusion of (UC2).
[Normalized gcd tails and direction counts (NGT1–NGT3)](ABC_SECANT_NORMALIZED_GCD_TAIL.md)
· [Primitive direction with an unrecorded cofactor (UC2)](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md#primitive-secant-coordinates-with-an-arbitrary-scaling-gcd)

## Limit of the transfer

These implications supply upper counts conditional on
caps for the two splice outputs or for the three physical
differences. The known shallow lifting identities force
neither of these caps on a power-dense configuration.
Polynomial transfer estimates ordinarily price the
radicals of the evaluated numerators and denominators.
Identity (ESL1) gives an exact way to price a polynomial
transfer when its physical signed output has a large
common factor, but it still needs a source of small
compensated radical mass.

The weighted same-target collision estimates in the
external abc bounds use exact common-power parts of
two fixed-shift terms after recording complementary
valuations. Browning–Verzobio's coefficient-uniform
generalized Fermat estimate requires prescribed
powers on the variable roots. The small-Pell
contribution for a balanced affine line uses a
prescribed square right-hand side of height one in
its norm equation. None of those hypotheses is
supplied by the shallow residue algebra or by the
support restoration identity alone.
[Zeng, complete modulus locking, §§3–5](https://hexagonmath.org/pdf/2611.00621v1)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)
· [Zeng, small Pell discriminants, Theorem 1.1 and §6](https://hexagonmath.org/pdf/2611.00611v1)
· [Radical thresholds for polynomial transfers](ABC_POLYNOMIAL_TRANSFER_RADICAL_THRESHOLD.md)
