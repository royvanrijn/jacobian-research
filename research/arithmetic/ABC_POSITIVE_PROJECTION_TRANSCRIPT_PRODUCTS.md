# Positive projection quotients and transcript products

A quotient of projected roots at two increasing orders has
nonnegative prime valuations. Two complementary ratios of these
quotients then carry three possible directions at a prime. A
product of two such records cannot cancel a prime between the
three directions. On balanced abc boxes this gives essentially
diagonal multiplicative energy at an arbitrarily small prescribed
power of the radical modulus. The proof uses actual addition
to recover the remaining atomic factors after the prime support
has been exposed.

## Bounded lists of short positive quotients

Use a height-first sequence of balanced original valuation boxes:
\(n=X^L\), at each fixed \(L>0\) take heights to infinity, and then
let \(L\to0\). The complete letters are \(A+B=C_0\) in the original
height shell with total radical cap \(n^{1+o(1)}\). The three
limiting degree measures have mass \(1/3\). Write \(c_A,c_B,c_C\)
for their nonatomic masses and \(S=c_A+c_B+c_C\). The pairwise
nonatomic mass condition gives \(1/2\le S\le1\). Here \(c_C\)
refers to the measure on \(C_0\).

Fix a prescribed pool \(P\) of polynomial prime product containing
2 and the fixed exceptional factors of the chosen boxes. Impose
that its contribution to the original normalized radical tends
to zero. For a prescribed pool of polynomial prime product one
can enforce this in a full-exponent near-endpoint sequence by
the fixed-support discount (EB5), taking height first.
For \(W=A,B,C_0\) let
\[
 x_j(W)=\prod_{p\notin P}p^{\lfloor v_p(W)/j\rfloor}.
\]
With \(1\le q<q+h\) put
\[
 M_q(W)=\frac{x_q(W)}{x_{q+h}(W)},\qquad
 \phi_q(e)=\lfloor e/q\rfloor-\lfloor e/(q+h)\rfloor .
\]
The \(M_q(W)\) are positive integers and their prime valuations
on their original terms are \(\phi_q(e)\ge0\).

**PP1.** Fix \(b,\zeta>0\). Along each limiting profile there
is a list of at most a fixed number of pairs \((q,q+h)\),
depending on the box and the two tolerances, with each order
\(O(L^{-\alpha})\) for any chosen \(1/2<\alpha<1\), such that

- at all selected probes,
  \(M_q(A)M_q(B)M_q(C_0)\le n^b\) for every row in the box,
  up to an arbitrarily small inner exponent allowance
- the union of their supports contains at least \(S-\zeta\)
  normalized mass on compact nonatomic carriers of the
  original letters, up to the outer box allowances.

Choose the carriers in a single positive compact degree band
\([r,R]\) and away from the countable union of the atoms of
the three limiting measures. Approximate them by continuity
neighborhoods in the exact slot boxing. Use the sum of the
three slot measures as the full weights \(w_e\), and the
carrier weights \(u_e\le w_e\). This makes the limiting
carrier marginal nonatomic, even in the presence of atoms
in the full weights.

To prove (PP1), repeat the bounded-probe argument of (MR4)
with \(Q=\lceil L^{-\alpha}\rceil\),
\(h=\lfloor \xi L Q^2\rfloor\) and \(q\) uniform on
\([Q,2Q)\). The fixed \(\xi>0\) is chosen small below.
Replace its two mask weights by
\[
      G_q=\sum_e u_e\,\mathbf 1_{\phi_q(e)>0},\qquad
      F_q=\sum_e w_e\,\phi_q(e).
\]
In a fixed compact degree range \(0\le Le\le R'\), choose
\(\xi R'<1\). The floors differ by zero or one, so the hit
condition is the same phase condition as for a parity contrast
in (MR4). On the carrier each hit probability is
\(\gg_r \xi\) and \(\mathbb E G_q\ll_R \xi\). The weighted
mixed compact expectation is \(O(\xi^2RR')+o(1)\).
Here the same rational lines in the Fourier calculation
have zero limiting mass against the nonatomic carrier.
The partner range \(Lf<r'\) contributes
\(O(\xi r')+o(1)\) by the crossing bound in (MR4).

For \(Lf>R'\) the exponent is bounded by \(\xi Lf+1\).
The uniform first moment and tail mass bounds therefore make
the part of \(F_q\) from this range \(O(\xi+1/R')\),
uniformly in \(q\). Its mixed contribution is
\(O_R(\xi^2+\xi/R')\). This argument also works with exact
slot boxing only up to a fixed cutoff at each \(L\), by
using this uniform tail bound beyond the cutoff. Choose
\(r'\) small and \(R'\) large, then \(\xi\) small.
The bound for \(\mathbb E(F_q G_q)\) can thereby be made
arbitrarily small relative to \(b\xi\). Discard indices
with \(F_q>b/2\), using a uniform upper estimate for
unboxed tail contributions and absorbing inner allowances.
The argument of (MR4) then leaves all but arbitrarily small
carrier mass on slots with a fixed positive good hit
frequency. Choose a bounded number of good probes to
cover this mass, by random sampling and averaging over
the miss probabilities. The order lists and any dyadic
refinements are fixed at each budget before height tends to
infinity. At a chosen probe \(F_q\) is exactly an upper
bound for the normalized logarithm of the three mask
product. This proves (PP1).
[Phase and mixed tail estimates (MR4)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md#a-bounded-number-of-small-modulus-probes)

## Two-row products have small fibres

Fix \(\delta>0\) and a positive integer weight list
\(w_1,\ldots,w_m\). Use (PP1) with
\(b=\delta/(4\sum_i w_i)\). Choose a common modulus
\(D>n^\delta\) supported on the prescribed pool, for instance
the first power of two exceeding that bound. Compute for each
row the vector
\[
 {\cal T}(A,B,C_0)=
 \left(\frac{M_q(A)}{M_q(C_0)},\
       \frac{M_q(B)}{M_q(C_0)}\right)_{q}
                      \pmod D.                             \tag{PP2a}
\]
When \(D\) is a power of two all coordinates are units.
Each ratio is read from the two corresponding projected
ratios at \(q,q+h\), and a common unit multiplier of
those projected residues cancels. This works separately on
each of the two complementary ratio lists. We will count
original row tuples with a specified product
\(\prod_i {\cal T}(A_i,B_i,C_i)^{w_i}\).

Within each coordinate of that product the exact rational
numerator times denominator is at most \(n^{\delta/4+o(1)}\).
Any two such positive rational numbers with the same
reduction modulo \(D\) have determinant of absolute value
less than \(D\), and so coincide. Thus the actual rational
product vector and its prime support are fixed by the
modular product record.

Choose cheap full divisors at fixed positive numerical
weight on every original letter with positive nonatomic
mass. At least two letters qualify. This choice uses the
compact block construction of (MD3). For any fixed
small \(\eta>0\) the joint dictionary of these divisors
can be required to have exponent less than \(\eta\),
with every numerical lower weight positive after fixing
\(\eta\) and the profile. For convenience count inside a
chart where these divisors are fixed. Pigeonholing the
small dictionary preserves the full endpoint exponent
with arbitrarily small fixed loss. All counts below
can also be summed over those charts with that small cost.

**PP2.** For \(m=2\), each fixed product in (PP2) has at most
\(O(n^\eta X^\varepsilon)\) ordered original pairs on a
fixed chart, at every sufficiently late budget. Here the
profile allowances and \(\zeta\) are chosen sufficiently
small in terms of \(\eta\). This is uniform in the product.

To see this, a prime appears on at most one original letter
of each row. At any probe its contribution to the two
valuations of (PP2) is a nonnegative multiple of one of
\[
                         (1,0),\qquad(0,1),\qquad(-1,-1).
\]
No sum of at most two nonzero positive multiples of these
vectors can vanish. The exact product support therefore
contains every prime detected by the probes on both
original rows. In particular it contains their chosen
nonatomic carriers except for the allowed missed mass.
Now apply the atomic completion part of (MD3) to each
row under any of its fixed available full divisors on the
opposite letter. Fix the two free letters' exact parts
on the supplied pool at subpower cost, freeze at small
price what remains outside the atom windows, and release
the pooled atomic groups by the original ratio modulo the
full divisor. The modular product thereby gives nearly
subpower completion of each row. Summing the inner
allowances proves (PP2). The argument includes repeated
rows, shared primes and varying complete targets.
[Atomic completion (MD3)](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)

For any subset of \(T\) rows in the chart, (PP2) for
weights \(1,1\) gives at least \(T^2 n^{-\eta}X^{-\varepsilon}\)
distinct product records, up to a constant. The
multiplicative energy counted on ordered original rows is
at most \(O(T^2 n^\eta X^\varepsilon)\), since the
sum of the squares of the product multiplicities is at
most their maximum times \(T^2\). A full-exponent subset
has essentially quadratic image at this small modulus
even if its original primes are correlated.

## Fixed positive products of more rows

**PP3.** For \(m\ge3\), at a fixed product and inside the same
kind of full-divisor chart, the ordered tuple count is
\(O(n^{\beta_m+\eta}X^\varepsilon)\), where
\[
 p_j=\frac{\binom{m-j}{2}}{\binom{m-1}{2}},\qquad
 \beta_m=\max_{1/2\le S\le1}\frac13
             \sum_{j=1}^m \min\{1-S+Sp_j,\,2Sp_j\}.
                                                               \tag{PP3}
\]
In particular \(\beta_3=1/3\), \(\beta_4=1/2\) and
\(\beta_m\le m(m-2)/(3(2m-3))\).
The choices of tolerances precede the fixed-budget height
limits. The estimate includes repeated rows, counted with
their slots, and arbitrary fixed positive weights.

Let \(P_0\) be the union support of the rational product
vector already recovered. At a detected prime outside
\(P_0\) all nonzero contributions at that probe sum to zero.
With positive weights, cancellation must use all three
directions displayed above. Thus the prime occurs in at
least three distinct row slots. For row \(i\) let \(\mu_i\)
be the normalized mass of its detected carrier in the
union of \(P_0\) and the raw prime supports of the
predecessors of that slot. Take a uniformly random ordering
and condition on the row occupying position \(j\). At least
two other slots contain every detected prime outside \(P_0\).
The chance that all of those slots occur after position \(j\)
is at most \(p_j\), by choosing any two of them. Thus the
expected known mass at position \(j\) is at least
\(S(1-p_j)-O(\zeta)\). The last two positions see every
detected prime for any order.

Given the previous rows, the prefix pool is fixed of
polynomial prime product. The zero-price effective radical
count bounds the number of possible next rows by
\(n^{(1-\mu_i)/3+\eta'}X^\varepsilon\), with an arbitrarily
small fixed allowance \(\eta'>0\). There is a second bound.
On any two free letters the already known prime support
can be fixed and the unknown diffuse residue frozen at
its combined raw cost. Complete the atomic groups by the
original ratio modulo the fixed full divisor on the opposite
letter as in (MD3). The two cheapest remaining diffuse
masses sum to at most \(2(S-\mu_i)/3+O(\zeta)\).
If one letter has zero nonatomic mass, a cheapest pair
may be chosen containing it, so the lock lies on one of
the two positive-mass letters. The same bound applies
after absorbing arbitrarily small costs outside the finite
atom windows. The last two rows consequently cost only
\(n^{O(\zeta)+\eta'}X^\varepsilon\) each.
[Fixed-pool effective radical discount (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The smaller exponent is the concave decreasing function
\(\mu\mapsto\min\{(1-\mu)/3,2(S-\mu)/3\}\) up to the
chosen allowances. Jensen's inequality bounds the expected
sum of these costs in the random order by
\(\frac13\sum_j\min\{1-S+Sp_j,2Sp_j\}+O_m(\zeta)\).
Thus each qualifying tuple has an order with at most this
total cost. Take a union over the finitely many orders and
refine the overlaps by a fixed finite mesh in the sequential
count. This gives (PP3) with its stated allowance. The
values for \(m=3,4\) follow by substitution. For the simpler
general bound, average \(p_j\) over the first \(m-2\)
positions. Its mean is \(m/(3(m-2))\), by summing the
binomial coefficients. Concavity in \(p_j\) and intersection
of the two affine functions of \(S\) give
\(m(m-2)/(3(2m-3))\). In particular a
positive product of three rows has exponent at most
\(1/3+\eta\) at a fixed record. The alternative arbitrary
signed products do not have the three-direction
cancellation property used here.

These bounds constrain the modular images of actual
balanced abc rows. To turn them into an unrestricted
one-third saving would require an independent upper bound
on the two-row product image, or an arithmetic lower bound
exceeding the essentially diagonal multiplicative energy
allowed by this count. The three-direction mechanism
alone supplies the small fibres and product growth.

## Separated integer weights

**PP4.** For any fixed number \(m\ge2\) of slots with weights
\(w_i=3^{i-1}\), a bounded list as in (PP1) gives
\(O(n^\eta X^\varepsilon)\) completions per fixed weighted
product on the same kind of full-divisor chart. In particular
its \(m\)-fold image has size at least \(T^m n^{-\eta}X^{-\varepsilon}\)
up to a constant, for any subset of \(T\) rows on the chart.

Fix a large compact cutoff \(H\) before choosing \(\xi\) in
(PP1). Enumerate the union of the primes from slots with
\(Lv_p(W)>H\) on the whole \(m\)-tuple, at raw cost
\(O(m/H)+o(1)\) by the numerical height bound.
Then fix the parts of every row on that union at subpower
cost. Choose \(\xi H<1\). Outside the union all mask
exponents are zero or one. For each coordinate and prime
the products supply a sum of the form
\(\sum_i3^{i-1}a_i\) with \(a_i\in\{-1,0,1\}\).
The signed digits are uniquely determined: in a comparison
of two different digit vectors the leading nonzero
difference has magnitude greater than the maximum
total magnitude of all the earlier differences. Hence
the fixed rational products recover every active prime
coordinate of every row outside the enumerated tail.
In particular each carrier is supplied to the atomic
completion. Take \(H\) large enough and all other
allowances small enough relative to \(\eta\). The
modular-to-rational step is unchanged because (PP1)
already prices the full masks, including the tail.
