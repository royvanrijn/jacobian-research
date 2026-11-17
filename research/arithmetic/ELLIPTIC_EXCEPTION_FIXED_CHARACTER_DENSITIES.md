# Arithmetic decay of each fixed coefficient character

Every fixed character of sufficiently large auxiliary order factors
through the rational reduction image at a small, precisely bounded
proportion of primes. For a character of prime order `ell`, that
proportion is `ell^(-2)+O(ell^(-3))`. The proof uses the actual full
torsion/division image and a finite Burnside count. It controls each
fixed Fourier mode arithmetically, but does not bound a Fourier tail
whose orders grow with the residue-prime cutoff.

## 1. The simultaneous affine image and the fixed-character event

Retain the fixed non-CM elliptic curve, the abstract free basis
`P_1,P_2` of `Gamma=E(Q)` modulo its order-two torsion, and
`Pi=widehat(Gamma)`. At a good prime let

\[
 \phi_p:\Pi\longrightarrow R_p,\qquad r_p=|R_p|,
 \qquad A_p=\{\phi_p=O\}.
\]

For a finite-order character `chi` of `Pi`, say that it factors
through `R_p` when it is trivial on `ker(phi_p)`. This is an exact
condition on the actual rational reduction image, not just on its
order or on the full group `E(F_p)`.

There is a fixed finite prime set `Sigma`, containing two and three,
for which every integer `m` supported outside `Sigma` has full
simultaneous affine image

\[
 \operatorname{Gal}\bigl(\mathbb Q(E[m],[m]^{-1}P_1,
                        [m]^{-1}P_2)/\mathbb Q\bigr)
       =E[m]^2\rtimes\operatorname{GL}_2(\mathbb Z/m).          \tag{1.1}
\]

Here is the extension of
[EOI1--EOI3](ELLIPTIC_EXCEPTION_ADELIC_OPEN_IMAGE.md) from prime
powers to arbitrary such `m`. Serre's open-image theorem supplies
the whole torsion factors outside a fixed finite set. The retained
uniform Kummer theorem bounds the translation-kernel index by a
constant `c` at every `m`. That index divides `m^4`. Include in
`Sigma` all primes at most `c`: any nontrivial index would have a
prime divisor larger than `c`, contradicting the bound. Thus the
translation kernel and torsion projection are both full, proving
(1.1). For any fixed finite set of auxiliary primes, inverse limits
give the full product of their affine Tate-module groups.

This uses precisely the Serre and Corvaja--Demeio--Javanpeykar--
Lombardo--Zannier interfaces retained and directly read in EOI.
No numerical exceptional set or computed Mordell--Weil basis is
asserted. In particular `Sigma` is fixed before the character order
is chosen.

Write

\[
 J_2(n)=n^2\prod_{\ell\mid n}(1-\ell^{-2}).
\]

There are exactly `J_2(n)` characters of exact order `n` of
`Gamma/nGamma=(Z/n)^2` when `n` is supported outside `Sigma`.

**Theorem EFC1 (fixed-order character density).** For every such
`n>=2` and every character `chi` of exact order `n`, the set of
good primes where `chi` factors through `R_p` has a natural density
`delta_chi` among primes. All exact-order-`n` characters have the
same density, and

\[
                   \boxed{0\le\delta_\chi\le J_2(n)^{-1}.}      \tag{1.2}
\]

The same density and bounds hold conditional on any fixed nonempty
finite Galois Frobenius condition, after enlarging `Sigma` for that
fixed condition. No uniformity as the Galois condition varies is
claimed.

**Existence of the density.** Put `S=supp(n)` and use the full
affine group from (1.1) at all powers of the primes in `S`. Its
Haar element has linear part `A=(A_ell)` and two independent
translation coordinates `t_1,t_2`, with linear and translation
coordinates mutually independent. Away from
`det(A_ell-1)=0`, the finite ambient primary group is

\[
 C_A=\prod_{\ell\in S}T_\ell(E)/(A_\ell-1)T_\ell(E).
                                                               \tag{1.3}
\]

Pink's Frobenius formula, retained as CDJLZ Theorem 6.7 in
[EAM1](ELLIPTIC_EXCEPTION_AFFINE_COSET_POLE_MEANS.md), identifies
the two rational primary reductions with
`(A_ell-1)^(-1)t_i mod T_ell(E)`. Their generated group `R` is
the `S`-primary part of `R_p`; the rational two-torsion contributes
nothing here.

On this nonsingular matrix locus, each of the two torsion targets
is locally constant: inverse matrices are continuous, and passage
to `V_ell/T_ell` has discrete target. Consequently their generated
group and the kernel of the coefficient map to it are locally
constant. Thus the condition that `chi` factor is locally constant.
The exceptional determinant-zero locus has Haar measure zero.
Indeed a nonzero polynomial over `Q_ell` has a zero set of measure
zero on `Z_ell^d`, by induction and Fubini: outside the zero set of
its nonzero coefficient polynomial, a one-variable fibre has
finitely many roots.

For completeness this only invokes finite Chebotarev: remove the
sets `v_ell(det(A_ell-1))>=b` for `ell in S`. Their measures
decrease to zero. On the remaining compact open set the event is
determined by one finite quotient, by a finite subcover of its
local constancy neighbourhoods. Finite-quotient Chebotarev gives
its density. The discarded set is also a finite-quotient event
and bounds the discrepancy above and below by a quantity tending
to zero. This proves existence of `delta_chi`, equal to the Haar
probability. The finite Chebotarev input is Milne ANT300 Theorem
8.31, as retained in EAM and EOI.

**Coefficient symmetry and the bound.** Right multiplication of
the ordered translation pair by `GL_2(prod_(ell in S) Z_ell)`
preserves Haar measure, commutes with the linear torsion action,
and acts transitively on the `J_2(n)` characters of exact order
`n`. Hence these characters have equal probabilities of factoring
through the generated group `R`.

For a finite abelian group `B`, let `a_n(B)` count its elements of
exact order `n`. It also counts its characters of exact order `n`.
Therefore the number of the coefficient characters that factor
through `R` is `a_n(R)`. Since `R` is a subgroup of `C_A`,

\[
                         a_n(R)\le a_n(C_A).                  \tag{1.4}
\]

This compares element counts. It does not assert that a subgroup
character extends to the ambient group with the same order.

For every `d|n`, Smith normal form gives

\[
 |C_A[d]|=\left|\ker(A-1:(\mathbb Z/d)^2\to(\mathbb Z/d)^2)\right|.
\]

The uniform group `GL_2(Z/d)` has one orbit on vectors of each
additive order `e|d`: at a prime power the common coordinate
valuation determines the orbit, and the Chinese remainder theorem
combines them. Its number of orbits is consequently `tau(d)`.
Counting pairs `(g,v)` with `gv=v`, first by vectors and then by
group elements, is Burnside's elementary identity and gives

\[
                    \mathbb E_A|C_A[d]|=\tau(d).               \tag{1.5}
\]

The determinant-zero null set does not change this finite-quotient
calculation. Mobius inversion now yields

\[
 \mathbb E_A a_n(C_A)
       =\sum_{d\mid n}\mu(n/d)\tau(d)=1.                      \tag{1.6}
\]

The last identity is multiplicative, and at a prime power it is
`(a+1)-a=1`. Average (1.4), and use coefficient symmetry, to get
`J_2(n)delta_chi=E a_n(R)<=1`, proving (1.2).

For the conditional version, let `L/Q` be the fixed finite Galois
field. Apply the same simultaneous affine argument over `L`,
enlarging `Sigma` using EOI and its Kummer constant over `L`.
At every relevant finite level the affine division field over
`Q` then has its full group even after base extension to `L`.
It is therefore linearly disjoint from `L`. Chebotarev in the
product, followed by the same zero-boundary approximation, makes
the auxiliary affine Haar distribution independent of the fixed
Frobenius condition. \(\square\)

## 2. A matching lower order at prime character orders

**Theorem EFC2 (sharp scale for a prime-order mode).** For every
prime `ell` outside `Sigma` and every character of exact order
`ell`,

\[
 \boxed{\frac{\ell-2}{\ell^3}\le\delta_\chi
                      \le\frac1{\ell^2-1}.}                   \tag{2.1}
\]

In particular `delta_chi=ell^(-2)+O(ell^(-3))`, with an absolute
implied constant. The conditional version of EFC1 has the same
lower bound after its fixed enlargement of `Sigma`.

**Proof.** Put `q=ell`. A uniform matrix in `GL_2(F_q)` has
distinct eigenvalues `1,lambda`, with `lambda!=1`, with probability

\[
                           \frac{q-2}{(q-1)^2}.               \tag{2.2}
\]

There are `q-2` possible `lambda`, and each conjugacy class has
centralizer of order `(q-1)^2`. Conditional on such a residue
matrix, `det(A-1)/q mod q` is a nonzero linear functional of the
first lift digit. Thus `v_q(det(A-1))=1` with probability
`1-1/q`. On these lifts (1.3) is the cyclic group `C_q`.
The two independent uniform translations generate a nonzero
subgroup with probability `1-q^(-2)`. Such a subgroup contributes
exactly `q-1` characters of exact order `q`. Dividing by the
`J_2(q)=q^2-1` equiprobable coefficient characters, this event alone
contributes

\[
 \frac{q-2}{(q-1)^2}\left(1-\frac1q\right)
       \left(1-\frac1{q^2}\right)\frac{q-1}{q^2-1}
                         =\frac{q-2}{q^3}.
\]

All omitted events contribute nonnegative amounts. Combine this
lower bound with EFC1. \(\square\)

For example, active-prime gating does not spoil these bounds.
The auxiliary multiplication levels are prime to three. Since
the curve is good at three, their division fields are unramified
there. By
[ETD2--ETD3](ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md),
the absolute density of primes at which `chi` factors and
`p=2 mod3`, `3` not dividing `#E(F_p)`, is `delta_chi/4`.
The corresponding density conditional on that gate is `delta_chi`.
This follows at every finite auxiliary level and then by the
zero-boundary limit above; it is not an independence assertion
about reductions at two varying residue primes.

## 3. What this proves about the Fourier spectrum

**Corollary EFC3 (square-summable fixed density spectrum).** The
nontrivial characters whose orders are supported outside `Sigma`
satisfy

\[
 \boxed{\sum_\chi\delta_\chi^2
       \le\sum_{\substack{n\ge2\\\operatorname{supp}(n)\cap\Sigma=\varnothing}}
                            \frac1{J_2(n)}<\infty.}           \tag{3.1}
\]

Define the reweighted pole sum and its exact mean by

\[
 T_Y=\sum_{\substack{p\le Y\\p\text{ good}}}
                  \frac{r_p}{p}\mathbf1_{A_p},\qquad
 H_Y=\mathbb E T_Y=\sum_{\substack{p\le Y\\p\text{ good}}}\frac1p.
\]

For every fixed such character,

\[
 \boxed{\frac{\mathbb E(T_Y\overline\chi)}{H_Y}
                                     \longrightarrow\delta_\chi.} \tag{3.2}
\]

Thus every fixed finite collection of these Fourier modes has
limiting squared energy bounded by its corresponding sum in
(3.1). The same assertions hold when the residue primes are
restricted by any fixed finite Frobenius condition, with the
enlargement of `Sigma` from EFC1 and its corresponding exact mean.

**Proof.** There are `J_2(n)` characters of each exact order `n`.
Equation (1.2) therefore gives the first inequality. Moreover
`J_2(n)>=n^2/zeta(2)`, so the series converges.

Character orthogonality on the pole subgroup gives exactly

\[
 \mathbb E(T_Y\overline\chi)
       =\sum_{\substack{p\le Y\\\chi\text{ factors through }R_p}}\frac1p.
                                                               \tag{3.3}
\]

Natural relative prime density `delta_chi` implies the same ratio
for these reciprocal-prime sums: partial summation of the
counting-function difference bounds the tail error by an arbitrary
fixed multiple of `H_Y`, while its initial finite contribution is
`o(H_Y)`. The denominator diverges, also under a positive-density
fixed Frobenius condition by the argument in EAM3. This proves
(3.2). Summing its squared limits over a fixed finite family proves
the last assertion. \(\square\)

The weights `r_p/p` are bounded by two at sufficiently large good
primes. These are still pole counts, not the formal-depth hazards
of ETP. No information about `e_p` is asserted.

## 4. An equivalent-Haar weak limit, without uniform integrability

Let `G_good=prod_(ell notin Sigma) Z_ell^2`, with normalized Haar
measure `m_good`, and let `pi_good:Pi->G_good` forget the exceptional
primary coordinates and rational torsion.

**Theorem EFC4 (the averaged pole measures have a positive weak
limit).** The probability measures

\[
 \nu_Y=(\pi_{\rm good})_*
               \left(\frac{T_Y}{H_Y}\,d\mathbb P\right)
\]

converge weakly on `G_good` to `g dm_good`, where

\[
 g>0\quad\text{almost everywhere},\qquad
 \int g\,dm_{\rm good}=1,\qquad g\in L^2(m_{\rm good}).        \tag{4.1}
\]

Thus the weak limit is equivalent to Haar measure. This is a
statement about averaged measures, not `L^2` convergence or uniform
integrability of the finite-time pole counts. The same assertion
holds with a fixed nonempty finite Frobenius condition and its
enlarged exceptional set from EFC1.

**Local densities.** At a good auxiliary prime `ell`, let a Haar
affine element give the finite group `R` generated by its two
primary reductions. For `x in Z_ell^2`, define

\[
 g_\ell(x)=\mathbb E_{A,t_1,t_2}
                 \left(|R|\,\mathbf1_{\{x\mapsto0\text{ in }R\}}\right).
                                                               \tag{4.2}
\]

Tonelli's theorem gives integral one: the coefficient map to `R`
is onto, so its kernel has Haar measure `1/|R|`. Its Fourier
coefficient at any local character is exactly its factoring
density. Parseval and EFC1 consequently give

\[
 \|g_\ell-1\|_2^2
       \le\sum_{a\ge1}\frac1{J_2(\ell^a)}
       =\frac{\ell^{-2}}{(1-\ell^{-2})^2}=:v_\ell.             \tag{4.3}
\]

One may justify the use of Parseval without assuming this density
is initially square integrable: its conditional expectations on
the finite coefficient quotients have exactly the corresponding
finite Fourier coefficients and uniformly bounded squared norms.
They are `L^2`-Cauchy by the summable coefficient bound, and their
limit equals the integral-one density (4.2), since all finite
quotient integrals agree. Thus (4.3) holds for that actual density.

There is also a pointwise lower bound, independent of translations:
when `det(A-1)` is a unit, `C_A=R=0` and the integrand in (4.2)
is identically one. Hence

\[
 g_\ell\ge h_\ell
    :=1-\frac{\ell^2-2}{(\ell-1)^2(\ell+1)}\ge\frac12
                     \qquad(\ell\ge3).                       \tag{4.4}
\]

For the displayed probability, a matrix with a nonzero fixed
vector has `ell-1` such vectors unless it is the identity, which
has `ell^2-1`. Burnside's mean number of nonzero fixed vectors
is one. Solving for the probability of any nonzero fixed vector
gives `(ell^2-2)/((ell-1)^2(ell+1))`, proving (4.4).

**The infinite product is positive.** The functions
`g_ell(x_ell)` are independent on `G_good`, have mean one, and
satisfy `sum_ell v_ell<infinity`. Their finite products are
`L^2`-Cauchy: independence gives squared norms
`prod_ell(1+||g_ell-1||_2^2)`, and the difference of two nested
products has squared norm equal to the earlier norm squared times
the later norm-square product minus one. These tails tend to zero.
Denote the resulting `L^2` limit by `g`.

To check positivity, for `z>=1/2` use

\[
 |\log z|\le2|z-1|,\qquad
 |\log z-(z-1)|\le2(z-1)^2.
\]

Equations (4.3)--(4.4) make the means of `log g_ell` absolutely
summable and their variances summable. The independent centered
logarithms therefore have an almost surely convergent series.
Here is an elementary justification of this last step. For
independent mean-zero variables, the probability that a partial
sum of a finite tail ever exceeds `epsilon` in absolute value is
at most `epsilon^(-2)` times the tail variance. Partition by the
first crossing; the future increment has zero conditional mean,
so the second moment at the final index bounds the crossing sum.
Choose tail indices with total variance at most `2^(-3j)` and
threshold `2^(-j)`. The resulting exceptional probabilities have
finite sum. The union bound on their tails proves the Cauchy
property almost surely. Adding the convergent mean series proves
convergence of the full logarithm series.

Its exponential is strictly positive and is also the almost sure
product limit. It agrees with the `L^2` limit by uniqueness of
limits in probability. This proves (4.1), including integral one
by `L^2` convergence.

**Identification of the weak limit.** Full simultaneous affine
Haar measure is a product over the auxiliary primes. A character
of composite order factors through the rational primary image
exactly when each of its primary components does. Consequently
its factoring density is the product of the local densities.
The Fourier coefficients of the product `g` are therefore exactly
the `delta_chi` of EFC1. Equation (3.2) proves convergence of each
Fourier coefficient of `nu_Y` to that of `g dm_good`.

Every continuous function on this compact profinite group is
uniformly approximable by a function on a finite quotient, by
uniform continuity and a finite open cover. Such a finite-quotient
function has a finite character expansion. Since all measures
have mass one, the fixed-mode convergence therefore proves weak
convergence. The conditional prime-gate version has the same
full affine distribution after its fixed enlargement of `Sigma`,
so the argument is unchanged. \(\square\)

The square summability in (3.1) concerns the **limiting density of
each fixed mode**. It is not a bound on the complete finite-time
Fourier energy of `T_Y/H_Y`. Taking the limit for each character
and then summing cannot bound the opposite order of operations;
energy may move to orders growing with `Y`. In particular neither
`L^2` convergence of those counts nor a bounded second moment is
proved by EFC4.

Even an equivalent-Haar weak limit does not force recurrence.
In [EPS2](ELLIPTIC_EXCEPTION_FINITE_PRIMARY_SEPARATION.md), every
fixed nontrivial coefficient character factors through only
finitely many Fermat blocks. The normalized reweighted measures
therefore converge weakly to Haar itself, while almost every
vector still has only finitely many returns. Uniform integrability
or another bound on the actual large-quotient clustering remains
an additional arithmetic obligation. No elliptic recurrence,
zero all-primes measure or Hasse infinitude follows from EFC1--EFC4.

The new arithmetic statements use the already retained open-image,
Kummer, Pink and finite Chebotarev interfaces, plus the written
finite-group count above. No numerical exceptional threshold,
new finite checker or formal verification is claimed.
