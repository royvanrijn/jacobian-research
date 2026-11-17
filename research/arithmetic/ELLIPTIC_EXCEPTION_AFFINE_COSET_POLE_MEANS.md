# Positive-density pole opportunities in every rational congruence coset

Every fixed affine congruence coset of the elliptic rational group
has a positive-density set of active primes where its reduction
contains the identity. At these primes its conditional pole-return
probability is at least `1/(2p)`. This includes cosets translated by
the rational two-torsion point. Consequently every fixed feasible
finite-place prefix has a divergent expected number of unweighted
pole returns. No almost-sure or weighted-hazard recurrence follows
from this mean statement alone.

## 1. Finite-level prescription away from three

Retain the fixed elliptic curve `E`, the rational group
`Gamma=E(Q)` of rank two and torsion `<T>` of order two, and
`Pi=widehat(Gamma)`. Give `Pi` normalized Haar probability `P`.
Choose a full integral free basis `P_1,P_2` for the proof only.
Existence of that basis is not a claim that displayed points are
already saturated.

The curve has `End_(Qbar)(E)=Z`, as proved from its nonintegral
`j` invariant in
[EKX1](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md).
Independence of the basis therefore makes
`Z=(P_1,P_2)` nondegenerate in `E^2`: its cyclic group is Zariski
dense, by the independence and Poincare-reducibility argument given
there.

Use `T_ell(E)` for the Tate module, `V_ell=T_ell(E) tensor Q_ell`,
and `pi_ell:V_ell/T_ell(E)->E[ell^infinity]`. A finite good-reduction
point has its usual `ell`-primary component.

**Lemma EAM1 (bounded torsion targets at fixed auxiliary primes).**
Let `L` be any finite set of primes different from three. For
`ell in L` and `i=1,2`, fix a torsion target
`zeta_(i,ell) in E[ell^infinity]`. There are a finite Galois
extension `C/Q`, unramified at three, and a nonempty conjugacy class
`Omega subset Gal(C/Q)` with this property: at every sufficiently
large good prime whose Frobenius lies in `Omega`, the two basis
points have these `ell`-primary targets after a suitable choice
of the place used to identify torsion under reduction. Invariant
integral relations involving the targets and rational torsion
therefore hold independently of that choice.

**Kummer and Frobenius inputs.** Corvaja--Demeio--Javanpeykar--
Lombardo--Zannier, *On the distribution of rational points on ramified
covers of abelian varieties*, Compositio Mathematica **158** (2022),
2109--2155, Theorem 6.9, printed p.2144 (PDF p.37), gives a constant
`c` such that for every positive integer `m`,

\[
 N_m=\operatorname{Gal}
  \bigl(\mathbb Q(E[m],[m]^{-1}P_1,[m]^{-1}P_2)/
                                     \mathbb Q(E[m])\bigr)
                 \subseteq E[m]^2
\]

has index at most `c`. Divisibility of the points is allowed by
Remark 6.10. The finite affine-action convention is equation (13),
pp.2142--2143. Theorem 6.7 and Remark 6.8 on p.2144 give the
Frobenius reduction formula

\[
 [\bar P_i]_\ell
      =\pi_\ell\bigl((A_\ell-1)^{-1}t_{i,\ell}\bigr),          \tag{1.1}
\]

where `A_ell` is the torsion action and `t_(i,ell)` its Kummer
translation. The operator `A_ell-1` is invertible over `Q_ell`
at a good Frobenius prime different from `ell`:
[published primary](https://doi.org/10.1112/S0010437X22007746),
[retained PDF](sources/cdjlz-2022-ramified-covers.pdf).

**Proof of the lemma.** Choose an integer `c_0>=c` and put
`c_*=c_0!`. The finite quotient `E[m]^2/N_m` has order at most
`c_0`, so

\[
                         c_*E[m]^2\subseteq N_m               \tag{1.2}
\]

for every `m`. This is what permits the small auxiliary primes,
including two; no assumption `ell>c` is made.

Let `s_ell` kill both targets, so they have order dividing
`ell^(s_ell)`, and put `a_ell=v_ell(c_*)`. Choose a good base
prime `p_0` outside `L union {3}` splitting the finite field

\[
 \mathbb Q\bigl(E[\prod_{\ell\in L}\ell^{a_\ell+s_\ell}],
       [\prod_{\ell\in L}\ell^{a_\ell}]^{-1}P_1,
       [\prod_{\ell\in L}\ell^{a_\ell}]^{-1}P_2\bigr).
\]

Such primes exist by Chebotarev. Take a compatible arithmetic
Frobenius `h` there. At every `ell in L`, its linear action
`A_ell` is the identity modulo `ell^(a_ell+s_ell)` and its two
translation coordinates belong to `ell^(a_ell)T_ell(E)`.
The matrices `A_ell-1` are invertible by the Frobenius input above.

Lift each target to `z_(i,ell) in V_ell`. The desired translation

\[
 t^*_{i,\ell}=(A_\ell-1)z_{i,\ell}
                               \in\ell^{a_\ell}T_\ell(E)      \tag{1.3}
\]

is integral. Its difference from the translation of `h` is in
the same subgroup. Take
`m=prod_(ell in L) ell^(n_ell)` with all `n_ell` sufficiently
large. By (1.2), a translation in `N_m` adjusts both coordinates
to (1.3) modulo `ell^(n_ell)`, while leaving the linear actions
unchanged modulo those powers. Let `g` be the resulting element
of the finite Galois group of

\[
 C=\mathbb Q(E[m],[m]^{-1}P_1,[m]^{-1}P_2).
\]

Here is the finite precision needed in this adjustment. Write
`d_ell=v_ell(det(A_ell-1))`. Choose
`n_ell>d_ell+s_ell+a_ell+2`. If a lift has actions `A'_ell,t'_i`
congruent to the chosen ones, then

\[
 t'_i-(A'_\ell-1)z_{i,\ell}
                         \in\ell^{n_\ell-s_\ell}T_\ell(E).
\]

The determinant valuation remains `d_ell`, and the inverse matrix
has entries of valuation at least `-d_ell`. Thus (1.1) is exactly
the prescribed target modulo `T_ell(E)`. This proves that the
target condition is determined by the finite element `g`.

Take `Omega` to be its conjugacy class. At a prime in that class
one can choose a place where Frobenius is `g`. Conjugating changes
the targets by the same torsion action, which fixes rational torsion
and preserves their integral relations. Finally `3` does not divide
`m`, and `E` has good reduction at three. The field `C` is
unramified there by finite etale multiplication and properness,
specifically Milne, *Abelian Varieties* (January 2, 2022),
Proposition 20.7 and Corollary 20.8, printed/PDF p.44:
[retained primary](sources/milne-abelian-varieties-2022.pdf).
If `L` is empty the lemma is immediate with `C=Q`.
\(\square\)

## 2. Pole opportunities in every affine coset

**Theorem EAM2 (positive-density affine-coset pole opportunities).**
For every `gamma in Gamma` and positive integer `M`, there is a
positive-density Chebotarev set `P_(gamma,M)` of good active primes
`p=2 mod3` such that

\[
 \boxed{\bar\gamma\in[M]R_p.}                       \tag{2.1}
\]

For `C_gamma,M=gamma+[M]Pi`, its exact conditional pole probability
at every sufficiently large prime of this set is

\[
 \boxed{\mathbb P(\bar P_p=O\mid C_{\gamma,M})
                   =\frac1{|[M]R_p|}\ge\frac1{2p}.}            \tag{2.2}
\]

In particular the expected number of unweighted pole returns,
conditional on that fixed coset, diverges.

**Proof.** It suffices to obtain (2.1) for an even multiple of
`M`, since its coset is contained in the original one. We may
therefore assume `M` is even during the construction. Write

\[
 \gamma=aP_1+bP_2+cT,\qquad c\in\{0,1\},\qquad s=v_2(M).
\]

At every odd prime `ell|M`, `ell!=3`, prescribe zero targets for
both basis points. At two choose the following targets, always
with `zeta_(2,2)=0`:

| Case | Target `zeta=zeta_(1,2)` | Consequence |
| --- | --- | --- |
| `c=0` | `zeta=0` | the two-primary part of `gamma` is zero |
| `c=1`, `t=v_2(a)<s` | `2^t zeta=T`, of order `2^(t+1)` | `a zeta+T=0` |
| `c=1`, `v_2(a)>=s`, including `a=0` | `2^s zeta=T`, of order `2^(s+1)` | `a zeta+T` belongs to `<T>=[M]<zeta>` |

These torsion targets exist in `E[2^infinity]`. In the second
row the odd part of `a/2^t` acts trivially on the order-two point
`T`. In the third row `M zeta=T`, since `M/2^s` is odd.
Thus every row makes the two-primary part of `gamma` divisible
by `M` within the subgroup generated by the targets and `T`.
This explicitly handles a nonzero torsion translate.

Apply EAM1 to all primes dividing `M` except three. Its field
`C` is unramified at three. The disjointness and joint-density
theorems [ETD2--ETD3](ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md)
allow this finite Frobenius condition simultaneously with

\[
 p=2\pmod3,\qquad3\nmid N_p.
\]

Explicitly, `Gal(Q(E[3])/Q)=GL_2(F_3)`, and its twelve matrices
of determinant `-1` with no eigenvalue one have proportion `1/4`.
They impose precisely these two conditions and are independent of
the field `C`. The resulting Chebotarev set has positive density
`|Omega|/(4[C:Q])`, after removing finitely many primes.

At an odd prime dividing `M` other than three, the rational image
has trivial primary part by the zero targets. Its three-primary
part is trivial because `3` does not divide `N_p`. The target
table makes `gamma` divisible by `M` in its two-primary image.
At primes not dividing `M`, multiplication by `M` is an
automorphism of the corresponding finite primary group. These
observations prove (2.1) in every Sylow subgroup, hence in `R_p`.
The integral relations in the table are unchanged by the common
Galois action on the targets.

Reduction pushes the normalized Haar measure of `gamma+[M]Pi`
uniformly onto the coset `bar gamma+[M]R_p`. Once that coset
contains the identity its probability is exactly `1/|[M]R_p|`.
Since `|[M]R_p|<=N_p<2p` at large good primes, (2.2) follows
also for the original `M`. A positive-density prime set has
divergent reciprocal-prime sum, proving the last assertion.
\(\square\)

## 3. Every fixed feasible prefix has a divergent unweighted mean

**Corollary EAM3 (fixed-prefix pole means).** Let
`C_0={alpha=a} cap intersection_(q in S) G_q` be any positive
finite-place event. There are a constant `c_(S,a)>0` and a
positive-density Chebotarev set of active primes such that

\[
 \mathbb P(\bar P_p=O\mid C_0)\ge c_{S,a}/p.         \tag{3.1}
\]

Consequently, for all sufficiently large `Y`,

\[
 \mathbb E\!\left(\sum_{\substack{p\le Y\\p\text{ active}}}
             1_{\{\bar P_p=O\}}\,\middle|\,C_0\right)
                          \gg_{S,a}\log\log Y.                \tag{3.2}
\]

**Proof.** By local openness and density, the full profinite
preimage of `C_0` contains a coset `C_gamma,M`, chosen with even
`M` if desired. This is the clopen coset construction of
[EOP2](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md), not a closure
of an arbitrary rational subset. For the prime set in EAM2,

\[
 \mathbb P(\bar P_p=O\mid C_0)
 \ge\frac{\mathbb P(C_{\gamma,M})}{\mathbb P(C_0)}
                 \frac1{2p}.
\]

The ratio is a fixed positive number. The Chebotarev set has
positive natural density among primes. Its counting function is
therefore bounded below by a positive multiple of `Y/log Y`, using
the elementary prime-count lower bound proved in
[ESR2](ELLIPTIC_EXCEPTION_SMALL_REDUCTION_IMAGES.md). Partial
summation gives a reciprocal-prime sum bounded below by a positive
multiple of `log log Y`, proving (3.2).
\(\square\)

The coset, the prime set and the constant depend on the fixed
prefix. This is not a uniform return bound conditioned on the
growing soluble prefix in ETP2. An infinite expected count does
not by itself imply infinitely many returns almost surely.

The same contained-coset argument works after intersecting the
prefix with any positive-measure clopen cell. Thus it settles the
divergent unweighted-mean hypothesis in the uniform-cell criterion
of ECM section 4; the uniform second-moment bound remains open.

Moreover these are **unweighted** pole means. At the new primes
the actual formal depth `e_p` is not controlled. The hazards `b_p`
can be small when `3|e_p`, so divergence of their weighted mean
has not been proved. The variance and recurrence obligations of
[ECM3](ELLIPTIC_EXCEPTION_COMMON_QUOTIENT_MOMENTS.md) and
[ETP5](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md) remain separate.
No value of the all-primes measure or Hasse infinitude follows.

## 4. The remaining depth condition for a weighted first moment

**Conditional corollary EAM4 (weighted means under a depth-density
hypothesis).** Suppose the active primes with `3|e_p` have natural
density zero among all primes. Then for every fixed feasible prefix
`C_0`,

\[
 \mathbb E\!\left(\sum_{\substack{p\le Y\\p\text{ active}}}
                         h_p\,\middle|\,C_0\right)
                          \gg_{S,a}\log\log Y.                \tag{4.1}
\]

Here `e_p` and the nonnegative hazards `h_p` are those of
[ETP2--ETP5](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md).
Set `h_p=0` at the finitely many primes preceding that sieve's
cutoff.
The density hypothesis is unproved. It is weaker than the
conditional density-one assertion `e_p=1` in ETP6.

**Proof.** Delete the density-zero exceptional set from EAM3's
positive-density prime set. The remaining set still has positive
natural density and the same lower-order reciprocal-prime bound.
For `3` not dividing `e_p`, the pole hazard is

\[
 b_p=1-\frac{p^{2-\delta_p}}{p^2+p+1},\qquad
 \delta_p=3\lceil e_p/3\rceil-e_p\in\{1,2\},
\]

so `b_p>=1/2`. Since `h_p>=b_p 1_{\{\bar P_p=O\}}`,
(3.1) and summation prove (4.1). This gives no second-moment or
almost-sure recurrence estimate. \(\square\)

This is a written deduction from the stated Kummer, Frobenius,
disjointness and Haar inputs. The deep Kummer theorem is imported,
not replayed. No computed free basis, effective Kummer constant,
numerical division field, new finite checker, formal verification
or external review is claimed.
