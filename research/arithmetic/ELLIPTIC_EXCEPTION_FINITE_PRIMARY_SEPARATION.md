# Finite-primary separation and the remaining quotient clustering

Every fixed rational congruence cell has a refinement on which a
positive-density set of elliptic pole returns is independent of any
prescribed finite set of primary coefficient coordinates. This removes
those correlations exactly. It does not bound the remaining large
common quotients. An explicit divisibility model shows that the latter
gap persists even with large marginal images and divergent return
means in every fixed cell.

## 1. Separating any fixed finite set of primary coordinates

Retain `Gamma=E(Q)=Z P_1 directsum Z P_2 directsum <T>` and
`Pi=widehat(Gamma)` with normalized Haar probability. The basis is
abstract; no displayed points are asserted to be saturated. At a good
prime write `phi_p:Pi->R_p` for the actual rational reduction, and
`A_p={phi_p=O}` for the pole event.

**Theorem EPS1 (finite-primary separation inside an affine cell).**
Let `C=gamma+[M]Pi`, and let `L` be any finite set of rational
primes. There are a subcell

\[
 C'=\gamma'+H\subseteq C,\qquad H=[N]\Pi,
 \qquad N=2\operatorname{lcm}(2,M),
\]

and a positive-natural-density Chebotarev set `P` of good active
primes `p=2 mod3` such that, for every sufficiently large `p in P`,

\[
 \bar\gamma'\in[N]R_p,\qquad
 \gcd\left(|[N]R_p|,\prod_{\ell\in L\cup\operatorname{supp}(N)}
                                          \ell\right)=1.       \tag{1.1}
\]

Consequently the conditional pole probability on `C'` is
`1/|[N]R_p|>=1/(2p)`. All selected pole events, jointly, are
independent of the pro-`L` coefficient coordinates on `H`.
Every relative common quotient `D^H_(p,q)` of two selected primes
has order prime to all primes in `L union supp(N)`.

**Proof.** Put `M_0=lcm(2,M)`, so `N=2M_0` is even and
`[M_0]Pi subseteq[M]Pi`. Write

\[
 \gamma=aP_1+bP_2+cT,\qquad c\in\{0,1\}.
\]

When `c=0`, take `gamma'=gamma`. When `c=1`, choose
`k in {0,1}` and set `gamma'=gamma+[kM_0]P_1` so that its
first coefficient `a'=a+kM_0` satisfies
`t=v_2(a')<v_2(N)`. Such a choice exists: if
`v_2(a)<v_2(M_0)` either choice works; otherwise choose the
bit making `a'/2^(v_2(M_0))` odd. This also proves `C' subset C`.

Apply [EAM1](ELLIPTIC_EXCEPTION_AFFINE_COSET_POLE_MEANS.md) at
all primes in `L union supp(N)` other than three. At every odd
auxiliary prime prescribe zero primary targets for both basis
points. At two prescribe both targets zero if `c=0`. If `c=1`,
prescribe

\[
 \zeta_{1,2}=\zeta,\quad 2^t\zeta=T,\quad
 \operatorname{ord}(\zeta)=2^{t+1},\qquad \zeta_{2,2}=0.
\]

Then `a' zeta+T=0`, and `[N]zeta=0` because
`t+1<=v_2(N)`. Thus the two-primary component of the whole
reduced coset `gamma'+[N]R_p` is zero. Its odd auxiliary-primary
components are zero as well. These relations are unchanged by the
common Galois action on the targets; the rational point `T` is fixed.

EAM1 supplies a nonempty finite Frobenius class in a Galois field
unramified at three. Impose in addition `p=2 mod3` and
`3` not dividing `#E(F_p)` using the disjointness theorem
[ETD2--ETD3](ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md).
The resulting Chebotarev set has positive density. Hence the
three-primary component is also zero. This proves the coprimality
in (1.1). At a prime not dividing `N`, multiplication by `N` is
invertible on the corresponding primary component of `R_p`.
At every prime dividing `N`, the reduced coset has zero primary
component by construction. Therefore `bar gamma' in[N]R_p`.
The Haar reduction and Hasse-bound argument of EAM2 gives the
stated conditional probability.

Since `N` is even, `H` is torsionfree and isomorphic to
`widehat(Z)^2`. It is the direct product of its pro-`L` part
`H_L=prod_(ell in L) Z_ell^2` and the complementary primary
coordinates. Each selected image `phi_p(H)` has order prime to
`L`; its restriction to `H_L` is therefore zero. All selected
return indicators, after translating `C'` to `H`, depend only
on the complementary coordinates. This proves joint independence
of the entire pro-`L` part, not merely vanishing of pairwise
Fourier coefficients. Finally `D^H_(p,q)` is a quotient of both
selected images by [EVO2](ELLIPTIC_EXCEPTION_POLE_VARIANCE_OBSTRUCTION.md),
so it has the asserted coprimality. \(\square\)

In particular every nontrivial finite-order character of `H`
whose order is supported on `L` has zero Fourier coefficient in
every selected pole sum. Conditioning further on any finite
pro-`L` congruence event leaves its full joint pole distribution
unchanged. The cell and prime set depend on `L`; no assertion
about their density uniformly as `L` grows is made.

## 2. Large common quotients can still prevent recurrence

**Proposition EPS2 (a divisibility countermodel).** There is a
family of surjective finite quotient maps
`psi_p:widehat(Z)^2->G_p`, indexed by primes `p>=3`, with all
of the following properties:

* `sqrt(p)<|G_p|<=p` at every index;
* every fixed finite set of primary primes is absent from `|G_p|`
  for all sufficiently large indices;
* in every fixed positive-measure clopen cell, the pole probability
  is at least `1/p` for every sufficiently large index;
* nevertheless almost every coefficient vector has only finitely
  many pole returns.

The same conclusions hold on restricting the indices to active
primes `p=2 mod3`. Thus the listed properties, including the
fixed-cell divergent means, cannot by themselves prove recurrence.

**Construction and proof.** Put

\[
 F_k=2^{2^k}+1\quad(k\ge0).
\]

These integers are pairwise coprime, since
`prod_(j<k) F_j=F_k-2` and all `F_j` are odd. Also

\[
 F_{k+1}=(F_k-1)^2+1<F_k^2,\qquad
                         \sum_k F_k^{-1}<\infty.               \tag{2.1}
\]

For each prime `p>=3`, take the unique `k` with
`F_k<=p<F_(k+1)`, and define

\[
 G_p=\mathbb Z/F_k\mathbb Z,\qquad
 \psi_p(x,y)=x\bmod F_k.
\]

The size assertion follows from (2.1). Pairwise coprimality also
means that every fixed prime divides at most one `F_k`, proving
eventual separation from every finite primary set.

Any fixed clopen cell is measurable modulo some integer `M`.
For all sufficiently large `k`, `gcd(F_k,M)=1`. The Chinese
remainder theorem then makes `x mod F_k` uniform and independent
of both coefficient coordinates modulo `M`. Consequently the
conditional pole probability in that cell is exactly

\[
                         1/F_k\ge1/p.                         \tag{2.2}
\]

In particular its expected count diverges, with lower order
`log log Y` up to the prime cutoff `Y`. For all primes this
follows from the retained elementary prime-count argument in
[ESR2](ELLIPTIC_EXCEPTION_SMALL_REDUCTION_IMAGES.md) and partial
summation. For active primes the same argument uses their
positive Chebotarev density, as in EAM3.

However the event for an entire block of indices is the single
event `B_k={x=0 mod F_k}`, of Haar probability `1/F_k`.
For each `K`, the union bound gives

\[
 \mathbb P\left(\bigcup_{k\ge K}B_k\right)
                               \le\sum_{k\ge K}F_k^{-1}
                                      \longrightarrow0.
\]

Almost every vector therefore hits only finitely many blocks.
Each block contains finitely many primes, so its total pole-return
count is finite. This remains true conditional on every
positive-measure clopen cell.

For two indices in one block, the common quotient is the full
group `Z/F_k Z`; for indices in different blocks it is trivial,
by the same Chinese remainder theorem. Thus repetition of a
large common quotient, rather than any fixed small-primary
entanglement, accounts for the failure of recurrence. \(\square\)

EPS2 is an abstract finite-group construction, not an elliptic
curve or a counterexample to an arithmetic recurrence conjecture.
It shows why EPS1 and EAM's first moments still require a bound
on the remaining common-quotient clustering, such as an appropriate
conditioned moment estimate. No weighted-hazard recurrence,
all-primes Haar value or Hasse-infinitude conclusion is inferred.

EPS1 imports only the previously stated Kummer, disjointness and
Chebotarev interfaces. EPS2 uses elementary profinite quotients,
the Chinese remainder theorem and the union bound. No new finite
checker or formal verification is claimed.
