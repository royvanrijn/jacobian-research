# Exact common quotients and second moments for elliptic returns

The dependence between two reduction events is measured by their exact
common quotient. Reduction-order estimates control the individual
probabilities but do not control this quotient or the distribution of
an earlier soluble set inside it. This note gives finite formulas for
those correlations and an elementary sufficient criterion for the
weighted recurrence required by the all-primes sieve. The needed
asymptotic moment bounds remain open.

## 1. Exact finite correlations

Put `Pi=widehat(Gamma)`, with its normalized Haar probability `P`,
as in [ETP1](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md).
All statements below push forward to the actual rational closure.
At a good prime let

\[
 \phi_p:\Pi\longrightarrow R_p\subseteq E_0(\mathbb F_p),
 \qquad K_p=\ker\phi_p,\quad r_p=|R_p|.
\]

The group `R_p` is the actual rational reduction image. For a set
`B_p subset R_p`, write `A_p=phi_p^(-1)(B_p)`. The pole return is
the special case `B_p={O}`.

**Proposition ECM1 (the common-quotient formula).** For two primes,
put

\[
 D_{p,q}=\Pi/(K_p+K_q),\qquad d_{p,q}=|D_{p,q}|.
\]

For pole returns,

\[
 \boxed{\mathbb P(A_p\cap A_q)=\frac{d_{p,q}}{r_pr_q},
 \qquad
 \operatorname{Cov}(\mathbf1_{A_p},\mathbf1_{A_q})
                    =\frac{d_{p,q}-1}{r_pr_q}.}     \tag{1.1}
\]

More generally, let `psi:Pi->G` be an actual finite quotient, let
`C subset G` be an event of positive probability `mu`, and set

\[
 D=\Pi/(\ker\psi+K_p).
\]

For `d in D`, define

\[
 c_d=\mathbb P(\psi\in C\mid d),\qquad
 b_d=|B_p\cap(R_p)_d|,
\]

where `(R_p)_d` is the fibre of `R_p->D`. Then

\[
 \boxed{\mathbb P(A_p\mid\psi\in C)
                   =\frac{\sum_{d\in D}c_db_d}{r_p\mu}.}       \tag{1.2}
\]

In the pole case this is `c_0/(r_p mu)`.

**Proof.** The actual joint image in `G times R_p` is exactly the
fibre product over `D`. To verify this, choose lifts of a compatible
pair. Their difference lies in `ker(psi)+K_p`, so adjusting the
lifts in those two kernels gives a common lift. Thus, conditional
on `d`, the two coordinates are independent and uniform in their
respective fibres. Since `|D|` is finite and all its classes have
equal measure, integrating gives (1.2).

Taking the two pole events, or using the subgroup index identity,
gives

\[
 [\Pi:K_p\cap K_q]\,[\Pi:K_p+K_q]=r_pr_q.
\]

This proves (1.1). All groups involved in this calculation are
finite quotients, so the index identity is the elementary finite
abelian subgroup cardinality formula. \(\square\)

There is an equivalent Fourier description. Let `X_p` be the
characters of `Pi` factoring through `R_p`, and put

\[
 a_p(\chi)=\frac1{r_p}\sum_{b\in B_p}\overline{\chi(b)}.
\]

Character orthogonality gives

\[
 \mathbb P(A_p\cap A_q)
        =\sum_{\chi\in X_p\cap X_q}a_p(\chi)
                                      \overline{a_q(\chi)}.    \tag{1.3}
\]

The common character group has order `d_(p,q)`. In particular the
absolute covariance is at most
`(d_(p,q)-1)|B_p||B_q|/(r_p r_q)`. In the pole case all the
Fourier coefficients are `1/r_p`, recovering (1.1). The common
quotient and its characters can be computed from a certified actual
joint finite image: if that image is `J subset R_p times R_q`,
then quotient `R_p` by `{r:(r,0) in J}`. This needs no full
Mordell--Weil basis. The result is not determined by `r_p,r_q` alone.

## 2. Conditioning on a full finite-place prefix

A full local root condition need not be measurable in any one finite
quotient. Formula (1.2) must therefore not be applied by treating such
a prefix as a finite subset without further justification.

**Proposition ECM2 (exact prefix-weighted moments).** Fix a positive
finite-place event

\[
 C_0=\{\alpha=a\}\cap\bigcap_{q\in S_0}\mathscr G_q,
                 \qquad \mathbb P(C_0)>0.
\]

Let `I` be a finite set of additional primes and let each `f_p` be
a rational-valued function of the reduction in `R_p`. Given a
certified joint finite quotient for `S_0 union I` as in EQA2 and
FQS3, all conditional moments of these finitely many functions
are exact rational calculations.

Specifically, choose a certified open product `W` inside the actual
joint closure, small enough that the reductions at `I` are constant
on its cosets. Let `J=H/W` be that exact finite joint image. For
`j in J`, let `w_j` be the conditional Haar probability of `C_0`
on that coset, and write `r_(j,p)` for its reduction at `p in I`.
Then

\[
 \boxed{\mathbb E\!\left(\prod_{p\in I}f_p\mid C_0\right)
   =\frac{\sum_{j\in J}w_j\prod_{p\in I}f_p(r_{j,p})}
                {\sum_{j\in J}w_j}.}              \tag{2.1}
\]

**Proof.** Conditional Haar measure on a coset of the product `W`
is product Haar measure. Hence `w_j` is the class indicator times
the product of the conditional local soluble weights for `q in S_0`.
These weights are rational by
[EQA2](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md),
including their infinite periodic annuli. Each `f_p` is constant
on that coset. Both numerator and denominator before normalization
are sums with the common factor `1/|J|`; it cancels to give (2.1).
The finite image `J` is computed by
[FQS2--FQS3](ELLIPTIC_EXCEPTION_FINITE_QUOTIENT_SATURATION.md) under
their supplied-quotient hypotheses. \(\square\)

This applies to return indicators and to the exact hazards of ETP4,
when the fixed local numbers `e_p,k_p` are supplied. It also gives
their exact conditional covariances. The weights in (2.1) retain the
prefix; unconditional covariances such as (1.1) cannot replace them.

## 3. A sufficient second-moment criterion for the all-primes measure

Let `0<=Z_p<=1` be any fixed sequence of return weights, for example
the hazards `h_p` of ETP4. For a prime cutoff `Y`, put

\[
 T_Y=\sum_{p\le Y}Z_p,\qquad
 m_Y=\mathbb E(T_Y\mid C_0),\qquad
 v_Y=\operatorname{Var}(T_Y\mid C_0).
\]

Finitely many omitted initial primes have no effect on divergence.

**Theorem ECM3 (conditional recurrence criterion).** If

\[
 m_Y\longrightarrow\infty,\qquad v_Y=o(m_Y^2),      \tag{3.1}
\]

then `sum_p Z_p=infinity` almost surely conditional on `C_0`.
In particular, if `C_0` is a fixed finite-place prefix in the
certified class and `Z_p=h_p`, then

\[
 \boxed{\mu_\infty=0,\qquad A_{K_0}(B)=o(\log B).}  \tag{3.2}
\]

The same implication holds for `Z_p=b_p 1_(bar P_p=O)`, or for
return indicators at the depth-one primes of ETP3. No claim that
(3.1) holds for this curve is made.

**Proof.** For any fixed real `M>=0` and large `Y`,

\[
 \mathbb P(T_Y\le M\mid C_0)
                  \le \frac{v_Y}{(m_Y-M)^2}\longrightarrow0.
\]

This inequality follows directly by integrating the squared
deviation on the indicated event. Since `T_Y` is nondecreasing,
`{sum Z_p<=M}` is contained in every `{T_Y<=M}`. It has conditional
probability zero. Taking the union over integer `M` proves almost
sure divergence.

For `Z_p=h_p`, ETP5 says that the intersection of the all-places
soluble set with divergent total hazard is Haar-null. Every
all-places soluble vector in the chosen class lies in `C_0`.
Thus its intrinsic measure is zero, and EFA3 gives (3.2).
Pole weights are bounded above by `h_p`; at a depth-one branch
return the hazard is bounded below by `1/2`. The other two versions
therefore also force divergent total hazard. \(\square\)

The moments in (3.1) can be expanded as

\[
 m_Y=\sum_{p\le Y}\mathbb E(Z_p\mid C_0),
\]
\[
 v_Y=\sum_{p\le Y}\operatorname{Var}(Z_p\mid C_0)
      +2\sum_{p<q\le Y}\operatorname{Cov}(Z_p,Z_q\mid C_0).
                                                               \tag{3.3}
\]

ECM2 makes every fixed truncation exact under its input hypotheses.
The uniform estimates on these sums as `Y` grows remain additional
arithmetic obligations. For comparison, for unconditioned pole
weights `Z_p=b_p 1_(bar P_p=O)`, set the weights to zero outside
the active primes of ETP4. ECM1 gives exactly

\[
 m_Y=\sum_{p\le Y}\frac{b_p}{r_p},\qquad
 v_Y=\sum_{p\le Y}b_p^2\left(\frac1{r_p}-\frac1{r_p^2}\right)
       +2\sum_{p<q\le Y}\frac{b_pb_q(d_{p,q}-1)}{r_pr_q}.
                                                               \tag{3.4}
\]

Equation (3.4) is not the conditional formula for a class or prefix.

## 4. What a bounded second moment does and does not imply

A weaker estimate `E(T_Y^2|C_0)<=K m_Y^2` with `m_Y->infinity`
only gives divergence on a subset of conditional measure at least
`1/K`. Indeed Cauchy--Schwarz gives

\[
 \mathbb P(T_Y>M\mid C_0)
       \ge\frac{(m_Y-M)^2}{\mathbb E(T_Y^2\mid C_0)}.
\]

First let `Y` grow and then `M` grow. Positive-probability recurrence
alone does not rule out a positive-measure all-soluble set in its
complement.

There is a useful uniform alternative to (3.1). Suppose the same
finite constant `K` gives this bounded second-moment conclusion
after conditioning on every positive-measure intersection `C_0 cap D`,
where `D` is a clopen finite-quotient cell, and the corresponding
means all tend to infinity. Then divergence is almost sure in
`C_0`. To prove it, a positive-measure finite-sum subset of `C_0`
could be approximated arbitrarily well in measure by finite unions
of clopen cells. One such cell would have conditional finite-sum
probability greater than `1-1/K`, contradicting the bound just
proved in that cell. This uses elementary regularity of Haar measure
on the compact profinite group, not an independence or zero-one
assumption. The uniformity in `D` is essential.

## 5. The missing correlation estimate

Growth of `r_p` supplies an upper bound on a marginal branch
probability, `|B_p|/r_p`, but gives no control over `d_(p,q)` or
the densities `c_d` in (1.2). For example, in the abstract quotient
`Z->Z/(2n)Z`, an earlier condition that the integer is odd makes
the conditional probability of the identity residue zero, however
large `2n` is. Identical quotient maps also have common quotient
of that full order, giving much larger pair intersections than the
product of their marginals. These are finite-group obstructions to
an inference from order alone, not counterexamples asserted for the
elliptic curve.

For the depth-one ETP3 route, a positive-density set of primes with
`e_p=1` does not by itself supply (3.1). One would still need
conditional return means and their correlations. A fixed-prefix
lower bound of order `1/p` on a positive-density prime set, together
with a variance negligible compared with the squared accumulated
mean, would satisfy ECM3. Alternatively a lower bound of order
`1/p` for the returns conditioned on the *growing* soluble prefix
would feed ETP2 and EFA4 directly. These are different obligations.

Finite-level Kummer or division-field information can determine
common quotients and fixed congruence patterns. Positive density
of such primes does not by itself bound the conditioned common
quotients or the moment sum for an increasing prime range. The
present exact interface makes that missing estimate explicit; it
does not establish it or a vanishing all-primes measure.

The proofs here use finite abelian groups, character orthogonality,
elementary second moments and the previously proved EQA, FQS,
ETP and EFA inputs. No new finite checker is needed for these
written identities. No formal verification or external review
is claimed.
