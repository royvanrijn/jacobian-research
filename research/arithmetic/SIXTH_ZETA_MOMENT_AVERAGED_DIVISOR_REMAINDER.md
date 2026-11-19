# The sixth zeta moment and an averaged divisor remainder

**Mathematical status.** The sixth-moment asymptotic and the signed
averaged estimate (SM.15) below remain **open**. What is established
here is a conditional implication: that one particular averaged
estimate for the full ternary-divisor remainder implies the expected
sixth-moment leading term. The implication is an application of
Nathan Ng's inspected approximate-functional-equation and explicit
main-term calculations. It is not a new solution of the moment problem
or a novelty claim for its reduction to additive divisor sums.

The new target is weaker than a pointwise additive-divisor conjecture:
it permits cancellation over shifts and dyadic scales, uses only the
specific moment kernels, and can be stated using the unshifted
coefficient \(d_3\) alone. The graph comparison estimates in the current
arithmetic theory do not by themselves establish this target.

## 1. The moment and the classical reduction

Write
\[
 d_3(n)=\sum_{n_1n_2n_3=n}1,\qquad
 a_3=\prod_p(1-p^{-1})^4(1+4p^{-1}+p^{-2}).
\]
The desired leading asymptotic is
\[
 \int_0^T\left|\zeta\!\left(\tfrac12+it\right)\right|^6dt
       \sim\frac{42a_3}{9!}T(\log T)^9.                  \tag{SM.1}
\]
For a fixed \(\phi\in C_c^\infty((0,\infty))\), set
\[
 \omega_T(t)=\phi(t/T),\qquad
 I(T;\phi)=\int_0^\infty\omega_T(t)
                 |\zeta(\tfrac12+it)|^6dt.              \tag{SM.2}
\]
Every derivative of this weight has scale \(T\); thus the parameter
\(T_0\) in Ng's notation is comparable to \(T\). No shrinking time
transition is required for the leading asymptotic.

The precise primary source is Nathan Ng,
[*The sixth moment of the Riemann zeta function and ternary additive
divisor sums*, arXiv:1610.04977v3, 6 July 2021](https://arxiv.org/html/1610.04977v3),
published in *Discrete Analysis* 2021:6. Its Conjecture 2 assumes a
complete smooth additive-divisor main term with error
\(O(P^C X^{\vartheta+\varepsilon})\),
\(\tfrac12\le\vartheta<\tfrac23\), uniformly in the stated
shift and derivative ranges. Its Theorem 1.1 and Corollary 1.2 then
give a sixth-moment asymptotic. Here that pointwise hypothesis is
replaced by a single signed aggregate of its remainders. The necessary
explicit main is retained in full.

## 2. The exact moment kernel

Fix \(G(s)=e^{s^2}\). Put
\[
 g_0(s,t)=
 \left\{
 \frac{\Gamma((\tfrac12+s+it)/2)
       \Gamma((\tfrac12+s-it)/2)}
      {\Gamma((\tfrac12+it)/2)\Gamma((\tfrac12-it)/2)}
 \right\}^{3},\qquad
 V_t(z)=\frac1{2\pi i}\int_{(1)}\frac{G(s)}s
                                      g_0(s,t)z^{-s}\,ds.
                                                               \tag{SM.3}
\]
Ng's Proposition 2.1, equations (2.3)--(2.6), applies without an
additive-divisor hypothesis. At zero complex shifts its two halves
are identical: the factor \(X_{\mathcal I,\mathcal J;t}\) is one.
Consequently, for every fixed \(A>0\),
\[
 I(T;\phi)=
 2\sum_{m,n\ge1}\frac{d_3(m)d_3(n)}{\sqrt{mn}}
 \int_0^\infty\omega_T(t)(m/n)^{-it}V_t(\pi^3mn)\,dt
                                      +O_{A,\phi}(T^{-A}).
                                                               \tag{SM.4}
\]
The exact diagonal is
\[
 D(T;\phi)=2\sum_{n\ge1}\frac{d_3(n)^2}{n}
                   \int_0^\infty\omega_T(t)V_t(\pi^3n^2)\,dt.
                                                               \tag{SM.5}
\]
Both the factor two and the smooth \(t\)-dependent weight matter.
The off-diagonal moment is not an unweighted positive divisor count.

Choose a smooth \(W_0\), supported in \([1,2]\), with
\[
 \sum_{M\in\mathcal D}W_0(x/M)=1\quad(x\ge1),\qquad
 \mathcal D=\{2^{k/2}:k\in\mathbb Z,\ k\ge-1\},
 \qquad W(u)=u^{-1/2}W_0(u).
\]
For a box \((M,N)\), define the globally smooth kernel on the positive
quadrant by
\[
 F_{M,N}(x,y)=W(x/M)W(y/N)\frac1T
       \int_0^\infty\omega_T(t)(x/y)^{-it}
                                             V_t(\pi^3xy)\,dt.
                                                               \tag{SM.6}
\]
It has compact box support. On the additive line \(x-y=r\),
the phase \((x/y)^{-it}\) equals \((1+r/y)^{-it}\).
Thus (SM.6) has exactly the restriction used by Ng's
\(f_{r;M,N}\) in (5.18). His off-line extension can be used when
applying the derivative bounds of Lemma 5.3. The global extension
(SM.6) avoids an undefined off-line phase at large negative \(r\).
All functionals used below depend only on the restriction to the
additive line, so this choice changes none of them.

For a compactly supported smooth \(f\) on the positive quadrant,
write
\[
 A_r[f]=\sum_{m-n=r}d_3(m)d_3(n)f(m,n),\qquad r\ne0.
                                                               \tag{SM.7}
\]
The factor \(u^{-1/2}\) in \(W\) is essential: the contribution
of a box is \(T/\sqrt{MN}\) times its additive-divisor sum.

## 3. Dyadic ranges and shift tails

Fix a sufficiently small \(\eta>0\), independently of \(T\).
An explicit enlarged set of retained boxes is
\[
 \mathcal B_T=\left\{(M,N)\in\mathcal D^2:
     \tfrac14\le M/N\le4,\quad
     \min(M,N)\ge T^{1-2\eta},\quad MN\le T^{3+2\eta}\right\}.
                                                               \tag{SM.8}
\]
The constants and enlarged exponents include the nonnegligible boxes
of Ng's Proposition 5.2 and equations (5.10)--(5.19). There are
\(O_\eta(\log T)\) comparable pairs. From that proposition, or directly
from the decay of (SM.3) and integration by parts in \(t\),
\[
 I(T;\phi)=D(T;\phi)+
 2\sum_{(M,N)\in\mathcal B_T}\frac{T}{\sqrt{MN}}
                         \sum_{r\ne0}A_r[F_{M,N}]
                                  +O_{A,\phi,\eta}(T^{-A}).
                                                               \tag{SM.9}
\]
The inner sum can include all nonzero integers. Both its discrete
and continuous main-term versions have only finitely many nonzero
terms because the box restricts \(|r|\ll M+N\).

The negligible tails have three distinct causes. If \(mn\) exceeds
\(T^{3+2\eta}\), the Mellin integral in (SM.3) can be moved to a
large positive line. If the endpoints are not comparable, their
time phase is nonstationary. If their scale \(X\asymp M\asymp N\)
is below \(T^{1-2\eta}\), every nonzero integer difference has
\(T|\log(m/n)|\gg T^{2\eta}\). Repeated integration by parts
therefore gives arbitrary power decay. Divisor bounds and the
available number of boxes consume only polynomial losses, which
are paid by choosing the number of integrations last.

Within a retained box, the same argument gives negligible total
contribution from
\[
                 |r|>C\frac XT\,T^\eta                 \tag{SM.10}
\]
for a sufficiently large fixed \(C\), after as many integrations
as the desired error requires. The bound also applies to the
explicit main in Section 4. Its Ramanujan and residue coefficients
cost at most fixed powers of \(|r|+1\) and \(\log X\); the
arbitrarily strong smooth decay pays those losses as well.

Thus the shift length is \(H\asymp X/T\), with a harmless
\(T^\eta\) enlargement. It becomes \(X^{1/3}\) at the top scale
\(X\asymp T^{3/2}\). At lower scales this equality must not be
used. Boxes near or below \(X=T\) have few or no significant
nonzero shifts.

## 4. The full unshifted additive-divisor main

The following definitions specialize Ng's Definition 4 to \(d_3\).
They concern the arithmetic main term, not a graph ideal kernel.
For \(z\) near one, set
\[
 \begin{split}
 g_3(z,n)&=\prod_{p^\alpha\parallel n}
   \frac{\sum_{j\ge0}d_3(p^{j+\alpha})p^{-jz}}
        {\sum_{j\ge0}d_3(p^j)p^{-jz}},\\
 G_3(z,n)&=\sum_{d\mid n}\frac{\mu(d)d^z}{\varphi(d)}
                  \sum_{e\mid d}\frac{\mu(e)}{e^z}
                                         g_3(z,ne/d),\\
 S_r(z,w)&=\sum_{q\ge1}\frac{c_q(r)G_3(z,q)G_3(w,q)}{q^{z+w}}.
 \end{split}                                                   \tag{SM.11}
\]
Here \(c_q(r)\) is the Ramanujan sum. Explicitly, with
\(u=p^{-z}\),
\[
 \begin{split}
 g_3(z,p^\alpha)&=
       \frac{(\alpha+1)(\alpha+2)}2
       -\alpha(\alpha+2)u+\frac{\alpha(\alpha+1)}2u^2,\\
 G_3(z,p^\alpha)&=
       \frac{p\,g_3(z,p^\alpha)-p^z g_3(z,p^{\alpha-1})}{p-1}
                                  \quad(\alpha\ge1).
 \end{split}                                                   \tag{SM.12}
\]
For example \(G_3(1,p)=2-p^{-1}\).

These formulas also justify the analytic domain of (SM.11).
For \(p\nmid r\), its Euler factor is
\(1-G_3(z,p)G_3(w,p)p^{-z-w}\), which differs from one by
\(O(p^{-2+4\delta})\) on
\(|z-1|,|w-1|\le\delta\), with any sufficiently small fixed
\(\delta>0\). At \(p\mid r\), the Ramanujan coefficients vanish
above exponent \(v_p(r)+1\), leaving a finite local expression.
Thus \(S_r\) is holomorphic near \((1,1)\), with normal
convergence of its defining series for each fixed \(r\).
Equivalently, (SM.12) gives
\(G_3(z,q)\ll d(q)^C q^{\delta}\) on this neighborhood.
The elementary bound \(|c_q(r)|\le(q,r)\le|r|\) then
majorizes the series by
\(\sum_q(q,r)d(q)^{2C}q^{-2+4\delta}\), which converges
when \(\delta\) is sufficiently small. Fixed-order
derivatives only insert logarithmic factors and retain this
convergence. This also gives a polynomial-in-\(|r|\) bound
adequate for the arbitrarily strong smooth tail estimate.
Its growth when \(r\) varies is retained in the tail estimates;
normal convergence for fixed \(r\) is not a claim of uniform
boundedness for all \(r\).

Define the main functional by the double residue
\[
 \mathcal M_r[f]=
 \underset{z=1}{\operatorname{Res}}
 \underset{w=1}{\operatorname{Res}}
 \left\{\zeta(z)^3\zeta(w)^3 S_r(z,w)
  \int_{\max(0,r)}^\infty
       f(u,u-r)u^{z-1}(u-r)^{w-1}\,du\right\}.            \tag{SM.13}
\]
This is an explicit functional, defined independently of any
asymptotic for \(A_r[f]\). It retains a polynomial of degree at
most four in \(\log u\) and \(\log(u-r)\), with coefficients
involving the derivatives of \(S_r\) at \((1,1)\). Only two
derivatives in each variable can contribute to the residue.

To check the exact source match, replace \(\zeta(z)^3\) by
\(\prod_{i=1}^3\zeta(z+a_i)\), and similarly in \(w\), using
Ng's \(G_{\mathcal I},G_{\mathcal J}\). The simple poles at
\(z=1-a_i\) and \(w=1-b_j\) yield the nine terms in his
Conjecture 2: the remaining zeta products, the denominator
\(q^{2-a_i-b_j}\), and the integral factors
\(u^{-a_i}(u-r)^{-b_j}\) agree exactly. Letting all shifts tend
to zero gives (SM.13). Hence no singular individual nine-term
expression has to be evaluated at coincident shifts.

The leading term alone would be
\(S_r(1,1)\int f(u,u-r)\log^2u\log^2(u-r)\,du/4\).
It is insufficient for the moment reduction. The other residue
coefficients depend on the growing integer \(r\) and must be
kept until the signed aggregate has been evaluated.

## 5. The single signed aggregate and conditional sufficiency

Put
\[
 \mathcal E_r[f]=A_r[f]-\mathcal M_r[f],\qquad
 \mathscr R(T;\phi)=
 2\sum_{(M,N)\in\mathcal B_T}\frac{T}{\sqrt{MN}}
                     \sum_{r\ne0}\mathcal E_r[F_{M,N}].  \tag{SM.14}
\]
There are no absolute values inside either sum. In particular,
\(\mathscr R\) is not a sum of positive errors. For real
\(\phi\), exchanging the endpoints pairs conjugate terms, so
the aggregate is real. The new estimate to be proved is
\[
       \boxed{\mathscr R(T;\phi)=o_\phi(T(\log T)^9)}
       \quad\text{for every fixed }
                       \phi\in C_c^\infty((0,\infty)).  \tag{SM.15}
\]

**Theorem SixthMomentAveragedCriterion.** Assume (SM.15), with the
explicit kernel and main functional (SM.3), (SM.6), and (SM.13).
Then
\[
 I(T;\phi)=\frac{42a_3}{9!}T(\log T)^9
                         \int_0^\infty\phi(v)\,dv
                             +o_\phi(T(\log T)^9),      \tag{SM.16}
\]
and the unsmoothed assertion (SM.1) follows.

**Proof.** Substitute (SM.14) into the unconditional identity
(SM.9). This gives the exact arithmetic remainder and an explicit
diagonal-plus-main expression. The evaluation of that explicit
expression is the main-calculation portion of Ng's Sections 4,
5.3--5.5 and 6. Its application here is as follows.

The only use of the additive-divisor hypothesis in the original
replacement step, Section 5.2, is to replace each actual
\(A_r[f]\) by its specified main and then bound the difference.
Here the difference is kept exactly as (SM.14). All subsequent
operations on the main use explicit smooth integrals, Ramanujan
sums, Euler products and contour residues. The extension of the
shift range and the removal of dyadic restrictions use the
independent smooth tails from Section 3.

For distinct tiny complex shifts, the explicit-main calculation
and the diagonal residues give the twenty-swap expression
\[
 \int_0^\infty\omega_T(t)
 \sum_{\substack{\mathcal S\subset\mathcal I,
                  \mathcal U\subset\mathcal J\\
                  |\mathcal S|=|\mathcal U|}}
     \mathcal Z_{\mathcal I_{\mathcal S},\mathcal J_{\mathcal U}}(0)
       (t/(2\pi))^{-\sum\mathcal S-\sum\mathcal U}\,dt.
                                                               \tag{SM.17}
\]
Here \(\mathcal Z_{\mathcal A,\mathcal B}(s)\) is initially
the series \(\sum_n\sigma_{\mathcal A}(n)\sigma_{\mathcal B}(n)n^{-1-s}\)
in its convergence half-plane; its value at zero in (SM.17) is
defined by meromorphic continuation. The paired swapped sets are
\(\mathcal I_{\mathcal S}=(\mathcal I\setminus\mathcal S)\cup(-\mathcal U)\)
and
\(\mathcal J_{\mathcal U}=(\mathcal J\setminus\mathcal U)\cup(-\mathcal S)\),
as in Ng's Definition 3. The notation for the first set therefore
also depends on the paired selection \(\mathcal U\).
The intermediate unwanted residues cancel between the two halves
of the functional equation by his Proposition 6.6. This
cancellation is an identity among the explicit main terms.

The contour and smooth-tail errors in this explicit calculation
are smaller than \(T(\log T)^9\) for fixed \(\phi\) and small
fixed \(\eta\). The explicit-main and diagonal expressions are
holomorphic in the six shifts; the apparent singularities in
the individual terms of (SM.17) cancel. Cauchy's formula on
the small separated contours used in Ng's proof coalesces them
to zero, without invoking any estimate for actual additive
divisor sums. The leading coefficient of (SM.17) at zero is
\(42a_3/9!\), as checked explicitly in Section 6. Lower residue
degrees contribute \(O_\phi(T(\log T)^8)\).

It follows that the explicit diagonal-plus-main expression is
the right side of (SM.16), up to \(o(T(\log T)^9)\).
Adding the exact remainder (SM.14) and using (SM.15) proves
(SM.16). This imports the inspected main-calculation portion
of Ng's proof; it does not invoke his whole conditional
Theorem 1.1 under an unstated weaker arithmetic hypothesis.

For unsmoothing, fix nonnegative smooth upper and lower
approximations to any compact interval in \((0,\infty)\),
first let \(T\) grow and then let their transition widths
shrink. A fixed smooth majorant of \([1,2]\) supplies
\[
 \int_U^{2U}|\zeta(\tfrac12+it)|^6dt
                              \ll U(\log U)^9
\]
for every sufficiently large \(U\). Summing over dyadic
intervals below \(\varepsilon T\) bounds the low-time tail
by \(O(\varepsilon T(\log T)^9)+o(T(\log T)^9)\).
Let \(\varepsilon\) tend to zero after \(T\). This proves
(SM.1). No power-saving moment remainder is asserted. \(\square\)

## 6. Leading-coefficient diagnostics

There is a short exact calculation of the coefficient in (SM.17).
At the leading logarithmic degree,
\[
 \mathcal Z_{\mathcal A,\mathcal B}(0)
       \sim a_3\prod_{a\in\mathcal A,b\in\mathcal B}(a+b)^{-1}.
\]
Take the generic scaled shifts
\(\boldsymbol\alpha=(1,2,3)z\) and
\(\boldsymbol\beta=(4,5,6)z\). For subsets \(S,U\) of the
unscaled vectors of the same size, put
\[
 c_{S,U}=\sum S+\sum U,\qquad
 d_{S,U}=\prod_{a\in\boldsymbol\alpha_S,
                    b\in\boldsymbol\beta_U}(a+b),
\]
where the vectors in this product are unscaled and have the
selected entries replaced by the negatives of the selected
entries from the other vector. The twenty exact rational
terms obey
\[
 \sum_{|S|=|U|}\frac{(-c_{S,U})^k}{d_{S,U}}=0
                      \quad(0\le k<9),\qquad
 \sum_{|S|=|U|}\frac{(-c_{S,U})^9}{d_{S,U}}=42.           \tag{SM.18}
\]
Taylor expansion of \(e^{-zc_{S,U}}\) therefore gives
\(42/9!\) for the constant after the nine polar powers
cancel. This is an exact twenty-term arithmetic identity.
It verifies the leading swap coefficient, not (SM.15).

The diagonal gives a second useful check. Since
\[
 \sum_{n\ge1}\frac{d_3(n)^2}{n^{1+2s}}
          =\zeta(1+2s)^9\mathcal A(1+2s),\qquad
                         \mathcal A(1)=a_3,
\]
the diagonal contour at \(s=0\) gives
\[
 D(T;\phi)=
 \frac{a_3}{9!}\,2(3/2)^9 T(\log T)^9\int\phi
                          +O_\phi(T(\log T)^8)
                          +o(T(\log T)^9).
                                                               \tag{SM.19}
\]
The full diagonal coefficient, scaled by \(a_3/9!\), is
\(19683/256\), whereas the desired total coefficient is
\(42\). The combined explicit off-diagonal main consequently
has coefficient
\[
                42-\frac{19683}{256}
                         =-\frac{8931}{256}.             \tag{SM.20}
\]
Its negativity is consistent with the oscillatory kernel and
rules out a positive-count surrogate that discards that phase.

## 7. Error estimates sufficient for the aggregate

At a dyadic scale \(X\), the weight is \(T/X\) and there
are about \(X/T\) significant shifts. A pointwise full-main
error \(O(X^\vartheta)\) would therefore contribute
\(O(X^\vartheta)\), before derivative and logarithmic losses.
At \(X=T^{3/2}\), this is \(T^{3\vartheta/2}\), explaining
the sufficient threshold \(\vartheta<2/3\).

An estimate \(O_A(X/(\log X)^A)\), even for every fixed
\(A\), instead gives a top-scale absolute error of order
\(T^{3/2}/(\log T)^A\). It does not establish (SM.15).
One cannot choose \(A=A(T)\) from an arbitrary-fixed-\(A\)
statement. A signed average preserving (SM.6) may require
less than such individual error bounds; that is the precise
analytic opportunity retained by (SM.15).

## 8. Scope of the deduction

The source statements and proof interfaces listed in Section 1
were inspected in Ng's exact v3 edition. The main functional
is a coalesced residue formula, not an assumed asymptotic.
The conditional implication is a classical application with
the arithmetic replacement error isolated in one signed
aggregate. Its proof does not establish that aggregate.

The exact coefficient identities and normalization diagnostics
are finite arithmetic checks. They supply no numerical evidence
for an infinite error estimate. Independent mathematical review,
formal verification, an effective onset, and a proof of the
sixth-moment conjecture are not claimed.

## 9. Graph precision and parameter obligations

At the top scale, set \(H=X/T=X^{1/3}\). Before the prefactor
\(1/H\), the target raw aggregate is \(o(X(\log X)^9)\).
If individual remainders have size \(X\), square-root cancellation
over \(H\) shifts gives only \(X\sqrt H=X^{7/6}\), hence
the moment scale \(X^{5/6}=T^{5/4}\). Thus ordinary square-root
shift cancellation alone does not reach the target.
This interpretation agrees with the explicit discussion in
[de Faveri and Pandey, Section 1.1, equation (3)](https://arxiv.org/html/2610.01035v1).
Their moment upper bound is not a uniform error theorem for the
individual dyadic blocks used here.

### 9.1. A local Fourier formulation and its limitation

For finite sequences \(a,b\), put
\(A(\theta)=\sum a_n e(n\theta)\) and
\(B(\theta)=\sum b_n e(n\theta)\), where
\(e(u)=e^{2\pi iu}\). If \(K\) is Schwartz and
\(\widehat K(v)=\int K(u)e(-uv)\,du\), Fourier inversion
and Poisson summation give the exact identity
\[
 \frac1H\sum_hK(h/H)\sum_n a_n\overline{b_{n+h}}
 =\int_0^1 A(\theta)\overline{B(\theta)}
       \sum_{j\in\mathbb Z}\widehat K(H(j-\theta))\,d\theta.
                                                               \tag{SM.21}
\]
For \(K(v)=\int\phi(u)e^{iuv}\,du\), with \(\phi\)
supported in a fixed positive interval, its Fourier multiplier
is supported in a band with \(|\theta|\asymp1/H\) on the
centered torus. This identifies a possible local-energy route.

Formula (SM.21) is a translation-invariant prototype. The actual
kernel (SM.6) depends on \(n\) through both the logarithmic
phase and the gamma cutoff. Freezing \(n\) over a factor-two
interval is not an established replacement. A proof must retain
that dependence, for example by a controlled Mellin separation.
Furthermore a per-box energy bound of order
\((X/H)(\log X)^8=T(\log X)^8\) could accumulate to
\(T(\log T)^9\) over \(O(\log T)\) boxes. The required
estimate must include an additional saving or control the signed
multiscale sum directly, with its cross terms.

### 9.2. Exact invariance does not supply moment precision

**Proposition InvariantMomentPrecisionObstruction.** Let
\(X=T^{3/2}\), \(L=\log X\), \(y=\exp(L^{2/5})\),
and let \(\mathcal P\) be any subset of the primes at most
\(y\). Write \(n_{\mathcal P}\) for the largest divisor of
\(n\) supported on \(\mathcal P\), and
\(n_{\mathcal P'}=n/n_{\mathcal P}\). Fix nonzero
\(W\ge0\) in \(C_c^\infty((1,2))\). With
\(t_0=3T/2\), \(\delta=e^{-\sqrt L}\), define
\[
 a_X(n)=\delta(n_{\mathcal P'})^{it_0},\qquad
 A_X(t)=\sum_n a_X(n)W(n/X)n^{-1/2-it}.
\]
Then \(a_X(pn)=a_X(n)\) for every \(p\in\mathcal P\),
and \(|a_X(n)|=\delta\). Every fixed-logarithmic absolute
mean or correlation bound holds, even with tests bounded by any
fixed power of \(L\). Nevertheless, uniformly in \(\mathcal P\),
\[
 \int_T^{2T}|A_X(t)|^2dt
       \gg_W \delta^2X(\log y)^{-6}
       \gg T(\log T)^9.                                 \tag{SM.22}
\]

**Proof.** Put \(f(n)=(n_{\mathcal P})^{-it_0}\). For
\(\Re s>1\), its exact series is
\[
 \sum_n f(n)n^{-s}=\zeta(s)C_{\mathcal P}(s),\qquad
 C_{\mathcal P}(s)=\prod_{p\in\mathcal P}
                  \frac{1-p^{-s}}{1-p^{-s-it_0}}.
\]
The finite Euler product is holomorphic on \(\Re s>0\).
On \(\Re s\ge1-c/\log y\), it is bounded by
\((\log y)^{C}\), uniformly in its imaginary argument and
in \(t_0\). Indeed its absolute value is at most
\(\prod_{p\le y}(1+p^{-\Re s})/(1-p^{-\Re s})\), and
\(\sum_{p\le y}p^{-\Re s}\ll\log\log y\).
At one, Mertens' product gives
\[
 |C_{\mathcal P}(1)|\ge
 \prod_{p\le y}\frac{1-p^{-1}}{1+p^{-1}}
                                      \gg(\log y)^{-2}.
\]
For fixed smooth \(U\), shift its Mellin inversion line to
\(\Re s=1-c/\log y\). The only crossed pole is that of
\(\zeta\) at one. Classical strip growth of \(\zeta\) and
the rapid Mellin decay of \(U\) bound the new line, uniformly,
by \(X\exp(-cL^{3/5})L^C\). Consequently
\[
 \sum_n f(n)U(n/X)
 =X C_{\mathcal P}(1)\int U
                   +O_U(Xe^{-cL^{3/5}}L^C).
\]
Apply this with \(U(u)=W(u)u^{-1/2}\). At \(t=t_0\),
\(a_X(n)n^{-it_0}=\delta f(n)\), so
\(|A_X(t_0)|\gg_W\delta\sqrt X(\log y)^{-2}\).
Removing the common \(X\)-phase and differentiating gives
\[
 \left|\frac d{dt}\{X^{i(t-t_0)}A_X(t)\}\right|
 \le\delta\sum_nW(n/X)n^{-1/2}|\log(n/X)|
                                      \ll_W\delta\sqrt X.
\]
The resonance retains at least half its magnitude on an interval
of width \(\gg_W(\log y)^{-2}\). This proves the first
inequality in (SM.22). Its ratio to the desired moment scale is
\(\gg T^{1/2}e^{-2\sqrt{(3/2)\log T}}(\log T)^{-57/5}\),
which tends to infinity. Exact invariance and all fixed-log
smallness assertions follow directly from the definition. \(\square\)

This is a proved information obstruction for general coefficients,
not a counterexample for \(d_3\) or for (SM.1). It uses no ZF
or graph input. Even its nonnegative envelope is the constant
\(\delta\); \(a_X/\delta\) is completely multiplicative.
The example shows that a new proof must use arithmetic information
stronger than invariance, bounded moments and arbitrary fixed-log
smallness alone.

There is a related obstruction at the rare-bank discard step.
Let \(Q\) be the product of selected bank primes, and
\(\eta_X=\prod_{p\mid Q}(1-1/p)\gg L^{-C}\). The ordinary
fundamental lemma at level \(X^\epsilon\), since
\(\log y=o(\log X)\), gives
\[
 \sum_{(n,Q)=1}W(n/X)n^{-1/2}
 \sim\eta_X\sqrt X\int W(u)u^{-1/2}\,du.
\]
The bounded abstract sequence
\(c_n=\mathbf1_{(n,Q)=1}W(n/X)n^{it_0}\) has every
successfully banked correlation identically zero. Yet its
Dirichlet-polynomial energy over a fixed small interval about
\(t_0\) is \(\gg_W\eta_X^2X\), by the positive cosine
sum after removing the common \(X\)-phase. This exceeds
\(T(\log T)^9\). Thus a density-only estimate for discarded
bank profiles cannot replace an estimate in the actual moment
norm. This auxiliary example is neither an invariant nor an
arithmetic \(d_3\) counterexample.

### 9.3. Power-size prime banks do not repair the existing probes

The existing [imported graph source](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md)
uses \(\asymp\log\log X\) banks and fixed mark parameters.
If primes instead lie in \([X^a,X^b]\), with fixed
\(0<a<b\), their total reciprocal mass is
\(\log(b/a)+o(1)\). Requiring each bank to have a fixed
positive mass permits only finitely many such banks. This already
changes the source hypotheses.

The same fixed-probe construction also has an intrinsic lower
bound, independent of prime sizes, at zero frequency. The next
statement concerns that construction only.

**Proposition FixedProbePowerBarrier.** In the exact-replacement
limit of the source's symmetric-probe ideal construction, fix
\(M=J+\ell\) slots per group, \(1\le m\le J\) designated
probes, and \(s_1,s_2\ge1\) groups in its two blocks. At
zero additive frequency and Mellin phase the residual norm is at
least
\[
 r^{s_1+s_2}>0,\qquad
 r=\binom{M}{m}^{-1}.                                  \tag{SM.23}
\]
The lower bound persists, with a factor one half, for replacement
operators whose uniform tuple-\(L^2\) error tends to zero.

**Proof.** The exceptional-pair formula after source (4.26)
is the residual
\(R=Q_1+Q_2-Q_1Q_2\), where \(Q_i\) is the product of
\(I-E_j\) over its designated probes and \(E_j\) averages
coordinate \(j\). In each group's label probability space
choose a mean-zero norm-one \(v\), and on its \(M\)
independent slots take
\[
 F_g=\binom{M}{m}^{-1/2}\sum_{|S|=m}\prod_{j\in S}v(p_j).
\]
Orthogonality gives \(\|F_g\|_2=1\), and symmetrization
fixes this vector. After averaging the last \(\ell\)
coordinates its squared norm is
\(c=\binom{J}{m}/\binom{M}{m}\). The squared mass surviving
all its \(m\) probes is exactly \(r\), since only their
one designated subset remains. For the tensor product over all
groups, the residual quadratic form is
\[
 \langle F',RF'\rangle
 =r^{s_1}c^{s_2}+c^{s_1}r^{s_2}-r^{s_1+s_2}
                                      \ge r^{s_1+s_2}.
\]
The unaveraged input has norm one, so this is a lower bound for
the corresponding ideal pairing operator norm. Uniform operator
approximation preserves half the bound. \(\square\)

Consequently larger primes alone cannot turn this fixed finite
construction into an \(X^{-\delta}\) contraction. Growing
its marks is also constrained: if every mark is at least
\(X^a\), a physical \(M\)-list dividing an integer of
size \(O(X^C)\) requires \(aM\le C+o(1)\). Thus fixed
positive-power prime sizes permit only bounded mark length.
The source CRT and residue-replacement arguments impose further
budgets; this calculation does not exclude a redesigned kernel
that uses the restricted moment weights instead of a norm bound
over arbitrary endpoint vectors.

## 10. The remaining analytic statement and verification

The precise remaining theorem is (SM.15), with the full main
(SM.13), the exact kernel (SM.6), and cancellation preserved
over all additive shifts and retained scales. None of the
obstructions above proves this estimate false. They show why
the current graph estimates do not establish it by a change of
parameters or by an absolute error sum.

A potentially stronger approach must control all bank profiles
in the actual time-moment norm, or prove a suitable arithmetic
local Mellin-energy estimate for the centered \(d_3\) coefficients.
The real kernel, cross terms and the sum of dyadic scales must
remain in that statement. These are open obligations here.

The [exact checker](verify_sixth_zeta_moment_main.py) verifies
80 swap terms on four generic rays, cancellation of 36 polar
coefficients, and two polynomial generating-function identities.
It returns the scaled coefficients \(42\), \(19683/256\),
and \(-8931/256\). Run
`python3 research/arithmetic/verify_sixth_zeta_moment_main.py`.
Its output explicitly records that the sixth-moment remainder
is not proved. The information obstruction and the fixed-probe
calculation have written internal application review; no external
human review or formal verification is asserted.
