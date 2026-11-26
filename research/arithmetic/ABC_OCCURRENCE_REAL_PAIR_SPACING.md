# Real pair spacing at the balanced occurrence endpoint

Two small-radical integers can be counted inside a fixed relative
power strip without requiring the pair to be primitive. The count has
half the sum of the individual radical exponents. At a fixed height
gap, this makes it possible to extract a full-exponent occurrence
matching induced on short real bins without a prime test.

## A uniform relative strip bound for actual pairs

**RSI1.** Fix \(A_0,\eta,\tau,K_0>0\). For all sufficiently small
fixed \(L>0\), put \(n=X^L\). For \(a,b\in[0,A_0]\) and
\(K>0\), the number of ordered positive integer pairs \(U,V\le X\)
satisfying
\[
 \operatorname{rad}(U)\le n^a,\qquad
 \operatorname{rad}(V)\le n^b,\qquad
 \left|\log(U/V)-\log K\right|\le K_0X^{-\tau}
\]
is
\[
 \ll_{A_0,\eta,\tau,K_0,L,\varepsilon}
               n^{(a+b)/2+\eta}X^\varepsilon
 \qquad(\varepsilon>0).
 \tag{RSI1}
\]
The estimate is uniform in \(a,b,K\). More generally, given
a fixed \(\zeta>0\) it holds at every fixed \(L>0\) for which
\[
               \frac{\tau}{L}
                       \ge \frac{a+b}{2}+\max(a,b)+\zeta .
 \tag{RSI1a}
\]
The constant may also depend on \(\zeta\). Here \(\tau\) may
depend on \(L\), and the margin in (RSI1a) is imposed before
the height limit. The same count with \(V=1\) gives the coefficient \(a/2+\eta\) for a single
small-radical integer in the strip. A relative condition
\(|(U/V)/K-1|\le K_0X^{-\tau}\) gives the same bound,
after increasing \(K_0\) in the logarithmic condition.

Fix an integer \(k\) at this budget with
\(1/(k+1)\le\eta L/8\). Split each integer into the product \(F\)
of its exact prime powers with exponent exceeding \(k\) and
the remaining exact-degree root factors,
\[
 U=F_U\prod_{e=1}^k x_e^e,\qquad
 V=F_V\prod_{e=1}^k y_e^e .
\]
Every \(F\) has \(\operatorname{rad}(F)\le X^{1/(k+1)}\).
Enumerating the two radicals and then their realizations at
height at most \(X\) costs \(n^{\eta/4}X^{o_X(1)}\) by
uniform fixed-support counting. The root factors \(x_e,y_e\)
are squarefree. Within each complete value their product is
bounded by its radical cap.

For each of the \(2k\) roots round its \(\log_n\)-height upward
to a positive point of a fixed mesh of size at most
\(\min(\eta,\zeta)/(16(k+1))\), charging a mesh step also for
height zero. Denote the resulting positive caps by \(b_j\)
and the signed exponents of the roots by \(e_j'\). Then
\(\sum_j b_j\le a+b+\min(\eta,\zeta)/8\) and
\(\max_j b_j\le\max(a,b)+\min(\eta,\zeta)/16\).
The number of cap vectors is
bounded at fixed \(A_0,k,\eta\). After \(F_U,F_V\) and a cap
vector are fixed, enlarge the possible roots to the independent
integer ranges \(1\le z_j\le n^{b_j}\). The original pairs
inject into these tuples for fixed cofactors, even if \(U,V\)
have common primes or the quotient admits more than one pair.

Let \(R=\frac12\sum_j b_j\) and \(H=X^\tau/(2K_0)\),
and use the triangular kernel
\[
 {\cal K}_H(s)=(1-H|s|)_+=\int e^{2\pi its}\,d\mu_H(t),
 \qquad
 d\mu_H(t)=H^{-1}
          \left(\frac{\sin(\pi t/H)}{\pi t/H}\right)^2dt .
\]
The quotient in the density at zero is defined by continuity.
The Fourier identity follows by convolving the indicator
of \([-1/2,1/2]\) with itself and rescaling. The measure
\(\mu_H\) is a probability measure. The kernel is at least
\(1/2\) on the target strip. Apply it with \(s=\sum_j e_j'\log z_j-\log(KF_V/F_U)\).
After Fourier expansion and taking absolute values, the count
in the fixed box is at most a constant times
\[
  \int \prod_j |S_j(t)|\,d\mu_H(t),
  \qquad S_j(t)=\sum_{z\le n^{b_j}}
                   \exp(2\pi it e_j'\log z).
\]

Take \(h_j=\lceil R/b_j\rceil\). Then
\(\sum_j(2h_j)^{-1}\le1\) and \(h_j b_j\le R+b_j\).
In the \(2h_j\)-th moment of \(S_j\), collapse the two halves to
products \(P,Q\le n^{h_j b_j}\). For \(P\ne Q\),
\(|e_j'\log(P/Q)|\ge n^{-h_j b_j}\), since \(|e_j'|\ge1\)
and the logarithm difference is at least \(1/\max(P,Q)\).
The upper bounds on the rounded caps and (RSI1a) imply
\(\tau/L\ge R+b_j+7\zeta/8\). Thus \(n^{-h_j b_j}>1/H\)
for large \(X\), and the kernel vanishes on every such term.
For fixed \(\tau>0\), the condition (RSI1a) follows uniformly
in \(a,b\) by taking \(L\) sufficiently small.
If \(P=Q\), choose the first half at cost \(n^{h_j b_j}\);
the ordered factorizations in the second half number \(X^{o_X(1)}\)
by the fixed-order divisor bound. The orders are fixed before
the height limit. Thus
\[
       \int |S_j(t)|^{2h_j}d\mu_H(t)
                    \ll n^{h_j b_j}X^{o_X(1)} .
\]
Hölder's inequality, with the constant function one supplying
any unused reciprocal exponent, gives
\(n^{\sum_j b_j/2}X^{o_X(1)}\) for the fixed box.
Restoring the cofactors and summing the finitely many cap
vectors proves (RSI1).

The logarithmic averaging and Hölder step are the singleton
case of (RG1). Applying them to full root tuples makes the
bound available for ordered pairs with arbitrary common prime
support and for a prescribed real quotient target. The
single-letter consequence improves the coefficient in (SI1)
from \(1-\tau/2\) to \(1/2\) for \(0<\tau<1\), with an
arbitrarily small additive allowance in base \(n\).
[Logarithmic moments (RG1)](ABC_REAL_STRIP_GROUP_DICTIONARIES.md#logarithmic-moments)
· [Relative short-interval bound (SI1)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#an-elementary-bound-in-relative-intervals)
· [Classical mean values of Dirichlet polynomials](https://doi.org/10.1112/jlms/s2-8.1.73)

If every exact-valuation layer of both inputs has radical
at most \(n^\vartheta\), the same estimate follows from
\[
 \frac{\tau}{L}\ge
 \frac{a+b}{2}+\min\{\max(a,b),\vartheta\}+\zeta .
                                                               \tag{RSI1b}
\]
Indeed the rounding now gives
\(\max_j b_j\le\min\{\max(a,b),\vartheta\}+\zeta/16\);
the identical exact-diagonal argument applies. In a
tall actual box with a wholly nonatomic weak profile,
the individual normalized exact-layer caps tend to
zero. Otherwise a layer of mass bounded below has
bounded scaled degree by the numerical height and
gives a limiting atom on a subsequence. This applies
to boxes formed before the height limit, with fixed
positive tolerances and the usual tail truncation.

At the level of arbitrary pairs this scale is sharp
for diffuse dictionaries. A dictionary of
\(n^{d-o(1)}\) distinct values at height comparable
to \(X\), divided into \(O(n^c)\) real bins of width
\(Xn^{-c}\), has at least \(n^{2d-c-o(1)}\) distinct
ordered near pairs when \(c<d\). The bound follows
by Cauchy–Schwarz on bin occupancies and survives
deleting the diagonal. The free dictionaries in
(IM4) without positive escaped weight give actual
small-radical examples with \(d=1/3\). The near-pair
bound in (RSI1) counts pairs regardless of coprimality,
as does this scale comparison.
[Diffuse free dictionaries (IM4)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#a-primitive-two-cap-induced-matching-with-diffuse-degrees)

## Real-bin inducedness under a fixed numerical gap

**RSI2.** Fix \(0<h_0\le h<1\), \(A_0>0\), and \(\eta>0\).
At a sufficiently small fixed budget \(L>0\), suppose a
complete-label matching of \(M\) actual rows satisfies
\[
 A_i+B_i=C_i,\qquad X/2<C_i\le X,
 \qquad X^{h_0}\le A_i\le X^h .
\]
Require \(\operatorname{rad}(B_i)\le n^a\) with
\(0\le a\le A_0\) and \(\operatorname{rad}(A_i)\le n^{A_0}\).
Let \(D=\lfloor X^\theta\rfloor\), where \(\theta\in(0,h)\)
is sufficiently small in terms of \(h_0,h,A_0,\eta\).
There is a subset of indices of size
\[
 \gg_{h_0,h,A_0,\eta,L,\varepsilon}
       \frac{M^2}{X^\varepsilon
                     (n^{a+\eta}+M n^\eta)}
 \tag{RSI2}
\]
and a single carry \(\kappa\in\{0,1\}\) on that subset
such that, writing \(b_D(x)=\lfloor x/D\rfloor\),
\[
  b_D(A_i)+b_D(B_j)+\kappa=b_D(C_k)
        \quad\Longleftrightarrow\quad i=j=k .
\]
Any original radical cap or repair cutoff on the target
continues to hold on the retained rows.

The same assertion holds when \(h=h(L)\in[h_0,1)\) has
\[
                     \frac{1-h(L)}{L}\ge 2a+\zeta
 \tag{RSI2a}
\]
for a fixed \(\zeta>0\), with \(L\) sufficiently small
for (SI2). The threshold on \(L\) and the choice of
\(\theta<h_0\) can be made independently of \(h(L)\).
The height tends to infinity at each fixed \(L,h(L)\).

All \(B_i\) are at least \(X/3\) for large \(X\).
Apply (RSI1) with \(K=1\) and \(\tau=1-h\) to bound
the number of distinct-index pairs with
\(\lvert B_i-B_j\rvert\le4X^h\) by
\(O(n^{a+\eta}X^\varepsilon)\). Both pair caps are \(a\),
so (RSI2a) is exactly (RSI1a) for this application.
The corresponding
logarithm difference is at most \(12 X^{-(1-h)}\).
Choose \(\theta<h\) so that (SI2), used with the
positive size floor \(h_0\) and the given short-role
radical cap, bounds the number of distinct-index pairs
with \(\lvert A_i-A_j\rvert\le4D\) by
\(O(Mn^\eta X^\varepsilon)\). Form the graph in which
an edge comes from either type of pair.
The independent set bound \(M^2/(M+2|E|)\) follows
by retaining each vertex that precedes all of its
neighbors in a uniform random ordering and applying
Cauchy's inequality to the resulting expected size.
On such an independent set retain the more frequent
of the two diagonal carries. This proves the size
bound (RSI2), including the loss of at most a factor
two for the carry.

On those indices the bin identity implies
\(\lvert A_i+B_j-C_k\rvert<2D\) by writing each
integer as its bin times \(D\) plus a remainder.
Using the actual \(k\)-th row,
\[
 |B_j-B_k|\le |A_i-A_k|+2D\le X^h+2D<4X^h
\]
for large \(X\), so \(j=k\) by independence of the
chosen indices. The same identity then implies
\(\lvert A_i-A_k\rvert<2D\), giving \(i=k\).
The diagonal is retained by the carry selection.
[Intervals at a positive size floor (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)

For a balanced occurrence cohort with individual caps
\(n^{1/3+o(1)}\) and row exponent \(1/3-o(1)\),
the complete-label matching extraction in (IM1)
permits \(a=1/3+o(1)\) and \(M=n^{1/3-o(1)}\).
For a fixed \(h<1\) this gives an induced matching with
the same limiting exponent on the real bins alone.
For any fixed \(c>2/3\) it also applies to
\(A_i\le X^{1-cL}\) when the cap and row exponent
allowances are sufficiently small in terms of \(c-2/3\).
If the tall matching dictionary used in the spacing step
has a wholly nonatomic exact-valuation profile, (RSI1b)
instead gives the same conclusion for any fixed \(c>1/3\),
using a sufficiently small fixed bound on its individual
layer masses. The positive lower height floor can follow
from (SH1).
Choose the normalized allowances and the short-interval
exponent before the fixed budget and its height limit.
The bin labels can be large, and the original radical
caps concern the complete integers. A power saving from
this induced configuration and control of summands
above the shrinking-gap threshold remain open.
[Exact and mixed induced matchings (IM1–IM3)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md)
· [Strict short-height sector (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)
