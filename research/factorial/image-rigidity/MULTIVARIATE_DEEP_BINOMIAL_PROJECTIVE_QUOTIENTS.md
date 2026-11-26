# Depth-weighted projective restrictions from binomial Newton faces

Consider ordinary two-monomial consecutive factorial moments \(P_n,P_N\)
with \(N=n+1\) and \(n\ge2\). A common rational monic factor \(H\)
has even degree \(d>0\), because the even moment is strictly positive
on the real axis. Write its roots, counted with multiplicity, as
\(\alpha_1,\ldots,\alpha_d\). None is zero.

Binomial Newton faces give higher prime-power information about the
quotients \(\alpha_i/\alpha_j\). The information can be used
simultaneously at several root valuations and at primes that change
with the profile. The resulting estimate joins the Newton face depths
to the denominator of a scale-invariant factor polynomial. Near a
positive Gaussian contact there is an additional contribution from
the root quotients that tend to roots of unity.

## Exact quotient denominator

For a prime \(p\), extend \(v_p\) to an algebraic closure of
\(\mathbb Q_p\), with \(v_p(p)=1\). Let \(s_i=v_p(\alpha_i)\) and
put
\[
  A_p=\sum_{i<j}|s_i-s_j|,\qquad D_H=\prod_p p^{A_p}.
\]
The number \(D_H\) is an integer: it is the least denominator of
\[
  Q_H(Y)=\prod_{i\ne j}(Y-\alpha_i/\alpha_j)\in\mathbb Q[Y].
\]
Indeed, for a rational monic polynomial, the minimum valuation of
its coefficients is the sum of the negative root valuations. The
coefficient selecting all negative root valuations has a unique
lowest-valuation term, unless there are none, in which case the
leading coefficient suffices. Apply this to the paired valuations
\(s_i-s_j,s_j-s_i\). The integer \(D_H\) is invariant under a common
rescaling of the roots.

Let \(L_H\) be the denominator from the [projective factor sieve](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md),
defined using the roots \(\alpha_i^d/\prod_j\alpha_j\). Then
\[
 v_p(L_H)\le A_p\le \frac{2(d-1)}{d}v_p(L_H).                    \tag{1}
\]
In fact, with \(\bar s=\sum_i s_i/d\), the prior formula is
\(v_p(L_H)=\tfrac d2\sum_i|s_i-\bar s|\). The lower estimate follows
by expressing \(s_i-\bar s\) as the average of the \(s_i-s_j\)
and using the triangle inequality. For the upper estimate use
\(|s_i-s_j|\le |s_i-\bar s|+|s_j-\bar s|\) for each unordered
pair. In degree two the two denominators agree.

Fix a positive integer \(\ell\), and define
\[
 C_{ij}=(\alpha_i/\alpha_j)^\ell+(\alpha_j/\alpha_i)^\ell-2,
 \qquad V_{\ell,H}(T)=\prod_{i<j}(T-C_{ij}).
\]
This monic polynomial is rational and scale-invariant by symmetry.
Its least denominator is exactly \(D_H^\ell\). When \(s_i\ne s_j\),
one term of \(C_{ij}\) has strictly smaller valuation than the
other terms, so \(v_p(C_{ij})=-\ell|s_i-s_j|\). When \(s_i=s_j\),
the valuation is nonnegative. The same minimum-coefficient rule
proves the denominator identity.

Let \(t_\ell\) count the unordered pairs with \(C_{ij}=0\), or,
equivalently, \((\alpha_i/\alpha_j)^\ell=1\). Form the nonzero integer
\[
 I_{\ell,H}=D_H^\ell [T^{t_\ell}]V_{\ell,H}(T).
\]
The coefficient is the product of the nonzero \(C_{ij}\), up to
a sign, with the empty product equal to one. The denominator
cancellation gives the exact local formula
\[
   v_p(I_{\ell,H})=
        \sum_{\substack{i<j\\s_i=s_j\\ C_{ij}\ne0}}v_p(C_{ij}).    \tag{2}
\]
In particular the negative valuations from distinct slopes have
already been accounted for before testing the equal-slope pairs.

## Face depth and the all-degree inequality

At a prime \(p\) and a shared root valuation \(s\), write either
parent as \(F(z)=\sum_k f_k z^k\). Suppose that the minimum of
\(v_p(f_k)+ks\) is attained at exactly two indices \(a<b\), and
that \(b-a\) divides \(\ell\). The face depth is
\[
 \min_{k\notin\{a,b\}}\bigl(v_p(f_k)+ks-v_p(f_a)-as\bigr)>0.
\]
All ordinary parent coefficients are nonzero, and the parent
degrees here are at least two, so this minimum is finite. Define
\(\delta_p(s)\) as the largest depth supplied by either parent
under these conditions, and set it to zero if neither supplies
such a face.

For a root \(\alpha\) of valuation \(s\), divide \(F(\alpha)=0\)
by \(f_b\alpha^a\). The remaining terms give
\[
      v_p\!\left(\alpha^{b-a}/(-f_a/f_b)-1\right)\ge\delta_p(s)
\]
when using the parent of maximal depth. For two roots of valuation
\(s\), divide the two quantities \(\alpha^{b-a}/(-f_a/f_b)\).
Both are units, so the quotient differs from one to at least
the same depth. Raising to \(\ell/(b-a)\) preserves that depth
by the integral binomial expansion. Therefore
\[
 v_p(C_{ij})\ge 2\delta_p(s)
   \qquad(s_i=s_j=s,\ C_{ij}\ne0),                               \tag{3}
\]
because \(C_{ij}=(\rho^\ell-1)^2/\rho^\ell\) for
\(\rho=\alpha_i/\alpha_j\). The argument requires no bound
on \(p\), and no separation or residue multiplicity assumption.

Write \(m_{p,s}=\#\{i:s_i=s\}\) and let \(t_{p,s}\) count the
\(\ell\)-torsion pairs among those roots. Every torsion pair has
equal valuations at each prime, so \(\sum_s t_{p,s}=t_\ell\).
By (2) and (3) one has
\[
 v_p(I_{\ell,H})\ge
 U_{p,\ell}:=
 \left\lceil 2\sum_s
       \left(\binom{m_{p,s}}2-t_{p,s}\right)\delta_p(s)
 \right\rceil.
\]
Consequently every degree-\(d\) common factor of every ordinary
profile satisfies
\[
   \sum_p U_{p,\ell}\log p
     \le \ell\log D_H+\binom d2(2n\ell+2)\log 2 .                \tag{4}
\]
To see the absolute Archimedean bound, divide the lower moment by
its constant coefficient and use a positive argument scale so
its constant and leading coefficients are both one. Log
convexity of the gamma factors bounds each resulting coefficient
by the corresponding binomial coefficient. The monic root bound
and its reciprocal give normalized moduli in
\([2^{-n},2^n]\). Thus each \(\alpha_i/\alpha_j\) has modulus in
\([2^{-2n},2^{2n}]\), and \(|C_{ij}|\le 2^{2n\ell+2}\).
Taking logarithms of \(|I_{\ell,H}|\) proves (4). It applies at
every order and every positive profile independently of any
absolute angular cutoff.

The estimate is a finite test on the parent arrays even if \(H\)
is unknown. For each possible \(d,t_\ell\) choose shared Newton
root multisets \((m_{p,s})_s\) with total \(d\) and counts bounded
by the smaller parent widths. Counts at a slope with reduced
denominator \(b_s\) must be divisible by \(b_s\), since each
nonzero face of a rational factor has lattice endpoints. This
uses the [Newton polygon factor theorem](https://arxiv.org/pdf/0807.2620)
in its ordinary order-one case. For each choice, distribute
\(t_\ell\) pairs among the groups with
\(0\le t_{p,s}\le \binom{m_{p,s}}2\).
At a fixed prime take the maximum of
\[
 \ell\sum_{s<s'}m_{p,s}m_{p,s'}(s'-s)-U_{p,\ell}
\]
over these local choices. The sum of these maxima weighted by
\(\log p\), plus the error term in (4), must be nonnegative.
This keeps the valuation dispersion and the attained depth
linked to the same local choice. Primes outside the endpoint
and coefficient support contribute zero. The logs can be
avoided by comparison of products of positive rationals
after clearing the finitely many slope denominators.

## Improvement near a finite positive diffuse contact

Fix the retained core and a positive strength interval as in
the [physical diffuse comparison](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
Use the finite projective sieve's matched reference factors
\(H_*\) and its root error \(O(\eta^\gamma)\) in the common
rational normalization. The general multiplicity-safe exponent
is \(\gamma=1/[N(nN)^4]\), with the sharper local exponent
\(1/m_*\) when the contacting roots are simple. Root quotients
and \(C_{ij}\) are unchanged by this normalization.

For a particular reference factor, let \(t_{\ell,*}\) be the
number of its unordered \(\ell\)-torsion pairs. For sufficiently
small \(\eta\), each actual \(\ell\)-torsion pair belongs to a
matched reference torsion pair, by positive separation from
the finitely many non-torsion reference values. Since
\(C_{ij}=(\rho^\ell-1)^2/\rho^\ell\), every remaining pair
that tends to \(\ell\)-torsion contributes \(O(\eta^{2\gamma})\)
to the nonzero product. Hence
\[
  \sum_p U_{p,\ell}\log p+
           2\gamma(t_{\ell,*}-t_\ell)\log(1/\eta)
     \le \ell\log D_H+O(1),                                    \tag{5}
\]
where the constant depends on the core, the reference data,
the paired orders and \(\ell\). The case \(t_\ell=t_{\ell,*}\)
includes exact attainment of every limiting torsion relation.

In particular, at a branch with no limiting \(\ell\)-torsion,
even a single prime with face depth tending to infinity
forces growth of \(D_H\) if its eligible valuation group has
at least two common roots. A varying set of such primes
contributes their full depth-weighted product. The bound
applies to several valuation groups at once. Combining (1)
and (5) gives a necessary growth bound on the projective
denominator previously controlled by the gcd of adjacent
factorial coefficient ratios.

These bounds still allow profiles with growing denominator
dispersion and insufficient binomial face depth. Excluding
that regime requires an additional argument.

## Exact arithmetic fixture

Take an odd prime \(p\) and positive integers \(A,K\), and set
\(u=p^A/(1+p^A)\), \(b=p^K/(1+p^K)\). Use
\[
 H(z)=(z^2+bz+1)(z^2+ubz+u^2),\qquad P_4=H,\quad P_5=H(z+1).
\]
This is a fixture for the local identity, with positive parent
coefficients and no claim that the parents are factorial moments.
The four root valuations at \(p\) are \(0,0,A,A\). The coefficient
ordinates of \(H\) are \(2A,A+K,0,K,0\), so both width-two
faces have depth \(\min(K,2A)\). For \(\ell=2\), each within-group
root quotient has
\(v_p(\rho^2-1)=K\): the two unscaled unit roots have sum
\(-b\) and unit difference. There are two such pairs, and
all four cross-group pairs have different valuations.
Thus \(v_p(D_H)=4A\) and (2) gives \(v_p(I_{2,H})=4K\).
When \(K\le2A\), the bound by the two face depths is attained.
