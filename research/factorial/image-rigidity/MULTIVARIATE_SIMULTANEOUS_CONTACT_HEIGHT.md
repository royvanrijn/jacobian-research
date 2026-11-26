# Simultaneous height bounds at algebraic diffuse contacts

Fix an ordinary two-monomial core and adjacent orders \(n,N=n+1\),
with \(n\ge2\). Suppose a sequence of exact ordinary common factors
approaches a fixed positive diffuse reference contact. Use the
[common rational argument normalization](MULTIVARIATE_DIFFUSE_CLOCK_PRIME_HEIGHT.md#the-two-eliminants-in-rational-coordinates).
Write \(H\) for a monic rational common factor of degree \(d>0\)
in that variable. Its roots \(\alpha_i\) are nonzero and nonreal,
and \(c=\prod_i\alpha_i>0\). Define
\[
       B_H(Y)=\prod_i(Y-\alpha_i^d/c).
\]
The [projective factor comparison](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#consequences-at-finite-physical-diffuse-contacts)
gives a finite list of limiting algebraic factor selections \(H_*\).
After selecting one of them, coefficients of \(H,B_H\) and the norm
\(c\) approach the respective limiting coordinates with error
\(O(\eta^\gamma)\). Here \(\eta\) is the largest added relative gap
and \(\gamma>0\) is the reference comparison exponent. The argument
below admits different error rates for different coordinates.

For a positive integer \(m\) and a finite prime set \(S\), write
\(m_S=\prod_{p\in S}p^{v_p(m)}\) and \(m_{S^c}=m/m_S\).
Select rational coordinates \(y_1,\ldots,y_r\) from the factor data
with real algebraic limits \(\theta_1,\ldots,\theta_r\) such that
\(1,\theta_1,\ldots,\theta_r\) are linearly independent over
\(\mathbb Q\). Let \(D_i\) be the positive reduced denominator
of \(y_i\), put \(D=\operatorname{lcm}_{i}D_i\), and write
\(e_i=|y_i-\theta_i|\).

**Theorem.** For every fixed finite prime set \(S\) and every
\(\epsilon>0\), sufficiently close to the reference factor one has
\[
 D^\epsilon D_{S^c}^{\,r+1}
               \prod_{i=1}^r (D_i)_S e_i\ \ge C>0.                 \tag{1}
\]
This bound is independent of the tail's row counts, gaps and heights.
The exclusion neighborhood from this theorem is qualitative. The same
conclusion holds for arbitrary rational vectors approaching the
selected algebraic coordinates.

In particular, when \(e_i=O(\eta^{\gamma_i})\) with \(\gamma_i>0\),
every such collision requires
\[
 (r+1)\log D_{S^c}+\sum_i\log (D_i)_S+\epsilon\log D
          \ \ge\ \left(\sum_i\gamma_i\right)\log(1/\eta)-O(1).
                                                                    \tag{2}
\]
One can use any independent subset, and test all such subsets
separately. Taking \(S\) empty permits freely moving denominator
primes. Then the exponent in a bound for \(D\) is
\(\sum_i\gamma_i/(r+1+\epsilon)\).
For one irrational coordinate the resulting exponent is
\(\gamma_1/(2+\epsilon)\), independent of its algebraic degree.
When the selected denominators are supported on \(S\), (2)
instead uses the sum of their logarithms and the term
\(\epsilon\log D\).

## Proof of the simultaneous bound

Use the primitive integer vector
\[
                       x=(D,Dy_1,\ldots,Dy_r).
\]
It is primitive because at any prime of \(D\) some reduced
denominator \(D_i\) attains the full valuation of \(D\); the
corresponding coordinate \(Dy_i\) is a unit there.
At the real place take the independent linear forms
\(x_0,x_1-\theta_1x_0,\ldots,x_r-\theta_rx_0\).
At each prime in \(S\) take \(x_0,\ldots,x_r\).
The
[\(p\)-adic Subspace Theorem, Schlickewei's Theorem 4.1](https://matwbn.icm.edu.pl/ksiazki/aa/aa31/aa3148.pdf)
puts solutions for which the product of the absolute values of these
forms is at most \(\|x\|_\infty^{-\epsilon}\) into finitely many proper
rational linear subspaces. Its hypotheses apply to the \(r+1\)
independent forms at each place. The same square-system statement
is given in
[Evertse's exposition, Theorem 1.1](https://pub.math.leidenuniv.nl/~evertsejh/dio2011-padicsubspace.pdf).

The real vector \((1,\theta_1,\ldots,\theta_r)\) lies in none
of the exceptional subspaces. Indeed membership in a proper
rational subspace implies a nontrivial rational linear relation
among its entries. Dividing \(x\) by \(D\), the vector
\((1,y_1,\ldots,y_r)\) approaches this point, so all sufficiently
close vectors avoid the finitely many subspaces.
On the same bounded neighborhood \(\|x\|_\infty\le C_0D\).
The real contribution to the product is \(D^{r+1}\prod_i e_i\).
The \(S\)-part of \(x_0\) is \(D_S\), and the \(S\)-part of
each \(x_i\) is at least \((D/D_i)_S=D_S/(D_i)_S\).
Hence the product of the forms is at most
\[
                      D_{S^c}^{\,r+1}\prod_i (D_i)_S e_i .
\]
It is greater than \((C_0D)^{-\epsilon}\) on the chosen
neighborhood, proving (1). If exact numerators are available,
write \(y_i=A_i/D_i\) in lowest terms. The stronger (1) uses
\((D_i)_S/(|A_i|)_S\) in place of \((D_i)_S\); this follows
from the exact \(S\)-part of \(x_i\). The limits \(\theta_i\)
are irrational, so each numerator is nonzero sufficiently close
to the reference.

## Denominator budgets from the moment arrays

The arithmetic data for (1) can come from common roots or from
upper bounds on their possible valuations. For a prime \(p\)
sort the common factor's root valuations in the common rational
variable as \(s_1\le\cdots\le s_d\) and put \(U=\sum_i s_i\).
The valuation of the reduced denominator of the coefficient
at distance \(k\) from the leading coefficient in \(H\) is at most
\[
                    \left(-\sum_{i\le k}s_i\right)_+.
\]
When the displayed sum is fractional it can be rounded down,
since the coefficient has integral valuation. The norm
denominator has exact valuation \((-U)_+\). For the coefficient
at distance \(k\) from the leading coefficient in \(B_H\) the
analogous bound is
\[
                      \left(kU-d\sum_{i\le k}s_i\right)_+ .       \tag{3}
\]
These follow by taking minimum valuations in the elementary
symmetric sums. For the common denominator of all the coefficients
of \(B_H\), the maximum of (3) is exact, as proved in the
[projective factor comparison](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators).
The individual bounds are no larger than that common bound.
The [oriented factorial block bounds](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md)
give explicit alternatives for the projective denominator.

For example, let \(L_H\) be the common projective denominator
and choose only projective coordinates, all with rate \(\gamma\).
A direct consequence of (2) is
\[
 (r+\epsilon)\log L_H+\log (L_H)_{S^c}
                 \ \ge\ r\gamma\log(1/\eta)-O(1).               \tag{4}
\]
This uses \(D\mid L_H\). An upper bound on the selected
denominators may also be substituted directly in (2);
the individual \(S\)-parts can give a stronger test than (4).
A uniform violation of the rate in (2) for each algebraic
reference selection with an irrational coordinate excludes
collisions accumulating at those contacts. The estimates allow
large denominator growth, which remains to be controlled for an
unrestricted tail.

## Distinct projective roots and rational relations

For irreducible \(H\), if the projective roots have \(e\)
distinct values, write \(B_H=B_0^{d/e}\) with \(B_0\) monic
irreducible. The
[distinct-root denominator identity](MULTIVARIATE_PROJECTIVE_BINOMIAL_PRIME_RESULTANTS.md#distinct-conjugates-and-denominators)
gives \(L_0=L_H^{e/d}\) for the least common denominator of
the coefficients of \(B_0\). Selecting a sequence with fixed
\(e\), the coefficients of \(B_0\) approach those of a monic
root \(B_{*,0}\) of \(B_{H_*}=B_{*,0}^{d/e}\), with the same
rate. To verify this last assertion, solve successively for
the coefficients of a monic polynomial from the coefficients
of its given power. At each stage the new coefficient is
\(e/d\) times the coefficient of the power minus a polynomial
in preceding coefficients. This is a fixed polynomial
coordinate map and so preserves the error rate on bounded sets.
The limit coefficients are algebraic. Thus one can use their
independent coordinates in (4) with \(L_0\), sometimes obtaining
a larger required exponent for \(L_H\).

Rational affine relations among all reference coordinates give
a separate elementary test. Choose such a relation
\(A_0+\sum_i A_i\theta_i=0\) with integer coefficients, and
let \(D'\) be the common denominator of the participating
actual coordinates. Either the relation holds exactly or
\[
                    1\ \le\ D'\sum_i |A_i|e_i .
\]
Indeed \(D'(A_0+\sum_i A_i y_i)\) is an integer with
absolute value bounded by the right side. In particular, if the
product on the right tends to zero, every rational affine
relation holds eventually. It suffices to check a finite rational
basis of the relations. This also handles a fully rational
reference vector under the indicated denominator growth bound.

## Joint approximation using deep binomial faces

The same theorem can also use genuine \(p\)-adic closeness of
some of the coordinates. Fix a subset \(I_p\) of their indices
for each \(p\in S\). For \(i\in I_p\) fix a finite list of algebraic
targets in \(\overline{\mathbb Q}_p\) and suppose \(y_i\) is
within \(p^{-\kappa_{p,i}}\) of one target in that list, where
\(\kappa_{p,i}\ge0\) can grow as the tail varies. Then the bound
(2) strengthens to
\[
\begin{split}
 (r+1)\log D_{S^c}+\epsilon\log D+
      \sum_{p\in S}\sum_{i\notin I_p} v_p(D_i)\log p
      \ \ge\ &\left(\sum_i\gamma_i\right)\log(1/\eta)\\
             &+\sum_{p\in S}\sum_{i\in I_p}\kappa_{p,i}\log p-O(1).
                                                               \tag{5}
\end{split}
\]
At \(p\), take \(x_0\) and replace \(x_i\) by
\(x_i-\theta_{p,i}x_0\) for the selected indices and targets.
The forms stay linearly independent. There are finitely many
choices of systems, so their union of proper rational exceptional
subspaces is still finite. The real target vector avoids this union
as before. At this prime the new product is
\[
 |D|_p^{\,r+1}\prod_{i\notin I_p}|y_i|_p
                  \prod_{i\in I_p}|y_i-\theta_{p,i}|_p.
\]
Using \(|y_i|_p\le p^{v_p(D_i)}\) and the certified closeness
gives (5). This proof does not require the real and local targets
to be conjugate to one another.

For projective coordinates, one source of local targets is a
fixed prime with a uniform shared slope and a deep binomial
parent face. If the width divides \(\ell\), the
[deep uniform-face identity](MULTIVARIATE_PROJECTIVE_BINOMIAL_PRIME_RESULTANTS.md#a-prime-power-bound-at-a-deep-uniform-face)
gives \(v_p(\beta^\ell-1)\ge\delta_p\) for every projective root
\(\beta\), where \(\delta_p\) is the face depth. Here \(\ell\)
can be fixed as a common multiple of the possible widths. Every
root of unity and every projective root in this identity is a
\(p\)-adic unit. Factoring \(Y^\ell-1\), the sum of the
\(\ell\) nonnegative valuations \(v_p(\beta-\zeta)\) is at
least \(\delta_p\). Choose \(\zeta\) with maximum valuation
\(a=v_p(\beta-\zeta)\), so \(a\ge\delta_p/\ell\). If
\(\beta\) itself is a root of unity the required approximation
is exact. Otherwise, for every \(\zeta'\ne\zeta\), the
ultrametric inequality and maximality give
\(v_p(\beta-\zeta')\le v_p(\zeta-\zeta')\). Indeed when the
two quantities \(a,v_p(\zeta-\zeta')\) are different, the
valuation on the left is their minimum. When they are equal,
maximality bounds the left by \(a\). Summing over \(\zeta'\)
now gives \(\delta_p\le a+v_p(\ell)\), since the product of
the other differences is the derivative \(\ell\zeta^{\ell-1}\).
It follows that one can take
\[
                 \delta_p^*=\max\{\delta_p/\ell,\,
                                          \delta_p-v_p(\ell)\}.
\]
In particular \(\delta_p^*=\delta_p\) for \(p\nmid\ell\).

Select such a root of unity for each projective root with
multiplicity. Each coefficient of the monic polynomial of
the selections is within \(p^{-\delta_p^*}\) of the actual
projective coefficient. This follows by expanding each
elementary symmetric sum and changing one unit factor at a time. The set
of possible torsion polynomials is finite. This argument applies
to \(B_H\) or to \(B_0\) using its distinct roots.
Consequently (5) charges every certified fixed-prime depth
together with the real contact error. In the case of \(r\)
independent projective coordinates with one common error rate
\(\gamma\), and with such a face at every \(p\in S\), it gives
\[
 (r+1+\epsilon)\log L_H
   \ \ge\ r\gamma\log(1/\eta)
     +r\sum_{p\in S}\delta_p^* \log p-O(1).                          \tag{6}
\]
Here \(D\mid L_H\) and all selected projective denominators are
units at the primes of \(S\). The primes in the denominators
may still vary. The sets of places used for (5) and (6) must
be fixed. For moving face primes the separate
[projective resultant bounds](MULTIVARIATE_PROJECTIVE_BINOMIAL_PRIME_RESULTANTS.md)
apply.
