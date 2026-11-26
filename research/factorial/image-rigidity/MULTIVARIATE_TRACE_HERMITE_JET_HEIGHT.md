# Hermite jet height bounds for moving trace fields

When several reciprocal traces of a possible common factor approach
the same algebraic value, the traces force successively small jets
of their rational polynomial at that value. Applying the Subspace
Theorem to those jets gives a denominator bound even if the root
fields and all denominator primes change with the profile.

Fix consecutive ordinary factorial moments and an irreducible
common factor \(H\) with roots \(\alpha_1,\ldots,\alpha_d\).
For one unordered Galois orbit \(O\) of pairs put
\[
 t_{ij}=\alpha_i/\alpha_j+\alpha_j/\alpha_i,\qquad
 T(Y)=\prod_{\{i,j\}\in O}(Y-t_{ij}),\qquad m=|O|.
\]
One can instead take any rational Galois-invariant multiset of
traces, including the full trace polynomial. Write \(D\) for the
least positive denominator of the monic rational polynomial \(T\).
The [trace denominator identity](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md)
gives
\[
             v_p(D)=\sum_{\{i,j\}\in O}
                          |v_p(\alpha_i)-v_p(\alpha_j)|           \tag{1}
\]
when using the orbit without additional multiplicities. The same
formula counts multiplicities in the more general case.

Take a sequence approaching a fixed positive algebraic diffuse
reference selection, and use \(\eta\to0\) for its maximum added
relative gap. By finite extraction we can fix the orbit labels.
Suppose the traces are matched to algebraic limits and there are
\(\nu_a\) traces approaching each distinct limit \(\theta_a\).
Write their positive error exponents in increasing order as
\(\gamma_{a,1}\le\cdots\le\gamma_{a,\nu_a}\), so for the corresponding
trace \(t\) one has \(|t-\theta_a|=O(\eta^{\gamma_{a,i}})\).
All traces stay bounded. These hypotheses follow with some effective
exponents from the
[algebraic contact comparison](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#consequences-at-finite-physical-diffuse-contacts).
The theorem also applies when the trace rates are given directly.
In particular the exponent for a matched quotient tending to
\(1\) or \(-1\) doubles in passing to its reciprocal trace.

Set
\[
 T_*(Y)=\prod_a(Y-\theta_a)^{\nu_a},\qquad
 w_{a,k}=\sum_{i=1}^{\nu_a-k}\gamma_{a,i}\quad(0\le k<\nu_a).
                                                               \tag{2}
\]
There are \(m\) weights \(w_{a,k}\). Write them in increasing order as
\(w_1^\uparrow,\ldots,w_m^\uparrow\). Let \(r_0\) be the dimension
over \(\mathbb Q\) of the span of the \(m+1\) coefficients of \(T_*\),
including its leading coefficient, and let \(r_*=\max(2,r_0)\).

**Theorem.** For every \(\epsilon>0\), sufficiently close to this
reference, either \(T=T_*\) exactly or
\[
             \log D\ \ge\
             \left(\frac{\sum_{\ell=1}^{r_*-1}w_\ell^\uparrow}{r_*}
                                  -\epsilon\right)\log(1/\eta). \tag{3}
\]
No actual root or quotient field is held fixed. When \(r_0=m+1\)
the exponent before \(\epsilon\) is the sum of all the jet weights
divided by \(m+1\). For each cluster of size \(\nu\) with equal rate
\(\gamma\), that cluster contributes
\(\nu(\nu+1)\gamma/2\) to the numerator.

Here is a refined version which can retain large weights even when
\(T_*\) has rational coefficient relations. For a rational linear
subspace \(W\) of the coefficient space containing \(T_*\) after
extension of scalars, let
\[
  S(W)=\max\left\{\sum_{(a,k)\in I}w_{a,k}:
       |I|=\dim W-1,\quad
       \{f\mapsto f^{(k)}(\theta_a)\}_{(a,k)\in I}
       \text{ are independent on }W\right\}.                    \tag{4}
\]
Require \(\dim W\ge2\). The exponent in (3) may be replaced by
\[
                     \min_W \frac{S(W)}{\dim W}-\epsilon .     \tag{5}
\]
The possible ratios form a finite list of weights divided by
integers, although the definition allows infinitely many rational
subspaces. Formula (3) is an explicit bound that avoids computing
the set of possible subspaces.

## Small jets and descent through rational subspaces

For \(0\le k<\nu_a\), expand
\[
 \frac{T^{(k)}(\theta_a)}{k!}
       =\sum_{\substack{I\subset\{1,\dots,m\}\\|I|=m-k}}
                     \prod_{i\in I}(\theta_a-t_i).
\]
Every term includes at least \(\nu_a-k\) factors from the cluster at
\(\theta_a\). Their combined rate is at least \(w_{a,k}\), since that
weight uses the least exponents in the cluster. Other factors are
bounded. Thus for the primitive integer polynomial \(P=DT\),
\[
                         |P^{(k)}(\theta_a)|
                            \le C D \eta^{w_{a,k}}.              \tag{6}
\]
All the jet forms are algebraic linear forms in the coefficients
of \(P\). Their common kernel on the polynomials of degree at most
\(m\) is exactly the line spanned by \(T_*\), because the vanishing
conditions impose divisibility by \(T_*\). This is also the
independence assertion for the \(m\) Hermite evaluations.

The projective coefficient vector of \(P\) has logarithmic height
\(\log D+O(1)\), by primitivity and bounded traces. Its ratio to
the leading coordinate tends to \(T_*\). A rational subspace
containing infinitely many of these vectors therefore contains
\(T_*\) after extension of scalars. Its dimension is at least
\(r_0\). Indeed a rational matrix whose columns form a basis of the
subspace writes every coefficient of \(T_*\) as a rational linear
combination of at most \(\dim W\) algebraic scalars. A subspace of
dimension one containing \(T_*\) consists only of exact projective
attainments, which can be discarded.

Fix such a subspace \(W\) of dimension \(r\) and choose a rational
coordinate isomorphism with \(\mathbb Q^r\). The \(m\) jet forms
have rank \(r-1\) on \(W\). Choose a basis of independent forms
with total weight \(S(W)\). The leading coefficient is independent
of them, since it is one on \(T_*\) and all jets vanish there.
At the real place these \(r\) linear forms, divided by the local
projective norm, have product \(O(\eta^{S(W)})\), by (6) and
the fixed coordinate isomorphism. These coordinates preserve
projective logarithmic height up to a bounded additive constant.

Use the [Subspace Theorem in projective form, equation (1.1)](https://arxiv.org/pdf/math/0408381)
over \(\mathbb Q\) with the real place and these algebraic forms.
They may have complex algebraic coefficients, with the ordinary
complex modulus at the chosen extension of the real place.
If infinitely many nonattainment points satisfy
\(\log D <(S(W)/r-\delta)\log(1/\eta)\) for some \(\delta>0\),
their product of normalized forms is at most
\(H(P)^{-r-\epsilon'}\) for a sufficiently small fixed
\(\epsilon'>0\), eventually absorbing the constants. The Subspace
Theorem places them in finitely many proper subspaces of \(W\).
For one given over an algebraic extension, take the rational span
of its rational points, which remains proper. One such rational
subspace contains infinitely many points, hence also \(T_*\).

If (5) fails infinitely often, start with the full space and apply
this procedure until the dimension can no longer decrease while
containing the limit and infinitely many nonattainment points.
The theorem at that dimension is a contradiction. This proves
(5). Moreover every basis of \(r-1\) distinct weights has total
at least \(\sum_{\ell=1}^{r-1}w_\ell^\uparrow\), whose ratio to
\(r\) is nondecreasing in \(r\) because the weights are ordered
and positive. This proves (3). The exclusion neighborhood from
the Subspace Theorem is qualitative.

A prescribed finite set of local Hermite jet certificates can
be included in the same argument. At each additional finite place
use a basis of independent small linear forms on \(W\), completed
by bounded coordinate forms. For example, if \(D\) is a unit
there and \(T\) has integral traces in clusters with algebraic
targets and certified powers of that prime as errors, the
product expansion in (6) gives the corresponding jet valuations.
Use the \(p\)-adic forms in the same projective theorem. After
division by the local norm, their certified depths add to
\(S(W)\log(1/\eta)\) in the smallness test. A finite number of
target systems may be used by partitioning into subsequences.
Only the places at which such extra forms are supplied must be
fixed. The other denominator primes remain unrestricted.

## A single merged irrational cluster

Suppose \(T_*=(Y-\theta)^m\), and all the rates in this cluster
are equal to \(\gamma\). If \([\mathbb Q(\theta):\mathbb Q]>m\),
the powers through \(\theta^m\) are linearly independent over
\(\mathbb Q\), so \(r_0=m+1\). Formula (3) gives
\[
                        D\ge \eta^{-m\gamma/2+\epsilon}.        \tag{7}
\]
This bound uses no field assumption on the varying trace roots.
There is also no exact attainment in degree \(m\) here. For
\(m>1\) the exponent is strictly larger than
\(m\gamma/(m+1)\), the exponent from using only rate \(\gamma\)
on the normalized ordinary coefficients.

For a general irrational \(\theta\) put
\(l=\min([\mathbb Q(\theta):\mathbb Q],m+1)\). A rational subspace
containing \(T_*\) contains each Galois conjugate
\((Y-\theta')^m\). The jets of order \(0,\ldots,l-2\) at
\(\theta\) are independent on these conjugates. In fact evaluate
the \(l-1\) jets on \(l-1\) distinct conjugates different from
\(\theta\), divide out the nonzero falling factorials in each row,
and factor \((\theta-\theta')^{m-l+2}\) in each column.
The remaining determinant is a Vandermonde determinant.
Thus these \(l-1\) largest weights can always be included in the
basis for \(W\). Completing to any independent basis and using
the smallest possible extra weights gives the explicit exponent
\[
  \gamma\min_{l\le r\le m+1}
     \frac{(l-1)(m-(l-2)/2)+(r-l)(r-l+1)/2}{r}                 \tag{8}
\]
before \(\epsilon\) in (5). If additionally \(P(\theta)\ne0\),
the direct norm bound gives the complementary exponent
\(m\gamma/[\mathbb Q(\theta):\mathbb Q]\): multiply \(P(\theta)\)
by a fixed algebraic denominator to make it integral, use (6)
at \(\theta\), and bound the other conjugate values by \(O(D)\).
One can take the larger of the two rates.

For a rational \(\theta_a\), each jet \(P^{(k)}(\theta_a)\) has
bounded denominator independent of the profile. Formula (6)
then either forces it to vanish exactly or gives the stronger
individual inequality \(\log D\ge w_{a,k}\log(1/\eta)-O(1)\).
The evaluation at \(2\) cannot vanish for the trace polynomial
of an irreducible common factor, because a trace of distinct
nonzero roots is \(2\) only when the two roots are equal.

## Comparison with the elementary discriminant bound

For squarefree \(T\) and a common trace rate \(\gamma\), the nonzero
integer discriminant of \(DT\) itself supplies a lower exponent
\[
                  \frac{\sum_a \nu_a(\nu_a-1)}{2(m-1)}\gamma
                                                               \tag{9}
\]
when \(m\ge2\). This follows because the product of the differences
inside each reference cluster is small, while the other differences
are bounded. The single-cluster discriminant estimate is also given in
[Alain, Section 2](https://arxiv.org/pdf/0706.1022). For a single
cluster (9) gives \(m\gamma/2\).
For several clusters the full-dimension jet exponent exceeds (9)
by \(\gamma(m^2-\sum_a\nu_a^2)/(m^2-1)\).
Thus the new gain over (9) is clearest for clustered traces at
several distinct targets. The jet bounds also apply without assuming
that \(T\) has distinct roots.

As an illustration of the comparison, take a limiting coefficient
vector from \(T_*=(Y-\theta)^3(Y-\theta^4)\) where
\([\mathbb Q(\theta):\mathbb Q]=8\) and \(\theta^4\ne\theta\).
The five coefficients are independent because their largest
powers of \(\theta\) have exponents \(0,4,5,6,7\).
Uniform trace rate \(\gamma\) would give the jet rate
\(7\gamma/5\), compared to \(\gamma\) in (9).
This example is a limiting polynomial shape for the height
comparison.

## Use with factorial blocks

For each trace orbit, the denominator comes from the exact
spread in (1). At a moving prime it can be bounded by maximizing
over shared slope choices after the
[truncated Newton face tests](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).
For the full trace denominator one can use the
[oriented factorial block bounds](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md).
An upper bound below the rate in (3) or (5) excludes accumulation
on the corresponding nonattained reference orbit. In particular
a budget that is too small to permit merged traces can exclude
that branch without first restricting the actual factor fields.
The needed uniform upper bound for arbitrary ordinary tails
remains an arithmetic obligation.
