# Independent-height windows from separated total-degree clusters

Fix distinct ordinary monomial offsets in any finite dimension
and increasing positive orders \(m_1<\cdots<m_s\), with one order
per offset. Partition the total-degree classes, in increasing
degree order, into consecutive nonempty clusters. Write
\(s_b\) for the number of offsets in cluster \(b\), and
\(p_b=s_1+\cdots+s_b\), with \(p_0=0\). Its assigned orders are
\[
                     I_b=\{m_{p_{b-1}+1},\ldots,m_{p_b}\}.
\]
Within a cluster, assume that all the within-class exponent
differences lie on one real line. The line may differ between
clusters. In two variables the condition is automatic.

For a degree class \(g\) of size \(r_g\), use an arbitrary
polynomial \(Q_g(z)\) of degree at most \(r_g-1\). Write
\(d_g\) for the class degree and \(d_{\min,b}\) for the
smallest degree in cluster \(b\). The **cluster condition**
is that the equations
\[
 \mathbb E\left(\sum_{g\in b}m^{d_g-d_{\min,b}}Q_g(Z_m)\right)^m=0,
 \qquad m\in I_b,\qquad Z_m\sim\mathcal N(0,1/m),
 \tag{1}
\]
have no common nonzero coefficient direction. This uses only the
assigned orders, internal degree differences and class sizes.

**Theorem.** Suppose every cluster satisfies (1). There is a
computable \(D\), depending only on these finite systems and
the prescribed orders, with the following property. If each
gap from the last degree of a cluster to the first degree of
the next is at least \(D\), then every sufficiently large
independent common translation excludes a common factorial
zero at the prescribed orders. This holds over every
characteristic-zero field. The required common-height cutoff
is computable for each fixed offset support. In two variables
there is a finite exact decision for all independent
nonnegative integer translations of that support.

Only the degrees at the cluster boundaries need to separate.
For example, put a single tripled degree class between any
number of singleton and doubled degree classes. Put the lower
classes in one cluster, the triple in its own cluster, and
the upper classes in a third cluster, omitting empty clusters.
For a consecutive window starting at \(N\), the order assigned
to the first row of the triple is \(N+s_{\rm low}\), where
\(s_{\rm low}\) is the number of lower-degree offsets.
The theorem applies at every prescribed start for which that
number is even, and at any odd value through 155. The only
required separations are at the boundaries of the triple.
The threshold may depend on the prescribed orders and on
the internal degree differences of the other clusters.

## 1. The ordered assignment of moment rows

Let \(W_b(m_i)\ge0\) be weights, positive at every row for each
nonzero cluster. Assume that if \(a<b\) and both clusters
are nonzero, then for every adjacent pair of rows,
\[
             \frac{W_b(m_{i+1})}{W_a(m_{i+1})}
                  > R^2\frac{W_b(m_i)}{W_a(m_i)}
 \tag{2}
\]
for some \(R>1\). Let \(S_i\) consist of the cluster labels
with weight at least \(1/R\) of the largest weight at row
\(i\). Assume there is at least one nonzero cluster. If
\(b\in S_i\) and \(a\in S_{i+1}\), then \(b\le a\):
otherwise the ratio on the right of (2) is at least
\(1/R\) before multiplication by \(R^2\). At row \(i+1\)
this would put \(a\) below the allowed fraction of the
maximum. Thus \(\max S_i\le\min S_{i+1}\).

For a boundary after cluster \(b\), let \(T_b\) be the
number of rows with \(\max S_i\le b\), and \(U_b\) the
number with \(\min S_i\le b\). Both sets of rows are
initial intervals. There is at most one row straddling
that boundary, so \(T_b\le U_b\le T_b+1\).
Take the first \(b\) with \(T_b\ge p_b\), which exists
since the last cluster has \(T_b=s\). If \(b>1\), then
\(U_{b-1}\le T_{b-1}+1\le p_{b-1}\), and for \(b=1\)
use \(U_0=0\). Hence \(S_i=\{b\}\) on every row assigned
to this cluster. In particular the same cluster exceeds
all other weights by a factor \(R\) throughout its own
assigned window.

For coefficient vectors \(a_g\), use
\[
             W_b(m)=\max_{g\in b} m^{d_g}\|a_g\|_\infty .
 \tag{3}
\]
When \(m'>m\), the growth ratio for a nonzero cluster lies
between \((m'/m)^{d_{\min,b}}\) and
\((m'/m)^{d_{\max,b}}\). Thus (2) holds if
\((m_{i+1}/m_i)^D>R^2\) for each adjacent pair and the gaps
are at least \(D\). If there is only one row, take \(D=0\).

## 2. Quantitative use of the cluster condition

At each \(m\), divide the phase of a cluster by its weight.
The resulting Gaussian polynomial has every coefficient
within each summand bounded by one, using (3). Therefore
its \(L^m\) norm is bounded by a common \(M_0\) depending
only on the sizes and orders. By (1) and compactness, for
each cluster there is \(\delta_b>0\) such that at least
one of its assigned orders has normalized moment of
absolute value at least \(\delta_b\). To see compactness,
normalize \(\max_{g\in b}\|a_g\|_\infty=1\). In the ratio
defining the normalized phase, subtract \(d_{\min,b}\)
from every degree in the numerator and denominator.
The denominator is positive on the compact coefficient
set. The resulting margin is independent of the absolute
cluster base degree. It can be bounded below by a positive
rational number using real quantifier elimination on the
real and imaginary parts of the rational Gaussian forms
and the finitely many maximum conditions.

Use any common rational lower bound \(\delta>0\).
The same margins with \(\delta/2\), and a uniform norm
bound \(M>M_0\), hold for the exact normalized cluster
phases at all sufficiently large heights, by Section 3.
Choose a rational \(R>1\) so large that for every order \(m\),
\[
              m\bigl(M+(k-1)M/R\bigr)^{m-1}
                         (k-1)M/R < \delta/2 ,
 \tag{4}
\]
where \(k\) is the number of clusters. This and (2)
determine the computable gap threshold. At the cluster
supplied by Section 1, divide the entire phase by its
weight at each assigned row. The other clusters together
have \(L^m\) norm at most \((k-1)M/R\). The power
difference identity and Hölder's inequality bound the
change of the normalized moment by the left side of (4).
At the order where the chosen cluster has margin at least
\(\delta/2\), the entire moment is nonzero. This argument
needs no estimates on joint correlations between clusters.

## 3. Exact normalization and effective cutoff

Write \(L_i>0\) for the shared heights. At order \(m\),
tilt the product-exponential density by the translated
common power. After removing the nonzero common
gamma factors, \(Y_i=X_i/(mL_i)\) are independent
normalized gamma coordinates of shape \(mL_i+1\).
Absorb \(\prod_i L_i^{b_{ji}}\) into the coefficient of
offset \(b_j\); the residual row factor from that column
is \(m^{d_g}Y^{b_j}\). This is the same
[exact normalization for independent heights](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#1-exact-normalization-and-the-case-of-distinct-degrees).

Choose an integral vector \(v_b\ne0\) along the common
direction for nontrivial classes in cluster \(b\).
If all its classes are singletons, any nonzero coordinate
vector works. Write the offsets in one class as
\(b_{g,0}+\lambda_j v_b\), and define
\[
 \sigma_b^2=\sum_i v_{b,i}^2/L_i,\qquad
 Z_{b,m}=\sigma_b^{-1}\sum_i v_{b,i}\log Y_i .
\]
Use the Vandermonde coefficients \(w_{gj\ell}\) satisfying
\(\sum_j w_{gj\ell}\lambda_j^t/t!=\delta_{\ell t}\) for
\(0\le t<r_g\). The exact basis elements in the phase are
\[
      Y^{b_{g,0}}\sigma_b^{-\ell}
                      \sum_j w_{gj\ell}e^{\lambda_j\sigma_b Z_{b,m}},
               \qquad 0\le\ell<r_g .
 \tag{5}
\]
This is an invertible change of column coefficients,
independent of \(m\). The weights in (3) refer to the
coefficients in (5).

For fixed \(m\), the variables \(\sqrt{L_i}\log Y_i\)
converge in moments to independent centered Gaussians
of variance \(1/m\). Their exponential moments at any
fixed real argument are uniformly bounded when all
heights are large, as in the
[log-gamma estimate for centered lifting](SFC_HOMOGENEOUS_CENTERED_HIGHER_LIFT.md#1-centered-phase-and-change-of-law).
The coefficients \(v_{b,i}/(\sigma_b\sqrt{L_i})\)
of these variables in \(Z_{b,m}\) have squares summing
to one. Therefore \(Z_{b,m}\) has uniform moment
convergence to \(Z_m\), by expanding each fixed moment,
and it has uniformly bounded exponential moments.
The exponential Taylor remainder in (5) shows that the
divided sum differs from \(Z_{b,m}^\ell\) by a term
tending to zero in every fixed moment norm.
The prefactor \(Y^{b_{g,0}}\) tends to one in those norms,
again by the exponential-moment bound. Expanding the
\(m\)-th power proves convergence of the normalized
cluster moments uniformly on the bounded row coefficients
of Section 2, with an \(L^m\) bound \(M>M_0\) once
the common-height cutoff is large enough. No constraint on
relative heights is used.

For effectiveness, the original forms after removing
\(\prod_i\Gamma(mL_i+1)\) have polynomial coefficients in
the heights, since every offset in the moment is a
nonnegative integer and the gamma ratios are rising
factorials. They represent the same equations for real
nonnegative heights. Search for an integer corner cutoff
by real quantifier elimination on their complex coefficient
equations. The proof of the theorem ensures termination.
In two variables the remaining integer heights lie on
finitely many horizontal or vertical fibers. The
[univariate resultant decision on these fibers](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#4-effective-search-in-two-independent-heights)
then applies. At integer heights all forms are rational.
If a common projective direction exists over a
characteristic-zero extension, the coefficient equations
in a chart also have a solution over
\(\overline{\mathbb Q}\), and hence over \(\mathbb C\).

## 4. Applications to repeated degree classes

A cluster of singleton and doubled classes satisfies (1)
on every assigned consecutive window. This is the
[Gaussian matrix test for any number of doubled classes](SFC_MULTIPLE_DOUBLED_DEGREE_INDEPENDENT_HEIGHTS.md#2-the-matrix-at-the-independent-height-limit)
since a window of its size contains at least as many
even orders as there are pairs. A homogeneous tripled
class satisfies (1) at every even assigned start by the
[even-start quadratic Gaussian theorem](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md#even-orders-and-distinct-valuations),
and at each odd assigned start through 155 by the
[exact odd-start certificates](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#4-a-verified-prefix-of-odd-windows).
This proves the stated example with a triple surrounded
by classes of size at most two. Several separated tripled
classes are allowed if each has a certified assigned window.
Larger multiplicities can be used whenever their assigned
cluster systems satisfy (1). For an unsupported cluster,
the theorem leaves its own shorter Gaussian system as
the remaining obligation.

## 5. A finite decision with no bound on the degrees

**Theorem.** Fix a finite increasing list of positive orders.
There is a finite exact decision for the following property:
for every contiguous sublist of \(r\) orders, the bivariate
graded Gaussian system is nonvanishing on projective
coefficient space for **every** ordered degree profile of
\(r\) columns. On success one obtains a positive computable
margin, uniform over all these degree profiles. On failure
one obtains a degree profile, a sublist and a Gaussian
coefficient direction witnessing failure.

Here a degree profile specifies the strictly increasing
degrees \(d_g\) and their positive class sizes \(r_g\).
The Gaussian system is (1) with every class in the
sublist included. Subtract the smallest degree, which
only removes a nonzero order factor. Any resulting abstract
profile can be realized by bivariate offsets after adding a
common degree to every class, so that each degree has enough
distinct monomials for its class size.

Proceed by increasing sublist length. For length one the
normalized Gaussian phase is a constant of absolute value
one, so the uniform margin is one. Suppose that every
shorter sublist has passed and its margin is available.
Take a positive common lower bound \(\delta\) for the shorter
margins. Use the pure Gaussian versions of (3) and (4) and
choose \(R\) with at most one other cluster in (4). Determine
\(D\) from the adjacent ratios of orders. These bounds can be
chosen for each sublist, or simultaneously for every sublist
of the next length.

If any adjacent degree gap is at least \(D\), split the
profile at that gap. The assigned orders of both parts
are shorter contiguous sublists. Sections 1–2 imply that
at an assigned order the whole normalized phase has moment
of absolute value at least \(\delta/2\). Indeed at the
selected row the weight of its selected cluster is the
global maximum, so the normalization for the whole profile
is exactly the one used in the perturbation estimate.
The lower bound is independent of the within-cluster
degree gaps.

The other profiles have all gaps less than \(D\). They
form a finite set after subtracting the smallest degree.
Test their rational Gaussian equations by exact algebraic
elimination. If a gate has a root, a point can be recovered
over the algebraic numbers. Otherwise bound the normalized
margin for each such profile by the compactness calculation
of Section 2. The minimum of these finitely many positive
bounds and \(\delta/2\) is a uniform margin for the new
length. This completes the decision. Its positive outcome
gives the finite independent-height decision for every
individual bivariate offset support on the tested sublists,
with no total-degree bound on the input.

This decision concerns Gaussian nonvanishing. A Gaussian
root alone does not give a factorial counterexample. The
odd quadratic problem already occurs in the single-degree
profiles. The reduction shows that certification across
arbitrarily large gaps can be obtained from shorter
sublists, leaving finitely many bounded gaps at each
new support size and fixed order list.

## 6. Initial windows through four columns

**Theorem.** For every bivariate offset support of \(s\le4\)
columns, the graded Gaussian gate at orders \(1,\ldots,s\)
has no common projective direction. Thus for every such
fixed support there is a finite exact decision for existence
of a complex common zero at the initial window under any
independent nonnegative integer translations. A nonvanishing
certificate holds over every characteristic-zero field.

Profiles with only singleton and doubled classes are
covered by the doubled-class test of Section 4. The
homogeneous triple is also covered there. For four columns
the other cases are a triple and a singleton, and the
homogeneous size-four class.

For a triple and a singleton, the assigned triple window
is \(1,2,3\) or \(2,3,4\). The exact checker below bounds
their margins using homogeneous ideal identities expressing
each pure power of degree \(\sum_{m\in I_b}(m-1)+1\)
in the ideal generated by the corresponding three
Gaussian equations. For each pure power, expand the
multipliers in monomials. On \(\|a\|_\infty=1\),
one of the pure powers has absolute value one.
The reciprocal of the largest sum of absolute multiplier
coefficients is therefore a margin for the three moments.
The calculation gives \(4/121\) and
\(2939770633/1033819864741\), respectively.

In the normalized Gaussian phases, the \(L^m\) norms of
\(Z_m\) and \(Z_m^2\) are at most \(2\) for \(1\le m\le4\),
since \(\mathbb E Z_m^{2m}\le2^m\). The normalized
cluster norms are at most five. With \(R>1\), the bound
on the left of (4) for two clusters and \(m\le4\) is at most
\(20000/R\). Using the smaller margin and
\(R=14066674\) gives a bound strictly below half the margin.
The smallest ratio of adjacent orders is \(4/3\), and
\((4/3)^{115}>R^2\). Thus every gap at least \(115\)
is covered by applying Sections 1–2 directly to the Gaussian system.

It remains to check both orders of the triple and singleton
at gaps \(1,\ldots,114\), and the homogeneous four-column
gate. The [exact checker](verify_sfc_initial_four_cluster_windows.py)
computes the rational ideal identities above and checks
all 229 remaining homogeneous systems modulo \(65537\).
For each, a standard basis gives affine dimension zero,
so the homogeneous common locus is only the origin.
A characteristic-zero projective common zero would
reduce to a projective one at this prime. This verifies
the required exclusion. Replay from the repository root
with Python 3, SymPy and Singular using

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_initial_four_cluster_windows.py
~~~
