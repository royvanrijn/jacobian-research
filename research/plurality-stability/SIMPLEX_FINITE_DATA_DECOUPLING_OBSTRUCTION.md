# Finite-data obstruction for independent simplex cells

Fix \(k\ge3\), \(d\ge k-1\), and \(0<\rho<1\). Let

\[
 T_i=\{x:\langle v_i-v_j,x\rangle>0\text{ for all }j\ne i\}
\]

be the cells of a centered regular simplex. Write \(P_r\) for the
Gaussian noise operator and \(s_r(E)=\langle1_E,P_r1_E\rangle\) for
the contribution of one set to stability. Fix any open annulus
\(D=\{a<|x|<b\}\) with \(0<a<b<\infty\).

For any finite collection of Gaussian-integrable real-analytic tests,
each cell \(T_i\) has arbitrarily close competitors \(E_i\) in Gaussian
measure that have the same test integrals and mass but strictly greater
\(s_\rho\). Each change consists of adding and removing finitely many
small balls in \(D\). We can preserve the integral of the cell's own
simplex score \(P_\rho1_{T_i}\) as well. In a second version, given
any finite list of correlations \(r_1,\ldots,r_m\in(0,1)\setminus\{\rho\}\)
and any nonnegative integer \(N\), we can preserve each exact \(s_{r_j}(T_i)\)
while preserving all polynomial moments of degree at most \(N\) and the
same simplex score integral.

The independently constructed sets have both overlaps and gaps. This
obstructs sharp upper bounds based only on independent cell constraints,
even if balance, finitely many integrals of the label sum, and finite
samples of stability at other correlations are imposed.

## Finite balls with prescribed linear data

Fix one cell \(T=T_i\). Let \(q=(q_1,\ldots,q_n)\) be a basis for the
span of the analytic tests, \(1\), and \(U=P_\rho1_T\). The score is
real analytic since it is a convolution of a bounded function against
the Gaussian, evaluated at \(\rho x\). For instance its Gaussian
Laplace transform representation converges locally uniformly for
complex \(x\). For \(x\in D\) off the cell boundary, put
\(v(x)=q(x)\) outside \(T\) and \(v(x)=-q(x)\) inside \(T\).

For each nonzero \(\ell\in\mathbb R^n\), some such \(x\) satisfies
\(\ell\cdot v(x)>0\). Otherwise \(F=\ell\cdot q\) is nonnegative in
\(D\cap T\) and nonpositive in \(D\setminus\overline T\). Choose a
facet with supporting hyperplane \(H\). There is a relatively open
patch of the facet in \(D\cap H\) on which all other defining
inequalities for \(T\) are strict. For the facet separating \(i,j\),
one can scale \(v_i+v_j\) into \(D\). The negative of this patch lies
in the interior of the complement, since there is at least one other
strict defining inequality and it reverses. Continuity forces \(F=0\)
on the facet patch, and analyticity forces \(F|_H=0\). Since \(F\)
is nonzero, there is a least order \(h\ge1\) for which its normal
Taylor coefficient along \(H\) is not identically zero. That coefficient
is analytic and nonzero on a dense open subset of \(H\). At a point of
the facet patch where it is nonzero, the sign conditions on opposite
sides of \(H\) force \(h\) to be odd. At a point of the negative patch
where it is nonzero, both sides belong to the complement, so \(h\)
cannot be odd. This proves the assertion.

By compactness of the unit sphere, choose distinct \(x_1,\ldots,x_J\)
so that the open hemispheres \(\{\ell:\ell\cdot v(x_j)>0\}\) cover
it. Then zero is interior to the convex hull of the \(v(x_j)\), and
there are \(w_j>0\) with \(\sum_j w_jv(x_j)=0\). To see strict
positivity, for each \(j\) take a convex representation of a sufficiently
small negative multiple of \(v(x_j)\), add a positive multiple of
\(v(x_j)\) to obtain a representation of zero, and average the \(J\)
representations.

Give the ball centered at \(x_j\) Gaussian volume \(t_j\ge0\) and
either add it to \(T\) or remove it, according to the sign \(e_j\) in
\(v(x_j)=e_jq(x_j)\). The balls are disjoint and contained strictly
inside their respective regions in \(D\) for small volumes. Allow
negative parameters by extending each signed ball integral oddly in
its volume. Let \(G(t)\) be the vector of resulting changes in the
\(q\)-integrals. The derivative of a ball integral with respect to its
Gaussian volume is the Gaussian surface average of the integrand on
the sphere, with limiting value at the center as the radius goes to
zero. Thus \(G\) is \(C^1\) near zero, \(G(0)=0\), and its derivative
has columns \(v(x_j)\), so it has full row rank. The implicit function
theorem gives a curve of exact zeros \(t(u)=uw+o(u)\), \(u\downarrow0\),
since \(w\) is in the kernel of the derivative. All its entries are
positive for small \(u>0\). This gives ordinary sets \(E(u)\) with
\(\gamma(E(u)\triangle T)=\sum_jt_j(u)\).

Put \(g=1_{E(u)}-1_T\). The exact score constraint gives

\[
 s_\rho(E(u))-s_\rho(T)
 =2\langle U,g\rangle+\|P_{\sqrt\rho}g\|_2^2
 =\|P_{\sqrt\rho}g\|_2^2>0.
\]

Every Hermite multiplier of \(P_{\sqrt\rho}\) is positive, which
gives strictness for nonzero \(g\). The gain is at least
\(c\gamma(E(u)\triangle T)^2\) for some \(c>0\) and all small \(u\).
Indeed take a polynomial \(p\) equal to one at \(x_1\) and zero at
all other centers, for example a suitably scaled product of the squared
distances to the other centers. Then
\(\langle p,g\rangle=e_1w_1u+o(u)\), and
\(|\langle p,g\rangle|\le\|P_{\sqrt\rho}^{-1}p\|_2
\|P_{\sqrt\rho}g\|_2\). Here the inverse is defined on the finite
Hermite expansion of \(p\).

## Matching finite samples of stability

Take the linear tests to be polynomials of degree at most \(N\) and
\(U=P_\rho1_T\). Remove repeated sample correlations. Augment the
moment-change map \(G\) by the equations \(s_{r_j}(E)-s_{r_j}(T)=0\).
At zero volume each new row of the derivative is the signed evaluation
of \(2P_{r_j}1_T\) at the ball centers. For the \(C^1\) assertion, use
the signed ball changes \(g_t\) also for negative volume parameters
and the formula

\[
 s_{r_j}(1_T+g_t)-s_{r_j}(T)
 =2\langle P_{r_j}1_T,g_t\rangle+\langle g_t,P_{r_j}g_t\rangle,
\]

where the quadratic form on the left is extended to real functions.
For each ball the derivative of a quadratic term is obtained by
averaging the Mehler kernel in the differentiated variable over its
Gaussian sphere, with the center as the zero-volume limit. In a
self-interaction there are two identical derivative contributions.
Continuity of the kernel on the compact \(\overline D\times\overline D\)
gives continuity of those derivatives, and they vanish at zero total
volume. In particular the quadratic terms do not contribute to the
derivative used to select the centers.

The scores at the distinct correlations \(\rho,r_1,\ldots,r_m\) are
linearly independent modulo polynomials of degree at most \(N\).
Indeed \(1_T\) has nonzero Hermite projections in unbounded degrees.
Otherwise it would equal a polynomial almost everywhere, which is
impossible for a nonconstant indicator. At any such degree \(n>N\)
a relation among the scores would give \(\sum_j a_jr_j^n=0\), with
the term for \(\rho\) included. Dividing by the largest active base
to the power \(n\) and letting \(n\) grow proves that its coefficient
is zero. Thus apply the same analytic sign argument and choice of
centers to the augmented derivative, using the polynomials and all
these scores. The full-rank implicit function argument then preserves
the sampled stabilities exactly. The integral against \(U\) remains
fixed, so the preceding gain estimate at \(\rho\) still holds.

## Pointwise partition constraint

The construction for each cell includes its own score integral. For
the reference simplex, \(U_i=P_\rho1_{T_i}>P_\rho1_{T_j}\) on the
interior of \(T_i\) for \(j\ne i\). Reflection in the hyperplane
orthogonal to \(v_i-v_j\) exchanges \(T_i,T_j\). For \(x,y\) on the
\(T_i\) side of that hyperplane the conditional Gaussian density at
\(y\) is strictly greater than at its reflection, which gives the
score inequality after integration over \(T_i\).

Consequently \(\sum_i\int 1_{T_i}U_i\,d\gamma\) is the integral of
the maximum score. A partition \(E_i\) with the same score integrals
would have to use the unique maximizing label almost everywhere,
forcing \(E_i=T_i\) up to null sets. This contradicts every strict
gain above. Since the constructed sets have total mass one, their
failure to partition gives both an overlap and a gap of positive
measure. The simplex comparison for genuine partitions remains a
separate question.
