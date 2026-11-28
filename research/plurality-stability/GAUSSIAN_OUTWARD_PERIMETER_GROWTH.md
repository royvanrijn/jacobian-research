# Perimeter growth from Gaussian threshold dissipation

Let \(A\) be a finite Gaussian partition in \(\mathbb R^d\)
with locally finite perimeter. Assume it satisfies the
stationary winner rule with constant prices at \(0<\rho<1\).
Write \(P(A;B)=\sum_i |D1_{A_i}|(B)\) for its unweighted
perimeter in \(B\). Use the discrete exterior limit \(m\),
defined almost everywhere by the physical-size windows
at \(\rho^{-n}x\) in the
[outward threshold theorem](GAUSSIAN_OUTWARD_THRESHOLD_LIMIT.md).
If \(m\) is fractional on a set of positive Lebesgue
measure, there exist \(C>0\) and \(R_0>0\) such that
\[
 \sum_{n\ge 1}\exp\!\left[
   -C\bigl(1+\rho^{nd}P(A;B_{R_0\rho^{-n}})\bigr)^2
 \right] <\infty.                                           \tag{OP1}
\]
In particular, the condition
\[
            P(A;B_R)=o\!\left(R^d\sqrt{\log\log R}\right)
                  \qquad (R\longrightarrow\infty)           \tag{OP2}
\]
forces the exterior limit to be ordinary. The dilated
indicators then converge strongly in local \(L^1\).

## A frequency interval carrying the fractional mass

Write \(u_n(x)=1_A(\rho^{-n}x)\), as a vector of label
indicators. Fix a smooth nonnegative cutoff \(\chi\le1\),
supported inside \(B_{R_0}\), with
\[
                 \delta=\int\chi^2(1-|m|^2)\,dx>0.
\]
Set \(v_n=\chi u_n\) and \(h=\sigma^2\rho^{2n}\), with
\(\sigma^2=1-\rho^2\). Thus \(h=h_{n+1}\) in (OT4).
By definition of the exterior limit,
\(H_h u_n(x)=H_{\sigma^2}1_A(\rho^{-n}x)\to m(x)\) almost
everywhere, where \(H_h\) is Gaussian convolution of variance
\(h\). The commutator \(H_h(\chi u_n)-\chi H_hu_n\) tends
to zero in \(L^2\). Indeed its \(L^1\) norm is \(O(\sqrt h)\)
by the translation estimate for \(\chi\), and the commutator
is uniformly bounded. Dominated convergence on the support
of \(\chi\) gives \(H_hv_n\to\chi m\) in \(L^2\). Since
\(\|v_n\|_2^2=M:=\int\chi^2\), for all large \(n\)
\[
 \sum_i\int (1-e^{-h|\xi|^2})|\widehat v_{n,i}(\xi)|^2\,d\xi
   =\|v_n\|_2^2-\|H_hv_n\|_2^2\ge\delta/2.
\]
Choose \(\alpha>0\) such that \(\alpha^2 M\le\delta/4\).
The part \(|\xi|\le\alpha/\sqrt h\) contributes at most
\(\delta/4\), so there is at least \(\delta/4\) Fourier
mass at \(|\xi|>\alpha/\sqrt h\).

Write \(V_n=\sum_i|Dv_{n,i}|(\mathbb R^d)\). The boundedness
and BV translation inequalities give, for \(L>0\),
\[
 \begin{aligned}
 \sum_i\int_{|\xi|>L}|\widehat v_{n,i}|^2
 &\le (1-e^{-1/2})^{-1}
                   \sum_i\langle v_{n,i},(I-H_{L^{-2}})v_{n,i}\rangle\\
 &\le C_d V_n/L .                                           \tag{OP3}
 \end{aligned}
\]
For the second inequality write the interaction loss as
one half of the averaged squared \(L^2\) translation
difference. Since the functions take values in \([0,1]\),
the squared \(L^2\) difference is at most the \(L^1\)
difference, which is bounded by the translation distance
times the total variation. Take
\(L_n=\max(2\alpha/\sqrt h,8C_dV_n/\delta)\).
At least \(\delta/8\) Fourier mass lies between
\(\alpha/\sqrt h\) and \(L_n\). Over this interval the
multiplier \(h|\xi|^2e^{-h|\xi|^2}\) is at least
\(\alpha^2 e^{-h L_n^2}\). The first spatial derivative
case of
[(OT4)](GAUSSIAN_OUTWARD_THRESHOLD_LIMIT.md#energy-along-the-outward-orbit)
therefore gives
\[
 \sum_n e^{-h L_n^2}<\infty .
\]
By the BV product rule and scaling,
\(V_n\le \rho^{n(d-1)}P(A;B_{R_0\rho^{-n}})
+\int |\nabla\chi|\,dx\). Substituting this into
\(hL_n^2\) proves (OP1).

If (OP2) holds, the summands in (OP1) are at least
\(n^{-1/2}\) for all large \(n\), a contradiction.
Thus \(m\) is ordinary in this case. Weak-* convergence
of indicator tuples to an indicator tuple gives local
\(L^1\) convergence, by testing each component against
the corresponding limiting indicator and its complement
on bounded sets.

A fractional limit also forces perimeter at least
\(c R_n^d\sqrt{\log\log R_n}\) at a set of geometric
scales \(R_n=R_0\rho^{-n}\) of natural density one.
Indeed, for every \(0<\eta<1\), on the indices with
\(C(1+\rho^{nd}P(A;B_{R_n}))^2<(1-\eta)\log n\)
the summands in (OP1) are at least \(n^{-1+\eta}\).
The number of these indices up to \(N\) is \(O(N^{1-\eta})\).
The constants are allowed to depend on the partition and
the cutoff witnessing the fractional phase.

## Cubic length bound in the active plane

In dimension two, the active network of a stationary
partition of positive masses and locally finite perimeter
has length \(O_{A,\rho}(R^3)\) in \(B_R\). For actual
balanced maximizers the constant is uniform on each compact
positive correlation interval. This sharpens the elementary
analytic upper bound
[in the global interface estimate](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#polynomial-growth-of-the-network).
It still leaves a gap to (OP2).

For each pair let \(F_{ij}\) be its adjusted score contrast,
and put \(F=\prod_{i<j}F_{ij}\). Choose a point \(p\)
off the zero set of \(F\). The Gaussian integral extends
each contrast holomorphically, and in any real unit
direction \(e\) it gives \(|F(p+ze)|\le M_0 e^{B|z|^2}\)
for constants \(M_0,B\) independent of \(e\).
There is a point \(t\in[R,2R]\) with
\[
                  |F(p+te)|\ge \exp(-C R^2)                \tag{OP4}
\]
for \(R\ge1\), uniformly in \(e\). To prove this, truncate
the Taylor series at degree \(m\). The remainder on
\(|z|\le2R\), by Cauchy at radius \(K R\) with \(K=100\),
is at most
\[
 M_0e^{BK^2R^2}\frac{(2/K)^{m+1}}{1-2/K}.
\]
At the \(m+1\) nodes \(t_j=R(1+j/m)\),
the Lagrange formula at zero has sum of absolute
coefficients at most
\[
 (2m)^m\sum_{j=0}^m\frac1{j!(m-j)!}
 \le (4e)^m.
\]
Choose \(m=\lceil C_1 R^2\rceil\) with \(C_1\) large
enough that the Taylor remainder times \((4e)^m\) is
at most \(|F(p)|/2\). The Lagrange bound then proves
(OP4) at one node.

Pick three directions \(e_\ell\) separated by \(120^\circ\),
and choose the corresponding base points
\(p_\ell=p+t_\ell e_\ell\) satisfying (OP4). Their convex hull
contains the triangle with vertices \(p+R e_\ell\):
it contains \(p\), since the offsets have a positive linear
dependence, and then contains the radial segments from
\(p\) to \(p_\ell\). Thus it contains a ball of radius \(R/2\)
and has width at least \(R\) in every direction.
The real score contrasts are bounded uniformly above,
so (OP4) gives the same lower bound, up to constants
in the exponent, for each individual contrast at the bases.
Jensen's formula on complex lines through these bases,
with outer and inner radii \(8R\) and \(4R\), bounds
the number of zeros of each contrast by \(CR^2\) within
the inner radius. For all large \(R\), this includes
every intersection with \(B_R\).

At any regular point \(x\) of an interface in \(B_R\)
and any unit tangent \(\tau\), at least one of the
angular projections \(x\mapsto\arg(x-p_\ell)\) has
Jacobian at least \(c/R\). Indeed
\(\max_\ell|(x-p_\ell)\cdot J\tau|\ge R/2\) by the
width bound, and \(|x-p_\ell|\le4R\) for large \(R\).
The angular projection is smooth away from its base.
Integrate on the rectifiable interface, omitting the
base points, and use coarea and the zero count per
angular fiber. This yields the cubic length bound.
For a sequence of maximizers on a compact correlation
interval use strong maximizer compactness and convergence
of the adjusted scores. The initial \(p\) can be chosen
off all score ties of the limiting partition and gives
a uniform nonzero product bound on the sequence.
The remaining constants are uniform on the correlation
interval.
