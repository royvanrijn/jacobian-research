# First-order asymptotics for unrestricted Gaussian cones

Fix positive masses \(m=(m_i)_{i=1}^k\) and \(0<r<1\). Write
\(p=k-1\) and \(V=V_p(r,m)\) for the unrestricted Gaussian
noise-stability optimum. Let \(C_{p+n}\) be the optimum over
conical partitions in dimension \(p+n\) and \(A_{p,n}\) the
optimum with an \(n\)-dimensional radial spectator block. Then
\[
 C_{p+n}=A_{p,n}+o(n^{-1}).                              \tag{DC1}
\]
In particular the radial expansion (AX4) is also the exact
first-order expansion of the unrestricted conical optimum.
More precisely, if balanced three-cell cones \(F_n\) converge in
Gaussian label distance to a fixed cylindrical planar maximizer
\(A\), then
\[
 \liminf_{n\to\infty} n\bigl(V-S_r(F_n)\bigr)\ge C(A),   \tag{DC2}
\]
where \(C(A)\) is the radial dilation loss from (CH2). It is
strictly positive for a nonsimplex \(A\), and a radial lift of
\(A\) attains the bound at first order. No conditional angular
mass constraint is needed.
[Radial expansion (AX4)](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#exact-leading-loss-for-the-entire-radial-profile-class)
[Dilation loss (CH2)](THREE_CELL_CONICAL_HOMOGENIZATION.md#stability-loss-of-the-balanced-lift-of-a-maximizer)

## Removing the spectator radius

Use the independent correlated pairs \((X,Y)\in(\mathbb R^p)^2\)
and \((Z,Z')\in(\mathbb R^n)^2\). Write
\(L=|Z|/\sqrt n\), \(L'=|Z'|/\sqrt n\), \(\omega=Z/|Z|\),
\(\omega'=Z'/|Z'|\), and \(t=\omega\cdot\omega'\). Represent the
cone by an ordinary profile tuple \(h_\omega(X/L)\). Define an
ordinary partition \(G_n(X,Z)=h_\omega(X)\) and write
\(\bar h(x)=\mathbb E_\omega h_\omega(x)\) and
\(q_\omega(x)=h_\omega(x)-\bar h(x)\). The profiles and mean are
allowed to vary with \(n\). Assume for now that \(F_n\) converges
in label distance to the cylinder over a fixed ordinary active
partition \(A\). Then
\[
 \int \mathbb E_\omega\,|h_\omega(x)-1_A(x)|^2\,d\gamma_p(x)
                                  \longrightarrow0.     \tag{DC3}
\]
Indeed the same integral in the law of \(X/L\) is bounded by a
constant times the Gaussian label distances from \(F_n\) to
\(A(X)\) and from \(A(X)\) to \(A(X/L)\). The second tends to
zero by Gaussian dilation continuity in measure and \(L\to1\)
in probability. The marginal expansion (AX2) changes the
integral by \(O(1/n)\). In particular \(G_n\) also converges to
the cylinder over \(A\).

For an integrable signed density \(D(x,y)\) on the active pair,
define
\[
 \mathcal H_D(t)=\iint D(x,y)
      \left\langle q_\bullet(x),R_tq_\bullet(y)\right\rangle
                                  \,dx\,dy,             \tag{DC4}
\]
where \(R_t\) is fixed-angle averaging on the spectator sphere
and the inner product sums labels and integrates the uniform
angle. Put \(a_n(x)=\|q_\bullet(x)\|_2\). By (DC3), this bounded
function tends to zero in Gaussian measure. Thus for each fixed
\(D\) the quantity
\(\epsilon_n(D)=\iint |D(x,y)|a_n(x)a_n(y)\,dx\,dy\) tends to
zero, by absolute continuity and boundedness. Since \(R_t\) is
a contraction, \(|\mathcal H_D(t)|\le\epsilon_n(D)\) for all \(t\).

On every compact subinterval of \((-1,1)\), \(R_t\) is uniformly
Lipschitz in operator norm for all \(n\ge5\). For degree \(j\ge1\)
its eigenvalue is the expectation of
\(z^j\), with \(z=t+i\sqrt{1-t^2}\,U\) and \(U\) a coordinate
of a uniform vector on \(S^{n-2}\). Here
\(|\partial_t z|\) is uniformly bounded and
\[
 j|z|^{j-1}\le(1-|z|)^{-1}
     \le \frac{2}{(1-t^2)(1-U^2)}.
\]
The expected reciprocal of \(1-U^2\) is \((n-3)/(n-4)\).
Differentiation under the integral gives a uniform eigenvalue
derivative bound, and the orthogonal harmonic decomposition
gives the operator bound. Degree zero is constant. Consequently
\(\mathcal H_D\) on that interval has Lipschitz constant at most
\(C\epsilon_n(D)\).
[Fixed-angle eigenvalue formula](GAUSSIAN_CONE_ANGULAR_SLICE_COMPARISON.md#conditional-angular-comparison)

Conditional on \((t,L,L')\), the angular pair is uniform at fixed
angle \(t\). The active pair has density \(j_r^{L,L'}\) for the
cone and \(j_r\) for \(G_n\), in the profile coordinates. The
angle mean terms separate from \(q\) because each angular marginal
is uniform. Thus the difference of the two stabilities is the
radial correction of \(\bar h\) plus
\(\mathbb E\mathcal H_{j_r^{L,L'}-j_r}(t)\).

Expand the density difference at \((1,1)\) as in (AX2), now
before averaging. Its remainder has expected \(L^1\) norm
\(O(n^{-3/2})\). Every quadratic term paired with \(q\) is
\(o(n^{-1})\), since it has a fixed derivative density and
the associated product of scale deviations has absolute mean
\(O(n^{-1})\). For a linear term with \(L-1\) and a derivative
density \(D\), write
\[
 \mathbb E[(L-1)\mathcal H_D(t)]
  =\mathbb E[L-1]\mathcal H_D(r)
   +\mathbb E[(L-1)(\mathcal H_D(t)-\mathcal H_D(r))].
\]
The first term is \(o(n^{-1})\) by the signed mean scale change
in (AX2). For the second term, on \(|t-r|\le\eta\) use the
Lipschitz bound and Cauchy–Schwarz. Both \(L-1\) and \(t-r\)
have second moments \(O(n^{-1})\). For \(t-r\) this follows
by expressing \(t\) as the empirical cross product divided by
the empirical lengths, expanding on the event where the lengths
are bounded away from zero and using the exponentially small
exception. More generally \(\Pr(|t-r|>\eta)\) is exponentially
small for fixed \(\eta>0\), by concentration of the two squared
lengths and the cross product. This and the contraction bound
handle the complement. The same proof applies to \(L'-1\).
We obtain the comparison
\[
 S_r(F_n)-S_r(G_n)=B_n(\bar h)-S_r(\bar h)+o(n^{-1})
                 =n^{-1}\mathfrak s(A)+o(n^{-1}),        \tag{DC5}
\]
where the last equality uses (AX2) and (DC3). Here \(B_n\)
is the radial lift stability and \(\mathfrak s\) the active
coefficient in (AX4). In this argument the marginal masses of
individual profiles do not enter.

## Correcting the average masses and optimizing

Assume now that \(A\) is an actual maximizer of masses \(m\).
The Gaussian masses of \(G_n\) are those of \(\bar h\) in
\(\gamma_p\), whereas the masses of \(F_n\) are those of \(\bar h\)
in the ratio law. Thus the required correction to \(G_n\) is
\(\Delta m=n^{-1}\mathfrak m(A)+o(n^{-1})\) by (AX2) and (DC3).
For any fixed \(\varepsilon>0\), choose bounded active patches on
each side of contacts along a spanning tree of the essential
closure graph of \(A\), on which the adjusted contact scores
differ by at most \(\varepsilon\). The patches can be chosen to
contain positive Gaussian mass of the appropriate source labels.
In the cylindrical patches for \(G_n\), a reservoir of uniformly
positive mass has the same source label as \(A\) and adjusted
score gap at most \(2\varepsilon\). Indeed \(G_n\) converges in
label distance to the cylinder, and its scores converge in
Gaussian \(L^2\) by contraction of \(T_r\). Route the prescribed
mass change along this tree by disjoint transfers of total mass
\(O(1/n)\). The linear stability change differs from
\(2\beta(A)\cdot\Delta m\) by \(O(\varepsilon/n)\). The quadratic
interaction of the transfers is nonnegative by positivity of
\(T_r\). The balanced repaired partition has stability at most
the unrestricted optimum \(V\). Letting \(\varepsilon\to0\)
after taking a limit and applying (DC5) gives
\[
 \limsup_{n\to\infty} n(S_r(F_n)-V)
       \le \mathfrak s(A)-2\beta(A)\cdot\mathfrak m(A).   \tag{DC6}
\]
For the balanced three-cell planar problem the right side is
\(-C(A)\) by (AX5). The radial lower bound in (AX4), used at
this fixed \(A\), also shows sharpness of (DC2).
[Small mass repair via contact patches](GAUSSIAN_NOISE_PROFILE_DIMENSION_REDUCTION.md#prices-and-small-exact-mass-repairs)

Finally take a sequence of conical maximizers along a subsequence
attaining the upper limit of \(n(C_{p+n}-V)\). The radial
comparison (AX3) makes them unrestricted near-maximizers. By
dimension-free compactness, after passage to a subsequence and
rotating the coordinates they converge in label distance to an
actual maximizer \(A\) in at most \(p=k-1\) active coordinates.
Pad this active space to dimension \(p\) if needed. Apply (DC6)
and take the maximum of its right side over all active
maximizers. Formula (AX4) gives the reverse bound by radial
profiles. This proves (DC1) and identifies the exact first-order
coefficient for unrestricted conical optima. In the balanced
three-cell case it is the minimum of \(C(A)\) over the planar
maximizers for the deficit \(V-C_{2+n}\).
[Compactness modulo orthogonal embeddings](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md#graded-limits-and-original-linear-coordinates)

## Correlation changes on the same scale

In the balanced three-cell problem, fix \(r_0\in(0,1)\) and
\(t\in\mathbb R\). Write \(C_{2+n}(r)\) for the full conical
value and \(V(r)\) for the unrestricted value. The same argument
identifies the conical value at \(r_0+t/n\):
\[
 \lim_{n\to\infty} n\bigl(C_{2+n}(r_0+t/n)-V(r_0)\bigr)
   =\max_{A\in\operatorname{Opt}(r_0)}
          \bigl[t S'_{r_0}(A)-C_{r_0}(A)\bigr].           \tag{DC7}
\]
The derivative here is taken with the partition fixed. The
Taylor remainder in the correlation is uniformly \(O(n^{-2})\)
over all indicator tuples and dimensions, since their total
squared Hermite norm is one and \(\sup_{j\ge2}j(j-1)r^{j-2}\)
is bounded on compact interior correlation intervals. For the
upper bound select maximizing cones on a subsequence giving the
upper limit. They are near-maximizing at \(r_0\) by the uniform
Taylor bound and the radial approximation. Apply compactness
there and then (DC2) at \(r_0\). Their derivatives converge to
the derivative of the limiting maximizer by the uniform norm
bound for the differentiated Ornstein–Uhlenbeck operator and
strong label convergence. For the lower bound use the balanced
radial lift of any fixed maximizer at \(r_0\) and the same Taylor
formula. The maximum is attained by compactness of the planar
maximizer set and continuity of the coefficient in (AX4), with
the price vector determined by its full-rank translation mass
derivative.

At a first tie, a fixed nonsimplex maximizer therefore incurs
its strict dilation loss in (DC7) even if it has a larger
correlation derivative than the simplex. The expression can
select the simplex at \(t=0\) and the nonsimplex at a larger
positive \(t\). A conical plateau at the tie correlation is
compatible with this first-order limit.
