# An optimal split with a fixed third cell can require nonaffine sets

The [affine separation reduction](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md)
starts with a line separating two cells of an actual three-cell
maximizer. Maximizing the stability of those two labels while holding
the third cell fixed cannot supply this premise in general. A bounded
convex third cell with no reflection symmetry excludes stationary
balanced affine splits at every \(0<r<1\). There is also
a radially symmetric union of Gaussian mass \(2/3\) for which every
maximizing split into sets of Gaussian mass \(1/3\) is nonaffine on an
explicit interval of correlations adjacent to one.

## A convex fixed cell and the reflection criterion

Let \(K\) be a compact convex body of planar Gaussian mass \(1/3\),
and split \(U=\mathbb R^2\setminus K\) into two cells of Gaussian
mass \(1/3\). For any \(0<r<1\), if this split is maximizing and
affine, its line passes through the Gaussian origin and is a
reflection axis of \(K\). The same necessity holds for a split
stationary under balanced compact variations of its open interface
in \(U\).

Choose orthonormal coordinates \((s,u)\) so the line is \(u=c\),
and let \(d=1_{A_1}-1_{A_2}=1_U\operatorname{sgn}(u-c)\). The
first variation of the sum of the two cells' stabilities gives
a constant value of \(T_r d\) on each open interface ray outside
\(K\). In detail, any smooth compact speed on one such ray with
zero Gaussian surface integral can be realized by a local normal
variation with exact balance, using a second compact speed of
nonzero integral to repair volume at order two. The stability
derivative, up to sign, is twice the integral of the speed times
\(T_r d\). Thus the score is constant on the ray. Gaussian smoothing
makes it real analytic, so it is constant on the entire line. As
\(s\to+\infty\) at \(u=c\), the conditional probability of \(K\)
tends to zero since \(K\) is bounded. The constant therefore equals
the full half-plane contrast \(T_r\operatorname{sgn}(u-c)\) on
the line. Consequently
\[
                T_r\bigl[1_K\operatorname{sgn}(u-c)\bigr](s,c)=0
                \qquad(s\in\mathbb R).                         \tag{FU3}
\]

Write \(\sigma=\sqrt{1-r^2}\) and let \(p_\sigma\) be the centered
normal density of standard deviation \(\sigma\). At observation
\((s,c)\) the conditional mean is \((rs,rc)\) and the coordinates
are independent. The left side of (FU3) is the convolution of
\(p_\sigma\) with
\[
       G(t)=\int 1_K(t,u)\operatorname{sgn}(u-c)
                                    p_\sigma(u-rc)\,du
\]
evaluated at \(rs\). Gaussian convolution is injective on this
compactly supported integrable function, by its nonzero Fourier
multiplier. Hence \(G=0\) almost everywhere. Each nondegenerate
section of the convex body is an interval \([a(t),b(t)]\).
For almost every such section, \(a(t)<c<b(t)\) and
\[
   \int_0^{c-a(t)}p_\sigma((1-r)c-v)\,dv
       =\int_0^{b(t)-c}p_\sigma((1-r)c+v)\,dv.                \tag{FU4}
\]

If \(c>0\), the integrand on the left is strictly larger
than the one on the right for each \(v>0\). Thus
\(b(t)-c>c-a(t)\) and \(a(t)+b(t)>2c>0\). On this section,
the lower label of \(U\) has ordinary Gaussian fiber mass
\(\Phi(a(t))\), which is strictly larger than the upper mass
\(\bar\Phi(b(t))\). The same strict ordering holds where the
section of \(K\) is empty, because the split of that fiber is
at \(c>0\). Singleton sections are immaterial. Integration
contradicts balance. Reflecting \(u\) excludes \(c<0\). With
\(c=0\), strict positivity and evenness of \(p_\sigma\) in
(FU4) give \(b(t)=-a(t)\) almost everywhere. This is equality
of \(K\) and its reflection up to null sets, and hence exact
equality since a convex body is the closure of its interior.

For example, start with the scalene triangle with vertices
\((-1,-1)\), \((2,-1)\), and \((-1,1)\). The origin is in its
interior, so scaling it concentrically varies its Gaussian mass
continuously and strictly from zero to one. There is a unique
positive scale giving mass \(1/3\). The scaled triangle has no
reflection axis, since a reflection of a nondegenerate triangle
must fix one vertex and switch the other two, making two side
lengths equal. The maximizing split of its complement exists by
the compactness argument below, so every maximizing split is
nonaffine. This holds at every \(0<r<1\).

## A radially symmetric example near one

Let \(\gamma_2\) be standard planar Gaussian measure and set
\[
 a=\sqrt{2\log(3/2)},\qquad b=\sqrt{2\log2},\qquad
 c=\sqrt{2\log6},\qquad \delta=b-a>0.
\]
Define \(D=B_a(0)\), \(B=\{b<|x|<c\}\), and \(U=D\cup B\).
The radial Gaussian tail \(\gamma_2(|x|>s)=e^{-s^2/2}\) gives
\(\gamma_2(D)=\gamma_2(B)=1/3\). For a split \(U=A_1\sqcup A_2\)
of these same masses, set \(A_3=U^c\). Write
\(X,Z\) for independent standard planar Gaussians,
\(Y=rX+\sigma Z\) for \(\sigma=\sqrt{1-r^2}\), and
\(C_r(E,F)=\Pr[X\in E,Y\in F]\). Then
\[
 S_r(A)=C_r(A_3,A_3)+C_r(U,U)-2C_r(A_1,A_2),              \tag{FU1}
\]
where symmetry of \((X,Y)\) accounts for the factor two. Thus
the third label and the union contribute fixed terms.

For the split into \(D\) and \(B\), whenever \(X\in D\) and
\(Y\in B\) we have \(|Y-rX|\ge b-ra>\delta\). Consequently
\[
 C_r(D,B)\le \tfrac13 e^{-\delta^2/(2\sigma^2)}.           \tag{FU2}
\]
Any balanced affine split of \(U\) is centered. The projection
of \(1_U\,d\gamma_2\) onto any unit direction has a symmetric
cumulative distribution, which is strictly increasing on
\((-c,c)\). The annulus already gives positive projection
mass on every subinterval of \((-c,c)\). Hence the unique
median is zero, and rotation invariance lets us take
\(A_1=U\cap\{x_1>0\}\) and \(A_2=U\cap\{x_1<0\}\).

For \(\sigma\le a/4\) consider the event given by
\[
 0<X_1<\sigma/2,\quad |X_2|<a/4,\quad
 -1<Z_1<-3/4,\quad |Z_2|<1.
\]
Here \(|X_1|<a/8\) and \(|X_2|<a/4\), so \(X\in D\).
The first coordinate of \(Y\) lies between \(-\sigma\) and
\(-\sigma/4\), and \(|Y_2|<a/4+\sigma\le a/2\), so \(Y\in D\).
This event therefore contributes to \(C_r(A_1,A_2)\). Let
\(\Phi\) and \(\phi\) be the one-dimensional standard Gaussian
cumulative distribution and density, and put
\[
 k_0=\tfrac12 \phi(a/8)\bigl(2\Phi(a/4)-1\bigr)
           \bigl(\Phi(-3/4)-\Phi(-1)\bigr)\bigl(2\Phi(1)-1\bigr)>0.
\]
Independence of the four coordinates gives a lower bound
\(k_0\sigma\) on the event probability. For
\[\sigma\le\min\{a/4,\,3k_0\delta^2/4\},\]
(FU2) and \(e^{-t}<1/t\) for \(t>0\) give
\[
 C_r(D,B)< \frac{2\sigma^2}{3\delta^2}
         \le \tfrac12 k_0\sigma < C_r(A_1,A_2).
\]
Equation (FU1) now makes the nonaffine split strictly better
than every balanced affine split.

A maximizing split of this fixed union exists. Relax the first
indicator to \(f\) with \(0\le f\le1_U\) and \(\int f\,d\gamma_2=1/3\),
and let the second assignment be \(1_U-f\). The feasible set
is weakly compact in \(L^2(\gamma_2)\), and the functional is
weakly continuous by compactness of \(T_r\). If a maximizer
has a positive-measure fractional region, restrict to a subset
where both assignments are bounded away from zero, split it
into equal-mass pieces and perturb by opposite signed amounts.
The functional is strictly convex along this line because all
Hermite eigenvalues of \(T_r\) are positive. This contradicts
maximality at the midpoint. Hence the maximum is attained by an
ordinary split, and the strict comparison excludes every affine
split of the given union from the maximum.

The example concerns the fixed-union problem (FU1). The common
choice of the third cell in a full three-label maximization is
an additional optimization variable.
