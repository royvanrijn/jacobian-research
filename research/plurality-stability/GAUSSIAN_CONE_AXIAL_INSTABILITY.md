# Axial instability and radial block limits for Gaussian cones

For every \(0<r<1\) and every choice of positive masses for at least
three labels, an optimal Gaussian conical partition cannot depend on
one latitude alone. A single rotation of a latitude boundary gives an
exactly mass-preserving improvement. For three labels, the best latitude
partition is separated from the unrestricted conical optimum by a
positive gap uniform in the ambient dimension and on compact sets of
positive correlations and mass triples.

For a fixed \(p\)-dimensional block, partitions depending on that block
divided by the normalized radius of an \(n\)-dimensional spectator block
have a uniform first-order stability expansion as \(n\to\infty\).
In the balanced planar three-cell case this gives the exact leading
loss for the optimum over this entire conical class.

## A single latitude rotation

Write \(P\) for the Gaussian angular operator on \(S^{d-1}\), \(d\ge2\),
at correlation \(r\). Its kernel \(K(u\cdot v)\) is a positive mixture
of \(e^{a(u\cdot v)}\) with \(a>0\), obtained by integrating the two
radii in the joint Gaussian density. In particular \(K'>0\) and
\(K''>0\). Its kernel is analytic on a neighborhood of \([-1,1]\),
since the real part of the radial quadratic form stays positive definite
for complex inner products close to this interval. It is strictly
positive definite as an operator, by identifying its quadratic form
with the Gaussian noise form on homogeneous extensions.
[Angular operator](GAUSSIAN_CONE_CONFORMAL_TRACE.md#angular-operator-and-conformal-derivative)

Suppose an optimal conical partition depends on the latitude
\(t=u\cdot v\) for a fixed unit vector \(v\). Its indicators choose
winners of the angular scores \(P u_i-\beta_i\), with constant mass
prices. Each score is analytic as a function of \(t\in(-1,1)\). Two
adjusted scores cannot coincide identically, by injectivity of \(P\),
the fact that \(P1=1\), and the positive masses of the cells. Thus the
labels form bands with transitions locally discrete in \((-1,1)\).
[Angular score rule](THREE_CELL_CONICAL_HOMOGENIZATION.md#equality-of-consecutive-conical-optima)

Choose a transition at \(x=0\) if present. Otherwise reverse the axis
if necessary and choose the least positive transition \(x\). The left
and right bands at \(x\) have labels \(i,j\), respectively. If \(x>0\),
the whole interval \((0,x)\) has label \(i\) up to null sets. Put
\(C=\{t>x\}\), \(g=u_j-u_i\), and \(s=2\,1_C-1\).

The axial integral kernel is
\[
 p(t,z)=\mathbb E K\bigl(tz+\sqrt{1-t^2}\sqrt{1-z^2}\,U\bigr),
\]
where \(U\) is the inner product of two independent uniform
directions in the \((d-1)\)-dimensional space perpendicular to \(v\).
Its law is symmetric and supported on \([-1,1]\). Integration in
\(z\) uses the uniform spherical latitude law. For \(z>x\),
differentiation in \(t\) at \(t=x\) gives a positive result.
Indeed the multiplier of \(K'\) is
\[
 z-\frac{x\sqrt{1-z^2}}{\sqrt{1-x^2}}U>0
\]
even at \(U=1\). For \(z<0\), the term \(z\mathbb E K'\) is strictly
negative, and the remaining term is nonpositive because
\(\mathbb E[U K'(b+aU)]\ge0\) for \(a\ge0\) by symmetry and
monotonicity of \(K'\). On \(z>x\) we have \(g-s\le0\), on \(z<0\)
we have \(g-s\ge0\), and on \((0,x)\) we have \(g=s\). At least
one third label has positive mass outside \((0,x)\). This gives
\[
                         (P(g-s))'(x)<0.                 \tag{AX1}
\]

Rotate \(C\) by a small nontrivial rotation taking \(v\) to a nearby
axis. Write \(\Delta=1_{C_{\mathrm{rot}}}-1_C\). Its support is in the
two bands meeting at \(x\), so \(u_j+\Delta,u_i-\Delta\) are feasible
indicators of the same masses. Their stability change is
\(2\langle\Delta,Pg\rangle+2\langle\Delta,P\Delta\rangle\). The same
change with \(g\) replaced by \(s\) is zero by rotational invariance
of the binary cap partition. The actual change therefore equals
\(2\langle\Delta,P(g-s)\rangle\). The function \(P(g-s)\) is
strictly decreasing near \(x\) by (AX1). The rotation adds mass below
the original latitude and removes an equal mass above it, so the
change is strictly positive.

The derivative comparison is the angular counterpart of the
nearest-transition argument for ordinary one-dimensional Gaussian
partitions. It uses positivity and monotonicity of \(K'\) and local
discreteness of the bands. The Gaussian angular operator supplies
all these hypotheses.
[Ordinary one-dimensional obstruction](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

## A uniform expansion for radial spectator blocks

Take integers \(p,n\ge1\). Let \(X,Y\in\mathbb R^p\) and
\(Z,Z'\in\mathbb R^n\) be independent Gaussian pairs, each with
standard marginals and correlation \(r\). Put \(L=|Z|/\sqrt n\),
\(L'=|Z'|/\sqrt n\). A conical partition invariant under the orthogonal
group on the second block has profile \(h(X/L)\) for an ordinary
partition \(h\) in \(\mathbb R^p\). Conversely every such profile
defines a partition in this class. For \(p=1\) these are exactly the
latitude partitions about the first axis.

Denote the standard marginal density in \(\mathbb R^p\) by \(\phi\)
and the joint Gaussian density at correlation \(r\) by \(j_r\).
Define \(\phi^a(x)=a^p\phi(ax)\) and
\(j_r^{a,b}(x,y)=(ab)^p j_r(ax,by)\). The densities of \(X/L\) and
\((X/L,Y/L')\) have the \(L^1\) expansions
\[
\begin{split}
 \mathbb E\phi^L&=\phi+n^{-1}M+O(n^{-3/2}),\\
 \mathbb E j_r^{L,L'}&=j_r+n^{-1}J_r+O(n^{-3/2}),\\
 M&=\tfrac14(\partial_{aa}-\partial_a)\phi^a\big|_{a=1},\\
 J_r&=\tfrac14(\partial_{aa}+\partial_{bb}-\partial_a-\partial_b
                 +2r^2\partial_{ab})j_r^{a,b}\big|_{a=b=1}.
\end{split}                                                   \tag{AX2}
\]
The bounds are uniform in \(r\) on compact subintervals of \((0,1)\)
for fixed \(p\). In fact with \(Q=L^2-1\) and \(Q'=(L')^2-1\), one
has \(\mathbb E Q^2=2/n\), \(\mathbb E QQ'=2r^2/n\), and
\(\mathbb E Q^4=O(n^{-2})\). The identity
\(L-1=Q/2-(L-1)^2/2\), with \(|L-1|\le |Q|\), gives
\[
 \mathbb E(L-1)=-\frac1{4n}+O(n^{-3/2}),\qquad
 \mathbb E(L-1)^2=\frac1{2n}+O(n^{-3/2}),\qquad
 \mathbb E[(L-1)(L'-1)]=\frac{r^2}{2n}+O(n^{-3/2}).
\]
Taylor expansion in \(L,L'\in[1/2,3/2]\) has an \(L^1\) remainder
bounded by the third absolute centered scale moments, which are
\(O(n^{-3/2})\) by the fourth moment bound. The needed density
derivatives are polynomials times Gaussian densities, uniformly
integrable over this square and the compact correlation interval.
The event outside the square has probability \(O(n^{-2})\). There
the densities themselves have unit integral, and the omitted first
and second centered scale moments are \(O(n^{-2})\) by the fourth
moment bound and Hölder's inequality. This proves (AX2).

Let \(A_{p,n}(r,m)\) be the conical optimum in this class at mass vector
\(m\), and \(V_p(r,m)\) the ordinary \(p\)-dimensional Gaussian optimum.
The relaxed invariant class is weakly compact on the sphere and the
angular operator is compact. A fractional maximum can be split in two
opposite balanced directions in the nonatomic profile law. Strict
positive definiteness therefore gives an ordinary maximum.
The uniform comparison
\[
                    A_{p,n}(r,m)=V_p(r,m)+O(n^{-1})         \tag{AX3}
\]
follows by applying (AX2) to the agreement event and to all cell
volumes of an arbitrary profile. For either inequality, transfer a
total mass equal to half the sum of the absolute label discrepancies
under the target marginal law. The marginal laws are nonatomic, and
the change in agreement is at most twice the transferred mass.
This also gives a uniform \(O(d^{-1})\) approximation of \(V_2(r,m)\)
by the unrestricted conical value in ambient dimension \(d\) for
three labels, using planar scalar dimension reduction for the upper
bound.
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)

For three positive masses \(V_1(r,m)<V_2(r,m)\), since a one-dimensional
maximizer exists, and its extension cannot be an unrestricted planar
maximizer by (G8). This gap is uniformly positive on compact
parameter sets by continuity of \(V_p\). The comparison in (AX3)
for \(p=1,2\) yields a positive gap between the latitude class and
the full conical optimum in all sufficiently large dimensions,
uniformly on the compact set. For each smaller dimension the latitude
class is weakly closed in the relaxed indicator space, and \(P\) is
compact and strictly positive definite. It attains an ordinary
maximum, by rounding a fractional optimum through a two-sided
balanced split in the nonatomic latitude law. The same is true for
the full conical optimum with the spherical law. The values are
continuous in masses by direct transfers and in correlation by
operator-norm continuity of the Gaussian noise operator. The rotation
above then gives a strictly positive minimum gap over the remaining
dimensions and the compact set.

## Exact leading loss for the entire radial profile class

Fix \(p,r,m\). For an ordinary Gaussian maximizer \(h=(h_i)\) of
masses \(m_i>0\), write \(\beta_i(h)\) for its score prices, so that
\(h\) chooses the winners of \(T_r h_i-\beta_i(h)\). Define
\[
 \mathfrak m_i(h)=\int h_i(x)M(x)\,dx,\qquad
 \mathfrak s(h)=\sum_i\iint h_i(x)h_i(y)J_r(x,y)\,dx\,dy .
\]
The exact asymptotic optimum in (AX3) is
\[
 \lim_{n\to\infty} n\bigl(A_{p,n}(r,m)-V_p(r,m)\bigr)
 =\max_{h\in\operatorname{Opt}_p(r,m)}
       \left[\mathfrak s(h)-2\sum_i\beta_i(h)\mathfrak m_i(h)\right].
                                                               \tag{AX4}
\]
The expression is independent of the common additive constant in the
prices because \(\sum_i\mathfrak m_i(h)=0\).

For the upper bound choose a maximizing profile \(h^{(n)}\) in the
ratio law. By Gaussian weak compactness pass to a limit \(h\) in a
subsequence attaining the limsup. By (AX2), its ordinary Gaussian
masses are \(m\). By compactness of \(T_r\) on the fixed-dimensional
Gaussian space and (AX3), its stability is \(V_p(r,m)\). Strict
positive definiteness and a two-sided balanced splitting exclude
fractional maximizers. Consequently the weak convergence to \(h\)
is strong in \(L^2\), since all tuples are ordinary and have squared
norm one. The functionals in (AX2) converge along this strong
convergence by integrability of the coefficient densities and
boundedness of the indicators.

The Gaussian mass correction from \(h^{(n)}\) to \(m\) is
\(n^{-1}\mathfrak m(h)+o(n^{-1})\). Route it along a spanning tree
of the essential closure graph of \(h\). The graph is connected,
and adjacent labels have tied adjusted scores at a contact point.
For each fixed \(\eta>0\) use bounded patches around these contacts
where the adjusted score differences have absolute value below
\(\eta\). Strong convergence supplies a uniformly positive reservoir
of each source label of \(h^{(n)}\) in its chosen patch. Gaussian
smoothing gives uniform convergence of the scores on these patches.
The linear stability change differs from twice the price of the mass
correction by \(O(\eta/n)+o(n^{-1})\), and the interaction of the
transfers is \(O(n^{-2})\) since the supports are bounded. The
corrected partition is feasible for \(V_p\). Taking the limit and
then \(\eta\downarrow0\) proves the upper bound in (AX4).
[Price-exact mass repairs](GAUSSIAN_NOISE_PROFILE_DIMENSION_REDUCTION.md#prices-and-small-exact-mass-repairs)

For the lower bound start with any fixed \(h\in\operatorname{Opt}_p\)
and use the same patches to correct its ratio masses. By (AX2) the
modified profile converges strongly to \(h\), and its Gaussian mass
change is \(-n^{-1}\mathfrak m(h)+o(n^{-1})\). Its Gaussian agreement
change therefore differs from \(-2n^{-1}\beta(h)\cdot\mathfrak m(h)\)
by \(O(\eta/n)+o(n^{-1})\). Formula (AX2) adds
\(n^{-1}\mathfrak s(h)+o(n^{-1})\) to pass to ratio agreement.
Again let \(\eta\downarrow0\) after taking the limit. The maximizer
selected in the upper bound attains the right side of (AX4).

When \(p=2\) and the three masses are equal, the translated lift
in (CH2) realizes the coefficient inside the maximum in (AX4)
for each fixed \(h\). In fact its \(O(n^{-1})\) translation changes
Gaussian stability to first order by twice the price of its marginal
mass change, which is \(-n^{-1}\mathfrak m(h)+o(n^{-1})\). Hence
\[
 \lim_{n\to\infty} n\bigl(V_2(r,m)-A_{2,n}(r,m)\bigr)
 =\tfrac12\min_{h\in\operatorname{Opt}_2(r,m)}
       \bigl(\|b_h\|_{a_h}^2-r^2\langle b_h,\mathcal K_h b_h\rangle\bigr).
                                                               \tag{AX5}
\]
The minimum is positive if the simplex falls strictly below the
global value, since it is attained and (CH2) makes the coefficient
strictly positive at every nonsimplex maximizer. If the simplex is
optimal the minimum is zero. Formula (AX5) controls partitions with
radial invariance in the spectator block. The unrestricted conical
class can use the spectator directions separately.
[Exact loss of the translated lift](THREE_CELL_CONICAL_HOMOGENIZATION.md#stability-loss-of-the-balanced-lift-of-a-maximizer)
