# Rigidity of isolated three-cell interface ends

Fix a balanced planar three-cell Gaussian maximizer at \(0<r<1\).
Use orthonormal coordinates about the Gaussian origin in which an
isolated end is a continuous graph \(y=h(x)\) for \(x>S\). In the
whole half-plane \(x>S\), one label fills the region below the
graph and a second fills the region above it, up to null sets.

- A bounded end has a finite transverse limit determined by its pair
  pressure, and \(h'(x)\to0\).
- A sublinear end with uniformly bounded oscillation on unit intervals
  is bounded.
- An end converging to a finite transverse limit from one side is a
  straight ray, as is an eventually monotone sublinear end or an
  eventually periodic end.
- A parallel group can be treated together when every graph approaches
  its finite offset from one side and the limiting strip slopes satisfy
  the conditions below.
- A maximizer with connected strict-winner cells and three pairwise
  obtuse ends of these straightening types is the centered simplex,
  under the network hypotheses specified below.

The boundedness and straightening arguments use the Gaussian winner
equations. The network consequence also uses regularity and complete
adjacency for actual maximizers.
[Winner and interface structure](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)
[Related affine limits for reflected graphs](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#local-affine-limits-at-infinity)
[Shifted-contour straightening of every bounded isolated graph](THREE_CELL_TAIL_FLATNESS.md)

## Conditional equation

Let \(f=1_{A_i}-1_{A_j}\) for the lower and upper labels, respectively.
The pair pressure \(\beta=\beta_i-\beta_j\) satisfies \(|\beta|<1\).
Indeed the pair wins with opposite strict inequalities in its two
positive-mass cells and \(T_rf\) lies strictly between \(-1\) and \(1\)
at every point. On the wall \(T_r f(x,h(x))=\beta\).

Write \(\sigma=\sqrt{1-r^2}\), let \(p_\sigma\) be the centered normal
density of variance \(\sigma^2\), and put

\[
 \psi(u)=2\Phi(u/\sigma)-1.
\]

For an arbitrary observation \((x,y)\), the conditional source
coordinates are independent normals with means \(rx,ry\) and variance
\(\sigma^2\). Integrating in the second coordinate gives

\[
 U(x,y)=T_rf(x,y)
   =\int_{s>S}\psi(h(s)-ry)p_\sigma(s-rx)\,ds+L(x,y),       \tag{IE1}
\]

where \(L\) is the contribution from \(s\le S\). At large \(x\) it is
\(O(e^{-c x^2})\), uniformly in \(y\). The same bound holds for its first
derivatives because differentiating the independent normal densities
inserts factors linear in their centered variables. The Gaussian tail
then absorbs the horizontal factor, and the absolute first moment of
the vertical factor is uniform.

## Bounded ends converge to the pressure offset

Suppose \(h\) is bounded on its tail and increase \(S\) so the bound
holds there. On the graph, (IE1) gives

\[
 \partial_y U(x,h(x))=-r A(x)+O(e^{-c x^2}),\qquad
 A(x)=\int_{s>S}\psi'(h(s)-rh(x))p_\sigma(s-rx)\,ds .
\]

Here \(A\) is bounded below for large \(x\) because \(\psi'>0\)
on the compact interval of arguments. The implicit function theorem
therefore makes \(h\) analytic on a sufficiently late tail. The
uniform bounds on the first two derivatives of a Gaussian convolution
of a bounded function, together with the lower bound on
\(|\partial_y U|\), give uniform bounds on \(h'\) and \(h''\) by
implicit differentiation.

Increase \(S\) to be within that smooth tail. Differentiate (IE1)
on the wall and integrate by parts in \(s\) in its horizontal
derivative. The endpoint term at \(S\) is exponentially small, and
the term at infinity vanishes. The result is

\[
 h'(x)=\frac{1}{A(x)}\int_{s>S}h'(s)\psi'(h(s)-rh(x))
                    p_\sigma(s-rx)\,ds+O(e^{-c x^2}).     \tag{IE2}
\]

Let \(M=\limsup_{x\to\infty}h'(x)\) and choose \(x_n\to\infty\)
with \(h'(x_n)\to M\). The derivative bounds and a diagonal
Arzelà–Ascoli argument give a subsequence for which
\(h(rx_n+u)\to H(u)\) in \(C^1_{mathrm{loc}}(\mathbb R)\) and
\(h(x_n)\to v\). The limit function is bounded and \(H'(u)\le M\)
for every \(u\). In (IE2), centered at \(rx_n\), Gaussian domination
gives \(M\) as the weighted mean of \(H'\) for the strictly positive
density proportional to \(\psi'(H(u)-rv)p_\sigma(u)\). This forces
\(H'\equiv M\), by continuity. Boundedness of \(H\) gives \(M=0\).
Using a sequence realizing the lower limit proves that it too is zero.

For each fixed \(u\), \(h(rx+u)-h(rx)\to0\) because \(h'\to0\).
By (IE1) and dominated convergence,

\[
  \psi(h(rx)-rh(x))\longrightarrow\beta .
\]

Set \(b=\psi^{-1}(\beta)/(1-r)\) and \(e=h-b\), so that
\(e(rx)-re(x)\to0\). Replacing \(x\) by \(x/r\) and taking the upper
limit in absolute value gives \(\limsup|e|\le r\limsup|e|\).
Boundedness then implies \(h(x)\to b\). In particular a periodic
tail, which is bounded by continuity, must be constant.

## Sublinear ends with bounded local oscillation or monotonicity

Here sublinearity means \(h(x)=o(x)\). For sufficiently large \(x\),
choose fixed \(R,C>0\) such that the conditional rectangle
\(|s-rx|\le R\), \(|t-rh(x)|\le C\) lies in \(s>S\) and has probability
greater than \((1+|\beta|)/2\). The rectangle cannot lie entirely on
one label side of the graph. Otherwise the expectation of \(f\) would
have absolute value greater than \(|\beta|\), contrary to stationarity.
Consequently \(rh(x)\) lies between the infimum of \(h\) on the
horizontal interval minus \(C\) and its supremum plus \(C\).

If the oscillation of \(h\) on each unit interval is uniformly bounded,
chaining across intervals on \([rx-R,rx+R]\) gives
\(|h(rx)-rh(x)|\le C_1\). Thus for all sufficiently large \(x\) and
every integer \(n\ge1\),

\[
 |h(x)|\le r^n|h(x/r^n)|+\frac{C_1}{1-r}.
\]

Sublinearity makes the first term tend to zero. This proves boundedness
and hence convergence to \(b\) with \(h'\to0\).

For an eventually monotone \(h\), reflect the vertical coordinate and
exchange the two labels if necessary to take \(h\) nondecreasing.
The rectangle bound then gives \(h(rx-R)\le rh(x)+C\). With
\(t=x+R/(1-r)\) and \(H(t)=h(t-R/(1-r))-C/(1-r)\), this reads
\(H(rt)\le r H(t)\) for all sufficiently large \(t\). If \(h\to+\infty\),
choose \(t_0\) sufficiently large that \(H(t_0)>0\). Iterating the
inequality at \(t_0 r^{-n}\) yields
\(H(t_0r^{-n})/(t_0r^{-n})\ge H(t_0)/t_0>0\) for all \(n\). This
contradicts sublinearity. A nondecreasing tail is bounded below, so
its limit is finite. It approaches that limit from one side, to
which the following argument applies.

## A one-sided approach forces exact straightness

Suppose \(h(x)\to b\) and \(e(x)=h(x)-b\) eventually has one sign.
By the bounded-end result \(\psi((1-r)b)=\beta\). Write \(g=e\) or
\(g=-e\), so \(g\ge0\) on a tail \(s>S\). The mean-value formula in
(IE1) at \(y=h(x)\) gives

\[
 \begin{split}
  r g(x)A_0(x)&=\int_{s>S}K(x,s)g(s)p_\sigma(s-rx)\,ds+B(x),\\
  K(x,s)&=\int_0^1\psi'((1-r)b+\theta(e(s)-re(x)))\,d\theta,\\
  A_0(x)&=\int_{s>S}K(x,s)p_\sigma(s-rx)\,ds ,
 \end{split}                                               \tag{IE3}
\]

where \(|B(x)|\le2\Pr[N(rx,\sigma^2)\le S]\). On \(x,s>S\) the
kernel \(K\) is bounded above and below by positive constants.

If \(g\) vanishes identically on a tail the result follows. Otherwise
there is a compact interval \(I=[a,a']\subset(S,\infty)\) of positive
length on which \(g\ge\delta>0\). Its contribution to the integral
in (IE3) is bounded below by a positive multiple of
\(p_\sigma(a-rx)\) for \(rx>a'\). This dominates the error \(B(x)\)
for all sufficiently large \(x\). The Gaussian upper bound for that
error uses \(\exp(-(rx-S)^2/(2\sigma^2))\), and the ratio grows
exponentially in \(x\) since \(a>S\). In particular \(g(x)>0\)
for all sufficiently large \(x\).

Fix \(0<\ell<r<u<1\). For large \(x\), \(I\) is disjoint from
\([\ell x,ux]\). Use its contribution to absorb the error and keep
only the remaining contribution on that moving interval. With
\(P_x(ds)=A_0(x)^{-1}K(x,s)p_\sigma(s-rx)1_{s>S}\,ds\), this gives

\[
  r g(x)\ge \int_{\ell x}^{u x}g(s)\,P_x(ds),
  \qquad \mu_x=P_x([\ell x,u x])\ge1-C_2e^{-c_2x^2}.        \tag{IE4}
\]

The mass estimate follows from the Gaussian tail and the uniform
kernel bounds. Choose \(R_0\) large enough that (IE4) holds for
\(x>R_0\), all its mass lower bounds there are positive, and \(g>0\)
on \([\ell R_0,R_0]\). Iterate (IE4) with the normalized probability
on \([\ell x,ux]\) until first reaching \(x_N\le R_0\). The terminal
point lies in \([\ell R_0,R_0]\), where \(g\) has a positive minimum
\(m\). The stopping depth is finite because each step decreases \(x\)
by a factor at most \(u\), and \(x_N\ge\ell^N x\) implies
\(N\ge\log(x/R_0)/\log(1/\ell)\). Along every path before stopping,
\(x_{N-1-j}>R_0u^{-j}\). The product of the retained masses is
therefore bounded below by the strictly positive infinite product
\(\prod_{j\ge0}(1-C_2\exp(-c_2R_0^2u^{-2j}))\).
Iteration bounds \(g(x)\) below by this product times \(m\) times
\(r^{-\log(x/R_0)/\log(1/\ell)}\), which diverges. This contradicts
\(g(x)\to0\), proving that \(h=b\) for all sufficiently large \(x\).

## Parallel ends approaching their offsets from one side

Suppose the entire jump graph in a half-plane \(x>S\) consists of
finitely many ordered continuous graphs \(y=h_j(x)\), with a constant
label in each successive strip. Write \(d_j\) for the difference of
the lower and upper label unit vectors at graph \(j\), and suppose
each \(h_j\) converges to a finite offset \(b_j\), with \(h_j-b_j\)
eventually of one sign. Assume

\[
 d_j\cdot d_l\le0\quad(j\ne l),\qquad
 D_j=\sum_l K^0_{jl}>0,\qquad
 K^0_{jl}=(d_j\cdot d_l)\,p_\sigma(b_l-rb_j).
\]

Then every graph is eventually constant. The sign and row-sum
conditions hold for the connected-network label patterns with finite
offsets checked in the algebraic-end argument, including two
alternating switches with distinct offsets.
[Limiting strip slopes](THREE_CELL_SEMIALGEBRAIC_NETWORK_RIGIDITY.md#no-sublinear-divergence-or-decaying-puiseux-term)

Put \(e_j=h_j-b_j\). At source horizontal coordinate \(s>S\) the
expected label vector as a function of the observation height \(y\)
is the top strip's label vector plus
\(\sum_l d_l\Phi((h_l(s)-ry)/\sigma)\). Contract it with \(d_j\)
and set \(y=h_j(x)\) to obtain the conditional pair contrast.
Its limit at the stationary wall uses \(h_l=b_l\) and \(y=b_j\),
so the limiting expression equals the pressure on that graph.
Subtract it inside the horizontal integral and use the mean-value
formula to obtain

\[
 \begin{split}
 r A_j(x)e_j(x)
   &=\int_{s>S}\sum_l K_{jl}(x,s)e_l(s)p_\sigma(s-rx)\,ds+B_j(x),\\
 K_{jl}(x,s)
   &=(d_j\cdot d_l)\int_0^1p_\sigma\bigl(b_l-rb_j+
                 \theta(e_l(s)-re_j(x))\bigr)\,d\theta,\\
 A_j(x)&=\int_{s>S}\sum_l K_{jl}(x,s)p_\sigma(s-rx)\,ds .
 \end{split}                                               \tag{IE5}
\]

The term \(B_j\) integrates the original pair contrast minus its
pressure over \(s\le S\), so \(|B_j(x)|\le2\Pr[N(rx,\sigma^2)\le S]\).
By taking \(S\) large, the kernels \(K_{jl}\) can be made uniformly
close to \(K_{jl}^0\) for all \(s,x>S\).

Let \(\eta_j\in\{-1,1\}\) satisfy \(g_j=\eta_je_j\ge0\) on the tail,
and put \(w_j=p_1(b_j)\), where \(p_1\) is the standard normal density.
Gaussian detailed balance gives \(w_jK_{jl}^0=w_lK_{lj}^0\).
The nonpositive off-diagonal entries then yield

\[
 \sum_j w_j\eta_j\eta_lK_{jl}^0\ge\sum_j w_jK_{jl}^0=w_lD_l .
\]

Choose \(\epsilon>0\) such that \(q=r(1+\epsilon)/(1-\epsilon)<1\).
Increasing \(S\) further gives \(A_j(x)\le(1+\epsilon)D_j\) and
\(\sum_j w_j\eta_j\eta_lK_{jl}(x,s)\ge(1-\epsilon)w_lD_l\), by
uniform approximation and positivity of all \(D_l\).
Sum (IE5) with weights \(w_j\eta_j\). The nonnegative function
\(H=\sum_j w_jD_jg_j\) satisfies

\[
 qH(x)\ge\int_{s>S}H(s)p_\sigma(s-rx)\,ds
                 -C\Pr[N(rx,\sigma^2)\le S].
\]

If \(H\) is nonzero arbitrarily far out, a fixed positive interval
of \(H\) absorbs the error exactly as in (IE3)–(IE4). Retain the
Gaussian mass on \([\ell x,ux]\) with \(0<\ell<r<u<1\). The same
first-passage iteration with multiplier \(q<1\) forces \(H\) to
grow, contrary to \(H(x)\to0\). Hence \(H\) vanishes on a tail.
Its coefficients are positive, so every \(e_j\) vanishes there.

## Three pairwise obtuse ends

Suppose every strict-winner cell is connected and unbounded, and
outside a compact set the jump graph consists of exactly three
proper ends. Assume end \(j\) is a transverse graph \(y=h_j(x)\)
with \(h_j(x)=o(x)\) in orthonormal coordinates whose longitudinal
direction is \(v_j\). Suppose \(v_i\cdot v_j<0\) for \(i\ne j\) and
each graph is eventually monotone, eventually periodic, or converges
to a finite limit from one side.

For fixed \(j\) the longitudinal projections of the other ends are
\(x(v_i\cdot v_j)+o(x)\to-\infty\), so far enough along the
\(v_j\) half-plane its graph is the only jump interface. Each side
has one almost-everywhere label, since it is connected and contains
no jump, and the two labels are distinct. Thus all three ends are
isolated in the sense used above and all three tails are straight.

Local analytic stratification and the finite number of straight
tails give a finite graph after one-point compactification. The
connected strict-winner cells give exactly three faces. There is
at least one genuine triple junction. Without one, a regular jump
component is a simple closed curve or a proper embedded line by
the active-interface regularity. It separates the plane, with its
incident connected labels on opposite sides. The third connected
label belongs to one side and cannot share a regular interface
with the label on the other, contrary to complete adjacency.

The vertex at infinity has degree three. Euler's formula with
three faces gives
\(\sum_z(\deg z-2)=4-2C\), where \(C\) is the number of graph
components. Actual interface vertices have no degree-one ends,
and every genuine triple has degree three or four. Thus there
is just one finite triple, of degree three, and \(C=1\).
A bounded loop would trap a bounded face, whereas all three
cells are connected and unbounded. All three arms therefore
join this triple to infinity. Each straight tail extends along
its regular analytic arm by the identity theorem applied to the
linear coordinate normal to its line. The resulting partition is
conical about the triple and is the centered regular simplex by
translated-cone rigidity.
[Arc and junction regularity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
[Complete adjacency](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#complete-adjacency)
[Graph counting](THREE_CELL_SEMIALGEBRAIC_NETWORK_RIGIDITY.md#counting-vertices-and-ends)
[Translated-cone rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

Shifted-contour straightening rules out curvature arbitrarily far
along any bounded isolated end. Deriving the required isolated
graph geometry for the ends of a general maximizer remains open.
[Exact straightness from boundedness](THREE_CELL_TAIL_FLATNESS.md#reaching-the-straightness-threshold)
