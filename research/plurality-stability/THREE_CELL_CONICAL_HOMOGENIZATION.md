# Conical approximation of balanced Gaussian partitions

Fix a positive scalar correlation \(0<r<1\). The supremum of
balanced Gaussian three-cell stability over conical partitions of
all finite dimensions equals the unrestricted optimum. It is attained
by a conical partition in a finite dimension exactly when the
regular simplex attains the unrestricted optimum. A nonsimplex
maximizer, if it exists, also has a sequence of balanced conical
approximants. An exact asymptotic formula below relates the stability
loss of an explicit such sequence to the quartic agreement moment
in (AG5).

## Exact equivalence of the suprema

For an ordinary planar partition \(A\) with indicator tuple \(f\),
write \(R_n(z)=|z|/\sqrt n\) in an auxiliary Gaussian block
\(\mathbb R^n\). Outside its null set \(z=0\), define
\[
                 F_n(x,z)=f(x/R_n(z)).                       \tag{CH1}
\]
The tuple is conical about zero in \(\mathbb R^{2+n}\). Its
Gaussian label distance from the cylinder \(f(x)\) tends to zero.
Indeed \(R_n\to1\) in probability under the standard Gaussian
law. For each bounded measurable indicator, rescaling the argument
is continuous at one in Gaussian \(L^1\). This follows by
approximating it by continuous compactly supported functions and
using the common Gaussian density bounds on compact sets for
scale factors near one, with uniformly small tails. Boundedness
then allows averaging the label discrepancy over \(z\).

If \(f\) is balanced, the deviations of the masses of \(F_n\) tend
to zero. To repair them exactly, view a cone as a label on the
unit sphere. Standard Gaussian direction is uniform and nonatomic
there. Transfer angular mass directly from surplus labels to
deficit labels, with total changed mass half the sum of the
absolute deviations. This keeps the partition conical. A change
on Gaussian mass \(\delta\) changes stability by at most
\(2\delta\), by the union bound on the two arguments. The resulting
balanced conical partitions have stability tending to \(S_r(A)\).

The unrestricted optimum is attained in the plane by scalar
dimension reduction. The construction proves equality of the
suprema, since conical partitions are themselves admissible for
the unrestricted problem. If the simplex is globally optimal,
it attains both optima. Conversely, any conical partition attaining
the supremum in a finite dimension is an unrestricted maximizer,
so it is cylindrical over a two-dimensional active space by
dimension reduction and full active rank. Its planar projection
is still conical: for every fixed positive dilation factor, pass
the ambient homogeneity to the active indicator using Fubini.
The resulting invariance gives zero radial derivative in
distributions away from the origin. The planar conical
classification therefore identifies the maximizer as a regular
simplex.
[Active rank and scalar conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).

## Stability loss of the balanced lift of a maximizer

Suppose now \(A\) is an actual balanced planar maximizer. Use the
interface conventions of the maximizer analysis, and write \(b=x\cdot N\)
for the outward dilation speed, \(a\) for the normal score derivative,
and \(\mathcal K\) for the interface interaction operator. There
are translations \(v_n=O(n^{-1})\) in the active plane for which
\[
              F_n^*(x,z)=f(x/R_n(z)-v_n)
\]
is exactly balanced. Their stability satisfies
\[
 S_r(F_n^*)=S_r(A)-
       \frac{\|b\|_a^2-r^2\langle b,\mathcal K b\rangle}{2n}
       +o(n^{-1}).                                         \tag{CH2}
\]
By (G1) the numerator is at least \((1-r)\|b\|_a^2\), which
is strictly positive unless \(A\) is a simplex. If \(A\) ties the
simplex \(T\), use \(H,K,E_T,E_A\) as in (AG5) and put
\(d_A=S'_r(A)-S'_r(T)\). Then the coefficient in (CH2) is
equivalently
\[
 \frac{r(1-r^2)}4
 \left(\mathbb E[HK(1_{E_T}-1_{E_A})]-2d_A\right).           \tag{CH3}
\]
[First extra-coordinate form, translation and dilation identities](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds).
[Quartic agreement comparison, AG5](THREE_CELL_AGREEMENT_MOMENT_RIGIDITY.md#a-single-quartic-comparison-without-matched-first-responses).

For the proof, let \(Z,Z'\) be the auxiliary Gaussian pair at
correlation \(r\), and put \(U_n=\log R_n(Z)\), \(V_n=\log R_n(Z')\).
These two log norms satisfy
\[
 \begin{split}
  &\mathbb E U_n=O(n^{-1}),\quad \mathbb E|U_n|^3=O(n^{-3/2}),\\
  &\mathbb E U_n^2=\frac1{2n}+O(n^{-3/2}),\qquad
    \mathbb E[U_nV_n]=\frac{r^2}{2n}+O(n^{-3/2}).             \tag{CH4}
 \end{split}
\]
In fact \(X_n=(|Z|^2-n)/n\) has zero mean, second moment \(2/n\),
and fourth moment \(O(n^{-2})\). The covariance of \(X_n\)
with the corresponding variable from \(Z'\) is \(2r^2/n\),
by the fourth moment of one correlated scalar Gaussian pair.
On a fixed small log-norm neighborhood,
\(U_n=X_n/2+O(X_n^2)\). Here the cubic errors in second
moments are \(O(n^{-3/2})\) by Hölder's inequality. The
complements of such neighborhoods have exponentially small
probability and log-norm moments. Indeed the chi-square moment
generating function and the Chernoff bound give
\[
 \mathbb P(\pm U_n>u)\le
        \exp\{-\tfrac n2(e^{\pm2u}-1\mp 2u)\},\qquad u>0.
\]
For \(u\) bounded away from zero the exponent is bounded above
by \(-c n u\), which also controls the tail moments. These
bounds prove (CH4).

Gaussian differentiation of the fixed cells under scale and
translation shows that
\(m(u,v)=(\gamma(e^u(A_i+v)))_i\) is smooth near \((0,0)\).
The averaged mass map \(m_n(v)=\mathbb E m(U_n,v)\) has
\(m_n(0)=m_*+O(n^{-1})\) by (CH4) and the exponential tails.
Its derivative at zero tends to the isomorphism \(C_Av\)
whose components are \(-v\cdot\int_{A_i}x,d\gamma\). Its
second derivatives in \(v\) stay uniformly bounded in a small
neighborhood. To justify differentiation also on the log-norm
tails, the \(L^1\) bounds on first and second translation
derivatives of a scaled Gaussian density are bounded by
constants times \(e^u\) and \(e^{2u}\), respectively. The
fourth moment of \(R_n\) is bounded, so Cauchy–Schwarz controls
these exceptional integrals. The inverse function theorem with
these derivative bounds gives \(v_n=O(n^{-1})\) for which
\(m_n(v_n)=m_*\). This is exact balance for \(F_n^*\).

For each pair of log scales \(u,v\) and translation \(p\), let
\(J(u,v,p)\) be the bilateral agreement at cross-covariance \(rI\)
of the partitions \(e^u(A+p)\) and \(e^v(A+p)\), less the sum
of their two mass vectors paired with the stationary prices of
\(A\). The first derivatives of \(J\) at zero all vanish by
stationarity, even though individual Gaussian masses may vary.
The needed second derivatives are
\[
 J_{uu}(0)=J_{vv}(0)=-\|b\|_a^2,
 \qquad J_{uv}(0)=\langle b,\mathcal K b\rangle .             \tag{CH5}
\]
For the mixed derivative, take the first label displacement
on each endpoint. For the diagonal derivative, write
\(E=x\cdot\nabla\) as a distributional operator on indicators.
The second derivative on the first endpoint is
\(\sum_i\langle E^2 f_i,T_r f_i-\beta_i\rangle_{\gamma_2}\).
The adjoint of \(E\) in Gaussian measure is \(|x|^2-2-E\).
Its multiplication term cancels since the adjusted score
differences vanish on the interfaces. The remaining term is
\(-\sum_i\langle E f_i,E T_r f_i\rangle_{\gamma_2}=-\|b\|_a^2\)
by the normal gradient formula. Cutoffs are legitimate by
the polynomial interface estimates of (G3) and (G6).

The priced integral \(J\) is smooth through order three in a
neighborhood of zero, by changing variables in the Gaussian
density integrals for the fixed measurable partition. The
third derivatives have an integrable common Gaussian majorant
times a polynomial. It is bounded even outside this parameter
neighborhood because agreement probabilities and masses are
bounded. Taylor expand \(J(U_n,V_n,v_n)\) on the small log-norm
neighborhood and use the exponential tails on its complement.
Terms quadratic in \(v_n\) and terms mixing \(v_n\) with either
log scale are \(o(n^{-1})\). By (CH4)–(CH5), the leading
change is
\[
 \frac1{2n}\big(-\|b\|_a^2+r^2\langle b,\mathcal K b\rangle\big).
\]
The averaged prices cancel at exact balance, proving (CH2).
Finally, the dilation equation of (G6) gives
\(\|b\|_a^2-r^2\langle b,\mathcal K b\rangle
=-r(1-r^2)\langle g,b\rangle\). Equation (AG5) identifies
the latter pairing with the quartic contrast in (CH3).

## The remaining finite-dimensional question

Each finite-dimensional conical problem in dimension at least two
attains its optimum,
by weak compactness of the relaxed angular assignment and
strict positivity of the Gaussian operator, which rules out a
fractional optimizer by a balanced split. The conical optima
are nondecreasing with dimension by cylindrical extension.
The equivalence above reduces scalar global optimality to
whether this increasing sequence attains its limit at a finite
dimension. It does not supply that attainment. Moreover, a
nonsimplex global maximizer would give balanced conical
partitions with loss of order \(1/n\) by (CH2) and with a
fixed positive Gaussian distance from all tripod orbits.
For the distance assertion, \(F_n^*\) tends in label distance
to the cylinder over \(A\). Any ambient tripod can be rotated
by an orthogonal map fixing the active plane of \(A\) to a
tripod using at most two other directions. Its distance from
that cylinder can therefore be tested on a compact orbit in
four dimensions. The minimum is positive for a nonsimplex
\(A\). Thus a dimension-uniform conical rigidity estimate
would have to control these lifts as well as exact conical
maximizers.

## Equality of consecutive conical optima

Let \(s_d^{\rm cone}\) be the maximum balanced three-cell cone
stability in dimension \(d\ge2\), at the fixed correlation \(r\).
If \(s_d^{\rm cone}=s_{d+1}^{\rm cone}\), every optimizer in
dimension \(d\) is a cylindrical extension of the regular planar
simplex. Consequently, once the conical optima exceed the simplex
value, they increase strictly at every subsequent dimension.

Let \(A\) attain \(s_d^{\rm cone}\), with indicator tuple \(f\),
and extend it cylindrically by a scalar coordinate \(z\). By the
assumed equality it is a conical optimizer in \(\mathbb R^{d+1}\).
Let \(U_i=T_r f_i\) be the scores in the original space, and let
\(R\) have the chi distribution of a standard Gaussian norm in
dimension \(d+1\). A direction of the lifted sphere has first
\(d\) coordinates \(a\omega\), where \(0<a<1\) and
\(\omega\in S^{d-1}\). The first variation with conical competitors
gives the angular winner rule with scores
\[
          V_i(a,\omega)=\mathbb E U_i(aR\omega)
\]
and constant prices. Indeed differentiation along a feasible
segment in the relaxed angular assignment makes the linear gain
nonpositive at the optimizer. The mass supergradient of the
linear assignment gives the prices as in (AM3). Joint spherical
area disintegration gives this rule for almost every \(a,\omega\).
Continuity in \(a\) then makes it valid for all \(0<a<1\) at
almost every \(\omega\). Taking \(a\downarrow0\), conicality and
balance give \(U_i(0)=1/3\) for all three cells. Each label occurs
on a positive-area angular set, so the prices must agree. Choose
them all zero.
[Angular linear assignment uses the price argument before AM3](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#balanced-translations-with-a-matrix-correlation).

We next recover the ordinary score identity needed for translation
variation, without requiring the original scores to order the cells
pointwise. On the original sphere put \(u_i(\omega)=f_i(\omega)\),
using the conical representatives. For each \(a\), the angular
winner rule gives the distributional equality
\[
                  \sum_i V_i(a,\omega)D_S u_i=0,              \tag{CH6}
\]
where \(D_S\) is the angular gradient. To see it, differentiate
\(\max_i V_i=\sum_i u_i V_i\). At differentiability points of
the maximum the gradient of each active score equals its gradient.
The product rule gives (CH6). The scores are smooth on the sphere
because bounded Gaussian score derivatives can be differentiated
under radial averaging, using the finite moments of \(R\).

For each smooth tangent test field \(W\), define
\(F_W(s)=\sum_i\langle D_S u_i,U_i(s\cdot)W\rangle\).
This function is continuous with growth \(O(1+s)\), by the
distributional product rule and bounded first derivatives of
the Gaussian scores. Pairing (CH6) with \(W\) and changing radial
variables gives
\[
       \int_0^\infty F_W(s)s^d e^{-s^2/(2a^2)},ds=0
                                  \qquad (0<a<1).
\]
By the substitution \(t=s^2\) this is a Laplace transform of a
locally integrable function of polynomial growth, vanishing on a
real interval of positive parameters. Analyticity and uniqueness
of the Laplace transform imply \(F_W(s)=0\) for every \(s>0\).
Uniqueness here also follows by multiplying by a fixed exponential
and using Fourier uniqueness after analytic continuation. Thus
\(\sum_i U_i(s\omega) D_S u_i=0\) for every radius \(s>0\).
Since \(f_i(x)=u_i(x/|x|)\), the weak Euclidean gradient away from
zero is \(|x|^{-1} D_S u_i(x/|x|)\). This follows by polar
integration by parts with tests supported away from zero. It gives
\(\sum_i U_i D f_i=0\) away from the origin and hence everywhere:
inserting a cutoff within radius \(\varepsilon\) creates error
\(O(\varepsilon^{d-1})\) for bounded indicators and smooth tests.
We have proved the needed ordinary score identity in distributions.

Let \(v\) be orthogonal to every cell centroid. The sheared family
\(f^t(x,z)=f(x-t|z|v)\) is still conical. Its first volume changes
vanish, since a translated cell has first mass change proportional
to its centroid projected on \(v\). Hence the volume discrepancies
are \(O(t^2)\) by differentiation of Gaussian densities and the
finite moments of \(|z|\). For the stability variation define
\(M_f(v,v)=\sum_i\langle D_v f_i,T_r D_v f_i\rangle\).
Differentiating the joint density of the translated original
coordinates, conditionally on auxiliary correlated scalars \(Z,Z'\),
gives zero first derivative and coefficient of \(t^2\) equal to
\[
    \mathbb E[|Z||Z'|] M_f(v,v)
        +\mathbb E|Z|^2\sum_i\langle D_v^2 f_i,U_i\rangle.
\]
The differentiation is justified by integrable bounds on each
Gaussian density derivative in translation and the finite auxiliary
moments. For the last term, Gaussian integration by parts uses
\(\sum_i U_i D_v f_i=0\) to cancel the factor \(v\cdot x\).
Since \(D_v T_r=rT_r D_v\), the term is \(-r M_f(v,v)\).
Thus the coefficient equals
\[
                 (\mathbb E|ZZ'|-r) M_f(v,v).                \tag{CH7}
\]
The parenthesis is strictly positive: \(r=\mathbb E[ZZ']\) and
the nonsingular Gaussian pair has positive probability of opposite
signs. Also \(M_f(v,v)>0\) unless all indicators are invariant
under translations along \(v\), by their positive Hermite
multipliers after smoothing.
[Gaussian distributional translation calculation](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#balanced-translations-with-a-matrix-correlation).

The cone shear can be repaired exactly in mass at an arbitrarily
small quadratic cost. For any fixed \(\eta>0\),
the essential angular cell closures of the unsheared lift have
connected intersection graph, since finitely many disconnected
unions of these closed sets would separate the sphere. At any
intersecting closure the two adjusted angular scores coincide
by continuity and the winner rule. For each direction of each edge
of a spanning tree choose a source patch of positive angular measure
where the destination penalty is at most \(\eta\). The sheared
indicators converge in measure to the unsheared ones by Gaussian
translation continuity, so source patches retain positive mass.
The sheared scores averaged on angular rays converge uniformly
to the original angular scores: Gaussian smoothing converges
uniformly on compact balls by \(L^2\) kernel bounds, and the large
radial tail is uniformly negligible for bounded indicators.
The exact total volume discrepancy is \(O(t^2)\), so it can be
removed by transfers along the chosen tree. Take angular transfers
to keep conicality. Their linear energy change is at least
\(-C\eta t^2+o_\eta(t^2)\) with \(C\) independent of \(\eta\).
Their quadratic interaction is nonnegative. For any positive
coefficient in (CH7), first take \(\eta\) small and then \(|t|\)
small, contradicting conical optimality in dimension \(d+1\).

Therefore the original shape is cylindrical over its centroid
span, of dimension at most two. Its active dimension is exactly
two, since a conical partition on a line has only two half-rays
and all three labels have positive mass. In this plane the three
centroids span the plane, are distinct, and sum to zero. Expand
the angular winner inequality as \(a\downarrow0\). Since
\(\nabla U_i(0)=r\,m_i/\sqrt{1-r^2}\), where \(m_i\) is the
centroid, the leading nonconstant term of \(V_i\) is a positive
common multiple of \(a\langle m_i,\omega\rangle\). This gradient
formula follows by Gaussian differentiation and rescaling the
cones. Thus the labels are the winners of these distinct linear
functions, whose pair ties have angular measure zero. Their
centroids are noncollinear, so the cells are three convex sectors.
Equal Gaussian mass gives sector angles \(2\pi/3\), proving the
classification and excluding repeated disjoint sectors.

By (CH1), scalar simplex optimality is equivalent to the existence
of adjacent plateaus in arbitrarily large dimensions. Indeed each
plateau has the simplex value, and the optima are nondecreasing in
dimension. If scalar optimality fails, there is a finite first
dimension whose value exceeds the simplex value, and the sequence
increases strictly at every dimension from that point onward.
Equality at a single adjacent pair of dimensions does not control
possible increases after that pair.
