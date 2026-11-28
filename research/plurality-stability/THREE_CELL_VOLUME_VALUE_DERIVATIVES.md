# Volume derivatives and label prices at Gaussian maximizers

Fix \(0<\rho<1\) and a finite Gaussian dimension \(d\ge2\). Let
\(\mathcal V_\rho(m)\) be the largest three-cell stability at Gaussian
masses \(m=(m_1,m_2,m_3)\) in the interior of the probability simplex.
For a maximizer \(A\), center its label prices \(\beta(A)\) so that
\(\sum_i\beta_i(A)=0\), with \(A_i\) the strict winners of
\(T_\rho1_{A_i}-\beta_i(A)\). Then for every \(v\) with \(\sum_i v_i=0\),
\[
  \lim_{t\downarrow0}\frac{\mathcal V_\rho(m+tv)-\mathcal V_\rho(m)}t
       =2\max_{A\text{ maximizing at }m}\beta(A)\cdot v.             \tag{VV1}
\]
The limit also holds for directions tending to \(v\) as \(t\downarrow0\).
In particular, at \(m_*=(1/3,1/3,1/3)\), differentiability of the optimal
value is equivalent to every maximizer having equal prices.
If one maximizer at balance has a nonzero centered price, then
\[
  \lim_{t\downarrow0}\frac{\mathcal V_\rho(m_*+tv)-\mathcal V_\rho(m_*)}t
        \ \ge\ |\beta(A)|\,|v| .                                     \tag{VV2}
\]
Thus a tied nonsimplex maximizer with unequal prices would force
strictly positive outgoing slopes of the volume-value function in every
direction. The statements here do not establish differentiability or
locate an adjusted-score junction.

## Compactness and the price vector

Work in the weakly compact class of simplex-valued functions at masses
\(m\), which includes ordinary partitions. The Ornstein–Uhlenbeck
operator \(T_\rho\) is compact on Gaussian \(L^2\) in fixed dimension:
its degree-\(n\) Hermite eigenspaces are finite-dimensional and the
eigenvalues \(\rho^n\) tend to zero. Consequently stability is weakly
continuous. A maximizer is ordinary. Indeed, a nonordinary function
admits a nonzero bounded perturbation in opposite directions, supported
where two coordinates are uniformly positive, with zero mean in every
coordinate. Strict positivity of the quadratic form of \(T_\rho\) makes
the average of the two perturbed stabilities strictly larger.

For a maximizer \(A\), the linear functional with scores
\(U_i=T_\rho1_{A_i}\) is maximal at \(A\) among simplex-valued functions
of the same masses. Expand the stability at \(A\) toward such a
competitor and use the nonnegative quadratic term. As a function of
the prescribed masses, the optimal value of that linear functional is
concave by taking mixtures. Its supergradient at the interior point
\(m\) yields prices \(\beta_i\). The pointwise maximization of
\(U_i-\beta_i\) is then attained by \(A\) almost everywhere. The winners
are strict almost everywhere by analyticity and injectivity of
\(T_\rho\): a constant pair contrast would imply that two cell
indicators differ by a constant, contrary to the positive masses of all
three labels. The positive mass of each winner and \(0\le U_i\le1\)
give \(|\beta_i-\beta_j|\le1\).

## Ambient transport and uniform mass correction

For a compactly supported smooth vector field \(X\), let \(A_s\) be the
image of \(A\) under its flow. The volume velocity is
\[
  (\mathcal M_A X)_i
       =\int_{A_i}\operatorname{div}(\phi_d X)\,dx ,
\]
where \(\phi_d\) is the standard Gaussian density. A change of variables
in the joint Gaussian integral gives the stability velocity
\[
  \left.\frac{d}{ds}S_\rho(A_s)\right|_{s=0}
      =2\sum_i\int_{A_i}\operatorname{div}(\phi_d U_i X)\,dx
      =2\beta(A)\cdot\mathcal M_A X .                                \tag{VV3}
\]
For the second equality, write \(W=\max_i(U_i-\beta_i)\). It is locally
Lipschitz and agrees together with its weak gradient with the winning
score on each strict winner set. Subtracting the price term therefore
leaves \(2\int\operatorname{div}(\phi_d W X)\,dx=0\).

The flux map onto the label sum-zero space is surjective. If a vector
\(c\) annihilates its range, the locally integrable function
\(\sum_i c_i1_{A_i}\) has zero distributional gradient, by testing with
\(\phi_dX\). It is constant and each \(A_i\) has positive mass, so all
\(c_i\) agree. Equation (VV3) also proves uniqueness of the centered
prices, since the same equation holds for any winning-score prices.

Choose two compact smooth fields whose fluxes at \(A\) are independent.
The joint flow of these fields, in two small parameters, has an
invertible volume derivative at zero. For partitions converging to
\(A\) in \(L^1(\gamma_d)\), this derivative remains uniformly
invertible, since it depends on integrals against bounded test
functions. The mass and stability maps under the two flows have
uniformly bounded second derivatives in a fixed small parameter ball.
This follows by the change of variables formulas, because the smooth
flow is the identity off a compact set and derivatives of the joint
Gaussian density there are bounded in absolute value by an integrable
Gaussian times a polynomial. The inverse function theorem with
uniform derivative bounds supplies an exact prescribed small mass
correction of size \(h\) using flow parameters of size \(O(|h|)\).
Equation (VV3) and the two Taylor bounds give a corresponding
stability change of \(2\beta\cdot h+O(|h|^2)\) at a maximizer.
The estimates apply uniformly along a strongly convergent sequence
of maximizers using the two fields at its limit.

## Directional derivative at interior masses

Transport any maximizer at \(m\) to \(m+tv\), using the preceding
correction. This gives the lower bound in (VV1) for that maximizer's
price. For the upper bound, let \(t_n\downarrow0\), \(v_n\to v\),
and take maximizers \(A_n\) at \(m+t_nv_n\). A weakly convergent
subsequence has a limit at masses \(m\) and its energies converge.
The lower bound from a fixed maximizer at \(m\) implies that the
limit is also maximizing. It is ordinary, so the weak convergence is
strong: the sum of squared \(L^2\) norms of the indicator vector is one
both along the sequence and at the limit.

The bounded centered prices have further convergent subsequences.
Along each one the score functions converge locally uniformly to the
limit partition's scores. For example, each Mehler kernel row is
bounded in \(L^2(\gamma_d)\) uniformly for its first argument on a
compact set, so strong \(L^2\) convergence suffices. Passing to
almost-everywhere convergence of the indicators shows that the limiting
prices are winning-score prices for the limit partition. By uniqueness,
these are its centered prices. The same argument at fixed masses shows
that the set of maximizing prices there is compact, so the maximum in
(VV1) is attained.

Use the two fields of the limit partition to correct each \(A_n\)
exactly to masses \(m\). The corrected partition has stability at most
\(\mathcal V_\rho(m)\) and stability change
\(-2t_n\beta(A_n)\cdot v_n+O(t_n^2)\). On a subsequence achieving the upper
limit, this proves the required upper bound in (VV1).

At balance, cyclic relabeling rotates the centered prices through
angles \(0,2\pi/3,4\pi/3\) in their two-dimensional label plane.
One of their dot products with \(v\) is at least
\(|\beta(A)|\,|v|/2\). Equation (VV1) gives (VV2). A function
differentiable at balance and invariant under all label permutations
has zero derivative on that plane, so its maximizing prices must all
vanish. Conversely, if all prices vanish, the version of (VV1) with
varying directions and compactness of the unit circle give
\(\mathcal V_\rho(m_*+h)-\mathcal V_\rho(m_*)=o(|h|)\).

## Strong convexity across positive masses

For \(d\ge2\), \(0<r<1\), and three positive masses,
\(\mathcal V_r(m)-(1-r)|m|^2\) is convex on the mass simplex.
In particular, for every maximizer \(A\) at \(m\), its centered
price vector gives a global bound at every other interior mass
vector \(m'\):
\[
 \mathcal V_r(m')\ge \mathcal V_r(m)+2\beta(A)\cdot(m'-m)
                               +(1-r)|m'-m|^2.                \tag{VV8}
\]

To prove this, write \(f_i=1_{A_i}\) in the two-dimensional
active plane and \(z_i=\int x f_i(x)\,d\gamma_2(x)\).
The [full active rank theorem](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)
holds for all three positive masses. Hence the translation
volume derivative \(C_Av=(-z_i\cdot v)_i\) is an
isomorphism onto the mass plane \(\sum_i h_i=0\).
For a translation by \(v\), the exact volume and stability
are smooth in \(v\), by translating the Gaussian densities
in their integral formulas. For prescribed small \(h\),
the inverse function theorem gives a translation
\(v=C_A^{-1}h+O(|h|^2)\) with masses exactly \(m+h\).

For the priced functional \(S_r-2\beta(A)\cdot m\), its first
translation derivative is zero by the winner identity.
Its second translation derivative divided by two is
\((1-r)M_A(v,v)\) by the bulk translation identity (AM3).
That identity is an equality for every translation direction.
The restriction to mass-tangent directions is only needed
when using maximality to sign the derivative.
The constant Hermite mode of \(D_v f_i\) is \(z_i\cdot v\)
by Gaussian integration by parts, so the positive Hermite
multipliers give \(M_A(v,v)\ge |C_Av|^2\).
The translated partition therefore gives
\[
 \mathcal V_r(m+h)
   \ge \mathcal V_r(m)+2\beta(A)\cdot h
       +(1-r)M_A(C_A^{-1}h,C_A^{-1}h)+o(|h|^2).
                                                               \tag{VV9}
\]
[Bulk translation Hessian, AM3](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#balanced-translations-with-a-matrix-correlation).

Put \(U(m)=\mathcal V_r(m)-(1-r)|m|^2\). This function is
continuous by the volume correction and compactness argument
for (VV1). Applying (VV9) with both signs of a displacement
shows on each line in the mass simplex that the lower limit
of the symmetric second difference quotient of \(U\) is
nonnegative at every point. To see that this implies
convexity, suppose the restriction of \(U\) exceeds an
endpoint chord somewhere on a closed subsegment. Subtract the
chord and add a sufficiently small positive multiple of the
squared line coordinate. The resulting function still has
an interior maximum strictly above its endpoint values.
At that maximum the symmetric second difference quotient of
\(U\) has upper limit bounded by minus twice that positive
multiple, a contradiction.

By (VV9), \(2\beta(A)-2(1-r)m\) is a lower supporting slope
for \(U\) to first order at \(m\). Convexity makes this a global
supporting slope, by the monotonicity of the chord slopes
from \(m\) along any direction. This proves (VV8). Applying
it in both directions at \(m,m'\), for arbitrary maximizers
\(A,B\) there, gives
\[
  \bigl(\beta(B)-\beta(A)\bigr)\cdot(m'-m)
                       \ge (1-r)|m'-m|^2.                    \tag{VV10}
\]
Taking \(m'\) and \(B\) to be obtained by transposing two
labels \(i,j\) shows
\((\beta_i(A)-\beta_j(A))(m_i-m_j)\ge(1-r)(m_i-m_j)^2\).
This controls the ordering of prices for unequal masses.
At balance, (VV8) is consistent with multiple maximizing
price vectors and therefore leaves the differentiability
question in (VV1) open.

## Quadratic volume sensitivity wherever the simplex is unique

Suppose at a fixed positive correlation the simplex is the unique
balanced maximizer up to its orthogonal orbit and relabeling. For
\(m=m_*+h\) sufficiently close to balance, every maximizer \(A\) at
masses \(m\) satisfies
\[
 \big|\mathcal V_\rho(m)-\mathcal V_\rho(m_*)\big|\le C|h|^2,
 \qquad
 \operatorname{dist}_{L^1}(A,\operatorname{orb} T)+|\beta(A)|
          \le C|h| .                                                \tag{VV4}
\]
This uses the [measurable local growth bound (B11)](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
from the three-cell interface gap. At any fixed finite dimension the
constants can be chosen uniformly on compact sets of positive
correlations throughout which the balanced simplex is unique.
At the interior masses in question, every maximizer is planar
after rotation and passage to a cylinder. Its stability and
centered prices are unchanged under extension to that cylinder,
and the full-space distance is at most the distance on the plane.
Thus planar constants suffice in every dimension \(d\ge2\).
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)

The Gaussian \(L^1\) distance of partitions here is the sum of
the distances of the three indicator functions. Take the two mass
correction fields from \(T\), rotated when using another simplex in
its orbit. For any \(A\) close to that simplex, the associated two
stability derivatives differ from the simplex derivatives by at most
\(C\|A-T\|_1\). In the first expression in (VV3), the change in
the indicator is bounded in \(L^1\) and the change in scores is
bounded in \(C^1\) on the supports of the fields. The latter bound
follows from uniform bounds on the Mehler kernel and its first
derivatives in the first variable over all second variables when
the first stays in a compact set. The simplex derivatives vanish because its label
prices coincide. The two flux columns stay invertible near \(T\), so
\[
                       |\beta(A)|\le C\|A-T\|_1 .                   \tag{VV5}
\]

As \(h\to0\), compactness and uniqueness imply that all the
maximizers in (VV4) are close to the simplex orbit. Apply an ambient
flow of size \(O(|h|)\) to correct \(A\) to a balanced \(B\). A
distance estimate independent of \(A\)'s perimeter is
\[
  \operatorname{dist}_{L^1}(A,\operatorname{orb} T)
      \le C\operatorname{dist}_{L^1}(B,\operatorname{orb} T)+C|h| .
\]
Indeed the inverse flow distorts \(L^1\) distances by a bounded
factor, and it changes any simplex indicator only in an
\(O(|h|)\)-neighborhood of its facets on the compact support.
With \(s=\operatorname{dist}_{L^1}(B,\operatorname{orb} T)\),
correcting the energy by (VV3)–(VV5) and applying (B11) gives
\[
  S_\rho(A)\le S_\rho(T)-cs^2+C|h|s+C|h|^2
              \le S_\rho(T)-\tfrac{c}{2}s^2+C'|h|^2 .
\]
Flowing the simplex to masses \(m\) gives
\(S_\rho(A)\ge S_\rho(T)-C''|h|^2\) because its first derivative
vanishes. Therefore \(s=O(|h|)\), which proves (VV4).

For the uniform statement, use compactness of the correlation set
and operator-norm continuity of \(T_\rho\) there. A sequence of
varying-mass maximizers with \(h\to0\) still has strongly convergent
subsequences with balanced maximizing limits. The local growth
constants are uniform by (B11) and a finite subcover, and the
transport bounds are uniform on the compact correlation set.

## Joint volume and correlation slopes at a tie

Let \(0<\rho<1\), \(m\) be interior and \(v\) sum to zero. With
\(\dot S_\rho\) denoting the correlation derivative at a fixed
partition, the one-sided joint derivative is
\[
 \lim_{t\downarrow0}
 \frac{\mathcal V_{\rho+t\sigma}(m+tv)-\mathcal V_\rho(m)}t
   =\max_{A\text{ maximizing at }(\rho,m)}
        \bigl[\sigma\dot S_\rho(A)+2\beta(A)\cdot v\bigr].          \tag{VV6}
\]
On compact correlation intervals the OU operator has bounded first
and second operator derivatives, from the Hermite multipliers
\(n\rho^{n-1}\) and \(n(n-1)\rho^{n-2}\). Taylor expansion in
correlation is therefore uniform over partitions. For the lower
bound, transport a maximizing \(A\) at the base point to masses
\(m+tv\), then change correlation. The transported indicator
converges strongly to \(A\), so the derivative term tends to
\(\sigma\dot S_\rho(A)\).

For the upper bound choose maximizers at the varying pairs of
parameters and pass to a strong limit as in (VV1), using
operator-norm continuity in correlation. Their price vectors
converge along subsequences to the limit maximizer's prices, since
the Mehler rows converge locally uniformly in Gaussian \(L^2\)
as correlation varies. Perform the exact mass correction at the
varying correlation, then move the corrected indicators back to
correlation \(\rho\). The correction converges strongly to the
identity on the sequence of indicators because the flows tend
uniformly to the identity, their Gaussian density distortions are
bounded, and the original indicators converge strongly. The two
Taylor expansions give the upper bound in (VV6).

For example, suppose \(A\) ties \(T\) globally at a first transition
\(\rho_*\) at balance. Put
\(d_A=\dot S_{\rho_*}(A)-\dot S_{\rho_*}(T)\ge0\), where the sign
follows from simplex optimality immediately to the left. Cyclic
relabeling \(A\) preserves its correlation derivative. Subtracting
the simplex's correlation derivative in (VV6) with \(\sigma=-1\)
yields
\[
 \lim_{t\downarrow0}
  \frac{\mathcal V_{\rho_*-t}(m_*+tv)-S_{\rho_*-t}(T)}t
       \ge |\beta(A)|\,|v|-d_A .                                 \tag{VV7}
\]
Consequently a tie with unequal prices gives a family of quota
directions with a linear gain over the balanced simplex immediately
below the transition. For a tangential tie with unequal prices,
the gain holds in every nonzero quota direction.

## Joint covariance and volume response

Let \(\mathcal V(R,m)\) be the optimal Gaussian stability with
positive prescribed masses \(m\) and symmetric positive
cross-covariance \(0<R<I\). This section allows any fixed numbers
of cells and ambient coordinates. For symmetric \(B\) and
\(\sum_i h_i=0\), the joint right directional derivative is
\[
 \lim_{t\downarrow0}
       \frac{\mathcal V(R+tB,m+th)-\mathcal V(R,m)}{t}
 =\max_{A\text{ maximizing at }(R,m)}
       [\operatorname{tr}(B M_A)+2\beta(A)\cdot h].          \tag{VV11}
\]
Here \(M_A\) is the matrix of price-free covariance derivatives
from (AM1)–(AM2), and \(\beta(A)\) is the centered vector of
winning-score prices. The limit is unchanged if the directions
\(B,h\) vary and converge to the displayed directions.
[Covariance differentiation, AM1–AM2](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#covariance-derivatives-and-rotation-stationarity).

The compactness and price proof for (VV1) applies to \(T_R\):
its positive Hermite multipliers give an ordinary maximizing
limit, and the locally uniform Gaussian kernel convergence
makes the centered prices continuous along strongly convergent
maximizers. The mass-flux proof used for (VV3) gives a choice
of smooth compact transports realizing all small interior
mass corrections, since its annihilator consists of constant
label vectors at any positive masses. Smooth volume
correction and (AM2) give the lower bound in (VV11) for a fixed
base maximizer, first transporting it to \(m+th\) and then
changing \(R\). The first-order mass term is \(2\beta(A)\cdot h\).

For the upper bound, take varying-parameter maximizers and
extract a strong limit after using the lower bound and kernel
continuity to identify it as a base maximizer. Correct their
volumes back to \(m\) using fixed transport fields at the limit,
then move the resulting partitions to \(R\). The mass correction
has first-order effect \(-2\beta(A)\cdot h\) at the limit by
(VV3), and the covariance change has first-order effect
\(-\operatorname{tr}(B M_A)\). Taylor remainders are uniform for the
fixed compact transport fields and a compact neighborhood of
\(R\) in the positive covariance domain. The rules for varying
directions and attainment of the maximum follow from the same
subsequence argument. This proves (VV11), of which (VV6) is
the scalar case.

### Which directions separate a tied shape from the simplex branch

Return to the balanced three-cell problem in the plane at
\(R=rI\), and suppose the simplex \(T\) is globally maximizing.
Fix a small closed neighborhood of its orthogonal orbit in
Gaussian label distance, contained in the strict local
comparison neighborhood. Denote by \(\mathcal V_T\) the
supremum over ordinary partitions in that neighborhood. Then
\[
 \lim_{t\downarrow0}
 \frac{(\mathcal V-\mathcal V_T)(rI+tB,m_*+th)}{t}
  =\max_{A\text{ maximizing at }(rI,m_*)}
        [\operatorname{tr}(B(M_A-M_T))+2\beta(A)\cdot h].
                                                               \tag{VV12}
\]
The same formula applies to converging directions.
[Uniform measurable simplex local comparison, B11](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison).

For the lower bound on \(\mathcal V_T\), transport \(T\)
to the perturbed masses. Its priced first derivative vanishes
since its centered prices are zero. For the upper bound take
ordinary partitions within \(o(t)\) of the local supremum. The
same compactness argument as for (VV11) shows their strong limits
are balanced scalar maximizers in the closed neighborhood,
hence are simplex partitions. The transport fields of a limit
correct the perturbed volumes exactly. The first-order effect
of this correction tends to zero by continuity of the transport
velocity in (VV3), evaluated at the limiting simplex. This use
of the transport derivative does not require stationary prices at
the approximating partitions. The first covariance derivative
is that of the limiting simplex, with \(M_T=\dot S_r(T)I/2\)
by the simplex rotation symmetry. The corrected value is
bounded by \(S_r(T)\). These estimates give the local derivative
\(\operatorname{tr}(B M_T)\), proving (VV12).

If any tied maximizer has \(\beta(A)\ne0\) or \(M_A\ne M_T\),
some direction has strictly positive right side in (VV12).
Every sequence attaining the global supremum to \(o(t)\) along
such a fixed direction eventually leaves the chosen simplex
neighborhood. The joint value is differentiable at the base
point exactly when every tied maximizer has zero centered
prices and \(M_A=M_T\). Indeed its directional derivative is a
support function by (VV11). A support function is linear
exactly when its slopes agree, and convergence of directions
in (VV11) upgrades this linear derivative to differentiability
by compactness of the finite-dimensional unit sphere. At a
first scalar transition, matching the full covariance derivative
matrix includes matching the scalar slope
\(\dot S_r(A)=\dot S_r(T)\).

For a first scalar transition, let
\(d_A=\dot S_r(A)-\dot S_r(T)\), let \(g_A\) be the difference
between the eigenvalues of \(M_A\), and take \(B=-I+B_0\)
with \(B_0\) traceless. The right side of (VV12) is at least
\[
          -d_A+g_A\|B_0\|_{\rm op}+|\beta(A)|\,|h| .          \tag{VV14}
\]
To see this, rotate \(A\) to align the ordered eigenvectors of
\(M_A\) and \(B_0\). Independently, cyclically relabel \(A\) so that
the price vector has angle at most \(\pi/3\) with \(h\).
The centered scalar covariance does not change either optimizer
after these actions. Thus a positive right side gives
directions with distant maximizers even when the mean of
the two correlations decreases.

### Gaussian midpoint and difference moments at a tie

Take a maximizer \(A\) at scalar correlation \(r\) and define
the independent standard Gaussian vectors
\[
    Z=\frac{X+Y}{\sqrt{2(1+r)}},\qquad
    W=\frac{X-Y}{\sqrt{2(1-r)}} .
\]
Let \(E\) be the event that \(X,Y\) receive the same label,
\(z_i=\int x1_{A_i}\,d\gamma\), and
\(P=\sum_i\beta_i(A)\int (xx^{\mathsf T}-I)1_{A_i}\,d\gamma\).
The following unnormalized conditional moment identities hold:
\[
\begin{split}
  \mathbb E[Z1_E]&=\sqrt{2(1+r)}\sum_i\beta_i(A)z_i,\\
  \mathbb E[(ZZ^{\mathsf T}-I)1_E]
      &=(1-r^2)M_A+(1+r)P,\\
  \mathbb E[(WW^{\mathsf T}-I)1_E]
      &=-(1-r^2)M_A+(1-r)P.                               \tag{VV13}
\end{split}
\]
Also \(\mathbb E[W1_E]=0\) and \(\mathbb E[ZW^{\mathsf T}1_E]=0\)
by interchanging \(X,Y\).

For the first identity translate both labels by \(v\), or
equivalently replace the density of \(Z\) in the integral for
stability by its translate with mean \(-\sqrt{2/(1+r)}v\).
The first variation of stability is \(2\beta(A)\) dotted with
the mass change, whose derivative is \((-z_i\cdot v)_i\).
For the second identity the translation Hessian divided
by two is \(\mathbb E[(ZZ^{\mathsf T}-I)1_E]/(1+r)\).
The price term in the Hessian divided by two is \(P\).
Their difference is \((1-r)M_A\) by (AM3) in arbitrary
translation directions. Finally, a covariance perturbation
\(B\) changes the covariances of \((X+Y)/\sqrt2\) and
\((X-Y)/\sqrt2\) by \(B\) and \(-B\), respectively. Density
differentiation and (AM2) give
\[
 M_A=\frac{\mathbb E[(ZZ^{\mathsf T}-I)1_E]}{2(1+r)}
          -\frac{\mathbb E[(WW^{\mathsf T}-I)1_E]}{2(1-r)}.
\]
Substitution gives the last line of (VV13).

Thus a tied shape with no first-order joint advantage in any
direction has the same conditional first and second moments
of \((Z,W)\) given agreement as the simplex. Equations
(VV11)–(VV13) are necessary first-order tests. They do not
classify tied shapes with equal centered prices and covariance
derivative matrices.
