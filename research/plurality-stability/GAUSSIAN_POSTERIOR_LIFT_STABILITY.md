# Stability of the posterior simplex lift at a maximizer

Let \(A\) be an actual balanced three-cell Gaussian
noise-stability maximizer at \(0<r<1\). Work in its
active plane. For \(0<s<1\) let \(F_s=\mathcal L_s A\)
be the posterior simplex lift. At \((x,z)\), with \(z\)
standard Gaussian in the centered label plane \(H\), this
assigns the maximizing label of
\(z_i-\beta_i(p_s(x))\), where \(p_s=T_s1_A\) and
\(\beta\) is the centered model price map. Its conditional
label probabilities given \(x\) are \(p_s(x)\), so the
lift is balanced.
[The posterior lift and its first moments, WC6](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-posterior-entropy-and-centroid-tradeoff)

For every \(0<s<1\),
\[
 S_r(F_s)=S_r(A)\quad\Longleftrightarrow\quad
             A\text{ is a regular simplex}.                 \tag{PL1}
\]
In particular the lift strictly decreases full stability
at any hypothetical nonsimplex maximizer, even though it
strictly improves its first moments. The near-identity
expansion, in the units of (JT4), is
\[
 S_r(A)-S_r(F_s)
   = r(1-s^2)\bigl[S'_r(A)-\mathcal E_r(A)\bigr]
                      +o(1-s^2),\qquad s\uparrow1 .           \tag{PL2}
\]
The coefficient is strictly positive for a nonsimplex
maximizer. These statements use maximality of the fixed
partition \(A\). They do not establish whether a
nonsimplex maximizer exists.
[The aligned label-direction interaction, JT4](THREE_CELL_COUPLED_TRANSLATION_TEST.md)

## Equality at one lift scale

Suppose \(S_r(F_s)=S_r(A)\). Dimension independence
makes \(F_s\) a maximizer in the lifted space. Its
cell first-moment matrix has \(x\)-block \(s M_A\)
and \(z\)-block \(\overline L/\sqrt2\), with
\(\overline L=\int L_{p_s(x)}\,d\gamma_2(x)\). The
latter matrix is positive definite on \(H\), since
every posterior mass is positive and the model
facet Laplacian is positive definite there. Put
\(B=\sqrt2\,\overline L^{-1}sM_A\). Each direction
\((u,-Bu)\) is orthogonal to all first moments of
the lift. Dimension reduction makes the lifted labels
invariant, modulo null sets, under translation
by each such direction.
[Dimension reduction over the centroid span](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md)

Fix \(u\). By Fubini the two model assignments with
prices \(\beta(p_s(x+u))+Bu\) and \(\beta(p_s(x))\)
agree for almost every \(z\), for almost every \(x\).
Hence they have the same masses. The uniqueness of
centered model prices gives equality of those price
vectors. Their continuity in \(x\) gives
\[
                    \beta(p_s(x+u))=\beta(p_s(x))-Bu
                    \qquad\text{for all }x,u.
\]
Thus the lifted winning scores are affine, with
affinely independent sites because their \(z\)
components are the coordinate sites in \(H\).
Stationary affine-fan rigidity makes \(F_s\) a
centered regular simplex. Equality in the last
first-moment comparison of (WC6) consequently makes
\(A\) a centered regular simplex. Conversely, for a
regular simplex \(A\) the lift is an orthogonal
cylinder of the same simplex by (WC6), which proves
(PL1).
[Stationary full-rank affine fans](AFFINE_SIMPLEX_FAN_RIGIDITY.md)

## A flip estimate for the quadratic interaction

Use the conventions of (JT4). Orient each wall
\(e=(i,j)\) out of \(A_i\), put \(d_e=e_i-e_j\) and
\(q_e=d_e/\sqrt2\), and write
\(d\nu_e=\phi_2\,d\mathcal H^1|_{\Sigma_e}\).
Set \(\epsilon=\sqrt{1-s^2}\). Write \(f=1_A\)
also for its cylinder, and put \(e_s=1_{F_s}-f\).
The exact conditional marginal identity is
\[
 \int (e_s)_i(x,z)\,d\gamma_H(z)=p_{s,i}(x)-f_i(x),
 \qquad
 \int |(e_s)_i(x,z)|\,d\gamma_H(z)
                         =|p_{s,i}(x)-f_i(x)| .             \tag{PL3}
\]
The second identity holds since \(f_i\) is an indicator.
By self-adjointness and the Ornstein–Uhlenbeck
semigroup identity, the mixed energy of \(1_{F_s}\)
and \(f\) at correlation \(r\) is exactly \(S_{rs}(A)\).
Consequently, with \(T_r\) now acting in the lifted space,
\[
 S_r(A)-S_r(F_s)
       =2[S_r(A)-S_{rs}(A)]-
                       \sum_i\langle(e_s)_i,T_r(e_s)_i\rangle .
                                                               \tag{PL4}
\]

The normalized discrepancy measures
\(\epsilon^{-1}(e_s)_i\,d\gamma_2\,d\gamma_H\) have
bounded total variation with weight
\(e^{\eta(|x|^2+|z|^2)}\), for any fixed
\(\eta\) satisfying
\[
                      \frac{r}{2(1+r)}<\eta<\frac14 .
\]
To verify this, choose \(1/2<q<1-2\eta\). By Hölder
in \(z\), the weighted integral of \(|(e_s)_i|\)
over a fiber is at most
\(C_{\eta,q}|p_{s,i}(x)-f_i(x)|^q\), using (PL3).
It remains to bound its \(x\)-integral against the
Gaussian weight \(e^{-(1/2-\eta)|x|^2}\) by \(O(\epsilon)\).

For an indicator \(g\) of a cell, split its flip
probability into a dilation by \(s\) and a Gaussian
translation of variance \(\epsilon^2\). The inequality
\((u+v)^q\le u^q+v^q\) permits separate bounds. For
the dilation the indicator difference is unchanged
by the \(q\)-power. Integrating the affine flow from
\(s\) to \(1\) gives \(O(1-s)\) times
\(\int(1+|x|)e^{-c|x|^2}\,d|Dg|(x)\) for some \(c>0\).
This follows first for smooth approximations by the
fundamental theorem of calculus and change of
variables, and then for the indicator by convergence
in weighted \(L^1\).

For the Gaussian translation, write
\[
 b_h(x)=\dashint_{B_h(x)} |g(y)-g(x)|\,dy .
\]
For \(h\le1\), cover the plane by squares of side \(h\)
and enlarge each square so that it contains \(B_h(x)\)
for every \(x\) in the original square. If \(p\) is
the fraction of the enlarged square in the cell, the
integral of \(b_h^q\) on the original square is at
most \(C h^2[p(1-p)^q+(1-p)p^q]\). This is at most
\(C h^2\min(p,1-p)^{1/2}\), since \(q\ge1/2\).
The relative planar isoperimetric inequality bounds
it by \(C h\) times the perimeter in the enlarged
square. With the Gaussian weight, the same argument
uses a smaller Gaussian exponent on the perimeter,
since each enlargement has bounded diameter and
the squares have bounded overlap.

Decompose the translation Gaussian into balls of
radii \(2^j\epsilon\) centered at the observation
point. Its flip probability is bounded above by
\(\sum_{j\ge0} C4^j e^{-c_0 4^j}b_{2^j\epsilon}\)
for a numerical \(c_0>0\), by bounding the density
on each annulus. The \(q\)-power estimate and the
square bound sum to \(O(\epsilon)\) on scales
\(2^j\epsilon\le1\). Bound \(b_h\) by one at larger
scales, where the sum is exponentially small.
Changing variables from \(sx\) to \(x\) only weakens
the fixed Gaussian exponent by a bounded factor.
All the weighted perimeter integrals are finite
by polynomial length growth of the active network.
This proves the weighted total variation bound.
[Polynomial length growth of the active planar network](GAUSSIAN_OUTWARD_PERIMETER_GROWTH.md#cubic-length-bound-in-the-active-plane)

## The normal-flux limit and the leading loss

As vector measures on the lifted space, the normalized
discrepancies converge locally weakly to
\[
 \frac{e_s}{\epsilon}\,d\gamma_2\,d\gamma_H
       \ \rightharpoonup\
          \sum_e d_e(q_e\cdot z)\,d\nu_e(x)d\gamma_H(z).      \tag{PL5}
\]
First take a compact chart in the interior of a regular
pair wall and write \(x=\xi+\epsilon tN_e(\xi)\).
The conditional source is
\(\xi+\epsilon(tN_e(\xi)+Z)+O(\epsilon^2)\) for fixed
\(\xi,t,Z\). The original cell \(i\) is on the negative
normal side, so the two incident posterior probabilities
tend to \(\Phi(-t)\) and \(\Phi(t)\). The spectator
probability tends to zero.

For the model price map these limits make
\((\beta_i-\beta_j)/\sqrt2\to t\), and the spectator
price tends above both incident prices by an unbounded
amount. Indeed a positive lower bound on either
incident probability bounds its positive price
differences by the Gaussian pair tail. A spectator
price bounded above relative to the two incident
prices would give a positive probability of winning.
Once the spectator disappears the binary Gaussian
threshold identifies the limiting price difference.
Thus for almost every fixed \(t,z\) the model
indicator contrast relative to \(f\) tends to
\(d_e(1_{\{t<q_e\cdot z\}}-1_{\{t<0\}})\), whose
integral in \(t\) is \(d_e(q_e\cdot z)\).

For domination on the chart, shrink the normal tube
so a label change from \(x\) requires a source
displacement of at least its normal distance to the
wall. Since \(x-sx=O(\epsilon^2)\) on the compact
tube, (PL3) bounds the integrated flip probability
by a Gaussian tail at \(|t|-O(\epsilon)\). This
is integrable in \(t\) and justifies the normal
coordinate limit. On compact subsets off the walls
the flips are exponentially small. Near a singular
point \(p\), one has the localized bound
\[
 \limsup_{\epsilon\downarrow0}\frac1\epsilon
       \int_{B_\delta(p)} |p_{s,i}-f_i|\,dx
                 \le C |Df_i|(B_{3\delta}(p)).
\]
Indeed when \(\epsilon|Z|<\delta\), the straight
affine flow from \(x\) to \(sx+\epsilon Z\) stays
in \(B_{3\delta}(p)\) for sufficiently small
\(\epsilon\). Integrating the indicator difference
there costs at most
\(C[(1-s)(|p|+\delta)+\epsilon|Z|]\) times the
perimeter of that disk, first by smooth
approximation and then in \(BV\). The complementary
Gaussian displacement gives an exponentially small
error for fixed \(\delta\). Genuine junctions are
locally finite and have zero perimeter measure, so
letting \(\delta\downarrow0\) removes their
contribution. These chart and singular bounds give
(PL5).
[Regular edges and isolated junctions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

The lifted relative Mehler kernel satisfies
\(M_r(w,w')\le C_r e^{\lambda_r(|w|^2+|w'|^2)}\),
where \(\lambda_r=r/[2(1+r)]<\eta\).
The weighted total variation estimate therefore makes
the contribution of the tails uniformly small in the
normalized quadratic energy. On compact sets the
kernel is continuous, so (PL5) passes to the twofold
integrals by approximation with sums of separated
continuous functions. For correlated standard
\(Z,Z'\) in \(H\) one has
\(\mathbb E(q_e\cdot Z)(q_f\cdot Z')=r\,q_e\cdot q_f\).
The resulting limit of the normalized energy in
(PL4) is \(r\mathcal E_r(A)\), since
\[
 \mathcal E_r(A)
    =\sum_{e,f}(d_e\cdot d_f)(q_e\cdot q_f)
            \iint M_r(x,y)\,d\nu_e(x)d\nu_f(y).
\]
Also \(2[S_r(A)-S_{rs}(A)]/(1-s^2)\to rS'_r(A)\).
This proves (PL2).

By (JT4), the coefficient in (PL2) is at least
\((1-r)\mathcal R_r(A)\). If this remainder were
zero, its matrix \(L\) would satisfy \(LN_e=q_e\)
on every wall, surface-almost everywhere. Complete
adjacency makes \(L\) invertible because the label
normals span \(H\). Put \(u_i=-L^{-1}P e_i\), where
\(P\) projects onto \(H\). The distributional
Jacobian of the piecewise constant field
\(x\mapsto u_{A(x)}\) is symmetric positive
semidefinite, since its density on each wall is
\((u_j-u_i)\otimes N_e=\sqrt2\,N_e\otimes N_e\).
Smoothing gives a bounded gradient field of a convex
function. Normalizing those convex potentials at
the origin and taking a locally uniform limit
shows that the original field is also a convex
gradient. A convex function whose gradient takes
these finitely many values is the maximum of affine
functions with those slopes. Indeed each attained
slope supplies a supporting affine function, and
convexity forces its intercept to be the same at
any two differentiability points of that slope.
The partition is therefore a full-rank affine fan,
and stationarity forces it to be a regular simplex.
Thus the coefficient in (PL2) is strictly positive
away from the simplex orbit.

## Uniform entrance into the local simplex regime

Fix a compact correlation band \(J\subset(0,1)\). There is
\(s_J>0\) such that every ordinary balanced three-cell partition
\(A\), in every finite dimension, satisfies
\[
        S_r(\mathcal L_s A)\le S_r(T)
                  \qquad (r\in J,\ 0<s<s_J),                    \tag{PL6}
\]
with equality exactly for a simplex input.

In fact (WC11) bounds the lifted first-moment deficit by
\(s^2W_1(T)\), uniformly in \(A\). Dimension-independent
first-moment stability makes the lift's Gaussian label
distance to the simplex orbit at most \(Cs\). For small \(s\)
this lies in the radius of the uniform measurable local
comparison (DL1) throughout \(J\). Equality then makes the
lift a simplex, and equality in (WC6) makes \(A\) a simplex.
The uniform conclusion (PL6) relies only on local stability.
It applies even at correlations where the unrestricted
global inequality for \(A\) is open.
[First-moment stability](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#gaussian-margin-and-partition-distance)
[Uniform measurable local stability](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)
