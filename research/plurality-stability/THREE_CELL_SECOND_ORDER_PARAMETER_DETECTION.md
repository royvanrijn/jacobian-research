# Second-order parameter tests at a balanced simplex tie

Consider a balanced planar three-cell Gaussian maximizer \(A\) at
\(0<r<1\). The first-order responses to unequal volumes and
anisotropic covariance can agree with those of a tied simplex.
The trace of the covariance Hessian in the two trace-free matrix
directions gives a further test. It differs from the scalar
correlation curvature by an exact expression in the dilation
speed. A different test applies when the volume-constrained
interface form has a zero mode that couples to a volume change
or to covariance.

## A covariance and dilation identity

Use the oriented interfaces, their Gaussian surface measure
\(\nu\), the stationary weight \(a\), and the interaction operator
\(K\) of the [maximizer variation formulas](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
Write \(f_i=1_{A_i}\). For a planar symmetric matrix \(B\), let
\[
 \Delta_B=\operatorname{tr}(B\nabla^2),\qquad
 L_B=\Delta_B-(Bx)\cdot\nabla,\qquad
 b_B(x)=(Bx)\cdot N(x).
\]
Here \(N\) is the outward normal of the first cell in the
oriented pair. Denote by \(S_R(A)\) the stability with
cross-covariance \(R\), at fixed standard Gaussian marginals,
and let \(M_A\) be its covariance gradient at \(R=rI\), in the
trace pairing.

For symmetric \(B,C\), there is the exact identity
\[
\begin{split}
 \sum_i\langle T_{\sqrt r}\Delta_B f_i,T_{\sqrt r}\Delta_C f_i\rangle
   ={}&r^2 d_R^2 S_{rI}(A)[B,C]
           +r\,\operatorname{tr}\!\left(\frac{BC+CB}{2} M_A\right)\\
     &+\frac{(1+r^2)\langle b_B,Kb_C\rangle
                       -2\langle b_B,b_C\rangle_a}{1-r^2},
\end{split}                                                   \tag{PD1}
\]
where \(\langle b_B,b_C\rangle_a=\int a b_B b_C\,d\nu\).
The derivatives on the left are Euclidean distributions before
Gaussian smoothing.

Indeed \(L_B\) is minus the sum of insertions of \(B\) on the
tensor factors of each Hermite chaos. The \(T_r\) cross pairing
of \(L_B f_i,L_C f_i\) separates insertions on distinct factors
and insertions on the same factor. The distinct factors give
\(r^2 d_R^2 S_{rI}[B,C]\). The same-factor term gives \(r\)
times the covariance derivative in direction \((BC+CB)/2\).
The chaos multiplier bounds justify summation.
[Covariance derivatives and Hermite lowering](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#covariance-derivatives-and-rotation-stationarity)

In distributions \(\Delta_B f_i=L_B f_i+(Bx)\cdot\nabla f_i\).
The last term is the negative outward speed \(b_B\), with the
signed incidence vector, on the interface measure before
Gaussian weighting. Its \(T_r\) pairing with the corresponding
term for \(C\) is \(\langle b_B,Kb_C\rangle\). Also
\(T_r L_C f_i=-r\,d_R(T_R f_i)|_{rI}[C]\), so the cross
pairing is \(r\langle b_B,g_C\rangle\), where \(g_C\) is the
derivative of the pair score on the interface.
The stationary affine dilation identity is
\[
                   K b_C=r^{-2}a b_C+(1-r^2)r^{-1}g_C .
\]
Substitute it for both cross pairings to obtain (PD1).
The polynomial-weighted cutoff bounds of (G6) give finite
weighted norms and absolutely convergent kernel interactions
for affine speeds. The differentiated scores have polynomial
growth at fixed \(r\), so the same bounds justify their interface
pairings by exhaustion. In Gaussian Hilbert space these
distribution pairings are the pairings after smoothing;
the Hermite coefficients of each fixed-order smoothed
derivative are square summable.
[Dilation identity and integrability, G6](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Anisotropic dilation identity, CT6](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions)

## Trace-free covariance curvature and quota curvature

Define
\[
\begin{split}
  \mathcal H_r(A)&=\sum_i\|T_{\sqrt r}\Delta f_i\|_2^2,\\
  b&=x\cdot N,\qquad B_a=\|b\|_a^2,\qquad
  E_b=\langle b,Kb\rangle,\qquad d=C_A^{\rm vol}b.
\end{split}
\]
Here \(C_A^{\rm vol}\) denotes the mass flux of a normal speed.
The quantity \(\mathcal H_r\) is the trace of the covariance
Hessian in the coordinates
\(R(x,y)=\left(\begin{smallmatrix}r+x&y\\y&r-x\end{smallmatrix}\right)\).
Taking \(B=C=I\) in (PD1) gives
\[
  \mathcal H_r(A)=r^2S_r''(A)+rS_r'(A)
                       +\frac{(1+r^2)E_b-2B_a}{1-r^2}.
                                                               \tag{PD2}
\]
[Trace-free covariance Laplacian, AM7](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#comparison-of-two-arbitrary-spreads-at-fixed-trace)

Suppose \(A\) is nonsimplex and ties \(T\) at \(r\), distinct
simplex-optimal scalar correlations accumulate at \(r\), and the
fixed-partition deficit \(D_\rho=S_\rho(T)-S_\rho(A)\) satisfies
\(D'(r)=0\). Then \(D''(r)\ge0\). At \(T\) the speed \(b\)
vanishes, so
\[
 \mathcal H_r(A)-\mathcal H_r(T)
   =-r^2 D''(r)- B_a
                         -\frac{1+r^2}{1-r^2}Q_0(b),    \tag{PD3}
\]
where \(Q_0(b)=B_a-E_b\).
If \(d=0\), balance gives \(Q_0(b)\ge0\).
Also \(B_a>0\) for a nonsimplex maximizer because \(b=0\)
would make the partition conical. In this case \(A\) has
strictly lower trace-free covariance curvature than \(T\).
This compares the two fixed shapes and their
rotation-averaged anisotropic perturbations to second order.
[Centered and translated cone rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

For a general mass flux, write \(P_{\rm tr}\) for weighted
projection onto translations and put
\(B_\parallel=\|P_{\rm tr}b\|_a^2\) and
\(B_\perp=\|b-P_{\rm tr}b\|_a^2\).
The interaction is \(1/r\) times the weighted squared norm
on translations and at most the weighted squared norm
on their orthogonal complement, with zero mixed pairing.
Thus
\[
 \mathcal H_r(A)-\mathcal H_r(T)
   \le -r^2D''(r)-B_\perp+
                      \frac{1-r}{r(1+r)}B_\parallel .     \tag{PD4}
\]
For a nonsimplex balanced maximizer, \(B_\perp>0\);
equality would give a cone about some point.
[Translation eigenspace and negative index, G17–G19](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

Assume now that the centered prices of \(A\) vanish and
that \(r\) is a first scalar transition with a unique balanced
simplex immediately to its left. If
\(\mathcal H_r(A)\ge\mathcal H_r(T)\), (PD3) gives
\[
                -Q_0(b)\ge
   \frac{1-r^2}{1+r^2}\big(B_a+r^2D''(r)\big).            \tag{PD5}
\]
If the right side exceeds \(\Lambda_r|d|^2\), where
\(\Lambda_r\) is the simplex-local mass coefficient in
(LM1), the tie gives separated maximizers at nearby
unequal quotas. The comparison uses a fixed compact smooth
speed approximating \(b\) in weighted norm and flux,
retaining the strict inequality. Exact volume corrections
of size \(O(t^2)\) after a flow of amplitude \(t\) do not
affect its quadratic gain because the base prices vanish.

In detail, let \(h\) be the chosen speed's mass flux and
let its quadratic gain exceed \(\Lambda_r|h|^2\) by \(\kappa>0\).
At \(r-\epsilon\) use the perturbed partition of masses
\(m_*+th\), with \(t=C\sqrt\epsilon\). Joint Taylor expansion
shows that it beats the simplex-local value by
\(C^2\kappa\epsilon+o(\epsilon)\), since \(D'(r)=0\).
The absence of maximizing partitions in a closed annulus
about the tripod orbit at the tie persists for nearby
parameters by strong maximizing compactness. On the quota
segment from \(m_*\) to \(m_*+th\), the sets admitting
an inner or an outer maximizer are closed and cover the
segment. The balanced endpoint admits only the inner orbit
and the perturbed endpoint admits only outer maximizers.
The two types therefore coexist at an unequal quota.
[Local mass curvature and coexistence, LM1 and LM8](THREE_CELL_LOCAL_MASS_PROFILE.md#consequences-at-a-first-balanced-tie)

If the fixed-shape trace-free curvature is at least the
simplex's and the dilation curvature fails this quota
comparison, then necessarily
\[
           \Lambda_r |d|^2\ge
                \frac{1-r^2}{1+r^2}(B_a+r^2D''(r)).
                                                               \tag{PD6}
\]
This requires a nonzero second radial moment imbalance,
even if the centered prices and covariance derivatives match.
The corrected-dilation bound (G22) can also be substituted
into (PD6); for \(w\) the translation with mass flux \(d\)
and \(h=b-w\) it gives
\[
 D''(r)\ge \frac{8(1-r)}{(1-r^2)^2}
             \left[(1+r)\|h\|_a^2+\langle w,h\rangle_a\right]_+ .
\]
[Balanced dilation and the mixed contact inequality, G21–G22](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#balanced-dilation-at-an-endpoint-tie)

## Zero modes coupled to a volume change

At a pressure-free first tie with the same scalar uniqueness
assumption, suppose a balanced weighted speed \(u\) satisfies
\(Q_0(u)=0\) but does not lie in the unrestricted operator
kernel. Constrained stability gives \(Q_0(u,f)=0\) for every
balanced \(f\), by polarization of the nonnegative balanced
form. Since translations span the volume changes, some
translation \(s_v\) has \(Q_0(u,s_v)\ne0\). For the nonzero
quota direction \(h=C_A^{\rm vol}s_v\), the speeds \(s_v+c u\)
all have flux \(h\), and their quadratic costs are
\[
                  Q_0(s_v)+2c\,Q_0(s_v,u).
\]
These costs have infimum \(-\infty\).
For any specified finite coefficient, choose a fixed \(c\)
and approximate that speed by compact regular speeds with
exact first flux \(h\), using two fixed bumps. Continuity of
the weighted form and flux preserves a strict coefficient
inequality. Realizing the flow gives the quota comparison
above with a quadratic gain exceeding any specified multiple
of \(|h|^2\). This forces separated maximizing orbits at
nearby unequal quotas, also at a first tie with nonzero scalar
slope gap, by choosing \(C\) large enough to dominate that
gap in the scalar Taylor expansion.
[Translation eigenspace and volume-flux continuity](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

## Unrestricted zero modes coupled to covariance

Consider a first tie invisible to the joint first-order
tests, so the centered prices vanish and \(M_A=M_T\).
Suppose a nonzero unrestricted zero mode \(u\) has
\(\langle u,b_B\rangle_a\ne0\) for a symmetric matrix \(B\).
The unrestricted kernel is balanced by the translation
and volume-flux identities. Testing the affine dilation
identity against \(u\) gives
\[
                    \langle g_B,u\rangle
                               =-\langle b_B,u\rangle_a/r .
                                                               \tag{PD7}
\]
Equality of the covariance derivatives gives zero scalar
deficit slope. The mixed contact inequality (G10) then
implies \(\langle g_I,u\rangle=0\). Here (G10) extends to
this weighted balanced speed by compact approximation,
fixed flux corrections, and the continuous pairing in (CT8).
We can replace \(B\) by its trace-free part and retain
a nonzero coupling in (PD7).
[Mixed contact inequality, G10](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)
[Continuous affine covariance pairing, CT8](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions)

There is a finite simplex-local quadratic upper bound for
this covariance change, uniformly in scalar \(\rho\)
near \(r\). For a balanced \(P\) in a fixed small neighborhood
of the simplex orbit, put \(\mathcal L=L_\rho(P;T_P)\),
its pointwise regret with the orientation \(T_P\) in (B11).
The global weighted estimate underlying (B3) gives
\(\sum_i\int w_\alpha |1_{P_i}-1_{T_{P,i}}|d\gamma
\le C\sqrt{\mathcal L}\). On a compact positive correlation
interval the covariance derivative of a Mehler kernel at
\(\rho I\) is that kernel times a quadratic polynomial.
It is bounded by \(C w_\alpha(x)w_\alpha(y)\) as in (B4),
with \(\alpha\) chosen there to absorb the polynomial.
Expand the difference of the two quadratic functionals
to obtain
\[
 |\operatorname{tr}(B[M_P-M_{T_P}])|\le C\sqrt{\mathcal L}.
\]
The simplex covariance derivative is scalar and pairs to
zero with \(B\), independently of the rotation. The
uniform second-order covariance remainder of (AM2), together
with (B11), therefore gives
\[
 S_{\rho I+tB}(P)\le S_\rho(T)-c\mathcal L
                            +C|t|\sqrt{\mathcal L}+Ct^2
                       \le S_\rho(T)+C't^2 .              \tag{PD8}
\]
[Weighted regret and local growth, B3 and B11](POSITIVE_CORRELATION_CONTINUATION.md)

Choose \(c\) so that \(2c\langle g_B,u\rangle\) exceeds the
finite upper coefficient in (PD8) after adding the fixed-shape
covariance curvature of \(A\). Approximate \(u\) by compact
smooth balanced speeds with fixed flux bumps. Since \(Q_0(u)=0\)
and (PD7) is continuous in the weighted space, for the chosen
fixed \(c\) one preserves the strict gain after paying the
shape cost \(c^2 Q_0\) of the approximating speed. Use its
exactly balanced flow of amplitude \(ct\) at covariance
\((r-\epsilon)I+tB\), with \(t=\sqrt\epsilon\).
The first scalar correlation derivatives are equal at the tie,
the first trace-free derivatives vanish, and the mixed
shape-covariance Taylor term is the chosen gain. Thus the
ordinary competitor beats (PD8) for small \(\epsilon>0\).

At fixed \(r-\epsilon\), continue from zero trace-free
perturbation to \(tB\). With a unique balanced simplex
immediately to the left of the original scalar tie, the
same annulus and compactness argument gives coexistence of
separated balanced maximizing orbits at an interior
anisotropic covariance of mean correlation \(r-\epsilon\).

These tests leave open invisible ties whose covariance
curvature and ordinary interface spectrum satisfy the
displayed necessary restrictions.

## Higher-order ties with three isolated ends

There is a reflection-free exclusion when all three pairs have
well-separated ends. Suppose \(A\) ties the simplex at \(r\),
distinct simplex-optimal correlations accumulate at \(r\), and
\(D'(r)=D''(r)=0\). For each pair of labels assume there are
orthonormal coordinates \((t,z)\) about the origin with the
following two properties.

- In the full half-plane \(t>S_0\), its interface is one finite
  bounded graph \(z=h(t)\), with its lower and upper labels
  occupying the two sides up to null sets.
- Some open angular cone about the negative longitudinal
  direction, outside a bounded set, belongs to the third
  label up to null sets.

Under these assumptions \(A\) is the simplex. No reflection or
eventual monotonicity of the bounded heights is required.

The conditions hold, for example, when the jump network outside a
compact set has exactly three proper arms, one for each pair, and
the arms are bounded transverse graphs in pairwise obtuse outward
directions. In fact, the other two arms have longitudinal
projection tending to \(-\infty\) in the direction of a fixed
arm, so its far half-plane contains only that arm. On a large
circle the separated directional neighborhoods of the three arms
have the prescribed pair transitions. The sector between the two
other directions and opposite the fixed arm therefore carries
the third label. The reverse of the fixed direction is strictly
between those other directions, by pairwise obtuseness. A narrow
cone about the reverse direction is disjoint from all three arm
neighborhoods and stays in that sector.

To prove the assertion, apply (G10) to arbitrary compact
zero-flux speeds on each individual pair. Because \(D''(r)=0\),
the correlation derivative \(g_e\) of each unadjusted pair score
is constant along its regular interface. The constants agree
across components of the same pair, by testing with two bumps
on different components whose integrals cancel.
[Mixed contact inequality, G10](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

In the coordinates of one end, orient the contrast \(F\) as the
score of the upper label minus that of the lower label. Write
\(\beta\) for its stationary value, \(\sigma^2=1-r^2\),
\(\psi(v)=2\Phi(v/\sigma)-1\), and \(G(z)=\psi(rz-b)\)
for the binary score of the affine pressure line with
\(G(b)=\beta\). The bounded isolated-end theorem gives
\(h(t)=b+u(t)\) with \(u,u'\to0\) and bounded derivatives.
Gaussian flattening strengthens the convergence to
\(|u(t)|\le C e^{-c t^2}\).
[Bounded isolated-end convergence](THREE_CELL_ISOLATED_END_RIGIDITY.md#bounded-ends-converge-to-the-pressure-offset)
[Gaussian flattening](THREE_CELL_TAIL_FLATNESS.md)

At the graph, the fixed-observation correlation derivative of
\(F\) tends to that of the affine score at \(z=b\). Indeed the
source difference from the affine halfspace consists of an
exponentially thin strip for large positive longitudinal
coordinates and a bounded contrast on the remaining left
half-plane. Differentiating the conditional Gaussian density
inserts a polynomial in source and observation coordinates.
Integrating over the strip or the left half-plane makes
the differentiated difference tend to zero. The difference
between evaluating the affine derivative at \(b+u(t)\) and
at \(b\) tends to zero by smoothness. Constancy of \(g_e\)
therefore gives \(LF=L_zG(b)\) along the graph, where
\(L=\partial_{tt}+\partial_{zz}-t\partial_t-z\partial_z\)
and \(L_z=\partial_{zz}-z\partial_z\).

Suppose \(\beta\ne0\), so \(b\ne0\). Set \(U=F-G\).
In the half-plane of the given end the decomposition is
\[
 U(t,z)=\int_{s>S}
       [\psi(rz-b-u(s))-\psi(rz-b)]p_\sigma(s-rt)\,ds + R(t,z),
\]
where \(S\) is large and \(p_\sigma\) is the centered normal
density of standard deviation \(\sigma\).
The remainder and its first two vertical derivatives are
\(O(\exp[-(rt-S)^2/(2\sigma^2)])\) for bounded \(z\), by
the conditional left tail and the bounded vertical Gaussian
derivative moments. This is the binary case of (CR2), with
no opposite arm in the half-plane. The computation in
(CR1), (CR4), and (CR5) uses only this decomposition, the
constant value of \(LF\), and the graph equation \(F=\beta\).
In particular, with the coefficient \(k_0(z)\) from (CR4),
the operator \(-\partial_{zz}+z\partial_z-k_0(z)\) cancels
the linear term in \(u(s)\) in the strip integrand.
The quadratic term at \(z=b\) is
\(-br\psi'((r-1)b)u(s)^2/[2(1+r)^2]\).
[Strip calculation and differential equation, CR1–CR5](THREE_CELL_REFLECTED_CONSTANT_CORRELATION_RIGIDITY.md#the-signed-quadratic-term)

If \(u\) vanishes on a tail, the end is straight.
Otherwise choose \(S\) large enough to fix the strict
quadratic sign, then keep a fixed interval \([a,a']\)
with \(a>S\) where \(u(s)^2\) is bounded below by a
positive constant. Its contribution for large \(t\)
dominates the left error by a Gaussian shift, since
\(p_\sigma(a-rt)/\exp[-(rt-S)^2/(2\sigma^2)]\to\infty\).
Thus the forcing in (CR5) has one strict sign on the tail.
Its first derivative coefficient is \(-t+O(1)\) with
bounded derivative, and its zeroth-order coefficient is
bounded. Removing the first derivative gives negative
zeroth-order coefficient on a tail. The maximum principle
on sign components, as in (CR5), makes \(u\) eventually of
one weak sign. The isolated one-sided end theorem then gives
exact straightness.
[One-sided isolated ends](THREE_CELL_ISOLATED_END_RIGIDITY.md#a-one-sided-approach-forces-exact-straightness)

Analyticity of the pair score continues the constant value
\(\beta\) along the whole supporting line \(z=b\) of this
straight end. When \(t\to-\infty\) on that line,
the conditional sample lies in the reverse third-label
cone with probability tending to one. Its pair contrast
tends to zero, a contradiction to \(\beta\ne0\).
Thus all three pair pressures vanish. The corresponding
affine halfspaces have \(b=0\), so their fixed-observation
correlation derivatives at the pressure line are zero.
Constancy gives \(g_e=0\) on every pair. The dilation
bound (G6) implies \(x\cdot N=0\) almost everywhere
on the interface, making \(A\) a cone about the origin.
Conical rigidity completes the proof.
[Dilation bound and its geometric consequence, G6 and G19](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

## Eliminating the reverse-cone hypothesis for \(r>1/2\)

In the preceding exclusion, the reverse-cone assumption can
be dropped when \(r>1/2\). One bounded isolated graph end in
a full forward half-plane for each pair suffices. To see this,
write \(e_{ij}\) for the unit longitudinal direction of one
such arm of each pair. Distinct directions have nonpositive
dot product. Otherwise the second arm would eventually
enter the first pair's isolated half-plane despite involving
its absent label. Boundedness of the transverse coordinates
and the two positive-mass sides of the second arm make that
impossible.

If no two directions are antipodal, the reverse direction
\(-e_{ij}\) has a strictly positive projection on at least
one of the other arms. If both projections are positive,
the overlap of their two forward half-planes supplies the
required third-label cone. If just one projection is
positive, \(-e_{ij}\) lies inside that arm's forward
half-plane and is distinct from its longitudinal direction.
Its transverse projection within that half-plane is
therefore nonzero. Boundedness of that arm's graph makes a
small cone about \(-e_{ij}\) a single label there. Tilt the
direction slightly into the other arm's forward half-plane
to identify that label as the third one. Hence the preceding
reverse-cone proof applies.

It remains to consider an antipodal pair. The third direction
must be perpendicular to both. After rotation and relabeling,
take the \(1/2\) direction to be rightward, the \(2/3\)
direction upward, and the \(1/3\) direction downward.
In the upper forward half-plane, label 2 is on the right
of the bounded graph and label 3 on the left. The orientation
follows from its overlap on the far right with the \(1/2\)
half-plane, where only the common label 2 is permitted.
Likewise label 1 is on the right and label 3 on the left in
the lower half-plane.

Let \(p,q\) be the limiting horizontal offsets of the
upper and lower graphs. Orient the two contrasts as \(U_2-U_3\)
and \(U_1-U_3\), respectively. Their pressure values are
\(\psi((r-1)p)\) and \(\psi((r-1)q)\).
As in the preceding proof, the constant correlation derivative
and the isolated end make any nonzero-pressure end exactly
straight. If both pressures are nonzero, evaluate the upper
straight pair score on \(x=p\) at heights tending to
\(-\infty\). Its limit is \(-\Phi((q-rp)/\sigma)\), since
only labels 1 and 3 are present in the lower half-plane.
Whole-line stationarity therefore gives \(p>0\) and
\[
        \Phi((q-rp)/\sigma)
               =2\Phi((1-r)p/\sigma)-1
               <\Phi((1-r)p/\sigma).
\]
This forces \(q<p\), whereas the same argument along
the lower straight line forces \(p<q\).

If exactly one of the two pressures is nonzero, take it
to be the upper one. Then \(q=0\) and the upper graph
is straight. The same comparison requires \(p>0\), and
with \(u=(1-r)p/\sigma\) and \(s=r/(1-r)>1\) gives
\[
                       2\Phi(u)-1=\Phi(-su).              \tag{PD9}
\]
The fixed-observation correlation derivative on \(x=p\)
at positive infinite height tends to \(2u\phi(u)/\sigma^2\),
the derivative of its binary pressure score. At negative
infinite height it tends to \(p\phi(su)/\sigma^3\),
the derivative of \(-\Phi(-\rho p/\sqrt{1-\rho^2})\) at
\(\rho=r\) with \(p\) fixed. For this reverse limit,
Gaussian flattening of the lower end makes the differentiated
strip error tend to zero even though the lower arm has
zero pressure and need not be straight. Constancy of the
derivative along the original straight arc extends to
the whole line by analyticity. Thus
\(\phi(su)=2\phi(u)/(1+s)\), which gives
\[
          u^2=\frac{2\log((1+s)/2)}{s^2-1}
                           \ge\frac{1}{s(1+s)}.
\]
For the inequality use
\(2\log((1+s)/2)\ge(s-1)/s\) for \(s>1\), obtained
by differentiation from equality at \(s=1\). But
\[
 2\Phi(u)-1\ge 2u\phi(u)
           \ge\frac{\phi(su)}{su}>\Phi(-su),
\]
where the last inequality is the strict Gaussian tail
bound. This contradicts (PD9).

The only remaining case has both vertical pressures zero,
so all pair pressures are zero. Then all \(g_e\) are zero
by their affine tail limits, and (G6) again gives the
centered simplex. This proves the strengthened exclusion.
