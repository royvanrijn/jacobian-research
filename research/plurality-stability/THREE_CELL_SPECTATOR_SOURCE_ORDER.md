# Uniform disorder of the spectator source at right contact

Fix a compact interval \(I\subset(0,1)\) of correlations. Let \(A\)
be a balanced planar Gaussian maximizer at \(r\in I\) outside the
centered simplex orbit, with equal stationary prices and the three
exact simplex pair-score lines. Write \(T\) for the reference
simplex, and assume the fixed-partition correlation derivatives
satisfy \(S'_r(A)\le S'_r(T)\). In particular these hypotheses hold
for a nonsimplex maximizing compact discrepancy at a tie approached
from the right by simplex-optimal correlations.

At least one reference line has a spectator source whose tilted
Gaussian \(L^1\) distance from every bounded one-sign-change profile,
with values in \([-2,1]\), is at least \(\varepsilon_I>0\).
The obstructing score signs occur in a fixed bounded window
depending only on \(I\). This separation does not require a bound on
the compact discrepancy radius. It also applies when the germ at the
origin is antipodal.

## A uniform excursion on one of the three germ arms

Let \(\tau_k\) be the outward reference ray direction of the pair
\(i,j\) opposite spectator \(k\). Write \(U_i=T_r1_{A_i}\),
\(\sigma^2=1-r^2\), and
\[
 h_k(s)=U_i(s\tau_k)-U_k(s\tau_k)=U_j(s\tau_k)-U_k(s\tau_k).
\]
Because the three exact pair equalities make the scores tie at zero,
the first derivatives of the pair contrasts have a common
nonzero proportionality factor to the reference pair normals.
Let \(\eta\) be its sign. As in (CG1), the active germ half-rays
are those parametrized by \(s=\eta u\), \(u>0\).
Uniform central gradient separation and a uniform upper bound on
the Gaussian-smoothed score Hessian give constants \(s_0,m_0>0\)
depending only on \(I\) such that for each spectator
\[
                    h_k(-\eta s_0)<-m_0<m_0<h_k(\eta s_0).
                                                               \tag{SO1}
\]
The proof of (CG1) uses the exact lines and equal prices at this
step. Exterior agreement is one way to obtain these hypotheses.
[Centered germ and the score gradient (CG1)](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius)
[Reference score lines and the tilted moments (RM1)–(RM3)](THREE_CELL_REFERENCE_MIRROR_MISMATCHES.md#the-line-equalities-determine-the-tilted-transverse-first-moments)

There are \(L>s_0\) and \(m_1>0\), also depending only on \(I\),
such that some \(k\) and \(s_1\in[s_0,L]\) satisfy
\[
                             h_k(\eta s_1)\le -m_1 .         \tag{SO2}
\]
Choose \(L\) to meet the uniform finite-arm comparison in
(OA6)–(OA8). The analogous comparison holds for any three straight
pair segments issuing from zero that belong to the respective
interfaces up to arclength-null sets: the proof retains only the
surface integrals over the three segments of length \(L\).
It forces \(S'_r(A)>S'_r(T)\) for a nonsimplex maximizer
with those segments.

To prove (SO2), if it fails, take a sequence for which all three
minima over \([s_0,L]\) tend to nonnegative values. Fix the reference
lines by rotation and labels and fix the sign \(\eta\) by a subsequence.
Strong maximizer compactness gives a limiting maximizer at a
correlation in \(I\). The limit remains nonsimplex by uniform
measurable local isolation of the simplex. The full-line equalities,
the zero-price winner rule, and the derivative inequality pass
to the limit by strong convergence. For the derivative this uses
the bounded OU multipliers \(j r^{j-1}\), uniformly on \(I\).
The uniform germ sign and the limiting inequalities give
\(h_k(\eta u)\ge0\) on all three segments from zero to \(L\).
The analytic functions \(h_k\) are nonzero by (SO1). Their zeros
on each compact segment are isolated. Off those zeros, the pair
contrast has an identically vanishing tangential restriction and a
nonzero analytic normal derivative except at isolated points, because
at zero that normal derivative is nonzero. Thus each segment lies
on its actual pair interface up to arclength-null sets. This
contradicts the finite-arm comparison. The argument needs only
the exact score lines for the limiting partition. It does not
require the limit to have a compact discrepancy or to satisfy
a right-accumulation assumption.
[Finite-arm comparison](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#arcs-close-to-the-origin-and-endpoint-implications)
[Maximizer compactness and local isolation](POSITIVE_CORRELATION_CONTINUATION.md)

## A quantitative obstruction to one source crossing

Use coordinates \(t,v\) along and normal to the chosen reference
line, and define
\[
 q_k(t)=\int 1_{A_k}(t\tau_k+v n_k)\,p_\sigma(v)\,dv,\qquad
 f(t)=1-3q_k(t),
\]
where \(n_k\) is a unit normal and \(p_\sigma\) is the centered
Gaussian density of variance \(\sigma^2\). Because
\(U_i=U_j=(1-U_k)/2\) on the line, the exact identity is
\[
 h_k(s)=\tfrac12(p_\sigma*f)(rs)=
 \tfrac12 e^{-r^2s^2/(2\sigma^2)}
    \int e^{rst/\sigma^2} f(t)\,d\gamma_{\sigma^2}(t).
                                                               \tag{SO3}
\]
Here \(\gamma_{\sigma^2}\) is the one-dimensional centered Gaussian
probability of variance \(\sigma^2\). This convolution is the
Gaussian case of Schoenberg's variation-diminishing theorem.
[Schoenberg, *Variation-diminishing integral operators of the convolution type*](https://acta.bibl.u-szeged.hu/13626/1/math_012_pars_b_097-106.pdf)

The one-crossing consequence used here has a direct proof. Say \(g\) has at most one
essential change if for some threshold \(b\) it has one sign to the
left and the opposite sign to the right, allowing zero values on
either side. Include single-sign functions in this class. For a
threshold profile, \(\int e^{u(t-b)}g(t)d\gamma_{\sigma^2}(t)\)
is monotone in \(u\) since its derivative is the transform of
\((t-b)g(t)\), which has one sign. Thus the convolution of \(g\)
cannot have alternating strict signs at three ordered points.
The conclusion is immediate for a single-sign \(g\). It does
not require the comparison profile \(g\) to be centered.

Choose \(k,s_1\) as in (SO2). The values of \(h_k\) at the three
points \(-\eta s_0,\eta s_0,\eta s_1\), read in either direction
along the line, have signs negative, positive, negative, with
magnitudes at least \(m=\min(m_0,m_1)>0\).
For \(f,g\) taking values in \([-2,1]\) and \(|s|\le L\),
Cauchy–Schwarz in (SO3) gives
\[
 \tfrac12 |(p_\sigma*(f-g))(rs)|
   \le \tfrac{\sqrt3}{2}e^{r^2 L^2/(2\sigma^2)}
                    \|f-g\|_{L^1(\gamma_{\sigma^2})}^{1/2}.
\]
The factor on the right has a uniform upper bound \(C_I\).
If \(g\) has at most one essential sign change, its transform
must differ by at least \(m\) from that of \(f\) at one of the
three points. Hence
\[
                       \|f-g\|_{L^1(\gamma_{\sigma^2})}
                                  \ge (m/C_I)^2 .           \tag{SO4}
\]
In particular the obstruction cannot consist of sign oscillations
with arbitrarily small tilted mass. The same estimate divided
by three bounds the distance of \(q_k\) from profiles in \([0,1]\)
crossing \(1/3\) at most once.

## Three separated source intervals

Let \(\varepsilon=(m/C_I)^2\) and \(\delta=\varepsilon/4\), and
write \(\mu_+(E)=\int_E f_+\,d\gamma_{\sigma^2}\) and
\(\mu_-(E)=\int_E f_-\,d\gamma_{\sigma^2}\), where
\(f_+=\max(f,0)\) and \(f_-=\max(-f,0)\). For each threshold
and each order of the two signs, the closest threshold profile
agrees with \(f\) wherever its sign is permitted and is zero
elsewhere. Thus the sum of the two wrong-sign masses at every
threshold is at least \(\varepsilon\). Comparison with single-sign
profiles also shows \(\mu_+(\mathbb R),\mu_-(\mathbb R)\ge\varepsilon\).

Choose \(b_1<b_2\) such that the left interval \((-\infty,b_1)\)
and the right interval \((b_2,\infty)\) each have positive-part
mass \(\delta\). The middle interval has positive-part mass at
least \(2\delta\). If it has negative-part mass at least \(\delta\),
the three intervals have sign masses at least \(\delta\) in the
order positive, negative, positive. If both outer intervals have
negative-part mass at least \(\delta\), use the reverse order.
Otherwise the middle and one outer interval have negative-part
mass less than \(\delta\) each. Require the positive sign on
their union and the negative sign on the opposite outer interval.
The wrong-sign mass is then less than \(3\delta\), a contradiction.

Both measures have uniformly bounded Lebesgue densities on the
correlation band and uniformly small tails. Intersect the three
intervals with a common bounded window and remove a fixed small
neighborhood of their two cuts. The window and removal width can
be chosen using only \(I\) so that each interval retains at least
\(\delta/2\) of the required mass. This gives three bounded source
intervals with a uniform positive separation and alternating
sign masses bounded below. Whether a balanced transfer can turn
that disorder into strict gain remains open.
