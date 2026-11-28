# Fixed-measure moments and X-ray constraints for compact simplex discrepancies

Fix \(0<r<1\) and a balanced planar Gaussian maximizer \(A\)
that agrees almost everywhere outside a disk with a centered regular
simplex \(T\). Let \(\sigma^2=1-r^2\) and \(\mu=N(0,\sigma^2 I_2)\).
The three reference pair lines of \(T\) are full score-equality
lines at \(A\), and the centered prices are zero. Thus on each
reference pair line, deconvolution gives a vanishing transverse
\(\mu\)-weighted contrast on every longitudinal slice, almost
everywhere. The same statements hold for \(T\).
[Full reference-line equalities and weighted projections](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

These identities give a fixed-measure formula for the scalar
correlation slope at a compact tie. They also give constraints
on a differentiable family of compact stationary shapes. In
particular a nonsimplex branch needs a nonzero speed after
rescaling to a fixed conditional variance. That speed must
simultaneously satisfy the linearized stationary equations and
all three homogeneous weighted projection equations.

## A fixed-measure moment series at a tie

For a partition \(P\) and integer \(n\ge0\), define
\[
 \begin{split}
 M_{n,i}(P)&=\int_{P_i}x^{\otimes n}\,d\mu(x),\\
 H_P(s)&=\sum_{n\ge0}\frac{(s/\sigma^2)^n}{n!}
                          \sum_i\|M_{n,i}(P)\|^2 ,
               \qquad |s|<1 .
 \end{split}
\]
Here \(x^{\otimes0}=1\), and the tensor norm is Euclidean.
Put \(\lambda(s)=\sqrt{1-s^2}/\sigma\). Changing variables
\(x=\lambda(s)x'\), \(y=\lambda(s)y'\) in the planar joint
Gaussian gives
\[
      S_s(\lambda(s)P)
          =(1-s^2)\iint 1_{\{P(x)=P(y)\}}
                          e^{s x\cdot y/\sigma^2}\,d\mu(x)d\mu(y)
          =(1-s^2)H_P(s).                                \tag{XJ1}
\]
The series follows from the tensor expansion of the exponential.
Absolute convergence and differentiation on compact subintervals
of \((-1,1)\) follow by bounding \(|x\cdot y|\) by
\((|x|^2+|y|^2)/2\), with a polynomial factor when differentiating.

For \(P=A,T\), the first shape derivative of
\(S_r(\lambda P)\) at \(\lambda=1\) vanishes. Along an
oriented active interface \(i,j\) it is the integral
of twice the pair score difference against the dilation speed
and Gaussian surface measure. The difference is zero since
the stationary prices are zero. The speeds away from the
exceptional disk vanish because the exterior rays are radial.
In particular this shape derivative vanishes even if the
dilation changes the Gaussian masses.
[First variation and the dilation speed](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Let \(D(s)=S_s(T)-S_s(A)\) for the fixed partitions and
assume \(D(r)=0\). Differentiating (XJ1) at \(r\) and
using the vanishing shape derivatives yields
\[
 D'(r)=\sigma^2[H_T'(r)-H_A'(r)] .                       \tag{XJ2}
\]
Consequently a tie with simplex-optimal correlations
approaching from the right requires \(H_T'(r)\ge H_A'(r)\).
Equivalently, let \((X,Y)\) be a standard correlated
Gaussian pair at correlation \(r\) and let \(E_P\) be
their same-label event for \(P\). The same differentiation
gives
\[
             S'_r(P)=
      \frac{\mathbb E[1_{E_P} X\cdot Y]-2rS_r(P)}{\sigma^2}.
                                                               \tag{XJ3}
\]
At the tie the right-side slope condition compares
the conditional mean of \(X\cdot Y\) on the two agreement
events and requires it to be at most the reference mean for \(A\).

The constant terms in \(H_T-H_A\) agree because all three
posteriors tie at the origin and \(U_i(0)=\mu(A_i)=1/3\).
The linear coefficient of \(H_T-H_A\) is strictly positive
if \(A\ne T\). To see this, let \(v_i\) be the unit centers
of the three reference sectors. The first test moment in each
longitudinal projection makes \(M_{1,i}(A)-M_{1,j}(A)\)
parallel to \(v_i-v_j\). The three edge vectors sum to
zero and any two are independent, so the three scalar
multipliers coincide. The total centroid is zero, hence
\(M_{1,i}(A)=c v_i\) for a common scalar \(c\).
The corresponding reference scalar is \(c_T>0\).
Pointwise maximization of \(v_i\cdot x\) gives \(3c\le3c_T\)
with equality only for \(T\). Pointwise minimization gives
\(3c\ge-3c_T\), with equality only for its antipode, by
central symmetry of \(\mu\). A compact discrepancy with
the stated exterior cannot be the antipode, so \(|c|<c_T\)
unless \(A=T\). The positive linear coefficient alone
does not determine the sign of the derivative at a zero
of the full series.

## The projection and interface equations for a moving branch

Suppose \(A_s\) is a family of ordinary balanced stationary
partitions at correlations \(s\) near \(r\), with \(A_r=A\).
Assume every \(A_s\) agrees with the same centered simplex
outside one fixed disk. Assume regular active interfaces
at the base are analytic arcs with nonvanishing pair-score
gradient. Use a differentiable motion of these arcs in
the stationary equations. Assume its bounded normal
speed \(u_e\) is supported in a common compact set and its
indicator derivative is the locally finite signed
interface measure with density \((d_e)_i u_e\) for label \(i\),
where \(d_e=e_i-e_j\) and the normal \(N_e\) points out
of label \(i\). For example these assumptions hold for
a finite regular network in the disk whose arcs and
junctions move differentiably with bounded normal speeds
and whose rays outside the disk stay fixed. The stated
indicator derivative excludes additional singular
contributions at junctions. These hypotheses also allow
a one-sided derivative in the correlation parameter.

Set \(B_s=\lambda(s)^{-1} A_s\), with the scale from (XJ1).
The centered prices of every \(A_s\) are zero by exterior
agreement. In the coordinates \(x=t\tau_{ij}+z n_{ij}\)
of the reference line for labels \(i,j\), the transverse
constraint for \(B_s\) then has the fixed centered Gaussian
variance \(\sigma^2\), because
\(\lambda(s)g_{\sqrt{1-s^2}}(\lambda(s)z)=g_\sigma(z)\).
Here \(g_\sigma\) is the scalar centered Gaussian density of
standard deviation \(\sigma\). Let \(b_e=x\cdot N_e\) on
the interfaces of \(A\). The rescaled normal speed is
\[
                         v_e=u_e+\frac r{\sigma^2}b_e .
\]
Differentiate the projection constraint against a compactly
supported longitudinal test \(\psi\). For each of the three
reference pairs \(i,j\) this gives
\[
 \sum_e[(d_e)_i-(d_e)_j]
   \int_{\Sigma_e}\psi(x\cdot\tau_{ij})g_\sigma(x\cdot n_{ij})
                                  v_e(x)\,ds(x)=0.        \tag{XJ4}
\]
The integral uses Euclidean arclength. It is the distributional
projection of the compactly supported indicator derivative,
with the transverse density as weight.

Write \(K_a=a^{-1}K\) in the weighted interface notation of the
maximizer variation, \(Q_0=I-K_a\),
\(a_e=|\nabla(T_r1_{A_i}-T_r1_{A_j})|\), and
\[
          g_e=\left.\partial_s T_s(1_{A_i}-1_{A_j})\right|_{s=r}
                                \quad\text{on }\Sigma_e.
\]
Differentiating the zero-pressure equality on a moving interface
gives \(g_e+(K u)_e-a_e u_e=0\). The convolution derivative
contributes \(K u\) and moving the observation in the outward
normal contributes \(-a_e u_e\). The dilation identity gives
\(g/a=-(I-r^2K_a)b/(r\sigma^2)\). Substitution yields the
interface equation for the same rescaled speed as in (XJ4),
\[
                               Q_0 v=-b/r .              \tag{XJ5}
\]
Each speed \(v_e\) is analytic along its connected regular
arc. Indeed \(u_e=(g_e+(K u)_e)/a_e\) there, where
\(g_e+(K u)_e\) is the restriction of an analytic Gaussian
convolution expression, \(a_e>0\) is analytic along the arc,
and \(b_e\) is analytic there. Equations (XJ4) and (XJ5)
therefore constrain the analytic speeds on the individual arcs.
[The mixed dilation identity, MR1](THREE_CELL_MIXED_DILATION_RESOLVENT.md#the-mixed-variation-and-the-mass-flux)

## A visibility condition that excludes motion

A regular actual arc on a reference line has \(v_e=0\).
At a regular point with a strict spectator gap, the full-line
equalities persist for all members of the family. Nonvanishing
of the transverse score derivative persists locally, so the
regular moving interface stays on that reference line and
its normal speed is zero. Dilation is also tangent there.

For the remaining regular arcs, use the following sufficient
condition on the base geometry. Suppose they form a finite
collection of connected analytic arcs. At each step, after
ignoring arcs already proved to have zero speed, choose a
reference pair and an open interval in its longitudinal
projection such that the projected locations in that interval
have exactly one remaining transverse crossing with nonzero
pair incidence coefficient. Require this crossing to be on
one selected arc. If successive steps remove every remaining
arc, no such moving stationary branch with a nonsimplex base
is possible.

Indeed, restrict (XJ4) to the selected projection interval.
Change variables on the single transverse crossing. The
coefficient, transverse Gaussian density and projection
Jacobian are nonzero, so the speed on that arc vanishes.
Analytic continuation along the arc lets it be removed
from all later projection equations. If every arc is removed,
\(v=0\), and (XJ5) gives \(b=0\). The indicators then have
zero distributional radial derivative and are constant
along almost every radial direction. The exterior fixes
their labels as \(T\), contrary to a nonsimplex base.
The visibility requirement is a sufficient geometry
hypothesis for this branch exclusion. It is open for
arbitrary compact stationary shapes, and the moving-branch
argument makes no exclusion of an isolated tie.
