# Finite-amplitude domination of the adverse radial mode

This note proves a direct comparison for **every** ordinary balanced Gaussian
three-cell partition in every dimension at
\(0<\rho\le221/500=0.442\). The new interval is
\([11/25,221/500]\). It isolates one adverse radial direction, uses the
simplex constraints to limit its amplitude, and pays for its positive
contribution with unavoidable angular loss. There is no assumption of
closeness to a simplex, stationarity, or continuation of maximizers.

The numerical inputs are a full-domain directed certificate, with analytic
tails, from [the checker](check_three_cell_finite_amplitude.py). Its
[verification record](three_cell_finite_amplitude_verification.json) retains
the entire cover and source pins. The result has internal checking only;
independent replay, formal verification and external review are unperformed.
The all-positive-correlation target remains open.

## FA1. Statement and imported comparison

For every \(d\ge2\), every measurable partition
\(A=(A_1,A_2,A_3)\) of \(\mathbb R^d\) with
\(\gamma_d(A_i)=1/3\), and every \(0<\rho\le221/500\),

\[
 S_\rho(A)\le S_\rho(T),\qquad
 \widetilde Q_\rho(1_A-1_T)\le2L_\rho(A;T).
 \tag{FA1}
\]

Here \(T\) is a regular three-cell simplex partition. DG10 of the
[direct comparison](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md) supplies
\(0<\rho\le11/25\). We prove the planar comparison on the new interval.
The same maximizer existence and dimension reduction as in DG4 then prove
the comparison for every competitor in each \(d\ge2\). These are
Heilman--Tarter, Lemma 2.3 and Theorem 1.9,
[Forum of Mathematics, Sigma 9 (2021), e65](https://doi.org/10.1017/fms.2021.56);
their hypotheses are checked in DG4. Only an actual maximizer is reduced
to a plane. The exact expansion (B1) in the
[continuation note](POSITIVE_CORRELATION_CONTINUATION.md) gives the second
inequality; no continuation theorem is used.

Write \(t=r^2/2\), \(d\nu(t)=e^{-t}dt\),
\(p_i(t)=\sigma\{\theta:(\sqrt{2t},\theta)\in A_i\}\), and
\(z_i=p_i-1/3\). Thus \(p(t)\in\Delta_3\),
\(\sum_i z_i=0\), and \(\int z_i\,d\nu=0\).
Use exactly \(K_+,J_+,j_1\) of (DG15), and the nonnegative elementary
moment penalty \(\Gamma\) of (DG26), including its boundary convention.
Put \(d=1-\rho^2\), \(\kappa=\rho^2/d\), and

\[
 V(t,s)=\frac{a}{d}\exp[-\kappa(t+s)-a],\qquad
 a=2\rho\sqrt{ts}/d,\qquad
 q(z)=\langle z,\mathcal K_+z\rangle-2\int J_+z^2d\nu.
\]

The unrestricted imported comparison (DG29) is

\[
 S_\rho(A)-S_\rho(T)
 \le\sum_i q(z_i)-\mathcal L(p),\qquad
 \mathcal L(p)=\iint V(t,s)\Gamma(p(t),p(s))d\nu(t)d\nu(s)\ge0.
 \tag{FA2}
\]

DG6 proves moment closure for arbitrary angular shapes, then applies
Friedberg--Luttinger Theorem I to the nonnegative monotone remainder of
the circular kernel. Heilman, [arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1),
Corollary A.3, supplies the scalar envelope. The exact theorem statements
and hypotheses were inspected for DG4/DG6. We reuse that argument and
do not assume an independently optimized arc profile is a partition.
This proof needs no new external spectral theorem.

## FA2. One direction and a uniformly negative remainder

Let \(\ell_\rho\) be the exact binary64 polynomial `CENTER` in the
retained positive-certificate source, with no approximation property
assumed. On \(L^2(\nu)\) define

\[
 \mathcal A(t,s)=
 \frac{K_+(t,s)-\ell_\rho(\sqrt t)-\ell_\rho(\sqrt s)}
      {\sqrt{J_+(t)J_+(s)}}.
 \tag{FA3}
\]

It is compact and self-adjoint by the Hilbert--Schmidt bound below.
Define an auxiliary trial, not a feasible mass profile, by

\[
 f(t)=t^{-1/2}\sum_{j=0}^7 b_j e^{-\lambda_j t}-M,\qquad
 M=\sqrt\pi\sum_j b_j/\sqrt{1+\lambda_j},
 \tag{FA4}
\]

where

\[
\begin{split}
 (\lambda_j)&=(0,1/16,1/4,1,4,16,64,256),\\
 (b_j)&=10^{-6}(-19606173,19029196,1042060,528898,
                                  3055,7376,-4335,-76).
\end{split}
\]

Then \(\int f\,d\nu=0\) exactly. Set
\(D=\int J_+f^2d\nu\), \(\phi=f/\sqrt D\), and
\(u=\sqrt{J_+}\phi\). The pole of \(f\) is permitted because
\(J_+(t)\asymp\sqrt t\) near zero, making \(u\in L^2(\nu)\).
The actual masses remain bounded. Uniformly on
\(11/25\le\rho\le221/500\), FA7 certifies

\[
\begin{split}
 268/1000\le D\le271/1000,&\qquad
 399/200\le a_u:=\langle u,\mathcal Au\rangle\le20031/10000,\\
 \|(\mathcal A-2)u\|^2\le1/50000,&\qquad
 \|\mathcal A\|_{HS}^2\le217/50.
\end{split}
 \tag{FA5}
\]

Let \(P\) project onto \(u^\perp\). Orthogonal block decomposition
gives \(\|P\mathcal AP\|_{HS}^2\le\|\mathcal A\|_{HS}^2-a_u^2\), so

\[
 \|P\mathcal AP\|\le\sqrt{217/50-(399/200)^2}<3/5.
\]

For each label decompose
\(y_i=\sqrt{J_+}z_i=\alpha_i u+g_i\), \(g_i\perp u\).
Put \(A=|\alpha|\) and \(G^2=\sum_i\|g_i\|^2\).
Since \(\sum_i y_i=0\), we have \(\sum_i\alpha_i=0\).
Mean-zero balance removes additive centering from each quadratic form.
The cross term satisfies
\(|\langle g_i,\mathcal Au\rangle|
 =|\langle g_i,(\mathcal A-2)u\rangle|\le\|g_i\|/\sqrt{50000}\).
Consequently

\[
 \sum_i q(z_i)\le eA^2+2RAG-\delta G^2,\qquad
 e=31/10000,\quad R^2=1/50000,\quad\delta=7/5.
 \tag{FA6}
\]

This isolates one direction without asserting it is an exact eigenvector.
If the scalar sum is negative, (FA2) already proves the comparison.
If it is nonnegative and \(A=0\), then \(G=0\) and all masses are
balanced on almost every circle. In the remaining case \(A>0\), (FA6)
implies

\[
 G^2/A^2\le H^2,
 \qquad H=\frac{R+\sqrt{R^2+\delta e}}\delta,
 \qquad H^2<0.002537.
 \tag{FA7}
\]

## FA3. Feasibility limits the adverse amplitude

The auxiliary pole forces a large correction if its amplitude is too
large. The certified local estimates, for \(0<x\le1/10\), are

\[
 e^{-x^2}\frac{J_+(x^2)}x\ge j_*:=107/1000,\qquad
 \phi(x^2)\ge c_*/x-C_*,\qquad c_*=47/25,\ C_*=14/25.
 \tag{FA8}
\]

Here is an analytic lower bound underlying the first estimate, rather
than a numerical quadrature near the pole. Integration by parts in
\(w\) of (DG15), using the nonnegative primitive
\(v-3v^2/2\), proves
\(w(a)\ge j_1 a e^{-a/2}/I_0(a)\). Therefore, with
\(q=\rho x/\sqrt d\),

\[
 \frac{J_+(x^2)}x\ge
 j_1\frac\rho{\sqrt d}e^{-\kappa x^2}
    [S(q/2)+S(q)/2],\quad
 S(q)=\sqrt\pi(1+2q^2)e^{q^2}\operatorname{erfc}(q)-2q
     =4\int_0^\infty v^2e^{-v^2-2qv}dv.
 \tag{FA9}
\]

The function \(S\) is positive and decreasing. Replace \(x\) by
\(1/10\) in its arguments and in the decreasing exponential, retain
the marginal factor \(e^{-x^2}\), and enclose the correlation interval
by Arb. This gives the first bound in (FA8).
For the second, use
\(\sum b_j e^{-\lambda_jx^2}\ge\sum b_j
 -x^2\sum_{b_j>0}b_j\lambda_j\), (FA5), and the explicit \(M\).
The checker verifies the resulting lower pole coefficient exceeds
\(47/25\), and \(|M|/\sqrt{268/1000}<14/25\).

Since \(\sum_i\alpha_i=0\), some \(\alpha_i\le-A/\sqrt6\).
Writing \(\epsilon_i=g_i/\sqrt{J_+}=z_i-\alpha_i\phi\), and using
\(z_i\ge-1/3\), yields

\[
 \frac{G^2}{A^2}\ge
 \int J_+(t)\left(\phi(t)/\sqrt6-1/(3A)\right)_+^2d\nu(t).
 \tag{FA10}
\]

Set \(A_0=27/1000\),
\(B_0=1/3+A_0C_*/\sqrt6\), and
\(X_0=A_0c_*/(\sqrt6B_0)<1/10\).
The right side of (FA10) increases with \(A\).
If \(A\ge A_0\), integrate (FA8) on \(0<x<X_0\), with
\(d\nu(x^2)=2xe^{-x^2}dx\), to get

\[
 \frac{G^2}{A^2}\ge
 2j_*\int_0^{X_0}(c_*/\sqrt6-B_0x/A_0)^2dx
 =\frac{j_*c_*^3 A_0}{9\sqrt6B_0}
 =:\Phi_0>0.0025648>H^2.
 \tag{FA11}
\]

This contradicts (FA7). Thus every remaining profile satisfies
\(A<A_0\). This conclusion imposes no bound on its pointwise error
\(\epsilon\), or on its angular shape.

## FA4. Uniform angular coercivity on a mass cube

Let \(\mathcal C=\{p\in\Delta_3:1/6\le p_i\le1/2\}\).
The full directed box cover in FA7 proves

\[
             \Gamma(p,q)\ge\frac1{200}|p-q|^2,
                      \qquad p,q\in\mathcal C.
 \tag{FA12}
\]

The reduction to a matrix inequality is exact. Let
\(d_i=q_i-p_i\),
\(c_i=\cos[\pi(p_i+q_i)/2]
             \operatorname{sinc}[\pi(q_i-p_i)/2]\),
where \(\operatorname{sinc}(0)=1\). Thus
\(b_i-a_i=c_i d_i\) for the moment lengths of DG6.
Let \(C(p)\) have diagonal one and off-diagonal entry
\(-\cos(\pi p_k)\), and let \(D_c=\operatorname{diag}(c_i)\).
The arc-moment Gram matrix is
\(G(a)=\operatorname{diag}(a_i)C(p)\operatorname{diag}(a_i)\).
Its kernel contains \(\mathbf1\), hence

\[
 r^TG(a)r=d^TD_cC(p)D_c d.
\]

The analogous reverse formula uses \(C(q)\).
Since a maximum is at least the average of its two entries, it suffices
to prove the average of the two matrices
\(D_cC(p)D_c/[2(\sum r_i+M(r))]\) and its reverse is at least
\(I/200\) on \(\mathbf1^\perp\). Positive upper enclosures of each
denominator give Loewner lower bounds because the actual numerators
are positive semidefinite. The checker compresses onto the orthonormal
basis \((1,-1,0)/\sqrt2,(1,1,-2)/\sqrt6\), subtracts \(I/200\),
and certifies its leading entry and determinant positive on every box.
The sinc enclosure includes its full Taylor remainder. The diagonal
\(p=q\) causes no division by \(|p-q|\).

## FA5. Arbitrary profiles on an annulus

Use the fixed band \(I=[5/16,2]\) in the variable \(x=\sqrt t\).
FA7 certifies

\[
 \sup_{x\in I}|\phi(x^2)|<P_*:=97/20,\qquad
 \iint_{\sqrt t,\sqrt s\in I}V(t,s)(\phi(t)-\phi(s))^2d\nu(t)d\nu(s)
       >3/2.
 \tag{FA13}
\]

Define reference vectors \(p^0(t)=\mathbf1/3+\alpha\phi(t)\).
Since \(|\alpha_i|\le\sqrt{2/3}A\), they lie inside \(\mathcal C\)
with distance from each coordinate boundary at least
\(\zeta=1/6-\sqrt{2/3}A_0P_*>0\).
For **all actual** \(p(t),p(s)\in\Delta_3\), even outside this cube,

\[
 \Gamma(p(t),p(s))\ge
 \frac{A^2}{400}(\phi(t)-\phi(s))^2
 -\frac1{20}(|\epsilon(t)|^2+|\epsilon(s)|^2),
                       \qquad \sqrt t,\sqrt s\in I.
 \tag{FA14}
\]

If both actual vectors belong to the cube, apply (FA12),
\(|v+w|^2\ge|v|^2/2-|w|^2\), and
\(|\epsilon(t)-\epsilon(s)|^2\le2(|\epsilon(t)|^2+|\epsilon(s)|^2)\).
The error coefficient \(1/20\) is larger than \(2/200\).
Otherwise at least one error norm exceeds \(\zeta\), while the
reference difference has norm at most \(2P_*A_0\).
The checker verifies
\(\zeta^2/20-(2P_*A_0)^2/400>0\).
The right side of (FA14) is then nonpositive, and \(\Gamma\ge0\)
proves the inequality. This is the step extending the local mass-cube
coercivity to all amplitudes of an actual profile's remaining error.

From (FA9), or directly its weaker exponential bound,
\(w(a)+\theta(a)\ge(3j_1/2)a e^{-a}/I_0(a)\), so

\[
 \int V(t,s)d\nu(s)\le\frac{2J_+(t)}{3j_1}
                         <\frac{200}{27}J_+(t).
 \tag{FA15}
\]

Integrate (FA14) on the band, discard the nonnegative loss elsewhere,
and use (FA13)--(FA15) and \(G^2=\int J_+|\epsilon|^2d\nu\):

\[
                 \mathcal L(p)\ge\frac3{800}A^2-\frac{20}{27}G^2.
 \tag{FA16}
\]

## FA6. Assembly

Combining (FA2), (FA6) and (FA16) gives

\[
 S_\rho(A)-S_\rho(T)
 \le-\frac{13}{20000}A^2+2RAG-\frac{89}{135}G^2\le0.
 \tag{FA17}
\]

The last quadratic form is strictly negative unless \(A=G=0\), since
\((13/20000)(89/135)-1/50000>0\).
Cases excluded before (FA7) already satisfy the comparison.
All integrals exist: the actual \(z_i\) are bounded, \(J_+\) is
positive and bounded, the conjugated kernel is Hilbert--Schmidt, and
\(\Gamma\) is dominated by the finite moment deficit as in DG6.
Thus (FA17) proves the claimed planar comparison, including profiles
with arbitrarily large pointwise residuals and disconnected angular cells.
FA1 completes the ambient-dimension and lower-correlation reductions.

## FA7. Directed computation and coverage

Run from the canonical repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_finite_amplitude.py --output .scratch/finite-amplitude-production.json
```

The checker rebuilds every table at 100-bit Arb precision and uses the
retained, unchanged `plur_cert.py`, `plur_cert_positive.py`, and
`check_positive_spectral_extension.py` for outward binary64 primitives,
exact one-dimensional radial reduction, interpolation curvature bounds,
and Hilbert--Schmidt tails. It changes the auxiliary trial only inside
its process. No exploratory cache or floating eigenvector is an input.

For the trial integrals it evaluates exactly the quantities
\(D\), \(I=\langle\sqrt{J_+}f,(\mathcal A-2)\sqrt{J_+}f\rangle\),
and \(B=\|(\mathcal A-2)\sqrt{J_+}f\|^2\), by the one-dimensional
formulas (DG44)--(DG46) with the coefficients (FA4). The analytic infinite
tails (DG46) remain valid because every expanded correlation is below
\(1/2\), and \(1-\kappa>0\). Then \(a_u=2+I/D\) and the raw residual
bound is \(B/D\).

There are 16384 exact dyadic radial panels on \([0,7]\), with midpoint
rule error \(h^3\sup|F''|/24\) on each panel. Interval jets enclose the
true first and second derivatives, rather than differentiating table
interpolants. If \(F\) has validated table nodes and curvature bound
\(M_2\), its true derivative is enclosed by the node secant plus/minus
\(M_2\max(|q-q_l|,|q-q_r|)\), by the mean value theorem. Its second
derivative is bounded by \(M_2\) on their full convex hull. The source's
increasing curvature bounds supply this hull bound even across cells.
Product, reciprocal, square root, and exponential jet rules then propagate
the enclosures through the analytic integrands. Positive lower clamps
use proved function bounds; they do not replace derivative enclosures.

Four correlation parent intervals each have 25 children. True parameter
derivative bounds over the entire parent widen each child's midpoint
integral by its full correlation radius. A second mean-value enclosure
for the small residual, first in the radial variable and then in the
correlation, reduces interval dependence. Thus the computation covers
all parameters, not just the 100 child midpoints. The exact Fraction
audit verifies both rational endpoints and each child-parent containment.

The trial interpolation table has 8192 graded bins up to \(q=4\).
Hilbert--Schmidt bounds use 20 expanded correlation intervals, 6400 radial
panels per coordinate, 64 upstream angular midpoint panels, 1024 q bins,
4096 angular-parameter bins up to 64, and the original analytic Gaussian
tails. The shared exponential table has 131072 bins up to 64. These
domains enclose every query on the new range; out-of-domain queries
raise an error. The exact rational audit checks every cap in (FA5).

For (FA12), simultaneously permute labels to sort the first profile.
Then \(p_1\in[1/3,1/2]\), \(p_2\in[1/4,5/12]\),
\(p_3=1-p_1-p_2\in[1/6,1/3]\); the second profile has arbitrary
\(q_1,q_2\in[1/6,1/2]\) and \(q_3=1-q_1-q_2\in[1/6,1/2]\).
Bisect the longest of these four rational coordinates until the two
matrix minors are positive. Only boxes proved infeasible by the third
coordinate or ordering constraints are discarded. The complete root
box is exhausted: 10431 attempted boxes, 4700 accepted and 516 discarded
as infeasible. The full binary-tree count is checked. Sinc uses its
Taylor polynomial through degree 14 and remainder at most
\((\pi/6)^{16}/17!\), valid since \(|p_i-q_i|\le1/3\).

For (FA13), 2048 exact dyadic intervals per coordinate cover
\([5/16,2]^2\). The positive integral is evaluated without a pole using
\(h(x)=xf(x^2)\):

\[
 D\,\mathcal D_I(\phi)=\frac{8\rho}{d^2}
 \iint_{I^2}e^{-(x^2+y^2+2\rho xy)/d}
                   (yh(x)-xh(y))^2dxdy.
\]

Every rectangular cell is enclosed over the whole correlation interval;
the lower sum and \(D\le271/1000\) give a normalized lower bound
greater than 1.6745. The same mesh bounds \(|\phi|\) by 4.8049,
using directed square-root normalization. FA8--FA11 and the final
quadratic determinant are separately enclosed with Arb. Negative
controls reject an excessive coercivity constant, an insufficient
midpoint remainder, and deliberately incomplete correlation covers.
These controls supplement the complete written reductions.

## Remaining all-correlation obligation

This is a finite-amplitude realization of the mode-isolation proposal,
on an explicit interval. It does not prove that the same one-mode
decomposition works at larger correlations. The
[mode-isolation note](THREE_CELL_RADIAL_ANGULAR_MODE_ISOLATION.md) proves
finite-dimensional scalar obstructions and finite angular domination
of quadratic directions on each compact correlation interval, using
the other agent's joint angular and near-one machinery. Passing from
those directional estimates to a uniform nonlinear amplitude certificate
for every remaining correlation is still open. In particular DG7 already
shows that the elementary penalty (FA2) alone fails for some arbitrarily
small bounded profiles at \(\rho=3/5\). Higher angular information is
necessary there. The present finite-amplitude result makes no claim of
overlap with the qualitative near-one interval.
