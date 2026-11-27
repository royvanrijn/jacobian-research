# Finite-amplitude domination of the adverse radial mode

This note proves a direct comparison for **every** ordinary balanced Gaussian
three-cell partition in every dimension at
\(0<\rho\le447/1000=0.447\). Sections FA1–FA7 give
\([11/25,221/500]\) and the [projected extension](#fa8-projected-extension-to-044295)
covers the next slab. The full-cap and radius-weighted extensions in FA9–FA10
cover the remaining two slabs. The proof isolates one adverse radial direction, uses the
simplex constraints to limit its amplitude, and pays for its positive
contribution with unavoidable angular loss. There is no assumption of
closeness to a simplex, stationarity, or continuation of maximizers.

The numerical inputs are a full-domain directed certificate, with analytic
tails, from [the checker](check_three_cell_finite_amplitude.py). Its
[verification record](three_cell_finite_amplitude_verification.json) retains
the entire cover and source pins. The result has internal checking only;
independent replay, formal verification and external review are unperformed.
FA11 gives a relative score margin. FA12 combines the same argument with the
[separate scalar comparison](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md)
to reach 0.535. FA13 carries the margin past the sharp scalar cutoff 0.58. The all-positive-correlation
target remains open.

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

## FA8. Projected extension to 0.44295

The direct comparison holds on the adjoining slab
\(221/500\le\rho\le8859/20000\). The
[directed checker](check_three_cell_finite_amplitude_projected_extension.py)
and [retained cover](three_cell_finite_amplitude_projected_extension_verification.json)
give the following replacements for (FA5):

\[
 \frac{268}{1000}\le D\le\frac{272}{1000},\qquad
 \frac{399}{200}\le a_u\le\frac{1003}{500},\qquad
 \|(\mathcal A-2)u\|^2\le\frac1{15000},\qquad
 \|\mathcal A\|_{\rm HS}^2\le\frac{22}{5}.
 \tag{FA18}
\]

The trial (FA4) and its exact centering are unchanged. Since
\(22/5-(399/200)^2<(13/20)^2\), the proof of (FA6) now gives
\[
 \sum_i q(z_i)\le \frac3{500}A^2+2\mathsf r AG-\frac{27}{20}G^2,
 \qquad \mathsf r^2=1/15000.
 \tag{FA19}
\]
When the sum is nonnegative and \(A>0\), this bounds \(G^2/A^2\)
by \(H^2<0.005328\), using (FA7) with these constants.
If \(A=0\), nonnegativity instead forces \(G=0\).

Use the full simplex circumradius in the
[amplitude obstruction (AP3)](THREE_CELL_FINITE_AMPLITUDE_PROJECTION.md#a-stronger-analytic-version-of-the-one-mode-calculation).
The current cover verifies (FA8) with \(c_*=15/8\),
\(C_*=14/25\), \(j_*=107/1000\), and origin cutoff \(1/10\).
Indeed \(|p-\mathbf 1/3|\le\sqrt{2/3}\). For \(A\ge A_0=1/100\)
the residual divided by \(A\) is therefore at least
\((c_*/x-C_*-\sqrt{2/3}/A_0)_+\) on the origin band.
Integrating its square against
\(J_+(x^2)e^{-x^2}2x\,dx\ge2j_*x^2\,dx\) over the positive part gives
\[
 \frac{G^2}{A^2}
 \ge\frac{2j_*c_*^3 A_0}{3(\sqrt{2/3}+A_0C_*)}
 >0.005719 .
 \tag{FA20}
\]
The integration endpoint is below \(1/10\), as required.
This contradicts the bound from (FA19), so any profile with a
nonnegative scalar sum and nonzero trial component has \(A<A_0\).

The cover verifies the same cube inequality (FA12) and annular
inputs (FA13). For \(A<A_0\) the references \(p^0\) have the same
coordinate margin formula \(\zeta=1/6-\sqrt{2/3}A_0 P_*\).
Inside the cube use
\(|v+w|^2\ge(5/6)|v|^2-5|w|^2\) in place of the bound in (FA14).
Then the coefficient of \(A^2(\phi(t)-\phi(s))^2\) is \(1/240\)
and the coefficient of the two squared residuals stays \(1/20\).
If one actual vector leaves the cube, the bound still holds because
\(\zeta^2/20>(2P_*A_0)^2/240\) and \(\Gamma\ge0\).
The resulting integrated loss is at least \(A^2/160-(20/27)G^2\)
by (FA13) and (FA15). Subtract it in (FA2) to obtain
\[
 S_\rho(A)-S_\rho(T)
 \le-\frac{A^2}{4000}+2\mathsf r AG-\frac{329}{540}G^2.
 \tag{FA21}
\]
The determinant is \((1/4000)(329/540)-1/15000>0\).
Together with the strictly losing case of a negative scalar sum,
this proves the planar comparison, and maximizer dimension
reduction gives the result in every \(d\ge2\).

Every tie on the added slab has \(A=G=0\) in the radial
decomposition, hence equal masses at almost every radius. By (DG2)
the angular self-deficit in (DG3) is nonnegative, so equality
forces its angular difference term to vanish. The strictly positive
Fourier multipliers of the circular kernel then force the slices to
agree across almost every pair of radii. For this single angular
partition, equality of the angular self-comparison holds at every
positive kernel parameter by continuity and full support of the radial
pair measure. Differentiating at zero gives the sum of squared first
circular Fourier moments of three equal arcs. For each cell, the
first Fourier magnitude is maximized uniquely by an arc of its mass:
choose the phase, then integrate cosine on its upper level set.
Thus the three cells are equal arcs tiling the circle. Equality cases
in higher dimensions are actual maximizers and reduce by the same
dimension theorem.

Run the supplemental replay from the repository root:

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_finite_amplitude_projected_extension.py --output .scratch/projected-extension-replay.json

The driver uses the SHA-pinned base implementation of FA7. Its
two audited substitutions in the angular integration retarget the
correlation enclosure and widen the trial upper normalization
from \(271/1000\) to \(272/1000\). The spectral calculations reuse
FA7's interpolation, midpoint remainder, and analytic infinite
tails on 76 child panels in two parents. Four Hilbert–Schmidt slabs
reuse the original kernels and tail bounds. Both covers are audited
against the exact rational endpoints. The full angular cube and
annular comparisons are recomputed, and Arb checks the constants
in (FA19)–(FA21), including (FA20) and the outside-cube margin.
The replay certifies this slab only. Its lower endpoint is already
covered by FA1.

## FA9. Full-cap extension to 0.445

On \([8859/20000,89/200]\) the
[directed checker](check_three_cell_finite_amplitude_cap_extension.py)
and [cover](three_cell_finite_amplitude_cap_extension_verification.json)
give
\[
 \frac{268}{1000}\le D\le\frac{273}{1000},\quad
 \frac{4011}{2000}\le a_u\le\frac{503}{250},\quad
 \|\mathcal A\|_{\rm HS}^2\le\frac{111}{25},\quad
 \|P\mathcal A u\|^2\le\frac1{60000},\qquad P=I-u u^{\mathsf T}.
 \tag{FA22}
\]
On each spectral child the last bound subtracts
\((a_u-2)^2\) from \(\|(\mathcal A-2)u\|^2\), with an exact
rational lower bound on the subtracted square. The orthogonal
block has squared norm at most \(111/25-(4011/2000)^2<(13/20)^2\).
Thus
\[
 \sum_i q(z_i)\le \frac3{250}A^2+2r AG-\frac{27}{20}G^2,
 \qquad r^2=1/60000.                                               \tag{FA23}
\]
For \(A\ge A_0=9/500\), (FA20) with the certified inputs
\(c_*=187/100,\ C_*=14/25,\ j_*=107/1000\) gives
\(G^2/A^2>0.010157\). A nonnegative value in (FA23) instead
requires \(G^2/A^2<0.009478\).

For \(A<A_0\), apply the
[exact angular cap formula](THREE_CELL_EXACT_ANGULAR_MOMENT_ENVELOPE.md#an-interior-range-with-a-closed-cap-formula)
on the cube \([5/24,11/24]^3\), which lies in its domain
\([\arccos(4/5)/\pi,1/2]^3\). Write \(\Gamma_{\rm ex}\) for the
deficit of its optimal pairing from \(\sum_i a_i(p)a_i(q)\).
As in (AM5), set \(w_i=a_i(p)a_i(q)\). Then
\[
 \Gamma_{\rm ex}\ge
 \frac{\sum_{i<j}w_iw_j(1-\cos(\pi(p_k-q_k)))}{\sum_i w_i}
 \ge \frac7{250}|p-q|^2,
 \quad k\text{ the label other than }i,j.
\]
Use \(w_i\ge\sin^2(5\pi/24)/\pi^2\),
\(\sum_iw_i\le9/(4\pi^2)\), and
\(1-\cos(\pi d)\ge(\pi^2/2-\pi^4/384)d^2\) for \(|d|\le1/4\).
The cover verifies the resulting coefficient bound.
With \(|v+w|^2\ge(2/7)|v|^2-(2/5)|w|^2\), this yields the
annular penalty
\(\frac1{125} A^2(\phi(t)-\phi(s))^2
 -\frac{43}{500}(|\epsilon_t|^2+|\epsilon_s|^2)\).
If an endpoint leaves the cube, (FA13) gives a residual at least
\(\zeta=1/8-\sqrt{2/3}A_0(97/20)\), and the certificate checks
\(\frac{43}{500}\zeta^2>(2(97/20)A_0)^2/125\).
The penalty is then nonpositive and holds by nonnegativity of the
actual angular deficit. The annular Dirichlet integral exceeds
\(8/5\). Integrating and using (FA15) bounds the angular deficit
below by \((8/625)A^2-(172/135)G^2\).
Subtract this in the comparison (DG29), retaining the actual
angular deficit in place of its lower bound \(\Gamma\).
The determinant is positive since
\((8/625-3/250)(27/20-172/135)>1/60000\).
Together with the large-amplitude comparison this gives strict
loss for nonzero radial defects and the same tie characterization
as FA8. Dimension reduction gives the global comparison.

The directed replay uses 160 spectral children, eight
Hilbert–Schmidt slabs, the same analytic tails as FA7, and the
certified cap and annular constants:

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_finite_amplitude_cap_extension.py --output .scratch/three-cell-cap-extension-replay.json

## FA10. Radius-weighted angular errors through 0.447

On \([89/200,447/1000]\), the
[weighted checker](check_three_cell_finite_amplitude_weighted_extension.py)
and [cover](three_cell_finite_amplitude_weighted_extension_verification.json)
give
\[
 D\in[268/1000,274/1000],\qquad
 a_u\in[4023/2000,8071/4000],\quad
 \|P\mathcal A u\|^2\le\frac1{60000},\quad
 \|\mathcal A\|_{\rm HS}^2\le\frac{89}{20}.
\]
Since \(89/20-(4023/2000)^2<(16/25)^2\),
\(\sum_i q(z_i)\le(71/4000)A^2+2r AG-(34/25)G^2\) with
\(r^2=1/60000\). Use \(A_0=1/40\). The pole clipping (FA20)
with \(c_*=373/200,\ C_*=14/25,\ j_*=107/1000\) gives
\(G^2/A^2>0.013929\) for \(A\ge A_0\). The scalar bound
is negative since its nonnegative region would require
\(G^2/A^2<0.013756\).

On the annulus \(I=\{t:5/16\le\sqrt t\le2\}\) let \(P_*=97/20\),
\(c=7/625\) and
\[
 E(t)=\max\!\left\{\frac{7499}{100000},
 \frac{c A_0^2(P_*+|\phi(t)|)^2}
 {(1/6-\sqrt{2/3}A_0|\phi(t)|)^2}\right\}.
\]
The denominator before squaring is positive on \(I\).
For \(A<A_0\) the following holds on \(I^2\):
\[
 \Gamma(p(t),p(s))\ge c A^2(\phi(t)-\phi(s))^2
                   -E(t)|\epsilon_t|^2-E(s)|\epsilon_s|^2.       \tag{FA24}
\]
In the cube \([1/6,1/2]^3\) the directed checker verifies
\(\Gamma\ge (2/125)|p-q|^2\). Split the square with
\(|v+w|^2\ge(7/10)|v|^2-(7/3)|w|^2\), and bound
\(|\epsilon_t-\epsilon_s|^2\) by twice the sum of the squares.
The required endpoint penalty is \(28/375<7499/100000\).
If an endpoint \(t\) leaves the cube, its residual has size
at least \(1/6-\sqrt{2/3}A_0|\phi(t)|\) and its own
\(E(t)|\epsilon_t|^2\) cancels the full positive term of
(FA24). The argument is identical if \(s\) leaves the cube.

The checker verifies the **weighted row** inequality
\[
 \frac{E(t)}{J_+(t)}\int_I V(t,s)d\nu(s)\le\frac{11}{20}
 \quad(t\in I).
\]
It bounds the row on 512 radial boxes by directed Riemann sums
after removing the common \(x e^{-\kappa x^2}\) factor, where
\(x=\sqrt t\). The separate annular integral divided by the
trial normalization is at least \(163/100\). Integrating (FA24)
therefore gives \(\mathcal L(p)\ge (7/625)(163/100)A^2-(11/10)G^2\).
The final negative quadratic has positive determinant since
\[
 \big((7/625)(163/100)-71/4000\big)(34/25-11/10)>1/60000 .
\]
This proves the comparison and tie rigidity as in FA8. This
extension certifies the entire larger cube by adaptive splitting
and checks 160 spectral children and eight Hilbert–Schmidt slabs,
with the FA7 analytic tails:

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_finite_amplitude_weighted_extension.py --output .scratch/three-cell-weighted-extension-replay.json

## FA11. A global planar margin relative to transfer cost

On \(11/25\le\rho\le447/1000\), every balanced ordinary planar
partition \(B\) admits an oriented simplex \(T_B\) such that
\[
 S_\rho(T_B)-S_\rho(B)\ge c L_\rho(B;T_B)
 \tag{FA25}
\]
with a uniform \(c>0\). Use the actual angular deficit
\(\mathcal D(t,s)=\sum_i(a_i(t)a_i(s)-u_i(t)\cdot u_i(s))\)
from (DG27), where \(u_i\) are the first angular vector moments and
\(a_i=\sin(\pi p_i)/\pi\). The derivation of (DG29) retains
\(\mathcal D\) before applying the lower bound \(\Gamma\).
With \(E=S_\rho(T)-S_\rho(B)\), all four slab arguments give
\[
 E\ge10^{-5}\int J_+|p-\mathbf1/3|^2d\nu
           +\frac1{50}\iint V(t,s)\mathcal D(t,s)d\nu(t)d\nu(s).
 \tag{FA26}
\]
Indeed, for small \(A\) apply only \(49/50\) of the angular
bound. For each slab, write the scalar estimate as
\(e_0 A^2+2rAG-\delta G^2\) and the angular bound as
\(a A^2-bG^2\). The four pairs \((a,b)\) are
\((3/800,20/27)\), \((1/160,20/27)\), \((8/625,172/135)\)
and \(((7/625)(163/100),11/10)\).
In each case \(49a/50-e_0-10^{-5}\) and
\(\delta-49b/50-10^{-5}\) are positive with product exceeding
\(r^2\). For large \(A\) the verified clipping lower bounds for
\(z=G^2/A^2\) exceed, respectively,
\(.0025648,.0057,.01015,.0139\). At these four rational lower
bounds, \(\delta z-2r\sqrt z-e_0>10^{-5}(1+z)\), using
the respective upper estimates \(r<.0045,.0082,.0041,.0041\)
and \(\sqrt z<.051,.076,.101,.118\) (at the four endpoints).
The derivative of the left side is at least \(1\) thereafter.
This proves (FA26) in the large-amplitude regime while keeping
the entire angular deficit.

Let \(J=[1,2]\), \(w(t)=\sqrt t\,e^{-t/3-2\sqrt t}\).
For \(s\in J\), the displayed formula for \(V\) gives
\(V(t,s)\ge c_0 w(t)\) uniformly on these slabs, using
\(\kappa\le1/3\) and \(2\rho/(1-\rho^2)\le4/3\).
Integration and (FA15) give \(J_+(t)\ge c_1 w(t)\).
Put \(d_i(t)=a_i(t)^2-|u_i(t)|^2\). The identity
\[
 2\mathcal D(t,s)=\sum_i\big[d_i(t)+d_i(s)+|u_i(t)-u_i(s)|^2
                                      -(a_i(t)-a_i(s))^2\big]
\]
bounds the sum of its three nonnegative terms in the measure
\(w(t)d\nu(t)1_J(s)d\nu(s)\) by \(C E\).
For the subtracted square use the 1-Lipschitz bound for \(a_i\),
the first term of (FA26), and \(\inf_J w>0\). Averaging over \(J\)
gives a slice \(s\) with \(|p(s)-\mathbf1/3|^2+\sum_i d_i(s)\le CE\)
and \(\int w\sum_i|u_i(t)-u_i(s)|^2d\nu\le CE\).
The three vectors \(u_i(s)\) close and their lengths differ from
\(a_*=\sin(\pi/3)/\pi\) by \(O(\sqrt E)\).
For small \(E\), standardize the direction of \(u_1(s)\).
The cosine law and the nonzero sine at \(2\pi/3\) give a regular
ordered triple of unit vectors \(v_i\), with an allowed reflection,
such that \(\sum_i|u_i(s)-a_*v_i|^2\le CE\).
For \(E\) bounded away from zero, (FA25) follows from boundedness
of \(L_\rho\), so only small \(E\) is needed.

At radius \(t\), define the angular score loss
\(R(t)=\sum_i(a_*-v_i\cdot u_i(t))\ge0\), using the
winning-direction simplex \(T_B\) with sites \(v_i\).
Taylor's theorem and \(\sum_i p_i=1\) bound
\(\sum_i(a_*-a_i(t))\) by \(C|p(t)-\mathbf1/3|^2\).
On \(p_i(t)\in[1/6,1/2]\) for every label, use
\[
 a_i(t)-v_i\cdot u_i(t)
 =\frac{d_i(t)+|u_i(t)-a_i(t)v_i|^2}{2a_i(t)} .
\]
Off this cube \(R\) is bounded and the \(w\,d\nu\)-measure is
\(O(E)\) by (FA26). Thus \(\int w R\,d\nu\le CE\).

On a fixed physical ball of radius \(R_0\), reflection bounds the
OU score loss of each transfer by
\(2\Phi(\rho x_N/\sqrt{1-\rho^2})-1\le C x_N\), where
\(x_N\) is its distance to the pair plane. The first inequality
integrates the nonnegative reflected kernel difference on the
winning cell instead of the full positive half-line. Here
\(x_N=\sqrt{2t}(v_i-v_j)\cdot\omega/\sqrt3\) for the unit angle
\(\omega\). Therefore \(L_{\rm core}\le C_{R_0}\int w R\,d\nu\le C_{R_0}E\).
By (B3) and its global version in the
[continuation note](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound),
the centered interaction of the core defect is at most
\(C L_{\rm core}\) and that of the tail at most
\(\varepsilon(R_0)L_{\rm tail}\), with \(\varepsilon(R_0)\to0\).
Positivity of the centered operator bounds the full interaction
by twice the sum of these bounds. Fix \(R_0\) so that
\(2\varepsilon(R_0)\le1\). The exact identity (B1) then gives
\(E\ge L_{\rm tail}-C L_{\rm core}\), hence (FA25).

Finally, (CT5) in the
[correlation contact note](THREE_CELL_CORRELATION_CONTACT.md#a-holomorphic-estimate-relative-to-score-regret)
gives, with this same choice of \(T_B\) at \(r=447/1000\),
\[
 |\Delta_{r+h}(B;T_B)-\Delta_r(B;T_B)|
 \le C|h|L_r(B;T_B)
\]
for sufficiently small \(|h|\), by summing the Taylor series.
Here \(\Delta_\rho=S_\rho(T_B)-S_\rho(B)\).
Thus the margin gives an additional open right neighborhood
of \(447/1000\) with strict planar comparison away from rotated
simplices. Maximizer dimension reduction gives the global comparison
and rigidity there. This argument gives a positive, unevaluated length
for the extra neighborhood.

## FA12. Scalar extension and transfer-cost margin through 0.535

The separate [global scalar comparison](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md)
supplies the interval through \(53/100\). A
[supplemental directed cover](three_cell_global_scalar_narrow_extension_verification.json)
extends its criterion (GH6), (GH9)–(GH10) to \([53/100,107/200]\).
It uses 100-bit Arb, 80 rational slabs of width \(1/16000\),
16384 radial panels up to 12, \(m=16\) row terms and \(N=96\).
The centered finite norm plus the trace tail is below \(399/200\),
leaving a radial deficit at least \(1/200\) against the full row \(J\).
The planar comparison and equality case follow as in the scalar proof.
Replay with the [retained supplemental checker](check_three_cell_global_scalar_narrow_extension.py):

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_global_scalar_narrow_extension.py --start 8480 --stop 8560 --denominator 16000 --modes 96 --cap 399/200

The scalar comparison also admits the relative margin (FA25) on
\([11/25,107/200]\), including a positive extension beyond \(107/200\).

In fact (GH2) estimates the full aligned-arc envelope, so the
derivation of (DG29) can subtract \(\iint V\mathcal D\,d\nu d\nu\)
from (GH4). For each same-label angular pairing split
\(e^{a\cos\theta}\) into \(ae^{-a}\cos\theta\) and a remainder
whose derivative as a function of \(\cos\theta\) is nonnegative.
Circular rearrangement bounds the remainder by the aligned arcs
of the same masses, with a constant added to the kernel if needed.
The removed first-mode interaction is exactly the difference
\(\mathcal D\) from the aligned arcs, with the radial factor \(V\).
The directed negative radial margin in (GH6)–(GH10) therefore gives
\[
 E\ \ge\ \frac1{200}\int J_1(t)|p-\mathbf1/3|^2d\nu
                   +\iint V(t,s)\mathcal D(t,s)d\nu(t)d\nu(s).
 \tag{FA27}
\]
Here \(J_1\) is the first term in (GH5), and the scalar
subintervals control \(J_m\) and \(J\), respectively, both at
least \(J_1\). For \(w(t)=\sqrt t\,e^{-(2/3)t-3\sqrt t}\)
one has \(J_1(t)\ge c w(t)\) by its explicit positive first term.
For \(s\in[1,2]\) also \(V(t,s)\ge c w(t)\), since on the whole
band \(\kappa\le2/3\) and \((2\rho/(1-\rho^2))\sqrt2\le3\).
This gives precisely the two estimates used in the reference-slice
and core–tail proof of (FA25). Thus the margin holds uniformly through
\(107/200\). Apply (CT5) at this endpoint and reduce actual maximizers
to the plane as in FA11. The resulting global comparison includes
a right neighborhood of \(107/200\) with simplex equality rigidity.

### Scalar trace tails on arbitrary compact positive bands

The convergence restriction in (GH10) can be removed by using
a shorter part of the same positive angular measure. For
\(0<v<1/3\) set \(q=2v-3v^2/2\) and \(a=\cos^2(\pi v)\).
In the even folded difference measure of (GH5), the positive
part on \([0,v]\) has mass \(q\) and \(\cos^2(\pi u)\ge a\).
The rest of the positive part has squared cosine at least \(1/4\),
and the negative part has squared cosine at most \(1/4\).
Thus \(M_{2j}-P_{2j}\ge q(a^j-4^{-j})\). As in (GH9), this gives,
with \(x=\sqrt t\),
\[
 J(t)\ge\frac q4 e^{-\kappa x^2}
                     (e^{a\kappa x^2}-e^{\kappa x^2/4}).
\]
On any compact band \([\rho_-,b]\subset(0,1)\), choose
\(a>(2b-1)/b^2\) by shrinking \(v\).
For sufficiently large fixed \(x\), this lower bound is at least
\((q/8)e^{-(1-a)\kappa x^2}\). The same Gaussian integration as
in (GH10) now gives a geometric trace tail on the outer range
with ratio \(r=2\rho/(1+a\rho^2)<1\) uniformly. The bounded
inner range is controlled by the first row term and the
factorial tail. This makes the tail summation available on each
compact positive band. It leaves the sign of the centered
finite block as a separate comparison.

## FA13. Transfer-cost margin past the sharp scalar endpoint

The [sharp aligned-arc scalar comparison](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)
extends the explicit cutoff to \(29/50\) with a negative radial
margin at least \(1/1000\) against its full row weight \(J\) on
\([53/100,29/50]\). Its bound estimates the same full aligned-arc
envelope as in FA12, so the entire first angular deficit
\(\iint V\mathcal D\) is still available. The sharp row has a
positive first term proportional to \(\sqrt t e^{-\kappa t}\)
by (SP4), with positive coefficient uniformly on the compact
band. The proof of (FA27) gives the same two weighted estimates
with \(1/1000\) as radial coefficient and
\(w(t)=\sqrt t e^{-2t/3-3\sqrt t}\) on the full combined band
\([11/25,29/50]\). Indeed on this band
\(\kappa\le2/3\) and \((2\rho/(1-\rho^2))\sqrt2\le3\).
The reference-slice and core–tail argument then gives a uniform
relative margin (FA25) up to \(29/50\). By (CT5) a positive right
neighborhood of \(29/50\) has the global inequality and simplex
equality rigidity. The length of that neighborhood is unevaluated.

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
