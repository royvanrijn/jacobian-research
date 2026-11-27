# Angular compatibility shared across radii

A single circular slice participates in every radial pairing. Optimizing
each pairing separately discards this constraint. This note retains it in
two ways: an exact rank-one first-moment optimization at arbitrary masses,
and a global finite-angular-order envelope with an explicit uniform error.
Neither result closes the all-positive-correlation simplex inequality.

The rank-one bound detects a finite-amplitude obstruction invisible to
every separately optimized pair envelope. For slices that split the circle
equally between two labels, with all three label pairs present, pairwise
alignment is possible for each pair but impossible simultaneously. The
weighted loss has a closed formula, and is stable under integrated mass
errors. The higher-order envelope converges uniformly over **all** mass
profiles, rather than only quadratically near balanced masses. Its global
moment optimization still has to be bounded to finish the direct proof.

These are written analytic results with internal checking only. There is
no claim of independent replay, formal verification or external review.

## GAC1. The first-mode comparison before pairwise relaxation

Use \(d\nu(t)=e^{-t}dt\), \(r=\sqrt{2t}\),
\(d\sigma=d\theta/(2\pi)\), and the definitions of DG15/DG26 in
the [direct comparison](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md).
For an ordinary planar partition write \(p_i(t)\) for its circular
masses and

\[
 u_i(t)=\int_{A_i(t)}(\cos\theta,\sin\theta)d\sigma,
 \qquad a_i(t)=\frac{\sin(\pi p_i(t))}{\pi}.
\]

Then \(\sum_i u_i=0\), \(|u_i|\le a_i\), and

\[
 \Delta(t,s)=\sum_i a_i(t)a_i(s)-\sum_i u_i(t)\cdot u_i(s)\ge0.
 \tag{GAC1}
\]

The proof of DG6, before replacing the actual moment deficit by
\(\Gamma\), gives

\[
 S_\rho(A)-S_\rho(T)\le U_\rho(p)-\iint V_\rho(t,s)\Delta(t,s)d\nu(t)d\nu(s),
 \tag{GAC2}
\]

where \(U_\rho(p)=\sum_i q_\rho(p_i-1/3)\), and, with
\(d=1-\rho^2\), \(\kappa=\rho^2/d\), \(v=2\rho\sqrt{ts}/d\),

\[
 V_\rho(t,s)=d^{-1}v e^{-\kappa(t+s)-v}.
\]

This uses the inspected Corollary A.3 of
[Heilman, arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1)
and the nonnegative monotone remainder in the circular rearrangement
argument of DG6. The latter imports
[Friedberg--Luttinger, Theorem I, printed pp. 1--2, (2)--(5)](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf).
No arc assumption on the actual slices is imposed.

There is a global rank-one minorant, without a radial cutoff:

\[
 V_\rho(t,s)\ge c_\rho h_\rho(t)h_\rho(s),\qquad
 c_\rho=\frac{2\rho}{d^2},\qquad
 h_\rho(t)=\sqrt t\,e^{-\rho t/(1-\rho)}.
 \tag{GAC3}
\]

Indeed \(2\sqrt{ts}\le t+s\), and
\(\kappa+\rho/d=\rho/(1-\rho)\). The nonnegativity of
\(\Delta\), not an operator-order assertion about this minorant, is
what permits its use. The weight is integrable, with
\(\int h_\rho d\nu=\sqrt\pi(1-\rho)^{3/2}/2\).

## GAC2. Exact global rank-one moment loss

Let \(w\ge0\) be any integrable radial weight, and put
\(W=\int w\,d\nu\), \(A_i=\int w a_i\,d\nu\). Define

\[
 \mathfrak D_w(p)=\sum_i A_i^2-
 \sup_{F:\,\sigma(F_i(t))=p_i(t)}
            \sum_i\left|\int w(t)u_{F,i}(t)d\nu(t)\right|^2.
 \tag{GAC4}
\]

The supremum ranges over ordinary measurable circular partitions at all
radii with the prescribed masses. It uses the **same** slice in every
pairing. It is finite and nonnegative by (GAC1). It admits the following
two-dimensional compact optimization, at arbitrary measurable masses.

Let \(H=\mathbf1^\perp\subset\mathbb R^3\), with orthonormal basis
\((1,-1,0)/\sqrt2,(1,1,-2)/\sqrt6\). For a positive semidefinite
\(2\times2\) matrix \(B\) of trace one, choose vectors
\(v_1,v_2,v_3\in\mathbb R^2\) whose Gram matrix is the extension of
\(B\) to \(H\). Thus \(\sum_i v_i=0\) and \(\sum_i|v_i|^2=1\).
Set

\[
 s_p(B)=\min_{\lambda\in\mathbb R^2}
                    \sum_i a_i(p)|v_i-\lambda|.
\]

The value is independent of the orthogonal factor chosen for the vectors.
Then

\[
 \boxed{\quad
 \mathfrak D_w(p)=\sum_i A_i^2-
 \left[\max_{B\succeq0,\,\operatorname{tr}B=1}
                    \int w(t)s_{p(t)}(B)d\nu(t)\right]^2.
 \quad}                                                        \tag{GAC5}
\]

The matrix region is the closed disk, using
\(B=\tfrac12\left(\begin{smallmatrix}1+x&y\\y&1-x\end{smallmatrix}\right)\),
\(x^2+y^2\le1\). The inner minimization is the convex weighted
geometric-median problem for three points. It may be restricted to their
convex hull. Boundary masses and singular matrices are included; if all
caps vanish, the support value is zero.

Here is a proof, including realizability of the support. For fixed masses
let

\[
 \mathcal K_p=\{(u_1,u_2,u_3):\sum_i u_i=0,\ |u_i|\le a_i(p)\}.
\]

The support of this compact convex body is attained by actual circular
partitions. If all masses are positive, every extreme point of
\(\mathcal K_p\) has at least two saturated caps: otherwise keep the
single saturated vector fixed and perturb the other two by opposite
small vectors, or perturb any two when no cap is saturated. For two
saturated caps, say \(i,j\), closure and the third cap imply that their
smaller angular separation is at least \(\pi(p_i+p_j)\). This follows
from the cosine law
\((a_k^2-a_i^2-a_j^2)/(2a_i a_j)=-\cos(\pi p_k)\).
Two arcs of masses \(p_i,p_j\) centered in these directions are therefore
disjoint. Their complement has the prescribed third mass and exactly
the required third moment. If one mass is zero, the body is a disk of
opposite vectors and its extreme points are complementary arcs. If one
mass is one, every cap is zero. Thus all support maxima are physically
realizable. This argument needs only the extreme points, not an assertion
that the entire norm-cap body is an ordinary moment set.

For fixed vectors \(v_i\), its support equals
\(\min_\lambda\sum_i a_i|v_i-\lambda|\). The upper bound follows
by subtracting \(\lambda\cdot\sum_i u_i=0\). At a minimizing
\(\lambda\), the subgradient condition gives vectors
\(u_i=a_i(v_i-\lambda)/|v_i-\lambda|\) for nonzero differences,
and \(|u_i|\le a_i\) for zero differences, whose sum is zero.
They attain the upper bound. This also proves the formula at vanishing
caps and coincident points.

Dualize the Frobenius norm of the integrated moment matrix. The dual
matrix can be centered on \(H\), have unit norm, and be factored through
the Gram matrix \(B\) above. Since the radial slices can be chosen
independently for this linear support problem, its supremum is
\(\int w s_{p(t)}(B)d\nu\). The measurable step can be checked without
a selection theorem: two arcs and their complement are described by a
finite choice of the complementary label and orientation, one rotation,
and a separation in \([\pi(p_i+p_j),\pi]\). Their moments depend
continuously on these parameters and on \(p\). Enumerate rational
rotations and normalized separations and select the first whose support
value is within \(\epsilon\) of the maximum. The resulting endpoints
are measurable in \(t\), preserve the masses exactly, and lose at most
\(\epsilon W\). Let \(\epsilon\downarrow0\).
The support value is continuous in \(p,B\) by this compact arc
parametrization, and bounded uniformly. Dominated convergence therefore
gives the maximum on the compact matrix disk. Squaring the resulting
norm supremum proves (GAC5).

The unsquared objective on the matrix disk is concave, so this is a convex
optimization problem. In fact, with \(\mathcal G_a\) the positive
semidefinite Gram matrices on \(H\) whose extended diagonal is bounded
by \(a_i^2\),
\(s_p(B)=\max_{C\in\mathcal G_a}
 \operatorname{tr}[(C^{1/2}BC^{1/2})^{1/2}]\).
The function in this maximum is jointly concave: it equals
\(\tfrac12\inf_{X\succ0}[\operatorname{tr}(XC)
                         +\operatorname{tr}(X^{-1}B)]\),
an infimum of linear functions of \((C,B)\).
The identity follows by completing the matrix square for positive definite
matrices and then by approximation for singular matrices. Maximizing a
jointly concave function over the fixed convex region \(\mathcal G_a\)
preserves concavity in \(B\), by testing convex combinations of its
maximizers. Integration against \(w\ge0\) preserves it as well.

There is also a direct certificate requiring no optimization over actual
angular sets. For any measurable positive definite matrix field \(X(t)\)
on \(H\), extend it by zero on \(\mathbf1\) and put
\(\lambda_i(X)=2X_{ii}+X_{jk}\), with distinct \(i,j,k\). Define

\[
 C_X=\int w\sum_i\lambda_i(X(t))_+a_i(t)^2d\nu(t),\qquad
 D_X=\int wX(t)^{-1}d\nu(t).
\]

If these integrals are finite, then

\[
 \mathfrak D_w(p)\ge
              [\sum_i A_i^2-C_X\lambda_{\max}(D_X)]_+.
 \tag{GAC5a}
\]

Indeed the squared matrix norm inequality gives, for every unit dual
matrix with Gram \(B\),
\(2\langle V,u(t)\rangle
 \le\operatorname{tr}(X(t)G_u(t))+
                             \operatorname{tr}(X(t)^{-1}B)\).
For a closing triple the first trace is
\(\sum_i\lambda_i(X)|u_i|^2
 \le\sum_i\lambda_i(X)_+a_i^2\), as in DG32--DG34.
Integrating and maximizing over \(\operatorname{tr}B=1\) bounds the
norm supremum by \((C_X+\lambda_{\max}(D_X))/2\).
Replace \(X\) by \(zX\) and minimize over \(z>0\) to get its
square at most \(C_X\lambda_{\max}(D_X)\). The zero case follows by
letting \(z\) tend to the appropriate endpoint. Exact cap-Gram supports,
when available, can replace the positive-part upper support in \(C_X\).
Piecewise constant positive matrices therefore give finite explicit
certificates for a shared-radius loss, including at arbitrary amplitudes.

One such certificate has a closed formula involving only radial integrals.
Choose a measurable label \(k(t)\), and let \(E_e\) be the three
sets where its complementary label pair is \(e\). Write
\(m_e=\int_{E_e}w\,d\nu\), \(W=\sum m_e\), and define

\[
 H_w=\int w\sum_i a_i^2d\nu,\quad
 E_w=\frac32\int w a_{k(t)}^2d\nu,\quad
 \Omega=\sum_e m_e^2-\sum_{e<f}m_em_f,\quad
 d_\pm=(W\pm\sqrt\Omega)/2.
\]

Then

\[
 \mathfrak D_w(p)\ge
 \left[\sum_i A_i^2-
       \inf_{q\ge1}(H_w-E_w+qE_w)(d_++d_-/q)\right]_+.
 \tag{GAC5b}
\]

For the proof, let \(P_e\) project onto
\((e_i-e_j)/\sqrt2\) in \(H\), and choose
\(X(t)=P_e+q(I_H-P_e)\) on \(E_e\). Its three trace coefficients
are \(1,1,(3q-1)/2\), the last at label \(k\); all are positive.
Consequently \(C_X=H_w+(q-1)E_w\), and
\(D_X=(1-1/q)\sum m_eP_e+(W/q)I_H\).
The two eigenvalues of \(\sum m_eP_e\) are \(d_+,d_-\): the
three doubled projector angles are separated by \(2\pi/3\), whose
cosine is \(-1/2\). Thus \(\lambda_{\max}(D_X)=d_++d_-/q\),
and (GAC5a) proves the claim.

Choosing \(k(t)\) to minimize its cap gives
\(0\le E_w\le H_w/2\). With \(H_0=H_w-E_w\), the infimum is

\[
 \begin{cases}
 (\sqrt{H_0d_+}+\sqrt{E_wd_-})^2,
                           &H_0d_-\ge E_wd_+,\\
 H_wW,                     &H_0d_-<E_wd_+ .
 \end{cases}
 \tag{GAC5c}
\]

This is the minimum of \(H_0d_++E_wd_-+E_wd_+q+H_0d_-/q\)
on \(q\ge1\), including its endpoint limits when a factor is zero.
On the half-mass family below, equal weighted type masses make (GAC5b)
sharp. In general
it is an explicit common-slice correction at every amplitude, without
solving the disk optimization. Its positivity or size relative to the
radial scalar form must still be checked on the profile being bounded.

The [two-slice moment analysis](THREE_CELL_EXACT_ANGULAR_MOMENT_ENVELOPE.md)
also identifies the norm-cap supports and an interior three-radius
transport obstruction. Formula (GAC5) keeps all radii in one optimization
and includes boundary masses, rather than optimizing their pairings
independently.

For any nonnegative valid pair penalty \(\gamma(p,q)\le\Delta\),
including DG26, one has

\[
 \mathfrak D_w(p)\ge\iint w(t)w(s)\gamma(p(t),p(s))d\nu(t)d\nu(s).
 \tag{GAC6}
\]

This follows by integrating the actual deficit and then taking the infimum
over the common slices. Combining (GAC2)--(GAC6), with \(w=h_\rho\), gives
the globally valid improvement

\[
\begin{split}
 S_\rho(A)-S_\rho(T)\le U_\rho(p)
 &-\iint(V_\rho-c_\rho h_\rho\otimes h_\rho)(t,s)
                         \gamma(p(t),p(s))d\nu(t)d\nu(s)\\
 &-c_\rho\mathfrak D_{h_\rho}(p).
\end{split}
 \tag{GAC7}
\]

It is at least as strong as subtracting \(\iint V_\rho\gamma\),
at every feasible amplitude. For a localized version, replace
\(h_\rho\) by its restriction to any measurable radial set in
(GAC3)--(GAC7). There is no stationarity or angular-shape hypothesis.

## GAC3. A sharp loss invisible to every pair envelope

Suppose each circular mass vector is one of
\((1/2,1/2,0),(0,1/2,1/2),(1/2,0,1/2)\).
Let \(E_{12},E_{23},E_{31}\) denote these three types, and put
\(m_e=\int_{E_e}w\,d\nu\), \(M=\sum_e m_e\),
\(m_{\max}=\max_e m_e\). Then

\[
 \boxed{\quad
 \pi^2\mathfrak D_w(p)
    =M^2-2\sum_e m_e^2+(2m_{\max}-M)_+^2
    \ge3\min_e m_e^2.
 \quad}                                                       \tag{GAC8}
\]

To prove it, orient the three type vectors cyclically as
\(b_{12}=e_1-e_2,b_{23}=e_2-e_3,b_{31}=e_3-e_1\).
Every slice moment matrix has form \(b_e\otimes z(t)\),
\(|z(t)|\le1/\pi\). Its integrated type vectors \(Z_e\) satisfy
\(|Z_e|\le m_e/\pi\), and its squared norm is
\(3\sum_e|Z_e|^2-|\sum_e Z_e|^2\).
This convex quadratic has a maximizer with each norm saturated: with
the other two vectors fixed, maximize over the closed disk and choose
a boundary point. For fixed lengths the smallest resultant in a plane
is \((2m_{\max}-M)_+/\pi\). The maximizing value is consequently
\([3\sum m_e^2-(2m_{\max}-M)_+^2]/\pi^2\).
It is attained by constant halfplane slices of each type with the
corresponding moment directions. The independent cap norm is
\((M^2+\sum m_e^2)/\pi^2\), giving equality in (GAC8).

For the lower bound, order \(m_1\ge m_2\ge m_3=s\).
If \(m_1\ge m_2+s\), the expression is \(4m_2s\ge4s^2\).
Otherwise write \(m_2=s+x,m_1=s+y\), with
\(0\le y-x\le s\). The expression becomes
\(3s^2+2s(x+y)-(x-y)^2\ge3s^2\).

In contrast, every pair of these profiles attains the independently
aligned arc envelope for any nonnegative symmetric decreasing circle
kernel. Identical types use identical complementary halfcircles.
Different types have exactly one common nonempty label; align that
halfcell and all other matched terms vanish. Thus every separately
optimized pair deficit is zero, including the higher polynomial pair
deficits of RMA6. Nevertheless the common-slice deficit in (GAC8) is
positive whenever all three weighted types occur. For
\(m_e=M/3\), it is \(M^2/(3\pi^2)\), one quarter of the
independent integrated cap envelope. This is a limitation of pairwise
optimization, not a counterexample to simplex optimality.

The bound is robust at arbitrary amplitudes. For any mass profile \(p\)
and any three-type reference \(p^0\), let
\(\epsilon=\int w|p-p^0|_1d\nu\), and compute \(m_e\) from the
reference. Then

\[
 \mathfrak D_w(p)\ge
 \left[\frac{M^2-2\sum m_e^2+(2m_{\max}-M)_+^2}{\pi^2}
                         -\frac{6W}{\pi}\epsilon\right]_+ .
 \tag{GAC9}
\]

For a proof, transfer excess angular mass directly to deficient labels
at each radius. This changes the total indicator \(L^1\) norm by exactly
\(|p-p^0|_1\), hence the Frobenius norm of the moment change by at most
that number. Measurable transfers follow by angular cumulative
quantiles, also when labels vanish. The cap vector has the same
Lipschitz bound, since \(|(\sin\pi p)'/\pi|\le1\).
The per-circle cap norm is at most \(3/(2\pi)\), by JA7, and the
actual moment norm is no larger. Each of the four integrated cap or
moment norms is therefore at most \(3W/(2\pi)\).
The difference of the two squared-cap terms and of the two squared
actual-moment terms is at most \(6W\epsilon/\pi\).
Apply (GAC8) to the transferred reference slices and then take the
infimum over the original actual slices. Nonnegativity gives (GAC9).

For a balanced three-type profile, each \(\nu(E_e)=1/3\).
Since \(h_\rho>0\) almost everywhere, all its weighted type masses
are positive. A uniform explicit bound is available: take
\(a=\log(12/11)\), \(b=\log12\), and
\(h_* =\min(h_\rho(a),h_\rho(b))>0\).
The radial measure outside \([a,b]\) is \(1/6\), and
\(\log h_\rho\) is concave, so \(m_e\ge h_*/6\).
Equations (GAC3) and (GAC8) give a strictly positive actual angular
loss at least \(\rho h_*^2/[6\pi^2(1-\rho^2)^2]\), independently
of the radial arrangement. No comparison with the entire scalar
contribution is inferred from positivity alone.

## GAC4. Global polynomial envelopes at all amplitudes

The pair hierarchy RMA6 optimizes a fresh pair for each \((t,s)\).
Instead define the angular polynomial

\[
 P_{n,v}(z)=e^{-v}\sum_{k=0}^{n}\frac{v^k(1+z)^k}{k!},\qquad n\ge0.
 \tag{GAC10}
\]

Both this polynomial and the remainder \(e^{vz}-P_{n,v}(z)\) are
nonnegative and nondecreasing on \([-1,1]\). Its angular Fourier
degree is at most \(n\), and its Fourier coefficients are nonnegative.
Indeed
\((1+\cos\theta)^k=2^{-k}\binom{2k}{k}
 +2^{1-k}\sum_{m=1}^k\binom{2k}{k-m}\cos(m\theta)\),
by the binomial expansion in \(e^{i\theta/2}\).
Here \(Q_\rho(t,s)=d^{-1}e^{-\kappa(t+s)}\).
Let \(\mathcal C_n(p)\) be the integral of the three independently
aligned-arc pairings for the kernel \(Q_\rho P_{n,v}(\cos\theta)\),
and let

\[
 \mathcal M_n(p)=\sup_{F,G:\,\sigma(F_i(t))=\sigma(G_i(t))=p_i(t)}
       \sum_i\iint Q_\rho P_{n,v}(\cos(\theta-\varphi))
       1_{F_i(t)}(\theta)1_{G_i(s)}(\varphi)
                    d\sigma(\theta)d\sigma(\varphi)d\nu(t)d\nu(s).
 \tag{GAC11}
\]

The same two radial fields \(F,G\) are used throughout the integral.
Both are ordinary; allowing soft partitions has the same supremum.
Circular rearrangement implies \(\mathcal M_n\le\mathcal C_n\).
Writing \(\mathcal C_\infty\) for the full exponential arc bound,
one obtains the global upper bound

\[
 S_\rho(A)-S_\rho(T)\le
 \mathcal E_n(p):=\mathcal C_\infty(p)-\mathcal C_n(p)
                         +\mathcal M_n(p)-S_\rho(T).
 \tag{GAC12}
\]

Only the remainder is rearranged independently. The polynomial part
keeps all angular compatibility across radii.

There is an explicit uniform approximation theorem. Let
\(\mathcal S_\rho(p)\) be the supremum of the actual full Gaussian
stability over ordinary partitions with circular masses \(p\), without
assuming balance. For every \(0<b<1\), put
\(R_b=(1+3b)/(4b)>1\). Then, simultaneously for all measurable
\(p:(0,\infty)\to\Delta_3\) and all \(0<\rho\le b\),

\[
 0\le\mathcal E_n(p)-[\mathcal S_\rho(p)-S_\rho(T)]
          \le\frac4{(1-b)^2}R_b^{-(n+1)}.
 \tag{GAC13}
\]

This is an all-amplitude **error bound**, not a sign theorem for
\(\mathcal E_n\). In particular its upper endpoint can still be
positive, and its global optimization is a genuine remaining problem.

To prove (GAC13), the full positive-correlation Gaussian operator is
positive semidefinite. Its bilinear maximum over two ordinary fields
with the same profile equals \(\mathcal S_\rho(p)\): use
\(2\langle F,T_\rho G\rangle\le S_\rho(F)+S_\rho(G)\),
and test \(F=G\). As \(0\le P_{n,v}\le e^{vz}\), the polynomial
bilinear maximum is no larger than this full maximum.
Conversely the remainder rearrangement in (GAC12) bounds every actual
competitor, so its right side is at least
\(\mathcal S_\rho(p)-S_\rho(T)\).

For the quantitative upper error, the aligned remainder sum is at most
its maximum kernel value times \(\sum_i p_i(t)p_i(s)\le1\).
At \(z=1\), the exponential tail satisfies

\[
 e^{-v}\sum_{k>n}\frac{(2v)^k}{k!}
       \le R_b^{-(n+1)}e^{(2R_b-1)v}.
\]

Put \(c=2R_b-1=(1+b)/(2b)\), so \(c\rho\le(1+b)/2<1\).
Using \(2\sqrt{ts}\le t+s\) again,

\[
 \iint Q_\rho(t,s)e^{cv}d\nu(t)d\nu(s)
 \le\frac{d}{(1-c\rho)^2}\le\frac4{(1-b)^2}.
\]

This proves the stated uniform error, including the whole radial axis
and all angular shapes. No Taylor interchange is used outside a
nonnegative dominated sum.

For completeness these optimization values are well defined. The
prescribed-profile soft feasible sets are weak-star compact; the
conditional mass constraints are closed when tested against arbitrary
bounded radial functions. Both the full Gaussian kernel and the
polynomial kernel are Hilbert--Schmidt at fixed \(\rho<1\): the
latter is pointwise bounded by the nonnegative Mehler kernel, whose
squared Hilbert--Schmidt norm in dimension two is
\((1-\rho^2)^{-2}\). Compactness gives joint weak continuity of
the bilinear form on bounded sets. Fixing one argument, choose an
extreme point of its linear maximizing face and then do the same for
the other. Such extreme points are ordinary. Indeed a fractional
pair of labels on a positive-measure set admits a bounded opposite
transfer with zero angular integral at each radius, by splitting
each atomless angular section into equal halves using its cumulative
measure. This contradicts extremality. The same argument identifies
the full soft self supremum with the ordinary one, or one may use
bilinear maximization and positivity as above.

## GAC5. Finite angular complexity and the unresolved sign

The supremum in (GAC11) can be approached by ordinary pairs whose
slice at almost every radius has at most \(6(n+1)\) switching points.
Add \(\epsilon\sum_i i\int\cos((n+1)\theta)1_{F_i}\) and the same
term for \(G\), integrating over the radial measure. At a maximizing
pair, each slice separately solves its prescribed-mass linear
assignment. The scores before this perturbation are trigonometric
polynomials of degree at most \(n\). Mass prices exist for the active
labels, as in JA2: use a supergradient of the finite concave assignment
value at the relative interior of its mass simplex. Every active pair
contrast after the perturbation has a nonzero harmonic of degree
\(n+1\), hence at most \(2(n+1)\) zeros. Three pairs give the stated
count. The linear perturbation is uniformly bounded by a constant
times \(\epsilon\); letting \(\epsilon\downarrow0\) proves the
approximation assertion. No compactness of the radial endpoint
functions or attainment of this complexity bound in the unperturbed
problem is assumed.

Thus the global comparison has an angular hierarchy with finite
slice complexity and a uniform all-amplitude error, retaining
compatibility lost by pairwise envelopes. The radial optimization
and the sign of (GAC12) remain unresolved. In the finite adverse-mode
decomposition of RMA4, (GAC7) or a higher-order global envelope can
replace the pairwise loss before comparing it with the positive
finite-mode contribution. Positivity of an angular correction,
quadratic directional domination, and uniform approximation error
do not by themselves establish that comparison. The direct global
inequality for all \(0<\rho<1\) remains open.
