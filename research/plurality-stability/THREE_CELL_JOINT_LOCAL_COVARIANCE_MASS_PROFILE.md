# Joint local covariance and mass response of the Gaussian tripod

Fix \(0<r<1\) and work in the Gaussian plane. Write \(T\) for
the regular three-sector partition, \(m_*=(1/3,1/3,1/3)\), and
\(T_O=OT\) for \(O\in O(2)\). Let \(V_T(R,m)\) be the supremum
over ordinary partitions with masses \(m\) in a sufficiently
small closed Gaussian \(L^1\) neighborhood of the orbit of \(T\).
No global optimality assumption on \(T\) is needed.

Use \(C\) for the facet mass flux and \(Q_r\) for the full
simplex interface form, including constant Hermite interaction
as in [(LM1)–(LM4)](THREE_CELL_LOCAL_MASS_PROFILE.md#exact-expansion-over-measurable-competitors).
Write \(K_B=\partial_R T_R|_{R=rI}[B]\) for the covariance derivative
in a symmetric matrix direction \(B\). For a facet speed \(f\)
on \(T_O\), put
\[
 G_{r,O,B}(f)=\sum_i\int K_B1_{T_{O,i}}\,d\lambda_i^f,\qquad
 H_{r,O}(B)=\partial_R^2 S_R(T_O)\big|_{R=rI}[B,B].
\]
The orientation is included in the facet measures
\(\lambda_i^f\), \(C\), and \(Q_r\). The matrix \(M_T\) in (AM2)
is scalar in the plane, so it does not depend on \(O\).
[Covariance differentiation](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#covariance-derivatives-and-rotation-stationarity)

Then for any symmetric \(B\) and \(\sum_i h_i=0\),
\[
\begin{split}
 V_T(rI+tB,m_*+th)
   &=S_r(T)+t\operatorname{tr}(BM_T)+t^2\Psi_r(B,h)+o(t^2),\\
 \Psi_r(B,h)
   &=\max_{O\in O(2)}
       \left\{\frac12 H_{r,O}(B)
          +\sup_{Cf=h}\bigl[2G_{r,O,B}(f)-Q_r(f)\bigr]\right\}.
                                                               \tag{JM1}
\end{split}
\]
The limit holds sequentially for converging directions \(B,h\)
and for \(r\) varying on compact positive correlation bands.
The inner supremum is attained modulo rotations. The formula
allows arbitrary measurable competitors.

## Relative covariance estimates for measurable transfers

For clarity, consider sequences \(r_j\to r\), \(B_j\to B\),
\(h_j\to h\), and \(t_j\to0\) positive. All matrices in the
segments from \(r_j I\) to \(r_j I+t_j B_j\) have eigenvalues
in a compact subinterval of \((0,1)\). Pick \(b<1\) strictly
above those eigenvalues and \(b/(2(1+b))<\alpha<1/4\).
For each fixed order through three, a covariance derivative
of the Mehler kernel is the kernel times a polynomial in its
two arguments whose coefficients are uniformly bounded
on those segments. This follows by differentiating the joint
Gaussian density with respect to the cross-covariance.
Diagonalization and the product of the one-dimensional
Mehler estimates bound the kernel by
\(C\exp(b(|x|^2+|y|^2)/(2(1+b)))\). Thus the derivatives
are bounded by \(Cw_\alpha(x)w_\alpha(y)\), with
\(w_\alpha(x)=e^{\alpha|x|^2}(1+|x|^2)\), by the strict
exponent margin.

For a competitor \(A\) and an oriented reference \(T_O\),
write \(e_i=1_{A_i}-1_{T_{O,i}}\) and \(L=L_{r_j}(A;T_O)\)
for its scalar score regret. By (B3)–(B5) the sum of the
weighted transfer masses is \(O(\sqrt L)\), with a uniformly
vanishing factor on the exterior of a growing ball. The
proof uses the saturated pair margin of the scalar simplex and
does not require balance. Since \(r_j\) is bounded away from zero,
the factor \(r_j^{-1/2}\) in (B3) is harmless here.
[Weighted exterior estimate and pair margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

Writing \(R_j=r_jI+t_jB_j\), the exact expansion at \(T_O\)
and these kernel estimates give
\[
\begin{split}
 S_{R_j}(A)-S_{R_j}(T_O)
 &=\sum_i\langle e_i,T_{r_j}e_i\rangle-2L\\
 &\quad +2t_j\sum_i\langle e_i,K_{B_j}1_{T_{O,i}}\rangle
               +O(t_j L+t_j^2\sqrt L).
                                                               \tag{JM2}
\end{split}
\]
In this display \(K_{B_j}\) is evaluated at \(r_j I\).
The linear pairing in the second line is \(O(\sqrt L)\).
The constant Hermite term \(\sum_i(\int e_i)^2\) is part of the
interaction in the first line.

## Flux limit and recovery

Take ordinary \(A_j\) within \(o(t_j^2)\) of the orbit-local
supremum at \((R_j,m_*+t_jh_j)\). Translating \(T\) to the
prescribed masses gives the lower bound
\(S_{r_j}(T)+t_j\operatorname{tr}(B_j M_{T,r_j})-O(t_j^2)\).
The first variation of \(T\) at \(r_jI\) vanishes without a
mass restriction, and the translation mass derivative has rank two.
The scalar energy of \(A_j\) is consequently at least
\(S_{r_j}(T)-O(t_j)\) by operator-norm covariance continuity.

Repair \(A_j\) to balance by the fixed local transports of
(LM3), rotated with a closest simplex. Their mass derivative
remains invertible uniformly in a small orbit neighborhood.
The repair changes the scalar energy by \(O(t_j)\) and sends
the competitor into the (B11) comparison neighborhood when
the original orbit radius is chosen sufficiently small.
The repaired partition is at distance \(O(\sqrt{t_j})\) from
the orbit by (B11). The uniform flow distortion and the
\(O(t_j)\) change of a simplex under the small flow give the
same \(O(\sqrt{t_j})\) distance for \(A_j\).

Use the closest-moment orientation \(O_j\) from the scalar
comparison and pass to a subsequence with \(O_j\to O\).
The normalized moment defect then satisfies the rotation
gauge. In particular the regrets \(L_j=L_{r_j}(A_j;T_{O_j})\)
tend to zero. In fact \(L_j=O(t_j^2)\). Otherwise normalize
the transfers by \(\sqrt{L_j/r_j}\) on a subsequence with
\(L_j/t_j^2\to\infty\). In (B13)–(B14) the limiting mass flux
is zero since \(t_j h_j/\sqrt{L_j/r_j}\to0\).
All additional terms in (JM2) are \(o(L_j)\), while
the lower bound for the nearly optimizing sequence gives
\(S_{R_j}(A_j)-S_{R_j}(T_{O_j})\ge -O(t_j^2)\).
The normalized mass cost and interaction therefore violate
the strict gap (B12), exactly as in the \(L_j=O(t_j^2)\)
argument for (LM1).
[Normal transfer compactness and the interaction limit](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)

Now normalize by \(t_j\). As in (LM2), the projected one-sided
transfer densities have weak \(L^2\) limits and the normalized
signed transfer measures converge to \(\lambda^f\) for a facet
speed with \(Cf=h\). The normal cost and the full interaction
give the upper bound \(-Q_r(f)\) for the scalar energy
coefficient. The functions \(K_{B_j}1_{T_{O_j,i}}\) converge
locally uniformly to \(K_B1_{T_{O,i}}\) by the continuous
covariance kernel derivatives and the strong convergence of
the rotated simplex indicators. The exterior weighted bound
on the normalized transfers makes the mixed pairing converge
to \(G_{r,O,B}(f)\). Taylor expansion of \(S_{R_j}(T_{O_j})\)
then proves the claimed upper bound.

For the reverse bound fix an orientation \(O\) and a smooth
speed of flux \(h\) compactly supported within facet interiors.
Change the coefficients of two compact mass bumps by \(o(1)\)
to realize flux \(h_j\), use normal displacements \(t_jf\),
and repair the remaining \(O(t_j^2)\) mass errors exactly.
The unpriced first variation at the simplex is zero. The
normal cost and full interaction give \(-Q_r(f)\) and the
mixed pairing gives \(2G_{r,O,B}(f)\). Approximation by such
speeds preserves the flux and proves the bound for all
weighted \(L^2\) speeds as in (LM2).

The null space of \(Q_r\) on balanced speeds consists of
rotations. The functional \(G_{r,O,B}\) vanishes on a rotation
speed, because twice its pairing with that speed is the
rotation derivative of \(\operatorname{tr}(B M_{T_O})\),
which is constant in the plane. Therefore coercivity on
the balanced quotient makes the supremum finite and attained.
The kernel bounds and the uniform quotient gap give
continuity of the maximized quadratic form in all parameters,
after identifying the facet spaces by rotations. The remaining
maximum is over the compact orientation group. This proves
the sequential and compact-band assertions in (JM1).

## The orientation envelope and a strict covariance response

Write \(B=bI+B_0\) with \(b=\operatorname{tr}(B)/2\) and
\(\operatorname{tr}(B_0)=0\). Scalar perturbations have
\(G_{r,O,I}=0\), since the scalar score differences vanish
on the simplex facets for every correlation. There are
constants \(c_r,\kappa_r\in\mathbb R\) such that
\[
 \Psi_r(B,h)=\tfrac12 S_r''(T)b^2
          +c_r\|B_0\|_{\rm F}^2+\Lambda_r|h|^2
          +2|\kappa_r|\|B_0\|_{\rm F}|h| .                   \tag{JM3}
\]
Indeed in any fixed orientation (JM1) is quadratic in
\((B,h)\), by the strictly positive Hessian on the balanced
quotient. The two-dimensional traceless covariance space and
the mass plane carry equivalent irreducible representations
of the tripod symmetries. A rotation by \(2\pi/3\) acts on
the former by angle \(4\pi/3\) and on the latter by a cyclic
permutation. A transposition acts by a reflection on both.
Thus the pure covariance quadratic terms are isotropic on
the traceless space and the mixed term is a multiple of one
isometric intertwiner. The scalar covariance space cannot
mix with the two nontrivial representations. Rotating the
reference simplex sweeps the full circle of traceless
directions, and maximizing the mixed term gives (JM3).

At \(h=0\) the coefficient in (JM1) is strictly above the
fixed-tripod coefficient \(\tfrac12H_{r,O}(B)\) whenever
\(B_0\ne0\). To see this, use outward facet normal \(n\),
unit ray direction \(\tau\), and arclength \(s\). The score
contrast is odd under the pair reflection and its outward
normal derivative is \(-a_r(s)\). The mixed derivative
formula (CT6) gives on this ray
\[
 (K_B1_{T_{O,i}}-K_B1_{T_{O,j}})(s\tau)
       =\frac{n^{\mathsf T}B\tau}{r}
                \bigl[2a_r'(s)-s a_r(s)\bigr].
                                                               \tag{JM4}
\]
For a nonscalar \(B\), at least one ray has
\(n^{\mathsf T}B\tau\ne0\), since three distinct tangent
lines cannot all be eigenlines. The bracket is nonconstant:
\(a_r(s)\) tends to a positive limit and \(a_r'(s)\to0\) by
the ray integral in (LM2). Choose a smooth speed on just
this ray with zero total flux and a nonzero pairing with
(JM4). A small multiple with the right sign makes
\(2G_{r,O,B}(f)-Q_r(f)>0\). This proves strict response.
[Mixed derivative on affine speeds](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions)

### The mixed coefficient is nonzero

Orient the three facets cyclically. In the standard sector,
the fixed mass-flux minimization for \(Q_r\) is governed by
the one-ray form
\[
 q(\varphi)=\int a_r\varphi^2\,d\nu
             -\iint K_+(s,u)\varphi(s)\varphi(u)\,d\nu(s)d\nu(u).
\]
Let \(\varphi\) minimize this form under \(\int\varphi\,d\nu=1\).
Here \(K_+\) is the standard-sector kernel of (S23) and
\(K_+1=a_r/r\) by (S17). Testing the normalized constant
shows \(q(\varphi)<0\). Its Euler equation is therefore
\(a_r\varphi-K_+\varphi=-\mu\) with \(\mu=-q(\varphi)>0\).
[Standard sector and translation row identity](SIMPLEX_INTERFACE_SYMMETRY.md#the-standard-sector-and-an-explicit-full-gap)

The operator \(P=r a_r^{-1}K_+\) is a compact positive
semidefinite self-adjoint Markov operator in \(L^2(a_r\,d\nu)\).
The eigenvalue one belongs to the constant function.
All eigenvalues on its orthogonal complement are strictly
less than \(r\). Otherwise, taking an eigenfunction there
with eigenvalue at least \(r\) and subtracting its \(\nu\)-mean
would give a nonzero mean-zero speed with nonpositive
standard form, contrary to (S23). Compactness makes the
spectral bound on that complement strictly less than \(r\).
If \(\overline{\phantom{f}}\) denotes the \(a_r\,d\nu\)-mean,
the Euler equation consequently gives a convergent Neumann
series in this Hilbert space:
\[
 \varphi-\bar\varphi
   =-\mu\sum_{n=0}^{\infty}r^{-n}P^n
                       \bigl(a_r^{-1}-\overline{a_r^{-1}}\bigr).
                                                               \tag{JM6}
\]

The Markov operator preserves nonincreasing functions.
For \(s_2>s_1\), the ratio between its two row densities
is nondecreasing in the integration variable \(u\).
In fact the nonseparable factor of \(K_+(s,u)\) is
\(F(rs u/(1-r^2))\), with \(F(z)=2e^z+e^{-z/2}\).
Both \((\log F)'(z)>0\) and \((\log F)''(z)\ge0\) for
\(z\ge0\), which gives positive mixed log derivative in
\(s,u\). Row normalization preserves the likelihood-ratio
ordering. Integrating a nonincreasing function against
the larger row decreases its expectation by the covariance
inequality with the increasing likelihood ratio.
Since \(a_r\) increases, every partial sum on the right
of (JM6) is nondecreasing. Convergence gives a nondecreasing
representative of \(\varphi\).

Set \(g(s)=2a_r'(s)-s a_r(s)\). The ray integral in (LM2)
makes \(a_r\) positive, increasing and concave on
\([0,\infty)\), so \(g\) is strictly decreasing. Also
\(a_r(s)\le 4a_r(0)/3\) by (LM6). Integration by parts
on the Gaussian ray, whose density is denoted by \(\nu'(s)\),
gives
\[
 \int g\,d\nu=\int s a_r(s)\,d\nu-2a_r(0)\nu'(0)
              \le-\tfrac23 a_r(0)\nu'(0)<0 .
\]
The opposite monotonicity of \(g\) and \(\varphi\) now gives
\(\int g\varphi\,d\nu<0\) by the covariance inequality and
\(\int\varphi\,d\nu=1\).

For \(B_0\ne0\) let \(d_e=n_e^{\mathsf T}B_0\tau_e\)
on the cyclically oriented rays. These numbers sum to
zero and are not all zero, by the \(120^\circ\) spacing
of the rays. Take \(h=C((d_e\varphi)_e)\ne0\).
The mass-constrained minimizer for this \(h\) is
\((d_e\varphi)_e\) modulo rotations by the standard-sector
decomposition and the Euler equation. Its pairing with
the covariance forcing is
\(r^{-1}\sum_e d_e^2\int g\varphi\,d\nu\ne0\).
After writing any other speed of that flux as the minimizer
plus a balanced speed, the mixed quadratic term in (JM1)
is exactly twice this pairing. In fact
\(\sum_e d_e^2=3\|B_0\|_{\rm F}^2/4\) and
\(|h|^2=3\sum_e d_e^2\) by the cyclic incidence matrix.
Since \(\nu'(0)/\nu([0,\infty))=2\phi_1(0)\)
and \(a_r(0)=3r\phi_1(0)/(2\sqrt{1-r^2})\), the bound on
\(\int g\,d\nu\) gives
\[
                 |\kappa_r|\ge\frac{1}{2\pi\sqrt{1-r^2}} .
                                                               \tag{JM7}
\]
The orbit-local second-order
coefficient therefore has distinct quadratic responses
as the relative orientation of covariance and quotas varies.
If the tripod is the unique globally maximizing orbit at
\((rI,m_*)\), the same expansion holds for the global value
\(\mathcal V\). Compactness of maximizers with varying quotas
and covariance puts all nearby optimizers in the chosen
neighborhood. In this case there is no joint second Fréchet
derivative at the base point. Such a derivative would make
\(\Psi_r\) a quadratic form, whereas its mixed term is positive
at both \(h\) and \(-h\) for every nonzero \(B_0,h\).

## A second-order test at a tie

Suppose \(A\) ties the tripod at \(rI,m_*\) and has
\(M_A=M_T\) and \(\beta(A)=0\), so that the first-order
discriminator (VV12) gives zero in every direction.
On the regular interfaces of \(A\), write \(Q_{r,A}\),
\(C_A\), and \(G_{r,A,B}\) for the corresponding second
variation form, mass flux and mixed covariance pairing.
For smooth compact interface speeds \(f\) away from
junctions with \(C_A f=h\), actual deformations with exact
mass correction give
\[
\begin{split}
 \liminf_{t\downarrow0}
   \frac{\mathcal V(rI+tB,m_*+th)-V_T(rI+tB,m_*+th)}{t^2}
  \ \ge{}
    \tfrac12\partial_R^2 S_R(A)\big|_{R=rI}[B,B]
        +2G_{r,A,B}(f)-Q_{r,A}(f)-\Psi_r(B,h).
                                                               \tag{JM5}
\end{split}
\]
Here \(\mathcal V\) is the unrestricted optimal value.
This follows by comparing the single deformed \(A\) with
(JM1). The exact mass correction is second order and
contributes nothing at that order because the initial
prices vanish. Thus any positive right side forces
maximizers away from the simplex neighborhood for the
corresponding small joint perturbations. Formula (JM5)
also applies in purely scalar covariance directions and
interpolates between the scalar shape test and the
local quota comparison.
[First-order joint contact test](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#joint-covariance-and-volume-response)
[Second variation at an actual maximizer](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

The comparison (JM5) leaves the original balanced scalar
tie open. A full exclusion would require a sign constraint
on the distant partition or an obstruction to maximizing
branches in the corresponding joint directions.

## Anisotropy creates unequal prices at balance

Suppose \(T\) is the unique maximizing orbit at \((rI,m_*)\).
Fix \(B_0\ne0\) and a mass direction \(h\). For any sequence
of actual maximizers \(A_j\) at
\((rI+t_j B,m_*+t_j h)\) with \(t_j\downarrow0\), let
\(\beta_j\) denote their centered winning-score prices.
If \(h\ne0\), then
\[
   \frac{\beta_j}{t_j}\longrightarrow
        \Lambda_r h+|\kappa_r|\|B_0\|_{\rm F}\frac{h}{|h|}.
                                                               \tag{JM8}
\]
If \(h=0\), then
\[
   \frac{|\beta_j|}{t_j}\longrightarrow |\kappa_r|\|B_0\|_{\rm F}>0 .
                                                               \tag{JM9}
\]
In particular the global value as a function of masses is
nondifferentiable at balance for all sufficiently small
positive perturbations \(tB\) with this fixed nonscalar
direction. This is so even though the scalar base has just
one maximizing orbit.

To prove these claims, compactness and uniqueness at the base
put all maximizing sequences in the orbit neighborhood of
(JM1). Its proof then gives \(L_j=O(t_j^2)\) and
\(\sum_i\int w_\alpha|e_{j,i}|\,d\gamma=O(t_j)\) after rotating
the reference into the moment gauge. Use two compact test
fields with independent mass flux at the reference simplex,
rotated with it. The first-variation identity (VV3) applied
at \(R=rI+t_j B\) gives \(|\beta_j|=O(t_j)\). Indeed the
mass-flux matrix stays invertible, the unpriced energy
derivatives at the reference scalar simplex vanish, and
their perturbation is \(O(t_j)\) by the Gaussian kernel
derivative bounds on the two compact supports. The volume
formula (VV3) remains valid for anisotropic positive
covariance by the same joint density change of variables.
[Transport and centered price formula](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#ambient-transport-and-uniform-mass-correction)

After passing to a subsequence, the orientations tend to
\(O\), the scaled prices to \(p\), and the normalized transfers
to a facet flux \(f\). The proof of (JM1) forces \(O\) and
\(f\) to attain the two optima in that formula. In addition,
stationarity identifies the limiting price. On a compact
relative-interior facet segment between \(i,j\) the scores
at the reference scalar simplex have a strict gap over the
third label. The pair contrast of the perturbed partition
differs from the reference contrast, divided by \(t_j\),
by a function converging locally uniformly to
\[
 K_B(1_{T_{O,i}}-1_{T_{O,j}})+
                  T_r(\lambda_i^f-\lambda_j^f)-(p_i-p_j).
\]
Here \(T_r\) applied to a surface measure denotes integration
against the Mehler kernel. This convergence follows from
the normal-flux limit and its exterior weighted control.
Spatial differentiation of the unscaled scores converges
locally uniformly to that of the reference scores, by the
same smooth kernel bounds. Thus only the two labels can
win on a small tube about the segment, and their interface
is the root of a contrast with strictly negative outward
normal derivative. On the normal scale \(t_j\) the signed
displacement is the displayed limit divided by \(a_r\),
by the first-order expansion of that root. Gaussian
integration in the normal direction identifies this
displacement with \(f\) on the segment. Exhausting the
facet interiors gives the weak equation
\[
           Q_r(f,u)=G_{r,O,B}(u)-p\cdot Cu
                                                               \tag{JM10}
\]
for compact regular test speeds and hence for all speeds
in the facet Hilbert space.

Equation (JM10) makes \(2p\) the mass derivative of the
optimized fixed-orientation quadratic coefficient. This
also follows by differentiating in flux using the linear
solution operator of the balanced quotient problem. Thus
by (JM3) before orientation maximization, \(p\) equals
\(\Lambda_r h+\kappa_r J(O^{\mathsf T}B_0O)\), with the
isometric intertwiner \(J\) used there. If \(h\ne0\), the
maximizing orientation aligns the second term with \(h\).
If \(h=0\), that term always has magnitude
\(|\kappa_r|\|B_0\|_{\rm F}\). These conclusions hold on
every convergent subsequence, proving (JM8)–(JM9).
Finally, at balance, relabeling any maximizer permutes its
centered prices. A nonzero centered price has at least two
distinct permutations, and the exact first quota derivative
(VV11) is their support function together with any other
maximizing prices. It cannot be linear. This proves the
nondifferentiability assertion.
