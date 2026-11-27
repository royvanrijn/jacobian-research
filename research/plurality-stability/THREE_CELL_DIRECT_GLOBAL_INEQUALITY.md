# A direct three-cell comparison and the angular information it requires

The unrestricted inequality
\(\widetilde Q_\rho(e)\leq2L_\rho(A;T)\) remains open here.
The exact expansion (B1) in the [continuation note](POSITIVE_CORRELATION_CONTINUATION.md)
makes it equivalent to \(S_\rho(A)\leq S_\rho(T)\), including for
partitions far from the simplex. The argument below proves the comparison
at every positive correlation when the three labels occupy equal portions
of almost every circle. For general balanced partitions it isolates an
exact missing inequality and proves that discarding the angular
Cauchy--Schwarz remainder cannot establish it near correlation one.
The counterexample is to an upper-bound strategy, not to simplex optimality.

## Primary results and scope

The closest inspected primary text is Heilman,
[arXiv:2306.03312v3](https://arxiv.org/html/2306.03312v3), Theorem 1.5,
Sections 1.4, 5, 6 and 8. Its current theorem statement gives the
unrestricted positive comparison for \(0<\rho\leq0.1082\), not every
positive correlation. Section 1.4 identifies the competing radial-mean
and angular terms. The proof of Lemma 6.1 uses the angular
Cauchy--Schwarz step studied below.

The circular rearrangement input needed for the restricted comparison is
only this fixed-measure statement: if \(E\subset S^1\) has normalized
length \(p\), and \(I_p\) is an arc of length \(p\), then

\[
 \iint 1_E(\theta)1_E(\varphi)e^{\kappa\cos(\theta-\varphi)}
             d\sigma(\theta)d\sigma(\varphi)
 \leq
 \iint 1_{I_p}(\theta)1_{I_p}(\varphi)e^{\kappa\cos(\theta-\varphi)}
             d\sigma(\theta)d\sigma(\varphi),\qquad \kappa>0.
 \tag{DG0}
\]

The statement, hypotheses and its application were inspected in the
rearrangement step of Heilman's proof of Lemma 6.1. That paper attributes
the input to Baernstein--Taylor, Theorem 2; the 1976 original was not
inspected. We use (DG0) as the standard imported circular rearrangement
inequality. No numerical constant or variable-mass optimization from
Heilman's Corollary 4.2 is needed below.

For general ambient dimension, existence and reduction of an actual
maximizer to a plane follow from Heilman--Tarter, Lemma 2.3 and Theorem 1.9,
[Forum of Mathematics, Sigma 9 (2021), e65](https://doi.org/10.1017/fms.2021.56).
The theorem concerns maximizers, not a dimension reduction of an arbitrary
competitor. The polar identities below are planar; the equal-circle
comparison also has a direct cylindrical extension. These comparisons
are not a literature-exhaustiveness or priority claim.

## Exact separation of radial masses and angular changes

Fix \(0<\rho<1\), put \(s_\rho=\sqrt{1-\rho^2}\), and use
\(d\sigma=d\theta/(2\pi)\) and \(d\mu(r)=r e^{-r^2/2}dr\).
For an ordinary measurable three-cell partition of the plane write

\[
 a_i(r,\theta)=1_{A_i}(r\cos\theta,r\sin\theta),\quad
 p_i(r)=\int a_i(r,\theta)d\sigma(\theta),\quad
 b_i=p_i-\tfrac13,\quad g_i(r,\theta)=a_i(r,\theta)-p_i(r).
\]

Balance says \(\int b_i\,d\mu=0\); it does not say \(b_i=0\).
Define

\[
\begin{split}
 \kappa&=\frac{\rho r u}{s_\rho^2},\qquad
 I_0(\kappa)=\int e^{\kappa\cos\theta}d\sigma(\theta),\qquad
 k_\kappa(\theta)=\frac{e^{\kappa\cos\theta}}{I_0(\kappa)},\\
 (U_\kappa f)(\theta)&=\int k_\kappa(\theta-\varphi)f(\varphi)d\sigma(\varphi),\\
 P_\rho(r,u)&=s_\rho^{-2}
  \exp\!\left[-\frac{\rho^2(r^2+u^2)}{2s_\rho^2}\right]I_0(\kappa).
\end{split}
 \tag{DG1}
\]

Then \(P_\rho\) is the reversible Markov kernel for the radius under
the planar OU operator, relative to \(\mu\). Write
\(dW_\rho=P_\rho(r,u)d\mu(r)d\mu(u)\), a probability measure.
Let \(g^*\) denote the centered indicator vector of a fixed partition
of the circle into three equal arcs, and set

\[
\begin{split}
 H_\kappa(r,u)&=\sum_i\langle g_i(r),U_\kappa g_i(u)\rangle_{\sigma},
 &H_\kappa^*&=\sum_i\langle g_i^*,U_\kappa g_i^*\rangle_{\sigma},\\
 \mathcal R_\rho(A)&=\sum_i\langle b_i,P_\rho b_i\rangle_{\mu},\\
 \mathcal D_\rho(A)&=\int\left[H_\kappa^*
                -\frac{H_\kappa(r,r)+H_\kappa(u,u)}2\right]dW_\rho,\\
 \mathcal C_\rho(A)&=\frac12\int\sum_i
       \langle g_i(r)-g_i(u),U_\kappa[g_i(r)-g_i(u)]\rangle_{\sigma}
                                                            dW_\rho.
\end{split}
 \tag{DG2}
\]

In \(H_\kappa(r,r)\), the angular operator still has parameter
\(\kappa=\rho r u/s_\rho^2\), not \(\rho r^2/s_\rho^2\).
The diagonal notation refers only to its two angular arguments.

**Proposition DG1 (exact polar deficit).** For every ordinary balanced
three-cell partition of \(\mathbb R^2\),

\[
 S_\rho(T)-S_\rho(A)
    =\mathcal D_\rho(A)+\mathcal C_\rho(A)-\mathcal R_\rho(A).
 \tag{DG3}
\]

Both \(\mathcal C_\rho\) and \(\mathcal R_\rho\) are nonnegative.
Consequently the unrestricted direct comparison is exactly the remaining
inequality

\[
 \mathcal R_\rho(A)\leq\mathcal D_\rho(A)+\mathcal C_\rho(A).
 \tag{DG4}
\]

To prove the identity, write the planar Mehler density in polar
coordinates. The angular conditional density is \(k_\kappa\), so

\[
 S_\rho(A)=\int\left[\sum_i p_i(r)p_i(u)+H_\kappa(r,u)\right]dW_\rho
          =\tfrac13+\mathcal R_\rho(A)+\int H_\kappa(r,u)dW_\rho.
\]

The second equality uses balance and the Markov property. For the
simplex the corresponding expression is
\(1/3+\int H_\kappa^*dW_\rho\). Polarizing the self-adjoint
quadratic form of \(U_\kappa\) proves (DG3).
Its Fourier eigenvalues are \(I_m(\kappa)/I_0(\kappa)>0\), as follows
directly by expanding
\(e^{(\kappa/2)(z+z^{-1})}\):

\[
 I_m(\kappa)=\sum_{j\geq0}
   \frac{(\kappa/2)^{2j+|m|}}{j!(j+|m|)!}>0.
\]

This proves \(\mathcal C_\rho\geq0\); positivity of the OU quadratic
form proves \(\mathcal R_\rho\geq0\). All angular functions are
bounded and all operators are Markov contractions, so all terms are
finite. Fubini and the polarization identity apply without a limiting
regularity assumption. This proof needs no rearrangement theorem.

The qualitative self-comparison in Heilman's Lemma 6.1 also gives
\(\mathcal D_\rho\geq0\). That additional imported fact alone does
not compare its size with \(\mathcal R_\rho\).

## The comparison under equal circular masses

**Proposition DG2 (equal-circle comparison).** Suppose that
\(p_i(r)=1/3\) for every label and almost every radius. Then
\(S_\rho(A)\leq S_\rho(T)\), equivalently
\(\widetilde Q_\rho(e)\leq2L_\rho(A;T)\), for every \(0<\rho<1\).
The same holds in \(\mathbb R^d=V\oplus V^\perp\), \(\dim V=2\),
if the three labels occupy one third of almost every circle
\(\{(r\cos\theta,r\sin\theta,z):\theta\in S^1\}\).

Indeed (DG0), separately for each label, bounds
\(H_\kappa(r,r)\) by \(H_\kappa^*\), since all three arc lengths
are \(1/3\). Hence \(\mathcal D_\rho\geq0\), while
\(\mathcal R_\rho=0\); now apply (DG3).
For the higher-dimensional assertion, apply the same angular argument
to the two slices \((r,z),(u,w)\). The remaining Mehler factor in
\((z,w)\) is positive and integrates to one. No reduction of arbitrary
partitions or stationarity assumption is used.

In particular this settles partitions equivariant under a planar
rotation by \(2\pi/3\) that cyclically permutes the three labels and
fixes \(V^\perp\). It permits arbitrary measurable angular changes
and does not require invariance under reflection. It also settles every
equal-angle sector partition with a measurable radius-dependent rotation.
These are consequences of circular rearrangement, not a solution for
all balanced partitions. The three-seed cycles inside a four-cell
tetrahedron in the [swap note](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md)
have a different slice geometry and are not covered by this assertion.

## Why dropping the angular remainder fails

The angular Cauchy--Schwarz bound replaces
\(H_\kappa(r,u)\) by
\([H_\kappa(r,r)+H_\kappa(u,u)]/2\). Thus it gives the valid upper bound

\[
 S_\rho(A)\leq\mathcal B_\rho(A)
 :=\tfrac13+\mathcal R_\rho(A)
       +\frac12\int[H_\kappa(r,r)+H_\kappa(u,u)]dW_\rho,
 \qquad \mathcal B_\rho(A)-S_\rho(A)=\mathcal C_\rho(A).
 \tag{DG5}
\]

**Proposition DG3 (failure of the angular Cauchy relaxation).** There
is an exactly balanced family of ordinary planar partitions, with each
label a single arc on every circle, for which
\(\mathcal B_\rho(A)>S_\rho(T)\) for every sufficiently large fixed
\(\rho<1\) and all sufficiently small nonzero family parameters.
Thus the stronger assertion
\(\mathcal R_\rho\leq\mathcal D_\rho\) is false, even in this class.

Put

\[
 b(r)=e^{-r^2/2}-\tfrac12,\qquad c=(1,-1,0),\qquad
 p_i(r)=\tfrac13+t c_i b(r).
 \tag{DG6}
\]

Take consecutive arcs of these lengths, with endpoints
\(0,2\pi p_1(r),4\pi/3,2\pi\). For \(|t|<2/3\) they define an
ordinary partition. Since \(\int b\,d\mu=0\), every cell has mass
exactly \(1/3\). At \(t=0\) the partition is \(T\).

For an arc of normalized length \(p\), its centered self-energy is

\[
 F_\kappa(p)=2\int_0^p(p-v)k_\kappa(2\pi v)dv-p^2,
 \qquad F_\kappa''(p)=2[k_\kappa(2\pi p)-1].
 \tag{DG7}
\]

Therefore, at each fixed \(\rho<1\),

\[
 \mathcal B_\rho(A(t))-S_\rho(T)
   =t^2\|c\|^2 E_\rho+o(t^2),\qquad
 E_\rho=\langle b,P_\rho b\rangle
                  -\int[1-m_\rho(r)]b(r)^2d\mu(r),
 \tag{DG8}
\]

where

\[
 m_\rho(r)=\int P_\rho(r,u)k_\kappa(2\pi/3)d\mu(u).
\]

The first derivatives cancel because \(\sum_i c_i=0\).
Differentiation under the integral is valid: near \(p=1/3\),
\(k_\kappa(2\pi p)\) and its needed angular derivatives are bounded
uniformly in \(\kappa\geq0\), and \(b\) is bounded.

Gaussian integration gives

\[
 \|b\|_\mu^2=\tfrac1{12},\qquad
 \langle b,P_\rho b\rangle=\frac1{4-\rho^2}-\tfrac14,\qquad
 \langle b,(I-P_\rho)b\rangle=\frac{s_\rho^2}{3(4-\rho^2)}.
 \tag{DG9}
\]

For example the joint Laplace transform of the two correlated planar
Gaussians gives
\(\mathbb E e^{-(|X|^2+|Y|^2)/2}=1/(4-\rho^2)\).
With \(a=\rho r/s_\rho\), the remaining multiplier is

\[
 m_\rho(r)=M(a),\qquad
 M(a)=e^{-a^2/2}\int_0^\infty v e^{-v^2/2-av/2}dv.
 \tag{DG10}
\]

In particular \(0<M(a)\leq e^{-a^2/2}\). Substitution
\(r=s_\rho a/\rho\) and dominated convergence show that

\[
 \lim_{\rho\uparrow1}\frac1{s_\rho^2}\int m_\rho b^2d\mu
   =\tfrac14 C,\qquad
 C=\iint_{(0,\infty)^2}av\,e^{-(a^2+v^2+av)/2}da\,dv
   =\frac43-\frac{4\pi}{9\sqrt3}.
 \tag{DG11}
\]

The last integral can be evaluated by differentiating the elementary
quadrant Gaussian integral in its cross coefficient, or by polar
coordinates. For completeness the general formula is

\[
 J(A,c):=\iint uv\,e^{-(Au^2+v^2+2cuv)/2}du\,dv
 =\frac1{A-c^2}-\frac{c\arccos(c/\sqrt A)}{(A-c^2)^{3/2}},
 \quad A>c^2,\ c\geq0.
\]

Here both integrals are on the positive quadrant and \(C=J(1,1/2)\).
Combining (DG8)--(DG11) proves

\[
 \boxed{\displaystyle
 \lim_{\rho\uparrow1}\frac{E_\rho}{1-\rho}
    =\frac{4-2\pi/\sqrt3}{9}>0.}
 \tag{DG12}
\]

Positivity is the elementary inequality \(\pi<2\sqrt3\).
First fix any sufficiently large \(\rho<1\), and then take sufficiently
small \(t\ne0\) in (DG8). This proves Proposition DG3 with the stated
order of quantifiers. It asserts nothing adverse about the actual
stability of \(A(t)\), since (DG5) has discarded \(\mathcal C_\rho\).
Unlike the independent-arc relaxation in the separate
[global obstruction note](THREE_CELL_GLOBAL_OBSTRUCTIONS.md), the
partition here remains disjoint on every circle. The loss occurs in
replacing interactions between different circles by angular self-energies.

## Remaining obligation and reproducible controls

An unrestricted direct proof must establish (DG4), retaining enough of
\(\mathcal C_\rho\) or using a different mechanism. Global balance
alone does not remove \(\mathcal R_\rho\); equal circular masses do.
Proposition DG3 rules out proving the stronger mass-only inequality after
the angular Cauchy step throughout \((0,1)\). It does not rule out
stronger coupled angular estimates, a different direct proof, or the
classification of actual maximizers.

The supplementary [checker](check_three_cell_direct_global.py) checks the
quadrant formula against independent numerical integration, exact radial
moments, and finite-correlation values of \(E_\rho\). It also checks the
near-one limit with the factor \(s_\rho^2/(1-\rho)\to2\).
Run from the repository root:

```sh
python3 research/plurality-stability/check_three_cell_direct_global.py
```

The retained [verification record](three_cell_direct_global_verification.json)
states the environment, inputs and actual internal outcome. Numerical
quadrature is a normalization diagnostic, not an interval certificate or
a proof of (DG12). The unrestricted identities and obstruction use the
written arguments above. The circular comparison imports (DG0) with the
source-inspection limits stated above. No independent replay, formal
verification or external mathematical review is asserted.
