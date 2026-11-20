# Unbounded right-half-plane zero counts at each integer Gamma cusp

This is a consequence of the audited complete-cycle bounded-window expansion,
not a global localization or consecutive-coprimality theorem.

Fix s>=1, a=s^2, d=2s+1, and

\[
 A_m(x)=\sum_{k=0}^m\binom mk(am+1)_{dk}x^{m-k},\qquad
 w_c=a^d\left(\frac{d+1}{d-1}\right)^{d+1}>0.
\]

**Corollary CriticalGammaRightHalfPlaneZeroGrowth.** For every integer K>=1 there exists M(s,K) such that every
integer m>=M(s,K) has at least 2K distinct roots of A_m in Re(x)>0.
Consequently the number of right-half-plane roots tends to infinity.

The model Z(h)=integral_R exp(-u^4/12+h*u)du has infinitely many simple
imaginary roots, in conjugate pairs. Choose K such pairs and a fixed compact
h-disc enclosing them. The complete original-cycle expansion (B5) and the
finite-window Rouche argument in
`BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md` give one distinct actual
zero for every selected model root, with

\[
 x_j(m)/m^d=w_c+K_w h_j m^{-3/4}+O_{s,K}(m^{-1}).
\]

The finite set of actual normalized roots therefore all has real part
greater than w_c/2 for sufficiently large m. Multiplication by the positive
real number m^d preserves the half-plane, proving the assertion. All these
roots are nonreal, since A_m has positive coefficients and hence no positive
real root.

This rules out a global Hurwitz or left-half-plane theorem for the actual
critical families, rather than merely showing that a proposed proof lacks
its hypotheses. It gives no positive limiting proportion of such roots,
no estimate for M(s,K), no information on roots whose scaled index tends
to infinity, and no common root of consecutive polynomials. The normalized
positive real pinch point w_c is compatible with each finite polynomial
being strictly positive on the positive real line.

Inputs are the actual complete-cycle asymptotic/Rouche proof and the
quartic model infinitude/simplicity proof; the statement is not inferred
from a weak zero-measure limit or the stationary cover alone. The statement does not upgrade publication or verification assurance.
