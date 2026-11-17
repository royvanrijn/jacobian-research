# Zero linear harmonic and rank-one quadratic boundary

17 November 2026. This is a new exhaustive boundary calculation under the
explicit overarching goal of solving the field-valued two-pair(3,3) Image
problem. It is not a coefficient-support census.

**SIC2B33RANK1QUAD20261117.** Let
\(F=H_2+H_3\) on \(q=2zw+t^2=1\), where \(H_2\ne0\) has trace-free
self-adjoint quadratic operator of rank1, and \(H_3\) is arbitrary harmonic
cubic. After a complex determinant-one linear isometry normalizing
\(H_2=w^2\), the following are equivalent:

1. Every positive moment of \(F\) vanishes.
2. The first ten moments vanish.
3. The cubic has only torus charges \(-3,-2,-1\):
   \[
   F=w^2+a w^3+bw^2t+c w(t^2-zw/2).
   \]

Ten is the **least consecutive set-theoretic purity prefix** on this
stratum. No equality of ordinary moment ideals or nilpotent scheme
stabilization is asserted.

All such pure phases are one-sided, hence lie in the SO3 nullcone and are
safe for every original multiplier. The full Image cutoff is

\[
 \mathcal E(P^mQ)=0\qquad(m>\deg_{Z,Y}Q).
\]

The strict cutoff is attained. Including \(H_2=0\) via the retained
[odd-slice theorem](ODD_HOPF_CUBIC_IMAGE_TRANSFER.md) therefore gives an
explicit safe container: \(H_1=0\), \(\operatorname{rank}B\le1\), with
\(H_3\) unrestricted before imposing purity. It is not asserted that
every coefficient point in that container is pure.

The nonzero-rank1 calculation below uses exact rational polynomial
identities and no msolve, numerical solver or Groebner-unit oracle.
It leaves the other \(D=0\) strata and the full(3,3) problem open.

## 1. A complete normal form

Self-adjoint rank1 means \(B=\kappa n\otimes n^\flat\). Trace zero gives
\(q(n)=0\). Over an algebraic closure, a determinant-one isometry and a
nonzero isotropic rescaling put \(H_2(x)=q(Bx,x)\) in the form \(w^2\).
All normalization steps preserve vanishing of moments and full all-marker
safety. Field extension introduces no loss for a characteristic-zero
field-valued identity.

Write the complete cubic as

\[
\begin{aligned}
H_3={}&\gamma_{-3}w^3+\gamma_{-2}w^2t
 +\gamma_{-1}w(t^2-zw/2)+\gamma_0(t^3-3ztw)\\
 &+\gamma_1z(t^2-zw/2)+\gamma_2z^2t+\gamma_3z^3.
\end{aligned}
\]

The indices are the charges under \((z,t,w)\mapsto(\rho z,t,\rho^{-1}w)\).
These seven harmonics span the entire cubic space.

The unipotent isometry

\[
 (z,t,w)\mapsto(z-st-s^2w/2,\ t+sw,\ w)
\]

fixes \(w^2\). It fixes the highest cubic charge, and changes the next
coefficient by

\[
 \gamma_2\mapsto\gamma_2-3s\gamma_3,\qquad
 \gamma_1\mapsto\gamma_1-2s\gamma_2\quad(\gamma_3=0).
\]

Combining the torus with phase amplitude \(\rho^2\) fixes the quadratic
part and sends \(\gamma_j\mapsto\rho^{j+2}\gamma_j\).
Thus every nonzero positive highest charge has one of the exhaustive forms:

- highest3: \(\gamma_3=1,\gamma_2=0\), with all five lower coefficients free;
- highest2: \(\gamma_3=0,\gamma_2=1,\gamma_1=0\), with all four lower coefficients free;
- highest1: \(\gamma_3=\gamma_2=0,\gamma_1\ne0\).

Only a coefficient already declared nonzero is divided by or normalized.
No cubic coefficient beyond the next coefficient is silently discarded.
These triangular charts are a complete orbit normalization of this
algebraically defined boundary family.

The raw first relevant moments are

\[
\begin{aligned}
\mu_2&=\gamma_0^2/7+\gamma_1\gamma_{-1}/21
       +4\gamma_2\gamma_{-2}/105+4\gamma_3\gamma_{-3}/35,\\
\mu_3&=(-8\gamma_0\gamma_2+2\gamma_1^2-4\gamma_3\gamma_{-1})/105.
\end{aligned}
\]

If the highest positive charge is1, \(\mu_3=2\gamma_1^2/105\ne0\).
If there is no positive charge, \(\mu_2=\gamma_0^2/7\), forcing
\(\gamma_0=0\). It remains to exclude exactly the two higher charts.

## 2. Highest charge two: first-five unit identity

Set \(\gamma_2=1,\gamma_1=\gamma_3=0\). Use the primitive integer moment
polynomials, which differ from the raw moments only by nonzero rational
scalars. The first two are

\[
 M_2=15\gamma_0^2+4\gamma_{-2},\qquad M_3=-\gamma_0.
\]

They set \(\gamma_0=\gamma_{-2}=0\) with constant pivots. With
\(a=\gamma_{-3},b=\gamma_{-1}\), the next two reduced equations are

\[
 f_4=16ab+208,\qquad f_5=b^2.
\]

The exact rational identity

\[
 \frac{13-ab}{2704}f_4+\frac{a^2}{169}f_5=1
\]

excludes this entire chart, including every value of both remaining
cubic parameters. No localization is used in this unit identity.

## 3. Highest charge three: explicit first-ten certificate

Set \(\gamma_3=1,\gamma_2=0\), and write \(y=\gamma_0,x=\gamma_1\).
The first two equations are

\[
 M_2=15y^2+5x\gamma_{-1}+12\gamma_{-3},\qquad
 M_3=x^2-2\gamma_{-1}.
\]

They solve, with constant nonzero pivots,

\[
 \gamma_{-1}=x^2/2,\qquad
 \gamma_{-3}=-5y^2/4-5x^3/24.
\]

Put \(u=\gamma_{-2}+yx/2\), an invertible polynomial coordinate change,
and \(v=27y^2+4x^3\). After nonzero rational rescaling, moments4,5,6,7
generate exactly the same ideal as

\[
\begin{aligned}
A_0&=u^2+5,\\
A_1&=v^2-2112x,\\
A_2&=287xv-2358yu,\\
A_3&=x^2v^2-36yxuv-1728y^2-16x^3.
\end{aligned}
\]

The certificate exports the complete invertible four-by-four polynomial
generator change. Its determinant is a nonzero rational constant.

Here is a direct ideal calculation. In the ideal \((A_0,A_1,A_2,A_3)\),
the following combinations are

\[
\begin{aligned}
B&=-\frac{2358A_3+7974x^2A_1-36xvA_2}{288}
   =14148y^2+58607x^3,\\
C&=(287xv+2358yu)A_2-287^2x^2A_1+2358^2y^2A_0\\
 &=5\cdot2358^2y^2+287^2\cdot2112x^3.
\end{aligned}
\]

Their constant coefficient determinant is831910506804, so both \(y^2\)
and \(x^3\), hence \(v\), belong to the ideal. Then \(A_1\) gives
\(x\), \(A_2\) gives \(yu\), and \(yA_0-u(yu)=5y\) gives \(y\).
The reverse containment is immediate. Thus the first-seven ideal on
this chart is exactly

\[
 (u^2+5,x,y).
\]

Moments8 and9 reduce to zero in this quotient. The primitive reduced
tenth moment has remainder \(-771876864\), a nonzero rational unit.
The exported certificate contains explicit polynomial multipliers
\(H_0,\ldots,H_3\) satisfying

\[
 \sum_{i=0}^3H_iA_i-\frac{f_{10}}{771876864}=1.
\]

The checker verifies this identity by polynomial expansion, together
with the generator change and all monic eliminations. It does not infer
the exclusion merely from a CAS output saying the ideal is one.

## 4. Sharp set-theoretic prefix and full marked safety

The two exact phases

\[
 F=w^2+z^3+\alpha w^2t,\qquad\alpha^2=-5,
\]

satisfy every moment through9, but

\[
 h(F^{10})=-\frac{2310144}{434113615}\ne0.
\]

They prove that prefix10 cannot be replaced by an earlier consecutive
prefix on the full stratum. They are finite-prefix false positives,
not new all-order pure phases. Since the only remaining field points
have negative charges, they are pure to every order automatically.
This proves the set-theoretic equivalence claimed at the start.

For a spherical polynomial marker of degree at most \(b\), its maximum
charge is at most \(b\). Each phase factor has charge at most \(-1\), so
every marked moment vanishes for \(m>b\). Linear rotations preserve
marker degree. The retained Hopf coefficient-balancing identity converts
every original bihomogeneous multiplier of dual/coordinate degrees
\((a,b)\) to such balanced markers; \(a>b\) gives zero by degree.
Finite sums handle arbitrary original multipliers. Hence this is the
full Image conclusion with cutoff \(m>\deg_{Z,Y}Q\).

For \(F=w^2+w(t^2-zw/2)\),
\(h(F^2z^2)=2/315\ne0\), so the strict boundary is attained. More
generally for positive even \(b\), the surviving charge-zero contribution
is \(h((zw)^b(t^2-zw/2)^b)>0\) in standard real-sphere coordinates.

## Evidence and assurance

- [Exact moment generator](boundary_rank1_quadratic.py).
- [Complete checker](verify_boundary_rank1_quadratic.py).
- [Polynomial unit identities and generator changes](boundary_rank1_quadratic_certificate.json).
- [Verification receipt](boundary_rank1_quadratic_verification.json).
- [Explicit safe-container equations](boundary_rank1_safe_container.json).

Run `python verify_boundary_rank1_quadratic.py`. The moment engine uses
torus balance and an exact finite binomial sum. Its formula is cross-checked
against direct spatial integration through6, while the false candidates
are independently integrated through10. Every polynomial unit identity
and coordinate change is checked exactly over \(\mathbb Q\).

The core nonzero-rank1 proof has no msolve dependency. The zero-quadratic
extension of the safe container uses the existing odd-harmonic theorem
and retains that theorem's assurance boundary. No ordinary coefficient-ring
stabilization, other \(D=0\) stratum, unrestricted cyclic-chart exhaustion
or full(3,3) Image theorem is claimed. Independent specialist, external and
formal verification remain outstanding.
