# The finite LR Rees/SAGBI module computation

This note implements the finite module problem in `OP-LR-REES` and
`OP-LR-II` for the degree-five base map \(F_2\).  It does not extend the
arbitrary jet expansion.  The calculation is exact over \(\mathbf Q\), uses
the ordinary-source-degree face compatible with

\[
 \deg u=2,\qquad \deg\gamma=3,
\]

and enlarges target-invariant coefficients to
\(R=\mathbf Q[u,\gamma]\) when forming the saturated normal quotient.
This is the same enlargement used in
[the torus-module note](TORUS_FILTERED_LR_MODULE.md).  Nonvanishing or
independence after this enlargement survives in the actual target quotient.

The complete matrices and all generator labels are stored in
[`lr_rees_sagbi_module_computation.json`](../artifacts/generated-results/lr_rees_sagbi_module_computation.json).

## 1. The target-field initial algebra

The target coordinate weights are

\[
 \operatorname{wt}(A,B,C)=(-2,-1,1).
\]

The invariant ring is the polynomial ring

\[
 T=\mathbf Q[P,Q],\qquad P=AC^2,\quad Q=BC.
\]

Under pullback by \(F_2\),

\[
 P\longmapsto a\gamma^2,\qquad Q\longmapsto b\gamma.
\]

Normalize

\[
 f_1=\frac{a\gamma^2}{6},\qquad f_2=4b\gamma.
\]

Their initial monomials are \(u^5\gamma^5\) and
\(u^4\gamma^4\).  The unique toric relation between these two monomials is
\(T_1^4-T_2^5\).  Exact subduction produces a third element \(f_3\), defined
by

\[
\begin{aligned}
 -\frac1{90}f_3={}&
 f_1^4-f_2^5-\frac{14}{3}f_1^3f_2
 \frac{461}{90}f_1^2f_2^2-\frac{739}{450}f_1f_2^3\\
 &-\frac{527}{720}f_2^4
 -\frac{1636136}{759375}f_1^3
 \frac{675947}{202500}f_1^2f_2 .
\end{aligned}
\]

Its initial monomial is

\[
 \operatorname{in}(f_3)=u^{12}\gamma^{14}.
\]

This exponent is off the diagonal containing the first two exponent vectors,
so it creates no new toric relation.  Hence

\[
 \boxed{\{f_1,f_2,f_3\}\text{ is a finite SAGBI basis}.}
\]

For an integer weight \(q\), the coefficient semi-invariants have the
following minimal \(T\)-generators:

\[
\begin{array}{c|c}
q&\text{generators}\\ \hline
q>0&C^q\\
q=0&1\\
q=-n<0&A^iB^{n-2i}\ (0\le i\le\lfloor n/2\rfloor),
 \quad A^{(n+1)/2}C\ \text{if \(n\) is odd}.
\end{array}
\]

The target-field module of field weight \(p\) is the direct sum of the
coefficient modules of weights \(p-2,p-1,p+1\) in the
\(\partial_A,\partial_B,\partial_C\) components.  This constructs every
weight module without a target-degree box.

## 2. Linear lifting does not commute with the initial face

The artifact records the initial lift of every module generator in weights
\(p=\pm1,\pm2\).  Already in weight \(p=1\), two generator lifts have degree
39 and proportional initial forms:

\[
 3\,\operatorname{in}\ell_{F_2}(AC\,\partial_A)
 =
 4\,\operatorname{in}\ell_{F_2}(C^2\,\partial_C).
\]

Their exact combination drops five degrees:

\[
 \deg\left(
 3\ell_{F_2}(AC\,\partial_A)
 -4\ell_{F_2}(C^2\,\partial_C)
 \right)=34,
\]

with new initial lift

\[
 -210x^{20}y^7z^6(x,0,-3z).
\]

This intermediate initial form is \( -28\) times the initial lift of
\(B\partial_A\).  One more exact subduction therefore gives

\[
 H=(28B+3AC)\partial_A-4C^2\partial_C
\]

with

\[
 \boxed{
 \deg\ell_{F_2}(H)=29,\qquad
 \operatorname{in}\ell_{F_2}(H)
 =-274x^{17}y^6z^5(x,0,-3z).
 }
\]

The four original weight-one generator lifts have degrees
\(19,34,39,39\).  The positive generators of the invariant initial algebra
have degrees \(20,25,66\).  Consequently no invariant multiple of an
original initial lift can have degree \(29\).  The displayed initial form is
therefore a genuinely new generator of the lifted initial submodule.

Thus the associated-graded linear lift acquires a kernel although the
unfiltered target lift is injective:

\[
 \boxed{\text{taking the initial face does not commute with linear target
 lifting}.}
\]

This is an explicit linear Rees-torsion class.  In particular,
`OP-LR-REES` has a negative answer for this face; adding the scalar SAGBI
generator \(f_3\) does not remove the module-level degree drop.

## 3. Saturated weight-zero normal quotient

Transport weight-zero source fields through the logarithmic differential
matrix \(J\).  After extending coefficients to \(R\), the target image and
normal quotient are

\[
 H_R=(a,b^2)e_A\oplus(b,a\gamma)e_B\oplus(\gamma)e_C,
\]

\[
 \boxed{
 N_R=
 R/(a,b^2)\oplus R/(b,a\gamma)\oplus R/(\gamma).
 }
\]

For the selected weighted face, exact standard bases give

\[
\begin{aligned}
 \operatorname{in}(a,b^2)
   &=(u^3,u^2\gamma,u\gamma^2,\gamma^3),\\
 \operatorname{in}(b,a\gamma)
   &=(u^3,u\gamma,\gamma^2),\\
 \operatorname{in}(\gamma)&=(\gamma).
\end{aligned}
\]

These ideals are also recorded in the certificate.  The quadratic matrices
below are reduced in the full quotient \(N_R\), which is stronger than
comparing only selected leading coefficients.

## 4. The finite quadratic matrices

The second fundamental form is \(T\)-bilinear in coefficient multipliers:

\[
 \operatorname{II}_{F_2}(fY,gZ)
 =(f\circ F_2)(g\circ F_2)\operatorname{II}_{F_2}(Y,Z).
\]

It is therefore enough to compute the six constant-component tensors and
multiply them by the finite semi-invariant generators above.

There are two surviving positive weights:

\[
\begin{array}{c|c|c}
p&\text{matrix size}&\text{nonzero columns in }N_R\\ \hline
1&3\times24&11\\
2&3\times24&8.
\end{array}
\]

For \(p\ge3\), every product of a weight-\(p\) generator and a
weight-\(-p\) generator has multiplier \(a^\alpha b^\beta\gamma^\chi\)
satisfying simultaneously

\[
 \alpha\ge1\text{ or }\beta\ge2,\qquad
 \beta\ge1\text{ or }(\alpha\ge1\text{ and }\chi\ge1),\qquad
 \chi\ge1.
\]

These are precisely the three elementary annihilator tests for the summands
of \(N_R\).  Consequently

\[
 \boxed{\operatorname{II}_{(F_2,p,-p)}=0\text{ in }N_R
 \quad\text{for every }|p|\ge3.}
\]

This is a semigroup cutoff, not the outcome of the regression range
\(3\le p\le8\) included in the checker.

## 5. A new independent \(p=2\) invariant

Let \(I_1\subset N_R\) be the \(R\)-submodule generated by all columns of the
\(p=1\) matrix.  Exact Gröbner-module reduction shows that seven of the eight
nonzero \(p=2\) columns lie in \(I_1\).  The exception is

\[
 \operatorname{II}_{F_2,2,-2}(\partial_A,A^2\partial_A).
\]

Modulo \(I_1\), together with the five defining relation vectors of \(N_R\),
its exact normal form is

\[
 \boxed{-\frac{987}{395}\,e_C\ne0.}
\]

There is also a compact separator independent of the module term order.
Evaluate at

\[
 (u,\gamma)=\left(\frac16,0\right)
\]

and apply the covector

\[
 \lambda=\left(0,-\frac{144}{79},1\right).
\]

The first entry of \(\lambda\) kills the first summand relations.  At the
selected point \(b=0\), \(a\gamma=0\), and \(\gamma=0\), so \(\lambda\)
descends through the remaining normal relations.  Direct rational evaluation
gives

\[
 \lambda(\operatorname{im}\operatorname{II}_{F_2,1,-1})=0,
\]

on all 24 matrix columns, whereas

\[
 \boxed{
 \lambda\!\left(
 \operatorname{II}_{F_2,2,-2}(\partial_A,A^2\partial_A)
 \right)=-\frac{987}{395}.
 }
\]

Every other \(p=2\) column evaluates to zero.  Thus the same scalar appears
both as the Gröbner-module remainder and as a one-point separating
functional.  The dependency-free checker
[`audit_lr_rees_sagbi_module_certificate.py`](../scripts/audit_lr_rees_sagbi_module_certificate.py)
replays this evaluation with standard-library `Fraction` arithmetic.

The exact module quotient gives the complete scheme-theoretic refinement.
The annihilator of the exceptional column modulo the full \(p=1\) image is

\[
 \operatorname{Ann}_R([\operatorname{II}_{2,-2}
  (\partial_A,A^2\partial_A)])=(\gamma,6u-1).
\]

Since every other \(p=2\) column belongs to the \(p=1\) image, it follows that

\[
 \boxed{
 \frac{\operatorname{im}\operatorname{II}_{p=1}
       +\operatorname{im}\operatorname{II}_{p=2}}
      {\operatorname{im}\operatorname{II}_{p=1}}
 \simeq R/(\gamma,6u-1)\simeq\mathbf Q.
 }
\]

Thus the extra quadratic image is one reduced length-one skyscraper at
\((u,\gamma)=(1/6,0)\); it has no hidden nilpotent thickness.

Therefore the answer to the decisive question is negative:

\[
 \boxed{\text{the known \(p=1\) class does not generate the entire
 quadratic normal image}.}
\]

Indeed, even the full \(p=1\) matrix does not generate the \(p=2\) image.
The class above is a new independent quadratic LR invariant.  The single
known column
\(\operatorname{II}_{F_2}(\partial_B,\partial_C)\) is also noncyclic: it
fails to generate several \(p=1\) columns as well as the new \(p=2\) class.

This conclusion is an exact saturated-module computation for the displayed
face.  It does not by itself classify higher mixed BCH forms or prove a
global statement for other filtrations.
