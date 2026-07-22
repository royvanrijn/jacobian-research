# Rank-two symplectic descent in every degree

This note proves rank-two completion for every normalized admissible seed in
every degree `N>=5`.  The weighted Laurent theorem on the two adapted charts
is the canonical source.  The degree-five and degree-six calculations are
generated specializations retained as readable formulas and independent
regressions; they are not separate theorem sources.  Stable moduli are
transferred by the coarse decorated-normalization map, not by full-cover
faithfulness.
These are internal repository theorems with exact certificates and no
recorded external specialist review.

## 1. Normalized seeds and fixed-third-component slices

A normalized degree-`N` weighted seed is

\[
 H(W)=\sum_{j=3}^N h_j(W^j-W^2),\qquad H'(1)=-1.
\]

It has dimension `N-3`.  Put

\[
 \kappa=H''(1),\qquad
 a=-\frac{1+\kappa}{2+\kappa}.
\]

The third component of the weighted map is

\[
 C=x(1+a,xy+x^2z).
\]

Thus fixing `gamma`, after fixing the harmless nonzero coefficient of `x^2z`,
is exactly the linear equation `H''(1)=kappa`.  It cuts the normalized seed
space from dimension `N-3` to `N-4`, provided this functional is nonconstant.
This proves the proposed dimension count for the *seed slice*; it does not by
itself prove completion or inequivalence of the resulting four-dimensional
maps.

For `a != 0`, every such slice has the same adapted quotient coordinate.  Set

\[
 W=Z+\psi(Q),\qquad Y=Q-\frac{XW}{3},
\]

and make the seed-dependent diagonal source change

\[
 (x,y,z)=\left(X,-\frac{3}{2a}Y,-\frac12W\right).
\]

Then, identically and independently of the seed and of `psi`,

\[
 2C=R=2X-3X^2Q.
\]

This is the common rank-two base needed in every degree.  The value `a=0`
(`kappa=-1`) requires another adapted chart; `kappa=-2` is outside the
weighted construction.

## 2. A rank-two flux obstruction

Let `K` be the rational function field of a bounded seed family and put

\[
 B=K[X,Q,Z],\qquad B_X=B[X^{-1}].
\]

On `B` use the quotient Poisson bracket

\[
 \{X,Z\}=-3X^2,\qquad \{Q,Z\}=6XQ-2,\qquad \{X,Q\}=0.
\]

For an adapted triple `Phi=(S,T,R)` with

\[
 \{S,T\}=1,\qquad \{R,S\}=\{R,T\}=0,
\]

let `w_Phi` be the horizontal derivation characterized by

\[
 w_\Phi(S)=w_\Phi(T)=0,\qquad w_\Phi(R)=1.
\]

Let `w_E={E,-}` be the fixed derivation in the common four-dimensional
coordinate system.  The difference

\[
 \delta_\Phi=w_\Phi-w_E
\]

has a normalized Hamiltonian `f_Phi in B_X`: integrate its `X`-component in
`Z` with zero `Z`-constant, and then integrate the remaining component after

\[
 v=X^{-1},\qquad \rho=R,qquad Q=\frac{v(2-\rho v)}3.
\]

This is the finite homotopy already used in the degree-five proof.  Define the
raw residue

\[
 \operatorname{Res}_X(\Phi)=\operatorname{PP}_{X=0}(f_\Phi),
\]

the sum of the negative powers of `X`.  For seeds and shears of bounded
degree, denominator and numerator degrees are uniformly bounded, so these
residues lie in a finite-dimensional Laurent space `V_{N,d}`.

Let `Gamma_{R,d}` be the bounded-degree group of base automorphisms preserving
`R`; it contains `Z -> Z+psi(Q)`.  The proposed obstruction is the orbit

\[
 \boxed{\operatorname{Flux}_{N,d}(H)
 = [\operatorname{Res}_X(\Phi_H)]\in V_{N,d}/\Gamma_{R,d}.}
\]

This quotient should initially be treated as an orbit set (or quotient
stack), not automatically as a vector-space quotient.  On an ansatz where
the shear response is linear, it can be replaced by the cokernel of the
linear response map.  The zero-orbit criterion is exact:

\[
 \operatorname{Flux}_{N,d}(H)=0
 \quad\Longleftrightarrow\quad
 \text{an allowed base change makes }f_\Phi\text{ polynomial}.
\]

Constants and linear terms of `psi` belong to the elementary base gauge, so
the first genuine coefficient is the coefficient of `Q^2`.

## 3. Canonical all-degree theorem

### Theorem — rank-two descent and stable moduli in every degree

For every `N>=5`, every normalized admissible seed (`kappa!=-2`) has an exact
symplectic rank-two completion in `A^4`, polynomially left--right equivalent
to the weighted map times `id_(A^1)`.  On `kappa!=-1,-2`, the complete Laurent
principal part lies on the universal four-residue line

\[
\Theta=X^{-4}+6QX^{-3}+\frac{45}{2}Q^2X^{-2}
+\frac{135}{2}Q^3X^{-1},
\]

and is canceled by a unique `Q^2` shear modulo elementary gauge.  On the
replacement chart `kappa=-1`, it lies on

\[
\Theta_{\rm exc}=X^{-4}-4QX^{-2}
\]

and is canceled by an `XQ` shear.

On the ordinary boundary-clean seed open, the
[coarse decorated-normalization map](DECORATED_NORMALIZATION_INVARIANT.md)
is generically etale of exact degree `N-2` onto the normalization of its
image, and that image has the full source dimension `N-3`.  Because rank-two
completion is an identity stabilization up to polynomial left--right
equivalence, it follows directly that the completed maps contain an
`(N-3)`-dimensional family of stable classes.  This conclusion is independent
of the affine-marked refinement below.

### Corollary — parameter-faithful rank-two moduli

After shrinking to the marked ordinary boundary-clean open, the assignment

\[
 H\longmapsto [G_H]_{\rm stable\ LR}
\]

is injective.  Indeed, `G_H` is polynomially left--right equivalent to
`F_H times id_(A^1)`, while the
[generic affine-mark faithfulness theorem](../papers/decorated-discriminant-normalization/main.tex)
recovers `H` from the stable class of `F_H`.  Thus the coarse theorem supplies
the `N-3`-dimensional result without faithfulness, and the affine mark upgrades
it to exact parameter separation on a nonempty open.

The proof of the two uniform residue statements is in Section 6.  Sections 4
and 5 instantiate them in degrees five and six as exact regression examples.

## 4. Degree-five regression formulas

The normalized degree-five surface has coordinates `(kappa,tau)`:

\[
 H_{\kappa,\tau}(W)=W^2(W-1)
 \left(
 \tau W^2+\left(\frac\kappa2-2\tau+2\right)W
 -\frac\kappa2+\tau-3
 \right).
\]

It satisfies `H'(1)=-1` and `H''(1)=kappa`.  It is useful in the exact
calculation to replace `kappa` by

\[
 a=-\frac{1+\kappa}{2+\kappa},\qquad
 \kappa=-\frac{1+2a}{1+a}.
\]

With `psi(Q)=s_2Q^2`, direct calculation over `Q(a,tau,s_2)` proves that the
localized Hamiltonian has denominator

\[
 15482880a^8(a+1)^2X^4
\]

and that its four low numerator coefficients are respectively

\[
 4096a^8\mathcal O,\quad
 24576Qa^8\mathcal O,\quad
 92160Q^2a^8\mathcal O,\quad
 276480Q^3a^8\mathcal O,
\]

where

\[
\begin{aligned}
\mathcal O={}&28(a+1)^2s_2+12a(a+1)^2\tau^2\\
&-18a(a+1)(4a+5)\tau
-216a^3-648a^2-738a-315.
\end{aligned}
\]

There are no lower powers of `X`.  Solving `mathcal O=0` and returning to
`kappa` gives

\[
 \boxed{
 s_2=K(\kappa,\tau):=
 \frac{
 12(\kappa+1)\tau^2
 -18(\kappa+1)(\kappa+6)\tau
 +9(\kappa^3+16\kappa^2+52\kappa+72)
 }{28(\kappa+2)}.}
\]

Up to a nonzero parameter unit, the raw residue is

\[
 (s_2-K(\kappa,\tau))
 \left(
 X^{-4}+6QX^{-3}+\frac{45}{2}Q^2X^{-2}
 +\frac{135}{2}Q^3X^{-1}
 \right).
\]

### Regression — generic-chart degree-five completion

On the chart `kappa != -2,-1`, the normalized Hamiltonian is polynomial
exactly for `s_2=K(kappa,tau)`.  The resulting four coordinates have the six
canonical brackets, and the source change is a polynomial automorphism.
Therefore every admissible weighted degree-five map on this chart has an
exact symplectic completion in `A^4`, polynomially left--right equivalent to
the original map times `id_(A^1)`.

Consequently:

* the completable locus on `kappa != -2,-1` is the entire dense chart of the
  degree-five seed surface, not a divisor;
* the completing quadratic shear is unique modulo constant and linear base
  gauge; and
* allowing higher-degree `psi` is unnecessary in degree five.

This correctly recovers the published line.  There `kappa=-9` and

\[
 \tau=\frac{3}{2(\lambda-1)},
\]

so the completing shear becomes

\[
 K\left(-9,\frac{3}{2(\lambda-1)}\right)
 =-\frac{27(57\lambda^2-138\lambda+73)}
 {196(\lambda-1)^2}.
\]

The published theorem is therefore the generic fixed-`gamma` line inside this
surface.  On the ordinary boundary-clean open, the coarse decorated-
normalization map is generically etale of exact degree three and completion
transfers its full two-dimensional image to stable moduli of exact symplectic
maps of `A^4`.  Thus degree
five is stronger than the fixed-`gamma` target `N-4=1`: it gives `N-3=2`
dimensions on this chart.

The exact certificate is
[`verify_degree_five_flux_surface.py`](../scripts/verify_degree_five_flux_surface.py).
The earlier slice-oriented exploratory replay remains
[`explore_degree_five_flux_surface.py`](../scripts/explore_degree_five_flux_surface.py).

### Exceptional divisors and the moduli open

The exact calculation separates four kinds of exceptional parameters:

1. `tau=0` is the degree-drop divisor: the coefficient of `W^5` in `H` is
   `tau`.
2. `kappa=-2` is genuinely outside the weighted marked-root construction,
   because `a=-(1+kappa)/(2+kappa)` is undefined.
3. `kappa=-1` gives `a=0`.  The generic-chart shear remains regular,
   `K(-1,tau)=45/4`, but that adapted source change divides by `a`.  The
   replacement chart below completes the entire divisor by a different
   `R`-preserving shear.
4. The ordinary-cusp condition removes the divisor

\[
\begin{aligned}
\Delta_{\rm cusp}(\kappa,\tau)={}&
\kappa^4-4\kappa^3\tau+18\kappa^3
+8\kappa^2\tau^2-36\kappa^2\tau+123\kappa^2\\
&-8\kappa\tau^3+24\kappa\tau^2-72\kappa\tau+376\kappa
+8\tau^4+24\tau^3-24\tau^2+72\tau+432.
\end{aligned}
\]

Indeed, `disc(H'')=432 Delta_cusp`.  The decorated-normalization theorem
removes additional proper incidence divisors for cusp-image collisions and
nonordinary bitangencies.  Their complement is already known to be a nonempty
open.  Intersecting it with

\[
 \tau(\kappa+2)\Delta_{\rm cusp}\ne0
\]

is therefore a nonempty two-dimensional open on which the maps have degree
five, the rank-two completion is defined, and the stable invariant is
generically finite.

### The `kappa=-1` replacement chart

At `kappa=-1`, the third component is

\[
 C=x(1+x^2z).
\]

Let `(X,Q,p,zeta)` have `{p,X}={zeta,Q}=1` and put

\[
\begin{aligned}
R&=2X+2X^3Q,\\
Z&=-2X^3p+(2+6X^2Q)\zeta,\\
E&=\frac{1-3X^2Q}{2}p+\frac92XQ^2\zeta.
\end{aligned}
\]

The momentum matrix has determinant `-1`, and

\[
\begin{gathered}
\{E,R\}=1,\qquad \{X,Z\}=2X^3,
\qquad \{Q,Z\}=-2(1+3X^2Q),\\
\{E,X\}=\frac{1-3X^2Q}{2},\qquad
\{E,Q\}=\frac92XQ^2,\qquad
\{E,Z\}=-6XQZ.
\end{gathered}
\]

Every polynomial automorphism of `A^3_(X,Q,Z)` preserving `R` fixes the two
components of `R=0`: their coordinate rings are `K[Q,Z]` and
`K[X,X^-1,Z]`, with different unit groups, so they cannot be interchanged.
It follows successively that

\[
 X\longmapsto X,\qquad Q\longmapsto Q,
 \qquad Z\longmapsto cZ+h(X,Q).
\]

Preservation of the quotient Poisson bracket (equivalently, preservation of
volume together with `R`) forces `c=1`.  Hence the full `R`-preserving Poisson
base group is

\[
 \boxed{Z\longmapsto Z+h(X,Q),\qquad h\in K[X,Q].}
\]

Use `W=Z+h(X,Q)` and the source change `(x,y,z)=(X,W,Q)`.  Then `2C=R`.
With `S=A/2` and `T=B`, direct calculation gives

\[
 \det\frac{\partial(S,T,R)}{\partial(X,Q,Z)}=-1,
 \qquad \{S,T\}=1,
 \qquad \{R,S\}=\{R,T\}=0.
\]

Write `A_tau=2tau^2-15tau-18`.  For the unsheared chart, the complete
negative principal part of the normalized Hamiltonian is

\[
 A_\tau\left(-\frac1{840}X^{-4}
 +\frac{Q}{210}X^{-2}\right).
\]

Restricting to the smaller ansatz `h=cQ^m`, `0<=m<=3`, gives a strict pole
filtration.  The highest new poles are

\[
 -\frac c4X^{-3},\quad
 \frac c{20}X^{-5},\quad
 -\frac c{112}X^{-7},\quad
 \frac c{660}X^{-9}
\]

for `m=0,1,2,3`, respectively.  In the associated pole filtration,
lower-degree shears cannot cancel the leading pole of a higher-degree shear.
Thus the filtered cubic ansatz successively forces

\[
 s_3=s_2=s_1=s_0=0,
\]

and polynomiality would then force

\[
 \boxed{2\tau^2-15\tau-18=0.}
\]

This explains why the former `psi(Q)` search saw only the two seed parameters

\[
 \tau=\frac{15\pm3\sqrt{41}}4
\]

within the cubic shear ansatz.

The full group contains the missing monomial `XQ`.  For

\[
 h(X,Q)=cXQ,
\]

the **complete** principal part is

\[
 \frac{105c-2A_\tau}{1680}
 \left(X^{-4}-4QX^{-2}\right).
\]

It is therefore canceled uniquely in this one-monomial family by

\[
 \boxed{c_\tau=\frac{2A_\tau}{105}
 =\frac{2(2\tau^2-15\tau-18)}{105}.}
\]

The source change `W=Z+c_tau XQ` is triangular and polynomially invertible.
The normalized Hamiltonian is polynomial, so the four output coordinates
have all six canonical brackets and are polynomially left--right equivalent
to the original weighted map times an identity.

### Corollary/regression — complete degree-five surface descent

Every normalized admissible degree-five seed (`kappa != -2`) has a rank-two
exact symplectic completion in `A^4`.  For `kappa != -1,-2` use the quadratic
`Q`-shear `K(kappa,tau)Q^2`; for `kappa=-1` use the replacement chart and the
shear `c_tau XQ` above.  Thus `kappa=-1` is only a chart divisor, not a
component of the flux obstruction locus.

The exact monomial replay is
[`explore_kappa_minus_one_flux.py`](../scripts/explore_kappa_minus_one_flux.py).
The command `--x-degree 1 --shear-degree 1` is the exceptional-chart
completion certificate.

## 5. Degree-six regression formulas

The normalized degree-six seed space has coordinates `(kappa,sigma,tau)`:

\[
H=W^2(W-1)\left(
\sigma W^3+\tau W^2
+\left(\frac\kappa2-3\sigma-2\tau+2\right)W
-\frac\kappa2+2\sigma+\tau-3
\right).
\]

On the generic chart put \(a=-(1+\kappa)/(2+\kappa)\).  With
\(W=Z+s_2Q^2\), the complete Laurent principal part lies in the universal
four-residue direction and has scalar obstruction

\[
\begin{aligned}
\mathcal O_{6,a}={}&308(a+1)^2s_2
+a(a+1)^2(753\sigma^2+627\sigma\tau+132\tau^2)\\
&-a(a+1)((2079a+2574)\sigma+(792a+990)\tau)\\
&-99(24a^3+72a^2+82a+35).
\end{aligned}
\]

The exact denominator is \(681246720a^{10}(a+1)^2X^4\), and the four low
numerator coefficients are

\[
16384a^{10}\mathcal O_{6,a},\quad
98304Qa^{10}\mathcal O_{6,a},\quad
368640Q^2a^{10}\mathcal O_{6,a},\quad
1105920Q^3a^{10}\mathcal O_{6,a}.
\]

Thus the unique completing quadratic shear is

\[
\boxed{
s_2=\frac{
-a(a+1)^2(753\sigma^2+627\sigma\tau+132\tau^2)
+a(a+1)((2079a+2574)\sigma+(792a+990)\tau)
+99(24a^3+72a^2+82a+35)
}{308(a+1)^2}.}
\]

This covers \(\kappa\ne-1,-2\).  On the missing admissible divisor
\(\kappa=-1\), use the replacement system of Section 3 and \(W=Z+cXQ\).
The complete denominator is \(18480X^4\).  Its only nonzero low numerator
coefficients are \(\mathcal E_6\) and \(-4Q\mathcal E_6\), in powers \(X^0\)
and \(X^2\), where

\[
\mathcal E_6=1155c-251\sigma^2-209\sigma\tau+858\sigma
-44\tau^2+330\tau+396.
\]

The exceptional divisor is therefore completed by

\[
\boxed{
c=\frac{251\sigma^2+209\sigma\tau-858\sigma
+44\tau^2-330\tau-396}{1155}.}
\]

### Corollary/regression — complete degree-six descent

Every normalized admissible degree-six seed \((\kappa\ne-2)\) has a
polynomial rank-two Hamiltonian after one explicit shear.  Use the displayed
\(Q^2\)-shear on \(\kappa\ne-1,-2\) and the displayed \(XQ\)-shear on
\(\kappa=-1\).  The resulting exact symplectic maps
\(\mathbb A^4\to\mathbb A^4\) are polynomially left--right equivalent to the
weighted maps times \(\operatorname{id}_{\mathbb A^1}\).

Fix `kappa=-9`, so again `a=-8/7` and the degree-five adapted coordinate
system applies.  With `W=Z+s_2Q^2`, the complete degree-six Laurent principal
part has the same universal four-residue direction as in degree five.  Its
scalar obstruction is

\[
\begin{aligned}
\mathcal O_6={}&2156s_2-6024\sigma^2-5016\sigma\tau
-11088\sigma\\
&-1056\tau^2-4752\tau+16929.
\end{aligned}
\]

Therefore the unique completing shear is

\[
\boxed{
s_2=\frac{
6024\sigma^2+5016\sigma\tau+11088\sigma
+1056\tau^2+4752\tau-16929}{2156}.}
\]

### Regression — degree-six fixed-`gamma` descent

Every normalized degree-six seed on the `kappa=-9` slice has a polynomial
rank-two Hamiltonian after the displayed quadratic shear.  For `sigma!=0`,
the resulting exact symplectic maps `A^4 -> A^4` have generic degree six and
form a two-dimensional algebraic family.  They are polynomially left--right
equivalent to the weighted maps times `id_(A^1)`.

This family contains a nonempty exact-double-zero boundary-clean open.  An
explicit witness is

\[
 (\sigma,\tau)=(1,0),\qquad
 H=\frac12W^2(W-1)(2W^3-11W+7),
\]

for which

\[
 \gcd(H,H')=W,qquad H''(0)=-7,qquad
 \operatorname{disc}(H/W^2)=\frac{1339}{4},qquad
 \operatorname{disc}(H'')=19184247360.
\]

On this open, the affine-mark theorem recovers the seed exactly.  Since
rank-two completion is polynomially left--right equivalent to adjoining one
identity coordinate, distinct `(sigma,tau)` parameters give distinct stable
classes after shrinking to the marked open.

### Corollary — two-dimensional degree-six stable moduli in `A^4`

There is a nonempty two-dimensional algebraic family of pairwise stably
polynomially left--right inequivalent degree-six exact symplectic Keller maps

\[
 \boxed{G_{\sigma,\tau}:\mathbb A^4\longrightarrow\mathbb A^4.}
\]

The exact certificate is
[`verify_degree_six_fixed_gamma_descent.py`](../scripts/verify_degree_six_fixed_gamma_descent.py).

The same witness proves that the marked ordinary boundary-clean open in the full
\((\kappa,\sigma,\tau)\)-space is nonempty.  The coarse decorated map is
generically etale of exact degree four onto a three-dimensional image, while
the generic-chart completion covers a dense open.  Hence degree six also has
a three-dimensional family of pairwise stably inequivalent exact symplectic
maps in \(\mathbb A^4\).  The dimension statement uses only the coarse map;
pairwise parameter separation uses its affine-marked refinement.

The full-space and exceptional-divisor certificates are
[`verify_degree_six_flux_surface.py`](../scripts/verify_degree_six_flux_surface.py)
and
[`verify_degree_six_kappa_minus_one_descent.py`](../scripts/verify_degree_six_kappa_minus_one_descent.py).

## 6. Uniform residue proof

The repeated four-residue direction is forced by the quotient grading.  Give
the localized base ring the weights

\[
\operatorname{wt}(X,Q,Z)=(1,-1,-2).
\]

For every normalized seed and quadratic shear, \((S,T,R)\) has weights
\((-2,-1,1)\).  The quotient bracket has degree \(3\), and the relative
horizontal vector field has component weights \((0,-2,-3)\).  The normalized
finite homotopy therefore gives

\[
f_H\in K[X,Q,Z,X^{-1}],\qquad \operatorname{wt}(f_H)=-4.
\]

Its negative-\(X\) part is independent of \(Z\): the second Hamiltonian
component is \((2-6XQ)\partial_Zf_H\), whose coefficient is a unit at \(X=0\),
so a negative power in \(\partial_Zf_H\) would contradict polynomiality of
the horizontal vector field.  Weight \(-4\) now forces

\[
(f_H)_-=c_4X^{-4}+c_3QX^{-3}+c_2Q^2X^{-2}+c_1Q^3X^{-1}.
\]

Polynomiality of the third component

\[
-3X^2\partial_Xf_H+(6XQ-2)\partial_Qf_H
\]

then gives

\[
c_3=6c_4,\qquad c_2=\frac{45}{2}c_4,
\qquad c_1=\frac{135}{2}c_4.
\]

Every raw residue on the generic chart is consequently a scalar multiple of

\[
\Theta=X^{-4}+6QX^{-3}+\frac{45}{2}Q^2X^{-2}
+\frac{135}{2}Q^3X^{-1}.                         \tag{A1}
\]

The \(R\)-preserving Poisson shear \(Z\mapsto Z+s_2Q^2\) has a
seed-independent flux cocycle.  Direct substitution gives

\[
\operatorname{Res}_X(s_2Q^2)=\frac{s_2}{135}\Theta. \tag{A2}
\]

If the unsheared seed has residue \(\mu(H)\Theta\), its unique completing
quadratic coefficient is \(s_2=-135\mu(H)\).

The exceptional chart \(\kappa=-1\) has the parallel grading

\[
\operatorname{wt}(X,Q,Z)=(1,-2,-1).
\]

Here \((S,T,R)\) again has weights \((-2,-1,1)\), the quotient bracket again
has degree \(3\), and the normalized Hamiltonian again has weight \(-4\).
The same unit-coefficient argument removes \(Z\) from its negative part, so

\[
(f_H)_-=d_4X^{-4}+d_2QX^{-2}.
\]

Polynomiality of

\[
2X^3\partial_Xf_H-(2+6X^2Q)\partial_Qf_H
\]

forces \(d_2=-4d_4\).  Thus the exceptional residue line is

\[
\Theta_{\rm exc}=X^{-4}-4QX^{-2}.                 \tag{A3}
\]

The seed-independent flux cocycle of \(Z\mapsto Z+cXQ\) is

\[
\operatorname{Res}_X(cXQ)=\frac{c}{16}\Theta_{\rm exc}. \tag{A4}
\]

Hence an unsheared exceptional residue \(\nu(H)\Theta_{\rm exc}\) is canceled
by \(c=-16\nu(H)\).  The degree-five and degree-six regression formulas above
are generated by specializing these two uniform identities.

### Proof of the canonical all-degree theorem

For every \(N\ge5\), every normalized admissible seed
\((\kappa\ne-2)\) has an exact symplectic rank-two completion in
\(\mathbb A^4\).  Use one quadratic \(Q^2\)-shear on
\(\kappa\ne-1,-2\) and one \(XQ\)-shear on \(\kappa=-1\).  In particular,
every fixed-\(\kappa=-9\) slice completes.  The proof is uniform in the
coefficients and degree of \(H\); degrees seven and eight have also been
replayed directly as exact probes.

The fixed-\(\kappa=-9\) slice has dimension \(N-4\) and meets the ordinary
exact-double-zero boundary-clean open.  Write \(H=W^2(W-1)P\), start from

\[
P_0(W)=\frac32-\frac52W,
\]

and perturb by \(\epsilon(W-1)^2W^{N-5}\).  This preserves \(P(1)=-1\) and
\(P'(1)=-5/2\), hence \(H''(1)=-9\), and gives exact degree \(N\) for generic
nonzero \(\epsilon\).  The bad root, critical-discriminant, and boundary
incidence conditions remove proper closed subsets; at \(\epsilon=0\), the
critical discriminant is \(216\).

The coarse decorated-normalization map has finite fibers on this open, so its
restriction to the fixed slice has image dimension \(N-4\).  On the full seed
open it is generically etale of exact degree \(N-2\) and has image dimension
\(N-3\).  Since completion is polynomially left--right equivalent to identity
stabilization, these are respectively \(N-4\) and \(N-3\) dimensions of stable
exact symplectic classes in \(\mathbb A^4\), without a faithfulness input.
After adding the affine root mark and shrinking once more, distinct seed
parameters give pairwise distinct stable classes on both opens.

The Laurent recurrence is checked by
[`verify_four_residue_recurrence.py`](../scripts/verify_four_residue_recurrence.py).
The degree-seven and degree-eight probes are generated by
[`explore_all_degree_fixed_gamma.py`](../scripts/explore_all_degree_fixed_gamma.py).

Thus the classical rank-two completion problem is closed on every admissible
normalized seed.  What remains classical is to express the scalar functions
\(\mu(H)\) and \(\nu(H)\) in a compact coefficient-invariant form and to
compare completions across the two charts intrinsically.

## 7. Regression-generation protocol

The uniform theorem makes new degree-by-degree searches unnecessary.  For an
independent fixed-degree replay or regression at `N>=5`:

1. Parametrize normalized seeds by `N-3` affine coordinates and use `kappa`
   as one coordinate.
2. Make the common diagonal change above on `a != 0` and compute the adapted
   triple `(S,T,R)` without expanding the source map first.
3. Choose `psi(Q)=sum_{j=2}^d s_jQ^j`, initially with `d=N-3`, and construct
   the normalized localized Hamiltonian by the finite homotopy.
4. Extract every coefficient of `X^{-1},...,X^{-M}`.  Clear only parameter
   denominators, not powers of `X`.
5. Eliminate the shear coefficients `s_j`.  The result must reproduce the
   universal one-dimensional residue response and its completing shear on the
   relevant chart; an extra seed condition signals a failed replay.
6. Repeat after imposing `kappa=kappa_0`.  Compare the rank of the shear
   response and the dimension of the resulting locus with `N-4`.
7. Transfer stable separation through the coarse decorated-normalization map;
   its exact generic degree is `N-2`, so parameter recovery is unnecessary.

The degree-five result warns against assuming that fixed `gamma` contributes
an additional completion equation: at least in degree five, `kappa` varies and
the shear absorbs the flux throughout the surface.

## 8. Related quantization result

The separate [filtered `A_2` quantization obstruction](RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md) studies the standard parity-preserving Weyl-symbol lift of one exact degree-five completion.  That specialized quantum calculation is not an input to the classical all-degree theorem.
