# Double-conic normal layers

## Status

This note proves `HC4NHM15`.  It continues the nonzero-restriction row left
by `HC4NHM13`; it does not close the complete double-conic packet.
The subsequent theorem `HC4NHM18` removes the finite exceptional locus left
here in \((3,3,2,2)\).  Therefore the current double-conic frontier starts
with restrictions supported on at least five points; see
[`HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md`](HC4_DOUBLE_CONIC_BALANCED_FOUR_ROOT_CLOSURE.md).

Let

\[
 q=xz-y^2,
 \qquad h_5=\operatorname{lift}(f_{10})+qG_3,
 \qquad f_{10}\ne0.
\tag{0.1}
\]

The results are:

1. the four conditions for \(q^4\mid\det\operatorname{Hess}(h_5)\) are
   computed as explicit binary covariants;
2. every decic with at most three distinct roots is excluded from the clean
   target \(\det\operatorname{Hess}(h_5)=q^4\ell\), \(\ell\ne0\);
3. every four-root decic at harmonic cross-ratio is excluded as well.
4. eight of the nine four-root partitions are excluded for arbitrary
   cross-ratio;
5. the ninth partition \((3,3,2,2)\) is empty over the function field
   \(\mathbf Q(\lambda)\), leaving only a finite exceptional cross-ratio
   locus.

Thus only a finite exceptional locus in \((3,3,2,2)\), together with the
partitions having at least five distinct roots, remains.  These are
candidates, not constructed survivors.  No Schur equation is imposed in
this note.

The bounded exact replay is

```bash
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-18
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-14
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-10
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group layer-6
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group support-one-two
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group support-three
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group harmonic-four
# long generic-cross-ratio replays
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-7111
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-6211
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-5311
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-5221
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4411
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4321
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-3331
# long one-chart polynomial replay and function-field replay
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group generic-four-4222-root-chart
.venv/bin/python scripts/verify_hc4_double_conic_normal_layers.py --group balanced-function-fields
```

SymPy performs exact rational harmonic decomposition and covariant
interpolation.  Singular 4.4.1 performs the exact rational modular standard
bases, including its final characteristic-zero verification.

## 1. Equivariant splitting

Put

\[
 \Box=\partial_x\partial_z-\frac14\partial_y^2.
\tag{1.1}
\]

For a binary form \(A_{2m}\), let \(H_m(A_{2m})\) be the unique ternary form
of degree \(m\) satisfying

\[
 \Box H_m(A_{2m})=0,
 \qquad
 H_m(A_{2m})(s^2,st,t^2)=A_{2m}(s,t).
\tag{1.2}
\]

This fixes the normalization of `lift`.  The harmonic decomposition gives

\[
 S^5(S^2U)=S^{10}U\oplus qS^6U\oplus q^2S^2U.
\tag{1.3}
\]

Consequently there are unique binary forms \(f=f_{10}\), \(g=g_6\), and
\(k=k_2\) such that

\[
 h_5=H_5(f)+qH_3(g)+q^2H_1(k).
\tag{1.4}
\]

Equivalently, \(G_3=H_3(g)+qH_1(k)\).  The Hessian determinant has the unique
decomposition

\[
 \det\operatorname{Hess}(h_5)
 =H_9(\Phi_{18})+qH_7(\Phi_{14})+q^2H_5(\Phi_{10})
   +q^3H_3(\Phi_6)+q^4H_1(\Phi_2).
\tag{1.5}
\]

It follows immediately that

\[
 \boxed{
 q^4\mid\det\operatorname{Hess}(h_5)
 \iff
 \Phi_{18}=\Phi_{14}=\Phi_{10}=\Phi_6=0.
 }
\tag{1.6}
\]

The remaining quotient is the line \(H_1(\Phi_2)\).  The clean packet also
requires \(\Phi_2\ne0\).

## 2. Transvectant convention

For binary forms \(A,B\), use the unnormalized transvectant

\[
 (A,B)_r=
 \sum_{i=0}^r(-1)^i\binom ri
 \partial_s^{r-i}\partial_t^iA\,
 \partial_s^i\partial_t^{r-i}B.
\tag{2.1}
\]

Write

\[
 [AB;r;C;s]=((A,B)_r,C)_s.
\tag{2.2}
\]

The coefficients below depend on (2.1); the zero loci do not depend on a
simultaneous change to another standard transvectant normalization.

## 3. The four normal layers

The first layer is

\[
\begin{aligned}
\Phi_{18}={}&
-\frac{[ff;4;f;2]}{56582064}
+\frac{[ff;6;f;0]}{18860688}\\
&+\frac{[ff;2;g;2]}{12150}
-\frac{19[ff;4;g;0]}{317520}
+\frac{16[ff;2;k;0]}{81}\\
&-\frac{8[gg;0;f;2]}{99}
+\frac{6[gg;2;f;0]}{55}
-80[fg;0;k;0]+32[gg;0;g;0].
\end{aligned}
\tag{3.1}
\]

The second layer is

\[
\begin{aligned}
\Phi_{14}={}&
-\frac{[ff;0;f;8]}{64818903663360}
+\frac{191[ff;2;f;6]}{66292060564800}\\
&+\frac{[ff;2;g;4]}{10024560}
+\frac{[ff;4;g;2]}{39584160}
-\frac{19[ff;6;g;0]}{227026800}\\
&+\frac{5[ff;2;k;2]}{4131}
+\frac{[ff;4;k;0]}{95256}
-\frac{[gg;0;f;4]}{33660}\\
&-\frac{3[gg;2;f;2]}{4675}
+\frac{[gg;4;f;0]}{360}
-\frac{25[fg;0;k;2]}{51}\\
&-\frac{9[fg;1;k;1]}{14}
-\frac{16[fg;2;k;0]}{945}
-160[kk;0;f;0]\\
&-\frac{27[gg;0;g;2]}{680}
+112[gg;0;k;0].
\end{aligned}
\tag{3.2}
\]

The third layer is

\[
\begin{aligned}
\Phi_{10}={}&
-\frac{[ff;0;f;10]}{608083813148797440}
+\frac{4399[ff;2;f;8]}{1303036742461708800}\\
&+\frac{29[ff;2;g;6]}{113837724000}
+\frac{53[ff;4;g;4]}{367783416000}
-\frac{31[ff;6;g;2]}{619783164000}\\
&-\frac{[ff;8;g;0]}{3960744480}
+\frac{109[ff;4;k;2]}{136216080}
-\frac{17[ff;6;k;0]}{227026800}\\
&-\frac{[gg;0;f;6]}{68108040}
-\frac{[gg;2;f;4]}{126126000}
-\frac{467[gg;4;f;2]}{26535600}\\
&+\frac{29[gg;6;f;0]}{75600}
+\frac{62[fg;2;k;2]}{19305}
-\frac{[fg;3;k;1]}{975}\\
&-\frac{263[fg;4;k;0]}{300300}
-\frac{32[kk;0;f;2]}{39}
-40[kk;2;f;0]\\
&-\frac{31[gg;0;g;4]}{409500}
+\frac{4[gg;0;k;2]}{143}
-\frac{18[gg;2;k;0]}{275}\\
&+144[kk;0;g;0].
\end{aligned}
\tag{3.3}
\]

The fourth layer is

\[
\begin{aligned}
\Phi_6={}&
\frac{[ff;2;f;10]}{77203623553056000}
+\frac{277[ff;4;f;8]}{19455313135370112000}\\
&+\frac{101[ff;4;g;6]}{101949562915200}
+\frac{19[ff;6;g;4]}{35399153790000}
-\frac{[ff;8;g;2]}{668375631000}\\
&-\frac{[ff;10;g;0]}{108020304000}
+\frac{19[ff;6;k;2]}{28605376800}
-\frac{53[ff;8;k;0]}{118822334400}\\
&-\frac{41[gg;0;f;8]}{2022808788000}
+\frac{397[gg;2;f;6]}{374594220000}
+\frac{377[gg;4;f;4]}{2829103200}\\
&+\frac{23[fg;4;k;2]}{7567560}
-\frac{[fg;5;k;1]}{249480}
-\frac{829[fg;6;k;0]}{78586200}\\
&-\frac{13[kk;0;f;4]}{4158}
-\frac{1747[gg;0;g;6]}{5096520000}
+\frac{211[gg;2;g;4]}{217800000}\\
&+\frac{8[gg;2;k;2]}{1925}
-\frac{307[gg;4;k;0]}{113400}
-\frac{4[kk;0;g;2]}{15}\\
&+8[kk;2;g;0]+96[kk;0;k;0].
\end{aligned}
\tag{3.4}
\]

### Why the interpolation is exact

Each \(\Phi_j\) is an \(SL_2\)-equivariant cubic polynomial in
\((f,g,k)\).  Clebsch--Gordan decomposition gives the following dimensions
for its multihomogeneous blocks.  The columns are ordered as

\[
 f^3, f^2g, f^2k, fg^2, fgk, fk^2, g^3, g^2k, gk^2, k^3.
\]

\[
\begin{array}{c|rrrrrrrrrr}
 &f^3&f^2g&f^2k&fg^2&fgk&fk^2&g^3&g^2k&gk^2&k^3\\ \hline
\Phi_{18}&2&3&2&2&1&0&1&0&0&0\\
\Phi_{14}&2&4&2&3&3&1&1&1&0&0\\
\Phi_{10}&2&4&2&4&3&2&1&2&1&0\\
\Phi_6&2&4&2&3&3&1&2&2&2&1
\end{array}
\tag{3.5}
\]

The checker constructs the natural iterated-transvectant spanning family in
every block.  Its exact evaluation matrices have precisely the ranks in
(3.5).  Evaluation is therefore injective on every possible covariant block,
and the exact Hessian values determine (3.1)--(3.4) uniquely.  This is finite
interpolation inside a proved covariant space, not agreement on a bounded
search box.

## 4. Root partitions with at most three support points

For the clean double-conic target, introduce a general ternary cubic \(G_3\)
and a general residual line

\[
 \ell=Ax+By+Cz.
\tag{4.1}
\]

For a fixed nonzero decic \(f\), let \(I_f\) be the 55-coefficient ideal of

\[
 \det\operatorname{Hess}(\operatorname{lift}(f)+qG_3)-q^4\ell.
\tag{4.2}
\]

Scaling normalizes the coefficient of \(f\).  The \(PGL_2\)-normal forms for
one, two, and three distinct roots are

\[
 s^{10},
 \qquad s^at^b,
 \qquad s^at^b(s-t)^c.
\tag{4.3}
\]

There are fourteen partitions:

\[
\begin{gathered}
(10);\\
(9,1),(8,2),(7,3),(6,4),(5,5);\\
(8,1,1),(7,2,1),(6,3,1),(6,2,2),\\
(5,4,1),(5,3,2),(4,4,2),(4,3,3).
\end{gathered}
\tag{4.4}
\]

For every partition \(\lambda\) in (4.4), exact rational modular standard
bases give

\[
 \boxed{
 I_\lambda+(uA-1)=(1),\quad
 I_\lambda+(uB-1)=(1),\quad
 I_\lambda+(uC-1)=(1).
 }
\tag{4.5}
\]

The three Rabinowitsch charts cover \(\ell\ne0\).  Hence all fourteen strata
are empty.  In particular none proceeds to the Schur equations.

## 5. Four roots: the first modulus and the harmonic orbit

Four support points introduce a cross-ratio.  The nine partitions are

\[
\begin{gathered}
(7,1,1,1),(6,2,1,1),(5,3,1,1),(5,2,2,1),(4,4,1,1),\\
(4,3,2,1),(4,2,2,2),(3,3,3,1),(3,3,2,2).
\end{gathered}
\tag{5.1}
\]

Use representatives

\[
 f=s^at^b(s-t)^c(s-\lambda t)^d.
\tag{5.2}
\]

At the harmonic cross-ratio orbit one may take \(\lambda=2\).  Repeating the
three exact line-chart calculations (4.5) gives the unit ideal for every
partition in (5.1).  Thus the complete harmonic four-root row is empty.

Retain \(\lambda\) as a polynomial variable and localize each residual-line
chart by

\[
 uA\lambda(\lambda-1)-1,
 \quad uB\lambda(\lambda-1)-1,
 \quad uC\lambda(\lambda-1)-1.
\tag{5.3}
\]

Exact characteristic-zero modular reconstruction gives the unit ideal in all
three charts for

\[
\boxed{
\begin{gathered}
(7,1,1,1),(6,2,1,1),(5,3,1,1),(5,2,2,1),\\
(4,4,1,1),(4,3,2,1),(3,3,3,1).
\end{gathered}}
\tag{5.4}
\]

Thus these seven complete cross-ratio families, not only their harmonic
members, are empty.

For \((4,2,2,2)\), the \(A\)-chart in (5.3) reconstructs exactly to the unit
ideal.  Here \(A\) is the value of the residual binary quadratic at one of
the three double roots.  Permuting the three double roots preserves the
partition and only changes the cross-ratio parameter.  A nonzero binary
quadratic cannot vanish at all three distinct double roots.  Hence one of
these three normalized \(A\)-charts is active, and the single chart theorem
closes the complete partition \((4,2,2,2)\).

The last row \((3,3,2,2)\) behaves differently because its roots have two
multiplicity types.  Exact function-field reconstruction gives

\[
 I_{3322,A}\mathbf Q(\lambda)[\mathbf u]
 =I_{3322,B}\mathbf Q(\lambda)[\mathbf u]
 =I_{3322,C}\mathbf Q(\lambda)[\mathbf u]=(1).
\tag{5.5}
\]

Therefore the generic cross-ratio fiber is empty and every possible
survivor lies in a finite exceptional denominator locus.  The checker also
verifies the same function-field statement for \((4,2,2,2)\) as a
calibration.  An attempted polynomial-ring marked-double-root certificate
for \((3,3,2,2)\) reached the 1200-second cap.  All three polynomial-ring
charts are unit over \(\mathbf F_{32003}\), and the harmonic fiber is exactly
empty, but neither fact identifies or eliminates the characteristic-zero
exceptional locus.  No complete \((3,3,2,2)\) closure is claimed.

## 6. Result and next action

> **Theorem `HC4NHM15` -- Double-conic normal-layer frontend.**  For
> \(q=xz-y^2\), write a ternary quintic uniquely as
> \(h_5=H_5(f_{10})+qH_3(g_6)+q^2H_1(k_2)\).  Then
> \(q^4\mid\det\operatorname{Hess}(h_5)\) if and only if the four binary
> covariants (3.1)--(3.4) vanish.  No clean target
> \(\det\operatorname{Hess}(h_5)=q^4\ell\), \(\ell\ne0\), exists when
> \(f_{10}\) has at most three distinct roots.  Nor does one exist for a
> four-root decic of harmonic cross-ratio.  In addition, the seven complete
> arbitrary-cross-ratio partitions in (5.4) and the partition
> \((4,2,2,2)\) are empty.  The remaining \((3,3,2,2)\) row is empty over
> \(\mathbf Q(\lambda)\), so only a finite exceptional cross-ratio locus can
> survive.

The next exact target is the exceptional denominator locus of
\((3,3,2,2)\).  It needs a smaller cross-ratio eliminant, a blockwise
normal-layer solve, or a denominator-tracking characteristic-zero
certificate over \(\mathbf Q(\lambda)\).  Only a characteristic-zero
survivor of the four layers should be passed to

\[
 C e=q^2\nabla s_3,
 \qquad
 (\nabla s_3)^{\mathsf T}e=q^2a.
\tag{6.1}
\]

The forms with at least five distinct roots remain a separate higher-moduli
frontier.

The exceptional-locus target stated above was subsequently completed by
`HC4NHM18`: the endpoint normal layers force the endpoint residual-line
coefficients to vanish, and two middle coefficients exclude the remaining
line direction for every cross-ratio.

For the many-root invariant continuation, `HC4NHM21` identifies the exact
clean saturation.  The ideal of the four displayed layers must first be
saturated by the three coefficients of `Phi_2`.  An unsaturated
`Disc(f)^N` certificate is impossible because the lower-Smith quintic
`x^5+z^5` has squarefree restriction `s^10+t^10` and all Hessian layers
zero.  After clean saturation, discriminant membership addresses only the
squarefree open; nullcone containment in the binary-decic invariant ring is
the correct target for all stable decics.  See
[`HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md`](HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md).
