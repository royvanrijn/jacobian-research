# The root-incidence derivative split for \(A_4\)

## 1. Scope

The cone construction has already isolated the determinant
\(4W^2K^3L\), and the subsequent ledger calculations show that a Keller
completion inside that model must alter at least two outputs by
source-dependent masks.  This note changes the starting point.  It returns
to the two-parameter generic \(A_4\) quartic and asks what the derivative
unit itself says in the rank-four root algebra.

The outcome is an exact partial construction:

- \(1/P'(T)\) has a compact four-term representative;
- the square discriminant can be split evenly between two incidence
  coordinates;
- the first new coordinate remains a primitive element, so the generic
  four-sheet field is not lost;
- one irreducible orientation pole survives polynomial pullback over the
  target coefficient ring.

Thus the quotient-algebra part of the proposed construction works.  It is
not yet an affine polynomial chart and does not give a Keller map.

## 2. Compact rank-four inverse

Use temporary coefficient variables \(A,B,C\) and put

\[
\begin{aligned}
P(T)={}&T^4-6ABT^2-8B^3T+B^2(9A^2-12CB),\\
R={}&4A^3B-3A^2C^2-6AB^2C+B^4+4BC^3.
\end{aligned}
\]

Then

\[
\operatorname{Disc}_T(P)=-110592B^8R.                \tag{2.1}
\]

In the free rank-four algebra

\[
\mathcal A=\mathbb Q(A,B,C)[T]/(P),
\]

define

\[
\begin{aligned}
n(T)={}&(B^2-AC)T^3+2B(A^2-BC)T^2\\
 &+B(3A^2C-7AB^2+4BC^2)T\\
 &-6B^2(A^3-2ABC+B^3).
\end{aligned}
\]

Exact reduction modulo \(P\) gives

\[
\boxed{\frac1{P'(T)}=\frac{n(T)}{48B^4R}.}            \tag{2.2}
\]

Equivalently, with

\[
N(T)=-2304B^4n(T),
\]

one has the unreduced discriminant identity

\[
\boxed{P'(T)N(T)\equiv\operatorname{Disc}(P)\pmod P.} \tag{2.3}
\]

The unreduced form is the useful one: it remembers that the denominator is
the square of one orientation equation.

## 3. Two-coordinate derivative-unit split

The polynomial \(N(T)\) has the explicit primitive

\[
\begin{aligned}
I(T)={}&576B^4(AC-B^2)T^4
       -1536B^5(A^2-BC)T^3\\
 &-1152B^5(3A^2C-7AB^2+4BC^2)T^2\\
 &+13824B^6(A^3-2ABC+B^3)T.
\end{aligned}                                       \tag{3.1}
\]

Thus \(I'(T)=N(T)\).  On the oriented chart

\[
\Omega^2=\operatorname{Disc}(P),
\]

introduce one incidence variable \(Q\) and set

\[
X=\frac{I(T)}{\Omega},\qquad Y=\frac Q\Omega.         \tag{3.2}
\]

Since \(\Omega\) is independent of \(T,Q\),

\[
\boxed{
\det\frac{\partial(X,Y)}{\partial(T,Q)}
=\frac{N(T)}{\Omega^2}
=\frac1{P'(T)}
}
\quad\text{in }\mathcal A.                           \tag{3.3}
\]

This is the desired algebraic cancellation: the two coordinates each
carry one copy of the orientation pole, and their combined Jacobian is the
derivative unit.

This split does not collapse the inverse cover.  Reduce \(I\) modulo \(P\)
and form the change-of-basis matrix from

\[
(1,T,T^2,T^3)
\quad\text{to}\quad
(1,I,I^2,I^3).
\]

Its determinant is a nonzero polynomial.  For example, at
\((A,B,C)=(0,1,0)\) one has

\[
I\equiv9216T\pmod P,
\]

so the determinant specializes to \(9216^6\ne0\).  Hence

\[
\mathbb Q(A,B,C)(I)=\mathbb Q(A,B,C)(T)              \tag{3.4}
\]

generically.  The root \(T\), and therefore the original quartic relation,
is rationally recoverable from \(I\).  This proves field recovery, not yet
a literal polynomial inverse equation in the new coordinates.

## 4. The surviving polynomiality obstruction

Now specialize to the Jensen--Ledet--Yui parameters

\[
\begin{aligned}
A_0={}&\alpha^3-\beta^3-9\beta^2-27\beta-54,\\
B_0={}&\alpha^3-3\alpha\beta^2+2\beta^3-9\alpha\beta
       +9\beta^2-27\alpha+27\beta+27,\\
C_0={}&\alpha^3-\beta^3+27.
\end{aligned}
\]

Put

\[
\rho=\beta^2+3\beta+9
\]

and

\[
\begin{aligned}
\sigma={}&2\alpha^3\beta+3\alpha^3
-3\alpha^2\beta^2-9\alpha^2\beta-27\alpha^2\\
&+\beta^4+6\beta^3+27\beta^2+54\beta+81.
\end{aligned}
\]

The orientation equation is the polynomial

\[
\boxed{\Omega=1728\rho B_0^4\sigma,}
\qquad
\operatorname{Disc}(P_{\alpha,\beta})=\Omega^2.      \tag{4.1}
\]

There is substantial cancellation in (2.2).  In fact

\[
R(A_0,B_0,C_0)=-27\rho^2\sigma^2,
\]

and all four coefficients of \(n(T)\) acquire a common factor
\(-3\rho\).  The reduced denominator of \(1/P'\) is consequently

\[
432\rho B_0^4\sigma^2.                               \tag{4.2}
\]

This cancellation is not enough to make (3.2) polynomial.  Modulo \(P\),
the \(T^3\)-coefficient of \(I\) is

\[
-1536B^5(A^2-BC).
\]

After specialization,

\[
A_0^2-B_0C_0=3\rho E,
\]

where

\[
\begin{aligned}
E={}&\alpha^4-\alpha^3\beta-6\alpha^3-\alpha\beta^3
+27\alpha\\
&+\beta^4+6\beta^3+27\beta^2+54\beta+81.
\end{aligned}
\]

Therefore the corresponding coefficient of \(I/\Omega\) is

\[
-\frac83\,\frac{B_0E}{\sigma}.                       \tag{4.3}
\]

The exact gcd calculation

\[
\gcd(\sigma,B_0E)=1
\]

shows that (4.3) has a genuine \(\sigma\)-pole.  Adding a
\(\mathbb Q[\alpha,\beta]\)-constant to the primitive changes only its
constant coefficient, so it cannot remove this obstruction.  The second
coordinate \(Q/\Omega\) retains the full orientation pole unless \(Q\) is
given its own source-dependent mask.

Hence the canonical triangular split proves precisely:

\[
\boxed{
\begin{gathered}
\text{the derivative unit admits a two-coordinate orientation split,}\\
\text{but target-only polynomial pullback fails already on }\sigma=0.
\end{gathered}}
                                                               \tag{4.4}
\]

This does not exclude nontriangular pairs linear in \(Q\), nor does it
exclude source charts in which two independent masks cancel the two poles.

## 5. Comparison with the ordinary root-incidence chart

The derivative-unit split must still be compared with the Jacobian of an
ordinary affine chart on the root incidence.  For the rational
\((U,V)\)-chart of the existing \(A_4\) construction, recall

\[
\alpha=\frac{N_1}{H},\qquad
\beta=\frac{N_2}{H},\qquad
\det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
=\frac{4K^3L}{H^3}.                                  \tag{5.1}
\]

Put

\[
\mathsf A=U^3-V^3-9V^2-27V-54.
\]

There are polynomials \(A_6,C_6\) defined by

\[
A_0(\alpha,\beta)=\frac{K^3A_6}{H^3},\qquad
C_0(\alpha,\beta)=\frac{K^3C_6}{H^3},
\]

while

\[
B_0(\alpha,\beta)=\frac{K^3L^2}{H^3}.
\]

The selected scaled root of \(P_{\alpha,\beta}\) on this chart is

\[
\boxed{
T_*=\frac{3\mathsf A K^3L}{H^3}.
}                                                     \tag{5.2}
\]

The checker verifies (5.2) by the denominator-free identity

\[
81\mathsf A^4-54A_6\mathsf A^2-24\mathsf A L^3
+9A_6^2-12C_6L^2=0.
\]

Define

\[
\Theta=
\frac{27\mathsf A^3-9A_6\mathsf A-2L^3}{2L}.
\]

Exact division shows that \(\Theta\in\mathbb Q[U,V]\), and direct
differentiation gives

\[
\boxed{
P'_{\alpha,\beta}(T_*)
=\frac{8K^9L^4\Theta}{H^9}.
}                                                     \tag{5.3}
\]

The square-discriminant orientation also pulls back completely.  With

\[
\rho_s=V^2+3V+9,\qquad
\sigma_s=\sigma(U,V),
\]

one has

\[
\boxed{
\Omega(\alpha,\beta)
=\frac{1728\rho_s\sigma_s\Theta K^{18}L^8}{H^{18}}.
}                                                     \tag{5.4}
\]

In particular,

\[
\frac{\Omega(\alpha,\beta)}
     {P'_{\alpha,\beta}(T_*)}
=\frac{216\rho_s\sigma_sK^9L^4}{H^9}.                \tag{5.5}
\]

This exposes the missing compatibility.  The ordinary suspension

\[
(U,V,Q)\longmapsto
\left(\alpha,\beta,\frac{Q}{P'_{\alpha,\beta}(T_*)}\right)
\]

has Jacobian

\[
\boxed{
\frac{4K^3L}{H^3}\,
\frac1{P'_{\alpha,\beta}(T_*)}
=\frac{H^6}{2\Theta K^6L^3},
}                                                     \tag{5.6}
\]

which is not constant.  The reciprocal unit actually required by the
ordinary \((U,V)\)-chart is

\[
\left(
\det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
\right)^{-1}
=\frac{H^3}{4K^3L}.                                  \tag{5.7}
\]

Thus the current root-incidence parametrization returns exactly to the old
\(K^3L/H^3\) cone ledger.  The rank-four derivative split is not a
standalone escape from that ledger: a successful construction must also
replace the affine root-incidence chart or realize (5.7) through two
source-dependent polynomial masks.

## 6. The correct chart unit in the root basis

The comparison in Section 5 identifies the function that must actually be
realized:

\[
\mathcal J^{-1}
=\frac{H^3}{4K^3L}.
\]

Four-branch interpolation in the root basis gives an exact formula

\[
\boxed{
\mathcal J^{-1}
=\frac{
3B_0p_0+B_0p_1T+p_2T^2-p_3T^3
}{
72B_0^2\rho\sigma
}.
}                                                     \tag{6.1}
\]

Here

\[
\deg(p_0,p_1,p_2,p_3)=(9,6,6,3);
\]

their explicit integer coefficients are recorded in
[`verify_a4_chart_unit_rank_four.py`](../scripts/verify_a4_chart_unit_rank_four.py).
The checker derives the four conjugate values independently on the
\((s,t)\)-presentation and verifies every coefficient of (6.1) by exact
polynomial identities.

The denominator is structured:

\[
\boxed{B_0^2\rho\sigma.}                              \tag{6.2}
\]

Moreover \(B_0,\rho,\sigma\) are irreducible and pairwise coprime over
\(\mathbb Q[\alpha,\beta]\).  Thus (6.2) has three genuine target boundary
components, but it admits the formal two-mask grouping

\[
B_0^2\mid\rho\sigma.
\]

Write

\[
\mathcal N(T)=3B_0p_0+B_0p_1T+p_2T^2-p_3T^3.
\]

On the localized fourfold chart, introduce two mask variables \(Q_1,Q_2\)
and use

\[
\boxed{
X=\frac{Q_1}{B_0^2},\qquad
Y=\frac{\mathcal N(T)Q_2}{72\rho\sigma}.
}                                                     \tag{6.3}
\]

Then

\[
\det\frac{\partial(X,Y)}{\partial(Q_1,Q_2)}
=\mathcal J^{-1}.                                    \tag{6.4}
\]

Consequently the rational map

\[
(U,V,Q_1,Q_2)\longmapsto(\alpha,\beta,X,Y)
\]

has determinant one and retains the original degree-four \(A_4\) function
field: the two mask variables are recovered uniquely at the generic point.
This is the exact localized two-mask Keller suspension sought by the
rank-four calculation.

It is still not polynomial.  The first mask has a genuine \(B_0^2\)-pole.
The coefficient \(p_3\) is coprime to both \(\rho\) and \(\sigma\), so the
second mask also has genuine poles along both remaining components.
Equivalently, integrating \(\mathcal N(T)\) to make a triangular incidence
coordinate leaves a \(T^3\)-coefficient \(p_2/3\), with

\[
\gcd(p_2,B_0)=1.
\]

No target-only redistribution removes these poles.

### 6.1. Quick local viability screen

The reduced boundary arrangement is not simple normal crossings.  Exact
Gröbner calculations give

\[
\operatorname{Sing}(\sigma)=V(a^2,\rho)
\]

and, for every pair chosen from \(B_0,\rho,\sigma\), the nontransverse
intersection scheme is again

\[
V(a^2,\rho).
\]

All three components meet, with triple-intersection ideal

\[
(B_0,\rho,\sigma)=(a^3,\rho).                        \tag{6.5}
\]

This scheme consists of the two conjugate points

\[
a=0,\qquad \beta^2+3\beta+9=0.
\]

The local structure is explicit.  Put

\[
z=\rho,\qquad c=2\beta+3.
\]

Since \(c^2+27=4z\), the function \(c\) is a unit along \(z=0\).
Reduction by \(\beta^2+3\beta+9-z\) gives

\[
\begin{aligned}
B_0&=a^3-3az+cz,\\
\sigma&=ca^3-3a^2z+z^2.                              \tag{6.6}
\end{aligned}
\]

The coefficient \(p_3\) controlling the \(T^3\)-term of
\(\mathcal N(T)\) has local form

\[
p_3=4ca^2-16az+63a+4cz.                              \tag{6.7}
\]

At the common cluster \((a,z)=(0,0)\),

\[
\operatorname{ord}(B_0^2)=2,\qquad
\operatorname{ord}(\rho\sigma)=3,\qquad
\operatorname{ord}(p_3)=1.                           \tag{6.8}
\]

Thus the numerator supplies only one of the three local orders required by
the second mask.  The simplest two-normal-crossings reciprocal chart fails
this viability test.

This is not a theorem excluding every resolved affine modification.  Such a
modification could blow up the common cusp/tangency cluster and distribute
the resulting exceptional valuations.  It does show that the apparent
two-factor grouping in (6.2) hides a higher-contact three-component
boundary package; realizing it polynomially would require a nontrivial
resolution-aware construction.

This is the final conclusion of the present construction:

\[
\boxed{
\begin{gathered}
\text{the correct unit has an exact two-mask localized realization,}\\
\text{but all three target boundary components remain genuine.}
\end{gathered}}
                                                               \tag{6.9}
\]

A polynomial counterexample would therefore require an affine source
modification in which the two source masks vanish with divisors
\(B_0^2\) and \(\rho\sigma\), respectively.  Constructing that modification
is additional geometry, not another quotient-ring syzygy.

## 7. Structured next search

The next calculation should no longer vary arbitrary shears of the cone.
Work in the basis \(1,T,T^2,T^3\) of

\[
\mathbb Q[\alpha,\beta,T]/(P_{\alpha,\beta})
\]

and write two general incidence coordinates

\[
\begin{aligned}
X&=x_0(\alpha,\beta,T)+Qx_1(\alpha,\beta,T),\\
Y&=y_0(\alpha,\beta,T)+Qy_1(\alpha,\beta,T).
\end{aligned}
\]

Reduction modulo \(P\) turns every condition into four coefficient
equations.  A bounded search must impose simultaneously:

1. the relative Jacobian equals the derivative unit;
2. two specified source masks make \(X,Y\) polynomial;
3. the power-basis determinant is nonzero, so \(T\) is recoverable;
4. elimination recovers the intended \(P_{\alpha,\beta}(T)\), not merely
   an unrelated quartic primitive.

The canonical pair (3.2) is a base point for that system.  Equation (4.3)
identifies the first divisor that a nontriangular deformation must move.
This is a finite free-algebra syzygy problem; it is not another search over
ambient cone automorphisms.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_root_incidence_derivative_split.py
.venv/bin/python scripts/verify_a4_chart_unit_rank_four.py
.venv/bin/python scripts/verify_a4_two_mask_local_viability.py
```

The checker verifies the compact inverse basis, the two-coordinate
Jacobian identity, generic root-field recovery, the square-discriminant
specialization, the residual \(\sigma\)-pole, the selected rational root,
and the complete comparison with the ordinary \((U,V)\)-chart ledger.  The
second checker verifies the exact rank-four expansion of the correct chart
unit, its three pairwise-coprime boundary factors, and the localized
two-mask determinant-one suspension.  The third checker computes the common
nontransverse cluster, its local normal forms, and the resulting order
deficit for the simple two-mask chart.
