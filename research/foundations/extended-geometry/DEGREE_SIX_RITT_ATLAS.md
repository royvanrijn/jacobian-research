# Degree-six Ritt atlas

Work over an algebraically closed field of characteristic zero.  This note
refines the degree-six Hessian-incidence solution in
[HESSIAN_RITT_INCIDENCE.md](HESSIAN_RITT_INCIDENCE.md).  It computes the
vertical decomposition loci only after imposing the normalized
Hessian-divisor constraints, then gives the requested comparison with the
omitted-value and affine-sheet boundary strata.  The exceptional-component
intersection multiplicities agree with the independent moment-coordinate
calculation in
[GAUSSIAN_EXCEPTIONAL_MOMENT_GEOMETRY.md](GAUSSIAN_EXCEPTIONAL_MOMENT_GEOMETRY.md).

## 1. Hessian-first normalized sextics

Write

\[
 H(W)=h_2W^2+h_3W^3+h_4W^4+h_5W^5+h_6W^6,
 \qquad h_6\ne0.
\]

Starting instead with the quartic Hessian

\[
 K(W)=k_0+k_1W+k_2W^2+k_3W^3+k_4W^4
\]

and integrating twice with zero constant and linear terms gives

\[
 (k_0,k_1,k_2,k_3,k_4)
 =(2h_2,6h_3,12h_4,20h_5,30h_6).
\]

The conditions `H(1)=0` and `H'(1)=-1` are

\[
 h_2+h_3+h_4+h_5+h_6=0,
 \qquad
 2h_2+3h_3+4h_4+5h_5+6h_6=-1.
\]

Thus

\[
 h_2=-(h_3+h_4+h_5+h_6),
 \qquad
 L:=h_3+2h_4+3h_5+4h_6+1=0.                 \tag{1.1}
\]

All eliminations below take place after (1.1).  The admissible open further
removes `H''(1)=-2`; the Hessian-clean open removes
`Disc(H'')=0`.  Neither removal changes the closed equations below.

## 2. The two decomposition surfaces

For `2 o 3`, normalize the inner cubic by

\[
 B_3=W^3+pW^2+qW,
 \qquad f=A(B_3)=aB_3^2+bB_3.
\]

Subtracting the linear term of `f` produces `H`.  Matching `H''=f''` and then
using (1.1) gives the hypersurface

\[
\begin{aligned}
 \Phi_{23}={}&32h_3h_5h_6^2+64h_3h_6^3+16h_4^2h_6^2
 -24h_4h_5^2h_6+64h_4h_6^3\\
 &+5h_5^4+64h_5h_6^3+64h_6^4=0.             \tag{2.1}
\end{aligned}
\]

For `3 o 2`, use

\[
 B_2=W^2+pW,
 \qquad f=A(B_2)=aB_2^3+bB_2^2+cB_2.
\]

The corresponding hypersurface is

\[
 \Phi_{32}=27h_3h_6^2-18h_4h_5h_6+5h_5^3=0. \tag{2.2}
\]

These are exact on `D(h_6)`, not merely necessary equations.  Rational
reconstruction is

\[
\begin{array}{c|c}
2\circ3 &
p={h_5\over2h_6},\quad
q={1\over2}\left({h_4\over h_6}-p^2\right),\quad
b=h_3-2h_6pq,\quad a=h_6,\\[4pt]
3\circ2 &
p={h_5\over3h_6},\quad
b=h_4-3h_6p^2,\quad
c=h_2-bp^2,\quad a=h_6.
\end{array}                                                   \tag{2.3}
\]

Consequently the saturated coefficient ideals are

\[
 I_{23}=(L,\Phi_{23}):h_6^\infty,
 \qquad
 I_{32}=(L,\Phi_{32}):h_6^\infty.                 \tag{2.4}
\]

A dense normalization chart of `D_23` is

\[
 a=-{p+1\over(p+q+1)D},\qquad
 b={p+q+1\over D},\qquad
 D=2p^2+5p-q+3.                                    \tag{2.5}
\]

A dense chart of `D_32` is

\[
 b=-{ap^3+6ap^2+9ap+4a+1\over2(p+1)},              \tag{2.6}
\]

\[
 c={ap^4+5ap^3+9ap^2+7ap+2a+p+1\over2}.            \tag{2.7}
\]

## 3. Their Ritt intersection

The complete coefficient-space intersection is

\[
 (L,\Phi_{23},\Phi_{32}):h_6^\infty.               \tag{3.1}
\]

On the dense chart (2.5), and away from its displayed denominators and
`p=-1`, it is the plane curve

\[
\begin{aligned}
 R(p,q)={}&4p^4+4p^3-18p^2q-27p^2-72pq-54p\\
          &-27q^2-54q-27=0.                         \tag{3.2}
\end{aligned}
\]

The coefficient ideal (3.1), rather than (3.2), retains the exceptional
chart points automatically.

The intersection is not forced onto any bad Hessian or affine-sheet
boundary.  A clean rational point is `p=-3,q=-2`, giving

\[
 \boxed{H=-{1\over16}W^2(W-1)(W^3-5W^2+20).}       \tag{3.3}
\]

It has

\[
 \operatorname{Disc}(H'')={103359375\over1024}\ne0,
 \qquad
 \gcd(H/W^2,H'/W)=1,
 \qquad H''(1)=-{25\over8}.
\]

Its two decompositions use

\[
 (a,b,p,q)_{23}=(-1/16,-1/2,-3,-2)
\]

and

\[
 (a,b,c,p)_{32}=(-1/16,7/16,-1/2,-2).
\]

## 4. Omitted-value components

The relation with omitted values is sharper than a mere intersection count.
Completing the outer quadratic gives

\[
 H-sW+t
 =a\left(B_3+{b\over2a}\right)^2,
 \qquad
 s=-bq,quad t={b^2\over4a}.                        \tag{4.1}
\]

Therefore

\[
 \boxed{\overline{\mathcal O_{222}}=\mathcal D_{23}.} \tag{4.2}
\]

The exact `(2,2,2)` stratum is the open where the cubic in (4.1) is
squarefree; its cubic-discriminant boundary contains types `(4,2)` and `(6)`.

For a cubic outer polynomial, adding a constant produces a cube precisely
when `3ac=b^2`.  In seed coordinates this is

\[
 \Psi_{33}:=9h_2h_6^2-3h_4^2h_6+h_4h_5^2=0.       \tag{4.3}
\]

Thus

\[
 \overline{\mathcal O_{33}}
 =V(L,\Phi_{32},\Psi_{33})\cap D(h_6).              \tag{4.4}
\]

When (4.3) holds,

\[
 H-sW+t
 =a\left(B_2+{b\over3a}\right)^3,
 \qquad
 s=-cp,quad t={b^3\over27a^2}.                    \tag{4.5}
\]

The exact `(3,3)` stratum is the open where the quadratic in (4.5) has two
distinct roots.  One rational clean witness is

\[
 \boxed{H={1\over27}W^2(W-1)(8W^3+20W^2-10W-45),} \tag{4.6}
\]

for which

\[
 H+W-1={8\over27}\left(W^2+{W\over2}-{3\over2}\right)^3.
\]

The cross-intersections are now explicit:

\[
 \mathcal D_{23}\cap\overline{\mathcal O_{222}}=\mathcal D_{23},
 \qquad
 \mathcal D_{32}\cap\overline{\mathcal O_{222}}
 =\mathcal D_{23}\cap\mathcal D_{32},              \tag{4.7}
\]

and `D_32 intersect O_33=O_33`.  Finally,

\[
 \mathcal D_{23}\cap\overline{\mathcal O_{33}}
\]

consists of type-`(6)` degenerations.  Write

\[
 H_\rho(W)=-{(W-\rho)^6-\rho^6+6\rho^5W
 \over6(5\rho^4-10\rho^3+10\rho^2-5\rho+1)}.       \tag{4.8}
\]

The four support points are cut out by

\[
 15\rho^4-20\rho^3+15\rho^2-6\rho+1=0.             \tag{4.9}
\]

The quartic has discriminant `10800`; its normalization and admissibility
denominators are coprime to it.  Thus the reduced support consists of four
admissible algebraic seeds.  Scheme-theoretically each point has intersection
length two: on the all-triple normalization, the all-double equation pulls
back as the square of this collision divisor.  In the `p` coordinate of
(2.5), equation (4.9) is

\[
 5p^4+20p^3+45p^2+54p+27=0,
 \qquad p=-3\rho.                                   \tag{4.10}
\]

There is no exact seed lying simultaneously in the open `(2,2,2)` and
`(3,3)` strata: uniqueness of the omitted value forces the common polynomial
to be both a square and a cube, hence a sixth power.

## 5. Degeneration of the affine-sheet distinction

The root-one sheet itself never ramifies on the normalized slice, since
`H'(1)=-1`.  What can degenerate is the exact-double, boundary-clean
distinction used to recognize it in the regular-reconstruction open.  Put

\[
 P=H/W^2.
\]

Its exact coefficient equation is

\[
 \boxed{\Delta_{\rm aff}
 =\operatorname{Res}_W(P,2P+WP')=0.}                \tag{5.1}
\]

The factor `P(0)=h_2` is the zero-cluster degeneration; the remaining factor
detects a nonzero common root of `H` and `H'`.  On the normalized exact-degree
slice, (5.1) factors as `h_6 h_2` times a quartic coefficient polynomial, up
to a nonzero scalar.

On the dense `2 o 3` chart, after deleting the displayed units and
denominators, the pullback is

\[
 \Delta_{23}=U_{23}V_{23},                           \tag{5.2}
\]

where

\[
 U_{23}=p^3+2p^2q+2p^2+2pq+p-q^2,                  \tag{5.3}
\]

\[
 V_{23}=4p^2q+p^2-10pq-2p-27q^2-14q-3.             \tag{5.4}
\]

Here `U_23=0` is the zero-cluster factor and `V_23=0` is the additional
multiple-root factor.

On the dense `3 o 2` chart the pullback, again up to exact-degree units and
denominators, is

\[
 \Delta_{32}=X_{32}Y_{32}Z_{32},                    \tag{5.5}
\]

with

\[
 X_{32}=ap^3+4ap^2+5ap+2a+1,                       \tag{5.6}
\]

\[
 Y_{32}=5ap^3+12ap^2+9ap+2a+2p+1,                  \tag{5.7}
\]

and

\[
\begin{aligned}
 Z_{32}={}&25a^2p^6+155a^2p^5+379a^2p^4+457a^2p^3
 +272a^2p^2+64a^2p\\
 &+23ap^3+47ap^2+24ap-2.                            \tag{5.8}
\end{aligned}
\]

The middle factor is `2(p+1)h_2`; the other two detect nonzero multiple
primitive roots.  Intersections with any atlas stratum are obtained simply
by adjoining the appropriate factors from (5.2) or (5.5) to its ideal.  In
particular the clean point (3.3) proves that the common Ritt curve is not
contained in this degeneration locus.

For completeness, the affine-sheet boundary on the common Ritt curve (3.2)
has the following projected equations on the same dense chart.  The
zero-cluster cut satisfies

\[
\begin{aligned}
\operatorname{Res}_q(R,U_{23})=-(p+1)^2(&560p^6+2520p^5+3699p^4
 +540p^3\\
 &-3402p^2-2916p-729),                              \tag{5.9}
\end{aligned}
\]

while the additional-multiple-root cut satisfies

\[
\begin{aligned}
\operatorname{Res}_q(R,V_{23})={}&108(p-3)(p+1)^3\\
 &\cdot(20p^4+100p^3-49p^2-624p-576).              \tag{5.10}
\end{aligned}
\]

The factor `p+1` belongs to the exceptional parameter chart and must be
resolved using the coefficient ideal (3.1); away from it, (5.9)--(5.10)
enumerate the finite boundary cuts exactly.

## 6. Concrete clean witnesses

The following four rational seeds lie on the Hessian-clean,
boundary-clean, admissible open and separate the main strata:

\[
\begin{array}{c|l}
\mathcal D_{23}\setminus\mathcal D_{32}
&{1\over27}W^2(-2W^4-4W^3-6W^2+5W+7),\\[2pt]
\mathcal D_{32}\setminus\mathcal D_{23}
&{1\over2}W^2(W-1)(W+1)(2W^2-3),\\[2pt]
\mathcal D_{23}\cap\mathcal D_{32}
&-{1\over16}W^2(W-1)(W^3-5W^2+20),\\[2pt]
\mathcal O_{33}\setminus\mathcal D_{23}
&{1\over27}W^2(W-1)(8W^3+20W^2-10W-45).
\end{array}                                         \tag{6.1}
\]

The executable certificate checks all coefficient identities,
reconstructions, omitted factorizations, Hessian discriminants, boundary
resultants, and witness exclusions:

```bash
.venv/bin/python scripts/verify_degree_six_ritt_atlas.py
```

This completes the requested degree-six atlas.  `OP-RITT` remains open in
general composite degree.

The exact primary decomposition, rational component witnesses, and finite
boundary cuts on the common Ritt curve are computed in the companion
[degree-six Ritt boundary atlas](DEGREE_SIX_RITT_BOUNDARY_ATLAS.md).
