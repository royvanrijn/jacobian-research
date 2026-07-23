# Degree-six Ritt boundary atlas

This note completes the boundary calculation left implicit in the
[degree-six Ritt atlas](DEGREE_SIX_RITT_ATLAS.md).  It separates three notions:

1. failure of the exact-double zero cluster;
2. collision of an additional primitive root; and
3. collision of the Hessian divisor.

The first two form the affine-sheet reconstruction boundary.  The third is
the deleted D1 Hessian divisor.  The distinguished root-one sheet itself
never ramifies, because every normalized seed satisfies `H'(1)=-1`.

## 1. Coefficient boundary and primary components

Use normalized coordinates

\[
 h_4=x,\qquad h_5=y,\qquad h_6=z,
\]

so that

\[
 h_3=-1-2x-3y-4z,qquad h_2=1+x+2y+3z.             \tag{1.1}
\]

Let `Phi_23` and `Phi_32` be the two Ritt equations from the main atlas.  Put

\[
 P=H/W^2,qquad
 \Delta_{\rm aff}=\operatorname{Res}_W(P,2P+WP').   \tag{1.2}
\]

Exact calculation gives

\[
 \Delta_{\rm aff}=z\,h_2\,Q_4(x,y,z),               \tag{1.3}
\]

where `z=0` is degree drop, `h_2=0` is the zero-cluster boundary, and
`Q_4=0` is the additional-multiple-root boundary.  The quartic is

\[
\begin{aligned}
Q_4={}&-4x^3z+x^2y^2-28x^2yz-56x^2z^2+6xy^3-48xy^2z\\
&-240xyz^2-18xyz-240xz^3-72xz^2+9y^4+4y^3\\
&-180y^2z^2-6y^2z-432yz^3-132yz^2-288z^4-176z^3-27z^2.
                                                               \tag{1.4}
\end{aligned}
\]

After saturation by `z`, exact primary decomposition gives

\[
 \mathcal D_{23}\cap V(\Delta_{\rm aff})
 =V(\Phi_{23},h_2)\ \cup\ V(\Phi_{23},Q_4),          \tag{1.5}
\]

with two irreducible curves, and

\[
 \mathcal D_{32}\cap V(\Delta_{\rm aff})
 =V(\Phi_{32},h_2)\ \cup\ C_X\ \cup\ C_Z,          \tag{1.6}
\]

with three irreducible curves.  On the dense `3 o 2` chart, `C_X` and `C_Z`
are the factors

\[
 X=ap^3+4ap^2+5ap+2a+1=0,                          \tag{1.7}
\]

and

\[
\begin{aligned}
Z={}&25a^2p^6+155a^2p^5+379a^2p^4+457a^2p^3
 +272a^2p^2+64a^2p\\
&+23ap^3+47ap^2+24ap-2=0.                           \tag{1.8}
\end{aligned}
\]

The zero-cluster factor on this chart is

\[
 Y=5ap^3+12ap^2+9ap+2a+2p+1=2(p+1)h_2.            \tag{1.9}
\]

Thus the five curves in (1.5)--(1.6) exhaust the exact-degree affine-sheet
boundary of the two Ritt surfaces.

## 2. Rational witnesses for all five curves

Every component has a rational Hessian-clean, weighted-admissible witness:

\[
\begin{array}{c|l|c}
\text{component}&H(W)&\text{degeneration}\\ \hline
D_{23}^{0}
&-\frac13W^3(W-1)(W^2+W+1)&\operatorname{ord}_0H=3\\[2pt]
D_{23}^{+}
&-\frac1{1024}W^2(W-1)(9W-17)^2(9W+7)&W=17/9\text{ double}\\[2pt]
D_{32}^{X}
&\frac1{100}W^2(W-6)^2(W-5)(W-1)&W=6\text{ double}\\[2pt]
D_{32}^{0}
&-\frac1{700}W^3(W-1)(11W^2-187W+876)&\operatorname{ord}_0H=3\\[2pt]
D_{32}^{Z}
&\frac1{147}W^2(W-1)(2W-5)(2W+5)^2&W=-5/2\text{ double}.
\end{array}                                                       \tag{2.1}
\]

Their decomposition parameters are

\[
\begin{array}{c|c}
D_{23}^{0}&(a,b,p,q)=(-1/3,1/3,0,0),\\
D_{23}^{+}&(a,b,p,q)=(-729/1024,-1,-2,-5/27),\\
D_{32}^{X}&(a,b,c,p)=(1/100,1/20,0,-6),\\
D_{32}^{0}&(a,b,c,p)=(-11/700,5/28,-45/7,-6),\\
D_{32}^{Z}&(a,b,c,p)=(8/147,-76/147,48/49,1/2).
\end{array}                                                       \tag{2.2}
\]

The first two witnesses remain on the all-double omitted component.  Their
omitted factorizations are

\[
 H_{23}^{0}-{1\over12}
 =-{1\over3}\left(W^3-{1\over2}\right)^2,          \tag{2.3}
\]

and, with `B=W^3-2W^2-5W/27`,

\[
 H_{23}^{+}+{5\over27}W-{256\over729}
 =-{729\over1024}\left(B+{512\over729}\right)^2.  \tag{2.4}
\]

The three `3 o 2`-only witnesses have no omitted value: in degree six every
omitted type lies in the all-double surface or all-triple curve, and these
witnesses avoid both.  They therefore give explicit surjective seeds with an
imprimitive vertical specialization on the affine-sheet boundary.

## 3. Global chart for the common Ritt curve

The centered common composition is

\[
 g_{c,q}(W)=\big((W-c)^3+q(W-c)\big)^2.              \tag{3.1}
\]

The normalized endpoint condition is

\[
\begin{aligned}
 E(c,q)={}&15c^4-20c^3+12c^2q+15c^2-8cq-6c\\
          &+q^2+2q+1=0.                              \tag{3.2}
\end{aligned}
\]

Writing

\[
 D=g'_{c,q}(1)-g'_{c,q}(0),
\]

the normalized seed is

\[
 H_{c,q}(W)=-{g_{c,q}(W)-g_{c,q}(0)-g'_{c,q}(0)W\over D}. \tag{3.3}
\]

The residual primitive polynomial

\[
 Q_{c,q}(W)={g_{c,q}(W)-g_{c,q}(0)-g'_{c,q}(0)W\over W^2}
\]

satisfies

\[
 Q_{c,q}(0)=15c^4+12c^2q+q^2,                       \tag{3.4}
\]

and

\[
\begin{aligned}
\operatorname{Disc}_W(Q_{c,q})={}&16c^2(c^2+q)^2(3c^2+q)\\
&\cdot(225c^4+285c^2q+64q^2).                       \tag{3.5}
\end{aligned}
\]

Finally,

\[
 \operatorname{Disc}_W(g''_{c,q})=108380160q^6.     \tag{3.6}
\]

Equations (3.4)--(3.6) separate the zero-cluster, extra-root, and Hessian
boundaries without any full coefficient elimination.

## 4. The six zero-cluster points

Eliminating `q` from (3.2) and (3.4) gives

\[
\boxed{
 Z_6(c)=560c^6-840c^5+411c^4-20c^3-42c^2+12c-1=0.} \tag{4.1}
\]

It is irreducible over `Q`, squarefree, and has discriminant

\[
 \operatorname{Disc}(Z_6)=83901850583040.            \tag{4.2}
\]

The corresponding coordinate is

\[
 q={1120c^5-1400c^4+472c^3+73c^2-62c+7\over2}.     \tag{4.3}
\]

All six points have `qD!=0`, are Hessian-clean and weighted-admissible, and
form one reduced degree-six orbit over `Q`.  In coefficient space,

\[
 (\Phi_{23},\Phi_{32},h_2):z^\infty
\]

has vector-space length six and equals its radical.

## 5. The six extra-root points

After deleting the normalization-pole solutions of (3.5), the exact extra
boundary consists of two rational points

\[
 (c,q)=(-1,-3),\qquad (c,q)=(1/3,-1/3),             \tag{5.1}
\]

and one quartic orbit

\[
\boxed{
 R_4(c)=180c^4-300c^3-49c^2+208c-64=0,}            \tag{5.2}
\]

with

\[
 q={60c^3-80c^2-43c+28\over12}.                    \tag{5.3}
\]

The quartic is irreducible and squarefree, with

\[
 \operatorname{Disc}(R_4)=13168189440000.           \tag{5.4}
\]

The rational seeds are

\[
 -{1\over36}W^2(W-1)(W+2)^2(W+3)                  \tag{5.5}
\]

and

\[
 -{1\over4}W^2(W-1)(3W-2)^2(3W+1).                \tag{5.6}
\]

The saturated coefficient ideal

\[
 (\Phi_{23},\Phi_{32},Q_4):z^\infty
\]

is reduced of total length six and splits over `Q` as `1+1+4`.

The discarded factors are genuine chart poles, not missing seeds.  For
example, `c=0,q=-1` and the two solutions of `q=-c^2` on (3.2) all satisfy
`D=0`.  The double solution `c=1/2,q=-3/4` does as well.

## 6. The four Hessian collisions

By (3.6), the Hessian divisor collides exactly at `q=0`.  Equation (3.2)
then becomes

\[
\boxed{
 T_4(c)=15c^4-20c^3+15c^2-6c+1=0.}                \tag{6.1}
\]

This polynomial is squarefree with discriminant `10800`.  These four points
are precisely the type-`(6)` support where the all-double and all-triple
omitted components meet.  Their omitted-component intersection length is
two, as independently checked in the Gaussian moment chart.

Importantly, `D`, weighted admissibility, and `Delta_aff` are all units at
the four points.  Thus they are Hessian collisions but not affine-sheet
boundary points.  This proves that the two boundary mechanisms must remain
separate in any higher-degree atlas.

## 7. Resulting boundary diagram

On the normalized common Ritt curve, the finite deleted support is therefore

\[
 \underbrace{6}_{\text{zero cluster}}
 \;\sqcup\;
 \underbrace{(1+1+4)}_{\text{extra root}}
 \;\sqcup\;
 \underbrace{4}_{\text{Hessian/type-(6)}}.          \tag{7.1}
\]

Every support is reduced as a boundary cut.  Only the intersection of the
two omitted-component schemes is doubled at the last four points.

The exact certificates are

```bash
.venv/bin/python scripts/verify_degree_six_ritt_boundary_atlas.py
Singular -q scripts/verify_degree_six_ritt_boundary_atlas.sing
```

The Python checker verifies the factorizations, fields, witnesses, clean-open
conditions, and excluded normalization poles.  The Singular checker performs
the exact-degree saturation, primary decompositions, component counts, and
local length calculations.

