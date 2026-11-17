# Complete affine singularity atlas for the F2 `k=1` target

> **Status.** Exact algebraic classification.  Every normalized degree
> `(3,5)` target curve
> `p=t^3+a*t`, `q=t^5+b*t^4+c*t^2+d*t` has exactly one of the twenty-one
> affine ADE packets in Tables 1 and 2 below.  Off the merger hypersurface
> `M=a^2+a*b^2-b*c+d=0`, the packet is one of twelve `A`-type rows.  On
> `M=0`, the collision quartic factors as a line times a depressed cubic,
> and its factor multiplicities give exactly nine `D/E` rows.  In
> particular the previously omitted generic `D4+A2` merger is a real
> stratum, while `A2+A5`, `2A2+A3`, `2A4`, and `A8` cannot occur as
> distinct-image packets.  This is a target-curve classification, not by
> itself a Keller-map exclusion.

The polynomial identities and a rational witness for every row are replayed
by
[`verify_f2_affine_k1_complete_singularity_atlas.py`](../scripts/verify_f2_affine_k1_complete_singularity_atlas.py).

## 1. Collision roots and the merger hypersurface

Write

\[
 p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt.            \tag{1.1}
\]

The collision theorem associates to an unordered normalization pair one
root `u=s+t` of

\[
 R(u)=u^4+bu^3+au^2+(2ab-c)u-(a^2+d).          \tag{1.2}
\]

The pair is diagonal precisely when

\[
 D(u)=3u^2+4a=0,                                \tag{1.3}
\]

and its target value is

\[
 x(u)=-u(u^2+a),\qquad
 y(u)=(u^2+a)(u^3+2au+ab-c).                   \tag{1.4}
\]

There is one and only one mechanism by which different collision roots can
have the same target image.

### Proposition 1.1 -- exact image-merger equation

If two distinct roots of `R` have the same pair `(x,y)`, then

\[
 \boxed{M=a^2+ab^2-bc+d=0}.                     \tag{1.5}
\]

Conversely, `M=0` is the closure of the equal-image locus, and on it

\[
 \boxed{
 R(u)=(u+b)C(u),\qquad
 C(u)=u^3+au+ab-c.}                             \tag{1.6}
\]

Every root of `C`, counted with multiplicity, maps to the single target point

\[
 (x,y)=(ab-c,-a(ab-c)).                         \tag{1.7}
\]

#### Proof

For distinct `u,v`, equality of the first coordinates in (1.4) gives

\[
 a=-(u^2+uv+v^2).                               \tag{1.8}
\]

After this substitution, `(R(u)-R(v))/(u-v)=0` gives

\[
 c=-b(u^2+uv+v^2)-u^2v-uv^2.                  \tag{1.9}
\]

The divided difference of `y` is `-(u+v)` times the left side of (1.9),
so the second target coordinate then agrees automatically.  Solving
`R(u)=0` for `d` and substituting (1.8)--(1.9) gives `M=0` identically.
Conversely, (1.6) is direct division.  Modulo `C`, equations (1.4) reduce
to (1.7).  This proves both directions. \(\square\)

Whenever `C` has at least two distinct roots this is a genuine merger of
distinct normalization pairs.  When its roots coalesce, (1.5) retains the
`D/E` limit strata.  Thus no additional equal-image discriminant is missing
from the atlas.

## 2. The twelve nonmerger packets

Assume first `M!=0`.  Different collision roots have different target
values.  An off-diagonal root of multiplicity `m` gives two smooth branches
of intersection multiplicity `m`, hence `A_(2m-1)`.  A diagonal root at a
nonzero critical parameter gives, after an analytic target shear, the
parametrization `(z^2,z^(2m+1))`, hence `A_(2m)`.  At the zero critical
parameter the only exceptional possibilities are the `E6/E8` rows of
Section 3, and both lie on `M=0`.

A nonzero common critical point can be scaled to `t=1`.  Then

\[
 a=-3,\qquad d=-5-4b-2c,                       \tag{2.1}
\]

and

\[
 R=(u-2)B(u),\quad
 B=u^3+(b+2)u^2+(2b+1)u+2-2b-c.                \tag{2.2}
\]

The controlling factorizations are

\[
\begin{aligned}
 \operatorname{Disc}(B)
  &=(3b+c-2)(4b^3-12b^2-69b-27c+50),\\
 M&=-(b+2)(3b+c-2),                             \tag{2.3}\\
 B(2)&=6b-c+20,\qquad B(-2)=-2b-c.
\end{aligned}
\]

If `B(2)=0`, then

\[
 B=(u-2)(u^2+(b+4)u+4b+9),\qquad
 \operatorname{Disc}_{u}(B/(u-2))=(b-10)(b+2). \tag{2.4}
\]

The marked diagonal root has multiplicity three only at `b=-7/2`; it can
never have multiplicity four.  The other diagonal root `u=-2` is common
when `c=-2b`, and can be double only at `b=5/2`.  On the two-critical locus

\[
 R=(u^2-4)(u^2+bu+1),\qquad M=4-b^2.           \tag{2.5}
\]

Finally, forcing `B` to have one triple root gives uniquely
`(b,c)=(1,-1)`, which has `M=0`.  Hence a separate `A2+A5` row is
impossible.  Equations (2.2)--(2.5) similarly rule out separate
`2A2+A3`, `2A4`, `A6+A2`, and `A8` rows.

The complete nonmerger list is therefore:

\[
\begin{array}{c|c|c}
&\text{affine packet}&(a,b,c,d)\\ \hline
1&4A_1&(1,0,0,0)\\
2&A_3+2A_1&(0,0,4,-3)\\
3&2A_3&(2,0,0,-5)\\
4&A_5+A_1&(3,-3,-17,-9)\\
5&A_7&(6,-4,-44,-37)\\
6&A_2+3A_1&(-3,0,1,-7)\\
7&A_2+A_3+A_1&(-3,-1/2,3,-9)\\
8&A_4+2A_1&(-3,-9/4,13/2,-9)\\
9&A_4+A_3&(-3,10,80,-205)\\
10&A_6+A_1&(-3,-7/2,-1,11)\\
11&2A_2+2A_1&(-3,0,0,-5)\\
12&A_4+A_2+A_1&(-3,5/2,-5,-5).
\end{array}                                                    \tag{2.6}
\]

These rows include the open immersed partitions, one- and two-cusp
strata, and every higher distinct-image specialization.

## 3. The nine merger packets

Now impose `M=0` and use (1.6).  The cubic discriminant is

\[
 \Delta_C=-4a^3-27(ab-c)^2.                    \tag{3.1}
\]

This gives a short exhaustive decision tree.

1. If `C` is squarefree, its three normalization points are smooth and
   have one common image.  If the line root `-b` is separate and
   off-diagonal, the packet is `D4+A1`; if it is separate and diagonal, it
   is `D4+A2`.  If `-b` is also a root of `C`, equivalently `c=-b^3`, one
   pair at the triple point has intersection multiplicity two and the
   packet is `D6`.
2. If `C=(u-rho)^2(u+2rho)` with `rho!=0`, then the simple root `-2rho`
   is diagonal and the double root `rho` records intersection multiplicity
   two between a cusp and a smooth branch.  This is `D5`.  A separate
   off-diagonal line root adds `A1`; the other diagonal root adds `A2`.
   If the line root joins `rho`, the intersection rises to three and gives
   `E7`; if it joins `-2rho`, the cusp rises to `A4` and gives the second
   `E7` realization.
3. A triple cubic forces `a=c=0`.  Then
   `R=u^3(u+b)`: for `b!=0` the packet is `E6+A1`, and for `b=0` it is
   `E8`.

Thus the merger list is exactly:

\[
\begin{array}{c|c|c}
&\text{affine packet}&(a,b,c,d)\\ \hline
13&D_4+A_1&(1,1,0,-2)\\
14&D_4+A_2&(-3,2,0,3)\\
15&D_5+A_1&(-3,0,-2,-9)\\
16&D_6&(1,1,-1,-3)\\
17&D_5+A_2&(-3,-2,4,-5)\\
18&E_7\text{ (intersection-three)}&(-3,-1,1,-7)\\
19&E_7\text{ (A4-cusp)}&(-3,2,-8,-13)\\
20&E_6+A_1&(0,1,0,0)\\
21&E_8&(0,0,0,0).
\end{array}                                                    \tag{3.2}
\]

The `D4+A2` row is essential.  It occurs generically when the three roots
of `C` are distinct and the residual line root is diagonal.  The simpler
point `(a,b,c,d)=(0,0,1,0)` is the same packet: its ordinary triple point
and ordinary cusp spend affine delta `3+1=4`.

## 4. Completeness and use in the Keller audit

The collision polynomial always has degree four, and its root
multiplicities are the four affine delta units of the rational quintic.
Proposition 1.1 separates all possibilities into `M!=0` and `M=0` without
losing an image-merger locus.  Section 2 exhausts every possible marking of
the roots by the two diagonal values; Section 3 exhausts the multiplicity
partitions of a cubic and its residual line.  Tables (2.6) and (3.2) are
therefore complete.

This repairs the earlier severe-stratum list.  It also changes the role of
the witness `(t^3+t,t^5-t)`: that curve is a `D6` merger with cyclic affine
complement, not a new unclassified escape.

The topological calculations have the following consequences, kept
separate from the algebraic classification:

- every row with cyclic affine-complement group is impossible as the sole
  ramified nonproperness component because its local meridian fixes affine
  sheets;
- the `E6/E8` rows are impossible as the complete exceptional set by
  Chau's common-coordinate theorem;
- at geometric degree six, the only remaining cubic-inertia action occurs
  on `A4+A2+A1`, and its filled cover has nonzero first homology.

Together these facts close the **one-component, degree-six, `k=1`** atlas.
They do not force the global F2 cover to have degree six, exclude a second
nonproperness component, or decide `k=2,...,24`.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_complete_singularity_atlas.py
```
