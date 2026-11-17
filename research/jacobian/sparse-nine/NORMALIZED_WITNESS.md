# A 49 term normalized Keller collision and its frame minimum

There is an explicit rational three-variable map H with

\[
\det JH=1,\qquad JH(0)=I,\qquad H(0)=H(e_1)=0,
\]

coordinate degrees (4,6,7), and nonlinear support counts (5,18,26).
It has 49 nonlinear monomial occurrences.

This is the known foundational counterexample in a sparser collision frame,
not a new mechanism. The minimum 49 proved below applies to the specified
frames based on the boundary of that map. It does not establish a global
support minimum.

## Explicit map

Put

\[
Y=y-2,\qquad Z=z+24x-16,\qquad U=1+xY,
\]

and define

\[
\begin{aligned}
H_1&=x-\frac32x^2Y-\frac12x^3Z,\\
H_2&=y+3xU^2Z+3xY^2(4+3xY),\\
H_3&=U^3Z+Y^2U(4+3xY)+16H_2-64H_1.
\end{aligned}
\]

Expansion gives the counts above; all rational coefficients are retained in
[the witness](normalized-49-witness.json). Substitution immediately gives
H(0)=H(e1)=0. Its derivative at zero is the identity.

For the determinant certificate, write F for the
[foundational map](../../foundations/verified/FOUNDATIONAL_GEOMETRY.md).
Its two points

\[
a=(0,-2,-16),\qquad b=(1,-2,8)
\]

both map to (0,-2,0). Let

\[
L=\begin{pmatrix}1&0&0\\0&1&0\\24&0&1\end{pmatrix},
\qquad A=JF(a).
\]

Then

\[
H(X)=(AL)^{-1}\bigl(F(a+LX)-F(a)\bigr).
\]

Since det JF=-2, det A=-2 and det L=1, the determinant of JH is one.
The exact checker also expands the determinant directly.

## All frames of the specified boundary class

Consider any distinct collision pair a,b for F with a.x=0. The restriction
of F to x=0 is (z+4y^2,y,0), hence injective, so b.x=s is nonzero.
Complete the collision direction by the original transverse axes:

\[
L_{a,b}=[b-a\mid e_y\mid e_z].
\]

The frame is, by definition,

\[
H_{a,b}(X)=(JF(a)L_{a,b})^{-1}
          (F(a+L_{a,b}X)-F(a)).
\]

The verified torus equivariance

\[
F(\lambda x,\lambda^{-1}y,\lambda^{-2}z)
=\operatorname{diag}(\lambda^{-2},\lambda^{-1},\lambda)\,F(x,y,z)
\]

lets us set s=1. This changes H only by diagonal source-target conjugation,
which preserves support cardinality. Specifically, for the torus source
matrix T and transformed natural completion L', one has
TL=L' diag(1,lambda^{-1},lambda^{-2}); the normalized maps are diagonally
conjugate.

Write the transformed b.y as w. The equations F(b)=F(a) then force

\[
\begin{aligned}
a(w)&=(0,4w+6,-63w^2-189w-142),\\
b(w)&=(1,w,2-3w).
\end{aligned}
\]

Indeed, F3(b)=0 fixes b.z=2-3w; then F2(b)=4w+6 and
F1(b)=(1+w)(2+w). Comparing with (a.z+4*a.y^2,a.y,0) gives a(w).
Thus the displayed one-parameter family exhausts this boundary frame class,
up to support-preserving diagonal changes.

No arbitrary transverse linear completion, nonlinear change of coordinates,
or collision pair with both x-coordinates nonzero is included.

## Exact sharpness within this class

Let H_w be the normalized map built from a(w),b(w) with transverse columns
ey,ez. Its generic nonlinear support has 67 occurrences. Every coefficient
is a polynomial in Q[w].

The [complete coefficient certificate](frame-family-certificate.json) retains
all 67 coefficients and the 30 distinct monic irreducible factors occurring
in them. For each such factor f, it counts a coefficient as zero precisely
when its polynomial remainder modulo f is zero. This determines support
at every root of f, including nonrational algebraic parameters.

The independent verifier reconstructs every coefficient in Q[x,y,z,w],
then checks complete factor products and pairwise polynomial gcds using
standard-library rational arithmetic. It does not need to trust the
symbolic irreducibility labels: the pairwise coprime factor basis and
complete product identities already imply identical zero patterns at every
root of each factor. A parameter that is not a root of any listed factor
has the generic support count. Therefore this
finite factor calculation covers all parameters over every characteristic-zero
field, including transcendental parameters.

The smallest strata are:

| Parameter condition | Nonlinear occurrences |
|---|---:|
| w=-2 | 49 |
| w=-3/2 | 54 |
| 21*w^2+62*w+48=0 | 59 |
| every other exceptional factor | at least 65 |
| generic parameter | 67 |

Consequently **49 is the sharp minimum in this boundary frame class**,
attained only at w=-2 after the torus normalization.

The map has the boundary profile (7,6,4). The ten-term lower bound concerns
all normalized raw-degree-at-most-seven collisions. The later
[expanded transverse-frame theorem](../collision-frame-search/PROOF.md)
constructs 48 terms outside this fixed-axis class, improving the current
normalized interval to [10,48]. The class minimum49 above remains valid.

