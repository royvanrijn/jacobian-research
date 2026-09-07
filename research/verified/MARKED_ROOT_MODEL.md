# The cubic map as a marked-root space

This note gives only the marked-root interpretation of the foundational map
and the two reconstruction charts needed to make that interpretation global.
The second chart is essential: it identifies the source divisor `x=0` with a
simple root at infinity.  The normalized factorization proof belongs to
[NORMALIZED_FACTORIZATION_MODEL.md](NORMALIZED_FACTORIZATION_MODEL.md), the
projective hyperplane geometry to
[FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md](FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md),
and the image and fiber analysis to
[IMAGE_AND_NONPROPERNESS.md](IMAGE_AND_NONPROPERNESS.md).

Work over a field of characteristic zero.  Write
`F(x,y,z)=(a,b,c)` for the map in
[FOUNDATIONAL_GEOMETRY.md](FOUNDATIONAL_GEOMETRY.md), and attach to a target
the binary cubic

\[
 Q_{a,b,c}(U,V)=cU^3-2U^2V+bUV^2-2aV^3.             \tag{1}
\]

Let

\[
 I=\{((a,b,c),[U:V])\in\mathbb A^3\times\mathbb P^1:
 Q_{a,b,c}(U,V)=0\}.                                \tag{2}
\]

The projection `I -> A^3` is the finite flat degree-three incidence cover:
each fiber is the projective root divisor of (1), and the fixed coefficient
`-2` prevents an entire projective-line fiber.  Let `I^simp` be the locus on
which the marked projective root is simple.

Equivalently, under

\[
 \operatorname{Sym}^3(\mathbb P^1)\simeq\mathbb P^3,
\]

this is the incidence of a binary cubic with one root marked, restricted away
from the repeated-marked-root divisor and from the tangent-nonosculating
hyperplane used to make the cubic coefficient chart affine.  The projective
geometry and the exceptional `(2,1)` contact type are proved in the incidence
note cited above.

## The affine-root reconstruction chart

On `V!=0`, put `u=U/V` and dehomogenize (1):

\[
 P(T)=cT^3-2T^2+bT-2a.                              \tag{3}
\]

On the source open `x!=0`, set

\[
 u=y+{1\over x},\qquad v={1\over x}.
\]

Direct substitution in `F` gives

\[
 a=u^2+uv-cu^3,qquad b=4u+2v-3cu^2.                \tag{4}
\]

Consequently

\[
 P(u)=0,qquad P'(u)=3cu^2-4u+b=2v.                 \tag{5}
\]

Thus an affine root `u` is simple exactly when `v!=0`.  Conversely a simple
affine root reconstructs a unique source point:

\[
 v={P'(u)\over2},\qquad x={1\over v},\qquad y=u-v,
\]

\[
 z=5v^2-3uv-cv^3.                                   \tag{6}
\]

Substitution of (6) into `F` returns `(a,b,c)`.  These are the same formulas
as the `t,r` formulas in
[IMAGE_AND_NONPROPERNESS.md](IMAGE_AND_NONPROPERNESS.md), with `t=u` and
`r=2v`.

## The global marked projective root

Direct expansion of the coordinates of `F` proves the polynomial identity

\[
 Q_{F(x,y,z)}(1+xy,x)=0.                             \tag{7}
\]

Every source point therefore marks the projective root

\[
 \boxed{[U:V]=[1+xy:x].}                            \tag{8}
\]

For `x!=0`, this is the affine root `U/V=u=y+1/x`.  Formula (8), unlike that
affine coordinate, remains defined across `x=0` because `x` and `1+xy` cannot
vanish simultaneously.

## The root-at-infinity reconstruction chart

On `U!=0`, set `s=V/U`.  Equation (1) becomes

\[
 R(s)=c-2s+bs^2-2as^3=0.                            \tag{9}
\]

Put

\[
 d=1-bs+3as^2=-{1\over2}R'(s).                      \tag{10}
\]

The root is simple on this chart exactly when `d!=0`.  The first two source
coordinates reconstruct as

\[
 x={s\over d},\qquad y=b-3as.                       \tag{11}
\]

For `s!=0`, the third coordinate initially has the form

\[
 z={2sd^2-3s^2d(b-3as)-cd^3\over s^3}.
\]

Using (9), hence `c=2s-bs^2+2as^3`, cancels the apparent pole and gives the
expression regular at `s=0`:

\[
\begin{aligned}
z=d\big(&a-4b^2+(b^3+22ab)s-(30a^2+8ab^2)s^2\\
        &+21a^2bs^3-18a^3s^4\big).                 \tag{12}
\end{aligned}
\]

Substitution of (9)--(12) into `F` returns `(a,b,c)`.  At `s=0`, the root
equation forces `c=0` and `d=1`, so

\[
 \boxed{(x,y,z)=(0,b,a-4b^2).}                      \tag{13}
\]

Thus the source divisor `x=0` is exactly the chart of the projective root
`[1:0]` at infinity.  It is always simple when it occurs, since
`R'(0)=-2`.  On the overlap `s!=0`,

\[
 u={1\over s},\qquad v={d\over s},
\]

and (11)--(12) agree with the affine reconstruction (6).

## Marked-root theorem

**Theorem.**  The morphism

\[
 \Phi:\mathbb A^3_{x,y,z}\longrightarrow I^{\rm simp},\qquad
 (x,y,z)\longmapsto(F(x,y,z),[1+xy:x])              \tag{14}
\]

is an isomorphism.  Under this isomorphism, `F` is the projection that
forgets which simple projective root was marked.

**Proof.**  Identity (7) puts the image in `I`.  If `x!=0`, equations (5) and
(8) give `P'(u)=2/x!=0`.  On the source open `1+xy!=0`, substitution in (10)
gives

\[
 d={1\over1+xy},
\]

so the marked root is simple there as well.  These two source opens cover.

On `V!=0` in `I^simp`, formulas (5)--(6) give a regular inverse because
`P'(u)` is invertible.  On `U!=0`, formulas (10)--(12) give a regular inverse
because `d` is invertible and (12) is polynomial at `s=0`.  The overlap
relations `u=1/s` and `v=d/s` show that the inverses agree.  Substitution
verifies both compositions.  The target component of (14) is `F`, so it is
the forgetful projection.  \(\square\)

The regular formula (12), rather than merely its limit at `s=0`, is the
additional argument that makes the marked-root interpretation global.

This projective-root organization is adopted with attribution to Andy
Jiang's public geometric-interpretation post; the provenance record is in
[SOURCES.md](../archive/legacy-notes/SOURCES.md) and
[PROVENANCE_AUDIT.md](../archive/legacy-notes/PROVENANCE_AUDIT.md).
