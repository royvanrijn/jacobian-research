# The cubic map as a marked-root space

This note repackages the affine cubic and fiber calculations in
[IMAGE_AND_NONPROPERNESS.md](IMAGE_AND_NONPROPERNESS.md) as one geometric
statement.  The affine reconstruction and the `3/1/0` fiber theorem were
already recorded there; the new organizing idea is to retain the inverse root
projectively and regard it as a marked root of a binary cubic.

Work over a field of characteristic zero.  Write `F(x,y,z)=(a,b,c)` for the
map in the README.

## The affine-root chart

On `x!=0`, put

\[
u=y+{1\over x},\qquad v={1\over x}.
\]

Substitution in the three coordinates of `F`, followed by elimination of
`x,y,z`, gives

\[
a=u^2+uv-cu^3,\qquad b=4u+2v-3cu^2.                 \tag{1}
\]

For

\[
P(T)=cT^3-2T^2+bT-2a,
\]

equation (1) immediately yields

\[
P(u)=0,\qquad P'(u)=3cu^2-4u+b=2v.                 \tag{2}
\]

Thus an affine root `u` is simple exactly when `v!=0`.  Conversely, a simple
root reconstructs a unique source point:

\[
v={P'(u)\over2},\qquad x={1\over v},\qquad y=u-v,
\qquad z=5v^2-3uv-cv^3.                            \tag{3}
\]

Substitution of (3) in `F` gives `(a,b,c)`.  These are the same formulas as
the `t,r` formulas in [IMAGE_AND_NONPROPERNESS.md](IMAGE_AND_NONPROPERNESS.md),
with `t=u` and `r=2v`.

## The projective root

For a target `(a,b,c)`, define the binary cubic

\[
Q_{a,b,c}(U,V)
 =cU^3-2U^2V+bUV^2-2aV^3.                          \tag{4}
\]

Its `V!=0` dehomogenization is `P(U/V)`.  Direct expansion of the coordinates
of `F` proves the polynomial identity

\[
Q_{F(x,y,z)}(1+xy,x)=0.                             \tag{5}
\]

Consequently every source point marks the projective root

\[
[U:V]=[1+xy:x].                                     \tag{6}
\]

For `x!=0`, (6) is the affine root `U/V=u=y+1/x`.  Unlike that affine
coordinate, (6) continues across `x=0`.

## The root-at-infinity chart

On `U!=0`, set `s=V/U`.  Equation (4) becomes

\[
R(s)=c-2s+bs^2-2as^3=0.                             \tag{7}
\]

Put

\[
d=1-bs+3as^2=-{1\over2}R'(s).                       \tag{8}
\]

The root is simple on this chart exactly when `d!=0`.  Solving the source
coordinates gives

\[
x={s\over d},\qquad y=b-3as.                        \tag{9}
\]

The third coordinate can first be written, for `s!=0`, as

\[
z={2sd^2-3s^2d(b-3as)-cd^3\over s^3}.
\]

Using (7), so that `c=2s-bs^2+2as^3`, the numerator is divisible by `s^3`.
After cancellation one obtains the expression regular at `s=0`

\[
\begin{aligned}
z=d\big(&a-4b^2+(b^3+22ab)s-(30a^2+8ab^2)s^2\\
        &+21a^2bs^3-18a^3s^4\big).                 \tag{10}
\end{aligned}
\]

Substitution of (7)--(10) in `F` returns `(a,b,c)`.  At `s=0`, the root
equation forces `c=0`, while `d=1`; hence

\[
(x,y,z)=(0,b,a-4b^2).                               \tag{11}
\]

Thus the existing `x=0` source chart is precisely the root `[1:0]` at
infinity.  On the overlap `s!=0`, one has

\[
u={1\over s},\qquad v={d\over s},
\]

and (9)--(10) agree with (3).

## The marked-root theorem

Let

\[
I=\{((a,b,c),[U:V])\in\mathbb A^3\mathbin\times\mathbb P^1:
Q_{a,b,c}(U,V)=0\}.
\]

The projection `I -> A^3` is the finite flat degree-three incidence cover:
the fixed coefficient `-2` prevents an entire projective-line fiber, and each
fiber is a degree-three Cartier divisor.  Let `I^simp` be its simple-root
locus.  Concretely:

- on `V!=0`, with `u=U/V`, simplicity means `P'(u)!=0`;
- on `U!=0`, with `s=V/U`, simplicity means `R'(s)!=0`, equivalently `d!=0`.

These conditions agree on the overlap.  They include the root at infinity:
when `s=0`, equation (7) says `c=0` and `R'(0)=-2`, so `[1:0]` is always
simple whenever it is a root.

**Theorem.**  The morphism

\[
\Phi:\mathbb A^3_{x,y,z}\longrightarrow I^{\rm simp},\qquad
(x,y,z)\longmapsto(F(x,y,z),[1+xy:x])               \tag{12}
\]

is an isomorphism.  Under this isomorphism, `F` is the projection that
forgets which simple projective root was marked.

**Proof.**  Identity (5) puts the image in `I`.  If `x!=0`, equations (2)
and (6) give `P'(u)=2/x!=0`.  If `1+xy!=0`, direct substitution in (8) gives

\[
d={1\over1+xy},
\]

so the marked root is again simple.  These source opens cover because `x`
and `1+xy` cannot vanish simultaneously.

On the `V!=0` part of `I^simp`, formulas (2)--(3) give a regular inverse:
`P'(u)` is invertible there.  On the `U!=0` part, formulas (8)--(10) give a
regular inverse: their only denominator is a power of the invertible function
`d`, and (10) is polynomial at `s=0`.  The overlap identities
`u=1/s`, `v=d/s` show that the inverses agree.  Substitution verifies both
compositions, proving (12).  The target component of (12) is `F`, so it is
the forgetful projection.  \(\square\)

The theorem is global: the regular expression (10), not merely its limit, is
the extra regularity argument needed at the projective root at infinity.

## Determinant and the removed ramification divisor

On the affine-root chart, (1) gives

\[
\det {\partial(a,b,c)\over\partial(u,v,c)}=2v.
\]

Meanwhile, from `u=y+1/x`, `v=1/x`, and
`c=2x-3x^2y-x^3z`,

\[
\det {\partial(u,v,c)\over\partial(x,y,z)}=-{1\over v}.
\]

Their product is `-2`.  The equation `v=0` is the affine-chart part of the
ramification divisor of the finite marked-root cover `I -> A^3`; globally it
is the repeated-marked-root locus removed in passing to `I^simp`.  It is not a
branch divisor of the affine Keller map `F`.  The latter is the restriction to
the simple-root locus and is everywhere etale, with determinant `-2`.

## Fibers, discriminant, and the omitted locus

Over an algebraically closed field, the projective roots of a nonzero binary
cubic have one of three multiplicity types.  Because a fiber of `F` is the set
of simple roots, the theorem gives

| Binary-cubic root type | Number of source points |
|---|---:|
| three distinct roots | 3 |
| one double and one simple root | 1 |
| one triple root | 0 |

In the repository's sign convention,

\[
Q_{\rm disc}=27a^2c^2-18abc+16a+b^3c-b^2.           \tag{13}
\]

Thus `Qdisc!=0` is the first row and `Qdisc=0` gives the last two rows.
The raw polynomial discriminant of `P` is `-4 Qdisc`; after suppressing the
inessential factor `4`, the ordinary cubic-discriminant convention used in
some sources is `-Qdisc`.  In particular, their sign is the opposite of the
repository's `Qdisc` convention.

The cubic is a cube exactly on

\[
3bc=4,\qquad 27ac^2=4.                              \tag{14}
\]

Here `c!=0`; equations (14) are equivalent to the repository's
`3bc=4, 12a=b^2`.  This is the triple-root locus `Gamma`, hence the omitted
locus.  The exact image and nonproperness arguments remain those in
[IMAGE_AND_NONPROPERNESS.md](IMAGE_AND_NONPROPERNESS.md); the marked-root
theorem packages their fiber count but does not replace their escaping-path
and boundedness proofs.

## The displayed collision

At `(a,b,c)=(-1/4,0,0)`,

\[
Q_{-1/4,0,0}(U,V)
=-2U^2V+{1\over2}V^3
=-2V(U-V/2)(U+V/2).
\]

Its three projective roots are

\[
[1:0],\qquad [1/2:1],\qquad [-1/2:1].
\]

They correspond respectively to

\[
(0,0,-1/4),\qquad
(-1,3/2,13/2),\qquad
(1,-3/2,13/2),
\]

under the marking `[1+xy:x]`.

## Relation with the weighted family

For this cubic, `I -> A^3` is a finite marked-root cover and the source is the
complement of its repeated-marked-root divisor.  The general weighted
inverse-pencil construction has the same pattern: form the finite cover that
marks a root of `H(W)-sW+t`, normalize it over the degenerate `C=0` fiber,
and retain its regular-reconstruction open.  The precise theorem is
[WEIGHTED_MARKED_ROOT_MODEL.md](WEIGHTED_MARKED_ROOT_MODEL.md).
[DICRITICAL_COMPACTIFICATION.md](DICRITICAL_COMPACTIFICATION.md) (Claim C16)
classifies the reconstruction-pole divisors on the corresponding normalized
inverse-graph compactification.

This interpretation is adopted with attribution to Andy Jiang's public
geometric-interpretation post; see [SOURCES.md](SOURCES.md) and
[PROVENANCE_AUDIT.md](PROVENANCE_AUDIT.md).  The overlap with the repository's
earlier C02--C03 computations is recorded there, and no priority claim is made.
