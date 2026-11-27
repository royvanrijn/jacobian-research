# The tangent norm class on a quadratic companion twist

The fixed elliptic double cover of the tangent quartic
gives a primitive point on a quadratic twist whose
coefficient depends only on the field. Its degree-three
Kummer class is the displayed edge norm class up to
inversion. Both isogeny Selmer groups of the companion
have exact class-field descriptions. This yields
rank, conductor and Tate–Shafarevich consequences on
the same prime configurations.

## The primitive point realizes the norm class

Use [CT1](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#the-tangent-identity-and-the-lost-cubic-direction)
and write \(-R(p,h)=d\delta^2\), with \(d\) signed
squarefree and \(\delta>0\).
Thus \(K=\mathbb Q(\sqrt d)\).
Use the double-cover model of
[FC1](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#the-fixed-elliptic-double-cover):
\[
 t=p/h,\quad U=\delta\sqrt d/h^2,\quad
 E:y^2=x^3+1,\quad T=(-1,0),\quad
 x(P)=\frac{U+1+t^2}{2t^2},\quad y(P)=\frac{x(P)+1}{t}.
\]
The point \(P\) has trace \(T\).
Write \(G_d:y^2=x^3+d^3\) and
use \((x,y)\mapsto(dx,d\sqrt d\,y)\)
on the anti-invariant subgroup to
define \(Q_d\) from \(Q=2P-T\).

**Theorem CD1.** For every genuine
tangent prime progression, \(Q_d\)
is primitive in \(G_d(\mathbb Q)\)
modulo torsion. Its \(x\)-coordinate
and its degree-three Kummer class
are given by
\[
 x(Q_d)=12(p h/\delta)^2-d,\qquad
 \bigl[y(Q_d)+d\sqrt d\bigr]=[\alpha_-]^2
              \quad\text{in }K^\times/K^{\times3}.
\]
For a fixed \(d\), all displayed points
\(Q_d\) lie in one coset modulo
\(2G_d(\mathbb Q)\).

In the notation of
[FC2](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#the-three-marked-fibres-force-primitivity),
\(H(K)\) has trace \(O\) or \(T\),
has torsion \(\{O,T\}\), and \(P\)
is primitive modulo torsion.
If \(2P-T\) is divisible by an integer
\(n>1\) in the free anti-invariant
group, its image \(2[P]\) in the
free group of \(H(K)\) is divisible
by \(n\). Primitivity forces \(n=2\),
and the proposed half has the same
free image as \(P\). Its trace would
be \(T\) since the torsion has trace
zero, contradicting anti-invariance.
This proves primitivity of \(Q_d\).
Taking differences for two trace-\(T\)
points \(P\) gives the common coset
modulo twice the anti-invariant group.

Put \(u=x(P)+1\) and \(y=y(P)\).
The quadratic for \(u\) is
\[
 t^2 u^2-(1+3t^2)u+3t^2=0,\qquad
 u-\frac3u=\frac{U}{t^2}.
\]
The doubling formula and then
translation by \(T\) give
\[
 x(2P)+1=\frac{(u^2-3)^2}{4y^2},\qquad
 x(2P-T)+1=\frac{3}{x(2P)+1}
          =\frac{12t^2}{U^2}.
\]
Scaling to \(G_d\) gives the asserted
coordinate.

For the degree-three statement there
is an identity of rational functions
on the double cover:
\[
 \frac{2t^3-3t^2-1+U}{2}
       =(y(P)+1)\left(\frac{1-t}{x(P)}\right)^3.
\]
It follows on substituting
\(t=(x+1)/y\), \(U=2t^2x-1-t^2\)
and \(y^2=x^3+1\).
The left side is \(\alpha_-/h^3\).
Choose the plus convention for the
[BES Kummer map](https://arxiv.org/html/1610.05759#S3.SS1)
on \(E/K\), represented by \(y+1\)
on the indicated affine points.
Its value at \(T\) is trivial.
Therefore the class at \(Q=2P-T\)
is \([\alpha_-]^2\).
The twist map multiplies \(y+1\)
by \(d\sqrt d=(\sqrt d)^3\)
and preserves the class.

## The two local Selmer conditions

Write \(m=-d/3\) and \(K_m=\mathbb Q(\sqrt m)\).
For every genuine tangent progression,
\(d\ne-3\), \(d\equiv5\bmod8\),
\(3\mid d\) and \(d/(-3)\) is a
square over \(\mathbb Q_3\).
In particular \(m\equiv1\bmod8\)
and \(m\equiv1\bmod3\). These local
conditions below can also be used
for signed squarefree \(d\) satisfying
them without a witness.
The natural isogeny goes from
\(G_d\) to \(G_m\), since
\(-27d^3=3^6m^3\).
Denote it by \(\phi\), and its
dual by \(\widehat\phi\).

**Proposition CD2.** Consider
\(\operatorname{Sel}_{\widehat\phi}(G_m)\)
in \(K^\times/K^{\times3}\)
and \(\operatorname{Sel}_\phi(G_d)\)
in \(K_m^\times/K_m^{\times3}\),
with norm-cube condition on
each side. Their local images
are unramified at every prime
\(q>3\). At \(2\) the first
image is zero and the second
is the full valuation line
in the split algebra. At \(3\)
the first image is the
unramified Kummer line in
\(\mathbb Q_3(\zeta_3)\), and the
second is the unit line in
the split algebra.
The global ratio is
\[
 c(\phi)=3^{2-\mathbf{1}_{d>0}}.
\]

At \(2,3\) the two curves are
locally the constant pair
\(E_{-27}, E_1\) in that order.
For \(E_{-27}\) at \(2\), the
unramified quadratic algebra
has residue field \(\mathbb F_4\).
A nonintegral \(x\) makes the
rational ordinate dominate the
Kummer constant, so after
removing a rational cube
the residue is rational.
An even integral \(x\)
would give \(y^2=5\bmod8\).
An odd \(x\) gives even \(y\)
and the residue of
\(y+3\sqrt{-3}\) is again rational.
Both cases give the trivial class.
On \(E_1\) over \(\mathbb Q_2\)
the point \((0,1)\) gives class
\(2\), which generates the
split local group.

At \(3\) put \(\lambda=\sqrt{-3}\),
so \(3\lambda=-\lambda^3\).
For \(E_{-27}\), a unit \(x\)
gives unit \(y\) and the class of
\(y+3\lambda\) lies on the
unramified line represented by
\(1+\lambda^3\mathcal O\);
a rational unit in the
anti-invariant calculation
does not affect this conclusion.
More explicitly, \(N(y+3\lambda)=x^3\)
and the class of \(y\) itself is
already trivial locally: its
rational norm-square is congruent
to a cube modulo \(9\), since
\(y^2=x^3-27\).
The nonintegral case follows
after normalizing by a rational
cube, with the irrational term
of valuation at least four and
the rational unit again a cube
by the norm identity modulo nine.
If \(v_3(x)=1\), writing \(x=3u\)
shows \(v_3(y)\ge3\) when
\(y\ne0\): \(u\equiv-1\bmod3\)
would give odd valuation of
\(y^2\), and \(u\equiv1\bmod3\)
gives \(v_3(u^3-1)\ge2\).
After dividing by the cube
\(3\lambda\), the irrational
Kummer calculation again lies
in \(1+\lambda^3\mathcal O\).
The case \(y=0\) here is trivial
in the cube quotient.
Valuation \(v_3(x)\ge2\) is
impossible. The point with
\(x=4,y^2=37\) gives a nonzero
class on the unramified line,
using the \(\lambda^3,\lambda^4\)
thresholds of
[RT1](ABC_COHERENT_REFLECTION_AT_THREE.md#the-two-local-thresholds).

For \(E_1\) at \(3\), the
class of \(y+1\) always has
valuation divisible by three.
This follows from
\((y+1)(y-1)=x^3\) and
coprimality at \(3\) when \(y\)
is integral. When \(y\) is
nonintegral its valuation is
already divisible by three.
The point \((0,1)\) gives the
nontrivial unit class \(2\).
Thus the image is the unit line.

At \(q>3\) dividing \(d\),
both quadratic algebras are
ramified and have zero
anti-invariant cube quotient.
Conjugation is trivial on
valuations modulo three and
on residue units, and every
principal unit is a cube.
At the other \(q>3\) the
curves have good reduction,
and the prime-to-\(q\)
isogeny image is the full
unramified local group.
At infinity the cohomology
groups are zero. The kernels
of \(\phi\) at \(2,3\) have
one rational point, so the
local \(c_q(\phi)\) are
both \(3\) by the second
image computations.
Elsewhere the finite local
ratios are \(1\).
At the real place, the
kernel has size \(3\) for
\(d>0\) and \(1\) for
\(d<0\); multiplication
through the isogeny on the
connected real elliptic
group is surjective. This
gives the stated product.

## Class-field descriptions of both Selmer groups

**Theorem CD3.** Let \(\mathfrak p_2\)
be either prime of \(K_m\) above \(2\).
For signed squarefree \(d\ne-3\)
with the local conditions of CD2,
the two Selmer groups identify with
\[
 \begin{split}
 \operatorname{Sel}_{\widehat\phi}(G_m)
   &\simeq \operatorname{Hom}
       (\mathrm{Cl}(K_m)/\langle[\mathfrak p_2]\rangle,\mathbb F_3),\\
 \operatorname{Sel}_{\phi}(G_d)
   &\simeq \operatorname{Hom}
       (\mathrm{Cl}_{(6)}(K),\mathbb F_3),
 \end{split}
\]
where \(\mathrm{Cl}_{(6)}(K)\)
is the ordinary ray class group
of modulus \(6\mathcal O_K\).
Their dimensions differ
by \(2-\mathbf{1}_{d>0}\),
with the ray dimension larger.

To see the reflection, extend a
Kummer class to \(K(\zeta_3)\),
whose three quadratic characters
correspond to \(K,K_m,\mathbb Q(\zeta_3)\).
The norm-cube condition in
\(K\) puts the class in the
\(K\)-eigenspace. Turning it
into a cubic character also
multiplies the eigensign by
the cyclotomic character.
The result descends to \(K_m\)
and is anti-invariant under its
conjugation. In the other
direction the same argument
reflects norm-cube classes
from \(K_m\) to anti-invariant
characters of \(K\).
Restriction to a quadratic
extension and the norm
are inverse up to a factor
of two on the appropriate
eigenspaces, so this is
an equivalence on the
cubic character groups.

On the \(K\)-Kummer side, the
condition at \(3\) is exactly
unramified after adjoining
\(\zeta_3\). The unit
conditions at other odd
primes are also exactly
unramified. The zero
condition at \(2\) requires
complete splitting of either
prime above \(2\) in the
reflected field \(K_m\).
The other prime then splits
as well by conjugation.
These conditions give the
first class-field description.

In the converse direction
the mirror is split at \(2,3\).
At \(2\) the split valuation
line reflects to all the
tame anti-invariant cubic
characters of the unramified
quadratic \(K_2\). The
conductor exponent is at
most one there.
At \(3\) the split unit
line is generated by \(2\).
Over \(K_3=\mathbb Q_3(\zeta_3)\)
the cubic root of \(-2=1+\lambda^2\)
generates a cubic extension
with conductor exponent two.
In fact the Newton polygon
for \(\xi^3=1+\lambda^2\)
gives \(v(\xi-1)=2/3\)
for \(v(\lambda)=1\).
In the extension
\(\pi=\lambda/(\xi-1)\)
is a uniformizer and
\(v_{\rm ext}(\sigma\pi-\pi)=2\)
for nontrivial cubic
\(\sigma\), by writing
\(\sigma\xi=\zeta_3\xi\).
The lower break is one.
The anti-invariant characters
of conductor at most two
constitute one line: the
valuation character is
invariant, while the
conjugation sign on
\(U_1/U_2\) is negative.
This is exactly the displayed
local character line.
The reflected conductor
therefore divides
\(2\mathfrak p_3^2=6\mathcal O_K\).

All ordinary ray class
characters modulo \(6\) are
anti-invariant in their
cubic part. An ideal
times its conjugate is
generated by its rational
norm, which is a unit
modulo \(6\) and can be
made \(1\bmod6\) by
changing sign. Thus
conjugation acts by
inversion. Infinite
places impose no extra
cubic condition.
This gives the second
description. Finally
[BES, equation (21)](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf)
and CD2 show the dimension
difference, since neither
rational isogeny kernel
has a nonzero point.

The first line is a
class-field reflection with
a dyadic splitting condition.
In the squarefree prime-gap
setup of
[CS1](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#three-local-images-in-a-common-quadratic-algebra),
intersect this companion
Selmer group with any
of the three cubic-twist
groups \(S_e\).
The only additional
condition is that its local
class at \(k=h/6\) vanish.
At a split \(k\), the
valuation line of \(S_e\)
has trivial intersection
with the companion unit
line, and at an inert
\(k\) each cubic-twist
local image is zero.
The conditions at \(2,3\)
on the companion side
are contained in all
the conditions on the
cubic-twist side. This
identifies the intersection
exactly. In particular
the displayed edge class
is in the intersection
and is rational on the
companion by CD1.

## Defined conductors and exact companion Sha groups

**Corollary CD4.** For the fields of
genuine tangent progressions the
companion curve has conductor
\(4d^2\), root number
\(-\operatorname{sgn}(d)\), and a
primitive nontorsion point.
In the imaginary case its
\(3\)-Selmer dimension is at
least two and even. For
each sufficiently large \(Y\)
there are
\(\gg Y^{1/4}/(\log Y)^3\)
distinct conductors at most
\(Y\) of these imaginary
companions with the
displayed edge class of
order three in the
maximal quadratic class
group. The lower bound
with prime \(h/6\)
has log power four.

At \(2,3\) the local model
is \(E_{-27}\), with
conductor exponents \(2,2\).
At the remaining primes
dividing \(d\) the
coefficient has valuation
three and the reduction
has type \(I_0^*\), with
conductor exponent two.
The other primes have good
reduction. Thus the
conductor is
\(36m^2=4d^2\).
Using
[Várilly-Alvarado, Lemma 4.1 and Proposition 4.4](https://msp.org/ant/2011/5-5/ant-v5-n5-p04-p.pdf),
the local root factors at
\(2,3\) are \(-1,+1\).
At \(q\ge5\) dividing
\(m\) the factor is
\((-1/q)\), so the
global root number is
\((-1/|m|)\), equal to
\(-\operatorname{sgn}(d)\)
since \(m\equiv1\bmod4\).
The finite \(3\)-Selmer
parity in
[BES, Proposition 49](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf)
is \(2-\mathbf{1}_{d>0}\)
modulo two by CD2,
with no rational
three-torsion.
The exhibited point
has positive rank,
so the imaginary
Selmer dimension is
at least two.

For the lower count one
can use the imaginary
squarefree fields of
CT4 with \(d\asymp \sqrt Y\)
and the constant chosen
so \(4d^2\le Y\).
They have order-three
edge classes and
distinct conductors.
Alternatively, in the
balanced imaginary
rectangle before its
squarefree sieve, FM1
removes the
\(O(Z^{2-\alpha+\varepsilon})\)
large-scale directions
with \(s>Z^\alpha\)
for any fixed
\(0<\alpha\le1/3\).
On the remaining ones
\(|d|\gg Z^{4-2\alpha}\)
exceeds \(4p(p-h)\).
The ideal of that
norm cannot be principal
by its one-sided supports
and the lower bound on
a nonrational generator's
norm. FM3 gives the
prime-form asymptotic
for distinct fields in
this larger fixed
rectangle sample, with
all but a power-saving
exceptional number
having the specified
edge direction.

**Proposition CD5.** For the
four indicated companion
curves, the exact ranks
and \(3\)-torsion in
Tate–Shafarevich are
\[
\begin{array}{c|r|c|c|c|c}
 (p,h)&d&r&\dim S&\dim T&|\Sha(G_d)[3]|\\ \hline
 (17,6)&-186843&4&1&3&1\\
 (101,30)&-256285803&2&1&3&9\\
 (167,114)&10179717&1&1&2&9\\
 (37,30)&2580117&1&1&2&9 .
\end{array}
\]
Here \(S=\operatorname{Sel}_{\widehat\phi}(G_m)\)
and \(T=\operatorname{Sel}_\phi(G_d)\)
in the table only.
In the two rank-one rows,
the stated prime progression
is the only unweighted one
in its field. The edge
ideal classes in the
first three rows have
order three. The last
row's noncube edge
class is a unit direction.

The ordinary quadratic class
groups, the class of a prime
over \(2\) in the mirror
and the ray group modulo
\(6\) are computed by
PARI ideal and ray
arithmetic with certified
ordinary quadratic class
groups. The two-descent
rank bounds coincide.
In each row the displayed
edge class is nontrivial
by the exact cube test.
By CD1 the one-dimensional
group \(S\) is realized by
rational points, so
\(\Sha(G_m)[\widehat\phi]=0\).
With no rational
three-torsion, the other
rational isogeny quotient
has dimension \(r-1\).
It follows that
\(\dim\Sha(G_d)[\phi]=\dim T-r+1\).
The exact sequence for
\(\widehat\phi\phi=[3]\)
then gives the stated
\(3\)-torsion on \(G_d\).
In the rank-one cases
use FC2 for uniqueness.
The ideal classes are
computed from the two
one-sided prime factors
of the edge generator.
These calculations,
including the base local
conductor exponents,
are reproducible with
[the Sage checker](verify_cubic_tangent_companion_descent.py).
The stated torsion
subgroups are determined
without a finiteness
hypothesis for full
Tate–Shafarevich.

The original rank-two
Mordell curves of CT4
are cubic twists at a
different conductor scale.
The companion construction
supplies a rational
primitive point and
two exact class-field
interfaces. A bound
for the companion rank
on average over this
marked subset is still
required for an
unbounded family of
nonzero companion Sha
from the isogeny
Selmer dimensions.

## Generic rank on the parameter line

**Proposition CD6.** Over \(\mathbb Q(t)\) the
companion \(Y^2=X^3+f(t)^3\) has torsion
of order two and free rank one. The
section of CD1 generates its free
part. The same rank-one conclusion
holds for the quadratic companion
of every fixed weighted tangent
quartic in WM1. In the weighted case
the section \(2P-T\) is primitive
precisely when \(T\) is not twice
a rational point on \(E_w\).

Over the function field of the
quartic double cover the twist
is the fixed elliptic curve.
A rational point over this function
field is a morphism from the smooth
projective cover to the fixed
curve. On the elliptic model each
such morphism over \(\mathbb Q\)
has the form \([n]P+B\), since
the rational endomorphism ring
is \(\mathbb Z\); use the
rational differential and period
lattice argument of
[FC3](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#field-preserving-maps-of-the-parameter-line).
The anti-invariance condition
for the quadratic twist is
\[
 [n](T-P)+B=-([n]P+B),
 \qquad\text{equivalently}\qquad
 2B=-[n]T .
\]
The map sending the section to
\(n\) therefore has kernel the
rational two-torsion.
All even \(n\) are admitted.
An odd \(n\) is possible
exactly when \(T\) is twice a
rational point. In the
unweighted case that is
impossible by FC1, so the
image is \(2\mathbb Z\) and
the coefficient \(n=2\) of
the section in CD1 generates.
For a weighted cover the
same calculation proves the
claims. Thus the imaginary
specializations with their
even finite \(3\)-Selmer
rank cannot be analyzed by
seeking a second independent
section on this same
parameter line.
