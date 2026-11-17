# Global Cox and boundary-deletion attack on the quartic JC(2) packets

## Status

This note is a reduction programme, not a proof of JC(2) or of the
geometric-degree-four case.

The local finite-flat, conductor-length, Euler, and sheet-graph attacks are
already exhausted.  The universal
[unibranch spectator models](UNIBRANCH_SPECTATOR_COUNTERMODELS.md) realize
the clean \(3+1\) cusp and a \(2+2\) self-collision in one finite-free
quartic cover.  They fail the Keller-normalization conditions only
globally: their deleted ramification complement is
\(\mathbb A^1\times\mathbb G_m\), not \(\mathbb A^2\).

The useful new conclusion is an exact Cox-lattice reduction:

> For every surviving quartic target grouping, the target boundary
> pullbacks miss one primitive saturation/complement character.  It can be
> represented by the class of the ramified boundary prime \(E\),
> equivalently the canonical or different class of the finite
> normalization.

Thus the proposed global obstruction is not class-group torsion by itself.
It is the incompatibility of deleting this primitive component while
retaining the affine companion cylinders in the canonical
\(\mathbb A^2\)-source open.

The later multi-factor attacks sharpen this point.  A unit-rank mismatch can
disappear after passing to a total family, and a complete exceptional
\(\mathbb P^1\) is forced only by particular two-chart completions.  Nonlinear
affine modifications can move or fill such a curve.  Accordingly, neither
"a unit appears" nor "a \(\mathbb P^1\) appears" is retained below as the
primary theorem.  They are possible witnesses for the canonical
boundary-deletion obstruction isolated in Section 6.

The exact integer audit is implemented by
`quartic_cox_boundary_lattice_atlas` in
[`cas/plane_boundary_exclusion.py`](cas/plane_boundary_exclusion.py) and
regressed by
[`cas/test_plane_boundary_exclusion.py`](cas/test_plane_boundary_exclusion.py).

## 1. Canonical finite-normalization setup

Let

\[
F=(P,Q):U=\mathbb A^2\longrightarrow Y=\mathbb A^2
\]

be a hypothetical noninvertible plane Keller map of geometric degree four.
Put

\[
A=k[P,Q]\simeq k[u,v],\qquad
B=\operatorname{Norm}_A k(x,y),\qquad
X=\operatorname{Spec}B.
\tag{1.1}
\]

The canonical factorization is

\[
U\lhook\joinrel\longrightarrow X
\mathop{\longrightarrow}^{\pi}Y.
\tag{1.2}
\]

The finite-normalization theorem gives:

1. \(B\) is finite free of rank four over \(A\);
2. \(X\setminus U\) has a nonempty divisorial part;
3. if its prime components are \(E_1,\ldots,E_r\), then
   \[
   B^\times=k^\times,\qquad
   \operatorname{Cl}(X)
   =\bigoplus_{i=1}^r\mathbb Z[E_i];
   \tag{1.3}
   \]
4. every ramification prime of \(\pi\) belongs to \(X\setminus U\); and
5. every target nonproperness component has a positive affine-sheet
   contribution.

Orevkov's global identity leaves the two rows in
[`JC2_QUARTIC_PACKET_FRONTIER.md`](JC2_QUARTIC_PACKET_FRONTIER.md):

- one ramified boundary \(E\), with one \(2\to3\) jump;
- one ramified boundary \(E\) and one unramified boundary \(D\), with no
  jump.

The two-boundary row splits further according to whether \(E\) and \(D\)
have the same target image.  These are the three target groupings audited
below.

## 2. What the coarse global invariants do not prove

The global hypotheses remove the spectator countermodel, but none of the
listed invariants gives an immediate contradiction alone.

### 2.1 Units and the free class group

The localization sequence is

\[
0\longrightarrow
U^\times/k^\times
\longrightarrow
\bigoplus_i\mathbb Z[E_i]
\longrightarrow
\operatorname{Cl}(X)
\longrightarrow
\operatorname{Cl}(U)
\longrightarrow0.
\tag{2.1}
\]

Both end groups vanish and the middle arrow is an isomorphism.  Hence no
nonzero boundary combination is principal and no nonconstant unit on
\(U\) is forced merely by listing the packet.

If \(C=(g=0)\) is the target image of the ramified boundary, the
one-boundary divisor relation is

\[
\operatorname{div}_X(g)=2E+A
\tag{2.2}
\]

or

\[
\operatorname{div}_X(g)=2E+A_1+A_2.
\tag{2.3}
\]

Here the \(A_i\) are closures of affine divisors meeting \(U\), not missing
boundary primes.  Equations (2.2)--(2.3) say only

\[
[A]=-2[E]
\quad\text{or}\quad
[A_1]+[A_2]=-2[E].
\tag{2.4}
\]

Both relations are compatible with the free boundary basis.

### 2.2 Canonical and different classes

The target has trivial canonical class.  In codimension one, the quartic
packets have one simple ramification prime and no other ramification
divisor.  Therefore

\[
K_X=\pi^*K_Y+R_\pi,\qquad
[K_X]=[R_\pi]=[E]
\quad\text{in }\operatorname{Cl}(X).
\tag{2.5}
\]

Equivalently, the dualizing module

\[
\omega_B\simeq\operatorname{Hom}_A(B,A)
\tag{2.6}
\]

has the primitive class \([E]\).  Finite freeness over \(A\) does not make
this rank-one \(B\)-module free.  Thus triviality of \(K_U\) is compatible
with (2.5): restriction kills the boundary class.

For reduced boundary \(E\) in the one-boundary row,

\[
[K_X+E]=2[E].
\tag{2.7}
\]

For reduced boundary \(E+D\) in the two-boundary row,

\[
[K_X+E+D]=2[E]+[D].
\tag{2.8}
\]

Again, restriction to \(U\) kills these classes.  There is no
canonical-class contradiction at this level.

### 2.3 Conductors are invisible to the divisor class group

The cusp and the \(2+2\) self-collisions contribute conductor at closed
points of the boundary curves.  On the normal surface \(X\), these are
codimension-two data.  They do not change \(\operatorname{Cl}(X)\) or the
codimension-one formula (2.5).

Consequently a conductor obstruction must be lifted to one of:

- an intersection or zero-cycle identity on a resolved completion;
- a Rees/Cox modification in which the two conductor branches become
  divisorial; or
- a global meridian relation retaining the conductor pairing.

No contradiction can follow from the class vector of the conductor alone.

### 2.4 Local monodromy is exhausted

A nondegenerate cusp gives two transpositions satisfying the braid relation
and generates \(S_3\) on three sheets.  One \(2+2\) collision gives a
perfect matching of the four sheets.  All \(24\cdot3=72\) combined packets
generate \(S_4\).

Thus sheet transitivity and the abstract local generators do not close the
packet.  A monodromy proof must use actual global relations in the target
curve complement, together with the primitive boundary character isolated
below.

## 3. Exact quartic Cox exponent lattices

For each target curve \(C_j=(g_j=0)\), write the boundary part of its
pullback as

\[
\sum_i b_{ji}E_i.
\tag{3.1}
\]

The vector \(b_j=(b_{j1},\ldots,b_{jr})\) is a boundary exponent, not a
relation asserting that the \(E_i\) vanish in \(\operatorname{Cl}(X)\).
The affine companion has class \(-b_j\).

The three quartic groupings give:

| grouping | boundary basis | target pullback vectors | affine companion class sums | target-lattice defect |
| --- | --- | --- | --- | --- |
| one boundary | \(E\) | \((2)\) | \((-2)\) | \(\mathbb Z/2\) |
| two boundaries, same target | \(E,D\) | \((2,1)\) | \((-2,-1)\) | free cokernel rank one |
| two boundaries, different targets | \(E,D\) | \((2,0),(0,1)\) | \((-2,0),(0,-1)\) | \(\mathbb Z/2\) |

In every row, append

\[
\epsilon_E=(1,0,\ldots,0).
\tag{3.2}
\]

Then the augmented exponent lattice is saturated:

\[
\begin{array}{c|c}
\text{grouping}&\text{unimodular/saturated certificate}\\ \hline
\text{one boundary}&\langle2,1\rangle=\mathbb Z,\\[1mm]
\text{same target}&
\det\begin{pmatrix}2&1\\1&0\end{pmatrix}=-1,\\[4mm]
\text{different targets}&
\gcd\left\{
\det\begin{pmatrix}2&0\\0&1\end{pmatrix},
\det\begin{pmatrix}0&1\\1&0\end{pmatrix}
\right\}=1.
\end{array}
\tag{3.3}
\]

This proves:

### Proposition 3.1 -- the unique missing primitive character

For every quartic target grouping:

1. target pullbacks alone do not give the complete primitive boundary
   lattice;
2. one additional primitive character suffices;
3. it can be chosen to be the ramified boundary character
   \(\epsilon_E\); and
4. this character is the canonical/different class \([K_X]=[E]\).

There is no unavoidable torsion after this character is retained.  The
Smith defects in the table are defects of the target-pullback exponent
lattice, not torsion in \(\operatorname{Cl}(X)\).

This is the direct interface with the multi-factor Cox calculations.  As in
those calculations, a rational rank count is not enough: the primitive
integral character must be realized geometrically.

## 4. The intrinsic Cox modification

Because the missing-boundary classes form a free basis, no auxiliary
shifting family is needed to define the divisorial Cox algebra

\[
\mathcal R_X
=
\bigoplus_{m\in\mathbb Z^r}
\mathcal O_X\left(\sum_i m_iE_i\right).
\tag{4.1}
\]

Let

\[
\widehat X=\operatorname{Spec}_X\mathcal R_X.
\tag{4.2}
\]

This relative spectrum is affine over the affine surface \(X\).  Finite
generation over \(k\) is not automatic and must not be silently assumed.

Let \(s_i\) denote the canonical homogeneous section of degree
\(\epsilon_i\).  Over \(U\), every boundary divisor is absent, so

\[
\mathcal R_X|_U
\simeq
\mathcal O_U[s_1^{\pm1},\ldots,s_r^{\pm1}]
\tag{4.3}
\]

and

\[
\widehat U
\simeq
\mathbb A^2\times\mathbb G_m^r.
\tag{4.4}
\]

The character units in (4.4) are expected Cox units.  They do not descend
to nonconstant units on \(U\).

For a target equation \(g_j\), put

\[
a_j=g_j\prod_i s_i^{-b_{ji}}.
\tag{4.5}
\]

The divisor identity

\[
\operatorname{div}_X(g_j)
=
\sum_i b_{ji}E_i+A_j
\tag{4.6}
\]

shows that \(a_j\) is a regular homogeneous Cox section with zero divisor
the affine companion.  Inside \(\mathcal R_X\) one has the exact relation

\[
\boxed{
g_j=a_j\prod_i s_i^{b_{ji}}.
}
\tag{4.7}
\]

For the rank-one packet this is the three-dimensional Cox equation

\[
\boxed{g=a\,s_E^2.}
\tag{4.8}
\]

Equation (4.8) is the correct global replacement for the principal
ramification equation in the spectator model.  The local model has a
global function cutting out \(E\); the Keller normalization has only the
primitive homogeneous section \(s_E\).

## 5. Three provisional obstruction witnesses

### 5.1 Equal-degree section and unit gate

Suppose the cusp plus its required \(2+2\) collisions produce two
nonproportional homogeneous Cox sections \(\alpha,\beta\) of the same
degree, and suppose their quotient is regular and nowhere zero on
\(\widehat U/\mathbb G_m^r=U\).  Then

\[
\alpha/\beta\in\mathcal O(U)^\times=k^\times,
\tag{5.1}
\]

a contradiction.

The missing lemma is not (5.1); that is immediate.  The missing lemma is
to construct \(\alpha,\beta\) canonically from the two conductor branches
and prove that their quotient has no affine zero or pole.  A packet ledger
which records only divisor classes cannot do this.

### 5.2 Exceptional-curve gate

The
[Davenport node-separation calculation](../extended-geometry/DAVENPORT_NODE_SEPARATING_AFFINE_MODIFICATION.md)
gives the relevant model:
two affine charts can restore both branches of a normalized boundary, but
their exceptional affine lines glue by reciprocal coordinates and form a
complete \(\mathbb P^1\).

If every Cox realization of the quartic conductor pairing forced the same
closed complete curve inside \(\widehat X\), the packet would be impossible.
Indeed an affine scheme contains no closed positive-dimensional proper
subscheme.

This calculation suggests the following strong lemma.

> **Strong Cox separation lemma.**  In the rank-one relation
> \(g=a\,s_E^2\), a clean \(3+1\) cusp together with a \(2+2\) boundary
> self-collision forces any normal primitive-\(E\) separation retaining
> both conductor branches to contain a closed exceptional
> \(\mathbb P^1\).

The existing Davenport calculation proves this only for its explicit
natural two-chart gluing.  It does not prove the universal quartic lemma.
Moreover, the later affine-modification chain shows that nonlinear
modifications can move or fill exceptional components.  The displayed
lemma is therefore a possible sufficient theorem, not the next canonical
goalpost.

### 5.3 Oriented monodromy gate

The primitive class \(E\) is the orientation missed by the doubled target
pullback.  On the target complement it is the integral refinement of the
sign/discriminant character of the \(S_4\)-monodromy.

The local cusp and matching generators are compatible with that character,
so parity alone is not an obstruction.  The target is to combine:

1. the global meridian product relation at infinity;
2. the conductor pairing of the two points in every \(2+2\) collision; and
3. extension of the primitive \(E\)-character across (4.8).

A contradiction must occur in this oriented extension problem, not in the
unoriented sheet graph.

## 6. Where three-dimensional normalization machinery applies

For the one-boundary row, \(\widehat X\) has expected dimension three when
the relevant Cox algebra is of finite type.  This is the point at which the
three-dimensional normalization and affine-modification machinery can feed
JC(2).

### 6.1 The finite target-side normalization

The Cox relation defines the finite-type affine hypersurface

\[
Z_g
=
\operatorname{Spec}
k[u,v,s,a]/(g(u,v)-a s^2).
\tag{6.1}
\]

Base-change the known finite-free surface normalization:

\[
T_g
=
B\otimes_A\mathcal O(Z_g)
\simeq
B[s,a]/(g-a s^2).
\tag{6.2}
\]

This is a finite-free \(\mathcal O(Z_g)\)-algebra of rank four.  It is a
domain with fraction field \(K(s)\): after inverting \(s\), eliminate
\(a=g/s^2\).  Since \(T_g\) is integral over \(\mathcal O(Z_g)\), both
rings have the same integral closure in \(K(s)\).

Form that canonical normalization:

\[
\overline Z_g
=
\operatorname{Norm}_{Z_g} K(s),
\qquad K=k(x,y).
\tag{6.3}
\]

Because \(Z_g\) is an affine finite-type \(k\)-scheme, this normalization is
finite.  Its generic degree is four.  After inverting \(s\),

\[
(Z_g)_s
\simeq
Y\times\mathbb G_m,
\qquad
(\overline Z_g)_s
\simeq
X\times\mathbb G_m.
\tag{6.4}
\]

The second identity follows because \(B\) is the integral closure of
\(k[u,v]\) in \(K\).  The distinguished Keller open gives

\[
U\times\mathbb G_m
\lhook\joinrel\longrightarrow
X\times\mathbb G_m
=
(\overline Z_g)_s.
\tag{6.5}
\]

This is the precise three-dimensional normalization object.  Unlike the
surface \(X\), the normal threefold \(\overline Z_g\) need not be
Cohen--Macaulay or flat over \(Z_g\); a closed-point Fitting defect on
\((s=0)\) can remain.

At the generic point of \(E\), choose a DVR parameter \(t\).  Tame simple
ramification gives

\[
g=\varepsilon t^2,\qquad \varepsilon\in B_E^\times.
\tag{6.6}
\]

The order \(T_g\) then has the local equation

\[
a s^2=\varepsilon t^2.
\tag{6.7}
\]

Its normalization adjoins the ratio \(w=t/s\); on the normalized generic
chart,

\[
t=sw,\qquad a=\varepsilon w^2.
\tag{6.8}
\]

After absorbing the unit and passing to the residue field
\(\kappa(E)\), the generic order and its normalization are the explicit
semigroup pair

\[
\kappa(E)[s,sw,w^2]
\ \subset\
\kappa(E)[s,w],
\tag{6.9}
\]

with conductor

\[
(s,t)=(s,sw)
\quad\text{in the order},\qquad
s\,\kappa(E)[s,w]
\quad\text{in the normalization}.
\tag{6.10}
\]

Thus the primitive \(E\)-character is exactly the normalization coordinate
missing from the doubled target pullback.  At the cusp and at the global
\(2+2\) pairing, the affine companion and the identification of two
distinct boundary points prevent this generic calculation from being the
whole normalization.  Those are precisely the special loci to which the
threefold conductor/Fitting analysis must be applied.

This also explains the exceptional-curve route in Section 5.2.  The
generic normalization uses the affine ratio \(w=t/s\).  A complementary
chart uses the reciprocal ratio \(q=s/t\).  If retaining the affine
companion at the cusp and both points of a \(2+2\) fiber forces both charts,
then \(wq=1\) glues their exceptional affine lines to a
\(\mathbb P^1\).  What remains unproved is exactly that both charts are
unavoidable for every normal Cox realization, rather than only for the
natural blowup.

The natural Cox morphism itself is not finite.  On \(s\ne0\) it restricts
to

\[
U\times\mathbb G_m
\longrightarrow
Y\times\mathbb G_m,
\tag{6.11}
\]

namely \(F\times\operatorname{id}\), which is still nonproper.  Instead,
the divisorial Cox algebra is the intersection of the corresponding
homogeneous valuation half-spaces in \(K(s)\), hence is integrally closed.
The integral closure defining (6.3) therefore embeds in it and gives a
comparison morphism

\[
\widehat X\longrightarrow\overline Z_g
\tag{6.12}
\]

whenever the divisorial Cox algebra is viewed inside \(K(s)\).  After
inverting \(s\), (6.12) is the open immersion in (6.5).  The essential new
problem is its behavior over \(s=0\).

A valid transfer must now prove at least one of:

1. a finite-type description of the relevant Cox/Rees chart and its map
   to \(\overline Z_g\);
2. a boundary-local comparison along \(s=0\), without claiming that
   \(\widehat X\to Z_g\) is finite; or
3. an affine modification description whose exceptional and conductor
   fibers are complete.

Once such an input is supplied, the threefold machinery can test
normality, depth, conductor, Fitting defects, affine modifications, and the
appearance of a complete exceptional curve.  The finite object to which
that machinery applies is \(\overline Z_g\), not the raw Cox morphism.

### 6.2 The canonical source-side bridge forced by \(\mathbb A^2\)

The distinguished open supplies a second finite-type threefold which is
more concrete than the intrinsic Cox algebra.  Put

\[
h=g(P,Q)\in k[x,y]
\tag{6.13}
\]

and define

\[
\begin{aligned}
C_{F,g}
&=
k[x,y,s,a]/(h-a s^2),\\
W_{F,g}
&=
\operatorname{Spec}C_{F,g}.
\end{aligned}
\tag{6.14}
\]

This is the Cartesian base change

\[
W_{F,g}
=
U\times_Y Z_g.
\tag{6.15}
\]

Hence \(W_{F,g}\to Z_g\) is étale and quasi-finite of generic degree four.
Give \(s\) degree \(1\), \(a\) degree \(-2\), and \(x,y\) degree zero.
Its degree-zero ring is exactly \(k[x,y]\).  Then (6.14) is a finite-type
graded model with

\[
(W_{F,g})_s
\simeq
\mathbb A^2\times\mathbb G_m.
\tag{6.16}
\]

The Keller condition makes \(F\) étale, so the pullback \(h\) of the
reduced target curve \(g=0\) is squarefree.  The hypersurface (6.14) is
Cohen--Macaulay.  Its singular locus is contained in

\[
\{s=0,\ h=h_x=h_y=0\},
\tag{6.17}
\]

which is a finite union of \(a\)-lines and hence has codimension at least
two.  Serre's criterion therefore makes \(W_{F,g}\) normal.  As a
hypersurface in affine four-space it is Gorenstein with trivial dualizing
module.

There are canonical inclusions in the common fraction field \(K(s)\):

\[
T_g
\lhook\joinrel\longrightarrow
\overline T_g
\lhook\joinrel\longrightarrow
C_{F,g}.
\tag{6.18}
\]

The first is normalization.  For the second, \(T_g\subset C_{F,g}\), and
every element integral over \(T_g\) is integral over the normal ring
\(C_{F,g}\), hence belongs to it.  Consequently (6.18) gives a finite-type,
\(\mathbb G_m\)-equivariant birational morphism.  It is quasi-finite because
its composition with the finite map \(\overline Z_g\to Z_g\) is the
quasi-finite base change (6.15).  A quasi-finite birational morphism to a
normal target is an open immersion.  Thus Zariski Main gives

\[
\boxed{
\Phi_{F,g}:W_{F,g}\lhook\joinrel\longrightarrow\overline Z_g.
}
\tag{6.19}
\]

After inverting \(s\), this is exactly

\[
U\times\mathbb G_m
\lhook\joinrel\longrightarrow
X\times\mathbb G_m.
\tag{6.20}
\]

The special fiber on the source side is

\[
(W_{F,g})_0
=
\operatorname{Spec}\bigl(k[x,y]/(h)\bigr)[a].
\tag{6.21}
\]

It is reduced and consists only of the affine-companion cylinders.  By
contrast, the special fiber of \(\overline Z_g\) also has the primitive
ramification component obtained generically by adjoining \(w=t/s\).
Thus (6.19) is a canonical **boundary deletion**: its special fiber omits
the primitive ramification component while retaining exactly the affine
companions.  There is no arbitrary affine modification left to choose.

This construction rules out an overbroad formulation of the Cox
separation attack.  An affine normal threefold with
\(\mathbb A^2\times\mathbb G_m\) as its \(s\ne0\) locus is already forced
by the hypothetical Keller map.  What is not automatic is compatibility
of the open boundary deletion (6.19) with the normalization
conductor and the finite-cover dualizing class at both the cusp and the
\(2+2\) connector.

### 6.3 Completed packet charts: the local conductor route is exhausted

The completed calculation can be made uniformly.  Let

\[
A_0=k[[t,r]]
\]

be the completed regular source surface at a ramification point, and let
\(\ell\in A_0\) be the smooth local equation of the other reduced pullback
branch.  Assume \(r\nmid\ell\).  The completed base-change order is

\[
O
=
A_0[[s,a]]/(a s^2-r^2\ell).
\tag{6.22}
\]

For the clean \(3+1\) cusp this equation is not merely formal bookkeeping.
The standard cubic block

\[
v=T^3+uT
\]

satisfies the exact factorization

\[
4u^3+27v^2
=(u+3T^2)^2(4u+3T^2).
\tag{6.23}
\]

Putting

\[
r=u+3T^2,\qquad \ell=4r-9T^2
\tag{6.24}
\]

gives (6.22).  The curves \(r=0\) and \(\ell=0\) have contact two.  At a
\(2+2\) connector there are two distinct source points.  At each point the
same equation holds, with \(\ell\) the other smooth branch.  Its contact
with \(r=0\) may be any finite \(m\ge1\); the calculation below is
independent of \(m\).

Adjoin

\[
z=\frac{r\ell}{s}=\frac{as}{r}.
\tag{6.25}
\]

Then the integral closure of \(O\) is

\[
\begin{aligned}
N
&=O[z]\\
&\simeq
A_0[[s,a,z]]/
\bigl(a\ell-z^2,\ as-rz,\ sz-r\ell\bigr).
\end{aligned}
\tag{6.26}
\]

The three equations are the \(2\times2\) minors of

\[
\begin{pmatrix}
a&z&r\\
z&\ell&s
\end{pmatrix}.
\tag{6.27}
\]

Here is a direct proof that (6.26) is the full normalization.  The element
\(z\) is integral and \(N\) has the same fraction field as \(O\).  The
determinantal presentation is Cohen--Macaulay by the Hilbert--Burch
theorem.  If \(\ell\) is invertible, putting \(w=z/\ell\) gives the regular
chart

\[
r=sw,\qquad z=\ell w,\qquad a=\ell w^2.
\tag{6.28}
\]

If \(r\) is invertible, then \(z=as/r\), and the smooth equation
\(\ell=as^2/r^2\) gives the other regular chart.  Since \(r=\ell=0\) has
codimension two in \(A_0\), every height-one point lies in one of these
charts.  Thus \(N\) satisfies \(R_1\); Cohen--Macaulayness gives \(S_2\),
so \(N\) is normal.

The conductor is equally explicit:

\[
\boxed{\mathfrak c=(O:N)=(r,s).}
\tag{6.29}
\]

Indeed, \(N=O+Oz\), and

\[
\operatorname{Ann}_O(z\bmod O)
=(s):(r\ell).
\tag{6.30}
\]

Modulo \(s\), the defining ideal is \((r^2\ell)\), whence

\[
(r^2\ell):(r\ell)=(r)
\]

because \(r\) and \(\ell\) are coprime.  This proves (6.29) and also

\[
N/O\simeq O/(r,s).
\tag{6.31}
\]

Since \(O\) is a Gorenstein hypersurface and \(N\) is finite
Cohen--Macaulay over it, finite duality gives

\[
\omega_N
\simeq
\operatorname{Hom}_O(N,\omega_O)
\simeq
\operatorname{Hom}_O(N,O)
\simeq
\mathfrak c
=(r,s)N.
\tag{6.32}
\]

At the packet endpoint this ideal needs two generators, so \(N\) is
non-Gorenstein there.  That is not an obstruction: the Zariski--Main source
open is locally \(D(r)\), where \((r,s)N=N\) and the canonical module is
trivial.  On the special fiber,

\[
N/(s)
=
A_0[[a,z]]/
\bigl(z^2-a\ell,\ rz,\ r\ell\bigr),
\tag{6.33}
\]

while localization at \(r\) gives

\[
(N/(s))_r
\simeq
(A_0/(\ell))_r[a].
\tag{6.34}
\]

Thus the canonical open deletes the ramification/conductor endpoint and
retains exactly the affine-companion cylinder with that endpoint removed.
The cusp and both completed branches of a \(2+2\) connector are compatible
with the same deletion.

The quartic member of the universal spectator family realizes the
\(3+1\) and \(2+2\) fibers in one finite smooth cover, so even simultaneous
completed compatibility is not a contradiction.  In that family \(r\) is
the global principal Jacobian coordinate and the distinguished open is
\(\mathbb A^1\times\mathbb G_m\).  For a Keller normalization, the local
equations \(r_i\) instead patch by units representing the primitive,
nonprincipal class \([E]\), while the degree-zero distinguished open must
be \(\mathbb A^2\).  Even this transition law is compatible.  If

\[
r_i=u_{ij}r_j,
\tag{6.35}
\]

then the global identity \(g=r_i^2\ell_i=r_j^2\ell_j\) gives

\[
\ell_i=u_{ij}^{-2}\ell_j,\qquad
z_i=u_{ij}^{-1}z_j.
\tag{6.36}
\]

In particular

\[
(r_i,s)=(r_j,s)
\tag{6.37}
\]

on overlaps.  Thus the normalizations and conductor ideals glue
canonically for every boundary class; \(z_i\) simply carries the inverse
primitive Cox character.  Nonprincipality by itself is not an obstruction.
The first unresolved datum is the global graded section algebra of the
deleted open, especially its degree-zero identity \(k[x,y]\), together
with the pairing of its cusp and connector endpoints.

### 6.4 The odd-square multiplication divisor

The grading compresses the remaining global datum further.  Give

\[
\deg s=1,\qquad \deg a=-2,\qquad \deg z=-1.
\tag{6.38}
\]

Normalization adds no degree-zero functions:

\[
\boxed{(\overline T_g)_0=B.}
\tag{6.39}
\]

Indeed, the integral closure of a graded domain is graded.  A homogeneous
degree-zero element of \(\overline T_g\) is integral over \(B=(T_g)_0\);
normality of \(B\) puts it back in \(B\).  In the completed chart the same
fact is visible monomial by monomial:

\[
a^jz^ks^{2j+k}
=(as^2)^j(sz)^k
=(r^2\ell)^j(r\ell)^k\in A_0.
\tag{6.40}
\]

By contrast, if \(R=k[x,y]\), the canonical source bridge has rank-one
free graded pieces

\[
(C_{F,g})_n=
\begin{cases}
R\,s^n,&n\ge0,\\
R\,a^m,&n=-2m,\quad m\ge1,\\
R\,a^{m+1}s,&n=-(2m+1),\quad m\ge0.
\end{cases}
\tag{6.41}
\]

Thus \((C_{F,g})_0=R\) is genuinely supplied by taking global sections on
the deleted open, not by normalization.

The first two negative pieces retain exactly the affine-companion divisor.
Locally,

\[
N_{-1}=A_0z,\qquad N_{-2}=A_0a,
\tag{6.42}
\]

and square multiplication is

\[
\mu_N:N_{-1}^{\otimes2}\longrightarrow N_{-2},
\qquad z^2=\ell a.
\tag{6.43}
\]

Hence

\[
\operatorname{coker}\mu_N\simeq A_0/(\ell).
\tag{6.44}
\]

On the source bridge,

\[
C_{-1}=R(as),\qquad C_{-2}=Ra,
\tag{6.45}
\]

and

\[
(as)^2=a^2s^2=ha,
\qquad
\operatorname{coker}\mu_C\simeq R/(h).
\tag{6.46}
\]

The global affine companion \(h=0\) is therefore not auxiliary geometry:
it is exactly the degeneracy divisor of the odd-square multiplication
between two rank-one graded pieces.  Its points at infinity carry the
contact-two cusp endpoint and the paired connector endpoints.  This is the
smallest global module datum that simultaneously remembers the
nonprincipal Cox character, the two-generated plane ring, and the packet
incidence.

### 6.5 The four degree-zero filters do not give a finite semigroup problem

There is a necessary finiteness gate before enumerating endpoint
semigroups.  The data proved in Sections 6.3--6.4 do not bound the pole
orders of \(x,y\), and the connector chart explicitly permits arbitrary
finite contact.  More strongly, factoriality and the rank-one module
conditions do not add such a bound.

> **Proposition 6.1 -- degree-zero semigroup insufficiency.**
> The following four inputs do not define a finite endpoint enumeration:
>
> 1. the completed models \(as^2=r^2\ell\), with cusp contact two and
>    arbitrary finite connector contact;
> 2. factoriality of \(R=k[x,y]\);
> 3. rank-one freeness of every graded piece of \(C_{F,g}\); and
> 4. \(\operatorname{coker}(C_{-1}^{\otimes2}\to C_{-2})=R/(h)\).
>
> At the level seen by these four inputs, infinitely many distinct
> connector pole pairs survive.  Consequently they cannot by themselves
> force a descended unit, an extra boundary component, or a monodromy
> contradiction.

#### Proof

First, the last three inputs are formal for every reduced
\(h\in k[x,y]\).  The ring \(R\) is a UFD independently of \(h\), while

\[
C_h=R[s,a]/(h-as^2)
\tag{6.47}
\]

has the normal forms

\[
(C_h)_d=
\begin{cases}
R s^d,&d\ge0,\\
R a^m,&d=-2m<0,\\
R a^{m+1}s,&d=-(2m+1)<0.
\end{cases}
\tag{6.48}
\]

Thus all pieces are free of rank one.  Moreover

\[
(as)^2=a^2s^2=ha
\tag{6.49}
\]

again gives cokernel \(R/(h)\).  None of these statements bounds a pole
order.

Here is an explicit infinite degree-zero test family.  For \(n\ge1\), put

\[
\begin{aligned}
c&=y^2-x^3,\\
d_n&=xy-y^{n+1}-1,\\
h_n&=c\,d_n.
\end{aligned}
\tag{6.50}
\]

The two factors are coprime: on the cusp normalization
\((x,y)=(t^2,t^3)\), the second factor becomes
\(t^5-t^{3n+3}-1\), which is nonzero.  The connector factor is smooth,
because \((d_n)_x=y\), and \(y=0\) is incompatible with \(d_n=0\).
Hence \(h_n\) is reduced.

The normalization of \(c=0\) has its standard endpoint pole pair

\[
(\operatorname{pole}(x),\operatorname{pole}(y))=(2,3).
\tag{6.51}
\]

The connector component is

\[
R/(d_n)
\simeq k[t,t^{-1}],
\qquad
y=t,\quad x=t^n+t^{-1}.
\tag{6.52}
\]

Indeed \(y^{-1}=x-y^n\) in \(R/(d_n)\), which gives the inverse
isomorphism.
At \(t=0\) and \(t=\infty\), its two endpoint pole pairs are respectively

\[
(1,0),\qquad(n,1).
\tag{6.53}
\]

They are distinct for all \(n\), and their row matrix is unimodular:

\[
\det
\begin{pmatrix}
1&0\\
n&1
\end{pmatrix}
=1.
\tag{6.54}
\]

Thus even the stronger *provisional* filter that asks the connector
generator vectors to span a saturated free lattice retains every \(n\).
The completed calculation does not remove them: its connector contact is
allowed to be any finite positive integer, and its formal equations contain
no pole-order constraint on the global functions \(x,y\).

This family is not asserted to lift to a quartic finite normalization or to
a Keller map.  It proves the narrower and necessary point: every condition
in the proposed four-filter degree-zero abstraction survives an unbounded
parametric family.  Therefore that abstraction is not yet a finite
semigroup problem. \(\square\)

There is a second loss of information.  The divisor \(h=0\), and hence the
odd-square cokernel, records its normalized endpoints but not the
equivalence relation saying which two endpoints have the same target value.
That pairing belongs to the finite map and its conductor, not to the
abstract \(R\)-module \(R/(h)\).  In particular, no sheet-monodromy test can
be recovered from the cokernel alone.

Factoriality must also be used with care.  A UFD need not have a saturated
numerical value semigroup for an arbitrary filtration, and normality need
not pass to its associated graded ring.  A semigroup saturation test is
therefore only a prefilter unless one first proves that the chosen
multivaluation has one-dimensional leaves and that \(x,y\) form a
Khovanskii/SAGBI basis for the relevant associated graded algebra.

There is a further coordinate-invariance obstruction which persists even
after the conductor pairing is restored.

> **Proposition 6.2 -- triangular-shear unboundedness.**
> Suppose a component of \(h=0\) has an endpoint \(p\), and let
> \[
> \alpha=\operatorname{pole}_p(x),\qquad
> \beta=\operatorname{pole}_p(y).
> \tag{6.55}
> \]
> Then the orbit of the ordered generator pair \((x,y)\) under polynomial
> automorphisms of \(k[x,y]\) has unbounded endpoint pole pairs.  This
> remains true after fixing the geometric degree, completed packet,
> conductor pairing, and sheet monodromy of a hypothetical Keller map.

#### Proof

After harmless constant translations, neither coordinate restricts to zero
on the component; such translations do not change its endpoint poles.  At
least one of \(\alpha,\beta\) is positive, since a missing point of the
normalization of an affine plane curve cannot have both affine coordinates
regular.  If \(\alpha>0\), use the triangular coordinates

\[
x_m=x,\qquad y_m=y+x^m.
\tag{6.56}
\]

For every \(m\) with \(m\alpha>\beta\), the two summands have unequal
orders at \(p\), so cancellation is impossible and

\[
\operatorname{pole}_p(y_m)=m\alpha.
\tag{6.57}
\]

If \(\alpha=0\), then \(\beta>0\), and the opposite shear
\((x+y^m,y)\) gives the same conclusion.  Both shears have Jacobian one
and polynomial inverses.

The shear extends to the bridge by fixing \(s,a\) and sending
\(h\) to its pullback \(h_m\):

\[
k[x,y,s,a]/(h-as^2)
\simeq
k[x,y,s,a]/(h_m-as^2).
\tag{6.58}
\]

If \(F\) were a quartic Keller map, precomposition by the same source
automorphism would preserve its constant Jacobian, function-field degree,
finite normalization up to isomorphism, completed conductor packets, and
monodromy.  Only the chosen polynomial generators and their numerical pole
pairs change.  Hence no bound on raw generator pole pairs can follow from
those invariant data. \(\square\)

For the family (6.50), the calculation is visible without any hypothetical
cover.  On the connector normalization,

\[
x=t^n+t^{-1},\qquad
y_m=t+(t^n+t^{-1})^m.
\tag{6.59}
\]

The endpoint pole pairs become

\[
(1,m)\quad\text{at }t=0,\qquad
(n,mn)\quad\text{at }t=\infty.
\tag{6.60}
\]

Thus the earlier determinant-one matrix is a useful witness, but not an
intrinsic factoriality condition: a polynomial coordinate change turns it
into a determinant-zero pole matrix while the ambient ring remains the
same UFD.  This is the numerical refinement of the
[stable affine-membership lemma](../extended-geometry/INTRINSIC_SELECTOR_ATTACK.md#4-stable-mixing-cannot-repair-the-failure):
the existence of a pole is intrinsic, whereas its magnitude in a chosen
coordinate generator is not.

The coordinate-invariant replacement is the marked multivaluation of the
whole ring.  Componentwise, for endpoint valuations
\(\nu_i=\operatorname{ord}_{p_i}\), set

\[
\Gamma(R;\nu_1,\ldots,\nu_s)
=
\left\{
(\nu_1(f),\ldots,\nu_s(f)):
0\ne f\in R/(h_j)
\right\}
\subset\mathbb Z^s,
\tag{6.61}
\]

where \(h_j\) is the relevant irreducible component.  Changing polynomial
coordinate generators permutes the functions being measured but not
\(\Gamma\).  A finite generator-pair compiler can equivalently use the
automorphism-minimized pole height

\[
H(\nu)
=
\min_{\substack{R=k[u,v]\\ \bar u,\bar v\ne0\text{ in }R/(h_j)}}
\sum_{i=1}^s
\left(
\max\{0,-\nu_i(u)\}
+
\max\{0,-\nu_i(v)\}
\right).
\tag{6.62}
\]

The minimum exists because the values are nonnegative integers and the
original pair \((x,y)\) is admissible.  It is an invariant of the marked
valuation orbit.  What is not yet proved is a degree-four bound on
\(H(\nu)\), finite generation of the relevant marked part of \(\Gamma\),
or an algorithm producing a minimizing coordinate pair.

For the connector component in (6.52), the minimizing calculation is
already exact:

\[
u=x-y^n=t^{-1},\qquad v=y=t.
\tag{6.63}
\]

Its two pole rows are \((1,0)\) and \((0,1)\), so \(H=2\).  The opposite
inequality is automatic because each of the two missing normalization
points must be a pole of at least one coordinate.  Thus the parameter \(n\)
is completely removed from the marked connector component by
automorphism minimization.  This demonstrates both why raw enumeration
overcounts and what a successful reduction step should output.

Numerical values alone still do not implement that reduction.

> **Proposition 6.3 -- the endpoint residues are necessary.**
> Two marked two-ended plane curves can have the same pole matrix for their
> displayed generators while a given triangular shear cancels the leading
> pole at both endpoints for one and not for the other.

Take \(c\in k^\times\) and, on \(k[t,t^{-1}]\), put

\[
X=t+t^{-1},\qquad
Y_c=t^2+c\,t^{-2}+t.
\tag{6.64}
\]

These functions define a birational plane curve with normalization
\(\mathbb G_m\).  Indeed \(t\) is integral over \(k[X]\) by
\(t^2-Xt+1=0\), while

\[
t=
\frac{Y_c-cX^2+1+c}{1+(1-c)X}
\tag{6.65}
\]

proves equality of fraction fields.  At both \(t=0\) and \(t=\infty\), the
pole row of \((X,Y_c)\) is \((1,2)\), independently of \(c\).
Nevertheless

\[
Y_c-\lambda X^2
\tag{6.66}
\]

loses its order-two pole at \(t=\infty\) exactly when \(\lambda=1\), and
loses it at \(t=0\) exactly when \(\lambda=c\).  A single quadratic shear
reduces both endpoints precisely for \(c=1\).

In general, whenever
\(\operatorname{pole}_{p_i}(y)=m\operatorname{pole}_{p_i}(x)\), the
relevant datum is the initial-residue ratio

\[
\rho_i=
\frac{\operatorname{in}_{\nu_i}(y)}
     {\operatorname{in}_{\nu_i}(x)^m}
\in\kappa(p_i)^\times.
\tag{6.67}
\]

The shear \(y-\lambda x^m\) cancels at every tied endpoint if and only if
all \(\rho_i\) are the image of the same \(\lambda\in k\).  Therefore an
effective automorphism reduction needs the residue-decorated associated
graded algebra, or at least these simultaneous initial-form ratios, not
only the numerical semigroup \(\Gamma\).

The residue gate now has a finite exact front end.

> **Proposition 6.4 -- finite monotone shear compiler.**
> For a fixed pair of nonzero Laurent polynomials \(u(t),v(t)\), the set of
> monomial triangular shears
> \[
> v\longmapsto v-c\,u^m
> \tag{6.68}
> \]
> which strictly decrease the sum of the pole orders at \(0,\infty\) is
> finite and is determined exactly by the pole orders and leading residues.
> The same holds in the opposite orientation.  Iterating only strict
> decreases terminates.

Indeed, a decreasing shear must cancel a current leading pole at at least
one endpoint.  At such an endpoint,

\[
\operatorname{pole}(v)=m\operatorname{pole}(u),
\qquad
c=\frac{\operatorname{in}(v)}{\operatorname{in}(u)^m}.
\tag{6.69}
\]

There are only two endpoints, so these equalities produce a finite
candidate list.  Direct Laurent expansion tests the total height after
each candidate shear, including any increase at the other endpoint.
Every accepted step lowers a nonnegative integer, proving termination.

The implementation
[`cas/endpoint_valuation_compiler.py`](cas/endpoint_valuation_compiler.py)
enumerates all such moves in both orientations and returns a complete
strict-descent certificate.  On (6.52) it finds
\(x\mapsto x-y^n=t^{-1}\) and height two.  Starting from the expanded
coordinates (6.59), it first removes \(x^m\) from \(y_m\) and then removes
\(y^n\) from \(x\).  On (6.64), its terminal monotone heights are two for
\(c=1\) and four for the checked \(c\ne1\) rows, exactly reflecting the
residue mismatch.

This is not yet a canonical normal-form theorem.  A general polynomial
shear can contain several powers whose combined effect is not represented
by a strictly decreasing monomial prefix.  The following extension closes
that single-shear gap.

> **Proposition 6.5 -- finite complete-polynomial shear compiler.**
> Modulo terms which do not affect either endpoint pole, all triangular
> polynomials \(P\) for which
> \[
> (u,v)\longmapsto(u,v-P(u))
> \tag{6.70}
> \]
> strictly lowers the two-endpoint pole height form a finite,
> residue-determined list.  They can be enumerated even when a proper
> leading prefix of \(P\) preserves or increases the height.

Process the nonconstant terms of \(P\) in descending degree.  Any term
which contributes to a final decrease must, at the stage when it is
processed, cancel a current leading pole at some endpoint.  Its degree and
coefficient are therefore forced by (6.69).  After that cancellation, the
next relevant degree is strictly smaller.  Consequently the search tree is
finite.  Terms which never cancel a current leading pole can be deleted:
they either leave both endpoint poles unchanged or introduce a pole which
lower-degree terms cannot remove.  This proves completeness after removing
endpoint-irrelevant terms.

The need to retain nondecreasing prefixes is real.  Put

\[
u=t^{-1}+t^2,\qquad
v=u^3-t^6+t^5+u^2.
\tag{6.71}
\]

The pole rows of \((u,v)\) are \((1,3)\) and \((2,5)\), of total height
eleven.  No monomial shear strictly lowers this height.  The first
cancellation is neutral:

\[
v-u^3=-t^6+t^5+u^2,
\tag{6.72}
\]

whose target pole sum is still eight.  The lower term then gives

\[
v-u^3-u^2=-t^6+t^5,
\tag{6.73}
\]

so the total height drops to nine.  The expanded checker verifies that the
polynomial compiler finds exactly this descending-degree certificate.

The implementation now iterates the best complete lowering polynomial
shear in either orientation.  Its terminal form admits no
height-decreasing triangular polynomial in either coordinate.  Alternating
directions do not create a further peak obstruction.

> **Proposition 6.6 -- marked multi-pole peak reduction for reduced Jung
> words.**
> Let \(p_1,\ldots,p_s\) be finitely many marked boundary valuations and put
> \[
> H(u,v)=\sum_i\bigl(\operatorname{pole}_{p_i}(u)
>                  +\operatorname{pole}_{p_i}(v)\bigr).
> \]
> Consider an alternating word of complete triangular factors
> \[
> z_{\epsilon_j}\longmapsto
> z_{\epsilon_j}+P_j(z_{1-\epsilon_j}),
> \qquad
> \epsilon_{j+1}=1-\epsilon_j,
> \tag{6.74}
> \]
> with every \(d_j=\deg P_j\ge2\).  If the whole word lowers \(H\), then its
> first complete triangular factor lowers \(H\).  In particular, a pair
> terminal for complete triangular reductions cannot be lowered by a
> reduced alternating Jung word.

#### Proof

It is enough to compare two consecutive factors at one marked valuation.
Write \(a,b\) for the pole orders before the first factor,

\[
x'=x+P(y),\qquad y'=y+Q(x'),
\tag{6.75}
\]

and write \(a',b'\) for the pole orders of \(x',y'\).  Set
\(d=\deg P\ge2\), \(e=\deg Q\ge2\), and

\[
\alpha=a'-a,\qquad \beta=b'-b.
\tag{6.76}
\]

We claim

\[
\boxed{\ \beta\ge\alpha.\ }
\tag{6.77}
\]

If \(\beta<0\), the leading pole of \(Q(x')\) must cancel that of \(y\).
Hence \(b=e a'\).  Since \(db=de a'>a'\), obtaining the smaller pole
\(a'\) in the first equation forces \(a=db\).  Therefore

\[
a-a'=(de-1)a',
\qquad
b-b'\le e a'\le(de-1)a',
\tag{6.78}
\]

which is (6.77).

Suppose now that \(\beta\ge0\).  If \(\alpha\le0\), (6.77) is immediate.
If \(\alpha>0\), the leading term of \(P(y)\) dominates \(x\), so
\(a'=db>a\).  Then \(ea'=edb>b\), the leading term of \(Q(x')\) dominates
\(y\), and \(b'=ea'\).  Thus

\[
\beta=(ed-1)b\ge db\ge db-a=\alpha,
\tag{6.79}
\]

because \(ed-1\ge d\).  This proves (6.77) in every case.  Notice that tied
initial residues determine how far a pole drops, but the inequality itself
does not require their values.

For the full alternating word, let \(\Delta_{j,i}\) be the change of the
pole order of the coordinate modified by factor \(j\) at \(p_i\).
Applying (6.77) to every consecutive pair gives

\[
\Delta_{1,i}\le\Delta_{2,i}\le\cdots\le\Delta_{r,i}.
\tag{6.80}
\]

After summing over the marked valuations, the factor-height increments
\(\Delta_j=\sum_i\Delta_{j,i}\) are nondecreasing.  If
\(\sum_j\Delta_j<0\), then necessarily \(\Delta_1<0\).  The first complete
factor therefore lowers the current marked multi-pole peak. \(\square\)

The conductor equivalence relation and residue characters are carried
unchanged by the coordinate automorphism.  Proposition 6.6 shows that they
are not needed for the inter-factor peak inequality; they remain necessary
for compiling the actual complete cancellation polynomial and for the
subsequent conductor/monodromy filters.

Affine factors create no separate local peak obstruction.  Translations,
nonzero diagonal scalings, and coordinate exchange preserve \(H\).  For a
linear change
\[
(u,v)\longmapsto(au+bv,cu+dv),
\qquad ad-bc\ne0,
\]
suppose no linear shear of \(u\) by \(v\), or of \(v\) by \(u\), lowers the
corresponding pole degree.  If \(a,d\ne0\), the two new coordinates have
pole degrees at least those of \(u,v\), respectively.  If \(b,c\ne0\), use
the off-diagonal pivots and obtain the same lower bound in the opposite
order.  One of these two pivot pairs is nonzero because the determinant is
nonzero.  Hence a lowering affine basis change already exposes a lowering
complete linear shear.

As bounded evidence only, the checker exhausts all 49 ordered Laurent pairs
whose two supports are nonempty subsets of
\(\{-1,0,1\}\), followed by all alternating two-step monomial shears of
degrees one or two and coefficients \(\pm1\).  Sixteen paths lower the
final height after a nondecreasing first step.  In every one, the initial
pair already admits a lowering complete polynomial shear, so none is a
terminal peak counterexample.  This finite grid is now a regression for
Proposition 6.6, not part of its proof.

The endpoint search becomes genuinely finite only after adjoining the
following global input:

1. a fixed compactification and the finite endpoint set;
2. the marked multivaluation semigroup \(\Gamma\) together with its
   initial-residue data, or a proved reduction to an automorphism-minimal
   coordinate pair;
3. a degree-four bound on the invariant height \(H(\nu)\), rather than on
   an arbitrary presentation of \(x,y\);
4. the conductor equivalence relation pairing connector endpoints; and
5. the global degree-four monodromy/meridian relation carried by that
   pairing.

With an invariant height bound \(B\) and a fixed minimizing coordinate
pair, the vectors
\[
\alpha=(\operatorname{pole}_{p_i}x)_i,\qquad
\beta=(\operatorname{pole}_{p_i}y)_i
\tag{6.81}
\]
lie in a finite box.  Local contact, principal-divisor, graded-module, and
monodromy tests can then be compiled as finite filters.  Proposition 6.6
closes the reduced alternating-direction peak step.  Deriving \(B\) and the
pairing map on the automorphism-reduced valuation data, rather than adding
more Laurent layers or running an unbounded raw-pair search, is the next
missing theorem.

The exact regression is
[`cas/test_degree_zero_endpoint_pairing.py`](cas/test_degree_zero_endpoint_pairing.py).
It checks the family (6.50)--(6.54) for representative \(n\) and the graded
normal forms (6.48).  Infinitude is the displayed uniform proof, not the
bounded loop.

### 6.6 What the later Cox attacks change

The following lessons from the multi-factor programme now apply directly:

1. [affine or vector-bundle stabilization](../extended-geometry/COX_VECTOR_BUNDLE_STABILIZATION_OBSTRUCTION.md)
   cannot remove a genuine unit or motivic defect, but nonlinear affine
   modification can change the boundary model;
2. a unit-rank mismatch on a fixed center can
   [disappear in a total family](../extended-geometry/DAVENPORT_POST_COORDINATE_ATTACKS.md);
3. reciprocal conductor charts produce an exceptional \(\mathbb P^1\) in
   their
   [natural completion](../extended-geometry/DAVENPORT_NODE_SEPARATING_AFFINE_MODIFICATION.md),
   but that curve is not invariant under every higher-dimensional surgery;
   and
4. an affine-space source with constant residue Jacobian can coexist with
   dicritical components when the target is a non-affine
   [Cox](../extended-geometry/AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md) or
   [Danielewski](../extended-geometry/SMOOTH_DANIELEWSKI_TRIPLE_ROOT_MAP.md)
   threefold.

Therefore a valid quartic obstruction must retain all four pieces which
those countermodels drop:

- the finite degree-four target normalization \(\overline Z_g\);
- the canonical graded overring \(C_{F,g}\) with degree-zero ring
  \(k[x,y]\);
- the ramification different/canonical character on the omitted component;
  and
- the completed conductor incidence at the \(3+1\) and \(2+2\) points.

The two-boundary rows do not add a second missing primitive character.
Their target vectors already account for the \(D\)-direction up to a
rank/saturation completion, while \(\epsilon_E\) completes the lattice in
both cases.  Geometrically eliminating or quotienting the unramified
direction still requires proof, but the integral obstruction remains the
same primitive ramification character.

## 7. Revised closure target

The global quartic target can now be stated without a false local premise:

> **Global quartic closure target.**  No normal affine surface \(X\), finite
> free of rank four over \(\mathbb A^2\), can have a distinguished open
> \(U=X\setminus(E\cup D)\simeq\mathbb A^2\) whose boundary and special
> fibers realize either quartic Orevkov packet.

Here \(D\) is omitted in the one-boundary row.  Section 6.3 shows that a
completed special-fiber deletion theorem is false: the cusp chart and both
connector charts have the same normal determinantal overring, conductor,
and dualizing module, and the source open is exactly their Gorenstein
locus \(D(r)\).  Equations (6.35)--(6.37) further show that the
nonprincipal transition cocycle glues without a defect.  The first sharp
lemma must therefore see the global sections of the deleted open:

> **Degree-zero endpoint-pairing theorem.**  For a degree-four
> finite-normalization packet with one ramified boundary \(E\), the
> canonical graded Zariski--Main open
> \(\Phi_{F,g}:W_{F,g}\to\overline Z_g\), locally obtained by the
> compatible deletions \(D(r_i)\), cannot have
> \((C_{F,g})_0=k[x,y]\) while its boundary valuation data realize both a
> clean \(3+1\) endpoint and a \(2+2\) connector in the degeneracy divisor
> of the odd-square multiplication map (6.46).

This formulation uses every global input which the spectator models lack:

- \(U^\times=k^\times\);
- \(\operatorname{Cl}(X)=\mathbb Z[E]\);
- the primitive canonical/different character \([E]\);
- the exact Cox relation \(g=a\,s_E^2\);
- the required \(2+2\) connector;
- the normal Gorenstein source bridge with degree-zero ring \(k[x,y]\);
  and
- the finite degree-four target normalization.

The completed rings, normalization conductor, reflexive dualizing
comparison, and transition law are now computed in (6.22)--(6.37), and all
are compatible.  Equations (6.39)--(6.46) further identify the first
global module carrying packet data: the square of the degree \(-1\)
rank-one piece maps to degree \(-2\), with cokernel \(k[x,y]/(h)\).  The
negative audit in Proposition 6.1 shows that this does not yet give a
finite semigroup enumeration.  The immediate task is to derive, from the
degree-four cover and its log boundary, a bound on the
automorphism-minimized multivaluation height (6.62), together with the
conductor equivalence relation pairing the connector endpoints.
Proposition 6.2 rules out a bound on arbitrary coordinate generators.
Only after adjoining the invariant data can one test whether
factoriality and \((C_{F,g})_0=k[x,y]\) force either:

1. trivialization of \([E]\), hence a nonconstant degree-zero unit on the
   purported \(\mathbb A^2\);
2. an additional homogeneous boundary component; or
3. a global monodromy/conductor pairing incompatible with the connector.

These are possible conclusions of the endpoint-pairing theorem, not
assumptions.
The coarse class of the conductor and the local non-Gorenstein endpoint
alone cannot supply the contradiction.

If the degree-zero endpoint-pairing theorem is proved, the clean one-boundary
packet is closed.  The same primitive-character audit then supplies the
correct integral input for the ramified-plus-unramified row.

The conductor-decorated semigroup continuation is
[`QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md`](QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md).
It computes the cusp and connector semigroups, restores arbitrary endpoint
pairings, and formulates the exact bounded lattice problem.  Its
no-finiteness theorem shows that conductor decoration alone still does not
bound the connector count, the automorphism-minimized marked height, or the
connector contact orders.  Thus it is a strengthened negative gate, not the
degree-zero endpoint-pairing theorem or a quartic exclusion.

## 8. Reproduction

Run

```bash
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
.venv/bin/python plane-jc/cas/test_degree_zero_endpoint_pairing.py
.venv/bin/python plane-jc/cas/experiment_quartic_endpoint_semigroups.py \
  --max-connectors 4 --max-pole 8 --max-contact 8 --cutoff 12 \
  --output artifacts/generated-results/quartic_endpoint_semigroups.json
Singular -q plane-jc/cas/quartic_completed_deletion.sing
```

In addition to the existing Orevkov, conductor, Euler, and monodromy
checks, the regression verifies:

1. the three target groupings;
2. their exact free/torsion cokernels;
3. the affine companion class vectors;
4. saturation after adjoining \(\epsilon_E\);
5. equality of \(\epsilon_E\) with the canonical class; and
6. the negative conclusions that neither a unit nor a class-group
   conductor obstruction is present at the coarse level; and
7. the cusp factorization, determinantal normalization relations, and
   monomial colon calculation giving conductor \((r,s)\) in both completed
   packet charts; and
8. the compatible transition exponents
   \((r,\ell,z,s,a)=(1,-2,-1,0,0)\); and
9. the degree \(0,-1,-2\) bridge and its affine-companion square
   cokernels \(A_0/(\ell)\) and \(k[x,y]/(h)\); and
10. the infinite packet-shaped degree-zero family with endpoint data
    \((2,3)\) and \((1,0),(n,1)\), including the unimodular connector
    matrix; and
11. the triangular-shear orbits
    \((1,m),(n,mn)\), which prove that raw generator pole pairs are not
    coordinate-invariant; and
12. equal numerical pole matrices with unequal endpoint residue ratios,
    which prove that simultaneous shear reduction needs initial-form data;
    and
13. exhaustive strict-height monomial shear reduction, including its
    terminating step certificates and residue-sensitive terminal rows; and
14. exhaustive lowering polynomial shears with nondecreasing prefixes,
    including the height-eleven to height-nine two-term witness (6.71);
    and
15. the bounded 49-pair alternating-peak grid and its sixteen
    non-counterexample lowering paths.

The normality and finite-duality conclusions are the written proofs in
Section 6.3.  The independent Singular replay computes the cusp
normalization and conductor from the hypersurface order, compares the
result with (6.26), and checks that the determinantal overring is normal.

The checker is pinned in `MATH_STATUS.json` with

```text
sha256:57b6cb95310985e30b48afacefc0e803253f37b2aa5d441d208930484c91f38a
```

computed by

```bash
shasum -a 256 plane-jc/cas/test_plane_boundary_exclusion.py
```

The Singular replay is pinned separately in `MATH_STATUS.json` under
`QLD4` with

```text
sha256:cb3cc6a7906b415ab73424398202883ff03c416cda63307f540295d9be29ab4d
```
