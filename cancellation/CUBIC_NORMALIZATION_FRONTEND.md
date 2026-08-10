# Cubic normalization and the affine-linear frontend

This note moves the degree-three classification problem one step upstream.
Instead of assuming a coordinate-preserving suspension, it starts with the
canonical finite normalization of a cubic Keller map.  It proves that the
only possible failure of flatness is zero-dimensional, extracts the
Deligne--Faddeev binary cubic when that defect is absent, and classifies the
case in which the resulting coefficient morphism is affine-linear.

The result does **not** make the proposed minimal-boundary gateway
operational or prove the resulting classification conjecture.  It replaces
one part of the undifferentiated phrase “extract a suspension” by two
concrete intrinsic obligations:

1. remove the zero-dimensional flatness defect of the finite normalization;
2. prove that its binary-cubic coefficient morphism is affine-linear (or
   prove a nonlinear replacement for the hyperplane-slice theorem).

The local defect search is now stratified by the formal-gauge theorem in
[`UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md`](UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md):

1. for a **smooth cubic symbol**, universal cotangent saturation of the full
   24-parameter quartic family is proved.  No further quartic saturation
   search belongs to this stratum; the next gates are global algebraization
   of the formal gauge and compatibility with the marked boundary and Keller
   open;
2. for the six **singular squarefree symbols**, genuine quartic nongauge
   moduli survive, of dimensions \(2,4,4,6,6,8\).  The complete
   deterministic nongauge complements now pass cotangent saturation, but
   retain the same parameter-independent multiplicity-six support-hull
   obstruction.  Thus quartic `C2` is closed on these complements while
   `C1` fails there; Keller geometry must exclude the models rather than
   deepen their collision ideal.  The rows are organized by their exact
   annihilators
   \[
   (x),\ (x^2),\ (yz),\ (y^3),\ (xyz),\ (x^3);
   \]
3. for the **double-line, triple-line, and zero symbols**, the gauge
   cokernels are faithful.  Generically étale and Keller compatibility must
   be tested before any further module-saturation calculation.

Work over an algebraically closed field `k` of characteristic zero.

## 1. The canonical cubic algebra

Let

\[
 F:U=\mathbb A^3\longrightarrow Y=\mathbb A^3
\]

be a dominant quasi-finite Keller map of geometric degree three.  Put

\[
 A=k[Y],\qquad L=k(U),\qquad
 B=\operatorname{Norm}_A(L),
\tag{1.1}
\]

and write

\[
 \pi:\bar X=\operatorname{Spec}B\longrightarrow Y
\]

for the canonical finite normalization.  The distinguished map
`U -> bar X` is the Zariski--Main open immersion.

The `A`-module `B` is finite, torsion-free, and has generic rank three.
It need not be locally free merely because `B` is normal: normality is an
`S_2` condition, while flatness over the regular threefold `Y` is a
maximal-Cohen--Macaulay condition.

### Proposition 1.1 -- flatness is a codimension-three question

The finite morphism `pi` is flat over every point of `Y` of codimension at
most two.  Consequently its nonflat locus is a finite zero-dimensional
closed subset of `Y`.

Equivalently, since `B` has generic rank three,

\[
 Z_{\mathrm{flat}}
 =
 V\!\left(\operatorname{Fitt}^{A}_{3}(B)\right)
\tag{1.2}
\]

is zero-dimensional, and `pi` is finite flat exactly when
`Z_flat` is empty.

#### Proof

Let `p in Spec A` have height at most two and let `q` be a prime of `B`
above it.  Integrality and equality of fraction-field dimensions give

\[
 \dim B_q=\dim A_p\le2.
\]

The local ring `B_q` is normal.  A normal local ring of dimension at most
two is Cohen--Macaulay.  Since `A_p` is regular and the map is finite,
miracle flatness makes `B_q` flat over `A_p`.  A finite flat module is
locally free, proving the first assertion.

For a finite module of generic rank three, the locally free rank-three
locus is the complement of `V(Fitt_3(B))`; the lower Fitting ideal
`Fitt_2(B)` is zero because it vanishes after tensoring with `Frac(A)` and
`A` is a domain.  The first assertion puts the support of (1.2) in
codimension three.  Since `A` has dimension three, that support is finite.
QED

This is important for the minimal-boundary program.  Height-one valuations,
unit lattices, and divisor ledgers cannot by themselves see
`Z_flat`.  The finite normalization itself does see it through (1.2), so
flatness is an intrinsic finite-stratum condition rather than a new
coordinate marking.

### Proposition 1.2 -- a Cartier boundary removes the point defect

Suppose every point of `bar X\setminus U` has a neighborhood on which the
boundary is an effective Cartier divisor `(t=0)` and `B/(t)` is
Cohen--Macaulay.  Then `Z_flat` is empty.

In particular, it is enough that the canonical boundary be Cartier and
normal as a surface.

#### Proof

The map is already étale, hence flat, at every point of `U`.  Let `q` be a
boundary point.  The Cartier equation `t` is a nonzerodivisor in the
three-dimensional local domain `B_q`.  Since `B_q/tB_q` is
Cohen--Macaulay of dimension two, a regular sequence of length two modulo
`t` lifts to a regular sequence of length three beginning with `t`.
Therefore `B_q` is Cohen--Macaulay.  Proposition 1.1, or miracle flatness
directly, makes it flat over the corresponding regular target local ring.

A normal surface is Cohen--Macaulay, giving the last assertion.  QED

This criterion shows exactly what a strengthened saturation theorem would
need to prove.  Generic DVR saturation is not enough: it must extend the
boundary equation as a Cartier class through the finitely many closed
points and control the scheme-theoretic surface there.

### Proposition 1.3 -- the binary-cubic package

Assume `Z_flat` is empty.  Then:

1. `B` is a finite locally free `A`-algebra of rank three;
2. the trace splitting
   \[
   B=A\oplus M,\qquad
   M=\ker(\operatorname{Tr}_{B/A}),
   \tag{1.3}
   \]
   makes `M` a projective module of rank two;
3. `M` is free by Quillen--Suslin;
4. after choosing a basis of `M`, the cubic algebra is represented by a
   binary cubic
   \[
   f_Y(U,V)
   =aU^3+bU^2V+cUV^2+dV^3,
   \qquad a,b,c,d\in A,
   \tag{1.4}
   \]
   uniquely up to the natural `GL_2(A)` action.

Thus the normalization canonically supplies, up to change of
Tschirnhausen basis, a coefficient morphism

\[
 \kappa_F:Y\longrightarrow
 \operatorname{Sym}^3(k^2)\simeq\mathbb A^4.
\tag{1.5}
\]

#### Proof

Finite flatness gives (1).  Since `3` is invertible, the unit map
`A -> B` is split by `(1/3)Tr`, giving (1.3).  Direct summands of a
finite locally free module are projective, so `M` has rank two.
Quillen--Suslin makes it free over the polynomial ring `A`.

The Deligne--Faddeev correspondence over an arbitrary base identifies
rank-three locally free commutative algebras with twisted binary cubic
forms on the rank-two quotient by the unit.  A basis of `M` untwists the
form and gives (1.4); changing the basis gives exactly the stated
`GL_2(A)` action.  QED

The choice of basis in (1.4) is not intrinsic.  The pair `(M,f_Y)` and its
polynomial `GL_2(A)` orbit, the flatness-defect scheme (1.2), and the
discriminant divisor are intrinsic.  The affine span of the four displayed
coefficients in one basis is **not** intrinsic: a target-dependent
Tschirnhausen shear can make a linear representative look nonlinear.  The
invariant condition used below is the existence of an affine-linear
full-rank representative somewhere in the polynomial `GL_2(A)` orbit.

### Proposition 1.4 -- the cubic critical-divisor degree budget

Let `E` be a boundary prime of the canonical normalization which is ramified
over its target prime divisor `D`.  Then the generic DVR decomposition over
`D` consists of exactly:

1. `E`, with ramification index two and residue degree one;
2. one prime contained in the Keller open `U`, with ramification index and
   residue degree both one.

In particular, no additional boundary prime lies over `D`.

#### Proof

Let `g` be a prime equation of `D` in `k[Y]`.  Dominance makes
`g\circ F` a nonconstant polynomial on `U=A^3`, hence a nonunit.  Any
height-one prime factor of `(g\circ F)` maps dominantly to `D`: its source
divisor has dimension two and the restriction of the quasi-finite map still
has zero-dimensional fibers.  This supplies an affine prime above `D`.
The Keller condition makes it unramified, so its contribution `ef` to the
generic DVR degree sum is at least one.

The critical boundary prime contributes at least two.  Since the function
field degree is three,

\[
 \sum_{\mathfrak q\mid D}e(\mathfrak q/D)f(\mathfrak q/D)=3.
\]

Both lower bounds are therefore equalities.  The critical contribution is
`2*1`, the affine contribution is `1*1`, and there is no room for another
prime.  QED

Consequently, under the hypothesis that there is only one critical boundary
prime, every additional boundary prime is unramified and maps to a target
divisor distinct from `D`.  The extra-simple-boundary obstruction in degree
three is therefore exactly the existence of a second, purely unramified
nonproperness divisor; it cannot hide over the critical discriminant.

### Proposition 1.4a -- a point defect lies over the non-SNC branch locus

Let \(D_{\mathrm{br}}\) be the reduced divisorial branch locus of the
canonical cubic normalization.  Then

\[
 \boxed{Z_{\mathrm{flat}}\subseteq
 \operatorname{NonSNC}(D_{\mathrm{br}})
 \subseteq\operatorname{Sing}(D_{\mathrm{br}}).}
\tag{1.5a}
\]

Equivalently, the normalization is finite flat over every smooth point of
the branch divisor and over every simple-normal-crossing point.  In
particular, an SNC branch divisor closes Certificate E outright.

#### Proof

Let \(p\) be an SNC point and pass to the strict henselization \(R\) of the
regular local ring \(\mathcal O_{Y,p}\).  There are regular parameters
\(t_1,\ldots,t_r\), with \(1\le r\le3\), such that the branch equation is
\(t_1\cdots t_r=0\).  Away from that equation, purity of the branch locus
makes the normalization finite étale.  Since the characteristic is zero,
the cover is tame along every \(t_i=0\).

The tame fundamental group of
\(\operatorname{Spec}R[1/(t_1\cdots t_r)]\) is the product of \(r\)
procyclic tame inertia groups.  Proposition 1.4 says that every nontrivial
inertia generator acts as a transposition on the three sheets.  These
generators commute, while two distinct transpositions in \(S_3\) do not.
They therefore act by the same transposition.  The cover of the complement
is the disjoint union of one quadratic orbit and one trivial orbit.  The
tame local structure theorem, equivalently Abhyankar's lemma, identifies
their normalizations over \(R\) with

\[
 R[s]/(s^2-t_1\cdots t_r)\qquad\text{and}\qquad R
\]

after absorbing a unit into one parameter.  Both are finite free over \(R\).
Flatness descends through the faithfully flat strict-henselization map, so
the original normalization is flat at \(p\).  QED

This is a genuinely codimension-three restriction, not a divisor-ledger
proof of global flatness.  For a singular discriminant it only says that
the two saturation quotients of Propositions 1.15--1.17 can be supported at
closed points of its non-SNC locus.  Those remaining points still require
the Keller open or an equivalent sheet-intersection argument.

### Proposition 1.4b -- an ordinary cusp also removes the point defect

Let \(p\in D_{\mathrm{br}}\) be an ordinary cusp point: after strict
henselization and completion, the pair \((Y,D_{\mathrm{br}})\) is

\[
 R=k[[u,v,w]],\qquad
 D_{\mathrm{br}}=(4u^3+27v^2=0),
\tag{1.5b}
\]

with no other branch component through \(p\).  Then the canonical cubic
normalization is finite flat at \(p\).

Consequently

\[
 \boxed{
 Z_{\mathrm{flat}}\subseteq
 D_{\mathrm{br}}\setminus
 \bigl(D_{\mathrm{br}}^{\mathrm{SNC}}
       \cup D_{\mathrm{br}}^{\mathrm{oc}}\bigr),
 }
\tag{1.5c}
\]

where \(D_{\mathrm{br}}^{\mathrm{oc}}\) is the ordinary-cusp locus.  Thus a
branch divisor having only SNC and ordinary-cusp singularities closes
Certificate E.

#### Proof

The local complement of the ordinary cusp has profinite fundamental group
the profinite completion of

\[
 B_3=\langle \sigma_1,\sigma_2\mid
 \sigma_1\sigma_2\sigma_1=
 \sigma_2\sigma_1\sigma_2\rangle,
\tag{1.5d}
\]

and both generators are meridians.  The degree-three valuation budget of
Proposition 1.4 makes the image of each meridian a transposition in \(S_3\).
There are only two cases up to relabeling.

If the two transpositions agree, the permutation representation has orbits
\(2+1\).  The corresponding normalization is

\[
 R[Z]/(Z^2-4u^3-27v^2)\ \oplus\ R.
\tag{1.5e}
\]

The quadratic summand is normal: after a linear change over the
algebraically closed coefficient field it has the \(A_2\) form
\(XY=4u^3\).  It is free of rank two over \(R\).

If the transpositions are distinct, they generate \(S_3\) and act
transitively.  The corresponding three-sheet normalization is the root
cover

\[
 R[T]/(T^3+uT-v).
\tag{1.5f}
\]

It is free of rank three, and it is regular because eliminating \(v\)
identifies its completed source ring with \(k[[u,T,w]]\).

These are all homomorphisms (1.5d) to \(S_3\) which send both meridians to
transpositions.  Finite étale covers of the complement are classified by
the resulting finite permutation sets, and a normal finite extension over
\(R\) is the unique integral closure extending its complement cover.
Hence the completed strict-henselian normalization is one of (1.5e) or
(1.5f), and is finite free.  Completion and strict henselization are
faithfully flat, so flatness descends to \(\mathcal O_{Y,p}\).  QED

The ordinary-cusp hypothesis is essential.  At a more complicated
non-SNC point, the local complement can admit additional three-sheet
permutation representations, and neither the divisor ledger nor the
codimension-one different determines their integral closures.

### Proposition 1.5 -- point-flatness is cubic fiber-minimality

For every `p in Y`, let

\[
 \lambda(p)=
 \dim_{\kappa(p)}\!\left(B\otimes_A\kappa(p)\right),
\tag{1.6}
\]

the scheme-theoretic length of the finite normalization fiber.  Then

\[
 \pi\text{ is flat at }p
 \quad\Longleftrightarrow\quad
 \lambda(p)=3.
\tag{1.7}
\]

In particular,

\[
 |Z_{\mathrm{flat}}|
 =
 \{p\in Y:\lambda(p)\ge4\}.
\tag{1.8}
\]

Thus the point-flatness obligation can be replaced by the intrinsic
**cubic fiber-minimality** condition `lambda(p)<=3` at the finitely many
closed target points under the boundary.

#### Proof

Localize at `p` and write `R=A_p`, `N=B_p`, and `K=Frac(R)`.  The
torsion-free `R`-module `N` has generic rank three.  Hence every generating
set of `N` has at least three elements, while Nakayama's lemma identifies
the least number of generators with

\[
 \mu_R(N)=\dim_{\kappa(p)}N/\mathfrak m_pN=\lambda(p).
\]

If `lambda(p)=3`, choose a minimal surjection `R^3 -> N`.  After tensoring
with `K` it is a surjection between three-dimensional vector spaces and is
therefore an isomorphism.  Its kernel has generic rank zero, but it is a
submodule of the torsion-free module `R^3`; consequently the kernel is zero.
Thus `N` is free of rank three.  The converse is immediate.

Since `lambda(p)` is always at least three, failure of flatness is
equivalent to `lambda(p)>=4`.  Proposition 1.1 shows that these points are
closed and finite.  QED

Proposition 1.5 is strictly scheme-theoretic.  A bound on the number of
geometric places, or even on the number of reduced points in a fiber, does
not bound `lambda(p)`: multiplicity can create an excess-length fiber
without creating another place.  What suffices is a length-three statement
for the canonical finite fiber itself.  This is weaker than the
Cartier--Cohen--Macaulay atlas in Proposition 1.2 and is the smallest direct
closed-point input that removes the Fitting obstruction.

### Proposition 1.6 -- cubic algebra structure alone does not force flatness

There are integral normal degree-three covers of a smooth factorial
threefold whose trace-free module is reflexive but not locally free.
Consequently neither normality nor the existence of the commutative
rank-three algebra structure eliminates the point defect in Proposition
1.1.

#### Proof

On `X=A^3=Spec k[x,y,z]`, consider the module `M` with minimal presentation

\[
 0\longrightarrow A
 \mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}
 A^3\longrightarrow M\longrightarrow0.
\tag{1.9}
\]

It is the second Koszul syzygy of `(x,y,z)`.  Away from the origin one entry
of the relation is a unit, so `M` is free of rank two.  At the origin the
Auslander--Buchsbaum formula gives depth two; hence `M` is torsion-free and
`S_2`, and therefore reflexive.  Its determinant is trivial.  The relation
vanishes in the origin fiber, so

\[
 \dim_k(M/\mathfrak mM)=3.
\tag{1.10}
\]

To apply Tan's nonflat triple-cover existence theorem without suppressing
its genericity hypothesis, compactify this module.  On
`\mathbb P^3=Proj k[w,x,y,z]`, let `\mathcal M` be the cokernel of

\[
 0\longrightarrow\mathcal O_{\mathbb P^3}(-1)
 \mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}
 \mathcal O_{\mathbb P^3}^{\,3}\longrightarrow\mathcal M\longrightarrow0.
\]

The map vanishes only at `[1:0:0:0]`.  The same depth computation as above
shows that `\mathcal M` is rank-two reflexive, and its restriction to the
chart `w\ne0` is `M`.  Tan's Theorem 7.6 constructs a reduced irreducible
triple cover from a **general** member of the basepoint-free cubic system on
the projectivization of a rank-two reflexive sheaf.  Its projective
corollary says that, after tensoring by a sufficiently negative power of a
very ample line bundle, any rank-two reflexive sheaf on a projective
factorial variety occurs as the trace-free sheaf of such a cover.  Apply
that corollary to `\mathcal M`.  Restriction to `w\ne0` removes the harmless
line-bundle twist and gives an integral normal triple cover of `X` with

\[
 \pi_*\mathcal O_{\bar X}\simeq\mathcal O_X\oplus M.
\tag{1.11}
\]

Its origin fiber has length `1+3=4`, so Proposition 1.5 makes it nonflat
there.  Equivalently, the presentation of the right side of (1.11) has
relation `(0,z,-y,x)`, and its rank-three Fitting ideal is exactly
`(x,y,z)`.  QED

The word **general** is essential here.  Global generation supplies a
basepoint-free linear system; it does not make an arbitrarily prescribed
cubic tensor reduced, irreducible, or normal.  Nor may one normalize an
arbitrary cubic construction and assume without proof that normalization
preserves the displayed trace module.  The existence assertion above uses
Tan's generic reduced-irreducible divisor and his reflexive-extension
argument.  Any search with a prescribed initial cubic must reprove the
corresponding genericity statement inside that constrained subsystem.

This example is not asserted to contain a Keller open isomorphic to `A^3`
or to satisfy the minimal boundary ledger.  It has a narrower purpose: it
rules out any proof of point-flatness that uses only normal finite
degree-three algebra structure.  A successful proof must use the
distinguished Keller open, the scheme intersections of its boundary, or an
equivalent fiber-minimality input.

### Proposition 1.7 -- determinantal classification of every point defect

Let `p` be a point of `Z_flat`, put `R=A_p`, and write the trace splitting

\[
 B_p=R\oplus M.
\tag{1.12}
\]

There is a unique integer `s>=1` and a minimal free resolution

\[
 0\longrightarrow R^s
 \mathop{\longrightarrow}^{\Phi}
 R^{s+2}\longrightarrow M\longrightarrow0
\tag{1.13}
\]

such that

\[
 \lambda(p)=s+3,\qquad
 \operatorname{Fitt}_3^R(B_p)
 =\operatorname{Fitt}_2^R(M)
 =I_s(\Phi).
\tag{1.14}
\]

The maximal-minor ideal `I_s(Phi)` is primary to the maximal ideal of `R`.
In particular, the defect has its smallest possible fiber length four
exactly when

\[
 0\longrightarrow R
 \mathop{\longrightarrow}^{(a,b,c)^{\mathsf T}}
 R^3\longrightarrow M\longrightarrow0,
\tag{1.15}
\]

where `(a,b,c)` is a parameter ideal.  If the Fitting defect is reduced at
`p`, then `(a,b,c)` is the maximal ideal; after completion and a regular
change of parameters, (1.15) is the Koszul model `(x,y,z)^T`.

In fact, **every reduced point defect is automatically this minimal
length-four case**.  For `s>=2`, minimality puts the entries of `Phi` in the
maximal ideal and therefore

\[
 I_s(\Phi)\subseteq\mathfrak m^s\subseteq\mathfrak m^2.
\]

An ideal primary to `mathfrak m` and contained in `mathfrak m^2` cannot
equal its radical `mathfrak m`.  Thus all higher determinantal rungs, and
all non-Koszul `s=1` parameter ideals, are necessarily nonreduced defects.

#### Proof

The finite normal ring `B_p` is torsion-free over `R` and has depth at least
two: at every prime above `p`, normality supplies `S_2`, and depth for a
finite module is the minimum of the depths at those primes.  The trace
summand `M` is therefore a reflexive rank-two `R`-module.

The regular local ring `R` has dimension three and finite global dimension.
Auslander--Buchsbaum gives

\[
 \operatorname{pd}_R(M)+\operatorname{depth}_R(M)=3.
\]

If `M` had projective dimension zero it would be free and `p` would be
flat.  Thus `depth_R(M)=2` and `pd_R(M)=1`.  Its minimal resolution has the
form (1.13), because the difference of the two free ranks is
`rank_R(M)=2`.

Minimality puts every entry of `Phi` in the maximal ideal.  Hence

\[
 \dim_{\kappa(p)}M/\mathfrak mM=s+2
\]

and the unit summand gives `lambda(p)=s+3`.  The direct-sum rule for Fitting
ideals and the presentation (1.13) give (1.14).  Proposition 1.1 says its
vanishing locus is supported only at the closed point, so `I_s(Phi)` is
maximal-ideal-primary.

For `s=1`, this ideal is `(a,b,c)`.  Three generators of a height-three
ideal in the Cohen--Macaulay ring `R` form a parameter regular sequence.
If the ideal is reduced and primary to the maximal ideal, it equals that
maximal ideal.  Its generators are then a regular system of parameters,
giving the completed Koszul normal form.  The preceding containment proves
that a reduced defect cannot have `s>=2`.  QED

This reduces the closed-point part of the conjecture to a concrete
determinantal exclusion problem.  There is no hidden infinite homological
complex: the excess fiber length `lambda(p)-3` is exactly the number of
relations `s`.  The first unresolved local case is therefore the single
Koszul relation (1.15), with compatibility imposed by the distinguished
Keller open and the one-boundary scheme intersections.

### Proposition 1.8 -- every reduced defect has a square-zero fiber

Assume the Fitting defect at `p` is reduced.  Then it automatically has
length four, and its canonical fiber algebra is uniquely

\[
 B\otimes_A\kappa(p)
 \simeq
 \kappa(p)\oplus V,
 \qquad
 \dim_{\kappa(p)}V=3,
 \qquad
 V^2=0.
\tag{1.16}
\]

Consequently the fiber has one geometric point, that point lies in the
boundary, and

\[
 F^{-1}(p)=\varnothing.
\tag{1.17}
\]

#### Proof

By Proposition 1.7, reducedness forces `s=1`; after completion and a regular
change of parameters,

\[
 M=\operatorname{coker}
 \left(R\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}R^3\right).
\]

Write the product of two trace-free elements as a scalar component and an
`M` component.

First consider an `R`-linear functional `M -> R`.  It lifts to a row
`w=(w_1,w_2,w_3)` satisfying

\[
 w_1z-w_2y+w_3x=0.
\]

Reducing modulo the square of the maximal ideal, the linear independence of
`x,y,z` forces the constant row of `w` to vanish.  Thus every such
functional has image in the maximal ideal, and the scalar component of
`V*V` vanishes in the closed fiber.

Next let `T` be an endomorphism of `M`.  It lifts to a matrix
`C in Mat_3(R)` satisfying

\[
 C(z,-y,x)^{\mathsf T}
 =r(z,-y,x)^{\mathsf T}
\]

for some `r in R`.  Reduction modulo the square of the maximal ideal forces
the constant matrix of `C` to be the scalar matrix `r(0)I`.  Therefore, for
each `u in V`, the `M` component of multiplication by `u` has the form

\[
 v\longmapsto\ell(u)v
\]

for a linear functional `ell:V -> kappa(p)`.  Commutativity gives

\[
 \ell(u)v=\ell(v)u
\quad\text{for all }u,v\in V.
\]

Since `dim V=3`, this identity forces `ell=0`.  Both components of `V^2`
therefore vanish, proving (1.16).

The square-zero algebra in (1.16) is local and has a single geometric
point.  If that point lay in the Keller open, the fiber local algebra would
be reduced of length one because `F` is étale.  Hence it lies in the
boundary, and there is no other fiber point in `U`, proving (1.17).  QED

### Proposition 1.8a -- the first multiplication symbol is a ternary cubic

Retain the completed reduced-defect model

\[
 R=k[[x,y,z]],\qquad
 M=\operatorname{coker}\left(
 R\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}R^3
 \right).
\tag{1.16a}
\]

A nonflat triple-cover multiplication with trace module `M` is encoded by
the generalized Miranda--Tan tensor

\[
 \phi:\operatorname{Sym}^3(M)\longrightarrow\det(M)\simeq R.
\tag{1.16b}
\]

Choose the displayed generators `e_1,e_2,e_3` of `M` and write

\[
 c_{ijk}=\phi(e_i,e_j,e_k).
\]

Then every `c_ijk` lies in `m^3`.  The degree-three initial symbols form a
ten-dimensional vector space.  More precisely, put

\[
 r=(z,-y,x)\in R_1\otimes k^3.
\]

After choosing a volume form on `k^3`, every degree-three symbol is uniquely
of the form

\[
 \boxed{
 \operatorname{in}_3(\phi)(u,v,w)
 =
 h(r\times u,r\times v,r\times w)
 }
\tag{1.16c}
\]

for a ternary cubic `h`, where the right side means the symmetric
trilinear polarization of `h`.  Thus, up to the harmless determinant twist,
the first multiplication-symbol space is

\[
 \operatorname{Sym}^3(k^3)^\vee.
\tag{1.16d}
\]

#### Proof

The relation in (1.16a) makes well-definedness of `phi` equivalent to

\[
 zc_{1jk}-yc_{2jk}+xc_{3jk}=0
\tag{1.16e}
\]

for every symmetric pair `(j,k)`.  Taking homogeneous degree `d` symbols,
the solution space is the kernel of the natural contraction map

\[
 K_d=
 \ker\left[
 \operatorname{Sym}^d(V^\vee)\otimes\operatorname{Sym}^3(V^\vee)
 \longrightarrow
 \operatorname{Sym}^{d+1}(V^\vee)
 \otimes\operatorname{Sym}^2(V^\vee)
 \right],
\tag{1.16f}
\]

where `V=k^3`.  The Pieri decomposition shows that the displayed map has
zero kernel for `d=0,1,2`.  For `d>=3`, its kernel is the Schur module

\[
 K_d\simeq\mathbb S_{(d,3)}(V^\vee).
\tag{1.16g}
\]

In particular,

\[
 \dim K_3=\dim\mathbb S_{(3,3)}(V^\vee)=10.
\]

For every ternary cubic `h`, formula (1.16c) annihilates `r` in each
argument and therefore satisfies (1.16e).  Polarizing the ten ternary
monomials gives ten linearly independent elements of `K_3`.  Since `K_3`
has dimension ten, these exhaust it and prove (1.16c)--(1.16d).  QED

This removes a large false search space.  At a reduced defect, one should
not enumerate arbitrary bilinear multiplication tables and then impose
associativity: the generalized triple-cover correspondence already packages
associativity, and its first nonzero local datum is just a ternary cubic.
The initial counterexample search therefore begins with the classical
plane-cubic types

\[
 \text{smooth},\quad
 \text{nodal},\quad
 \text{cuspidal},\quad
 \text{reducible}
\tag{1.16h}
\]

of plane cubics, followed by higher-order lifts of (1.16c).

These are geometric types, not literally four `GL_3` orbits.  The smooth
type has its `j`-modulus, line-plus-conic splits into transverse and tangent
types, and three distinct lines split into a triangle and three concurrent
lines.  A useful working list of normal forms is

| type | normal form or family |
|---|---|
| smooth | `X^3+Y^3+Z^3-3 lambda XYZ`, `lambda^3 ne 1` |
| nodal | `Y^2Z-X^2(X+Z)` |
| cuspidal | `Y^2Z-X^3` |
| line plus transverse conic | `Z(XY-Z^2)` |
| line tangent to conic | `Z(YZ-X^2)` |
| three lines forming a triangle | `XYZ` |
| three concurrent lines | `XY(X-Y)` |
| double line plus line | `X^2Y` |
| triple line | `X^3` |

### Proposition 1.8b -- the cubic is the exceptional incidence curve

Blow up the defect point in `Spec R` and let
`E\simeq\mathbb P^2` be the exceptional plane.  After removing pullback
torsion, the transform of `M` on `E` is the universal quotient

\[
 0\longrightarrow\mathcal O_E(-1)
 \longrightarrow\mathcal O_E^{\,3}
 \longrightarrow Q\longrightarrow0.
\tag{1.16i}
\]

For `[r]\in E`, cross product identifies

\[
 Q_{[r]}=k^3/\langle r\rangle
 \mathop{\longrightarrow}^{\sim}r^\perp,
 \qquad [u]\longmapsto r\times u.
\tag{1.16j}
\]

Consequently the exceptional cubic divisor supplied by the order-three
symbol `h` is

\[
 Z_h=
 \left\{
 ([r],[q])\in\mathbb P^2\times\mathbb P^2:
 r\cdot q=0,\ h(q)=0
 \right\}.
\tag{1.16k}
\]

It is the pullback of the plane cubic `C_h=V(h)` along the second projection
of the incidence flag variety.  In particular, scheme-theoretically,
`Z_h\to C_h` is a `\mathbb P^1`-bundle.  Thus reducedness, irreducibility,
and the number of components of this first exceptional divisor are exactly
those of `C_h`.

#### Proof

On the blowup chart the Koszul relation has a common exceptional parameter.
Saturating removes that parameter and leaves the tautological line
`\langle r\rangle\subset k^3`, giving (1.16i).  Cross product has kernel
`\langle r\rangle` and image `r^\perp`, proving (1.16j).  On the diagonal
`u=v=w`, formula (1.16c) becomes `h(r\times u)`.  Under (1.16j), its zero
scheme is therefore (1.16k).  The incidence projection has fiber the line
`r\cdot q=0` over every `[q]`, proving the last assertion.  QED

This gives a first local ranking, but not a classification.  Smooth `h`
produces a ruled exceptional surface over an elliptic curve.  Nodal and
cuspidal `h` are irreducible with rational normalization.  Reducible `h`
produces several exceptional components, and a double or triple line gives
a nonreduced initial exceptional divisor.  None of these facts alone decides
normality of the original finite threefold: an exceptional divisor on a
modification may be nonnormal even when the unmodified cover is normal.

Higher lifts have the form

\[
 \phi=\phi_h+\phi_4+\phi_5+\cdots,\qquad
 \phi_d\in K_d\simeq\mathbb S_{(d,3)}(k^3)^\vee.
\tag{1.16l}
\]

Tan's unconstrained generic-section theorem does not by itself say what
happens after `h` is fixed.  In the present punctured affine model, however,
the constrained genericity statement follows from the same method.

### Proposition 1.8c -- every nonzero cubic type has normal integral lifts

For every nonzero ternary cubic `h`, there is a generalized cubic tensor
`\phi` with initial symbol `h` whose reflexive triple-cover extension is
integral and normal and has trace module `M`.

#### Proof

Work first with the algebraic Koszul module `M_A` over
`A=k[x,y,z]`; completion at the origin recovers (1.16a).  Put

\[
 T_A=
 \operatorname{Hom}_A(\operatorname{Sym}^3M_A,\det M_A)
\]

and let `m=(x,y,z)`.  Formula (1.16c) supplies a polynomial tensor
`\phi_h` with initial symbol `h`.  Every perturbation in `mT_A`
has order at least four and therefore preserves that initial symbol.

On the punctured affine threefold `U=Spec A\setminus\{m\}`, one has
`m\mathcal O_U=\mathcal O_U`.  Hence

\[
 \widetilde{mT_A}|_U=\widetilde{T_A}|_U.
\tag{1.16m}
\]

The module `M_A|_U` is locally free, and `\widetilde{T_A}|_U` is the rank-four
bundle of ordinary binary cubics.  Quasi-compactness permits a
finite-dimensional subspace `W\subset mT_A` which generates this
bundle on `U`.  Enlarging `W` by coordinate multiples makes its induced
linear system on `\mathbb P(M_A|_U)` separate the base directions as well as
restrict to the complete cubic system on every projective-line fiber.

For general `\psi\in W`, the section `\phi_h+\psi` is nowhere the zero
binary cubic on a fiber: vanishing imposes four independent linear
conditions over a three-dimensional base.  Its divisor in
`\mathbb P(M_A|_U)` is therefore finite flat of degree three over `U`.
Characteristic-zero Bertini, applied to this basepoint-free system which is
not composed with a pencil, makes the general divisor smooth and
irreducible.

Tan's codimension-two extension construction now extends this punctured
triple cover to the normalization over `Spec A`.  The extension is integral
and normal.  Its trace-free module and `M_A` agree on `U`; both are
reflexive, so they agree on `Spec A`.  Finally, the perturbation lies in
`mT_A`, so its order-three symbol remains `h`.  Localizing and completing
at the origin gives the asserted cover of the completed model; normality is
preserved because these finite-type local rings are excellent.  QED

Thus initial cubic type alone cannot be excluded by abstract existence of
a normal integral cover.  What remains nontrivial is compatibility of such
a lift with a single global Keller open and its marked boundary.  A
bounded-degree computational ansatz must also verify that it contains
enough of the perturbation space used above; genericity in the full module
does not imply genericity in an arbitrarily small ansatz.

### Proposition 1.8d -- a reduced defect forces branch multiplicity six

Let \(p\) be a reduced point defect and let \(h\), possibly zero, be the
degree-three ternary-cubic symbol of Proposition 1.8a.  On the exceptional plane
\(E=\mathbb P(T_pY)\), define

\[
 \mathcal H_h([r])
 =
 \operatorname{Disc}\!\left(h|_{r^\perp}\right).
\tag{1.16n}
\]

This is a well-defined section of \(\mathcal O_E(6)\).  It is the
degree-six initial form of the branch discriminant at \(p\).  Consequently:

\[
 \boxed{
 \begin{array}{ll}
 h\ \text{squarefree}
   &\Longrightarrow\
     \operatorname{mult}_p(D_{\mathrm{br}})=6,\\[1mm]
 h\ \text{non-squarefree or zero}
   &\Longrightarrow\
     \operatorname{mult}_p(D_{\mathrm{br}})\ge7.
 \end{array}}
\tag{1.16o}
\]

In the squarefree case the projectivized tangent cone is the discriminant
of line sections of \(C_h=V(h)\).  For a smooth cubic it is the dual sextic.
For a singular reduced cubic it is the dual curve together with the
pencils through its singular points, counted with their discriminant
multiplicities.

#### Proof

Proposition 1.8b identifies the exceptional cubic over a direction
\([r]\) with

\[
 C_h\cap r^\perp.
\]

The transformed degree-three cover is ramified at \([r]\) exactly when this
length-three line section is nonreduced.  Its binary-cubic discriminant is
therefore (1.16n).

Choose a chart \(r_3\ne0\) and write the line as
\[
 q_3=-\frac{r_1}{r_3}q_1-\frac{r_2}{r_3}q_2.
\]
The discriminant of the resulting binary cubic, multiplied by \(r_3^6\),
is homogeneous of degree six in \(r_1,r_2,r_3\).  The transition under a
different basis of \(r^\perp\) is the sixth power of the determinant, so
these chart expressions glue to a section of \(\mathcal O_E(6)\).

The section is nonzero exactly when a general line meets \(C_h\) in three
distinct points, equivalently when \(h\) is squarefree.  Since no
multiplication symbol exists in orders zero, one, or two, higher-order
terms of the generalized cubic tensor cannot contribute below this
degree-six discriminant.  Proposition 1.4 gives simple generic
ramification, so the algebra discriminant and the reduced branch equation
have the same codimension-one divisor.  Thus a nonzero
\(\mathcal H_h\) is the initial branch equation and gives multiplicity six.
If it vanishes identically, the first possible branch term has strictly
larger integral order.  This proves (1.16o).  QED

For the standard squarefree degenerations, the degree-six tangent cone
factors as follows after choosing dual coordinates \(A,B,C\):

| \(C_h\) | \(\mathcal H_h\), up to a nonzero scalar |
|---|---|
| smooth cubic | irreducible dual sextic |
| nodal cubic | \(C^2\) times the dual quartic |
| cuspidal cubic | \(C^3\) times the dual cubic |
| line plus transverse conic | \(A^2B^2(4AB-C^2)\) |
| line tangent to conic | \(B^4(A^2-4BC)\) |
| triangle | \(A^2B^2C^2\) |
| three concurrent lines | \(C^6\) |

Here the factor \(C=0\) is the pencil of lines through the displayed
singular point.  A double or triple component makes
\(\mathcal H_h=0\), in agreement with the second row of (1.16o).  The same
row includes a zero degree-three symbol, when the multiplication tensor
begins in still higher order.

Proposition 1.4b now has a numerical converse for reduced defects: not only
must their branch point be worse than an ordinary cusp, whose multiplicity
is two, but its branch multiplicity is at least six.  This still does not
exclude the defect from a boundary-minimal Keller map, because the present
boundary invariant records divisorial components and different orders, not
closed-point branch multiplicity.

### Exact computation 1.8e -- symbol-stratified double saturation

The order-three symbol can now be tested against the two canonical
saturation modules, provided the chosen full multiplication tensor is kept
explicit.  Put

\[
 A=\mathbb Q[x,y,z],\qquad
 M=\operatorname{coker}\!\left(
 A\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}A^3
 \right),\qquad B=A\oplus M.
\tag{1.16p}
\]

For a ternary cubic \(h\), let \(\phi_h\) be the homogeneous generalized
triple-cover tensor from (1.16c).  The multiplication is recovered without
choosing a punctured chart.  If \(r=(z,-y,x)\), its trace-free part
\(\mu:\operatorname{Sym}^2M\to M\) is characterized by

\[
 \det(r,\mu(u,v),w)=3\phi_h(u,v,w),
\tag{1.16q}
\]

and the scalar part is

\[
 s(u,v)=\operatorname{tr}_M(\mu_u\mu_v).
\tag{1.16r}
\]

These formulas give a twelve-generator \(A\)-presentation of
\(Q=\Omega_{B/A}\).  The kernel of

\[
 B\longrightarrow Q^3,\qquad
 b\longmapsto(b\,de_1,b\,de_2,b\,de_3)
\tag{1.16s}
\]

is \(\operatorname{Ann}_B(Q)\), so its first four syzygy coordinates give
an \(A\)-presentation of
\(T=B/\operatorname{Ann}_B(Q)\).  Exact Singular calculations give:

| symbol stratum | \(\dim T\) | \(H^0_{(x,y,z)}(Q)\) | \(\operatorname{Ext}^2_A(T,A)\) |
|---|---:|---:|---:|
| smooth | \(2\) | \(0\) | length \(6\) |
| nodal | \(2\) | \(0\) | length \(6\) |
| cuspidal | \(2\) | \(0\) | length \(6\) |
| line plus transverse conic | \(2\) | \(0\) | length \(6\) |
| line tangent to conic | \(2\) | \(0\) | length \(6\) |
| triangle | \(2\) | \(0\) | length \(6\) |
| three concurrent lines | \(2\) | \(0\) | length \(6\) |
| double line | \(3\) | \(0\) | dimension \(1\) |
| triple line | \(3\) | \(0\) | dimension \(1\) |
| zero tensor | \(3\) | \(0\) | \(0\) |

Thus every squarefree homogeneous model passes conormal saturation but
fails support saturation by a length-six canonical dual of \(C/T\).  The
length-six module has three-dimensional top, is annihilated by
\(\mathfrak m^2\), and therefore has Hilbert function
\(3+3t\).  The repeated coefficient six is consequently a uniform
two-layer invariant, not only an equality of total lengths.  The
double- and triple-line homogeneous models have the wrong support
dimension, so in particular they fail the purity hypothesis: their support
obstruction is not finite length.  The zero homogeneous
tensor passes both displayed module tests only because its square-zero
algebra is nowhere generically étale; it is not a degree-three cover of the
required generic kind.

This last behavior is a property of the lift, not of the cubic symbol
alone.  The order-four constraint space has dimension \(24\).  For one
explicit integral linear combination of its exact kernel basis, adding the
same order-four tensor to each of the nine nonzero orbit representatives
and to the zero symbol
gives

\[
 \dim T=2,\qquad H^0_{(x,y,z)}(Q)=0,\qquad
 \dim_{\mathbb Q}\operatorname{Ext}^2_A(T,A)=6
\tag{1.16t}
\]

in every row, including the double line, triple line, and zero symbol.
Consequently no generically valid ternary-cubic stratum passes double
saturation in these two exact leading models: the surviving obstruction is
always support saturation once purity and generic étaleness are restored.

There is also an exact deformation statement for the seven squarefree
rows.  Let \(\psi _4\) be the order-four tensor used above and introduce a
parameter \(t\).  Over
\(\mathbb Q[t,x,y,z]\), the family

\[
 \phi_h+t\psi _4
\tag{1.16ta}
\]

has uniformly saturated cotangent presentation.  Its relative
\(\operatorname{Ext}^2\) module is supported scheme-theoretically on the
\(t\)-axis, has no \(t\)-torsion, and has multiplicity six.  More strongly,
the computed presentation is equal to the scalar extension of its
specialization at \(t=0\).  Thus the length-six defect is flat and constant
along each of these seven chosen deformation lines; it is not merely an
agreement between their two endpoints.

The line test can be enlarged without introducing a multivariate
coefficient ring.  Choose the exact nullspace basis
\(\psi_{4,0},\ldots,\psi_{4,23}\) of the order-four constraint space.  For
every squarefree \(h\) and every basis index \(i\), the family

\[
 \phi_h+t\psi_{4,i}
\tag{1.16tb}
\]

has uniformly saturated cotangent presentation, relative
\(\operatorname{Ext}^2\) multiplicity six, no \(t\)-torsion, and radical
support equal to the collision axis.  Thus all \(7\cdot24=168\) basis-axis
families retain the flat length-six support defect.

Literal equality of the Gröbner presentation with its central
specialization is not invariant.  It changes in exactly four basis-axis
pairs: concurrent lines in directions \(13,17\), line tangent to a conic
in direction \(10\), and the triangle in direction \(17\).  None of these
four changes affects flatness, support, or multiplicity.  Accordingly,
presentation equality is useful as a strong coordinate check on the
original seven lines, but it is not part of the geometric
double-saturation criterion.

The following coordinate-subspace calculations predate the universal
formal-gauge proof.  They remain useful independent regressions, but they
are not the active smooth-symbol frontier.  For every pair
\(0\leq i<j<24\), compute the full coordinate
plane over \(\mathbb Q[u,v,x,y,z]\):

\[
 \phi_h+u\psi_{4,i}+v\psi_{4,j}.
\tag{1.16tc}
\]

On all \(\binom{24}{2}=276\) coordinate planes, the cotangent presentation
is saturated and the relative \(\operatorname{Ext}^2\) presentation is
equal to the scalar extension of its value at \(u=v=0\).  Its support is
the parameter plane and its relative multiplicity is six.  This includes
every specialization on every two-basis plane, and excludes every
order-four escape supported on at most two basis tensors for the smooth
cubic.

The same polynomial-base test remains tractable on coordinate
three-spaces.  For every \(0\leq i<j<k<24\), the family

\[
 \phi_h+p_0\psi_{4,i}+p_1\psi_{4,j}+p_2\psi_{4,k}
\tag{1.16td}
\]

over \(\mathbb Q[p_0,p_1,p_2,x,y,z]\) has saturated cotangent presentation,
support equal to the parameter three-space, and relative multiplicity six.
After pruning contractible free summands, its rank-three
\(\operatorname{Ext}^2\) presentation is pulled back from the parameter
origin in all \(\binom{24}{3}=2024\) rows.  Two raw resolutions change
ambient presentation, but pruning identifies both changes as nonminimal
artifacts.  Thus every smooth order-four direction supported on at most
three basis tensors retains the length-six defect.

The coordinate-plane calculation also closes for all six singular
squarefree symbols.  Among their
\(6\binom{24}{2}=1656\) planes, 1652 have the same pruned rank-three
presentation as the parameter origin.  Four planes have a nonconstant
ambient presentation:

\[
\begin{split}
 &(h_{\mathrm{cusp}},\{3,8\}),\quad
 (h_{\mathrm{cusp}},\{4,9\}),\\
 &(h_{\mathrm{line\text{-}tangent\text{-}conic}},\{2,9\}),\quad
 (h_{\mathrm{node}},\{5,11\}).
\end{split}
\tag{1.16te}
\]

This presentation jump is not a flatness defect.  In each exceptional
row the relative module is annihilated by
\((x,y,z)^2\), so it has an exact finite presentation over
\(S=\mathbb Q[p_0,p_1]\).  Direct calculation gives

\[
 \operatorname{Fitt}_6^S(M)=S,\qquad
 \operatorname{Fitt}_5^S(M)=0.
\tag{1.16tf}
\]

Hence \(M\) is locally free of rank six over the full parameter plane,
and Quillen--Suslin makes it free.  Together with (1.16tc), all
\(7\binom{24}{2}=1932\) coordinate planes for squarefree cubic symbols
retain the flat length-six support defect.

There is now also a plane whose generic tensor has full support in the
fixed 24-element kernel basis.  Put

\[
 \psi_+=\sum_{i=0}^{23}\psi_{4,i},\qquad
 \psi_-=\sum_{i=0}^{23}(-1)^i\psi_{4,i}.
\]

For each of the seven squarefree symbols, compute the complete family

\[
 \phi_h+u\psi_++v\psi_-
\tag{1.16tg}
\]

over \(\mathbb Q[u,v,x,y,z]\).  In all seven rows the cotangent
presentation is saturated, the relative
\(\operatorname{Ext}^2\) support is exactly the parameter plane, its
relative multiplicity is six, and its pruned rank-three presentation is
the scalar extension of the presentation at \(u=v=0\).  On
\(D(u^2-v^2)\), its coefficient in every one of the 24 fixed basis
directions is nonzero.  Thus (1.16tg) is an exact two-parameter family of
full-support quartic tensors retaining the length-six defect for every
squarefree cubic symbol.

This removes basis sparsity as a necessary escape mechanism: a quartic
tensor can use all 24 directions and still retain the same obstruction.
It does not give a Zariski-open statement in the full 24-dimensional
kernel, because (1.16tg) is only one plane.

The purity-restoring lift for the three degenerate symbols is stable on a
translated version of the same plane.  Let \(\psi_g\) be the deterministic
full-support quartic tensor used in (1.16t).  For all nine nonzero cubic
orbit representatives and the zero symbol, compute

\[
 \phi_h+\psi_g+u\psi_++v\psi_- .
\tag{1.16th}
\]

Over \(\mathbb Q[u,v,x,y,z]\), every one of these ten affine planes has
saturated cotangent presentation, support exactly equal to the parameter
plane, relative \(\operatorname{Ext}^2\) multiplicity six, and pruned
rank-three presentation pulled back from \(u=v=0\).  The coordinates of
\(\psi_g\) in the fixed primitive basis are all nonzero, so the generic
member of (1.16th) also has full support.

For the double-line, triple-line, and zero cubic symbols, this is stronger
than the central endpoint calculation (1.16t): the order-four translation
restores pure two-dimensional ramification support and the finite
length-six defect uniformly on an affine parameter plane.  It still does
not prove normality of these nonhomogeneous algebras or compatibility with
a distinguished Keller open.

The plane computations do **not** by themselves prove lift-independence.
The later formal-gauge theorem in
[`UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md`](UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md)
covers every linear combination in the full smooth-symbol 24-parameter
order-four space and proves its cotangent saturation.  Consequently the
smooth row is closed as a quartic saturation search.  Its remaining problem
is global: algebraize the formal equivalence and test the marked boundary
and Keller-open compatibility.

For the six singular squarefree rows, the gauge-cokernel atlas gives the
exact quartic nongauge dimensions and annihilators

\[
\begin{array}{c|c|c}
h&\dim (Q_h)_4&\operatorname{Ann}(Q_h)\\ \hline
\text{nodal}&2&(x)\\
\text{cuspidal}&4&(x^2)\\
\text{line + transverse conic}&4&(yz)\\
\text{line + tangent conic}&6&(y^3)\\
\text{triangle}&6&(xyz)\\
\text{concurrent lines}&8&(x^3).
\end{array}
\]

The routing certificate generated by
[`compile_support_saturation_cases.py`](../scripts/compile_support_saturation_cases.py)
as `support_saturation_cubic_annihilator_frontier.json` first isolated this
queue.  The double-saturation computation below now resolves the complete
chosen complement in every row.  It does not classify higher-order formal
deformations: removing a gauge component at quartic order can create terms
of order five and above.

### Exact computation 1.8f -- singular-squarefree nongauge double saturation

Let \(K_4\) be the 24-dimensional space of compatible quartic tensors in
the fixed primitive basis \(\psi_0,\ldots,\psi_{23}\).  For each singular
squarefree cubic symbol \(h\), choose the following deterministic complement
of the degree-four image of the determinant-twisted formal gauge map:

\[
\begin{array}{c|c|c}
h&\dim K_4/\operatorname{im}(G_h)_4&\text{basis indices}\\ \hline
\text{nodal}&2&0,1\\
\text{cuspidal}&4&0,1,2,5\\
\text{line + transverse conic}&4&15,19,20,23\\
\text{line + tangent conic}&6&8,12,15,16,20,21\\
\text{triangle}&6&0,1,15,19,20,23\\
\text{concurrent lines}&8&0,1,2,3,4,5,6,7.
\end{array}
\tag{1.16ti}
\]

Write \(S_h=\mathbb Q[p_0,\ldots,p_{r-1}]\),
\(A_h=S_h[x,y,z]\), and form the complete displayed complement family

\[
 \phi_h+\sum_{i=0}^{r-1}p_i\psi_{j_i}.
\tag{1.16tj}
\]

Let \(Q_h=\Omega_{B_h/A_h}\), let
\(T_h=B_h/\operatorname{Ann}_{B_h}(Q_h)\), and put
\(\mathfrak m=(x,y,z)\).  Exact characteristic-zero module calculation
gives, in every one of the six rows,

\[
 H^0_{\mathfrak m}(Q_h)=0.
\tag{1.16tk}
\]

On the other hand, \(\operatorname{Ext}^2_{A_h}(T_h,A_h)\) has radical
support exactly \(V(\mathfrak m)\), multiplicity six, and a pruned
rank-three presentation equal literally to the scalar extension of its
central presentation at \(p_0=\cdots=p_{r-1}=0\).  Moreover

\[
 \mathfrak m^2\operatorname{Ext}^2_{A_h}(T_h,A_h)=0,
 \qquad
 \operatorname{Ext}^3_{A_h}(T_h,A_h)=0.
\tag{1.16tl}
\]

These statements commute with every geometric specialization of $S_h$.
More precisely, give the twelve cotangent generators $de_i,e_jde_i$
filtered degrees two and four.  Let $\mathcal R(Q_h)$ be the weighted
Rees presentation, and let $W_h$ be the cokernel of the action map

\[
 B_h\longrightarrow Q_h^{\oplus3},\qquad
 b\longmapsto (bde_1,bde_2,bde_3).
\tag{1.16tm}
\]

The same component filtration defines $\mathcal R(W_h)$.  Exact module
saturation gives

\[
\begin{aligned}
 (N_{\mathcal R(Q_h)}:t^\infty)&=N_{\mathcal R(Q_h)},&
 \mathcal R(Q_h)/(t)&\simeq S_h\otimes_{\mathbb Q}Q_{h,0},\\
 (N_{\mathcal R(W_h)}:t^\infty)&=N_{\mathcal R(W_h)},&
 \mathcal R(W_h)/(t)&\simeq S_h\otimes_{\mathbb Q}W_{h,0}.
\end{aligned}
\tag{1.16tn}
\]

Thus both filtered modules are strict and $S_h$-flat.  Flatness of
$Q_h$ and $W_h$ makes the kernel of (1.16tm) commute with arbitrary base
change, so the specialization of $T_h$ is the intrinsic module

\[
 T_{h,s}=B_{h,s}/\operatorname{Ann}_{B_{h,s}}(Q_{h,s})
\tag{1.16to}
\]

for every geometric point $s\to\operatorname{Spec}S_h$, rather than only
the specialization of a chosen relative presentation.  It is also
$S_h$-flat.  The constant relative `Ext^2` presentation is $S_h$-flat,
and relative `Ext^3` vanishes, so cohomology and base change applied to a
finite free resolution of $T_h$ gives

\[
 \operatorname{Ext}^2_{A_h}(T_h,A_h)\otimes_{S_h}k(s)
 \simeq
 \operatorname{Ext}^2_{k(s)[x,y,z]}
     (T_{h,s},k(s)[x,y,z]).
\tag{1.16tp}
\]

Finally, a closed-point-torsion class in $Q_{h,s}$ would have a
closed-point-torsion initial class in its associated graded module, but
(1.16tn) identifies that module with the saturated central cotangent
module.  Therefore

\[
 H^0_{(x,y,z)}(Q_{h,s})=0,
 \qquad
 \operatorname{length}_{k(s)}
 \operatorname{Ext}^2(T_{h,s},k(s)[x,y,z])=6
\tag{1.16tq}
\]

for every geometric $s$.  Hence `C2` passes and `C1` fails on every fiber
of each complete displayed quartic nongauge family, with the same
square-zero multiplicity-six support-hull defect.  Explicitly, Proposition
1.15 and finite-length duality give, with $(-)^\vee$ denoting the canonical
finite-length dual,

\[
 L_{h,s}^{\vee}\simeq
 \operatorname{Ext}^2(T_{h,s},k(s)[x,y,z]),
 \qquad
 \operatorname{length}L_{h,s}=6,
 \qquad
 (x,y,z)^2L_{h,s}=0.
\tag{1.16tr}
\]

There is a second, independent local consequence.  Let

\[
 J_{h,s}=\operatorname{Ann}_{B_{h,s}}(Q_{h,s}),
 \qquad
 \mathfrak n_s=(x,y,z,e_1,e_2,e_3)\subset B_{h,s}
\]

be the Kähler different and the collision maximal ideal.  Reducing the
same exact annihilator presentation modulo its products by
`x,y,z,e_1,e_2,e_3` gives a relative Nakayama module

\[
 \mathcal G_h=J_h/\mathfrak n_hJ_h\simeq S_h^{\oplus6}.
\tag{1.16tr1}
\]

More precisely, its radical support is the parameter axis, its multiplicity
and pruned presentation rank are six, and its pruned presentation is
literally the scalar extension of the central presentation.  The strict
Rees base-change certificate for the annihilator identifies every geometric
fiber with the intrinsic different.  Consequently

\[
 \boxed{\dim_{k(s)}J_{h,s}/\mathfrak n_sJ_{h,s}=6}
\tag{1.16tr2}
\]

for every geometric parameter value in all six rows.  Nakayama's lemma says
that this dimension is the minimal number of generators of the localized
different.  A nonzero principal ideal has one minimal generator.  Hence the
Kähler different is not locally principal, and in particular is not
Cartier, at any collision in these six quartic families.

This is an exact incompatibility between the persistent length-six quartic
models and the Cartier-different hypothesis.  Its scope stops at the full
quartic nongauge layer: order-five and higher formal corrections could
change the annihilator presentation, and no Keller-open or normality claim
is made here.

The certificate is
[`cubic_double_saturation_stratification.json`](../artifacts/generated-results/cubic_double_saturation_stratification.json),
replayed by
[`verify_cubic_double_saturation_stratification.py`](../scripts/verify_cubic_double_saturation_stratification.py).
It records all representative tensors, complement ranks, exact presentation
tests, and hashes.  The conclusion, including the non-Cartier calculation,
concerns the first quartic nongauge layer;
it does not prove that these algebras contain a Keller open or that every
higher formal deformation is equivalent to one of (1.16tj).

<!-- status-consumer: KDSQ6 cd423f625f1f3cd2 -->

### Exact computation 1.8f.1 -- nodal sextic persistence of the non-Cartier different

Let `eta` generate the nodal gauge cokernel

\[
 \ker(C)/\operatorname{im}(G_{\mathrm{nod}})
 \simeq \mathbb Q[y,z](-3).
\]

The complete homogeneous compatible spaces in collision degrees four,
five, and six have the exact decompositions

\[
 24=22+2,
 \qquad
 42=39+3,
 \qquad
 64=60+4,
\tag{1.16tr3}
\]

where the first summands are the formal-gauge images and complements are

\[
 \langle y\eta,z\eta\rangle,
 \qquad
 \langle y^2\eta,yz\eta,z^2\eta\rangle,
 \qquad
 \langle y^3\eta,y^2z\eta,yz^2\eta,z^3\eta\rangle.
\tag{1.16tr4}
\]

Form the exact nine-parameter polynomial tensor family containing these two
quartic, three quintic, and four sextic directions.  Over
`S=Q[p_0,...,p_8]`, the cotangent presentation is saturated and the
support-hull obstruction remains the scalar extension of the central
length-six module.  The collision Nakayama quotient of the intrinsic
Kähler different satisfies

\[
 J/\mathfrak nJ\simeq S^{\oplus6}.
\tag{1.16tr5}
\]

Weighted Rees presentations for `Omega` and for
`coker(B -> Omega^3)` are `t`-saturated and have their literal central
initial presentations.  Thus the annihilator commutes with every geometric
specialization, and Proposition 1.15c gives

\[
 \boxed{
 \dim_{k(s)}J_s/\mathfrak n_sJ_s=6
 }
\tag{1.16tr6}
\]

on every fiber.  Hence the Kähler different is non-Cartier throughout the
complete nodal order-six normal-form family.

This finite-jet certificate closes the nodal generator-persistence problem
through the complete sextic gauge quotient.  By itself it does not control
the order-seven terms created when sextic gauge directions are removed; the
all-orders theorem immediately below supplies that missing argument.
Normality and Keller-open compatibility remain open.  The pinned certificate
is
[`nodal_sextic_different_persistence.json`](../artifacts/generated-results/nodal_sextic_different_persistence.json),
replayed by
[`verify_nodal_sextic_different_persistence.py`](../scripts/verify_nodal_sextic_different_persistence.py).

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

### Theorem 1.8f.2 -- all-orders nodal non-Cartier persistence

Let \(k\) be a characteristic-zero field, let
\(\widehat A=k[[x,y,z]]\), and fix the nodal cubic tensor
\(\phi_{\mathrm{nod}}\) attached to

\[
 h_{\mathrm{nod}}=Y^2Z-X^2(X+Z).
\]

Let \(\Phi\) be any compatible formal tensor with
\(\Phi-\phi_{\mathrm{nod}}\) of collision order at least four.  Then:

1. after a determinant-twisted formal gauge transformation, \(\Phi\) has
   the form
   \[
    \phi_f=\phi_{\mathrm{nod}}+f(y,z)\eta,
    \qquad f\in(y,z)k[[y,z]],
   \tag{1.16tr7}
   \]
   where \(\eta\) is the compatible tensor attached to \(Z^3\);
2. for the resulting completed cubic algebra \(B_f\),
   \[
    H^0_{(x,y,z)}(\Omega_{B_f/\widehat A})=0;
   \tag{1.16tr8}
   \]
3. if
   \(J_f=\operatorname{Ann}_{B_f}(\Omega_{B_f/\widehat A})\) and
   \(\mathfrak n_f=(x,y,z,e_1,e_2,e_3)\), then
   \[
    \boxed{\dim_k J_f/\mathfrak n_fJ_f=6.}
   \tag{1.16tr9}
   \]

Consequently the intrinsic Kähler different is non-Cartier for every
compatible formal correction with fixed nodal leading symbol.  For an
algebraic or localized tensor, the same conclusion holds at the collision:
completion is faithfully flat, so a cyclic local different could not acquire
six minimal generators after completion.

#### Proof

The exact graded calculation (5.15)--(5.17) gives

\[
 \ker C=\operatorname{im}G_{\mathrm{nod}}+A\eta,
 \qquad
 \ker C/\operatorname{im}G_{\mathrm{nod}}
 \simeq A/(x)(-3)=k[y,z](-3).
\tag{1.16tr10}
\]

Suppose terms below collision degree \(d\) have already been put in normal
form.  Decompose the compatible degree-\(d\) remainder, using (1.16tr10), as

\[
 G_{\mathrm{nod}}(D_{d-3})+f_{d-3}(y,z)\eta.
\]

The gauge transformation with first nonidentity term \(D_{d-3}\) removes
the first summand.  Its action on terms of degree at least four starts in
degree \(d+1\), so it does not alter any earlier normal-form term.  Iterating
in the collision-adic topology gives (1.16tr7).  This is an existence
statement; stabilizers can change the representative \(f\) in higher order.

It remains to show that no tail \(f\) can make the different cyclic.  Introduce
one universal coefficient \(u\) and form the algebra \(B_u\) over

\[
 R=k[u,x,y,z],
 \qquad \phi_u=\phi_{\mathrm{nod}}+u\eta.
\tag{1.16tr11}
\]

The exact multiplication reconstruction is affine-linear in \(u\).  In the
trace-free products the only changes are

\[
 \delta(e_1^2)=(0,0,-3uy^2),\quad
 \delta(e_1e_2)=(0,0,-3uyz),\quad
 \delta(e_2^2)=(0,0,-3uz^2),
\tag{1.16tr12}
\]

and the other three trace-free products do not change.  The six scalar
corrections are also linear in \(u\); explicitly they are

\[
\begin{array}{c|c}
11&6u(y^4-y^2z^2)\\
12&u(-9xy^3+6y^3z-6yz^3)\\
13&-9uy^4\\
22&u(-18xy^2z+6y^2z^2-6z^4)\\
23&-9uy^3z\\
33&0.
\end{array}
\tag{1.16tr13}
\]

In particular there are no hidden \(u^2\)-terms.  Exact annihilator
reduction over \(R\) gives

\[
 J_u/\mathfrak n_{\mathrm{rel}}J_u\simeq k[u]^{\oplus6},
 \qquad
 \mathfrak n_{\mathrm{rel}}=(x,y,z,e_1,e_2,e_3).
\tag{1.16tr14}
\]

More concretely, the pruned presentation is independent of \(u\) and is the
direct sum of six copies of the row \((z,y,x)\).  The same calculation gives
cotangent saturation over \(R\).

Give \(u,x,y,z\) collision weight one, \(de_i\) weight two, and \(e_jde_i\)
weight four.  Put

\[
 Q_u=\Omega_{B_u/R},\qquad
 W_u=\operatorname{coker}
 \left(B_u\longrightarrow Q_u^{\oplus3}\right).
\]

The exact weighted-Rees calculation gives

\[
\begin{aligned}
 (N_{\mathcal R(Q_u)}:t^\infty)&=N_{\mathcal R(Q_u)},&
 \operatorname{gr}Q_u&\simeq Q_0\otimes_k k[u],\\
 (N_{\mathcal R(W_u)}:t^\infty)&=N_{\mathcal R(W_u)},&
 \operatorname{gr}W_u&\simeq W_0\otimes_k k[u].
\end{aligned}
\tag{1.16tr15}
\]

Now let \(f\in(x,y,z)k[[x,y,z]]\); the normal form only needs the smaller
ideal \((y,z)\).  In the graph specialization \(u\mapsto f\), the initial
form of \(u-f\) is either \(u\), if \(\operatorname{ord}(f)>1\), or
\(u-f_1(x,y,z)\), if \(f_1\) is the linear part.  It is monic in \(u\).
Multiplication by a monic polynomial in \(u\) is injective on \(M[u]\) for
every module \(M\).  Hence (1.16tr15) and the filtered-regularity lemma show
that \(u-f\) is regular on both \(Q_u\) and \(W_u\), after completion as
well.

Apply this to the exact sequence

\[
 0\longrightarrow J_u\longrightarrow B_u
 \longrightarrow Q_u^{\oplus3}\longrightarrow W_u\longrightarrow0.
\tag{1.16tr16}
\]

The vanishing of
\(\operatorname{Tor}_1^R(R/(u-f),W_u)\) identifies the base change of \(J_u\)
with the intrinsic annihilator \(J_f\).  Reducing (1.16tr14) modulo the graph
and then modulo the collision maximal ideal sends \(u\) to zero and gives
(1.16tr9).  Likewise the graph associated graded of \(Q_u\) is the saturated
central cotangent module, proving (1.16tr8).  This proves the theorem.

The pinned certificate is
[`nodal_all_orders_different_persistence.json`](../artifacts/generated-results/nodal_all_orders_different_persistence.json),
replayed by
[`verify_nodal_all_orders_different_persistence.py`](../scripts/verify_nodal_all_orders_different_persistence.py).
It independently replays (1.16tr10), (1.16tr12)--(1.16tr15), and the
six-generator calculation.  The recursive normal-form and monic-graph
arguments are the formal deductions above.  Normality, algebraization of the
infinite gauge, and compatibility with a distinguished Keller open remain
open.

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

### Theorem 1.8f.3 -- all singular-squarefree formal tails are non-Cartier

Let \(h\) be any of the six singular squarefree ternary-cubic symbols in
Exact computation 1.8f, and let \(\phi_h\) be its compatible cubic tensor.
Every compatible formal tensor

\[
 \Phi=\phi_h+\text{terms of collision order at least four}
\]

is determinant-twisted formally gauge-equivalent to

\[
 \phi_h+\sum_{i=1}^{r_h} f_i(x,y,z)\eta_{h,i},
 \qquad f_i\in(x,y,z)k[[x,y,z]],
\tag{1.16tr17}
\]

where the minimal generator counts and one exact choice of source cubics
for the tensors \(\eta_{h,i}\) are

\[
\begin{array}{c|c|c}
h&r_h&\text{source cubics}\\ \hline
\text{nodal}&1&Z^3\\
\text{cuspidal}&2&Z^3,\ XZ^2\\
\text{line + transverse conic}&2&Y^3,\ X^3\\
\text{line + tangent conic}&3&Y^2Z,\ Y^3,\ XY^2\\
\text{triangle}&3&Z^3,\ Y^3,\ X^3\\
\text{concurrent lines}&4&Z^3,\ YZ^2,\ Y^2Z,\ XZ^2.
\end{array}
\tag{1.16tr18}
\]

For the completed cubic algebra \(B_\Phi\), put

\[
 Q_\Phi=\Omega_{B_\Phi/\widehat A},\qquad
 J_\Phi=\operatorname{Ann}_{B_\Phi}(Q_\Phi),\qquad
 \mathfrak n_\Phi=(x,y,z,e_1,e_2,e_3).
\]

Then

\[
 H^0_{(x,y,z)}(Q_\Phi)=0,\qquad
 \boxed{\dim_kJ_\Phi/\mathfrak n_\Phi J_\Phi=6.}
\tag{1.16tr19}
\]

In particular, no compatible formal higher tail in any singular-squarefree
row can make the intrinsic Kähler different Cartier.

#### Proof

The exact pruned presentations of
\(\ker C/\operatorname{im}G_h\), after exact changes of minimal quotient
generators, are

\[
\begin{array}{c|c}
h&P_h\\ \hline
\text{nodal}&(x)\\[2pt]
\text{cuspidal}&
\begin{pmatrix}3x&0\\-z&x\end{pmatrix}\\[8pt]
\text{line + transverse conic}&
\begin{pmatrix}0&y\\z&0\end{pmatrix}\\[8pt]
\text{line + tangent conic}&
\begin{pmatrix}
0&0&-6y&3y^2\\
y&0&-2z&0\\
z&y&x&-z^2
\end{pmatrix}\\[12pt]
\text{triangle}&
\begin{pmatrix}0&0&x\\0&y&0\\z&0&0\end{pmatrix}\\[8pt]
\text{concurrent lines}&
\begin{pmatrix}
3x&0&0&0\\
-y&2x&0&0\\
-z&0&2x&0\\
0&-2y+z&y-2z&x
\end{pmatrix}.
\end{array}
\tag{1.16tr20}
\]

Direct reduction verifies both

\[
 \ker C=\operatorname{im}G_h+\sum_iA\eta_{h,i}
\tag{1.16tr21}
\]

and equality of the quotient with the displayed presentation.  The
annihilators read from the same exact modules are respectively

\[
 (x),\ (x^2),\ (yz),\ (y^3),\ (xyz),\ (x^3),
\tag{1.16tr22}
\]

agreeing with the gauge-cokernel atlas.  Apply (1.16tr21) to the first
unremoved homogeneous remainder and iterate as in the proof of Theorem
1.8f.2.  This gives (1.16tr17).  Relations among the generators make the
coefficients \(f_i\) nonunique, but do not obstruct existence.

For each row introduce independent coefficients \(u_1,\ldots,u_{r_h}\) and
form the universal normal-coefficient algebra over

\[
 R_h=k[u_1,\ldots,u_{r_h},x,y,z].
\]

Exact annihilator reduction gives, in every row,

\[
 J_h/\mathfrak n_{\mathrm{rel}}J_h
 \simeq k[u_1,\ldots,u_{r_h}]^{\oplus6}.
\tag{1.16tr23}
\]

The pruned presentation is literally the direct sum of six copies of
\((x,y,z)\), independently of all \(u_i\).  Give every \(u_i,x,y,z\)
collision weight one and use the cotangent component weights of (1.16tr15).
For both \(Q_h\) and
\(W_h=\operatorname{coker}(B_h\to Q_h^{\oplus3})\), exact Rees saturation
gives

\[
\begin{aligned}
 (N_{\mathcal R(Q_h)}:t^\infty)&=N_{\mathcal R(Q_h)},&
 \operatorname{gr}Q_h&\simeq Q_{h,0}\otimes_k k[\mathbf u],\\
 (N_{\mathcal R(W_h)}:t^\infty)&=N_{\mathcal R(W_h)},&
 \operatorname{gr}W_h&\simeq W_{h,0}\otimes_k k[\mathbf u].
\end{aligned}
\tag{1.16tr24}
\]

The support-hull \(\operatorname{Ext}^2\) presentation itself need not be
constant: for the tangent line--conic row its pruned presentation has four
noncentral reductions.  This does not enter the annihilator argument;
(1.16tr23)--(1.16tr24) are the required intrinsic packet.

Specialize successively along
\(u_i-f_i(x,y,z)\).  After the preceding graph equations have been imposed,
the initial form of the next equation is still monic in the unused variable
\(u_i\).  It is therefore regular on both associated-graded packets in
(1.16tr24).  Filtered regularity makes the graph equations a regular
sequence on \(Q_h\) and \(W_h\).  Applying the resulting Tor vanishing to
the analogue of (1.16tr16) commutes the intrinsic annihilator with the full
graph specialization.  Equation (1.16tr23) then gives the second equality
of (1.16tr19), and the first follows from the saturated central associated
graded cotangent.  This proves the theorem.

The exact data are pinned in
[singular_squarefree_all_orders_different_persistence.json](../artifacts/generated-results/singular_squarefree_all_orders_different_persistence.json)
and replayed by
[verify_singular_squarefree_all_orders_different_persistence.py](../scripts/verify_singular_squarefree_all_orders_different_persistence.py).
Normality, algebraization of an infinite formal gauge, and compatibility
with a distinguished Keller open remain outside the theorem.

### Corollary 1.8f.4 -- every squarefree cubic formal collision is non-Cartier

Let \(h\) be any squarefree ternary-cubic symbol, smooth or singular.  For
every compatible formal tensor with leading symbol \(h\), the intrinsic
Kähler different at the collision has six minimal generators and is not
Cartier.

#### Proof

The six singular rows are Theorem 1.8f.3.  In the smooth row the exact
gauge-cokernel theorem gives no positive-order normal coefficients:
successive formal gauge transformations identify every compatible tail with
the central smooth algebra.  Direct central annihilator reduction gives

\[
 H^0_{(x,y,z)}(Q_h)=0,\qquad
 \dim_kJ_h/\mathfrak n_hJ_h=6.
\tag{1.16tr25}
\]

Formal gauge isomorphism and faithful flatness of completion preserve the
minimal generator count.  This proves the corollary.

Consequently, if boundary-minimal Keller geometry supplies a Cartier Kähler
different at a cubic collision, its leading ternary-cubic symbol cannot be
squarefree.  The only remaining atlas rows are double line, triple line, and
zero; those rows must still pass the independent generically-étale and
Keller-compatibility gate before their translated quartic models are
relevant.

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

For the double-line, triple-line, and zero rows the quotient annihilator is
zero and its generic rank is respectively one, two, and four.  These rows
must first pass generically étale and Keller compatibility.  The translated
plane calculation above is not a substitute for that gate.

### Proposition 1.8g -- the smooth homogeneous defect is normal but not affine-open compatible

Let \(h\) be a smooth ternary cubic and let \(B_h=A\oplus M\) be the
homogeneous algebra of Exact computation 1.8e.  Give \(x,y,z\) degree one
and the three displayed generators of \(M\) degree two.  Then:

1. \(B_h\) is a normal integral graded domain, finite of generic rank three
   over \(A\);
2. its only nonfree \(A\)-module fiber is the reduced Koszul defect at the
   homogeneous vertex;
3. \(\operatorname{Spec}B_h\) contains no open subset isomorphic to
   \(\mathbb A^3\).

Thus \(B_h\) is an explicit global normal cubic defect model, but it cannot
be the canonical normalization of a Keller map with distinguished source
\(\mathbb A^3\).

#### Proof

The relation \(ze_1-ye_2+xe_3\) has degree three.  Formula (1.16q) makes
the trace-free part of \(e_i e_j\) quadratic in \(x,y,z\) times an \(e_k\),
while (1.16r) makes its scalar part quartic.  Hence \(B_h\) has the stated
positive grading.

Because \(B_h\) is finite over \(A\), every homogeneous prime in
\(\operatorname{Proj}B_h\) has at least one of \(x,y,z\) nonzero.  On this
locus Proposition 1.8b identifies the projective cover with

\[
 Z_h=\{([r],[q])\in\mathbb P^2\times\mathbb P^2:
       r\cdot q=0,\ h(q)=0\}.
\tag{1.16u}
\]

The second projection makes \(Z_h\) a \(\mathbb P^1\)-bundle over the
smooth plane cubic \(C_h\).  It is therefore smooth and integral.  The
punctured homogeneous spectrum is the corresponding
\(\mathbb G_m\)-torsor over \(Z_h\), so it too is smooth and integral.

As an \(A\)-module,

\[
 B_h=A\oplus M,
\]

and \(M\) is reflexive.  Thus \(B_h\) is \(S_2\), and in particular the
restriction map from \(B_h\) to the ring of sections on the punctured
spectrum is injective.  The punctured spectrum is integral, so this
injection excludes both zero divisors and nilpotents in \(B_h\).  Hence
\(B_h\) is a domain.
It is regular in codimension one because its punctured spectrum is smooth
and the vertex has codimension three.  Serre's criterion now proves
normality.  The fiber and Fitting assertions follow from Propositions
1.5--1.8 and the fixed presentation of \(M\).

It remains to exclude an affine-space open.  Base change to an algebraic
closure \(\bar k\); an affine-space open over the original field would
remain one after this base change.  Write
\(X_{\bar k}=\operatorname{Spec}(B_h\otimes\bar k)\) and remove its
codimension-three vertex.  The grading action is free there because at
least one of the degree-one coordinates \(x,y,z\) is nonzero.  Hence the
punctured spectrum is the \(\mathbb G_m\)-torsor associated with
\(\mathcal O_{Z_h}(-1)\).  Removing the zero section of a line bundle gives

\[
 \operatorname{Pic}(X_{\bar k}\setminus\{0\})
 \simeq
 \operatorname{Pic}(Z_{h,\bar k})/
 \mathbb Z[\mathcal O_{Z_h}(1)].
\tag{1.16v}
\]

Normality and the codimension-three complement identify the left side with
\(\operatorname{Cl}(X_{\bar k})\).  Since \(Z_h\to C_h\) is a
projective-line bundle,
\(\operatorname{Pic}^0(Z_{h,\bar k})\simeq
\operatorname{Pic}^0(C_{h,\bar k})\), an elliptic curve.  The group
\(\operatorname{Pic}^0(C_{h,\bar k})(\bar k)\) is not finitely generated,
and quotienting \(\operatorname{Pic}(Z_{h,\bar k})\) by the one cyclic
class in (1.16v) does not make it finitely generated.  Thus
\(\operatorname{Cl}(X_{\bar k})\) is not finitely generated.

If an open \(U\simeq\mathbb A^3_{\bar k}\) existed, the localization
sequence for divisor class groups would make
\(\operatorname{Cl}(X_{\bar k})\) generated by the finitely many
divisorial components of its complement, because
\(\operatorname{Cl}(\mathbb A^3_{\bar k})=0\).  This contradicts the
preceding non-finite-generation.  QED

The cone argument is special to the homogeneous lift.  For a general
higher tensor the exceptional \(\operatorname{Pic}^0\) is only formal
blowup data and need not algebraize in the local class group; this is
exactly the warning in E6 of the closure protocol.

Only after this constrained-lifting gate should a candidate be tested for:

1. preservation of normality inside the chosen global or bounded-degree
   ansatz (`S_2` is already forced by the underlying module `R\oplus M`);
2. a localization isomorphic to a localization of `k[u,v,w]`;
3. vanishing of `Omega_{B/A}` on that affine open;
4. the exact DVR splitting `(e,f)=(2,1)+(1,1)` over the critical divisor;
5. absence of every other boundary prime whose image is a target divisor.

For constructing an example, the nodal and cuspidal types are the first
rational candidates; the smooth family is the cleanest constrained-Bertini
test and may instead expose an elliptic-boundary obstruction.  Reducible
and nonreduced types should follow only after the irreducible cases, because
their extra exceptional components enlarge rather than simplify the
boundary-prime ledger.

There is also a geometric warning.  Ordinary trace and discriminant
identities are built into (1.16b), and the tame different/log-crepancy
ledger is codimension-one automatic.  They cannot eliminate the Koszul
defect.  The genuinely Keller-specific gates begin only after (1.16c):
rationality of the cubic function field, the colored
`(2,1)+(1,1)` boundary decomposition, and an open complement isomorphic to
`A^3`.  If such an `A^3` complement exists and the finite map is étale on
it, its Jacobian is automatically a nonzero constant because every unit of
`k[A^3]` is constant.

### Corollary 1.9 -- reduced point-flatness becomes sheet separation

Let `D` be the target divisor of the critical boundary prime and let
`A_D` be the closure in `bar X` of the unique affine prime over `D` from
Proposition 1.4.  A reduced length-four defect over a closed point `p in D`
forces the critical boundary prime and `A_D` to meet at the unique point of
the fiber over `p`.

Hence any one of the following excludes reduced minimal defects over `D`:

1. `F(A^3)` contains every closed point of `D`;
2. the critical boundary and `A_D` are separated over closed points of `D`;
3. their intrinsic scheme intersection has no fiber whose local algebra is
   the square-zero algebra in (1.16).

#### Proof

Both the critical boundary and `A_D` are finite and dominant over `D`, hence
surjective.  Proposition 1.8 says the fiber over `p` has only one point.
Both closures must contain that point.  Each of the three stated conditions
prevents the square-zero collision.  QED

This is the first direct bridge from the point-flatness defect to the
intrinsic boundary markings.  Generic height-one saturation supplies the
two sheets but cannot say whether they collide at a closed point.  The
scheme-intersection part of the Zariski--Main package can say exactly that.
The correct target is not unconditional sheet separation.  In the flat
foundational model the two sheets meet over a triple-root cubic, but the
fiber there is the curvilinear length-three algebra

\[
 \kappa(p)[\epsilon]/(\epsilon^3),
\tag{1.18}
\]

not the embedding-dimension-three square-zero algebra (1.16).  The remaining
local task is therefore to prove that the intrinsic intersection is
curvilinear/fiber-minimal.  Proposition 1.8 handles every reduced defect;
only nonreduced Fitting defects remain outside the square-zero analysis.

### Proposition 1.10 -- curvilinear fibers force point-flatness

Let `p in Y`.  If the finite algebra

\[
 B\otimes_A\kappa(p)
\tag{1.19}
\]

is generated by one element as a `kappa(p)`-algebra, then `B_p` is free of
rank three over `A_p`.  Consequently, if every closed canonical fiber is
curvilinear, then `Z_flat` is empty.

Here a finite fiber is called curvilinear when each of its local Artin
factors has embedding dimension at most one.  At the closed points relevant
to Proposition 1.1, `kappa(p)=k` is algebraically closed and infinite; a
finite product of such factors is generated by one element, so this
geometric definition implies the hypothesis above.

#### Proof

Put `R=A_p` and choose `bar t` generating `B_p/\mathfrak m_pB_p` as a
residue algebra.  Lift it to `t in B_p` and set `S=R[t]`.  The element `t`
is integral, so `S` is a finite `R`-module.  The generator assumption says

\[
 B_p=S+\mathfrak m_pB_p.
\]

Apply Nakayama's lemma to the finite module `B_p/S`; it gives `B_p=S`.

Let `K=Frac(R)`.  Since `B_p` has generic rank three and is a domain,
the minimal polynomial of `t` over `K` has degree three.  Its coefficients
are integral over `R` and lie in `K`; the regular local ring `R` is
integrally closed, so they lie in `R`.  If this monic cubic is `f(T)`, then
division by `f` shows

\[
 B_p=R[t]\simeq R[T]/(f),
\]

which is free with basis `1,t,t^2`.  The last assertion follows from
Proposition 1.1, since only closed points could be nonflat.  QED

This is the cleanest current point-flatness target.  The stronger
Cartier--Cohen--Macaulay atlas of Proposition 1.2 is sufficient but not
necessary.  It is enough to extract from the intrinsic scheme intersections
that every collision fiber is curvilinear.  The foundational triple-root
fiber (1.18) satisfies this condition, whereas every reduced defect fails it
maximally by Proposition 1.8.  Proposition 1.10 also excludes all
nonreduced and higher determinantal defects without classifying their
individual algebra structures.

### Proposition 1.11 -- the intrinsic cotangent-cyclicity test

For a closed point `p in Y`, put

\[
 C_p=B\otimes_A\kappa(p).
\]

The following conditions are equivalent:

1. the finite fiber `Spec C_p` is curvilinear;
2. `Omega_{C_p/kappa(p)}` is locally generated by at most one element;
3.
   \[
   \operatorname{Fitt}^{C_p}_1
   \left(\Omega_{C_p/\kappa(p)}\right)=C_p;
   \tag{1.20}
   \]
4.
   \[
   \bigwedge\nolimits^2\Omega_{C_p/\kappa(p)}=0.
   \tag{1.21}
   \]

Moreover,

\[
 \Omega_{C_p/\kappa(p)}
 \simeq
 \Omega_{B/A}\otimes_B C_p,
\tag{1.22}
\]

so (1.20)--(1.21) are intrinsic Fitting conditions already contained in the
scheme-theoretic Zariski--Main package.  If they hold at every closed
collision fiber, the canonical normalization is flat.

#### Proof

Decompose `C_p` into its Artin local factors.  For a factor with maximal
ideal `n` and residue field `k`, the standard cotangent identification gives

\[
 \Omega_{C_p/k}\otimes k\simeq n/n^2.
\]

Thus its embedding dimension is at most one exactly when its cotangent
module is cyclic, by Nakayama.  A finite module is locally generated by at
most one element exactly when its first Fitting ideal is the unit ideal,
proving the equivalence of 1--3.  Over each Artin local factor, a module is
cyclic exactly when its second exterior power vanishes, proving equivalence
with 4.  Formula (1.22) is base change for Kähler differentials.  The last
assertion is Proposition 1.10.  QED

The test distinguishes the two collisions without normal forms.  For the
foundational fiber `k[epsilon]/(epsilon^3)`,

\[
 \Omega\simeq k[\epsilon]/(\epsilon^2)\,d\epsilon
\]

is cyclic.  For the reduced defect `k plus V`, `V^2=0` and `dim V=3`, so

\[
 \Omega\otimes k\simeq V
\]

requires three generators.  Intrinsic curvilinearity is therefore exactly
cotangent cyclicity, not an additional chosen root coordinate.

### Proposition 1.12 -- primitive nilradical generation is equivalent

Let `N_p` be the nilradical of the finite collision algebra `C_p`.  Then
`Spec C_p` is curvilinear if and only if `N_p` is locally principal.
Consequently, a primitive conormal element which generates `N_p` on every
closed collision fiber implies point-flatness.

#### Proof

Over the algebraically closed residue field, each Artin local factor of
`C_p` has residue field `k`.  Its maximal ideal is therefore exactly its
nilradical.  The minimal number of generators of that ideal is

\[
 \dim_k N_p/N_p^2,
\]

which is the embedding dimension.  Hence the local factor is curvilinear
exactly when its nilradical is principal.  Proposition 1.10 gives the final
assertion.  QED

For the foundational collision, `N=(epsilon)` and `N^3=0`; for the reduced
defect, `N=V`, `N^2=0`, and `N` needs three generators.  Thus a generic
primitive conormal class is not by itself enough: the closed-point
saturation theorem must say that its specialization generates the whole
nilradical.  Nilradicals, their powers, and their exact nilpotency indices
are already strata data in the formal intrinsic package.

### Proposition 1.13 -- `S_2` extension of the primitive class

Set

\[
 Q=\Omega_{B/A},
\qquad
 T=B/\operatorname{Ann}_B(Q).
\tag{1.23}
\]

Suppose globally on `Spec T` (or on open neighborhoods covering every
closed collision point) that:

1. `T` is pure of dimension two and satisfies Serre's condition `S_2`;
2. `Q` has rank one, full support, and satisfies `S_1` as a `T`-module;
3. the intrinsic primitive conormal class `tau in Q` generates `Q` at every
   point of `Spec T` of codimension at most one.

Then, on those neighborhoods,

\[
 Q=T\tau
\tag{1.24}
\]

through every closed collision point.  If the neighborhoods cover every
collision point, all fiber cotangent modules are cyclic and the canonical
normalization is flat.

#### Proof

Multiplication by `tau` gives a map

\[
 T\longrightarrow Q.
\]

Condition `S_2` implies `S_1`; together with purity this says that `T` has
no embedded associated primes.  The rank-one and codimension-zero part of
assumption 3 therefore makes this map injective.  Its cokernel `K` is
supported in codimension at least two by the rest of assumption 3.

Localize at a point of that support.  The module `T` has depth two, `Q` has
depth at least one, and `K` has finite length.  In

\[
 0\longrightarrow T\longrightarrow Q\longrightarrow K\longrightarrow0,
\]

the first local-cohomology sequence begins

\[
 0=H^0_{\mathfrak m}(Q)\longrightarrow
 K\longrightarrow H^1_{\mathfrak m}(T)=0.
\]

It forces `K=0`.  Thus (1.24) holds.  At points outside `Spec T` the
relative cotangent module is zero.  Base change therefore makes every
fiber cotangent module cyclic, provided the stated neighborhoods cover
every collision point, and Propositions 1.11 and 1.10 give flatness.  QED

The hypotheses are exactly a codimension-two Hartogs package.  The more
concrete conditions

\[
 T\text{ Cohen--Macaulay of pure dimension two},\qquad
 Q\text{ rank-one and without embedded associated points over }T
\tag{1.25}
\]

imply conditions 1--2, but are not necessary.

This theorem is compatible with the foundational triple-root collision.
There `Q` is a cyclic module over its possibly nonreduced hypersurface
support, so both are maximal Cohen--Macaulay; no reduced-support assumption
is imposed.  The remaining extraction problem has now become a precise
depth statement: prove that saturation and boundary monotonicity make the
ramification support `S_2`, make its cotangent module `S_1`, and prevent a
codimension-one zero of the primitive class.  Proposition 1.14 packages the
two closed-point failures as exact `Ext` modules.

### Proposition 1.14 -- the two-`Ext` obstruction ledger

Retain (1.23).  Assume that `T` is pure of dimension two and satisfies
`S_1`, that `Q` has rank one and full support over `T`, and that `tau`
generates `Q` in codimension at most one.  Then multiplication by `tau` is
injective and its cokernel

\[
 K=Q/T\tau
\tag{1.26}
\]

has finite length.  Moreover:

1. the failure of `T` to be `S_2` at closed points is measured by
   \[
   \operatorname{Ext}^2_A(T,A);
   \tag{1.27}
   \]
2. if (1.27) vanishes, then
   \[
   \operatorname{Ext}^3_A(Q,A)
   \simeq
   \operatorname{Ext}^3_A(K,A),
   \tag{1.28}
   \]
   so `Ext^3_A(Q,A)` is the canonical dual of the primitive-generation
   defect `K`;
3. consequently,
   \[
   \boxed{
   \operatorname{Ext}^2_A(T,A)=0,\qquad
   \operatorname{Ext}^3_A(Q,A)=0
   }
   \tag{1.29}
   \]
   imply cotangent cyclicity at every collision and finite flatness of the
   canonical normalization.

#### Proof

The `S_1` hypothesis removes embedded associated primes; as in Proposition
1.13, codimension-zero generation makes `T -> Q` injective, and
codimension-one generation makes `K` zero-dimensional.

Because `A` is a regular ring of dimension three and `T` is a pure
two-dimensional finite `A`-module, local duality identifies its only
possible closed-point `S_2` deficiency with `Ext^2_A(T,A)`.  Purity already
and `S_1` remove `Ext^3_A(T,A)`.  Thus (1.27) vanishes exactly when `T`
is `S_2`.
In that case `T` is Cohen--Macaulay of codimension one over `A`, so

\[
 \operatorname{Ext}^i_A(T,A)=0\qquad(i\ne1).
\tag{1.30}
\]

Apply `Hom_A(-,A)` to (1.26).  A finite-length module over the regular
threefold has `Ext^i_A(K,A)=0` for `i<3`.  Using (1.30), the end of the
long exact sequence gives precisely (1.28).  The right side of (1.28)
vanishes if and only if `K=0`.  Conditions (1.29) therefore give
`Q=T tau`; Propositions 1.11 and 1.10 finish the proof.  QED

This is the smallest current algebraic attack on the collision gap.  It
replaces two qualitative depth assertions by two explicit modules computed
from finite presentations already present in the intrinsic package.

### Proposition 1.15 -- double saturation computes both `Ext` obstructions

Retain the hypotheses and notation of Proposition 1.14.  Let

\[
 C=T^{[2]}
 :=\operatorname{Ext}^1_A(
      \operatorname{Ext}^1_A(T,A),A)
                                                               \tag{1.31}
\]

be the canonical `S_2` hull of `T`, with its evaluation map
`iota:T -> C`.  Then `iota` is injective, is an isomorphism at every point
of `Spec T` of codimension at most one, and `C` is pure of dimension two
and satisfies `S_2`.  Put

\[
 L=C/T,\qquad K=Q/T\tau.                                      \tag{1.32}
\]

Then `L` and `K` have finite length and

\[
 \boxed{
 \operatorname{Ext}^2_A(T,A)
 \simeq \operatorname{Ext}^3_A(L,A).
 }                                                            \tag{1.33}
\]

Consequently,

\[
 \operatorname{Ext}^2_A(T,A)=0
 \quad\Longleftrightarrow\quad
 L=0
 \quad\Longleftrightarrow\quad
 T=C.                                                         \tag{1.34}
\]

After these equivalent conditions hold,

\[
 \boxed{
 \operatorname{Ext}^3_A(Q,A)
 \simeq \operatorname{Ext}^3_A(K,A),
 }                                                            \tag{1.35}
\]

and hence the two-`Ext` certificate of Proposition 1.14 vanishes exactly
when

\[
 \boxed{L=0,\qquad K=0.}                                      \tag{1.36}
\]

In words, the closed-point problem is the conjunction of two literal
saturation statements:

1. **support saturation:** the ramification support already equals its
   `S_2` hull `C`;
2. **conormal saturation:** the primitive class already generates the
   relative cotangent module.

#### Proof

For a pure `S_1` module of codimension one over a Gorenstein ring, the
codimension-one canonical bidual in (1.31) is its `S_2`-ification: the
evaluation map is injective, is an isomorphism in codimension at most one,
and its target is `S_2`.  This can also be read from the biduality spectral
sequence for `RHom_A(RHom_A(T,A),A)`.  Thus `L` is zero-dimensional.
Since `C` is a pure `S_2` module of dimension two over the regular
threefold `A`, it is Cohen--Macaulay of codimension one.  Therefore

\[
 \operatorname{Ext}^i_A(C,A)=0\qquad(i\ne1).            \tag{1.37}
\]

Apply `Hom_A(-,A)` to

\[
 0\longrightarrow T\longrightarrow C\longrightarrow L
 \longrightarrow0.
\]

The segment in cohomological degrees two and three, together with (1.37),
gives (1.33).  A finite-length module over the regular threefold has only
one nonzero `Ext`, in degree three, and its degree-three `Ext` is its
canonical dual.  It vanishes exactly when the module does.  This proves
(1.34).

When `L=0`, the module `T=C` is Cohen--Macaulay.  Applying
`Hom_A(-,A)` to

\[
 0\longrightarrow T\mathop{\longrightarrow}^{\tau}Q
 \longrightarrow K\longrightarrow0
\]

gives (1.35), exactly as in Proposition 1.14.  Canonical duality for the
finite-length module `K` proves (1.36).  QED

Proposition 1.15 is stronger operationally than merely naming the two
`Ext` modules.  Both obstruction modules are now duals of finite quotients
already visible in the formal intrinsic package.  The first quotient is
computed by the canonical bidual of the scheme-theoretic ramification
stratum; geometrically it can be compared with the finite `S_2` model
carried by the normalized critical boundary.  The second is computed by a
single cokernel of the primitive conormal class.

### Proposition 1.15a -- Cartier different versus the normalization deficiency

Work locally at a closed collision.  Let $(A,\mathfrak m)$ be a regular
local ring of dimension three and let $(B,\mathfrak n)$ be a finite local
normal $A$-algebra of dimension three.  Let

\[
 Q=\Omega_{B/A},\qquad J=\operatorname{Ann}_B(Q),\qquad T=B/J.
\]

Assume that the Kähler different is Cartier at the collision:

\[
 J=dB
\tag{1.37a}
\]

for a nonunit nonzerodivisor $d\in B$.  Then $T$ is pure and `S_1`,
is `S_2` away from the closed point, and its canonical support-hull defect
satisfies

\[
 \boxed{
 L=T^{[2]}/T
 \simeq H^1_{\mathfrak n}(T)
 \simeq
 (0:_{H^2_{\mathfrak n}(B)}d).
 }
\tag{1.37b}
\]

Local duality over $A$ consequently gives

\[
 \operatorname{Ext}^2_A(T,A)
 \simeq
 \operatorname{Hom}_A(L,E_A(k)),
\tag{1.37c}
\]

with the evident residue-field extension when
$B/\mathfrak n\ne A/\mathfrak m$.  In particular,

\[
 \boxed{
 L=0
 \quad\Longleftrightarrow\quad
 B\text{ is Cohen--Macaulay at }\mathfrak n.
 }
\tag{1.37d}
\]

Thus a Cartier codimension-one different does not by itself close `C1`.
It identifies the remaining class exactly as the $d$-torsion in the
degree-two deficiency module of the normalization.  Under the additional
CM hypothesis it does close `C1`; since $B$ is then a maximal
Cohen--Macaulay module over the regular local ring $A$, Auslander--Buchsbaum
makes $B$ finite free over $A$, so the local cubic normalization is
already flat.

#### Proof

Normality gives `S_2`, hence

\[
 H^0_{\mathfrak n}(B)=H^1_{\mathfrak n}(B)=0.
\]

The quotient of an `S_2` ring by a nonzerodivisor is `S_1`.  At every
nonclosed point of $V(d)$, the relevant localization of $B$ has
dimension at most two and is Cohen--Macaulay; its quotient by $d$ is
therefore Cohen--Macaulay.  Hence $T$ is `S_2` on the punctured spectrum.
The punctured-spectrum exact sequence for an `S_1` two-dimensional module
identifies its `S_2`-ification quotient with $H^1_{\mathfrak n}(T)$.

Apply local cohomology to

\[
 0\longrightarrow B
 \mathop{\longrightarrow}^{d}B
 \longrightarrow T\longrightarrow0.
\]

The preceding vanishing gives

\[
 H^1_{\mathfrak n}(T)
 =\ker\!\left(
 d:H^2_{\mathfrak n}(B)\longrightarrow H^2_{\mathfrak n}(B)
 \right),
\]

which proves (1.37b).  Formula (1.37c) is local duality and agrees with
Proposition 1.15.

If $H^2_{\mathfrak n}(B)=0$, normality already supplies depth at least
two, so $B$ has depth three and is Cohen--Macaulay.  Conversely, suppose
$H^2_{\mathfrak n}(B)\ne0$.  Its Matlis dual is a nonzero finite module;
Nakayama therefore gives a nonzero socle in
$H^2_{\mathfrak n}(B)$.  Because $d\in\mathfrak n$, that socle is killed
by $d$, so the last kernel is nonzero.  This proves (1.37d).  The final
flatness statement follows from Auslander--Buchsbaum.  QED

### Corollary 1.15b -- the local CM/Cartier collision fork

Under the hypotheses of Proposition 1.15a, the following are equivalent:

\[
 B\text{ is Cohen--Macaulay},\qquad
 T=B/dB\text{ is }S_2,\qquad
 T=T^{[2]},\qquad
 \operatorname{Ext}^2_A(T,A)=0.
\tag{1.37e}
\]

Thus a persistent length-six support-hull defect together with a Cartier
different implies that `B` is non-Cohen--Macaulay and that the ramification
surface `T` fails `S_2`; it is not, by itself, a contradiction.  There are
exactly two ways to close this fork from the present data:

1. prove from boundary-minimal Keller geometry that `T` is `S_2` (or prove
   local Cohen--Macaulayness of `B` directly);
2. prove that a collision carrying the persistent defect cannot have a
   locally principal Kähler different.

Exact computation 1.8f proves the second alternative on every geometric
fiber of all six complete quartic nongauge families, where the different
requires six generators.  Exact computation 1.8f.1 extends it through the
complete quintic and sextic normal-form quotients for the nodal symbol.
Theorem 1.8f.2 then closes the nodal row through every compatible formal
tail.  The higher-order rows for the other five singular-squarefree symbols
remain open.

#### Proof

Because `T` is `S_1`, is pure of dimension two, and is `S_2` away from the
closed point, it is `S_2` exactly when its canonical hull quotient `L`
vanishes.  Proposition 1.15a identifies this with Cohen--Macaulayness of
`B`, and local duality identifies it with vanishing of the displayed
`Ext^2`.  The final assertions are Exact computation 1.8f and Theorem
1.8f.2 within their stated scopes.  QED

### Proposition 1.15c -- relative persistence of a non-Cartier different

Let `S` be a connected Noetherian base and let `B` be an `S`-flat algebra
with a collision section defined by an ideal `n` such that `B/n=S`.  Let

\[
 Q=\Omega_{B/A},\qquad
 J=\operatorname{Ann}_B(Q),\qquad
 W=\operatorname{coker}\big(B\longrightarrow Q^{\oplus3}\big),
\]

where the displayed map sends `b` to `(b de_1,b de_2,b de_3)`.  Assume
that `de_1,de_2,de_3` generate `Q` as a `B`-module and that `Q` and `W`
are `S`-flat.  Then `J` commutes with arbitrary base change on `S`.  If in
addition

\[
 \mathcal G=J/\mathfrak nJ
\tag{1.37f}
\]

is finite locally free of rank `r` over `S`, then at every geometric point
`s` of `S`,

\[
 \boxed{
 \mu_{B_s,\mathfrak n_s}(J_s)
 =\dim_{k(s)}J_s/\mathfrak n_sJ_s=r.
 }
\tag{1.37g}
\]

In particular, `r>1` excludes a locally principal, hence a Cartier,
Kähler different on every fiber.  For a finite exhaustive separated good
filtration, it is enough to certify local freeness of `G` by a strict
filtered/Rees presentation whose initial module is a free rank-`r` module
over `S`.

#### Proof

Let `I` be the image of `B -> Q^3`.  Flatness of `Q` and `W` makes
`0 -> I -> Q^3 -> W -> 0` remain exact after every base change and makes
`I` flat over `S`.  Flatness of `B` and `I` then makes
`0 -> J -> B -> I -> 0` remain exact.  Thus the fiber of `J` is the
annihilator kernel in the fiber algebra, not merely the specialization of a
chosen presentation.  Right exactness identifies the fiber of (1.37f) with
`J_s/n_sJ_s`.  Local freeness gives its dimension `r`, and Nakayama's lemma
identifies that dimension with the minimal number of local generators.
The filtered criterion is the standard strict-Rees flatness criterion.
QED

The new geometric target is therefore the following **double-saturation
theorem**:

> Saturated minimal-boundary intersections force `T=C` and `Q=T tau`.

It is enough to prove this at the finite set of closed collisions.  No
global choice of a root coordinate, Tschirnhausen basis, or suspension
chart occurs in the statement.

### Proposition 1.16 -- the two saturation defects are coupled

Retain Proposition 1.15, and let `Z` be the finite closed collision locus.
Put

\[
 P=H^0_Z(Q),
\tag{1.38}
\]

the maximal submodule of the relative cotangent module supported at closed
collisions.  Then there is a canonical exact sequence

\[
 \boxed{
 0\longrightarrow P\longrightarrow K\longrightarrow L
 \longrightarrow H^1_Z(Q)\longrightarrow0.
 }
\tag{1.39}
\]

In particular,

\[
 L=0\quad\Longrightarrow\quad K\simeq P.              \tag{1.40}
\]

Consequently the two-`Ext` certificate is equivalent to

\[
 \boxed{
 L=0,\qquad P=0.
 }
\tag{1.41}
\]

Thus the closed-point theorem can be stated without an independent
generation condition after support saturation:

> The ramification support equals its canonical `S_2` hull, and the
> relative cotangent module has no closed-point torsion.

#### Proof

The exact sequence

\[
 0\longrightarrow T\longrightarrow C\longrightarrow L
 \longrightarrow0
\]

and the `S_1` condition on `T` give `H^0_Z(T)=0`; the `S_2` condition on
`C` gives `H^0_Z(C)=H^1_Z(C)=0`.  Since `L` is supported on `Z`, its
local-cohomology sequence gives a canonical isomorphism

\[
 H^1_Z(T)\simeq L.                                     \tag{1.42}
\]

Apply local cohomology to

\[
 0\longrightarrow T\mathop{\longrightarrow}^{\tau}Q
 \longrightarrow K\longrightarrow0.
\]

Here `H^0_Z(T)=0`, `H^0_Z(Q)=P`, `H^0_Z(K)=K`, and
`H^1_Z(K)=0` because `K` has finite length.  The beginning of the long
exact sequence, followed by (1.42), is exactly (1.39).

If `L=0`, (1.39) gives (1.40).  Proposition 1.15 says that vanishing of the
two `Ext` modules is equivalent to `L=K=0`; by (1.40), this is equivalent
to `L=P=0`.  QED

Proposition 1.16 changes the preferred proof strategy.  It is no longer
necessary to prove closed-point generation of `tau` directly.  Prove:

1. the canonical bidual map `T -> T^[2]` is onto;
2. `Omega_{B/A}` has no associated prime supported at a closed collision.

Then (1.39) forces `Q=T tau`, and the cubic normalization is flat.
Both statements are presentation-theoretic: the first is a canonical
bidual cokernel, while the second is the zeroth local cohomology or,
equivalently, the finite-length torsion of a Jacobian cokernel.

### Proposition 1.17 -- the finite presentation certificate

Let

\[
 F_1\mathop{\longrightarrow}^{\Psi}F_0\longrightarrow Q
 \longrightarrow0
\tag{1.43}
\]

be a finite free `A`-presentation, put `N=im(Psi)`, and let

\[
 I=\operatorname{Fitt}^A_3(B).
\tag{1.44}
\]

The support of `I` is the finite point-flatness locus.  The closed-point
cotangent torsion is

\[
 \boxed{
 H^0_I(Q)\simeq (N:_{F_0}I^\infty)/N.
 }
\tag{1.45}
\]

Therefore Certificate E vanishes if and only if the following two finite
presentation tests pass:

\[
\boxed{
\begin{aligned}
 T&\longrightarrow
 \operatorname{Ext}^1_A(
   \operatorname{Ext}^1_A(T,A),A)
 &&\text{is onto},\\
 N:_{F_0}I^\infty&=N.
\end{aligned}}
\tag{1.46}
\]

#### Proof

An element of `Q=F_0/N` belongs to `H^0_I(Q)` exactly when it is
annihilated by some power of `I`.  A representative `v in F_0` has this
property exactly when

\[
 I^n v\subset N
\]

for some `n`, which is the definition of
`v in N:_{F_0}I^\infty`.  This proves (1.45).  Proposition 1.15 identifies
surjectivity of the first map in (1.46) with `L=0`; Proposition 1.16 then
identifies the second equality with `P=0` and proves the equivalence.  QED

Formula (1.45) is the module form of the general
[support-saturation principle](../verified/SUPPORT_SATURATION_PRINCIPLE.md).
For a finite module \(M\) and an ideal \(\mathfrak a\), that principle
identifies \(H^0_{\mathfrak a}(M)=0\) simultaneously with positive
grade, avoidance of associated primes, and presentation saturation.  Its
[`S_1` boundary theorem](../verified/SUPPORT_SATURATION_PRINCIPLE.md#the-geometric-depth-theorem)
now gives a restrictive structural shortcut here.  After the first row of
(1.46) gives `T=T^[2]`, it is enough to prove that `Q=Omega_{B/A}` is
`S_1`: the finite collision locus has relative height two on the pure
two-dimensional full support, so \(H_I^0(Q)=0\).  Normality of `B` does not
itself prove this module-depth condition; a Cohen--Macaulay/perfect
presentation would certify it.

This is the current smallest machine-checkable flatness certificate.  It
requires no enumeration of collision algebras: compute one canonical
bidual and one module saturation.  The exact Singular calibration in
`scripts/verify_cubic_double_saturation.sing` separates a pure surface
summand from a single closed-point cotangent summand.

#### Integrated support-saturation gate

The remaining point-flatness path is now frozen as the `C0`--`C3` row of
[`SUPPORT_SATURATION_PATHS.json`](../verified/SUPPORT_SATURATION_PATHS.json):

| Stage | Repository obligation | Accepted certificate | Current status |
|---|---|---|---|
| `C0` / `G0` | Construct \(Q=\Omega_{B/A}=F_0/N\) and \(I=\operatorname{Fitt}_3^A(B)\) | Proposition 1.17 and the checked finite presentation | proved |
| `C1` / `G2` | Prove \(T\to T^{[2]}\) is onto, hence the finite collision locus has positive relative height on the pure full support of \(Q\) | canonical bidual cokernel, or equivalent \(S_2\)-hull certificate | open |
| `C2` / `G3` | Exclude embedded collision support in \(Q\) | \(S_1\), CM/perfectness, associated-prime avoidance, a regular element in \(I\), or direct colon equality | open |
| `C3` / `G4` | Kill the cotangent defect | apply `SST1` to get \(H_I^0(Q)=0\), then use Proposition 1.16 | conditional on `C1`--`C2` |

There is no separate `G1` stage here: the target is vanishing of the whole
module (H_I^0(Q)), not a distinguished class whose support must first be
localized.

The failure output is also prescribed.  A nonzero `C1` cokernel is retained
as \(L=T^{[2]}/T\).  A minimal component of `Supp(Q)` inside `V(I)` fails
relative height and cannot be saturated away.  An embedded associated prime
containing `I` must be reported with its torsion witness.  None of these
outcomes may be replaced by generic local freeness or by raising `I` to a
higher power.

<!-- status-consumer: SST1 12c5cb15e8b6de26 -->

### Proposition 1.18 -- the grade-one perfect shortcut

Retain Proposition 1.14.  Suppose `T` and `Q` admit balanced finite free
presentations

\[
\begin{aligned}
 0&\longrightarrow A^r\mathop{\longrightarrow}^{\Phi}A^r
   \longrightarrow T\longrightarrow0,\\
 0&\longrightarrow A^s\mathop{\longrightarrow}^{\Psi}A^s
   \longrightarrow Q\longrightarrow0,
\end{aligned}
\tag{1.47}
\]

with `det(Phi)` and `det(Psi)` nonzero.  Then Certificate E vanishes and
the canonical cubic normalization is flat.

It is enough, more invariantly, that `T` and `Q` be perfect `A`-modules of
grade one.

#### Proof

The nonzero determinants make both displayed left maps injective over the
domain `A`.  Thus `T` and `Q` have projective dimension one.  Dualizing
(1.47) gives

\[
 \operatorname{Ext}^i_A(T,A)=
 \operatorname{Ext}^i_A(Q,A)=0
 \qquad(i\ge2).
\]

In particular the two modules in (1.29) vanish, and Proposition 1.14 gives
flatness.  The perfect grade-one formulation is equivalent locally to such
balanced resolutions.  QED

Under the one-critical-prime hypothesis, the radicals of the two
determinants in (1.47) are supported on the single branch equation.
Adjugate matrices then turn these presentations into matrix factorizations
of suitable powers of that equation.  This is the strongest convenient
target for the different/conductor attack; Proposition 1.17 remains
available when a perfect presentation cannot be extracted globally.

## 2. From a binary cubic to normalized factorization

### Proposition 2.1 -- affine Hartogs maximality

Let `X` be a normal separated integral variety and let `U=A^n` be a dense
open.  If `X\setminus U` has codimension at least two, then `U=X`.

#### Proof

Normal Hartogs extension gives

\[
 \Gamma(X,\mathcal O_X)=\Gamma(U,\mathcal O_U)
 =k[x_1,\ldots,x_n].
\]

The coordinate functions define a morphism `r:X -> U` whose restriction to
`U` is the identity.  If `j:U -> X` is the open immersion, the morphisms
`jr,id_X:X -> X` agree on dense `U`.  Their equalizer is closed because `X`
is separated, so they agree everywhere.  Thus `j` and `r` are inverse
isomorphisms.  QED

Apply this to the étale locus `V` of the canonical finite normalization.
The Keller open `U=A^3` lies in `V`.  If the intrinsic divisor ledger records
every boundary valuation and every recorded boundary divisor is critical,
then `V\setminus U` has no codimension-one component.  Proposition 2.1
forces `U=V`.

For a finite flat cubic written as a binary cubic, `V` is exactly the full
simple-marked-root locus: the relative differential at the marked root
vanishes precisely when that root is repeated.  Consequently the “no extra
simple boundary” assumption below is automatic from the exhaustive minimal
divisor ledger; it is retained in the theorem so the statement also applies
without that ledger.

### Proposition 2.1a -- the phantom-boundary factor

Let `S_F` be the nonproperness set of `F` and let

\[
 D_F=\pi\!\left(\operatorname{Supp}\Omega_{B/A}\right)
\tag{2.0a}
\]

be the branch set of the canonical finite normalization.  Then

\[
 S_F=\pi(\partial_F),\qquad D_F\subseteq S_F.
\tag{2.0b}
\]

In codimension one, the irreducible components of
`S_F\setminus D_F` are in bijection with the boundary primes at whose
generic points `pi` is unramified.

For a cubic package with one critical boundary prime, `D_F` has one
irreducible divisorial component.  If `delta_F` is its reduced equation and
`j_F` is a reduced equation for `S_F`, then

\[
 j_F=\delta_F\,u_F.
\tag{2.0c}
\]

The following are equivalent:

1. there is no second unramified boundary divisor;
2. the **phantom-boundary factor** `u_F` is a unit;
3. `S_F` and `D_F` agree in codimension one.

In particular, irreducibility of `S_F` closes the unramified-boundary gap.

#### Proof

The restriction `F=pi|_U` fails the valuative criterion for properness
exactly over the finite image of `bar X\setminus U`, giving the first
identity in (2.0b).  The Keller condition makes `pi` étale on `U`, so the
support of its relative cotangent module lies in the boundary and gives the
second inclusion.

Let `E` be a boundary prime mapping to a target divisor.  Its image belongs
to `D_F` exactly when `Omega_{B/A}` is nonzero at the generic point of
`E`, equivalently when the corresponding DVR extension is ramified.  This
proves the codimension-one correspondence.

In geometric degree three, Proposition 1.4 shows that the unique critical
prime has `(e,f)=(2,1)` and exhausts the ramification over its image.
Hence its image is the unique irreducible branch divisor.  Since
`A=k[Y]` is factorial, reduced equations `j_F,delta_F` exist, and
`delta_F` divides `j_F`.  Their quotient records exactly the remaining
divisorial images.  This proves all three equivalences.  QED

Both factors in (2.0c) are intrinsic and computable.  The branch equation is
obtained from the discriminant/different or the Fitting support of
`Omega_{B/A}`; the nonproperness equation is obtained by eliminating the
graph at infinity.  Thus the second-divisor problem is no longer a search
over compactifications: it is the single unit test `u_F in k^*`.

### Corollary 2.1b -- the phantom factor is absent under boundary minimality

If \(F\) is boundary-minimal among nonproper geometric-degree-three Keller
maps in the sense of Definition 2.3 of
[`MINIMAL_BOUNDARY_CLASSIFICATION.md`](MINIMAL_BOUNDARY_CLASSIFICATION.md),
then \(S_F\) is irreducible and \(u_F\in k^\times\).

Indeed, the foundational cubic is a degree-three competitor whose
nonproperness equation is the irreducible branch equation verified above.
The first lexicographic entry of the boundary-minimality invariant is
therefore at most one, while nonproperness makes it at least one.  Hence
\(S_F\) has one irreducible divisorial component.  Since \(D_F\) is a
nonempty divisor contained in \(S_F\), Proposition 2.1a gives the claim.

### Proposition 2.2 -- there is no global monogenic shortcut

Assume the canonical normalization is finite flat of degree three and its
étale locus is exactly `U=A^3`.  Then `B` is not generated by one element as
an `A`-algebra.

#### Proof

Suppose `B=A[t]`.  Since `t` has degree three over `Frac(A)`, integrality and
normality of `A` give a monic cubic `f(T) in A[T]` with

\[
 B\simeq A[T]/(f).
\]

The relative differential module is

\[
 \Omega_{B/A}\simeq B/(f'(t))\,dt,
\]

so the étale locus is the principal open `D(f'(t))`.  Its restriction to
`U` is therefore a unit.  But

\[
 \Gamma(U,\mathcal O_U)^*=k^*,
\]

hence `f'(t)=c` for some `c in k^*` in the common function field.  The
nonzero polynomial `f'(T)-c` has degree two and annihilates `t`, contradicting
that its minimal polynomial has degree three.  QED

Thus Proposition 1.10 is deliberately local.  Curvilinear generators prove
flatness point by point, but they cannot be patched into one global affine
root coordinate once the Keller open is the full étale locus.  The
binary-cubic `P^1` of root directions, and therefore its transition
markings, are genuinely necessary.

Call the intrinsic binary-cubic orbit **affine-linear of full rank** if some
Tschirnhausen basis makes `a,b,c,d` affine-linear functions on `Y` and their
affine span has dimension three.  Equivalently, one representative
`kappa_F` is an isomorphism from `Y` onto an affine hyperplane

\[
 H_\ell=\{C\in\operatorname{Sym}^3(k^2):\ell(C)=1\}
\tag{2.1}
\]

for a nonzero `ell in Sym^3(k^2)^*`.

Assume in addition that the distinguished affine open `U` is the full
simple-marked-root locus of the cubic algebra.  This condition is intrinsic:
the omitted divisor is exactly the Fitting support of relative
differentials of the finite cubic cover, with no extra unramified boundary
prime.

On that simple-root locus the cubic factors projectively as

\[
 C=LQ,\qquad \deg L=1,\quad\deg Q=2,
\tag{2.2}
\]

with `Res(L,Q)` nonzero.  The line of marked factors is a line bundle on
`U`.  Since

\[
 \operatorname{Pic}(\mathbb A^3)=0,
\qquad
 k[\mathbb A^3]^*=k^*,
\tag{2.3}
\]

one may choose `L,Q` globally and use their relative scaling to normalize

\[
 \operatorname{Res}(L,Q)=1.
\tag{2.4}
\]

The transformation law

\[
 \operatorname{Res}(\lambda L,\lambda^{-1}Q)
=\lambda\operatorname{Res}(L,Q)
\tag{2.5}
\]

shows both existence and uniqueness up to the harmless normalized scalar.
Thus `U` becomes the normalized factorization slice

\[
 X_\ell=
\{(L,Q):\operatorname{Res}(L,Q)=1,\ \ell(LQ)=1\}.
\tag{2.6}
\]

This extracts the primitive resultant marking from the finite cubic
normalization; it is not supplied as a suspension coordinate.

## 3. Affine-linear cubic classification

### Theorem 3.1 -- unique affine-linear cubic normalization

Let `F:A^3 -> A^3` be a geometric-degree-three Keller map.  Assume:

1. its canonical finite normalization has empty flatness defect (1.2);
2. its Deligne--Faddeev coefficient morphism is affine-linear of full rank;
3. its distinguished affine open is the full simple-marked-root locus.

Then `F` is polynomially left--right equivalent to the foundational map.

#### Proof

Proposition 1.3 gives the binary cubic and (2.1).  Section 2 identifies the
source map with normalized multiplication

\[
 X_\ell\longrightarrow H_\ell,\qquad(L,Q)\longmapsto LQ.
\tag{3.1}
\]

The restriction of `ell` to the twisted cubic has one of the three contact
types

\[
 (1,1,1),\qquad(2,1),\qquad(3).
\]

The cubic hyperplane-orbit theorem computes

\[
\begin{array}{c|c}
\text{contact type}&[X_\ell]\\ \hline
(1,1,1)&\mathbb L^3-\mathbb L,\\
(2,1)&\mathbb L^3,\\
(3)&\mathbb L^3-\mathbb L^2.
\end{array}
\tag{3.2}
\]

Since the distinguished source is `A^3`, only the tangent-nonosculating
type `(2,1)` is possible.  All hyperplanes of this type are in one
`PGL_2` orbit, and normalized factor rescaling restores both equations in
(2.6).  Hence all of their multiplication maps are linearly left--right
equivalent.

For the representative

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2,
\]

the equations are

\[
 a^2e-abd+b^2c=1,\qquad ad+bc=1.
\tag{3.3}
\]

The global coordinates

\[
 b=1+ay,\qquad
 c=1-\frac32ay+a^2z
\tag{3.4}
\]

together with the forced polynomial formulas for `d,e` identify (3.3)
with `A^3`; normalized multiplication is the foundational map after
diagonal source and target changes.  QED

### Corollary 3.2 -- positive labels are consequences in this frontend

Inside the hypotheses of Theorem 3.1, the positive cubic quotient and
conormal labels are automatic:

\[
 y=\frac{b-1}{a},\qquad
 z=\frac{c-1+\frac32ay}{a^2}.
\tag{3.5}
\]

They are the first two saturated coefficients in the unique
tangent-nonosculating factorization slice.  Thus the separate positive-chart
straightening assumptions in the suspension approach are unnecessary in
the affine-linear finite-normalization frontend.

### Corollary 3.3 -- curvilinear normalization gateway

Let `F:A^3 -> A^3` have geometric degree three.  Assume:

1. every closed fiber of its canonical finite normalization is
   curvilinear, equivalently cotangent-cyclic in the sense of Proposition
   1.11; it is enough that the ramification support is `S_2`, the cotangent
   module is `S_1`, and the primitive class generates in codimension one as
   in Proposition 1.13.  Equivalently, the two intrinsic obstruction
   modules in Proposition 1.14 may be required to vanish;
2. the intrinsic binary-cubic orbit, once extracted, has an affine-linear
   full-rank representative;
3. the canonical boundary has no unramified prime over a second target
   divisor.

Then `F` is polynomially left--right equivalent to the foundational map.

#### Proof

Proposition 1.10 gives finite flatness and hence the binary-cubic package of
Proposition 1.3.  The third assumption and affine Hartogs maximality identify
the Keller open with the full simple-marked-root locus.  Theorem 3.1 then
applies.  QED

This is strictly more geometric than assuming an empty Fitting defect:
curvilinearity is read directly from the finite collision algebras in the
scheme-intersection package.  It is also calibrated to keep the
foundational triple-root fiber.

## 4. What remains open

Theorem 3.1 does not construct the finite-normalization witness of
`MINIMAL_BOUNDARY_CLASSIFICATION.md` or prove its eight predicates from an
unmarked boundary-minimal map.  Within this finite-normalization frontend,
the cubic extraction problem is split into exact tests:

1. **point-flatness:** prove `Fitt_3^A(B)=A`, excluding a defect supported at
   finitely many target points.  Equivalently, prove cubic fiber-minimality
   `lambda(p)<=3` there (Proposition 1.5); Proposition 1.2 gives the
   stronger sufficient route through a Cartier--Cohen--Macaulay canonical
   boundary.  For a reduced minimal defect, Propositions 1.7--1.8 reduce
   this further to excluding the square-zero collision of the ramified and
   affine sheets while retaining the foundational curvilinear triple-root
   collision.  Most efficiently, Proposition 1.10 shows that intrinsic
   curvilinearity of all collision fibers eliminates every point defect at
   once.  The local symbol calculation inside this test is stratified:
   smooth symbols proceed directly to global algebraization and boundary
   compatibility; every displayed singular-squarefree quartic nongauge
   complement passes `C2` but fails `C1` on every geometric parameter
   fiber, so Keller geometry must exclude its constant multiplicity-six
   support defect; under a Cartier different Proposition 1.15a identifies
   this exclusion exactly with local Cohen--Macaulayness of the
   normalization; and double-line,
   triple-line, and zero symbols must pass the generically étale/Keller gate
   before saturation;
2. **coefficient linearity:** prove that the intrinsic binary-cubic orbit
   has a full-rank affine hyperplane representative modulo polynomial
   Tschirnhausen gauge;
3. **unramified boundary:** for an arbitrary cubic, exclude a recorded
   height-one boundary prime inside the étale locus.  By Proposition 1.4 it
   would have to map to a second target nonproperness divisor, distinct from
   the critical discriminant.  For a boundary-minimal cubic this is already
   excluded by Corollary 2.1b.

Proposition 2.1 shows that there is no additional codimension-two or
codimension-three version of the third obstruction: once unramified boundary
divisors are excluded, the full simple-root locus equals `U`.  If
“divisor-minimal” is strengthened to say that every recorded boundary prime
is critical, the third test is automatic.  If all three tests hold, the
one-place/two-place suspension distinction is no longer needed to prove
cubic uniqueness: Theorem 3.1 goes directly from the canonical normalization
to the foundational map.

Proposition 2.2 rules out a tempting over-simplification after flatness:
the locally monogenic collision charts cannot be replaced by a single
global monic cubic coordinate.  Their `P^1` transition is precisely where
the remaining coefficient-gauge and intrinsic-marking information lives.

For the boundary-minimal problem, the second item is the main geometric
obstruction after the closed-point saturation test.  A general finite flat
cubic algebra over `A^3` gives a nonlinear morphism
`A^3 -> A^4`; neither flatness nor the discriminant divisor alone makes its
image a hyperplane.

Here affine-linearity is always understood modulo the polynomial
`GL_2(A)` Tschirnhausen gauge.  The
[cubic gauge-straightening theorem](CUBIC_GAUGE_STRAIGHTENING.md) proves that
every invariant-parameter nonlinear slice

\[
 C_1=q-3C_0h,\qquad q\ne0,\quad h\in\ker D_+,
\]

where `D_+` is the translation locally nilpotent derivation, is
gauge-equivalent by explicit source and target polynomial automorphisms to
the tangent slice `C_1=q`; the opposite shear has the symmetric theorem.
Moreover, a variable-time upper or lower shear is a polynomial automorphism
only when its time is invariant.  Thus these theorems exhaust the
single-unipotent-shear automorphisms, rather than merely giving a broad
family of examples.  The open case is transverse to finite compositions of
these polynomial unipotent gauge orbits or fails to admit such a
Tschirnhausen reduction.

## 5. Reproduction and external inputs

Run

```bash
.venv/bin/python scripts/verify_cubic_normalization_frontend.py
.venv/bin/python scripts/verify_cubic_symbol_double_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_deformation_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_quartic_tangent_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_singular_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_three_space_saturation.py
.venv/bin/python scripts/verify_universal_cubic_cotangent_saturation.py
.venv/bin/python scripts/verify_cubic_formal_gauge_cokernel_atlas.py
.venv/bin/python scripts/compile_support_saturation_cases.py --case cubic-frontier
.venv/bin/python scripts/verify_support_saturation_compiler.py
.venv/bin/python plane-jc/cas/test_cubic_cusp_local_model.py
Singular -q scripts/verify_cubic_double_saturation.sing
```

The checker verifies the universal cubic-algebra multiplication table,
trace-zero splitting, trace discriminant, the codimension-three reflexive
module warning with its excess-length-four special fiber, the canonical
`S_2`-hull calibration, the coupling of conormal failure with point torsion,
the degree-six line-section discriminants for every reduced ternary-cubic
type, and the exact tangent-hyperplane quotient coordinates.  The Singular
checker verifies the module-saturation formula (1.45) on a pure surface
module with one closed-point cotangent summand.
The symbol-stratified checker reconstructs the generalized triple-cover
multiplication and the presentations of `T` and `Omega` for all nine
nonzero ternary-cubic orbit representatives and the zero symbol.  It runs
both saturation tests for the homogeneous tensor and for one explicit
order-four lift.  These are exact leading-model computations, not a proof
for arbitrary higher lifts.
The deformation checker works over `Q[t,x,y,z]`.  For the seven squarefree
symbols it verifies that the family obtained by scaling that order-four
tensor has uniformly saturated cotangent presentation and a
parameter-independent relative `Ext^2` presentation of multiplicity six
on the collision axis.  It does not test arbitrary order-four directions.
The quartic tangent checker tests all 24 nullspace-basis axes for every
squarefree symbol.  It verifies the four invariant family conditions in
all 168 rows and records the four rows where literal Gröbner-presentation
equality changes.  This is a spanning collection of lines, not the
universal 24-parameter family or a proof for every linear combination.
The smooth-plane checker is a four-worker calculation over
`Q[u,v,x,y,z]`.  It verifies all 276 full coordinate two-planes, including
every specialization, and proves that each relative `Ext^2` presentation
is pulled back from the parameter-plane origin.  It does not control
mixtures supported on three or more basis directions.
The singular-plane checker treats the other six squarefree symbols.  It
verifies 1,652 central pruned presentations and proves the remaining four
planes flat of rank six from their exact finite parameter-ring
presentations via `Fitt_6=(1)` and `Fitt_5=(0)`.  Thus no coordinate plane
for a squarefree cubic symbol is left open.
The smooth three-space checker is the longer four-worker continuation over
`Q[p0,p1,p2,x,y,z]`.  It verifies all 2,024 coordinate three-spaces after
pruning nonminimal free summands.  It does not control mixtures supported
on four or more basis directions.
The universal formal-gauge checker closes precisely that smooth-symbol
gap without a 27-variable saturation.  It derives the nine gauge columns
from the finite determinant-twisted action over the dual numbers, verifies
`ker(C)=im(G)+A*eta` and an explicit matrix identity
`G*L=[x*eta,y*eta,z*eta]`, constructs rational gauge lifts for all 24
quartic basis tensors, and replays central cotangent saturation.  The
result is
`H^0_(x,y,z)(Omega)=0` over
`Q[u1,...,u24,x,y,z]` for the smooth symbol.  It does not assert universal
saturation for a singular cubic symbol.
The gauge-cokernel atlas checker computes the six singular squarefree
quartic quotient dimensions and their exact annihilators.  The
`cubic-frontier` compiler adapter imports that proved atlas and emits the
annihilator-indexed search queue; it performs no new singular saturation
and asserts no Keller compatibility.
The cusp checker enumerates all nine transposition-valued \(B_3\)
representations on three letters and verifies the monic cubic root model,
the \(2+1\) Kummer model, and the curvilinear length-three cusp fiber.

External structural inputs:

- the cubic-algebra/binary-cubic equivalence over an arbitrary base:
  Wood, *Parametrizing quartic algebras over an arbitrary base*,
  [Theorem 2.1](https://msp.org/ant/2011/5-8/ant-v5-n8-p05-p.pdf);
- nonflat triple covers with prescribed reflexive trace-free sheaf,
  including the generic reduced-irreducible divisor construction in
  Theorem 7.6 and the generalized cubic tensor used in Proposition 1.8a:
  Tan, *Triple covers on smooth algebraic varieties*,
  [Theorem 7.6, pp. 162--163](https://math.ecnu.edu.cn/~sltan/Triple-Cover.pdf);
- the canonical bidual characterization of the `S_2` condition for modules:
  Hashimoto, *Canonical and n-canonical modules of a Noetherian algebra*,
  [Theorem 8.4](https://www.math.okayama-u.ac.jp/~hashimoto/paper/ncan7.pdf);
- projective modules over a polynomial ring are free:
  Suslin, *Projective modules over a polynomial ring are free*,
  [MathNet](https://www.mathnet.ru/eng/dan40545);
- maximal Cohen--Macaulay modules over regular local rings are free:
  [Stacks Project, Tag 00NT](https://stacks.math.columbia.edu/tag/00NT);
- complements of affine opens and the normal Hartogs argument used in
  Proposition 2.1:
  [Stacks Project, Tag 0BCQ](https://stacks.math.columbia.edu/tag/0BCQ).
- purity and uniqueness of the normal finite extension of a finite étale
  complement cover:
  Stacks Project, Tags
  [0BMB](https://stacks.math.columbia.edu/tag/0BMB) and
  [0EY6](https://stacks.math.columbia.edu/tag/0EY6);
- the ordinary-cusp complement braid presentation used in Proposition 1.4b
  is also recorded, with its Zariski--van Kampen calculation, in
  [`../plane-jc/JC2_QUARTIC_PACKET_FRONTIER.md`](../plane-jc/JC2_QUARTIC_PACKET_FRONTIER.md).

The normalized hyperplane-orbit theorem and its motivic exclusions are
internal:
[`../verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md`](../verified/FOUNDATIONAL_INCIDENCE_CONSTRUCTION.md).
