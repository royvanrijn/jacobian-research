# Log geometry of controlled-boundary suspensions

This note recasts the determinant ledger as a relative-canonical-divisor
identity, packages the canonical finite normalization as a marked log-crepant
cover, and proves a marked `A^1/G_m` dichotomy for the plane core.  It also
separates that plane theorem from the three-dimensional lift.  In the
reciprocal `G_m` branch, valuation straightening and a canonical locally
nilpotent derivation reduce the lift to the number of components in one
general boundary fiber.  The full `A`-adic Keller congruence forces the unique
cancellation jet in the whole source ring, producing a global slice before
any Stein or puncture marking is chosen.  The
polynomial `A^1` branch retains a separate chart-straightening problem.

The point of the reformulation is organizational.  It does not assert that
every one-boundary suspension is weighted or cancellation-type.  Rather, it
identifies three layers:

1. a formal log-Jacobian identity, valid for every suspension square;
2. a plane-core classification once the primitive boundary coordinates are
   marked;
3. a vertical-chart problem which, in the reciprocal branch, reduces to the
   Stein degree of a canonical `G_a`-quotient boundary, with affine
   modifications used only after choosing graph charts.

The completion-theoretic part of the reciprocal third layer is settled here
without a primitive residue or puncture-completeness assumption.  What
remains is the upstream extraction of the reciprocal height-one markings
from an arbitrary divisor-minimal suspension, and a polynomial-chart problem
in the weighted branch.

The upstream audit in
[`ONE_BOUNDARY_CHART_CLASSIFICATION.md`](ONE_BOUNDARY_CHART_CLASSIFICATION.md)
shows that this last sentence needs a precise definition of divisor
minimality.  If it means only one affine pole divisor, the proposed global
dichotomy is false: a positive-oriented rational chart can cancel the
Jacobian of a one-boundary core and give the identity Keller map.  The audit
proves the signed chart theorem under the additional boundary-monotonicity
hypothesis and separates normal height-one primitivity from the still-open
tangential residue markings.

Work over a characteristic-zero field `k`.  All varieties are integral and
all generically finite maps are generically separable.

## 1. The determinant ledger is a relative-canonical identity

The affine determinant statement itself is the
[boundary-cancelled incidence lemma](CONTROLLED_BOUNDARY_SUSPENSIONS.md#1-boundary-cancelled-incidence-lemma).
This section does not re-prove the three family Jacobians; it upgrades that
common ledger to arbitrary rational top forms and birational models.

Consider a commutative square of smooth equidimensional varieties and
dominant rational maps

```text
 X  -------- F -------->  Y
 |                         |
 | alpha                   | beta
 v                         v
 Z  -------- Phi ------->  T
```

with

\[
 \beta\circ F=\Phi\circ\alpha.                         \tag{1.1}
\]

For a dominant rational map `g:V dashrightarrow W`, choose nonzero rational
top forms `omega_V,omega_W` and write

\[
 g^*\omega_W=J_g\omega_V.
\]

The divisor

\[
 R_g:=\operatorname{div}(J_g)
\]

is the relative canonical divisor associated with these trivializations.
Changing either top form changes both sides by the corresponding principal
canonical divisor.  Equivalently, the collection of these divisors on all
normal birational models is the relative Jacobian `b`-divisor of the induced
function-field inclusion.  This interpretation is preferable when `g` is
rational: negative coefficients then record poles of the chart rather than
pathology.

### Proposition 1.1 -- logarithmic chain rule

For the square (1.1),

\[
 \boxed{R_\alpha+\alpha^*R_\Phi
       =R_F+F^*R_\beta.}                                \tag{1.2}
\]

If `R_Phi=rE` and `F` has constant nonzero Jacobian, then

\[
 \boxed{R_\alpha+r\alpha^*E=F^*R_\beta.}                \tag{1.3}
\]

Thus (1.3) is exactly the determinant ledger.

#### Proof

Pull back a rational top form on `T` around the two sides of (1.1).  The
ordinary chain rule gives

\[
 (J_\beta\circ F)J_F=(J_\Phi\circ\alpha)J_\alpha
\]

in `k(X)^*`.  Taking principal divisors gives (1.2).  If `J_F` is a nonzero
constant and `div(J_Phi)=rE`, this reduces to (1.3).  QED

Put

\[
 \Lambda_X:=R_\alpha+r\alpha^*E,
 \qquad
 \Lambda_Y:=R_\beta.
\]

Then (1.3) says

\[
 K_X+\Lambda_X=F^*(K_Y+\Lambda_Y)                       \tag{1.4}
\]

after compatible choices of canonical divisors.  Formula (1.4) is a crepancy
identity for marked sub-pairs.  The word "sub-pair" matters: the coefficients
of `Lambda_X,Lambda_Y` may be negative or greater than one.  For cancellation,
the pole of the rational source chart contributes a negative coefficient and
cancels the positive core ramification.  For the weighted square, the ledger
has positive coefficient three on both sides.

This chart-level statement should be distinguished from the reduced
log-crepant normalization constructed next.

## 2. The canonical normalization is a reduced log-crepant cover

Let `F:U->Y` be dominant and quasi-finite, with `U,Y` normal affine, and set

\[
 \pi:\bar X_F=\operatorname{Norm}_Y k(U)\longrightarrow Y,
 \qquad
 \partial_F=(\bar X_F\setminus U)_{\mathrm{red}}.        \tag{2.1}
\]

Let `B_F` be the reduced union of

1. the codimension-one branch divisors of `pi`; and
2. the images in `Y` of the prime components of `partial_F`.

Define the full reduced divisor over this target boundary by

\[
 D_F=(\pi^{-1}B_F)_{\mathrm{red}}.                       \tag{2.2}
\]

It is important that `D_F` contains both sorts of primes over `B_F`:

- primes contained in the distinguished affine open `U`;
- missing reconstruction primes contained in `partial_F`.

The open immersion `U -> bar X_F` colors these two sorts and is retained as
part of the object.

### Proposition 2.1 -- tame log crepancy in codimension one

Assume `Y` is smooth.  As an equality of Weil divisors, or equivalently of
rank-one reflexive sheaves,

\[
 \boxed{
 K_{\bar X_F}+D_F=\pi^*(K_Y+B_F),}                      \tag{2.3}
\]

\[
 \boxed{
 \omega_{\bar X_F}(D_F)
 \simeq
 \bigl(\pi^*\omega_Y(B_F)\bigr)^{**}.}                 \tag{2.4}
\]

Here `omega_(bar X_F)(D_F)` denotes the reflexive hull of the corresponding
rank-one sheaf.  This formulation does not require `bar X_F` to be smooth,
factorial, or `Q`-factorial.

Hence the canonical finite normalization is a reduced log-crepant cover in
codimension one.

#### Proof

It suffices to compare coefficients at every height-one prime `E` of
`bar X_F`.  Let `Z=pi(E)` and let `e=e(E/Z)`.  If `Z` is not contained in
`B_F`, then `pi` is unramified at the generic point of `E`, `e=1`, and both
boundary coefficients are zero.

Suppose `Z` is a component of `B_F`.  The extension of the generic DVR of
`Z` by the generic DVR of `E` is tame: both residue fields have
characteristic zero.  Its different exponent is therefore `e-1`.  Thus the
coefficient at `E` of

\[
 K_{\bar X_F}-\pi^*K_Y+D_F-\pi^*B_F
\]

is

\[
 (e-1)+1-e=0.                                           \tag{2.5}
\]

The divisor equality follows.  Both sides of (2.4) are rank-one reflexive
and agree at every height-one point, so they agree globally.  QED

At a generic boundary point the characteristic-monoid map is

\[
 \mathbb N\longrightarrow\mathbb N,
 \qquad 1\longmapsto e.                                \tag{2.6}
\]

This records the ramification index logarithmically.  It does not record the
residue-field extension, whether `E` lies in the reconstruction open, or a
nonreduced intersection with another boundary component.  Those are genuine
decorations, not consequences of the reduced divisorial log structure.

### Proposition 2.2 -- the logarithmic Jacobian is a unit

At the generic point of every component `E` of `D_F`, the determinant of the
map on logarithmic differentials

\[
 \pi^*\Omega_Y^1(\log B_F)
 \longrightarrow
 \Omega_{\bar X_F}^1(\log D_F)                         \tag{2.7}
\]

is a unit.  Equivalently, the logarithmic ramification divisor of the finite
map of pairs defined by (2.1)--(2.2) vanishes in codimension one.

#### Proof

Let `Z=pi(E)`, choose a target uniformizer `t` and an upstairs uniformizer
`u`, and put `e=e(E/Z)`.  In the upstairs DVR,

\[
 t=a u^e,
 \qquad a\in\mathcal O_{\bar X_F,E}^*.                 \tag{2.8}
\]

Therefore

\[
 \frac{dt}{t}
 =e\frac{du}{u}+\frac{da}{a}.                          \tag{2.9}
\]

The coefficient of `du/u` is congruent to the nonzero scalar `e` modulo
`u`.  The residue extension is separable in characteristic zero, so the
tangential differential determinant is also a unit at the generic point.
Thus the full log determinant is a unit.  Formula (2.3) is its determinant
line-bundle identity.  QED

This is the literal logarithmic version of the ledger: the ordinary
Jacobian contributes the zero of order `e-1`, while the source and target
boundary generators contribute `1-e`, leaving a logarithmic unit.

## 3. The decorated log-normalization object

Define `L(F)` to consist of the finite map of reduced pairs

\[
 (\bar X_F,D_F)\longrightarrow(Y,B_F),                  \tag{3.1}
\]

together with the following markings.

### Divisorial vertices

For every prime `E` over a prime `Z`, retain

\[
 E\longrightarrow Z,
 \quad e(E/Z),
 \quad k(Z)\subset k(E),
 \quad f(E/Z),                                          \tag{3.2}
\]

the order of the different, and the color

\[
 \epsilon(E)=
 \begin{cases}
 0,&E\cap U\ne\varnothing,\\
 1,&E\subset\partial_F.
 \end{cases}                                            \tag{3.3}
\]

The complete residue morphism `E->Z` is stronger than the number `f`.

### Strata

For every collection of upstairs or downstairs boundary components, retain
their reduced and scheme-theoretic intersections, the maps between those
intersections, and, at the formal level, their completed local maps.  Contact
multiplicity and nilpotency index belong here.  They are not labels on a
single valuation.

This gives the forgetful ladder

\[
 \mathcal L^{\mathrm{formal}}(F)
 \longrightarrow
 \mathcal L^{\mathrm{sch}}(F)
 \longrightarrow
 \mathcal L^{\mathrm{red}}(F).                          \tag{3.4}
\]

The reduced incidence or dual complex is a combinatorial shadow of
`L^red(F)` and should not be added as independent data.

### Proposition 3.1 -- stable base change

For every `s>=0`,

\[
 \mathcal L(F\times\operatorname{id}_{\mathbb A^s})
 \simeq
 \mathcal L(F)\times(\mathbb A^s,\varnothing).          \tag{3.5}
\]

Under this identification:

1. every prime becomes `E times A^s` or `Z times A^s`;
2. every valuation is its Gauss extension and has value zero on the new
   variables;
3. `e`, `f`, the residue morphism, different order, and color `epsilon` are
   unchanged;
4. every stratum ring `R` becomes `R[t_1,...,t_s]`;
5. reducedness and the exact nilpotency index are unchanged;
6. the reduced incidence complex is canonically identical;
7. unit lattices of principal boundary complements and their orientation maps
   are unchanged.

#### Proof

Assertions 1--6 are direct applications of the
[stable normalization functoriality theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md):
normalization, height-one Gauss valuations, relative differentials and
Fitting ideals, scheme intersections, and nilradicals all commute with the
polynomial base change.  In particular,

\[
 \operatorname{Nil}(R[t_1,\ldots,t_s])
 =\operatorname{Nil}(R)[t_1,\ldots,t_s],                \tag{3.6}
\]

and the same equality holds for every power of the nilradical.  The
distinguished open also base-changes, so the color is unchanged.

For assertion 7, if `R` is a domain then every unit of `R[t]` already lies
in `R^*`.  Thus adjoining affine variables neither creates a boundary unit
nor changes the induced automorphism of a rank-one unit lattice.
This proves assertion 7 and the proposition.  QED

The last assertion concerns the canonical affine boundary system.  If one
instead compactifies the stabilizing `A^1` to `P^1`, an additional infinity
component appears and changes the ordinary proper boundary complex.

## 4. Rational critical curves and their punctures

The `A^1/G_m` split can be expressed intrinsically in logarithmic terms.

### Lemma 4.1 -- units and punctures

Let `k` be algebraically closed and let `E` be a smooth rational affine curve.
Write

\[
 \bar E\simeq\mathbb P^1,
 \qquad
 S=\bar E\setminus E.
\]

Then

\[
 \operatorname{rank}\bigl(\mathcal O(E)^*/k^*\bigr)
 =|S|-1.                                                \tag{4.1}
\]

If `|S|<=2`, then

\[
 |S|=1\Longleftrightarrow E\simeq\mathbb A^1,
 \qquad
 |S|=2\Longleftrightarrow E\simeq\mathbb G_m.          \tag{4.2}
\]

#### Proof

Taking divisors embeds `O(E)^*/k^*` into the degree-zero divisors supported
on `S`.  Conversely, every degree-zero divisor on `P^1` is principal.  This
identifies the unit group modulo scalars with the lattice

\[
 \{(a_p)_{p\in S}\in\mathbb Z^S:\sum_pa_p=0\},
\]

which has rank `|S|-1`.  An automorphism of `P^1` carries one puncture to
infinity and two punctures to zero and infinity, proving (4.2).  QED

The log-canonical degrees are

\[
 \deg(K_{\bar E}+S)=
 \begin{cases}
 -1,&E\simeq\mathbb A^1,\\
 0,&E\simeq\mathbb G_m.
 \end{cases}                                            \tag{4.3}
\]

Thus the weighted/cancellation split is also a split between a log-Fano
critical normalization and a log-Calabi--Yau critical normalization.

Unit rank and number of punctures are consequently redundant under the
smooth-rational hypothesis.  The useful additional information is the vector
of valuations of the distinguished reconstruction functions at the
punctures.

## 5. A marked plane-core dichotomy

The affine type of the normalized critical curve alone does not determine a
plane core.  The following theorem records exactly the primitive-function
markings which force the two cores in the repository.

### Theorem 5.1 -- marked `A^1/G_m` core normal forms

Let `K` be a characteristic-zero field.

#### The one-place case

Let

\[
 \phi:\mathbb A^2_{w,q}\longrightarrow\mathbb A^2_{q,t},
 \qquad
 \phi(w,q)=(q,T(w,q)).                                  \tag{5.1}
\]

Suppose the reduced critical divisor `E` is isomorphic to `A^1`, the
restriction of `w` is a coordinate on `E`, and

\[
 q|_E=h(w),
 \qquad
 \det D\phi=u(h(w)-q),quad u\in K^*.                  \tag{5.2}
\]

Then, after a target shear preserving `q` and a nonzero scaling of `t`,

\[
 \boxed{\phi(w,q)=(q,wq-H(w)),\qquad H'=h.}             \tag{5.3}
\]

#### The two-place case

Fix `p in K`, let

\[
 \chi:\mathbb A^2_{s,q}\longrightarrow\mathbb A^2_{q,t},
 \qquad
 \chi(s,q)=(q,T(s,q)),                                  \tag{5.4}
\]

and let `E` be its reduced irreducible critical divisor.  Suppose there is an
isomorphism `nu:G_m,Y -> E` such that, for some `m>=1`,

\[
 \nu^*s=Y^{-m},
 \qquad
 \nu^*(q-ps)=Y.                                        \tag{5.5}
\]

If

\[
 \det D\chi=-cD(s,q)^r,
 \qquad c\in K^*,\quad r\ge1,                          \tag{5.6}
\]

where `D` is a reduced equation of `E`, then, after multiplying `D` by a unit
and applying a target shear preserving `q`,

\[
 \boxed{D=1-s(q-ps)^m,}                                \tag{5.7}
\]

\[
 \boxed{
 T(s,q)=c\int_0^s\{1-v(q-pv)^m\}^r\,dv.}              \tag{5.8}
\]

Thus the marked one-place core is weighted and the marked two-place core is
cancellation-type.

#### Proof

The one-place statement is the simple-section normal form.  Since
`det Dphi=-T_w`, (5.2) integrates to

\[
 T=u(wq-H(w))+g(q),
\]

and a target shear and scaling remove `g` and `u`.

For the two-place statement, equations (5.5) give

\[
 s=Y^{-m},
 \qquad
 q=Y+pY^{-m}.                                          \tag{5.9}
\]

They identify `G_m` with the plane curve

\[
 1-s(q-ps)^m=0:                                        \tag{5.10}
\]

on that curve `q-ps` is automatically invertible and is the inverse
parameter `Y`.  Hence (5.10) is irreducible, reduced, and has the same prime
ideal as `E`, proving (5.7) up to a unit.  Finally,

\[
 \det D\chi=-T_s.
\]

Equations (5.6)--(5.7) therefore give

\[
 T_s=c\{1-s(q-ps)^m\}^r.
\]

Integration at fixed `q` gives (5.8) plus a polynomial in `q`, which a target
shear removes.  QED

The hypotheses (5.5) can be stated without choosing `Y`: the normalized
critical curve has two punctures, `q-ps` is a primitive unit with valuation
vector `(1,-1)`, and `s` has valuation vector `(-m,m)`.  Since the units of
`G_m` are scalar monomials, these intrinsic markings recover (5.5) after a
scaling and possibly interchanging the punctures.

## 6. What the curve type does not prove

The markings in Theorem 5.1 are necessary.  Two elementary enlargements show
why the bare isomorphism type of `E` cannot imply the proposed dichotomy.

### Higher-power one-place cores

For every `r>=2` and every graph `D=q-h(w)`, the coordinate-preserving map

\[
 (w,q)\longmapsto
 \left(q,-\int_0^w(q-h(v))^r\,dv\right)                 \tag{6.1}
\]

has critical normalization `A^1` and Jacobian `D^r`.  Whether a particular
such core admits a divisor-minimal polynomial Keller suspension is an extra
threefold question.  The `A^1` type alone does not force simple weighted
ramification.

### Unmarked Laurent cores

On `G_m`, one may replace the primitive function `Y` in (5.9) by a more
general Laurent polynomial `L(Y)` and consider

\[
 s=Y^{-m},
 \qquad
 q=pY^{-m}+L(Y).                                       \tag{6.2}
\]

When this parameterization is an embedding, its image still has normalized
critical curve `G_m`, but its plane equation need not be (5.7).  Unit rank and
two punctures distinguish the multiplicative branch; the primitive valuation
vector is what selects the cancellation core inside that branch.

This is not merely a hypothetical warning.  The affine-plane complement
problem admits infinite characteristic-zero families of pairwise
nonequivalent embeddings of `G_m` with isomorphic complements.  Therefore an
isomorphism of complements, the abstract curve type, and the rank-one unit
lattice cannot by themselves imply cancellation normal form.  The marked
primitive functions and the suspension ledger are essential extra data.

## 7. Numerical compression for cancellation cores

For cancellation type `(m,r)`, the known labels satisfy

\[
 N=r(m+1)+1,
 \qquad
 e_\Delta=r+1,
 \qquad
 \mu=mr(m+1).                                          \tag{7.1}
\]

Consequently, inside the cancellation family,

\[
 r=e_\Delta-1,
 \qquad
 m=\frac{N-1}{e_\Delta-1}-1,                           \tag{7.2}
\]

and

\[
 \boxed{
 \mu=(N-1)
 \left(\frac{N-1}{e_\Delta-1}-1\right).}              \tag{7.3}
\]

Thus `mu` is not independent of `(N,e_Delta)` after a map is known to be
cancellation-type.  It remains essential before classification: for
`e_Delta=2`, reduced weighted contact has `mu=1`, while cancellation contact
has

\[
 \mu=(N-1)(N-2)>1.                                     \tag{7.4}
\]

The scheme-theoretic contact is therefore a branch detector, even though it
becomes numerically redundant within the cancellation branch.

## 8. The reciprocal-boundary bridge

Theorem 5.1 classifies the marked plane cores.  It does not prove that every
divisor-minimal threefold suspension supplies the primitive markings used
there, nor that every one-boundary birational chart is triangular.

There is an important distinction between the two known families.  The
weighted vertical charts are polynomial birational morphisms.  The
cancellation source chart is not a morphism across its boundary and should
not itself be called an affine modification.

### Proposition 8.1 -- the cancellation chart is a reciprocal boundary link

For the cancellation construction put

\[
 A=1+xy^m,
 \quad
 B=A^{r+1}z+y^{m+1}h(A),
 \quad
 P=AB,
 \quad
 Q=y+xB,
 \quad
 s=x/A,                                                 \tag{8.1}
\]

and on `A^3_(s,P,Q)` put

\[
 D=1-s(Q-sP)^m.                                        \tag{8.2}
\]

Then the rational chart `alpha:(x,y,z) dashrightarrow (s,P,Q)` restricts to
an isomorphism

\[
 \boxed{
 \mathbb A^3_{x,y,z}\setminus V(A)
 \simeq
 \mathbb A^3_{s,P,Q}\setminus V(D),}                  \tag{8.3}
\]

and on this common open

\[
 \boxed{D=A^{-1}.}                                     \tag{8.4}
\]

The inverse is

\[
 y=Q-sP,
 \qquad
 x=s/D,
 \qquad
 B=PD,                                                 \tag{8.5}
\]

\[
 z=PD^{r+2}
 -(Q-sP)^{m+1}D^{r+1}h(D^{-1}).                        \tag{8.6}
\]

#### Proof

On the source,

\[
 Q-sP=y+xB-(x/A)(AB)=y,
\]

so

\[
 D=1-(x/A)y^m
   =(A-xy^m)/A=A^{-1}.                                 \tag{8.7}
\]

Equations (8.5)--(8.6) follow by solving the definitions of `s,P,Q,B` and
then `B=A^(r+1)z+y^(m+1)h(A)`.  They are regular on `D!=0`, and direct
substitution gives the identity in both directions.  QED

Thus the two boundary divisors are opposite ends of the common open, not one
exceptional divisor viewed in two affine coordinate systems.  Indeed, the
valuation of `V(A)` satisfies `v_A(D)=-1`, so it has no center on the affine
`(s,P,Q)` chart; conversely `v_D(A)=-1`.  The determinant identity

\[
 \det D\alpha=-A^r,
 \qquad
 D\circ\alpha=A^{-1}                                   \tag{8.8}
\]

is the discrepancy balance between these reciprocal ends.

This observation changes the missing classification problem.  One should
classify affine completions of the common principal open and the elementary
birational links between them, rather than assume that the rational chart is
already a one-exceptional-divisor affine morphism.

### Lemma 8.2 -- the boundary orientation is forced by the unit lattice

Let `R,S` be factorial affine `k`-domains with `R^*=S^*=k^*`, and let
`a in R`, `d in S` be irreducible.  Every isomorphism

\[
 \varphi:S[d^{-1}]\xrightarrow{\sim}R[a^{-1}]          \tag{8.9}
\]

satisfies

\[
 \boxed{\varphi(d)=c a^\epsilon,
 \qquad c\in k^*,\quad\epsilon\in\{1,-1\}.}           \tag{8.10}
\]

#### Proof

Unique factorization gives

\[
 R[a^{-1}]^*/k^*=\mathbb Z[a],
 \qquad
 S[d^{-1}]^*/k^*=\mathbb Z[d].                         \tag{8.11}
\]

Indeed, if `b/a^i` has inverse `c/a^j`, then `bc=a^(i+j)`, so every prime
factor of `b` and `c` is associated to `a`.  An isomorphism in (8.9) induces
an automorphism of the lattice `Z`, hence multiplication by `epsilon=1` or
`-1`.  This proves (8.10).  QED

Call `epsilon` the **boundary orientation**.  It is independent of rescaling
the prime equations.  The cancellation link has orientation `-1` by (8.4).
The weighted target chart restricted to `C!=0` has orientation `+1`, since
its boundary coordinate pulls back to `C` itself.  Proposition 3.1 shows that
this sign survives stabilization.

Consequently the reciprocal equation in the cancellation bridge is not an
extra coordinate hypothesis.  Once the two affine completions and their
effective boundary generators are specified, it is exactly the
orientation-reversing case of the intrinsic rank-one unit-lattice map.
The sign is only a first prefilter: exceptional complement isomorphisms show
that many inequivalent boundary embeddings can have the same unit lattice and
the same orientation.

### Proposition 8.3 -- valuation straightening produces a canonical LND

Let `k` be algebraically closed of characteristic zero, let

\[
 R=k[u,v,w],\qquad S=k[s,P,Q],
\]

and identify `S[D^(-1)]` with `R[A^(-1)]`.  Suppose

\[
 D=1-sY^m,\qquad Y=Q-sP,\qquad D=A^{-1},                \tag{8.12}
\]

where `A` and `D` are prime boundary equations.  Assume that the two primitive
boundary markings are

\[
 v_A(Y)=0,\qquad v_A(P)=1.                              \tag{8.13}
\]

Then the three elements

\[
 x=As,\qquad y=Y,\qquad B=P/A                           \tag{8.14}
\]

belong to `R` and satisfy

\[
 \boxed{A=1+xy^m,\qquad P=AB,\qquad Q=y+xB.}            \tag{8.15}
\]

Moreover,

\[
 R[A^{-1}]=k[x,y,B,A^{-1}],
 \qquad
 \det {\partial(s,P,Q)\over\partial(x,y,B)}=-A^{-1}.   \tag{8.16}
\]

If the original chart determinant is

\[
 \det {\partial(s,P,Q)\over\partial(u,v,w)}=-cA^r,
 \qquad c\in k^*,                                      \tag{8.17}
\]

then

\[
 \det {\partial(x,y,B)\over\partial(u,v,w)}=cA^{r+1}. \tag{8.18}
\]

Consequently the Jacobian derivation

\[
 \partial=\operatorname{Jac}_{u,v,w}(x,y,-)             \tag{8.19}
\]

is a nonzero locally nilpotent derivation of `R`, and

\[
 \partial(B)=cA^{r+1},\qquad
 (\ker\partial)[A^{-1}]=k[x,y,A^{-1}].                  \tag{8.20}
\]

#### Proof

Every element of `R[A^(-1)]` can have a pole only at `A`.  The valuations in
(8.13), together with `v_A(s)=-1` forced by (8.12), therefore show that the
three expressions in (8.14) are regular on `Spec R`; normality is enough for
this step.  Multiplying `1-sY^m=A^(-1)` by `A` gives

\[
 A=1+(As)Y^m=1+xy^m,
\]

and the other identities in (8.15) are immediate.  The forward and inverse
formulas in (8.14)--(8.15) prove the equality of localized rings.  A direct
Jacobian calculation gives the second identity in (8.16), whence (8.18).

On the localization, (8.18) says simply

\[
 \partial=cA^{r+1}{\partial\over\partial B}.             \tag{8.21}
\]

It is locally nilpotent there.  For every `f in R`, a sufficiently high
iterate of `partial` vanishes after localization and hence already vanishes
in the domain `R`.  Thus `partial` is locally nilpotent on `R`.  Since `A` is
in its kernel, kernels commute with this localization, which proves (8.20).
QED

This proposition extracts the cancellation skeleton without assuming a
triangular source coordinate.  The derivation (8.19) is canonical up to the
choice of source volume form.  Its fixed-point scheme is supported on `A=0`,
because `partial(B)=cA^(r+1)` is nonzero off that divisor.  Its divisorial
content, if nontrivial, is therefore a power `A^a`; after dividing by this
content one obtains the primitive LND with the same kernel.  The residual
fixed scheme of that primitive action has codimension at least two.  Thus the
ledger reduces the failure of straightening to one integer `a` and one
boundary-supported codimension-two scheme.

The integer `a` has a differential interpretation.  The coefficients of
`partial` are the `2 by 2` minors of the Jacobian matrix of `(x,y)`, so their
ideal is

\[
 \operatorname{Fitt}_1\Omega_{R/k[x,y]}.                \tag{8.21a}
\]

Its divisorial part is `A^a`.  Hence `a=0` is equivalent to generic
smoothness of the quotient map `(x,y):A^3 -> A^2` along the marked boundary.
If reconstruction-open status is required to include this generic
smoothness, primitivity already forces unit content and only the residual
codimension-two fixed scheme remains.

### Proposition 8.4 -- the marked boundary recovers the plane coordinates

Retain the hypotheses of Proposition 8.3 and put `K=ker(partial)`.  Suppose
the unique divisor `E=V_K(A)` maps dominantly to `V(A) in A^2_(x,y)`; in the
decorated-normalization language, the recorded boundary place is not
contracted and has a residue-degree label.  Then

\[
 \boxed{K=k[x,y].}                                      \tag{8.22}
\]

#### Proof

Miyanishi's kernel theorem gives `K isomorphic to k^[2]` for every nonzero
locally nilpotent derivation of `k^[3]`.  Equation (8.20) defines a birational
morphism

\[
 \phi:\operatorname{Spec}K\longrightarrow\mathbb A^2_{x,y}              \tag{8.23}
\]

which is an isomorphism off `A=0`.  The polynomial `A`, irreducible in `R`,
is also irreducible in `K`, so `E` is the only curve that could occur in a
positive-dimensional fibre.  It is not contracted by hypothesis.  Hence
`phi` is quasi-finite.

Zariski's Main Theorem, applied to the normal target with the same function
field, makes `phi` an open immersion.  Its image contains `D(A)` and the
generic point of `V(A)`, so its complement has codimension two.  Regular
functions extend across such a subset of the normal affine plane.  Therefore

\[
 K=\Gamma(\operatorname{Spec}K,\mathcal O)
  =\Gamma(\operatorname{im}\phi,\mathcal O)
  =k[x,y].                                               \tag{8.24}
\]

QED

This is exactly where the boundary labels remove the exceptional-complement
ambiguity.  An abstract isomorphism of two `G_m` complements need not extend
in either direction.  Here one direction is already the morphism (8.23), and
the noncontraction marking makes it quasi-finite.

### Proposition 8.5 -- a one-component boundary fiber forces a slice

Retain Propositions 8.3--8.4 and let

\[
 \pi=(x,y):\operatorname{Spec}R\longrightarrow\operatorname{Spec}K.
                                                               \tag{8.25}
\]

Suppose the general closed fiber of `pi` over `V_K(A)` is irreducible.  Then
the canonical derivation (8.19) is fixed-point-free and

\[
 \boxed{R=K[z]=k[x,y,z]}                                 \tag{8.26}
\]

for a polynomial slice `z`.

#### Proof

Replace `partial` by its irreducible factor; this does not change its kernel.
The quotient is surjective by the surjectivity theorem for locally nilpotent
derivations of `k^[3]`.  It is equidimensional: off `A=0` this follows from
(8.20), while a two-dimensional special fiber would contain the divisor
`V_R(A)` and contract it to a point.  This is impossible because
`(AR) intersect K=AK`, so `V_R(A)` dominates `V_K(A)`.

The source and quotient are smooth factorial affine varieties, and (8.20)
trivializes `pi` as an `A^1`-bundle over `D(A)`.  Masuda's principal-plinth
theorem says that if the plinth ideal were nontrivial, the general closed
fiber over each component of its zero divisor would be a disjoint union of at
least two affine lines.  Every plinth component must lie in `A=0`, again by
the trivialization.  The assumed irreducibility of the general boundary fiber
therefore forces the plinth ideal to be the unit ideal.  The irreducible
derivation has a slice and `R=K[z]`.

Finally `(x,y,z)` is a polynomial coordinate system.  Since polynomial
coordinate volume forms differ only by a nonzero scalar, the original
Jacobian derivation (8.19) is itself a nonzero scalar multiple of
`partial/partial z`.  It too is fixed-point-free.  QED

The hypothesis has an intrinsic boundary formulation.  Let `E=V_R(A)`, let
`L_A` be the algebraic closure of `k(V_K(A))` inside `k(E)`, and put

\[
 \sigma_A=[L_A:k(V_K(A))].                              \tag{8.27}
\]

This is the generic Stein degree of `E -> V_K(A)`.  In the present factorial
`G_a` setting, Masuda's theorem identifies it with the number of affine-line
components in a general closed boundary fiber.  Thus the hypothesis of
Proposition 8.5 is precisely

\[
 \boxed{\sigma_A=1.}                                    \tag{8.28}
\]

This **quotient Stein degree** is the relative-dimension-one analogue of the
residue-degree label for a generically finite boundary map.  It should be
included in reconstruction-open status.  It is stable under a product with a
trivial affine factor.

In fact it is already present when the marked reconstruction generator is
interpreted literally.  Let `partial_0` be the irreducible factor of
`partial`.  It descends to a nonzero LND on `E`: an irreducible LND on the
smooth threefold cannot vanish along the whole divisor.  If
`partial_0(B)` is still divisible by `A`, the residue `bar(B) in k(E)` is
invariant.  The invariant field of the induced boundary action is the
relative algebraic closure `L_A`; if `bar(B)` is the marked primitive
generator, then

\[
 L_A=k(V_K(A))(\bar B),\qquad
 \boxed{\sigma_A=f_{\bar B}}.                           \tag{8.28a}
\]

If `partial_0(B)` is not divisible by `A`, it is already a boundary slice and
`sigma_A=1`.  Consequently a primitive reconstruction residue with recorded
degree `f_(bar B)=1` implies the hypothesis of Proposition 8.5.  No additional
Stein label is needed; `sigma_A` is the correct interpretation of the
existing residue-degree/reconstruction-open decoration in relative dimension
one.

### Proposition 8.6 -- the Keller core collapses the displayed residue

Retain Propositions 8.3--8.4 and suppose that the plane-core primitive

\[
 T=C_0\int_0^s\{1-t(Q-Pt)^m\}^r\,dt                    \tag{8.28b}
\]

is regular on `Spec R`.  Then there is a constant root `q in k` of
`J_(mr,r)` such that

\[
 \boxed{\bar B=qy^{m+1}\in k(V_K(A)).}                 \tag{8.28c}
\]

If, in addition, `bar(B)` generates `L_A` over `k(V_K(A))`, then
`f_(bar B)=sigma_A=1`.

#### Proof

Put `t=su`, use `s=x/A`, `P=AB`, and `Q=y+xB`.  Since
`x=(A-1)/y^m` after localizing at `y`, equation (8.28b) becomes

\[
 {C_0x\over A^{r+1}}
 \int_0^1
 \left[A-xu\{y+xB(1-u)\}^m\right]^rdu.                \tag{8.28d}
\]

The functions `A` and `y` are coprime.  Regularity therefore forces the
constant term of the numerator at `A=0` to vanish in `k(E)`.  Writing
`b=bar(B)` and using `x=-y^(-m)` there gives

\[
 \boxed{
 J_{mr,r}\left({b\over y^{m+1}}\right)
 =\int_0^1u^r
   \left\{1-{b\over y^{m+1}}(1-u)\right\}^{mr}du=0.}   \tag{8.28e}
\]

The nonzero polynomial `J_(mr,r)` has coefficients in the algebraically
closed ground field.  Hence `b/y^(m+1)` is one of its constant roots, so
`b=qy^(m+1)` belongs to `k(V_K(A))`.  If `b` generates `L_A`, this equality
also gives `L_A=k(V_K(A))` and degree one.  QED

Equivalently, if `H(T)` is the boundary minimal polynomial of `bar(B)`, the
polynomials

\[
 H(T),\qquad
 y^{(m+1)mr}J_{mr,r}(T/y^{m+1})                        \tag{8.28f}
\]

have no common factor of degree greater than one over `k(y)`.  This is one
exact symbolic obstruction implemented by the reciprocal-link classifier.
It rules out the general `f_(bar B)>1` plinth model when the displayed
residue is primitive, not merely the quadratic example.  It does not by
itself exclude a larger relative constant field invisible to `bar(B)`.

### Proposition 8.6a -- unsliced Hensel rigidity produces the slice

Retain Proposition 8.6.  No primitivity or puncture-completeness hypothesis
on `bar(B)` is required.  Let `q` be the constant obtained there and let
`h_q(A)` be the unique degree-`r` cancellation jet with constant term `q`.
Then

\[
 \boxed{B-y^{m+1}h_q(A)\in A^{r+1}R.}                 \tag{8.28g}
\]

Consequently

\[
 z={B-y^{m+1}h_q(A)\over A^{r+1}}\in R,qquad
 \partial(z)=c,                                       \tag{8.28h}
\]

and hence

\[
 \boxed{R=k[x,y,z],\qquad\sigma_A=1.}                 \tag{8.28i}
\]

#### Proof

Work in the UFD `S=R[y^(-1)]` and put `H=B/y^(m+1)`.  Substituting
`x=(A-1)/y^m` into (8.28d) gives

\[
 T={C_0x\over A^{r+1}}\Phi_{m,r}(A,H),                 \tag{8.28j}
\]

where

\[
 \Phi_{m,r}(A,H)=\int_0^1
 \left[A+(1-A)u\{1-(1-A)H(1-u)\}^m\right]^rdu.        \tag{8.28k}
\]

The prime `A` does not divide `x`, so regularity of `T` implies

\[
 \Phi_{m,r}(A,H)\in A^{r+1}S.                         \tag{8.28l}
\]

At `A=0`, this is the spectral equation
`J_(mr,r)(bar H)=0`.  Proposition 8.6 gives `bar H=q`.  The construction of
the finite cancellation jet gives

\[
 \Phi_{m,r}(A,h_q(A))\in A^{r+1}k[A].                 \tag{8.28m}
\]

Take the polynomial divided difference

\[
 U={\Phi_{m,r}(A,H)-\Phi_{m,r}(A,h_q(A))
       \over H-h_q(A)}\in S.                          \tag{8.28n}
\]

Modulo `A`, both arguments are `q`, and therefore

\[
 \bar U=J'_{mr,r}(q)\in k^*.                          \tag{8.28o}
\]

Here squarefreeness of the cancellation parameter polynomial is essential.
Thus `v_A(U)=0`.  Equations (8.28l)--(8.28m) show that
`A^(r+1)` divides `(H-h_q(A))U`; the prime valuation now gives

\[
 H-h_q(A)\in A^{r+1}S.                                \tag{8.28p}
\]

Multiplying by `y^(m+1)` yields (8.28g) after localization.  Since
`A=1+xy^m` is coprime to `y`, `A`-adic divisibility descends from `S` to
`R`, proving (8.28g) in the polynomial source ring.

Equation (8.20) and the fact that `x,y,A` lie in `ker(partial)` now give
`partial(z)=c`.  Hence `z/c` is a slice, so the slice theorem gives
`R=(ker partial)[z]`; Proposition 8.4 identifies the kernel with `k[x,y]`.
Finally `R/AR=(k[x,y]/(A))[z]`, so the generic boundary fiber is one affine
line and `sigma_A=1`.  QED

The point is that Hensel uniqueness works in the whole nilpotent algebra
`S/A^(r+1)S`, not only in the base subring `k[y,y^(-1),A]/(A^(r+1))`.
Hidden boundary constants therefore cannot alter any coefficient of the
finite jet.  This is stronger than detecting a hidden cover after it has
formed.

### Corollary 8.6b -- rank-two boundary flags are automatic

Under Propositions 8.3--8.4 and polynomiality of the Keller core, the full
Stein normalization of the quotient boundary is

\[
 C_A=V_K(A)=\operatorname{Spec}k[y,y^{-1}].            \tag{8.28q}
\]

Its smooth completion has exactly the two flags `y=0` and `y=infinity`, and
their valuation vector is `(1,-1)`.  Thus an initially unmarked
divisor-minimal reciprocal suspension has no missing rank-two boundary flag:
the polynomial core forces the full Stein extension to be trivial before any
flag decoration is chosen.

#### Proof

Proposition 8.6a gives
`R/AR=(k[x,y]/(1+xy^m))[z]=k[y,y^(-1),z]`.  The relative algebraic closure
of `k(y)` in `k(y)(z)` is `k(y)`, so `C_A=G_m`.  Its two standard end
valuations give `(1,-1)`.  QED

### Theorem 8.7 -- completion of the reciprocal cancellation branch

In addition to Propositions 8.3--8.4, suppose that the plane-core primitive
is the polynomial

\[
 T=C_0\int_0^s\{1-t(Q-Pt)^m\}^r\,dt\in R,
 \qquad C_0\in k^*.                                    \tag{8.29}
\]

Then, after a polynomial source change and harmless nonzero scalings,

\[
 \boxed{
 A=1+xy^m,\quad
 B=A^{r+1}z+y^{m+1}h_q(A),\quad
 P=AB,\quad Q=y+xB,}                                   \tag{8.30}
\]

where `q` is a root of the cancellation parameter polynomial and `h_q` is
its unique finite cancellation jet.  Thus the reciprocal link is a
cancellation construction, not merely a link with the same critical curve.

#### Proof

Proposition 8.6a gives the polynomial slice (8.28h) directly and identifies
`R=k[x,y,z]`.  Equations (8.15) then give `P=AB` and `Q=y+xB`, so (8.30)
follows.  QED

The theorem now closes the arbitrary-rational-chart, hidden-Stein-cover, and
rank-two-flag gaps simultaneously in the reciprocal branch.  Neither
primitivity of `bar(B)`, Stein degree one, nor puncture-completeness is an
assumption.  The full `A^(r+1)`-adic Keller congruence produces the slice and
makes all three conclusions automatic.  The remaining global problem is
upstream: deriving the reciprocal link and its height-one valuation markings
from a completely arbitrary divisor-minimal polynomial suspension.

### Affine modifications after choosing a completion

An affine modification of `Spec A` with denominator `f` and center ideal `I`
is

\[
 \operatorname{Spec}A[I/f]\longrightarrow\operatorname{Spec}A. \tag{8.36}
\]

It is an affine chart of the blow-up of `(I,f)`.  The prime exceptional
divisors of the normalized blow-up are controlled by the Rees valuations of
that ideal.  Such modifications can describe affine charts in a graph
completion of (8.3), or the operation of replacing one completion boundary
by another.  They remain the natural machinery, but only after this
completion data has been chosen.

### Residual supply problem

Propositions 8.3--8.6a, Corollary 8.6b, and Theorem 8.7 close the reciprocal
bridge once its height-one link and noncontraction marking have been supplied.
More explicitly, start with an
elementary reciprocal-boundary link satisfying:

1. there are prime boundary equations `A,D` and an isomorphism
   `X\V(A) = Z\V(D)`;
2. the induced unit-lattice map has boundary orientation `-1`, equivalently
   `D=cA^(-1)` by Lemma 8.2;
3. the Jacobian `b`-divisor and core ramification satisfy the one-boundary
   ledger at the two reciprocal valuations;
4. one coordinate is preserved as the family parameter;
5. the marked boundary divisor is not contracted by the invariant-plane map;
6. both affine completions are polynomial threefolds;
7. the graph completion is divisor-minimal, with no unrecorded divisorial
   valuations or target ledger.

The valuation markings already force, before any source-coordinate
straightening,

\[
 A=1+xf(y),
 \qquad
 s=x/A,
 \qquad
 Q=y+x(P/A).                                            \tag{8.37}
\]

In the normalized `G_m` case Theorem 5.1 makes `f` a translated/scaled pure
power.  Boundary noncontraction recovers the invariant plane by Proposition
8.4.  Proposition 8.6 collapses the displayed residue `bar(B)` to the base,
and Proposition 8.6a upgrades this one coefficient to the complete
`A^(r+1)`-adic cancellation jet and a global slice.  Corollary 8.6b then
recovers the two primitive end flags automatically.

Thus complicated or nonreduced graph centers are no longer an independent
gap in the reciprocal branch: any hidden boundary splitting is incompatible
with the full Keller congruence.  The weighted branch retains its separate
polynomial-chart straightening problem.  The two branches unify at the level
of the suspension square and its Jacobian `b`-divisors, not as one category of affine birational
morphisms.

## 9. Falsification tests for the bridge

A proposed straightening theorem should be tested against the following
classes before an exhaustiveness claim is made.

1. **Boundary Stein splitting.**  Compute `sigma_A` directly for proposed
   reciprocal links.  Proposition 8.6 proves that a displayed primitive value
   greater than one is incompatible with a polynomial Keller core;
   Proposition 8.6a excludes even an extension hidden from that residue by
   forcing the complete cancellation jet and slice.
2. **The standard plinth countermodel.**  The LND
   `partial(x)=0, partial(y)=-2z, partial(z)=x^2` has kernel
   `k[x,x^2y+z^2]`; its general special fiber is two affine lines and its
   fixed locus lies over their collision.  The classifier recovers its
   quadratic residue polynomial and full Stein field exactly; Proposition
   8.6 explains why it cannot carry the cancellation Keller core with that
   residue primitive.  Displaying only a degree-one base residue leaves the
   quadratic field hidden from the spectral constant term, which the
   full-Stein stage detects and the unsliced Hensel congruence excludes.
3. **Boundary contraction.**  Drop the residue-degree/noncontraction marking
   in Proposition 8.4 and construct birational morphisms of affine planes
   which are isomorphisms off one curve but contract that curve.  This tests
   the exact necessity of the marking.
4. **Nonreduced completion centers.**  Test whether a graph-chart center with
   one radical prime but nontrivial infinitesimal structure produces precisely
   a boundary-supported fixed scheme for (8.19).
5. **Laurent deformations.**  Use (6.2) to enumerate two-place critical
   embeddings and test which admit a polynomial threefold lift.  Theorem 8.7
   shows that none survive the full `A^(r+1)`-adic congruence.
6. **Higher-power `A^1` cores.**  Test (6.1) for `r>=2` against low-weight
   polynomial and birational vertical ledgers.
7. **Nontrivial target ledgers.**  Allow `R_beta` to absorb a divisor which
   cannot be canceled by the source chart alone.  This is the most direct
   escape from the current triangular rigidity theorem.
8. **Two primitive reconstruction variables.**  The existing third-divisor
   obstruction applies only when two source valuations share one primitive
   variable.

Any surviving example would refine the hypotheses of the bridge problem;
failure across bounded families would supply evidence but not a proof.

## 10. Further structural connections

### Root-stack uniformization

The map (2.6) suggests taking an `e`-th root stack along each ramified target
component.  Tame ramification is then absorbed into the stack stabilizer and
the normalized cover becomes etale at the generic boundary points after the
appropriate base change.  A global statement requires compatible roots at
boundary intersections.  In this language the weighted discriminant carries
order two, while cancellation type `(m,r)` carries order `r+1`.  Residue
covers and the distinguished-open coloring must still be retained.

### Tropical boundary slopes

The puncture valuations of the distinguished functions form a one-dimensional
tropical curve.  In the cancellation branch,

\[
 v(s)=(-m,m),
 \qquad
 v(Y)=(1,-1).                                           \tag{10.1}
\]

The exponent `m` is therefore a tropical slope.  The determinant ledger is
an equality of valuation vectors at every boundary prime.  A tropical
prefilter may classify the possible one- and two-ended cores before affine
modification theory is invoked.

### Logarithmic and derived contact

The reduced divisorial log structure records the support and monoid maps but
forgets the thick cancellation intersection.  The scheme-theoretic stratum
ideal, or equivalently its derived intersection algebra, retains this defect.
It is natural to ask whether the nilpotency index `mu` is the length of a
canonical nonsaturated logarithmic intersection.  A positive answer would
replace an ad hoc contact label by a standard log-intersection construction.

## 11. Resulting classification target

The realistic theorem now has a proved plane statement, a proved marked
cancellation lift, and a separate weighted lift problem.

The central intrinsic formulation and its exact degree-three branch collapse
are stated in
[`MINIMAL_BOUNDARY_CLASSIFICATION.md`](MINIMAL_BOUNDARY_CLASSIFICATION.md).
The present section supplies the log-geometric proof skeleton for that
conjecture.
The one-place cubic core marking is made automatic, and the remaining
two-place defect is isolated, in
[`CUBIC_MARKING_EXTRACTION.md`](CUBIC_MARKING_EXTRACTION.md).

### Proved core statement

A coordinate-preserving one-boundary plane core with normalized critical
curve having at most two punctures is weighted or cancellation-type once the
primitive boundary-coordinate valuation vectors of Theorem 5.1 are imposed.

### Proved reciprocal cancellation lift

A divisor-minimal orientation-reversing suspension with the primitive
height-one `G_m` markings is cancellation-type provided its marked boundary
divisor is not contracted.  No primitive reconstruction residue, Stein
degree, or puncture marking is required.

Noncontraction is already part of the residue-degree/reconstruction-open
decoration.  Proposition 8.6a applies Hensel uniqueness in the entire source
algebra modulo `A^(r+1)`, produces the slice directly, and makes full Stein
degree one and the primitive two-place flags consequences.

### Remaining open suspension statement

The corrected open statement is that every **nontrivial, height-one
saturated** divisor-minimal one-boundary Keller suspension satisfies the
boundary-monotonicity and residue-saturation properties needed to admit,
after polynomial left--right equivalence, one of two straightened chart
mechanisms:

1. an orientation-preserving polynomial chart supplying the one-place
   markings of Theorem 5.1;
2. an orientation-reversing reciprocal boundary link supplying its primitive
   height-one `G_m` markings and boundary noncontraction.

Without the emphasized strengthening in this formulation, there is
a third elementary mechanism.  For `r>=1`, put

\[
 \alpha(a,y,z)=(a^{-r}y,a,z),
 \qquad
 \Phi(w,q,p)=(q,q^rw,p).
\]

Their composition is the identity.  The sole boundary is smooth, the
orientation is positive, and the two Jacobian powers cancel, but `alpha` is
rational rather than polynomial.  This example is divisor-minimal in the naïve
one-pole-support sense.  It is excluded by boundary monotonicity or by a
precise nontriviality condition, not by factoriality or smoothness alone.

Proving the open suspension statement would yield the desired dichotomy:

\[
 \begin{array}{c|c|c|c}
 \text{critical normalization}&\text{log type}&\text{chart mechanism}
 &\text{suspension type}\\ \hline
 \mathbb A^1&\text{log Fano}&\text{polynomial, orientation }+1
 &\text{weighted}\\
 \mathbb G_m&\text{log Calabi--Yau}&\text{reciprocal, orientation }-1
 &\text{cancellation}.
 \end{array}                                            \tag{11.1}
\]

The current theory therefore no longer stops at an unspecified rational-chart
classification inside the reciprocal branch: that branch is complete once
the height-one reciprocal link is present.  Rank-two flags, the full Stein
degree, and the slice are forced by the Keller core.  The remaining suspension
problem is the upstream classification of arbitrary divisor-minimal charts
into the two displayed height-one mechanisms.  The orientation-preserving
weighted branch still has its separate polynomial-chart straightening
problem.

The upstream chart note proves the following partial replacement.  If the
two affine models are factorial, their height-one link is saturated, and
positive orientation is boundary-monotone, then the unit lattice forces
exactly:

\[
 D=cA\quad\text{with a polynomial chart},
 \qquad\text{or}\qquad
 D=cA^{-1}\quad\text{with a reciprocal link}.         \tag{11.2}
\]

On a factorial graph completion, displayed exceptional primes which consume
the full relative Picard rank are the complete exceptional list.  In the
reciprocal row, a polynomial effective target ledger must vanish once the
source ledger cancels.  These statements settle the normal height-one and
hidden-divisor parts under explicit hypotheses.  They do not produce
`v_A(Q-sP)=0`, `v_A(P)=1`, or the weighted vertical formulas; those are the
residue-saturation and positive straightening frontiers respectively.

## References and repository links

- The original determinant ledger and coordinate-preserving normal forms are
  in [`CONTROLLED_BOUNDARY_SUSPENSIONS.md`](CONTROLLED_BOUNDARY_SUSPENSIONS.md).
- The intrinsic pole-support criterion, relative Picard-rank audit,
  boundary-monotonicity theorem, and the same-orientation countermodels are
  in
  [`ONE_BOUNDARY_CHART_CLASSIFICATION.md`](ONE_BOUNDARY_CHART_CLASSIFICATION.md).
- The canonical finite normalization and invariant ladder are in
  [`BOUNDARY_GEOMETRY.md`](BOUNDARY_GEOMETRY.md).
- The construction-independent stable base-change theorem is
  [`../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md`](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md).
- Rigidity after triangularization is in [`RIGIDITY.md`](RIGIDITY.md).
- The finite-cancellation operator, its leading equation, and the Laurent
  Hensel argument used in Theorem 8.7 are in
  [`../archive/cancellation-components/GENERALIZED_CANCELLATION_MECHANISM.md`](../archive/cancellation-components/GENERALIZED_CANCELLATION_MECHANISM.md).
- The exact valuation, LND, boundary-elimination, automatic full-Stein-field,
  hidden-cover, spectral-gcd, and unsliced Hensel-multiplier certificates are
  implemented in
  [`../jcsearch/reciprocal.py`](../jcsearch/reciprocal.py) and exercised by
  [`../scripts/verify_reciprocal_link_classifier.py`](../scripts/verify_reciprocal_link_classifier.py).
  Usage and certificate semantics are documented in
  [`RECIPROCAL_LINK_CLASSIFIER.md`](RECIPROCAL_LINK_CLASSIFIER.md).
- The tame logarithmic Jacobian, both marked plane cores, the reciprocal
  cancellation chart, and the numerical compression are checked exactly by
  [`../scripts/verify_log_geometry_of_suspensions.py`](../scripts/verify_log_geometry_of_suspensions.py).
- The affine-modification framework is due to
  [Kaliman--Zaidenberg](https://arxiv.org/abs/math/9801076),
  *Affine modifications and affine hypersurfaces with a very transitive
  automorphism group*, and its global form to
  [Dubouloz](https://arxiv.org/abs/math/0503142),
  *Quelques remarques sur la notion de modification affine*.
- The need for more than complement type and unit rank is reinforced by
  [Blanc--Furter--Hemmig](https://arxiv.org/abs/1609.06682),
  *Exceptional isomorphisms between complements of affine plane curves*;
  in characteristic zero they construct infinite families of nonequivalent
  `G_m` embeddings with isomorphic complements.
- Miyanishi's kernel theorem and current structure results for locally
  nilpotent derivations of `k^[3]` are summarized in
  [Dasgupta--Gaifullin](https://arxiv.org/abs/2306.00510),
  *On locally nilpotent derivations of polynomial algebra in three variables*.
  The fixed-point-free translation theorem is due to
  [Kaliman](https://arxiv.org/abs/math/0207156), and its characteristic-zero
  field form appears in
  [Daigle--Kaliman](https://doi.org/10.4153/CMB-2009-054-5).
- The principal-plinth theorem used in Proposition 8.5 is
  [Masuda](https://arxiv.org/abs/2312.05455),
  *Factorial affine `G_a`-varieties with principal plinth ideals*; its
  Corollary 3.11 turns irreducibility of the general boundary fiber into a
  trivial `A^1`-bundle.
- The quasi-finite step in Proposition 8.4 is the algebraic form of
  [Zariski's Main Theorem, Stacks Project Tag 03GS](https://stacks.math.columbia.edu/tag/03GS).
- The tame local model and the identity between tame different exponent and
  `e-1` are standard DVR facts; convenient references are the Stacks Project
  Tags [0EYF](https://stacks.math.columbia.edu/tag/0EYF),
  [0EYG](https://stacks.math.columbia.edu/tag/0EYG), and
  [0C1F](https://stacks.math.columbia.edu/tag/0C1F).
