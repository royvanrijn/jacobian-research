# Log geometry of controlled-boundary suspensions

This note recasts the determinant ledger as a relative-canonical-divisor
identity, packages the canonical finite normalization as a marked log-crepant
cover, and proves a marked `A^1/G_m` dichotomy for the plane core.  It also
separates that plane theorem from the genuinely open problem of classifying
the three-dimensional birational charts which lift the core to a polynomial
Keller map.

The point of the reformulation is organizational.  It does not assert that
every one-boundary suspension is weighted or cancellation-type.  Rather, it
identifies three layers:

1. a formal log-Jacobian identity, valid for every suspension square;
2. a plane-core classification once the primitive boundary coordinates are
   marked;
3. a reciprocal-completion straightening problem for the vertical charts,
   with affine modifications used only after choosing graph charts.

Only the third layer remains substantially open.

Work over a characteristic-zero field `k`.  All varieties are integral and
all generically finite maps are generically separable.

## 1. The determinant ledger is a relative-canonical identity

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

Normalization commutes with adjoining polynomial variables because the
polynomial ring over a normal domain is normal.  The assertions about
height-one primes and `(e,f)` follow from the corresponding Gauss valuations.
The relative differential module and its zeroth Fitting ideal commute with
this flat base change, so the different and all of its divisorial orders are
preserved as well.
Flat base change sends every sum of boundary ideals to its polynomial
extension.  Finally,

\[
 \operatorname{Nil}(R[t_1,\ldots,t_s])
 =\operatorname{Nil}(R)[t_1,\ldots,t_s],                \tag{3.6}
\]

and the same equality holds for every power of the nilradical.  The
distinguished open also base-changes, so the color is unchanged.  These are
exactly the assertions proved in the existing stable-normalization
proposition.  QED

For assertion 7, if `R` is a domain then every unit of `R[t]` already lies
in `R^*`.  Thus adjoining affine variables neither creates a boundary unit
nor changes the induced automorphism of a rank-one unit lattice.

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

### Affine modifications after choosing a completion

An affine modification of `Spec A` with denominator `f` and center ideal `I`
is

\[
 \operatorname{Spec}A[I/f]\longrightarrow\operatorname{Spec}A. \tag{8.12}
\]

It is an affine chart of the blow-up of `(I,f)`.  The prime exceptional
divisors of the normalized blow-up are controlled by the Rees valuations of
that ideal.  Such modifications can describe affine charts in a graph
completion of (8.3), or the operation of replacing one completion boundary
by another.  They remain the natural machinery, but only after this
completion data has been chosen.

### Bridge problem

Classify elementary reciprocal-boundary links between affine three-spaces
satisfying:

1. there are prime boundary equations `A,D` and an isomorphism
   `X\V(A) = Z\V(D)`;
2. the induced unit-lattice map has boundary orientation `-1`, equivalently
   `D=cA^(-1)` by Lemma 8.2;
3. the Jacobian `b`-divisor and core ramification satisfy the one-boundary
   ledger at the two reciprocal valuations;
4. one coordinate is preserved as the family parameter;
5. the reconstruction field has one marked primitive generator whose
   puncture valuation vector is primitive;
6. both affine completions are polynomial threefolds;
7. the graph completion is divisor-minimal, with no unrecorded divisorial
   valuations or target ledger.

The desired cancellation straightening conclusion is that, after source and
target automorphisms, the link has

\[
 A=1+xf(y),
 \qquad
 s=x/A,
 \qquad
 Q=y+sP.                                                \tag{8.13}
\]

Once (8.13) is reached, the existing rigidity theorem forces `f` to be a
translated/scaled pure power and forces the finite cancellation jet.  The
weighted branch has a separate polynomial-chart straightening problem.  The
two branches unify at the level of the suspension square and its Jacobian
`b`-divisors, not as one category of affine birational morphisms.

It should initially be treated as a problem, not as a theorem or unrestricted
conjecture.  Even after choosing a graph completion, affine modifications can
have complicated or nonreduced centers.  The assumptions that both affine
completions are `A^3`, that a family coordinate is preserved, and that
reconstruction is primitive must all be used.

## 9. Falsification tests for the bridge

A proposed straightening theorem should be tested against the following
classes before an exhaustiveness claim is made.

1. **Nonreduced completion centers.**  Let a graph-chart modification center
   have one radical prime but nontrivial infinitesimal structure.  Test
   whether the ledger sees only its divisorial valuation while the completed
   link retains inequivalent contact data.
2. **One Rees valuation but a nontriangular graph chart.**  Search ideals in
   `k[x,y,z]` whose normalized blow-up has one exceptional prime and whose
   relevant affine charts give two `A^3` completions of the same principal
   open.
3. **Laurent deformations.**  Use (6.2) to enumerate two-place critical
   embeddings and test which admit a polynomial threefold lift.
4. **Higher-power `A^1` cores.**  Test (6.1) for `r>=2` against low-weight
   polynomial and birational vertical ledgers.
5. **Nontrivial target ledgers.**  Allow `R_beta` to absorb a divisor which
   cannot be canceled by the source chart alone.  This is the most direct
   escape from the current triangular rigidity theorem.
6. **Two primitive reconstruction variables.**  The existing third-divisor
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

The realistic theorem can now be divided into a proved plane statement and
one open lift statement.

### Proved core statement

A coordinate-preserving one-boundary plane core with normalized critical
curve having at most two punctures is weighted or cancellation-type once the
primitive boundary-coordinate valuation vectors of Theorem 5.1 are imposed.

### Open suspension statement

Every divisor-minimal one-boundary Keller suspension with one primitive
reconstruction generator admits, after polynomial left--right equivalence,
one of two straightened chart mechanisms:

1. an orientation-preserving polynomial chart supplying the one-place
   markings of Theorem 5.1;
2. an orientation-reversing reciprocal boundary link supplying its two-place
   markings and covered by the bridge problem.

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

The current theory therefore stops at a sharply located bridge: the
logarithmic ledger, stable decorated normalization, and marked plane-core
dichotomy are formal or elementary; straightening the reciprocal completion
link, together with the separate polynomial weighted chart, is the remaining
global problem.

## References and repository links

- The original determinant ledger and coordinate-preserving normal forms are
  in [`CONTROLLED_BOUNDARY_SUSPENSIONS.md`](CONTROLLED_BOUNDARY_SUSPENSIONS.md).
- The canonical finite normalization and invariant ladder are in
  [`BOUNDARY_GEOMETRY.md`](BOUNDARY_GEOMETRY.md).
- The exact stable base-change proposition is in
  [`../papers/marked-root-multiplicity/stable-functoriality.tex`](../papers/marked-root-multiplicity/stable-functoriality.tex).
- Rigidity after triangularization is in [`RIGIDITY.md`](RIGIDITY.md).
- The tame logarithmic Jacobian, both marked plane cores, the reciprocal
  cancellation chart, and the numerical compression are checked exactly by
  [`../scripts/verify_log_geometry_of_suspensions.py`](../scripts/verify_log_geometry_of_suspensions.py).
- The affine-modification framework is due to
  [Kaliman--Zaidenberg](https://arxiv.org/abs/math/9801076),
  *Affine modifications and affine hypersurfaces with a very transitive
  automorphism group*, and its global form to
  [Dubouloz](https://arxiv.org/abs/math/0503142),
  *Quelques remarques sur la notion de modification affine*.
- The tame local model and the identity between tame different exponent and
  `e-1` are standard DVR facts; convenient references are the Stacks Project
  Tags [0EYF](https://stacks.math.columbia.edu/tag/0EYF),
  [0EYG](https://stacks.math.columbia.edu/tag/0EYG), and
  [0C1F](https://stacks.math.columbia.edu/tag/0C1F).
