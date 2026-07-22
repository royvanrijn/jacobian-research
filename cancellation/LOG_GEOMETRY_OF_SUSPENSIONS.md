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
3. an affine-modification straightening problem for the vertical charts.

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
6. the reduced incidence complex is canonically identical.

#### Proof

Normalization commutes with adjoining polynomial variables because the
polynomial ring over a normal domain is normal.  The assertions about
height-one primes and `(e,f)` follow from the corresponding Gauss valuations.
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

## 8. The affine-modification bridge

Theorem 5.1 classifies the marked plane cores.  It does not prove that every
divisor-minimal threefold suspension supplies the primitive markings used
there, nor that every one-valuation birational chart is triangular.

An affine modification of `Spec A` with denominator `f` and center ideal `I`
is

\[
 \operatorname{Spec}A[I/f]\longrightarrow\operatorname{Spec}A. \tag{8.1}
\]

It is an affine chart of the blow-up of `(I,f)`.  The prime exceptional
divisors of the normalized blow-up are controlled by the Rees valuations of
that ideal.  This converts the missing rational-chart classification into a
problem about centers, normalized blow-ups, and discrepancies.

### Bridge problem

Classify birational suspension charts between affine three-spaces satisfying:

1. the chart is an isomorphism away from one prime denominator divisor;
2. the normalized graph has one relevant Rees valuation;
3. the Jacobian `b`-divisor is supported on that valuation;
4. one target coordinate is preserved as the family parameter;
5. the reconstruction field has one marked primitive generator;
6. the source and target total spaces are both polynomial threefolds;
7. the suspension square is divisor-minimal.

The desired straightening conclusion is that, after source and target
automorphisms, the chart belongs to one of the following forms:

- the polynomial weighted charts;
- a rank-one triangular modification with
  `A=1+xf(y)` and `Q=y+sP`.

Once the second form is reached, the existing rigidity theorem forces `f` to
be a translated/scaled pure power and forces the finite cancellation jet.
Thus the bridge problem, rather than the plane determinant calculation, is
the actual missing classification theorem.

It should initially be treated as a problem, not as a theorem or unrestricted
conjecture.  Affine modifications with one exceptional valuation can still
have complicated centers.  The assumptions that both total spaces are
`A^3`, that a family coordinate is preserved, and that reconstruction is
primitive must all be used.

## 9. Falsification tests for the bridge

A proposed straightening theorem should be tested against the following
classes before an exhaustiveness claim is made.

1. **Nonreduced centers.**  Let `I` have one radical prime but nontrivial
   infinitesimal structure.  Test whether the ledger sees only its Rees
   valuation while the modification retains inequivalent contact data.
2. **One Rees valuation but a nontriangular center.**  Search ideals in
   `k[x,y,z]` whose normalized blow-up has one exceptional prime and whose
   affine modification is again factorial or isomorphic to `A^3`.
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
the normalized cover becomes etale in codimension one after the appropriate
base change.  In this language the weighted discriminant carries order two,
while cancellation type `(m,r)` carries order `r+1`.  Residue covers and the
distinguished-open coloring must still be retained.

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
reconstruction generator supplies, after polynomial left--right equivalence,
the marked coordinates of Theorem 5.1 and a vertical chart covered by the
bridge problem.

Proving the open suspension statement would yield the desired dichotomy:

\[
 \begin{array}{c|c|c}
 \text{critical normalization}&\text{log type}&\text{suspension type}\\ \hline
 \mathbb A^1&\text{log Fano}&\text{weighted}\\
 \mathbb G_m&\text{log Calabi--Yau}&\text{cancellation}.
 \end{array}                                            \tag{11.1}
\]

The current theory therefore stops at a sharply located bridge: the
logarithmic ledger, stable decorated normalization, and marked plane-core
dichotomy are formal or elementary; straightening the threefold modification
is the remaining global problem.

## References and repository links

- The original determinant ledger and coordinate-preserving normal forms are
  in [`CONTROLLED_BOUNDARY_SUSPENSIONS.md`](CONTROLLED_BOUNDARY_SUSPENSIONS.md).
- The canonical finite normalization and invariant ladder are in
  [`BOUNDARY_GEOMETRY.md`](BOUNDARY_GEOMETRY.md).
- The exact stable base-change proposition is in
  [`../papers/marked-root-multiplicity/stable-functoriality.tex`](../papers/marked-root-multiplicity/stable-functoriality.tex).
- Rigidity after triangularization is in [`RIGIDITY.md`](RIGIDITY.md).
- The affine-modification framework is due to Kaliman--Zaidenberg,
  *Affine modifications and affine hypersurfaces with a very transitive
  automorphism group*, and its global form to Dubouloz,
  *Quelques remarques sur la notion de modification affine*.
- The tame local model and the identity between tame different exponent and
  `e-1` are standard DVR facts; convenient references are the Stacks Project,
  Tags `0EYF`, `0EYG`, and `0C1F`.
