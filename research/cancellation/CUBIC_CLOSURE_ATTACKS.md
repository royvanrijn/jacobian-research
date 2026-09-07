# Cubic closure attacks

This note is the active closure protocol for the geometric-degree-three
minimal-boundary conjecture.  It replaces three qualitative gaps by three
intrinsic certificates built from the same cubic discriminant package.
The former global `OP-UG3` direction is now false: the universal flat cubic
testbed has a genuinely ungraded specialization with two phantom boundary
components.  The certificates below remain the active protocol for the
narrower boundary-minimal cubic classification.

Work over an algebraically closed field of characteristic zero.  Let

\[
 A=k[Y],\qquad
 B=\operatorname{Norm}_A k(\mathbb A^3),\qquad
 Q=\Omega_{B/A},\qquad
 T=B/\operatorname{Ann}_B(Q).
\]

Let `tau in Q` be the primitive conormal class extracted in codimension one.

## 1. The three-certificate closure criterion

For a minimal-boundary cubic Keller map, each of the following certificates
is intrinsic.

### Certificate E -- closed-point extension

\[
 \mathcal E_2=\operatorname{Ext}^2_A(T,A),\qquad
 \mathcal E_3=\operatorname{Ext}^3_A(Q,A).
\tag{1.1}
\]

Assuming `S_1` purity, rank one, and codimension-one generation:

- `E_2=0` says that the two-dimensional ramification support is `S_2`;
- once `E_2=0`, `E_3` is the canonical dual of
  `Q/T tau`;
- `E_2=E_3=0` therefore gives `Q=T tau`, cotangent-cyclic collision
  fibers, and finite flatness.

This is Proposition 1.14 of
[`CUBIC_NORMALIZATION_FRONTEND.md`](CUBIC_NORMALIZATION_FRONTEND.md).

There is now a more geometric form.  Let `C=T^[2]` be the canonical `S_2`
hull, equivalently the codimension-one canonical bidual, and set

\[
 L=C/T,\qquad K=Q/T\tau.
\tag{1.5}
\]

Under the preceding `S_1` and codimension-one hypotheses, both quotients
automatically have finite length and

\[
 \mathcal E_2\simeq\operatorname{Ext}^3_A(L,A),\qquad
 L=0\Longrightarrow
 \mathcal E_3\simeq\operatorname{Ext}^3_A(K,A).
\tag{1.6}
\]

Thus Certificate E vanishes exactly when `L=K=0`.  This is Proposition
1.15 of the frontend.  It converts E from a depth problem into a
**double-saturation problem** for two explicit scheme-theoretic cokernels.

Proposition 1.16 couples those cokernels.  If

\[
 P=H^0_Z(Q)
\tag{1.7}
\]

is the closed-point torsion in the cotangent module, then

\[
 0\to P\to K\to L\to H^1_Z(Q)\to0.
\tag{1.8}
\]

Hence `L=0` makes `K=P`, and Certificate E is equivalently

\[
 \boxed{L=0,\qquad P=0.}
\tag{1.9}
\]

The preferred E attack is therefore support saturation plus exclusion of
embedded closed points in `Q`; direct specialization of `tau` is no longer
an independent obligation.

### E0. SNC/cusp localization -- closed

Proposition 1.4a of the frontend removes the entire simple-normal-crossing
locus of the critical discriminant from the flatness-defect scheme:

\[
 Z_{\mathrm{flat}}\subseteq\operatorname{NonSNC}V(\delta_F).
\tag{1.10}
\]

Indeed, after strict henselization, the tame inertia generators around the
SNC components commute.  Degree three makes every nontrivial one a
transposition, hence all are the same transposition.  The normalization is
the disjoint union of the free quadratic Kummer algebra
`R[s]/(s^2-t_1...t_r)` and the free trivial sheet `R`.  The two saturation
tests in (1.9) only need to be run at closed non-SNC points of the
discriminant.  This closes E whenever the branch divisor is SNC and leaves
the non-SNC local models.

Proposition 1.4b closes the ordinary-cusp case as well.  The two cusp
meridians satisfy the \(B_3\) braid relation and act by transpositions on
three sheets.  Equal transpositions give the free
`R[Z]/(Z^2-delta) plus R` Kummer model; distinct transpositions give the
free monic root cover `R[T]/(T^3+uT-v)`.  Hence

\[
 Z_{\mathrm{flat}}\subseteq
 V(\delta_F)\setminus
 \bigl(V(\delta_F)^{\mathrm{SNC}}\cup
       V(\delta_F)^{\mathrm{oc}}\bigr).
\tag{1.11}
\]

Certificate E is therefore closed for a discriminant with only SNC and
ordinary-cusp singularities.  The remaining local target is a
worse-than-ordinary-cusp point.

Proposition 1.8d quantifies that target for every reduced defect.  If \(h\)
is its ternary-cubic multiplication symbol, the degree-six initial branch
equation is

\[
 \mathcal H_h([r])=\operatorname{Disc}(h|_{r^\perp}).
\tag{1.12}
\]

It is nonzero exactly for squarefree \(h\), giving branch multiplicity six;
for non-squarefree \(h\) the multiplicity is at least seven.  Hence the
surviving E attack has two explicit rows:

1. exclude a reduced Koszul defect at a branch point of multiplicity at
   least six;
2. exclude a nonreduced higher determinantal defect.

The present boundary invariant does not record closed-point branch
multiplicity, so this numerical jump is a new obstruction target rather
than a completed contradiction.

The exact symbol-stratified audit in Proposition 1.8e of the frontend tests
the double-saturation modules on all nine nonzero ternary-cubic orbit
representatives and the zero symbol.  For the homogeneous tensor, every
squarefree row has

\[
 H^0_{\mathfrak m}(Q)=0,\qquad
 \operatorname{length}\operatorname{Ext}^2_A(T,A)=6.
\tag{1.13}
\]

This Ext module is killed by \(\mathfrak m^2\), has a three-dimensional
top, and hence has Hilbert function \(3+3t\) in every squarefree row.
The double- and triple-line rows have a one-dimensional support obstruction
and hence fail purity.  The zero homogeneous tensor passes the module tests
but is nowhere generically étale.  Adding one explicit order-four kernel
tensor restores a finite length-six support obstruction in all ten rows
while leaving the cotangent presentation saturated.  This does not prove
the result for every higher lift.  It does show that neither cubic orbit
type nor conormal saturation is the missing discriminator in these leading
models: the live gate is support saturation \(T=T^{[2]}\), together with
normality and the global Keller-open compatibility of the chosen lift.

The seven squarefree rows admit a stronger exact comparison.  Along the
families \(\phi_h+t\psi _4\), computed over
\(\mathbb Q[t,x,y,z]\), the cotangent saturation vanishes uniformly and
the relative \(\operatorname{Ext}^2\) presentation is pulled back from
\(t=0\).  It has no \(t\)-torsion, is supported on the collision axis, and
has multiplicity six.  Thus the support defect is flat along these seven
lines.  An enlarged audit tests all 24 nullspace-basis axes for every
squarefree symbol.  All 168 families retain cotangent saturation,
collision-axis support, no parameter torsion, and relative multiplicity
six.  Four rows change their literal Gröbner presentation, showing that
presentation equality is stronger than, and should not be confused with,
the invariant flatness statement.  The calculation still does not cover
arbitrary linear combinations in the full order-four kernel or arbitrary
higher orders.

For the smooth symbol, all 276 full coordinate planes
\(u\psi_{4,i}+v\psi_{4,j}\) have also been tested over
\(\mathbb Q[u,v]\).  On every plane the relative Ext presentation is
pulled back from the origin, with support equal to the parameter plane and
relative multiplicity six.  This closes every specialization supported on
at most two basis vectors.

The other six squarefree symbols have now been tested on all 1,656
coordinate planes.  In 1,652 rows the pruned presentation is pulled back
from the origin.  Four rows have a nonconstant ambient presentation, but
their exact finite \(\mathbb Q[p_0,p_1]\)-presentations satisfy
\(\operatorname{Fitt}_6=(1)\) and \(\operatorname{Fitt}_5=0\), proving
local freeness of rank six, hence freeness by Quillen--Suslin.  Therefore
all 1,932 squarefree-symbol coordinate planes retain the length-six defect;
any remaining singular-symbol escape uses at least three basis tensors.

Basis sparsity is not itself an escape criterion.  With
\(\psi_+=\sum_i\psi_{4,i}\) and
\(\psi_-=\sum_i(-1)^i\psi_{4,i}\), the exact polynomial plane

\[
 \phi_h+u\psi_++v\psi_-
\]

has saturated cotangent presentation and a relative length-six support
defect for every squarefree symbol.  Its pruned rank-three Ext presentation
is pulled back from \(u=v=0\), and on \(u^2-v^2\ne0\) all 24 fixed kernel
coordinates are nonzero.  This is the first parameter family of
full-support quartic tensors in the audit.  It is one plane, not a
Zariski-open subset of the universal 24-space.

After translation by the deterministic generic quartic lift \(\psi_g\),
the same statement holds for all ten cubic-symbol strata:

\[
 \phi_h+\psi_g+u\psi_++v\psi_- .
\]

Every affine plane has saturated cotangent presentation and a relative
length-six support defect pulled back from its center.  In particular, the
double-line, triple-line, and zero symbols retain the pure finite defect
uniformly after the order-four lift, rather than only at one tested
endpoint.  The generic tensor again has all 24 kernel coordinates nonzero.
Normality and Keller-open compatibility remain separate and are not
implied by this module calculation.

All 2,024 smooth coordinate three-spaces have now also been tested.  Their
pruned rank-three relative Ext presentations are pulled back from the
parameter origin, again with relative multiplicity six and saturated
cotangent presentation.  Two raw resolution presentations vary, but the
variation disappears under pruning.  Consequently the first untested
smooth quartic direction in that coordinate-subspace audit had to use at
least four basis tensors.

The separate universal formal-gauge certificate
[`UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md`](UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md)
now closes those remaining smooth directions.  If \(K\) is the graded
module of compatible tensor corrections and \(G\) is the
determinant-twisted simultaneous collision-coordinate/coefficient-module
gauge differential, the exact identities
\[
 K=\operatorname{im}G+A\eta,\qquad
 G L=[x\eta,\ y\eta,\ z\eta]
\]
give \(K/\operatorname{im}G\simeq\mathbb Q(-3)\).  Hence every compatible
term of collision degree at least four is gauge.  Formal equivalence with
the saturated smooth central fiber and faithful detection of
\((x,y,z)\)-power torsion by completion prove
\[
 H^0_{(x,y,z)}(\Omega)=0
\]
over the full 24-parameter smooth family.  This is the universal
cotangent-saturation theorem, not another finite coordinate sweep.

The same gauge proof has an exact boundary.  For the six singular
squarefree symbols, \(K/\operatorname{im}G_h\) has positive-dimensional
support and nonzero quartic part, so formal reduction to the central fiber
is unavailable.  That failure of the smooth proof is not a remaining
saturation problem: `KDSQ6` treats each complete quartic nongauge
complement, and `SSADPALL` puts every compatible higher tail into its
minimal multi-coefficient formal normal form and proves that the cotangent
presentation stays saturated while the intrinsic different remains
six-generated.  What remains open is normality/algebraization and
compatibility with a distinguished Keller open, together with the separate
boundary-geometric step that would force the different to be Cartier.

For a smooth homogeneous symbol, the obstruction can also be made global.
The corresponding graded algebra is a normal integral threefold.  Its
punctured spectrum is a \(\mathbb G_m\)-torsor over the
\(\mathbb P^1\)-bundle \(Z_h\to C_h\), where \(C_h\) is the smooth plane
cubic.  After algebraic closure, the elliptic
\(\operatorname{Pic}^0(C_h)\) survives in the divisor class group modulo
the grading class.  That class group is not finitely generated, whereas
the class group of a normal variety containing an \(\mathbb A^3\) open is
generated by the finitely many divisorial components of its boundary.
Hence this homogeneous algebra is an explicit globally defined normal
length-six defect model, but it cannot contain the distinguished Keller
open.

The coefficient-six problem is now closed intrinsically.  The fixed
six-row tail is tautological symmetric multiplication on the exceptional
plane.  If \(\mathcal Q\) is its universal rank-two quotient, its projective
charge is

\[
 c_2(\operatorname{Sym}^2\mathcal Q)=6,
\]

while the transpose vertex module is
\(V^\vee\oplus\bigwedge^2V^\vee(1)\), has square-zero maximal-ideal action,
and has length \(3+3=6\).  The binary-cubic discriminant simultaneously has
class \(c_1((\det\mathcal Q)^{\otimes6})=6H\).  Thus the branch
multiplicity and support length are the first- and second-Chern forms of
the same triple-cover quotient bundle, not orbit-by-orbit coincidences.
[Section 6.1 of the universal cubic quartic-kernel note](UNIVERSAL_CUBIC_QUARTIC_KERNEL_SATURATION.md#61-the-intrinsic-localized-six-charge)
contains the proof and the exact scope.  It also records the terminology
warning that this colength/localized-Chern charge is not automatically the
classical Buchsbaum--Rim multiplicity.

<!-- status-consumer: KLC6 88e899a9645c4a70 -->

The remaining continuation gates are:

1. derive from boundary-minimal Keller geometry either the support-hull and
   conormal equalities \(C=T\) and \(Q=T\tau\), or a Cartier-different
   theorem strong enough to invoke the six-generator obstruction; `KDCD3`
   shows that Cartierness closes the support hull only with the still-unproved
   local-CM input;
2. analyze the double-line, triple-line, and zero leading-symbol rows, which
   remain behind the independent generic-etaleness and Keller-compatibility
   gates after `SSADPALL` excludes every squarefree row under a Cartier
   boundary theorem;
3. only after a local model survives, prove normality, algebraize the formal
   gauge, construct the actual \(\mathbb A^3\) Keller open, and then address
   the remaining coefficient/base-change rigidity.

### Certificate P -- no phantom boundary (closed for boundary-minimal cubics)

Let `delta_F` be the reduced branch equation and `j_F` the reduced
nonproperness equation.  The unique cubic critical prime gives

\[
 j_F=\delta_F u_F.
\tag{1.2}
\]

The quotient `u_F` is a unit exactly when there is no second unramified
boundary divisor.  Equivalently, the branch and nonproperness hypersurfaces
agree in codimension one.  It is enough to prove that the nonproperness
hypersurface is irreducible.

Proposition 2.5 of
[`MINIMAL_BOUNDARY_CLASSIFICATION.md`](MINIMAL_BOUNDARY_CLASSIFICATION.md)
now closes this certificate.  The foundational cubic is a degree-three
competitor with exact invariant `mu=(1,1,1)`, while nonproperness and the
nonempty branch divisor give the componentwise lower bounds after each
lexicographic equality.  Boundary minimality therefore forces the same
invariant.  In particular `#Irr(B_F)=1`; the sole boundary prime is critical
of different order one, and its branch image is the unique component of
`B_F`.  Thus `u_F` is a unit by Proposition 2.1a of the cubic normalization
frontend.

This closure has exactly the stated scope.  For an arbitrary cubic,
Proposition 1.4 already proves that an extra boundary prime cannot lie over
the critical discriminant: the ramified `(2,1)` boundary sheet and affine
`(1,1)` sheet exhaust the degree.  The only remaining possibility is a
separate unramified target divisor, precisely a nonconstant factor of
`u_F`.  Boundary minimality rules it out.  The specialization
`A=0, gamma=1` in
[`UNIVERSAL_CUBIC_UNGRADED_TESTBED.md`](UNIVERSAL_CUBIC_UNGRADED_TESTBED.md)
shows that this alternative actually occurs and can force genuine
ungradedness, so no universal reduction to the boundary-minimal stratum is
possible.

### Certificate G -- coefficient-gauge rigidity

For alternating upper/lower times, the homogeneous correction equation in
degree `n` is

\[
 \mathcal L_n(p_n,q_n)=-\mathcal R_n,\qquad
 \mathcal L_n=D_+\oplus D_-.
\tag{1.3}
\]

Its cokernel is

\[
 \operatorname{coker}\mathcal L_n=
 \begin{cases}
 0,&4\nmid n,\\
 k\Delta^{n/4},&4\mid n.
 \end{cases}
\tag{1.4}
\]

Thus the entire two-shear formal problem has one scalar compatibility
condition every four degrees.  Single shears, polynomial Borel gauges,
transported two-shear gauges, all homogeneous-linear cancellations, and all
two-monomial cancellations are already classified or excluded.
Proposition 2.8 adds a decisive warning: the first bilinear interaction
always has zero discriminant projection.  A raw “first scalar is nonzero”
argument cannot work.

This is Propositions 2.2--2.8 of
[`CUBIC_GAUGE_STRAIGHTENING.md`](CUBIC_GAUGE_STRAIGHTENING.md).

### Closure implication

If Certificates E and P vanish and the binary-cubic orbit is reduced to the
tangent-nonosculating slice—either directly or by Certificate G
rigidity—then the canonical normalization gateway applies and the map is
polynomially left--right equivalent to the foundational cubic.

## 2. Attack E: prove two Ext vanishings

This is the highest-priority attack because it simultaneously extracts the
closed-point conormal marking and removes the finite-flatness defect.

### E1. Support-hull attack

Compute the canonical finite map

\[
 T\longrightarrow C
\tag{2.1}
\]

from the scheme-theoretic ramification support to its canonical `S_2`
hull.  Algebraically,

\[
 C=\operatorname{Ext}^1_A(
      \operatorname{Ext}^1_A(T,A),A).
\]

Purity and `S_1` already make this map injective and an isomorphism in
codimension at most one.  The only remaining cokernel is

\[
 L=C/T,
\tag{2.2}
\]

supported at closed collisions.  Proposition 1.15 proves the exact
identity

\[
 \operatorname{Ext}^2_A(T,A)
 \simeq\operatorname{Ext}^3_A(L,A).
\tag{2.3}
\]

The target is no longer an unspecified depth theorem: prove that no
regular function on the normalized critical boundary is missing from the
scheme intersection at a closed collision.

Three concrete ways to kill `L` are now available:

1. **conductor:** show that the conductor of `T subset C` is not supported
   at a closed point;
2. **length:** show that `L nonzero` adds excess length to the corresponding
   cubic collision fiber, contradicting degree three;
3. **intersection:** identify `T` directly with the pushforward of the
   saturated ramified-sheet/affine-sheet intersection.

### E2. Conormal-cokernel attack

Once `L=0`, compute

\[
 K=\operatorname{coker}(T\mathop{\longrightarrow}^{\tau}Q).
\tag{2.4}
\]

It is finite length and its canonical dual is `Ext^3_A(Q,A)`.  The formal
boundary package already retains the nilradical and its powers in every
collision fiber.  Proposition 1.12 of the frontend identifies `K=0`
fiberwise with principal generation of that nilradical.  The target theorem
is therefore:

> The primitive codimension-one conormal class remains a generator after
> every saturated ramified-sheet collision.

This is a Nakayama statement on one marked element, not a classification
of all Artin cubic algebras.

Proposition 1.16 makes a shorter attack available.  Under `L=0`,

\[
 K=H^0_Z(Q).
\]

Thus it is enough to prove that the Jacobian presentation of `Q` has no
finite-length submodule.  Equivalently, every associated prime of `Q` has
dimension two.  This is the exact “absence of point torsion” statement
suggested by saturated scheme intersections.

Proposition 1.17 makes the test algorithmic.  If

\[
 F_1\mathop{\longrightarrow}^{\Psi}F_0\to Q\to0,
\qquad N=\operatorname{im}\Psi,\qquad
 I=\operatorname{Fitt}_3^A(B),
\]

then

\[
 H^0_I(Q)=(N:_{F_0}I^\infty)/N.
\]

Certificate E is therefore exactly:

1. the canonical bidual map `T -> T^[2]` is onto;
2. the presentation submodule `N` is `I`-saturated.

Both are finite syzygy calculations in the intrinsic normalization algebra.

### E3. Perfect-presentation attack

Compute finite presentations of `T` and `Q` from:

1. a presentation of the finite normalization algebra `B`;
2. the Jacobian presentation of `Omega_{B/A}`;
3. the intrinsic primitive class `tau`.

The target statements are the exact projective-dimension bounds

\[
 \operatorname{pd}_A(T)=1,\qquad
 \operatorname{Ext}^3_A(Q,A)=0.
\tag{2.5}
\]

Unlike a direct classification of collision algebras, (2.5) can be checked
by syzygies and saturation.  Proposition 1.15 also permits a smaller
calculation: present only the two finite cokernels `L` and `K`.

Proposition 1.18 gives the clean sufficient endpoint: balanced square
presentations

\[
 0\to A^r\to A^r\to T\to0,\qquad
 0\to A^s\to A^s\to Q\to0
\]

make both modules perfect of grade one and kill Certificate E immediately.
With one critical prime, their determinants are supported on the branch
equation, so the different/conductor attack can seek matrix factorizations
of its powers.

### E4. Boundary-intersection attack

The codimension-one branch theorem already supplies `T tau=Q` away from
closed collisions.  It should also identify `T=C` away from those points.
The scheme intersection of the ramified boundary sheet and affine sheet
should therefore identify both finite modules in (1.5):

- a curvilinear collision has `L=K=0`;
- the reduced Koszul defect contributes a nonzero closed-point deficiency;
- a nonreduced higher defect must contribute to the same two-cokernel
  ledger.

The required new theorem is:

> Saturated boundary monotonicity makes the ramification-support
> presentation perfect of grade one and leaves no finite-length submodule
> in the relative cotangent module.

This statement is exactly (2.5), with no coordinate marking added.

### E5. Different/conductor attack

Compare `T` with the Kähler different and compare `Q` with the inverse
different.  The unique tame `(e,f)=(2,1)` critical prime fixes these modules
in codimension one.  Reflexive saturation now has a named remainder:
the support conductor quotient is `L` and the primitive-generation quotient
is `K`.  Principal generation eliminates the latter by Proposition 1.12 of
the frontend.

Proposition 1.15a now closes the local algebra under the precise additional
hypothesis that was previously hidden.  If the Kähler different is Cartier,
`Ann_B(Q)=dB`, then at a closed collision

\[
 L\simeq (0:_{H^2_{\mathfrak n}(B)}d).
\]

Consequently `L=0` if and only if the normal threefold $B$ is
Cohen--Macaulay there.  Codimension-one principality plus normality/`S_2`
does **not** suffice: the remaining support-hull class is exactly the
Cartier equation's torsion in the degree-two normalization deficiency.
Thus this route has one sharply stated geometric obligation—derive local
CM from the saturated minimal-boundary intersection—or else it must retain
the nonzero deficiency class.

The quartic double-saturation stratification supplies the complementary
attack.  For each of the six singular-squarefree nongauge families, exact
annihilator reduction gives at every geometric collision

\[
 \mu(J_{\mathfrak n})
 =\dim_k J/\mathfrak nJ
 =6.
\]

Hence none of these persistent length-six quartic models has a Cartier
Kähler different.  This creates a genuine fork:

- prove the boundary-minimal support is `S_2`, hence `B` is locally CM under
  a Cartier different;
- or prove that the six-generator initial different cannot become cyclic
  under compatible higher corrections.

The second route is now a filtered Nakayama problem, not another `Ext`
calculation.  Proposition 1.15c of the frontend supplies the lifting lemma:
flatness of the cotangent and annihilator-cokernel packets commutes `J` with
base change, and local freeness of `J/nJ` of rank six preserves the exact
generator count.  Upper semicontinuity alone is insufficient because the
number of generators can drop under deformation.

For the nodal symbol, the filtered Nakayama route is now closed to all
formal orders.  The cyclic gauge cokernel gives the recursive normal form
`h_nod+f(y,z)*eta`.  The universal family `h_nod+u*eta` has
`J/nJ=Q[u]^6`; after assigning `u` collision weight one, the cotangent and
annihilator-control Rees packets are strict with central initial form.
Because the initial form of every graph equation `u-f` is monic in `u`, the
annihilator commutes with every polynomial or formal graph specialization.
Thus the nodal different remains six-generated and non-Cartier after every
compatible higher correction.  Minimal multi-coefficient normal forms and
the same successive monic-graph argument now close the other five
singular-squarefree rows.  The entire higher-order filtered Nakayama queue
is therefore closed.  The formally rigid smooth row has the same
six-generator central different, so no squarefree cubic formal collision
can acquire a Cartier Kähler different from compatible higher terms.  A
boundary-geometric Cartier theorem would force the leading symbol into the
three non-squarefree rows before their generic-étaleness gate.

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

<!-- status-consumer: KDCD3 c6e56b39bbb498d8 -->

<!-- status-consumer: KDSQ6 cd423f625f1f3cd2 -->

### E6. Exact one-boundary class group and the affine Chern no-go

The one-boundary hypothesis determines the global divisor class group more
strongly than mere cyclic generation.  Put

\[
 X=\operatorname{Spec}B,\qquad U\simeq\mathbb A^3,
 \qquad (X\setminus U)_{\mathrm{red}}=D,
\]

where \(X\) is normal and \(D\) is its unique boundary prime.  Restriction of
Weil divisors gives the exact sequence

\[
 B^*\longrightarrow \mathcal O(U)^*
 \xrightarrow{\operatorname{ord}_D}\mathbb Z[D]
 \longrightarrow \operatorname{Cl}(X)
 \longrightarrow \operatorname{Cl}(U)\longrightarrow0.
\tag{2.6}
\]

Indeed, every prime divisor other than \(D\) meets \(U\), so the kernel of
restriction on divisor groups is exactly \(\mathbb Z[D]\); exactness then
follows directly from the principal-divisor presentation of the Weil class
group.  This is also the elementary one-prime specialization of divisor
localization; see [Stacks Project, Section 31.27](https://stacks.math.columbia.edu/tag/0BE0).

Now \(\mathcal O(U)^*=k^*\) and
\(\operatorname{Cl}(U)=0\).  Also \(B^*=k^*\): a unit of \(B\) restricts to
a unit on \(U\), hence to a constant, and every nonzero constant is already
a unit of \(B\).  Thus the first map in (2.6) is an isomorphism and

\[
 \boxed{\operatorname{Cl}(X)\simeq\mathbb Z[D].}
\tag{2.7}
\]

In particular \(D\) is primitive, has infinite order, and is not principal.
More explicitly, if a rational function had divisor \(nD\), its restriction
to \(U\) would be a unit and therefore a constant; hence \(n=0\).  The former
fallback through a primitive rational boundary function is consequently
unavailable.  This does **not** prove that \(D\) is non-Cartier: an affine
normal scheme can have nontrivial Picard group.  It therefore does not by
itself settle the Cartier--Cohen--Macaulay hypothesis of Proposition 1.2 of
the frontend.

Normality identifies Cartier classes with a subgroup of the Weil class
group, so there is a unique \(m\geq0\) such that

\[
 \operatorname{Pic}(X)=m\mathbb Z[D]\subseteq\mathbb Z[D],
\tag{2.8}
\]

where \(m=0\) denotes the zero subgroup.  The boundary \(D\) is Cartier
exactly when \(m=1\).  Thus the global Cartier question has become a single
primitive-generator test, not a class-group-vanishing problem.  At the same
time, this only controls the reflexive divisorial different.  It does not
make the possibly nonsaturated Kähler different
\(\operatorname{Ann}_B\Omega_{B/A}\) principal; the quotients \(C/T\) and
\(H_Z^0(\Omega_{B/A})\) retain precisely that closed-point information.
The rank-one-reflexive interpretation of \(\operatorname{Cl}(X)\) is
[Stacks Project, Lemma 31.30.2](https://stacks.math.columbia.edu/tag/0EBM).

This generator test is exactly the Kähler-different test under the standing
purity hypothesis of Certificate E.  Let

\[
 J=\operatorname{Ann}_B\Omega_{B/A},\qquad
 \mathfrak p_D=I(D),\qquad T=B/J,
\]

and assume \(T\) is \(S_1\) with irreducible support \(D\).  The tame
different order one gives
\(J_{\mathfrak p_D}=\mathfrak p_DB_{\mathfrak p_D}\), while
\(\sqrt J=\mathfrak p_D\).  Thus
\(\mathfrak p_D/J\) is a submodule of \(T\) which vanishes at the generic
point of \(D\).  If it were nonzero, one of its associated primes would be
an associated prime of \(T\).  The \(S_1\) hypothesis says that the only
such prime is \(\mathfrak p_D\), contradicting generic vanishing.  Hence

\[
 \boxed{J=\mathfrak p_D.}
\tag{2.9}
\]

Consequently the following are equivalent:

\[
 D\text{ is Cartier}
 \quad\Longleftrightarrow\quad
 m=1
 \quad\Longleftrightarrow\quad
 J\text{ is locally principal}.
\tag{2.10}
\]

Moreover \(C/T\) is literally the \(S_2\)-hull deficiency of the reduced
boundary surface \(B/\mathfrak p_D\).  At a certified squarefree formal
collision, the six-generator theorem for \(J\) therefore says that the
local image of \([D]\) is nonzero and that \(m\ne1\).  The exact global
squarefree-exclusion target is now

\[
 \boxed{\operatorname{Pic}(X)\longrightarrow
        \operatorname{Cl}(X)=\mathbb Z[D]\text{ is surjective}.}
\tag{2.11}
\]

This remains a boundary-geometric theorem to prove; neither (2.7) nor
boundary minimality alone supplies the surjectivity.

The codimension-one cubic ledger supplies two useful refinements.  Since the
unique ramification prime is \(D\) and its tame different order is one, the
relative dualizing class, equivalently the different-divisor class, is the
generator \([D]\); the reflexive different ideal has class \(-[D]\).  If
\(\delta\) is a prime equation of the target branch and \(E\) is the closure
of the unique unramified affine sheet above it, Proposition 1.4 gives

\[
 \operatorname{div}_X(\pi^*\delta)=2D+E,
 \qquad [E]=-2[D].
\tag{2.12}
\]

Thus even the exact height-one ledger contains compensation rather than a
vanishing boundary class.

There is a parallel obstruction to using the localized six-charge directly
on the affine target.  For a collision point \(p\in Y=\mathbb A^3\), the
forget-support map

\[
 K_0^{\{p\}}(Y)\longrightarrow K_0(Y)
\]

kills \([k(p)]\): the Koszul resolution on three regular parameters gives
\([k(p)]=(1-3+3-1)[\mathcal O_Y]=0\).  Consequently the intrinsic class
\(6[k(p)]\) of KLC6 also maps to zero, regardless of whether the supported
module itself vanishes.  Ordinary affine \(K_0\), ordinary affine Chern
classes, and (2.7) therefore cannot supply the missing contradiction.

The viable global target is now precise: extend the perfect different
complex to a compactification, or retain it in a relative/logarithmic
localization group, and prove a conservation identity in which the boundary
term is computed.  Only such an identity still sees both \(6[p]\) and its
possible compensation along \(D\) and \(E\).  This is a narrowing theorem,
not the missing conservation identity.

<!-- status-consumer: OBCL3 b4b3a3f27b8a97b1 -->

The smooth-symbol exceptional divisor does not make the required Cartier or
compactified conclusion automatic.  If \(h\) is smooth, Proposition 1.8b gives a
\(\mathbb P^1\)-bundle over an elliptic cubic and hence
\(\operatorname{Pic}^0(Z_h)\ne0\).  The one-boundary localization sequence
now makes \(\operatorname{Cl}(\bar X)=\mathbb Z[D]\), while every local
class group \(\operatorname{Cl}(\mathcal O_{\bar X,q})\) is a cyclic quotient
of that boundary class.  But \(\operatorname{Pic}^0(Z_h)\) initially belongs to the
formal Picard group of the completed blowup.  Completion need not preserve
surjectivity on divisor class groups, and formal divisor classes need not
algebraize on the punctured local scheme.  Therefore

\[
 \operatorname{Pic}^0(Z_h)\ne0
 \quad\not\Longrightarrow\quad
 \operatorname{Cl}(\mathcal O_{\bar X,q})\ \text{is noncyclic}
\tag{2.13}
\]

without a separate algebraization theorem.  Any revival of this route must
prove that the elliptic formal classes descend to the actual punctured
neighborhood; cyclicity alone is not a contradiction.  This distinction is
genuine: Brevik--Nollet,
[*Local Picard Groups*](https://arxiv.org/abs/1110.1867), Theorem 1.2,
construct algebraic UFD local rings whose completions are prescribed normal
hypersurface singularities with larger formal class groups.

## 3. Certificate P: closed by divisor minimality

No separate elimination or compactification argument is needed for a
boundary-minimal cubic.

### The irreducibility argument

The first entry of the boundary-minimality invariant is the number of
irreducible target boundary components.  The foundational cubic realizes
the value one, and nonproperness rules out zero.  Hence `S_F` is irreducible.
Since it contains the branch divisor, Proposition 2.1a forces

\[
 S_F=V(\delta_F),\qquad u_F\in k^*.
\tag{3.1}
\]

For cubic maps not assumed boundary-minimal, graph-at-infinity elimination
or the former degree-three sheet argument remain possible tests.  They are
outside the minimality closure path.

## 4. Attack G: finish the coefficient orbit

### G1. Bilinear obstruction attack -- closed negatively

The proposed implication

\[
 \operatorname{proj}_{k\Delta}(\mathcal R_4)\ne0
\tag{4.1}
\]

cannot hold.  Proposition 2.8 proves, in every pair of homogeneous degrees,
that the full bilinear interaction lies in
`D_+(R)+D_-(R)`.  In particular, the degree-four projection is identically
zero for the ten-dimensional quadratic cancellation kernel.  This route
should not receive further effort.

### G2. Cubic-recursion attack

After the bilinear cancellation, the exact remaining Jacobian term is

\[
 q\bigl(D_-(p)H(q)-H(p)D_-(q)-qD_-(p)\bigr).
\tag{4.2}
\]

Compute its invariant projection together with the corrections forced in
lower degrees.  For quadratic leading times, exact basis recursions are
soluble through degree eight, so the next useful calculation must be
generic on the coupled-kernel quotient, not another monomial enumeration.

This attack is now a fallback.  Formal constant-Jacobian solvability is
too close to a general formal volume-preserving problem and does not by
itself give polynomial invertibility.

### G3. Newton-face attack

Apply the lowest-face argument of Proposition 2.6 to every exposed face of
the Newton polytopes of `p` and `q`.  A face containing one monomial is
already impossible.  Hence any unresolved pair must have:

1. equal lowest degrees;
2. at least two monomials on both interacting faces;
3. zero discriminant projection at every degree `4r`.

This is a sharply constrained finite support problem.

### G4. Avoid full gauge classification -- preferred

It is enough to prove that the particular coefficient morphism extracted
from the minimal-boundary package reaches the tangent slice.  The full
group `SL_2(k[Y])` need not be classified.  Intrinsic primitive conormal
generation from Certificate E supplies a preferred root direction; use it
to reduce the Tschirnhausen bundle to a Borel subgroup, where Corollary 2.3
already finishes the gauge.

Proposition 2.9 of the gauge note makes the required output exact.  Extract:

1. a saturated rank-one flag
   \[
    0\to L\to M\to N\to0
   \tag{4.3}
   \]
   in the trace-free rank-two module, with `L,N` projective;
2. primitive generators for which the mixed cubic coefficient is a scalar
   `q in k^*`;
3. the coordinate equality
   \[
    A=k[C_0,C_2,C_3].
   \tag{4.4}
   \]

Quillen--Suslin and projectivity split the flag, the mixed coefficient
becomes `C_1=q`, and all residual basis changes are Borel.  Equation (4.4)
then identifies the target with the full tangent hyperplane and invokes
the cubic hyperplane theorem.

Thus the real E-to-G question is no longer “is the coefficient map
linear?”  It is:

> Does the double-saturated conormal direction descend to a projective
> Tschirnhausen flag whose three complementary coefficients generate the
> target ring?

This formulation matches the localized chart equality and primitive
quotient conditions already appearing in the positive cubic chart theorem.

Proposition 2.10 gives a less assumption-heavy version of the third step.
After the flag and `C_1=q` are extracted, the original Keller map is the
cartesian pullback of the foundational tangent-hyperplane map

\[
 \mu_q:V_q\to H_q
\tag{4.5}
\]

along the remaining coefficient map

\[
 G=(C_0,C_2,C_3):Y\to H_q.
\tag{4.6}
\]

Thus (4.4) is equivalent to proving **foundational base-change rigidity**:

> If `Y times_(H_q) V_q` has simple-root open `A^3` and the pulled-back
> cover has the minimal-boundary package, then `G` is a polynomial
> automorphism.

The first attacks on this narrower statement should be:

1. compare the scheme-theoretic ramified/affine-sheet intersection before
   and after base change; double saturation should exclude nontransverse
   fibers of `G`;
2. compare the pulled-back discriminant `Delta o G` with the phantom-factor
   equation; a divisor of nonproperness of the base change not lying over
   `Delta` must be the forbidden unramified component;
3. use `Y times_(H_q) V_q isomorphic to A^3` to compare class groups and
   units of the two normalized factor spaces.

Geometric degree does not by itself prove `deg G=1`: the upper horizontal
map in the cartesian square has the same degree as `G`.

### G5. Diagonal toric orbit limits -- closed

One direct way to seek a nonflat counterexample is to degenerate the
foundational map inside its polynomial left--right orbit.  In the displayed
foundational coordinates, no diagonal one-parameter subgroup produces a new
Keller map.

Write the foundational map as

\[
\begin{aligned}
 F_1&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
 F_2&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
 F_3&=2x-3x^2y-x^3z.
\end{aligned}
\tag{4.7}
\]

Give `x,y,z` weights `a,b,c`, and let `m_i` be the least monomial
weight in `F_i`.  A diagonal output renormalization with weights `r_i`
has a polynomial limit only when `r_i <= m_i`.  If that limit still has
nonzero constant Jacobian, then

\[
 r_1+r_2+r_3=a+b+c.
\tag{4.8}
\]

The linear terms `z,y,2x` give
`m_i <= (c,b,a)_i`.  Equality in (4.8) therefore forces all six
inequalities to be equalities:

\[
 (r_1,r_2,r_3)=(m_1,m_2,m_3)=(c,b,a).
\tag{4.9}
\]

Put

\[
 p=2b-c,\qquad q=a-b+c.
\tag{4.10}
\]

Every monomial-weight difference from the corresponding linear term is
`alpha*p+beta*q`, where the nonzero pairs `(alpha,beta)` are exactly

\[
 (0,1),(1,0),(1,1),(1,2),(2,1),(2,2),(2,3),(3,2),(3,3).
\tag{4.11}
\]

Since both `(1,0)` and `(0,1)` occur, the admissible cone is exactly
`p>=0, q>=0`, with lineality direction `(a,b,c)=(-1,1,2)`.  Its four
relatively open faces have initial maps

\[
 \begin{array}{c|c}
 p>0,\ q>0&(z,y,2x)\\
 p=0,\ q>0&(z+4y^2,y,2x)\\
 p>0,\ q=0&(z,y+3xz,2x)\\
 p=q=0&F.
 \end{array}
\tag{4.12}
\]

The first three maps are triangular polynomial automorphisms and the last
is the foundational map.  Thus the displayed diagonal toric boundary of
the foundational orbit contains no geometric-degree-three counterexample.
The exact checker is
[`verify_foundational_toric_degenerations.py`](../scripts/verify_foundational_toric_degenerations.py).

This statement does **not** classify constant-conjugate (oblique)
one-parameter subgroups or nonlinear polynomial orbit limits.  Those are
genuinely larger counterexample searches, not consequences of (4.12).

<!-- status-consumer: FTD3 f4b5bf44c04dba69 -->

## 5. Execution order

The shortest credible closure path is:

1. prove E1/E2, obtaining flatness and the preferred primitive root
   direction;
2. use that direction in G4 to reduce coefficient gauge to the proved Borel
   theorem.

If G4 fails, run G2 or G3 as the exact fallback.

No broader classification theorem is needed on this path.
