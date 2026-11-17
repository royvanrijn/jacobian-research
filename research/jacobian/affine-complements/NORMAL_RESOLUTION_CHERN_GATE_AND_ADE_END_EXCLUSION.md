# The resolved Chern class and the complete normal ADE end exclusion

PF2RULED109, 16 November 2026. This removes the SMOOTH-N hypothesis
from the complete ADE-end exclusion in
[PF108](ALL_ADE_END_CANONICAL_ENERGY_EXCLUSION.md). Every actual normal
finite Keller normalization is included, with its WHOLE marked source
U=A2 and ALL omitted primes retained. No field-degree, index, residue,
prime-count or finite-image-singularity bound is imposed.

**Normal-source Chern theorem.** Let F:A2->A2 be a nonautomorphic complex
Keller map, pi:N->A2 its original finite normalization, and
rho:Y->N the MINIMAL resolution. Write D_i for the strict transforms
of its r omitted primes and E_j for the m exceptional components.
A compact core of Y admits a small complex deformation to a simply
connected Stein filling V of the same oriented end, with

    c1(TV) != 0 in H^2(V;Q),    b2(V)>=r+m.                (1)

The complex structure changes; neither the finite target map nor the
marked plane is asserted to remain holomorphic after deformation.
Only topology, the homotopy class of the complex tangent bundle and
the boundary contact isotopy class are carried across.

**Complete normal ADE end exclusion.** The oriented end of N is NEVER
the negative-Cartan link A_l (any l>=1), D_l (any l>=4), E6, E7 or E8.
Equivalently it is never the oriented boundary of any connected
negative-definite all-minus-two sphere-plumbing TREE. These include
SINGULAR actual normalizations. Orientation is the outward complex
core orientation; arbitrary spherical space forms and reversed
orientations are not being excluded.

More generally (1) excludes an oriented end for which every simply
connected Stein filling of every positive tight contact structure has
rationally zero first Chern class. This is a sufficient comparison
criterion, not an assertion that arbitrary ends have that property.

## 1. Boundary divisor classes remain independent on the resolution

The following elementary lemma keeps the WHOLE plane hypothesis visible.
For ANY smooth algebraic surface Y with dense open U=A2 and finitely
many irreducible boundary divisors B_1,...,B_s, localization gives

    Pic(Y) = direct-sum_(i=1)^s Z*[B_i].                 (2)

Indeed Pic(U)=Cl(U)=0, so subtracting a rational principal divisor
makes any divisor supported on the boundary. A principal boundary
combination restricts to an invertible regular function on U. Since
O(U)^*=C*, that function is constant in the whole rational function
field; every coefficient of its divisor is zero. Thus there are NO
relations in (2). Intersections among boundary components and compact
exceptional curves cause no change to this argument. No extra interior
blowups to make the boundary a good resolution are required.

There is also an injective rational Chern-class map

    Pic(Y) tensor Q -> H^2(Y;Q).                         (3)

Here is a proof that does NOT assume rational-singularity vanishing.
Complete Y to a smooth projective rational surface S and arrange the
boundary T=S minus Y to be simple normal crossing, by blowups outside
Y only. On the rational S, the exponential sequence gives
Pic(S) isomorphic to H^2(S;Z), because H^1(O_S)=H^2(O_S)=0.
Every algebraic line bundle on Y extends to S by extending a divisor.
The kernel of H^2(S;Q)->H^2(Y;Q) is the span of the classes of the
irreducible T components. This follows from the pair sequence and the
degree-two Thom classes; crossing points have real codimension four
and add no degree-two generator.

Consequently a line bundle with rationally zero restricted Chern class
has an extension whose positive multiple is a divisor supported on T.
That multiple becomes trivial on Y. Since (2) is free, the original
line bundle is trivial. This proves (3). In particular a nonzero
boundary divisor class in (2) cannot become a topologically trivial
canonical class. This proof uses rationality of S, factoriality of U
and its constant units separately; it is not a general statement for
arbitrary open surfaces.

## 2. The pulled-back Jacobian gives a nonzero canonical class

For the actual minimal resolution, Y is rational because its function
field is C(x,y). The map h=pi*rho:Y->A2 is holomorphic and algebraic.
Pull back the target's standard volume form:

    tau=h*(du wedge dv),
    div(tau)=sum_i(e_i-1)D_i + sum_j a_j E_j.            (4)

It is a NONZERO regular section of K_Y. On the whole marked U it is
c*dx wedge dy, so its zero divisor is supported on precisely the
boundary components in (2). The coefficients on D_i are the original
finite transverse indices minus one. Every a_j is an integer >=1:
h contracts E_j to a point, so its differential has zero determinant
along the generic smooth point of E_j. The a_j are canonical orders
of a contracted curve; they are NOT finite ramification indices of
new omitted primes.

If there is an exceptional curve, (4) is nonzero immediately. If there
is none, the retained smooth argument in PF103 applies: all e_i=1
would make the finite pi etale and hence an isomorphism over A2.
An omitted divisor on that A2 would define a nonconstant polynomial
which restricts to a unit on U, impossible. Isolated omitted points
are ruled out by normal Hartogs extension and affineness of U, as in
PF103. Hence a nonautomorphism still has a nonzero coefficient.

Thus (4) is a nonzero effective boundary vector. By (2) and (3),

    c1(K_Y) != 0 in H^2(Y;Q).                            (5)

This is proved without clearing Laurent keys, extending a source
symmetry to the target field, or computing cyclic discrepancy chains.
It also agrees with the all-degree divisor localization input of
[PF57](CANONICAL_FORMATS_AND_BINARY_CONDUCTOR.md), now keeping the
resolution's exceptional divisors and the topological Chern class.

## 3. Use a minimal complex resolution, then deform it to Stein

The original radial end is retained from PF98. Resolution changes
nothing outside a compact set. A large core X in Y retracts Y: the
original finite radial function has no critical points sufficiently
far out, and all resolution curves are contained inside. A small
affine-embedding-norm perturbation near its boundary makes that
boundary strictly pseudoconvex, preserving its oriented isotopy type.
The function is not required to be strictly plurisubharmonic on the
compact exceptional curves. It could not be so there.

This core is holomorphically MINIMAL. Every compact complex curve in
Y maps to a point in the affine N, because all affine coordinate
functions are constant on a compact irreducible curve. Thus it is
exceptional for rho. Minimal resolution has no exceptional smooth
rational curve of square -1. This excludes every holomorphic (-1)
curve in the core, not merely ones in a selected boundary chain.

Now import Bogomolov--de Oliveira Theorem(2'): a minimal relatively
compact strictly pseudoconvex smooth surface has a small Stein
deformation. Proposition(2) supplies C-infinity trivial deformation,
and Theorem(5)'s proof keeps the boundary contact structure in its
isotopy class by a small deformation and Gray's theorem. Epstein's
2012 addendum applies exactly this theorem to MINIMAL RESOLUTIONS OF
NORMAL STEIN SPACES, confirming that the input is not restricted to
an already smooth affine source or to a local single-curve resolution.
These deformation theorems are imported, not independently reproved.

After identifying the underlying compact manifolds, the complex
tangent bundles form a continuous family. Therefore c1(TV)=c1(TX)
under this identification. Equation(5) survives, because X retracts Y.
The new Stein structure has an exact symplectic form, hence is
symplectically MINIMAL. This is the missing bridge: the ORIGINAL
resolved Kahler form was NOT assumed exact or symplectically minimal.

General position moves every loop in Y off the complex boundary
divisors and their finite singular strata, into the simply connected
U=A2. Thus pi1(Y)=pi1(X)=pi1(V)=0. Equations(2)--(3) also give
b2(V)>=s=r+m. This proves (1).

Minimality is essential here. A nonminimal resolution obtained by
blowing up a smooth point has a holomorphic exceptional (-1) curve
and is not an input to Theorem(2'). We neither apply it to that
blowup nor assert that an arbitrary singularity has a smoothing whose
Milnor fibre is the minimal resolution. A Stein deformation of the
smooth resolution is a different assertion.

## 4. Close the whole normal ADE end class

Suppose the actual end is an oriented negative ADE link of rank L.
Use exactly the primary contact and filling inputs audited in PF108.
Honda covers A type; spherical L-space plus Ghiggini's all-minus-two
count covers all D/E types. There is ONE positive tight contact
structure, the canonical one. It is therefore also the contact
structure on the Stein filling V just constructed.

Exactness makes V minimal. Ohta--Ono's theorem identifies its oriented
smooth topology with that of the ADE Milnor fibre. In particular

    b1(V)=0, b2(V)=L, chi(V)=1+L, sigma(V)=-L,
    the intersection form is NEGATIVE DEFINITE.          (6)

The reference canonical contact structure has theta=L-2, from the
zero-Chern Milnor fibre. For V itself, Gompf's formula gives

    theta=c1(TV)^2-2(1+L)-3(-L)=c1(TV)^2+L-2.

Contact identification forces c1(TV)^2=0. The negative-definite
inverse intersection form on H^2(V;Q), whose boundary is a rational
homology sphere, then forces c1(TV)=0 in rational cohomology. This
contradicts (1). ALL actual NORMAL ADE end packets are excluded,
including all singular cyclic quotient sources, arbitrary residue
covers, every degree and every number of omitted primes.

The unbounded graph classification from PF108 makes this precisely
the negative-definite all-minus-two TREE plumbing end class. No end
is identified merely by its homology order or its intersection matrix;
the oriented boundary identification remains an explicit premise.

## 5. Genuine finite countercontrol: why the plane and its units matter

The q=2 homogeneous model from PF107 descends through the scalar
order-two action to a finite degree12 map from N=C2/{+/-1}:

    P=Y^4+X^3*Y,
    Q=Y^6+(3/2)*X^3*Y^3+(3/8)*X^6,
    Jac_upstairs=-(9/4)*X^8.

Its end really is the A1 link RP3 and its transverse index is9.
Remove the critical boundary X=0. The quotient open is C* times C,
NOT A2. In the minimal-resolution chart v=X/Y, t=Y^2,

    P=t^2*(1+v^3),
    Q=t^3*(1+(3/2)*v^3+(3/8)*v^6),
    Jac_(v,t)=-(9/8)*v^8*t^4.

Here D={v=0}, E={t=0}. The apparent canonical vector (8,4) is FOUR
times the principal relation (2,1), because X^2=v^2*t is a unit on
the chosen open. Thus K_Y is trivial for this A1 resolution, exactly
as it should be. Equations(2)--(5) cannot be imported when the open
has nonconstant units. The finite cover and its ADE end exist; the
theorem excludes their realization with the WHOLE marked plane,
not the finite cover itself. This control also prevents promoting
positive Jacobian coefficients alone to a nonzero Chern class.

The unramified affine quadric from PF103 has K=0 and is still allowed
as a source-only extension, with no compatible finite Keller target.
The blowup controls from PF108 still show why minimality cannot be
inferred from theta and homology alone. Neither source is an actual
counterexample to the theorem.

## 6. Scope, provenance and the remaining goal

The new primary import is Bogomolov--de Oliveira's Stein deformation
theorem, with C-infinity triviality and boundary-contact preservation.
Their author PDF has a draft header Volume00,1998; the publication is
Contemporary Mathematics207(1997),25-41, DOI10.1090/conm/207/02717,
as pinned in Epstein's bibliography. Do not silently relabel draft
pagination as published pagination. The
[manifest](sources/normal-resolution-chern/manifest.json) records full
pages actually rendered and read, the independently written Picard/
Chern and Jacobian arguments, and the retained PF108 classification
audits. No independent whole deformation/classification proof, formal
verification or external review is claimed.

Other oriented ends, common-leading-direction residue germs and the
ACTUAL Laurent finite-algebra finish line remain OPEN. This closes
a complete normal higher-degree end class; it does not prove B=A
for every Keller map or complete the full Wright goal. Published
Domrina2000 quartic retirement and closure before review expense
remain intact. Prior sealed proofs and finite covers are preserved.

See the [checker](verify_normal_resolution_chern.py) and
[receipt](normal_resolution_chern_receipt.json) for exact controls.
