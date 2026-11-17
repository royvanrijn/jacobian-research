# Global Haar-orbit splitting and the boundary radical

16 November 2026. Internal written proofs with explicitly audited classical
inputs. Independent replay, formal verification and external review are not
claimed. Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**;
review stays deferred.

There is an exact GLOBAL replacement for the two local obstructions. Under
eventual purity, safety holds iff the actual Haar-monodromy orbit contains
a globally fixed ambient lift. Equivalently the rational unit splits off
the actual Cauchy FUNCTION connection, before quotienting by polynomials.
Global semisimplicity of this connection is therefore sufficient, even
when individual local monodromies have nontrivial Jordan blocks.

Every failure must use a genuine mixed-weight GLOBAL monodromy extension
from boundary weight3 or4 of the fibre. Pure compact-middle cohomology,
of arbitrary rank and with arbitrary allowed local unipotent monodromy,
cannot by itself cause a failure. In particular if generic fibre H2 is
pure of weight2, EVERY eventually pure polynomial composite is safe.

See the [checker](verify_haar_global_radical.py),
[receipt](haar_global_radical_verification.json) and
[source audit](haar_global_radical_source_audit.json).

## 1. Keep actual Cauchy functions before the polynomial quotient

Let X=SL2/Gamma, ANY finite Gamma, and f any nonconstant regular function.
Normalize the algebraic volume Omega by its nonzero period Z on the compact
core. On the FULL regular base S in the affine parameter line put

    U_s=X minus {f=s},   B=H3(U_s;C),
    C_Q(s)=Z^-1 integral_gamma Q Omega/(s-f).

H3 denotes cohomology. Continue the physical Haar lift gamma from infinity,
and let W be the complex span of its ENTIRE monodromy orbit in H_3(U_s).
The ambient map q:W -> C is normalized by q(gamma)=1. All finite
bifurcations caused by infinity are retained. Loops are global loops in S.

Define a finite-rank differential module of ACTUAL FUNCTION germs

    M_f=sum_Q C(s)*C_Q(s),                               (1)

BEFORE taking any quotient by C[s]. The single-pole spanning theorem in
the [quadratic-log note](CUBIC_NILPOTENCE_AND_QUADRATIC_HAAR_LOGARITHMS.md)
identifies this with a quotient of the FULL complement de Rham connection.
It is closed under differentiation. It always contains the literal
rational function1, since

    s*C_1-C_f=1.                                        (2)

In fact its local system is canonically W^vee. To check the kernel,
represent a section by a rational cohomology section alpha. If its gamma
period is identically zero, continuation around every loop makes its
period zero on every translate of gamma, since alpha is single-valued.
It therefore annihilates W at every regular point. The converse is
immediate. The analytic annihilator is an algebraic subconnection:
regular-singular algebraic connections on an algebraic curve retain their
monodromy subobjects. Alternatively, its horizontal maps have moderate
growth at all P1 punctures, hence are meromorphic there and algebraic.
Thus the kernel is exactly Ann(W), not an unrelated generic-fibre kernel.

The line C(s)*1 in (1) corresponds to the ambient line L=Q(-2) in B
and is dual to q. Do NOT first pass to C_Q modulo C[s]: that would kill
the function1 and localize away the point module of the pure unit.

## 2. A global invariant lift gives a rational retraction

Suppose that some v in W is fixed by EVERY global monodromy and q(v)=1.
Define R on (1) by replacing gamma with v in every period representation.
This is well-defined precisely because v belongs to W: every relation
between actual gamma-periods annihilates its full orbit. It is C(s)-linear,
commutes with differentiation, and sends the literal function1 to1.

The resulting periods are single-valued on S. Algebraic regular-singular
periods have only power/log growth at ALL punctures of P1, including
infinity. Single-valuedness leaves a meromorphic germ with finite pole
order. Thus every value of R is a GLOBAL rational function. Consequently

    R:M_f -> C(s),   R|C(s)=identity.                    (3)

Conversely, a differential retraction (3) gives an invariant dual vector
v in W with q(v)=1. This is the ordinary flat duality of this actual
connection. No Hodge property of the actual period kernel is assumed.

Assume now EVENTUAL PURITY: the original unit exponential period is a
nonzero polynomial, with constant term1. The exact actual entire image
V_f contains C[t] and has dim(V_f/tV_f)=1, hence is indecomposable.
These are the retained [intrinsic state facts](INTRINSIC_HAAR_STATES_AND_NONSPLIT_PERIOD_EXTENSIONS.md).

Under the coefficientwise Borel/Fourier identification,

    B(sum a_n t^n/n!)=[sum a_n s^(-n-1)] modulo C[s],
    t=-partial_s,   partial_t=s,

the polynomial submodule C[t] becomes the point module
D_0=span{s^-1,s^-2,...}. Since R fixes rational functions and kills no
nonzero principal part, composing it with PP_0 gives

    B(V_f) -> D_0,   [g] -> PP_0(R(g)).                  (4)

It is well-defined modulo polynomials. It commutes with differentiation
and multiplication by s modulo polynomials, and is identity on D_0.
Thus (4) splits C[t] off V_f. Indecomposability forces V_f=C[t], which
is exactly full ALL-original-marker safety. The principal-part argument
is the same exact one as in the [algebraic trace theorem](ALGEBRAIC_RESOLVENT_SPLITTING_AND_ZERO_VALUE_LOGARITHMS.md);
the new ingredient is the global invariant-lift retraction instead of
algebraicity of every period.

Conversely safety makes every C_Q a finite Laurent polynomial with its
only finite pole at0. Hence every loop fixes its period. The original
simple-pole classes span full complement de Rham, so every loop fixes
gamma itself. Therefore W=C*gamma and gamma is the required v.

**Global splitting theorem.** For EVERY eventually pure phase on X,

    ALL original markers safe
      iff C(s)*1 splits off M_f as a differential module
      iff W contains a GLOBAL invariant v with q(v)=1
      iff gamma itself is globally invariant
      iff dim W=1.                                     (5)

A globally fixed lift in the FULL homology extension is insufficient
unless it lies in W. Likewise separate fixed lifts for individual loops
need not be one simultaneous fixed lift. These distinctions preserve
the residual actual beta from the local theorem. Once the GLOBAL lift
exists inside W, (5) proves that purity kills that residual beta too.

If global monodromy on M_f (equivalently W) is SEMISIMPLE, its trivial
unit submodule has a complement, so (5) proves safety. In particular
global semisimplicity of the FULL complement B suffices, because M_f is
a quotient. This is global representation semisimplicity, not the
semisimplicity of each local matrix. A reductive global representation
may have nontrivial local unipotent matrices. On the eventually pure
locus, actual global semisimplicity is also necessary: safety makes
M_f=C(s) and W one-dimensional trivial.

## 3. Global complement semisimplicity reflects ALL polynomial recomposition

For ANY nonconstant polynomial H of degree d, choose the parameter open
S' so that H:T -> S' is a finite unramified cover and all primitive and
composite fibre/complement systems are regular. Removing extra regular
values does not change their monodromy image. Pullback from the d
primitive complements to the complement of their union gives an exact
sequence of GLOBAL local systems

    0 -> L tensor ker(pi_* C -> C)
      -> pi_* B_f -> B_(H(f)) -> 0,                     (6a)

where pi:T -> S' is that cover and the map on pi_* C is SUM. Indeed
all d ambient classes restrict to the SAME ambient class, so their
sum-zero relations form the displayed kernel. On residues the map is
the isomorphism from the direct sum of every primitive branch H2(F)(-1)
to the composite fibre cohomology. Exactness of the complement sequences
then proves (6a). All branches and disconnected fibres are retained.

Induction and restriction across a finite-index subgroup preserve
semisimplicity for finite-dimensional characteristic-zero representations.
For restriction, take a normal finite-index core; the sum of the finitely
many conjugates of a simple core submodule is semisimple and invariant,
so gives all of a simple representation. For the converse average a
core-equivariant projection over the finite quotient. Induction restricts
to a finite direct sum of such conjugate core representations, so the
same averaging applies.

Thus if B_f has semisimple global monodromy, so does pi_* B_f and its
quotient B_(H(f)). Conversely, restrict a semisimple B_(H(f)) to T.
Pullback from the distinguished primitive complement injects B_f into
that restriction: it is injective on its ambient line and on its
individual residue summand. A subrepresentation of a semisimple
representation is semisimple. Hence, BEFORE PURITY,

    B_f globally semisimple iff B_(H(f)) globally semisimple. (6b)

It follows that global semisimplicity for ONE fixed primitive f proves
ALL-original-marker safety of EVERY eventually pure H(f), in EVERY
degree. No purity of f or f-a is assumed. Full source semisimplicity
is sufficient rather than necessary for safety; the actual Haar period
can kill nonsplit source components.

## 4. All global nonsemisimplicity comes from higher boundary weights

The complement mixed Hodge system has original weights4,5,6:

    0 -> L=Q(-2) -> B -> K=H2(F_s)(-1) -> 0.            (6)

Every GrW^k B is a geometric polarizable pure variation. On the smooth
algebraic base S, its global monodromy representation is semisimple by
Deligne's Hodge II Theorem4.2.6, with Proposition4.2.5 and its geometric
fixed-part hypotheses. The proof and hypotheses are read in full below.
This theorem is NOT applied to an arbitrary analytic punctured-disc VHS.

Let A be the finite-dimensional complex algebra generated by global
monodromy on a fibre of B, and J its Jacobson radical. Since J kills
every semisimple module, it kills each GrW^k B. Consequently

    J W4 B=0,   J W5 B subset W4 B,   J B subset W5 B,
    J^3=0.                                            (7)

Faithfulness of A on B turns the last vanishing on B into an algebra
identity. More exactly J is the ideal of A acting trivially on all three
weight grades. That ideal has cube zero and hence kills every simple
module; conversely J kills the semisimple grades. This also directly
identifies the radical without confusing it with a coefficient-ring
purity ideal.

Let M be a fibre of the actual quotient (1) and u its nonzero rational
unit vector. It is an A-module and u is a trivial submodule. Elementary
finite-dimensional radical theory gives

    u has an A-linear retraction to C
      iff u NOT in J M.                                (8)

Indeed all maps to the trivial simple module kill J M. Conversely M/JM
is semisimple, so if the image of u is nonzero, its trivial line has a
retraction, which pulls back to M. In combination with (5),

    eventually pure full failure iff u lies in J M.      (9)

Here J^3 M=0 for EVERY phase. Formula (9) locates the exact actual unit
obstruction in a GLOBAL boundary-lowering radical, not merely in an
arbitrary local nilpotent block. J kills the entire pure compact-middle
weight4 system. A failure therefore needs higher boundary weight5 or6
and an actual nonsplit mixed-weight GLOBAL monodromy extension. Nonzero
boundary ranks alone do not prove that extension or the failure.

The radical J is not generated merely by the logarithms of local
unipotent matrices: a globally semisimple representation can contain
such matrices. Nor is J the ordinary coefficient-state ideal of earlier
all-ring results. No arbitrary-ring topological promotion is made.

## 5. Complete pure-middle closure, including ALL polynomial composites

Suppose that the GENERIC smooth fibre of a fixed f has

    H2(F_s;Q) PURE of weight2.                           (10)

No bound on its rank is assumed. Then (6) is entirely pure of weight4.
It is a geometric polarizable pure variation, so its global monodromy
is semisimple. By (5), EVERY eventually pure f is safe. No local fibre
semisimplicity, vanishing of H2, tameness, isolated critical points or
properness of f is required.

For ANY nonconstant polynomial H the generic fibre of H(f) is the
disjoint union of all primitive branch fibres. Its H2 remains pure of
weight2, because it is a direct sum of those pure groups. The same
argument therefore proves

    If (10), EVERY eventually pure H(f) is safe
    for ALL original markers, in EVERY H degree.         (11)

More generally the conclusion follows whenever the FULL complement
weight filtration splits as GLOBAL local systems, because each pure
grade is semisimple. That splitting is NOT asserted universally.

In a smooth proper SNC compactification of the generic fibre, the
weight2 piece is the compact-middle image; weight3 comes from boundary
curve H1 and weight4 from the boundary graph. If every boundary curve
is rational and the graph is a forest, both higher pieces vanish, so
(10)--(11) apply WITHOUT the earlier local semisimplicity hypothesis.
The intrinsic condition is actual weight purity, not just a count of
visible curves in an unresolved compactification.

Thus nonzero generic H2 alone is no longer the cohomological obstruction:
every full exceptional failure needs NONZERO GrW3 or GrW4 of fibre H2,
and the corresponding mixed-weight global monodromy must fail to split
in a way putting the actual unit in (9). Arbitrary compact-middle
degeneration alone, even size-three local blocks, cannot cause failure.

## 6. Scope of the remaining whole-quotient problem

The complete exceptional target is now (5), or equivalently exclusion of
u in JM in (9), for EVERY actual eventually pure regular phase. Neither
global source semisimplicity nor the pure-middle hypothesis is proved
for all remaining phases. Actual pure SAFE exceptional biform phases
already have large higher-weight boundary spaces; those spaces cannot
be universally eliminated. Their coupling to the actual unit is what
must vanish, or a genuine original-phase failure must realize it.

All finite bifurcations, the actual Haar orbit, all original markers,
all primitive branches and every recomposition degree remain. The local
intrinsic xi and residual beta tests, the N^3/log-degree bounds and the
canonical original-form lattice defects remain valid. The new global
retraction closes a class of arbitrary-rank compact-middle systems and
removes the need to test beta separately IF its global orbit criterion
is established. It does not prove universal MZ or construct a failure.

Classical inputs: Deligne Hodge II, complete PDF39--43/printed42--46,
especially the hypotheses4.2.2, Proposition4.2.5 and full Theorem4.2.6
proof; MIT18.706 Spring2023 Lecture5, complete PDF2--3, Proposition5.6
and Corollary5.7 with proofs; retained geometric regularity, exact
complement MHS, original de Rham detection and actual indecomposability.
The orbit retraction, its Borel/principal-part application, finite-cover
semisimplicity reflection, exact global-radical failure criterion and all-composite pure-middle closure
are the written deductions above. Finite controls are not universal
machine verification or actual exceptional phase realizations.
