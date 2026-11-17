# Full intrinsic local Tate criterion

16 November 2026. Internal written proof. Independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**; review deferred.

For EVERY eventually pure phase on SL2/Gamma, ANY finite Gamma, full
all-original-marker safety is equivalent to vanishing of the intrinsic
LOCAL unit obstruction at zero in the canonical geometric BOUNDARY source.
There is NO global semisimplicity hypothesis on its reduced fibre quotient.
The possible hidden genus--Tate branch of the preceding
[graph theorem](TATE_GRAPH_NECESSITY_AND_SEMISIMPLE_BOUNDARY_QUOTIENTS.md)
is now closed. The remaining obstruction is an intrinsic rational Tate
class, not an unknown actual Haar-lift or actual relation kernel.

See the [checker](verify_haar_local_source.py),
[receipt](haar_local_source_verification.json) and
[source audit](haar_local_source_source_audit.json).

## 1. Canonically remove the lifted nontrivial genus socle

Use the full regular algebraic base S and canonical admissible source

    0 -> L=Q(-2) -> E -> V -> 0,
    A=W5 E, P=GrW5 E, T=GrW6 E.                        (1)

W4 E=L, P is pure polarizable weight5, and T has finite global monodromy.
The nearby and limiting Hodge types have p,q between2 and3. In particular,
at ANY P1 puncture, after unipotent base change,

    M6 E=E, M3 E=0, N(E) subset M4 E,
    N acts as zero on the original GrW6 E.             (2)

These are the established geometric level-one bounds. They do NOT assert
M=W when genus degenerates. All quotients below are admissible VMHS only
when explicitly constructed as Hodge quotients.

Let P0 be the pure Hodge complement of P^global in P, obtained from a
polarization and pure fixed part. It has NO trivial subquotient. Form

    D=Hom(P0,A), H0=H0(S,D).                           (3)

D is an admissible VMHS with original weights -1 and0. Mixed fixed part
makes H0 a CONSTANT MHS. Its weight -1 part is

    W-1 H0=Hom_local_systems(P0,L)=0.                   (4)

Hence H0 is PURE of weight0. Evaluation is a VMHS morphism, so

    C=image(P0 tensor H0 -> A)                         (5)

is a PURE weight5 subvariation. As an underlying local system it is
exactly the sum of ALL nontrivial simple subsystems of A: any such
subsystem injects into P, hence is reached by a flat map from P0; the
converse follows from semisimplicity and absence of trivial constituents.
Thus C is the canonical nontrivial socle, even when L does NOT split
globally in A. No source-unit splitting has been assumed.

Put Ebar=E/C and Abar=W5 Ebar=A/C. The local system Abar has NO nontrivial
simple submodule. To verify this last assertion rather than assume it,
the injection C -> P and global semisimplicity of P give a projection
A -> P -> image(C) -> C killing L. Hence A=C direct_sum A' as LOCAL
systems. A nontrivial simple submodule of A/C would lift into A', contrary
to maximality of (5). The projection need not itself be Hodge.

This Hom/fixed-part construction strengthens the preceding constant-genus
argument: it removes the LIFTED nontrivial socle in every source, including
globally nonsplit sources. It does not remove a genus quotient for which
no global lift exists.

## 2. Removing C preserves the intrinsic LOCAL unit test

At a fixed puncture let N denote the logarithm after unipotent base
change. For pure C of original weight5, with types(3,2),(2,3), its
monodromy filtration satisfies

    M4 C=image N_C.                                   (6)

Indeed N_C^2=0 by its Hodge strip; the only blocks have size1 centred
at5 or size2 with grades4 and6. All inclusions and quotient maps here
are strict for the actual LIMIT MHS, by the retained geometric limit
argument, not an arbitrary relatively filtered exactness assertion.

Suppose u is in image N on Ebar. Lift a preimage to x in E. Then
Nx-u belongs to C. Both Nx and u lie in M4 E by (2), so strictness puts
their difference in M4 C=image N_C. Subtract its preimage in C from x
to get u in image N_E. The reverse implication is immediate. Therefore

    u NOTin image N_E iff u NOTin image N_Ebar.         (7)

This is exactly reflection of LOCAL unit splitting. Averaging over the
finite local cover restores the original local monodromy and its
eigenvalue1 unit block. No assertion about the full unreduced compact-
middle source B is made: a pure weight4 quotient does not satisfy (6).

## 3. Every global unit-retaining quotient reflects the local unit

Let q:E -> Q be ANY global local-system quotient retaining u. It need
not be Hodge and can be an actual complex period quotient. The image
q(C) is semisimple with no trivial constituent; it cannot contain u.
Consequently

    Qbar=Q/q(C), Ebar -> Qbar                          (8)

still retain u. Put K=ker(Ebar -> Qbar). Its intersection with Abar
injects into the semisimple pure quotient Abar/L=P/image(C), so it is
globally semisimple. By Section1 it has no nontrivial simple constituent.
Thus

    K intersect Abar subset F=Abar^global.             (9)

Mixed fixed part makes F a CONSTANT mixed Hodge subvariation, and W4 F=L.
At the chosen puncture its limiting filtration is its original W. Strict
limit inclusion into Ebar therefore gives

    F intersect M4 Ebar=M4 F=W4 F=L.                  (10)

Suppose now u lies in image N_Qbar. Lift its preimage to x in Ebar, so
Nx-u=k in K. Since N is zero on GrW6 Ebar, k belongs to Abar. By (2),
it also belongs to M4 Ebar. Equations(9),(10) force k into L. But K
is disjoint from L, so k=0. Hence u already lies in image N_Ebar.

It follows that if E has a local unit retraction, so does Qbar, by (7),
and hence so does Q by composition with (8). The converse is composition
with q. We have proved, BEFORE purity, for EVERY puncture:

    unit splits LOCALLY in E iff it splits LOCALLY
    in ANY GLOBAL quotient retaining that unit.        (11)

The word GLOBAL in (11) is essential: the kernel must be invariant under
the entire algebraic-base monodromy, so (3),(9) apply. An arbitrary local
quotient can destroy a source local splitting, as the admissible local
control in the graph theorem shows. That local countercontrol is not a
global algebraic variation with a global quotient of the required type.
Equation(11) does NOT assert reflection of GLOBAL unit splitting for every
abstract quotient, or finite local/global monodromy of those quotients.

## 4. Full safety equals the intrinsic boundary obstruction at zero

The full actual scalar Cauchy FUNCTION module M surjects its actual
boundary quotient Q, a global quotient of E retaining u=1. Under (11)
a local source-unit retraction at zero passes to Q. Regular singularity
gives Q -> meromorphic germs fixing u. Compose M -> Q with this map and
PP0 modulo polynomials. It kills actual relations and C[s], and is
identity on delta0. Under EVENTUAL PURITY the actual Borel image contains
delta0; the retraction and actual entire indecomposability give safety
of EVERY original marker. This is the retained actual local-functional
argument; Q need not itself be realized as scalar functions or be Hodge.

Conversely full safety under eventual purity gives a GLOBAL source-unit
retraction by the intrinsic safety theorem, hence a local one. Define
xi_boundary,a as in the graph theorem, the coinvariant obstruction of the
dual extension of (1). We have the FULL all-phase equivalence

    under EVENTUAL PURITY,
    ALL original markers safe iff xi_boundary,0=0.     (12)

There is NO global semisimplicity assumption on V, no fibre eigenvalue1
semisimplicity, no zero-genus condition and no actual relation-kernel
calculation. The hidden-xi-zero branch is now EMPTY for genuine eventual-
pure phases. The prior source-global/constant-TYPE criterion remains
equivalent; under eventual purity (12) forces global constancy from this
ONE local condition, rather than from a blanket local-to-global theorem.

By the exact geometric Tate-factor theorem, xi_boundary,0 is detected
entirely by the connecting map on the IMAGE J0 of local invariants in
the top original Tate grade. It is not automatically the entire top
Tate space and is not the full-source intrinsic class. In particular
the graph-necessity theorem is retained, and strengthened: EVERY failure
requires a NONZERO intrinsic reduced LOCAL Tate obstruction at zero.
Genus and actual kernels no longer supply an independent residual gate.

## 5. Every polynomial H and finite intrinsic bad values

Let Bad_boundary(f) be the finite set of finite primitive parameters with
xi_boundary,a!=0. The geometric induction/ambient-SUM law for E and its
dual diagonal lattice, together with the cyclic coinvariant calculation,
sends the intrinsic class for H(f) at zero to

    (xi_boundary,a)_(H(a)=0), with factor1.             (13)

ALL roots, multiplicities and primitive components remain; no cross-root
cancellation or primitive purity is assumed. Combining (12),(13) gives
for EVERY f and EVERY nonconstant polynomial H, in EVERY degree:

    IF H(f) is eventually pure,
    full safety iff H^-1(0) misses Bad_boundary(f).     (14)

The global source-unit split status still reflects H before purity.
If the primitive source is globally nonsplit and H's roots miss the
finite bad set, H(f) CANNOT be eventually pure: (14) would make it safe
and force primitive source splitting. Thus every eventually-pure
composite of a nonsplit source must map SOME primitive bad finite value
to zero. If the finite bad set is empty, a globally nonsplit source has
NO eventually-pure polynomial composite in ANY degree. This includes
purely parabolic global obstructions not detected at finite punctures.
Hitting a bad value is not an existence theorem for pure composites.

The previous Long pure composite with scalar Norm1 remains a control:
its full tuple has NONZERO reduced local class at zero; no scalar norm
is substituted in (13). It is not exceptional invariant. Known actual
pure SAFE exceptional phases with unbounded top graph remain safe because
their intrinsic reduced local class is zero, not because their graph
space vanishes.

The complete exceptional target is now to prove xi_boundary,0=0 on the
ENTIRE exceptional eventually-pure locus, or realize an admissible
nonzero class with an actual pure invariant phase and original bad marker.
Equivalently prove the previously recorded constant-TYPE section theorem
on that full locus. No exceptional universal vanishing, pure producer,
coefficient-ring topological promotion, universal machine proof or review
gate is claimed. Both COMPLETE exceptional MZ cases remain OPEN.
