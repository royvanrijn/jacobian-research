# Intrinsic source splitting and full Haar safety

16 November 2026. Internal written proofs. Independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**; review stays deferred.

For EVERY eventually pure regular phase on SL2/Gamma, ANY finite Gamma,
full all-original-marker safety is equivalent to splitting the ambient
unit in the canonical GLOBAL boundary source. The actual relation space
is no longer an independent remaining gate. Constant boundary genus may
have arbitrary rank. All polynomial recompositions and all branches remain.

This strengthens the [constant-genus reduction](CONSTANT_GENUS_TRANSLATIONS_AND_ACTUAL_RECOMPOSITION.md).
Its abstract translation/actual-relation models did not impose geometric
admissibility. The mixed fixed-part theorem supplies that missing condition
for the geometric residual source. A local functional on its actual quotient
then suffices; a scalar-function realization of that quotient is unnecessary.

See the [checker](verify_haar_intrinsic_safety.py),
[receipt](haar_intrinsic_safety_verification.json) and
[source audit](haar_intrinsic_safety_source_audit.json).

## 1. Hypotheses and statement

Let X=SL2/Gamma and f be ANY nonconstant regular phase. On its FULL regular
base S, write E for the canonical boundary quotient of H3(X minus F_s),
constructed in the [boundary theorem](CANONICAL_BOUNDARY_REDUCTION_AND_LOWEST_WEIGHT_INVARIANTS.md).
This is an admissible graded-polarizable rational VMHS with

    W4 E=L=Q(-2),  P=GrW5 E,  T=GrW6 E.                 (1)

P is a polarizable pure genus variation; T has finite GLOBAL monodromy.
The full scalar Cauchy FUNCTION connection M is a quotient of the full
complement source. It contains its literal rational unit u=1, since
s*C_1-C_f=1. Removing the nontrivial compact-middle part preserves all
unit retractions. The subsequent genus quotient will also preserve them.
These constructions concern original regular markers, not an enlarged
amplitude algebra. Eventual purity means h(f^n)=0 for all sufficiently
large n; full safety means h(Q*f^n)=0 eventually for EVERY original Q.

The new theorem is

    under EVENTUAL PURITY:
    all original markers safe
      iff some global flat p:E -> Q has p(u)=1.          (2)

Here p is a LOCAL-SYSTEM retraction, not necessarily a Hodge morphism;
the target denotes the trivial local system, with no imposed Tate weight.
Rational and complex existence are equivalent: the splitting equations
are rational linear equations in a finite-dimensional space. This does
not assert that every eventually pure exceptional phase satisfies (2).

## 2. Mixed fixed part makes the second quotient admissible

Assume the right side of (2). Its restriction splits L in A=W5 E. Hence
A is semisimple as a rational LOCAL SYSTEM, because A/L=P is semisimple.
Set

    F=A^global,  F'=(A^dual)^global,
    C5=sum_g image(g-I)|A.                              (3)

The mixed fixed-part theorem, Steenbrink--Zucker Proposition4.19,
applies to BOTH A and A^dual: they are admissible graded-polarizable
variations over a smooth algebraic curve, including the required good
behavior at ALL boundary points. It identifies F and F' as CONSTANT
sub-VMHS. Their Hodge structures may be mixed and need not split by weight.

Semisimplicity of the underlying local system makes the evaluation
pairing F tensor F' -> Q perfect. It is a morphism of mixed Hodge
structures, so F' is isomorphic to F^dual as a MHS. Consequently

    A -> (F')^dual ~= F                                 (4)

is a VMHS morphism whose restriction to F is the identity. Its kernel
is precisely C5: it is the sum of the globally nontrivial isotypic
parts in the semisimple local system A. Therefore

    A=F direct_sum C5 AS VMHS,
    C5 is pure of original weight5.                     (5)

Indeed all of W4 A=L is invariant and lies in F, so W4 C5=0. This argument
uses the fixed part of the DUAL as well as the original variation;
an arbitrary local-system complement would not suffice.

Admissible graded-polarizable variations are closed under subquotients.
Equivalently here, the canonical Hodge subvariation C5 has its strict
geometric limit maps, so the quotient inherits the admissible limits.
Thus the residual source

    V=E/C5,   0 -> F -> V -> T -> 0                     (6)

IS an admissible VMHS. Its W4 is L and its W5 is the constant MHS F.
No algebraic curve/1-motive presentation or Hodge realization of the
ACTUAL period quotient is being asserted. The new geometric conclusion
is for V, rather than for an arbitrary kernel of scalar periods.

The source retraction p annihilates C5, which has no trivial quotient.
It descends to V and still has p(u)=1. As before, quotienting the actual
boundary period module by the image of C5 gives N, with u retained and

    Hom_D(M,C(s)) ~= Hom_D(N,C(s)),                      (7)

preserving evaluation on u. There is a differential quotient q:M -> N.
Neither q nor its kernel is required to respect a mixed Hodge structure.

## 3. Every LOCAL monodromy of V is finite

Fix ANY puncture of S in P1, INCLUDING INFINITY. Take a finite local
cover killing the finite monodromy on T and write Nloc=log(Tloc).
Nloc acts trivially on EVERY original GrW V: F is globally constant
and the top is now constant. Admissibility gives the relative monodromy
filtration M(Nloc,W). It follows that

    M(Nloc,W)=W,  Nloc(W_k) subset W_(k-2).              (8)

For clarity, the first equality does NOT follow merely by testing W
as a candidate before the second inclusion is known. On each GrW_k,
zero graded monodromy makes its induced M concentrated in degree k.
Two finite filtrations with this property coincide: simultaneous
splitting has nonzero components only on the diagonal. Hence M=W,
and the defining relative-monodromy inclusion gives the second statement.

By (1) and (6), (8) now says

    Nloc(V) subset L,  Nloc(F)=0.                       (9)

The descended local-system retraction p commutes with monodromy, so
p*Nloc=0. Since p is injective on L, (9) forces Nloc=0. Thus the original
local monodromy on V is finite. Every actual quotient N also has finite
local monodromy. This holds at phase zero even when zero is a regular
value; in that case its local monodromy is simply the identity.

This does NOT imply finite GLOBAL monodromy of V. After a finite global
cover the compactified base can have positive genus and nonzero handle
translations can remain. Global square-zero translation span R need not
vanish, and dim F need not be1. These facts do not obstruct the next step.

## 4. A local quotient functional splits the actual Borel unit

Here is the functional lemma used to pass through an actual, possibly
non-Hodge quotient. Let M be a regular-singular scalar-function differential
module containing its literal rational u=1. Suppose q:M -> N retains u
and local monodromy at s=0 on N has a semisimple eigenvalue1 block.
There is a monodromy-invariant covector ell on N with ell(u)=1.
Regular singularity identifies it with a horizontal local functional

    phi:N -> K0,   K0=convergent meromorphic germs at0.   (10)

Single-valued horizontal solutions of a regular-singular system have
at most a finite-order pole, so (10) lands in K0, not in a logarithmic
or essentially singular extension. It is a differential-module map;
phi*q fixes EVERY rational multiple of the unit. The finite-local-
monodromy conclusion of Section3 supplies precisely these hypotheses.

Compose with principal part, PP0:K0 -> delta0=K0/O0. This quotient is
the space of finite negative Laurent polynomials, as a Weyl module
for s and derivative_s. In particular PP0 is Weyl-linear and kills
C[s]. Therefore the restriction of PP0*phi*q to the actual Cauchy
functions descends through their polynomial quotient. Its restriction
to delta0 is the identity: rational negative powers of u are fixed.

The original exponential moment image Vtime is an indecomposable Weyl
module of entire functions. Its Borel map is

    B(sum a_n*t^n/n!)=[sum a_n*s^(-n-1)] modulo C[s],
    t acts as -derivative_s,  derivative_t acts as s.    (11)

All its elements have convergent Cauchy germs at infinity and admit the
regular-singular continuations used above. Under EVENTUAL PURITY its
nonzero polynomial unit moment generates Otime=C[t] inside Vtime.
The image of Otime in (11) is exactly delta0. Consequently the map
constructed from (10) retracts B(Vtime) onto delta0. Transporting across
the Borel correspondence gives a Weyl retraction

    Vtime -> Otime,  identity on Otime.                 (12)

Indecomposability forces Vtime=Otime. Every original marked moment
function is polynomial, proving full safety. A single compatible
continuation is enough: the resulting retraction on the actual image
is well-defined because it was obtained from a functional on M and
kills its actual relations and all polynomials.

Crucially, N has NOT been embedded into scalar functions. The FUNCTION
module is M; composition M -> N -> K0 is the needed operation. A global
retraction on N is unnecessary for constructing (12). Afterwards (12)
and full safety do imply the existing global actual-unit splitting theorem.

## 5. Necessity, the actual relation space, and every polynomial H

Necessity in (2) is the retained actual global-unit theorem: under eventual
purity, full safety gives a global unit retraction on M, hence on the
boundary actual quotient and on E by composition. Sections2--4 prove
the previously missing sufficiency. No fibre eigenvalue1 semisimplicity,
zero fixed genus, tame behavior or supplied relation-kernel calculation
is assumed.

The earlier constant-space criterion remains valid. If source splitting
holds, V has translation radical image R subset F and its actual quotient
has K=F intersect ker(V -> N). Under eventual purity it tested

    u NOTin R+K.                                        (13)

The present theorem now proves (13) AUTOMATICALLY for the genuine actual
kernel whenever eventual purity holds. This is not an unconditional
identity for arbitrary abstract K, and does not prove R=0 or K=0. In
particular, the prior abstract model translating top weight6 directly
into a weight5 fixed-genus vector fails (8) if used as LOCAL monodromy.
If such translations are handle monodromy on a finite cover instead,
the scalar actual-image/Borel argument still excludes its obstructing
K under eventual purity. No abstract model was an exceptional phase.

For EVERY nonconstant polynomial H, the intrinsic source-unit splitting
status is unchanged by recomposition, BEFORE purity, by the retained
finite induction/SUM and primitive-branch injection. Combining with (2):

    whenever H(f) is eventually pure,
    H(f) is fully safe iff the source unit of f splits.  (14)

All roots and all original markers are retained. No purity of f or f-a
is inferred. This removes the actual K_H test for EVERY such H at once.
It supplies no producer of eventual-pure composites.

There is a further exhaustive consequence. Recall the actual connection
injection M_f -> H^*M_H(f), so actual rank CANNOT DROP. If the source
unit splits and ANY eventually pure H(f) exists, (14) makes the composite
safe and its actual rank1. Hence r_f=1. The safe-centre theorem then gives
a UNIQUE c with f-c fully safe, and exactly the composites H(c)=0 are
eventually pure and fully safe. Thus for ANY f in the source-split branch:

* r_f=1: a unique safe centre c exists; eventually pure composites are
  exactly H(c)=0, and all of them are fully safe.
* r_f>1: NO eventually pure H(f) exists, in ANY degree, WHETHER OR NOT
  M_f is globally semisimple.

If the source unit is nonsplit, EVERY eventually pure H(f), IF ANY,
is a full failure. The source gate does not itself assert such H exist.
This strengthens the earlier semisimple-rank>1 exclusion to the entire
source-split branch.

## 6. Exact remaining whole-quotient target and source boundary

For EVERY eventually pure phase on EACH exceptional quotient, prove that
its canonical GLOBAL boundary source unit splits, or realize a genuine
admissible source-nonsplit eventually pure phase and original bad marker.
That one intrinsic condition is now necessary and sufficient. There is
no separate residual actual-relation or fixed-genus gate to close under
eventual purity. No universal source-splitting or pure-phase producer has
been established. The known Long pure failure has nonsplit source and
is not right2O/2I invariant, consistently with this theorem.

The deep new import is Steenbrink--Zucker Proposition4.19, with its good
variation hypotheses3.13, not an elementary theorem proved by the checker.
The original scanned primary article has been cached and its relevant
pages visually read. Its AppendixA and the published André Section5
provide the dual/fixed-part context; the prior canonical geometric
admissibility and strict-limit audits remain dependencies. Finite exact
linear-algebra and Laurent controls check transcription and distinguish
nonadmissible abstract data; they do not machine-prove the universal VMHS
or Borel arguments. No arbitrary-ring topological promotion is claimed.

Both COMPLETE exceptional MZ cases remain OPEN. Review remains deferred
until a substantial endpoint, ideally a whole quotient, is completed.
