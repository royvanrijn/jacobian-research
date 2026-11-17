# Unipotent Haar extensions and general recomposition

16 November 2026. Internal written proofs, with precise classical inputs.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**. Review stays
deferred.

The intrinsic connecting obstruction is Tate even WITHOUT semisimple
fibre monodromy. More precisely it factors through the image of local
fibre invariants in the top Tate weight; it need not use that whole
weight space. If it is nonzero, every ambient lift moves. If it vanishes,
the remaining actual Haar motion has a canonical class modulo fixed
fibre cycles. This gives two separate gates for EVERY phase.

Intrinsic extension classes also reflect through EVERY polynomial
recomposition without semisimplicity or purity. Under the natural cyclic
coinvariant isomorphism, the composite class is the tuple of primitive
classes at the roots over zero. No multiplicity factor and no cancellation
occur at this quotient level. The earlier 1/e formula concerns fixed-vector
representatives under semisimplicity and is recovered as a special case.

See the [checker](verify_haar_unipotent_tate.py),
[receipt](haar_unipotent_tate_verification.json) and
[source audit](haar_unipotent_tate_source_audit.json).

## 1. The intrinsic class before any semisimplicity hypothesis

Let X=SL2/Gamma for ANY finite Gamma, f a nonconstant regular function,
and a a finite puncture of the FULL fibre/complement bifurcation system,
including finite bifurcations caused by infinity. On rational homology,
write the actual tube extension

    0 -> F=H_2(F_s;Q) -> V=H_3(U_s;Q) -> Q -> 0,      (1)

with monodromies A on F and T on V. For any lift gamma of the ambient
generator write (T-I)gamma=Tube(d). Its intrinsic class is

    xi_(f,a)=[d] in F/(A-I)F.                         (2)

Changing the lift adds (A-I)b, so (2) is ALWAYS independent of the lift.
Over Q, A-I is invertible off generalized eigenvalue 1; only that part
contributes. The class is zero exactly when a rational fixed lift of
the ambient generator exists. We do not assert an integral fixed lift;
the integral cokernel can contain torsion invisible to this rational test.

In cohomology the extension is

    0 -> L=Q(-2) -> V^c=H3(U_s;Q)
                          -> K=H2(F_s;Q)(-1) -> 0.    (3)

After finite base change making monodromy unipotent, let N denote its
logarithm. The connecting map

    delta:ker N_K -> L(-1),
    delta(v)=N_(V^c)(a lift of v)                      (4)

is well-defined. The Tate twist records that N is of type (-1,-1).
Restricted to the ORIGINAL eigenvalue-1 part, its vanishing is dual to
xi_(f,a)=0. One may equivalently use the original unipotent logarithm;
base change scales the connecting map by a nonzero integer. Original
non-1 eigenspaces cannot map to the ambient eigenvalue-1 line.

The existing canonical residue test computes the same zero/nonzero
condition, with its full complement connection and canonical lattice
retained. In particular

    delta!=0 iff L is contained in image N_(V^c).      (5)

For the forward direction use a nonzero connecting value. Conversely,
if a lift maps nontrivially into L, its image in K is N_K-invariant,
so gives a nonzero value in (4). Formula (5) is finite linear algebra
on the actual canonical system, not a claim that such a system has
been constructed for every phase.

## 2. Exact limits and an elementary bound on invariant vectors

Use W for the ORIGINAL fibre weight filtration and M=M(N,W) for the
relative monodromy filtration. They are generally different. On (3),
the original weights of K are 4,5,6 and L is pure of weight 4.

We need exactness for the LIMIT mixed Hodge structures, not merely for
unfiltered vector spaces. Here is why it holds for the geometric sequence.
After unipotent base change, the canonical extension construction is exact
on local systems: in the exponential logarithmic frame, a flat map is a
constant linear map intertwining the nilpotents. The limiting Hodge maps
preserve F: they do so on the punctured disc and the extended filtrations
are subbundles. They preserve M by functoriality of the relative filtration.
That functoriality follows directly from Deligne's inductive kernel/image
identities in Proposition1.6.13. Thus these are morphisms of the limit
mixed Hodge structures supplied by geometric admissibility. Ordinary
mixed Hodge strictness then makes their underlying exact sequence strict
for M and F. The same applies to every W-truncated geometric sequence.

This does NOT say that arbitrary sequences of relatively filtered vector
spaces are M-strict, or that the entire category of limit data is abelian.
Brosnan--El Zein explicitly warn that their extra original filtration
need not be strict in arbitrary limit-data sequences. Here original
W-strictness comes from the nearby geometric mixed Hodge sequence and
its canonical flat extension; strictness for M comes from the limit
mixed Hodge maps just described.

**Invariant-vector bound.** If a finite filtered nilpotent space (B,N,W)
has relative filtration M and W_b B=B, then

    ker N subset M_b B.                              (6)

Proof: for ell>b and every original weight k<=b, N is injective from
gr_M^ell gr_W^k B to gr_M^(ell-2) gr_W^k B. This is the positive-degree
injectivity in a monodromy filtration centred at k. Hence N is injective
on gr_M^ell B itself: its maps on every associated W grade are injective.
The highest nonzero M grade of a vector killed by N therefore cannot
exceed b. This proves (6) without assuming that N is strict for original W.

The relative filtration restricts to W_5 K, by its defining properties
and uniqueness. Applying (6) gives

    ker N_K intersect W_5 K subset M_5 K.             (7)

## 3. The general intrinsic obstruction factors through Tate invariants

In the exact limit mixed Hodge sequence, (4) is the connecting morphism
of the two-term complexes defined by N. Its target L(-1)=Q(-3) is pure
of LIMIT weight 6. Therefore (7) and mixed Hodge strictness imply

    delta(ker N_K intersect W_5 K)=0.                 (8)

Equivalently choose an M_5 lift using the strict limit sequence. Its
image under N lies in M_3 V^c, whereas the ambient line is pure M_4;
that image must be zero. The twist in (4) and these two weight indices
give the same argument in two conventions.

Let J_a be the image of the ORIGINAL eigenvalue-1 local invariants in
the top original weight:

    J_a=image(ker N_(K,1) -> gr_W^6 K).                (9)

The top original piece is pure Tate, with finite monodromy. After the
chosen base change its monodromy logarithm is zero, so its limit is still
pure Tate of weight 6. The image in (9) is a rational Tate subspace.
Working first with the whole unipotent base change avoids assuming that
the original spectral projector preserves the moving nearby Hodge
filtration. Restriction to the original eigenvalue-1 Tate summand is
then harmless. By (8), there is a unique induced morphism

    delta_J:J_a -> Q(-3)                              (10)

carrying ALL information in the intrinsic obstruction xi_(f,a).

**General Tate-obstruction theorem.** Equations (8)--(10) hold for every
phase and finite puncture, WITHOUT semisimple fibre1, purity, tameness
or isolated singularities. Under eventual purity at0, delta_J!=0 is a
full all-original-marker FAILURE certificate: every ambient lift moves,
and the full de Rham/cycle theorem detects it by a fixed original marker.

The space J_a is the IMAGE of actual local invariants in the top Tate
weight, not automatically that entire Tate weight. That top weight is
the boundary graph space H_1(Delta(D))(-3), using K=H2(F)(-1).
Statement (10) does NOT identify the full actual variation with a class
in W_(-4)H_2(F) when fibre1 is nonsemisimple. It confines the intrinsic
COINVARIANT extension obstruction; additional fibre-internal motion
can remain. No canonical lift of the dual J_a class into fibre homology
is asserted.

## 4. The exact remaining actual-lift class when xi=0

If xi_(f,a)=0, choose a fixed lift gamma_fix in the generalized eigenvalue-1
extension of (1). Let gamma_1 be the projected actual physical Haar lift.
There is b in F_1 with

    gamma_1=gamma_fix+Tube(b).

Its class

    beta_gamma=[b] in F_1/ker(A-I)                     (11)

is independent of the chosen fixed lift, since two fixed lifts differ
by a FIXED fibre class. The actual projected variation is

    pi_1(A)d=(A-I)b.                                  (12)

The induced map from the quotient (11) to image(A-I) is an isomorphism.
Thus, for EVERY eventually pure phase at0, there are precisely two gates:

* If xi!=0, equivalently delta_J!=0, the phase fails, independently of
  a choice of ambient lift.
* If xi=0, all-original-marker safety is equivalent to beta_gamma=0.

This is a complete decomposition of the existing all-marker criterion,
not a proof that either obstruction vanishes universally. The quotient
in (11) is by FIXED cycles, and is NOT the coinvariant quotient in (2).
A split extension can have nonzero (11) for a specified moving lift.
When fibre1 is semisimple, (11) is zero automatically, recovering the
previous intrinsic complement criterion. All nonsemisimple cases remain
in scope; no arbitrary fixed lift replaces the physical Haar lift there.

## 5. Intrinsic extension classes reflect through ALL recompositions

Let P=H(f), H nonconstant, and let a run through its roots above0,
with multiplicities e_a. No semisimplicity is assumed. The composite
fibre is the disjoint union of the primitive fibres over its nearby
root branches. On an a-block, after compatible flat identifications,
its monodromy has the cyclic form

    A_P(v_0,...,v_(e-1))=(A_a v_(e-1),v_0,...,v_(e-2)).

There is a natural isomorphism of RATIONAL coinvariants

    Phi:F^P/(A_P-I)F^P -> direct_sum_a F^a/(A_a-I)F^a,
    Phi_a([(v_j)_j])=[sum_j v_j].                     (13)

Indeed the cyclic relations identify adjacent summands, and the final
wrap relation is exactly (A_a-I)v. This proves both descent and
bijectivity. A different cyclic starting branch only applies A_a,
which acts trivially on these coinvariants.

**General recomposition theorem.** Under (13),

    Phi(xi_(H(f),0))=(xi_(f,a))_a.                    (14)

Proof: choose a common positive multiple L of all e_a and set k_a=L/e_a.
The exact common-loop tube formula, valid without semisimplicity, is

    [S_L(A_P)d_P]_(a,j)=S_(k_a)(A_(a,j))d_(a,j).

In composite coinvariants S_L acts as L. After summing the a-branches
in primitive coinvariants, the right side is k_a times e_a copies of
the SAME intrinsic primitive class. Independence of the primitive lift
in (2) makes those classes identical after compatible transport.
Thus L*Phi_a(xi_P)=L*xi_a. Division by L over Q proves (14).

There is NO 1/e factor in (14): it is a sum on the coinvariant quotient.
Under semisimple fibre1, the fixed representative of the composite
class has e equal branch components c_a/e, whose sum is c_a. Hence (14)
recovers the earlier exact 1/e theorem rather than contradicting it.
Unipotent Jordan blocks are fully retained in (13)--(14).

Equivalently, the composite rational tube extension has a fixed ambient
lift iff every primitive extension over a root has one. Another direct
proof uses the fibre product of their tube extensions over the same
ambient Q. It is isomorphic to the composite extension by the maps to
every primitive complement. Choose primitive fixed lifts, obtain an
L-loop-fixed composite lift, then average its first L translates to
obtain a one-loop-fixed lift. This too uses no semisimplicity.

## 6. A finite intrinsic failure set in EVERY primitive fibration

For a fixed f define

    E_f^ext={a in the FULL finite bifurcation set : xi_(f,a)!=0}.

This is defined WITHOUT any fibre semisimplicity assumption. From (14),
for EVERY nonconstant polynomial H,

    xi_(H(f),0)!=0 iff H^-1(0) meets E_f^ext.           (15)

If H(f) is eventually pure, meeting this finite set is therefore a
full-original-marker failure certificate in ANY recomposition degree.
Avoiding it proves that some fixed lift exists; in the nonsemisimple
case the remaining beta of the ACTUAL Haar lift still must be tested.
An empty intrinsic failure set is NOT a universal positive MZ theorem.
Nor does a bad puncture by itself produce a pure polynomial composite.

All statements apply on the actual regular quotient of a canonical
primitive factor: 2I, 2O, or 2T for the octahedral sign case. Finite-cover
marker trace, the genuinely non-tame unit defect and all original
markers/branches remain. The complete remaining problem is to exclude
or realize delta_J, and when it vanishes to control beta_gamma, for
EVERY eventual-pure exceptional phase and composite. Both whole cases
remain OPEN. No general limit-MHS/canonical-lattice/cycle producer,
all-ring topological promotion, numerical bound, degree/support campaign
or external-review gate is introduced.

Classical imports are the geometric limit MHS and logarithmic extension,
Deligne mixed-Hodge strictness and the relative-filtration kernel/image
construction. The invariant-vector bound, the all-phase Tate factor of
the connecting obstruction, the residual-class decomposition and the
general coinvariant transfer are written proofs above. Selected finite
controls check normalization and logical boundaries, not the universal
geometric functoriality or existence arguments.
