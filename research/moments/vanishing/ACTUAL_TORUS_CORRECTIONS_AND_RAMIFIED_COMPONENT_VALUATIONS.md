# Actual torus corrections and ramified component valuations

16 November 2026. Internal written deductions from named corrected
primary theorems. Independent replay, formal verification and external
review are not claimed. Both COMPLETE exceptional Mathieu--Zhao cases
remain **OPEN**; review deferred.

The remaining Neron component has an exact valuation vector in the ACTUAL
toric subgroup, even when the semiabelian extension is nonsplit and its
abelian quotient degenerates. One first removes the abelian specialization
by an identity-component lift, then evaluates torus characters on the
resulting torus point. Characters need not extend to the whole group.

This applies to EVERY canonical geometric section. It replaces neither
the full exceptional pure locus nor the unproved implication from actual
unweighted purity to vanishing. It identifies precisely which orders of
the actual Picard gluing data that implication must control.

See the [checker](verify_haar_toric_correction.py),
[receipt](haar_toric_correction_verification.json) and
[source audit](haar_toric_correction_source_audit.json).

## 1. Corrected Neron exactness, including degenerating abelian quotients

Let K=C((z)), R=C[[z]], and let

    0 -> T -> G -> A -> 0                               (1)

be the ACTUAL toric/abelian decomposition of a semiabelian variety.
Use Neron LFT models and write Phi for their component groups.
Halle--Nicaise, [arXiv:1209.5556v1](https://arxiv.org/pdf/1209.5556v1),
Chapter4, Propositions1.1.3 and3.3.1, give

    G^0(R) -> A^0(R) SURJECTIVE,
    Phi_T -> Phi_G -> Phi_A -> 0 EXACT,
    ker(Phi_T -> Phi_G) TORSION.                       (2)

The preprint numbering is used throughout. The complete statements and
proofs are on PDF45--46 and56--57. Section1.1 expressly REPLACES the
incorrect perfect-residue-field case of Bosch--Xarles Lemma4.2 and
checks that the repaired input suffices for the subsequent results.
We use the CORRECTED propositions, not that old lemma or merely a
reported correction. Their cited standard Neron/rank inputs remain
named imports; this is not a reconstruction of the whole monograph.

No good-reduction hypothesis on A, splitting of (1), or extension of
torus characters to G occurs in (2). Phi_A is finite. If T is split of
rank r, Chapter2,(3.4.3)--(3.4.4) identify

    Phi_T=X_*(T)=Hom(X^*(T),Z),
    t |-> (chi |-> ord_z chi(t)).                     (3)

Since Phi_T is torsion-free, (2) implies the exact sequence

    0 -> X_*(T) -> Phi_G -> Phi_A -> 0.                (4)

In particular the first map becomes an isomorphism over Q. It need
not be a primitive integral inclusion, and Phi_G may retain finite
torsion. The torus in (3) is the generic-fibre subgroup T of G; it is
NOT the possibly larger torus in the special fibre created by abelian
degeneration. Such degeneration contributes only the finite Phi_A
in (4), not additional free valuation coordinates.

## 2. Correct the ACTUAL section into T and take every character order

Assume first that T is split and let sigma be ANY point of G(K).
Choose n>0 killing Phi_A, or just the component of pi(sigma). Then
[n]pi(sigma) belongs to A^0(R). By (2) choose a lift

    b in G^0(R), pi(b)=[n]pi(sigma).

Use additive notation on G and put

    t=[n]sigma-b in T(K),
    v_z(sigma)=val_T(t)/n in X_*(T) tensor Q.           (5)

Equivalently, in ANY character basis chi_1,...,chi_r, its coordinates
are ord_z chi_i(t)/n. This is evaluation on the torus point t, not
on sigma in G. A nonsplit extension can have no suitable characters
on G at all; that causes no difficulty for (5).

The vector (5) is independent of BOTH b and n. Indeed, two lifts of
the same abelian point differ by an element of

    T(K) intersect G^0(R).

Its torus component maps to zero in Phi_G, and the injectivity in
(4) makes all its character orders zero. For different n, pass to a
common multiple and multiply the chosen lifts; the same argument
then proves independence after division by n. This also proves
additivity of v_z by adding identity-component lifts.

Under the rational isomorphism in (4), v_z(sigma) is exactly the
component nu_G(sigma). Consequently, BEFORE ANY PURITY,

    nu_G(sigma)=0 iff v_z(sigma)=0
      iff t belongs to T(R)
      iff ALL chi_i(t) are UNITS in R.                 (6)

The n chosen above suffices: the equivalent component statement is
that [n]sigma meets the identity component. Sigma itself need not
do so; finite component torsion is allowed. Neither sigma nor its
abelian image is required to be a torsion point of its generic fibre.

This supplies an actual coordinate interpretation of the previous
component test, rather than replacing G by a fictitious split product
T times A. The existence of b is the corrected theorem's content.
A general algorithm producing the actual formal lift b or the Neron
models from every phase is NOT supplied here.

Chapter2,(3.1.5) also gives a useful coordinate-free version. The
maximal quasi-compact open subgroup G^qc of the Neron lft model has
exactly the torsion components. Thus (6) is equivalent to

    sigma in G^qc(R).                                 (7)

This boundedness condition retains all finite components and is not
the assertion sigma in G^0(R).

## 3. A nonsplit torus, normalized covers and the original invariant space

For general T take a finite Galois splitting extension L/K. Over C((z))
it has the form z=u^e after a change of uniformizer. Work with the
ACTUAL base-changed semiabelian group and its Neron lft model over
R_L=C[[u]]. Choose n and b there as in Section2. Define

    v_a(sigma)=val_u([n]sigma-b)/(e*n)
       in X_*(T_L) tensor Q.                          (8)

The division by e normalizes the valuation to the ORIGINAL parameter z.
For infinity use z=1/s. The vector is unchanged after further finite
base extension, with character lattices identified: a prior lift b
remains an identity-component lift after base change, and each character
valuation multiplies by the further ramification index. Independence
of n and b then applies. A common splitting extension compares any
two choices. A change of character basis acts by GL_r(Z) and preserves
zero; an isogeny gives the corresponding rational lattice isomorphism.

The vector in (8) is fixed by the local inertia action I. For gamma
in Gal(L/K), gamma(b) and b lift the same [n]pi(sigma). Their difference
lies in T(L) intersect G_L^0(R_L), whose torus valuation is zero by (4).
Since valuations are Galois invariant on the field, (8) is therefore
I-invariant as a cocharacter vector. No assumption that b descends
to K is made.

The retained finite-cover reflection of rational component zero shows
that v_a=0 iff the ORIGINAL component is torsion. More precisely,
(8) induces a canonical rational isomorphism

    Phi_G tensor Q ~= (X_*(T_L) tensor Q)^I.            (9)

Injectivity follows from that zero-reflection and (4). For surjectivity,
clear denominators in an invariant rational cocharacter to obtain an
I-invariant integral cocharacter lambda. It descends to a K-cocharacter
of T. The point lambda(z) in T(K) has normalized vector lambda under
(8), so these points span the entire right side. This proves (9), not
just a dimension equality. It is functorial for semiabelian homomorphisms.

The dimension is the rank of the ORIGINAL invariant cocharacter space,
not necessarily the full split rank r. In particular an anisotropic
local torus has zero rational component space even if r>0; arbitrary
abelian degeneration remains permitted.

For the canonical 1-motive of a phase,

    M_f=[Z -> G_f], 1 |-> sigma_f,
    T_H(M_f)=(E_f(2))^dual,

this invariant cocharacter space is, up to the displayed Tate twists,
the DUAL of the original eigenvalue-one top Tate space. The previous
identity J_boundary=T^1 remains intact. The vector (8) is a coordinate
form of the same reduced local class; it is not the full unreduced
compact-middle class or a new genus obstruction.

## 4. Full MZ safety is exactly the corrected-unit condition

Combine (6),(8) with the
[full local theorem](FULL_INTRINSIC_LOCAL_TATE_CRITERION.md) and the
[Neron component criterion](LOCAL_TATE_NORMAL_FORM_AND_NERON_COMPONENTS.md).
For EVERY nonconstant regular phase f on SL2/Gamma, ANY finite Gamma,

    IF f is eventually pure,
    ALL original markers are safe iff v_0(sigma_f)=0
      iff ALL character values of [n]sigma_f-b
          are units on a local torus-splitting cover.  (10)

No genus vanishing, local genus semisimplicity, global genus--Tate
splitting, good reduction of A_f, or split extension G_f is assumed.
No general logarithmic 1-motive theorem requiring POTENTIALLY GOOD
reduction is applied to a degenerating abelian quotient. In particular
we do not replace an algebraic group by a rigid analytic uniformizing
cover and declare it algebraic.

In the actual Picard/NS/trivialization construction, t in (5) is the
distinguished gluing section after correcting its abelian specialization
by b. The remaining WHOLE exceptional implication is now precisely

    unweighted eventual purity of f
      => ord_u chi_i([n]sigma_f-b)=0 for EVERY i,       (11)

for Gamma=2O,2I, on ANY splitting cover. It requires the ACTUAL gluing
section and its ACTUAL identity-component lift; arbitrary torus factors
or characters of an unrelated split model do not establish (11).
The corrected propositions guarantee existence, not vanishing.

The [pullback controls](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md)
remain actual invariant phases with v_0!=0 and v_-1!=0 for EVERY Bezout
parameter. They show why exceptional invariance alone cannot prove (11).
Their one fixed original marker still detects every nonzero polynomial
composite, but supplies no unweighted pure example. The entire exceptional
pure locus, rather than that sufficient family, remains the main target.

## 5. Every recomposition and its ramification factors

Let H be ANY nonconstant polynomial and retain the established algebraic
1-motive law with semiabelian Weil restriction, diagonal lattice and
ALL-BRANCH section. At a primitive root a of H of multiplicity e_a,
write locally r=H(s), z=s-a, so r=c*z^e_a times a unit, c!=0.
Take a common splitting cover r=u^E, with E divisible by the required
branch and torus splitting indices. In each of the e_a branches over a,
ord_u(z)=E/e_a. Applying (8) to each ACTUAL toric correction gives,
under the branchwise cocharacter identifications,

    v_(H(f)),0 = direct_sum_(H(a)=0)
                    (v_a/e_a,...,v_a/e_a),
                  with e_a entries for that root.     (12)

The abelian correction causes no extra terms: the tuple of primitive
identity-component lifts is an identity-component lift after the common
base change. Its orders are independent of all other choices by Section2.
Units in the local parameter relation contribute no valuation.

Formula(12) is the NORMALIZED split-cocharacter description. Summing
within ONE ramification cycle gives e_a*(v_a/e_a)=v_a. This is exactly
the factor1 cyclic-coinvariant/Shapiro law for the intrinsic class.
Thus the factors1/e_a in split coordinates and factor1 in primitive
coinvariants are compatible. They are not interchangeable assertions.

Different roots and different torus directions remain separate summands.
No scalar norm can replace (12): opposite entries can cancel after a
norm even though the full vector is nonzero. In particular

    v_(H(f)),0=0 iff v_a=0 for EVERY root a of H.        (13)

Whenever H(f) is eventually pure, (13) is exactly full marker safety.
No purity of f is inferred; no H or phase is produced. All original
markers, ramification indices and primitive components are retained.

The entire exceptional MZ target is (11), or an actual invariant
eventually pure phase violating it with an original bad marker.
No universal exceptional vanishing, general identity-lift producer,
degree/support campaign, coefficient-ring promotion, universal machine
proof or external-review gate is claimed. Both COMPLETE cases remain OPEN.
