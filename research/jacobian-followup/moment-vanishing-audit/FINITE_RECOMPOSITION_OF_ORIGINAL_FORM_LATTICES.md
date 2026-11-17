# Finite recomposition of original-form lattices

16 November 2026. Internal written proofs with precise classical inputs.
Independent replay, formal verification and external review are not claimed.
Review stays deferred. Both COMPLETE exceptional Mathieu--Zhao cases remain
**OPEN**.

Polynomial recomposition can destroy cohomological tameness by introducing
nonisolated critical fibres. It cannot manufacture the ORIGINAL-volume
lattice defect needed for an eventually pure normalized period from a
primitive map whose original forms have no such defect.

More generally the good original-form lattice condition is preserved under
EVERY nonconstant polynomial recomposition. The original-form image remains
finite free when the starting map is tame, with an exact ordinary-rank
formula, even if its composite is not tame. No assertion about finiteness
of a raw Brieskorn quotient is made. The resulting exclusion applies to
all original forms, all markers and every recomposition degree.

See the [checker](verify_haar_finite_recomposition.py),
[receipt](haar_finite_recomposition_verification.json),
[source audit](haar_finite_recomposition_source_audit.json) and the preceding
[tame-period theorem](TAME_ORIGINAL_PERIOD_RIGIDITY_AND_PURE_BOUNDARY_DEFECTS.md).

## 1. The exact original-form image

Let X be smooth complex affine of dimension n>=2 and f:X->A1_s
nonconstant. Use the FULL ordinary top Gauss--Manin module

    G_f=Omega^n(X)[v]/(d-v*df wedge)Omega^(n-1)(X)[v],
    v=partial_s,      s=f-partial_v.                     (1)

Let M_f be the image of the ORIGINAL top forms Omega^n(X) in G_f.
This is a C[s]-module. The word image matters: for a nonisolated phase
it need not equal a raw quotient Omega^n/(df wedge dOmega^(n-2)).
Neither a localization nor the actual compact-period quotient replaces G_f.

At s=a use the decreasing V-filtration, with s-a times partial_s
having eigenvalue alpha on gr_V^alpha. Point modules have indices
-1,-2,... . Write L_a=V_a^(>-1)(G_f)_a, an Euler-stable finite analytic
O_a-module with logarithmic residue eigenvalues in (-1,0]. The classical
tame theorem gives

    (M_f)_a subset L_a for every finite a,                (2)

and M_f finite free over C[s], if f is cohomologically tame in Sabbah's
sense. The same inputs for M-tame functions are as stated in Nice
Proposition 3.2.3. These exact imports were pinned in the preceding note.
The new local transfer below applies to ANY regular holonomic G_f and L_a;
tameness is only one sufficient input for (2).

## 2. A coherent Euler-lattice lemma

Let B be a holonomic D-module germ in a coordinate r, and A a finite
O_r-submodule stable under E_r=r partial_r. If every residue eigenvalue
of E_r on A/rA is strictly greater than -1, then

    A subset V_r^(>-1)B.                                (3)

Only real-index geometric V-filtrations are used here. To prove (3),
finite generation over O_r bounds A below in the V-filtration. If its
lowest nonzero V index were beta<=-1, multiplication by r would put
rA at index at least beta+1. The nonzero lowest graded piece of A would
therefore survive as an E_r-stable subquotient of A/rA, with eigenvalue
beta. This contradicts the residue hypothesis. This argument also
handles possible O_r-torsion; it assumes neither flatness of A nor a
splitting of punctual terms from B.

## 3. The finite-pushforward identity, with its shift retained

Let h:A1_s->A1_r be a nonconstant polynomial and put P=h(f). Derived
direct images compose. Finite pushforward is exact on regular holonomic
modules, by Riemann--Hilbert and exactness of finite pushforward on
perverse sheaves. Consequently

    G_P=h_+G_f.                                         (4)

In the local graph/relative-de-Rham model, h_+G_f is the cokernel of

    partial_s-h'(s)u on G_f[u],
    partial_r=u,          r=h(s)-partial_u.              (5)

Let q denote the cokernel map on u-independent sections. The original
top form omega for P corresponds under (4) to q([omega]_f). This follows
from the transfer-module composition map: its original unit sends
omega tensor1 tensor1 to omega tensor1. In left-module notation the
one-form transfer factor is already included in (5); no additional
Jacobian factor may be dropped or inserted.

Near a point a over c, choose analytic coordinates z at a and r-c at c
so that h(s)-c=z^e. Changing s to z multiplies the one-form transfer
factor by an analytic unit, which preserves L_a. For any m in L_a,
(5) gives the EXACT identity

    e E_(r-c) q(z^j m)
       =q(z^j E_z m)+(j+1-e)q(z^j m),   0<=j<e.        (6)

Indeed q(partial_z(z^(j+1)m))=e q(u z^(e+j)m), whereas
(r-c)u acts on an original section as z^e u-1. The subtraction of 1
is essential. Omitting it would use the wrong left-module transfer.

The image A=q(L_a) is finite over O_(r-c), since O_z is finite free
with basis 1,z,...,z^(e-1). It is Euler-stable by (6). If alpha ranges
over the source logarithmic residues in (-1,0], the eigenvalues before
taking a quotient of A/(r-c)A are

    beta_(alpha,j)=(alpha+j+1)/e-1,  0<=j<e.             (7)

To see this with a nonconstant connection matrix, filter L_a/z^eL_a
by powers of z. Its diagonal blocks in (6) give (7); higher Taylor
coefficients only raise the z-filtration. The same conclusion holds
for nonsemisimple residue matrices and for the image quotient.

Every number in (7) lies in (-1,0]. In particular none is -1.
The Euler-lattice lemma therefore proves

    q(L_a) subset V_c^(>-1)(h_+G_f)_c.                  (8)

All points a over c contribute a finite direct sum. Thus if every
original form for f satisfies (2) at these points, every original form
for h(f) satisfies it at c. Ramification of any order does not create
an original-form index <=-1 from these good lattices.

## 4. Global consequences in EVERY recomposition degree

**Finite-recomposition theorem.** Suppose (2) holds at every finite a.
Then for EVERY nonconstant polynomial h,

    (M_(h(f)))_c subset V_c^(>-1)(G_(h(f)))_c
                                           for all c.   (9)

Every formal state of the exponential source for h(f) sends each
ORIGINAL top form either to zero or to a function which is NOT a finite
exponential-polynomial sum. Its ordinary moment resolvent is therefore
zero or nonrational. Proof: (8) gives (9). The original-form CRT
projectors and the V-strict point-module argument in the tame-period
note now apply without assuming that h(f) itself is tame.

In particular, for cohomologically tame f, or M-tame f with the stated
original-lattice theorem, NO nonconstant h(f) has an eventually pure
normalized original volume period. This handles arbitrary polynomial
degrees and the nonisolated critical fibres created by h'.

There is also an exact finiteness statement. Suppose M_f is finite
free over C[s], is contained in the generic connection, and has generic
rank rho=rank_(C(s)) G_f. These hypotheses hold in the tame setting above.
Then q identifies M_f, restricted to C[r] by r=h(s), with M_(h(f)):

    M_(h(f)) is finite free of rank deg(h)*rho over C[r]. (10)

Every original form maps by q, so the image is exactly M_(h(f)). Over
the unramified generic base, q is the usual one-form pushforward and is
injective. Its kernel is therefore C[r]-torsion. M_f is free and C[s]
is finite free over C[r], so that kernel is zero. This proves (10).
The image still generates the full source over D_r by construction.

The generic-rank identification used here is the ordinary relative de Rham
comparison: over the generic smooth fibre, every middle relative form can
be wedged with df to give an original top form, after inverting a nonzero
polynomial in the base. The generic Jacobian ideal is a unit ideal, so
this supplies all middle cohomology classes. Hence the original-form
image spans the full generic ordinary connection. Containment in the
good lattices excludes its intersection with punctual torsion.

Here rho is the ordinary generic Gauss--Manin rank, equivalently the
generic middle-fibre Betti number. It is NOT the generic Fourier rank or
the total Milnor number. Equation (10) does not assert finiteness of a
raw Brieskorn quotient, does not infer tameness from finiteness, and
does not remove punctual source or parameter-zero extension data.

## 5. A compulsory primitive defect above phase value zero

The local conclusion does not require (2) for ALL forms or ALL a.
For a particular original volume form Omega, (8) gives

    [Omega]_f in V_a^(>-1)G_f at EVERY a with h(a)=0
       => [Omega]_(h(f)) in V_0^(>-1)G_(h(f)).           (11)

An eventually pure normalized volume period is a nonzero polynomial
in exponential time. The original-unit defect theorem excludes its
class from V_0^(>-1). Therefore EVERY factorization P=h(f) of such
a phase has

    [Omega]_f notin V_a^(>-1)G_f
                       for SOME root a of h.            (12)

This is necessary, not sufficient for purity or failure. It does NOT
say that f-a is pure. Bad source classes can be killed by ramification,
so no converse to (8), or general reflection of V indices, is asserted.
For example if e>1 and m=delta_0, (5) gives q(partial_z m)=0;
punctual derivatives must not be treated as a generic connection.

Apply (12) to a canonical primitive factor. The retained normalization
argument writes every nonconstant regular invariant P on SL2 as

    P=h(f), f regular and primitive, unique up to affine change. (13)

It normalizes C[P] in its relative algebraic closure in C(SL2): the curve
is rational and has only one point at infinity because SL2 is rational
with constant units. Thus its ring is C[f]. The finite right group acts
on f affinely; after centering, it acts through a cyclic character.

For 2I this character is trivial, so f is already regular on SL2/2I.
For 2O it is trivial or the sign character; in the latter case f is
regular on SL2/2T and h is even. These are the prior exact descent facts,
not new invariant-ring calculations.

**Primitive-wildness consequence.** For EVERY nonconstant eventually
pure exceptional phase, its canonical primitive factor is itself
non-tame in the preceding cohomological/M-tame senses, and has the
original-unit defect (12) above phase value zero. This holds on SL2,
and on the quotient where f is actually regular: 2I, 2O in the invariant
case, or 2T in the sign case. No equivalence of tameness under an etale
cover is assumed; apply the theorem separately to the same normalized
pure phase on each of these spaces.

For an actual failure, the generic middle cohomology of that primitive
factor must also be nonzero. Otherwise a generic fibre of h(f), a
disjoint union of its deg(h) primitive fibres, has zero middle
cohomology and the retained full-source acyclic-fibre theorem proves
all-original-marker safety. This is only a necessary failure condition.

The full exceptional MZ problem has NOT been reduced to pure primitive
phases: (12) does not give purity of f-a, and the exact scalar-centre
reflection theorem concerns SAFETY, not this unproved purity reflection.
The complete remaining implication still includes all polynomial
composites of the genuinely non-tame primitive maps.

## 6. Genuine controls for the image and tameness distinctions

First take X=SL2 and f=trace=x+w. In coordinates A=(x-w)/2, its fibres
satisfy A^2+yz=s^2/4-1. A proper partial compactification is

    A^2+yz=(s^2/4-1)h^2 in P3 times A1_s.                (14)

The boundary h=0 is a smooth conic and the map is locally a product
there. The only critical points are two finite Morse points at s=+-2.
Thus f is cohomologically tame. Its generic fibre is an affine quadric
surface with middle Betti number 1.

Set P=f^2. Its critical locus contains the ENTIRE noncompact surface
f=0. Standard cohomological tameness with finite critical support is
therefore lost. Nonetheless (9)--(10) apply and M_P is finite free of
ordinary rank 2, with no original polynomial period. Direct Haar gives

    h(P^m)=binom(2m,m)/(m+1)>0 for all m>=0.             (15)

Only the central term x^m w^m in (x+w)^(2m) survives, with Haar value
1/(m+1). The unit ordinary resolvent is

    C_1^P(r)=(1-sqrt(1-4/r))/2,                         (16)

on its physical infinity branch. At zero it has a half-power singularity,
not the 1/r pole of a pure normalized unit. This is an actual non-tame
composite of a tame primitive map; it is not an exceptional phase.

Second, the raw quotient warning is substantive. On A2 with f=x and
P=x^2, the raw module Omega^2/(dP wedge dO) is C[y] dx wedge dy, of
infinite dimension and annihilated by P. But the ORIGINAL-form image
in the full top Gauss--Manin source is zero: for any Q, take
beta=-(integral Q dy)dx. Then

    (d-u*dP wedge)beta=Q dx wedge dy.

This proves the distinction without a source-rank estimate or a sampled
truncation. There is no normalized top volume state on this A2 control.

## 7. Sources and remaining full target

The new finite transfer and Euler-lattice arguments are written above.
Classical inputs are direct-image composition (Schnell, Proposition 17.8,
with the transfer/Spencer constructions in Lecture 17), finite perverse
exactness (Sabbah Nice Section 2.1(a)(7)), and the previously pinned
original tame lattice and strict V-filtration facts. The
[source audit](haar_finite_recomposition_source_audit.json) records the
exact pages and hashes. Selected exact controls check the Euler shift,
nonsemisimple residue blocks and the two geometric examples; they do not
machine-verify the universal functoriality or lattice proof.

The remaining full cycle-vanishing theorem must concern genuinely non-tame
primitive fibrations and all their composites, with an original-unit defect
at a root above zero. Recomposition-induced nonisolated critical fibres,
generic rank, finite original-form image or an abstract splitting cannot
replace the actual Haar-cycle coupling. Both complete exceptional MZ cases
remain OPEN. No degree/support/prefix campaign or external-review gate is
introduced.
