# Canonical energy excludes every oriented ADE end

PF2RULED108, 16 November 2026. This is an ALL-DEGREE exclusion for the
ACTUAL smooth finite Keller normalization with its WHOLE marked plane
open and ALL omitted primes. Deep contact and filling classifications
are declared imports; their proofs are not independently reproduced.

**Theorem.** Let F:A2->A2 be a nonautomorphic complex Keller map and
pi:N->A2 its original finite normalization. If N is smooth, its oriented
end cannot be the link of ANY nontrivial rational double point:

    A_m (m>=1), D_m (m>=4), E6, E7 or E8.

Orientation means the outward boundary orientation of a compact core
in the complex surface. The excluded orientation is the boundary of
the NEGATIVE Cartan plumbing, equivalently the complex link orientation
on S3/Gamma for Gamma a nontrivial finite subgroup of SU(2).
Opposite orientations, arbitrary spherical space forms, singular N
and the full Wright goal are NOT covered by this theorem.

Equivalently the actual end cannot be the oriented boundary of ANY
connected NEGATIVE-DEFINITE TREE plumbing of spheres all of square -2.
This equivalent unbounded graph class is proved below without a finite
graph enumeration. There is no restriction on field degree, boundary
prime count or indices, infinity representations, or finite image
singularities. No actual omitted knot is identified with a plumbing core.

## 1. The source obstruction, retaining all primes

Use the exact actual-source results in
[PF103](MULTIPLE_BOUNDARY_CANONICAL_CONTACT_AND_A_TYPE_END_EXCLUSION.md).
The boundary N minus U is the disjoint union of r embedded affine lines.
Relative Thom and general position give H2(N;Z)=Z^r and pi1(N)=0.
In the compact basis dual to ALL omitted divisors the intersection
form Q is negative definite. Write B=-Q and k_i=e_i-1. The actual
Jacobian divisor is K_N=sum k_i D_i and k is NONZERO: otherwise the
original finite map is etale and the marked plane open is the whole
normalization, making F an automorphism. Thus

    E=k^T B^(-1) k>0.                                   (1)

A small strictly plurisubharmonic perturbation of the original radial
core X gives its actual positive tight boundary contact structure xi,
without changing the complex tangent bundle or compact topology. Then

    b2(X)=r, chi(X)=1+r, sigma(X)=-r,
    theta(W,xi)=c1(TX)^2-2chi(X)-3sigma(X)=r-2-E.         (2)

The symplectic form on this Stein core is EXACT. Stokes forbids every
closed symplectic surface, so X is symplectically MINIMAL. All these
claims concern the original finite normalization. Deleting other
boundary primes or resolving a singular source and treating its
exceptional curves as an exact Stein filling is not permitted.

Here is the reusable criterion. Suppose an oriented rational homology
sphere Y has a unique positive tight contact structure, a reference
filling with c1=0, b1=0 and negative rank L, and every minimal filling
of that structure has b2<=L. Its theta is L-2. If W=Y as ORIENTED
manifolds, uniqueness and (2) force E=r-L, while minimality gives
r<=L. This contradicts (1). In the ADE case the imported rank theorem
actually gives equality r=L.

## 2. Identify the actual contact structure on every D/E link

The A-type uniqueness and resulting exclusion were already proved in
PF103 using Honda's classification. For D and E write the negative
plumbing as a central -2 sphere with three all-minus-two arms of
lengths l1,l2,l3. In the surgery convention of Ghiggini's Figure2 this is

    M(-2; l1/(l1+1), l2/(l2+1), l3/(l3+1)).

The complete list is

| Type | Arm lengths | Rank L | H1 order |
|---|---|---:|---:|
| D_m, m>=4 | 1,1,m-3 | m | 4 |
| E6 | 1,2,2 | 6 | 3 |
| E7 | 1,2,3 | 7 | 2 |
| E8 | 1,2,4 | 8 | 1 |

These links are spherical, so Ozsvath--Szabo Proposition2.3 makes them
L-spaces. Ghiggini Theorem1.3 counts their positive tight structures
as the product of |a+1| over all coefficients a in the three negative
continued fractions -1/r_i. Every coefficient is -2, so the count is
ONE, for every m and each E type. In particular an arbitrary actual
Stein contact structure on this oriented end is the canonical one;
canonical identification is not being assumed from topology alone.

The sign is essential: the primary theorem is for e0=-2. Complete PDF
pages were rendered and read, including its surgery Figure2 and proof;
the text extractor drops minus signs and is not the orientation audit.
The L-space theorem and tight-structure count are imported. They do
not assert uniqueness of fillings by themselves.

## 3. Minimal filling rank supplies the contradiction

Ohta--Ono's Main Theorem states that a minimal symplectic filling of
the canonical contact link of a simple singularity is diffeomorphic
to its Milnor fibre. We use only its Betti-number consequence, not an
unsupported claim that the actual finite normalization is biholomorphic
to a model or that its omitted divisors are vanishing cycles.

The reference Milnor fibre has negative Cartan form and rank L. Its
first Chern class is zero: as a smooth affine hypersurface its tangent
bundle is stably trivial, by the exact tangent/normal sequence in C3.
Equivalently the minimal ADE resolution has trivial canonical bundle
and the same oriented filling topology. Hence

    theta_canonical = L-2.

Uniqueness of the contact structure identifies this with the actual
theta. Exactness makes X minimal, and Ohta--Ono gives r=L. Substitution
in (2) yields

    L-2-E=L-2,  so E=0,

contrary to (1). This proves the whole stated end class simultaneously,
including unbounded A and D ranks and arbitrary omitted-prime counts.
Together with PF103 it is the COMPLETE ADE end exclusion, not just an
E8 or one-prime consequence.

The retrieved primary Ohta--Ono source is the authors' MPIM2002/83
PostScript preprint, converted locally by Ghostscript solely as a reading
aid. Full pages1-4 were inspected, including Main Theorem on page2,
orientation definitions, and the ADE table. The published version is
J.DifferentialGeom69(2005),1-42, DOI10.4310/jdg/1121540338; its binary
was NOT obtained and is not falsely described as freshly audited.
The older Ohta--Ono1999 source separately supplies the negative Cartan
Milnor lattice and trivial canonical class; its full pages1-2 are audited.
The [manifest](sources/ade-end/manifest.json) pins original source bytes,
the derived reading PDF, source versions and exact imported scopes.

## 4. Why this is every negative-definite all-minus-two tree

Let C=2I-adjacency(T) for a finite connected tree T and suppose C is
positive definite. No vertex can have valency at least4: the principal
four-leaf star has a nonzero vector with value2 at its centre and1 at
its four leaves whose quadratic square is zero.

There cannot be TWO trivalent vertices. Take their connecting path,
with t vertices, and two off-path neighbors at each endpoint. On this
principal subgraph put value2 on the path and1 on the four leaves.
Its square is

    8t+8-8(t-1)-16=0.

Extra branches carry value zero and do not alter that square. This
again contradicts positive definiteness. Thus T is a path (A type)
or has precisely one trivalent vertex with three nonempty arms.

For an arm of length l, the chain determinant is l+1 and the inverse
entry at its attaching end is l/(l+1). These follow from the recurrence
D_l=2D_(l-1)-D_(l-2), D_0=1,D_1=2. The central Schur complement is

    S=2-sum l_i/(l_i+1)=sum 1/(l_i+1)-1>0.            (3)

Put p_i=l_i+1 and order 2<=p1<=p2<=p3. Equation(3) forces p1=2 and
p2=2 or3. For p2=2 every p3>=2 is possible, giving D_(p3+2).
For p2=3 it forces p3=3,4 or5, giving E6,E7,E8. No other case is
positive definite. This is an unbounded classification proof, not a
bounded graph search. It also gives the determinant and Seifert
orbifold Euler characteristic used in the table: detC=prod(p_i)*S,
and chi_orb=S>0, with rational Seifert Euler number -S.

## 5. Controls and exact limits

The reference ADE Milnor fibres are genuine Stein fillings with E=0;
they are not ruled out as fillings. The actual affine index-one quadric
control in PF103 also survives: nonzero actual Keller ramification was
essential to (1). For a formal ramified Cartan source packet, changing
one index from1 to3 makes E>0 and lowers theta below L-2, so it fails
the contact comparison if the oriented end is really the ADE link.

Conversely, blow up a reference filling once. Its abstract form becomes
C_ADE directsum[1], with characteristic vector(0,...,0,1), E=1 and
rank L+1. Its theta is STILL L-2 and |H1| is unchanged. Thus theta and
homology alone do not prove the exclusion; the MINIMAL filling theorem
is indispensable. The blowup carries a compact symplectic exceptional
sphere and is not an exact Stein source.

In particular the formal E8 rank8 packet retained in PF104, with index3
at the long-arm end, has E=8 and theta=-2. The oriented Poincare end
would require theta=6, so that end identification rules it out. Merely
having intersection matrix -E8 or |H1|=1 DOES NOT identify a boundary
with the Poincare sphere. The former numerical control remains valid
without that additional end hypothesis. No unrelated lattice or cover
is discarded because its determinant matches an ADE determinant.

Resolving singular N introduces compact exceptional curves and changes
the exactness/minimality argument. The theorem does NOT exclude singular
normalizations with these ends. No all-spherical-end or reversed-
orientation assertion is made. The common-tangent residue gap from
PF106/107, all f1 packets beyond the stated end exclusions, and ACTUAL
Laurent finite algebra remain OPEN. Published Domrina2000 quartic
retirement persists; external review stays deferred until closure.

[The checker](verify_ade_end_exclusion.py) and
[receipt](ade_end_exclusion_receipt.json) test exact lattice, contact,
Schur, rank and scope controls. The all-degree proof is the argument
above with explicitly imported classification theorems. Formal
verification, independent whole-classification proof and external
review remain false; the full research goal remains active.
