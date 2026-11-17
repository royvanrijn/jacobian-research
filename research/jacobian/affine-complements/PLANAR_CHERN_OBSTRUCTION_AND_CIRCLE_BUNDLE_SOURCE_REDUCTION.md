# Planar Chern obstruction and the complete circle-bundle source reduction

PF2RULED110, 16 November 2026. Throughout, N is the ORIGINAL finite
normalization of a nonautomorphic complex plane Keller map, its WHOLE
marked open U is A2, and every omitted prime is retained. The oriented
end is the outward boundary of its actual complex core. All statements
hold in every field degree, without bounds on indices or singularities.

**Planar obstruction.** If the actual positive contact end (W,xi) is
supported by a planar open book, then c1(xi) is NONZERO IN INTEGRAL
H^2(W;Z). In particular an actual integral-homology-sphere contact end
cannot be planar. The actual end cannot be S3. Together with PF109 this
excludes every oriented S3/Gamma for finite Gamma contained in SU(2),
including the trivial group, in the negative complex-link orientation.

**Circle-bundle source reduction.** If the actual oriented end is
L(p,1), the negative-Euler-p circle bundle over S2, p>=2, then N is
SMOOTH and has EXACTLY ONE omitted prime D, an embedded A1. Its index e
satisfies

    p>=3,   2<=e<=p-1,   e-1 = p modulo 2.              (1)

Its residue degree is one. The actual end knot K_D generates
H_1(W)=Z/p; its source intersection form is [-p], and its canonical
contact invariant is theta=-1-(e-1)^2/p. Thus EVERY singular-source or
multiple-prime packet in this entire oriented circle-bundle end class
is excluded, in ALL degrees. The surviving smooth one-prime class is
NOT excluded. No target knot type or field degree follows just from
this end identification.

**Sharp source controls.** For EVERY p>=3 there is a smooth affine
surface with a WHOLE plane open, one A1 complement, this oriented end,
and a globally exact plane volume with divisor (p-2)D. It saturates
e=p-1 in (1). No compatible finite Keller target map is supplied.
The remaining global finite-algebra implication stays OPEN.

## 1. The normal-source Chern bridge and planarity

Use [PF109](NORMAL_RESOLUTION_CHERN_GATE_AND_ADE_END_EXCLUSION.md)
with rho:Y->N the MINIMAL resolution, r>=1 strict omitted primes and
m>=0 exceptional components. Its compact resolved core admits a small
C-infinity-trivial Stein deformation V with the SAME boundary contact
isotopy class, and

    pi_1(V)=0,   b2(V)>=r+m,   c1(TV)!=0 in H^2(V;Q).   (2)

The nonzero Chern assertion uses the free algebraic Picard basis on
ALL resolution boundary components, constant units of the WHOLE A2,
and the rational Chern injection proved there. The deformation changes
the complex structure; neither the target map nor the marked plane is
claimed to remain holomorphic. We do not assume that the original
resolved Kahler form is exact.

Import Ozsvath--Stipsicz--Szabo, *Planar open books and Floer homology*,
Corollary1.5: a contact manifold with c1(xi)=0 that admits a Stein
filling of NONZERO integral first Chern class cannot be planar.
Rational nonvanishing in (2) implies integral nonvanishing, so its
contrapositive proves the stated obstruction. A rationally zero
boundary Chern class is NOT enough: nonzero torsion is allowed.
In an integral homology sphere H^2(W;Z)=0, giving the asserted special
case. Planarity means the genus of a CONTACT-SUPPORTING open book;
genus-zero bases of Seifert pieces, plane target curves, and planar
knot diagrams do not establish this hypothesis.

For S3, Stein fillability implies tightness. Eliashberg1992
Theorem2.1.1 says every positive tight contact structure on oriented
S3 is isotopic to the standard one. That structure is supported by
the disk open book: on |z1|^2+|z2|^2=1 take binding z1=0, fibration
arg(z1), and standard contact form alpha=sum r_j^2 dtheta_j.
On a page d(alpha)=2 dx2 wedge dy2 is positive, and on the binding
alpha=dtheta2 is positive. Thus it is planar and has c1=0, contradicting
(2) and Corollary1.5. Combining this with the nontrivial ADE quotient
groups already treated in PF109 gives the finite-SU(2) statement.
Other spherical space forms and opposite link orientations remain open.

## 2. All L(p,1) Stein fillings have b2 at most one

Use the PUBLISHED Plamenevskaya--Van Horn-Morris (PVHM2010),
*Planar open books, monodromy factorizations and symplectic fillings*.
Theorem1.2 and Corollary1.3, with the explicit reference filling in
Section2, give b2=1 for every Stein filling of every tight L(p,1)
when p!=4. Theorem1.2 also gives it for the virtually overtwisted
L(4,1) structure. The reference is D4 plus one two-handle along the
Legendrian unknot with tb=-p+1, so its intersection form is [-p].
The p-1 tight structures have rotation numbers

    R = -p+2, -p+4, ..., p-2.                          (3)

The two extreme structures are universally tight. The published
count p-1 is used; the earlier arXiv draft's different count is not.
These contact/filling classification theorems are imported, not
independently reproved here.

For completeness the exceptional p=4 rank bound needs NO uniqueness
claim for the rational-ball filling, and NO fresh import from an
unretrieved McDuff original or erratum. By PVHM Theorem1.1 (Wendl),
every Stein filling of either standard/conjugate contact structure
comes from an allowable positive factorization of its planar monodromy,
up to symplectic deformation. Its page is a disk with three holes,
and monodromy is one twist about each hole and one about all three
(PVHM Remark2.4). In a positive factorization let a_i count twists
enclosing just hole i, b_ij count twists enclosing exactly holes i,j,
and c count twists enclosing all three. Curves with no enclosed holes
are homologically trivial and do not occur in an allowable fibration.
Pairwise and capped-single-hole multiplicities give

    b_ij+c=1,  a_i+b_ij+b_ik+c=2.

Hence b_ij=1-c, a_i=c, and c is either 0 or 1. Factorization length
is respectively 3 or 4. This calculation counts enclosed subsets; it
does not identify isotopy classes of curves or claim uniqueness of the
three-twist filling. The Lefschetz fibration has chi=1-3+length.
For the simply connected V in (2), Stein two-complex topology gives
b1=b3=0; thus b2=length-3 is 0 or 1. That is all we need. The standard
four-twist case consists of boundary-parallel curves and has the
reference two-handle filling; the three-twist case has b2=0 and is
ineligible in (2). This explicitly handles the exceptional case
instead of silently applying Corollary1.3 at p=4.

## 3. Eliminate all singular and multiple-prime circle-bundle packets

Apply the preceding bound to (2):

    1 <= r+m <= b2(V) <= 1.

Therefore r=1, m=0 and b2=1. Since the minimal resolution of a singular
normal surface has a nonempty exceptional divisor, m=0 means N is
smooth. The omitted curve is an embedded A1 by the retained actual
normalization boundary theorem; no component-open deletion is used.
The b2=0 exceptional branch at p=4 is excluded already by this rank
inequality (also by the nonzero rational Chern class).

Now [PF102](PRIMITIVE_END_CANONICAL_CONTACT_AND_TRIPLE_TORUS_EXCLUSION.md)
applies without extra smoothness or sole-prime assumptions. It gives
the primitive marking, Q=[-p], and theta=-1-k^2/p with k=e-1>0 and
k=p modulo2. For the reference Stein surgery of rotation R in (3),
c1 evaluates as R, so c1^2=-R^2/p and theta=-1-R^2/p. Equality of
the boundary plane-field invariant forces k^2=R^2. Consequently
1<=k<=p-2 and k=p modulo2, exactly (1). There is no k at p=2.
For the surviving values the integral boundary Chern class is
-k[D] modulo p, nonzero, consistent with Section1.

For smooth ACTUAL N the residue-one theorem is already supplied by
the retained Borisov boundary result audited in
[PF105](TRIPLE_BOUNDARY_RESIDUE_ONE_AND_NORMAL_ITERATED_EXCLUSION.md).
It is credited as an import, not a new residue theorem. The finite
target compatibility remains an additional requirement. In particular
the end condition alone implies neither a torus branch-image knot,
an unknot, a monodromy group, nor the torus-only relation p|d.

## 4. An unbounded family of genuine affine, exact-volume source controls

For every integer p>=3 glue charts A2_(u,t) and A2_(v,s) over u*v=1 by

    s=u^p*t+u,       t=v^p*s-v^(p-1).                  (4)

This defines an A1 torsor over P1 under O(-p), hence a smooth surface
N_p. The FIRST ENTIRE chart is U=A2 and its complement is precisely
D={v=0}, an A1 with coordinate s. These are actual algebraic charts,
not only a formal boundary germ or a C* times C open.

Here is a direct affine compactification proof, without importing
the Danilov--Gizatullin isomorphism theorem. Homogenize (4) as
(s,z)=(u^p*t+u*z,z), giving an extension
0->O(-p)->E->O->0. In E(1) the chart transition is
(s',w')=(u^(p-1)*t'+z',u^(-1)*z'). The section (0,1) in the u chart
becomes (1,v) in the v chart and is nowhere zero. Thus O(-1) injects
into E with quotient O(-(p-1)). This extension splits, since
Ext^1(O(-(p-1)),O(-1))=H^1(P1,O(p-2))=0. Therefore

    P(E)=F_(p-2),   N_p=F_(p-2) minus S,
    S=C0+(p-1)F,    C0^2=-(p-2),   S^2=p,   S.C0=1.   (5)

The removed section is the subline O(-p) at z=0, whose normal bundle
is Hom(O(-p),O)=O(p); this identifies its positive self-intersection.
The displayed class follows from being a section and its square.
It is ample: a divisor a*C0+b*F is ample when a>0 and b>a*(p-2),
or directly use Nakai with the two ruled-cone generators C0,F.
Here a=1,b=p-1, so the complement is AFFINE. A tube around S is the
positive Euler-p disk bundle. The outward boundary of its complement
has the opposite orientation, giving the negative-Euler-p end L(p,1).
This is an actual oriented end identification, not recognition from
the determinant of a formal lattice.

The plane volume and its primitive extend globally:

    tau=du wedge dt=-v^(p-2) dv wedge ds,
    alpha=u*dt
         =v^(p-1) ds + p*v^(p-2)*s dv -(p-1)*v^(p-3) dv,
    d(alpha)=tau,      div(tau)=(p-2)D.                (6)

Every exponent in alpha is nonnegative for p>=3. Thus tau is a
globally EXACT regular two-form, standard on the WHOLE plane, and
its canonical coefficient is positive. The source has numerical
index e=p-1 and saturates (1). Its boundary Chern residue is 2 mod p,
nonzero for every p>=3. As a torsor with contractible fibers N_p has
the homotopy type of P1. Its rank-one intersection form has determinant
of absolute value p, by the boundary pair sequence and H_1(W)=Z/p.
Its compact generator is nonzero in the compactification and orthogonal
to the positive S; the Hodge index theorem makes its square negative.
Thus its intersection form is [-p].

The endpoint p=2 diagnoses the obstruction: (4) still gives the
affine quadric source, with zero canonical coefficient, but the
displayed primitive has a pole -dv/v and is not a regular global
primitive. This endpoint is outside the positive-canonical family;
failure of this particular primitive alone is not a proof that no
global primitive exists.

Equations (4)-(6) hold symbolically for arbitrary p, so the controls
are an UNBOUNDED construction, not extrapolation from sampled indices.
They supply no pair of regular functions (P,Q) giving a finite map
to A2 with dP wedge dQ a nonzero constant multiple of tau. Existence
of that compatible finite target would imply a Keller realization;
it is not asserted. The controls show exactly why source affineness,
the whole plane, exact volume, positive canonical coefficient and
this end topology by themselves do not finish the remaining class.

## 5. Evidence and remaining scope

The source pins and complete-page visual audits are in
[the manifest](sources/planar-circle-bundle/manifest.json): OSS full
PDF2/14; PUBLISHED PVHM full PDF2/4/5/12/13; Eliashberg full PDF6/8,
with Theorem2.1.1 on printed169/fullPDF6. The original McDuff PDF and
erratum returned HTTP403; neither was audited or made a new theorem
dependency. PVHM's original proofs and the written p=4 rank argument
above suffice for the needed bound. Deep OSS, Eliashberg and filling
classification results are imported, not independently reproved.

[The checker](verify_planar_circle_bundle.py) seals sources and
predecessors, verifies the arbitrary-p transition/primitive identities,
the exceptional factorization-length equations, contact theta and
parity identities, and scope countercontrols. Its finite sample
controls do NOT prove universal assertions; the written arguments do.
[The receipt](planar_circle_bundle_receipt.json) records the result.

Other graph/spherical ends, common-leading residue germs and the
[actual Laurent finite algebra](LAURENT_FINITE_ALGEBRA_FINISH_LINE.md)
remain OPEN. Contact planarity is not assumed from a graph or knot
diagram. Published Domrina2000 quartic retirement and shared trace
priorities are retained. Formal verification and external review stay
false; the full research goal remains active.
