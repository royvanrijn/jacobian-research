# Multiple-boundary canonical contact and the complete A-type end exclusion

16 November 2026. PF2RULED103. Work over C. The finite normalization and
the WHOLE marked source are essential throughout.

**Theorem.** Let F:A2 -> A2 be a nonautomorphic Keller map, and let
pi:N -> A2 be the normalization of its target in C(A2_source). If N is
smooth, its oriented end is not L(p,p-1), for ANY p>=2. There is no bound
on the field degree, the number of omitted primes, their ramification
indices, or the finite singularities of their images. Orientation is
the boundary orientation of a compact core in the complex surface N;
L(p,p-1) means negative p/(p-1) surgery, equivalently the boundary of
the negative A_(p-1) plumbing.

The universal input is a formula for the actual source lattice and
canonical class. This extends [PF102](PRIMITIVE_END_CANONICAL_CONTACT_AND_TRIPLE_TORUS_EXCLUSION.md)
to arbitrarily many smooth boundary primes. It does not exclude
singular N, other oriented lens spaces, general graph-manifold ends,
or all Keller maps. Published quartic retirement remains intact.

## 1. Retain every boundary prime

Use U=A2_(x,y) inside the ACTUAL finite N, not a component-open model
formed by deleting other primes. The imported Borisov boundary theorem,
with its quantifiers audited in [PF2BATT1](../jc2-boundary-attachment/PROOF.md),
makes the boundary curves pairwise disjoint, unibranch, with normalization
A1, and smooth wherever N is smooth. Thus here

    N minus U = disjoint union D_1,...,D_r,  D_i = A1 embedded in N.

There are no extra isolated omitted points. More generally an affine
open subset of a normal affine surface cannot omit an isolated point
away from its omitted divisors. On an affine neighborhood V of such a
point q, disjoint from the omitted divisors, the coordinates x,y of U
extend across q by normal Hartogs extension. They define V -> U. Its
composition with U -> N agrees with V -> N on V minus q, hence everywhere
by separatedness, forcing q into U, a contradiction.

Relative Thom isomorphism for this disjoint smooth divisor, followed by
the pair sequence and the contractibility of U and D_i, gives

    H_2(N;Z)=Z^r; H_1(N)=H_3(N)=H_4(N)=0.              (1)

Choose the compact-cycle basis z_i by z_i dot D_j=delta_ij. In particular
the classes c1(O(D_i)) are the DUAL basis, not an arbitrary basis of
divisor coefficients. Also pi_1(N)=0: general position moves any loop
off the real-codimension-two D_i, so pi_1(U) -> pi_1(N) is surjective.

Take the original finite radial compact core X from
[PF98](PROPER_FINITE_MARKED_END_AND_COMPONENT_OPEN_AUDIT.md), with boundary W.
There are no radial critical points sufficiently far out. Off D this
follows from the Keller local inverse. On each D_i, its finite map to
the target is a nonconstant polynomial parametrization, whose squared
norm has positive radial derivative for sufficiently large parameter.
There are finitely many D_i. Thus X retracts N and each proper D_i has
one end circle K_i in W.

Let Q=(z_i dot z_j). In the PF98 smooth projective rational
compactification S of the ORIGINAL finite map, resolve only infinity.
The nef pullback H of a target infinity line has H^2=degree(pi)>0 and
is supported outside N. The z_i inject into H_2(S): their intersections
with the closures of D_j remain delta_ij. They are orthogonal to H.
Hodge index makes Q negative definite. Put B=-Q, positive definite.
Pair duality now gives

    H_1(W)=coker Q,  chi(X)=1+r,  sigma(X)=-r.          (2)

The classes K_i are the images of the dual basis in coker Q and
COLLECTIVELY generate it. No individual K_i is asserted primitive or
to have order det B. Its exact order is the least common denominator
of the entries in column i of B^-1. This is a useful marked comparison
with [PF101](FULL_PERIPHERAL_RATIONAL_LINKING_AND_LONGITUDE_JOINING_RECURSION.md),
whose filled homology and core orders must match the actual source.

## 2. Canonical energy of the actual source

Let e_i be the generic transverse ramification index along D_i, and
write k_i=e_i-1. The target volume pulls back to c dx wedge dy on U.
Since pi is finite between smooth surfaces, its divisor on N is

    K_N = sum k_i [D_i],  k_i>=0,
    <c1(TN),z_i> = -k_i,  Q_ii = k_i mod 2.           (3)

Here (3) uses the actual generic Jacobian divisor, not a reconstructed
boundary model. No assumption that every omitted prime ramifies is made.
For a nonautomorphic Keller map, k is NONZERO. Otherwise purity makes
pi finite etale everywhere; the connected cover of A2 is trivial.
Then N=A2. An omitted divisor would be principal on N and its defining
polynomial would pull back to a nonzero unit on U=A2, hence be constant
on dense U and N, impossible. Thus U=N and F is an automorphism.

Set

    E = k^T B^-1 k.                                  (4)

Positive definiteness and k!=0 give E>0. This is the ordinary rational
quadratic form of the source intersection lattice; it is not a
field-trace construction or an additional positivity hypothesis.

As in PF102, the unperturbed radial function can have degenerate Levi
form at ramification. On the compact core and a collar add a sufficiently
small positive multiple of a squared affine-embedding norm. This is
strictly plurisubharmonic, and its regular sublevel is isotopic to X.
Its complex tangencies give a positive tight contact structure xi on W.
With the actual complex structure retained, Gompf's invariant is

    theta(W,xi)=c1(TN)^2-2chi(X)-3sigma(X)
               = r-2-E.                              (5)

The Chern square is the rational inverse-form evaluation -E, since
H^2(W) is torsion. The same Stein core is symplectically minimal:
its symplectic form is dd^c(phi), hence exact, so Stokes forbids ANY
closed symplectic surface, including a symplectic exceptional sphere.
This proves the needed minimality, rather than inferring it from a
smooth diffeomorphism type or treating a topological blowup label as a
symplectic blowup.

## 3. Exclude every oriented A-type lens end

Suppose W is orientation-preservingly diffeomorphic to L(p,p-1), p>=2.
The negative continued fraction is a chain of p-1 entries -2.
Honda Theorem2.1 gives exactly one positive tight contact structure:
the product of |(-2)+1| is one. In the standard negative A_(p-1)
Legendrian surgery chain every unknot has tb=-1 and rotation zero.
Gompf's Stein handle theorem and Chern-evaluation formula give a Stein
filling with c1=0, chi=p and sigma=-(p-1). Its boundary has

    theta_unique = p-3.                              (6)

Uniqueness makes (6) the invariant of the actual xi after the oriented
diffeomorphism. The comparison does not assume that the actual omitted
knots are the plumbing surgery cores, or that the actual X is that
plumbing.

Import Aceto--McCoy--Park Theorem1.1: a minimal symplectic filling of
(L(p,q),xi) with finite fundamental-group order delta has
b2 <= p/delta^2-1. Here pi_1(X)=0, so delta=1 (NOT the Keller field
degree), and exactness proved minimality. Therefore

    r <= p-1.                                       (7)

Equations (5) and (6) imply E=r-p+1<=0, contradicting E>0.
This proves the theorem simultaneously for every p, every field degree
and every number of omitted primes. In particular RP3 is now excluded
for smooth actual N with arbitrarily many omitted primes, removing the
sole-prime restriction in the RP3 consequence of PF102. The opposite
orientation is not automatically included: for p>2 it is generally a
different oriented lens space. RP3 itself admits orientation reversal.

The selected primary statements and their hypotheses were read in full
pages, rendered and visually inspected; see the
[source manifest](sources/multiple-boundary-contact/manifest.json).
Deep filling/contact classification is imported, not independently
reproved. In the pinned AMP v1 the later NON-simply-connected computation
(3.7) has an apparent b2-1 versus b2+1 typo; no part of that calculation
is used. The d=1 proof is the separate first paragraph of Section3.2.
The theorem used here has the stated minus-one bound. The bound is also
consistent with its standard negative A-chain equality example.

## 4. Preserve the zero-canonical controls

For every r>=1, choose distinct a_0=0,a_1,...,a_r with a_i!=0 for i>0,
and let P(W)=product_(j=0)^r(W-a_j). The smooth affine surface

    N_r: f V=P(W)

contains the WHOLE marked A2_(x,f) by

    W=fx,  V=x product_(i=1)^r(fx-a_i).

Its exact complement is the disjoint union of r affine lines
D_i={f=0,W=a_i}. The inverse is x=W/f for f!=0, and
x=V/product_(i=1)^r(W-a_i) near the retained f=0,W=0 line.
Smoothness follows from P having simple roots: a singular point would
require f=V=P=P'=0. Near D_i the regular coordinate
t=(W-a_i)/f=V/(P(W)/(W-a_i)) gives

    f=u, x=a_i/u+t, dx wedge df=dt wedge du.

Thus every boundary index is one and the canonical class is zero.
These are genuine affine plane extensions, not finite Keller targets.
No claim that their volume is globally exact is needed. The r=1
quadric/RP3 control of PF102 remains valid. The all-r negative A-chain
Stein filling used in (6) likewise has E=0 and attains r=p-1; the
contradiction requires k!=0. We do not infer that an arbitrary affine
plane extension with k=0 is impossible.

## 5. Reproducibility and next boundary

[The checker](verify_multiple_boundary_contact.py) verifies exact inverse
forms, characteristic parity, marked core orders, rational congruence
invariance, symbolic all-r Gram and energy identities, the contact
comparison, source pins, predecessor seals and explicit affine controls.
[The receipt](multiple_boundary_contact_receipt.json) separates these
finite controls from the written universal proof and imported theorems.
No finite sweep is the reason the all-degree exclusion holds.

The original finite-algebra/Laurent-boundary implication remains open.
Use (2)--(5) with full marked peripheral data when attacking other ends;
do not identify a lens space merely from cyclic homology, drop primes
while retaining finiteness, or transfer the smooth theorem across
singular points without discrepancy/topological correction terms.
External review remains deferred under the closure-first priority.
