# Canonical content, forbidden normal binary models, and the required conductor

15 November 2026. PF2RULED57 gives an all-field-degree format obstruction
for the actual finite normalization of a plane Keller map. It cannot be
a global affine complete intersection or a closed subscheme of a target
P1-bundle. Every admissible binary inverse must have a divisorial
normalization conductor. More generally, if the canonical boundary vector
is primitive, its only possible projective complete-intersection formats
are a cubic in P1 and two quadrics in P2; the former is already excluded
for a Keller map by the full boundary budget and Orevkov's maximal-index
remark. These assertions do not exclude the two-quadric quartic algebra
or higher-degree finite algebras with other resolutions.

## 1. The full boundary is a basis, in arbitrary degree

Keep the original normal finite normalization pi:N->A2 of a hypothetical
nonautomorphic plane Keller map on U=A2. Write E_1,...,E_m for ALL omitted
prime divisors, including unramified ones, and e_i for their transverse
indices. The retained finite-flatness theorem is PFN2. The complete
boundary Euler budget is PF2OMIT2:

    d-1 = sum_i e_i + sum_i E_i^ex,   E_i^ex>=0.       (1)

Residue degrees and coincident target supports are retained; (1) is
indexed by omitted primes, not by distinct image equations. The external
boundary hypotheses of that retained theorem remain premises here.

The divisor localization argument gives an isomorphism

    Cl(N) = direct-sum_i Z*[E_i].                     (2)

Indeed any Weil divisor becomes principal on factorial U; subtracting
that rational function leaves a boundary divisor. A principal boundary
combination would be the divisor of a unit on U, hence of a constant,
and is zero. This proof uses neither smoothness of N nor disjointness
of its omitted primes. Possible isolated boundary points do not change
the divisor class group.

In codimension one the divisor of pi^*(du wedge dv) is the tame different,
so its canonical class in (2) is

    [K_N] = sum_i (e_i-1)*[E_i].                      (3)

The target canonical bundle is trivial. All ramification is outside U.
At least one e_i>1: otherwise purity makes the full finite normal cover
etale, and connected finite etale triviality over complex A2 makes d=1.
These are the retained purity inputs of PF2PRIM1. In particular (3) is
nonzero and nontorsion. Define its positive integer content by

    g = gcd_i(e_i-1), ignoring zero entries.           (4)

The class is primitive exactly when g=1; one index-two prime suffices.
For the surviving quartic, the basis is [R],[D], the class is [R],
and g=1. No smoothness assumption is needed for (2)-(4).

## 2. Adjunction excludes whole classes of global presentations

First, N cannot be a global complete intersection cut out by a regular
sequence of polynomials in an affine-space bundle A2 x A^r, finite over
the target. Ambient canonical bundle and determinant of the normal
bundle are trivial, so Koszul adjunction would give [K_N]=0, contrary
to (3). This includes a monogenic normal algebra presented by one monic
equation. It does not prohibit local complete-intersection equations.

Now suppose N is a global relative complete intersection in a projective
bundle P(V)->A2 of relative dimension r, with fibrewise degrees
a_1,...,a_r. Equations can be twisted by target line bundles; all such
line bundles are trivial on A2. Remove fibrewise independent linear
equations by passing to their projective subbundle, so a_j>=2 and r>=1.
Let L denote the pullback of O_P(V)(1) to N. The determinant Euler
sequence for P(V), followed by Koszul adjunction, gives

    [K_N] = s*[L],   s=sum_j a_j-r-1,
    d = product_j a_j.                               (5)

Here s>=0. If s=0, (3) is contradicted. For s>0, the free boundary
basis forces s to divide every e_i-1, hence

    s divides g.                                    (6)

This is an equality of integral classes, not a numerical intersection
test. It is valid even if N was initially allowed to be singular:
the complete-intersection hypothesis supplies the relevant dualizing
line bundle, whose image in Cl(N) must still satisfy (3).

Since sum(a_j-1)=s+1 and a_j<=2^(a_j-1), (5)-(6) give the closed
degree bound for this presentation class

    r<=g+1,   d<=2^(g+1).                            (7)

If g=1, then s=1, so sum(a_j-1)=2. There are exactly two formats:

    r=1, (a_1)=(3), d=3;
    r=2, (a_1,a_2)=(2,2), d=4.                       (8)

In particular every degree>=5 packet with primitive canonical vector
requires a projective presentation beyond a global complete intersection
of hypersurfaces in one projective bundle. This applies in particular
whenever an index-two omitted prime occurs, with any other primes and
any residue degrees. It is not an exclusion of those finite covers.
For example (2,3) in P2 has s=2,d=6 and is not ruled out by this
argument when g=2. Treating g=1 as automatic in all packets would be
incorrect. The degree bound (7) is conditional on the stated global
complete-intersection format, not a universal Keller degree bound.

The adjunction input can be obtained from the ordinary conormal/Koszul
determinant calculation; its divisor version is Stacks Tag0B4A,
Lemmas48.14.1-2. No freeness of a vector bundle over A2 is required.

## 3. No normal binary model in ANY nontrivial Keller degree

For an embedding N->P(V) over the target with V of rank2, N is an
effective Cartier divisor of fibre degree d. Its class is O(d) times
a target line bundle, and (5) becomes

    [K_N]=(d-2)*[L].                                 (9)

For d>=3, at a ramified E_i, the integral boundary coefficient of L
is a positive integer. Thus e_i-1>0 and (9) imply e_i>=d-1.
The full budget (1) forces e_i=d-1, with just this one omitted prime
and zero excess. This maximal index is excluded by Orevkov's final
Remark in *On three-sheeted polynomial mappings of C2*, complete
PDF10, already retained in PF2RULED46. For d=2, (9) would make K_N
zero; for d=1 the proposed nonautomorphic packet is absent. Therefore

    the actual N has no closed target-relative P1-bundle embedding
    in any nonautomorphic Keller degree.              (10)

This proof removes the g=1 hypothesis specifically for binary models,
using the full boundary budget and the maximal-index exclusion.
It does not forbid a morphism N->P1 or a ruled completion: an inverse
coordinate provides such a morphism, but the joint map with pi can
identify points and has a nonnormal image. Replacing that image by N
without tracking its conductor would lose (10).

## 4. Exact class of the divisorial surface conductor

Let Z be any finite-flat degree-d binary root cover over A2, whose
normalization is the actual N, and let nu:N->Z be the finite birational
normalization. Set L=nu^*O_Z(1). The root cover is a relative Cartier
divisor, so omega_Z=O_Z(d-2). Let C_cond be the effective divisorial
part of its normalization conductor on N. Finite duality in codimension
one identifies the canonical module of N with the conductor times
nu^*omega_Z. Consequently

    [C_cond]=(d-2)*[L]-[K_N].                         (11)

This can also be checked by the local derivative/different formula
PF2OMIT4(D3)-(D4); it is an adjunction form of that exact local ledger.
We assert the Weil class and divisorial exponents, not absence of
additional codimension-two module questions in singular settings.

For an admissible source-finite primitive coordinate T from PF2OMIT3,
the whole infinity divisor is sum_i p_i E_i, with every p_i>0. Thus

    [C_cond]=sum_i((d-2)*p_i-(e_i-1))*[E_i].          (12)

These are class-vector coefficients. They need not be the actual
conductor multiplicities along E_i: C_cond can contain curves meeting
U, with no boundary component in its support. This distinction is
essential to applying the formula correctly.

For d>=4, (1) and Orevkov give e_i<=d-2 for every omitted prime,
including unramified ones. All coefficients in (12) are therefore
strictly positive, at least1. In particular C_cond is nonzero for
EVERY such coordinate, in EVERY field degree d>=4, at unrestricted
source pole orders and with all omitted primes retained.

Let J_Z be the target normalization-index determinant. The exact
trace ledger gives, at every target prime H,

    v_H(disc Z)=v_H(disc N)+2*length(B_H/O_H),
    sum_(E over H) f_E*c_E=2*length(B_H/O_H).         (13)

Thus Norm(C_cond)=div(J_Z^2), and J_Z is a nonconstant polynomial.
Neither the conductor nor this square factor can be omitted from
an inverse equation to obtain a normal model. Non-normality occurs
along a divisor; it cannot be confined to a finite set of surface
points. The finite-flat hypersurface order is Cohen-Macaulay/S2,
which also gives this last conclusion by Serre's normality criterion.

This C_cond is the conductor of the SURFACE map N->Z. The finite
normalization conductor of the TARGET CURVE Ctilde->C is a different
object, supported at finitely many curve points. QCL4's warning that
those curve-point corrections do not change Cl(N) remains valid.
No earlier curve conductor is silently made divisorial here.

## 5. The quartic cubic order has exactly the same index determinant

For the surviving quartic, write a binary inverse on an affine chart
as f(t)=a*t^4+b*t^3+c*t^2+d*t+e. Its monic cubic resolvent order has
the polynomial

    rho(s)=s^3-c*s^2+(b*d-4*a*e)*s
                    +(4*a*c*e-a*d^2-b^2*e).          (14)

Universally, disc(rho)=disc(f). The three generic roots correspond to
the three partitions of the four quartic roots into pairs; with the
actual S4 monodromy its fraction field is the actual cubic resolvent
field. This holds when a is nonunit: (14) defines a monic integral
cubic order and does not assert that this order is normal.

The normal cubic resolvent C3 of PF56 has discriminant unit*G_R,
as does the normal quartic algebra B. Comparing (13) for the two orders
gives, for every target divisor, equality of their normalization indices.
With target determinant-line trivializations,

    J_rho = constant*J_Z,
    disc(rho)=constant*G_R*J_Z^2.                    (15)

Thus the binary-order conductor and PF56's unavoidable cubic index
are the SAME target index, not two independent defects to clear
separately. Wood's Propositions6.1-3 give the geometric binary-quartic
and line-bundle-resolvent construction over any base; on A2 the line
bundle is trivial. The explicit identity (14) is also verified directly.

In the complementary r2q1 row of PF53, choose its coordinates so x
has pole1 along R and y is regular there, while y has pole p>0 along
the unramified D and x is regular there. Keep p unrestricted. For a suitable generic nonzero
constant gamma, T=x+gamma*y is admissible, separates every other
generic geometric sheet over the two boundary image curves, and has
infinity divisor R+pD. Equations (11)-(12) specialize to

    [C_cond]=[R]+2*p*[D].                             (16)

The minimal pole p=2 gives [R]+4[D]; this value is not assumed for
the whole r2q1 row.

At the generic point of C=pi(R), 1/T is a uniformizer of the ramified
quadratic factor; the other two sheets are finite and separated. At
the generic point of H=pi(D), the infinite D factor has degree1
and the other sheets are finite and separated. Hence the binary
order is normal at both generic target points:

    gcd(J_Z,G_R*H_D)=1, but J_Z is nonconstant.       (17)

The generic choice is justified by the finite number of geometric
sheet pairs: a given pair can agree on x+gamma*y for at most one
constant gamma. Intersect this finite-exception complement with the
admissibility and nonzero-gamma conditions of PF2OMIT3. The pole
orders distinguish the sole infinite factor from every finite factor.
Degree1 at D needs no primitive local generator beyond the base ring.

Consequently this entire r2q1 row requires at least one index-curve
component DISTINCT from both omitted image supports. At its generic
point pi is etale and the full fibre lies in U. The index records
identification of distinct source points by (pi,T), not new ramification
or another omitted prime. Distinct support does not mean disjoint
curves: this index curve may meet C and H at finitely many points.

## 6. The remaining two-quadric presentation and its global constraint

For g=1 the quartic (2,2) format satisfies adjunction with L=K_N;
it survives the format obstruction. If an actual full quartic
two-quadric presentation is used, with a globally defined conic pencil,
the binary cubic det(u*Q1+v*Q2) cuts out its cubic resolvent. The
standard oriented discriminant identity for quartic algebras and their
cubic resolvents is recorded in Wood (2011), Section3, PDF8. It must
be applied to the full quartic algebra, with its orientation, not an
arbitrary pair of matrices whose base locus has extra components.

Its cubic discriminant is the reduced G_R, so the resulting finite
cubic order is normal: at a DVR, a nonzero normalization index would
add an even positive discriminant exponent; finite flatness gives S2.
It is therefore the same actual C3 as in PF56. Its local product
description there is Gorenstein, with fibres2+1 along C. Hence its
binary cubic is nonzero in every target fibre, and its projective
zero scheme realizes the whole finite cubic cover; a zero form with
an entire P1 fibre is not being substituted for that algebra.

If the pencil had a globally everywhere-nonsingular member, the section
it defines in the pencil P1-bundle would avoid this finite cubic scheme.
The complement of that section is an affine-line bundle. Over affine
A2 its line-bundle and additive-torsor twists are trivial (Pic(A2)=0,
H1(A2,O)=0). It would give a global monogenic model of C3, contradicting
PF56. Thus every nonzero global polynomial pencil member, in any
matrix trivialization, must have nonzero NONCONSTANT determinant:

    det(lambda*Q1+mu*Q2) is never a unit.             (18)

The determinant is not identically zero for a nonzero pair: that would
give a rational point of the irreducible degree3 resolvent over the
target fraction field. Polynomial changes of pencil basis do not
remove this constraint. Locally a smooth conic can be chosen; no
global choice follows. This is a necessary gate on any actual full
two-quadric presentation, not a construction or exclusion of all such
pencils. We do not infer existence of a global matrix model from the
numerical format calculation (8).

## Evidence, controls, and next geometric question

[The checker](verify_canonical_formats.py) verifies exact class-vector
arithmetic, the complete primitive-format classification by the sum
identity, symbolic quartic/resolvent discriminants, local binary index
examples, and predecessor/source seals. It does not turn these finite
controls into a formal proof of normalization or adjunction.
[The manifest](sources/canonical-formats/manifest.json) pins the primary
adjunction statement, both Wood sources, the freshly inspected Orevkov
Remark, and the retained full boundary and conductor proofs.

The next full-algebra question is whether the marked source and trace
dual can coexist with the complete two-quadric pencil, its globally
unavoidable singular members, and the common binary/cubic conductor.
For higher field degrees with g=1, use presentations beyond the excluded
global complete-intersection formats and retain the same canonical and
conductor data. All results are necessary conditions or presentation
exclusions. No finite Keller pair, full quartic exclusion, higher-degree
exclusion, formal verification, independent full reproof, external
review, novelty priority or completion of the research goal is claimed.
