# Original log symbol and two-orbit boundary support

16 November 2026. Written uniform calculation of the ORIGINAL Haar
clutching class and its graph residues. Both FULL exceptional
Mathieu--Zhao cases remain **OPEN**. Review remains deferred.
No degree campaign or new holonomic calculation is used.

The [contact primitive](ORIGINAL_CONTACT_PRIMITIVE_AND_POLAR_DIVISOR_RESIDUES.md)
now has an explicit two-chart transition. It is a fixed two-logarithm
symbol in the minimum bottom-row Klein forms. Its node residues are
determinants of integer valuations. In the remaining dominant-left-root
region, nonzero residues can occur only above TWO right-root orbits.

This computes the clutching representative, NOT the intrinsic extension
class. An actual pure SAFE control below has nonzero symbol residues
on the graph with the extra Klein cuts, although its original fibre
has H^2=0. The cut graph and the intrinsic graph must not be identified.

## 1. The original Haar class has a fixed logarithmic transition

Use the original X=SL2/Gamma, Gamma=2O or 2I, and the bottom-row forms
F=F_n(z,w), G=G_m(z,w) of the retained
[Klein model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md). Here

    (n,m)=(8,12) or (12,20),
    nm=2|Gamma|=4N, N=24 or 60.

They have no common zero on the nonzero bottom row. Thus the ORIGINAL
opens X_F={F!=0}, X_G={G!=0} cover X. Their effective charges are
k_F=-n/2 and k_G=-m/2. The original contact primitives satisfy

    B_F=-lambda wedge dlog F/(2k_F),
    B_G=-lambda wedge dlog G/(2k_G),
    dB_F=dB_G=Omega.                                 (1)

On their overlap the exact transition is

    T=B_F-B_G=(2/(nm)) dlog F wedge dlog G
              =(1/|Gamma|) dlog F wedge dlog G.      (2)

To verify equality as FORMS, use a chart with bottom row proportional
to (1,u) and projective determinant r. Then

    lambda=dlog r-2v du,
    dlog F=-(n/2)dlog r+dlog F_n(1,u),
    dlog G=-(m/2)dlog r+dlog G_m(1,u).

Substitution gives (2), including its coefficient. The term v du
wedges to zero with either bottom-projective logarithmic derivative.
This is not a declaration that lambda is closed.

In the de Rham Mayer--Vietoris sequence of X_F,X_G, the connecting
image of [T] is [Omega], with the sign fixed by the difference in (2).
Indeed the two primitives in (1) are precisely its Cech descent data.
The image is the NONZERO original ambient Haar line, not a new
Gaussian phase or a weighted Haar functional.

There is a rational Betti lift as well: the two dlog classes divided
by 2*pi*i are the pullbacks of the integral circle classes on C*.
Their cup product is a rational class on the overlap, and (2) fixes
its coefficient and comparison factor (2*pi*i)^2. Thus the transition
has full ambient cohomological meaning, not merely a pointwise identity
of meromorphic functions.

The old additive-slice pole/transgression theorem is credited in
GHAARPOLERESIDUE20261115. The new calculation is this particularly
simple contact transition and its uniform graph valuation data.

## 2. A single coordinate on the RIGHT projective quotient

Put d=gcd(n,m)=4 and

    q=F^(m/d)/G^(n/d)
      =F^3/G^2 on 2O, F^5/G^3 on 2I.                (3)

It has charge zero and depends ONLY on the bottom row line. The
projective forms in its numerator and denominator have no common
root and have equal degree nm/d=N. Therefore q:P1->P1 has degree N.
It is right-G invariant, G=Gamma/{+1,-1}, of order N. Every generic
orbit already has N points, so each generic q-fibre is exactly one
orbit. Consequently

    C(P1)^G=C(q).                                    (4)

In particular q is a coordinate on P1/G, not a further nontrivial
cover of that quotient. In the same chart the coefficient calculation
in Section 1 gives

    dlog F wedge dlog G=2 lambda wedge dlog q,
    T=(1/N) lambda wedge dlog q.                     (5)

This form is closed, because d lambda is horizontal and q depends
only on the right projective coordinate u. Formula (5) does NOT make
the original three-form exact globally; the primitives remain on
the two different opens.

The zeros and poles of q are two distinct right-root orbits. For O,
F_8 has 8 simple roots and

    G_12=z^2 w^2(z^4-w^4)^2

has 6 distinct roots, each of multiplicity 2. The minimum O orbit
size 6 makes each set one orbit. For I, F_12 and G_20 have 12 and
20 simple roots, each one orbit. Simplicity follows also from the
minimum I orbit size 12: a multiple root of either form would use
at least 24 degrees. Two distinct root orbits would likewise need
at least 24 roots. These use the retained full Klein/orbit data,
not a finite phase-degree check.

Thus on the quotient there are just the two distinguished points
q=0 and q=infinity for this transition. This statement concerns its
zero/pole support, not the full bifurcation set of an arbitrary phase.

## 3. Exact graph residues and tame functions

Fix ANY nonconstant original phase f on a common regular parameter
open where its cut fibres have their generic components. Shrinking
that open does not discard branches in a later local comparison.
The transition (2) restricts to the CUT fibre

    F_s^*=F_s intersect {F G!=0}.

Choose a smooth proper SNC compactification of this surface, resolving
also the zero/pole divisors of the two original functions F,G. This
may require a finite parameter cover to label components, as in the
[algebraic residue theorem](ALGEBRAIC_HAAR_LOG_COEFFICIENTS_AND_ORIGINAL_RESIDUES.md).
At an oriented node D_i intersect D_j use local parameters u,v and
write

    F=u^(a_i) v^(a_j) A,
    G=u^(b_i) v^(b_j) B,

with A,B invertible. The exact double-log coefficient is

    r_(T,ij)=(1/|Gamma|)(a_i b_j-a_j b_i).           (6)

Only the integer divisorial valuations enter. Higher jets of the
units cannot change (6). In particular all these coefficients lie
in (1/48)Z or (1/120)Z, uniformly in every phase degree and support.

There is also an actual rational function on each component:

    tau_i=(-1)^(a_i b_i)
                 (F^(b_i)/G^(a_i))|_(D_i),          (7)

where restriction means cancellation of its order along D_i first.
This is the elementary tame symbol notation, not an import of a
regulator comparison theorem. Its node orders and the single residue
of T are

    ord_j tau_i=b_i a_j-a_i b_j,
    Res_i T=-(1/|Gamma|)dlog tau_i.                  (8)

All zero/pole points of tau_i occur at the included SNC boundary
nodes. The degree-zero divisor of this rational function therefore
gives the signed vertex balance for (6). Hence the residue vector
is a graph cycle, as the logarithmic comparison requires. Equations
(6)--(8) retain the actual functions as well as their integer orders;
the functions may contain essential parameter-dependent units.

In particular, at any node above a point INSIDE the original affine
F_s, the coefficient is zero. The two original regular functions
have disjoint zero sets there. Near a zero of F, G is a unit, so
both b_i and b_j vanish; the analogous argument holds for G. All
nonzero graph residues of (2) are consequently over original infinity.
This does not say they belong to the graph of the uncut F_s: added
Klein divisors can meet the old infinity boundary and create cycles.

## 4. Two-orbit support above the remaining dominant-left-root region

For an eventually pure phase align its unique strict leading left
root to ell={x=y=0}, using the retained
[leading-root calculation](DOMINANT_ROOT_RESIDUE_AND_FIBRE_GRAPH_GLUING.md).
The bottom row is nonzero everywhere on this line. At a right point
where both F_n and G_m are nonzero, choose a bottom-row projective
chart. The EXACT original functions near that point are

    F=Delta^(-n/2) F_n(1,u),
    G=Delta^(-m/2) G_m(1,u),                          (9)

and both right factors are units. For EVERY divisor over this open,

    a_i=-(n/2) v_i(Delta),
    b_i=-(m/2) v_i(Delta).                           (10)

Their valuation vectors are proportional. Every determinant (6)
at a node above this open is zero, after any phase resolution and
in every degree. Moreover the tame function (7) there is

    tau_i=constant_sign * q^(-2 v_i(Delta)).         (11)

Here q is a regular unit. These are identities in the actual functions,
not an assumption about the shape of an ambient resolution graph.

It follows that, ABOVE ell, nonzero symbol graph residues can occur
only at the 14 right points for O or 32 for I constituting the two
root orbits in Section 2. On the quotient these are just q=0 and
q=infinity. This is a uniform localization of the CLUTCHING SYMBOL's
double residues, independent of phase degree, coefficient field
extensions used to label the graph, or the graph's number of cycles.

It is NOT yet localization of xi_boundary. The leading-root theorem
localizes original volume POLES. Section 4 localizes a PARTICULAR
logarithmic transition's graph residues inside that region. Neither
assertion by itself proves that other terms in the actual relative
source or its gluing can be discarded.

## 5. Actual pure SAFE controls separate the cut graph from the intrinsic source

This distinction has a sharp original exceptional control. Take

    f=F_n(x,y), n=8 or 12.                            (12)

This is an already closed one-weight PURE SAFE phase. Every fixed
original marker vanishes for all sufficiently large powers by weight
selection. No new safe candidate family is claimed.

For s!=0 its original fibre F_s is an affine-line torsor over the
smooth affine curve {F_n(x,y)=s}/Gamma. Smoothness follows from
homogeneity and Euler's identity, and the finite action is free.
An affine curve has no H^2 and the affine-line fibre is contractible,
so

    H^2(F_s;Q)=0.                                    (13)

The full complement source is therefore just the ambient Haar line;
the intrinsic extension and intrinsic graph obstruction are ZERO.

In contrast T restricts NONTRIVIALLY to F_s^*. Let

    Z={bottom row: F G!=0}/Gamma.

The original projection X_(FG)->Z is the principal additive torsor
of the retained quotient theorem; its regular F-slice trivializes
it. In particular H^2(Z)->H^2(X_(FG)) is an isomorphism. The form
dlog F wedge dlog G comes from Z, and is nonzero there because its
Mayer--Vietoris connecting image in Section 1 is [Omega]!=0.

On each additive fibre the polynomial F_n(x+t z,y+t w)-s has degree
n and leading coefficient F_n(z,w), invertible on Z. Therefore

    F_s^* -> Z

is FINITE FLAT of degree n. Both source and target are smooth surfaces.
The proper degree trace satisfies trace pullback=n: define pushforward
by Poincare duality and proper Borel--Moore pushforward, then use the
projection formula and the degree of the map. Thus the pullback on
rational cohomology is injective, even at ramification.
It follows that [T|_(F_s^*)]!=0.

This class has weight 4. To check that it really gives nonzero graph
residues, upstairs Z is a finite quotient of the nonzero bottom rows
minus the finitely many F,G root lines. This arrangement complement
is a C* bundle over a punctured P1; the line bundle is trivial there.
Its H^2 is pure Tate of weight 4, as is the invariant quotient. The
injective finite pullback preserves this weight. Logarithmic comparison
therefore makes the graph residue vector (6) NONZERO on the cut fibre.

Thus actual original PURE SAFE phases have a nonzero clutching-symbol
graph class on F_s^*, although H^2 of the uncut F_s is ZERO. The added
cuts are essential to that class. Its nonvanishing cannot be identified
with a nonzero intrinsic failure vector.

## 6. The unresolved comparison is now explicit

The whole original unit has the fixed rational Betti/de Rham clutching
representative (2), the quotient coordinate (3), the integer node
determinants (6), and the actual tame functions (7). Compute their
image through the FULL relative support/Mayer--Vietoris source,
remove the prescribed compact-middle part, and make the actual
abelian identity-component correction. Only THEN can the remaining
cycle valuations be identified with xi_boundary,0.

One must retain the cut-graph kernel demonstrated by (12)--(13), the
parameter-dependent units of (7), all original markers, and all
primitive branches. Equality of (6) with the physical Haar variation
c_0, or with the corrected intrinsic cochain, has NOT been proved and
is not asserted. Ordinary marked resolvent algebraicity likewise does
not follow from a logarithmic representative of the ambient unit.

The positive target is still zero intrinsic corrected coupling on
the ENTIRE pure locus; the negative target is an actual invariant
pure phase with a fixed original bad marker. Capacity, automorphism,
primitive and completed family-transfer inputs remain intact. Both
FULL cases remain OPEN. No finite collection of degrees or additional
decidability result has been substituted for this obligation.
