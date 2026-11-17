# Minimal-degree pseudo-plane compression and the residual-ramification dichotomy

16 November 2026. PF2RULED128. Written all-degree deductions with explicit
published inputs and small exact controls; no novelty or independent-review
claim. The full Wright/JC2 and marked trace/energy implications remain OPEN.

## 1. Actual setup and two different kinds of ramification

Let F:X=A2->Y=A2 be an actual complex polynomial Keller map of generic
degree d>1. Write N for the full finite normalization over Y. For a
GENUINE block factorization write

    N --h--> M --p--> Y,    deg(h)=t, deg(p)=b>1, d=bt.

Put V equal to the ENTIRE p-etale locus and Z=h^-1(V). The retained
HC4MYKU0/1 construction, applied to the intermediate normalization embedded
in C[X], makes V smooth AFFINE. It contains h(X). Z is normal affine,
finite of degree t over V, and contains the WHOLE original X openly.
The complement Z\X is contained in the original omitted prime curves.

Every loop in Z can avoid its isolated singularities and these complex
curves, and then contracts in X. Thus pi1(Z)=1. This uses the original
plane; it does not assert that M contains a plane open.

Remove the branch locus of h|Z from V and take its inverse image in Z.
The resulting connected covering has degree t. Its subgroup H0 of index
t maps trivially to pi1(V), because the map factors through Z. The map
from the complement fundamental group onto pi1(V) is surjective. Its
kernel contains H0, so

    pi1(V) is finite, m=|pi1(V)| divides t.                 (1)

This argument also works if Z is singular. If h|Z is etale, Z is the
connected universal cover of V, hence is Galois and m=t. This conclusion
does NOT say that h is etale on all N: it may ramify over the critical
curves of p, which were removed from V.

In PF122 notation let sigma_j be an original branch meridian, with
support w_j on bt letters and support a_j on b blocks. Set

    delta_j=w_j-t*a_j,    Delta=sum_j delta_j.

Include ALL original branch supports, including those invisible on
blocks. A moved block contributes all t letters, so delta_j is exactly
the number of moved letters in fixed blocks. Consequently

    Delta=0  <=>  h|Z is etale.                            (2)

Indeed a fixed block is a generically unramified local sheet of p; any
movement inside it is precisely relative ramification over V. Every
divisorial branch of h|Z gives such movement. Conversely their absence
and purity over smooth V make h|Z etale, including over codimension two.

Having no invisible TARGET component alone is insufficient in general:
visible inertia can still move letters inside fixed blocks. In the
binary affine PF123 case the visible delta_C is already zero, so (2)
does reduce to absence of invisible original critical branch supports.

## 2. A proper primitive quotient detects the FULL generic deck group

Let G be original monodromy, H its original sheet stabilizer and J the
stabilizer of a block containing that sheet. A primitive block quotient
means J is maximal proper in G, with [J:H]=t and [G:J]=b. Write
D=Aut_C(Y)(C(X))=N_G(H)/H for the GENERIC rational deck group.

Suppose Delta=0. Section 1 makes C(X)/C(M) Galois, so H is normal in J.
Therefore J is contained in N_G(H). Maximality gives N_G(H)=J or G.
The latter makes H normal in the faithful original action, hence makes
the original extension Galois. The retained KMFI1/HC4MYIB0 obstruction
excludes this for a nonautomorphic Keller map. Thus

    Delta=0 => N_G(H)=J, D=J/H, |D|=t.                    (3)

In particular M is the finite normalization in the fixed field of the
FULL generic deck group. Equivalently, a proper primitive quotient with
trivial generic deck group and t>1 requires residual ramification over V.
Without primitivity, J need not be the full normalizer. Without the
Keller Galois obstruction, N_G(H)=G is possible. Both premises matter.

## 3. Universal-cover compression when chi(V)=1

Assume chi(V)=1. By (1), b1(V;Q)=0. A smooth affine complex surface has
the homotopy type of a real two-dimensional CW complex, so chi(V)=1
makes all positive rational Betti numbers zero. Thus V is a smooth
Q-homology plane. Its log Kodaira dimension is negative infinity:
pullback of log pluricanonical forms under the generically finite
dominant X=A2->V injects them into those of A2, which all vanish. This
is the same whole-plane argument already retained in HC4MYLA1.

Published input A: Dubouloz--Palka, arXiv:1701.01425v2, Lemma 2.9,
printed/PDF p6, gives an A1-fibration over A1, all reduced fibres A1,
and fundamental group the free product of the cyclic multiple-fibre
groups. By finiteness in (1), there is at most one multiple fibre. Hence
V is an affine pseudo-plane of type m, including V=A2 when m=1.

Published input B: Miyanishi, Affine pseudo-planes and affine
pseudo-coverings, Oberwolfach Report 19/2005, printed p1111/PDF p39,
paragraph after Definition 2.5, and Lemma 2.6(1)(i) on p1112/PDF p40.
The universal finite etale cyclic cover W->V of degree m contains an
open A2 which maps SURJECTIVELY to V. The statement is explicitly for
all affine pseudo-planes, with no torus-action hypothesis, and credits
Masuda--Miyanishi, Trans.AMS357(2005),4867--4883, Lemma 2.6. We adopt
the stated covering theorem; its classification/construction is not
claimed independently proved. For m=1, V=A2 supplies the same statement.
Both cached primary PDFs and their relevant complete pages were reread
and visually checked for this use. No new Q-factorial classification
from the separately searched 2007 paper is needed.

Take the plane chart from input B and compose with p|V. It gives an
ACTUAL polynomial map

    F':A2 -> W -> V -> A2,    generic degree bm.           (4)

All arrows are etale on their domains, so F' is etale. Its polynomial
Jacobian is nowhere zero over C and therefore is a nonzero constant.
Since b>1, F' is another nonautomorphic plane Keller map. By (1), its
degree divides bt and is at most d; if m<t, it is STRICTLY smaller.
The chart in (4) need not be the original source X.

Now assume F has MINIMAL GENERIC DEGREE among all nonautomorphic plane
Keller maps. Equations (1) and (4) force m=t. The finite map Z->V lifts
to W because pi1(Z)=1: pull back the finite etale cover and take the
component which is a degree-one cover of Z. This is an algebraic lift.
The induced Z->W is finite, of degree t/m=1. Normality makes it an
isomorphism. Therefore

    chi(V)=1 and degree-minimal F
      => Z=W, h|Z universal finite etale cyclic of degree t,
         Delta=0.                                       (5)

This is a minimal-counterexample reduction, not an exclusion for EVERY
higher-degree map. It neither assumes a plane intermediate M nor uses
degree induction on maps whose source is merely a smooth affine surface.

Moreover Z is the full original etale locus
U_N: an etale point of the composite finite map lies over an etale point
of p (use local finite degrees, including at normal singularities), and
the converse holds since h|Z is etale. The relative cyclic Galois group
of order t is a subgroup of the FULL generic deck group D. HC4MYLA0
makes D act freely on U_N and proves that |D| divides chi(U_N)=t.
Therefore |D|=t: NO PRIMITIVITY is needed for this chi1 conclusion.
The intermediate field is the fixed field of the full D, and

    D=C_t,  b0=t-1,  U_N=W.                              (6)

Here b0 counts ALL original omitted UNRAMIFIED primes. The group need
not preserve the chosen original A2 chart or its boundary basis.

## 4. Degree-four one-cusp binary quotient in a minimal counterexample

For PF123--124 binary affine degree b=4r and cusp count kappa=0 or1,
M is smooth, with r disjoint critical A1s, and

    chi(M)=r*(3-kappa), chi(V)=r*(2-kappa).              (7)

The latter subtracts exactly the r critical curves. Hence chi(V)=1
holds exactly for r=1,kappa=1: the DEGREE-FOUR, ONE-CUSP quotient.
This conclusion requires no hypothesis on the oriented end of M.

For this quotient in a degree-minimal counterexample, (5)--(6) give:

* There are NO invisible original critical branch supports. The sole
  original critical image is the visible curve C, and all its critical
  preimages are smooth by PF123. The entire original N is consequently
  smooth: its possible singularities lie on critical primes, and Z is
  smooth as a finite etale cover of V.
* The FULL generic deck group is C_t, and b0=t-1. A trivial deck group
  with t>1, noncyclic full deck group, or any residual ramification is
  incompatible with a MINIMAL counterexample having this quotient.
* Since a centralizer of a transitive permutation group acts freely,
  C_t has four regular orbits on the 4t sheets. Its centralizer is
  C_t wr S4, so original G is contained in C_t wr S4. In particular G
  is SOLVABLE. This is not a theorem excluding solvable monodromy.

There is a further exact local profile, without assuming pi1(M)=1.
The relative cover N->M is cyclic Galois, branched only over the single
critical E. Let k be its common inertia order, k|t. A finite etale
cover of E=A1 is trivial. The t/k critical primes above E each map
isomorphically to E and have original index 2k. Thus

    #critical primes=t/k,  all e_i=2k, f_i=1,
    pi1(M)=C_(t/k).                                    (8)

For the fundamental group, M=V union E; adding the smooth divisor kills
the normal subgroup generated by its meridian. Its image in pi1(V)=C_t
is exactly inertia C_k. The local cyclic smooth-axis model from PF123
justifies the prime/residue assertions. When k=1, N->M would be etale
everywhere; this is included in the formula, not silently excluded by
a simply-connectedness assumption on M.

## 5. Boundaries and the next global gate

No existence of such a finite Keller normalization is asserted. (5)
allows a cyclic universal cover and (8) allows k<t. End recognition,
primitivity of E in Pic(M), an automatic global LND extension, and
inversion for solvable monodromy are NOT proved here. For arbitrary
nonminimal maps, residual ramification can instead trigger the smaller
counterexample (4). The cases chi(V)>1 remain outside this compression.

The independent [quadratic-refinement criterion](../cross-programme-connections/BREAKTHROUGH_SEARCH_TOPOLOGICAL_COUPLING_2026_11_16.tex)
can obstruct a globally generating affine-transvection tuple when a
single quadratic refinement has all the required values. PF126 only
forces the local preferred longitude to be trivial under a common
alternating form; it does NOT supply the global compatibility equations.
No new closure is inferred from that handoff here.

The companion checker verifies source seals, unbounded degree/Euler
identities, small actual permutation groups, and countercontrols for
the hypotheses. Its finite examples support the written arguments;
they are not a proof of the classification or of JC2. The shared
full-source trace/energy and empty-original-fibre gates stay OPEN.
