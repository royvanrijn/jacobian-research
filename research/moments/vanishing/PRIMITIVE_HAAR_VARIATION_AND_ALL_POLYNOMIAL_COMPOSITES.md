# Primitive Haar variation and all polynomial composites

16 November 2026. Internal written proof; independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**. Review remains deferred.

There is no need to assume that purity reflects through a polynomial.
For an eventually pure composite H(f), all-original-marker safety is
equivalent to the vanishing of specified primitive Haar-cycle variations
over the roots of H. The primitive map f need not itself be pure, tame,
or have connected generic fibres. All recomposition degrees are covered.

A uniform consequence classifies EVERY eventually pure polynomial
composite of a fixed f under one geometric condition on its actual Haar
monodromy orbit. Either f has a unique safe centre c, and those composites
are exactly H(c)=0 and are all safe, or none exist. This condition has NOT
been established for every exceptional primitive map.

See the [checker](verify_haar_recomposition_cycle.py),
[receipt](haar_recomposition_cycle_verification.json) and
[source audit](haar_recomposition_cycle_source_audit.json).

## 1. The actual primitive cycles, including all branches

Let X=SL2/Gamma for any finite Gamma in SL2(C), with compact core
M=SU2/Gamma, normalized Haar h and descended volume Omega. In particular
this includes 2O, 2I and the 2T cover required in the octahedral sign case.
The integral facts H2(X)=H4(X)=0 and H3(X)=Z give the retained
[tube sequence](HAAR_CYCLE_VARIATION_AND_FIBRE_BOUNDARY_LINKING.md).

Fix a nonconstant regular f:X->A1_s, a nonconstant polynomial H of degree d,
and P=H(f). Let B_f contain the FULL finite bifurcation set of fibres and
complements, including those caused by infinity. Enlarge it when necessary.
Remove H(B_f) and all critical values of H from the r-plane. Over this open
set, H(s)=r has d distinct roots b_i(r), and

    F^P_r = disjoint_union_i F^f_(b_i(r)),
    U^P_r = X minus F^P_r subset U^f_(b_i(r)).          (1)

Choose r_infty sufficiently large that ALL roots are outside f(M).
The actual compact cycle M lifts to U^P_(r_infty). Continue it to a small
punctured disc around r=0 along ONE specified path in the remaining base.
Write gamma_P for this lift. The inclusions in (1) send it to cycles gamma_i
in the individual primitive complements. These are precisely the physical
primitive Haar cycles continued along the lifted root paths. They are not
arbitrary splittings of the primitive tube extensions.

Write H^-1(0)={a}, let e_a be the multiplicity of a, and label the nearby
roots b_(a,j), 0<=j<e_a. Choose any positive common multiple

    L of all e_a,             k_a=L/e_a.               (2)

One positive r-loop permutes the e_a roots over a cyclically. Its L-th
power fixes every label; the path of b_(a,j) is k_a positive primitive
loops about a. This follows in a local coordinate H(s)=z^e_a.
The small disc is chosen to avoid all other removed values.

## 2. Exact transfer of the integral variation

Let A be the one-loop monodromy on H2(F^P_r;Z), and let d_P be defined by

    (T^P_U-I)gamma_P=Tube(d_P).                         (3)

For the primitive complement at b_(a,j), let T_(a,j) be the one-loop
monodromy on H2(F^f_(b_(a,j));Z), and define d_(a,j) by

    (T^f_(U,a,j)-I)gamma_(a,j)=Tube(d_(a,j)).           (4)

No purity assumption is made in this section. All these cycles are
integral and uniquely determined by their specified lifts.
Under the direct sum (1), the exact formula is

    [sum_(v=0)^(L-1) A^v d_P]_(a,j)
       =sum_(v=0)^(k_a-1) T_(a,j)^v d_(a,j).           (5)

To prove it, apply each inclusion in (1) to (T^P_U^L-I)gamma_P.
A tube about a DIFFERENT fibre becomes zero: its normal circle bounds
its disc in the complement with only the selected fibre removed. A tube
about the selected fibre remains its primitive tube. Complex multiplication
by H'(b_i) preserves the normal orientation, so there is no extra integer
factor or sign. The individual inclusions therefore detect every summand
of the direct sum of tubes. Naturality of transport identifies the
projected motion with (T^f_(U,a,j)^k_a-I)gamma_(a,j). Telescoping gives (5).

Equivalently A^L acts on the selected summand as T_(a,j)^k_a. The cyclic
induced block of A over a has e_a summands, and A^e_a acts as a conjugate
of T_(a,j) on each. This retains the branch permutation; replacing the
composite fibre by one primitive fibre would lose it.

## 3. A complete primitive criterion on the eventually pure locus

Over Q, let pi_(a,j,1) be the generalized eigenvalue-1 projector of
T_(a,j), and put

    epsilon_(a,j)=pi_(a,j,1)d_(a,j).                   (6)

**Primitive-variation theorem.** If P=H(f) is eventually pure, then

    P is safe for EVERY original fixed marker
       iff epsilon_(a,j)=0 for EVERY a,j.              (7)

Here eventual purity means h(P^m)=0 for all sufficiently large m.
There is no assertion that f-a is pure or eventually pure.

If P is safe, its actual cycle is fixed by the retained cycle criterion,
so the left side of (5) is zero. Project its right side to eigenvalue 1.
The polynomial S_k(T)=1+T+...+T^(k-1) is invertible on the generalized
eigenvalue-1 space in characteristic zero: its value at 1 is k. Thus every
epsilon_(a,j) is zero.

Conversely suppose all (6) vanish. Let S be the direct sum of the primitive
generalized eigenvalue-1 spaces in (1). It is A-stable, since transport
between successive branches conjugates their full-return monodromies.
On a cyclic block, every generalized eigenvalue-1 vector of A belongs to S:
its full-return eigenvalue is 1^e_a=1. Applying the projection to S in (5)
therefore makes the generalized eigenvalue-1 part of its LEFT side zero.
But that part is S_L(A)pi_(P,1)d_P. Again S_L(A) is invertible there, so
pi_(P,1)d_P=0. The retained eventually-pure cycle theorem now gives d_P=0
and all-original-marker safety. This proves (7).

The converse uses eventual purity and the intrinsic rank-one quotient at
exponential time zero. Without that hypothesis a nontrivial root-of-unity
eigenvalue of A can have a unipotent block. After e_a turns it produces
primitive eigenvalue-1 motion although pi_(P,1)d_P=0. The checker retains
an exact four-dimensional cyclic-block countercontrol. The statement is
not an unconditional equality between these spectral projections.

The criterion is valid for ANY one compatible continuation from infinity.
If P is safe all epsilon vanish for every such continuation. An eventually
pure failure has a nonzero epsilon on every such choice. Do not discard
branches, choose a convenient fixed lift, or independently replace the
specified primitive continuation at each root.

## 4. The exact primitive logarithm interpretation

The condition epsilon_(a,j)=0 is equivalent to absence of INTEGER-exponent
logarithms at s=a in ALL original primitive Cauchy periods

    C_Q^f(s)=integral_(gamma_(a,j)) Q Omega/(s-f),
                       Q in C[X],                     (8)

with the chosen continuation. This equivalence does not require primitive
purity. Indeed the generalized eigenvalue-1 part of the complement cycle
has motion Tube(epsilon_(a,j)). Regular-singular normal form identifies its
nonconstant unipotent part with the integer-class logarithms in (8).
If all original periods have none, the same holds for their derivatives,
which supply every higher pole. Original rational top forms and those
higher poles span H3_dR of the primitive complement. The perfect de
Rham--Betti pairing therefore detects and kills precisely this projected
motion. The reverse implication follows from a fixed eigenvalue-1 cycle.
This is the same spanning argument as in the retained cycle theorem,
applied without its additional pure-to-safe step.

Thus a pure composite failure requires an ACTUAL primitive original marked
integer logarithm at some root of H, on an induced branch. Neither an
abstract bad source lattice, a large invariant boundary space, a generic
fibre eigenvalue alone, nor a noninteger logarithm supplies this condition.
Pure composite safety, conversely, forces all these primitive marked
integer logarithms to vanish even when the primitive unit is not pure.

## 5. A uniform theorem for EVERY polynomial composite of a fixed f

Fix a physical primitive basepoint b outside f(M) and B_f. Inside
H3(U^f_b;Q) form the actual Haar monodromy orbit span

    W_f=span_Q{rho(g)gamma_b : g in pi1(C minus B_f,b)}. (9)

This is a finite-dimensional monodromy-stable space. Transport it to each
finite puncture a. Assume that local complement monodromy on W_f is
SEMISIMPLE on its generalized eigenvalue-1 part at EVERY finite a in B_f.
Other eigenvalues may have unipotent blocks. Points outside B_f have
trivial local monodromy. This hypothesis concerns the actual Haar orbit,
not the entire source and not merely H2 of the primitive fibre.

Every cycle gamma_(a,j) induced by ANY H and any root path belongs to the
corresponding transported W_f. The hypothesis consequently kills (6).
Equation (7) and the retained
[safety reflection theorem](POLYNOMIAL_RECOMPOSITION_REFLECTION_AND_UNIVERSAL_FAILURE.md)
give the following COMPLETE classification, with NO degree bound on H:

* If there is no scalar c for which f-c is all-marker safe, NO nonconstant
  polynomial H(f) is eventually pure.
* Otherwise that centre c is unique. The eventually pure nonconstant
  composites are EXACTLY the H(f) with H(c)=0, and ALL of them are safe
  for EVERY original marker.

For completeness, necessity follows because an eventually pure H(f) is
safe by (7), and reflection gives a safe centre c with H(c)=0. Conversely
if f-c is safe and H(c)=0, expansion of H(c+z)^m has no terms of degree
less than m, so it preserves every eventual marked cutoff. No purity
reflection is used in either direction.

Even the remaining distinction between pure and eventually pure is a
finite UNIT-jet condition in this alternative. Choose N with
h((f-c)^k)=0 for k>N, and set a_k=h((f-c)^k), a_0=1. For H(c)=0,

    h(H(f)^m)=sum_(k=0)^N a_k [z^k]H(c+z)^m.          (10)

Every m>N vanishes. Thus purity is equivalent to (10) being zero for
1<=m<=N, independent of the degree of H. No universal N, no uniform
original-marker cutoff, and no implementation for producing W_f or these
jets is asserted. When N=0 all H(c)=0 are pure. A safe centre alone does
not silently identify eventual purity with purity.

## 6. Controls, descent and the complete remaining target

Actual one-weight invariant seeds f on either quotient are safe with
c=0 and N=0. Their primitive original Cauchy periods are finite Laurent
polynomials; the full spanning argument makes the physical Haar cycle
globally fixed. Hence W_f is the trivial line, and every zero-constant
H(f) satisfies the theorem in every degree. This is a control in an already
closed sector, not a new exceptional sector claim.

The genuine trace control f=x+w on SL2 has H(f)=f^2 with Catalan unit
moments and is not pure. At the only root 0 of H the primitive family
extends smoothly and topologically trivially, so (6) vanishes. This
shows why the variation criterion must be applied WITH eventual purity.
The [finite-lattice theorem](FINITE_RECOMPOSITION_OF_ORIGINAL_FORM_LATTICES.md)
already gives the stronger no-eventual-purity result for every composite
of this tame primitive map.

For an exceptional P, use its canonical regular primitive factor. On 2I
it is invariant. On 2O it is invariant or sign-semi-invariant; in the latter
case work on SL2/2T, where it is regular and H is even. The exact finite
cover marker trace makes safety upstairs and downstairs equivalent for
this pulled-back P. The primitive defect/tameness theorem remains intact.

This gives an all-degree way to control composites without assuming that
their primitive factors are pure. To finish a WHOLE exceptional quotient
by this route it would suffice to prove the W_f condition for EVERY
canonical primitive factor, or directly prove the branch conditions (6)
for every eventually pure H(f). Neither universal statement is proved.
Conversely any exceptional failure must violate that W_f condition at a
finite primitive puncture and exhibit actual marked integer-log coupling
on a compatible branch. Both full decisions remain OPEN. No expanding
degree/support/prefix campaign or external-review gate is introduced.
