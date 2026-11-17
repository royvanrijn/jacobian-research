# Quadratic source reconstruction and the published quartic exclusion

Claims HC4MYKR0/1/2. All-degree plane statements except for the explicitly
imported degree-four theorem. Written proofs and primary-source adoption;
neither a universal inverse nor an HC4 counterexample is supplied.

## 1. A necessary correction to the research frontier: HC4MYKR0

A. V. Domrina, *On four-sheeted polynomial mappings of C2. II. The general
case*, Izvestiya: Mathematics 64:1 (2000), 1--33,
[DOI](https://doi.org/10.1070/im2000v064n01ABEH000273), proves that a
polynomial self-map of the complex plane with generic sheet number four
cannot have a nonzero constant Jacobian. The definition and main theorem
are on printed page1 of the [primary full text](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=273&what=fullteng).
There is NO one-boundary, irreducibility, monodromy, polynomial-degree or
genericity restriction in this theorem. In characteristic zero the sheet
number is exactly [C(x,y):C(P,Q)].

The distinction between the two papers is essential. Domrina--Orevkov,
Part I, proves the one-dicritical case (Theorem0.2); Part II handles the
remaining cases. Part II Proposition1.2 reduces its remaining proof to two
dicritical components; that reduction is not an extra hypothesis of its
main theorem. The final paragraph on printed page32 finishes the cases.
The [source audit](sources/quadratic-reconstruction-context.json) records
what was inspected. This is an adoption of a published theorem, NOT an
independent reproduction of its entire splice-diagram proof.

Consequently ALL proposed complex plane Keller realizations of generic
degree four are excluded, including the residual S4 packets behind QCL4,
QLD4, the older endpoint-pairing gates, and the quartic specialization of
KQ2. Their conditional algebra, completed threefold normalization and
local compatibility identities remain valid. They are controls or possible
independent-proof tasks; they are no longer counterexample candidates.
This does not validate any unfinished internal implication chain.

Makar-Limanov's 2021 *On the Newton polyhedron of a Jacobian pair*,
Section8, also treats generic degree at least six for a counterexample as
established, citing the low-sheet literature including Zoladek2008. We did
not inspect Zoladek's actual Theorem6.12 in this turn. We therefore do NOT
register a fresh degree-five exclusion here or substitute a claim about
total polynomial degree, prime degree, or gcd of polynomial degrees.
Full HC4 and JC2 remain open.

## 2. Setup and the square of the trace-zero source

Use [KQ0/1](SOURCE_MODULE_LAYERS_AND_TRACE_ZERO_KERNEL.md): B=C[P,Q],
R=C[x,y], F=Frac(B), L=Frac(R), d=[L:F], W=D_B, and the actual finite
normalization with omitted primes E_i. Write D_j for their distinct target
images, e_i for ramification index and f_i for residue degree. Set

    N = ker(Tr_(L/F):R -> F),
    A = image(Sym^2_B N -> R),
    ell_j = sum_(E_i over D_j) e_i f_i.

A means finite B-linear sums of products ab with a,b in N. It is not
merely the set of individual squares or a symmetric tensor product without
its field multiplication map. The number ell_j counts missing generic
sheets over D_j, including ramification multiplicities and residue degree.

The Leibniz rule and trace compatibility with target derivations show that
A is a W-submodule of R: delta(ab)=(delta a)b+a(delta b), with delta N
contained in N. Thus no further differential closure of A is necessary.

Assume d>=3. Then A contains the ENTIRE semisimple IC socle I of R.
Indeed, put V=ker(Tr:L -> F). The trace pairing on the separable field L
is nondegenerate. If c is perpendicular to all products V V, then
cV is contained in V-perp=F. For c nonzero this is impossible because
dim_F(V)=d-1>1. Hence V V spans L and A_F=L. Now
I/(I intersect A) has generic rank zero. A quotient of the semisimple
full-support module I is a sum of full-support simples, so this quotient
is zero. This proves I contained in A without identifying I with the
finite normal algebra.

## 3. The exact quadratic defect: HC4MYKR1

Let T=R/I and T_j be its direct summand supported on D_j, as in KQ0.
The all-degree formula is

    A/I = direct_sum_(ell_j>=2) T_j,
    R/A = direct_sum_(ell_j=1) IC_(D_j).                 (1)

The first equality is an equality of embedded submodules of T. Thus the
source is recovered as A=N N over B EXACTLY when every nonproperness
image loses at least two sheets. A sole omitted prime automatically has
e>=2 by the retained ramification theorem, so one omitted prime suffices
for R=A, in ANY degree d>=3 and for ANY residue degree or monodromy.
Neither residue-degree-one nor double transitivity is needed.

Here is a local proof including all boundary characters. Pass at the
generic point of D_j to a faithfully flat completed strictly unramified
DVR extension O=k[[u]] with algebraically closed residue field. Tame
normalization splits the base-changed source algebra into retained O
factors and omitted k((t)) factors with t^e=u. Every retained factor is
unramified. At least one retained factor exists by KMFI1. Each omitted
prime gives f_i such factors; hence their degrees sum to ell_j.

For a retained factor the IC lattice is O. In an omitted factor it is

    O*1 plus direct_sum_(1<=a<e) k((u))*t^a.

This is the elementary transverse IC extension: the trivial character
extends as O, while nontrivial Kummer characters are already localized.
The quotient of the omitted factor by this lattice is k((u))/O in its
integer character. This describes the generic pieces of T_j, not a
claim that the full normal algebra is I. Flat base change commutes with
the trace kernel and with the image defining A.

For EACH omitted ramified factor e>=2, put t^-1 and t^-(e-1) in that
factor and zero in the others. Both have zero trace. Their product is
u^-1 in that factor alone, with trace e/u. This generates its entire
integer-character quotient under transverse differentiation.

Next consider EACH omitted unramified factor, when ell_j>=2. There is
another omitted factor of degree h>=1. A trace-zero element can have
u^-1 in the chosen unramified factor and -u^-1/h in the other omitted
factor, zero elsewhere. A second trace-zero element has 1 in the chosen
factor and -1 in one RETAINED O factor, zero elsewhere. Their product
is u^-1 in the chosen omitted factor alone. Again it generates that
factor's quotient. This covers all omitted factors, including unramified
ones sharing an image with ramified or other unramified primes.

If ell_j=1 there is exactly one omitted factor, with e=f=1. All other
factors are retained O factors. Trace zero then forces the omitted
coordinate of every element of N to be regular. Every product in A is
therefore regular there, so A has zero image in T_j.

Finally T is semisimple, and each of its simple factors is the IC
extension of its generic local system along one D_j. A submodule is
determined by those generic subspaces: there is no point-supported
summand in T. The local calculations recover every character over each
eligible D_j after faithful base change and recover none over the other
D_j. They therefore give (1) globally, including singularities and
intersections of target curves. This proves the stated gluing; a generic
rank comparison without semisimplicity would not suffice.

In particular, once the embedded N inside the ACTUAL field L and its
field multiplication are supplied, all covered boundary extensions are
fixed by the quadratic image. An arbitrary abstract extension class is
not extra construction freedom there. This does NOT reconstruct L, N,
the finite normalization or polynomial source generators from monodromy
alone, and does not make their existence automatic.

## 4. The quadratic trace obstruction: HC4MYKR2

Put M=Tr(R). For d>=3 the trace diagram in KQ1 and (1) give

    Tr(A) contains B,
    Tr(A)/B = direct_sum_(ell_j>=2) IC_(D_j)
              inside M/B = direct_sum_j IC_(D_j).        (2)

The trace on every boundary summand is nonzero on its trivial character;
therefore it sees each eligible image. In particular R=A implies
Tr(N N)=M. Equation(2) is an embedded differential-module statement,
not closure of M under ring multiplication or a numerical Gram matrix.

In ALL degrees, including d=1, one obtains the exact inverse criterion

    the plane Keller map is invertible
      iff Tr(ab) lies in B for EVERY a,b in N.           (3)

The forward implication has L=F and N=0. Conversely, a nonautomorphism
has an omitted ramified prime by the retained boundary theorem. Over its
generic image the first local construction above gives trace e/u in the
base-changed quadratic image. If all original quadratic traces lay in B,
their flat base change would lie in O, a contradiction. This argument
does not use d>=3 or automatic regularity of individual first traces.

Criterion(3) is necessary and sufficient, but its regularity condition is
NOT proved automatic from the Keller determinant. Regular holonomic
self-duality does not supply an ordinary B-valued trace pairing. For
example the rank-one nontrivial Kummer IC module for t^2=u contains
t^-1, whose square has trace 2/u. Confusing these dualities would create
a false proof of inversion.

## 5. Controls, limits and the next useful question

The algebraic check uses the split trace-zero hyperplane as a uniform
model: for v_i=e_i-e_d, two distinct v_i,v_j multiply to e_d, and
v_i^2=e_i+e_d. These identities prove full generic quadratic span for
d>=3; d=2 has only the diagonal square line. Local controls check each
of the three cases used above. These are controls of the written proof,
not a bounded-degree proof of it.

For the cyclic near-cover (t,v)->(t^d,v) on C* times C, the trace-zero
characters multiply to the whole Laurent source ring for every d>=3.
The source is not A2. Its finite normalization is a global complete
intersection, already excluded as an actual nontrivial Keller
normalization by KN0/1/2. This is only a control against claiming that
quadratic reconstruction forces polynomial source generators.

Thus the module-extension gate is reduced to ACTUAL embedded field
multiplication, with a precise exceptional quotient for single unramified
missing sheets. What is still needed is a global constraint on this
embedding or on the marked energy traces that forces inversion, or an
actual characteristic-zero polynomial source construction. Another
abstract Ext calculation or quartic realization search does not supply
that constraint. Full HC4/JC2 remains open. Existing energy, native-action,
arithmetic-lifting and rational-recovery obligations are preserved.
