# Exceptional primitive fibres and low-rank unit obstructions

17 November 2026. Written all-degree exclusion of fibre architectures,
including two Euler-zero cases. Both FULL exceptional Mathieu--Zhao
statements remain **OPEN**. External review is deferred. No independent
replay, formal verification or universal machine certificate is claimed.

No nonconstant invariant phase on SL2/2O or SL2/2I can have an algebraic
two-torus or Gm times A1 as an upstairs geometric generic fibre component.
More generally, its primitive upstairs fibre cannot have rank-two units
that generate its function field. This eliminates the entire automorphic
polynomial-recomposition architecture of xw, including the retained
unit-blind marked-resolvent control (xw)^2. It does not eliminate all
unit-invisible corrected classes or prove their coupling vanishes.

## 1. Keep the primitive base and its stabilizer

Use the retained
[primitive-base character theorem](PRIMITIVE_BASE_CHARACTERS_AND_FIBRE_EULER_CONSTRAINTS.md).
For every nonconstant original invariant phase f on SL2, there are a
primitive regular q, a shift q0, a finite character chi of Gamma and a
polynomial G such that

    gamma(q0)=chi(gamma)q0,
    K=ker chi contains [Gamma,Gamma],
    f=G(q0^e),       e=|image chi|.

The geometric generic fibre F of q0 is integral. K fixes the base and
acts FREELY, hence faithfully, on F. These assertions hold after passing
to the algebraic closure k of C(q0); the action fixes k. All geometric
generic components of f upstairs have this primitive fibre after the
corresponding base extension. The stabilizer K must not be replaced by
Gamma when the base character is nontrivial.

For 2I, K=2I because 2I is perfect. For 2O, K is either 2O or its
commutator 2T. The latter has

    2T = Q8 semidirect C3,       (2T)'=Q8,
    (2T)^ab=C3.

To see the commutator assertion, let the order-three generator cycle
the quaternion generators i,j,k. The commutators with i,j generate Q8,
and the quotient by Q8 is cyclic of order three. Thus 2T has no nontrivial
character to {1,-1}, and its derived subgroup is nonabelian. The same two
properties hold for 2I: its derived subgroup is itself. These are deck
groups on F, not the global boundary-graph monodromy group in the unit
detector. No comparison between those two actions is assumed.

## 2. Rank-two unit fields cannot carry these faithful actions

**Lemma.** Let H be a finite group with no nontrivial character to C2
and with nonabelian derived subgroup. Let F be an integral variety over
an algebraically closed field k of characteristic zero. Suppose

    Lambda=O(F)^*/k^* is free of rank at most two,
    k(F) is generated over k by the invertible regular functions of F.

Then H has no faithful algebraic action on F.

**Proof.** The action on units gives a finite integral representation

    alpha:H -> GL(Lambda).

In rank two its determinant is a character to C2, hence trivial. Average
a positive definite real inner product over this finite representation.
The image is then conjugate to a subgroup of SO(2), so it is abelian.
In ranks zero and one the same conclusion is immediate. Consequently
H' is contained in N=ker alpha.

For g in N and any unit u, write g(u)=c_g(u)u, with c_g(u) in k^*. This
defines a homomorphism from Lambda to k^*. It is independent of the unit
representative since g fixes k. For g,h in N, c_gh=c_g c_h, so

    N -> Hom(Lambda,k^*)

is a group homomorphism. It is injective: if g fixes every unit, it fixes
their generated function field, and therefore is the identity automorphism
of the integral variety. Thus N is abelian. Its subgroup H' would be
abelian, a contradiction. QED.

The hypothesis concerns the ENTIRE unit lattice. Two independent units
alone do not suffice if further units are moved into their expressions.
Likewise a rational surface or a birational torus chart is not enough:
the functions must be regular invertible functions on the entire F and
generate its function field. These distinctions retain boundary data
that a birational presentation could discard.

Applying the lemma to H=2T in the 2O case, or H=2I in the 2I case, proves
the unconditional necessary constraint

    an exceptional primitive upstairs fibre cannot have
    rank <=2 units generating its entire function field.       (1)

For a surface the field-generation hypothesis already requires rank at
least two, so the substantive excluded case has rank exactly two.

In particular an algebraic two-torus is impossible: its units are exactly
the constant multiples of u^a v^b and its field is k(u,v). Equivalently,
every automorphism of a two-torus is monomial with an invertible integral
exponent matrix and a constant translation. The translation kernel is
abelian; the proof above explains why the exceptional stabilizer cannot
fit into it. No classification of finite subgroups of GL2(Z) is needed.

The group obstruction can be made explicit. ANY 2T action on such an F
kills its central involution: its Q8 subgroup maps into the abelian unit
kernel, so [Q8,Q8]={1,-I} acts trivially. The same holds for ANY 2O action
by restriction to 2T. ANY 2I action on such an F is entirely trivial:
perfectness first kills the abelian unit-lattice image, then the remaining
abelian unit-kernel image. Thus the original free right action contradicts
the particular element -I already, not merely an unspecified kernel.

This also excludes any other entire fibre with the stated unit lattice
and generated function field. No torus chart is substituted for such a
fibre, and the conclusion is not inferred just from its rational cohomology.

## 2a. A unit-field degree constraint without field generation

There is a stronger consequence. Retain the ENTIRE unit lattice Lambda
of rank at most two, but now put U=k(O(F)^*) without assuming U=k(F).
For ANY 2T or 2O action on F, the central element -I acts trivially on U.
For ANY 2I action, the whole group acts trivially on U.

Indeed the unit-lattice image of 2T is abelian by the averaged-inner-product
argument, so Q8 acts on each unit by a constant scalar. Such scalar
actions commute, hence [Q8,Q8]={1,-I} fixes every unit. Restrict a 2O
action to 2T. For 2I, perfectness makes its abelian unit-lattice image
trivial. The remaining action by scalars on units is a homomorphism to
an abelian group, so perfectness kills it too. This proof does not require
the action on units to be faithful.

Suppose now trdeg_k U=2. Since the unit lattice is finite rank and F is
a surface, k(F)/U is a finite field extension. Let

    d=[k(F):U].

The original free action is faithful on k(F). In the 2O case its central
subgroup J={1,-I} fixes U; in the 2I case J=2I fixes U. The fixed-field
calculation therefore gives

    d=|J| [k(F)^J:U],
    2 divides d on 2O,       120 divides d on 2I.           (1a)

For completeness, a finite faithful group J of field automorphisms has
[k(F):k(F)^J]=|J|. Each element satisfies its J-orbit polynomial, of
degree at most |J|. In characteristic zero a finite subextension is
separable and simple, so this bounds every finite subextension by |J|
and hence bounds the entire extension. Conversely the |J| distinct
automorphisms are at most the degree of a finite extension. Equality
follows, and the tower law gives (1a).

Thus field generation, d=1, is only the first excluded case. More
generally any primitive upstairs surface whose full unit field has
transcendence degree two must have this divisibility. These are degrees
of the GENERIC UNIT MAP on the geometric generic fibre, not matrix or
Peter--Weyl degrees of the phase. If trdeg U is smaller, no finite degree
d is asserted. A pair of chosen units cannot replace the full unit field.
Nor is a boundary normal-function character automatically a regular
unit on F. No corrected coupling vanishing follows from (1a) alone.

## 3. Free finite actions on Gm times A1 are cyclic

**Lemma.** Every finite group acting freely on Gm times A1 over k is cyclic.

**Proof.** Its coordinate ring is k[t,t^-1,z], whose units are k^*t^Z.
Every automorphism therefore induces a base automorphism

    t -> a*t or a/t,       a in k^*.

If a nonidentity group element fixes this base, it restricts on every
fibre t=b to a finite-order affine automorphism of A1. Such an automorphism
has a fixed point: a nontrivial translation cannot have finite order in
characteristic zero, and an affine map with nonunit multiplier has its
usual fixed point. The identity fixes every point. Either case contradicts
freeness. Thus the finite group embeds in Aut(Gm).

A base inversion t->a/t has a fixed point b with b^2=a. The corresponding
group element restricts to a finite-order affine automorphism of that
fibre, again with a fixed point. Hence inversions are impossible. The
group embeds in k^* acting by multiplication on t, and every finite
subgroup of k^* is cyclic. QED.

The exceptional stabilizers K are noncyclic, so this rules out Gm times
A1 as a primitive upstairs fibre. The argument uses freeness, unlike
the unit-field lemma. It does not assert that every noncyclic group is
excluded from an invariant cylinder recomposition: a noncyclic Gamma
could have cyclic K in other quotient problems. The actual exceptional
character restriction is essential.

Both a two-torus and Gm times A1 have Euler characteristic zero. Thus
these exclusions supply information not provided by the earlier Euler
divisibility condition, which is vacuous for those fibres.

## 4. Entire automorphic recomposition families are excluded

Use the original matrix entries

    [[x,y],[z,w]],       xw-yz=1.

For p=xw and s outside {0,1}, its fibre has the exact presentation

    x,y in k^*,       w=s/x,       z=(s-1)/y.

It is an algebraic two-torus. The original function field is C(p)(x,y), so p is
primitive. For p=xy and s!=0, its fibre is

    x in k^*, z in k,       y=s/x,       w=(1+yz)/x,

namely Gm times A1; the original function field is C(p)(x,z), so it too is primitive.
These displays over C also prove the geometric generic assertions.

Let phi be ANY algebraic automorphism of SL2 and H any nonconstant
polynomial. Primitivity and the fibre isomorphism type are preserved
under composing p with phi. If H(p composed phi) were invariant under
2O or 2I, the primitive character theorem would give the freely acting
noncyclic exceptional stabilizer K on this same primitive fibre. The
preceding lemmas contradict it. Therefore, for BOTH exceptional groups,

    H((xw) composed phi) is never invariant,
    H((xy) composed phi) is never invariant.                 (2)

This holds in EVERY degree and for EVERY automorphism, with no tame
automorphism hypothesis. Affine shifts of the primitive coordinate do
not change the generic fibre argument.

In particular the diagnostic
[(xw)^2 unit-invisible block](RATIONAL_MONODROMY_BLOCKS_OF_THE_UNIT_INVISIBLE_OBSTRUCTION.md)
cannot be converted into an exceptional phase by ANY automorphism and
polynomial recomposition. That diagnostic is already NONPURE. This new
exclusion is about its entire invariant architecture, not a purity claim
for a finite list of examples. The earlier Euler-one obstruction still
separately excludes all automorphic recompositions of Long's small phase
and of a matrix coordinate.

## 5. A corresponding constraint on the original quotient fibre

Put r=q0^e, primitive on the ORIGINAL quotient. Over its geometric generic
fibre F_down, the pullback to SL2 has e components permuted transitively
by Gamma. Each component maps to F_down as a connected finite etale
Galois cover with group K. Consequently

    pi_1^et(F_down) surjects onto K.                        (3)

This is the defining finite-cover property of the geometric etale
fundamental group; no properness of F_down is used. It follows that
F_down cannot have abelian geometric etale fundamental group. On 2I it
cannot even have prosolvable geometric etale fundamental group, since
the quotient 2I maps onto the nonsolvable A5. Formula (3) concerns the
deck group of the actual cover, not boundary-graph monodromy.

This constraint does not force a particular boundary graph, unit residue,
genus, or corrected Neron component. In particular it gives no automatic
way to recover the full fundamental group from the rational Tate graph.

## 6. Source context and the remaining purity obligation

[Arzhantsev--Gaifullin, Math.Nachr.290 (2017),662--671](https://doi.org/10.1002/mana.201600295),
Example2.3, records the usual monomial automorphism structure of an
algebraic torus. Only this elementary comparison is used here; the unit
argument, finite-group obstruction and cylinder fixed-point proof above
are written directly. No classification theorem for rigid affine surfaces
or birational finite-group actions is imported.

The new result is a necessary constraint on actual exceptional phases,
independent of purity. It rules out entire Euler-zero prototypes for
realizing the unit-invisible obstruction. It does NOT show that every
remaining invisible class is zero, that all primitive fibres have the
excluded form, or that the unit sees a nonzero corrected component on
every surviving mixed phase. General rank-two cohomology and general
rational surfaces are not excluded by the stated unit hypotheses.

The full target remains to prove that the actual balanced ruling/Neron
class vanishes on the entire surviving pure locus, or to realize a
nonzero class by an actual pure exceptional phase with ONE fixed original
bad marker. The capacity, automorphism, primitive-character and full
family-transfer inputs are retained. No degree campaign, further
decidability result, ordinary/EGF identification or external review is
claimed. Both FULL cases remain OPEN.
