# Canonical trace-dual reconstruction of the entire etale locus

Claims HC4MYKU0/1/2. The global reconstruction theorem below is proved
completely, for finite locally free generically etale algebras in every
characteristic and dimension. Its plane application constructs a canonical
etale model for EVERY embedded intermediate field. It does not prove
that these models are polynomial planes or resolve HC4/JC2.

## 1. Full reconstruction theorem: HC4MYKU0

Let B be a Noetherian integral domain, F=Frac(B), and S a finite locally
free commutative B-algebra. Assume L=S tensor_B F is finite etale over F;
L may be a product of separable fields. Normality of S is NOT required.
Embed S in L and define its trace-dual module and generated algebra by

    Lambda = {q in L : Tr_(L/F)(q*S) is contained in B},
    A_can = B[Lambda] inside L.

Then Lambda is finite locally free over B, and

    Spec(A_can) -> Spec(S)

is an OPEN IMMERSION whose image is EXACTLY the full source etale locus
of Spec(S)->Spec(B). In particular that entire etale locus is AFFINE.
A_can is a finite-type etale B-algebra; it need not be finite over B.
If B is normal, A_can is normal. The construction is intrinsic and
independent of a primitive element or basis.

Equivalently, A_can is the universal S-algebra for S-algebra maps to rings
T satisfying Omega_(S/B) tensor_S T=0. Under the finite-flat hypotheses
this means that Spec(T) maps into the etale locus. No scalar target
discriminant is inverted in defining this universal algebra.

### Complete proof

The generic trace pairing is nondegenerate. It identifies Lambda with
Hom_B(S,B). Since S is finite locally free, so is Lambda. The trace of
multiplication by any element of S lies in B; consequently S is contained
in Lambda, which contains1 and is an S-module. Thus B[Lambda]=S[Lambda]
is of finite type over B and over S.

The finite-projective dual-basis element gives a canonical tensor

    e = sum_i s_i tensor q_i in S tensor_B Lambda,
    a = sum_i s_i*Tr(q_i*a) for every a in S.          (1)

A projective dual frame suffices; a global free basis is unnecessary.
After extending F to split L, this tensor is the diagonal idempotent
sum_j epsilon_j tensor epsilon_j. Hence generically

    e^2=e,   mu(e)=1,
    (s tensor1)*e=(1 tensor s)*e for all s in S.       (2)

These are the elementary separability-idempotent identities, not a
semisimplicity assertion for a ramified special fibre.

Put T=S tensor_B A_can. It embeds in L tensor_F L: S is B-flat, and
A_can embeds in L; equivalently T is finite locally free over the
B-torsion-free ring A_can. Thus the identities (2) descend to T without
an undetected torsion term. Its idempotent summand eT is EXACTLY A_can,
via multiplication mu:T->A_can, with inverse a->e*(1 tensor a).
Indeed (2) reduces e*(s tensor a) to e*(1 tensor sa), proving surjectivity;
mu(e)=1 proves injectivity.

Therefore the base change S tensor_B A_can has a clopen factor A_can
at the diagonal section. Its relative differentials vanish on that factor:

    Omega_(S/B) tensor_S A_can
      = Omega_(T/A_can) tensor_T A_can = 0.           (3)

Every point of Spec(A_can) consequently maps to the unramified locus of
Spec(S)->Spec(B). Flatness of the ORIGINAL finite map makes it the
etale locus U. This step cannot be replaced by arguing that the image
of Lambda after tensoring contains1: a pullback module can have torsion.
The actual idempotent proves (3).

Conversely, at every point p in U the trace section generates the
relative dualizing module, so Lambda_p=S_p. This is the usual local
perfect-trace statement for an etale factor. Precisely, the Dedekind
complementary module is the relative dualizing module by
[Stacks49.8.2](https://stacks.math.columbia.edu/tag/0BW2), and the unit
locus of its trace section is the etale locus by
[Stacks49.9.6](https://stacks.math.columbia.edu/tag/0BW9). The hypotheses
are finite, flat and generically etale, exactly as here.

Since Lambda/S is a finite S-module, each p has a principal neighbourhood
D(h) in Spec(S) with Lambda_h=S_h. On it (A_can)_h=S_h. These
neighbourhoods cover U, while (3) already puts the entire image in U.
Thus Spec(A_can) is isomorphic to U. This proves the global theorem,
including all codimension-two points. The universal property follows
because finite-flat etaleness is exactly vanishing of the relative
differentials at the image points. Normality, when B is normal, follows
from the normality of an etale B-algebra. No surface classification,
divisor-Cartier assumption or characteristic-zero argument entered.

## 2. A finite algebra construction rather than a source-dependent module

Choose ANY finite B-generating family q_1,...,q_N of Lambda. Then

    A_can=B[q_1,...,q_N].                              (4)

Thus a supplied finite-projective presentation of S gives a finite
canonical construction, by taking its module dual and identifying it
with its generic trace-dual. No unknown polynomial source coordinates
or infinite trace-zero module are inputs.

If a free B-basis s_1,...,s_d is supplied, set

    G_ij=Tr(s_i*s_j),
    q_i=sum_j (G^-1)_ji*s_j.

These q_i form a B-basis of Lambda and (4) uses d generators. The
algorithm computes field trace and a finite matrix inverse. It does NOT
assert that G^-1 has polynomial entries, that Lambda is an algebra, or
that the resulting etale algebra is a polynomial ring. One may instead
use projective dual frames and module presentations globally; no global
freeness theorem is needed for the reconstruction theorem.

## 3. Canonical etale towers for every plane intermediate field: HC4MYKU1

Let K:A2->A2 be an ACTUAL complex polynomial Keller map, B=C[P,Q],
R=C[x,y], F=C(P,Q), L=C(x,y), and let F subset E subset L be ANY
finite embedded intermediate field. Let S_E be the full normalization
of B in E. It is finite and normal of dimension2, hence Cohen--Macaulay.
Miracle flatness makes it finite locally free over the regular surface B;
see [Stacks00R4](https://stacks.math.columbia.edu/tag/00R4). Local
dimensions agree under the finite map, and its fibres have dimension0.
Thus the entire theorem in Section1 applies, in EVERY generic degree.

Write Lambda_E for its trace-dual and A_E=B[Lambda_E]. Then U_E=Spec(A_E)
is its canonical maximal etale open and is smooth affine over C.
Moreover Lambda_E is contained in R. To check this, take any height-one
source prime. Its actual ramification index over the target is1, so the
intermediate valuation tower has index1 at its centre on Spec(S_E).
The local inverse different is regular there. Every element of Lambda_E
therefore has nonnegative source valuation at every source height-one
prime. Normality of R makes it a source polynomial.

We obtain the canonical factorization by ETALE maps

    A2_source -> U_E -> A2_target.

The first map is etale because it is a morphism between etale schemes
over the same target. This constructs an etale model, NOT an affine-plane
model for E. Kulikov's known countercontrol shows why that distinction
must not be erased.

For embedded E_1 subset E_2, trace transitivity gives

    Lambda_(E_1) subset Lambda_(E_2), A_(E_1) subset A_(E_2).

Indeed Tr_(E_2/E_1)(S_(E_2)) is contained in S_(E_1), by integrality and
normality. Pairing with any q in Lambda_(E_1) then has regular target
trace. Thus the U_E form a choice-free tower of smooth affine etale
models. No intermediate field is given a polynomial reconstruction open
merely by this construction.

For E=L, U_L is the FULL maximal etale completion of the actual source:
it fills ALL omitted unramified primes. In the retained plane geometry
it is Z minus its ramified omitted primes. Consequently

    A_L subset B[N] subset R,
    A_L=R iff NO omitted unramified prime exists.      (5)

Here N=ker(Tr:R->F) and B[N] is the KS0 completion filling only the
single-unramified-lost-sheet primes. Those two completions must not be
confused: a shared-image unramified omitted prime can make the first
inclusion in (5) strict even when B[N]=R.

If b_0 is the number of omitted unramified primes, Pic(U_L)=Z^b_0 on
those added curves, by the same source localization argument as KS0.
Thus factoriality of this ACTUAL completion is equivalent to b_0=0 and
equality with the given polynomial source. This is not a plane-recognition
theorem for arbitrary etale surfaces. In particular, when every omitted
prime is ramified, (4) reconstructs the ENTIRE source from the finite
normalization. The one-boundary nonautomorphic case is included.

### Relation to the already retained two-generator row

PF2RULED53 already proves, under its entire simple-pole r2q1 hypotheses,
Lambda=S+x*S for an actual source coordinate x. That module equality and
its symmetric trace pencil are PRIOR results, not new here. The complete
global reconstruction now gives A_L=S[x] exactly, including normality and
all finite target points. It keeps every unramified omitted prime and
does not manufacture a descent of the source translation to B or F.
PF2RULED52 already excludes such descent in a noninvertible map.

This is a concrete finite-algebra compatibility statement relevant to
the remaining Laurent-key/finite-target work; it does not certify an
unproved key-chain necessity or exclude that entire higher-degree row.

## 4. The differential closure is finitely generated: HC4MYKU2

Keep the ACTUAL complex plane Keller source and ANY intermediate E of
Section3. Let I_L be the full semisimple IC socle of its source module R,
and T_L=R/I_L the semisimple curve-supported quotient from KQ0. Let I_E
be the IC extension of the generic permutation subrepresentation E inside
L. This is a canonical semisimple submodule of I_L. It contains B; put
I_(E,0)=ker(Tr:I_E->B) for its augmentation summand.

The complete identities are

    W_B*S_E=I_E,
    A_E=B[W_B*S_E]=image(Sym^2_B I_E->A_E),
    A_E=image(Sym^2_B I_(E,0)->A_E) if [E:F]>=3.      (6)

Thus the algebra generated by ALL target derivatives of the finite
normalization is already the FINITELY GENERATED canonical trace-dual
algebra. It is determined by the finite input S_E. This is an exact
differential-module/algebra bridge, not a proposed infinite derivative
calculation or a new assumption of differential stability of S_E.

Here are the global and local steps. A_E is etale over B, so the target
derivations preserve it and its fraction field E. Inside R, the
intersection A_E intersect I_L is exactly I_E: its generic subspace is E,
and a semisimple full-support IC module is determined by that subspace.
Hence A_E/I_E embeds in T_L and is SEMISIMPLE with no point-supported
summand. This transfers the actual source-layer theorem without assuming
that the intermediate normalization has its own plane open, disjoint
boundary, or only quotient singularities.

At a generic target divisor, split the complete tame algebra. A_E has
regular DVR factors at ALL unramified primes and Laurent field factors at
ALL ramified primes. In a ramified factor t^e=u, the IC lattice consists
of the regular integer character and all localized nontrivial characters.
The finite normal ring S_E lies in this lattice. Therefore the image of
S_E in A_E/I_E vanishes generically on every curve. Each such element
could only generate a point-supported differential submodule, which the
semisimple quotient does not contain. Thus S_E is contained in I_E
GLOBALLY. Now W_B*S_E is a full-generic-rank submodule of semisimple I_E,
and equals I_E. This proves the first equality in (6).

The quadratic image I_E*I_E is W_B-stable and contains I_E because1 is
in I_E. In every ramified factor the IC elements t^-1,t^-(e-1) have
product u^-1, which generates its integer-character boundary quotient.
There is no missing unramified factor in A_E. Semisimplicity therefore
upgrades the generic calculation to I_E*I_E=A_E. Since A_E is a ring,
this proves the second equality in (6), including finite generation of
the full derivative algebra.

For [E:F]>=3, the trace-zero hyperplane has full quadratic generic span
by the retained nondegenerate-trace argument. Its quadratic image hence
contains I_E. The same two fractional poles have trace zero, so they
recover every ramified boundary quotient again. This proves the final
augmentation identity. The bound is essential: in a quadratic split
field the trace-zero line has only a rank-one square. No unsupported
degree-two augmentation reconstruction is asserted.

The first theorem in Section1 does not require characteristic zero or
any polynomial source. Equation(6) DOES use the actual plane Keller
embedding and the entire KQ0 source-layer theorem. The finite algebra
S_E, its differential IC module I_E and its canonical etale algebra A_E
are three DIFFERENT objects. In particular S_E is not itself stable under
the target derivations, and I_E is generally not an algebra.

## 5. Hostile controls and sharp hypothesis boundaries

**Connected ramified cover with a retained sheet.** For
S=B[t]/(t^3+a*t+b), B=C[a,b], the inverse different is
Lambda=(3*t^2+a)^-1*S. Hence

    A_can=S[(3*t^2+a)^-1].

It retains unramified points ABOVE the target discriminant. For example
t=1,a=-3/4,b=-1/4 lies on the cover and has derivative9/4, while
4*a^3+27*b^2=0. Inverting the TARGET discriminant would wrongly remove
it. This source has nonconstant units and is not a whole polynomial
Keller source. The cubic ramification/companion distinction was already
retained in QLD4/PF2OMIT4; this is a trace-dual-ring control of it.

**Finite flat nonnormal orders.** For S_0=B[w]/(w*(w-x)), the theorem
correctly returns B[x^-1] times B[x^-1]. Its full normalization is B times
B, which is etale even at x=0. Therefore applying the construction to
an arbitrary primitive ORDER instead of the full normalization can lose
actual unramified points. The theorem allows orders; its HC4 application
must use the specified full normalization or prove source landing anew.

**Flatness is essential.** For B=C[x,y], I=(x,y), let
S={(f,g) in B times B : f-g in I}. It is B plus I as a B-module, generically
F times F, but its origin fibre has length3 instead of generic rank2.
Its trace-dual is B times B: x*u,y*u in B force u in B, since x,y are
coprime. The generated algebra is therefore B times B, which includes
the two origin points. The actual etale locus of Spec(S)->Spec(B) is
TWO PUNCTURED PLANES, which is nonaffine. The projective dual tensor in
(1) was unavailable, and the theorem would be FALSE without flatness.
Normality in higher dimension alone must not replace finite flatness.

**Characteristic is not hidden.** The separable characteristic-two
Artin--Schreier algebra F2[a,t]/(t^2+t+a) is finite etale and equals its
own trace-dual algebra. It is nontrivial. The reconstruction theorem is
characteristic-free; the later complex inverse conclusion uses affine
simple connectedness and does not extend to this example.

## 6. Full HC4 target and the last unsolved implication

The user target remains a proof or counterexample for FULL HC4, equivalently
JC2 through the retained Tang reproduction. Sections1--4 complete a
global construction theorem, not that target. No external-review request
or external review is claimed.

For the actual full finite normalization the precise missing implication
can now be stated using a FINITE canonical lattice:

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.             (OPEN)

The right side would force Lambda_L subset Lambda_L^vee=S_L, hence
Lambda_L=S_L. The finite trace pairing would be perfect and the finite
normalization etale. Over the complex affine plane it is then trivial,
giving the Keller inverse. In a supplied free basis this is exactly the
unproved regularity of G^-1. Conversely an inverse plainly has it.

This identifies the last implication; it is NOT proved automatic, nor is
the criterion itself reported as a new exclusion. Canonical module
generation is now explicit and complete. What remains is the actual
plane-source constraint on its multiplication, or the full marked
canonical-energy integrality criterion. Merely enlarging the etale algebra,
normalizing a presentation, or importing a coframe/UFD replacement cannot
finish that implication. Quartic realization remains retired by Domrina.

The inverse criterion in this section uses no differential IC/source-layer
claim from Section4. Its finite-etaleness step can also be checked
directly: after Lambda=S, the finite projective trace form of S is
perfect. In every geometric fibre, a nilpotent element is orthogonal to
the whole fibre algebra, since its product with any element is nilpotent
and has multiplication trace zero. Perfection therefore makes each fibre
reduced. A finite reduced algebra over an algebraically closed
characteristic-zero field is a product of that field, so all geometric
fibres are etale. Finite flatness then makes S/B etale. Affine-space
cover triviality and the actual birational Keller denominator-removal
argument give the inverse, as above. The stronger unresolved differential
input is unnecessary for this purely algebraic criterion.
