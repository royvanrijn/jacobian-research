# Alternating generic Keller monodromy and compiler obstructions

The polynomial profile compilers that realize alternating arithmetic fibres
have symmetric generic monodromy. This is a uniform obstruction to this
construction method. The normalization-unit obstruction excludes the original
Laurent extension, including its stabilizations. Neither argument excludes
every polynomial plane Keller extension with alternating geometric group.

The [fixed-map theorem](FIXED_KELLER_MAP_ALTERNATING_FAMILY.md) now supplies
one polynomial map with infinitely many A10 whole fibres and generic S10.
That settles the fixed-fibre-family question while leaving the generic
alternating realization question separate.

## Proposition 1. Polynomial coefficient profiles retain symmetric monodromy

Consider the quadratic-gauge compiler over a characteristic-zero field k,
allowing arbitrary polynomial coefficient profiles in its first output V.
Assume its generic inverse polynomial has degree N>=3 and, after normalization,
has the form

\[
 h_V(S)+YS^2-Z/2,\qquad h_V'(0)\ne0,
\]

with Y and Z independent target coordinates. Then both its arithmetic and
geometric generic groups are S_N. In particular no choice of these profiles
produces generic A_N for N>=3.

Proof. Apply the retained
[quadratic-tilt theorem UAQTILT1](../universal-keller-structure/PROOF.md)
over k(V). Its hypotheses are exactly the indicated degree, characteristic
and nonzero linear coefficient, with two independent parameters Y,Z.
It gives S_N before and after algebraic extension of the coefficient field.
The [fixed-map construction](FIXED_KELLER_MAP_ALTERNATING_FAMILY.md) has
h_V'(0)=1, so it lies in this class. Specializing Y,Z to a target curve
can reduce the fibre group without changing the ambient generic group.

## Proposition 2. Retained obstructions do not become a universal exclusion

The [Laurent feasibility theorem](ALTERNATING_LAURENT_FEASIBILITY.md) excludes
polynomial realization of the fixed extension defined by
s=(H(z)+a)/z, a!=0, with independent stabilized coordinates and the indicated
base coordinate retained. Its integral closure contains the nonconstant
unit z. This obstruction survives stabilization of that extension.

It does not apply to a different extension with the same abstract alternating
group. Similarly, the rational one-coordinate factor theorem there excludes
Keller constructions retaining a nontrivial rational factor F1=f(z); it
does not classify all Keller maps with alternating monodromy.

The smallest natural alternating actions illustrate the difference. The
regular degree-three action of A3 fails the retained point-stabilizer
normal-generation test: its stabilizer is trivial. In the degree-four action
of A4, a point stabilizer C3 normally generates A4, and three-cycle inertia
has a fixed sheet and normally generates A4. Thus these elementary tests
alone do not exclude this action. The degree-five A5 action likewise passes:
its stabilizer A4 normally generates the simple group A5, and the
double-transposition class has a fixed sheet and generates A5. Within the
original knot passport the first degree is eight; within the literal rational
square-pencil family it is ten. Passing these necessary tests is not a
realization theorem in any of these degrees.

Proof. The extension-specific assertions are precisely the two retained
normalization arguments. For the finite-group statements, conjugates of a
point-stabilizer three-cycle generate A4, and A5 is simple. Double
transpositions form a nonidentity conjugacy class in A5; their generated
subgroup is normal and hence all A5. The fixed-sheet assertions follow
directly from the displayed cycle types. These statements only test the
specified elementary conditions.

## Theorem 3. A primitive-monomial criterion cannot exclude odd prime degrees

Let k have characteristic zero and let ell be any odd prime. There are
algebraically independent u,v in k[x,y] such that

\[
 [k(x,y):k(u,v)]=\ell,
\]

and every nonconstant monomial x^i y^j, i,j>=0, generates k(x,y) over
k(u,v). The additional element x+y also generates the extension. One
uniform choice is

\[
 u=x+y^{\ell-1},\qquad v=xy+2y^\ell=y^\ell+uy.
\]

These polynomial maps are not Keller maps: their determinant is
x+(ell+1)y^(ell-1), which is nonconstant.

Proof. Work in k[u,y]=k[x,y]. The element v=y^ell+uy is transcendental
over k(u), so u,v are algebraically independent. The relation
y^ell+uy-v=0 is irreducible in k(u,v)[y]: in k[u,v,y] it is primitive
and linear in v with unit coefficient, and its quotient is the integral
domain k[u,y]. Gauss's lemma applies and proves the degree ell assertion.

Put R=k(u,v). The ring k[u,y] is integral over the normal polynomial
ring k[u,v], so R intersect k[u,y]=k[u,v]. Suppose a nonconstant
monomial belongs to R, and therefore

\[
 (u-y^{\ell-1})^i y^j=F(u,v),\qquad F\in k[u,v].
\]

Give y,u,v weights 1,ell-1,ell. The left side is homogeneous of weight
d=(ell-1)i+j>0. The injective graded substitution v=y^ell+uy shows
that F is homogeneous of the same weight. Setting u=0 gives
F(0,y^ell)=(-1)^i y^d, so d=ell r for a positive integer r, and the
u^0 term of F is (-1)^i v^r. Every other monomial u^a v^b of F has
(ell-1)a+ell b=ell r. Hence a is a positive multiple of ell. In
particular no such term contributes to the coefficient of u^1 after
substitution for v. Comparing that coefficient gives

\[
 (-1)^i r y^{d-\ell+1}
       =(-1)^i(-i)y^{d-\ell+1}.
\]

Characteristic zero implies r=-i, contrary to r>0 and i>=0. Thus
every nonconstant monomial lies outside R. By the prime degree, every
one is a primitive element.

Also x+y=u-y^(ell-1)+y has y-degree ell-1. Any polynomial F(u,v)
with positive v-degree has y-degree a positive multiple of ell under
v=y^ell+uy; a polynomial with v-degree zero has y-degree zero.
Thus x+y is outside R and also primitive. Direct differentiation gives
the stated nonconstant determinant. This proves a counterexample to the
general field criterion, with no counterexample to a Keller assertion.

## Source boundary and remaining question

Moskowicz's [arXiv v1 preprint](https://arxiv.org/html/2407.13795v1),
Section 2.8, Answer 2.21, is invoked in Theorem 3.2 as an assertion that
the primitive-monomial hypothesis forces degree two. Section 2.8 expressly
permits a nonconstant Jacobian in that general auxiliary question.
Theorem 3 above disproves that general degree-two forcing
assertion. It does not disprove the preprint's Keller conclusion, and
does not audit every other argument or version. The asserted auxiliary
criterion cannot be imported here as a checked proof premise.

The [source ledger](alternating_generic_keller_source_ledger.json)
pins that exact version and the inspected passages. The
[exact checker](verify_alternating_generic_keller_audit.py) controls the
uniform formulas and finite small-prime instances; the all-exponent
and all-prime argument is the proof above. Its
[receipt](alternating_generic_keller_audit_verification.json) records
internal replay, without formal verification or external review.

The remaining target is either a polynomial Keller extension outside
the retained quadratic-tilt and Laurent-factor architectures with generic
alternating geometric group, or an obstruction applying to every such
plane extension. Neither has been proved here.
