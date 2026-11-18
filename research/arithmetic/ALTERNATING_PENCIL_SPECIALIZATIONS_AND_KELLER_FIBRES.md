# Alternating-pencil specializations and full Keller fibres

The degree-eight pencils give a fixed quadratic discriminant resolvent in every
unbranched specialization. Hilbert irreducibility supplies full-group field
specializations. The retained finite-etale Keller compiler then realizes their
degree-eight root algebras as complete fibres of polynomial maps. The group of
a specialized fibre and the generic monodromy of its compiler map are different
objects.

The two inputs are the [quadratic-field model](DEGREE_EIGHT_ALTERNATING_PENCIL.md)
over K=Q(sqrt(-329)), and the
[rational model](RATIONAL_ALTERNATING_PENCIL_REALIZATION.md) over Q. For either
model, write k for its coefficient field and

\[
 P_s(T)=H(T)-sT+a,\qquad
 \Delta(s)=\operatorname{Disc}_T P_s=-7^7 Q_4(s)^2.
\]

Here H is monic, a is nonzero, and Q_4 has four distinct roots. Put
k^+=k(sqrt(-7)). In the first case this is
Q(sqrt(-329),sqrt(-7)); in the rational case it is Q(sqrt(-7)).
The element -7 is nonsquare in both coefficient fields.

## Proposition 1. Every unbranched splitting field contains the same quadratic extension

For every s_0 in k with Q_4(s_0) nonzero, let M_(s_0) be the splitting field
of P_(s_0) over k. Then

\[
                      k^+\subseteq M_{s_0}.
\]

This assertion does not require P_(s_0) to be irreducible or to have full
symmetric group.

Proof. Order its eight distinct roots r_1,...,r_8 in M_(s_0). The
Vandermonde product delta=product_(i<j)(r_i-r_j) belongs to that field and
has square Delta(s_0). Consequently

\[
             \left(\frac{\delta}{343Q_4(s_0)}\right)^2=-7.
\]

All denominators are nonzero. This constructs the fixed quadratic extension
inside the splitting field. The quotient is independent of root ordering
up to sign.

The sign action on delta is the discriminant character. If the specialized
group is S8, k^+ is its unique quadratic subfield. For a proper subgroup
it is still the specified discriminant resolvent, but uniqueness among all
quadratic subfields is not asserted. The quadratic field belongs to the
Galois closure; it is not thereby a subfield of the degree-eight root
algebra. In fact on an S8 specialization it is not contained in the root
field, whose stabilizer S7 contains odd permutations.

## Theorem 2. Hilbert specializations with S8 and A8

There are infinitely many s_0 in k for which P_(s_0) is irreducible and

\[
 \operatorname{Gal}(M_{s_0}/k)=S_8,\qquad
 \operatorname{Gal}(M_{s_0}/k^+)=A_8.
\]

For the same parameters, P_(s_0) remains irreducible over k^+. Thus the
degree-eight root field remains a degree-eight field after adjoining
sqrt(-7).

Proof. The retained monodromy proofs establish arithmetic S8 over k(s),
geometric A8, and algebraic constant field k^+ in the generic splitting
field. The polynomial is monic and irreducible over k(s): as a polynomial
linear in s, its coefficients -T and H(T)+a are relatively prime because
a is nonzero. Its bad specialization set is exactly the finite branch
set Q_4=0.

Apply Hilbert irreducibility in its full arithmetic-group form. The precise
imported formulation is Zywina, *Hilbert's irreducibility theorem and the larger
sieve*, Theorem 1.1, pages 1--2 of the retained author manuscript dated
30 November 2010. It supplies a Hilbert set of full S8 specializations;
in the height ordering stated there it has density one. No regularity
hypothesis is imposed in that theorem. It is therefore legitimate to apply
it over k even though the generic S8 splitting field has the nontrivial
constant subfield k^+.

For such a specialization, Proposition 1 and the discriminant character
identify Gal(M_(s_0)/k^+) with the sign kernel A8. Its natural action on the
eight roots is transitive, which proves irreducibility over k^+ as well.
Equivalently, the intersection of the S7 root-field stabilizer with A8 is
A7, giving degree eight after this base change.

Over k^+ one may also apply Hilbert irreducibility directly to the regular
A8 generic splitting field, obtaining infinitely many k^+-parameters with
group A8. The preceding argument is stronger for the stated purpose: it
already obtains infinitely many such A8 specializations with s_0 in the
smaller field k.

Only infinitely many parameter values are asserted here. No counting or
pairwise-distinctness theorem for the resulting number fields is used.
The imported theorem and its application are recorded in the
[source ledger](alternating_specialization_source_ledger.json).

## Proposition 3. A completely explicit S8 number field

For the rational model, take s_0=2. Then

\[
 P_2(T)=T^8+T^6+\frac9{28}T^4+\frac{233}{5292}T^2
                         -2T-\frac1{5488}
\]

is irreducible over Q and has splitting-field group S8. Its splitting
field contains Q(sqrt(-7)), and the group after adjoining this field is A8.

Proof. All coefficient denominators are units at 11, 23 and 47. The
following reductions have the indicated irreducible factors:

\[
\begin{aligned}
 P_2\bmod11
 &= (T^3+T^2-3T-4)
    (T^5-T^4+5T^3-4T^2-3),\\
 P_2\bmod23
 &= (T+1)
    (T^7-T^6+2T^5-2T^4-10T^3+10T^2+3T-5),\\
 P_2\bmod47
 &=T^8+T^6+2T^4+10T^2-2T-17.
\end{aligned}
\]

The last polynomial is irreducible, proving irreducibility over Q and
transitivity. Frobenius at 23 supplies a seven-cycle fixing one sheet,
so the transitive group is primitive. Frobenius at 11 has type (3)(5);
its cube is a five-cycle fixing three sheets. The inspected
[Jordan theorem](sources/jones-primitive-cycle-1209.5169v1.pdf), Theorem 1.1,
therefore puts A8 inside the group. The nonsquare discriminant has class -7,
so the group contains an odd permutation and is S8. Alternatively the
eight-cycle at 47 itself supplies an odd permutation.

The [exact checker](verify_alternating_specialization_fibres.py) verifies
every displayed factorization, factor irreducibility, squarefreeness and
the exact discriminant. This is an arithmetic certificate for one
number-field specialization, not an inference from generic inertia.

## Theorem 4. The complete degree-eight algebras become full Keller fibres

For either coefficient field k and every unbranched s_0 in k, the algebra

\[
                         A_{s_0}=k[T]/P_{s_0}(T)
\]

is realized as the entire scheme-theoretic fibre over a k-rational target of
a polynomial map F_(s_0):A3_k -> A3_k with determinant one and generic
degree eight. Its coordinate degrees are at most 50.

For the Hilbert parameters in Theorem 2, that fibre is a connected degree-eight
field scheme with splitting-field group S8 over k and A8 after base change
to k^+. The same construction after scalar extension realizes
A_(s_0) tensor_k k^+ as the entire degree-eight field fibre over k^+.

Proof and explicit compiler. Apply the retained
[characteristic-zero polynomial realization theorem](../../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex),
label thm:polynomial-realization, and its scalar-extension and arithmetic-transfer
results in [the following section](../../unpublished/common-arithmetic-fibers/sections/03-realization.tex).
Their exact hypotheses are separability and degree at least three. These hold
because Q_4(s_0) is nonzero and the degree is eight.

Here is the supplied-presentation construction, so that this consequence does
not rely on a name alone. Choose a translation tau in k with
P_(s_0)'(tau) P_(s_0)'''(tau) nonzero and write

\[
                     P_{s_0}(\tau+S)=\sum_{j=0}^8g_jS^j.
\]

Such a tau exists over the characteristic-zero field k. For these pencils
one can take tau=1 except at the single parameter s_0=H'(1):
the cubic coefficient is 76 in the quadratic-field model and 541/7 in
the rational model, both independent of s_0.

In the new source coordinates (x,y,w), define

\[
 \vartheta=1+xy,\qquad
 q=\vartheta^2w+\frac{g_1}{g_3}y^2(1+3\vartheta),\qquad
 \Pi=\vartheta q,
\]
\[
 B=y+3\frac{g_3}{g_1}xq+2\frac{g_2}{g_1}\vartheta q
       +\sum_{j=4}^8j\frac{g_j}{g_1}\vartheta^2x^{j-2}q^j,
\]
\[
 C=x(5-3\vartheta)-\frac{g_3}{g_1}x^3w
       -\sum_{j=4}^8(j-2)\frac{g_j}{g_1}(xq)^j.
\]

Set

\[
                    F_{s_0}=(\Pi,-B/2,C),\qquad
                    y_{s_0}=(1,0,-2g_0/g_1).
\tag{1}
\]

These are polynomial expressions: their denominators are the nonzero
coefficient scalars g_1 and g_3, not source functions. The coefficient
specialization s_0 is fixed in the map, not retained as a target variable.

For the Jacobian and fibre calculation, introduce the rational chart

\[
 S=x/\vartheta,\quad Q=y+xq,\quad D=1-SQ+\Pi S^2=1/\vartheta,
\]
\[
 G_\Pi(S)=g_1S+\Pi(g_2S^2+g_3S^3)
                         +\sum_{j=4}^8g_j\Pi^jS^j,\qquad
 \beta=\frac{G_\Pi'(S)/g_1-1-\Pi S^2}{S}.
\]

The numerator defining beta is divisible by S. The displayed polynomials
satisfy B=Q+beta and C=2G_Pi(S)/g_1-BS^2. The source chart has Jacobian
D^-1 and the plane map (S,Q) -> (B,C) at fixed Pi has Jacobian -2D.
The factor -1/2 in (1) therefore gives determinant one. The identity holds
on a dense chart and hence globally for these polynomial maps.

Every test-algebra point of the fibre has vartheta*q=1; thus vartheta is a
unit on the whole fibre. The marked S is therefore defined everywhere on
it. The inverse equation is

\[
 G_\Pi(S)-\frac{g_1}{2}(BS^2+C)=0.
\]

At the target in (1), it becomes P_(s_0)(tau+S)=0. Its derivative divided
by g_1 is exactly D. Separability makes this derivative a unit over every
commutative test algebra: evaluate a Bezout identity U P+V P'=1.

Conversely a root S in a test algebra reconstructs a source point. Put
d=P_(s_0)'(tau+S)/g_1, Q=-beta(1,S), and

\[
 \vartheta=d^{-1},\quad x=S/d,\quad y=Q-S,\quad q=d,\qquad
 w=d^2\left(d-\frac{g_1}{g_3}y^2(1+3/d)\right).
\]

The identities above verify the fibre equations, and extraction of S
recovers the same root. Conversely they recover the same source point
from every fibre point. These operations commute with test-algebra
homomorphisms, proving the full scheme isomorphism with Spec A_(s_0),
not only agreement of geometric points.

The rational chart is birational. Generically the inverse equation has degree
eight and is irreducible because it is linear in the independent constant
target coordinate C, with nonzero scalar C coefficient. It therefore gives
the actual degree-eight source extension. The largest polynomial degree in
the displayed B coordinate is 6*8+2=50. These are the degree and fullness
claims of the retained construction.

For the explicit s_0=2 field in Proposition 3 and tau=1, the complete
coefficient list is

\[
 (g_0,\ldots,g_8)=
 \left(\frac{54125}{148176},\frac{35387}{2646},
       \frac{237995}{5292},\frac{541}{7},\frac{2389}{28},
       62,29,8,1\right).
\]

Thus (1) is an explicit rational polynomial map, given without expanded
degree-50 output, whose full fibre over

\[
                         (1,0,-54125/990836)
\]

is precisely Spec Q[T]/P_2. The checker verifies derivative invertibility,
source reconstruction, target equations and extraction of S in the exact
degree-eight quotient, together with the two Jacobian factors.

## The generic-group distinction

The splitting field and Galois action of the finite fibre are preserved
because its algebra is unchanged. This does not prescribe the generic
monodromy of F_(s_0). Its generic inverse polynomial is the three-target
polynomial G_Pi(S)-(g_1/2)(BS^2+C), not the original Laurent pencil
H(T)-sT+a. No generic A8 conclusion for the polynomial Keller map is made
from the specialized A8 fibre.

There is also no conflict with the
[normalization-unit obstruction](../jacobian/monodromy-theory/ALTERNATING_LAURENT_FEASIBILITY.md).
That obstruction excludes polynomial realization of the exact Laurent
extension over an affine base retaining its original cover coordinate s.
Here s_0 is a fixed scalar and a different polynomial map is constructed
from its finite quotient algebra. The finite fibre is the transferred
object; the original function-field extension is not.

## Verification and source scope

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_alternating_specialization_fibres.py
```

The [receipt](alternating_specialization_fibres_verification.json) records
Python 3.12.3 and SymPy 1.14.0, exact modulo-11/23/47 factors, translated
coefficients, the normalized target, Jacobian factors, quotient reconstruction
and invalid-input controls. The Hilbert and scheme-theoretic quantified
conclusions additionally use the written proof and the specified imported
theorems. No formal verification of this new bridge, independent external
replay, external human review or priority is claimed.
