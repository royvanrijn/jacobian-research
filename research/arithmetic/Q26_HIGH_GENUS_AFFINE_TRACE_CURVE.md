# The degree-twenty-six affine-trace curve has genus at least 1051

The explicit accessory curve controlling normalized degree-twenty-six
paired pencils with affine canonical trace is geometrically integral,
has geometric splitting group S12 over the slope line, and has genus
at least1051. It has only finitely many rational points and admits no
rational parametrization. Consequently the normalized rational pencils
in this slice form a finite set.

This is a finiteness theorem without a point list. It neither constructs
a rational degree-twenty-six paired pencil nor excludes every such
pencil. The general non-affine-trace problem remains open. The
previous local Q29 and number-field points remain compatible with
the result.

## 1. The exact accessory equation has no vertical component

Use h6 in the
[uniform odd affine-trace theorem](ODD_AFFINE_CRITICAL_TRACE_CURVES.md).
Its normalized data are

\[
 f(z)=N(z^2)/z,\quad \deg N=13,\quad N(1)=N'(1)=0,
 \quad \theta=s+t\rho.
\]

The canonical scale rho separates the six nonzero critical pairs.
Away from explicitly specified exceptional slopes, the theorem
reconstructs its monic polynomial q and the numerator N rationally
from s,t. The remaining equation is Phi_6(s,t)=N(1)=0, of degree
twelve in s. Every rational good model belongs to that chart: its
only possible rational exceptional slope is t=-2, and that slope
is impossible for a good model.

Let D0(t) be the monic least common denominator in Q[t] of all
thirteen coefficients of Phi_6 as a polynomial in s. The exact
coefficient arithmetic gives

\[
 \deg D_0=144,\qquad
 \gcd_{\mathbb Q[t]}\{[s^i](D_0\Phi_6):0\le i\le12\}=1.
 \tag{1}
\]

The gcd in(1) is computed **before** removing any coefficient content.
Thus no vertical component has been discarded. Clear only the rational
numerical denominators and divide the common integer content, choosing
a positive final leading coefficient, to obtain F in Z[s,t]. Then

\[
 \deg_s F=12,\quad \deg_t F=144,\quad
 \deg_t([s^i]F)=144-i\quad(0\le i\le12).                 \tag{2}
\]

Its coefficient gcd in Q[t] and its integer content are both one.
The [input](q26_high_genus_affine_trace_curve_inputs.json) retains
every coefficient. The checker rebuilds(1)--(2) from the frozen
uniform-curve coefficient input, over Q[t], rather than inferring
primitivity from a reduction modulo a prime.

## 2. The arithmetic splitting group is S12

Specialize t=2, divide the integer coefficient content and choose
positive leading coefficient, obtaining f2 in Z[s]. Its leading
coefficient is

\[
                 a=699925604026943163215025\ne0.
\]

All factors in the following exact reductions are distinct and
irreducible, and their leading coefficients are units:

| Prime | Factor degrees of f2 |
| --- | --- |
| 17 | 12 |
| 47 | 1,11 |
| 71 | 1,2,9 |

The input contains the monic factors and their scalar products.
The checker verifies each product, separability and irreducibility
by the exact finite-field criterion

\[
 X^{p^d}=X\pmod g,\qquad
 \gcd(g,X^{p^{d/\ell}}-X)=1
                  \quad\text{for every prime }\ell\mid d.
 \tag{3}
\]

No factorization routine is used to establish the retained factors.
To apply the integral monic form of
[Milne, *Algebraic Number Theory*3.00, Theorem8.23,p138](sources/milne-ant300.pdf),
use the polynomial a^11*f2(X/a). It is monic integral and has the
same splitting field. Since a is a unit at17,47,71, this change
preserves the displayed factor degrees.

Thus the specialized group contains a12-cycle, an11-cycle, and a
permutation of type(1,2,9). The ninth power of the latter is a single
transposition. The first cycle makes the group transitive. The
stabilizer of the fixed letter of the11-cycle is transitive on the
other eleven letters, so the group is two-transitive. Its conjugates
of a transposition therefore contain every transposition, and it is
S12. This is also the criterion in Milne's Lemma8.26,p139.

The specialization t=2 has nonzero leading coefficient and
discriminant. Its splitting group embeds into the generic splitting
group: the generic splitting cover is unramified there, and the
specialized splitting group is the residue group of a decomposition
group. Consequently the arithmetic splitting group of F over Q(t)
is S12. In particular F is irreducible over Q(t), hence over Q[s,t]
by(1).

The prime29 is used below for an entirely different purpose. The
degree-twelve coefficient of f2 vanishes modulo29, so its reduction
there has degree eleven. No degree-twelve Frobenius inference at29
is made.

## 3. An exact modular certificate gives 2124 simple branch points

Write D(t)=Disc_s F in Z[t], and reduce at29. The full finite-field
polynomial calculation gives

\[
 \deg\overline D=3036,\qquad
 \deg\gcd(\overline D,\overline D')=782.                 \tag{4}
\]

Put G=gcd(Dbar,Dbar') and define

\[
       A=\frac{\overline D/G}{\gcd(\overline D/G,G)}.
 \tag{5}
\]

Up to a nonzero scalar, A is exactly the product of the irreducible
factors of Dbar having multiplicity one. This formula remains valid
when a factor's multiplicity is divisible by29: that factor is then
entirely retained in G and contributes nothing to Dbar/G. The exact
gcd calculations show

\[
 \deg A=2124,\qquad
 \gcd(A,A')=1,\qquad
 \gcd(A,\overline{\operatorname{LC}_s F})=1.              \tag{6}
\]

The entire polynomial Dbar and A are retained. Their binding to F
has two separate exact checks. First, the checker computes the
bivariate discriminant over F29. Second, it independently computes
ordinary univariate discriminants at3169 distinct slope values in

\[
 \mathbb F_{29^3}=\mathbb F_{29}[z]/(z^3+2z+27),          \tag{7}
\]

using a univariate Euclidean resultant algorithm. At each chosen
value the degree in s remains twelve. The extension polynomial is
verified irreducible. Five values with zero leading coefficient
are skipped; the3169 checked values are distinct.

This is a complete polynomial identity certificate. Indeed a
degree-twelve discriminant is homogeneous of degree22 in the
coefficients, so(2) implies the safe bound

\[
                     \deg_t D\le22\cdot144=3168.          \tag{8}
\]

The retained polynomial has degree3036, within this bound. Equality
at3169 distinct values proves equality as polynomials over the
extension field and hence over F29. These are not sample-based
estimates of the discriminant or branch count. No characteristic-zero
discriminant factorization is needed.

**Lemma 1.** The characteristic-zero discriminant D has at least2124
distinct simple geometric roots at which LC_s(F) is nonzero.

**Proof.** Each geometric root of A lies in a finite residue extension
of F29. Over the corresponding finite unramified extension of Q29,
the simple-root Hensel lemma lifts it to a root of D. Its derivative
and LC_s(F) remain units. Distinct residue roots have distinct lifts.
These are algebraic characteristic-zero roots of D and are simple.
The simple-root lemma is
[Milne, Proposition7.31,p113](sources/milne-ant300.pdf). ∎

At a simple discriminant root with unit leading coefficient, the
inertia in the degree-twelve root cover is a transposition. Locally
one may normalize the polynomial to be monic. A ramification cycle
of length e contributes at least e-1 to the root-difference
discriminant valuation, while collisions among wholly unramified
roots contribute an even valuation. Valuation one therefore forces
exactly one2-cycle and all other sheets unramified. Equivalently,
this is the usual simple-discriminant branch calculation in
characteristic zero.

## 4. Geometric S12, high genus and finiteness

**Theorem 2.** The arithmetic and geometric splitting groups of F
are S12. Its smooth projective normalization C is geometrically
integral and satisfies g(C)>=1051. The curve has no nonconstant
rational parametrization. For every fixed number field K, C(K)
and the affine point set F(K)=0 are finite.

**Proof.** The geometric group is normal in the arithmetic S12.
By Lemma1 it contains a transposition. All transpositions are
conjugate in S12 and generate S12, so that normal subgroup is the
whole S12. This also proves geometric irreducibility and makes
C->P1_t a degree-twelve cover.

Each of the2124 simple branch points contributes one to the
ramification of this cover. Riemann--Hurwitz gives

\[
 2g(C)-2=-24+\operatorname{Ram}(C/P^1)
                       \ge-24+2124=2100.
\]

Thus g(C)>=1051. Other branch points may increase this bound;
no exact genus is asserted. Riemann--Hurwitz also excludes a
nonconstant map P1->C, so there is no rational parametrization.
Finiteness of C(K) follows from
[Faltings, *Endlichkeitssätze für abelsche Varietäten über Zahlkörpern*,
Satz7, printedp365](sources/faltings-1983-endlichkeitssaetze.pdf).
Nonsingular affine K-points lift to C(K), and there are only finitely
many singular affine points, proving affine finiteness as well. ∎

**Corollary 3.** There are only finitely many normalized good rational
degree-twenty-six pencils N(z^2)-wz whose canonical trace is affine
in its canonical scale. No nonconstant rationally parametrized family
of these normalized coefficient models exists.

**Proof.** The complete rational-slope statement of the uniform
theorem puts every such model on F(Q)=0, away from the specified
denominators and good-open exclusions. Reconstruction gives a
unique N at each valid point. Theorem2 makes the point set finite.
A rationally parametrized family of normalized coefficient models
would recover a nonconstant rational map from P1 to this curve by
the rational canonical inverse, which is impossible. ∎

The conclusion is about normalized coefficients. Rational coordinate
rescalings or twists can produce further raw coefficient models and
are not counted as new points on this normalized chart. No claim is
made about finiteness of all Q-isomorphism classes without this
normalization. Over a larger fixed number field, F(K)=0 remains
finite, but exceptional slopes outside the triangular chart would
need separate analysis before asserting finiteness of every model
in that field.

The good Q29 point from the uniform theorem and its number-field
lift show that the good locus is nonempty geometrically. Finiteness
does not determine whether its set of Q-rational points is empty.
The non-affine-trace degree-twenty-six paired-pencil problem is
unaffected by this restriction.

The S12 group in this note is the monodromy of the accessory curve
over its slope parameter t. It is distinct from the A26 monodromy
of a reconstructed Laurent cover over its pencil parameter w.
Finiteness concerns normalized model parameters (s,t), not w.
Even one rational cover, if found, could still supply infinitely
many Hilbert specializations and corresponding Keller fibres.

## Exact evidence and scope

The [checker](verify_q26_high_genus_affine_trace_curve.py),
[inputs](q26_high_genus_affine_trace_curve_inputs.json),
[receipt](q26_high_genus_affine_trace_curve_verification.json) and
[source ledger](q26_high_genus_affine_trace_curve_source_ledger.json)
retain all primitive coefficients, modular polynomials, Frobenius
factors, exact arithmetic checks and imported theorem locators.
The frozen uniform operator theorem supplies the paired-pencil
reconstruction; this note proves the arithmetic and genus statement
for its exact accessory curve.

There is no exact genus computation, effective rational-point list,
rational existence or nonexistence theorem, arbitrary-pencil
classification, formal verification, external human review,
independent external replay or novelty claim.
