# A rational composition normal form for every pure odd-start phase

Date: 20 November 2026. This is a characteristic-zero theorem for
the exact two-height profile, not a classification of arbitrary
sphere phases or arbitrary Image polynomials.

Let K be ANY characteristic-zero field, and fix an odd integer a>=1.
Normalize the regular unit-center phase

 p=t^a(C(x)+t²D(x))/x, C,D in K[x], C(0)=-D(0)=1.

Purity means vanishing at EVERY positive exponent. The retained
rationality and endpoint results give rational functions

 X=sqrt(-D/C), X(0)=1,
 u=xX^a/C=x+O(x²), W=X^(a-1)/C,
 T=W/u'=1+O(x).

If v is the selected inverse of u at0, purity is exactly

 T(v(u))+T(v(-u))=2.                                  (N1)

The rationality theorem is an imported premise. The new composition
argument below uses only K(u,T); it introduces no rational primitive,
logarithmic decomposition, or polynomial-quotient hypothesis.

## 1. The exact iff

**Theorem OddStartRationalCompositionNormalForm.** A normalized
polynomial pair C,D is pure if and only if it is obtained as follows.

Choose rational functions w(x),U(z),V(z) OVER K such that

 w(0)=0, w'(0)=1,
 U(0)=0, U'(0)=1, V(0)=0,
 U(-z)=-U(z), V(-z)=-V(z).                            (N2)

All values in(N2) mean regularity at the indicated finite point.
Define rational functions

 u=U(w), T=1+V(w), X=u/(x T u'),
 C=x X^a/u, D=-C X².                                 (N3)

Require the TWO finite polynomiality conditions C,D in K[x].
Then C(0)=-D(0)=1 automatically, and the resulting phase is pure.

No polynomiality of w,U,V or of an outer phase is imposed.
The data form an exhaustive rational parametrization with exact
finite tests for each proposed triple. They do not give a finite
list of all triples or a purity decision from finitely many moments.

If one additionally requires

 K(U(z),V(z))=K(z),                                  (N4)

then w is the canonical generator of the ACTUAL field E=K(u,T).
Under(N4) the triple(w,U,V) is UNIQUE for the normalized pair C,D.
Without(N4), the construction is still sufficient and exhaustive,
but it permits redundant presentations of the same phase.

## 2. The field involution and the marked rational place

Embed K(x) into K((u)) by x=v(u). This is injective: x is
transcendental, u(x) is nonconstant rational, and formal substitution
respects rational identities. Let

 E=K(u,T) subset K(x).

In the formal embedding, replacement u by-u sends

 u -> -u, T -> 2-T

by(N1). Thus it preserves E and defines a nontrivial K-involution
sigma of E. Its square is the identity. This does NOT assert that
sigma extends to K(x).

Actual Lüroth over K gives E=K(f) for some rational f in K(x).
The restriction of the x=0 place is K-rational: evaluating rational
functions there has residue K. Since u in E has x-order1, this
place is unramified in K(x)/E.

Here is an elementary normalization without hidden place theory.
The value f(0) belongs to K union{infinity}. A K-Mobius change makes
it0. Then f has some positive x-order e. Express u=R(f).
Since u has order1, R has positive orderk at0 and e*k=1.
Hence e=k=1. After dividing f by f'(0) in K*, obtain a generator z
with z(0)=0,z'(0)=1.

Formal replacement u->-u preserves order and sends the expansion
z(v(u))=u+O(u²) to-u+O(u²). Consequently sigma fixes the marked
place and has derivative-1 there.

Every K-automorphism of K(z) is Mobius. Indeed the field degree of
a rational map equals its degree; an automorphism and its inverse
therefore both have degree1. Its normalized form here is

 sigma(z)=-z/(1+beta z), beta in K.                    (N5)

This explicitly shows the other fixed point is K-rational: it is
infinity for beta0 and -2/beta otherwise. No algebraic closure or
square-root choice is needed.

**Lemma CanonicalAntiInvariantGenerator.** Put

 w=z/(1+beta z/2).                                    (N6)

Then E=K(w), w(0)=0,w'(0)=1 and sigma(w)=-w. Moreover this w is
the UNIQUE generator with these three properties.

The displayed identities follow by direct substitution in(N5).
For uniqueness another generator is g(w), with g Mobius over K.
It fixes0, is regular there, has derivative1 and is odd. Write
g(z)=z/(1+gamma z). The identity g(-z)=-g(z) forces gamma0 in
characteristic zero, so g=z. QED.

This argument proves descent over the original K and the exact
marked-base normalization; it does not invoke a generic assertion
that a ramification index divides a field degree.

## 3. Necessity and sufficiency

Since u,T belong to K(w), write u=U(w), T=1+V(w).
The expansions at0 give the three normalizations for U,V in(N2).
The involution gives U(-z)=-U(z),V(-z)=-V(z). Because
K(u,T)=K(w), condition(N4) holds. Finally the definitions of u,T
give X=u/(x T u') and therefore(N3). This proves necessity.

Conversely choose any data(N2), without requiring(N4), and suppose
the rational C,D in(N3) are polynomials. We have u=x+O(x²),
T=1+O(x), X=1+O(x), C=1+O(x),D=-1+O(x). Thus the phase is regular,
X is the selected square root, and its reconstructed coordinate and
Jacobian are exactly u,T.

Let eta(s) be the formal inverse of U at0. Since U is odd with
derivative1, eta is odd. Compatibility of formal inverses gives

 w(v(s))=eta(s), T(v(s))=1+V(eta(s)).

The last summand is odd. Thus(N1) holds, proving ALL-ORDER purity
by the endpoint theorem. No rational return on the original x-line
or finite moment cutoff is assumed.

If(N4) holds, E=K(w) and w is the unique generator in Lemma(N6),
so the construction recovers its own canonical triple. This proves
the uniqueness/bijection assertion. If V=0, condition(N4) forces
degU1; oddness and the derivative normalization then give U(z)=z
and canonical w=u.

Odd rational functions may equivalently be written

 U(z)=z A(z²)/B(z²), V(z)=z C0(z²)/D0(z²),

with B(0),D0(0) nonzero and A(0)=B(0). Thus the parity and base
conditions are themselves finite rational coefficient identities.

## 4. Explicit finite polynomiality tests

Write u=A/B,T=M/N as rational functions in x, with polynomial
numerators/denominators and nonzero denominators. Put

 Delta=A'B-AB'.

The base conditions guarantee Delta and M are not identically zero.
Directly from(N3),

 X=A B N/(x M Delta),
 C=A^(a-1) B^(a+1) N^a /
                (x^(a-1) M^a Delta^a),
 D=-A^(a+1) B^(a+3) N^(a+2) /
                (x^(a+1) M^(a+2) Delta^(a+2)).        (N7)

Therefore the exact polynomiality tests are the TWO divisibilities

 x^(a-1) M^a Delta^a divides A^(a-1) B^(a+1) N^a,
 x^(a+1) M^(a+2) Delta^(a+2)
                   divides A^(a+1) B^(a+3) N^(a+2).  (N8)

They are valid with or without cancelling common factors first.
Polynomial division gives finitely many exact remainder equations.
Over a field with effective arithmetic they are finite algorithms;
over an abstract K they are finite algebraic conditions, without
assuming a computable presentation of K.

Primitivity(N4), when desired, is a separate function-field condition.
It is not needed for sufficiency. Here is an exact finite criterion. Write U=P0/Q0,V=R0/S0 in reduced
form. In K(z)[Y], form the two polynomials

 P0(Y)Q0(z)-P0(z)Q0(Y),
 R0(Y)S0(z)-R0(z)S0(Y).

Their monic gcd has degree[K(z):K(U,V)]. Indeed its common roots
are precisely the distinct images of z under embeddings over K(U,V):
the reduced denominators cannot vanish at a common root, and U(z)
is transcendental. Characteristic zero makes the generic fiber
separable. If V is constant the second polynomial is zero and the
first gives the same criterion. Thus(N4) is equivalent to gcd Y-z.
This is a finite Euclidean computation over K(z), not a shortcut
using the gcd of the two rational-map degrees.

## 5. Concrete canonical coordinates and scope controls

For the nontrivial Phi1 examples, E=K(u), so canonical w=u,
U(z)=z,V(z)=0. In particular the two-root family

 C=(1+2x/3)^(a-2)(1+x)^(a+3),
 D=-(1+2x/3)^a(1+x)^(a+5), a>=3 odd,

has canonical w=x(1+2x/3)²/(1+x)³, which is not polynomial.
The former cubic-height polynomial-quotient mechanism cannot be
inserted into this theorem.

A second simple example makes rational OUTER coefficients essential.
Take lambda!=0 and

 w=x/(1+lambda x), U(z)=z, V(z)=-lambda z.

The original pair is polynomial:

 X=(1+lambda x)²,
 C=(1+lambda x)^(2a+1),
 D=-(1+lambda x)^(2a+5).

But the corresponding outer pair, reconstructed from U,V with
source variablez, is

 Cbar=(1-lambda z)^(-a), Dbar=-(1-lambda z)^(-a-2).

It is rational with finite poles, so no polynomial-outer premise is
available even in this elementary canonical example.

For the retained hyperbolic family put H=1+lambda x and beta=n*lambda/2.
Then

 w=2(H^n-1)/(n*lambda*(H^n+1)),
 U(z)=z/(1-beta²z²),
 V(z)=-2beta z/(1+beta²z²).                           (N9)

Direct substitution reproduces its u and T. Here(N4) follows from

 z=2 U(z)V(z)/(V(z)-2beta U(z)).                      (N10)

These are rational identities. Composing w with an admissible
polynomial source map, including the retained degree59 counterexample,
does not alter U,V and preserves their primitivity.

## 6. Polynomial outer phases impose a different restriction

The canonical anti-invariant coordinate need not give polynomial outer
coefficients. The following restriction concerns attempts to use a
polynomial outer pair with the same admissible-pullback formulas.

For clarity, suppose a nontrivial polynomial pure unit-center
outer phase of the same odd height has degrees cbar,ebar, and a
rational map w normalized at0 produces polynomial C,D by the
same-height pullback formulas

 C=Cbar(w)(w/x)^(a-1)/(w')^a,
 D=Dbar(w)(w/x)^(a+1)/(w')^(a+2).

At a finite pole of w of orderk>=1 (necessarily away from0),
the exact local orders are

 ord C=a+k(1-cbar),
 ord D=a+2+k(1-ebar).                                (N11)

The retained sharp pole bound says ebar>=a+3. Polynomial D
therefore forces k1,ebar=a+3. The retained equality theorem
makes that outer phase the sharp Phi1 family. Consequently a
nonconstant-Phi POLYNOMIAL outer phase admits no finite pole of w;
its admissible rational source maps must actually be polynomial.
The outer rational phases of(N3) are not subject to this premise.

These examples explain why the reflection normal form allows both
rational w and rational outer coefficients.

## 7. Compatibility with the former rational-primitive stable field

**Lemma RationalPrimitiveStableFieldEquality.** Let R be the normalized
rational primitive part of omega=Tdu, including any polynomial part.
Then R belongs to E=K(u,T). Consequently

 K(u,T)=K(u,R,T).

This holds at every height once u,T are rational. It does not assert
that the new anti-invariant generator is the former polynomial
generator used in cubic-height descent.

To prove it, write E=K(w). The differential descends to
eta=T(w)U'(w)dw. Decompose eta into a rational exact part dRbar
and logarithmic parts. Pullback gives dRbar(w) and logarithmic
differentials of rational functions of x. Such logarithmic
differentials have only simple finite poles and O(1/x) coefficients
at infinity. A derivative of a nonconstant rational function cannot
have that form: a finite primitive pole gives a differential pole
of order at least2, and a polynomial primitive part gives a
polynomial derivative. Hence uniqueness of rational exact/log
separation gives R-Rbar(w) constant.

The descended differential is regular at the marked point because
w is unramified there and omega is regular. Thus Rbar is regular
there and may be normalized to value0, as R is. The constant is0.
The rational primitive parts are defined over K by their normalized
uniqueness, so the argument does not require coefficient extension.
This proves the field equality. Its use of exact/log separation is
confined to this compatibility lemma, not the core normal form.

## 8. Sources and verification scope

J.S. Milne, *Fields and Galois Theory*, v5.10, Theorem9.19
(printedp117) is the actual arbitrary-ground-field Lüroth statement.
Lemma9.24 (printedp118) gives the rational-map field degree formula,
and Example9.25 (printedp119) identifies degree-one generators as
Mobius functions. These actual statements and their hypotheses were
read in the retained primary PDF.

The canonical odd-start rationality theorem and all-height endpoint
criterion are imported as stated. The marked-place normalization,
split involution, generator uniqueness, composition converse and
finite divisibilities are supplied explicitly above. The optional
polynomial-outer pole restriction uses the retained sharp degree
and equality theorem. No historical novelty of Lüroth,
fractional-linear involutions or rational pullbacks is asserted.

The standard-library verifier checks the involution straightening,
full rational reconstruction and divisibilities on exact fixtures,
including nonpolynomial w, rational outer coefficients, hyperbolic
data and the degree59 pullback. It compares raw angular moments
and rejects malformed parity, base and polynomiality data. These
finite controls supplement the written all-field, all-order proof;
they do not provide a finite purity prefix or external verification.


The [checker](verify_odd_start_rational_composition_normal_form.py),
[inputs](odd_start_rational_composition_inputs.json),
[source ledger](odd_start_rational_composition_source_ledger.json) and
[receipt](odd_start_rational_composition_verification.json) retain the exact
reconstruction and polynomiality controls.
