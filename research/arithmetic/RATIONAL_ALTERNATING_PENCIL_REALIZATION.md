# Rational degree-eight alternating pencils

The rational-coefficient geometric realization question has a positive answer.
There are infinitely many pairwise nonisomorphic covers over Q with four finite
double-transposition branch values, pole profile (7,1), and geometric group A8.
Their discriminants have constant square class -7. This is compatible with the
[earlier universal obstruction](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md):
a literal square in Q[s] is impossible.

The earlier [quadratic-field example](DEGREE_EIGHT_ALTERNATING_PENCIL.md) still
has no rational form. The examples here are different covers. The coefficient
family below is not asserted to classify every cover with this passport.

## Theorem 1. An exact rational coefficient family

For c,t in Q put

\[
 d=\frac{47c^2+7t^2}{168},\qquad
 e=\frac{379c^3+357ct^2+196t^3}{21168},\qquad
 a=-\frac{(7t^2-3c^2)^2}{87808},
\]
\[
 H_{c,t}(z)=z^8+cz^6+dz^4+ez^2,\qquad
 f_{c,t}(z)=\frac{H_{c,t}(z)+a}{z}.
\]

Define Q_(c,t)(s)=s^4+B s^2+C, where

\[
 B=\frac{6912c^7-14112c^5t^2+5488c^4t^3+8575c^3t^4
             -7203c^2t^5+31213ct^6+16807t^7}{533655864},
\]
\[
 C=\frac{(7t^2-3c^2)^2(4c^2-7ct+7t^2)^2
                  (48c^2+56ct+35t^2)^3}{379718108242115328}.
\]

Then the following is an identity over Q[c,t,s]:

\[
 \operatorname{Disc}_z(H_{c,t}(z)-sz+a)
       =-7^7 Q_{c,t}(s)^2.                                      \tag{1}
\]

If

\[
                    t\ne0,\qquad t\ne2c,\qquad t\ne-c,          \tag{2}
\]

then a is nonzero, f has eight distinct simple finite critical points,
Q_(c,t) has four distinct roots, and precisely two critical points lie over
each of its roots. Thus the finite branch permutations are four double
transpositions, and infinity has profile (7,1).
Condition (2) is the exact nondegenerate domain for this rational family.

Proof. Write

\[
 T(u)=7u^4+5cu^3+3du^2+eu-a.
\]

Then z^2 f'(z)=T(z^2). Exact calculation gives

\[
 \operatorname{Disc}_u T
 =\frac{t^2(t-2c)^6(c+t)^2(48c^2+56ct+35t^2)}{1259712},          \tag{3}
\]
\[
 B^2-4C=-\frac{t^6(t-2c)^5(c+t)^3}{10806531246}.                \tag{4}
\]

The quadratic forms 4c^2-7ct+7t^2 and 48c^2+56ct+35t^2 are positive
definite over R. The equation 7t^2=3c^2 has no nonzero rational solution,
by the odd valuation of 3 or 7. Under (2), the pair (c,t) is nonzero, so
a and C are nonzero. Equations (3) and (4) give both required
squarefreeness statements. Conversely, t=0,2c,-c makes (3) and (4)
vanish, so it does not give this simple-branch passport.

Because T(0)=-a is nonzero, the four distinct u-roots lift to eight distinct
nonzero z-roots. Each is a simple zero of f'. The
[exact resultant identity](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md),
Lemma 3, counts discriminant order one for each such critical point.
Equation (1), with Q squarefree, therefore pairs the eight critical points
two at each of four values. Finally f has poles of orders seven and one.

The retained [checker](verify_rational_alternating_pencil.py) recomputes
(1), (3), and (4) as exact symbolic identities in all the displayed
indeterminates. These identities do not rest on a finite parameter sample.

One way to derive the parametrization is to begin with
H=z^8+cz^6+dz^4+ez^2. The coefficient equations for its normalized
discriminant to be (s^4+B s^2+C)^2 have a component with

\[
 a=-\frac{(25c^2-84d)^2}{21952},
\]
\[
 4148928e^2+(791056c^3-3358656cd)e+
 145375c^6-1474760c^4d+4806704c^2d^2-4917248d^3=0.
\]

The discriminant of this quadratic in e is

\[
                  2^{10}7^5(168d-47c^2)^3.
\]

Taking 168d-47c^2=7t^2 makes it a square and gives the formulas above;
changing t to -t exchanges the quadratic roots. The identity check proves
the displayed family, without asserting that this component exhausts
the whole coefficient locus.

## Proposition 2. A rational example with full geometric A8

The member c=t=1 is

\[
 H(z)=z^8+z^6+\frac9{28}z^4+\frac{233}{5292}z^2,\qquad
 a=-\frac1{5488},
\]
\[
 Q(s)=s^4+\frac{5960}{66706983}s^2+
                    \frac{2685619}{1483273860320763}.
\]

It has the passport of Theorem 1, geometric monodromy A8, arithmetic
monodromy S8, and Galois-closure constant field Q(sqrt(-7)).
Its discriminant is a square over Q(sqrt(-7))[s].

Proof. The critical discriminant is 139/314928, and (4) is nonzero.
For the group argument, H(z)-sz+a is irreducible over Q(s):
as a polynomial linear in s its two coefficients -z and H(z)+a are
coprime, and Gauss's lemma applies. The arithmetic group G is transitive.
The pole of order seven supplies a 7-cycle. A block containing its fixed
sheet is either that sheet or all eight sheets, so G is primitive.

Modulo 11 at s=2 the polynomial is

\[
 z^8+z^6+7z^4+2z^2-2z+1
 =(z^3+z^2-3z-4)(z^5-z^4+5z^3-4z^2-3).
\]

The two factors are distinct and irreducible over F_11. All denominators
are units at 11, and the reduced polynomial is squarefree. Frobenius
therefore supplies a (3)(5) permutation in G; its cube is a 5-cycle.

The imported Jordan theorem, in the inspected
[Jones formulation](sources/jones-primitive-cycle-1209.5169v1.pdf),
Theorem 1.1, says that a primitive group containing a prime cycle with
at least three fixed points contains the alternating group. It gives
A8 contained in G. Equation (1) has nonsquare constant class -7 over Q(s),
so G contains an odd permutation and is S8.

The geometric group is a normal subgroup of S8, is contained in A8
because the discriminant becomes a square over algebraic constants, and
contains the nontrivial 7-cycle inertia. Simplicity of A8 forces it to be
A8. The quotient S8/A8 has order two and is detected by sqrt(Delta);
since sqrt(Delta)=343 sqrt(-7) Q(s), the algebraic constant field is
exactly Q(sqrt(-7)).

This is a genuine rational coefficient realization of geometric A8.
It does not give an arithmetic A8 cover over Q(s), which the universal
leading-coefficient obstruction excludes.

## Corollary 3. Infinitely many rational covers with the same passport

For every integer k>=0, set c=1 and t=1+11k in Theorem 1. Each member has
the same four-double-transposition passport, geometric A8, arithmetic S8,
and constant field Q(sqrt(-7)). These covers are pairwise nonisomorphic,
even allowing geometric changes of both source and target coordinates.

Proof. These t are positive and avoid 0,2,-1, so (2) holds. Reducing
coefficients modulo 11 gives exactly the member c=t=1. At s=2 the same
irreducible cubic and quintic supply the same 5-cycle certificate. The
7-cycle and Jordan argument above applies to each member separately.

For inequivalence, the unique (7,1) branch value fixes target infinity.
Its two poles, distinguished by their orders, fix source infinity and zero.
Every equivalence between normalized Laurent forms is therefore a source
scaling z -> lambda z and target affine transformation
s -> mu s+nu. For coefficients b_j of the Laurent expression, the
nonzero-denominator ratio

\[
                     I=\frac{b_1 b_7^2}{b_5^3}
\]

is invariant: numerator and denominator both acquire mu^3 lambda^15.
Here b_7=b_5=1 and I=e(t), with

\[
                     e(t)=\frac{379+357t^2+196t^3}{21168}.
\]

This is strictly increasing for positive t, since its derivative is
(714t+588t^2)/21168>0. The chosen integer sequence thus gives distinct
geometric cover classes. Source or target conic forms introduce no extra
equivalences: the distinguished branch point and the two distinguished
poles are rational points and split the genus-zero curves.

## Discriminant geometry

For each cover in Corollary 3, the horizontal line at its rational value
a passes through exactly four ordinary nodes of the degree-eight reduced
plane discriminant of H(z)-sz+v. Here v denotes the second plane coordinate
to distinguish it from the coefficient parameter t.

At each branch value, the two distinct critical points r,u are nonzero
and have H''(r)H''(u) nonzero. The discriminant normalization has derivatives
H''(r)(1,r) and H''(u)(1,u), with distinct tangent slopes. The two branches
are smooth and transverse. Each node consumes two intersections with
v=a, and the four nodes exhaust degree eight. This is the sharp
collinear-node configuration established in the
[all-degree bridge](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md).
The assertion concerns nodes on this line; it imposes no additional
generic-singularity claim away from it.

## Verification and limits

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_rational_alternating_pencil.py
```

The [receipt](rational_alternating_pencil_verification.json) records
Python 3.12.3 and SymPy 1.14.0 from the
[pinned research environment](../requirements.txt). It checks the
all-parameter discriminant identity, both nondegeneracy factorizations,
the seed model, the irreducible (3,5) factorization modulo 11, preservation
of that reduction for every t=1+11k, and invalid-input controls.

The algebraic identities are exact computational evidence. Monodromy and
the infinite inequivalence assertion additionally use the written proofs
and the specified imported Jordan theorem. This is internal verification;
independent human review, formal verification, and exact priority are
not claimed. No statement here gives rational realizations in every
even degree or a polynomial affine-space Keller map.
