# A rational degree-ten alternating pencil with square discriminant

The degree-ten rational target is realized by an explicit integer-coefficient
polynomial. Its five finite branch values each have two distinct simple
critical points. Both arithmetic and geometric monodromy are A10.

## Theorem 1. The rational degree-ten example

Put

\[
 H(z)=9z^{10}-864z^8+28161z^6-350150z^4+1127500z^2,
 \qquad a=250000,
\]
\[
 W(s)=59049s^4-126938174857200000s^2+
                         68130645548641000000000000000.
\]

Then

\[
 \boxed{\operatorname{Disc}_z(H(z)-sz+a)
                    =\bigl(2187s\,W(s)\bigr)^2.}                 \tag{1}
\]

The square root has degree five and five distinct real roots. The rational
function f(z)=(H(z)+a)/z has exactly ten distinct simple finite critical
points, all real, paired two at each of those five branch values. Its two
poles have orders nine and one. Its geometric monodromy over Qbar(s) and
arithmetic monodromy over Q(s) are both A10; the algebraic constant field
of its splitting field is Q.

The [exact checker](verify_rational_degree_ten_pencil.py) recomputes (1) in
Z[s]. The [receipt](rational_degree_ten_pencil_verification.json) records the
input, arithmetic, software, exact identities, finite-field group certificate
and rejection controls. The monodromy proof below uses more than the square
discriminant.

## Distinct paired critical points

The numerator factors as

\[
 H(z)+a=(z^2-25)^2(9z^6-414z^4+1836z^2+400).
\]

Let

\[
 T(u)=81u^4-4023u^3+40230u^2-44700u+10000.
\]

Differentiation gives

\[
 z^2f'(z)=zH'(z)-H(z)-a=(z^2-25)T(z^2).                         \tag{2}
\]

The following exact values certify nondegeneracy:

\[
 \operatorname{Disc}_u T=1801825579550588085000000,\qquad
 T(25)=-7182500,\qquad T(0)=10000.
\]

Thus (u-25)T(u) has five distinct nonzero roots. Each lifts to two distinct
simple roots in z of (2), so the ten finite critical points are all simple.
A second direct check gives

\[
 \operatorname{Disc}_u((u-25)T(u))
                   =92953129806939475496581031250000000000\ne0.
\]

Write W(s)=A s^4+B s^2+C. Here A,C are positive, B is negative, and

\[
 B^2-4AC=21114280070272604391840000000000>0.
\]

The two roots of A y^2+B y+C are distinct and positive: their product is
C/A>0 and their sum is -B/A>0. Consequently W has four distinct nonzero
real roots, and sW(s) has five distinct real roots.

For completeness, the ten critical points are real as well. The signs of
T at u=0,1/2,1,25,50 are respectively +,-,+,-,+. The intermediate value
theorem supplies one root in each of the four disjoint positive intervals.
Since T has degree four, those are all its roots. Together with u=25 this
proves the claim about real critical points.

The [even-degree resultant identity](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md),
Lemma 3, expresses the discriminant, up to a nonzero constant, as the product
over the critical points r of s-f(r). Every simple critical point contributes
one to the discriminant order at its value. Equation (1), with its squarefree
degree-five square root, therefore places exactly two distinct critical
points over each of five values. The branch permutation there is a double
transposition. This rules out index-three ramification as the explanation
for the square discriminant.

Since H has degree ten and a is nonzero, f has poles of order nine at
infinity and one at zero. Its infinity inertia has type (9)(1).

## Full arithmetic and geometric monodromy

The polynomial F_s(z)=H(z)-sz+a is irreducible over Q(s), and also over
Qbar(s). Indeed it is linear in s in Qbar[z,s], with coprime coefficients
-z and H(z)+a because a is nonzero. Gauss's lemma gives irreducibility.
Thus both monodromy actions are transitive.

Let G be the arithmetic group over Q(s). The nine-cycle at infinity forces
primitivity: a block containing its fixed sheet is invariant under that
cycle, hence is either the singleton or the whole ten-point set. All blocks
in a transitive block system have equal size, so the block system is trivial.

An independent finite-field certificate gives the other required cycle.
Modulo 13, at s=12, the leading coefficient 9 is invertible, and the monic
normalization of F_12 factors as

\[
 9^{-1}F_{12}(z)=
 (z^3+2z^2+3z+1)
 (z^7-2z^6-4z^5-3z^3-3z^2+4z+4)
                  \quad\text{in }\mathbb F_{13}[z].              \tag{3}
\]

The factors are distinct and irreducible. The checker verifies their
irreducibility with Rabin's modular-power and gcd criterion, without calling
a group-classification package. The product is squarefree, so this is good
reduction for the specialization. Its Frobenius, lifted through the
unramified splitting algebra, supplies a permutation of cycle type (3)(7)
in G. Its cube is a seven-cycle fixing three points.

The precise imported Jordan theorem is Theorem 1.1 of Gareth A. Jones,
*Primitive permutation groups containing a cycle*, in the
[retained primary source](sources/jones-primitive-cycle-1209.5169v1.pdf).
A primitive group containing a prime cycle with at least three fixed points
contains the alternating group. Applying it to the seven-cycle gives
A10 contained in G. Equation (1) is a literal square over Q(s), so G is
contained in A10. Hence G=A10.

The geometric group is a normal subgroup of G=A10 and contains the
nontrivial nine-cycle inertia at infinity. Simplicity of A10 implies that
it is A10 as well. The arithmetic/geometric quotient is trivial, so the
splitting field has no algebraic constant extension of Q.

One could alternatively use the primitive double-transposition theorem in
Jones's retained Beauville paper, Proposition 6.1(3), since the degree is
at least nine. Certificate (3) instead connects this example directly to
the same prime-cycle Jordan theorem used for the degree-eight examples.

## Corollary 2. Five real collinear ordinary nodes

For the degree-ten reduced plane discriminant curve

\[
 \mathcal D_H=\{(s,t):\operatorname{Disc}_z(H(z)-sz+t)=0\},
\]

the rational horizontal line t=250000 contains exactly five ordinary nodes.
All five nodes are real and both branches at each node have real distinct
tangent slopes. This attains the maximum number of ordinary nodes on a line
for a degree-ten integral plane curve.

Indeed the [normalization and local-node argument](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md),
Theorem 4, applies here. The normalization is
r -> (H'(r),rH'(r)-H(r)); at a simple critical point of f its tangent is
H''(r)(1,r), with H''(r) nonzero. Two distinct real critical points over a
value therefore give two smooth transverse real branches. They are nonzero,
so the horizontal line is transverse to both. Each node consumes exactly two
of the ten intersections. No additional intersection remains by Bezout.
The theorem concerns the chosen line, without asserting generic singularities
elsewhere on the discriminant curve.

## A monic rational normalization

The same cover has the rational presentation

\[
 \begin{split}
 H_0(z)&=z^{10}-\frac{96}{25}z^8+\frac{3129}{625}z^6
                    -\frac{14006}{5625}z^4+\frac{1804}{5625}z^2,\\
 a_0&=\frac{16}{5625}.
 \end{split}
\]

If f_0=(H_0+a_0)/z, the integer model satisfies
f(z)=9\cdot5^9 f_0(z/5). These are rational source and target changes.
In the coordinate u=z^2 of the monic model,

\[
 H_0(z)+a_0=(u-1)^2
       \left(u^3-\frac{46}{25}u^2+\frac{204}{625}u+\frac{16}{5625}\right).
\]

This factorization separates the zero branch value from the four others.
The displayed model was obtained from the quartic critical-pair equations;
the retained certificate verifies the resulting coefficients directly.
No classification of rational degree-ten models is asserted.

## Verification and limits

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_rational_degree_ten_pencil.py
```

The retained receipt was produced with the additional argument
`--write-receipt research/arithmetic/rational_degree_ten_pencil_verification.json`.
The pinned environment is Python 3.12.3 with SymPy 1.14.0 from
[research requirements](../requirements.txt).

The checker verifies the full discriminant identity, both squarefreeness
conditions, real-root sign conditions, rational normalization, the exact
modulo-13 factorization and Rabin irreducibility. Its negative controls include
a square-discriminant degree-ten example whose critical polynomial is a
square and whose inertia has index three, demonstrating why that mechanism
does not answer this problem.

The proof imports Jordan's theorem and uses the stated good-reduction and
specialization argument. This is internal exact verification and internal
mathematical checking. No formal verification, independent external replay,
external human review or priority for the coefficient formula is claimed.
The result realizes arithmetic A10 over Q(s), rather than producing a
polynomial affine-space Keller map with generic monodromy A10.
