# An explicit degree-eight alternating pencil

This is an exact algebraic example with four finite double-transposition branch
values and geometric monodromy A8. Its coefficients lie in a quadratic number
field. The displayed cover has no form over Q; existence of a different rational
coefficient example up to a constant discriminant factor remains a separate
question. No identification with the braid orbit of the particular knot tuple
is claimed.

## The polynomial and its discriminant

Put

```text
w^2 = -329,
e = (-1009 + 94 w)/10584,
a = -625/21952,
H(z) = z^8 + z^6 + e z^2.
```

Either choice of w gives a conjugate example. The coefficient e has irreducible
polynomial

```text
p(e) = 4148928 e^2 + 791056 e + 145375.
```

Let

```text
B = (3246152 e - 22375)/19765032,
C = -625 (6143242504 e + 991745875)/3515908409649216,
Q(s) = s^4 + B s^2 + C.
```

Then the exact identity is

```text
Disc_z(H(z) - s z + a) = -7^7 Q(s)^2.
```

In particular it is a literal square over
L = Q(sqrt(-329),sqrt(-7)), with square root
343 sqrt(-7) Q(s). Over K = Q(sqrt(-329)), the constant square class is -7.

The [exact checker](verify_degree_eight_alternating_pencil.py) computes the
discriminant over Q[e] and reduces the identity modulo p. All arithmetic is exact.
The [internal replay receipt](degree_eight_alternating_pencil_verification.json)
records Python 3.12.3 and SymPy 1.14.0, pinned by
[the research requirements](../requirements.txt).

## The branch cycles are double transpositions

For f(z) = H(z)/z + a/z, differentiation gives

```text
z^2 f'(z) = 7 z^8 + 5 z^6 + e z^2 - a = T(z^2),
T(u) = 7 u^4 + 5 u^3 + e u - a.
```

Its critical polynomial is squarefree. Indeed,

```text
Disc_u T = -(32616905272 e + 10887485125)/784147392,
Norm_K/Q(Disc_u T) = 1119061279296875/7809875684352.
```

Since a is nonzero, T has no zero root. Thus T(z^2) has eight distinct
nonzero roots, and all finite critical points of f are simple.

The quartic Q is also squarefree. For Q = s^4+B s^2+C, it suffices that
C and B^2-4C are nonzero. Reduction modulo p gives

```text
B^2 - 4C = -2209 (171096632 e + 35613875)/329616413404614.
```

Both this expression and C are nonzero, since e has degree two over Q and each
displayed numerator has degree one. The checker additionally computes the
nonzero norm of Disc_s Q.

For completeness, the discriminant multiplicity counts the critical points here.
Writing N(z) = H(z)+a and R(z) = zN'(z)-N(z), one has R = zF_s'-F_s for
F_s(z)=N(z)-sz. For fixed nonzero a, the resultant identity gives

```text
Disc_z F_s = constant * product_{R(r)=0} (f(r)-s).
```

The constant is nonzero, and all roots r of R are distinct. This can also be seen
locally from a simple ramification point: two roots coalesce with square-root
separation, so the discriminant has a simple zero per critical point. Since Q
has four distinct roots, and the discriminant equals a nonzero constant times
Q squared, each finite branch value has exactly two distinct simple critical
points. Its local permutation is therefore a double transposition.

The remaining branch point is infinity. The two poles have orders seven at
infinity and one at zero, giving inertia (7)(1). This is the required passport:
four double transpositions and one (7)(1).

## Full geometric monodromy

The square discriminant alone would allow proper subgroups such as AGL3(2).
An additional exact certificate excludes that possibility.

Work over L = K(sqrt(-7)); the discriminant is then a square, so the arithmetic
Galois group G of H(z)-sz+a over L(s) is contained in A8. The polynomial is
irreducible: it is linear in s in L[z,s], and its s coefficient -z is relatively
prime to N(z), since N(0)=a is nonzero. Hence G is transitive.

The inertia at infinity supplies a 7-cycle fixing one point. In a transitive
degree-eight action this forces primitivity. To see this, the block containing
the fixed point is invariant under that 7-cycle, so it contains either only that
point or all eight points. Any block system is therefore trivial.

There is a prime of L above 11 with

```text
sqrt(-329) -> 1, sqrt(-7) -> 2, e -> 10, a -> 5.
```

At s=2 the reduction factors as

```text
z^8 + z^6 + 10 z^2 - 2 z + 5
 = (z^3 + 4z^2 - 4z - 3)
   (z^5 - 4z^4 - z^3 + 2z^2 - 2z + 2)    over F_11.
```

The two factors are distinct and irreducible over F_11; the checker verifies
these facts exactly. All coefficient denominators are units at 11, and the
reduced polynomial has nonzero discriminant. Its splitting algebra is finite
etale near this specialization. A decomposition group at the residue point
therefore lifts the residue Frobenius to a permutation of the eight generic
roots. Frobenius has cycle type (3)(5) in a
specialization subgroup of G. Its cube is a 5-cycle, fixing three sheets.

Jordan's theorem now gives A8 contained in G, so G=A8. The precise imported
statement is Theorem 1.1 of
[Gareth A. Jones, *Primitive permutation groups containing a cycle*](https://arxiv.org/pdf/1209.5169v1):
a primitive finite permutation group containing a prime-length cycle with at
least three fixed points contains the alternating group. The
[retained source](sources/jones-primitive-cycle-1209.5169v1.pdf) and its exact
source record are in the
[bridge source ledger](alternating_hurwitz_source_ledger.json).

Finally, the geometric Galois group is normal in G=A8. It is nontrivial because
it contains the inertia at infinity. Simplicity of A8 implies that it is A8.
This proof does not infer full monodromy from the passport alone.

Over K(s), the arithmetic group is S8: it contains the A8 found after extension
to L, and the discriminant has nonsquare constant class -7 in K. The field
of constants in its splitting field is exactly L. These arithmetic statements
are supplementary; the geometric group is A8 independently of the chosen
coefficient field containing K.

## Rational descent of this example

The displayed cover has no rational form, even allowing geometric changes of
source and target coordinates.

The branch value with a ramification point of index seven is unique; it must
therefore be rational in any Q-defined form. Its two preimages, of indices seven
and one, are individually rational as well. Over Q one may consequently put that
branch value at target infinity and these two source points at infinity and zero.
A rational form would then be a Laurent polynomial with rational coefficients,
with powers from -1 through 7.

An isomorphism between any two such normalized covers acts by z -> lambda z on
the source and s -> mu s+nu on the target. If b_j denotes the coefficient of z^j,
the ratio

```text
I = b_1 b_7^2 / b_5^3
```

is defined because b_5=1 for our model, and nonvanishing is preserved by these
coordinate changes. It is invariant: both numerator and denominator acquire the factor
mu^3 lambda^15. Target translation changes only b_0. For our cover I=e, which
is not rational. A rational form would give a rational value of I, a
contradiction. In fact the field of moduli is exactly K, since e generates K
and the cover is already defined over K.

This does not rule out other degree-eight Laurent polynomials over Q with
discriminant a constant times a square. A literal square in Q[s] is impossible
for every rational degree-eight pencil H(z)-sz+a: the leading coefficient in s
of its discriminant is -7^7 h_8^6, a negative nonsquare. The distinction between
a literal square and a constant times a square is essential.

## Consequence for the plane discriminant curve

For this H the horizontal line t=a passes through exactly four ordinary nodes
of the curve Disc_z(H(z)-sz+t)=0. Each of the four branch values above has two
distinct simple critical points. The local discriminant parametrization is

```text
s=H'(z),    t=zH'(z)-H(z).
```

At a simple critical point of the restricted pencil, its derivative is
H''(z)(1,z) with H''(z) nonzero. Distinct critical points therefore give distinct
tangent slopes. There are exactly two smooth transverse branches at each of
the four points. No additional intersection with t=a is possible, since the
restricted discriminant is the degree-eight polynomial -7^7 Q(s)^2.

## Verification and scope

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_degree_eight_alternating_pencil.py
```

The retained receipt was produced by that checker with the additional
`--write-receipt research/arithmetic/degree_eight_alternating_pencil_verification.json`
argument. It checks the discriminant identity, the critical and branch
squarefreeness, the prime-11 factorization and factor irreducibility, and three
negative controls. The monodromy and descent conclusions use the arguments
above and their stated imported theorems. This is internal exact replay and
internal mathematical checking, without a claim of formal verification,
external review, novelty, or a certified braid-orbit match to the knot tuple.
