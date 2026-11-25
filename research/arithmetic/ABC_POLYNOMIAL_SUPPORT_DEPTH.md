# Polynomial fake support occurrence cannot have deep simultaneous precision

Let \(A,B,C\) be pairwise coprime nonconstant polynomials over
a characteristic-zero field and put \(D=A+B-C\ne0\).
If \(A,B,C\) are linearly independent over the constants, then
\[
              \sum_{P\mid ABC} v_P(D)\deg P
                   \le 3\deg\operatorname{rad}(ABC)-3.          \tag{1}
\]
Here the sum is over monic irreducible polynomial divisors.
In particular a fake polynomial family with a small *formal*
radical cannot meet high-depth additive congruences on its
own colored supports.

Indeed \(W=W(A,B,C)\) is the nonzero order-zero-through-two
Wronskian. At a divisor \(P\) of \(A\), replace \(C\) by
\(D\) in its third column. This changes just the sign.
Put \(e=v_P(A)\), \(r=v_P(D)\). In every determinant term
the derivative orders landing in those two columns total
at most three. In characteristic zero each derivative can
lower a prime-order valuation by at most one. The remaining
polynomial column has nonnegative valuations. Thus
\(v_P(W)\ge e+r-3\). Apply the same argument on the other
two colored supports, and use disjointness to sum. Because
\(\deg W\le\deg A+\deg B+\deg C-3\), the weighted sum gives
(1). Factoring over an algebraic closure gives the equivalent
statement at individual roots.

If the constants rank is two, a nonzero \(D\) either
is proportional to one letter or has no common root with
any letter. In the latter case \(D,F\) are a basis for each
letter \(F\), and a common root would force every letter
to vanish there. In the former case \(D\) cannot meet
another colored support. Thus the rank-two exception cannot
give a fake passing positive precision on all three supports.

The same reasoning works after generic affine-line specialization
of a multivariate polynomial model, retaining the top degrees
and all support separations. Independence of the restrictions
need not be assumed: the low-rank case above handles divisibility
on all colors. The field of constants can include the generic
line parameters. This does not infer a bound for arbitrary
integer evaluations, since prime powers in evaluated values
need not come from formal polynomial factors.

More generally, suppose \(F_1,\ldots,F_m\) are pairwise coprime,
of constants rank \(r\ge2\), and \(D\ne0\) is any constant
linear combination independent of each individual \(F_i\).
All Wronskians of bases of their span agree up to constants.
At a divisor of \(F_i\), choose a basis containing \(D,F_i\).
The largest two derivative orders now sum to \(2r-3\).
Using one basis for the degree upper bound and all colored
supports for the degree lower bound gives
\[
 \sum_{P\mid \prod_i F_i}v_P(D)\deg P
 \le (2r-3)\deg\operatorname{rad}\!\left(\prod_i F_i\right)
                         -\binom r2 .
\]
The degree upper bound here may use the sum of all input
degrees, or a smaller sum for a particular basis.

## Actual depth-three passing pairs on degree-at-most-two factor families

There is an integer-specialization consequence that does not
assume evaluated radicals have their formal sizes. Let
\(A,B,C\) be homogeneous binary forms of the same fixed
positive degree, pairwise coprime over \(\mathbf Q\), with uniformly
subpower coefficient heights, and let \(D=A+B-C\ne0\).
Assume all distinct irreducible factors of \(ABC\) over
\(\mathbf Q\) have degree at most two, and write \(\ell\)
for their number. For primitive integer parameter pairs
\((s,t)\) of dyadic height \(H\), require
\(A(s,t)B(s,t)C(s,t)D(s,t)\ne0\) and the actual condition
\[
             \operatorname{rad}(A(s,t)B(s,t)C(s,t))^3
                  \ \mid\ \gcd(D(s,t),(ABC)(s,t)).              \tag{2}
\]
In particular (2) follows when the original values all
have prime exponents at least three and all the third-layer
tests pass. Then the number of such primitive pairs is
\[
                          H^{\,2-1/(2\ell+1)+o(1)}              \tag{3}
\]
as an upper bound. The subpower constants are uniform
under the coefficient assumptions at fixed degrees.

Suppose first the constants rank of \(A,B,C\) is three.
The homogeneous version of (1) follows by choosing a
projective chart whose point at infinity avoids all roots
of the forms, so every degree is preserved. Put
\(G=\gcd(D,ABC)\) as a primitive binary polynomial, and
\(S=\deg\operatorname{rad}(ABC)\). Then
\(\deg G\le3S-3\). The coprime forms left after removing
\(G\) give a fixed nonzero homogeneous resultant. Its
Bezout identities with powers of the two homogeneous
coordinates show that at a primitive pair, the gcd
of their evaluated values divides that resultant, with
scalar contents included. These contents and the
resultant have subpower size, by fixed-degree coefficient
and factor height bounds. Consequently
\[
 \gcd(D(s,t),(ABC)(s,t))\le H^{3S-3+o(1)},\qquad
 \operatorname{rad}(A(s,t)B(s,t)C(s,t))
                     \le H^{S-1+o(1)}                         \tag{4}
\]
by (2).

Write \(F_1,\ldots,F_\ell\) for the distinct primitive
formal factors. The product of their evaluated radicals
is at most the radical in (4) times a subpower factor:
overlapping primes are on the support of fixed pairwise
homogeneous resultants. Pick \(0<\theta<1\).
If \(\prod_i |F_i(s,t)|\ge H^{S-\theta}\), then
\[
             \prod_i\frac{|F_i(s,t)|}
                              {\operatorname{rad}(F_i(s,t))}
                         \ge H^{1-\theta-o(1)} .
\]
At least one factor has deficit \(H^{(1-\theta)/\ell-o(1)}\).
For any nonzero integer \(n\), its largest square divisor
is at least \(|n|/\operatorname{rad}(n)\), prime by prime.
Thus the selected factor value is divisible by \(r^2\)
for some \(r\ge H^{(1-\theta)/(2\ell)-o(1)}\).
Counting multiples of such squares up to its cap
\(H^{\deg(F_i)+o(1)}\) gives at most
\(H^{\deg(F_i)-(1-\theta)/(2\ell)+o(1)}\) values.
A linear value has \(O(H)\) input pairs. An irreducible
quadratic value has subpower many pairs by the
coefficient-uniform quadratic norm bound
[for fixed-factor recovery](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md). In either case
this branch has count
\(O(H^{2-(1-\theta)/(2\ell)+o(1)})\).

Otherwise \(\prod_i |F_i(s,t)|<H^{S-\theta}\).
Choose a real projective chart with denominator comparable
to \(H\). Its parameter approaches a real root of
the product within \(H^{-\theta+o(1)}\).
Indeed the root factorization first gives a root within
\(H^{-\theta/S+o(1)}\). All distinct roots are separated
by \(H^{-o(1)}\) using the fixed-degree squarefree
discriminant and coefficient bounds, so exactly one is
this close. A nonreal root could not be this close to
the real axis since its conjugate is separated. All
other factors are bounded below by subpower losses,
giving the stronger stated proximity. Counting integer
pairs in the finitely many resulting parameter intervals
costs \(O(H^{2-\theta+o(1)}+H)\). Taking
\(\theta=1/(2\ell+1)\) balances the two branches.

If the constants rank is two, either \(D\) is proportional
to one letter or it is coprime to all of them. In both
cases choose a letter coprime to \(D\). Its evaluated
radical under (2) lies on the fixed subpower resultant
support. Any one of its irreducible factors then has only
subpower choices for its value. Fixing a linear value
has \(O(H)\) input pairs and fixing a nonlinear quadratic
value has subpower many, as above. This gives (3) in
the lower-rank case also.

For example this proves a power loss in the number of
actual parameter pairs in any fixed-format linear and
quadratic polynomial countermodel attempting to pass
third layers. If its target radical scale is \(H^2\),
as in the usual primitive Pythagorean parametrization,
the passing target count is also a power below \(H^2\). The scope is fixed-format polynomial
parametrizations with subpower coefficient heights.

## Extension through degree four

The same exponent (3) holds when every irreducible factor
has degree at most four. Only the square-divisor step for
a chosen factor \(F\) changes. Set
\(\eta=(1-\theta)/(2\ell)\), so a square divisor \(r^2\)
has \(r\ge H^{\eta-o(1)}\).

For \(r\le H^{7/4}\), use projective congruence packing rather
than bounding the integer values of \(F\).
There are \(H^{o(1)}\) projective residue classes at modulus
\(r^2\) in which \(F(s,t)\equiv0\) for primitive pairs.
To verify this, at each prime choose a unit-coordinate
chart. The univariate defining polynomial in that chart
is separable over characteristic zero with subpower
derivative resultant \(J\). At a modulus prime of exponent
\(e>2v_p(J)\) every approximate root has derivative valuation
at most \(v_p(J)\) and lifts to one of the finitely many
exact roots at precision \(p^{e-v_p(J)}\). Below that
threshold the residue-class count is trivially at most
\(p^{2v_p(J)}\). Combining by CRT, including the two
chart choices per prime, gives the asserted subpower bound.

Distinct primitive parameters in one class have cross
determinants divisible by \(r^2\). Angular packing in the
height-\(H\) box gives \(O(1+H^2/r^2)\) such pairs up to
the bounded sign multiplicity. The sum over the indicated
\(r\) is
\[
               O(H^{2-\eta+o(1)}+H^{7/4+o(1)}).
\]
The latter term is smaller at the balancing choice
\(\eta=\theta=1/(2\ell+1)\), since \(\ell\ge3\).

A larger \(r\) is impossible for factors of degree at
most three. In degree four and \(r>H^{7/4}\), put
\(n=F(s,t)/r^2\), a nonzero integer with
\(|n|\le H^{1/2+o(1)}\). For each \(n\), the primitive pairs
embed, up to sign, as rational points
\([s^2:st:t^2:r]\) on the curve in projective three-space
defined by
\[
                UW=V^2,\qquad n Z^2=Q_F(U,V,W),
\]
where \(Q_F(s^2,st,t^2)=F(s,t)\). This curve is absolutely
irreducible of degree four. On \(W=1\) the equation
\(n Z^2=F(V,1)\) is irreducible since \(F(V,1)\) has
distinct roots. At \(W=0\) there are only finitely many
points, and the two quadric equations have no common
surface component. Heights in these coordinates are
\(H^{2+o(1)}\). The coefficient-uniform projective curve
bound \(O_{d,\varepsilon}(B^{2/d+\varepsilon})\) gives
\(H^{1+o(1)}\) points per \(n\), and hence
\(H^{3/2+o(1)}\) pairs in this large-divisor range.
[Heath-Brown, The density of rational points on curves and surfaces](https://arxiv.org/abs/math/0405392)

The other angular and rank branches are unaffected.


## Depth three via a formally deficient factor

There is a stronger estimate covering higher-degree forms.
Under the assumptions of (2), let \(F_i\) be any distinct
formal irreducible factor of \(ABC\) with \(r_i=v_{F_i}(D)<3\).
In constants rank three, (1) ensures that there is such a
factor. In rank two one can choose a factor coprime to
\(D\) by the rank-two argument above.

Outside a set \(\mathcal S\) whose prime product is
\(H^{o(1)}\), any prime of \(F_i(s,t)\) satisfies
\[
                   v_p(D(s,t))=r_i v_p(F_i(s,t)).
\]
Take \(\mathcal S\) to contain the primes of the
homogeneous resultant of \(F_i\) and \(D/F_i^{r_i}\)
and the relevant scalar contents. The nonzero resultant
has subpower size and bounds common prime divisors at a
primitive pair. By (2) the evaluated factor is supported
on \(\mathcal S\) if \(r_i=0\); it is cube-full away from
\(\mathcal S\) if \(r_i=1\), and square-full away from
\(\mathcal S\) if \(r_i=2\).

Consequently a deficient factor of degree \(d\ge5\)
gives the bound \(H^{d/5+o(1)}\) for \(r_i=2\), and
\(H^{d/6+o(1)}\) for \(r_i=1\). The argument uses the
prime-power kernel choices on the actual values of a
single factor and the uniform primitive fixed-power count.
A rank-two factor or any factor with \(r_i=0\) gives
at most \(H^{1+o(1)}\) pairs for all degrees, by the
fixed-support and fixed-value counts.
[Primitive powerful values of binary forms](ABC_BINARY_POWER_TWISTS.md#square-full-and-cube-full-values)

For completeness, a deficient factor of degree at most
three allows at most \(H^{3/2+o(1)}\) pairs. Its values
are square-full away from \(\mathcal S\). The integers
supported on \(\mathcal S\) in the polynomial height
range have subpower count, and the square-full count at
most \(U\) is \(O(U^{1/2})\) by the representation
\(u^2v^3\) with \(v\) squarefree. Fixing an irreducible
nonlinear value costs only subpower many pairs, while a
linear value costs \(O(H)\) and has cap \(O(H^{1+o(1)})\).

For a quartic deficient factor, \(|F_i(s,t)|<H^{15/4}\)
costs \(O(H^{7/4+o(1)})\) by the real root localization
used above. Otherwise the largest square divisor \(r^2\)
of its absolute value satisfies
\(r^2\ge |F_i(s,t)|^{2/3}H^{-o(1)}\).
Indeed an exponent \(e\ge2\) contributes at least
\(2e/3\) to the maximal square divisor, and the product
of the primes having exponent one divides
\(\prod_{p\in\mathcal S}p\). Thus \(r\ge H^{5/4-o(1)}\).
For \(r\le H^{7/4}\) the projective congruence packing
above sums to \(O(H^{7/4+o(1)})\). For larger \(r\),
the quartic curve argument above gives \(O(H^{3/2+o(1)})\).

In particular, if all formal factors have degree at
most nine, (2) has at most \(H^{9/5+o(1)}\) primitive
parameter pairs. The same conclusion holds whenever
there is a deficient factor of degree at most nine,
irrespective of all other factor degrees. A deficient
factor of degree ten or eleven with \(r_i=1\) also
gives a power saving, using the cube-full bound.
