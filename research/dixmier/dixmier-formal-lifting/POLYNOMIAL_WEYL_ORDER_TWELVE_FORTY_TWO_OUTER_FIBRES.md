# Small-root exclusion and the outer fibres of the order-\((12,42)\) row

The order-\(34\) resultant provides a third coprime leading
exponent. At a small root, this forces continuation of a
monomial face until both components are constant on the
corresponding position line. At position infinity it gives
a separate reduction to Laurent charts of momentum degree
two or one. Together with the finite high-root pole calculation,
these charts determine the geometric degree of the transferred
planar pair.

Use \(x\) for position and \(X\) for momentum, with
\(\{F,G\}=F_XG_x-F_xG_X\).
Let \(d=\deg h\), and let \(H\) denote the number of
geometric multiplicity-fifteen roots. We work over an
algebraically closed field of characteristic zero after
the finite transfer.

## The two Euler faces in the center

**Proposition PolynomialWeylFortyTwoTwoEulerTransfer.**
The [planar transfer](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#transfer-of-the-high-root-profile)
of a hypothetical proper polynomial Weyl pair in the reduced
\((12,42)\) row can preserve the following additional data.
At every multiplicity-three root \(b\), the planar initials at
weights \(w(X)=1,\ w(x-b)=-1\) are nonzero multiples of
\((x-b)^6 X^{12}\) and \((x-b)^{21}X^{42}\).
At the global weight \(w(x)=1,\ w(X)=-1\), the initials are
nonzero multiples of \(x^{2d}X^{12}\) and \(x^{7d}X^{42}\).
The other high-root faces and scalar successor caps in the
planar transfer can be kept simultaneously.

At a small root the
[Weyl lowest Euler polynomials](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_RESULTANT.md#a-third-indicial-shift-at-a-small-root)
have exponents \(-6,-21\) and degrees \(12,42\) in
\(E=(x-b)D_x\).
For the outer Euler weight first apply the
[infinity shear](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree).
Here \(d>6\). The termwise bound in the proof gives
\(\deg p_i-i\le2(d-6)\) and \(\deg q_i-i\le7(d-6)\)
for the coefficients of \(D_x^i\).
The leading differential coefficients give equality.
Thus the highest global Euler components have exponents
\(2(d-6),7(d-6)\) and degrees \(12,42\) in \(xD_x\).

The shear is \(D_x\mapsto D_x+f(x)\) with \(f\) polynomial.
It preserves the finite root data. In fact it is the
exponential of commutation with a polynomial in \(x\),
up to sign. At a translated finite position, subtract
its constant value from that polynomial before using
the commutator. The result has negative weight at the
small Euler face and at the first high face.
In the shifted ramified chart \(x-b=\eta+\omega/t\),
with weights \(w(t)=7,\ w(\eta)=-12\), it has weight
at most \(-7\). The contraction loss there is nine.
The exponential is finite on original Weyl polynomials
and each correction strictly lowers the indicated local
weight. Thus the local caps and leading symbols survive.

Impose these Euler supports in the finite system used
for the planar transfer. After a number-field solution
is chosen, enlarge the splitting field to include the roots
of the finitely many Euler coefficient polynomials.
For a component \(y^a f(E)\), \(E=yD_y\), with
\(a\ne0\) in \(\mathbb F_p\), its \(p\)-th power at
the Euler grade is exactly
\[
        y^{pa}\prod_{j=0}^{p-1}f(E+ja)
          =(\operatorname{lc} f)^p y^{pa}(E^p-E)^{\deg f}.
\]
Indeed the Euler grading is exact and the shifts of
each split root range over \(\mathbb F_p\).
Since \(E^p-E=y^pD_y^p\), these become the asserted
monomials in the center. Here \(y=x-b\) at a finite
position and \(y=x\) at the global weight.
The other weights use the positive-contraction calculation
of the planar transfer. The same finite system return,
using the
[center Poisson theorem of Lauritzen–Thomsen](https://arxiv.org/pdf/1508.02499v1),
gives the simultaneous characteristic-zero data.

## A common primitive face for three operators

**Lemma PlanarFortyTwoThreeWeightPropagation.**
Let \(P,Q,W\) be the transferred planar pair and its
[polynomial resultant](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-shared-planar-resultant).
Use either orientation \(w_s=(1,-s)\) on the two variables,
with \(s\ge1\). Suppose at the beginning of an interval
all three weights are positive, in ratio \(6:21:17\).
On the connected part continuing to the right while the
weight of \(P\) is positive, the three weights retain that
ratio, and their leading symbols are, up to scalars,
\[
                               R^6,\quad R^{21},\quad R^{17}.
\]
The same statement works after a Laurent shear of \(X\)
by a function of \(x\), in \(k[x,x^{-1},X]\).

Write \(p_s,q_s,w_s^*\) for the three weights.
The scalar bracket of \(P,Q\) forces their initial forms
to commute while both weights are positive, since
\(p_s+q_s+s-1>0\).
On each open interval without ties the leading monomials
then have proportional exponent vectors, so \(q_s/p_s\)
is constant there. By continuity the ratio \(21/6\)
persists. For \(W\), the leading bracket with \(P\)
vanishes whenever \(w_s^*+p_s+s-1>q_s\), by
\(\{W,P\}=2\kappa Q\).
The same ratio argument gives \(w_s^*/p_s=17/6\)
on the initial connected interval of that inequality.
The slack there is \(p_s/3+s-1>0\), so this interval
continues as long as \(p_s>0\).
Finally the weighted Euler identities convert the
commutation into a constant ratio of powers.
Factorization and \(\gcd(6,21,17)=1\) give the
single polynomial or Laurent polynomial \(R\).
The lattice exponents of the Laurent unit are also
divisible by the required powers since their gcd is one.

## Exclusion of the multiplicity-three roots

**Theorem PolynomialWeylFortyTwoNoSmallRoots.**
In the reduced order-\((12,42)\) row of a proper polynomial
Weyl pair, every geometric root of \(h\) has multiplicity
fifteen. In particular \(H\) is even and at least two,
and \(\deg h=15H\) is a positive multiple of thirty.

Suppose instead that a small root exists, and take the
two-Euler planar transfer, translating that root to zero.
At \(w(X)=1,\ w(x)=-1\), \(P,Q\) have the asserted
pure top monomials. The planar resultant \(W\) has
\(X\)-degree \(34\) with leading coefficient proportional
to \(h^{17/3}\). Its Euler weight is at least \(17\)
by the order \(17\) of that coefficient at the root.
If its Euler weight is at least \(17\), it exceeds the
threshold for a nonzero leading bracket with \(P\):
the sum of weights would exceed \(\operatorname{wt}Q=21\).
Every monomial commuting with \(x^6X^{12}\) is a
power of \(xX^2\). The \(X\)-degree bound on \(W\)
therefore gives Euler weight exactly \(17\) and a pure
top monomial. The three tops are powers \(6,21,17\)
of \(xX^2\).

Now increase the penalty on \(x\) using \(w(X)=1\),
\(w(x)=-s\), from \(s=1\) to \(2\).
The initial monomial of \(P\) gives positive weight
for \(s<2\), so the three-weight lemma applies.
At a first tie in \(1<s<2\), the common primitive face
would have the endpoint \(xX^2\).
Any other support point would have strictly smaller
\(x\)-degree because \(xX^2\) wins immediately before
the tie. Its \(x\)-degree therefore is zero. Its
\(X\)-degree would be \(2-s\), which is not an integer.
There is no such tie. At \(s=2\), continuity gives
nonpositive weights for both \(P,Q\), so both are
constant on \(x=0\). This contradicts the nonzero
Jacobian. The prior degree divisibility \(6\mid\deg h\)
now gives the assertion about \(H\).

## The outer Newton path

**Theorem PlanarFortyTwoOuterFibreProfile.**
Take the two-Euler planar transfer of the remaining
no-small-root profile, and put \(D=d/3=5H\).
The generic \(P\)-fibre has twelve unramified points
over \(x=\infty\), all regular for \(Q\).
They occur in exactly one of the following two patterns.

- The terminal Laurent chart has \(X\), after integer Laurent
  shifts, of scale \(x^{-k}\) with \(k=D/2\). The restrictions
  of \(P,Q,W\) at the boundary have degrees \(12,42,34\)
  in the scaled variable.
- There is a first split at an integer \(q\) with
  \(2\le q<D/2\). Two Laurent charts then have scale \(x^{-k}\)
  with \(k=D-q\). On each, the restrictions have degrees
  \(6,21,17\), and each accounts for six points.
  
In each pattern the ramification index of \(Q\) at every
one of the twelve outer points on the normalized generic
\(P\)-fibre is \(k-1\).

At the global Euler face, the transferred tops of \(P,Q\)
are powers \(6,21\) of \(x^D X^2\).
The leading \(X\)-coefficient of \(W\) gives weight at least
\(17(D-2)\). Any top of that weight or higher commutes
with \(P\), since the sum of weights is strictly larger
than the weight of \(Q\).
A monomial commuting with \(x^{6D}X^{12}\) has
\(x\)-degree \(D/2\) times its \(X\)-degree.
The bound \(\deg_X W=34\) makes its Euler top a
multiple of \((x^D X^2)^{17}\).
The preceding three-weight lemma now applies to
\(w(x)=1,\ w(X)=-s\).

At a first tie before \(D/2\) the primitive endpoint
is \(x^D X^2\). A tie of smaller \(X\)-degree is at
integer or half-integer \(s\).
We exclude a noninteger half-integer \(s=m/2\), with
\(m\) odd and \(1<s<D/2\), as follows.
The only possible other \(X\)-degree is zero,
so the primitive face is
\[
                      R=x^D X^2(1+z),\qquad
                      z=\gamma x^{-m}X^{-2},\quad \gamma\ne0.
\]
There is a homogeneous \(A\in k[x,x^{-1},X]\) with
\(\{A,R\}=cR^N\), \(c\ne0,\ N\ge0\).
For completeness, polynomials of bounded degree in \(P,Q\)
give quadratically many independent elements and only
linearly many possible initial weights, by bounding
both ends of their supports. At fixed weight the
homogeneous centralizer of \(R\) has dimension at most
one by the Euler identities. Thus some initial bracket
with \(R\) is nonzero. Iterate \(\{\cdot,P\}\) in
\(k[P,Q]\), where this derivation is locally nilpotent.
At the last noncommuting initial the next initial is
its leading bracket with \(R^6\) and commutes with \(R\).
It is a scalar times an integral power of \(R\),
because \(x^mX^2+\gamma\) has a simple non-coordinate
factor in the Laurent polynomial ring. Dividing by
\(6R^5\) gives the displayed bracket for \(A\).
The power \(N\) is nonnegative because that factor
cannot appear in a denominator of a Laurent polynomial
in \(x\) with polynomial \(X\)-dependence.

Put \(J=D-m>0\).
The weight equation and odd \(m\) give the form
\(A=x^{D(N-1)+1}X^{2N-1}F(z)\) with \(F\) Laurent.
Writing \(r=1+z\), direct differentiation gives, after
rescaling \(A\),
\[
       (D-2)F r+
       \bigl(2(N-1)J+2-m\bigr)zF r'
                  -2J zF' r=r^N .
\]
A negative lowest power \(j\) would have coefficient
multiplier \(D-2-2Jj\ne0\), so \(F\) is polynomial.
The coefficient at degree \(\deg F+1\) has multiplier
\(J(2N-1-2\deg F)\ne0\), giving \(\deg F\le N-1\).
At \(z=-1\), after successively dividing by \(r^j\)
for \(0\le j<N\), the surviving coefficient of
\(F/r^j\) is a nonzero multiple of
\(2(N-1-j)J+2-m\), which is odd.
This forces \(r^N\mid F\), a contradiction.
The argument extends the corresponding
[half-integer exclusion on the planar F7 left side](POLYNOMIAL_WEYL_FORTY_TWO_DIRECT_F7_CERTIFICATE.md#the-necessary-left-subring).

At an integer tie \(s=q<D/2\), the primitive form is
\(x^{D-2q}g(Xx^q)\) with \(g\) quadratic.
If it has a double root, shift by that root times \(x^{-q}\).
This is a Laurent shear preserving the bracket and retaining
the same primitive endpoint immediately to the right.
There are finitely many integer slopes below \(D/2\).
If there is no split, continuation to \(D/2\) gives
weight zero for all three elements. After the accumulated
Laurent shifts, substitute \(X=x^{-D/2}v\) in that chart.
The boundary polynomials have the asserted degrees
from the endpoint, and their coefficients are regular
in \(x^{-1}\).

If the quadratic has distinct roots, follow each of
them by shifting by that root times \(x^{-q}\).
Immediately to the right the primitive endpoint is
\(x^{D-q}X\) on both charts, since each root was simple.
Now a new tie before \(D-q\) must be at an integer
slope, with the other term of \(X\)-degree zero.
Shifting its unique root preserves this endpoint.
There are finitely many such slopes, and the chart
at \(D-q\) has weight zero and the asserted degrees.
The two charts have distinct initial values of
\(Xx^q\) before their shifts, so their six generic
points each are distinct. In the unsplit case the
degree-twelve boundary accounts for all generic points.
In either case, because \(D\) is even, the final scales
are integral, and the distinct roots of the boundary
polynomial at a generic level lift to unramified
points over infinity.

Write \(t=1/x\) in a terminal chart. The original \(X\)
has form of a finite Laurent shift in \(t\) plus \(t^kv\).
The change of the two-form is
\(dX\wedge dx=-t^{k-2}dv\wedge dt\).
Hence on the generic \(P\)-level, where the boundary
\(P_v\) is a unit, \(dQ\) has order \(k-2\) in \(t\).
It follows that \(Q\) is regular with the asserted
ramification index.

## Exact geometric degree and genus

**Corollary PlanarFortyTwoExactDegreeAndGenus.**
For the planar transfer of a hypothetical reduced
\((12,42)\) Weyl pair, the geometric degree is \(8H\).
The generic \(P\)-fibre has \(4H+12\) punctures
and its smooth completion has genus
\[
                              g=6k-5H-11 ,
\]
where \(k\) is the integer from the outer profile.
For example, the smallest remaining leading-factor
degree \(30\) has \(H=2\), geometric degree \(16\),
twenty punctures, and possible genera \(9,15,21,27\).

At every finite position outside the roots of \(h\),
the leading \(X\)-coefficient of \(P\) is a unit,
so all points of its generic fibre are affine and
\(Q\) is regular. The only finite roots are high.
The [high-root pole calculation](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#the-complete-pole-profile-over-a-high-root)
shows exactly two poles of \(Q\) there, of orders \(7,1\),
and two additional points at infinity with finite \(Q\)
from its nonzero simple characteristic root. The chart
at the zero characteristic root stays affine. Indeed
\(P(X,0)\) is quadratic by the first weight bound, so
its critical point has a finite \(X\)-coordinate and
lifts by the implicit function theorem; its critical
series in the first direction parameter is \(O(x^3)\).
This identifies the quadratic chart there with the
finite critical-point chart. The nonzero simple
characteristic and both quadruple-root branches have
\(X\) with a pole. The twelve outer points contribute
no poles. This counts the degree and all punctures.

On the smooth affine generic \(P\)-fibre the Keller
condition gives \(dQ\) nowhere zero.
The nonzero simple characteristic chart also gives no
ramification of \(Q\) on that fibre. Its coordinate
change has unit Jacobian there. The high-root poles
contribute \(6H\) to ramification. The twelve outer
points contribute \(12(k-2)\).
The [generic fibre is geometrically irreducible](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_PLANAR_FIBRES.md#monodromy-of-the-generic-position-projection),
so Riemann–Hurwitz for \(Q\) gives
\(2g-2=-16H+6H+12(k-2)\), as claimed.

The [target cylinder](POLYNOMIAL_WEYL_ORDER_TWELVE_FORTY_TWO_TARGET_CYLINDER.md#a-target-chart-with-the-residue-variable)
excludes every split terminal profile. It leaves
\(k=M=5H/2\) and genus \(10H-11\).
