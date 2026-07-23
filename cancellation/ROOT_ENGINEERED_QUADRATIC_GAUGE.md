# Root-engineered quadratic gauges

This note isolates a second polynomial gauge of the `(m,r)=(1,1)`
controlled-boundary cancellation core.  It explains the four-point map

\[
 t=1+xy,\qquad q=t^2z+y^2(1+3t),
\]

\[
 \left(
 tq,\;
 3y-11tq+9xq-2t^2x^2q^4,\;
 6x-9x^2y-3x^3z+x^4q^4
 \right),                                                \tag{1}
\]

proves that its displayed four-point collision is the complete fiber, and
extends it to a normalized all-degree family.  The main new feature is a
quadratic-tilt inverse pencil.  The weighted tangent construction has a
linear tilt `H(W)-sW+t`; here the inverse equation has the form

\[
 H_P(S)-\frac{g_1}{2}(BS^2+C).                           \tag{2}
\]

The derivative of (2) is again exactly the reconstruction coordinate.  Thus
the inverse pencil, ramification divisor, and Jacobian cancellation remain
three aspects of one incidence.

## What the second construction adds

This is not merely another polynomial map with a collision.

1. **A genuinely different construction mechanism.**  The weighted family
   is a polynomial weighted suspension of the linear-tilt plane incidence
   \[
    H(W)-sW+t=0,
   \]
   whose affine discriminant normalization is `A^1`.  The present family
   uses the quadratic-tilt incidence
   \[
    G_P(S)-\frac{g_1}{2}(BS^2+C)=0,
   \]
   together with a birational source chart whose Jacobian is the reciprocal
   of the plane Jacobian factor.  Its affine discriminant normalization is
   `G_m`.  Both constructions realize the same broad cancellation principle,
   but their plane incidences and suspension geometries are different.

2. **A programmable complete fiber.**  For any admissible squarefree rooted
   seed
   \[
    G(S)=g_1S+\cdots+g_NS^N,\qquad g_1g_3g_N\ne0,
   \]
   the complete fiber over `(1,0,0)` is exactly the root set of `G`, with one
   source point reconstructed from each root.  The reconstruction scales are
   the barycentric weights `g_1/G'(r_i)` and obey the standard Lagrange
   balance relations.  Choosing `G` therefore gives direct control over
   rational or number-field fibers, symmetries, affine balance relations,
   heights, and the arithmetic of the root configuration.  Section 8 explains
   why every split squarefree configuration admits an admissible marking
   after rerooting.

3. **The foundational map as the first seed.**  The choice
   `G(S)=S^3+S` recovers the foundational counterexample exactly.  The
   quartic is therefore the next member of one explicit all-degree extension,
   not an unrelated example.

4. **A stable discriminator.**  In fixed degree `N>=4`, the weighted
   and quadratic mechanisms have different affine discriminant
   normalizations, `A^1` and `G_m`.  The complete canonical-boundary
   calculation promotes this distinction to a stable left--right invariant:
   after deleting the intrinsically selected second boundary vertex, the
   full normalized strata are `A^1 x G_m` and `G_m^2`, whose unit ranks are
   one and two.  Thus no quadratic-gauge map is stably polynomially
   left--right equivalent to a boundary-clean weighted map of the same
   degree.

## 1. The short chart for the four-point map

On `t!=0`, put

\[
 P=tq,\qquad S=\frac{x}{t},\qquad Q=y+xq.
\]

Then

\[
 Q-PS=y,\qquad
 D:=1-S(Q-PS)=\frac1t,                                   \tag{3}
\]

and direct substitution turns (1) into

\[
 \begin{aligned}
 A&=P,\\
 B&=3Q-11P+6PS-2P^4S^2,\\
 C&=6S-3QS^2+P^4S^4.                                    \tag{4}
 \end{aligned}
\]

The reciprocal chart and controlled-divisor factors are

\[
 \det\frac{\partial(P,S,Q)}{\partial(x,y,z)}=t=D^{-1}
                                                                  \tag{5}
\]

and

\[
 \det\frac{\partial(B,C)}{\partial(S,Q)}
 =-18(1-SQ+PS^2)=-18D.                                  \tag{6}
\]

These are the degree-four quadratic-gauge inputs to the
[boundary-cancelled incidence lemma](CONTROLLED_BOUNDARY_SUSPENSIONS.md#33-root-engineered-quadratic-gauge).
Here `D=t^{-1}`, the polynomiality condition is already visible in (1), and
the lemma gives

\[
 \boxed{\det DF=-18}
\]

everywhere, without a separate three-variable Jacobian expansion.

Eliminating `Q` from (4) gives the inverse equation over a target `(P,B,C)`:

\[
 E_{P,B,C}(S)
 =P^4S^4-6PS^3+(B+11P)S^2-6S+C=0.                       \tag{7}
\]

On the incidence (4),

\[
 \boxed{\partial_SE_{P,B,C}(S)=-6D.}                    \tag{8}
\]

At `(P,B,C)=(1,0,0)`,

\[
 E_{1,0,0}(S)=S(S-1)(S-2)(S-3).                         \tag{9}
\]

The first target coordinate is one, so every point in this fiber has
`tq=1` and lies in the chart `tq!=0`.  Each of the four roots in (9) is
simple, (8) gives `D!=0`, and the reconstruction formulas

\[
 Q=\frac{B+11P-6PS+2P^4S^2}{3},\qquad
 t=D^{-1},\qquad x=\frac{S}{D},
\]

\[
 y=Q-PS,\qquad q=PD,\qquad
 z=\frac{q-y^2(1+3t)}{t^2}                               \tag{10}
\]

give exactly one source point per root.  They are

\[
 \left(0,\frac{11}{3},-\frac{475}{9}\right),\quad
 \left(-3,\frac43,\frac{125}{81}\right),\quad
 \left(6,\frac13,-\frac7{81}\right),\quad
 \left(-3,\frac23,-\frac19\right).                       \tag{11}
\]

Thus (11) is the complete fiber, not merely a collision certificate, and
the geometric degree is four.

## 2. The quadratic gauge of the cancellation core

Keep

\[
 D=1-SQ+PS^2
\]

and let `a!=0`.  For an arbitrary polynomial `b(P,S)`, define `R(P,S)` by

\[
 R(P,0)=0,\qquad
 R_S+S^2b_S=2a(1+PS^2).                                  \tag{12}
\]

Equivalently,

\[
 R(P,S)=2aS+\frac{2a}{3}PS^3
        -\int_0^S U^2b_U(P,U)\,dU.                       \tag{13}
\]

Set

\[
 \mathcal B=aQ+b(P,S),\qquad
 \mathcal C=R(P,S)-aQS^2.                                \tag{14}
\]

Then

\[
 \begin{aligned}
 \det\frac{\partial(\mathcal B,\mathcal C)}
          {\partial(S,Q)}
 &=-aS^2b_S-aR_S+2a^2QS\\
 &=-2a^2(1-SQ+PS^2)\\
 &=-2a^2D.                                               \tag{15}
 \end{aligned}
\]

Together with the reciprocal source chart (5), this is the ledger in
Lemma 1.2 and gives the constant three-dimensional determinant `-2a^2`.

This is an exact gauge of the standard `(m,r)=(1,1)` plane core.  Indeed put

\[
U=aQ,\qquad
V=2a\int_0^S(1-TQ+PT^2)\,dT
   =2aS-aQS^2+\frac{2a}{3}PS^3.
\]

If

\[
 K(P,S)=\int_0^S T^2b_T(P,T)\,dT,
\]

then (14) is simply

\[
 \boxed{\mathcal B=U+b,\qquad \mathcal C=V-K.}            \tag{16}
\]

The paired corrections in (16) preserve the same relative area form:

\[
 d\mathcal B\wedge d\mathcal C
 =dU\wedge dV
 =-2a^2D\,dS\wedge dQ
\]

at fixed `P`.  They are not generally a target automorphism because `b` and
`K` depend on the marked root `S`.  They are instead a new polynomial gauge
of the controlled-boundary plane map.

The exponent two in (14) is forced within this affine-linear-in-`Q` ansatz.
If one replaces `S^2` by a polynomial `phi(S)`, the coefficient of `Q` in
the plane Jacobian is `a^2 phi'(S)Q`.  Matching the `-SQ` term of `D` forces
`phi'(S)` to be proportional to `S`; after normalization and removal of an
irrelevant constant, `phi(S)=S^2`.

## 3. Engineering a prescribed root polynomial

Let

\[
 G(S)=g_1S+g_2S^2+g_3S^3+\cdots+g_NS^N
 \in k[S],\qquad g_1g_3g_N\ne0.                          \tag{17}
\]

Put

\[
 a=-\frac{g_1}{2},\qquad
 b_0=-g_2,\qquad
 b_1=\frac{g_1-3g_3}{2},\qquad
 b_{j}=-\frac{j+2}{2}g_{j+2}\quad(j\ge2),                \tag{18}
\]

or, at `P=1`,

\[
 b(1,S)
 =-\frac{G'(S)+2a(1+S^2)}{2S}.                           \tag{19}
\]

The numerator in (19) is divisible by `S` because `2a=-g_1`.  Lift (19) by

\[
 b(P,S)=b_0P+b_1PS+
 \sum_{j=2}^{N-2}b_jP^{j+2}S^j.                         \tag{20}
\]

With `R` from (12), one obtains the useful coefficient-weight identity

\[
 \boxed{
 -R(P,S)-b(P,S)S^2
 =g_1S+P(g_2S^2+g_3S^3)
   +\sum_{k=4}^Ng_kP^kS^k.
 }                                                        \tag{21}
\]

Denote the right side by `G_P(S)`.  Eliminating `Q` from (14) gives

\[
 G_P(S)+\mathcal B S^2+\mathcal C=0.                    \tag{22}
\]

After scaling the last two target coordinates by `1/a=-2/g_1`, the
normalized inverse pencil is

\[
 \boxed{
 E_{P,B,C}(S)
 =G_P(S)-\frac{g_1}{2}(BS^2+C).
 }                                                        \tag{23}
\]

On its source incidence,

\[
 \boxed{\partial_SE_{P,B,C}(S)=g_1D.}                   \tag{24}
\]

At `(P,B,C)=(1,0,0)`, equation (23) is exactly `G(S)=0`.

### Tangent-line duality: the conceptual reason

Define the parametrized plane curve

\[
 X(S)=S^2,\qquad
 Y_P(S)=\frac{2G_P(S)}{g_1}.                             \tag{24a}
\]

A target pair `(B,C)` represents the affine line

\[
 Y=BX+C.
\]

Equation (23) says precisely that this line passes through the marked point
`(X(S),Y_P(S))`.  Thus the second target coordinate is not an ad hoc
correction:

\[
 \boxed{C=Y_P(S)-BX(S).}                                \tag{24b}
\]

The first target coordinate is the slope after a triangular source gauge.
Indeed put

\[
 \beta(P,S)
 =\frac{G_P'(S)/g_1-1-PS^2}{S}.                         \tag{24c}
\]

The numerator vanishes at `S=0`, and (21) makes `beta` polynomial.  Formula
(26) in the rational chart is exactly

\[
 B=Q+\beta(P,S).
\]

Consequently

\[
 \frac{G_P'(S)}{g_1}-BS
 =1-SQ+PS^2
 =D.                                                     \tag{24d}
\]

At fixed `P`, use `(S,B)` as source coordinates.  Differentiating (24b)
gives

\[
 \boxed{dC+S^2\,dB=2D\,dS,}
\]

and therefore

\[
 \det\frac{\partial(B,C)}{\partial(S,B)}=-2D.            \tag{24e}
\]

Thus `D` has a geometric meaning: it is half the transversality of the
line to the parametrized curve.  The equation `D=0` says that the line is
tangent.  In particular, the repeated-root discriminant of (23) is the
affine dual curve of `(X(S),Y_P(S))`.  Restricting the boxed one-form
identity to that dual curve immediately gives `dC/dB=-S^2`, anticipating
(36).

The triangular change `(S,Q)->(S,B)` has determinant one.  Hence (24e) is
the normalized form of (15).  Combining it with (5) gives the entire Keller
identity as the zero--pole cancellation

\[
 dP\wedge dB\wedge dC
 =-2D\,dP\wedge dS\wedge dQ
 =-2\,dx\wedge dy\wedge dz.                             \tag{24f}
\]

More generally, for any parametrized curve `(X(S),Y(S))`, the marked-line
map

\[
 (S,B)\longmapsto(B,Y(S)-BX(S))
\]

has Jacobian `-(Y'-BX')`, the derivative of its line-incidence equation.
A constant-Jacobian suspension results whenever that derivative is a
controlled divisor `D` and a second chart contributes `D^{-1}`.  The
weighted tangent core uses the graph coordinate `X=W`; the present family
uses the ramified horizontal coordinate `X=S^2`.  This is the geometric
reason that the two constructions have parallel inverse, discriminant, and
reconstruction formulas.

The general implication in the preceding paragraph is now Lemma 1.2 rather
than a family-local criterion. For the present application its design data
are

\[
 D_0=1+PS^2,\qquad \kappa=S,\qquad
 \lambda=2,\qquad X=S^2,\qquad Y=2G_P/g_1.             \tag{24i}
\]

The common lemma separates construction into three independent design
problems: find a reciprocal-Jacobian chart, integrate its `Q`-coefficient to
obtain the horizontal curve coordinate `X`, and choose `Y` so that the
resulting slope correction and its pullback are polynomial. Equations
(18)--(21) solve the third problem while prescribing the fibre polynomial.

## 4. The denominator-free all-degree map

The rational plane gauge becomes polynomial on affine three-space after one
forced source choice.  Put

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t).                       \tag{25}
\]

Then define

\[
 \boxed{
 \begin{aligned}
 F_{G,1}={}&tq,\\
 F_{G,2}={}&y+3\frac{g_3}{g_1}xq
       +2\frac{g_2}{g_1}tq
       +\sum_{k=4}^N k\frac{g_k}{g_1}
          t^2x^{k-2}q^k,\\
 F_{G,3}={}&x(5-3t)-\frac{g_3}{g_1}x^3z
       -\sum_{k=4}^N(k-2)\frac{g_k}{g_1}(xq)^k.
 \end{aligned}
 }                                                        \tag{26}
\]

Formula (26) is the polynomial algebraization required by the common lemma.
Consequently

\[
 \boxed{\det DF_G=-2.}                                   \tag{27}
\]

It also follows directly from (23)--(24) that:

1. the generic degree of `F_G` is `N`;
2. simple roots of (23) reconstruct uniquely;
3. repeated roots are exactly the reconstruction poles `D=0`;
4. if `G` is squarefree, the complete fiber over `(1,0,0)` consists of the
   roots of `G`.

For a root `r` of a squarefree `G`, the reconstruction is especially short.
One has

\[
 D_r=\frac{G'(r)}{g_1},\qquad
 t_r=D_r^{-1},\qquad x_r=\frac r{D_r},\qquad q_r=D_r,     \tag{28}
\]

and in fact all remaining coordinates have the closed form

\[
 \boxed{
 y_r=\frac{1-D_r}{r},\qquad
 z_r=D_r^3-\frac{g_1}{g_3}D_r(D_r+3)y_r^2.
 }                                                        \tag{28a}
\]

The first quotient is regular at the marked root `r=0`; its value there is
`-2g_2/g_1`.  Formula (28a) follows from
`D=1-rQ+r^2`, `Q=y+r`, and (25).

These reconstruction scales are classical barycentric weights.  If

\[
 G(S)=g_N\prod_{i=1}^N(S-r_i),
\]

then `t_i=g_1/G'(r_i)`, and Lagrange interpolation gives the moment
relations

\[
 \boxed{
 \sum_i t_i r_i^m=0\quad(0\le m\le N-2),\qquad
 \sum_i t_i r_i^{N-1}=\frac{g_1}{g_N}.
 }                                                        \tag{28b}
\]

In particular, the reconstructed fiber satisfies
`\sum_i t_i=\sum_i x_i=0`.  Thus the source points are not an unrelated
list: their chart scales form the unique barycentric balance attached to
the prescribed root configuration.

Taking

\[
 G_N(S)=\prod_{j=0}^{N-1}(S-j)                           \tag{29}
\]

works over `Q` for every `N>=3`: its coefficients `g_1` and `g_3` are
nonzero and all roots are simple integers.  Hence (26) gives, uniformly in
every degree `N>=3`, a determinant-`-2` polynomial map with a complete
`N`-point rational fiber.  For `N>=4`, its component degrees are

\[
 \deg(F_{G,1},F_{G,2},F_{G,3})=(7,\,6N+2,\,6N).           \tag{30}
\]

For this seed the barycentric pattern becomes especially concrete:

\[
 \boxed{
 t_i=(-1)^i\binom{N-1}{i},\quad
 D_i=q_i=\frac{(-1)^i}{\binom{N-1}{i}},\quad
 x_i=i(-1)^i\binom{N-1}{i}.
 }                                                        \tag{30a}
\]

Also `y_0=2H_{N-1}`, where `H_{N-1}` is the harmonic number, and
`y_i=(1-D_i)/i` for `i>0`.  In degree four, (30a) gives
`t=(1,-3,3,-1)` and `x=(0,-3,6,-3)`, exactly the `x`-coordinate pattern
of the four source points in (2).

This is a second, denominator-free route to the existence statements in the
[all-degree rational-fiber theorem](../verified/ALL_DEGREE_RATIONAL_FIBERS.md)
and the
[geometric-degree spectrum](../verified/GEOMETRIC_DEGREE_SPECTRUM.md).
This construction favors a compact formula and freely prescribed fiber over
low coordinate degree.  The weighted family remains substantially smaller
in degree.

## 5. The foundational cubic is the first seed

Take

\[
 G(S)=S^3+S.
\]

Then `g_1=g_3=1`, `g_2=0`, the sums in (26) are empty, and

\[
 q=t^2z+y^2(1+3t).
\]

Formula (26) becomes

\[
 \begin{aligned}
 F_{G,1}&=tq
 =t^3z+y^2t(1+3t),\\
 F_{G,2}&=y+3xq,\\
 F_{G,3}&=x(5-3t)-x^3z.
 \end{aligned}
\]

Since `1+3t=4+3xy` and `x(5-3t)=2x-3x^2y`, this is exactly the
foundational map.  Thus the root-engineered family is an all-degree
extension of the foundational cubic, not merely a family sharing its
Jacobian mechanism.

The four-point map (1) is the unnormalized specialization

\[
 G(S)=S(S-1)(S-2)(S-3)
 =S^4-6S^3+11S^2-6S.                                    \tag{31}
\]

Its last two coordinates are three times the normalized coordinates in
(26), explaining `-18=3^2(-2)`.

## 6. A smaller-coefficient quartic

The seed

\[
 G(S)=S(S-1)(S+1)(S-2)
 =S^4-2S^3-S^2+2S                                       \tag{32}
\]

gives

\[
 t=1+xy,\qquad q=t^2z-y^2(1+3t),
\]

\[
 \boxed{
 \widetilde F=
 \left(
 tq,\;
 -y+3xq+tq-2t^2x^2q^4,\;
 -2x+3x^2y-x^3z+x^4q^4
 \right).
 }                                                        \tag{33}
\]

Here `det D\widetilde F=-2`, and the complete fiber over `(1,0,0)` is

\[
 (0,1,5),\quad
 (-1,2,-9),\quad
 \left(\frac13,-4,-27\right),\quad
 \left(\frac23,-1,45\right).                             \tag{34}
\]

The four source points correspond respectively to the marked roots
`0,1,-1,2`.

## 7. Discriminant normalization and monodromy

Fix `P!=0` and abbreviate `H(S)=G_P(S)`.  A repeated root `r` of (23)
satisfies

\[
 B(r)=\frac{H'(r)}{g_1r},\qquad
 C(r)=\frac{2H(r)-rH'(r)}{g_1}.                          \tag{35}
\]

These are exactly the slope and intercept of the tangent line to

\[
 \left(X(r),Y(r)\right)
 =\left(r^2,\frac{2H(r)}{g_1}\right).
\]

Since `H'(0)=g_1`, one has `r!=0` on the critical divisor.  Differentiating
(35) gives the weighted Legendre relation

\[
 \boxed{\frac{dC}{dB}=-r^2.}                             \tag{36}
\]

Thus (36) is the standard envelope identity `dC+X\,dB=0` for the affine
dual curve, rather than a separate coincidence.  It follows that `r^2`
lies in the function field of the discriminant.  That
field cannot be `k(r^2)`, because `B(r)` has a pole of odd order one at
`r=0`.  Hence `r` itself lies in the discriminant function field, and (35)
is birational.  Both missing points `r=0,infinity` map to infinity, so the
affine reduced discriminant has normalization `G_m`, rather than the `A^1`
normalization of the weighted linear-tilt pencil.

The inverse polynomial (23) is irreducible over `k(P,B,C)` because it has
degree one in `C` with unit coefficient.  The map `r -> B(r)` has degree
`N-1`: its pole orders are one at zero and `N-2` at infinity.  Away from
finitely many values of `B`, its `N-1` critical points are simple, and
birationality of (35) makes their critical values distinct.  A generic
vertical `C`-line is therefore the cover of a Morse polynomial.  The
classical branch-cycle theorem, in the same form used by the
[universal linear-tilt proof](../verified/UNIVERSAL_SYMMETRIC_MONODROMY.md),
gives

\[
 \boxed{\operatorname{Mon}_{\rm geom}(F_G)
       =\operatorname{Mon}_{\rm arith}(F_G)=S_N.}         \tag{37}
\]

Thus the quadratic-tilt family has the same universal symmetric monodromy as
the weighted linear-tilt family, but a different affine discriminant
normalization and a birational rather than polynomial suspension chart.

## 8. Further structural patterns

Several redundancies and linearities become transparent in (26).

### Projective seed

Multiplying `G` by a nonzero scalar changes none of the ratios in (25)--(26).
The normalized map depends on the projective class of the seed polynomial.

The condition `g_3!=0` is only a condition on the chosen marked root.  If a
squarefree degree-`N` polynomial is split over `k`, some root `rho` satisfies
`G'''(rho)!=0`: otherwise the degree-`N-3` polynomial `G'''` would vanish at
all `N` roots.  Rerooting by `S -> S+rho` then gives both
`g_1=G'(rho)!=0` and `g_3=G'''(rho)/6!=0`.  Thus every split squarefree root
configuration admits at least one quadratic-gauge marking.

### The quadratic coefficient is a target shear

The coefficient `g_2` appears only as

\[
 2(g_2/g_1)F_{G,1}
\]

in the second coordinate.  It is removed by a triangular target
automorphism.  Equivalently, in (23) the term `Pg_2S^2` is absorbed by
translating `B`.

### Higher coefficients are additive decorations

For `k>=4`, adding `lambda S^k` to the seed does not change `t`, `q`, or the
first coordinate.  It adds exactly the paired terms

\[
 \Delta F_2=
 k\frac{\lambda}{g_1}t^2x^{k-2}q^k,\qquad
 \Delta F_3=
 -(k-2)\frac{\lambda}{g_1}(xq)^k.                        \tag{38}
\]

The coefficient ratio `k:(-(k-2))` is forced by the differential equation
(12).  The `q^4` pair in (1) is the first decoration.

### The coefficient-weight jump is minimal term by term

Consider a prospective term `cP^alpha S^k` in `G_P`.  Its contribution to
the slope correction (24c) is

\[
 kcP^\alpha S^{k-2},
\]

while its net contribution to the intercept (24b) is

\[
 (2-k)cP^\alpha S^k.                                    \tag{38a}
\]

Under `P=tq`, `S=x/t`, these have `t`-orders
`alpha-k+2` and `alpha-k`.  For an isolated term with `k>=4`,
polynomiality therefore forces

\[
 \alpha\ge k.
\]

The choice `alpha=k` in (21) is the minimal diagonal monomial lift.  Degree
two is exceptional because its intercept contribution vanishes identically.
Degrees one and three form the coupled cubic skeleton whose apparent poles
cancel through (25).  Thus the jump

\[
 (w_1,w_2,w_3,w_4,w_5,\ldots)=(0,1,1,4,5,\ldots)
\]

is explained by tangent-line geometry and boundary regularity.  Lower-degree
families, if they exist in this chart, must use cancellations between
several coefficient terms rather than improving one monomial weight at a
time.

### The source jet is forced

If (25) is replaced by

\[
 q=t^2z+y^2h(t),
\]

polynomiality of the low part of the quadratic gauge forces

\[
 h(t)\equiv\frac{g_1}{g_3}(1+3t)\pmod{t^2}.              \tag{39}
\]

Adding `t^2k(t)` to `h` is removed by the polynomial source shear
`z -> z+y^2k(t)`.  Thus the ratio `g_1/g_3` and the jet `(1,3)` are the
essential source data.

### Relation to the two established cores

The three incidence mechanisms now line up as follows.

| construction | inverse tilt | discriminant normalization | suspension chart |
|---|---|---|---|
| weighted tangent core | `H(W)-sW+t` | `A^1` | polynomial vertical maps |
| standard cancellation core | integral in `S` | controlled boundary | birational |
| quadratic gauge | `H_P(S)-(g_1/2)(BS^2+C)` | `G_m` | birational |

The weighted and quadratic pencils both identify the derivative of the
inverse equation with the reconstruction coordinate.  The standard and
quadratic cancellation cores both cancel the plane Jacobian divisor against
the reciprocal Jacobian of the chart (5).

## 9. Stable separation from the weighted family

The fixed-`P` normalization `G_m` from Section 7 is only the first indication
of the distinction: a target automorphism need not preserve the displayed
slice `P=1`.  The intrinsic comparison uses the complete Zariski--Main
normalization boundary.

For the quadratic family, the only target images of canonical boundary
divisors are the reduced discriminant `Z_Delta` and `Z_0=V(P)`.  Over
`Z_Delta` there is one boundary prime with `(e,f)=(2,1)`.  Put

\[
 d=N-3,\qquad h=\gcd(d,2).
\]

Over `Z_0`, the inverse Newton polygon has slopes

\[
 0,\qquad 1,\qquad\frac{N-1}{N-3}.
\]

The first block is the residue-degree-two affine divisor `q=0`, the second
is the residue-degree-one affine divisor `t=0`, and the last gives `h`
boundary primes, each with

\[
 (e,f)=\left(\frac d h,1\right).
\]

The degree sum

\[
 2+1+h\frac d h=N
\]

proves exhaustion.  These labels intrinsically order `(Z_Delta,Z_0)`.
For `N>=6` the second-vertex ledger already differs from the weighted
ledger of `N-3` unramified boundary primes.

For every `N>=4`, delete the intrinsically selected second vertex and
normalize the ramified target stratum.  The weighted intrinsic-boundary
theorem and the quadratic calculation give

\[
 \begin{aligned}
 \operatorname{Norm}(Z_\Delta^{\rm wt}\setminus Z_0^{\rm wt})
   &\simeq\mathbb A^1\times\mathbb G_m,\\
 \operatorname{Norm}(Z_\Delta^{\rm quad}\setminus Z_0^{\rm quad})
   &\simeq\mathbb G_m^2.
 \end{aligned}
\]

After adjoining any number `s` of identity variables, their unit groups are

\[
 k^\times\xi^{\mathbb Z}
 \quad\text{and}\quad
 k^\times P^{\mathbb Z}r^{\mathbb Z},
\]

of ranks one and two modulo `k^\times`.  Stable normalization functoriality
preserves the ordered strata and stabilization adds only an `A^s` factor, so
the unit-rank gap cannot disappear.  Therefore:

\[
 \boxed{\text{No quadratic-gauge map is stably polynomially
 left--right equivalent to a boundary-clean weighted map.}}
\]

The complete boundary proof, Laurent normalization calculation, and precise
scope of the resulting “two components” statement are in the
[stable-separation theorem](../verified/QUADRATIC_WEIGHTED_STABLE_SEPARATION.md).

## 10. Stable moduli inside the quadratic family

The first moduli count suggested by the fixed-`P` pencil is `N-3`, but the
full three-dimensional family has one additional equivalence.  After dividing
by `g_1` and shearing away `g_2`, put `a_j=g_j/g_1`.  The source and target
scalings

\[
 (x,y,z)\longmapsto(\alpha x,\alpha^{-1}y,\beta z),
\qquad
 (P,B,C)\longmapsto(\beta P,\alpha^{-1}B,\alpha C)
\]

identify the seeds

\[
 a_3\longmapsto\alpha^{-2}\beta^{-1}a_3,\qquad
 a_j\longmapsto\alpha^{1-j}\beta^{-j}a_j\quad(j\ge4).
\]

The `beta` action is independent of root scaling because the cubic skeleton
has `P`-weight one while the higher decorations have `P`-weight `j`.

On the coefficient-torus locus `a_3a_4\cdots a_N!=0`, the intrinsic
normalization Fitting ideal is

\[
 \boxed{
 J=-1+3a_3Pr^2+\sum_{j=4}^Nj(j-2)a_jP^jr^{j-1}.
 }
\]

Its Laurent support orders the two punctures, excludes a residual
`r -> P^m r` ambiguity, and recovers all coefficients modulo exactly the
two displayed scalings.  Hence stable equivalence and ordinary
left--right equivalence have the same orbits on this locus, and

\[
 \boxed{\dim\mathcal M_N^{\rm quad}=N-4.}
\]

The degree-four coarse quotient is one point, with stack stabilizer `mu_5`;
for `N>=5` the quotient is a torus of dimension `N-4`.  Thus the tempting
`N-3` theorem is false: the weighted and quadratic loci are stably disjoint,
but their dimensions are `N-3` and `N-4`, respectively.  The exact
classification is proved in the
[quadratic-gauge stable-moduli theorem](../verified/QUADRATIC_GAUGE_STABLE_MODULI.md).

The comparison with the standard cancellation family is also complete.
Although the stationary-point ladder is the `m=1` cancellation column, only
its first rung belongs to the root-engineered quadratic-gauge family up to
stable polynomial left--right equivalence.  The foundational cubic is the
unique common stable class.  For every `N>=4`, the cancellation
ramified-stratum Fitting generator is binomial while the quadratic generator
has two-dimensional Laurent support; independently, their thick
boundary-intersection nilpotency indices are `mr(m+1)` and `2`.  See the
[quadratic-gauge/cancellation stable-intersection theorem](../verified/QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md).

## 11. Incidence suspensions through horizontal degree four

The marked-line specialization of the
[boundary-cancelled incidence lemma](CONTROLLED_BOUNDARY_SUSPENSIONS.md#33-root-engineered-quadratic-gauge)
is an equivalence inside the coordinate-preserving ansatz: for a fixed
reciprocal chart, `X_S=lambda*kappa` is forced, and the remaining equation
for `Y` is linear.

There is also a complete bounded rechart calculation.  Preserve the marked
root and the `P`-fibration, and consider every determinant-one affine
rechart

\[
 \widetilde P=A(S)P+E(S),\qquad
 \widetilde Q=A(S)^{-1}Q+H(\widetilde P,S).
\]

Then

\[
 \widetilde\kappa=SA(S),\qquad
 D_0=
 1+SAH+\frac{(\widetilde P-E)S^2}{A}.
\]

Polynomiality of `D_0` forces `A|S^2`.  Up to scale, the complete list is

\[
\begin{array}{c|c|c}
A&\widetilde\kappa&X\\ \hline
1&S&S^2\\
S&S^2&S^3\\
S^2&S^3&S^4.
\end{array}
\]

The cubic and quartic cases are genuine reciprocal rational charts, but the
first incidence coordinate contains respectively `Q/S` and `Q/S^2`.
Neither the chart shear `H` nor the incidence correction `beta` depends on
`Q`, so these poles cannot cancel.  On the source, `S=x/t` and `Q=y+xq`;
the poles occur generically along `x=0`.  Hence only `X=S^2` produces a
polynomial Keller map in this entire bounded class.

The weighted case `X=S` remains separate: it is a polynomial vertical
suspension, not a reciprocal rechart of the quadratic source chart.  The
full proof and the precise remaining chart classes are in the
[degree-four incidence-suspension classification](../verified/INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md).

## 12. Next questions

The present formulas suggest several bounded follow-ups.

1. **Non-fibre-affine chart search.**  A new higher-power tilt must change
   the marked root nontrivially, abandon the `P`-fibration, introduce a
   second compensating boundary ledger, or use more than one marked-line
   coefficient.  Search these classes separately; the diagonal reweighting
   route is exhausted.
2. **Degree reduction.**  Formula (30) is much larger than the weighted
   degree profile.  Search for alternate lifts of the coefficient weights
   in (21) that retain polynomiality but lower the `q` powers.
3. **Fiber arithmetic.**  Use the barycentric relations (28b) to choose
   rational root configurations with small reconstructed heights, or to
   impose extra linear relations on the complete fiber.

## 13. Exact verification

Run

```bash
.venv/bin/python scripts/verify_root_engineered_quadratic_gauge.py
.venv/bin/python scripts/verify_quadratic_weighted_stable_separation.py
.venv/bin/python scripts/verify_quadratic_gauge_stable_moduli.py
.venv/bin/python scripts/verify_incidence_suspension_degree_four.py
```

The checker verifies the universal plane determinant, the coefficient-weight
identity, the denominator-free pullback through degree six, the inverse
derivative relation, the foundational specialization, both quartic maps and
their complete fibers, the general discriminant differential identity, and
the consecutive-root seeds through degree twelve.  The stable-separation
checker additionally verifies the two affine `P=0` branches, the complete
Newton ledger and degree sum through degree 64, the high-branch source pole,
and the Laurent recovery identities for the intrinsic `G_m^2` normalization.
The stable-moduli checker verifies the independent two-torus left--right
action, the intrinsic Fitting polynomial and support rigidity, and the exact
`N-4` quotient dimension.  The incidence-suspension checker verifies the
universal determinant criterion, all three degree-four reciprocal recharts,
and the unavoidable cubic and quartic source poles.
