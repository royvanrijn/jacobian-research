# Davenport node-separating affine modifications

The conductor obstruction in the normalized Davenport boundary can be
separated inside an affine plane.  An explicit secant-slope graph is
\(\mathbb A^2\), its Jacobian is one linear boundary equation, and the
derivative divisor pulls back with exactly the expected nodal square.

This does not yet give the absolute Sunada--Keller pair.  One triangular
orientation is impossible in every degree.  The opposite orientation has
one unique torus-clean solution, the parabola \(T+Y^2=0\), whose affine-plane
graph separates the node but deletes the two branches above the coordinate
origin.  Moreover, the puncture-swapping involution cannot act over the node
contraction.  The remaining problem is consequently an alternating
polynomial-coordinate or multi-chart gluing problem, rather than the
previous abstract conductor problem.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0,
\]

and put \(J(T,Y)=g_T'(Y)\).

## 1. An affine-plane chart separating the node

The torus node is

\[
p=\left(-\frac14,-\frac12\right).
\]

Adjoin its secant slope

\[
W=\frac{4T+1}{2Y+1}.                                  \tag{1}
\]

The graph

\[
(2Y+1)W=4T+1                                          \tag{2}
\]

is an affine plane: eliminate \(T\) to obtain

\[
T=\frac{(2Y+1)W-1}{4}.
\]

Thus the graph map is the polynomial morphism

\[
\pi:\mathbb A^2_{Y,W}\longrightarrow\mathbb A^2_{T,Y},
\qquad
(Y,W)\longmapsto
\left(\frac{(2Y+1)W-1}{4},Y\right).                   \tag{3}
\]

Its Jacobian is

\[
\operatorname {Jac}(\pi)=\frac{2Y+1}{4}.              \tag{4}
\]

Exact substitution gives

\[
\pi^*J=(2Y+1)^2S(Y,W),                                \tag{5}
\]

where \(S=0\) is smooth on \(Y\ne0\), hence on the relevant torus.  On the
exceptional line \(Y=-1/2\),

\[
S\left(-\frac12,W\right)
=
\frac{W\bigl((3a+6)W-4a-16\bigr)}{64}.                \tag{6}
\]

The two points are therefore

\[
W_-=0,\qquad W_+=2-\frac{2a}{3},                      \tag{7}
\]

and are distinct.  They are precisely the two tangent directions of the
node.  Equations (4)--(5) bring two previously separate ledgers together:
the affine modification contributes one boundary Jacobian, while the
nodal derivative contributes its square.

This is the requested ambient separation of the two-point conductor, but
only on one affine chart.

## 2. The unavoidable extra point in the minimal chart

The image of (3) is the complement of the line \(Y=-1/2\), together with
the node \(p\).  Restricting \(J\) to that contracted line gives

\[
J\left(T,-\frac12\right)
=-\frac{(4T+1)^2\bigl((12a+20)T-1\bigr)}{64}.         \tag{8}
\]

Besides the double node, there is one further torus point,

\[
\left(
\frac1{12a+20},-\frac12
\right)
=
\left(
\frac1{56}-\frac{3a}{112},-\frac12
\right).                                             \tag{9}
\]

The graph (2) deletes (9).  Accordingly the strict transform is the
node-separated boundary with one extra puncture, not the original
three-puncture normalization used by the boundary involution.

This is the smallest possible defect among affine-linear charts.

## 3. First triangular orientation in every degree

Put

\[
n=4T+1,\qquad u=2Y+1.                                \tag{10}
\]

Every finite-slope triangular polynomial coordinate through the node has
the form

\[
d=u+f(n),\qquad f(0)=0.                              \tag{11}
\]

The graph

\[
dW=n                                                  \tag{12}
\]

is always an affine plane.  Indeed, in coordinates \((d,W)\),

\[
n=dW,\qquad u=d-f(dW).
\]

The Jacobian from \((d,W)\) to \((n,u)\) is \(-d\), hence the Jacobian to
\((T,Y)\) is \(-d/8\).  Thus this is the complete triangular family of
affine-plane node contractions, not merely a bounded-degree ansatz.

Use the normalization parameter \(v\) from the preceding audit.  Define

\[
\begin{aligned}
D(v)&=(v+a)
\left(v^2+(-2+3a)v-8-2a\right),\\
Y_0(v)&=(2-a)v^2-2-11a,\\
z(v)&=\frac{3a-2-v^2}{4a+16}.
\end{aligned}                                        \tag{13}
\]

Then

\[
Y=\frac{Y_0}{D},\qquad
u=\frac{U}{D},\qquad
n=\frac{N}{D^2},                                     \tag{14}
\]

where

\[
U=D+2Y_0,\qquad N=D^2+4zY_0^2.                       \tag{15}
\]

Let \(e=\deg f\ge1\), with leading coefficient \(c_e\).  On the normalized
boundary,

\[
u+f(n)=\frac{H(v)}{D(v)^{2e}},
\]

with

\[
H
=UD^{2e-1}
+\sum_{i=1}^{e}c_iN^iD^{2e-2i}.                     \tag{16}
\]

At each of the three roots of \(D\), this function has pole order \(2e\).
The two node branches are

\[
v_+=a+4,\qquad v_-=-a-4,
\]

and a transverse contraction gives a simple zero at each.  Away from the
torus, the only additional finite zeros are the two branches above the
origin \(Y_0=0\); they occur together and simply when \(1+f(1)=0\).  The
remaining zero multiplicity is at \(v=\infty\).

Consequently, if the graph introduced no residual torus point, its numerator
would have to be one of

\[
H=cR_2,\qquad H=cR_4,
\quad
R_2=(v-v_+)(v-v_-),\quad R_4=R_2Y_0.                 \tag{17}
\]

This converts the apparently unbounded coefficient problem into a
three-puncture norm calculation.  Reducing (16) modulo \(D\), only the
highest term survives, so (17) would force

\[
N^e\equiv c'R_j\pmod D,\qquad j=2\text{ or }4.        \tag{18}
\]

Write

\[
Q(v)=v^2+(-2+3a)v-8-2a
\]

and let \(\beta\) be either root of \(Q\).  Comparing (18) at the rational
puncture \(v=-a\) and at \(\beta\) gives

\[
\left(\frac{N(\beta)}{N(-a)}\right)^e
=\frac{R_j(\beta)}{R_j(-a)}.                          \tag{19}
\]

The quadratic norms to \(K\) are

\[
\begin{aligned}
\operatorname N\!\left(\frac{N(\beta)}{N(-a)}\right)
  &=\frac{343(5a-2)}{64},\\
\operatorname N\!\left(\frac{R_2(\beta)}{R_2(-a)}\right)
  &=\frac{7(a+2)}8,\\
\operatorname N\!\left(\frac{R_4(\beta)}{R_4(-a)}\right)
  &=\frac{49(3a+2)}{32}.
\end{aligned}
\]

Taking the absolute norm \(K/\mathbb Q\) yields

\[
\left(\frac{7^6}{2^6}\right)^e
\in
\left\{
\frac{7^2}{2^4},
\frac{7^4}{2^6}
\right\}.                                            \tag{20}
\]

The \(7\)-adic valuations would give \(6e=2\) or \(6e=4\), impossible for
every integer \(e\ge1\).  Therefore

\[
\boxed{\text{No graph }d=u+f(n)\text{, in any degree, separates the node
without adding a torus puncture.}}
\]

The infinite-slope coordinate \(n=0\) is a node tangent and loses one branch
from its affine chart, so it is not an omitted transverse case.

## 4. Opposite triangular orientation and the unique parabola

It remains to reverse the triangular direction:

\[
d=n+g(u),\qquad g(0)=0.                              \tag{21}
\]

The graph \(dW=u\) is again \(\mathbb A^2\), now with

\[
u=dW,\qquad n=d-g(dW),
\]

and Jacobian \(d/8\) to \((T,Y)\).

For \(e=\deg g\ge3\), the leading polar class at the three punctures is
\(U^e\).  A graph with no residual torus zero can have finite zero divisor
only

\[
R_2Y_0^r,\qquad r\ge0.
\]

The relevant absolute norms are

\[
\operatorname N_{\!K/\mathbb Q}
\operatorname N\!\left(\frac{U(\beta)}{U(-a)}\right)
=\frac{7^2}{2^2},
\]

and

\[
\operatorname N_{\!K/\mathbb Q}
\operatorname N\!\left(
\frac{R_2(\beta)Y_0(\beta)^r}
     {R_2(-a)Y_0(-a)^r}
\right)
=\frac{7^{2+2r}}{2^{4+2r}}.                          \tag{22}
\]

Equality would force \(e=r+1\) from the \(7\)-adic valuation and
\(e=r+2\) from the \(2\)-adic valuation.  Hence no degree \(e\ge3\) works.
Degree one is either covered by the first orientation or is a node tangent.

Degree two is exceptional because \(n\) and \(u^2\) have the same pole
order.  Write

\[
d=n+b_2u^2+b_1u.
\]

Reduction modulo \(D\) shows that the leading class \(N+b_2U^2\) can match
a torus-clean divisor only when

\[
b_2=1,\qquad r=2.
\]

The full polynomial identity then fixes \(b_1=-2\):

\[
N+U^2-2UD
=\frac{a-3}{14}R_2Y_0^2.                             \tag{23}
\]

In the original coordinates,

\[
d=n+u^2-2u=4(T+Y^2).                                 \tag{24}
\]

Thus, up to harmless scaling, the unique torus-clean opposite triangular
chart is

\[
\delta W=2Y+1,\qquad \delta=T+Y^2.                   \tag{25}
\]

It is explicitly \(\mathbb A^2_{\delta,W}\):

\[
Y=\frac{\delta W-1}{2},\qquad T=\delta-Y^2.
\]

Its Jacobian is \(\delta/2\), and exact substitution gives

\[
\pi_{\mathrm{par}}^*J=\delta^2S_{\mathrm{par}},
\]

where \(S_{\mathrm{par}}=0\) is smooth and

\[
S_{\mathrm{par}}(0,W)
=\frac{(W+2)\bigl((a-2)W+6a+12\bigr)}{16}.           \tag{26}
\]

The two exceptional points are distinct.  The contracted parabola meets
the derivative curve by

\[
J(-Y^2,Y)=(a-2)Y^4(2Y+1)^2.                         \tag{27}
\]

Therefore (25) has no residual torus intersection: besides the node, it
deletes only the coordinate origin, with its two normalization branches.
This is strictly better than the secant chart (2), but it does not preserve
the full affine normalized boundary used in the unit-rank comparison.

Together, Sections 3--4 classify both one-shear triangular orientations.
Any remaining one-chart coordinate must have alternating triangular length
at least two in a Jung--van der Kulk decomposition.

## 5. The complementary two-chart gluing

The secant and parabola charts repair one another on the derivative curve.
Use

\[
u=2Y+1,\qquad w=\frac{4T+1}{u}
\]

on the secant chart, and

\[
\delta=T+Y^2,\qquad z=\frac{u}{\delta}
\]

on the parabola chart.  Direct substitution gives

\[
\delta=\frac{u(w+u-2)}4,\qquad
z=\frac4{w+u-2},                                     \tag{28}
\]

with inverse

\[
u=\delta z,\qquad
w=\frac4z-\delta z+2.                                \tag{29}
\]

The extra torus point (9) omitted by the secant chart has
\(\delta=T+1/4\ne0\), so it belongs to the parabola chart.  Conversely, the
origin omitted by the parabola chart has \(u=1\), so it belongs to the
secant chart.  Hence the two strict transforms cover the full normalized
affine derivative boundary.

The ambient gluing is nevertheless not affine.  On the exceptional fiber,
(28) reduces to

\[
z=\frac4{w-2}.                                       \tag{30}
\]

Thus the exceptional affine line in the \(w\)-chart and the exceptional
affine line in the \(z\)-chart glue by reciprocal coordinates.  Together
they form a complete \(\mathbb P^1\).  This \(\mathbb P^1\) is the closed
fiber over the node, so the glued surface cannot be affine.

Therefore the most direct two-chart repair succeeds exactly on the boundary
curve but fails at the ambient level:

\[
\boxed{\text{full normalized boundary}+
\text{natural two-chart gluing}
\Longrightarrow\text{exceptional }\mathbb P^1.}
\]

Any affine two-chart replacement must contract or modify this complete
fiber while retaining its two distinguished tangent points and the Cox
Jacobian ledger.

## 6. Global no-one-chart theorem

The preceding classifications suggest a stronger formulation which no
longer refers to a Jung normal form.

Let \(C\subset\mathbb A^2_{T,Y}\) be any polynomial coordinate line through
the node, parametrized by

\[
t\longmapsto (T(t),Y(t)).
\]

If its affine-plane graph separated the two node branches while preserving
the full affine derivative boundary, then \(C\) would meet \(J=0\) only at
the node and transversely to both branches.  After translating \(t\), this
would force

\[
J(T(t),Y(t))=\kappa t^2,\qquad \kappa\ne0.            \tag{31}
\]

Put

\[
p=\deg T(t),\qquad q=\deg Y(t).
\]

The Newton vertices of \(J\) are \((3,0)\), \((2,0)\), and \((0,6)\).
If \(p>2q\), the \(T^3\) term is the unique highest-degree term after
substitution.  If \(p<2q\), the \(Y^6\) term is unique.  Both contradict
(31), so

\[
p=2q.                                                 \tag{32}
\]

The outer Newton polynomial is

\[
\Phi(z)
=(-3a-5)z^3+(6a-9)z^2+(5a+5)z+1.
\]

Over \(K\), it has one linear factor and one irreducible quadratic factor.
Since the leading coefficient ratio
\(\operatorname {lc}(T)/\operatorname {lc}(Y)^2\) lies in \(K\), it must
equal the unique \(K\)-root

\[
z_0=\frac{1+2a}{7}.                                  \tag{33}
\]

Write \(T=z_0Y^2+X\), so \(\deg X<2q\).  In the transformed polynomial,
the competing top terms are \(XY^4\) and \(Y^5\).  If
\(\deg X>q\), the first is unique; if \(\deg X<q\), the second is unique.
Thus \(\deg X=q\), and cancellation fixes the leading coefficient

\[
b=-\frac{a+4}{7}.
\]

Write \(X=bY+X_1\), with \(\deg X_1<q\).  The next competing terms are
\(X_1Y^4\) and \(Y^4\).  Hence \(X_1\) must be constant, and cancellation
fixes it.  The coordinate curve is forced to be

\[
T=
\frac{1+2a}{7}Y^2
-\frac{a+4}{7}Y
+\frac{a+4}{7}.                                      \tag{34}
\]

But exact substitution gives

\[
J(T,Y)
=\frac{
(12a-8)Y^3+(-30a+20)Y^2+(18a-12)Y-(5a+6)
}{7}.                                                 \tag{35}
\]

This has degree three in \(Y\), so along \(Y(t)\) it has degree \(3q\),
which cannot equal two.  This contradicts (31).

Therefore:

\[
\boxed{\text{No polynomial coordinate line gives a one-chart affine-plane
separation preserving the full derivative boundary.}}
\]

This removes every alternating polydegree, not just the one-shear families
of Sections 3--4.  The remaining construction must genuinely use a
multi-chart modification or a higher-dimensional non-coordinate surgery.

## 7. Why the direct normalization graph is worse

The rational parameter \(v\) from the
[boundary-involution audit](DAVENPORT_BOUNDARY_INVOLUTION.md) has an
ambient representative which simplifies modulo \(J\) to

\[
v=-\frac{L(T,Y)}{Y(T+Y^2)},                           \tag{36}
\]

where

\[
L=(4a-2)TY+(a-2)T+2(a+1)Y^3.
\]

The denominator divisor is the union of the line \(Y=0\) and the parabola
\(T=-Y^2\), of class \(2\mathbb L-1\).  Its reduced center consists of the
origin and the torus node.  The graph of (36) therefore has class

\[
\mathbb L^2-(2\mathbb L-1)+2\mathbb L
=\mathbb L^2+1.                                      \tag{37}
\]

After embedding \(K\) in \(\mathbb C\), its Hodge--Deligne polynomial
differs from that of \(\mathbb A^2\).  The direct normalization graph is
not an affine plane.  The secant construction (2) is genuinely better:
it preserves the ambient affine plane, at the price of the single point
(9).

## 8. Interaction with the puncture involution

The branch parameters from the normalized-boundary audit are

\[
v_+=a+4,\qquad v_-=-a-4.
\]

Under (3), they give the two exceptional points (7).  The unique
puncture-swapping involution fixes \(v_-\), but

\[
\iota(v_+)=-\frac{5a+6}{4}.
\]

The image point has

\[
(T,Y,W)
=
\left(
-\frac14,\frac{a+3}{2},0
\right),                                             \tag{38}
\]

which is not on the exceptional line.  Hence the involution does act on
the abstract normalized boundary after the conductor is separated, but it
cannot lift to an automorphism over the same contraction (3): such a lift
would have to preserve the exceptional fiber.  The same argument applies
to the parabola contraction (25), since its exceptional fiber also maps to
the node.

The old conductor obstruction is therefore replaced by two sharper facts:

1. every one-shear affine-plane contraction either introduces a torus
   puncture or is the unique parabola chart, which deletes the coordinate
   origin;
2. the complementary two-chart gluing restores the full boundary but
   creates a complete exceptional \(\mathbb P^1\);
3. the desired involution is not equivariant for either minimal
   contraction.

## 9. Remaining construction

The next search no longer needs an unspecified “ambient normalization.”
It has two concrete forms:

1. contract or affinely modify the explicit exceptional \(\mathbb P^1\)
   in the complementary gluing while retaining its two tangent points and
   allowing the involution to move the exceptional divisor; or
2. use a higher-dimensional non-coordinate affine modification whose
   exceptional fiber is affine but whose center still realizes the full
   three-puncture normalization and the doubled Cox ledger.

All polynomial-coordinate one-chart constructions are now excluded by
Section 6.  The first remaining option is no longer an unspecified gluing
problem: its transition is (28)--(29), and the precise obstruction is the
complete exceptional fiber (30).

## 10. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_node_separation.py
```

The checker verifies the affine-plane graph, its Jacobian and squared
derivative pullback, smoothness of the torus strict transform, the two
separated tangent directions, the extra point (9), the all-degree
one-shear norm classification, uniqueness and smoothness of the parabola
chart, the complementary transition and exceptional \(\mathbb P^1\), the
global Newton no-one-chart theorem, the class (37), and failure of
equivariance for the puncture involution.
