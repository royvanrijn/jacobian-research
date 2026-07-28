# The `PC(2)` graph route to `HC_4`

## Status

Meng--Yang define `HC_n` by asking whether every polynomial potential with
constant nonzero Hessian determinant has polynomial formal Legendre transform.
Their 24 July 2026 preprint proves `HC_5` false and records

\[
 HC_n\text{ true for }n\leq 3,\qquad
 HC_n\text{ false for }n\geq 5,\qquad
 HC_4\Longrightarrow JC_2.
\]

See [Meng--Yang, arXiv:2607.22198](https://arxiv.org/abs/2607.22198).

The explicit `PC(2)` map

\[
 G=(R,T,D,S):\mathbb A^4\longrightarrow\mathbb A^4
\]

is an exact polynomial symplectic Keller map with a complete rational
three-point fiber.  It is therefore a natural source for a four-variable
gradient Keller map.  The most direct proposal is to view its graph as a
Lagrangian fourfold, choose a linear symplectic polarization, and use the
complementary graph coordinates as the gradient of a generating function.

The exact outcome of the full linear search is negative:

> **Complete linear-polarization obstruction.** No linear Lagrangian
> projection of the graph of `G` has constant nonzero Jacobian.

This eliminates the proposed *linear* graph-polarization route, including
projections that might create collisions unrelated to the known fiber.

The nonlinear escape has also been tested in three all-degree shear classes
and in expanded sparse, dense-random, and two-step searches.  Exact jet ideals
and boundary classifications rule out unrestricted quadratic--cubic
potentials in all 16 coordinate charts.  The complete caustic and global
Schur analysis now also rules out every quartic potential in the principal
chart `1000`; the final 43-parameter affine-normal edge is an exact rational
unit ideal on \(Y=0\).  Quartic potentials in the other unresolved charts
and general higher-degree nonlinear polynomial polarizations remain open.

A direct imitation of the Meng--Yang one-variable Schur operation has also
been carried out.  Among all coordinate graph charts, omitted coordinates,
and linear auxiliary functions of the polynomial graph parameters, only two
polynomial-coordinate generating families survive.  Both preserve the
rational collision, but every polynomial quadratic pivot in either family
has zero or nonconstant descended Hessian determinant.  A second linear
Schur descent of the Meng--Yang five-variable polynomial is likewise
impossible.  Arbitrary polynomial auxiliary coordinates are now reduced to
the same two families; non-coordinate graph embeddings remain open.

The exact linear certificates are
[`scripts/verify_hc4_linear_polarization_obstruction.py`](scripts/verify_hc4_linear_polarization_obstruction.py)
and
[`scripts/verify_hc4_all_linear_projection_obstruction.py`](scripts/verify_hc4_all_linear_projection_obstruction.py).

## 1. The graph and its ambient Darboux coordinates

Use the source bracket convention

\[
 \{p,x\}=\{z,q\}=1
\]

and the target pairs

\[
 \{D,R\}=\{S,T\}=1.
\]

On the graph ambient space, with symplectic form
\(\omega_{\rm target}-\omega_{\rm source}\), take the four Darboux pairs

\[
 Q=(x,q,R,T),\qquad M=(-p,-z,D,S).
\]

The graph map

\[
 u=(x,q,p,z)\longmapsto (Q(u),M(u))
\]

is Lagrangian.  A linear symplectic polarization consists of two Lagrangian
coordinate systems `(A,B)` with

\[
 \Omega=\sum_{i=1}^4 dA_i\wedge dB_i.
\]

If the independent projection \(A|_{\Gamma_G}\) is a polynomial
automorphism, then

\[
 b(A)=B\circ(A|_{\Gamma_G})^{-1}
\]

has symmetric Jacobian because the graph is Lagrangian.  Over characteristic
zero it integrates to a polynomial \(\Psi\) with \(b=\nabla\Psi\).  Moreover

\[
 \det\operatorname{Hess}\Psi
 =\frac{\det d(B|_{\Gamma_G})}{\det d(A|_{\Gamma_G})}.
\]

Since a polynomial automorphism has constant nonzero Jacobian, constant
nonzero Hessian determinant requires

\[
 \det d(B|_{\Gamma_G})\in k^\times.                  \tag{1}
\]

Thus (1) is a necessary gate that can be checked before polynomial inversion
or integration.

## 2. The 16 Lagrangian charts

For each of the four Darboux pairs, either keep `(Q_i,M_i)` or swap it to
`(M_i,-Q_i)`.  A mask

\[
 \epsilon\in\{0,1\}^4
\]

therefore gives chart coordinates `(q_0,m_0)`.  Every linear Lagrangian
subspace lies in at least one of these 16 standard charts.  In a chart it has
the form

\[
 B=m_0+Kq_0,\qquad K=K^T.                            \tag{2}
\]

This is the finite reduction of the full linear-polarization problem: the
continuous part is the ten-parameter symmetric matrix \(K\), and the
collision condition imposes four linear equations on it.

The determinant calculations are made in the polynomial Darboux coordinates

\[
 (X,Y,W,D)
\]

already constructed for `PC(2)`.  The change
\((x,q,p,z)\leftrightarrow(X,Y,W,D)\) is a polynomial automorphism with
constant Jacobian, so constancy and nonvanishing in (1) are unchanged.

## 3. The certified fiber and collision equations

The fiber over \((R,T,D,S)=(0,0,0,-1/8)\) is

\[
\begin{aligned}
 P_0&=\left(0,0,\frac1{24},-\frac18\right),\\
 P_+&=\left(1,\frac23,\frac{247}{96},-\frac{89}{64}\right),\\
 P_-&=\left(-1,-\frac23,\frac{247}{96},-\frac{89}{64}\right).
\end{aligned}
\]

For a selected pair, let \((\Delta Q,\Delta M)\) be its ambient displacement.
The collision survives the complementary projection exactly when

\[
 \Delta B=\Delta m_0+K\Delta q_0=0.                 \tag{3}
\]

The three displacements used by the checker are

\[
\begin{array}{c|c|c}
\text{pair}&\Delta Q&\Delta M\\ \hline
P_+-P_-&(2,4/3,0,0)&(0,0,0,0)\\
P_+-P_0&(1,2/3,0,0)&(-81/32,81/64,0,0)\\
P_--P_0&(-1,-2/3,0,0)&(-81/32,81/64,0,0).
\end{array}
\]

Equation (3) is solved exactly over \(\mathbb Q\) in every chart.

## 4. A hand-checkable obstruction for \(P_+,P_-\)

For the symmetric pair, four charts have inconsistent collision equations.
The remaining 12 charts have six free parameters after solving (3).

Ten charts contain a nonconstant monomial of \(\det dB\) whose coefficient is
nonzero and independent of all six parameters:

\[
\begin{array}{c|c|c}
\epsilon&\text{monomial in }(X,Y,W,D)&\text{forced coefficient}\\ \hline
0000&Y^2D&15\\
0001&X^3D^2&-96\\
0010&X^5W^3&-27/2\\
0011&X^7Y^3D&-972\\
0101&W^2&-32/9\\
0110&X^4&-2\\
0111&X&32/3\\
1001&W^2&-8\\
1010&X^4&-9/2\\
1011&X&24.
\end{array}
\]

So \(\det dB\) cannot be constant in those charts.

In chart `0100`, four coefficient equations would successively force

\[
 f=0,\qquad a=0,\qquad c=1,
\]

after which the coefficient of \(YW^2\) is \(-160/9\), a contradiction.
In chart `1000` they force

\[
 f=0,\qquad a=0,\qquad c=-3/2,
\]

after which the coefficient of \(YW^2\) is \(-40\).

This gives a short exact proof for all 16 charts for the pair \(P_+,P_-\).

## 5. Exact coefficient ideals for the other two pairs

For each of \(P_+-P_0\) and \(P_--P_0\), equation (3) is consistent in all 16
charts and leaves six parameters.  It is unnecessary to expand the entire
determinant obstruction.

Restrict \(\det dB\) to the two slices

\[
 X=0,\qquad Y=W=D=0.
\]

If \(\det dB\) were constant, every nonconstant coefficient on both slices
would vanish.  After denominators and scalar contents are removed, these
coefficients generate an ideal in the six chart parameters.  In every chart
the exact rational Gröbner basis is the unit ideal.  Hence the constancy
conditions are inconsistent.

Together with Section 4, this proves the stated obstruction for every pair in
the known three-point fiber.

## 6. Coordinate polarizations and bounded shear regression

Two smaller searches are retained as transparent regressions.

[`scripts/search_hc4_graph_polarizations.py`](scripts/search_hc4_graph_polarizations.py)
checks the 16 coordinate polarizations directly.  Collision survival leaves
four cases.  Two independent projections have rank at most three; the other
two have nonconstant Jacobians:

\[
\begin{aligned}
\det d(x,q,D,T)
 &=\frac{2X}{3}\left(
 3WX^3Y+3WX^2+9X^2Y^2+3XY-4\right),\\
\det d(x,q,D,S)
 &=\frac16\left(
 3WX^4Y^2+6WX^3Y+3WX^2
 +9X^3Y^3+12X^2Y^2-XY-3\right).
\end{aligned}
\]

[`scripts/search_hc4_lagrangian_shears.py`](scripts/search_hc4_lagrangian_shears.py)
searches the symmetric shear chart

\[
 B=M+KQ,\qquad K(3,2,0,0)^T=0.
\]

It finds no Keller projection in the complete integral box
\([-3,3]^6\), totaling \(117649\) matrices.  The all-chart certificate makes
this bounded negative result logically unnecessary, but it remains useful as
a fast regression and as a template for later nonlinear ansatz searches.

## 7. Consequence for the proposed program

The gates occur in this order:

1. retain a certified collision under \(B\);
2. make \(B|_{\Gamma_G}\) Keller;
3. make \(A|_{\Gamma_G}\) a polynomial automorphism;
4. integrate \(B(A)\) to a polynomial potential.

The linear search fails at gate 2 for every possible linear Lagrangian
complement \(B\), even without imposing a collision.  Thus a linear
projection cannot escape by creating a new collision: it never becomes a
four-variable Keller map in the first place.

The only live extension is genuinely nonlinear.  This is the natural place
to reuse the mechanism that forced the unique shear
\(Z\mapsto Z-9Q^2\).

## 8. Unrestricted linear projections

Drop all collision equations.  In each of the 16 Lagrangian charts write

\[
 B=m_0+Kq_0,\qquad K=K^T,
\]

with ten independent entries in \(K\).  Constancy of \(\det dB\) requires
every nonconstant coefficient to vanish.

Exact coefficient ideals restricted to the three hyperplanes

\[
 X=0,\qquad Y=0,\qquad W=0
\]

are the unit ideal in 14 charts.  In the two remaining charts, `1110` and
`1111`, the radical conditions force

\[
 K=\operatorname{diag}(0,0,0,\ell).
\]

The resulting complementary projection has identically zero Jacobian for
every \(\ell\).  Consequently:

\[
\boxed{\text{No linear Lagrangian projection of }\Gamma_G
\text{ has constant nonzero Jacobian.}}
\]

This completely closes the proposed “linear projection with a new collision”
escape.

## 9. Nonlinear shear obstructions

### 9.1 The natural chart, in every degree

Consider the full polynomial momentum-shear class

\[
 B=M+\nabla V(Q),\qquad Q=(x,q,R,T),
\]

for an arbitrary polynomial \(V\).  Put

\[
\begin{aligned}
A_0={}&3WX^4Y^2+6WX^3Y+3WX^2+9X^3Y^3+12X^2Y^2-XY-3,\\
C_0={}&3WX^3Y+3WX^2+9X^2Y^2+3XY-4.
\end{aligned}
\]

The coefficient of \(D^2\) in \(\det dB\) is

\[
 -6X^2\left(A_0+4XC_0V_{TT}\right).
\]

Constancy therefore forces

\[
 V_{TT}=-\frac{A_0}{4XC_0}
        =-\frac{3}{16X}+O(1).
\]

But \(V_{TT}(Q)\) is polynomial in the graph coordinates.  The forced
negative-\(X\) principal part is an exact obstruction, directly analogous to
the calculation that selected \(Z-9Q^2\).  Unlike that earlier calculation,
there is no free shear coefficient that can cancel this pole.

### 9.2 The two exceptional charts

The linear calculation suggests charts `1110` and `1111`, so test their
smallest nonlinear all-degree subclasses.

In chart `1111`, take \(V=V(D,S)\).  If

\[
 A=V_{DD},\qquad \Delta=\det\operatorname{Hess}_{D,S}V,
\]

then

\[
 \det dB=\frac{\Delta A_0-4AXC_0}{6}.
\]

At \(X=0\), constancy forces \(\Delta\) to be a nonzero constant.  The
remaining equation forces the boundary value

\[
 A=\frac{\Delta Y}{16}.
\]

Yet on this boundary

\[
 S=\frac{W+4Y^2}{2}.
\]

Fixing \(S\) while varying \(Y\) contradicts \(A\) being a function of
\((D,S)\).

In chart `1110`, take \(V=V(D,T)\).  Here

\[
 \det dB=\frac{4XC_0\Delta+AA_0}{6}.
\]

Constancy first forces \(A\) to be constant and then forces
\(\Delta\) to be a multiple of

\[
 \frac{A_0+3}{4XC_0}.
\]

If this expression factored through \(T\), its restriction at \(X=0\), where
\(T=Y\), would force it to equal \(T/16\) identically.  Exact subtraction
gives a nonzero polynomial numerator, so no such factorization exists.

Thus both exceptional two-variable potential classes are obstructed in every
degree.

### 9.3 Sparse cubic search

Two bounded searches were run in the exceptional charts.  The first perturbed
the residual singular linear family:

\[
 V=\frac{\ell}{2}u_4^2+\sum_i c_im_i,
\]

where the \(m_i\) range over all 20 cubic monomials in the chart variables.
One- and two-monomial supports use

\[
 c_i\in\{-2,-1,1,2\},\qquad \ell\in[-2,2].
\]

Three-monomial supports are also included with \(c_i\in\{-1,1\}\).  This
gives \(61200\) potentials in each exceptional chart, \(122400\) total.
A nine-point modular determinant filter found no survivors.

There is a simple reason this family cannot work: every cubic Hessian vanishes
at the origin, while the residual quadratic family is singular there.

The genuinely nonlinear search therefore allowed

\[
 V=V_2+V_3,
\]

where the symmetric Hessian of \(V_2\) has at most two nonzero entries and
\(V_3\) contains one or two cubic monomials.  Every nonzero coefficient is
\(\pm1\).  This gives \(160800\) potentials per chart, \(321600\) total.
There were no modular survivors.

The exhaustive sparse search was then expanded in two independent directions:

\[
\begin{array}{c|c|c}
\text{quadratic support}&\text{cubic support}&\text{total candidates}\\ \hline
\leq3&\leq2&1\,857\,600\\
\leq2&\leq3&3\,987\,840.
\end{array}
\]

Again there were no modular survivors.  A deterministic dense search sampled
another \(320000\) quadratic--cubic potentials, with all 30 coefficients
drawn independently from \([-2,2]\), across all 16 charts.  It also found no
survivors.

### 9.4 Unrestricted quadratic--cubic jet ideals

The bounded searches can be replaced by an exact calculation in several
charts.  Let

\[
 V=V_2+V_3
\]

have all ten quadratic coefficients and all twenty homogeneous-cubic
coefficients free.  Expand \(\det dB\) by total degree at the graph origin.
Set every positive-degree coefficient to zero and adjoin

\[
 zc_0-1,
\]

where \(c_0\) is the constant determinant term.  This saturates by
\(c_0\ne0\).  Exact Gröbner reduction over \(\mathbb Q\) gives:

\[
\begin{array}{c|c}
\text{chart}&\text{unit ideal reached by}\\ \hline
0000&\text{all-degree principal-part obstruction}\\
0001&\text{degree-four jets}\\
0010&\text{boundary Schur classification and cubic slice}\\
0011&X=0,\ Y=0\text{ coefficient ideals}\\
0100&\text{degree-four jets}\\
0101&\text{degree-four jets}\\
0110&X=0,\ Y=0,\ W=0\text{ coefficient ideals}\\
0111&Y=0\text{ coefficient ideal}\\
1000&\text{boundary Schur classification, \(Y=0\) slice, and two points}\\
1001&X=0,\ Y=0,\ W=0\text{ coefficient ideals}\\
1010&\text{degree-four jets}\\
1011&\text{degree-five jets}\\
1100&X=0,\ Y=0,\ W=0\text{ coefficient ideals}\\
1101&X=0,\ Y=0,\ W=0\text{ coefficient ideals}\\
1110&\text{degree-five jets}\\
1111&\text{degree-four jets}.
\end{array}
\]

Thus unrestricted quadratic--cubic single shears are exactly ruled out in
all 16 charts.

The native matrix-first slice checker also accepts an unrestricted quartic
part, adding all 35 homogeneous quartic coefficients.  As an initial
higher-degree screen, chart `0111` on \(Y=0\) and chart `0011` on
\(X=0,Y=0\), both over \(\mathbb F_{32003}\), reached the 300-second runtime
limit.  These timeouts are not mathematical evidence.  They show that the
direct quartic coefficient ideal needs structural elimination before exact
characteristic-zero computation.

The chart-`1000` checker now implements the first such elimination on the
irreducible characteristic branch.  It can impose the five quadratic
coefficients of

\[
\det\operatorname{Hess}\bigl(V_4(a,b,0,b)\bigr)=0,
\]

or cover that zero-Hessian locus by

\[
V_4(a,b,0,b)=\lambda(a+tb)^4,\qquad
V_4(a,b,0,b)=\lambda(ta+b)^4.
\]

On either affine chart five quartic coefficients are eliminated linearly
and only \(\lambda,t\) are introduced.  Over \(\mathbb F_{32003}\), the
five-equation \(X=0\) screen and the \(a\)-normalized pure-power screens on
\(X=0\) and \(Y=0\) each still reached 300 seconds.  These timeouts are
again not evidence.  They locate the remaining complexity in the lower
quartic coefficients and motivate imposing the exact \(A,C,E\) divisor
branches before any further Gröbner run.

### 9.5 Two-step symplectic compositions

To leave the single-shear class, set

\[
 M_1=M+\nabla V(Q),\qquad
 Q_1=Q+\nabla W(M_1).
\]

The four components of \(Q_1\) Poisson-commute.  The search allowed

- \(V=0\) or one signed cubic monomial;
- a quadratic Hessian in \(W\) with at most two signed entries; and
- zero or one signed cubic monomial in \(W\).

All \(337881\) resulting two-step symplectic polarizations were tested at
nine modular points.  None survived the constant-nonzero determinant gate.

The smallest genuinely untested class now lies in degree-at-least-four
single-shear potentials not covered by the all-degree subclasses, denser
two-step compositions, or nonlinear symplectic transformations not
expressible as alternating elementary shears
\(m_0\mapsto m_0+\nabla V(q_0)\).

### 9.6 Caustic Schur reduction and collision conditioning

For a chart \((q_0,m_0)\), put

\[
A=\frac{\partial q_0}{\partial h},\qquad
C=\frac{\partial m_0}{\partial h},\qquad
\Delta=\det A,\qquad
S=CA^{-1}.
\]

The Lagrangian identity makes \(S\) symmetric.  Four charts have a simple
projection caustic through the graph origin.  If
\(N=C\operatorname{adj}(A)\), then \(S=N/\Delta\), and the exact local data
are

\[
\begin{array}{c|c|c}
\text{chart}&\nabla\Delta(0)&N(0)\\ \hline
0010&(-8/3,0,0,0)&-\frac12vv^t\\
0111&(0,5/2,0,0)&\frac12ww^t\\
1000&(-32/3,0,0,0)&-2vv^t\\
1101&(0,10,0,0)&2ww^t,
\end{array}
\]

where

\[
v=(0,1,0,-1)^t,\qquad w=(0,1,0,1)^t.
\]

In every case \(A(0)\) has rank three, so the caustic is smooth there and
the Schur polar numerator has rank one.  The two \(X\)-caustic charts admit
the stronger boundary identities

\[
\begin{aligned}
\left.XS_{0010}\right|_{X=0}
 &=\frac3{16}
 \begin{pmatrix}4r/3\\1\\0\\-1\end{pmatrix}
 \begin{pmatrix}4r/3&1&0&-1\end{pmatrix},\\
\left.XS_{1000}\right|_{X=0}
 &=\frac3{16}
 \begin{pmatrix}0\\1\\2r/3\\-1\end{pmatrix}
 \begin{pmatrix}0&1&2r/3&-1\end{pmatrix},
\end{aligned}
\qquad r=2W+9Y^2.
\]

These formulas reduce the leading Laurent coefficient of

\[
\det(S+\operatorname{Hess}V)=\kappa/\Delta
\]

to a scalar Schur-complement equation.  They do not by themselves obstruct
arbitrary polynomial Hessian data; normal derivatives and Hessian
integrability still have to be imposed.

For chart `0010`, the boundary equation has a sharper triangular form.  Use
adapted potential coordinates

\[
u_0=a,\qquad u_1=b,\qquad u_2=c,\qquad u_3=b+d,
\]

so the caustic image is \(a=d=0\).  Write

\[
f=V|_{a=d=0},\qquad
g=V_a|_{a=d=0},\qquad
h=V_d|_{a=d=0},
\]

and put \(A=V_{aa}\), \(B=V_{ad}\), \(C=V_{dd}\) on that plane.  Define

\[
L=5b f_{cc}-(f_{bb}f_{cc}-f_{bc}^2)
  =-\det\operatorname{Hess}_{b,c}\left(f-\frac56b^3\right).
\]

The restriction of \(\det dB\) to \(X=0\) has \(W\)-degree two.  Its
\(W^2\) coefficient is

\[
-\frac{32}{9}\det
\begin{pmatrix}
f_{bb}-5b&f_{bc}&h_b-b/2\\
f_{bc}&f_{cc}&h_c\\
h_b-b/2&h_c&C+13b/64
\end{pmatrix}.
\]

More importantly, the \(W^2,W^1,W^0\) equations are triangular in
\(C,B,A\), with respective pivots

\[
\frac{32}{9}L,\qquad \frac83L,\qquad \frac12L.
\]

Thus \(L\ne0\) uniquely forces rational boundary values of \(C,B,A\).
There is no pointwise obstruction on this branch; polynomial divisibility
of the forced expressions and compatibility with higher normal jets are
the remaining gates.

The degenerate branch \(L=0\) can be closed for potentials of degree at
most three.  The two-variable zero-Hessian normal form gives

\[
f=\frac56b^3+P(\alpha b+\beta c)
\]

up to affine terms.  Nonzero boundary determinant forces \(P''\ne0\).
The remaining equation contains

\[
4\beta b^3+\alpha g_c-2\alpha-\beta g_b.
\]

Since \(g\) has degree at most two for a cubic potential, constancy forces
\(\beta=0\).  After rescaling, the complete relevant cubic Cauchy data are

\[
\begin{aligned}
f&=\frac56b^3+\frac{\lambda}{2}b^2,\\
g&=(\mu+2)c+g_1b+g_2b^2,\\
h&=\frac14b^2+h_1b+h_2b^2,
\end{aligned}
\qquad \lambda\mu\ne0.
\]

Adjoin all 13 cubic-or-lower terms containing at least two of the normal
variables \(a,d\).  On the smaller graph slice \(Y=D=0\), constancy gives
49 coefficient equations in 19 potential parameters; saturation by
\(\lambda\mu\) adds one variable.  Exact `slimgb` reduction over
\(\mathbb Q\) gives the unit ideal.  Therefore:

\[
\boxed{\text{Chart `0010` has no quadratic--cubic solution on the branch }
L=0.}
\]

It remains to show that every cubic boundary solution lies on this branch.
Parameterize arbitrary degree-at-most-three boundary data by seven
coefficients of \(f\), five each of \(g,h\), and three each of \(A,B,C\).
The \(W^2,W^1,W^0-\kappa\) identities give 52 coefficient equations in
these 26 parameters.  After adjoining \(z\kappa-1\), exact Gröbner reduction
over \(\mathbb Q\) contains

\[
\begin{gathered}
f_{03},f_{12},f_{21},f_{02},f_{11},g_{02},
h_{02},h_{11},h_{01},\\
g_{11}^2,\qquad (f_{30}-5)^2,\qquad
f_{20}(g_{01}-2)^2-2\kappa.
\end{gathered}
\]

At every field-valued solution this is precisely the normalized \(L=0\)
family above, with

\[
\lambda=f_{20},\qquad \mu=g_{01}-2,\qquad
\lambda\mu\ne0.
\]

Combining the boundary classification with the 49-equation slice
certificate proves

\[
\boxed{\text{Chart `0010` has no unrestricted quadratic--cubic
single-shear solution.}}
\]

For potentials of higher degree, the generic branch \(L\ne0\) remains open:
the forced rational data can have degrees larger than those allowed in the
cubic classification.

Chart `1000` has a parallel but distinct fixed-image calculation.  At
\(X=0\), fixing the independent image to \((a,b,0,b)\) forces

\[
D=-\frac a2-\frac{W^2}{3}-3Wb^2-\frac{29}{4}b^4.
\]

Use the same adapted coordinates

\[
u_0=a,\qquad u_1=b,\qquad u_2=c,\qquad u_3=b+d,
\]

but now the caustic image is \(c=d=0\).  Put

\[
f=V|_{c=d=0},\qquad g=V_c|_{c=d=0},\qquad
h=V_d|_{c=d=0}
\]

and denote the normal Hessian entries by \(V_{cc},V_{cd},V_{dd}\).
The fixed-image determinant again has \(W\)-degree two.  With

\[
L=5b f_{aa}-(f_{aa}f_{bb}-f_{ab}^2)
  =-\det\operatorname{Hess}_{a,b}
       \left(f-\frac56b^3\right),
\]

its \(W^2\) coefficient is

\[
-\frac{32}{9}\det
\begin{pmatrix}
f_{aa}&f_{ab}&h_a\\
f_{ab}&f_{bb}-5b&h_b-b/2\\
h_a&h_b-b/2&V_{dd}-5b/64
\end{pmatrix}.
\]

The \(W^2,W^1,W^0\) equations are triangular in
\(V_{dd},V_{cd},V_{cc}\), with pivots

\[
\frac{32}{9}L,\qquad \frac{16}{3}L,\qquad 2L.
\]

For arbitrary cubic boundary data, the three determinant identities give
52 coefficient equations in 26 parameters.  After adjoining
\(z\kappa-1\), exact Gröbner reduction over \(\mathbb Q\) forces every
field-valued solution into

\[
\begin{aligned}
f&=\frac56b^3+\frac{\lambda}{2}b^2,\\
g&=\frac{\mu+1}{2}a+g_1b+g_2b^2,\\
h&=h_1b+h_2b^2,
\end{aligned}
\qquad
\kappa=\frac{\lambda\mu^2}{2},\qquad \lambda\mu\ne0.
\]

In particular every cubic boundary solution has \(L=0\).  Adjoin all
degree-at-most-three terms containing at least two of the normal variables
\(c,d\).  The restriction of the full determinant identity to \(Y=0\)
has 185 coefficients.  Together with saturation by \(\lambda\mu\), their
exact Gröbner basis has 22 elements and leaves a zero-dimensional residue.
Adding the two exact graph evaluations

\[
(X,Y,W,D)=(1,1,0,0),\qquad (1,1,1,0)
\]

makes the quotient ideal the unit ideal over \(\mathbb Q\).  Therefore

\[
\boxed{\text{Chart `1000` has no unrestricted quadratic--cubic
single-shear solution.}}
\]

As in chart `0010`, this does not close higher-degree potentials on the
generic branch \(L\ne0\).

Formal compatibility does not obstruct either generic \(X\)-caustic branch.
Let \(Q\) be a variation of the normal \(2\times2\) Hessian block and put

\[
\tau=2W+9b^2.
\]

The exact boundary linearizations in charts `0010` and `1000` are

\[
D_H(\det J)(Q)=\frac L2Q(v_\epsilon,v_\epsilon),
\qquad
v_{0010}=\left(1,\frac43\tau\right),\qquad
v_{1000}=\left(2,\frac43\tau\right).
\]

At normal prolongation order \(r\), the new pure-normal derivative is a
binary symmetric tensor \(T\) of order \(r+2\).  Its principal symbol is

\[
T\longmapsto \frac L2T(v_\epsilon,\ldots,v_\epsilon)
 \in k[W]_{\le r+2}.
\]

The second component of \(v_\epsilon\) is affine-linear and nonconstant in
\(W\).  In the monomial bases of binary forms and \(W\)-polynomials, this
map is triangular with nonzero diagonal after \(L\) is inverted.  It is
therefore an isomorphism at every order.  Tangential differentiation fixes
the mixed derivatives, and the symbol then uniquely solves the new
pure-normal derivatives.  Consequently the generic branch is formally
recursively solvable rather than formally inconsistent.

For chart `1000`, the complete first prolongation has \(W\)-degree three.
After clearing the already-forced Schur denominators, its pivots for
\(V_{ddd},V_{cdd},V_{ccd},V_{ccc}\) are

\[
819200L^4,\quad1843200L^4,\quad1382400L^4,\quad345600L^4.
\]

The second principal symbol has pivots

\[
1024L,\quad3072L,\quad3456L,\quad1728L,\quad324L.
\]

Thus prolonging the differential ideal cannot by itself close \(L\ne0\).
The remaining gates are global polynomial divisibility of the recursively
forced rational jets and termination of the normal Taylor series.

The complete first divisor-local chain is explicit in chart `1000`.  Put

\[
\begin{aligned}
A&=bf_{aa}-2f_{aa}h_b+2f_{ab}h_a,\\
C&=4b^3f_{aa}-2f_{aa}g_b+2f_{ab}g_a-f_{ab},\\
B&=f_{aa}^3f_{bbb}-5f_{aa}^3
-3f_{aa}^2f_{ab}f_{abb}
+3f_{aa}f_{aab}f_{ab}^2-f_{aaa}f_{ab}^3.
\end{aligned}
\]

Polynomiality of the already-forced normal Hessian is the first divisor
gate.  After clearing the displayed scalar denominators, the numerators of
\(V_{dd},V_{cd},V_{cc}\) restrict to \(L=0\) as

\[
-\frac{16A^2}{f_{aa}},\qquad
-\frac{4AC}{f_{aa}},\qquad
-\frac{20(C^2-2f_{aa}\kappa)}{f_{aa}}.
\]

Therefore every reduced irreducible component met by polynomial normal
Hessian data, with \(f_{aa}\ne0\) generically, must satisfy

\[
\boxed{A=0,\qquad C^2=2f_{aa}\kappa.}
\]

Since \(\kappa\ne0\), both \(f_{aa}\) and \(C\) are units at its generic
point.  This is the unique componentwise branch; it is stronger than the
factor alternatives obtained from the third jets alone.

The \(W^3\) equation in the first prolongation reduces to

\[
256L^3V_{ddd}+\mathcal R=0
\]

for a polynomial differential expression \(\mathcal R\) in the tangential
jets of \(f\) and \(h\).  The other three \(W\)-coefficients recursively
solve \(V_{cdd},V_{ccd},V_{ccc}\).  On \(L=0\), in the chart
\(f_{aa}\ne0\), their leading Laurent coefficients are respectively

\[
\begin{aligned}
V_{ddd}&:\ -\frac{A^3B}{4f_{aa}^3}L^{-3},\\
V_{cdd}&:\ -\frac{A^2CB}{4f_{aa}^3}L^{-3},\\
V_{ccd}&:\ -\frac{AB(3C^2-2f_{aa}\kappa)}
 {12f_{aa}^3}L^{-3},\\
V_{ccc}&:\ -\frac{CB(C^2-2f_{aa}\kappa)}
 {4f_{aa}^3}L^{-3}.
\end{aligned}
\]

Considered without the preceding normal-Hessian gate, these four leading
third-jet residues would allow \(B=0\), or
\(A=0\) together with \(C=0\) or
\(C^2=2f_{aa}\kappa\).  Those are not separate polynomial branches:
normal-Hessian polynomiality has already forced
\(A=0,\ C^2=2f_{aa}\kappa\), and excludes \(C=0\).
Consequently all four \(L^{-3}\) coefficients vanish automatically on the
actual branch.

The factor \(B\) is still useful for describing its characteristic
intersection.  Modulo \(L\),

\[
B=f_{aa}(f_{ab}\partial_a-f_{aa}\partial_b)L.
\]

Hence \(B=0\) says that the Hamiltonian derivation of \(f_a\) preserves the
caustic divisor.  More precisely, let \(p\) be a reduced irreducible
component of \(L=0\), with \(f_{aa}\ne0\) generically on \(p\).  If
\(B=0\) modulo \(p\), then the Hamiltonian field
\(f_{ab}\partial_a-f_{aa}\partial_b\) is tangent to \(p\).  It is nonzero
at the generic point, so \(f_a\) is constant in the one-variable function
field of \(p\).  Therefore

\[
p\mid(f_a-c)
\]

for some constant \(c\).  In particular, if \(L\) itself is irreducible and
\(\deg f=d\ge4\), then \(L\mid(f_a-c)\) forces
\(\deg L\le d-1\).  The generic top degree \(2d-4\) must disappear, so

\[
\det\operatorname{Hess}(f_d)=0.
\]

The binary zero-Hessian theorem then gives

\[
f_d=\lambda(\alpha a+\beta b)^d.
\]

Thus the irreducible characteristic branch reduces every higher-degree
search to a pure-power leading tangential form.

Because normal-Hessian polynomiality already gives \(A=0\), the next two
characteristic \(V_{ddd}\) residues, which contain \(A^2\) and \(A\),
vanish automatically.  The first nontrivial later condition comes from
\(V_{cdd}\).  Put

\[
H_2=f_{aa}^2h_{ab}-f_{aa}f_{aab}h_a
-f_{aa}f_{ab}h_{aa}+f_{aaa}f_{ab}h_a.
\]

The unique branch can now be prolonged componentwise.  Assume

\[
A=0,\qquad C^2=2f_{aa}\kappa
\]

on a reduced irreducible component \(p\) of \(L=0\), with
\(f_{aa}C\ne0\) generically.  The \(V_{cdd}\) pole first gives

\[
h_aH_2=0.
\]

Put \(F=f_{aa}\), \(P=f_{ab}\), and define

\[
S_a=2F\,\partial_aC-Cf_{aaa},\qquad
S_b=2F\,\partial_bC-Cf_{aab},
\]

\[
R=21bF^2-8Fg_{aa}+8f_{aaa}g_a-4f_{aaa}.
\]

After removing the nonzero factors \(C\) and \(F\), the next four Laurent
conditions are \(K_1=K_2=K_3=K_4=0\).  Exact reduction gives the compact
derivative identities

\[
\begin{aligned}
H_2&=-\frac F2\,\partial_aA,\\
K_1&=-F^2\partial_bA-2PH_2,\\
K_2&=C(Fh_{aa}-f_{aaa}h_a)
-6(2g_a-1)H_2+3h_aS_a,\\
K_3&=FS_b-PS_a,\\
K_4&=CR-12(2g_a-1)S_a.
\end{aligned}
\]

These identities are valid after substituting the component equations, but
do not differentiate them away.  They yield a useful multiplicity
obstruction.  If \(h_a\ne0\) generically on \(p\), then \(H_2=0\), hence
\(\partial_aA=0\); next \(K_1=0\) gives \(\partial_bA=0\).  Since \(p\mid A\)
and \(p\) is reduced,

\[
\boxed{p^2\mid A.}
\]

For a quartic potential, \(\deg A\le4\).  Therefore any component of degree
at least three on this \(h_a\ne0\) branch forces the stronger polynomial
identity \(A=0\).  Only line and conic components can support a nonzero
\(A\).

The \(h_a=0\) line branch also reduces exactly.  For a nonvertical line
\(a=\rho b+\sigma\), tangent compatibility gives

\[
h_{ab}=-\rho h_{aa},\qquad
h_{bb}=\frac12+\rho^2h_{aa}.
\]

Writing \(U=F\rho+P\), substitution into the Laurent equations yields

\[
H_2=-FU h_{aa},\qquad
K_1=2FU^2h_{aa},
\]

\[
K_2=Fh_{aa}\bigl(C+6(2g_a-1)U\bigr).
\]

Because \(FC\ne0\), \(K_1=K_2=0\) forces

\[
h_{aa}=h_{ab}=0,\qquad h_{bb}=\frac12.
\]

The factor \(U\) is the derivative of \(f_a\) along the line.  If \(U\ne0\),
tangential differentiation of \(C^2/F=2\kappa\), together with \(K_3=0\),
forces \(S_a=S_b=0\); then \(K_4=0\) forces \(R=0\).  If \(U=0\), the line is
characteristic for the Hamiltonian field of \(f_a\) and joins the
characteristic branch.  Vertical lines have the same rigid \(h\)-jet
conclusion by using \(\partial_a\) as tangent.

Characteristic lines have an all-degree normal form.  Put

\[
u=a-\rho b-\sigma,\qquad v=b.
\]

The condition \(U=0\) says that \(f_a=f_u\) is constant on \(u=0\).  Hence

\[
f=cu+\phi(v)+u^2q(u,v).
\]

Direct substitution into the caustic polynomial gives

\[
\left.L\right|_{u=0}
=f_{uu}(0,v)\bigl(5v-\phi''(v)\bigr).
\]

Since \(f_{uu}\ne0\) generically on the chosen component,

\[
\boxed{\left.f\right|_{u=0}
=\frac56v^3+\alpha v+\beta.}
\]

For a quartic potential, \(q\) has degree at most two.  This replaces the
arbitrary quartic boundary datum on every characteristic line by six
coefficients in \(q\), the line parameters, and three affine constants.

The next Laurent conditions reduce just as sharply.  Write

\[
h_0(v)=h(0,v),\qquad h_1(v)=h_u(0,v),\qquad F(v)=f_{uu}(0,v).
\]

On the characteristic line,

\[
A=F(v-2h_0'),\qquad
H_1=F^3(2h_0''-1),\qquad
H_2=F(Fh_1'-F'h_1).
\]

The normal-Hessian gate \(A=0\), followed by the later condition
\(h_1H_2=0\), gives

\[
h_0=\frac14v^2+\delta,\qquad h_1=\lambda F.
\]

Thus every surviving quartic characteristic line has the displayed rigid
quadratic \(h\)-restriction, and its normal first jet is one scalar multiple
of \(F\).

The remaining normal-Hessian condition now closes the line.  If
\(g_0(v)=g(0,v)\), then

\[
\left.C\right|_{u=0}
=F(v)\bigl(4v^3-2g_0'(v)+\rho\bigr).
\]

For a quartic potential, \(\deg g_0\le3\), so the parenthesized polynomial
has degree exactly three.  The equation \(C^2=2F\kappa\) would imply

\[
F(v)\bigl(4v^3-2g_0'(v)+\rho\bigr)^2=2\kappa,
\]

which is impossible because \(F\ne0\) and \(\kappa\ne0\).  Therefore

\[
\boxed{\text{No quartic characteristic line component survives.}}
\]

Noncharacteristic lines admit an equally explicit surviving normal form.
Again put \(u=a-\rho b-\sigma,\ v=b\), and write

\[
F=f_{aa}|_{u=0},\qquad
U=\partial_v(f_a|_{u=0}),\qquad
f_0=f|_{u=0}.
\]

The caustic and normal-Hessian equations restrict to

\[
Ff_0''-U^2=5vF,\qquad C^2=2F\kappa.
\]

Because \(F\) has degree at most two and \(C\) is polynomial, the second
equation forces \(\deg F\) even.  If \(\deg F=2\), the term \(4v^3F\) makes
\(\deg C=5\), impossible in \(C^2=2F\kappa\).  Hence \(F\) is a nonzero
constant.  The first equation then forces

\[
U=mv+n
\]

with \(\deg U\le1\).  The constant case \(m=0\) leaves the uncancellable
term \(4Fv^3\) in \(C\), so \(m\ne0\).  Integration gives

\[
\begin{aligned}
f_0={}&\frac{m^2}{12F}v^4
+\left(\frac56+\frac{mn}{3F}\right)v^3
+\frac{n^2}{2F}v^2+\alpha v+\beta,\\
f_a|_{u=0}={}&\frac m2v^2+nv+\eta.
\end{aligned}
\]

If \(g_1=g_a|_{u=0}=g_{12}v^2+g_{11}v+g_{10}\), constancy of \(C\) first
forces

\[
g_{12}=-\frac{2F}{m},
\]

and then determines the cubic and quadratic coefficients of
\(g_0=g|_{u=0}\):

\[
[v^3]g_0=\frac{mg_{11}+ng_{12}}{3F},\qquad
[v^2]g_0=\frac{2mg_{10}+2ng_{11}-m}{4F}.
\]

The remaining first-prolongation equations are triangular on this form.
Writing \(h_1=h_a|_{u=0}\) and
\(A_3=f_{aaa}|_{u=0}\), the equations \(A=0\) and \(h_1H_2=0\) force
\(h_1=\lambda\) to be constant and

\[
\begin{aligned}
h|_{u=0}&=\frac14v^2
+\frac{\lambda}{F}\left(\frac m2v^2+nv\right)+\delta,\\
[u^2]h&=\frac{\lambda A_3}{2F}.
\end{aligned}
\]

Indeed, \(A=0\) first gives \(\deg h_1\le1\); solving \(H_2=0\) for
\([u^2]h\) then requires the nonconstant linear polynomial \(U\) to divide
the constant \(F^2h_1'\), hence \(h_1'=0\).

Here the case \(\lambda=0\) is included.  Then \(K_1=0\) is automatic.
The noncharacteristic tangent equations and \(K_3=0\) give
\(S_a=S_b=0\).  Coefficient comparison in \(S_a\) successively forces

\[
\begin{aligned}
A_3&=0,\\
[v][u^2]g&=0,\qquad [u^2]g=-\frac{2F^2}{m^2},\\
[v]g_1&=-\frac{4Fn}{m^2}.
\end{aligned}
\]

After these substitutions, \(K_4=0\) becomes

\[
R=21F^2v+\frac{32F^3}{m^2}=0,
\]

which is impossible because \(F,m\ne0\).  Together with the characteristic
calculation,

\[
\boxed{\text{No quartic line component of }L=0\text{ survives.}}
\]

The characteristic conic branch also has a sharp transverse obstruction.
Suppose first that the conic has nonzero \(a^2\)-coefficient.  After an
affine shear preserving \(b\), write

\[
u=a-\rho b-\sigma,\qquad v=b,\qquad
p=u^2+q(v).
\]

Characteristicity means

\[
f_a-c=p\ell
\]

for an affine-linear polynomial \(\ell\).  On \(p=0\), the equation
\(C^2=2f_{aa}\kappa\) becomes

\[
\ell D^2=2\kappa p_u,\qquad
D=p_u(4v^3-2g_b)+p_b(2g_a-1).
\]

The line \(p_u=2u=0\) cuts the smooth projective conic in two distinct
points.  Since \(p_u/\ell\) must be a square in its function field, divisor
parity forces \(\ell=xu\).  Caustic divisibility then integrates exactly to

\[
\boxed{
f=cu+\frac{x}{4}p^2+\frac56v^3+\alpha v+\beta .
}
\]

The normal-Hessian equation reduces to \(xD^2=4\kappa\), so \(D\) must be
constant in the conic coordinate ring.  But for an arbitrary cubic \(g\),
division by \(p=u^2+q(v)\) leaves

\[
\operatorname{rem}_p(D)=8uv^3+\text{terms of total degree at most three}.
\]

The coefficient \(8\) is independent of every coefficient of \(g\) and
of the conic.  Therefore

\[
\boxed{\text{No transverse quartic characteristic conic survives.}}
\]

For completeness, conics with \(p_{aa}=0\) split into parabolic and
hyperbolic normal forms.  In the parabolic form

\[
p=a+q(b),\qquad \deg q=2,
\]

divisor parity forces \(\ell=z\ne0\).  Caustic divisibility gives

\[
f=ca+\frac z2p^2+\frac56b^3+\text{affine}.
\]

Writing \(C=zD\), the equation \(C^2=2z\kappa\) forces \(D\) to be
constant.  Exact coefficient comparison for a general cubic \(g\) then
makes \(g_{aa}\) constant.  If \(h_a\ne0\) generically on the conic, the
\(K_2,K_3,K_4\) chain would require

\[
R=21bz^2-8zg_{aa}=0,
\]

which is impossible.  If \(h_a=0\), mixed-derivative compatibility gives
\(h=b^2/4+\text{constant}\), and the already-forced normal Hessian entry is

\[
V_{cd}=\frac5{16}b^3.
\]

This also is impossible for a quartic potential, since every second
derivative of \(V\) has degree at most two.

In the hyperbolic form, use coordinates

\[
p=uv-1,\qquad
\partial_a=\partial_u,\quad
\partial_b=\partial_v-\rho\partial_u.
\]

The coordinate ring is \(k[v,v^{-1}]\).  Its unit parity leaves precisely
the two tangent multipliers \(\ell=xu\) and \(\ell=xv\).  For
\(\ell=xu\), caustic coefficient comparison leaves the nonzero residue

\[
\left.L\right|_p=\frac{x^2}{v^4}.
\]

For \(\ell=xv\), caustic divisibility succeeds, but the equation
\(C^2=2f_{aa}\kappa\) requires \(D\) to be constant while every cubic \(g\)
has

\[
\left.D\right|_p=4v^4+\text{lower Laurent powers}.
\]

Consequently

\[
\boxed{\text{No asymptotic quartic characteristic conic survives.}}
\]

The noncharacteristic conic now has only one branch left.  First suppose
\(h_a=0\) generically on a geometrically irreducible conic \(p=0\).  Since
\(\deg h_a\le2\) and \(\deg(h_b-b/2)\le2\), componentwise vanishing gives

\[
h_a=\mu p,\qquad h_b-\frac b2=\nu p
\]

for constants \(\mu,\nu\).  Symmetry of the mixed derivative gives
\(\mu p_b=\nu p_a\).  The two partial derivatives of a geometrically
irreducible conic are not proportional, so \(\mu=\nu=0\).  Hence

\[
h_a=0,\qquad h_b=\frac b2
\]

as polynomial identities.  But the Schur solution then has the universal
entry

\[
V_{cd}=\frac5{16}b^3,
\]

again contradicting the quartic degree bound.  Thus \(h_a=0\) is closed
for every conic, not only the characteristic normal forms above.

It remains to take \(h_a\ne0\) generically.  The multiplicity argument
gives \(p^2\mid A\), while the tangent derivative of
\(C^2/F=2\kappa\) and \(K_3=0\) form the system

\[
p_bS_a-p_aS_b=0,\qquad FS_b-PS_a=0.
\]

Its determinant

\[
U=Fp_b-Pp_a
\]

is the derivative of \(f_a\) along the conic.  On the
noncharacteristic branch \(U\ne0\), hence

\[
S_a=S_b=0.
\]

Then \(K_4=0\) forces \(R=0\), while \(K_2=H_2=0\) give

\[
\partial_a\!\left(\frac{h_a}{F}\right)=
\partial_b\!\left(\frac{h_a}{F}\right)=0
\qquad\text{modulo }p.
\]

In characteristic zero \(h_a/F\) is constant on the conic.  Since both
polynomials have degree at most two, there are constants \(\lambda,\mu\)
such that

\[
h_a=\lambda f_{aa}+\mu p.
\]

Put \(r=h-\lambda f_a\).  The equation \(A=0\) on \(p\) also says that
\(r_b-b/2\) vanishes on \(p\).  Both \(r_a\) and \(r_b-b/2\) have degree at
most two, so

\[
r_a=\mu p,\qquad r_b-\frac b2=\nu p
\]

for constants \(\mu,\nu\).  Equality of mixed derivatives gives
\(\mu p_b=\nu p_a\).  Since the two partial derivatives of an irreducible
conic are not proportional, \(\mu=\nu=0\).  Therefore

\[
h=\lambda f_a+\frac14b^2+\text{constant}
\]

as a polynomial identity.  Substitution into the Schur solution gives

\[
V_{cd}=\frac5{16}b^3+\lambda\left(g_a-\frac12\right).
\]

For a quartic potential, \(g\) has degree at most three and hence \(g_a\)
has degree at most two.  The cubic term \(5b^3/16\) cannot cancel.
Consequently

\[
\boxed{\text{No quartic noncharacteristic conic survives.}}
\]

Together with the line and characteristic calculations,

\[
\boxed{\text{No reduced line or conic component of the chart-1000
caustic can support quartic polynomial normal data.}}
\]

In fact, the same proof closes every noncharacteristic component.  If
\(\deg p\ge3\) and \(h_a=0\) modulo \(p\), the degree bounds make
\(h_a=0,\ h_b=b/2\) polynomial identities, giving the same impossible
\(V_{cd}=5b^3/16\).  If \(h_a\ne0\), then \(p^2\mid A\) and
\(\deg A\le4\) force \(A=0\) identically.  The tangent/\(K_3\) system
again gives \(S_a=S_b=0\), and \(K_2=H_2=0\) makes \(h_a/F\) constant on
\(p\).  Since \(\deg(h_a-\lambda F)\le2<\deg p\), this is the polynomial
identity \(h_a=\lambda F\); the identity \(A=0\) then gives
\(h_b=b/2+\lambda P\).  Thus the same uncancellable cubic occurs:

\[
\boxed{\text{No reduced noncharacteristic component of the chart-1000
caustic can support quartic polynomial normal data.}}
\]

A characteristic component satisfies \(p\mid(f_a-c)\), so
\(\deg p\le3\).  Lines and conics are closed above.  The only remaining
reduced caustic component is therefore a characteristic cubic, followed
by the unit-\(L\) branch if that cubic is excluded.

The characteristic cubic is also impossible.  Write
\(f_a-c=\alpha p\) and rescale \(p\) so \(\alpha=1\).  If \(f_4\) is the
quartic homogeneous part, comparison of the degree-four part of
\(p\mid L\) gives

\[
(f_4)_a\mid\det\operatorname{Hess}(f_4).
\]

The exact binary-quartic divisibility equations have, up to a shear
preserving \(b\), only three possible leading cubic forms:

\[
p_3=a^3,\qquad p_3=ab^2,\qquad p_3=b^3.
\]

Writing the general lower terms, integrating \(f_a=p\), and imposing the
full divisibility \(p\mid L\) gives respectively

\[
\begin{aligned}
p&=a^3+A a^2+D a+F,\\
p&=(a+C)(b+B/2)^2,\\
p&=(b+C/3)^3.
\end{aligned}
\]

The first depends only on \(a\), and the other two are visibly reducible.
Over an algebraically closed characteristic-zero field none is an
irreducible cubic.  Hence

\[
\boxed{\text{No reduced characteristic cubic component survives.}}
\]

Consequently, if \(L\) is not the zero polynomial, polynomial normal data
force

\[
\boxed{L\in k^\times.}
\]

There remains a distinct degenerate identity branch \(L\equiv0\), not
represented by an irreducible divisor.  On its \(f_{aa}\ne0\) part, the
binary zero-Hessian theorem writes

\[
f=\frac56b^3+P(\alpha a+\beta b)+\text{affine},
\qquad \alpha P''\ne0.
\]

The identities \(A=0,\ C^2=2f_{aa}\kappa\) reduce the second equation to

\[
P''(\alpha a+\beta b)
\bigl(4\alpha b^3-2\alpha g_b+2\beta g_a-\beta\bigr)^2
=2\kappa.
\]

Both factors on the left must be polynomial units.  But the second retains
the term \(4\alpha b^3\), while \(g_a,g_b\) have degree at most two.  Thus
this part is impossible.  The edge \(f_{aa}=0\) reduces directly to

\[
\begin{gathered}
f_a\in k,\qquad f_{bb}=5b-\delta,\\
h_a=0,\qquad g_a=\gamma,\\
\kappa=-\frac{\delta}{2}(2\gamma-1)^2,\qquad
\delta(2\gamma-1)\ne0.
\end{gathered}
\]

This affine-normal edge and the unit-\(L\) branch are the two surviving
global systems.

The unit branch closes at the normal-Hessian degree gate.  Put
\(\phi=f-5b^3/6\).  The binary constant-Hessian equations through degree
four put \(\phi\), after a shear preserving \(b\), in one of the two
triangular forms

\[
\phi=q ab+P(a),\qquad \phi=q ab+P(b),\qquad q\ne0.
\]

In the first form, write
\(P=p_4a^4+p_3a^3+p_2a^2+\text{affine}\).  For an arbitrary cubic \(g\),
the forced \(V_{cc}\) has

\[
[a^2b^6]V_{cc}=-\frac{48p_4}{q^2},\qquad
[ab^6]V_{cc}=-\frac{24p_3}{q^2},\qquad
[b^6]V_{cc}=-\frac{8p_2}{q^2}.
\]

Since \(\deg V_{cc}\le2\), \(P\) is affine.  This moves the data to the
second orientation.  Write

\[
P''=r_2b^2+r_1b+r_0.
\]

If \(r_2\ne0\), the \(b^6\) coefficient first forces the \(ab^2\)
coefficient of \(g\) to vanish, after which the \(b^5\) coefficient is the
nonzero scalar \(99/80\).  Hence \(r_2=0\).  If \(r_1\ne0\), the
\(a^4b,a^3b,a^2b,a b^3,b^5\) chain successively ends in
\(99q^2=0\); the \(r_1=0\) chain does the same.  Therefore

\[
\boxed{\text{The chart-1000 unit-\(L\) branch is impossible for quartic
potentials.}}
\]

It remains only to test the affine-normal zero-\(L\) edge above.  Adjoin
all quartic terms with at least two normal variables: arbitrary quadratic
\(V_{cc},V_{cd},V_{dd}\), arbitrary affine pure-normal third derivatives,
and arbitrary constant pure-normal fourth derivatives.  This is the full
43-parameter quartic completion of the reduced boundary data.  On the
single graph hyperplane \(Y=0\), native Singular expansion gives 1039
coefficient columns.  Saturating by the common constant determinant and
running `slimgb` over \(\mathbb Q\) gives

\[
\operatorname{reduce}(1,G)=1,\qquad |G|=44.
\]

Thus the edge is exactly inconsistent:

\[
\boxed{\text{The affine-normal zero-\(L\) edge has no quartic completion.}}
\]

Combining the caustic-component theorem, the characteristic-cubic
classification, the zero-\(L\) reduction, the unit-\(L\) degree
obstruction, and the final rational unit ideal gives the complete chart
result

\[
\boxed{\text{Chart `1000` admits no quartic polynomial potential with
constant nonzero determinant.}}
\]

No collision equations are used, so this is stronger than failure to
retain the known \(PC(2)\) collision.

The other six formerly unresolved quadratic--cubic charts can be closed
without classifying their boundary Cauchy data.  Let

\[
V=V_2+V_3
\]

have all ten quadratic and twenty homogeneous-cubic coefficients free.  Form

\[
J_\epsilon=d\left(m_\epsilon+\nabla V(q_\epsilon)\right).
\]

If \(\det J_\epsilon\) is a nonzero constant, then its restriction to every
coordinate hyperplane is the same constant.  Substitute a hyperplane into
the matrix \(J_\epsilon\) before taking its determinant, extract all
nonconstant coefficients in the remaining three source variables, and
adjoin \(z c_0-1\), where \(c_0=\det J_\epsilon(0)\).  Direct `slimgb`
reduction over \(\mathbb Q\) gives:

\[
\begin{array}{c|c|c}
\text{chart}&\text{hyperplanes}&
 \text{coefficient columns, including the constant}\\ \hline
0011&X=0,\ Y=0&100,\ 157\\
0110&X=0,\ Y=0,\ W=0&100,\ 263,\ 1074\\
0111&Y=0&192\\
1001&X=0,\ Y=0,\ W=0&190,\ 438,\ 931\\
1100&X=0,\ Y=0,\ W=0&190,\ 622,\ 1565\\
1101&X=0,\ Y=0,\ W=0&206,\ 536,\ 1602.
\end{array}
\]

Every saturated rational ideal in the table is the unit ideal.  These
calculations impose no collision equations and therefore prove the stronger
statement

\[
\boxed{\text{No coordinate chart admits an unrestricted
quadratic--cubic single-shear solution.}}
\]

The matrix-first substitution is computationally essential: expanding the
four-variable determinant before taking a hyperplane section produces severe
expression swell but no stronger equations.

There is also an all-degree shortcut excluding a tempting lower-dimensional
escape.  Suppose a three-dimensional Keller map is jointly affine-linear in
two source variables:

\[
F(s,u,v)=b(s)+u\,a_1(s)+v\,a_2(s).
\]

The coefficients of \(u\) and \(v\) in its Jacobian determinant give

\[
\det(a_1',a_1,a_2)=\det(a_2',a_1,a_2)=0.
\]

Hence the two-plane spanned by \(a_1,a_2\) is constant.  Write
\(a_1\wedge a_2=\rho(s)\omega_0\).  The constant Jacobian equation is

\[
\rho(s)\,\omega_0(b'(s))\in k^\times,
\]

so both polynomial factors are units.  After a linear target change, the
third output is affine-linear in \(s\), while the first two outputs are
affine-linear in \(u,v\) with a \(2\times2\) coefficient matrix of constant
nonzero determinant.  Solving first for \(s\) and then for \(u,v\) gives a
polynomial inverse.  Thus every such map is an automorphism; a rank-two
affine-linear Schur descent cannot produce an `HC_4` counterexample.

The certified collision supplies a separate all-degree incidence
obstruction.  In charts `1100` and `1101`, the symmetric points \(P_+\) and
\(P_-\) have identical \(q_0\)-values but distinct \(m_0\)-values.  Therefore

\[
m_0(P_+)+\nabla V(q_0(P_+))
\ne
m_0(P_-)+\nabla V(q_0(P_-))
\]

for every polynomial \(V\).  No single shear in either chart can retain that
pair, and in particular it cannot retain the complete three-point fiber.
This does not exclude retaining one of the pairs involving the third point
or creating a new collision.

The quadratic--cubic jet checker now accepts `--collision 01`, `02`, `12`,
or `all`.  It eliminates the resulting linear coefficient constraints before
the determinant ideal.  A selected pair leaves 26 shear parameters and the
complete three-point condition leaves 22 in each consistent chart.
Collision-conditioned degree-two and degree-three modular reductions in
chart `0010` still reached the 120-second runtime bound; those timeouts are
not mathematical evidence.

### 9.7 The caustic recursion as a graded meromorphic jet module

The uniform symbol above admits a useful logarithmic packaging, but one
distinction is essential.  The order-raising derivatives are taken in the
two image-normal variables \(c,d\), whereas \(L=0\) is a divisor in the
tangential Cauchy-data plane \(a,b\).  Thus the verified recursion is a
graded meromorphic jet module supported along \(L=0\); it is not yet a
finite-rank logarithmic connection in the usual sense.

Fix either \(X\)-caustic chart, fix polynomial tangential Cauchy data, and
let \(R\) be the resulting tangential coefficient ring, localized away from
the denominators declared nonzero on the branch.  At normal prolongation
order \(r\), put \(n=r+2\) and write

\[
 t_n=(t_0,\ldots,t_n),\qquad
 t_i=T(\underbrace{e_c,\ldots,e_c}_{n-i},
       \underbrace{e_d,\ldots,e_d}_{i})
\]

for the new pure-normal symmetric \(n\)-tensor.  If

\[
 v_\epsilon(W)=\left(\alpha_\epsilon,
       \frac43(2W+9b^2)\right),\qquad
 \alpha_{0010}=1,\quad\alpha_{1000}=2,
\]

then the coefficient of the new tensor in the prolonged determinant
equations is

\[
 \frac L2\,T(v_\epsilon(W),\ldots,v_\epsilon(W)).
                                                        \tag{9.7.1}
\]

Let \(P_n=R[W]_{\le n}\), and identify it with \(R^{n+1}\) by coefficient
extraction.  Equation (9.7.1) defines

\[
 \sigma_n=L M_n:\operatorname{Sym}^n(R^2)^*\longrightarrow P_n,
\]

where, in the displayed tensor basis and the basis \(1,W,\ldots,W^n\),

\[
 (M_n)_{ji}=
 \begin{cases}
 \displaystyle
 \frac12\binom ni\alpha_\epsilon^{\,n-i}
 \left(\frac43\right)^i
 \binom ij2^j(9b^2)^{i-j},&j\le i,\\[6pt]
 0,&j>i.
 \end{cases}                                           \tag{9.7.2}
\]

The matrix is triangular and

\[
 \det M_n=
 2^{-(n+1)}
 \left(\prod_{i=0}^n\binom ni\right)
 \alpha_\epsilon^{\,n(n+1)/2}
 \left(\frac83\right)^{n(n+1)/2}\in k^\times.          \tag{9.7.3}
\]

Consequently \(M_n\) is unimodular over \(R\) at every order.  Before
substituting previously solved rational jets, the full prolonged equation
therefore has the exact form

\[
 \boxed{\quad L M_n t_n=\Psi_n(t_0,\ldots,t_{n-1})\quad}                 \tag{9.7.4}
\]

for a polynomial differential expression \(\Psi_n\) in lower normal jets
and tangential derivatives.  Substitution of earlier solutions creates the
higher visible powers of \(L^{-1}\); it does not change the single factor
\(L\) in the new principal symbol.

This gives a precise two-term normal-symbol complex

\[
 \mathcal C_n^\bullet=
 \left[
 \operatorname{Sym}^n(R^2)^*
 \xrightarrow{\ L M_n\ }P_n
 \right].
\]

If \(R\) is a domain and \(L\ne0\), then

\[
 H^0(\mathcal C_n^\bullet)=0,\qquad
 H^1(\mathcal C_n^\bullet)\simeq (R/(L))^{n+1}.         \tag{9.7.5}
\]

After inverting \(L\), the complex is acyclic and (9.7.4) uniquely solves
every new pure-normal jet.  This is the symbol-level reason that ordinary
formal prolongation cannot obstruct the generic \(L\ne0\) branch.

Now let \(p\) be a reduced irreducible component of \(L=0\), and work at
its generic point.  Write \(L=u p^e\), with \(u\) a unit.  The obstruction
to extending the \(n\)-th solved jet across \(p\) is not presently an
endomorphism but the cokernel class

\[
 \rho_{p,n}=
 \left[u^{-1}M_n^{-1}\Psi_n\right]
 \in (R_{(p)}/p^e)^{n+1}.                             \tag{9.7.6}
\]

For a simple component, pole-free extension is equivalent to
\(\rho_{p,n}=0\).  Once that class vanishes, successive coefficients in the
\(p\)-adic filtration give higher residue conditions.  This recovers the
componentwise calculations already found in chart `1000`:

- at the normal-Hessian stage the three residue coordinates are, up to
  units,
  \[
  A^2,\qquad AC,\qquad C^2-2f_{aa}\kappa,
  \]
  hence \(A=0\) and \(C^2=2f_{aa}\kappa\) at the generic point;
- the four leading third-jet residues are the displayed
  \(A^3B,A^2CB,AB(3C^2-2f_{aa}\kappa)\), and
  \(CB(C^2-2f_{aa}\kappa)\), so they vanish on the actual
  normal-Hessian branch;
- the later conclusion \(p^2\mid A\) is a depth-two cancellation in this
  \(p\)-adic residue filtration, not a new pointwise branch.

This also identifies the exact Spencer statement currently proved.  The
normal two-term symbol complex is acyclic off \(L=0\), while its only
symbol cokernel is the \(L\)-torsion module (9.7.5).  The full Spencer
complex, including tangential compatibility and all nonlinear syzygies,
has not been computed.  The explicit first-prolongation identities verify
the first nontrivial part of that larger complex but do not prove its
finite generation.

#### Why a residue operator and \(b\)-function are not yet available

A logarithmic connection would require a finite-rank lattice
\(\Lambda\) with an operator

\[
 \theta_p=p\partial_p,\qquad \theta_p\Lambda\subseteq\Lambda,
\]

whose reduction modulo \(p\) is a residue endomorphism.  The present
recursion does not supply this object:

1. raising normal degree is not the derivation \(\partial_p\);
2. the rank of the new-jet space is \(n+1\), so it grows with \(n\);
3. \(\Psi_n\) is nonlinear in lower jets.

Accordingly there is currently no intrinsic residue spectrum, no notion of
integral eigenvalue resonance, and no justified Bernstein--Sato polynomial
for the nonlinear infinite jet recursion.  Calling (9.7.6) a residue
*class* rather than a residue *operator* keeps this distinction explicit.

A genuine \(V\)-filtration program can nevertheless be stated precisely.
For each surviving formal branch, linearize the full differential ideal,
form its Rees--Spencer module for the \(p\)-adic filtration, and prove:

1. coherence (preferably holonomicity) and finite generation under normal
   prolongation;
2. regular singularity for \(\theta_p\), producing a finite-rank residue
   operator and a Bernstein polynomial \(b_p(s)\);
3. stabilization of the successive residue/colon ideals, so only finitely
   many cancellation types occur;
4. a separate normal-Euler or nilpotence criterion forcing the uniquely
   solved normal tail either to terminate or to remain infinite.

The fourth item is logically independent of extension across \(p=0\): a
\(b\)-function can control poles without forcing a power series in \(c,d\)
to be a polynomial.

The resulting target statement is therefore recorded as a conjectural
finite-generation theorem, not as a consequence of the present
calculation:

> **Caustic resonance/termination conjecture.** For each coordinate chart
> and each reduced caustic component, the \(p\)-adic Rees--Spencer module of
> the normal recursion is differentially finitely generated.  Its
> pole-free locus has finitely many residue-cancellation strata, and on
> every such stratum the unique normal solution either terminates by a
> finite-order algebraic condition or has infinitely many nonzero normal
> jets.

If proved uniformly in the degree of the boundary data, this conjecture
would give the desired conclusion: outside finitely many resonant
cancellation types, formal solutions acquire a caustic pole or infinite
normal degree.  Equations \(A=0\),
\(C^2=2f_{aa}\kappa\), and \(p^2\mid A\) are the first three pieces of its
candidate \(V\)-filtration.

## 10. Direct one-variable Schur ascent and descent

The literal backward prescription has an immediate integrability gate.  For

\[
 \Xi(x,t)=\Phi(x)+tA(x)+\frac12t^2B(x),
\]

the Schur complement of the \(t,t\) entry is

\[
\begin{aligned}
\mathcal S={}&
\operatorname{Hess}\Phi+t\operatorname{Hess}A
+\frac12t^2\operatorname{Hess}B\\
&-\frac{1}{B}
(\nabla A+t\nabla B)(\nabla A+t\nabla B)^t .
\end{aligned}
\]

It is symmetric.  It therefore cannot literally equal the Jacobian of the
displayed `PC(2)` map: for \(G=(R,T,D,S)\), the skew entry

\[
 \frac{\partial T}{\partial x}-\frac{\partial R}{\partial q}
\]

equals \(16/3\) at \((x,q,p,z)=(0,0,0,1)\).  A generating-family
construction, rather than matrix equality, is needed to transfer the
symplectic graph.

### 10.1 The two polynomial-coordinate generating families

Fix a coordinate graph chart \((q_\epsilon,m_\epsilon)\).  Omit one of the
four \(q_\epsilon\)-coordinates and adjoin a linear auxiliary function

\[
 u=a_0X+a_1Y+a_2W+a_3D.
\]

Requiring the resulting four functions of \((X,Y,W,D)\) to have constant
nonzero Jacobian is a linear coefficient problem in the \(a_i\).  The exact
all-chart calculation leaves only

\[
 (\epsilon,j)=(0010,3),\qquad(0011,3),
\]

where \(j=3\) is the omitted fourth coordinate.  In both cases the
determinant is

\[
 \frac{Xa_1-3a_2}{3}.
\]

Thus \(a_1=0\), \(a_2\ne0\), and additions of \(X,D\) are base gauge.  After
normalization the essential auxiliary coordinate is \(u=W\), and

\[
 w=(x,q,d,u)=(X,Q,D,W)
\]

is a polynomial coordinate system.

The linear restriction on \(u\) can be removed completely when the goal is
to transfer a pair from the certified three-point fiber.  For an omitted
coordinate \(j\), define the Jacobian derivation

\[
 V_{\epsilon,j}(u)
 =\det\frac{\partial(q_{\epsilon,\widehat j},u)}
 {\partial(X,Y,W,D)}.
\]

A polynomial auxiliary coordinate requires
\(V_{\epsilon,j}(u)\in\mathbb Q^\times\).  Equality of the four
complementary graph coordinates at a certified pair first forces

\[
 \epsilon\in\{0000,0001,0010,0011\};
\]

only the symmetric pair \(P_+,P_-\) survives.  In charts `0000` and `0001`,
for every \(j\), the four coefficients
\(V_{\epsilon,j}(X),V_{\epsilon,j}(Y),
V_{\epsilon,j}(W),V_{\epsilon,j}(D)\) have a nonunit common factor or all
vanish.  Their ideal therefore cannot contain a nonzero constant.

For charts `0010` and `0011`, omissions \(j=0,2\) are also impossible.
For `0010`, \(j=0\) has a fixed point at the origin and \(j=2\) has a
nonunit coefficient gcd.  For `0011`, \(j=2\) again has a nonunit gcd, while
\(j=0\) has the common zero

\[
 (X,Y,W,D)=
 \left(\frac{\sqrt{34}-4}{6},\,1,\,
 \frac{84+3\sqrt{34}}{25},\,0\right).
\]

It remains to consider \(j=1,3\).  For \(j=1\), the retained triples are
\((X,D,T)\) and \((X,D,S)\).  Over \(K=\mathbb Q(X,D)\), the last retained
coordinate has the form

\[
 f(Y,W)=a(Y)W+b(Y),
\]

with

\[
 a_T=3X(1+XY)^2,\qquad
 a_S=\frac12(1+XY)^3.
\]

Neither admits a polynomial Jacobian mate.  Indeed, if

\[
 V=-a\partial_Y+(a'W+b')\partial_W,\qquad V(u)=c\in K^\times,
\]

and \(c_n(Y)W^n\) is the top \(W\)-term of \(u\), then

\[
 -ac_n'+na'c_n=0,
\]

so \(c_n\) is a scalar multiple of \(a^n\).  Subtracting the corresponding
multiple of \(f^n\) lowers the \(W\)-degree without changing \(V(u)\).
Iteration reduces to \(u=c_0(Y)\), where

\[
 -a(Y)c_0'(Y)=c.
\]

This is impossible because \(a\) is nonconstant.

For \(j=3\), both charts retain \((X,Q,D)\).  In the polynomial coordinates
\((X,Q,D,W)\),

\[
 V_{\epsilon,3}=-\partial_W.
\]

Consequently every polynomial slice is

\[
 u=-cW+f(X,Q,D),\qquad c\in\mathbb Q^\times.
\]

The function \(f\) is retained-coordinate gauge.  Thus \(W\) is the unique
essential auxiliary coordinate in every degree, and the two generating
families below exhaust the polynomial coordinate-auxiliary route that
transfers the certified collision.

Let \(t=O_\epsilon(w)\) be the omitted graph equation and let
\(A_\epsilon(w)\) be its complementary coordinate:

\[
\begin{array}{c|c|c|c}
\epsilon&O_\epsilon&A_\epsilon&s=A_\epsilon(P_\pm)\\ \hline
0010&T&S&-1/8\\
0011&S&-T&0.
\end{array}
\]

If \(P_\epsilon(w)\) is the exact primitive of the restricted canonical
one-form, the canonical one-auxiliary generating family is

\[
 F_\epsilon(w,t)
 =P_\epsilon(w)+A_\epsilon(w)(t-O_\epsilon(w)).
\]

On \(t=O_\epsilon(w)\), its gradient is the complementary graph coordinate
together with a zero auxiliary component.  In both charts the certified
collision lifts to

\[
 w_\pm=\left(\pm1,\mp\frac23,0,\frac{13}{2}\right)
\]

and gives equal gradients of \(F_\epsilon\).

Every quadratic modification preserving the same graph first jet is

\[
 F_{\epsilon,K}
 =F_\epsilon+K(w)(t-O_\epsilon(w))^2.
\]

Partial Legendre transform at the collision value \(s\) gives

\[
 \psi_{\epsilon,K}
 =P_\epsilon-sO_\epsilon
  -\frac{(A_\epsilon-s)^2}{4K}.                    \tag{10.1}
\]

In each chart \(A_\epsilon-s\) is irreducible over
\(\mathbb Q[x,q,d,u]\).  Hence (10.1) is polynomial only if

\[
 K=c(A_\epsilon-s)^e,\qquad
 c\in\mathbb Q^\times,\quad e\in\{0,1,2\}.
\]

Put \(r=1/(4c)\) and \(p=2-e\).  It remains to test

\[
 \psi_{\epsilon,p,r}
 =P_\epsilon-sO_\epsilon-r(A_\epsilon-s)^p,
 \qquad p=0,1,2.
\]

The following exact evaluations close every case.  Points in the table are
in the \(w=(x,q,d,u)\) coordinates.

\[
\begin{array}{c|c|c|c}
\epsilon&p&\text{Hessian determinants}&\text{consequence}\\ \hline
0010&0&\det H(0)=0&\text{not nonzero constant}\\
0010&1&\det H(0)=0&\text{not nonzero constant}\\
0010&2&\det H(0)=-2r^2,
 \det H(0,0,0,1)=-10r^2&r=0\text{ forced}\\ \hline
0011&0&\det H(0)=1,
 \det H(1,1,0,0)=-120632/3&\text{not constant}\\
0011&1&
\begin{array}{l}
\det H(0)=1,\\
\det H(1,0,0,0)=144r^2+120r+1,\\
\det H(1,1,0,0)=
(936r^2-19319r-120632)/3
\end{array}
&\text{difference polynomials coprime}\\
0011&2&\det H(0)=1,
\det H(1,0,0,0)=(32r+3)/3&r=0\text{ forced}.
\end{array}
\]

Since \(r\ne0\), no polynomial quadratic-pivot Schur descent in either
linear-auxiliary `PC(2)` generating family has constant nonzero Hessian
determinant.

### 10.2 No second linear descent of the Meng--Yang potential

There is also a clean obstruction to descending the Meng--Yang polynomial
itself along a constant linear direction.  Write

\[
 \Psi_{\rm MY}=A^2+11A+2B
\]

in variables \((x_1,x_2,y_1,y_2,y_3)\), with \(A,B\) linear in the three
\(y\)-variables as in their paper.  The coefficient ideal of the third
directional derivative \(D_v^3\Psi_{\rm MY}\) contains \(v_1^3,v_2^3\), and
the derivative vanishes identically after \(v_1=v_2=0\).  Thus the only
directions in which \(\Psi_{\rm MY}\) is at most quadratic are the three
dual-variable directions.

For such a direction \((a,b,c)\), put

\[
 L=D_vA,\qquad K=D_vB.
\]

The quadratic pivot is \(2L^2\).  A partial Legendre transform at a constant
dual value \(\sigma\) can be polynomial only if

\[
 L\mid 2K-\sigma .
\]

For every nonzero \((a,b,c)\), \(L\) and \(K\) have the same total degree, so
the quotient would be a scalar \(\ell\).  The complete coefficient ideal of

\[
 2K-\sigma-\ell L
\]

has Gröbner basis

\[
 \{\sigma,a,b,c\}.
\]

Only the zero direction remains.  Therefore the five-variable
Meng--Yang potential has no further polynomial partial Legendre transform
along any constant linear direction.

These results are exact obstructions, not evidence against every conceivable
Schur descent.  The remaining variants require a non-coordinate embedding
of the symplectic graph or a higher-degree generating family whose critical
equation still has a polynomial solution.  They no longer imitate the
Meng--Yang operation in a controlled finite ansatz.

## Reproduction

Run:

```bash
.venv/bin/python scripts/search_hc4_graph_polarizations.py
.venv/bin/python scripts/search_hc4_lagrangian_shears.py --bound 3
.venv/bin/python scripts/verify_hc4_linear_polarization_obstruction.py
.venv/bin/python scripts/verify_hc4_all_linear_projection_obstruction.py
.venv/bin/python scripts/verify_hc4_nonlinear_shear_obstructions.py
.venv/bin/python scripts/search_hc4_sparse_cubic_shears.py
.venv/bin/python scripts/search_hc4_sparse_quadratic_cubic_shears.py
.venv/bin/python scripts/search_hc4_dense_random_cubic_shears.py
.venv/bin/python scripts/search_hc4_two_step_symplectic_shears.py
.venv/bin/python scripts/verify_hc4_caustic_schur_reduction.py
.venv/bin/python scripts/verify_hc4_direct_schur_descent.py
.venv/bin/python scripts/verify_hc4_0010_boundary_schur_chain.py
.venv/bin/python scripts/verify_hc4_0010_cubic_boundary_classification.py
.venv/bin/python scripts/verify_hc4_0010_degenerate_cubic_branch.py
.venv/bin/python scripts/verify_hc4_1000_boundary_schur_chain.py
.venv/bin/python scripts/verify_hc4_1000_cubic_boundary_classification.py
.venv/bin/python scripts/verify_hc4_1000_degenerate_cubic_branch.py
.venv/bin/python scripts/verify_hc4_x_caustic_formal_compatibility.py
.venv/bin/python scripts/verify_hc4_logarithmic_normal_symbol.py
.venv/bin/python scripts/verify_hc4_1000_divisor_local_chain.py
.venv/bin/python scripts/verify_hc4_1000_characteristic_line_normal_form.py
.venv/bin/python scripts/verify_hc4_1000_noncharacteristic_line_normal_form.py
.venv/bin/python scripts/verify_hc4_1000_characteristic_conic_obstruction.py
.venv/bin/python scripts/verify_hc4_1000_asymptotic_characteristic_conics.py
.venv/bin/python scripts/verify_hc4_1000_noncharacteristic_conic_reduction.py
.venv/bin/python scripts/verify_hc4_1000_characteristic_cubic_obstruction.py
.venv/bin/python scripts/verify_hc4_1000_zero_caustic_identity.py
.venv/bin/python scripts/verify_hc4_1000_unit_caustic_obstruction.py
.venv/bin/python scripts/search_hc4_1000_affine_normal_edge.py \
  --slices Y --characteristic 0
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 0011 --slices X Y
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 0110 --slices X Y W
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 0111 --slices Y
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 1001 --slices X Y W
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 1100 --slices X Y W
.venv/bin/python scripts/verify_hc4_multislice_cubic_obstruction.py \
  --chart 1101 --slices X Y W
.venv/bin/python scripts/verify_hc4_quadratic_cubic_jet_obstruction.py \
  --chart 1111 --characteristic 0 --max-degree 4
.venv/bin/python scripts/verify_hc4_quadratic_cubic_jet_obstruction.py \
  --chart 1110 --characteristic 0 --algorithm slimgb --max-degree 5
```
