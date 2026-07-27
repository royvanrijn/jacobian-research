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
potentials in ten of the 16 charts.  General nonlinear polynomial
polarizations remain open.

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
0100&\text{degree-four jets}\\
0101&\text{degree-four jets}\\
1000&\text{boundary Schur classification, \(Y=0\) slice, and two points}\\
1010&\text{degree-four jets}\\
1011&\text{degree-five jets}\\
1110&\text{degree-five jets}\\
1111&\text{degree-four jets}.
\end{array}
\]

Thus unrestricted quadratic--cubic single shears are exactly ruled out in
ten of the 16 charts.  The remaining charts

\[
0011,\ 0110,\ 0111,\ 1001,\ 1100,\ 1101
\]

remain algebraically unresolved; several degree-four and degree-five
reductions hit the declared runtime bound and are not counted as evidence.

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

The smallest genuinely untested class now lies in the six unresolved
quadratic--cubic charts, denser two-step compositions, or nonlinear
symplectic transformations not expressible as alternating elementary shears
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
.venv/bin/python scripts/verify_hc4_0010_boundary_schur_chain.py
.venv/bin/python scripts/verify_hc4_0010_cubic_boundary_classification.py
.venv/bin/python scripts/verify_hc4_0010_degenerate_cubic_branch.py
.venv/bin/python scripts/verify_hc4_1000_boundary_schur_chain.py
.venv/bin/python scripts/verify_hc4_1000_cubic_boundary_classification.py
.venv/bin/python scripts/verify_hc4_1000_degenerate_cubic_branch.py
.venv/bin/python scripts/verify_hc4_quadratic_cubic_jet_obstruction.py \
  --chart 1111 --characteristic 0 --max-degree 4
.venv/bin/python scripts/verify_hc4_quadratic_cubic_jet_obstruction.py \
  --chart 1110 --characteristic 0 --algorithm slimgb --max-degree 5
```
