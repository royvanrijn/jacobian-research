# Double Schur descent for the quadratic-gauge families

## Status

This note tests the most direct two-step Meng--Yang descent on the full
root-engineered quadratic-gauge family.  No four-variable Hessian
counterexample is found.

The first, six-to-five-variable step is the general
[Meng--Yang Schur-descent bridge](MENG_YANG_SCHUR_DESCENT_BRIDGE.md);
its determinant identity is imported here.  The calculations below begin
with the additional pivot needed for a five-to-four-variable descent.
The separate [continuation note](SCHUR_DESCENT_CONTINUATIONS.md) gives the
exact reduced-determinant remainder and the simultaneous source-block
version.  Its codimension-one block-affine theorem globally excludes the
pure-source two-pivot route, including nonlinear source recharts; the cases
not covered here are therefore genuinely mixed source--dual or coisotropic.

There are two exact negative results.

1. Eliminating two constant linear combinations of the three dual variables
   can never leave a constant Hessian determinant.  This holds for every
   admissible cubic gauge parameter and in every higher seed degree.
2. A source-first Meng--Yang descent followed by a constant dual pivot also
   fails.  In degree at least four there is no nonzero constant source
   direction along which the gauge map is affine.  In degree three the only
   such direction is the original `z`-direction, but its coefficient row has
   no constant nonzero linear combination.

These are ansatz theorems, not a result about `HC(4)`.  Nonlinear symplectic
changes, nonlinear remaining dual coefficients, non-coordinate
coisotropic embeddings, and nonconstant quadratic blocks with exceptional
divisibility remain open.

## 1. The double-dual candidate

Let

\[
 F=(P,B,C):\mathbb A^3\longrightarrow\mathbb A^3
\]

be one of the determinant-minus-two quadratic-gauge maps and let

\[
 \Phi(x,u)=u_1P(x)+u_2B(x)+u_3C(x)
\]

be its six-variable Meng doubling.  Its Hessian determinant is the nonzero
constant `-4`.

Choose a constant target splitting

\[
 (P,B,C)\longmapsto(A_1,A_2,L)
\]

and put \(A=(A_1,A_2)^t\).  Let \(\Lambda\) be an invertible symmetric
constant `2 x 2` matrix.  Repair the two dual variables \(r=(r_1,r_2)^t\)
by

\[
 \widehat\Phi
 =r^tA+sL-\frac12(r-\mu)^t\Lambda^{-1}(r-\mu).
\]

The two critical equations are linear and have the polynomial solution

\[
 r=\mu+\Lambda A.
\]

Thus both successive Schur complements are polynomial.  If the selected
gauge collision lies over \(F_0\), write \(A_0=A(F_0)\) and take
\(\mu=-\Lambda A_0\).  Up to an irrelevant constant, the critical value is

\[
 \boxed{
 \psi(x,s)
 =sL(x)+\frac12(A(x)-A_0)^t\Lambda(A(x)-A_0).
 }                                                     \tag{1.1}
\]

At two source points \(x_+\ne x_-\) with \(F(x_+)=F(x_-)=F_0\), the points
\((x_+,0)\) and \((x_-,0)\) have equal gradients under (1.1).  Therefore
one instance of (1.1) with constant nonzero Hessian determinant would
immediately refute `HC(4)`.

The collision condition is consequently automatic.  Only the Hessian
determinant remains.

## 2. A necessary bordered-Hessian equation

The quadratic expression in (1.1) is only one possible polynomial
\(h(x)\).  For the more general form

\[
 \psi(x,s)=sL(x)+h(x),
\]

write \(g=\nabla L\), \(H=\operatorname{Hess}L\), and
\(R=\operatorname{Hess}h\).  Then

\[
 \operatorname{Hess}\psi
 =
 \begin{pmatrix}
 sH+R&g\\
 g^t&0
 \end{pmatrix}.
\]

The coefficient of \(s^2\) in its determinant is

\[
 \boxed{
 K(L)=-g^t\operatorname{adj}(H)g.
 }                                                     \tag{2.1}
\]

In particular, a constant Hessian determinant requires

\[
 K(L)=0.                                               \tag{2.2}
\]

This obstruction is independent of the repaired quadratic form, its
centering, and every lower term \(h\).

## 3. Exact cubic-family obstruction

Normalize the cubic seed by

\[
 \alpha=\frac{g_3}{g_1}\ne0,\qquad
 \beta=\frac{g_2}{g_1}.
\]

Put

\[
 t=1+xy,\qquad
 q=t^2z+\alpha^{-1}y^2(1+3t)
\]

and use

\[
\begin{aligned}
 P&=tq,\\
 B&=y+3\alpha xq+2\beta tq,\\
 C&=x(5-3t)-\alpha x^3z.
\end{aligned}                                         \tag{3.1}
\]

These formulas include the entire admissible two-parameter cubic gauge
family and satisfy \(\det D(P,B,C)=-2\).

For

\[
 L=pP+bB+cC,
\]

expand \(\alpha^2K(L)\) in \(x,y,z\).  Three coefficients are triangular:

\[
\begin{aligned}
 [x^8]\,\alpha^2K(L)&=9\alpha^4c^4,\\
 [z^2]\,\alpha^2K(L)&=9\alpha^2(p+2\beta b)^4,\\
 [1]\,\alpha^2K(-2\beta b,b,0)&=9\alpha^4b^4.
\end{aligned}                                         \tag{3.2}
\]

Over characteristic zero, with \(\alpha\ne0\), equation (2.2) therefore
forces

\[
 c=0,\qquad p=-2\beta b,\qquad b=0,\qquad p=0.
\]

Thus

\[
 \boxed{K(L)\ne0\quad\text{for every nonzero target linear form }L}
                                                               \tag{3.3}
\]

throughout the cubic gauge family.  No choice of \(\Lambda\) in (1.1) can
repair this obstruction.

## 4. All higher seed degrees

For a seed of exact degree \(N\ge4\), let
\(\gamma=g_N/g_1\ne0\).  The top \(z\)-degree of
\(L=pP+bB+cC\) comes only from the degree-\(N\) decoration.  Its coefficient
is \(\gamma f_N(x,y)\), where

\[
 f_N
 =x^{N-2}t^{2N}\bigl(Nbt^2-(N-2)cx^2\bigr).            \tag{4.1}
\]

For a polynomial

\[
 L=f(x,y)z^N+\text{lower powers of }z,
\]

the coefficient of \(z^{4N-2}\) in \(K(L)\) is

\[
\begin{aligned}
 I_N(f)=-Nf\bigl[
  &Nf(f_{xx}f_{yy}-f_{xy}^2)\\
  &-(N+1)(f_x^2f_{yy}-2f_xf_yf_{xy}+f_y^2f_{xx})
 \bigr].                                               \tag{4.2}
\end{aligned}
\]

The leading \(x\)-adic term already decides (4.2).  For
\(f=A x^m t^k\), its coefficient at \(y=0\) is

\[
 [x^{4m}]I_N(f)
 =-A^4kN\bigl(3km-kN+m^2+mN\bigr).                    \tag{4.3}
\]

If \(b\ne0\), use

\[
 (m,k)=(N-2,2N+2).
\]

The last factor in (4.3) is

\[
 2(3N^2-7N-4)>0\qquad(N\ge4).                          \tag{4.4}
\]

The \(c\)-term in (4.1) has two higher powers of \(x\), so it cannot cancel
this leading coefficient.  More explicitly, at \(y=0\), a term of
\(x\)-order \(r\) gives derivative orders

\[
 (f,f_x,f_y,f_{xx},f_{xy},f_{yy})
 =(r,r-1,r+1,r-2,r,r+2);
\]

every quartic summand of (4.2) has total order \(4r\).  Replacing any
\(b\)-factor by a \(c\)-factor raises that order by two.  Hence \(b=0\).

With \(b=0\), formula (4.1) is the single monomial with

\[
 (m,k)=(N,2N),
\]

and the last factor in (4.3) is \(6N^2\).  Thus \(c=0\).  It remains only
\(L=pP\), but direct restriction gives

\[
 K(P)(0,0,z)=9z^2.                                    \tag{4.5}
\]

Therefore \(p=0\), proving (3.3) in every exact degree \(N\ge4\) as well.

## 5. Source-first successive descent

The original Meng--Yang descent first eliminates a source coordinate along
which the Keller map is affine.  Let

\[
 v=a\partial_x+b\partial_y+c\partial_z
\]

be a constant source direction.  The common first gauge coordinate is
\(P=tq\).  Two coefficients of its second directional derivative are

\[
 [y^4]\,\alpha D_v^2P=6a^2,\qquad
 [1]\,\alpha D_v^2P=8b^2.                              \tag{5.1}
\]

Thus \(D_v^2F=0\) forces \(a=b=0\): only the \(z\)-direction can be affine.

For every exact seed degree \(N\ge4\), the second target coordinate contains

\[
 N\gamma t^2x^{N-2}q^N.
\]

Its second \(z\)-derivative is

\[
 N^2(N-1)\gamma t^6x^{N-2}q^{N-2}\ne0.                \tag{5.2}
\]

Hence no nonzero constant affine source direction exists in the
higher-degree maps, and a source-first Meng descent cannot start.

For the cubic family, \(z\) is affine and

\[
 \partial_z(P,B,C)
 =\bigl(t^3,\ 3\alpha xt^2+2\beta t^3,\ -\alpha x^3\bigr).
                                                               \tag{5.3}
\]

If a constant target covector \((p,b,c)\) made (5.3) a scalar unit, setting
\(y=0\) would force \(b=c=0\), and varying \(xy\) would then force \(p=0\).
There is no constant nonzero second pivot.  This recovers, uniformly in the
cubic parameters, the obstruction that the foundational coefficient row
is polynomially unimodular but not constant-linearly unimodular.

## 6. What remains worth searching

The calculation closes the most literal parameter-family extension of the
Meng--Yang construction:

- constant target splittings with two repaired dual variables;
- arbitrary constant symmetric quadratic repair on that dual block;
- every admissible root-engineered seed degree;
- constant linear source directions followed by a constant dual pivot.

The surviving search space is nonlinear.  The most plausible next tests
are:

1. a polynomial symplectic change in which the retained dual coefficient
   is not a constant linear form \(L=pP+bB+cC\);
2. a nonconstant `2 x 2` pivot whose determinant is a unit after a
   cancellation intrinsic to the gauge chart;
3. a coisotropic embedding adapted to the marked-line coordinates
   \((P,S,Q)\), while proving that both critical solutions pull back
   polynomially to \((x,y,z)\);
4. a search allowing the two Schur variables to mix source and dual
   directions nonlinearly.

For any resulting scalar affine pivot, the determinant target is the weaker
condition
\[
 \det M(\mu+\lambda A,w)\in K
\]
from Proposition 1.1 of the continuation note; the full pencil
\(\det M(s,w)\) need not vanish identically.  For a simultaneous pivot block,
the first gate is the corank bound of Theorem 3.1 there.

The existing nonlinear toric calculation for the foundational potential
shows why polynomial coordinates alone are insufficient: a nonlinear point
change can produce a unit pivot while changing the Hessian determinant and
losing the collision.  Any extension of items 1--4 must check all three
properties after the change: polynomial critical solutions, constant
Hessian determinant, and an equal-gradient pair.

## 7. First nonlinear test: triangular target coordinates

The first test in item 1 can be closed exactly through target degree three.
The cubic gauge parameters do not enlarge this calculation.  Indeed, for

\[
 \alpha=g_3/g_1\ne0,\qquad \beta=g_2/g_1,
\]

the linear source--target normalization

\[
 z'=\alpha z,\qquad
 P'=\alpha P,\qquad
 B'=B-2\beta P,\qquad
 C'=C                                             \tag{7.1}
\]

turns the cubic gauge map into the foundational representative
\((\alpha,\beta)=(1,0)\).  Thus use

\[
\begin{aligned}
 t&=1+xy,&q&=t^2z+y^2(1+3t),\\
 P&=tq,&B&=y+3xq,&
 C&=x(5-3t)-x^3z.                                  \tag{7.2}
\end{aligned}
\]

For each target permutation \((U,V,W)\), consider the triangular coordinate

\[
 L=W+H(U,V),                                       \tag{7.3}
\]

where

\[
\begin{aligned}
H={}&a_{10}U+a_{01}V+a_{20}U^2+a_{11}UV+a_{02}V^2\\
 &+a_{30}U^3+a_{21}U^2V+a_{12}UV^2+a_{03}V^3.
                                                               \tag{7.4}
\end{aligned}
\]

The necessary equation remains \(K(L)=0\).  Sparse coefficient extraction
gives short triangular contradictions in all three orientations.

For \(L=P+H(B,C)\), the coefficients at

\[
\begin{gathered}
(34,22,10),\ (32,0,8),\ (34,14,10),\ (34,6,10)
\end{gathered}
\]

in \((x,y,z)\), successively after the preceding forced zeros, are

\[
\begin{aligned}
-1549681956a_{30}^4,\qquad&
729a_{03}^4,\\
-11337408a_{21}^4,\qquad&
-61236a_{12}^4.
\end{aligned}                                      \tag{7.5}
\]

Thus the cubic part vanishes.  The two remaining coefficients

\[
 [x^8z^6]K=-1259712a_{20}^4,\qquad
 [z^2]K=-9(12a_{20}-1)                              \tag{7.6}
\]

then force simultaneously \(a_{20}=0\) and \(a_{20}=1/12\).

For \(L=B+H(P,C)\), four analogous top coefficients force

\[
 a_{30}=a_{03}=a_{21}=a_{12}=0.                     \tag{7.7}
\]

The next coefficients successively force

\[
 a_{20}=a_{02}=a_{11}=a_{10}=a_{01}=0,              \tag{7.8}
\]

after which the constant spatial coefficient of \(K\) is \(9\).

Finally, for the full cubic polynomial \(L=C+H(P,B)\), no parameter
elimination is needed:

\[
 \boxed{[x^8]K(L)=9.}                                \tag{7.9}
\]

Consequently:

> **Cubic triangular-target obstruction.**  After the canonical linear
> normalization (7.1), no target permutation followed by a triangular
> shear \(W\mapsto W+H(U,V)\) of total degree at most three can pass the
> necessary bordered-flatness equation \(K(L)=0\).  Hence none can produce
> a four-variable constant-Hessian double-dual descent.

This does not cover an arbitrary nonlinear target automorphism, a single
shear of degree at least four, a word with a higher-degree factor, or a
nonlinear source--dual symplectic change.  The quadratic length-two case is
treated next.

## 8. Two quadratic triangular target shears

Write \(X_0=P,X_1=B,X_2=C\).  For distinct \(i,j\), with \(k\) the remaining
index, consider the six ordered words

\[
\begin{aligned}
A&=X_i+Q(X_j,X_k),\\
L&=X_j+R(A,X_k),                                    \tag{8.1}
\end{aligned}
\]

where

\[
\begin{aligned}
Q(U,V)&=q_{10}U+q_{01}V+q_{20}U^2+q_{11}UV+q_{02}V^2,\\
R(U,V)&=r_{10}U+r_{01}V+r_{20}U^2+r_{11}UV+r_{02}V^2.
                                                               \tag{8.2}
\end{aligned}
\]

The second shear genuinely depends nonlinearly on the first precisely when
\(r_{20}\ne0\) or \(r_{11}\ne0\).  If both vanish, then for \(r_{10}\ne0\)

\[
r_{10}^{-1}L=X_i+
\left(r_{10}^{-1}X_j+Q(X_j,X_k)+r_{10}^{-1}R(0,X_k)\right),     \tag{8.3}
\]

which is one quadratic triangular shear.  The case \(r_{10}=0\) is even
more immediate, so Section 7 excludes this last stratum.

For \(r_{20}\ne0\), the first exact \(x\)-axis gates are

\[
\begin{array}{c|c}
(i,j,k)&\text{top coefficient of }K(x,0,0)\\ \hline
(0,1,2),(1,0,2)&9437184q_{02}^8r_{20}^4\\
(0,2,1),(1,2,0)&9437184q_{20}^8r_{20}^4\\
(2,0,1),(2,1,0)&2304r_{20}^4.
\end{array}                                                   \tag{8.4}
\]

The last two words are impossible.  In the two words retaining \(C\),
the next gates are

\[
[x^{12}]K=2304q_{10}^8r_{20}^4,\qquad [x^8]K=9.               \tag{8.5}
\]

For \((i,j,k)=(0,1,2)\), the residual equations force

\[
q_{02}=0,\quad r_{02}=-q_{01}^2r_{20}-q_{01}r_{11},\quad
q_{20}=-4,\quad q_{10}=0,                                    \tag{8.6}
\]

and then give the incompatible coefficients

\[
3r_{10}(23763r_{10}^3-920r_{20}),\qquad
-6(4005r_{10}^3-8r_{20}).                                   \tag{8.7}
\]

The other residual word reduces by exact elimination to

\[
L=P+r(B+uPC)^2,\qquad r\ne0.                                  \tag{8.8}
\]

Here

\[
[z^2]K=9-4r(8u+9)(2u+3),                                    \tag{8.9}
\]

and, after imposing (8.9),

\[
[y^{36}]K(1,y,0)=
-\frac{910050728661u^8}{4(2u+3)^4(8u+9)^4}.                  \tag{8.10}
\]

Thus \(u=0\), then \(r=1/12\), but the remaining \(y\)-axis coefficient is
\(35641/144\).

For the bilinear stratum \(r_{20}=0,\ r_{11}\ne0\), the first gates are

\[
\begin{array}{c|c}
(i,j,k)&\text{top coefficient of }K(x,0,0)\\ \hline
(0,1,2),(1,0,2)&186624q_{02}^4r_{11}^4\\
(0,2,1),(1,2,0)&2304q_{20}^4r_{10}^4\\
(2,0,1),(2,1,0)&9r_{10}^4.
\end{array}                                                   \tag{8.11}
\]

Both branches of the middle product in (8.11) are retained.  Exact axis
and transverse-line descent closes them using

\[
\begin{aligned}
48v^4+48v^3+836v^2-268v+99&>0,\\
9v^2+27v+53&>0                                               \tag{8.12}
\end{aligned}
\]

over the reals.  After the forced cancellations, the non-immediate normal
forms end with one of

\[
-\frac{1712421}{4},\qquad
-\frac{1712421}{4}r^4,\qquad
-\frac{291761109}{4}r^8,                                    \tag{8.13}
\]

with the displayed \(r\ne0\).  The checker records every branch and proves
the positivity statements by exact real-root counts.

Consequently:

> **Two-quadratic-shear obstruction.**  For the foundational normalized
> cubic gauge and real shear coefficients, none of the six ordered
> compositions (8.1)--(8.2) satisfies \(K(L)=0\).  Hence no length-two word
> in positive-degree quadratic triangular target automorphisms supplies
> the retained coordinate of a four-variable constant-Hessian double-dual
> descent.

This is an ansatz theorem, not a result about `HC(4)`.  It does not exclude
a word containing a cubic-or-higher shear, length at least three, a general
tame or wild target automorphism, or a nonlinear source--dual symplectic
transformation.  The degree-four single-shear case is treated next.

## 9. Quartic triangular target coordinates

Return to

\[
L=W+H(U,V),                                                 \tag{9.1}
\]

but now let \(H\) have arbitrary positive total degree at most four.  Write
the homogeneous quartic layer as

\[
H_4=a_{40}U^4+a_{31}U^3V+a_{22}U^2V^2+a_{13}UV^3+a_{04}V^4.
                                                               \tag{9.2}
\]

The full sparse calculation includes all nine lower-degree coefficients;
the following extreme spatial coefficients are not contaminated by them.
For \(L=P+H(B,C)\), successive substitution of the preceding forced zeros
gives

\[
\begin{array}{c|c}
\text{spatial coefficient}&\text{value}\\ \hline
[x^{46}y^{30}z^{14}]K&-396718580736a_{40}^4\\
[x^{44}z^{12}]K&2304a_{04}^4\\
[x^{46}y^{22}z^{14}]K&-3367210176a_{31}^4\\
[x^{46}y^6z^{14}]K&-139968a_{13}^4\\
[x^{46}y^{14}z^{14}]K&-25194240a_{22}^4.
\end{array}                                                   \tag{9.3}
\]

Thus \(H_4=0\).  For \(L=B+H(P,C)\), the corresponding five coefficients
are

\[
\begin{array}{c|c}
[x^{46}y^{46}z^{14}]K&-16128a_{40}^4\\
[x^{44}z^{12}]K&2304a_{04}^4\\
[x^{46}y^{34}z^{14}]K&-10800a_{31}^4\\
[x^{46}y^{10}z^{14}]K&-2736a_{13}^4\\
[x^{46}y^{22}z^{14}]K&-6336a_{22}^4.
\end{array}                                                   \tag{9.4}
\]

Again \(H_4=0\).  Section 7 then excludes both remaining cubic-or-lower
shears.  Finally, in the third orientation the quartic parameters do not
alter the earlier scalar obstruction:

\[
\boxed{[x^8]K(C+H(P,B))=9.}                                  \tag{9.5}
\]

Consequently:

> **Quartic triangular-target obstruction.**  Over every
> characteristic-zero coefficient field, no target permutation followed
> by one triangular shear of total degree at most four satisfies
> \(K(L)=0\).  Hence no such shear supplies a four-variable
> constant-Hessian double-dual descent.

Together with Section 8, this closes the two cheapest nonlinear
target-only continuations: one quartic shear and two quadratic shears.
It does not exclude a shear of degree at least five, a word containing a
cubic-or-higher factor, length at least three, a general target
automorphism, or a nonlinear source--dual transformation.

## 10. First nonconstant `2 x 2` block test

The foundational \(z\)-coefficient row is

\[
\partial_z(P,B,C)=\bigl(t^3,3xt^2,-x^3\bigr).          \tag{10.1}
\]

Although no constant combination of (10.1) is a unit, the row is
polynomially unimodular.  With \(v=xy\),

\[
(1-3v+6v^2)t^3
-x^3y^3(6v^2+15v+10)=1.                               \tag{10.2}
\]

Thus the first nonconstant unit pivot exists.  Completing (10.1) to an
\(\operatorname{SL}_3\) matrix is exactly the nonlinear source-first
reduction leading to the Meng--Yang five-variable potential analyzed in
[`HC5_NONLINEAR_TORIC_DESCENT.md`](HC5_NONLINEAR_TORIC_DESCENT.md).

After its next unit-pivot elimination, a relative polynomial correction
\(C(x,y)\in SL_2\) acts on

\[
\beta=(-yP(xy),xQ(xy)),\qquad G=\beta C,
\]

and the four-variable determinant is \(16(\det DG)^2\).  The first
factorization-independent non-toric search takes all four entries of \(C\)
to have total degree at most four.  Its sixty parameters give 45 equations
from \(\det C=1\) and 218 nonconstant-Jacobian equations.  Their exact ideal
over \(\mathbb Q\) is the unit ideal, even before collision equality is
imposed.  Starting instead from the known degree-ten toric correction,
whose coefficient map is \((y,2x)\), an arbitrary affine \(SL_2\)
perturbation also has unit ideal after the collision equations are added.

This closes raw non-toric relative corrections through degree four and
affine perturbations of the known toric solution.  It does not exclude raw
degree at least five, quadratic-or-higher perturbations of the toric
solution, transformations mixing base and dual variables, or non-coordinate
coisotropic embeddings.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_double_schur_gauge_obstruction.py
.venv/bin/python scripts/verify_hc4_triangular_target_shears.py
.venv/bin/python scripts/verify_hc4_two_quadratic_target_shears.py
.venv/bin/python scripts/verify_hc4_quartic_target_shears.py
.venv/bin/python scripts/verify_hc4_nontoric_sl2_correction_degree4.py
```

The checker verifies the universal bordered-determinant identity, the three
cubic coefficient equations, the all-degree leading-layer formula, the
source-direction classification, and the missing cubic unit pivot.  The
second checker verifies the degree-at-most-three triangular-target
obstruction by sparse exact coefficient extraction.  The third checker
verifies all six ordered length-two quadratic words, including exact
real-root counts for the two positivity gates.  The fourth checker verifies
that all five quartic coefficients vanish in the first two orientations
and that the third retains \([x^8]K=9\).  The fifth checker verifies the
bounded non-toric relative-\(SL_2\) obstruction and the affine perturbation
obstruction at the known toric correction.

## External source

The one-variable Schur-descent lemma and the five-variable foundational
example are from Guowu Meng and Liang Yang,
[*A five-variable counterexample to the Hessian conjecture, and the
low-dimensional status of the Jacobian and Hessian
conjectures*](https://arxiv.org/abs/2607.22198).
