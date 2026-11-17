# Davenport quadratic survivors

The full quadratic-vector calculation leaves two classes:

\[
\mathbf F(T,Y,U)
=\mathbf h(T,Y)+U\mathbf a(T,Y)+U^2\mathbf b(T,Y),
\]

with either the affine through-origin direction

\[
\mathbf b=(T,0,0),
\]

or a nonlinear projectively developable direction.  This note performs
the next calculation in both classes.

## 1. The through-origin equation is a Laurent zero-Jacobian equation

Write

\[
\mathbf h=(T,g,H),\qquad
\mathbf a=(a_1,a_2,a_3),\qquad
\mathbf b=(T,0,0).
\]

The highest surviving coefficient is

\[
2T\,J(a_2,a_3)-a_2(a_3)_Y+a_3(a_2)_Y=0.
\tag{1}
\]

Make the quadratic base change and rescale the auxiliary coordinate:

\[
T=s^2,\qquad u=sU,
\]

and put

\[
\alpha=\frac{a_1(s^2,Y)}s,\qquad
A=\frac{a_2(s^2,Y)}s,\qquad
B=\frac{a_3(s^2,Y)}s.
\]

The map becomes the Laurent-polynomial map

\[
\widetilde{\mathbf F}
=
\left(
s^2+\alpha u+u^2,\;
G+Au,\;
K+Bu
\right),
\tag{2}
\]

where \(G=g(s^2,Y)\) and \(K=H(s^2,Y)\).  The source change has
Jacobian two:

\[
\det\frac{\partial(T,Y,U)}{\partial(s,Y,u)}=2.
\]

The coefficient of \(u^3\) in (2) is

\[
\boxed{2J_{s,Y}(A,B).}
\tag{3}
\]

Consequently (1) is equivalent to

\[
J_{s,Y}(A,B)=0.
\tag{4}
\]

Thus \(A\) and \(B\) are algebraically dependent in
\(K[s,s^{-1},Y]\).  At the function-field level, Lüroth gives a rational
parameter \(w\) such that

\[
A=P(w),\qquad B=Q(w)
\tag{5}
\]

for rational one-variable functions \(P,Q\).  The original polynomiality
condition is the additional parity constraint

\[
A(-s,Y)=-A(s,Y),\qquad B(-s,Y)=-B(s,Y).
\tag{6}
\]

This is the structural classification of the projective PDE: its
solutions are one-parameter Laurent curves with an odd descent datum.

## 2. The next coefficient also integrates

Assume first that (5) is represented inside the Laurent ring.  Put

\[
\Delta=P Q'-Q P'.
\]

The coefficient of \(u^2\) in (2) is

\[
\Delta J(w,\alpha)
+2P'J(w,K)
-2Q'J(w,G).
\tag{7}
\]

Where \(P'\ne0\), equation (7) is

\[
\boxed{
J\left(
w,\;
K+\frac{\Delta}{2P'}\alpha-\frac{Q'}{P'}G
\right)=0.
}
\tag{8}
\]

Thus the two top determinant equations have become two successive
one-variable centralizer equations.

### 2.1 Primitive monomial parameters are impossible

Take a primitive monomial parameter

\[
w=s^rY^m,\qquad r,m>0,\qquad \gcd(r,m)=1,
\tag{9}
\]

with \(r\) odd.  The kernel of \(J(w,-)\) in
\(K[s,s^{-1},Y]\) is exactly \(K[w]\): on a monomial
\(s^iY^j\), the Jacobian vanishes precisely when

\[
rj=mi.
\]

Let \(P,Q\) be odd polynomials, as required by (6), and suppose their
projective ratio is nonconstant.  After ordering their initial odd orders,
the two coefficients in (8) have the expansions

\[
\frac{Q'}{P'}=\ell_0+O(w^2),\qquad
\frac{\Delta}{2P'}=O(w^3).
\tag{10}
\]

Equation (8) and parity therefore give

\[
H=\ell_0g+c+Y^2M(T,Y)
\tag{11}
\]

for a polynomial \(M\); in fact the residual vanishes to order at least
\(2m\) in \(Y\).  The constant determinant coefficient is

\[
J_0=-a_1J(H,g)-a_2H_Y+a_3g_Y.
\tag{12}
\]

Both \(a_2\) and \(a_3\) are divisible by \(Y\), while (11) makes
\(J(H,g)\) divisible by \(Y\).  Hence

\[
Y\mid J_0.
\]

It cannot be a nonzero constant.  Therefore:

\[
\boxed{
\text{Every nonconstant primitive monomial common-parameter branch
of (1) is impossible.}
}
\tag{13}
\]

The smallest example is

\[
w=sY,\quad A=w,\quad B=w^3,
\]

or

\[
a_2=TY,\qquad a_3=T^2Y^3.
\]

Equation (8) integrates to

\[
H=3TY^2g-TY^3a_1+R(TY^2),
\]

and direct substitution visibly gives \(Y\mid J_0\).

What remains in the through-origin class is therefore a non-monomial
Lüroth parameter, a rational parameter whose Laurent presentation has
nontrivial boundary, or the constant-projective-direction branch.

## 3. The nonlinear developable equation

Now write a general nonlinear projective direction locally as

\[
\mathbf b=\lambda(T,Y)\,\mathbf v(q(T,Y)),
\tag{14}
\]

where \(\mathbf v(q)\) is a polynomial or rational curve in target
three-space.  Put

\[
\phi=[\mathbf v,\mathbf v',\mathbf a].
\tag{15}
\]

The \(U^4\)-coefficient is exactly

\[
J_4
=\lambda\left(
-J(q,\lambda)\phi
+2\lambda J(q,\phi)
\right).
\tag{16}
\]

Away from \(\lambda\phi=0\), this is equivalent to

\[
\boxed{
J\left(q,\frac{\phi^2}{\lambda}\right)=0.
}
\tag{17}
\]

Thus the normal component of the linear direction is not free:

\[
\phi^2=\lambda R(q)
\tag{18}
\]

over the generic function field.  There are two branches:

1. the tangent branch \(\phi=0\), where
   \(\mathbf a\in\langle\mathbf v,\mathbf v'\rangle\);
2. the transverse branch (18), which forces the square class of the
   radial scale \(\lambda\) to come from the curve parameter \(q\).

Equation (17) is the nonlinear analogue of the Laurent dependence
equation (4).

## 4. The pure Davenport conic

The first nonlinear curve worth testing is

\[
\mathbf v(q)=(1,q,q^2),\qquad \lambda=1.
\tag{19}
\]

Here

\[
\phi=q^2a_1-2qa_2+a_3.
\]

Take the parameter to be the Davenport output itself:

\[
q=g(T,Y).
\]

Since every shifted Davenport polynomial is closed, its polynomial
centralizer is \(K[g]\).  Equation (17) gives

\[
\phi=R(g).
\]

Use the moving frame of the conic:

\[
\begin{aligned}
\mathbf a
={}&x(1,g,g^2)+v(0,1,2g)+R(g)(0,0,1),\\
a_2={}&gx+v,\\
a_3={}&g^2x+2gv+R(g).
\end{aligned}
\tag{20}
\]

### 4.1 The \(U^3\)-equation

Direct calculation gives

\[
J_3=2J(g,\Phi),
\]

where

\[
\Phi
=H+Tg^2-v^2-R'v-\frac R2x.
\]

Closedness of \(g\) therefore integrates this coefficient:

\[
\boxed{
H=-Tg^2+v^2+R'v+\frac R2x+S(g).
}
\tag{21}
\]

### 4.2 The \(U^2\)-equation

After (21), the next coefficient is again exact:

\[
J_2=-\frac12J(g,\Psi),
\tag{22}
\]

with

\[
\begin{aligned}
\Psi={}&
-\frac R2x^2-4v^2x-2R'vx\\
&+\left(-8Tg-8g+4S'\right)v
+2R''v^2\\
&+T(2R-4gR').
\end{aligned}
\tag{23}
\]

Hence

\[
\boxed{\Psi=L(g)}
\tag{24}
\]

for a third one-variable polynomial \(L\).

Equations (20), (21), and (24) solve

\[
J_5=J_4=J_3=J_2=0
\]

for this nonlinear conic direction.  The remaining equations are
\(J_1=0\) and \(J_0\in K^*\).

### 4.3 The \(R=0\) branch is impossible

When \(R=0\), equation (24) is

\[
-4v^2x+4\bigl(S'-2(T+1)g\bigr)v=L(g).
\tag{25}
\]

Put

\[
B_0=S'-2(T+1)g.
\]

For \(v\ne0\), equation (25) gives

\[
x=\frac{B_0}{v}-\frac{L}{4v^2}.
\tag{26}
\]

After this substitution, the \(U\)-coefficient is again exact:

\[
J_1=J(g,\Omega),
\tag{27}
\]

where

\[
\boxed{
\Omega=
\left(2T+4-S''\right)v^2
+\frac{L'}2v
+\frac L4x
+2TgS'
-2(T+1)^2g^2.
}
\tag{28}
\]

Thus \(J_1=0\) gives \(\Omega=M(g)\).  Eliminating \(x\) between
(25) and (28), and multiplying by \(16v^2\), gives the quartic

\[
\begin{aligned}
0={}&
16(2T+4-S'')v^4
+8L'v^3\\
&+16\left(
2TgS'-2(T+1)^2g^2-M
\right)v^2\\
&+4LB_0v-L^2.
\end{aligned}
\tag{29}
\]

If \(v\notin K(T,g)\), then \(v\) generates the prime degree-seven
Davenport extension, so its minimal polynomial over \(K(T,g)\) has
degree seven.  It cannot satisfy the nonzero quartic (29).  The quartic
cannot vanish coefficientwise because its leading coefficient contains
the term \(32T\), whereas \(S''\) depends only on \(g\).

If \(v\in K(T,g)\), equation (26) also puts \(x\) in the base field and
the map has the fatal \(g_Y\) factor.  Finally, the excluded division
case \(v=0\) has

\[
J_1=-2\bigl(S'-2(T+1)g\bigr)g\,g_Y,
\]

which cannot vanish identically.  Therefore

\[
\boxed{\text{the entire }R=0\text{ pure-Davenport conic branch fails}.}
\tag{30}
\]

## 5. What this opens

If \(x,v\in K[T,g]\), every output is a polynomial in \((T,g,U)\).
The chain rule then gives

\[
\det D\mathbf F
=g_Y\det\frac{\partial\mathbf F}{\partial(T,g,U)},
\]

so this subbranch is impossible.

Therefore a genuine solution of (24) must use at least one of \(x,v\)
outside \(K[T,g]\).  Because the Davenport extension

\[
K(T,g)\subset K(T,Y)
\]

has prime degree seven, any such element generates the full degree-seven
field.  This is the first surviving quadratic-vector gate that is
intrinsically aligned with the desired monodromy rather than merely
masking its derivative.

The next calculation is now precise: analyze polynomial points
\((x,v)\in K[T,Y]^2\) on the cubic relation (24), then impose \(J_1=0\).
The factor-through component is already excluded; every remaining
component must itself carry the full Davenport field.

The remaining conic case has \(R\ne0\).  Equation (24) is then quadratic
in \(x\).  If \(v\) were in \(K(T,g)\), that quadratic would put \(x\) in
an extension of degree at most two, impossible inside a degree-seven
extension unless \(x\) were already in the base.  Consequently every
genuine \(R\ne0\) survivor must have

\[
\boxed{K(T,g,v)=K(T,Y).}
\tag{31}
\]

More explicitly, the discriminant of (24) as a quadratic in \(x\) is

\[
\begin{aligned}
\mathcal D(v)={}&
16v^4+16R'v^3
+4(RR''+(R')^2)v^2\\
&+8R\bigl(S'-2(T+1)g\bigr)v\\
&+4R^2T-8RR'Tg-2LR.
\end{aligned}
\tag{32}
\]

It must be a square in \(K(T,Y)\):

\[
\mathcal D(v)
=\left(Rx+4v^2+2R'v\right)^2.
\tag{33}
\]

The next gate is therefore completely explicit: solve the degree-seven
polynomial-point problem (32)--(33), then impose the \(J_1=0\)
differential equation in the \((T,g)\) coordinates.  Unlike the closed
\(R=0\) branch, this discriminant is a genuine quartic and does not by
itself give a degree contradiction.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_quadratic_survivors.py
```

The checker verifies the Laurent base change, both through-origin top
coefficients, the primitive monomial divisor obstruction, the general
nonlinear equation (17), the exact integrations (21)--(24), and the
quartic discriminant (32), together with the degree-four obstruction
(29) for \(R=0\).
