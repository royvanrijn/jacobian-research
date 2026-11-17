# Davenport quadratic-in-\(U\) incidence

The independent marked-line programme leaves genuinely nonlinear dependence
on an auxiliary variable as its first change of mechanism.  The smallest
constant-direction quadratic ansatz has an exact obstruction.

For either Davenport polynomial \(g_T(Y)\), consider

\[
\boxed{
\Phi(T,Y,U)=
\left(
T+\alpha U^2,\;
g_T(Y)+\beta U+\delta U^2,\;
h(T,Y)+k(T,Y)U+\gamma U^2
\right),
}
\tag{1}
\]

where

\[
\alpha,\beta\ne0,\qquad
\delta,\gamma\in K,\qquad h,k\in K[T,Y].
\]

No output in (1) exposes \(T\), \(Y\), or \(U\).  The old parameter and
cover value are both replaced, and the third output carries a
sheet-dependent coupling.  Thus (1) lies beyond the affine-mask and
retained-parameter classes already excluded.

Nevertheless:

> **Quadratic-incidence no-go.**  The map (1) cannot have nonzero constant
> Jacobian.

The same statement holds for the conjugate line-cover polynomial.

## 1. The three determinant equations

Write subscripts for partial derivatives.  Direct expansion gives

\[
\det D\Phi=J_0+J_1U+J_2U^2,
\]

with

\[
J_0=g_Yk-\beta h_Y,
\tag{2}
\]

\[
J_1=
2\gamma g_Y-\beta k_Y-2\delta h_Y
+2\alpha(g_Th_Y-g_Yh_T),
\tag{3}
\]

\[
J_2=
-2\delta k_Y
+2\alpha(g_Tk_Y-g_Yk_T).
\tag{4}
\]

Assume \(\det D\Phi=C\in K^*\).  Equation (4) becomes

\[
J\left(g-\frac{\delta}{\alpha}T,k\right)=0.
\tag{5}
\]

## 2. The shifted Davenport polynomial is closed

For every \(\lambda\in\overline K\), put

\[
f_\lambda=g-\lambda T.
\]

As a polynomial in \(Y\),

\[
[Y^7]f_\lambda=\frac17,\qquad
[Y^6]f_\lambda=0,\qquad
[Y^5]f_\lambda=(1+a)T.
\tag{6}
\]

If \(f_\lambda=G(\varphi)\) were a nontrivial polynomial composite, then
the prime \(Y\)-degree seven would force

\[
\deg G=7,\qquad \deg_Y\varphi=1.
\]

The constant \(Y^7\)-coefficient forces the \(Y\)-slope of \(\varphi\) to
be constant.  The zero \(Y^6\)-coefficient then forces its translation to
be constant.  But the \(Y^5\)-coefficient of \(G(\varphi)\) would also be
constant, contradicting (6).

Hence every \(f_\lambda\) is closed.  Its polynomial Jacobian centralizer
is therefore

\[
\ker J(f_\lambda,-)=K[f_\lambda].
\tag{7}
\]

From (5), there is a polynomial \(K_0\) such that

\[
k=K_0(f_{\delta/\alpha}).
\tag{8}
\]

## 3. The remaining coefficient is impossible

Equation (2) gives

\[
h_Y=\frac{k g_Y-C}{\beta}.
\]

Choose \(L'=K_0\).  Integration in \(Y\) yields

\[
h=
\frac{L(f_{\delta/\alpha})-CY}{\beta}+c(T).
\tag{9}
\]

Substituting (8)--(9) into (3) cancels every undifferentiated \(K_0\)-term
and gives

\[
\boxed{
\Theta(T,f_{\delta/\alpha})g_Y
-2\alpha Cg_T
+2\delta C=0,
}
\tag{10}
\]

where

\[
\Theta=
2\beta\gamma-\beta^2K_0'(f_{\delta/\alpha})
-2\alpha\beta c'(T).
\]

If \(K_0'\) is nonconstant, the first term in (10) has \(Y\)-degree at
least

\[
7+6=13,
\]

whereas \(g_T\) has \(Y\)-degree five.  This is impossible.

If \(K_0'\) is constant, then \(\Theta\in K[T]\).  The monic
\(Y^6\)-term of \(g_Y\) forces \(\Theta=0\).  Equation (10) would then say

\[
g_T=\frac{\delta}{\alpha},
\]

but

\[
[Y^5]g_T=1+a\ne0.
\]

This is the final contradiction.

If \(\beta=0\), equation (2) already reads \(g_Yk=C\), which is impossible
in \(K[T,Y]\).  Thus the nonzero-\(\beta\) normalization loses no candidate
inside this constant-direction class.

## 4. Consequence

The first nonlinear auxiliary suspension is closed before monodromy,
boundary, or Hessian calculations are needed.  A surviving quadratic
construction must change at least one structural feature of (1):

1. use a nonconstant \(U^2\)-coefficient depending on \(T,Y\);
2. use two independent auxiliary variables;
3. move beyond a single quadratic direction in target space; or
4. combine the nonlinear suspension with a non-elementary reciprocal
   source modification.

The smallest next coefficient ansatz is obtained by replacing
\(\delta,\gamma\) with low-degree polynomials.  Its highest powers of \(U\)
should be tested first: they determine whether the shifted-centralizer
reduction survives before the lower Keller equations are expanded.

## 5. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_quadratic_u_incidence.py
```

The checker verifies the universal determinant coefficients, the exact
Davenport derivative degrees and primitivity data, and the final degree
contradiction.
