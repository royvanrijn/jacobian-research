# The affine \(A_4\) Keller frontier

## 1. Outcome

The oriented-quartic checkpoint gave a determinant-one \(A_4\) inverse cover
on a discriminant complement. Two further reductions now remove the
orientation denominator and put the remaining absolute problem into
polynomial affine coordinates.

1. Squaring on the oriented cubic quotient gives a polynomial degree-four
   map with generic \(A_4\) monodromy and constant hypersurface-residue
   Jacobian \(4\).
2. Rational invariant coordinates turn the same cover into a rational map
   of affine planes with a completely factored Jacobian ledger.
3. Homogenizing its common denominator produces a literal polynomial map
   \(\mathbb A^3\to\mathbb A^3\) with generic monodromy \(A_4\).

The polynomial affine map is not yet Keller: its determinant is
\(4W^2K^3L\). The construction therefore reaches prescribed monodromy and
affine-space algebraization simultaneously, and isolates constant-Jacobian
completion as the only remaining gate.

## 2. Oriented cubic quotient

Let

\[
 \delta(x,y)=x^2y^2-4x^3-4y^3+18xy-27
\]

and define the affine oriented cubic surface

\[
 \mathcal S=\{d^2=\delta(x,y)\}\subset\mathbb A^3.    \tag{2.1}
\]

If \(a,b,c\) are the roots of

\[
 Z^3-xZ^2+yZ-1,
\]

then \(abc=1\), \(x=a+b+c\), \(y=ab+bc+ca\), and \(d\) is the Vandermonde
orientation. Thus \(\mathcal S\) is the affine quotient of
\((\mathbb G_m)^2\) by the cyclic permutation

\[
 (a,b,c)\longmapsto(b,c,a).
\]

Squaring the three roots gives

\[
\boxed{
\begin{aligned}
X&=x^2-2y,\\
Y&=y^2-2x,\\
D&=d(xy-1).
\end{aligned}}                                      \tag{2.2}
\]

The discriminant identity is

\[
 \delta(X,Y)=\delta(x,y)(xy-1)^2,                    \tag{2.3}
\]

so (2.2) is a polynomial self-map of \(\mathcal S\).

The kernel of squaring on \(abc=1\) is \(V_4\). Cyclic permutation acts
transitively on its three nonidentity elements, so the normal closure group
is

\[
 V_4\rtimes C_3=A_4.
\]

After quotienting source and target by the common \(C_3\), (2.2) has degree
four and monodromy \(A_4\) in its natural four-point action.

On \(dD\ne0\), use the residue form

\[
 \omega_{\mathcal S}=\frac{dx\wedge dy}{2d}.
\]

Since

\[
 \det\frac{\partial(X,Y)}{\partial(x,y)}
 =4(xy-1),
\]

equation (2.2) gives

\[
 \boxed{F^*\omega_{\mathcal S}=4\omega_{\mathcal S}.} \tag{2.4}
\]

Thus the orientation absorbs the full Jacobian divisor polynomially, with no
derivative-unit suspension and no \(D^{-2}\) coordinate.

The inverse \(x\)-coordinate satisfies

\[
 T^4-2XT^2-8T+X^2-4Y=0,                              \tag{2.5}
\]

whose discriminant is

\[
 4096\,\delta(X,Y)=(64D)^2.                          \tag{2.6}
\]

For example, the target

\[
 (X,Y,D)=(-6,3,27)
\]

has the complete fiber polynomial

\[
 T^4+12T^2-8T+24.
\]

Both this quartic and its cubic resolvent are irreducible, while its
discriminant is \(1728^2\); hence it is an explicit \(A_4\) fiber.

The limitation is geometric. The surface \(\mathcal S\) has the three
quotient singularities corresponding to triple roots

\[
 (x,y,d)=(3\zeta,3\zeta^2,0),\qquad \zeta^3=1.
\]

Therefore (2.2) is a constant-residue-Jacobian Cox map, not a polynomial
Keller map of smooth affine space.

## 3. An affine two-parameter \(A_4\) polynomial

[Jensen--Ledet--Yui, Theorem
2.2.9](https://library.slmath.org/books/Book45/files/book45.pdf)
give a generic \(A_4\) polynomial over the rational
\((\alpha,\beta)\)-plane. Put

\[
\begin{aligned}
A={}&\alpha^3-\beta^3-9\beta^2-27\beta-54,\\
B={}&\alpha^3-3\alpha\beta^2+2\beta^3-9\alpha\beta
     +9\beta^2-27\alpha+27\beta+27,\\
C={}&\alpha^3-\beta^3+27.
\end{aligned}                                       \tag{3.1}
\]

After scaling their root coordinate by \(B\), the monic polynomial becomes

\[
\boxed{
\begin{aligned}
P_{\alpha,\beta}(Z)
={}&Z^4-6ABZ^2-8B^3Z\\
   &+B^2(9A^2-12CB).
\end{aligned}}                                      \tag{3.2}
\]

All coefficients are polynomial in \(\alpha,\beta\). Its discriminant is the
literal square

\[
\begin{aligned}
\operatorname{Disc}(P_{\alpha,\beta})
=\bigl[
 &1728(\beta^2+3\beta+9)B^4\\
 &\cdot(2\alpha^3\beta+3\alpha^3-3\alpha^2\beta^2
 -9\alpha^2\beta-27\alpha^2\\
 &\qquad+\beta^4+6\beta^3+27\beta^2+54\beta+81)
\bigr]^2.                                           \tag{3.3}
\end{aligned}
\]

This removes the oriented-discriminant hypersurface from the target: generic
\(A_4\) monodromy already lives over an affine plane.

## 4. Rational source coordinates

Let \(C_3\) act on \(\mathbb Q(s,t)\) by

\[
 (s,t)\longmapsto\left(t,\frac1{st}\right).
\]

The fixed field is rational. Write its standard generators as \(U,V\).
Applying the same generators to \((s^2,t^2)\) gives the target parameters
\(\alpha,\beta\). Exact elimination yields

\[
 \alpha=\frac{N_1(U,V)}{H(U,V)},\qquad
 \beta=\frac{N_2(U,V)}{H(U,V)},                      \tag{4.1}
\]

where

\[
\begin{aligned}
H={}&8U^3-6UV^2-18UV-54U\\
   &\quad-2V^3-9V^2-27V-27,\\
K={}&4U^2+4UV+6U+V^2+3V+9,\\
M={}&U^2+2V^2+6V+18,\\
L={}&U^3-3UV^2-9UV-27U\\
   &\quad+2V^3+9V^2+27V+27,
\end{aligned}                                       \tag{4.2}
\]

\[
 N_1=MK,                                             \tag{4.3}
\]

and

\[
\begin{aligned}
N_2={}&8U^3V+12U^2V^2+36U^2V+108U^2\\
     &+6UV^3+36UV^2+108UV+162U\\
     &+V^4+9V^3+27V^2+54V.
\end{aligned}                                       \tag{4.4}
\]

The rational Jacobian factors completely:

\[
\boxed{
\det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
=\frac{4K^3L}{H^3}.
}                                                    \tag{4.5}
\]

This is the exact divisor ledger for affine \(A_4\) Kellerization. The three
roles are separated:

- \(H=0\) is the common pole boundary;
- \(K=0\) is a multiplicity-three Jacobian divisor;
- \(L=0\) is the remaining reduced Jacobian divisor.

## 5. A polynomial affine-space map with \(A_4\) monodromy

Introduce \(W\) and homogenize the common denominator:

\[
\boxed{
\Phi(U,V,W)=
\bigl(WN_1(U,V),\,WN_2(U,V),\,WH(U,V)\bigr).
}                                                    \tag{5.1}
\]

On the target chart with third coordinate nonzero, the first two coordinate
ratios recover (4.1). Once \(U,V\) are chosen, \(W\) is reconstructed
uniquely. Hence

\[
 \operatorname{gdeg}(\Phi)=4
\]

and its generic inverse monodromy is the same natural \(A_4\) action.

The determinant is

\[
\boxed{
\det D\Phi=4W^2K^3L.
}                                                    \tag{5.2}
\]

Thus (5.1) is an absolute polynomial \(\mathbb A^3\)-map with the prescribed
generic monodromy, but it is not Keller. Compared with the original
fiber-transfer problem, the gap has been reduced to cancelling one explicit
four-column divisor ledger without altering the function-field extension.

## 6. Defect-lift obstruction

There is an obvious alternative: thicken the oriented surface inside
\(\mathbb A^3_{x,y,d}\). Put

\[
 R=d^2-\delta(x,y)
\]

and preserve (2.2) modulo \(R\) by considering

\[
\begin{aligned}
\widetilde X&=X+Rf,\\
\widetilde Y&=Y+Rg,\\
\widetilde D&=D+Rh.
\end{aligned}                                       \tag{6.1}
\]

The affine-linear coefficient equations give the unit ideal over
\(\mathbb Q\), even with the desired constant determinant free.  The
[follow-up ledger reduction](A4_LEDGER_REDUCTION_AND_RIGIDITY.md) strengthens
this to every polynomial degree.  At \((-1,-1,0)\), the derivative of the
unmodified map has rank one, and the derivative of every correction
divisible by \(R\) lies in the same one-dimensional row direction.
Therefore every map (6.1), for arbitrary polynomial \(f,g,h\), has zero
Jacobian at that point.

Thus the defect-multiple ambient route is closed.  A viable lift must work
directly with the cone ledger (5.2) or change the oriented map away from the
defect-preserving class (6.1).

## 7. Next attack

Only the cone-ledger route survives.  The
[exact target-ledger identity](A4_LEDGER_REDUCTION_AND_RIGIDITY.md)

\[
 \mathcal B(\Phi)=W^3K^3L^2,\qquad
 \frac{\mathcal B(\Phi)}{\det D\Phi}=\frac{WL}{4}
\]

absorbs the entire \(K^3\) column and one copy of \(L\).  The remaining
construction problem is a genuinely coupled affine modification for the
two-factor residual boundary \(WL\).  Polynomial target reparametrization
inside the present cone model and block-triangular stabilization are both
proved insufficient in that follow-up.
The subsequent
[pure-target lift](A4_PURE_TARGET_LEDGER_LIFT.md) adjoins one coordinate and
makes the fourfold determinant exactly the pullback of the target
\(\mathcal B\)-divisor while preserving the generic \(A_4\) extension.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_affine_keller_frontier.py
```

The checker verifies the oriented discriminant identity, constant residue
Jacobian, inverse quartic, scaled generic-polynomial discriminant square,
rational invariant formulas, the complete Jacobian ledger, the polynomial
cone determinant, and the affine-linear ambient-lift unit ideal.
