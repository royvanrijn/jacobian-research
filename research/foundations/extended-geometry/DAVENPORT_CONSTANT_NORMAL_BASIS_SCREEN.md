# Davenport constant-normal basis screen

The surviving pure-Davenport conic branch from
[the quadratic-survivor calculation](DAVENPORT_QUADRATIC_SURVIVORS.md)
has

\[
\mathbf b=(1,g,g^2)
\]

and a nonzero one-variable normal coefficient \(R(g)\).  This note begins
the residual calculation with

\[
R=c\in K^*
\]

and the natural primitive generator

\[
v=Y.
\]

The calculation is exact but does not eliminate reduced square roots with
\(Y^5\) or \(Y^6\) components for arbitrary \(c\).  The normalized slice
\(c=1\) is eliminated completely by the \(J_1\)-fiber certificate below.

## 1. The degree-seven algebra

Put \(z=g_T(Y)\) and work over

\[
F=K(T,z),\qquad
E=F[Y]/(g_T(Y)-z).
\]

The classes

\[
1,Y,\ldots,Y^6
\]

form the standard \(F\)-basis of the degree-seven Davenport field.

For constant \(R=c\), the discriminant equation from the conic calculation
becomes

\[
\boxed{
W^2
=16v^4+8cBv+4c^2T-2cL,
}
\tag{1}
\]

where

\[
B=S'(z)-2(T+1)z,\qquad
W=cx+4v^2.
\]

With \(v=Y\), equation (1) is

\[
W^2\equiv
16Y^4+8cB\,Y+C
\pmod{g_T(Y)-z},
\tag{2}
\]

where \(C=4c^2T-2cL(z)\).

Writing

\[
W=w_0+w_1Y+\cdots+w_6Y^6
\]

and reducing (2) gives seven explicit quadratic equations over \(F\).
This is the full unbounded Davenport-basis square-root system for the
choice \(v=Y\).

## 2. Basis degree at most three

Specialize the necessary identity to \(z=0\).  Then

\[
B=B(0)\in K,\qquad
C=4c^2T-2cL(0).
\tag{3}
\]

If the reduced degree of \(W\) is at most three, its square has degree at
most six and no reduction by the degree-seven equation occurs.  Comparing
coefficients in

\[
W^2=16Y^4+8cB(0)Y+C
\]

gives successively

\[
w_3=0,\qquad
w_2=\pm4,\qquad
w_1=w_0=0,\qquad
B(0)=0,\qquad C=0.
\]

But

\[
\frac{\partial C}{\partial T}=4c^2\ne0,
\]

whereas \(L(0)\) is constant.  Hence:

\[
\boxed{\deg_YW\le3\text{ is impossible}.}
\tag{4}
\]

## 3. Basis degree four

Suppose

\[
W=w_4Y^4+w_3Y^3+w_2Y^2+w_1Y+w_0,
\qquad w_4\ne0.
\]

Since \(W^2-(16Y^4+8cB(0)Y+C)\) is divisible by
\(7g_T(Y)\), write its quotient as \(h_1Y+h_0\).  Put

\[
p=w_4,\qquad t=\frac{w_3}{w_4},\qquad k=p^2.
\]

The coefficients of \(Y^8,Y^7,Y^6,Y^5,Y^4\) determine

\[
h_1,h_0,w_2,w_1,w_0
\]

in terms of \(p,t\).  The remaining \(Y^3\) and \(Y^2\) equations are
linear in \(k\).  Eliminating \(k\) gives, up to a nonzero factor, the
quintic

\[
\begin{aligned}
Q_5(t)={}&t^5+\frac13t^4
+T\left(2+\frac43a\right)t^3\\
&+T\left(-4-\frac23a\right)t^2\\
&+\left[
T^2\left(-\frac{23}{3}+a\right)
+T\left(1+\frac a3\right)
\right]t\\
&+T^2\left(7+\frac73a\right).
\end{aligned}
\tag{5}
\]

Exact factorization over

\[
K(T)=\mathbb Q(a)(T),\qquad a^2+a+2=0,
\]

shows that \(Q_5\) is irreducible.  In particular it has no root
\(t\in K(T)\).  Since \(t=w_3/w_4\) would be rational over \(K(T)\), this
is impossible:

\[
\boxed{\deg_YW=4\text{ is impossible}.}
\tag{6}
\]

The obstruction is independent of the lower coefficients \(B(0),C\) and
therefore independent of the choices of \(S,L,c\).

## 4. Exact residual gate

For the natural generator \(v=Y\), every remaining constant-normal
candidate must have

\[
\boxed{\deg_YW\in\{5,6\}}
\tag{7}
\]

in its reduced Davenport-basis representative.

These two cases are genuinely different from the discarded low-basis
branches.  Their squares use the degree-seven relation in several
coefficients, and the unrestricted seven-equation ideal is
positive-dimensional before the conditions

\[
B=S'(z)-2(T+1)z,\qquad
C=4c^2T-2cL(z)
\]

and \(J_1=0\) are imposed.  A bounded Gröbner search here would not be a
proof of the unrestricted statement.

The next exact calculation should therefore retain the derivation:
reduce the \(J_1=0\) equation in the same seven-element basis and combine
it with the degree-five and degree-six square systems.  The present result
removes all lower-basis components before that differential elimination.

## 5. The \(J_1\) fiber closes \(c=1\)

Set \(c=1\), but retain square roots of all reduced basis degrees.  It is
more efficient to use the function-field coordinates \((T,z)\).  The
implicit generator \(v=Y(T,z)\) satisfies

\[
v_T=-\frac{g_T}{g_Y},\qquad
v_z=\frac1{g_Y}.
\tag{8}
\]

After differentiating the square equation \(W^2=\mathcal D(Y)\), the
\(J_2\)-coefficient vanishes identically.  Substitution in \(J_1\) leaves
an equation linear in \(W\):

\[
\mathcal A\,W+\mathcal C=0
\tag{9}
\]

in the Davenport field.

Any global solution must survive the ordinary specialization

\[
(T,z)=(1,0).
\]

The degree-seven polynomial \(g_1(Y)\) is separable there.  Put

\[
s_1=S'(0),\quad s_2=S''(0),\quad
\ell_0=L(0),\quad\ell_1=L'(0).
\]

Eliminating \(W\) between (9) and

\[
W^2=16Y^4+8s_1Y+4-2\ell_0
\]

gives

\[
\mathcal C^2-\mathcal A^2
\left(16Y^4+8s_1Y+4-2\ell_0\right)=0
\]

in

\[
K[Y]/(g_1(Y)).
\]

Its seven basis coefficients are quartic equations in
\((s_1,s_2,\ell_0,\ell_1)\).  Exact Gröbner reduction over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0,
\]

returns

\[
\boxed{\{1\}.}
\tag{10}
\]

Thus the necessary fiber is empty.  This includes the previously
unresolved degree-five and degree-six square roots:

\[
\boxed{\text{the complete }R=1,\ v=Y\text{ branch is impossible}.}
\tag{11}
\]

Treating \(c\) as a symbolic coefficient makes the same Gröbner
calculation substantially larger.  Equation (11) is therefore a theorem
for the \(c=1\) slice, not a normalization of all \(c\in K^*\).

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_constant_normal_basis.py
.venv/bin/python scripts/verify_davenport_constant_normal_j1_fiber.py
```

The checker constructs all seven basis equations, proves the
degree-at-most-three obstruction, derives the degree-four resultant, and
verifies irreducibility of the quintic (5) over \(K(T)\).  The second
checker constructs the specialized \(J_1\) system and verifies the unit
Gröbner basis (10).
