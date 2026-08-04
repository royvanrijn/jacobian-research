# The plane characteristic-two map has no Keller lift modulo four

## 1. Statement

Let

\[
\begin{aligned}
 P&=x+x^2y+x^4+x^6y^2,\\
 Q&=y+x^5+x^6y+x^7y^2+x^8y^3
\end{aligned}
\]

over \(\mathbb F_2\).  This is the determinant-one, noninjective plane map
from the
[characteristic-two audit](HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md).

**Theorem.**  There are no polynomials
\(\widetilde P,\widetilde Q\in(\mathbb Z/4)[x,y]\) reducing to \(P,Q\)
whose Jacobian determinant is constant.  In particular this map has no
Keller lift through \(W_2(\mathbb F_2)=\mathbb Z/4\), and hence no compatible
lift through any \(W_n(\mathbb F_2)\), through \(\mathbb Z_2\), or to
characteristic zero.

The obstruction is unstable: after adjoining one identity coordinate, the
map \((P,Q,z)\) has an explicit compatible Keller lift through every finite
Witt ring \(W_n(\mathbb F_2)=\mathbb Z/2^n\).

This is an unrestricted polynomial statement: no degree bound is imposed on
the correction terms.

## 2. First-order Jacobian correction

Use the displayed zero-one formulas as integral representatives \(P_0,Q_0\).
Every lift modulo four is uniquely of the form

\[
 \widetilde P=P_0+2A,\qquad \widetilde Q=Q_0+2B
\]

for polynomials \(A,B\) modulo two.  Direct differentiation gives

\[
 \det D(P_0,Q_0)=1+2K,
\]

where, modulo two,

\[
\begin{aligned}
K={}&xy+x^7y+x^{10}y+x^5y^2+x^{11}y^2\\
   &+x^9y^3+x^{12}y^3+x^{13}y^4.                 \tag{2.1}
\end{aligned}
\]

In particular

\[
 [xy]K=1.                                         \tag{2.2}
\]

Modulo four, the determinant of an arbitrary lift is

\[
 \det D(\widetilde P,\widetilde Q)
 =1+2\bigl(K+L(A,B)\bigr),                        \tag{2.3}
\]

with

\[
 L(A,B)=A_xQ_y+P_xB_y+A_yQ_x+P_yB_x              \tag{2.4}
\]

over \(\mathbb F_2\).  The signs become plus signs in characteristic two.

## 3. The one-coefficient obstruction

The four derivatives of the reduced map are

\[
 P_x=1,\quad P_y=x^2,\quad
 Q_x=x^4+x^6y^2,\quad Q_y=1+x^6+x^8y^2.          \tag{3.1}
\]

Now inspect the coefficient of \(xy\) in (2.4).

* In \(A_xQ_y\), only the constant term of \(Q_y\) could contribute.  It
  would require differentiating \(x^2y\), but its \(x\)-derivative is zero
  in characteristic two.
* In \(P_xB_y=B_y\), the only candidate is \(xy^2\), whose \(y\)-derivative
  is zero.
* Every monomial of \(A_yQ_x\) has \(x\)-degree at least four.
* Every monomial of \(P_yB_x=x^2B_x\) has \(x\)-degree at least two.

Therefore, for arbitrary polynomial corrections of any degree,

\[
 [xy]L(A,B)=0.                                    \tag{3.2}
\]

Equations (2.2)--(3.2) show that the coefficient of \(xy\) in every lifted
Jacobian is \(2\) modulo four.  It cannot be a constant.  This proves the
theorem.

## 4. One stabilization gives a full finite-Witt tower

The same calculation gives a sharp contrast.  Since

\[
 \det D(P_0,Q_0)=1+2K,
\]

put \(h=2K\) and, for every \(n\ge2\), define

\[
 S_n=\sum_{j=0}^{n-1}(-h)^j,
 \qquad
 \widetilde F_n(x,y,z)=\bigl(P_0,Q_0,zS_n\bigr) \pmod {2^n}.  \tag{4.1}
\]

This reduces to \((P,Q,z)\).  Its Jacobian matrix is block lower triangular.
The finite geometric-series identity gives

\[
 \det D\widetilde F_n
  =(1+h)S_n=1-(-h)^n=1\pmod {2^n}.               \tag{4.2}
\]

Moreover \(S_{n+1}\equiv S_n\pmod {2^n}\), so these polynomial maps form a
compatible tower.  Thus the obstruction class is nonzero in the plane
correction complex but becomes exact after one identity stabilization.  This
is not merely a failure of the particular \([xy]\) detector: (4.1) cancels
the entire Jacobian error at every finite Witt level.

The degrees of \(S_n\) grow with \(n\).  Their inverse limit is the restricted
two-adic power series

\[
 (1+2K)^{-1}=\sum_{j\ge0}(-2K)^j,                              \tag{4.3}
\]

not a polynomial in \(\mathbb Z_2[x,y]\).  Accordingly, (4.1) is a compatible
formal/Witt lift; it is not a finite-degree characteristic-zero polynomial
Keller map.

## 5. Consequences and boundary of the result

The obstruction closes the direct mixed-characteristic route for this exact
plane map.  It is stronger than the earlier observation that the displayed
integer formulas themselves are not Keller: allowing arbitrary higher-degree
corrections divisible by two does not help.

It does **not** rule out lifting some source/target equivalent
characteristic-two plane representative, nor does it obstruct unrelated
positive-characteristic Keller maps.  Stabilization at all finite Witt levels
is decided by (4.1), but polynomial algebraization over \(\mathbb Z_2\) is
not.  The most useful continuations are to determine whether a compatible
stable tower can have uniformly bounded degree, and to package the plane
obstruction as an equivalence-aware first Bockstein class under polynomial
left-right changes.

## 6. Exact reproduction

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_w2_obstruction.py
```

The checker computes the integral Jacobian error, verifies (3.1), and checks
the \(xy\)-coefficient functional on each of the four all-degree summands in
(2.4).  It verifies the universal geometric-series induction step and replays
the first finite levels of the explicit stable Witt tower.  The all-level
claim is the exact identity (4.2), not an inference from that regression.  The
checker performs no bounded correction search.
