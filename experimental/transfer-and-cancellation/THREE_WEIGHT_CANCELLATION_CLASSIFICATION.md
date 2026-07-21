# Classification after one additional weight

This note carries out the first relaxation proposed beyond the completed
two-weight cancellation theorem.  The exponent in the `Q` coordinate is
allowed to vary independently.  The constant-Jacobian equation forces it
back to the old value, so this third weight produces no branch beyond C24.

Throughout, `k` is a characteristic-zero field, `f in k[y]` is nonconstant,
`a,b>=1`, and `c>=0` are integers.  Let

\[
 A=1+xf(y),\qquad B=A^bz+g(y,A),                              \tag{1}
\]

\[
 P=A^aB,\qquad Q=y+xA^cB.                                    \tag{2}
\]

The former two-weight skeleton imposed `c=a-1`.  Here `c` is independent.
Put

\[
 h=c-a,\qquad s=xA^h,\qquad W(t)=t f(Q-Pt).                   \tag{3}
\]

Then `Q=y+sP`.  For a nonzero polynomial `Theta in k[W]`, consider the most
general one-factor polynomial resolvent derivative in this inverse variable,

\[
 R=C\int_0^s\Theta\bigl(t f(Q-Pt)\bigr)dt,
 \qquad C\in k^*.                                             \tag{4}
\]

All expressions are initially interpreted in `k[x,y,z,A^(-1)]`.  The
question is which weights and which `Theta` can make the localized Jacobian
a nonzero constant.  Polynomial cancellation is a subsequent, stronger
condition.

## 1. The third-weight factor

At fixed `y`, differentiation of `s=xA^h` gives

\[
 {\partial s\over\partial x}
 =A^{h-1}\{(h+1)A-h\}.                                      \tag{5}
\]

At fixed `(s,y)`, one has `P=A^aB` and `Q=y+sP`, whence

\[
 \det {\partial(P,Q)\over\partial(y,B)}=-A^a.                 \tag{6}
\]

Since `partial B/partial z=A^b`, equations (5)--(6) give

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}
 =-A^{b+c-1}\{(h+1)A-h\}.                                   \tag{7}
\]

Holding `(P,Q)` fixed, (4) gives `R_s=C Theta(W(s))`.  On the
source,

\[
 Q-Ps=y,\qquad W(s)=s f(y)=(A-1)A^h.                         \tag{8}
\]

Therefore

\[
 \boxed{\quad
 \det {\partial(P,Q,R)\over\partial(x,y,z)}
 =-C A^n\{(h+1)A-h\}
       \Theta\bigl((A-1)A^h\bigr),\quad n=b+c-1\geq0.
 \quad}                                                       \tag{9}
\]

## 2. Complete weight classification

### Theorem 2.1

The expression in (9) is a nonzero constant if and only if

\[
 c=a-1,
 \qquad
 \Theta(W)=\lambda(1-W)^{a+b-2}
 \quad(\lambda\in k^*).                                     \tag{10}
\]

Thus allowing the additional weight `c` gives no new localized Keller
branch.  It recovers exactly the derivative power and weights of the
two-weight skeleton.

**Proof.**  Suppose first that the right side of (9) is a nonzero constant as
a Laurent polynomial in `A`.

If `h` is neither `-1` nor `0`, the factor

\[
 L_h(A)=(h+1)A-h
\]

vanishes at the nonzero point `A=h/(h+1)`.  The other factors in (9) are
finite there, including `(A-1)A^h`, even when `h<0`.  The full product would
therefore vanish, a contradiction.

If `h=0`, then `L_h(A)=A` and `(A-1)A^h=A-1`.  Constancy would require

\[
 A^{n+1}\Theta(A-1)\in k^*.
\]

But `n>=0` and `Theta(A-1)` is a nonzero polynomial in `A`, so the left side
has positive degree.  This is impossible.

Hence `h=-1`, equivalently `c=a-1`.  Now `L_h=1`,
`W=(A-1)A^(-1)=1-A^(-1)`, and (9) is constant exactly when

\[
 A^n\Theta(1-A^{-1})\in k^*.
\]

The substitution `W=1-A^(-1)` is invertible with `A=(1-W)^(-1)`, so this is
equivalent to

\[
 \Theta(W)=\lambda(1-W)^n.
\]

Finally `n=b+c-1=a+b-2`, proving necessity.  Direct substitution in (9)
proves sufficiency.  QED

### Corollary 2.2 (no branch beyond C24)

Assume in addition that all three coordinates `(P,Q,R)` are polynomial.
Then, up to translations, nonzero scalings, and polynomial left--right
equivalence, the map is a C24 map.

**Proof.**  Theorem 2.1 forces `c=a-1` and the derivative
`Theta(W)=lambda(1-W)^(a+b-2)`.  After absorbing `lambda` into `C`, this is
exactly equations (1)--(2) of
[the generalized two-weight theorem](GENERALIZED_CANCELLATION_MECHANISM.md).
If `a=b=1`, then `e=a+b-2=0` and (4) is `R=Cs=Cx/A`, which is not a
polynomial because `A=1+xf(y)` does not divide `x`.  Hence a polynomial member
has `e>=1`.  The two-weight weight-rigidity theorem then forces `a=1`,
`b=e+1`, and its complete spectral classification shows that every polynomial
member is left--right equivalent to C24.  QED

The conclusion is deliberately scoped to one additional monomial weight and
a polynomial derivative depending on the single natural inverse factor
`W=t f(Q-Pt)`.  A second independent resolvent factor is not covered.

## 3. Exact regression

Run

```bash
.venv/bin/python scripts/verify_three_weight_cancellation.py
```

The script independently differentiates representative symbolic instances,
checks formula (9), exhausts a bounded box of integer weights and polynomial
`Theta`, and confirms that every constant-Jacobian solution has precisely the
form (10).  The all-weight theorem is the Laurent-polynomial argument above,
not the bounded search.
