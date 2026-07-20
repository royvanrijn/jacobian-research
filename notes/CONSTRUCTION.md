# Construction and anatomy

## 1. Weighted invariants

The map respects

\[
(x,y,z)\mapsto(\lambda x,y/\lambda,z/\lambda^2).
\]

Use invariant coordinates `v=xy`, `s=x^2z`, and a scale coordinate `x`. For the
announced map define

\[
u=1+v,\quad \gamma=1-3v/2-s/2,\quad w=u\gamma,
\]

with seed `p(w)=2w-3w^2` and `q(w)=w^2-2w^3`. Then the map can be written

\[
F=(\alpha/x^2,\;2\beta/x,\;2x\gamma),
\quad \alpha=u+q(w)/\gamma^2,\quad
\beta=1+p(w)/\gamma.
\]

The apparent poles cancel after substituting `w=u gamma`. This presentation is
the construction mechanism: rational weighted coordinates are composed so that
their Jacobian factors and poles cancel in the final affine coordinates.

## 2. Fiber reduction

For target `(A,B,C)`, set

\[
P=BC/4,\qquad Q=AC^2/4.
\]

The inverse problem reduces to

\[
w^2-w^3=Pw-Q.
\]

For a root `w`, reconstruct

\[
\gamma=P-p(w),\quad u=w/\gamma,\quad v=u-1,
\quad x=C/(2\gamma),\quad y=v/x,
\]

and then `s=2(1-3v/2-gamma)`, `z=s/x^2`. Generically the cubic has three
distinct roots, hence three sheets. At its discriminant locus roots do not merge
at finite source points: `gamma` tends to zero and reconstruction sends `x` to
infinity.

For `C=0`, the source splits into `x=0` and
`x^2 z=2-3xy`. The former maps bijectively to the target plane; the latter gives
two more sheets and diverges as `x -> 0`.

## 3. General weighted-lift recipe

Choose a polynomial `p`, a nonzero constant `c`, satisfying

\[
p(0)=0,\qquad p(1)=-c,\qquad \int_0^1p(w)\,dw=0.
\]

Define `q(0)=0` by `q'(w)=w p'(w)/c`. Put
`kappa=p'(1)/c`, require `kappa != -2`, set
`a=-(1+kappa)/(2+kappa)`, and choose `b != 0`. With

\[
v=xy,\ s=x^2z,\ u=1+v,\ \gamma=1+av+bs,\ w=u\gamma,
\]

form

\[
F_p=(\alpha/x^2,\beta/x,x\gamma),\quad
\alpha=u+q(w)/\gamma^2,\quad \beta=c+p(w)/\gamma.
\]

The endpoint and integral constraints make both quotients polynomial. The
Jacobian is the constant `bc`, and the inverse equation is

\[
\int_0^w p(s)\,ds=wP-cQ.
\]

A degree-`d` seed therefore produces generic degree `d+1`. Degree one cannot
satisfy the constraints nontrivially; the quadratic seed is unique up to scale,
so the announced cubic-sheet example is minimal *within this family*.

The script `weighted_lift.py` checks these claims for the quadratic seed and can
accept further symbolic seeds.

