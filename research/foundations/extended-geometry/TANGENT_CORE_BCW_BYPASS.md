# Tangent-core BCW bypass: the exact structural seed

## Outcome

The tangent core does bypass the expanded degree-seven representation, but
not by simply declaring its five incidence relations to be new output
coordinates.  For the foundational seed, the relations for
`gamma,W,C,s,t` have degrees

\[
 (3,3,2,2,3),
\]

so the full suspension square is already a quadratic--cubic polynomial
presentation in eight variables.  However, the Jacobian of these five
relations has rank only three on each of the three lifted collision points.
No three further coordinate functions can complete five rank-three rows to
an invertible eight-by-eight Jacobian.  Thus the naive eight-variable
residual-coordinate construction is not a Keller map and cannot certify a
new route bound.

The exact audit is
[`audit_tangent_core_bcw_bypass.py`](../scripts/audit_tangent_core_bcw_bypass.py).

## 1. Cubic structural presentation

Use the weighted normalization

\[
 \begin{aligned}
 \gamma_0&=1-\frac32xy-\frac12x^2z,\\
 U&=1+xy,\\
 W_0&=U\gamma_0,\\
 C_0&=x\gamma_0,
 \end{aligned}
\]

and put `H(W)=W^2(1-W)`.  The tangent core is

\[
 s_0=H'(W)+\gamma=2W-3W^2+\gamma,
 \qquad
 t_0=Ws_0-H(W)=W\gamma+W^2-2W^3.
\]

For the weighted foundational target `(A,B,C)`, the suspension square is
equivalent to the five equations

\[
 \begin{aligned}
 r_\gamma&=\gamma-\gamma_0,\\
 r_W&=W-(1+xy)\gamma,\\
 r_C&=C-x\gamma,\\
 r_s&=BC-(2W-3W^2+\gamma),\\
 r_t&=AC^2-(W\gamma+W^2-2W^3).
 \end{aligned}                                      \tag{1}
\]

Every expression in (1) has total degree at most three.  Direct substitution
recovers

\[
 C=C_0,\qquad BC=s_0(W_0,\gamma_0),\qquad
 AC^2=t_0(W_0,\gamma_0),
\]

without expanding either degree-seven target coordinate.

The three source collision points lift to

\[
 (W,\gamma,C)=(1,1,0),(0,0,0),(0,0,0)
\]

over the common weighted target `(-1/4,0,0)`.  The structural core therefore
retains the collision as an incidence, although it intentionally identifies
the last two points in the `(W,gamma,C)` projection.

## 2. Why the five residuals cannot be outputs

Let `R=(r_gamma,r_W,r_C,r_s,r_t)` and differentiate with respect to

\[
 (x,y,z,A,B,C,W,\gamma).
\]

At every lifted collision point,

\[
 \operatorname{rank} DR=3.                         \tag{2}
\]

An eight-dimensional map retaining all five entries of `R` as coordinate
functions can add only three more Jacobian rows.  Its rank on the collision
is therefore at most `3+3=6`, not eight.  This excludes arbitrary nonlinear
choices of the remaining three coordinates, not merely coordinate
projections.  The script also checks all 56 completions by three input
coordinates as a finite regression.

Geometrically, (2) is the same ramification that the weighted-suspension
determinant calculation cancels only after comparing the two vertical maps.
Turning both vertical equations into simultaneous residual outputs retains
the ramification instead of cancelling it.

There is a stronger row-by-row consequence.  The differential of `r_t`
vanishes at the two lifts with `(W,gamma,C)=(0,0,0)`.  Hence `r_t` cannot be
a coordinate of any terminal Keller map, regardless of the other seven
coordinates.  Moreover, no four residuals are simultaneously independent at
all three lifts.  Of the ten residual triples, only

\[
 (r_\gamma,r_W,r_C),\qquad (r_\gamma,r_C,r_s)       \tag{3}
\]

have row rank three at every collision point.  Thus every valid trace must
release at least two structural relations, must use the `t` relation only
transiently, and has only the two triples in (3) as maximal terminal
residual sets.

## 3. Correct search state

The viable search is a **coupled stable-shear search**.  Its state should
contain the exact polynomials in (1), but a relation is used transiently by a
triangular source or target shear and then released; all five residual rows
must not survive together in the terminal map.

A state needs to record:

- the eight current coordinate polynomials and their exact inverse shear
  trace;
- live structural nodes `gamma,W,C,s,t` and which output currently owns each
  node;
- the residual form of every unused side of (1);
- the maximum degree, quadratic part, cubic-output rank, and constant kernel;
- the transported three collision points.

Admissible transitions are triangular source shears, triangular target
shears, and linear coordinate exchanges.  A terminal state must pass, in
this order:

1. degree at most three and identity linear part;
2. exact determinant one (not only the suspension determinant on the
   incidence);
3. three distinct lifted points with a literal common image;
4. exact cubic-output rank and constant-kernel quotient;
5. either the two-parameter determinant identity or rank-compressed
   homogenization.

The correct optimization target remains the final homogeneous quotient
dimension.  The structural count five is therefore a search seed, not a
claimed `3+5` degree-reduced witness.  Any improvement of the 22-variable
route must come from an explicit stable shear trace satisfying the five
terminal checks above.

The collision-rank pruning should be applied before symbolic expansion:
discard every state which freezes `r_t`, every state which freezes four
residuals, and every three-residual state other than the two sets in (3).
This reduces the structural branch search to the two meaningful terminal
charts without imposing a monomial beam width.
