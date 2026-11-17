# All-degree polynomial-graph obstruction for the Meng–Yang family

13 November 2026. The bounded graph frontier is replaced by an all-degree
result, **HC4MYGALL1**, using the reproduced external structural theorem
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md). The proof below does not
use a plane degree bound or assume JC2. The broader mixed/coisotropic
problem remains open. No novelty priority, formal verification or external
review is claimed. The elementary no-plane and formal-algebraic arguments
are registered separately as **HC4MYGAP1**; they do not depend on Tang.

## The all-degree theorem

Over any characteristic-zero field let LN be nonzero, M arbitrary, and put

\[
 u=1+xy,\quad A=u^3p+3xu^2q-x^3r,
\]
\[
 B=y^2u(4+3xy)p+(y+3xy^2(4+3xy))q+(2x-3x^2y)r,
 \qquad\Psi=L A^2+M A+N B.
\]

For **every** polynomial R(x,y,p,q), set f_R=Psi(x,y,p,q,R).

**Theorem.** The polynomial f_R admits no two-dimensional constant vector
plane along which it is jointly affine. Consequently, if det Hess(f_R)
is a nonzero constant, grad(f_R) is a polynomial automorphism.

In particular, for the v2 member (L,M,N)=(1,13,2), no polynomial graph
containing both marked points can have nonzero constant Hessian determinant,
**in any degree**. More generally no polynomial graph in this chart can
produce any gradient counterexample. This does not say that every unmarked
polynomial graph has nonconstant determinant.

### 1. A direction invisible to x is forced to be the q direction

Write D_v for constant directional differentiation. On x=0 the graph
disappears from the potential:

\[
 f_R(0,y,p,q)=L p^2+M p+N(4y^2p+yq).
\]

For v=(0,b,c,d), restriction commutes with D_v, and

\[
 (D_v^2 f_R)|_{x=0}
   =8N b^2p+16Nbc y+2L c^2+2Nbd.
\]

If D_v^2 f_R=0 identically, comparison of the p coefficient gives b=0,
and the constant coefficient then gives c=0. Thus any such direction
with zero x component must be a multiple of e_q.

### 2. The q direction is never affine for a polynomial graph

The complete defining formulas give the identity

\[
 x^2 B+(2-3xy)A=u(u+1)p+2x(2u+1)q.                 \tag{1}
\]

Work over K(x,y,p) and regard A_R as a polynomial in q. Its derivative is

\[
 (A_R)_q=x(3u^2-x^2 R_q)\ne0,
\]

because 3u^2-x^2 R_q restricts to 3 at x=0. Hence n=deg_q A_R>=1.
By (1),

\[
 f_R=L A_R^2+
 \left(M-\frac{N(2-3xy)}{x^2}\right)A_R
 +\frac{N u(u+1)}{x^2}p+\frac{2N(2u+1)}xq.
\]

The q degree of this expression is exactly 2n: the other summands have
q degrees at most n and one. Characteristic zero therefore gives
(f_R)_qq nonzero. Equivalently, its leading q coefficient after two
derivatives is L*a_n^2*(2n)*(2n-1), which cannot vanish.
This is a polynomial-degree argument for all n, not a bounded search.

### 3. Excluding every cotangent plane

Every two-dimensional constant vector plane has a nonzero vector with
zero x component. If f_R were jointly affine along that plane, its second
directional derivative in that vector would vanish. Step 1 makes it a
nonzero multiple of e_q, contradicting Step 2. Thus no such plane exists.

A nonautomorphic constant-Hessian gradient would, by HCTANG2, have a
constant-linear cotangent form. Its two fibre directions span just such
a jointly affine plane. This proves the theorem over C. For an arbitrary
characteristic-zero field, embed a finitely generated coefficient subfield
in C and use the unique formal inverse to descend polynomial invertibility.
The elementary no-plane argument itself already holds over every such field.

### 4. Transferring the known collision

For v2 the two marked points are
(1,-3/2,0,0,0) and (-1,3/2,0,0,0), with common gradient
g=(0,0,-1/2,0,0). If both belong to the graph, the chain rule gives
grad f_R=Di_R^T grad Psi at them. The last component of g is zero,
so the two pulled-back gradients coincide independently of all first
derivatives of R. This contradicts the theorem if the Hessian determinant
is a nonzero constant. It is precisely the marked-graph collision route
inside OP-HC4-SCHUR that is now closed, not every mixed reduction.

## An explicit nonrational formal graph with constant determinant

The polynomial restriction in the proof has content. Let kappa also be
nonzero and choose the cotangent polynomial

\[
 F=L p^2+M p+4Ny^2p+Nyq+\kappa xp,
 \qquad\det\operatorname{Hess}F=(N\kappa)^2.
\]

Let A0=A|r=0, Psi0=Psi|r=0 and define the polynomial

\[
 C_R=N(2-3xy)-Mx^2-2Lx^2 A0,\qquad
 D_R=(\Psi0-F)/x.
\]

D_R is polynomial because F and Psi0 agree on x=0. The equation
Psi(x,y,p,q,R)=F is exactly

\[
 Lx^5R^2+C_R R+D_R=0.                              \tag{2}
\]

At x=0 its derivative in R is 2N, a unit. It therefore has a unique
formal solution in K[y,p,q][[x]], with

\[
 R_0=\frac{\kappa p-3(2Lp+M)(yp+q)-N(7y^3p+12y^2q)}{2N}.
\]

This formal graph has exactly the nonzero constant Hessian determinant
of F. It is also the unique constant-determinant formal graph for this
trace and constant, by the retained all-normal theorem HC4MYGJ2.
At v2 with kappa=1, R0 has total degree four and R1 has total degree
five, so the complete first normal jet already has total degree six.
It is therefore outside the earlier degree-at-most-four first-jet
obstruction HC4MYGJ1, rather than contradicting it.

Its quadratic discriminant is

\[
 \Delta=C_R^2-4Lx^5D_R,
 \quad [q]\Delta=-12LNx^3(xy+2),\quad \deg_q\Delta=1.
\]

Thus Delta is not a square in K(x,y,p)(q): it has order one at its
linear prime. Equation (2) has no rational solution in K(x,y,p,q),
although its formal solution exists. In particular the formal solution
is not a polynomial graph.

There is an explicit all-order growth certificate. Choose the square-root
branch with sqrt(Delta)|x=0=2N, so

\[
 R=\frac{-C_R+\sqrt\Delta}{2Lx^5}.
\]

For every j>=2 the first q^j term occurs at x order 3j-5 and has coefficient

\[
 \frac NL\binom{1/2}{j}\left(-\frac{6L}{N}\right)^j\ne0.
\]

This follows by writing Delta=Delta0-12LNx^3(xy+2)q with
Delta0|x=0=4N^2 and expanding its square root in q. Neither -C_R nor A0
has q degree greater than one. The formula proves nontermination in
arbitrarily large q degree; checking a few orders is only a regression.

For v2 and kappa=1, the algebraic graph component (Psi-F)/x=0 even contains both
marked points, with a double root R=0 in each fibre. However Delta vanishes
at them and its q derivatives are -12 and 12 respectively. Thus its
projection to the four graph coordinates is ramified there. No branch
is a regular graph near either marked point: differentiating (2) there
would give 0*R_q+3=0. The formal solution near x=0 cannot be used to
transfer a global marked collision through those points.

This supplies an explicit algebraic local-to-global boundary: exact
constant determinant and arbitrarily many correct normal jets coexist
with nonrationality and failure of a regular marked graph. It does not
construct a Keller counterexample on affine four-space.

## Evidence and remaining scope

The [checker](verify.py) verifies (1), the universal slice derivative,
the complete quadratic equation/discriminant, the constant determinant
and marked branch obstruction. It also generates nine normal coefficients
at the v2 parameters and checks the first nonzero higher-q coefficients.
The all-degree q argument and binomial formula are written proofs.

```sh
timeout 60s .venv/bin/python research/hessian/hc4-polynomial-graphs/verify.py
```

The imported HCTANG2 structural theorem remains a written reproduction
with the scope and assurance stated in its canonical review. This result
uses neither a stronger plane-degree frontier nor a new proof of JC2.
Other eliminated-coordinate graph charts, non-graph embeddings, rational
Schur pivots and general coisotropic reductions require their own argument.
The graph coordinates are (x,y,p,q), up to constant affine changes.
Arbitrary nonlinear reparameterizations of the same geometric graph are
also outside this theorem: the gradient and Hessian conditions change
under such a reparameterization.
