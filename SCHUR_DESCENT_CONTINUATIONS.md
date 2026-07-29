# Exact Schur remainders and simultaneous descent

## Status

This note records two elementary strengthenings of the
[Meng--Yang Schur-descent bridge](MENG_YANG_SCHUR_DESCENT_BRIDGE.md) and
uses them to sharpen the remaining `HC(4)` search.

The proved additions are:

1. an exact one-pivot remainder formula, which weakens the published
   identically-singular-pencil hypothesis for any fixed descent;
2. a simultaneous \(r\)-pivot theorem, which gives a direct
   \(2n\)-to-\((2n-r)\) descent for Keller maps that are affine-linear in
   an \(r\)-dimensional source block;
3. an explicit collision lift showing that every value of the parameter
   \(\mu\), not only the \(y=0\) normalization, gives a collision in the
   doubled-Keller family;
4. a codimension-one block-affine theorem showing that the apparent direct
   \(6\)-to-\(4\) source-block route cannot start from a Keller
   counterexample;
5. a rigidity observation closing pure univariate higher-degree repairs.

These statements do not settle `HC(4)`.  They replace several broad search
phrases by narrower algebraic targets.

Throughout, \(K\) is a field of characteristic zero.

## 1. The exact one-pivot remainder

Use the notation

\[
 \Phi(t,w)=tA(w)+B(w),\qquad
 g=\nabla A,\qquad
 M(s,w)=\operatorname{Hess}_w(B+sA),
\]

and assume only

\[
 \det\operatorname{Hess}_{t,w}\Phi=c\in K^\times.
\]

Set

\[
 D(s,w)=\det M(s,w)
\]

and, for \(\lambda\ne0\),

\[
 \psi_{\lambda,\mu}
 =B+\frac{\lambda}{2}A^2+\mu A.
\]

> **Proposition 1.1 (exact Schur remainder).**
>
> \[
>  \boxed{
>  \det\operatorname{Hess}\psi_{\lambda,\mu}
>  =
>  D(\mu+\lambda A(w),w)-\lambda c.
>  }
>  \tag{1.1}
> \]

Indeed, the bordered determinant identity gives

\[
 g^T\operatorname{adj}(M(s,w))g=-c
\]

identically in \(s,w\), while

\[
 \operatorname{Hess}\psi_{\lambda,\mu}
 =M(\mu+\lambda A,w)+\lambda gg^T.
\]

The rank-one determinant identity now gives (1.1).

Meng--Yang's hypothesis \(D(s,w)\equiv0\) is therefore a clean structural
sufficient condition, but it is not necessary for a fixed descent.  The
exact necessary-and-sufficient condition is

\[
 D(\mu+\lambda A(w),w)\in K.                          \tag{1.2}
\]

If the constant in (1.2) is \(d\), the descended determinant is
\(d-\lambda c\), which must also be nonzero.

Two elementary examples show that both qualifications are real.  Take
\(w=(x,y)\) and \(A=x\).  For

\[
 \Phi=tx+x^3+\frac12y^2
\]

the bordered Hessian determinant is \(c=-1\), but \(D=6x\), so
\(\det\operatorname{Hess}\psi_{\lambda,\mu}=6x+\lambda\) is not constant.
Thus a constant bordered determinant alone does not suffice.  On the other
hand, for

\[
 \Phi=tx+\frac12x^2+\frac12y^2
\]

one has \(c=-1\) and \(D=1\), giving
\(\det\operatorname{Hess}\psi_{\lambda,\mu}=1+\lambda\).  Hence
\(D\equiv0\) is not necessary.

This is the first simplification to use in new searches: compute the
specialized reduced determinant in (1.2), not every homogeneous face of the
full descended Hessian.  A nonzero or nonconstant pencil \(D(s,w)\) may
still become constant on the graph \(s=\mu+\lambda A(w)\).

## 2. Every \(\mu\) carries a doubled-Keller collision

Let

\[
 F(x',t)=F_0(x')+tF_1(x')
\]

be a Keller map and suppose

\[
 p_\pm=(x'_\pm,t_0),\qquad
 F(p_+)=F(p_-),\qquad p_+\ne p_-.
\]

For the doubled potential \(\Phi=\langle y,F\rangle=tA+B\), fix arbitrary
\(\lambda\in K^\times\) and \(\mu\in K\).  Put

\[
 A_0=\frac{t_0-\mu}{\lambda}
\]

and choose any covector \(\eta\in K^n\) whose \(t\)-component is \(A_0\).
Because \(DF(p_\pm)\) is invertible, define

\[
 y_\pm=DF(p_\pm)^{-T}\eta.                            \tag{2.1}
\]

Then

\[
 \nabla\Phi(p_\pm,y_\pm)
 =\bigl(DF(p_\pm)^Ty_\pm,F(p_\pm)\bigr)
 =\bigl(\eta,F(p_\pm)\bigr),
\]

so the two doubled gradients agree.  Moreover,

\[
 A(x'_\pm,y_\pm)
 =\langle y_\pm,F_1(x'_\pm)\rangle
 =\eta_t=A_0.
\]

The collision-transfer normalization
\(\mu=t_0-\lambda A_0\) is therefore automatic by construction.

> **Corollary 2.1.**  For a doubled Keller collision over a common affine
> pivot coordinate, every member \(\psi_{\lambda,\mu}\), with
> \(\lambda\ne0\), has a gradient collision.  The choice \(\mu=t_0\)
> makes \(\eta=0\) and \(y_\pm=0\); other values of \(\mu\) move the
> collision to nonzero dual coordinates.

This explains uniformly why the Meng--Yang representatives
\(A^2+11A+2B\) and \(A^2+13A+2B\) are both counterexamples: the latter is
preferred only because its collision lies on \(y=0\).

## 3. Simultaneous \(r\)-variable descent

Let \(t=(t_1,\ldots,t_r)^T\), let \(w\) have length \(m\), and write

\[
 \Phi(t,w)=t^TA(w)+B(w),
\]

where \(A=(A_1,\ldots,A_r)^T\).  Put

\[
 J=D_wA,\qquad
 M(s,w)=\operatorname{Hess}_w(B+s^TA).
\]

Thus

\[
 \operatorname{Hess}_{t,w}\Phi
 =
 \begin{pmatrix}
  0&J\\
  J^T&M(t,w)
 \end{pmatrix}.
\]

Let \(\Lambda\) be an invertible symmetric \(r\times r\) matrix and
\(\mu\in K^r\).  Define

\[
 \psi_{\Lambda,\mu}
 =B+\mu^TA+\frac12A^T\Lambda A.                       \tag{3.1}
\]

There is an exact formula before imposing a corank condition.  For a
symmetric \(r\times r\) matrix \(Q\), put

\[
 \mathcal B(Q;s,w)
 =
 \det
 \begin{pmatrix}
  Q&J\\
  J^T&M(s,w)
 \end{pmatrix}.
\]

Schur complementation of the repaired pivot gives

\[
 \boxed{
 \det\operatorname{Hess}\psi_{\Lambda,\mu}
 =
 (-1)^r\det(\Lambda)\,
 \mathcal B\!\left(
  -\Lambda^{-1};\,\mu+\Lambda A(w),w
 \right).
 }
 \tag{3.2}
\]

The original constant-Hessian hypothesis says only
\(\mathcal B(0;s,w)=c\).  Thus the exact general gate is constancy and
nonvanishing of the right-hand side of (3.2).  The following rank condition
makes that constancy automatic by making \(\mathcal B\) independent of
its corner \(Q\).

> **Theorem 3.1 (simultaneous Schur descent).**  Assume
>
> \[
>  \det\operatorname{Hess}_{t,w}\Phi=c\in K^\times,
>  \qquad
>  \operatorname{rank}M(s,w)\le m-r
>  \tag{3.3}
> \]
>
> over \(K(s,w)\).  Then
>
> \[
>  \boxed{
>  \det\operatorname{Hess}\psi_{\Lambda,\mu}
>  =(-1)^r\det(\Lambda)c.
>  }
>  \tag{3.4}
> \]

### Proof

Repair all pivots at once:

\[
 \widehat\Phi
 =\Phi-\frac12(t-\mu)^T\Lambda^{-1}(t-\mu).
\]

Its critical equation has the polynomial solution

\[
 t^*=\mu+\Lambda A,
\]

and its critical value is (3.1).  Schur complementation of the pivot
\(-\Lambda^{-1}\) gives

\[
 \det\operatorname{Hess}\widehat\Phi
 =(-1)^r\det(\Lambda)^{-1}
  \det\operatorname{Hess}\psi_{\Lambda,\mu}.          \tag{3.5}
\]

It remains to see that the quadratic repair does not change the determinant.
Work over \(K(s,w)\).  Since the full bordered matrix with zero corner is
invertible, adjoining its \(r\) border rows and columns to \(M\) shows
\(\operatorname{rank}M\ge m-r\).  Hypothesis (3.3) therefore gives equality.
After invertible row and column operations on the \(w\)-block, eliminate
the nonsingular \((m-r)\)-dimensional part of \(M\).  The remaining block has
the form

\[
 \begin{pmatrix}
  Q&U\\
  V&0
 \end{pmatrix},
\qquad U,V\text{ of size }r\times r,
\]

whose determinant is independent of the corner \(Q\).  Replacing the zero
corner of \(\operatorname{Hess}\Phi\) by \(-\Lambda^{-1}\) consequently
leaves its determinant equal to \(c\).  Equation (3.5) gives (3.4).

For \(r=1\), the rank condition is exactly \(\det M=0\), so this recovers
Meng--Yang.

## 4. The block-affine Keller corollary

Let a Keller map be jointly affine-linear in an \(r\)-dimensional source
block:

\[
 F(x,t)=F_0(x)+\sum_{i=1}^r t_iF_i(x),
\qquad x\in K^{n-r}.                                  \tag{4.1}
\]

For \(\Phi=\langle y,F\rangle=t^TA+B\), the remaining variables are
\(w=(x,y)\), so \(m=2n-r\).  The reduced Hessian has block form

\[
 M(s,w)=
 \begin{pmatrix}
  *&C(s,x)^T\\
  C(s,x)&0
 \end{pmatrix},
\]

where \(C\) has \(n-r\) columns.  Hence

\[
 \operatorname{rank}M\le2(n-r)=m-r,
\]

which is exactly (3.3).  Since the doubling determinant is

\[
 c=(-1)^n(\operatorname{Jac}F)^2,
\]

Theorem 3.1 gives

\[
 \det\operatorname{Hess}\psi_{\Lambda,\mu}
 =(-1)^{n+r}\det(\Lambda)(\operatorname{Jac}F)^2.     \tag{4.2}
\]

If \(F\) has a collision

\[
 F(x_+,t_0)=F(x_-,t_0),\qquad x_+\ne x_-,
\]

the vector version of (2.1), with

\[
 \eta_t=\Lambda^{-1}(t_0-\mu),
\]

transfers that collision to every \(\psi_{\Lambda,\mu}\).

> **Corollary 4.1.**  A three-variable Keller counterexample that can be
> put in the form
>
> \[
>  F(x,t_1,t_2)=F_0(x)+t_1F_1(x)+t_2F_2(x)
> \]
>
> while retaining a collision over one common pair \((t_1,t_2)\) produces
> a four-variable Hessian counterexample in one simultaneous descent.

Formally, this is a one-shot route to `HC(4)`.  The next theorem shows that
its premise is impossible for a Keller counterexample.

## 5. Two useful obstructions

### 5.1 Codimension-one block-affine Keller maps are automorphisms

> **Theorem 5.1.**  Let \(F:K^n\to K^n\) be a Keller map.  If, in some
> polynomial source coordinate system \((x,t_1,\ldots,t_{n-1})\),
>
> \[
>  F(x,t)=a(x)+B(x)t,
>  \tag{5.1}
> \]
>
> then \(F\) is a polynomial automorphism.

The elementary proof is recorded for use in this descent problem; no
literature-wide novelty or priority claim is made.

Write the columns of \(B\) as \(b_1,\ldots,b_{n-1}\), and let
\(N(x)\in K[x]^n\) be their cofactor normal, so

\[
 N^Tb_i=0,\qquad
 \det(v,b_1,\ldots,b_{n-1})=N^Tv.
\]

The Jacobian identity expands as

\[
 \det\bigl(a'+B't,b_1,\ldots,b_{n-1}\bigr)=c\in K^\times.
\]

Every coefficient of \(t_i\) vanishes, hence

\[
 N^Tb_i'=0.
\]

Differentiating \(N^Tb_i=0\) gives \(N'^Tb_i=0\).  Over \(K(x)\), the
columns of \(B\) have rank \(n-1\), so their left kernel is one-dimensional;
therefore \(N'\) is proportional to \(N\).  Equivalently, every projective
ratio of two nonzero components of \(N\) has derivative zero.  The constants
of \(d/dx\) on \(K(x)\) are \(K\), so

\[
 N(x)=q(x)N_0
\]

for a fixed nonzero vector \(N_0\in K^n\) and \(q\in K[x]\).  The constant
term of the Jacobian identity is

\[
 c=N^Ta'=q(x)N_0^Ta'(x).
\]

Since the product of two polynomials is a unit, \(q\in K^\times\).  Thus
the normal \(N\) is constant.

Apply a constant target-linear change whose last row is \(N^T\).  The map
then has the triangular form

\[
 \widetilde F(x,t)
 =
 \bigl(g(x)+C(x)t,\ \alpha x+\beta\bigr),
\qquad
\alpha\in K^\times,\quad
\det C\in K^\times.
\]

Its inverse is polynomial:

\[
 x=\frac{z_n-\beta}{\alpha},\qquad
 t=C(x)^{-1}\bigl(z_{1,\ldots,n-1}-g(x)\bigr).
\]

This proves the theorem.

For \(n=3\), Theorem 5.1 closes Corollary 4.1 completely: no polynomial
source rechart of a Keller counterexample can produce the jointly affine
two-dimensional source block.  The family-specific constant-direction
calculation in the
[double-Schur gauge audit](HC4_DOUBLE_SCHUR_GAUGE_AUDIT.md) remains a useful
local diagnostic, but the obstruction is in fact global and allows nonlinear
source recharts.

### 5.2 Pure higher-degree pivot repairs do not enlarge the search

Suppose the one-variable repair is

\[
 \widehat\Phi(t,w)=\Phi(t,w)-R(t)
\]

and its critical equation \(R'(t)=A(w)\) has a solution
\(t=T(A(w))\) for a polynomial \(T\), with \(A\) nonconstant.  Since
substitution \(K[z]\to K[w]\), \(z\mapsto A(w)\), is injective,

\[
 R'(T(z))=z.
\]

Polynomial degrees force

\[
 \deg R'=\deg T=1.
\]

Therefore \(R\) is quadratic.  A pure univariate repair of degree at least
three cannot provide a new polynomial critical solution.  What remains open
under the phrase "higher-degree critical equation" must involve
\(w\)-dependent coefficients, a mixed source--dual transformation, or some
other special divisibility mechanism.

## 6. Prioritized continuation

The exact formulas suggest the following order.

1. **Relax the singular-pencil gate.**  In every nonlinear coordinate or
   coisotropic chart with an affine pivot, test the composite
   \(D(\mu+\lambda A,w)\) from (1.2).  Requiring \(D(s,w)\equiv0\) is
   unnecessarily strong.
2. **Move to mixed source--dual or coisotropic pivots.**  Theorem 5.1 rules
   out every pure-source codimension-one block, even after a nonlinear
   polynomial source rechart.  A viable two-pivot construction must use
   genuinely mixed canonical variables or a non-coordinate embedding.
3. **Use the corank budget before coefficient expansion.**  For such a
   proposed
   simultaneous \(r\)-pivot chart, first test
   \(\operatorname{rank}M\le m-r\).  Success makes the determinant automatic.
   Failure does not exclude descent; it sends the problem to the exact
   bordered-graph expression (3.2).
4. **Separate determinant and collision work.**  For doubled Keller maps,
   (2.1) makes collision transfer linear algebra once the pivot coordinates
   are fixed.  The real search is the determinant/corank condition.
5. **Classify the \((\lambda,\mu)\)-family.**  The parameters give genuine
   counterexamples, but their polynomial-equivalence classes and the minimum
   degree obtainable after allowed gauges remain open.  This is secondary to
   `HC(4)` but may expose invariants useful in the descent search.

The current obstructions remain substantial: constant linear two-pivot
routes are excluded in the root-engineered gauge families, the known toric
correction restores constant determinant only with an invertible descended
gradient, and the coordinate four-variable homogeneous chart is closed
through the displayed cubic--quartic--sextic range.  The live directions are
therefore the relaxed remainder (1.2) in genuinely mixed
source--dual/coisotropic pivots, and higher homogeneous layers not covered by
the existing coordinate chart.
