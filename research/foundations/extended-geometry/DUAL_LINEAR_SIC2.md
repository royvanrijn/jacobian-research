# The dual-linear `SIC(2)` theorem

## 1. Statement and scope

Let \(k\) be a field of characteristic zero and put

\[
 A_2=k[w_1,w_2,x,y].
\]

For \(\alpha=(\alpha_1,\alpha_2)\), define

\[
 \mathcal E_2(w^\alpha q(x,y))
 =\partial_x^{\alpha_1}\partial_y^{\alpha_2}q(x,y).       \tag{1.1}
\]

This note treats every polynomial that is homogeneous of degree one in the
dual variables:

\[
 p=w\mathbin{\cdot}H=w_1P(x,y)+w_2Q(x,y).                  \tag{1.2}
\]

No Keller hypothesis, homogeneity in \(x,y\), or normalization at the
origin is imposed.

> **Theorem 1.1 — Dual-linear `SIC(2)`.**
> Let \(p\) have the form (1.2).  If
> \[
>  \mathcal E_2(p)=\mathcal E_2(p^2)=0,                    \tag{1.3}
> \]
> then there are \(a,b,c_1,c_2\in k\), with \((a,b)\ne(0,0)\),
> and \(f\in k[t]\) such that
> \[
> \begin{aligned}
>  H(x,y)&=(c_1,c_2)+(b,-a)f(ax+by),\\
>  p&=c_1w_1+c_2w_2+(bw_1-aw_2)f(ax+by).                  \tag{1.4}
> \end{aligned}
> \]
> Conversely, every polynomial in the normal form (1.4) has
> \(\mathcal E_2(p^m)=0\) for every \(m\geq1\).  For every
> \(g\in A_2\),
> \[
>  \mathcal E_2(gp^m)=0\qquad(m\gg0).                      \tag{1.5}
> \]
> More explicitly, if
> \[
>  d=\max(0,\deg f),\qquad G=\deg_{x,y}g,
> \]
> then
> \[
>  \boxed{m>(d+2)G\quad\Longrightarrow\quad
>  \mathcal E_2(gp^m)=0.}                                  \tag{1.6}
> \]

The usual SIC premise assumes
\(\mathcal E_2(p^m)=0\) for every \(m\geq1\), so (1.3) is automatic.
Therefore \(\operatorname{SIC}(2)\) holds on the complete dual-linear
stratum, and its infinite pure-moment premise collapses there to the first
two moments.

This theorem is all-degree and nonhomogeneous in \(x,y\).  It does not
cover polynomials of dual degree at least two.

## 2. The first two contractions

Write

\[
 H=(H_1,H_2)=(P,Q),\qquad \delta=\operatorname{div}H=P_x+Q_y.
\]

The first contraction is

\[
 \mathcal E_2(p)=\delta.                                   \tag{2.1}
\]

For the second, expand

\[
 \mathcal E_2(p^2)
 =\sum_{i,j=1}^2\partial_i\partial_j(H_iH_j).
\]

Applying the product rule and collecting the four types of terms gives

\[
 \boxed{
 \mathcal E_2(p^2)
 =2H\mathbin{\cdot}\nabla\delta
  +\delta^2+\operatorname{tr}\big((JH)^2\big).
 }                                                         \tag{2.2}
\]

Under (1.3), \(\delta=0\), and hence

\[
 0=\mathcal E_2(p^2)=\operatorname{tr}\big((JH)^2\big).
                                                               \tag{2.3}
\]

For a two-by-two matrix,

\[
 \operatorname{tr}(A^2)=(\operatorname{tr}A)^2-2\det A.
\]

Equations (2.1)--(2.3) therefore imply

\[
 \operatorname{tr}JH=0,\qquad \det JH=0.                   \tag{2.4}
\]

Thus the first two contractions force \(JH\) to be nilpotent.  In the
earlier Keller-constrained theorem, the determinant equation supplied the
second equality in (2.4); here the second pure contraction supplies it
without any Keller input.

## 3. The binary zero-Hessian lemma

The trace equation in (2.4) supplies a polynomial Hamiltonian
\(h\in k[x,y]\):

\[
 P=h_y,\qquad Q=-h_x.                                      \tag{3.1}
\]

Indeed, integrate \(P\) in \(y\).  The \(y\)-derivative of the resulting
\(h_x+Q\) is zero, so a final antiderivative in \(x\) gives (3.1).

Moreover,

\[
 \det JH=h_{xx}h_{yy}-h_{xy}^2.                            \tag{3.2}
\]

We need the following elementary classification.

> **Lemma 3.1.**
> If \(h\in k[x,y]\) and
> \[
>  h_{xx}h_{yy}-h_{xy}^2=0,
> \]
> then
> \[
>  h=\phi(ax+by)+\ell(x,y)                                 \tag{3.3}
> \]
> for some \(\phi\in k[t]\), some \((a,b)\ne(0,0)\), and an
> affine-linear polynomial \(\ell\).

### Proof

Put

\[
 A=h_{xx},\qquad B=h_{xy},\qquad C=h_{yy}.
\]

If \(A=0\), then \(B^2=AC=0\), so \(B=0\).  Thus \(h_x\) is constant and
(3.3) follows with the one-variable part depending on \(y\).

Suppose \(A\ne0\), work in \(k(x,y)\), and set \(r=B/A\).  Since
\(AC=B^2\), we have \(B=rA\) and \(C=r^2A\).  Equality of the third mixed
derivatives gives

\[
 A_y=B_x,\qquad B_y=C_x,
\]

and substitution yields

\[
 r_y=r\,r_x.                                               \tag{3.4}
\]

Write \(r=R/S\) with coprime \(R,S\in k[x,y]\).  Clearing denominators in
(3.4) gives

\[
 S(R_yS-RS_y)=R(R_xS-RS_x).                               \tag{3.5}
\]

Reduction modulo \(R\) shows \(R\mid R_y\), hence \(R_y=0\).  Reduction
modulo \(S\) then gives \(S_x=0\), and (3.5) reduces to

\[
 -S_y=R_x.
\]

Consequently either \(r\) is constant or, for some \(\gamma\ne0\),

\[
 R=\gamma x+\rho,\qquad S=\sigma-\gamma y.                 \tag{3.6}
\]

The second case is incompatible with a polynomial Hessian.  Polynomiality
and coprimality give

\[
 A=S^2T,\qquad B=RST,\qquad C=R^2T
\]

for some \(T\in k[x,y]\).  The identity \(A_y=B_x\) becomes

\[
 S T_y-R T_x=3\gamma T.                                    \tag{3.7}
\]

In the independent affine coordinates \(R,S\), the left side is

\[
 -\gamma(R\partial_R+S\partial_S)T.
\]

Equation (3.7) would require every nonzero homogeneous part of the
polynomial \(T\) to have degree \(-3\).  Hence \(T=0\), contradicting
\(A\ne0\).

Thus \(r=\lambda\in k\).  Then

\[
 (A,B,C)=A(1,\lambda,\lambda^2).
\]

The derivation \(D=\partial_y-\lambda\partial_x\) annihilates both \(h_x\)
and \(h_y\), so \(Dh\) is constant.  After subtracting an affine-linear
polynomial, \(Dh=0\), and the polynomial kernel of \(D\) is
\(k[x+\lambda y]\).  This proves (3.3). \(\square\)

## 4. The dual-linear normal form

Apply Lemma 3.1 to (3.2).  Write

\[
 h=\phi(ax+by)+\alpha x+\beta y+\eta.
\]

Equation (3.1) gives

\[
 H=(b,-a)\phi'(ax+by)+(\beta,-\alpha).                    \tag{4.1}
\]

Taking \(f=\phi'\) and \(c=(\beta,-\alpha)\) proves (1.4).
Notice that no origin normalization is needed: it is exactly the constant
vector \(c\) that survives in the general normal form.

## 5. Pure and mixed contractions

Put

\[
 u=bw_1-aw_2,\qquad s=ax+by,\qquad r=c_1w_1+c_2w_2,
\]

and let

\[
 D=b\partial_x-a\partial_y,\qquad
 C=c_1\partial_x+c_2\partial_y.
\]

Then

\[
 p=r+uf(s),\qquad D(s)=0.                                  \tag{5.1}
\]

Expand a multiplier as

\[
 g=\sum_\alpha w^\alpha g_\alpha(x,y).
\]

The binomial theorem and the definition of \(\mathcal E_2\) give

\[
 \mathcal E_2(gp^m)
 =\sum_{\alpha}\sum_{k=0}^{m}\binom{m}{k}
 \partial^\alpha C^{m-k}D^k
 \left(g_\alpha f(s)^k\right).                             \tag{5.2}
\]

All constant-coefficient derivatives commute, and \(Df(s)=0\).  Therefore

\[
 D^k\left(g_\alpha f(s)^k\right)
 =(D^kg_\alpha)f(s)^k.                                    \tag{5.3}
\]

Let \(G=\deg_{x,y}g\).  If \(k>G\), equation (5.3) is zero.  If
\(k\leq G\), the polynomial inside \(C^{m-k}\) has degree at most

\[
 G+kd\leq(d+1)G.
\]

When \(m>(d+2)G\), we have

\[
 m-k\geq m-G>(d+1)G,
\]

so the \(C^{m-k}\) derivative also kills the term.  Every summand of
(5.2) is zero, proving (1.6).  Taking \(g=1\) proves all pure
contractions vanish.

## 6. Keller-constrained corollary

Suppose additionally that

\[
 H(0)=0,\qquad JH(0)=0,\qquad \det J((x,y)-H)=1.            \tag{6.1}
\]

If only \(\mathcal E_2(p)=0\), then

\[
 1=\det(I-JH)
  =1-\operatorname{tr}JH+\det JH
\]

forces \(\det JH=0\).  Thus the proof above applies without separately
assuming the second contraction.  The normalizations in (6.1) remove the
constant vector and the constant and linear parts of \(f\), giving

\[
 H=(b,-a)f(ax+by),\qquad f\in t^2k[t].                    \tag{6.2}
\]

Here \(p=uf(s)\) has no constant-dual summand \(r\), so (5.2) contains only
\(k=m\).  The sharper cutoff is

\[
 m>\deg_{x,y}g.                                            \tag{6.3}
\]

This recovers the previous Keller-constrained theorem as a strict
corollary of Theorem 1.1.

## 7. Relation to the complete frontier

The theorem closes every two-pair polynomial of dual degree one.  The
unrestricted conjecture nevertheless fails in bidegree \((4,4)\) by the
[full-rank counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
The safe-stratum classification beginning at dual degree two is organized by

\[
 \operatorname{End}(\operatorname{Sym}^d)
 \cong\bigoplus_{r=0}^{d}\operatorname{Sym}^{2r}.
\]

For \(d\geq2\), there is no vector field \(H\) and no two-by-two Jacobian
whose determinant is exposed by the second moment.  The formerly proposed
replacement target was the all-degree moment--nullcone theorem.  It is
already false at \(d=3\) by the full-rank semistable Rodrigues survivor and
at every \(d\geq4\) by the propagated bidegree-\((4,4)\) counterexample.
The full \(d=3\) **SIC** classification remains open: the Rodrigues orbit
itself is safe by its all-order integration-by-parts cutoff.

Thus this result moves the exact safe boundary from “Keller provenance” to
“all dual-linear polynomials,” while the known **SIC** counterexample begins
only at dual degree four.

## 8. Reproduction

Run

```bash
python3 scripts/verify_dual_linear_sic2.py
```

The dependency-free checker verifies (2.2) on a generic exact polynomial
pair and replays the general normal form, all pure contractions, the mixed
cutoff, and the sharper Keller subcase on exact integer examples.  The
checker is a regression; the all-degree proof is Sections 2--6.
