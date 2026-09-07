# The minimal coupled Danielewski ansatz

The residue-period obstruction eliminates every contraction that preserves
the suspension coordinate \(a\), and more generally every target triple in
which at most one coordinate depends on \(a\).  This note treats the first
finite ansatz that escapes that obstruction.

Use the two-dicritical arithmetic model

\[
P(x)=x(x^2+1),\qquad
x=bc,\qquad
w=b(x^2+1).
\]

Inside the source coordinate ring \(k[a,b,c]\), consider

\[
\mathcal V=
\operatorname{span}_k
\{1,a,c,x,w,ac,ax,aw\}.                            \tag{1}
\]

The terms \(ac,ax,aw\) allow two or all three target coordinates to have
independent \(a\)-dependence.  Thus (1) is the minimal genuinely coupled
bilinear extension of the affine-linear Danielewski coordinates.

## 1. Determinants as a Pluecker ledger

Discard the constant basis vector and write

\[
(f_1,\ldots,f_7)=(a,c,x,w,ac,ax,aw).
\]

Let

\[
M=
\frac{\partial(f_1,\ldots,f_7)}
     {\partial(a,b,c)}
\]

be the \(7\) by \(3\) source Jacobian.  If

\[
(Y_1,Y_2,Y_3)^T=C(f_1,\ldots,f_7)^T+\text{constants}
\]

for a \(3\) by \(7\) constant matrix \(C\), Cauchy--Binet gives

\[
\det\frac{\partial(Y_1,Y_2,Y_3)}
         {\partial(a,b,c)}
=
\sum_{\substack{I\subset\{1,\ldots,7\}\\|I|=3}}
p_I(C)\,\Delta_I,                                  \tag{2}
\]

where \(p_I(C)\) is a \(3\) by \(3\) Pluecker minor of \(C\) and
\(\Delta_I\) is the corresponding \(3\) by \(3\) minor of \(M\).

The determinant problem is therefore linear in the Pluecker coordinates,
even though decomposability of those coordinates is nonlinear.

## 2. The coefficient lock

For a polynomial \(H(a,b,c)\), write \([m]H\) for the coefficient of the
monomial \(m\).  Direct differentiation gives

\[
\boxed{
[b^2c^2]\Delta_I=3[1]\Delta_I
\quad\text{for all }I.
}                                                    \tag{3}
\]

There is only one minor with either coefficient nonzero:

\[
\Delta_{\{a,c,w\}}
=
\det\frac{\partial(a,c,w)}{\partial(a,b,c)}
=-(1+3b^2c^2).                                     \tag{4}
\]

Every other basic minor has zero constant coefficient and zero
\(b^2c^2\) coefficient.  By (2), the same relation holds for every target
triple in \(\mathcal V\):

\[
\boxed{
[b^2c^2]J=3[1]J.
}                                                    \tag{5}
\]

If \(J\) were a nonzero constant \(\lambda\), its left side would be zero
and its right side would be \(3\lambda\), a contradiction in
characteristic zero.

## 3. Classification theorem

**Theorem.**  Let \(k\) have characteristic zero.  No three polynomials

\[
Y_1,Y_2,Y_3\in
\operatorname{span}_k\{1,a,c,x,w,ac,ax,aw\},
\]

with \(x=bc\) and \(w=b((bc)^2+1)\), have nonzero constant Jacobian with
respect to \((a,b,c)\).

The proof is exactly (2)--(5).  Notice that it does not merely find no
decomposable Pluecker point.  The affine linear space cut out by the
constant-Jacobian equations is already empty before the Grassmannian
relations are imposed.

The same proof permits adjoining \(a^2\): every minor involving \(a^2\)
but not \(a\) has a factor \(a\), while minors involving both have
proportional \(da\)-rows and vanish.  Thus the extension
\(\mathcal V+k a^2\) is also excluded.

## 4. The maximal diagonal-quadratic enlargement

The coefficient ledger becomes stronger after adjoining all the
ambient-quadratic monomials except \(cx\) and \(xw\).  Put

\[
\begin{aligned}
\mathcal V_{\rm diag}=\operatorname{span}_k\{&
1,a,c,x,w,ac,ax,aw,\\
&a^2,c^2,x^2,w^2,cw\}.
\end{aligned}                                      \tag{6}
\]

There are \(\binom{12}{3}=220\) basic nonconstant Jacobian minors.  Every
one satisfies

\[
\boxed{
[x^2]J-[x^4]J+[x^6]J=3[1]J,
\qquad x=bc.
}                                                    \tag{7}
\]

Thus no triple in \(\mathcal V_{\rm diag}\) has nonzero constant
Jacobian.  Among the full ambient-quadratic basis in \(a,c,x,w\), the only
monomials not yet covered by (7) are

\[
cx,\qquad xw.                                       \tag{8}
\]

These are exactly the individual quadratic directions whose new minors
escape the linear coefficient functional in (7).

## 5. The minimal escape chart also fails

Add both directions in (8), but no diagonal quadratic:

\[
\mathcal V_{\rm esc}
=\operatorname{span}_k
\{1,a,c,x,w,ac,ax,aw,cx,xw\}.                      \tag{9}
\]

Any nonzero constant-Jacobian triple in this space has an invertible
linear part.  Indeed, the only source-linear terms among the generators
are \(a,c,w=b+\text{higher terms}\).  A constant target-linear change
therefore normalizes the triple to

\[
\begin{aligned}
A&=a+\sum_{j=1}^5 A_jN_j,\\
C&=c+\sum_{j=1}^5 B_jN_j,\\
W&=w+\sum_{j=1}^5 C_jN_j,
\end{aligned}
\qquad
(N_1,\ldots,N_5)=(ac,ax,aw,cx,xw).                 \tag{10}
\]

Its desired constant is then \(-1\).  Six linear coefficient equations
give

\[
\begin{gathered}
A_1=A_3=B_5=C_4=0,\qquad
B_1=-C_3,\\
A_2=-2(C_5+B_4).                                   \tag{11}
\end{gathered}
\]

Set

\[
r=C_5,\qquad s=B_4,\qquad
p=A_4B_2,\qquad q=C_2A_5.
\]

Five remaining coefficients of
\(\det\partial(A,C,W)/\partial(a,b,c)+1\) give

\[
\begin{aligned}
E_2&=2q+4r^2+5rs+2p+4s^2-3=0,\\
E_4&=4q+8r^2+11rs+4p+8s^2=0,\\
K&=qs+2r^2s+rp+2rs^2=0,\\
E_3&=3K+2(r+s)=0,\\
C_2A_4&=0.                                         \tag{12}
\end{aligned}
\]

The first two equations imply

\[
E_4-2E_2=rs+6=0.
\]

The next two imply \(r+s=0\), hence \(r^2=6\).  Substituting
\(s=-r\) into \(K=0\) gives \(p=q\).  The last equation forces
\(p=q=0\): if \(C_2=0\), then \(q=0\), while if \(A_4=0\), then
\(p=0\).

Finally \(E_2\) becomes

\[
4(6)+5(-6)+4(6)-3=15=0,
\]

a contradiction.  Therefore \(\mathcal V_{\rm esc}\) also contains no
constant-Jacobian triple.

Combining (7) and (12), any candidate in the full ambient-quadratic space
must use at least one escape direction \(cx\) or \(xw\) and at least one
diagonal direction from

\[
a^2,c^2,x^2,w^2,cw.                                \tag{13}
\]

The subsequent
[quadratic interaction theorem](DANIELEWSKI_QUADRATIC_INTERACTION_FRONTIER.md)
checks all 68 charts containing at most four extra directions.  It raises
this to a support lower bound of five: one escape plus four diagonal
directions, or both escapes plus three diagonal directions.

## 6. What the locks mean

The distinguished minor (4) is the ordinary contraction

\[
(a,c,w),
\]

whose determinant is \(-P'(bc)\).  The cubic correction \(3b^2c^2\) is
the finite-chart trace of the residue factor.  The three first mixed
coordinates \(ac,ax,aw\) create many new minors, but none creates a bare
\(b^2c^2\) term capable of cancelling it without also changing the
constant coefficient.

Merely activating \(a\) in several linear Danielewski directions is
insufficient.  Diagonal quadratic terms alone are also insufficient, and
the two minimal escape terms alone are insufficient.  Only their
interaction remains at ambient degree two.

Any larger search should first compute its analogue of the linear
functional

\[
\ell(H)=[b^2c^2]H-3[1]H.                           \tag{14}
\]

Only basis extensions whose new Jacobian minors escape
\(\ker\ell\) can possibly remove this obstruction.  This provides a
cheap exact prefilter before any Grassmannian or coefficient-ideal
calculation.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_minimal_coupled_ansatz.py
```

The checker constructs all 35 basic minors, verifies the first universal
coefficient lock, confirms that (4) is its unique nonzero carrier, checks
all 220 minors of the diagonal-quadratic extension and their stronger
lock, verifies the normalized minimal-escape equations and their unit
ideal, and performs deterministic random-matrix Cauchy--Binet regressions.
