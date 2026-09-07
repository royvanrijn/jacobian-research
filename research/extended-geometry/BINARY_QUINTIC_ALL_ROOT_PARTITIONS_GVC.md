# Binary sextic GVC for every quintic leading symbol

## 1. Theorem

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_5+\Lambda_6+\cdots
\]
be a constant-coefficient operator in two variables with lowest positive
order five, and let \(\deg P=6\).

> **Theorem 1.1 — complete quintic-leading sextic row.**
> If
> \[
>  \Lambda^m(P^m)=0\qquad(m\ge1),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Thus none of the seven binary-quintic root partitions
\[
 (5),\ (4+1),\ (3+2),\ (3+1+1),\ (2+2+1),\
 (2+1+1+1),\ (1+1+1+1+1)
\tag{1.3}
\]
supports a genuinely nonhomogeneous binary GVC counterexample with
\(\deg P=6\).  Arbitrary lower pieces of \(P\) and arbitrary higher
operator jets are included.

The all-order premise is used once, to classify the leading homogeneous
pair.  After that classification, the local correction calculations use
moments only through order six.  The quintuple-root component is the
unique component that needs moment six.

## 2. The leading locus is a Hall-matching locus

After scalar extension, split
\[
 \Lambda_5=\prod_{i=1}^5D_{v_i},
 \qquad
 P_6=\prod_{j=1}^6L_j.
\tag{2.1}
\]
Retain the translation variable in the split-symbol constant-term
construction:
\[
 H_z(t)=
 \frac{P_6(z+\sum_i t_iv_i)}{t_1t_2t_3t_4t_5}.
\tag{2.2}
\]
The top homogeneous part of (1.1) and the
Duistermaat--van der Kallen theorem give
\[
 0\notin\operatorname{Newt}(H_z)
\tag{2.3}
\]
at a generic \(z\).

The Newton polytope before division by \(t_1\cdots t_5\) is the
Minkowski sum of the six simplices
\[
 \operatorname{conv}\bigl(
 0,\ e_i:\ L_j(v_i)\ne0
 \bigr).
\tag{2.4}
\]
The vector \((1,1,1,1,1)\) belongs to this sum exactly when the five
derivative copies can be matched to five distinct factors \(L_j\), with
\(L_j(v_i)\ne0\).  This is the integral bipartite matching polytope.

Suppose a root direction \(v\) occurs in \(\Lambda_5\) with multiplicity
\(e\), and let \(c\) be the number of factors \(L_j\) annihilating \(v\).
A Hall-deficient subset cannot contain two nonparallel derivative
directions: every nonzero \(L_j\) interacts with at least one of them, so
its neighborhood then has all six vertices.  The only possible deficient
subset consists of the \(e\) copies of one direction.  Hall's inequality
fails exactly when
\[
 6-c<e,
 \qquad\text{or equivalently}\qquad
 c\ge7-e.
\tag{2.5}
\]

Consequently the leading pure-zero locus is precisely the union, over
the roots of multiplicity \(e\), of
\[
 P_6=y^{\,7-e}C_{e-1}(x,y)
\tag{2.6}
\]
after coordinates put that root direction at \(X\).  For multiplicities
\(e=1,2,3,4,5\), these are respectively
\[
 y^6,\quad y^5C_1,\quad y^4C_2,\quad
 y^3C_3,\quad y^2C_4.
\tag{2.7}
\]
This simultaneously classifies the leading loci on all seven
partitions (1.3).

## 3. Local normal form and weights

Fix one component (2.6).  Normalize the exposed monomial of the leading
symbol to
\[
 M_e=X^eY^{5-e}.
\tag{3.1}
\]
The other terms of \(\Lambda_5\) have the form
\[
 X^{e+j}Y^{5-e-j},\qquad j\ge1.
\tag{3.2}
\]
Local formal division gives
\[
 \lambda=U W,\qquad U(0,0)=1,
\tag{3.3}
\]
where later normalized terms of \(W\) are supported on
\[
 X^aY^{d-a},\qquad a<e\ \text{or}\ d-a<5-e.
\tag{3.4}
\]
The locally finite differential unit \(U(\partial)\) is invertible on
polynomials and changes neither premise nor conclusion.

On the projective chart whose highest \(x\)-term in (2.6) is
\(x^iy^{6-i}\), put
\[
 k=e-i,\qquad
 w(x)=k+1,\qquad w(y)=k,\qquad W_0=5k+e.
\tag{3.5}
\]
Then
\[
 w(M_e)=w(x^iy^{6-i})=W_0,
\tag{3.6}
\]
while every term (3.2) has weight \(W_0+j\).  Formula (3.4) shows that
only finitely many higher jets can lie on or below the threshold.

The remainder of the proof is therefore local in the root multiplicity
\(e\).  The checker retains generic coefficients in (3.2), so the
calculations below apply equally to monomial, repeated-root, and
squarefree quintic symbols.

## 4. Multiplicity four

Take \(M_4=X^4Y\) and
\[
 P_6=C y^6+Axy^5+Bx^2y^4+Dx^3y^3.
\tag{4.1}
\]
Write
\[
 W_6=\ell_0Y^6+\ell_1XY^5+\ell_2X^2Y^4
 \ell_3X^3Y^3+\ell_6X^6.
\tag{4.2}
\]
Moment one solves
\[
 p_4=-5A\ell_1-2B\ell_2-30C\ell_0-\frac32D\ell_3.
\tag{4.3}
\]
Exact defect-one moments through order four give
\[
 \sqrt{J_1}=
 (D\ell_1,D\ell_0,B\ell_0,D(15D\ell_2+7p_5)).
\tag{4.4}
\]
Its three minimal components are
\[
 (D,B),\qquad(D,\ell_0),\qquad
 (\ell_0,\ell_1,15D\ell_2+7p_5).
\tag{4.5}
\]

On \(D\ne0\), normalize \(D=1\).  Complete defect-two equations give
\[
 h_0=-\frac{151}{196}\ell_2^2,\qquad
 355\ell_2^2=0,\qquad h_1=0,\qquad
 h_2=\frac1{40}\ell_3^2.
\tag{4.6}
\]
Thus the migrating \(\ell_2/p_5\) pair vanishes.  The remaining
weight-nine terms form a one-sided equality face.

On \(D=0,B\ne0\), normalize \(B=1\).  Defect two contains
\[
 330\ell_1^2+42\ell_1p_5+11p_5^2.
\tag{4.7}
\]
If \(\ell_1=0\), then \(p_5=0\), and the remaining equation gives
\(h_0=-13\ell_2^2/42\).  If \(\ell_1\ne0\), put
\(s=p_5/\ell_1\).  A complete defect-three coefficient is
\[
 221s^3+546s^2+2970s+21840.
\tag{4.8}
\]
The resultant of (4.7) and (4.8) is
\[
 1909272615840\ne0.
\tag{4.9}
\]

On \(D=B=0,A\ne0\), the only migrating pair is
\((u,v)=(\ell_0,p_5)\).  Extremal coefficients of moments four and five
are, up to nonzero constants,
\[
\begin{aligned}
 G_3={}&13v^3+54uv^2+1638u^2v+302328u^3,\\
 G_4={}&323v^4+680uv^3+8580u^2v^2\\
 &+465120u^3v+98062800u^4.
\end{aligned}
\tag{4.10}
\]
At \(u=0\), \(G_3\) gives \(v=0\).  At \(u\ne0\), their ratio-polynomial
resultant is
\[
 340491886133329409922608332032\ne0.
\tag{4.11}
\]
Then the defect-two scalar forces \(\ell_1=0\).

On the pure-\(y\) chart, a defect-four coefficient of moment five is a
nonzero multiple of \(p_5^4\).  After \(p_5=0\), defect two is a nonzero
multiple of \(\ell_0^2\).

The isolated simple-root component \(P_6=x^6\) has equality chain
\[
\begin{aligned}
 P_{[6]}={}&x^6-30t x^4y-720t^2x^2y^2-154320t^3y^3,\\
 W_{[6]}={}&X^4Y+tX^6.
\end{aligned}
\tag{4.12}
\]
Moment four is a nonzero multiple of \(t^4\).

## 5. Multiplicities three and two

For multiplicity three, take \(M_3=X^3Y^2\) and
\[
 P_6=C y^6+Axy^5+Bx^2y^4.
\tag{5.1}
\]
Allow the generic strict cofactor terms
\(\alpha X^4Y+\beta X^5\).  They do not change the radical
\[
 \sqrt{J_1}=
 (Bp_5,Ap_5,B\ell_0,B(56B\ell_1+5p_4)).
\tag{5.2}
\]

On \(B\ne0\), a defect-two positive-degree coefficient is
\(\ell_1^2\); after \(\ell_1=0\), the only remaining non-strict equation
is
\[
 56h_1-4\ell_2^2+q_4=0.
\tag{5.3}
\]
All three terms lie on the weight-eight equality face.

On \(B=0,A\ne0\), the possible \(\ell_0/p_4\) ratio satisfies
\[
\begin{aligned}
 f(s)&=7s^2+330s+30030,\\
 g(s)&=13s^3+540s^2+32760s+3023280.
\end{aligned}
\tag{5.4}
\]
Their resultant is
\[
 5436596606290200\ne0.
\tag{5.5}
\]
On the pure-\(y\) chart, extremal coefficients successively give
\[
 p_5^2=0,\qquad p_4^3=0,\qquad \ell_0^2=0.
\tag{5.6}
\]

For multiplicity two, swap coordinates and work with \(M=X^3Y^2\) and
\[
 P_6=Ax^5y+Bx^6.
\tag{5.7}
\]
Generic strict cofactor terms are retained.  The defect-one radical is
\[
\begin{aligned}
 \sqrt{J_1}=(Bp_1,Ap_1,Bp_0,Ap_0,
 A(60A\ell_6+p_2)).
\end{aligned}
\tag{5.8}
\]
On \(A\ne0\), defect two exposes \(\ell_6^2\); the remaining equality
relation is
\[
 420h_7-20\ell_5^2+q_1=0.
\tag{5.9}
\]
On the pure-\(x\) chart, moment three first gives \(p_2^2=0\).  The
weight-twelve equality chain is
\[
\begin{aligned}
 P_{[12]}&=x^6-60t x^3y^2-10620t^2y^4,\\
 W_{[12]}&=X^3Y^2+tX^6,
\end{aligned}
\tag{5.10}
\]
and its third moment is a nonzero multiple of \(t^3\).

## 6. Multiplicities one and five

At a simple root, use \(M_1=X^4Y\) and \(P_6=x^6\).  Retaining all four
strict coefficients of the quintic cofactor, defect one kills the four
nonterminal coefficients of \(P_5\) and gives \(p_4=-30t\).  Defects two
and three give
\[
 q_2=-720t^2,\qquad r_0=-154320t^3,
\tag{6.1}
\]
and moment four gives \(t^4=0\).  Thus every sixth-power component
arising from any simple root is closed.

Multiplicity five is the
[quintuple-root theorem](BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md).
Its longest branch retains normalized operator jets through order ten:
moment five solves the last \(Y^{10}\) coefficient, and moment six is a
nonzero multiple of \(56908187t^6\).

## 7. All-order conclusion

For every component, equations (3.5)--(3.6) leave a positive common
weight and threshold.  Every polynomial term above the threshold is
eliminated by the displayed exact equations.  Every unused normalized
operator jet is strict; only finitely many earlier jets lie on the
equality face, and those are exactly the terms retained above.

The generic strict cofactor terms (3.2) were retained in the
multiplicity-three, multiplicity-two, and simple-root calculations.
They do not alter the exposed squares, ratio polynomials, or terminal
moments.

Every final equality face has either:

- a coordinate-derivative deficit linear in the number of selections; or
- equal positive input and operator weight, in which case a nonzero
  output of positive ordinary degree would have weight zero, which is
  impossible.

The common weight-defect lemma therefore proves (1.2).  Since the Hall
classification covers every root of every binary quintic, the seven
partitions (1.3) are exhausted.

## 8. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quintic_all_root_partitions_gvc.py
```

The checker uses exact sparse SymPy contraction and Singular over
\(\mathbb Q\).  It brute-forces the Hall matching criterion on every
partition of five, verifies the three new defect-one radicals, retains
generic strict cofactor coefficients, replays every localized
higher-defect identity, audits the weighted supports, and invokes the
existing quintuple-root checker.
