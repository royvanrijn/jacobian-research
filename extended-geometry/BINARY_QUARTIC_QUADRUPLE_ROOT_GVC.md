# Binary quartic-leading GVC with a quadruple root

## 1. Theorem and scope

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_4+\Lambda_5+\cdots
\]
be a constant-coefficient operator in two variables whose lowest nonzero
symbol is a fourth power, and let \(\deg P=5\).

> **Theorem 1.1 — quadruple-root quartic-leading quintic GVC.**
> If
> \[
>  \Lambda^m(P^m)=0\qquad(1\le m\le5),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Thus no genuinely nonhomogeneous binary GVC counterexample occurs on the
quartic root partition \((4)\), even with arbitrary higher operator jets and
arbitrary lower polynomial pieces.  The later
[triple-plus-simple](BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md) and
[double-root](BINARY_QUARTIC_DOUBLE_ROOT_GVC.md) theorems close the other
three repeated-root partitions.

## 2. Weierstrass and polynomial normalization

After scalar extension and a linear change, normalize
\(\Lambda_4=X^4\).  Formal Weierstrass division gives
\[
 \lambda(X,Y)=U(X,Y)W(X,Y),
\qquad
 W=X^4+a(Y)X^3+b(Y)X^2+c(Y)X+d(Y),
\tag{2.1}
\]
where \(U(0,0)=1\) and
\[
 \operatorname{ord}a\ge2,\quad
 \operatorname{ord}b\ge3,\quad
 \operatorname{ord}c\ge4,\quad
 \operatorname{ord}d\ge5.
\tag{2.2}
\]
The differential unit \(U(\partial)\) is injective on polynomials and
commutes with \(W(\partial)\), so it does not change either the pure
vanishing premise or the eventual mixed conclusion.

The leading pure equation gives
\[
 P_5=y^2C_3(x,y).
\tag{2.3}
\]
A shear \(x\mapsto x+ty\) preserves \(X^4\).  Choosing \(t\) away from the
finitely many roots of \(C_3(t,1)\), then rescaling, gives
\[
 P_5=y^5+Axy^4+Bx^2y^3+Cx^3y^2.
\tag{2.4}
\]
This chart therefore covers every nonzero \(P_5\) after scalar extension.

Write the first normalized operator jet as
\[
 W_5=\ell_0Y^5+\ell_1XY^4+\ell_2X^2Y^3+\ell_3X^3Y^2.
\tag{2.5}
\]

## 3. The defect-one radical

For \(m\ge2\), the defect-one coefficient is
\[
 X^{4(m-1)}W_5(P_5^m).
\tag{3.1}
\]
The possible \(P_4\) contribution vanishes by \(x\)-degree:
\[
 3(m-1)+4<4m.
\tag{3.2}
\]
Exact characteristic-zero elimination of the coefficients of (3.1) for
\(m=2,3\) gives
\[
\sqrt{J_1}
=
(C\ell_2,C\ell_1,C\ell_0,B\ell_0).
\tag{3.3}
\]
Its three minimal components are
\[
 (B,C),\qquad (C,\ell_0),\qquad
 (\ell_0,\ell_1,\ell_2).
\tag{3.4}
\]
They are respectively the \(A\ell_0\) crossing, the \(B\ell_1\) crossing,
and the \(C\ell_3\) threshold face, together with their boundary
specializations.

## 4. The two lower threshold branches

### 4.1 The \(A\ell_0\) branch

Put \(B=C=0\).  If \(A=0\), moment one solves the \(x^4\)-coefficient of
\(P_4\) as \(p_4=-5\ell_0\), and moment two is a nonzero scalar multiple of
\(\ell_0^2\).

If \(A\ne0\), normalize \(A=1\).  Moment one gives
\[
 p_4=-5\ell_0-\ell_1.
\tag{4.1}
\]
On \(\ell_0\ne0\), moment two then solves
\[
 p_3=
 -\frac{
 795\ell_0^2+310\ell_0\ell_1+28\ell_0\ell_2+19\ell_1^2
 }{2\ell_0}.
\tag{4.2}
\]
After (4.1)--(4.2), the primitive coefficient identity is
\[
 [y]\,W^4(P^4)\doteq\ell_0^3.
\tag{4.3}
\]
Here \(\doteq\) means equality up to a nonzero rational scalar.  The
complete \(W_6,W_7\) jets and the complete \(P_3,P_2\) pieces were retained
in (4.3); none can reach this output monomial.  Jets of still larger defect
cannot enter a degree-one output of moment four.  Hence \(\ell_0=0\).

With \(\ell_0=0\), moment one gives \(p_4=-\ell_1\), and moment two has
primitive scalar \(\ell_1^2\).  Thus \(\ell_1=p_4=0\).

Now use weights
\[
 w(x)=3,\qquad w(y)=2.
\tag{4.4}
\]
Every monomial of \(P\) has weight at most \(11\), while every remaining
monomial of \(W\) has weight at least \(12\).  Therefore
\[
 W^m(QP^m)=0\qquad(m>w(Q)).
\tag{4.5}
\]

### 4.2 The \(B\ell_1\) branch

Put \(C=\ell_0=0\).  The case \(B=0\) is the preceding branch, so normalize
\(B=1\).  Moment one gives
\[
 p_4=-A\ell_1-\frac12\ell_2.
\tag{4.6}
\]
Retaining the complete \(W_6\) and \(P_3\) data, moment three has
\[
 [y]\,W^3(P^3)\doteq\ell_1^2.
\tag{4.7}
\]
Thus \(\ell_1=0\).  The remaining scalar equations from moments two and
three are
\[
 24k_0+13\ell_2^2=0,\qquad
 \ell_2(12k_0+\ell_2^2)=0,
\tag{4.8}
\]
where \(k_0=[Y^6]W_6\).  Their nonzero ratios conflict, so
\[
 \ell_2=k_0=p_4=0.
\tag{4.9}
\]

Use weights \(w(x)=2,w(y)=1\), with threshold \(W=7\).  The equality faces
are
\[
\begin{aligned}
 P_{[7]}&\in\langle x^2y^3,x^3y\rangle,\\
 W_{[7]}&\in\langle XY^5,Y^7\rangle.
\end{aligned}
\tag{4.10}
\]
Every equality operator selection has \(y\)-order at least five, while
every equality polynomial selection has \(y\)-degree at most three.  The
weight-defect lemma therefore leaves a linear \(y\)-derivative deficit and
proves the mixed conclusion.

## 5. The terminal \(C\ell_3\) face

On the third component of (3.4), put
\[
 \ell_0=\ell_1=\ell_2=0.
\]
If \(C=0\), the preceding branches apply.  Normalize \(C=1\), and put
\(t=\ell_3\).  Moment one gives \(p_4=-t/2\).

Write \(k_i=[X^iY^{6-i}]W_6\).  Two positive-degree coefficients of moment
three give
\[
 k_0=0,\qquad 7Bk_0+3k_1=0,
\tag{5.1}
\]
so \(k_0=k_1=0\).  Moment two then gives
\[
 k_2=-\frac5{24}t^2.
\tag{5.2}
\]

Let \(h_i=[X^iY^{7-i}]W_7\).  The scalar third moment gives
\[
 h_1=\frac{47}{144}t^3-\frac73Bh_0.
\tag{5.3}
\]
After (5.1)--(5.3), the degree-one fourth-moment coefficient is
\[
 [y]\,W^4(P^4)\doteq h_0.
\tag{5.4}
\]
Hence \(h_0=0\), and (5.3) becomes \(h_1=47t^3/144\).

All operator terms of weight below eight have now vanished.  The weight-eight
faces are
\[
\begin{aligned}
 P_{[8]}&=x^3y^2-\frac t2x^4,\\
 W_{[8]}&=
 X^4+tX^3Y^2-\frac5{24}t^2X^2Y^4
 +\frac{47}{144}t^3XY^6+z_0Y^8.
\end{aligned}
\tag{5.5}
\]
Their fourth moment is
\[
 W_{[8]}^4(P_{[8]}^4)
 =
 958003200(40541t^4+80640z_0).
\tag{5.6}
\]
Thus
\[
 z_0=-\frac{40541}{80640}t^4.
\tag{5.7}
\]
The terminal fifth moment is
\[
 W_{[8]}^5(P_{[8]}^5)
 =
 -19931886558904320000\,t^5.
\tag{5.8}
\]
Consequently \(t=0\), and then (5.2), (5.3), and (5.7) kill every remaining
weight-eight correction.

Every later Weierstrass jet has weight strictly greater than eight.  The
only equality pair is now \((X^4,x^3y^2)\), which has a linear
\(x\)-derivative deficit.  The weight-defect lemma proves (1.2).

## 6. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quartic_quadruple_root_gvc.py
```

The checker uses exact sparse SymPy contraction and Singular over
\(\mathbb Q\).  It verifies the radical (3.3), the branch identities
(4.3), (4.7)--(4.8), (5.1)--(5.4), the terminal values (5.6)--(5.8), and
the final weight separators.

The earlier modular search in
[`BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md`](BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md)
is now only discovery history for this orbit.  Its negative samples are not
used in the proof.
