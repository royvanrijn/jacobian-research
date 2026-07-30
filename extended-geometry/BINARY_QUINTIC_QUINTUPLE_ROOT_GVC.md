# Binary sextic GVC with a quintuple-root leading symbol

## 1. Theorem and scope

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_5+\Lambda_6+\cdots
\]
be a constant-coefficient operator in two variables whose lowest nonzero
symbol is a fifth power, and let \(\deg P=6\).

> **Theorem 1.1 — quintuple-root quintic-leading sextic GVC.**
> If
> \[
>  \Lambda^m(P^m)=0\qquad(1\le m\le6),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Thus the first \(r=5,\deg P=6\) root orbit contains no genuinely
nonhomogeneous binary GVC counterexample.  Arbitrary lower pieces of \(P\)
and arbitrary higher operator jets are included.

The later
[`all-root-partitions theorem`](BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md)
uses this result as its multiplicity-five component and closes the entire
quintic-leading sextic row.

## 2. Normalization and defect-one radical

After scalar extension and a linear change, normalize
\(\Lambda_5=X^5\).  Formal Weierstrass division gives
\[
 \lambda(X,Y)=U(X,Y)W(X,Y),
\qquad
 W=X^5+\sum_{i=0}^4 a_i(Y)X^i,
\tag{2.1}
\]
where \(U(0,0)=1\) and \(\operatorname{ord}a_i\ge6-i\).  The locally
finite differential unit \(U(\partial)\) is invertible on polynomials, so
it changes neither the pure premise nor the eventual mixed conclusion.

The top-degree part of moment one gives
\[
 X^5(P_6)=0,
\]
hence \(P_6=y^2C_4(x,y)\).  A shear preserving \(X^5\), followed by
rescaling, gives
\[
 P_6=y^6+Axy^5+Bx^2y^4+Cx^3y^3+Dx^4y^2.
\tag{2.2}
\]
Write
\[
 W_6=\ell_0Y^6+\ell_1XY^5+\ell_2X^2Y^4
 +\ell_3X^3Y^3+\ell_4X^4Y^2.
\tag{2.3}
\]
Moment one solves the \(x^5\)-coefficient of \(P_5\):
\[
 p_5=-6\ell_0-A\ell_1-\frac25B\ell_2
 -\frac3{10}C\ell_3-\frac25D\ell_4.
\tag{2.4}
\]
For moments at least two, \(P_5\) cannot enter the defect-one
positive-degree coefficients by \(x\)-degree.  Exact moments two and
three give
\[
\sqrt{J_1}=
(D\ell_0,D\ell_1,D\ell_2,D\ell_3,C\ell_0,C\ell_1).
\tag{2.5}
\]
The three minimal components are
\[
 (D,C),\qquad(D,\ell_0,\ell_1),\qquad
 (\ell_0,\ell_1,\ell_2,\ell_3).
\tag{2.6}
\]

## 3. The \(Dx^4y^2\) component

Suppose \(D\ne0\), normalize \(D=1\), and use the third component of
(2.6):
\[
 \ell_0=\ell_1=\ell_2=\ell_3=0.
\]
Put \(t=\ell_4\).  Then \(p_5=-2t/5\).

Use weights \(w(x)=2,w(y)=1\), with threshold ten.  Unlike the
degree-five quadruple-root calculation, normalized operator terms below
or on the threshold continue through order ten.  Retaining the complete
jets at every relevant defect, moments determine them successively:
\[
\begin{array}{c|c}
\text{operator order}&\text{coefficients of weight at most ten}\\ \hline
6&tX^4Y^2,\\
7&-\frac15t^2X^3Y^4,\\
8&\frac{3971}{12375}t^3X^2Y^6,\\
9&-\frac{38939}{78750}t^4XY^8,\\
10&uY^{10}.
\end{array}
\tag{3.1}
\]
More explicitly, the positive-degree coefficients first force all
lower-weight terms in orders seven through nine to vanish.  The remaining
scalar equations give the four displayed equality coefficients.

The complete weight-ten faces are therefore
\[
\begin{aligned}
 P_{[10]}={}&x^4y^2-\frac25tx^5,\\
 W_{[10]}={}&X^5+tX^4Y^2-\frac15t^2X^3Y^4
 +\frac{3971}{12375}t^3X^2Y^6\\
 &-\frac{38939}{78750}t^4XY^8+uY^{10}.
\end{aligned}
\tag{3.2}
\]
Their fifth moment is
\[
 W_{[10]}^5(P_{[10]}^5)
 \doteq4872527t^5-5906250u.
\tag{3.3}
\]
After solving (3.3), the sixth moment is a nonzero rational multiple of
\[
 56908187\,t^6.
\tag{3.4}
\]
Thus \(t=u=0\).  Every operator term of order at least eleven has weight
strictly greater than ten.  The final equality pair
\((X^5,x^4y^2)\) has an \(x\)-derivative deficit.

## 4. The \(Cx^3y^3\) component

Put \(D=0\), \(C=1\), and \(\ell_0=\ell_1=0\).  Complete defect-two
moments give
\[
 \ell_2=0,\qquad
 1400B h_0+600h_1+321\ell_3^2=0,
\tag{4.1}
\]
where \(h_i=[X^iY^{7-i}]W_7\).  Hence
\[
 h_1=-\frac73Bh_0-\frac{107}{200}\ell_3^2.
\tag{4.2}
\]
The defect-three third moment contains
\[
 1484Bh_0\ell_3+1344h_0\ell_4+280h_0p_4
 -565\ell_3^3.
\tag{4.3}
\]
If \(h_0=0\), equation (4.3) gives \(\ell_3=0\), hence \(h_1=0\).
If \(h_0\ne0\), solve (4.3) for \(p_4\).  Retaining the complete
\(W_8,W_9,P_3,P_2\) data, a positive-degree coefficient of moment five
is a nonzero rational multiple of \(h_0^2\), a contradiction.

Weights \(w(x)=3,w(y)=2\), with threshold fifteen, now leave only
\((X^5,x^3y^3)\) on the equality face.

## 5. The \(Bx^2y^4\) component

Put \(D=C=0\) and normalize \(B=1\).  Complete defect-two moments first
give \(\ell_0=0\).  Their scalar remainder is
\[
\begin{aligned}
 E={}&1525A^2\ell_1^2+1590A\ell_1\ell_2+2100\ell_1^2
 +280\ell_1\ell_3\\
 &+50\ell_1p_4+196\ell_2^2.
\end{aligned}
\tag{5.1}
\]
If \(\ell_1=0\), then \(E=0\) gives \(\ell_2=0\).  If
\(\ell_1\ne0\), solve (5.1) for \(p_4\).  With complete
\(W_7,W_8,P_4,P_3\) data retained, the defect-three fourth moment has a
positive-degree coefficient proportional to \(\ell_1^3\), a
contradiction.

Weights \(w(x)=4,w(y)=3\), with threshold twenty, leave only
\((X^5,x^2y^4)\) on the equality face.

## 6. The \(Axy^5\) and \(y^6\) components

Put \(D=C=B=0\).  If \(A\ne0\), normalize \(A=1\).  The defect-two scalar
equation is
\[
 3516\ell_0^2+1162\ell_0\ell_1+84\ell_0\ell_2
 +2\ell_0p_4+61\ell_1^2=0.
\tag{6.1}
\]
If \(\ell_0=0\), this gives \(\ell_1=0\).  If \(\ell_0\ne0\), normalize
\(\ell_0=1\).  Moments two and three solve successively for
\(p_4,p_3\).  A defect-four positive-degree coefficient of moment five
is then already the nonzero constant
\[
 9306726025998591590400000000,
\tag{6.2}
\]
so this localized branch is empty.

Weights \(w(x)=5,w(y)=4\), with threshold twenty-five, leave only
\((X^5,xy^5)\) on the equality face.

Finally, if \(A=0\), normalize \(P_6=y^6\).  The complete weight-thirty
faces begin as
\[
 P_{[30]}=y^6-6\ell_0x^5,\qquad
 W_{[30]}=X^5+\ell_0Y^6.
\tag{6.3}
\]
Their second moment is a nonzero multiple of \(\ell_0^2\), so
\(\ell_0=0\).  Weights \(w(x)=6,w(y)=5\) then leave
\((X^5,y^6)\), with an \(x\)-derivative deficit.

## 7. All-order conclusion

Each projective branch ends with a positive common weight and threshold
\[
 w(P)\le W\le w(\Lambda).
\tag{7.1}
\]
All strict selections have positive weight surplus.  Every equality face
is the pair \((X^5,x^iy^{6-i})\) for some \(0\le i\le4\), so it has an
\(x\)-derivative deficit linear in the number of selections.  The
weight-defect lemma from the
[`degree-five frontier`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md)
proves (1.2).

## 8. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quintic_quintuple_root_gvc.py
```

The checker uses exact sparse SymPy contraction and Singular over
\(\mathbb Q\).  It verifies the radical (2.5), every localized branch
identity, the complete order-six-through-ten terminal chain, and the
final weighted support inequalities.
