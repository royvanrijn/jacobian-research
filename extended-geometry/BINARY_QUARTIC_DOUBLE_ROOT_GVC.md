# Binary quartic-leading GVC with maximal root multiplicity two

## 1. Theorem and consequence

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_4+\Lambda_5+\cdots
\]
be a constant-coefficient operator in two variables, let \(\deg P=5\),
and suppose the lowest nonzero symbol \(\Lambda_4\) has root partition
\((2+2)\) or \((2+1+1)\).

> **Theorem 1.1 — remaining repeated-root quartic GVC.**
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

Arbitrary lower pieces of \(P\) and arbitrary higher operator jets are
included.

Together with the
[`quadruple-root theorem`](BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md),
the
[`triple-plus-simple theorem`](BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md),
and the squarefree row of the
[`degree-five frontier theorem`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md),
this closes every quartic-leading binary quintic row.

> **Corollary 1.2 — binary GVC through polynomial degree five.**
> Every constant-coefficient operator in two variables satisfies the GVC
> conclusion for every polynomial \(P\) of degree at most five.

Indeed, the earlier lowest-order and degree-five frontier theorems already
close lowest positive orders other than the quartic-leading row.

## 2. The \((2+2)\) orbit

Normalize
\[
 \Lambda_4=X^2Y^2.
\tag{2.1}
\]
Formal division by the leading monomial, followed by removal of the
differential unit, gives
\[
 W=X^2Y^2+
 \sum_{d\ge5}
 \left(
  a_{d,0}Y^d+a_{d,1}XY^{d-1}
  +a_{d,d-1}X^{d-1}Y+a_{d,d}X^d
 \right).
\tag{2.2}
\]
The leading pure-zero locus is
\[
 P_5=x^4L(x,y)\qquad\text{or}\qquad P_5=y^4L(x,y).
\tag{2.3}
\]
The two families are exchanged by \(x\leftrightarrow y\), so write
\[
 P_5=Cx^5+Ax^4y.
\tag{2.4}
\]
Put
\[
 W_5=\ell_0Y^5+\ell_1XY^4+\ell_4X^4Y+\ell_5X^5.
\tag{2.5}
\]

### 2.1 The linear-tip chart \(A\ne0\)

Normalize \(A=1\).  The defect-one equations give
\[
 p_0=0,\qquad p_1=-28\ell_5,\qquad
 p_2=-6\ell_4-30C\ell_5.
\tag{2.6}
\]
Write
\[
 W_6=k_0Y^6+k_1XY^5+k_5X^5Y+k_6X^6
\]
and
\[
 P_3=r_0y^3+r_1xy^2+r_2x^2y+r_3x^3.
\]
The defect-two equations include
\[
\begin{aligned}
 [y]\,W^3(P^3)&\doteq\ell_5^2,\\
 W^2(P^2)&\doteq
 2550C^2\ell_5^2+360C\ell_4\ell_5+140k_6
 -2\ell_4^2+7\ell_5p_3+r_0.
\end{aligned}
\tag{2.7}
\]
Thus \(\ell_5=0\), and
\[
 r_0=2\ell_4^2-140k_6.
\tag{2.8}
\]

Use weights \(w(x)=1,w(y)=2\), with threshold six.  The complete equality
faces are
\[
\begin{aligned}
 P_{[6]}&=x^4y-6\ell_4x^2y^2
 +(2\ell_4^2-140k_6)y^3,\\
 W_{[6]}&=X^2Y^2+\ell_4X^4Y+k_6X^6.
\end{aligned}
\tag{2.9}
\]
Putting \(u=\ell_4,v=k_6\), their third and fourth moments are, up to
nonzero rational factors,
\[
\begin{aligned}
 F(u,v)&=u(u^2-28v),\\
 G(u,v)&=2789u^4-71680u^2v-62320v^2.
\end{aligned}
\tag{2.10}
\]
If \(u=0\), then \(G=-62320v^2\).  If \(v=u^2/28\), then
\[
 G=\frac{7326}{49}u^4.
\tag{2.11}
\]
Hence \(u=v=0\).

Every remaining operator term has weight greater than six.  The equality
pair \((X^2Y^2,x^4y)\) has a \(y\)-derivative deficit.

### 2.2 The fifth-power chart \(A=0\)

Normalize \(C=1\), so \(P_5=x^5\).  Moment one gives
\[
 p_2=-30\ell_5.
\tag{2.12}
\]
The defect-two equations include
\[
\begin{aligned}
 W^2(P^2)&\doteq
 15\ell_4p_1+2550\ell_5^2+p_0p_4+p_1p_3,\\
 [y]\,W^3(P^3)&\doteq p_0p_1,\\
 [x]\,W^3(P^3)&\doteq12\ell_5p_0+p_1^2,\\
 [x^2]\,W^4(P^4)&\doteq p_0^2.
\end{aligned}
\tag{2.13}
\]
They force
\[
 p_0=p_1=\ell_5=p_2=0.
\tag{2.14}
\]
Weights \(w(x)=2,w(y)=3\), with threshold ten, now leave only
\((X^2Y^2,x^5)\) on the equality face.  It has a \(y\)-derivative
deficit.  The \(y^5\) chart follows by symmetry.

## 3. The \((2+1+1)\) orbit

Normalize
\[
 \Lambda_4=X^2Y(X-Y)=X^3Y-X^2Y^2.
\tag{3.1}
\]
Using \(X^3Y\) as the division monomial gives the differential-unit normal
form
\[
 W=X^3Y-X^2Y^2+
 \sum_{d\ge5}
 \left(
  a_{d,0}Y^d+a_{d,1}XY^{d-1}
  +a_{d,2}X^2Y^{d-2}+a_{d,d}X^d
 \right).
\tag{3.2}
\]
The exact leading moment radical has the three components
\[
 P_5=x^5,\qquad P_5=(x+y)^5,\qquad
 P_5=y^4(Cy+Ax).
\tag{3.3}
\]

### 3.1 The linear-tip chart \(A\ne0\)

Normalize \(A=1\).  The defect-one equations give
\[
 p_4=0,\qquad p_3=28\ell_0,\qquad
 p_2=42\ell_0+6\ell_1+30C\ell_0.
\tag{3.4}
\]
The defect-two third moment contains
\[
 [x]\,W^3(P^3)\doteq\ell_0^2,
\tag{3.5}
\]
so \(\ell_0=0\).  The scalar second moment then gives
\[
 r_3=140k_0+2\ell_1^2.
\tag{3.6}
\]

With weights \(w(x)=2,w(y)=1\) and threshold six, the complete equality
faces are
\[
\begin{aligned}
 P_{[6]}&=xy^4+6\ell_1x^2y^2
 +(140k_0+2\ell_1^2)x^3,\\
 W_{[6]}&=-X^2Y^2+\ell_1XY^4+k_0Y^6.
\end{aligned}
\tag{3.7}
\]
For \(u=\ell_1,v=k_0\), moments three and four are proportional to
\[
\begin{aligned}
 F(u,v)&=u(u^2+28v),\\
 G(u,v)&=2789u^4+71680u^2v-62320v^2.
\end{aligned}
\tag{3.8}
\]
As in (2.10)--(2.11), the two branches of \(F=0\) reduce \(G\) to a
nonzero multiple of \(v^2\) or \(u^4\).  Hence \(u=v=0\).

The final equality pair \((-X^2Y^2,xy^4)\) has an \(x\)-derivative
deficit.

### 3.2 The pure-\(y\) chart \(A=0\)

Normalize \(C=1\).  Defect one gives
\[
 p_4=0,\qquad
 p_2=\frac{3p_3+60\ell_0}{2}.
\tag{3.9}
\]
The defect-two equations contain
\[
\begin{aligned}
 [y]\,W^3(P^3)&\doteq p_3^2,\\
 W^2(P^2)\big|_{p_3=0}&\doteq\ell_0^2.
\end{aligned}
\tag{3.10}
\]
Thus \(p_3=\ell_0=p_2=0\).

Weights \(w(x)=3,w(y)=2\), with threshold ten, leave only
\((-X^2Y^2,y^5)\) on the equality face.  It has an \(x\)-derivative
deficit.

### 3.3 The isolated fifth powers

For \(P_5=x^5\), defect one gives
\[
 p_0=p_1=p_2=0,\qquad p_3=-20\ell_5.
\tag{3.11}
\]
The next defect gives
\[
 r_0=0,\qquad r_1=-340\ell_5^2.
\tag{3.12}
\]
With weights \(w(x)=1,w(y)=2\), the complete equality faces are
\[
\begin{aligned}
 P_{[5]}&=x^5-20\ell_5x^3y-340\ell_5^2xy^2,\\
 W_{[5]}&=X^3Y+\ell_5X^5.
\end{aligned}
\tag{3.13}
\]
Their third moment is
\[
 W_{[5]}^3(P_{[5]}^3)
 =301335552000\,\ell_5^3.
\tag{3.14}
\]
Hence \(\ell_5=0\).  The remaining pair \((X^3Y,x^5)\) has a
\(y\)-derivative deficit.

The involution
\[
 (X,Y)\longmapsto(X,X-Y)
\tag{3.15}
\]
preserves \(X^2Y(X-Y)\) and exchanges the two simple-root annihilator
lines.  Its contragredient primal action sends the \(x^5\) component to
the \((x+y)^5\) component, so (3.11)--(3.14) close both.

## 4. All-order conclusion

After the displayed eliminations, every chart has a positive weight and
threshold satisfying
\[
 w(P)\le W\le w(\Lambda).
\tag{4.1}
\]
All strict selections have positive weight surplus.  Every equality
selection has a coordinate-derivative deficit linear in the number of
selections.  The weight-defect lemma from the
[`degree-five frontier`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md)
therefore proves (1.2).

## 5. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quartic_double_root_gvc.py
```

The checker performs exact sparse contractions over \(\mathbb Q\).  It
retains every normalized higher operator and lower polynomial jet capable
of entering each displayed defect coefficient, verifies the two
two-parameter equality systems, and checks all final weighted faces.
