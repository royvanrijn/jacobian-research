# Binary septic GVC for every quadratic leading symbol

## 1. Theorem and scope

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda _2+\Lambda _3+\cdots
\]
be a constant-coefficient operator in two variables whose lowest nonzero
positive homogeneous part has order two, and let \(\deg P=7\).

> **Theorem 1.1 (complete quadratic-leading septic row).**  If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Arbitrary lower terms of \(P\), arbitrary higher operator jets, and both
root partitions \((2)\) and \((1+1)\) of \(\Lambda _2\) are included.  This
is the first complete row on the polynomial-degree-seven frontier.  Together
with the lowest-order-one theorem and the result for \(\deg P\leq r\), it
gave the strict intermediate narrowing
\[
 \boxed{
 \text{a degree-seven binary counterexample must have lowest order }
 r\in\{3,4,5,6\}.}
\tag{1.3}
\]

The later cubic- and high-order septic theorems close degree seven, and the
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md) proves unrestricted
binary GVC.  This row is retained as an independent exact regression.

The proof uses the all-order premise only to enter each successive Newton
face.  Every exact face radical below is determined by moments of order at
most nine.  The final passage to (1.2) is the common weight-defect lemma from
the [degree-five frontier](BINARY_DEGREE_FIVE_GVC_FRONTIER.md#2-weight-defect-lemma).

## 2. Hall locus and exact local division

After scalar extension, split
\[
 \Lambda _2=D_{v_1}D_{v_2},
 \qquad P_7=\prod_{j=1}^7L_j.
 \tag{2.1}
\]
If a root direction of \(\Lambda _2\) has multiplicity \(e\), and \(c\) of
the seven factors of \(P_7\) annihilate it, Hall failure is
\[
 7-c<e,
 \qquad\text{or equivalently}\qquad c\geq8-e.
 \tag{2.2}
\]
A subset containing two nonparallel derivative directions sees every
polynomial factor, so there is no further deficient subset.  Up to a linear
change of variables, the leading locus is therefore
\[
\begin{array}{c|c|c}
e&\Lambda _2&P_7\\ \hline
1&XY&y^7,\\
2&X^2&y^6(Ax+Cy).
\end{array}
\tag{2.3}
\]
On \(A\ne0\), a shear preserving \(X^2\) leaves \(P_7=xy^6\); the boundary
is \(P_7=y^7\).

Formal division by the leading symbol gives the exact normal forms
\[
\begin{aligned}
 XY\text{ chart}:&\quad W=XY+H(X)+K(Y),\\
 X^2\text{ chart}:&\quad W=X^2+a(Y)X+b(Y),
\end{aligned}
\tag{2.4}
\]
where all displayed corrections start in total order three.  The removed
unit and its inverse act locally finitely on polynomials, so (1.1)--(1.2) are
unchanged.  The Newton ledgers below are complete for (2.4): at each stage
the integer exponent solutions of the displayed weight equality are listed,
and a channel between two equality values is recorded as a one-coefficient
crossing.  After the last displayed common threshold, every omitted channel
has strict positive defect.

The Hall matching and every exponent ledger are replayed exactly by
[`verify_binary_quadratic_septic_gvc.py`](../scripts/verify_binary_quadratic_septic_gvc.py).

## 3. The distinct-root chart

Swap variables so that \(P_7=x^7\), and write
\[
 W=XY+\sum_{j\ge3}a_jX^j+\sum_{j\ge3}b_jY^j.
 \tag{3.1}
\]
Since the coefficient of \(XY\) is a unit, the first equation \(W(P)=0\)
solves all twenty-one mixed coefficients of \(P\), from high total degree to
low total degree.  Retain the operator through order fourteen; all later
terms kill \(P^2\).  Put \(p_{0j}=[y^j]P\).  Successive coefficients of the
second moment, read after the preceding rows have vanished, are
\[
\begin{array}{c|c}
\text{output coefficient}&\text{value}\\ \hline
[x^5y^4]&2520p_{06}\\
[x^8]&635040a_3^2\\
[x^5y^3]&1680p_{05}\\
[x^5y^2]&1008p_{04}\\
[x^6]&80015040a_4^2\\
[x^5y]&504p_{03}\\
[x^5]&168p_{02}\\
[x^4]&3166732800a_5^2\\
[x^2]&41912640000a_6^2\\
[1]&86467046400a_7^2.
\end{array}
\tag{3.2}
\]
Thus
\[
 P=f(x)+cy,
 \qquad
 W=Y\Gamma(X,Y)+H(X),
 \qquad \operatorname {ord}H\geq8>\deg f.
 \tag{3.3}
\]
The transverse-linear/high-order theorem in
[the separable escape note](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md) now gives
(1.2).  This closes the distinct-root partition.

## 4. The double line over \(P_7=xy^6\)

Use \(W=X^2+a(Y)X+b(Y)\).  The two half-integral faces are
\[
\begin{array}{c|c|c}
\text{slope}&W&P\\ \hline
3/2&X^2+aY^3&xy^6+cx^3y^3\\
5/2&X^2+aY^5&xy^6+cx^3y.
\end{array}
\tag{4.1}
\]
Both exact radicals are \((a,c)\).  The remaining odd pure-\(Y\) crossings
are triangular; their second-moment pivots are
\[
\begin{array}{c|c}
\text{operator correction}&\text{extremal coefficient}\\ \hline
aY^7&15966720a\\
aY^9&319334400a\\
aY^{11}&1916006400a.
\end{array}
\tag{4.2}
\]

For an integral slope \(j=2,3,4,5,6\), the complete face is
\[
\begin{aligned}
 W_j&=X^2+hXY^j+zY^{2j},\\
 P_j&=xy^6+
 \sum_{\substack{0\leq i+k\leq7\\ji+k=j+6\\(i,k)\ne(1,6)}}
 c_{ik}x^iy^k.
\end{aligned}
\tag{4.3}
\]
The correction supports are respectively
\[
\begin{array}{c|l}
j& (i,k)\\ \hline
2&(2,4),(3,2),(4,0)\\
3&(2,3),(3,0)\\
4&(2,2)\\
5&(2,1)\\
6&(2,0).
\end{array}
\tag{4.4}
\]
Exact moments through order seven give
\[
 \sqrt{(\text{coefficients of }W_j^m(P_j^m),\ 1\leq m\leq7)}
 =(h,z,c_{ik})
 \tag{4.5}
\]
for every row of (4.4).

At the final common threshold use \(w(x)=6,w(y)=1\) and \(W=12\).  The sole
equality pair is
\[
 (X^2,xy^6).
 \tag{4.6}
\]
Its \(x\)-derivative demand is \(2m\), while its polynomial supply is \(m\).
The weight-defect lemma absorbs the bounded strict selections introduced by
\(Q\), proving (1.2) on the non-pure chart.

## 5. The endpoint \(P_7=y^7\)

This is the only branching chart.  Its first two complete faces are
\[
\begin{aligned}
 W_{3/2}&=X^2+vY^3,
 &P_{3/2}&=y^7+px^2y^4+qx^4y,\\
 W_2&=X^2+BXY^2+AY^4,
 &P_2&=y^7+zxy^5+qx^2y^3+\rho x^3y.
\end{aligned}
\tag{5.1}
\]
Their exact radicals are
\[
 \sqrt{I_{3/2}}=(v,p,q),
 \qquad
 \sqrt{I_2}=(\rho,q,A,Bz).
 \tag{5.2}
\]
The second equality already follows from the first two moments.  It leaves a
\(B\)-axis, a \(z\)-axis, and their intersection.

### 5.1 The \(B\)-axis

Normalize \(B=1\).  The three integral faces are
\[
\begin{array}{c|c|c}
\text{weight}&W&P\\ \hline
3&XY^2+AY^5&y^7+zxy^4+qx^2y\\
4&XY^2+AY^6&y^7+zxy^3\\
5&XY^2+AY^7&y^7+zxy^2.
\end{array}
\tag{5.3}
\]
Every displayed face has origin radical.  The three intervening polynomial
crossings \(cx^3\), \(cx^2y^2\), and \(cx^2\) have pivots \(10080c\),
\(4c\), and \(3360c\), in moments two, one, and two respectively.  The final pair
\[
 (XY^2,y^7)
 \tag{5.4}
\]
has common threshold seven for \(w=(5,1)\) and has positive linear
\(x\)-derivative demand against zero \(x\)-supply.

### 5.2 The \(z\)-axis

Normalize \(z=1\).  For \(j=3,4,5\), the complete first-output face is
\[
\begin{aligned}
 W_j&=X^2+uXY^j+vY^{j+2},\\
 P_j&=y^7+xy^5+cx^2y^{5-j}.
\end{aligned}
\tag{5.5}
\]
There is also one top half-integral face
\[
 (X^2+aY^5,\,xy^5+cx^3),
\tag{5.6}
\]
with radical \((a,c)\).  The three ordinary top faces are
\[
 W=X^2+uXY^j+vY^{2j},
 \qquad P=xy^5+cx^2y^{5-j},
 \qquad j=3,4,5,
 \tag{5.7}
\]
and again have radicals \((u,v,c)\).  These top faces and (5.5) are
different output layers: the latter retain the lower \(y^7\) channel and
the pure operator correction \(Y^{j+2}\).  The unpaired \(Y^8\) and \(Y^9\)
crossings have second-moment pivots \(7257600a\) and \(14515200a\).  Thus
only \((X^2,xy^5)\) remains.  It has common
threshold ten for \(w=(5,1)\), and its \(x\)-demand \(2m\) exceeds its
supply \(m\).

### 5.3 The intersection

When \(B=z=0\), the slope-\(5/2\) face
\[
 W=X^2+AY^5,
 \qquad P=y^7+qx^2y^2
 \tag{5.8}
\]
has radical \((A,q)\).  At slope three,
\[
\begin{aligned}
 W&=X^2+BXY^3+AY^6,\\
 P&=y^7+zxy^4+qx^2y,
\end{aligned}
\qquad
 \sqrt I=(q,A,Bz).
\tag{5.9}
\]
Unlike the sextic row, the two axes in (5.9) are not yet at a common
threshold; each migrates once more.

On the \(B\)-axis, the only half crossing is killed by the moment-two pivot
\(20160c\), and the final face
\[
 W=XY^3+AY^7,
 \qquad P=y^7+zxy^3
 \tag{5.10}
\]
has radical \((A,z)\).  The surviving pair \((XY^3,y^7)\) has common weight
seven for \(w=(4,1)\) and an \(x\)-deficit.

On the \(z\)-axis, both
\[
\begin{array}{c|c}
W&P\\ \hline
X^2+uXY^4+vY^7&y^7+xy^4+qx^2\\
X^2+uXY^4+vY^8&xy^4+qx^2
\end{array}
\tag{5.11}
\]
have radical \((u,v,q)\).  The last pair \((X^2,xy^4)\) has common weight
eight for \(w=(4,1)\), with demand \(2m\) and supply \(m\).

Equation (5.4), the \(j=5\) member of (5.7), (5.10), and the second row of
(5.11) are final common-threshold coordinate-deficit faces.  The
weight-defect lemma handles all strict lower
polynomial terms, all strict higher normalized operator jets, and the fixed
multiplier \(Q\).  This closes the pure endpoint.

## 6. Completion and the remaining frontier

Every nonzero binary quadratic has root partition \((2)\) or \((1+1)\).
Section 2 exhausts their Hall loci.  Section 3 closes the distinct-root
partition; Section 4 closes the non-pure double-line chart; and Section 5
closes the pure endpoint.  This proves Theorem 1.1.

The computation also identifies where the new degree really entered.  The
slope-three intersection radical in (5.9) is the same two-axis pattern as in
degree six, but its axes are one unit short of a common threshold.  The extra
migration (5.10)--(5.11), not an adelic collision, is the first new septic
shell.  Both branches terminate.  Thus the next degree-seven attack should
start at lowest order three, rather than enlarge projected carry-signature
censuses on the already closed quadratic row.

## 7. Exact replay

Run

```bash
.venv/bin/python scripts/verify_binary_quadratic_septic_gvc.py
```

The command uses SymPy exact sparse contraction and Singular exact radicals
over \(\mathbb Q\).  It checks the seven-factor Hall classification, all
twenty-one distinct-root eliminations and the ten-row second-moment ladder,
the complete half-integral and integral support ledgers, every displayed
crossing coefficient, and all radicals in the repeated-root branch tree.
The all-order implication is the written Hall, local-division, and final
weight-defect argument, not a bounded extrapolation from the checker.
