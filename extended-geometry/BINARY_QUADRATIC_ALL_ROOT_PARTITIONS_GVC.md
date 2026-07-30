# Binary sextic GVC for every quadratic leading symbol

## 1. Theorem

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda _2+\Lambda _3+\cdots
\]
be a constant-coefficient operator in two variables with lowest positive
order two, and let \(\deg P=6\).

> **Theorem 1.1 — complete quadratic-leading sextic row.**
> If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Arbitrary lower pieces of \(P\), arbitrary higher operator jets, and both
root partitions \((2)\) and \((1+1)\) of the binary quadratic are included.
Together with the already proved rows of lowest orders three, four, and
five, and the theorem for \(\deg P\le r\), this gives:

> **Corollary 1.2 — binary GVC through polynomial degree six.**
> Every constant-coefficient differential operator in two variables
> satisfies GVC for every polynomial of degree at most six.

The all-order premise in (1.1) is used first for the leading Hall
classification.  After that reduction, the distinct-root chart needs only
the first two moments.  The repeated-root Newton faces use moments through
eight.

## 2. Hall classification and local division

After scalar extension, factor
\[
 \Lambda _2=D_{v_1}D_{v_2},\qquad
 P_6=\prod_{j=1}^6L_j.
\tag{2.1}
\]
The translated split-symbol construction and the
Duistermaat--van der Kallen theorem show that the two derivative copies
cannot be matched to two distinct polynomial factors on which they act
nontrivially.

If a root direction of \(\Lambda _2\) has multiplicity \(e\), and \(c\)
factors of \(P_6\) annihilate it, Hall failure is
\[
 6-c<e,\qquad\text{equivalently}\qquad c\geq7-e.
\tag{2.2}
\]
A subset containing two nonparallel derivative directions sees all six
polynomial factors, so there is no other deficient subset.  The two local
models are therefore
\[
\begin{array}{c|c|c}
e&\Lambda _2&P_6\\ \hline
1&XY&y^6,\\
2&X^2&y^5(Ax+Cy).
\end{array}
\tag{2.3}
\]
On \(A\ne0\), a shear preserving \(X^2\) removes \(C\), leaving \(xy^5\).
The boundary \(A=0\) is \(y^6\).

Formal division by the leading symbol removes irrelevant multiples.  Thus
\[
\begin{aligned}
 XY\text{ chart}:&\quad
 W=XY+H(X)+K(Y),\\
 X^2\text{ chart}:&\quad
 W=X^2+a(Y)X+b(Y),
\end{aligned}
\tag{2.4}
\]
where the displayed series start in total order three.  The removed unit
and its inverse act locally finitely on polynomials, so this changes neither
(1.1) nor (1.2).

Every Newton face used below is complete for (2.4).  The
[common weight-defect lemma](BINARY_DEGREE_FIVE_GVC_FRONTIER.md#2-weight-defect-lemma)
absorbs all terms strict on the final face.  Intermediate faces only extract
the next coefficient; they are not asserted to be final separators.

## 3. The distinct-root chart

Swap the variables so that \(P_6=x^6\).  Write
\[
 W=XY+\sum_{j\ge3}a_jX^j+\sum_{j\ge3}b_jY^j.
\tag{3.1}
\]
Because the \(XY\) coefficient is a unit, \(W(P)=0\) solves every mixed
coefficient of \(P\) triangularly, from high total degree down.  Retain all
operator terms through order twelve, since later terms kill \(P^2\).

Put \(p_{0j}=[y^j]P\).  Successive coefficients of the full second moment
are
\[
\begin{aligned}
[x^4y^3]W^2(P^2)&=1200p_{05},\\
[x^4y^2]W^2(P^2)&=720p_{04},\\
[x^6]W^2(P^2)&=201600a_3^2,\\
[x^4y]W^2(P^2)&=360p_{03},\\
[x^4]W^2(P^2)&=120(113040a_4^2+p_{02}),\\
[xy]W^2(P^2)&=-2880a_4(15840a_4^2+p_{02}),\\
[x^2]W^2(P^2)&=213580800a_5^2,\\
[1]W^2(P^2)&=466560000a_6^2.
\end{aligned}
\tag{3.2}
\]
Each line is read after the preceding parameters have vanished.  The fifth
line solves \(p_{02}\); the sixth then forces \(a_4=0\).  Hence
\[
 P=f(x)+cy,\qquad
 W=Y\Gamma(X,Y)+H(X),\qquad \operatorname{ord}H\ge7>\deg f.
\tag{3.3}
\]
The transverse-linear/high-order theorem in
[the separable escape note](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
now proves (1.2).  This closes the whole distinct-root partition, not only
the displayed finite jet.

## 4. The double line over \(P_6=xy^5\)

Use \(W=X^2+a(Y)X+b(Y)\).  The half-integral Newton slopes
\[
 \frac32,\ \frac52
\tag{4.1}
\]
have complete faces
\[
\begin{array}{c|c}
W&P\\ \hline
X^2+aY^3&xy^5+cx^3y^2,\\
X^2+aY^5&xy^5+cx^3,
\end{array}
\tag{4.2}
\]
and both exact radicals are \((a,c)\).  The \(Y^7\) and \(Y^9\) crossings
have extremal second-moment coefficients
\[
 2419200a,\qquad 14515200a,
\tag{4.3}
\]
so they vanish as well.

For each integral slope \(j=2,3,4,5\), the complete equality face is
\[
\begin{aligned}
 W_j&=X^2+hXY^j+zY^{2j},\\
 P_j&=xy^5+
 \sum_{\substack{0\le i+k\le6\\ji+k=j+5\\(i,k)\ne(1,5)}}
 c_{ik}x^iy^k.
\end{aligned}
\tag{4.4}
\]
Exact moments through order six give
\[
 \sqrt{(\text{coefficients of }W_j^m(P_j^m),\,1\le m\le6)}
 =(h,z,c_{ik})
\tag{4.5}
\]
for every \(j=2,3,4,5\).  For orientation, after the second moment solves
the pure-\(Y\) companion, the last three terminal coefficients are
\[
\begin{array}{c|c}
j&\text{nonzero terminal coefficient}\\ \hline
3&-3888000h^3,\\
4&-2193454080000h^3,\\
5&-15435360000000h^3.
\end{array}
\tag{4.6}
\]

After (4.2)--(4.6), weights \(w(x)=5,w(y)=1\) give the common threshold
ten.  The equality pair is only
\[
 (X^2,xy^5).
\tag{4.7}
\]
It has \(x\)-derivative demand \(2m\) and polynomial \(x\)-supply \(m\).
The common weight-defect lemma absorbs the bounded number of strict
selections introduced by \(Q\).  Thus (1.2) holds on the non-pure
double-line chart.

## 5. The pure-sixth-power endpoint

The endpoint \(P_6=y^6\) contains the only nested branching.  The complete
slope-\(3/2\) face is
\[
 W=X^2+vY^3,\qquad
 P=y^6+px^2y^3+qx^4,
\tag{5.1}
\]
and its radical through moment five is \((v,p,q)\).

At slope two the complete face is
\[
\begin{aligned}
W&=X^2+BXY^2+AY^4,\\
P&=y^6+zxy^4+qx^2y^2+rx^3.
\end{aligned}
\tag{5.2}
\]
Moments through order eight give
\[
 \sqrt I=(r,q,A,Bz).
\tag{5.3}
\]
Thus there are two axes and their intersection.

### 5.1 The \(B\)-axis

Normalize \(B=1\).  The first secondary face is
\[
 W=XY^2+AY^5,\qquad
 P=y^6+zxy^3+qx^2,
\tag{5.4}
\]
and has radical \((A,z,q)\).  The last possible migration is the
common-threshold weight-four face
\[
 W=XY^2+AY^6,\qquad P=y^6+zxy^2,
\tag{5.5}
\]
whose radical is \((A,z)\).  The surviving equality pair
\[
 (XY^2,y^6)
\tag{5.6}
\]
has common weight six and a linear \(x\)-derivative deficit.

### 5.2 The \(z\)-axis

Normalize \(z=1\).  The first secondary face is
\[
\begin{aligned}
W&=X^2+uXY^3+vY^5,\\
P&=y^6+xy^4+px^2y,
\end{aligned}
\tag{5.7}
\]
and has radical \((u,v,p)\).  The final weight-four face is
\[
\begin{gathered}
 [1](X^2+aY^6)(y^6+xy^4)=720a,\\
 [y](X^2+aY^7)^2(y^6+xy^4)^2=161280a,
\end{gathered}
\tag{5.8}
\]
so the intervening pure-\(Y\) crossings vanish.  The final weight-four
face is
\[
 W=X^2+uXY^4+vY^8,\qquad P=xy^4+qx^2,
\tag{5.9}
\]
with radical \((u,v,q)\).  Its surviving equality pair is
\[
 (X^2,xy^4),
\tag{5.10}
\]
at common weight eight, again with a linear \(x\)-derivative deficit.
The term \(y^6\) is strict below this face.

### 5.3 The intersection

When \(B=z=0\), the slope-\(5/2\) face
\[
 W=X^2+AY^5,\qquad P=y^6+qx^2y
\tag{5.11}
\]
has radical \((A,q)\).  The final slope-three face is
\[
\begin{aligned}
W&=X^2+BXY^3+AY^6,\\
P&=y^6+zxy^3+qx^2.
\end{aligned}
\tag{5.12}
\]
Its exact radical is
\[
 (q,A,Bz).
\tag{5.13}
\]
The two components are
\[
 (X^2+BXY^3,y^6),\qquad
 (X^2,y^6+zxy^3).
\tag{5.14}
\]
Both sides have common weight six.  On the first component every operator
selection has positive \(x\)-derivative demand and the polynomial face has
none; on the second, the demand is \(2m\) and the supply is at most \(m\).

Equations (5.6), (5.10), and (5.14) are therefore final
common-threshold coordinate-deficit faces.  The weight-defect lemma covers
all strict lower polynomial pieces, all strict higher normalized operator
jets, and the fixed multiplier \(Q\).  This closes the pure-sixth-power
endpoint.

## 6. Completion

Every nonzero binary quadratic has root partition \((2)\) or \((1+1)\).
The Hall locus (2.3) is exhaustive.  Section 3 closes the distinct-root
partition, while Sections 4--5 close both projective charts of the double
line.  This proves Theorem 1.1.

For \(\deg P\le6\), a nonzero order-zero operator term forces \(P=0\);
lowest order one is covered in arbitrary degree by formal drift
straightening; lowest orders two through five are now closed; and lowest
order at least six is covered by the theorem \(\deg P\le r\).  This proves
Corollary 1.2.

## 7. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quadratic_all_root_partitions_gvc.py
```

The checker uses exact sparse contraction in SymPy and exact radicals over
\(\mathbb Q\) in Singular.  It verifies the Hall classification, the full
distinct-root first-equation reduction, every repeated-root primary and
secondary Newton face, and all final coordinate-axis radicals.  The
all-order step is the written Hall, local-division, and weight-defect
argument, not a bounded search.
