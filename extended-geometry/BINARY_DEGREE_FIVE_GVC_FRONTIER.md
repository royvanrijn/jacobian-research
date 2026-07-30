# Binary degree-five GVC frontier

## 1. Statement and scope

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_r+\Lambda_{r+1}+\cdots
\]
be a nonzero constant-coefficient operator in two variables with no
order-zero term, where \(\Lambda_r\ne0\), and let \(\deg P=5\).

This note closes the frontier requested after the quadratic-leading row:

> **Theorem 1.1 — cubic-leading quintic row.**  If \(r=3\), then
> \[
>  \Lambda^m(P^m)=0\quad(m\ge1)
> \]
> implies
> \[
>  \Lambda^m(QP^m)=0\quad(m\gg0)
> \]
> for every fixed polynomial \(Q\).  All eight leading normal forms in
> Proposition 3.17 of
> [the escape-obstruction note](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
> close.

> **Theorem 1.2 — squarefree quartic-leading quintic row.**  Suppose
> \(r=4\) and the binary quartic symbol of \(\Lambda_4\) is squarefree.
> The first four pure equations already imply the same mixed conclusion.
> The result is uniform in the quartic cross-ratio.

No binary counterexample occurs on these rows.  Together with the already
proved \(r=2\) row, this removes every case explicitly named in the requested
degree-five frontier.

There is one status caveat.  The repository did not previously contain a
theorem for the non-squarefree \(r=4,\deg P=5\) **nonhomogeneous** rows.
The balanced low-root and split-symbol theorems do not cover arbitrary
higher operator jets.  Thus Theorems 1.1–1.2 alone do not justify the
unqualified statement that every binary operator satisfies GVC through
degree five.  That universal corollary requires those omitted discrete
order-four partitions as an additional dependency.

## 2. Weight-defect lemma

The all-order step used below is elementary but its common threshold is
essential.

Let \(w=(w_x,w_y)\) be a positive integral weight and suppose that, for one
integer \(W\),
\[
 w(x^iy^j)\le W
 \quad\text{for every monomial of }P,
 \qquad
 w(\partial_x^a\partial_y^b)\ge W
 \quad\text{for every monomial of }\Lambda.
\tag{2.1}
\]
For a monomial term selected in the expansion of
\(\Lambda^m(QP^m)\), put
\[
\begin{aligned}
 \delta_\Lambda&=\sum_{\text{operator selections}}
   (w(\partial_x^a\partial_y^b)-W),\\
 \delta_P&=\sum_{\text{polynomial selections}}
   (W-w(x^iy^j)).
\end{aligned}
\tag{2.2}
\]
Its output weight is
\[
 w(\text{output})=w(Q_{\rm mon})-\delta_\Lambda-\delta_P.
\tag{2.3}
\]
Since output exponents are nonnegative,
\[
 \delta_\Lambda+\delta_P\le w(Q_{\rm mon})\le w(Q).
\tag{2.4}
\]
Every strict selection has positive integral defect.  Hence at most \(w(Q)\)
operator or polynomial factors can lie off the equality faces.

If equality is strict on one side of (2.1), (2.3) immediately gives
vanishing for \(m>w(Q)\).  Otherwise every equality-face pair below has an
\(x\)- or \(y\)-derivative deficit linear in \(m\); only \(O_Q(1)\) factors
can be exceptional by (2.4), so the mixed expression vanishes eventually.
A coarse bound used on the equality components is
\[
 m>\deg_xQ+6w(Q).
\tag{2.5}
\]

The intermediate weights below extract successive coefficients.  They are
not all common-threshold separators.  The lemma is invoked only after the
last elimination, with the final weights displayed explicitly.

## 3. Triple-root cubic symbol

### 3.1 Exact Weierstrass normalization

Normalize the leading symbol to \(X^3\).  Formal Weierstrass division gives
\[
 \lambda(X,Y)=U(X,Y)W(X,Y),
\qquad
 W=X^3+a(Y)X^2+b(Y)X+c(Y),
\tag{3.1}
\]
where
\[
 U(0,0)=1,\qquad
 \operatorname{ord}a\ge2,\quad
 \operatorname{ord}b\ge3,\quad
 \operatorname{ord}c\ge4.
\tag{3.2}
\]
The formal constant-coefficient operators \(U(\partial)\) and
\(U(\partial)^{-1}\) act finitely on each polynomial.  Since all
constant-coefficient operators commute,
\[
 \Lambda^m(F)=U(\partial)^mW(\partial)^m(F).
\tag{3.3}
\]
Thus \(\Lambda^m(F)=0\) if and only if \(W(\partial)^m(F)=0\), for
\(F=P^m\) and \(F=QP^m\).  We may work with \(W\) exactly, not merely with
a truncated jet.

### 3.2 The \(x^2y^3\) and \(x(x-y)y^3\) forms

Use \(w=(3,1)\) and \(W=9\).  The equality operator lattice is
\[
 (3,0),\ (2,3),\ (1,6),\ (0,9),
\tag{3.4}
\]
and the equality polynomial lattice is
\[
 (2,3),\ (3,0).
\tag{3.5}
\]
Write
\[
\begin{aligned}
 P_{[9]}&=x^2y^3+Rx^4+Sx^3y+Tx^3,\\
 W_{[9]}&=\partial_x^3
 +a_2\partial_x^2\partial_y^2
 +a_3\partial_x^2\partial_y^3\\
&\quad
 b_3\partial_x\partial_y^3+\cdots
 +b_6\partial_x\partial_y^6
 +c_4\partial_y^4+\cdots+c_9\partial_y^9.
\end{aligned}
\tag{3.6}
\]
The first equation gives
\[
 R=-\frac{b_3}{2},\qquad S=-2a_2,\qquad T=-2a_3.
\tag{3.7}
\]
The second moment successively gives
\[
\begin{gathered}
 c_4=b_3=b_4=c_5=c_6=0,\\
 b_5=-\frac{a_2a_3}{3},\qquad
 b_6=-\frac{5a_3^2}{12}.
\end{gathered}
\tag{3.8}
\]
The decisive third-moment coefficient is
\[
 [y^3]\,W^3(P^3)=-6531840a_2^3.
\tag{3.9}
\]
After \(a_2=0\), the remaining coefficients give
\[
 c_7=c_8=0,\qquad c_9=\frac{1169}{2160}a_3^3.
\tag{3.10}
\]
Finally
\[
 W^4(P^4)=67315784417280a_3^4.
\tag{3.11}
\]
The equality pair is therefore
\[
 (\partial_x^3,x^2y^3),
\tag{3.12}
\]
which has a linear \(x\)-derivative deficit.

For the second normal form, the shear \(x\mapsto x+y/2\) preserves
\(\partial_x^3\) and sends
\[
 x(x-y)y^3
 \longmapsto
 x^2y^3-\frac14y^5.
\tag{3.13}
\]
The added term has weight five, four below the equality face.  It cannot
alter (3.7)–(3.11), and (2.4) absorbs it in the mixed expression.

### 3.3 The \(xy^4\) form

Use \(w=(2,1)\) and \(W=6\).  The complete equality pair is
\[
\begin{aligned}
 P&=xy^4+ux^2y^2+vx^3,\\
 W&=\partial_x^3+A\partial_x^2\partial_y^2
  +B\partial_x\partial_y^4+C\partial_y^6.
\end{aligned}
\tag{3.14}
\]
The preceding transverse coefficients are killed by
\[
 [x^2]M_2=54720c_4^2,\qquad
 [y^2]M_2=28800b_3^2,\qquad
 [y^2]M_3=4311014400c_5^2.
\tag{3.15}
\]
Moment one on (3.14) gives
\[
 v=-\frac23Au-4B.
\tag{3.16}
\]
Let \(\mu_m\) be the scalar \(m\)-th moment after (3.16).  Exact
characteristic-zero elimination gives
\[
 \sqrt{(\mu_2,\mu_3,\mu_4,\mu_5)}=(B,C,Au)
 =(u,B,C)\cap(A,B,C).
\tag{3.17}
\]
The two components are
\[
\begin{array}{c|c}
 u=0&
 \bigl(\partial_x^2(\partial_x+A\partial_y^2),xy^4\bigr),\\
 A=0&
 \bigl(\partial_x^3,xy^4+ux^2y^2\bigr).
\end{array}
\tag{3.18}
\]
Both have a linear \(x\)-derivative deficit, so (2.4) proves the mixed
conclusion.

### 3.4 The \(y^5\) form

The \(w=(2,1)\) ladder first gives
\[
 [y]M_2=576Tb(2T+35b),
\qquad
 [x]M_3=311040T^2b^2(5T+84b).
\tag{3.19}
\]
The two nonzero ratios conflict, hence \(Tb=0\).

If \(T\ne0\), then \(b=0\), the next equation gives \(c_5=0\), and the
equality pair is the already proved quartic face
\[
\begin{aligned}
 P&=x^2y^2-\frac23Ax^3,\\
 W&=\partial_x^3+A\partial_x^2\partial_y^2
 -\frac29A^2\partial_x\partial_y^4
 +\frac{136}{405}A^3\partial_y^6.
\end{aligned}
\tag{3.20}
\]
Its fourth moment is
\[
 3361505280A^4.
\tag{3.21}
\]

If \(b\ne0\), then \(T=0\) and the next two equations give
\[
 c_5=-\frac{Zb}{20},\qquad H=-\frac{17Z^2}{20}.
\tag{3.22}
\]
The remaining face
\[
\begin{aligned}
 P&=y^5+Zxy^3-\frac{17}{20}Z^2x^2y,\\
 W&=b\left(\partial_x\partial_y^3
 -\frac Z{20}\partial_y^5\right)
\end{aligned}
\tag{3.23}
\]
has
\[
 W^3(P^3)=-37666944Z^3b^3.
\tag{3.24}
\]
Thus it reduces to the one-sided pair
\((b\partial_x\partial_y^3,y^5)\).

Finally, if \(T=b=0\), the residual equation is
\[
 3859200c_5^2=0.
\tag{3.25}
\]
All triple-root branches are now one-sided or strict.

## 4. Double-root and squarefree cubic symbols

The double-root representatives are
\[
 (\partial_x^2\partial_y,x^5),\quad
 (\partial_x^2\partial_y,xy^4),\quad
 (\partial_x^2\partial_y,y^5).
\tag{4.1}
\]
For the squarefree representative
\[
 \partial_x\partial_y(\partial_x+\partial_y)
 =\partial_x^2\partial_y+\partial_x\partial_y^2,
\tag{4.2}
\]
the top polynomial is \(x^5\).  The second summand in (4.2) is strictly
above every \(x^5\) extraction face below, so that row has the same
calculation.

### 4.1 The \(x^5\) row

Four exhaustive equality lattices occur:
\[
\begin{array}{c|c|c}
w&\Lambda_{\rm face}&P_{\rm face}\\ \hline
(3,5)&\partial_x^2\partial_y&x^5+Dy^3\\
(1,2)&\partial_x^2\partial_y+U\partial_x^4&
x^5+Ax^3y+Bxy^2\\
(2,5)&\partial_x^2\partial_y&x^5+Ey^2\\
(1,3)&\partial_x^2\partial_y+T\partial_x^5&
x^5+Cx^2y.
\end{array}
\tag{4.3}
\]
Their decisive equations are
\[
\begin{aligned}
 M_2(x^5+Dy^3)&=1440Dxy,\\
 M_1&=6x(A+20U),\\
 M_2|_{A=-20U}
 &=144\left(10(B+340U^2)x^2-3200U^3y\right),\\
 M_2(x^5+Ey^2)&=480Ex,\\
 M_1(x^5+Cx^2y)&=2(C+60T),\\
 M_2|_{C=-60T}&=2592000T^2.
\end{aligned}
\tag{4.4}
\]
After these eliminations, \(w=(3,10)\), \(W=15\) is a strict final
separator:
\[
 w(P)\le15,\qquad w(\Lambda)\ge16.
\tag{4.5}
\]
Hence the mixed expression vanishes for \(m>w(Q)\).  In the squarefree
case, \(\partial_x\partial_y^2\) has weight \(23\), so (4.5) is unchanged.

### 4.2 The \(xy^4\) row

The three extraction faces use weights \((2,1),(5,2),(3,1)\).  The first is
\[
\begin{aligned}
 \Lambda&=\partial_x^2\partial_y
 +H\partial_x\partial_y^3+J\partial_y^5,\\
 P&=xy^4-6Hx^2y^2+Bx^3.
\end{aligned}
\tag{4.6}
\]
The second moment gives
\[
 B-2H^2+140J=0,\qquad
 H(B+18H^2-100J)=0.
\tag{4.7}
\]
The two nonzero branches are killed by the third-moment coefficients
\[
 71884800H^3y^3,\qquad
 9638092800J^2xy.
\tag{4.8}
\]
The next face has
\[
 M_2(xy^4;\partial_x^2\partial_y+M\partial_y^6)
 =161280My.
\tag{4.9}
\]
The last face is
\[
\begin{aligned}
 \Lambda&=\partial_x^2\partial_y
 +K\partial_x\partial_y^4+L\partial_y^7,\\
 P&=xy^4-12Kx^2y,
\end{aligned}
\tag{4.10}
\]
where moment two gives \(L=-23K^2/70\), and moment three is
\[
 -3318921216K^3.
\tag{4.11}
\]
The strict final separator is \(w=(7,2)\), \(W=15\):
\[
 w(P)\le15,\qquad w(\Lambda)\ge16.
\tag{4.12}
\]

### 4.3 The \(y^5\) row

The first two faces kill the \(\partial_y^4\) coefficient and the \(x^3\)
polynomial correction:
\[
\begin{aligned}
 M_2(y^5-30E x^2y^2;
   \partial_x^2\partial_y+E\partial_y^4)
 &=1468800E^2y^2,\\
 M_3(y^5+Dx^3;\partial_x^2\partial_y)
 &=129600D^2y^2.
\end{aligned}
\tag{4.13}
\]
The final \(w=(2,1)\), \(W=5\) face is
\[
\begin{aligned}
 \Lambda&=\partial_x^2\partial_y
 +H\partial_x\partial_y^3+J\partial_y^5,\\
 P&=y^5+Axy^3+Bx^2y,\qquad B=-3AH-60J.
\end{aligned}
\tag{4.14}
\]
Exact moments through order four give
\[
 \sqrt{(M_1,M_2,M_3,M_4)}=(J,B,AH).
\tag{4.15}
\]
The two survivor lines are
\[
\begin{array}{c|c}
H=J=B=0&P=y^5+Axy^3,\\
A=J=B=0&\Lambda=\partial_x^2\partial_y
 +H\partial_x\partial_y^3.
\end{array}
\tag{4.16}
\]
Both have a one-unit \(x\)-degree deficit, so (2.4) applies.

For the requested mixed test, take \(Q=x^2\).  On the two lines, the first
two values are respectively
\[
\begin{aligned}
 &2y^2(9Ax+5y^2),\qquad 720A^2y^4,\\
 &10y^2(12Hx+y^2),\qquad 302400H^2y^4,
\end{aligned}
\tag{4.17}
\]
and every value from \(m=3\) onward is zero.  Thus neither stable component
is a counterexample.

## 5. Squarefree quartic-leading row

Normalize
\[
 D_4=\partial_x\partial_y(\partial_x-\partial_y)
 (\partial_x-\lambda\partial_y),
\qquad \lambda\ne0,1.
\tag{5.1}
\]
For
\[
 P_5=\sum_{i=0}^5a_i x^{5-i}y^i,
\]
the first equation gives
\[
\begin{aligned}
 a_1&=\frac{(\lambda+1)a_2-\lambda a_3}{2},\\
 a_4&=\frac{(\lambda+1)a_3-a_2}{2\lambda}.
\end{aligned}
\tag{5.2}
\]
After (5.2), exact projective saturation over
\(\mathbb Q[\lambda,1/(\lambda(\lambda-1))]\) of the coefficient ideal from
\[
 D_4^2(P_5^2),\qquad D_4^3(P_5^3)
\tag{5.3}
\]
is the union of the four sections
\[
 P_5=x^5,\quad y^5,\quad(x+y)^5,\quad(\lambda x+y)^5.
\tag{5.4}
\]
There are no exceptional harmonic or equianharmonic cross-ratios.  Root
permutation reduces the correction calculation to \(P_5=x^5\).

Use translation to kill \([x^4]P_4\).  A differential unit
\(1+\alpha\partial_x+\beta\partial_y\) kills the
\(\partial_x^4\partial_y\) and
\(\partial_x^3\partial_y^2\) coefficients of \(\Lambda_5\).  The
degree-\((m-1)\) pieces of moments \(m=4,3,2\) successively give
\[
 [y^4]P_4=[xy^3]P_4=[x^2y^2]P_4=0.
\tag{5.5}
\]
Thus
\[
\begin{aligned}
 P_4&=t x^3y,\\
 \Lambda_5&=a\partial_x^5
 +u\partial_x^2\partial_y^3
 +v\partial_x\partial_y^4+w\partial_y^5,
\end{aligned}
\tag{5.6}
\]
and moment one gives
\[
 a=-\frac t{20}.
\tag{5.7}
\]

With \(w=(3,5)\), the maximum polynomial face is
\[
 x^5+c\,y^3
\]
and the unique minimum operator face is
\(\partial_x^3\partial_y\).  The extremal third moment is
\[
 (\partial_x^3\partial_y)^3(x^5+cy^3)^3
 =65318400c\,x.
\tag{5.8}
\]
Hence \(c=0\).

Now use \(w=(1,2)\), \(W=5\).  The complete equality pair is
\[
\begin{aligned}
 L&=\partial_x^3\partial_y+a\partial_x^5,\\
 F&=x^5+t x^3y+sxy^2.
\end{aligned}
\tag{5.9}
\]
Its first three moments are
\[
\begin{aligned}
 M_1={}&6(20a+t),\\
 M_2={}&1440(2520a^2+112at+t^2+2s),\\
 M_3={}&2177280(600600a^3+25740a^2t+330at^2\\
 &\hspace{27mm}+330as+t^3+6ts).
\end{aligned}
\tag{5.10}
\]
Equations \(M_1=M_2=0\) give
\[
 a=-\frac t{20},\qquad s=-\frac{17}{20}t^2,
\tag{5.11}
\]
and then
\[
 M_3=-37666944t^3.
\tag{5.12}
\]
Thus \(a=t=s=0\).

At the terminal pair, equality in (2.1) occurs only at
\[
 (\partial_x^3\partial_y,x^5).
\tag{5.13}
\]
If \(K=w(Q)\), a contributing term has at most \(K\) off-face selections.
There are at least \(m-K\) equality \(\partial_y\)-derivatives, while its
input \(y\)-degree is at most
\(\deg_yQ+5K\).  Hence
\[
 \Lambda^m(QP^m)=0
 \qquad
 \text{for }m>\deg_yQ+6K.
\tag{5.14}
\]
For example, \(Q=y\) gives \(60x^2\) at \(m=1\) and zero for all \(m\ge2\).

Only the first four full pure equations were used: moments one through
three classify \(P_5\), moment four starts the triangular \(P_4\)
elimination, and the terminal faces use moments at most three.

## 6. Finite-field searches

The modular search is deliberately separate from the proof.

On the double-root \(xy^4\) face, after the first equation the parameters
are \((B,H,J)\).  Exhaustive enumeration at \(p=101,103,107\) gives
\[
\begin{array}{c|c|c}
p&\text{after }M_2&\text{after }M_3\\ \hline
101&201&1\\
103&205&1\\
107&213&1.
\end{array}
\tag{6.1}
\]
The single point is the origin.

On the \(y^5\) face, enumerate \((A,H,J)\) after
\(B=-3AH-60J\):
\[
\begin{array}{c|c|c|c}
p&\text{after }M_2&\text{after }M_3&\text{after }M_4\\ \hline
101&10201&201&201\\
103&9589&205&205\\
107&12721&425&213.
\end{array}
\tag{6.2}
\]
The final set is exactly
\[
 \{(A,0,0)\}\cup\{(0,H,0)\},
\tag{6.3}
\]
of size \(2p-1\).  The residual stabilizer acts by
\[
 (A,H,B,J)\longmapsto(tA,tH,t^2B,t^2J),
\tag{6.4}
\]
so the nonzero survivor set has two projective orbits.

These two face searches enumerate \(6{,}696{,}142\) raw triples.

For the quartic row, every \(\lambda\ne0,1\) and every projective point of
the four-dimensional first-moment kernel was checked at
\[
 p=11,13,17,19,23,29,31.
\tag{6.5}
\]
This is \(2{,}082{,}612\) projective top-form tuples.  Moments two and three
always leave exactly the four fifth powers in (5.4).

The modular calculations use primitive characteristic-zero coefficient
equations before reduction.  They are nullcone screens and regression
tests, not substitutes for the exact radical and face-separation proofs.

## 7. Reproduction

Run the exact moment and mixed-tail checker:

```bash
.venv/bin/python scripts/verify_binary_degree_five_gvc_frontier.py
```

If Singular is available, also verify the two residual face radicals:

```bash
.venv/bin/python scripts/verify_binary_degree_five_gvc_frontier.py --singular
```

The uniform squarefree-quartic saturation is the expensive exact command:

```bash
.venv/bin/python scripts/verify_binary_degree_five_gvc_frontier.py \
  --singular-top
```

Reproduce the \(8{,}778{,}754\) finite-field tuples with:

```bash
.venv/bin/python scripts/search_binary_degree_five_gvc_faces_mod_p.py
```

The generated modular search record is
[`binary_degree_five_gvc_face_search.json`](../artifacts/generated-results/binary_degree_five_gvc_face_search.json).

## 8. Why this is not yet the universal degree-five theorem

The squarefree row is only the partition \((1,1,1,1)\) of a quartic
symbol.  The repeated-root partitions have larger leading pure-zero loci.
For
\[
 D_4=\partial_x^a\partial_y^b,\qquad a+b=4,
\]
write \(P_5=x^5p(y/x)\).  The \(m\)-th leading moment tests exactly the
coefficient band
\[
 [t^k]p(t)^m,\qquad bm\le k\le (b+1)m.
\tag{8.1}
\]
Applying the one-variable constant-term theorem after every rational
rescaling shows that the support of \(p\) must lie wholly below or wholly
above the interval \([b,b+1]\).  Consequently the leading loci are
\[
\begin{array}{c|c}
\text{root partition}&\text{leading pure-zero locus}\\ \hline
(4)&P_5=y^2C_3(x,y),\\
(3,1)&P_5=x^5\ \text{or}\ y^3C_2(x,y),\\
(2,2)&P_5=x^4L(x,y)\ \text{or}\ y^4L(x,y).
\end{array}
\tag{8.2}
\]
For the remaining partition, normalize
\[
 D_4=\partial_x^2\partial_y(\partial_x-\partial_y).
\]
Exact moments one through four have radical
\[
 P_5=x^5,\qquad P_5=(x+y)^5,\qquad P_5=y^4L(x,y).
\tag{8.3}
\]
Thus the omitted rows have positive-dimensional leading families; they
cannot be inferred from the squarefree calculation by specialization.

The higher-jet problem is finite but still nontrivial.  For example, in
the \((4)\) row, Weierstrass normalization followed by the weight
\(w(x)=2,w(y)=1\), \(W=8\), leaves the below-threshold operator terms
\[
 X^2Y^3,\ XY^4,\ XY^5,\ Y^5,\ Y^6,\ Y^7,
\tag{8.4}
\]
and the equality operator chain
\[
 X^4,\ X^3Y^2,\ X^2Y^4,\ XY^6,\ Y^8.
\tag{8.5}
\]
An arbitrary-jet elimination theorem for all families in (8.2)--(8.3)
is not proved here.  This is the precise additional dependency needed
before upgrading Theorems 1.1--1.2 to “every binary operator through
polynomial degree five.”
