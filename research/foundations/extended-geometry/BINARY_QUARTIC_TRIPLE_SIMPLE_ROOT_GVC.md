# Binary quartic-leading GVC with root partition \((3+1)\)

## 1. Theorem and scope

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_4+\Lambda_5+\cdots
\]
be a constant-coefficient operator in two variables whose lowest nonzero
symbol has one triple root and one simple root, and let \(\deg P=5\).

> **Theorem 1.1 — triple-plus-simple quartic-leading quintic GVC.**
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

Thus no genuinely nonhomogeneous binary GVC counterexample occurs on the
quartic root partition \((3+1)\).  Arbitrary lower pieces of \(P\) and
arbitrary higher operator jets are included.

The proof uses the leading pure-zero classification already recorded in
the
[`degree-five frontier`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md):
after scalar extension and normalization of the leading symbol to
\(X^3Y\), the degree-five leading form is either
\[
 P_5=x^5
 \qquad\text{or}\qquad
 P_5=y^3(Cy^2+Axy+Bx^2).
\tag{1.3}
\]

## 2. Differential-unit normal form

Formal monomial division gives
\[
 \lambda(X,Y)=U(X,Y)W(X,Y),\qquad U(0,0)=1,
\tag{2.1}
\]
where
\[
 W=X^3Y+
 \sum_{d\ge5}
 \left(
  a_{d,0}Y^d+a_{d,1}XY^{d-1}
  +a_{d,2}X^2Y^{d-2}+a_{d,d}X^d
 \right).
\tag{2.2}
\]
Indeed, every omitted monomial is divisible by \(X^3Y\) and is absorbed
recursively into \(U\).  The locally finite differential operator
\(U(\partial)\) is invertible on polynomials and commutes with
\(W(\partial)\).  Hence it changes neither (1.1) nor (1.2).

Write
\[
 W_5=\ell_0Y^5+\ell_1XY^4+\ell_2X^2Y^3+\ell_5X^5
\tag{2.3}
\]
and
\[
 P_4=p_0y^4+p_1xy^3+p_2x^2y^2+p_3x^3y+p_4x^4.
\tag{2.4}
\]
The defect of a term is its loss in polynomial degree plus its excess in
operator order.  Terms of different defect land in different ordinary
output degrees, so the equations below are exact coefficient equations
of the full moments.

## 3. The \(y^3C_2\) family

### 3.1 Defect-one radical

For
\[
 P_5=Cy^5+Axy^4+Bx^2y^3,
\tag{3.1}
\]
moment one gives
\[
 p_3=-20C\ell_0-4A\ell_1-2B\ell_2.
\tag{3.2}
\]
After (3.2), the positive-degree parts of moments two and three have,
up to nonzero rational scalars, the equations
\[
 B(14A\ell_0+4B\ell_1+p_4),\qquad
 B^2\ell_0,\qquad B^3\ell_0.
\tag{3.3}
\]
Exact characteristic-zero elimination gives
\[
 \sqrt{J_1}
 =
 (B\ell_0,\;B(4B\ell_1+p_4)).
\tag{3.4}
\]
Its two minimal components are
\[
 (B),\qquad(\ell_0,4B\ell_1+p_4).
\tag{3.5}
\]

### 3.2 The quadratic-tip chart \(B\ne0\)

Normalize \(B=1\).  Equation (3.5) gives
\[
 \ell_0=0,\qquad p_4=-4\ell_1,\qquad
 p_3=-4A\ell_1-2\ell_2.
\tag{3.6}
\]
Write
\[
 W_6=k_0Y^6+k_1XY^5+k_2X^2Y^4+k_6X^6.
\]
The defect-two equations include
\[
\begin{aligned}
 [x]\,W^3(P^3)&\doteq\ell_1^2,\\
 [y]\,W^3(P^3)&\doteq
 186A\ell_1^2+60k_0+37\ell_1\ell_2,\\
 W^2(P^2)&\doteq
 10A^2\ell_1^2+21Ak_0+5A\ell_1\ell_2
 +20C\ell_1^2+6k_1.
\end{aligned}
\tag{3.7}
\]
Thus
\[
 \ell_1=k_0=k_1=0.
\tag{3.8}
\]

Use weights \(w(x)=2,w(y)=1\), with threshold seven.  The complete
weight-seven faces are
\[
\begin{aligned}
 P_{[7]}&=x^2y^3-2\ell_2x^3y,\\
 W_{[7]}&=X^3Y+\ell_2X^2Y^3+h_0Y^7.
\end{aligned}
\tag{3.9}
\]
Here \(h_0=[Y^7]W_7\); no other lower polynomial piece or higher normalized
operator monomial has weight seven.  Their third and fourth moments give
\[
\begin{aligned}
 W_{[7]}^3(P_{[7]}^3)&\doteq-20h_0+\ell_2^3,\\
 W_{[7]}^4(P_{[7]}^4)\big|_{h_0=\ell_2^3/20}
 &\doteq\ell_2^4.
\end{aligned}
\tag{3.10}
\]
Consequently \(\ell_2=h_0=0\).

Every remaining operator term has weight strictly greater than seven.
The only equality pair is \((X^3Y,x^2y^3)\), with an \(x\)-derivative
deficit.

### 3.3 The linear-tip chart \(B=0,A\ne0\)

Normalize \(A=1\), retaining arbitrary \(C\).  The defect-two scalar part
of moment two is
\[
\begin{aligned}
 E_2={}&1340C^2\ell_0^2+480C\ell_0\ell_1
 +4C\ell_1p_4+56\ell_0\ell_2+6\ell_0p_2\\
 &+20\ell_1^2+4\ell_2p_4+p_2p_4,
\end{aligned}
\tag{3.11}
\]
while a positive-degree coefficient of moment three is
\[
 E_3=p_4^2+10\ell_0p_4+330\ell_0^2.
\tag{3.12}
\]
If \(\ell_0=0\), then (3.12) gives \(p_4=0\), and (3.11) gives
\(\ell_1=0\).

Suppose instead that \(\ell_0\ne0\), and put
\[
 s=p_4/\ell_0.
\]
Then
\[
 f(s)=s^2+10s+330=0.
\tag{3.13}
\]
Since \(\gcd(f,s+6)=1\), equation (3.11) solves \(p_2\).  Retaining the
complete \(W_6,W_7,P_3,P_2\) data, the \(x\)-coefficient of the
defect-three fourth moment is
\[
 [x]\,W^4(P^4)\doteq
 g(s)=143s^3+840s^2+13860s+480480.
\tag{3.14}
\]
But
\[
 \operatorname{Res}_s(f,g)=889363523400\ne0.
\tag{3.15}
\]
Thus the branch \(\ell_0\ne0\) is empty.

With weights \(w(x)=3,w(y)=2\) and threshold eleven, every remaining
operator term has weight greater than eleven.  The equality pair is
\((X^3Y,xy^4)\), again with an \(x\)-derivative deficit.

### 3.4 The pure-\(y\) chart \(B=A=0\)

Normalize \(C=1\).  The defect-two scalar equation is
\[
 1340\ell_0^2+4\ell_1p_4+p_2p_4=0.
\tag{3.16}
\]
If \(p_4\ne0\), solve (3.16) for \(p_2\).  With complete
\(W_6,W_7,P_3,P_2\) data retained, the defect-three fourth moment has
\[
 [y]\,W^4(P^4)\doteq p_4^3.
\tag{3.17}
\]
Hence \(p_4=0\), and then (3.16) gives \(\ell_0=0\).

Weights \(w(x)=4,w(y)=3\), with threshold fifteen, now give a strict
separator except for \((X^3Y,y^5)\).  That pair has an
\(x\)-derivative deficit.

## 4. The isolated \(x^5\) family

Normalize \(P_5=x^5\).  The defect-one equations give successively
\[
 p_3=-20\ell_5,\qquad p_1=p_2=p_0=0.
\tag{4.1}
\]
Write
\[
 P_3=r_0y^3+r_1xy^2+r_2x^2y+r_3x^3.
\]
The defect-two equations are
\[
 r_1=-340\ell_5^2,\qquad r_0=0.
\tag{4.2}
\]

For weights \(w(x)=1,w(y)=2\) and threshold five, the complete equality
faces are therefore
\[
\begin{aligned}
 P_{[5]}&=x^5-20\ell_5x^3y-340\ell_5^2xy^2,\\
 W_{[5]}&=X^3Y+\ell_5X^5.
\end{aligned}
\tag{4.3}
\]
Their third moment is
\[
 W_{[5]}^3(P_{[5]}^3)
 =301335552000\,\ell_5^3.
\tag{4.4}
\]
Thus \(\ell_5=0\), and (4.1)--(4.2) kill every other equality correction.
The final pair \((X^3Y,x^5)\) has a \(y\)-derivative deficit.

## 5. All-order conclusion

In all four projective charts, the exact moment equations leave a common
positive weight and threshold such that
\[
 w(P)\le W\le w(\Lambda).
\tag{5.1}
\]
Every strict selection has positive weight surplus.  Every equality
selection has a fixed coordinate-derivative deficit linear in the number
of selections.  The weight-defect lemma from the
[`degree-five frontier`](BINARY_DEGREE_FIVE_GVC_FRONTIER.md)
therefore gives
\[
 \Lambda^m(QP^m)=0
\]
for every fixed \(Q\) and all sufficiently large \(m\).

## 6. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quartic_triple_simple_root_gvc.py
```

The checker uses exact sparse SymPy contraction and Singular over
\(\mathbb Q\).  It verifies the radical (3.4), the complete branch
identities (3.7), (3.10)--(3.17), (4.1)--(4.4), and the final support
separators.  The calculations labelled complete retain every normalized
operator and polynomial jet capable of entering the displayed defect and
output coefficient.
