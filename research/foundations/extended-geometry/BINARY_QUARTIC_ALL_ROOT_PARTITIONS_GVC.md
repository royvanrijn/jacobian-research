# Binary sextic GVC for every quartic leading symbol

## 1. Theorem

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_4+\Lambda_5+\cdots
\]
be a constant-coefficient operator in two variables with lowest positive
order four, and let \(\deg P=6\).

> **Theorem 1.1 — complete quartic-leading sextic row.**
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

Thus no genuinely nonhomogeneous binary GVC counterexample occurs in the
\((r,\deg P)=(4,6)\) row.  Arbitrary lower pieces of \(P\), arbitrary
higher operator jets, and all five root partitions of the binary quartic
are included.

## 2. Hall classification of the leading locus

After scalar extension, split
\[
 \Lambda_4=\prod_{i=1}^4D_{v_i},\qquad
 P_6=\prod_{j=1}^6L_j.
\tag{2.1}
\]
Retaining the translation variable in the split-symbol constant-term
construction and applying the Duistermaat--van der Kallen theorem shows
that the four derivative copies cannot be matched to four distinct
polynomial factors on which they act nontrivially.

If a root direction has multiplicity \(e\), and \(c\) polynomial factors
annihilate that direction, Hall failure is
\[
 6-c<e,\qquad\text{equivalently}\qquad c\ge7-e.
\tag{2.2}
\]
A deficient subset containing two nonparallel derivative directions sees
all six polynomial factors, so no other Hall failure is possible.

Consequently, at a root of multiplicity \(e\), coordinates may be chosen so
that
\[
 M_e=X^eY^{4-e},\qquad
 P_6=y^{7-e}C_{e-1}(x,y),\qquad 1\le e\le4.
\tag{2.3}
\]
The checker enumerates every distribution of the six factors for the five
quartic root partitions.

Local formal division removes multiples of \(M_e\).  At defect \(d\), a
normalized operator monomial \(X^aY^{4+d-a}\) therefore satisfies
\[
 a<e\quad\text{or}\quad 4+d-a<4-e.
\tag{2.4}
\]
The formal differential unit is invertible on polynomials and changes
neither (1.1) nor (1.2).

## 3. Quadruple root

Put \(\Lambda_4=X^4\) and
\[
 P_6=Cy^6+Axy^5+Bx^2y^4+Dx^3y^3.
\tag{3.1}
\]
The exact defect-one radical has three charts, according to the last
nonzero coefficient \(D,B,A\), plus the pure-\(y\) boundary.

### 3.1 The \(D\)-chart

All sub-threshold jets reduce to the weight-twelve face
\[
\begin{aligned}
 W&=X^4+tX^3Y^3-\frac{33}{80}t^2X^2Y^6
       +\frac{1029}{1920}t^3XY^9+uY^{12},\\
 P&=x^3y^3-\frac32t x^4.
\end{aligned}
\tag{3.2}
\]
Moment four is a nonzero multiple of
\[
 44113t^4+51200u.
\tag{3.3}
\]
After solving (3.3), moment five is a nonzero multiple of \(t^5\).
Thus \(t=u=0\).

### 3.2 The \(B\)-chart

The terminal weight-eight face is
\[
\begin{aligned}
 W&=X^4+\ell X^3Y^2+hX^2Y^4+kXY^6+zY^8,\\
 P&=x^2y^4+p x^3y^2+q x^4,\qquad
 q=-2h-\frac12\ell p.
\end{aligned}
\tag{3.4}
\]
Exact moments two through six have radical
\[
 \sqrt J=(z,k,h,\ell p).
\tag{3.5}
\]
The two components are
\[
 (W,P)=(X^4+\ell X^3Y^2,x^2y^4),\qquad
 (W,P)=(X^4,x^2y^4+p x^3y^2),
\tag{3.6}
\]
and both have a linear \(x\)-derivative deficit.

### 3.3 The remaining boundaries

On \(D=B=0,A\ne0\), moment two kills the tilted pair
\[
 W=X^4+aY^5,\qquad P=xy^5-a x^5.
\tag{3.7}
\]
On the pure-\(y\) boundary, moment two first kills \(Y^5\).  The remaining
\(XY^4/x^3y^2\) ratio satisfies
\[
 s+14=0,\qquad s^2+30s+2002=0,
\tag{3.8}
\]
which are coprime.  A common positive weight then leaves only
coordinate-deficit equality faces.  This closes the entire quadruple-root
orbit.

## 4. Triple root

Use \(M_3=X^3Y\), retaining the strict cofactor \(AX^4\).

### 4.1 Non-pure charts

On the \(x^2y^4\) chart, the terminal common-weight face is
\[
\begin{aligned}
 W&=X^3Y+\ell X^2Y^3+hXY^5+kY^7,\\
 P&=xy^5+p x^2y^3+q x^3y,\qquad
 q=-20h-2\ell p.
\end{aligned}
\tag{4.1}
\]
Moments two through five have radical
\[
 \sqrt J=(k,h,\ell p).
\tag{4.2}
\]
Both components are coordinate-deficit axes.

On the \(xy^5\) boundary, moments two kill the \(Y^5\) and \(XY^4\)
migrations.  The \(Y^6/x^4\) ratio is eliminated by
\[
\begin{aligned}
 180180u^2+180uv+v^2,\qquad
 11174042880u^3+7567560u^2v+15120uv^2+143v^3,
\end{aligned}
\tag{4.3}
\]
whose resultant is
\[
 257444582155895534400.
\tag{4.4}
\]

### 4.2 Pure-sixth-power endpoint

At \(P_6=y^6\), defect one first kills the \(Y^5/x^5\) pair.  The highest
remaining tilted face
\[
 W=X^3Y+aXY^4,\qquad
 P=y^6+b x^2y^3+q x^4
\tag{4.5}
\]
has radical \((q,ab)\) through moment four.

Now retain every adjacent equal-weight spectator.  With
\[
\begin{aligned}
 W={}&X^3Y+AX^4+\ell_1XY^4+\ell_2X^2Y^3
                 +h_0Y^6+h_1XY^5,\\
 P={}&y^6+p_1xy^4+p_2x^2y^3+q_3x^3y,
\end{aligned}
\tag{4.6}
\]
moments one through six have the exact radical
\[
\sqrt J=
(q_3,h_0,h_1p_2,\ell_2p_2,\ell_1p_2,\ell_1p_1).
\tag{4.7}
\]

If \(p_2\ne0\), then
\(\ell_1=\ell_2=h_1=0\).  The only balanced tail is
\[
 W=X^3Y+k_0Y^7,\qquad P=x^2y^3,
\tag{4.8}
\]
and its third moment is a nonzero multiple of \(k_0\).  Thus \(k_0=0\);
every later pure operator jet has defect at least four.

If \(p_2=0,p_1\ne0\), then \(\ell_1=0\).  Each non-pure operator selection,
together with the polynomial \(x\)-degree needed to absorb it, contributes
total defect at least three.

If \(p_1=p_2=0\), there are exactly two further cost-two balances.  They are
\[
\begin{aligned}
 (X^3Y+\ell_1XY^4,\ y^6+q_2x^2y^2),\\
 (X^3Y+\ell_1XY^4,\ y^6+r_3x^3).
\end{aligned}
\tag{4.9}
\]
Their first scalar coefficients are nonzero multiples of
\[
 \ell_1^2q_2,\qquad \ell_1^3r_3
\tag{4.10}
\]
in moments two and three.  After these products vanish, every non-pure
operator selection together with the polynomial \(x\)-degree needed to
absorb it contributes total defect at least three.  Pure operator selections
also have defect at least three.

More explicitly, on the \(p_2\ne0\) component, if \(N\) operator selections
carry \(x\)-derivatives and \(S=m-N\) are pure, then
\[
 D_\Lambda\ge4S,\qquad
 D_P\ge\frac{3N-\deg_xQ}{2}.
\tag{4.11}
\]
Since \(N\le(2m+\deg_xQ)/3\), (4.11) gives
\[
 D_\Lambda+D_P>2m+\deg Q
\tag{4.12}
\]
for all sufficiently large \(m\).  The other components have the stronger
bound \(D_\Lambda+D_P\ge3m-O_Q(1)\).  This contradicts the necessary degree
condition for a nonzero term and closes the triple-root endpoint.

## 5. Double root

Use \(M_2=X^2Y^2\), retaining the generic strict quadratic cofactor.

### 5.1 Non-pure chart

On the \(xy^5\) chart, the only possible increasing migration is
\[
 \ell_1\longleftrightarrow p_2,\qquad
 h_0\longleftrightarrow q_3,\qquad
 k_0\longleftrightarrow r_3.
\tag{5.1}
\]
Extremal coefficients of moments two, three, and four are nonzero
multiples of
\[
 \ell_1^2,\qquad h_0^2,\qquad k_0^2.
\tag{5.2}
\]
The terminal face
\[
 W=X^2Y^2+hXY^5+zY^8,\qquad
 P=xy^5-30h x^2y^2
\tag{5.3}
\]
has first scalar equations
\[
 2h^2+7z,\qquad h(55h^2+476z).
\tag{5.4}
\]
Hence \(h=z=0\).

### 5.2 Pure-sixth-power endpoint

At \(P_6=y^6\), defect one kills \(p_4,p_5\), and moment two kills the
remaining \(Y^5/x^2y^3\) pair.

The complete weight-relevant face is
\[
\begin{aligned}
 W={}&X^2Y^2+\ell XY^4+h_0Y^6
 +AX^3Y+h_1XY^5+k_0Y^7\\[-2mm]
 &+BX^4+k_1XY^6+z_0Y^8,\\
 P={}&y^6+p\,xy^4+q_2x^2y^2+q_3x^3y+\rho x^3.
\end{aligned}
\tag{5.5}
\]
Moments one through six have radical
\[
 \sqrt J=(\rho,q_3,q_2,h_0,\ell p).
\tag{5.6}
\]

If \(p=0\), every surviving pure operator jet has defect at least three.
Every \(x\)-bearing operator selection, together with the polynomial
\(x\)-degree needed to absorb it, also costs at least three total defect.
Thus \(D_\Lambda+D_P\ge3m-O_Q(1)\).

If \(p\ne0\), then \(\ell=0\).  A leading \(X^2Y^2\) selection needs two
\(p\,xy^4\) factors, while an \(XY^5\) selection already has operator
defect two.  Pure selections have defect at least three.  The resulting
integer linear bound is
\[
 D_\Lambda+D_P\ge\frac52m-O_Q(1).
\tag{5.7}
\]
Both bounds eventually exceed \(2m+\deg Q\), so the double-root endpoint
is closed.

## 6. Simple root and squarefree symbols

Use \(M_1=XY^3\) and retain the complete strict cubic cofactor
\[
 XY^3+AX^2Y^2+BX^3Y+GX^4.
\tag{6.1}
\]
Exact defect radicals give the following triangular support reduction:

1. defect one kills the \(x^3,x^4,x^5\) terms of \(P_5\); moment two then
   kills \(Y^5\leftrightarrow xy^4\);
2. defect two kills the \(x^2,x^3,x^4\) terms of \(P_4\), with
   \(q_1=-120s\); moment two kills \(Y^6\leftrightarrow xy^3\);
3. defect three kills the \(x^2,x^3\) terms of \(P_3\);
4. defect four kills the \(x^2\) term of \(P_2\);
5. defects five and six introduce no new dangerous terms.

Consequently every surviving \(x\)-bearing polynomial term has defect at
least three.  Every surviving pure-\(Y\) operator jet also has defect at
least three.

For a monomial selection in \(\Lambda^m(QP^m)\), let \(N\) be the number of
operator factors carrying an \(x\)-derivative, \(S=m-N\), and \(R\) the
number of selected \(P\)-factors carrying \(x\).  Nonvanishing requires
\[
 R+\deg_xQ\ge N.
\tag{6.2}
\]
Therefore
\[
 D_\Lambda+D_P
\ge3S+3R
\ge3m-3\deg_xQ.
\tag{6.3}
\]
On the other hand, ordinary degree compatibility requires
\[
 D_\Lambda+D_P\le2m+\deg Q.
\tag{6.4}
\]
Equations (6.3)--(6.4) are incompatible when
\[
 m>\deg Q+3\deg_xQ.
\tag{6.5}
\]
This closes every simple-root component, including all squarefree quartic
symbols.

## 7. Completion of the row

Every binary quartic has a root of multiplicity \(1,2,3\), or \(4\).
The Hall classification (2.3) lists every leading pure-zero component at
that root.  Sections 3--6 close every projective chart and every
pure-sixth-power intersection.  The finite-tail estimates include arbitrary
higher normalized operator jets, not merely the jets displayed in the
radical calculations.

Hence Theorem 1.1 holds for all five quartic root partitions.

## 8. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_binary_quartic_all_root_partitions_gvc.py
```

The checker uses exact sparse contraction in SymPy and radical computations
over \(\mathbb Q\) in Singular.  It verifies the Hall classification, every
terminal repeated-root face, both complete pure-sixth-power endpoint
radicals, the triple-root \(Y^7\) terminal coefficient, and the complete
simple-root defect layers through defect four.
