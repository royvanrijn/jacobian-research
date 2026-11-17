# Binary sextic GVC for every cubic leading symbol

## 1. Theorem

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_3+\Lambda_4+\cdots
\]
be a constant-coefficient operator in two variables with lowest positive
order three, and let \(\deg P=6\).

> **Theorem 1.1 — complete cubic-leading sextic row.**
> If
> \[
>  \Lambda^m(P^m)=0\qquad(m\ge1),
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \]

Arbitrary lower pieces of \(P\), arbitrary higher operator jets, and all
three root partitions of the binary cubic are included.  Consequently the
only row that remained at this stage had lowest operator order two; it is
now closed by the
[complete quadratic-leading theorem](BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md).

## 2. Leading Hall locus and local supports

Factor \(\Lambda_3\) and \(P_6\) after scalar extension.  The translated
split-symbol constant-term construction and Hall's theorem show that a
root direction of multiplicity \(e\) can be deficient only when at least
\(7-e\) polynomial factors lie on its annihilator line.  Thus
\[
 M_e=X^eY^{3-e},\qquad
 P_6=y^{7-e}C_{e-1}(x,y),\qquad 1\le e\le3.
\tag{2.1}
\]
A subset containing two nonparallel derivative directions sees every
polynomial factor, so (2.1) is exhaustive for the partitions
\((3),(2+1),(1+1+1)\).

Local division removes multiples of \(M_e\).  Every weighted face below is
therefore finite.  Once its equality pair has a linear coordinate deficit,
the common-threshold defect lemma from the degree-five theorem absorbs all
strict terms.  The checker separately enumerates the Hall locus.

## 3. Triple root

Put \(\Lambda_3=X^3\), so
\[
 P_6=Cy^6+Axy^5+Bx^2y^4.
\]

### 3.1 The \(B\)-chart

Defect one kills the weight-\((4,3)\) and \((3,2)\) tilts.  The
weight-\((2,1)\) face is
\[
\begin{aligned}
 W&=X^3+aX^2Y^2+bXY^4+cY^6,\\
 P&=x^2y^4-4ax^3y^2+(2a^2-2b)x^4.
\end{aligned}
\tag{3.1}
\]
Moment two gives
\[
 a^2-27b,\quad 7a^3-25ab-84c,\quad
 8a^4-33a^2b+24ac+43b^2.
\tag{3.2}
\]
The first two equations give \(b=a^2/27,\ c=41a^3/567\), and the last
becomes \(43744a^4/5103\).  Hence \(a=b=c=0\).

The last two slopes have weights \((3,1)\) and \((4,1)\).  Successive
moments solve their three trailing coefficients and leave respectively
nonzero multiples of \(a^4\).  The terminal pair
\((X^3,x^2y^4)\) has an \(x\)-deficit linear in \(m\).

### 3.2 The \(A\)-chart

The first three slopes \(4/3,3/2,5/3\) die in moment two.  At slope two,
moment one eliminates \(q\) and moments two through seven have radical
\[
 (h,k,q,\ell p).
\tag{3.3}
\]
The \(p\ne0\) component is the already proved degree-five
\((X^3,x^2y^3)\) chart, with \(xy^5\) strictly below its terminal face.
The \(\ell\ne0\) component is the quartic-leading
\((X^2Y^2,xy^5)\) face, with \(X^3\) strict.  On the coordinate boundary,
the \(XY^5/x^3\) face dies in moment two.  Weight \((3,1)\) is then a
strict separator.

### 3.3 Pure-sixth-power endpoint

The \(Y^4/x^3y^2\) tilt dies in moments one and two.  The next face
\[
 (X^3+aXY^3,\ y^6+bx^2y^3+cx^4)
\]
has \(c=-ab/2\), followed by the incompatible nonzero ratios
\[
 328a+7b=0,\qquad 12a-13b=0.
\tag{3.4}
\]

If \(a=0\), the \(Y^5/x^3y\) tilt dies and the complete weight-\((2,1)\)
face
\[
\begin{aligned}
W&=X^3+\ell X^2Y^2+hXY^4+kY^6,\\
P&=y^6+pxy^4+qx^2y^2+rx^3
\end{aligned}
\tag{3.5}
\]
has exactly the three support-separated components
\[
\begin{array}{c|c}
\min\deg_XW&\max\deg_xP\\ \hline
3&2\\
2&1\\
1&0.
\end{array}
\tag{3.6}
\]
Equivalently, after moment one eliminates \(r\), the affine chart
saturations give \(k=q=hp=0\) on \(\ell\ne0\), \(k=p=q=0\) on
\(h\ne0\), and no \(k\ne0\) chart.  The possible \(x^2y^3\) spectator is
the closed degree-five cubic chart.

If \(a\ne0\), normalize \(a=1\).  The successive \(Y^5\) and \(Y^6\)
faces are killed by nonzero cubic coefficients, and the only remaining
term above the final weight bound is \(x^2y\).  Its defect-five
coefficient in the second moment is \(20160[x^2y]P\).  Weight \((4,1)\)
is then strict.  This completes the triple-root orbit.

## 4. Double root

Retain the strict cofactor in
\[
 \Lambda_3=X^2Y+\alpha X^3.
\]

On \(P_6=xy^5+Cy^6\), the slopes
\[
 \frac32,\ 2,\ \frac52,\ 3,\ 4
\]
are triangular.  The slope-two equations first force the \(Y^5\) and
\(XY^3\) parameters to vanish.  The remaining scalar faces give nonzero
multiples of their square or cube.  Weight \((5,1)\) is strict.

At \(P_6=y^6\), the slope-\(3/2\) face dies.  The complete slope-two face
\[
\begin{aligned}
W&=X^2Y+bXY^3+cY^5,\\
P&=y^6+pxy^4+qx^2y^2+rx^3
\end{aligned}
\tag{4.1}
\]
has radical
\[
 (c,r,bp).
\tag{4.2}
\]
The \(b\ne0\) component is the \(XY^3/y^6\) ladder of Section 3.3; the
\(p\ne0\) component is the closed degree-five \(X^2Y/xy^4\) chart.  On
their intersection, the \(Y^6/x^2y\) face dies in moment two and weight
\((5,2)\) is final.  Thus both cubic partitions containing a double root
are closed.

## 5. Simple root and completion

Write the complete strict cofactor as
\[
 \Lambda_3=XY^2+AX^2Y+BX^3.
\tag{5.1}
\]
Local division leaves only pure \(Y\)-jets below the relevant weights.
The defect equations are triangular:

1. defect one kills the \(x^2,x^3,x^4,x^5\) terms of \(P_5\), and the
   remaining \(Y^4/xy^4\) pair dies in moment two;
2. defect two kills the \(x^2,x^3,x^4\) terms of \(P_4\), and the
   \(Y^5/xy^3\) pair dies in moment three;
3. defect three kills the \(x^2,x^3\) terms of \(P_3\), and the
   \(Y^6/xy^2\) pair dies in moment two;
4. defect four kills the \(x^2\) term of \(P_2\).

Consequently every surviving pure operator selection and every surviving
\(x\)-bearing polynomial selection has defect at least four.  If \(N\)
operator selections carry an \(x\)-derivative, \(S=m-N\) are pure, and
\(R\) polynomial selections carry \(x\), nonvanishing requires
\[
 R+\deg_xQ\ge N.
\]
Therefore
\[
 D_\Lambda+D_P\ge4S+4R\ge4m-4\deg_xQ.
\tag{5.2}
\]
Ordinary degree compatibility gives
\[
 D_\Lambda+D_P\le3m+\deg Q,
\tag{5.3}
\]
which contradicts (5.2) for \(m>\deg Q+4\deg_xQ\).

Every cubic has a root of multiplicity one, two, or three.  Sections 3--5
close every Hall chart at such a root, including arbitrary strict tails.
This proves Theorem 1.1.

## 6. Reproduction

Run
```bash
.venv/bin/python scripts/verify_binary_cubic_all_root_partitions_gvc.py
```
with SymPy, Singular, and msolve available.  The replay uses exact
characteristic-zero arithmetic.  The msolve calls test emptiness of
explicit affine saturations; no modular sample is promoted to a proof.
