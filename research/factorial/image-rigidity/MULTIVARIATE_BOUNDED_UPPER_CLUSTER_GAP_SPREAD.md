# All-degree coprimality with an arbitrary spread of gaps in a bounded upper cluster

For distinct nonnegative integer exponent columns \(a,b\), let
\[
 Q_j(z)=\sum_{k=0}^j\binom jk
                    \prod_i((j-k)a_i+kb_i)!\,z^k .
\]
Omit equal rows, and let \(M\) be the largest remaining exponent.
Select a nonempty group of at most \(K\) unequal rows and let \(X\)
be their largest gap. Fix an integer \(B\ge0\) and an integer
\(J\ge4\) at least the
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** There are computable constants
\(c,\eta>0\) and \(X_0\), depending only on \(K,B,J\), with
the following property. Sum the gaps outside the selected group,
omitting rows whose larger exponent is at most \(JX\) and rows
whose two exponents lie in \([M-\eta X,M]\). Call this sum \(G\).
If
\[
     X\ge X_0,\qquad G\le cX,\qquad
     M-\max(a_i,b_i)\le BX\quad\text{in the group},
\]
then \(Q_n(z)\) and \(Q_{n+1}(z)\) are coprime for every \(n\ge1\)
over every characteristic-zero field. The smallest selected gap
may stay fixed as \(X\) grows. The number and total gap of the
omitted rows are unrestricted.

Take \(\delta=1/(128J^3)\). Define a sequence of integers by
\[
 T_1=1,\qquad
 T_{j+1}=1+\left\lceil\frac{2KT_j}{\delta}\right\rceil,
                  \qquad 1\le j\le K.
\]
Set
\[
 c=\frac{\delta}{2T_K},\qquad \eta=\frac1{5T_K},\qquad
 X_0=\max\left\{24J,12B,\max_{1\le j\le K}
                                T_jD_0(J,T_j,BT_j)\right\},
\]
where \(D_0(J,A,B')\) is the explicit bound in the
[large factorial mark prime budget](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md).
The numerical value of \(J\) is unevaluated.

Suppose first that the total curvature is below \(12\),
since [the curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
covers the other case. Order the selected gaps as
\(d_1=X\ge d_2\ge\cdots\ge d_r>0\) and append \(d_{r+1}=0\).
Let \(i\) be the first index with \(d_{i+1}<X/T_{i+1}\).
It exists since \(i=r\) meets the inequality. By minimality,
\(D=d_i\ge X/T_i\), also for \(i=1\).
The first \(i\) rows have gaps in \([D,AD]\), where \(A=T_i\).
Write \(M'\) for their largest exponent and let \(B'=BT_i\).
Their larger exponents lie in \([M'-B'D,M']\), and
\(M-M'\le BX\le B'D\).

Apply the large factorial mark prime budget to these first \(i\)
rows, with \(D,A,B',M'\) and \(H=J^2(A+B')\).
Here \(D\ge X/T_i\ge D_0(J,A,B')\).
Rows of larger exponent at most \(JX\) contribute no prime
larger than \(HD\), since \(HD/J\ge JT_iD\ge JX\).

Likewise a row with both entries in \([M-\eta X,M]\) contributes
no large prime that divides a mark. For any marked order \(m\)
the marks \(x_u=mM'-u\) satisfy \(mD/5\le u<mAD\)
by the bounds on the window in that proof. Every rising
factor at order \(m\) is greater than
\(m(M-\eta X)\ge m(M'-D/5)\ge x_u\). Its difference from
\(x_u\) is positive and at most
\(m(M-M')+u\le m(A+B')D\).
At order \(m'=m-1\), let \(r\) be a rising factor of this row.
Then \(m'(M-\eta X)<r\le m'M\). Consequently
\[
 0<m(r-m'M')+m'u\le mm'(M-M')+m'u\le mm'(A+B')D .
\]
This is \(mr-m'x_u\), so both orders give units at every
prime \(p>HD\) dividing a mark. The proof of the prime budget
uses only the large-prime parts of retained factors at the
marks, so these rows can be omitted from their log-height
estimate.

The total retained gap contributing to that estimate is at most
\(G+KX/T_{i+1}\), where the second term bounds the gap sum of the
rows dropped from the original group. Our choices give
\[
 G+\frac{KX}{T_{i+1}}\le\frac{\delta X}{T_i}\le\delta D .
\]
The first row, of gap \(X\), gives \(M'>X^2/12\) from the
curvature bound. Hence \(M\le M'+BX\le2M'\), and
\(M'\ge 2J\) by the lower bounds on \(X\). Each contributing
retained row of upper exponent \(M_i\) has
\(\log(J M_i)\le2\log M'\). Its contribution is consequently
bounded as required by (1) in the prime budget,
\[
 \mathcal L\le 2\delta D\log M'=\frac{D\log M'}{64J^3}.
\]
This certifies physical coprimality in every degree for the
full original profile.
