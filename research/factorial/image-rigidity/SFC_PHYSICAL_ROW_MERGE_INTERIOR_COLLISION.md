# An interior collision on a physical Gaussian row-merger path

There is a three-column ordinary factorial support and a
single collision-free row merger such that the order-two,
order-three, and order-four moment forms are projectively
regular at both ordinary endpoints. Along the actual
Gaussian interpolation between those endpoints they have
a common nonzero projective zero at a physical squared
correlation strictly between \(21/25\) and \(17/20\).

Take the rows of the source exponent matrix to be
\[
 \begin{pmatrix}
 5&2&3\\
 2&5&3\\
 5&4&4\\
 4&5&4\\
 5&5&4
 \end{pmatrix}.
\]
Correlate the Gaussian coordinates in the first two rows
and leave the last three coordinates independent.
At correlation one, the first two rows merge to \((7,7,6)\).
The two remaining swapped rows and the last row distinguish
all three merged columns.

This disproves strict nonvanishing along the entire physical
row-merger path for actual integer occupation profiles.
The particular source and merged forms both satisfy the
consecutive condition for this window.

## Genuine physical moments

Let \(X_1=|Z_1|^2\) and
\(X_2=|\sqrt q Z_1+\sqrt{1-q}Z_2|^2\), with independent
standard circular complex Gaussians \(Z_1,Z_2\) and
\(0\le q\le1\). Let \(X_3,X_4,X_5\) be independent
unit-rate exponentials, independent of the pair.
Write \(U,V,W\) for the monomials with the displayed columns.
For all nonnegative integers \(i,j,k\), define
\[
 p_\ell=i A_{\ell1}+j A_{\ell2}+k A_{\ell3}
 \quad(1\le\ell\le5),
\]
where \(A\) is the displayed matrix. Their joint moments
are
\[
 {\mathbb E}_q[U^i V^j W^k]
 =\left(\prod_{\ell=1}^5 p_\ell!\right)
       \sum_{h=0}^{\min(p_1,p_2)}
                   \binom{p_1}{h}\binom{p_2}{h}q^h .      \tag{1}
\]
This is the
[correlated Gaussian identity](SFC_ROW_MERGING_OBSTRUCTIONS.md#1-beta-jacobi-and-confluent-divided-difference-identities),
with \(q=|\rho|^2\) in the convention where \(\rho\) is
the off-diagonal Gaussian covariance. At \(q=0\) it gives
the ordinary source factorials. At \(q=1\) Vandermonde's
identity gives the ordinary factorial for the summed row.

Swapping \(X_1,X_2\) and simultaneously \(X_3,X_4\)
preserves the law, interchanges \(U,V\), and fixes \(W\).
Hence every mixed moment of \(B=U-V\) to an odd power
with \(A_+=U+V\) and \(W\) vanishes. For a real field \(t\)
put
\[
 R=A_++tW,\qquad V_2={\mathbb E}_q B^2,\qquad
 V_4={\mathbb E}_q B^4,\qquad M_2={\mathbb E}_q R^2 .
\]
At every physical \(q\) the last three independent variables
distinguish \(U\) and \(V\), so \(V_2>0\). The real
polynomial \(R\) cannot vanish identically, as its
\(U,V\) coefficients are both one and the three columns
are distinct. Full positive-coordinate support gives
\(M_2>0\), including at \(q=1\) after identifying the first
two variables.

## Two scalar polynomials

Use the complex coefficient vector
\((1+i s,1-i s,t)\), so the phase is \(R+i s B\).
Its second moment vanishes exactly when \(s^2=M_2/V_2\).
After substituting this positive real value, the third
and fourth moments, multiplied by \(V_2\) and \(V_2^2\)
respectively, are
\[
 \begin{aligned}
 Q_3(q,t)&=V_2\,{\mathbb E}_q R^3
                  -3 M_2\,{\mathbb E}_q(R B^2),\\
 Q_4(q,t)&=V_2^2\,{\mathbb E}_q R^4
                  -6 V_2 M_2\,{\mathbb E}_q(R^2B^2)
                  +V_4 M_2^2 .
 \end{aligned}                                             \tag{2}
\]
Formula (1) makes these integer polynomials in \(q,t\).
The
[exact certificate](verify_physical_row_merge_interior_collision.py)
checks the following strict signs on the rational rectangle
\[
             [21/25,17/20]\times[-2540,-2520].
\]
On the left field edge \(Q_3(q,-2540)<0\), on the right
field edge \(Q_3(q,-2520)>0\), and on the entire rectangle
\(\partial_t Q_3(q,t)>0\).
Thus each \(q\) in the correlation interval has a unique
cubic root \(t(q)\) in the open field interval, and it
varies continuously. At the two correlation edges the
certificate encloses the root and gives
\[
\begin{array}{c|c|c}
q &\text{an enclosure for }t(q)& \operatorname{sign} Q_4(q,t(q))\\ \hline
21/25 &
[-2650734275/1048576,\ -1325367135/524288]&-\\
17/20 &
[-1328604315/524288,\ -2657208625/1048576]&+
\end{array}
\]
Therefore continuity gives some \(q\) strictly between
the two edges with \(Q_4(q,t(q))=0\). Taking the positive
square root \(s=\sqrt{M_2/V_2}\) then gives a nonzero
projective coefficient vector with three consecutive
zero moments.

## Projective regularity at the two endpoints

Use general complex coefficients \((r+z,r-z,t)\).
The phase is \(r A_++tW+zB\). If \(r\ne0\), divide the
vector by \(r\), and eliminate \(z^2=-M_2/V_2\) by the
second moment equation. The third and fourth equations
then reduce to the same two polynomials (2) at the
rescaled field. The checker gives
\(\gcd_{\mathbb Q[t]}(Q_3(q,t),Q_4(q,t))=1\) at both
\(q=0,1\), excluding a common complex root in this chart.

If \(r=0\), the second moment reads
\(t^2{\mathbb E}_q W^2+z^2 V_2=0\).
Any nonzero solution has \(t\ne0\).
Substituting into the third moment gives
\[
  \frac{t^3}{V_2}
      \bigl(V_2 {\mathbb E}_q W^3
             -3{\mathbb E}_q W^2\,{\mathbb E}_q(WB^2)\bigr).
\]
The parenthesis is the coefficient of \(t^3\) in
\(Q_3(q,t)\). The exact endpoint cubic degrees are three,
so this chart contains no common zero either.

## Reproduction

Run the checker from the repository root with
python-flint installed.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_physical_row_merge_interior_collision.py
~~~

It constructs (2) from the exact integer coefficients
in (1). To certify strict signs over an interval, it
substitutes an affine map from \([0,1]\) and converts
the resulting power coefficients \(c_j\) to Bernstein
coefficients
\(\sum_{j\le k}c_j\binom{k}{j}/\binom{d}{j}\)
for degree \(d\). For the rectangle it applies the
same conversion separately in both variables.
A strict common sign of these coefficients implies
the indicated sign everywhere. The two edge cubic
roots are enclosed by exact rational bisection.
All certificate tests use integer or rational arithmetic.
