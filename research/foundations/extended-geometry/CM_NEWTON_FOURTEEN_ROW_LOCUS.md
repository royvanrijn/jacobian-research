# The optimal fourteen-row conductor-18 window

**Theorem (optimal initial window).**
For the fixed top with two opposite fourth-order roots,
\(\deg B\leq5\), \(D=-z^{-1}+c+ez\), and correction degree at most two,
the formal response in shifts \(0,\ldots,13\) detects every nonzero
coherent input over every characteristic-zero field. Fourteen is the
least length of a uniform initial window.

The [three-pole reduction](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#three-extra-pole-inputs-for-odd-degree-five-and-even-degree-one)
detects the branch \(b=0\) in shifts 0 and 1. For \(b\ne0\), the first
seven rows eliminate the correction and angular coefficients. Blindness
in the first fourteen shifts gives a common zero of \(S_7,\ldots,S_{13}\)
in \((r,\ell,h,v,w)=(b^2,b_1/b,b_2/b,c,e)\). It suffices to exclude
such zeros over \(\overline{\mathbb Q}\), since a proper ideal in the
normalized rational polynomial ring has a geometric point over
\(\overline{\mathbb Q}\).

## Reduction at 331

The [affine integral-place computation](CM_NEWTON_FOURTEEN_ROW_INTEGRAL_PLACES.md)
gives the unit ideal for these seven rows modulo \(331\), so at any
place over \(331\) at least one normalized coordinate has negative
valuation. Set
\[
\delta=\min\{\nu(r)/2,\nu(\ell),\nu(h)/2,\nu(v),\nu(w)/2\}<0.
\]
After adjoining a scaling parameter of valuation \(-\delta\), scale the
coordinates with weights \((2,1,2,1,2)\). At least one scaled coordinate
is a unit, and the weighted homogenizer has positive valuation. Thus
the resulting residue lies on the weighted projective boundary.

The [projective boundary computation](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#the-projective-boundary-at-331)
gives the same three components with just \(S_7,\ldots,S_{13}\).
The exact normal computations for \(L_1\) and \(L_2\) still exclude
\(L_1\) and restrict \(\bar v/\bar\ell\) on \(L_2\) to
\(2/17,5/56,1/8\). These chart comparisons use \(\delta<0\).
The last residue gives the \(w=0\) endpoint of \(L_3\).

The [nonendpoint charts](CM_NEWTON_331_L3_NORMAL_CHARTS.md)
and the [\(v=0\) endpoint charts](CM_NEWTON_331_L3_VZERO_ENDPOINT.md)
of \(L_3\) pass on the seven-row prefix. In the latter's
\(\tau\)-least chart the leading equation \(xU+vV=0\) has \(U\) a unit.
If \(v=0\) exactly, this forces \(x=0\), hence \(\ell=0\), contrary to
the exact affine identity for \(\ell=v=0\) through shift 13.
Thus the chart can divide by \(v\). In the other charts the first
allows the parameter \(\eta=0\), and the second already requires \(v\)
nonzero. The [\(w=0\) endpoint charts](CM_NEWTON_331_L3_WZERO_ENDPOINT.md)
also pass on the prefix, using the weight-nine span in the second
normalization. The checker replays each exclusion on these same rows.
We are left with
\[
s=\ell/v \in \mathcal{O}_{\nu},\qquad
\bar s\in\{17/2,56/5\}\subset\mathbb F_{331}.
\]

## Two slope minors for the shorter prefix

Use the divided slope rows and ordinary homogenizations \(H_i\) of the
[seventeen-row slope computation](CM_NEWTON_SEVENTEEN_ROW_LOCUS.md#an-exact-univariate-eliminant).
They use \(s=\ell/v,\ x=v^2\) and homogenizing variable \(t\).
The degree-nine multiples of the full row set have rank 685 over
\(\mathbb Q(s)\), and their span contains \(t^9\). This generic rank and
membership follow from the two quotient restriction maps and the pivot
minor in that computation.

Quotient by \(H_0,H_2\) using the same monic reductions in \(r,h\).
Their degree-nine span has dimension 321. The full quotient row span
has dimension 364 and contains the reduction of \(t^9\).
The shorter prefix uses only \(H_0,\ldots,H_6\).
Its compressed matrix has 428 rows and 394 columns, with coefficients
in \(\mathbb Z_{(331)}[s]\). A nonzero \(364\)-minor over
\(\mathbb F_{331}[s]\) therefore proves that its row span over
\(\mathbb Q(s)\) equals the full quotient row span. For any such
minor \(\Delta(s)\), adjugate multiplication expresses
\(\Delta(s)t^9\) in the homogeneous prefix ideal. This follows by
solving for the reduction of \(t^9\) in the selected rows, clearing
the determinant, then lifting the relation through \((H_0,H_2)\).
Dehomogenizing at \(t=1\) makes \(\Delta(s)=0\) at an affine common zero
with \(v\ne0\).

For each candidate residue, the checker selects a minor using
successive modular Schur complements at \(s=\alpha\), then verifies the
same orders over \(\mathbb Q\) on that selected submatrix. The following
are the successive pivot counts for the shorter matrix.

| \(\alpha\) | Pivot counts | Order of the selected determinant |
|---|---|---:|
| \(17/2\) | \(359,2,0,0,0,2,0,1\) | \(19\) |
| \(56/5\) | \(361,1,0,1,0,0,1\) | \(10\) |

In particular
\(\Delta_\alpha(s)=(s-\alpha)^{m_\alpha}U_\alpha(s)\) with
\(U_\alpha\in\mathbb Z_{(331)}[s]\) and \(U_\alpha(\alpha)\) a unit.
If \(s\) has the indicated residue and \(\Delta_\alpha(s)=0\),
then \(s=\alpha\). On each of these exact slices the first seven
residuals already give an affine unit identity after the substitution
\(\ell=\alpha v\) and division of odd rows by \(v\). The exact coefficient
systems at degree nine each have 847 rows and 715 columns, with pivot
ranks 680 and 682 respectively. This rules out both slopes.

The [algebraic blind phase through shift 12](CM_NEWTON_SEVENTEEN_ROW_LOCUS.md#exact-coherent-blindness-through-shift-12)
shows that thirteen rows are insufficient. The sign rotation of the
top gives the same result in the other sign convention.

For fixed algebraic coefficients, the same good-prime persistence
argument detects each nonzero coherent input at one of the moments
indexed by \(m+kp\), \(0\leq k\leq13\), in the class \(p\equiv7\bmod36\)
at every sufficiently large good prime. In the class \(p\equiv25\bmod36\),
shift zero remains sufficient.

## Marker recovery for a known phase

**Corollary (fourteen-row marker rank).**
Fix a phase in the theorem over a characteristic-zero field and let
\(M_{0:13}\) be its formal response matrix in shifts \(0,\ldots,13\).
Use the two coherent coordinates in \(A=(a_0+a_1z)(z^2\pm1)^2\),
with the sign of the top, and the three correction coefficients in
\(C=C_0+C_1z+C_2z^2\). Then
\[
\operatorname{rank} M_{0:13} =
\begin{cases}
5,& B\not\equiv0,\\
2,& B\equiv0.
\end{cases}
\]

For \(B\ne0\), put \(q=\min\{j:[z^j]B\ne0\}\) and \(b_q=[z^q]B\).
The [correction block calculation](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#a-reconstruction-matrix-for-unmarked-moments)
gives column \(C_j\) its first nonzero entry in row \(q+j+1\),
with pivot a nonzero rational multiple of \(b_q\).
The rows \(q+1,q+2,q+3\) thus give a triangular minor of rank three.
They are all in the window since \(q+3\le8\).
If \(B=0\), the phase is even in the height variable \(t\).
Every correction column is zero because it uses only odd powers of \(t\).

The theorem makes the two coherent columns independent modulo the span
of the correction columns. The full rank is therefore two plus the
correction rank, giving the stated formula. At a known phase the
fourteen responses determine all five marker inputs if \(B\ne0\),
and determine the two coherent inputs if \(B=0\). The phase
coefficients are fixed throughout this linear recovery.

## Replay

    sage -python research/foundations/extended-geometry/verify_cm_newton_fourteen_row_locus.py
