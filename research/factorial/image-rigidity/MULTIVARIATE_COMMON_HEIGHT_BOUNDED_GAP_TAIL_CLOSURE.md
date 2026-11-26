# All-degree closure for common-height tails of bounded gaps

Fix any finite ordinary two-monomial core \(C\), and an integer
\(D\ge1\). Append a nonempty collection of unequal rows, all with
the same upper exponent \(M\), with gaps in \(\{1,\ldots,D\}\).
Their orientations and multiplicities may be chosen independently.

**Theorem.** There is a computable height \(M_0(C,D)\) such
that every indicated profile with \(M\ge M_0(C,D)\) satisfies
physical SFC(2) in every degree. There is no bound on the number
of appended rows. This gives exact exclusion through every
positive Gaussian exceptional strength of the fixed core for
such ordinary tails.

The proof uses any integer \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Its numerical value is unevaluated. In particular a finite
catalogue of fixed cores has a common computable height bound.
The divisibility proof allows unmatched rows, unequal orientation
counts, and every sign of the Gaussian phase response.

## A bounded signed gap forces substantial counts of each sign

Write \(r\) for the number of appended rows, \(r_\pm\) for their
orientation counts, and \(h\) for their total signed gap.
Let \(\Delta_C\) denote the signed degree difference of the core.
At a physical common root the
[absolute signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
gives \(|h|\le H=146J-74+|\Delta_C|\). Since every gap is at
least one and at most \(D\), this implies
\[
                        r_\pm\ge\frac{r-H}{D+1}.
\]
For example \(r_- -D r_+\le -h\le H\) yields the lower bound
for \(r_+\). At \(r\ge2H\), each orientation therefore has a gap
type of multiplicity at least \(r/[2D(D+1)]\).

Set \(b=8D\), \(t_0=r/[4D(D+1)]\), and \(t_s=t_0 b^{-s}\).
At most \(2D\) positive multiplicities occur among the distinct
oriented gap types. One of the \(2D+1\) disjoint intervals
\([t_s/b,t_s)\), \(s=0,\ldots,2D\), therefore contains none
of these multiplicities. Keep as the selected types those of
multiplicity at least \(t_s\). Both signs occur and every
selected type has multiplicity at least
\[
                    t_s\ge\frac{r}{Q_1},
                 \qquad Q_1=4D(D+1)(8D)^{2D}.
\]
The total number of omitted tail rows is at most \(t_s/4\).
Indeed every omitted type has multiplicity less than \(t_s/b\).

For \(r<2H\), if both signs occur, select all types with
positive counts and use \(t_s=1\). If only one sign occurs, reverse
the columns if necessary so it is the increasing one, and use
the first block of any appended row. In all cases put
\[
                           Q=\max(Q_1,2H).
\]

## A native literal factor with a uniform count reserve

Consider a terminal \(4\le m<J\), and set \(m'=m-1\).
If the selected types have both signs, take their smallest
gaps \(d_+\) and \(d_-\) in each orientation. Use the
[common upper central mark](MULTIVARIATE_COMMON_UPPER_EXPONENT_DEEP_MARKS.md#a-common-middle-interval)
at an offset \(u\) with \(0\le u<m\min(d_+,d_-)\) and
\(m\nmid u\). That mark gives an interval \(I\) of width at
most one on which no selected row contributes the literal
factor \(x=mM-u\); at every index outside \(I\) at least one
selected type contributes \(x\). Each selected type has
multiplicity at least \(t_s\). Inside \(I\), the omitted tail
rows can contribute at most \(t_s/4\) copies of \(x\).
Thus the difference in the total literal count between any
outside index and any inside index is at least \(3r/(4Q)\).
For small \(r\) with both signs there are no omitted rows.

If the tail has only increasing rows, choose any gap \(d\)
present and mark \(x=m(M-d)+1\), so \(u=md-1\).
It occurs at every positive index in the chosen row, and
no tail row contributes it at zero. The same count difference
bound holds with \(I=\{0\}\), since this case has \(r<2H\).
In all cases \(0\le u<mD\) and \(m\nmid u\).

## A prime-depth bound independent of multiplicity

Strip each row's smaller-end factorial in each degree.
For \(2\le j<J\) let \(B_j\) be \(j!\) times the product
of the full stripped ramps of the fixed core in degree \(j\).
They depend only on \(C,J\). For a prime \(p\mid x\) put
\(E=v_p(x)\) and
\[
                 L_p=\max_{1\le h\le J^2D}v_p(h).
\]
If \(E>L_p\), then at the marked order all tail factors
different from \(x\) have valuations at most \(L_p\), by
subtracting \(x\). At the other order use
\[
            m(m'M-h)=m'(mM-u)+(m'u-mh),\qquad
                              0\le h<m'd_i .
\]
The last offset is nonzero modulo \(m\), since \(m\nmid u\),
and has magnitude at most \(J^2D\). The valuation of the
left side is its valuation, so every tail factor at order
\(m'\) has valuation at most \(L_p\).

After removing the literal contributions from order \(m\),
the remaining coefficient ordinates lie in
\([0,v_p(B_m)+rmD L_p]\).
At order \(m'\), all slopes have magnitude at most
\(v_p(B_{m'})+rm'DL_p\).
Subtract any such slope times the coefficient index in order
\(m\). For an outside index and its nearest inside index the
resulting difference is at least
\[
           \frac{3r}{4Q}E-v_p(B_m B_{m'}^m)-rJ^2D L_p .
\]
This uses the count reserve and \(m+mm'=m^2<J^2\).
If
\[
                     E>2Q\bigl(v_p(B_m B_{m'}^m)+J^2D L_p\bigr),
\]
the difference is positive. Every edge of shared slope at
order \(m\) is then contained in \(I\), of total width at
most one. By the
[Newton polygon degree comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
a common rational factor has degree at most one. The even
ordinary moment is strictly positive on the real line, so a
common rational factor must have even degree.

It follows that at any physical common root the positive
integer \(mM-u\) divides the fixed bound
\[
 K=\bigl((J^2D)!\bigr)^{2QJ^2D}
           \left(\prod_{2\le j<J} B_j^J\right)^{2Q}.
\]
Indeed \(L_p\le v_p((J^2D)!)\), and the displayed degree
comparison applies whenever the claimed valuation bound fails.
This gives \(M<K+JD+1\), independently of \(r\) and the
multiplicity choices. The first two pairs use
[ordinary low-order separation](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
and \(m\ge J\) uses the absolute cutoff. This proves the result
for the same physical variable and transfers to every
characteristic-zero field by a rational Bézout identity.
