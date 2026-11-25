# Scaled Hermite interlacing for the full Cayley deficit ratio

Let \(n\ge2\), \(N=n+1\), and let \(\mathrm{He}_j\) denote the
monic probabilists' Hermite polynomials. If
\[
                       \frac{n-1}{n}\le \rho\le1,
\]
then the zeros of \(\mathrm{He}_n\) strictly interlace the zeros of
\(\mathrm{He}_N\) multiplied by \(\sqrt\rho\). Both endpoints
of the scale interval are included. The elementary matrix
argument below extends the same comparison from
[\(\rho\in[n/N,1]\)](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md#quantitative-comparison-of-the-hermite-pairs).

The three-term recurrence identifies the zeros of \(\mathrm{He}_j\)
with the eigenvalues of the real Jacobi matrix \(J_j\) with zero
diagonal and off-diagonal entries \(\sqrt{k}\), \(1\le k<j\).
The eigenvalues are simple. Deleting the first or the last vertex
gives strict interlacing, since an eigenvector of a path with
nonzero edges cannot vanish at an endpoint and the three-term
determinant recurrence precludes a shared eigenvalue with that
principal block.

The leading \(n\)-vertex block of \(\sqrt{\rho}\,J_N\) is
\(\sqrt{\rho}\,J_n\), with positive eigenvalues at most those
of \(J_n\). The trailing block has edge ratios against \(J_n\)
given by
\[
                         c_k=\sqrt{\rho(k+1)/k},
                       \qquad 1\le k\le n-1.
\]
They decrease with \(k\) and their minimum is at least one.
To compare positive eigenvalues, put \(v_n=0\) and define
\(v_k=\log c_k-v_{k+1}\) successively down to \(k=1\).
Each \(v_k\) is an alternating sum of a nonnegative decreasing
sequence, so \(v_k\ge0\). The trailing block is \(J_n\) multiplied
on both sides by \(\operatorname{diag}(e^{v_k})\).
In odd/even vertex order, the positive eigenvalues of either
matrix are the positive singular values of its off-diagonal
rectangular block. Scaling its rows or columns by numbers at least
one can only increase the singular values, by the variational
comparison of the corresponding Gram matrices. This proves
that the positive trailing eigenvalues are at least those
of \(J_n\).

Each positive eigenvalue of \(J_n\) thus lies between the
corresponding positive leading and trailing eigenvalues, which
both lie strictly between the same consecutive eigenvalues of
\(\sqrt{\rho}\,J_N\). Changing the signs on alternating vertices
conjugates each path matrix to its negative, so reflection
gives the negative half. The odd-degree matrix has its simple
middle root zero, and the even-degree matrix has an equal
number of positive and negative eigenvalues. This supplies
the middle comparison. At \(n=1\) the same interlacing conclusion
is immediate for every \(\rho>0\).

For the separately normalized near-flat Cayley anchors,
the order-\(j\) imaginary heights are proportional to the
Hermite roots with variance \(1/[2(j-1)]\). If \(p=1-u>0\) and
\(q=1-w\) are the two outer deficits, moving both anchors to
the common lower scale replaces the upper variance relative
to the lower by
\[
                        \rho=\frac q p\frac{n-1}{n}.
\]
The [linked outer bounds](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md#linked-factors-and-comparison-across-orders)
give \(u^{n/(n-1)}\le w\le u\), hence
\(1\le q/p\le n/(n-1)\). The upper inequality for the deficits
follows from \(1-u^a\le a(1-u)\) for \(a\ge1\).
Thus all the physical values of the deficit ratio are within the
strict anchor-interlacing range. The statement concerns the
Hermite anchors. A transfer to the actual zeros at nonzero
deficit still uses an error bound.

