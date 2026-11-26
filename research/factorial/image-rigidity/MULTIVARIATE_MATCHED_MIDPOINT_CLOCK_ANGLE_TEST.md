# A signed matched-moment test at the two exact midpoint clocks

Let \(a,b\) be distinct nonnegative integer exponent columns. Remove equal
rows and put \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\), and
\(\delta_i=|d_i|/m_i\le2\). Write \(L_j(s)\) for the
endpoint-normalized log-factorial symbol per degree. At order \(j\)
the exact polynomial, after positive normalization of its coefficient
and variable, is
\[
 P_j(z)=\sum_{k=0}^j\binom jk e^{jL_j(k/j)}z^k .
\]
Define the row and total midpoint curvatures and the total signed skew by
\[
 C_{j,i}=L_{j,i}''(1/2),\qquad C_j=\sum_i C_{j,i},\qquad
 B_j=L_j'''(1/2).
\]
Here \(L_{j,i}\) is the contribution of row \(i\) to \(L_j\).
The exact Gaussian reference at order \(j\) is
\[
 A_j(z)=\sum_{k=0}^j\binom jk e^{-C_j k(1-k/j)/2}z^k .
\]

**Theorem.** Fix \(n\ge2\) and \(N=n+1\).
With \(t=C_n/2>0\) one may take
\[
                      g=\frac{t}{72nN(1+t)}
\]
as a lower bound on the cross-degree root chord gaps of the exact
two Gaussian references. The comparison below also allows any
stronger lower bound for that pair. Put
\[
 e_j=\frac{|B_j|}{72\sqrt3}
       +\frac1{24}\sum_i C_{j,i}\delta_i^2,\qquad j=n,N .
\]
The exact factorial pair has disjoint root arguments if
\[
                  2(ne_n)^{1/n}+2(Ne_N)^{1/N}<g.             \tag{1}
\]
The comparison uses the signed skew and the curvature mass of each
row. It allows off-circle factorial roots and private rows.

**All-order consequence.** Put
\[
 W=\sum_i\frac{d_i^2}{m_i},\qquad
 B_\infty=-\sum_i\frac{d_i^3}{m_i^2},\qquad
 D=\sum_i\frac{d_i^4}{m_i^3}.
\]
For every positive rational \(w\) there is a computable positive
tolerance \(e_*(w)\) such that the exact ordinary moment polynomials
have disjoint consecutive root arguments in every degree whenever
\[
 W\ge w,\qquad
 \frac{|B_\infty|}{72\sqrt3}
   +\left(\frac1{24}+\frac1{144\sqrt3}\right)D<e_*(w).         \tag{2}
\]
One can use the exact errors \(e_j\) in place of the expression in
(2) on the finite set of orders below the absolute ordinary cutoff.
The test uses signed skew and row-curvature-weighted moments.
The tolerance depends on the absolute cutoff
and is unevaluated. The
[absolute weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
allows the lower restriction on \(W\) to be dropped by choosing \(w\)
at its weak cap, or any smaller positive rational value.

## Signed chord bound

Use the exact curvature integral from the
[matched tanh hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md).
Its rowwise expansion at \(1/2\), before taking the sum of the
linear terms, gives for \(0<u<1\)
\[
 \left| L_j''(u)-C_j-(u-1/2)B_j\right|
 \le \sum_i C_{j,i}
       \frac{\delta_i^2|u-1/2|^2}{1-\delta_i|u-1/2|}.
                                                               \tag{3}
\]
Let \(G(s,u)=\min(s,u)(1-\max(s,u))\) be the interval Green
kernel for minus the second derivative with zero endpoint values.
The zero-endpoint second primitive of \(u-1/2\), evaluated at
\(s\), is \(((s-1/2)^3-(s-1/2)/4)/6\). Its sup norm is
\(1/(72\sqrt3)\). Also
\[
 \frac{G(s,u)}{1-\delta_i |u-1/2|}
 \le\frac{\min(u,1-u)}{1-2|u-1/2|}=\frac12 .
\]
Since \(\int_0^1|u-1/2|^2\,du=1/12\), integration of (3)
therefore gives
\[
          |L_j(s)+C_j s(1-s)/2|\le e_j\quad(0\le s\le1).    \tag{4}
\]
The Green bound also justifies the integration when \(\delta_i=2\);
the model remainder times the kernel is then integrable at the
endpoints. The exact symbol is smooth there since the original Gamma
arguments are at least one.

## The Gaussian clocks and a bound on root rays

The exact row formula gives
\[
 C_j=\sum_i d_i^2 j\psi'(1+jm_i).
\]
The positive curvature increment in the matched hierarchy implies
\(C_N>C_n\), and the strict decrease of \(\psi'\) on positive
arguments gives \(C_N/N<C_n/n\). Thus the upper raw Gaussian
clock \(C_N/(2N)\) lies strictly between \(C_n/(2N)\) and
\(C_n/(2n)\). The
[polynomial physical gap theorem](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md)
gives circular angular distance greater than \(t/[36nN(1+t)]\)
throughout this clock interval, at all \(t>0\).
A chord is at least half its circular distance, which gives the
announced \(g\).

The exponents in \(P_j\) and \(A_j\) are nonpositive, since
\(L_j\) is convex with zero endpoints. The exponential has derivative
at most one there, so (4) gives
\[
                    |P_j(z)-A_j(z)|\le j e_j(1+|z|)^j.
\]
Write any root of \(P_j\) as \(r\omega\), with \(r>0\) and
\(|\omega|=1\). Each unit-circle root \(\zeta\) of \(A_j\)
satisfies
\[
                   \frac{|r\omega-\zeta|}{1+r}
                         \ge\frac{|\omega-\zeta|}{2}.
\]
This follows on squaring from \(4r/(1+r)^2\le1\).
Factor the monic \(A_j\) at the exact root to conclude
\[
      \prod_{A_j(\zeta)=0}|\omega-\zeta|\le 2^j j e_j .
\]
Some reference root is therefore within \(2(je_j)^{1/j}\) of the
exact root's unit projection, with zero error interpreted exactly.
A common ray would contradict (1) by the triangle inequality.

## One tolerance for all degrees

Choose a computable integer \(J\ge4\) at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The exact ordinary pair is already covered when \(N\ge J\).
At \(N=2\), the first moment has its root on the negative real
axis, whereas the quadratic is strictly positive on that axis.
Indeed its value at \(-r\) is the expectation of
\((X^a-rX^b)^2\) under independent mean-one exponential variables,
which is positive because the exponent columns are distinct.

For \(3\le N<J\), the matched hierarchy gives \(C_n\ge W/2\).
If \(W\ge w\), the lower normalized time is \(t=C_n/2\ge w/4\).
The polynomial physical gap bound gives a common chord gap at least
\[
                         g=\frac{w}{72J^2(w+4)} .
\]
There is no upper curvature restriction. This bound uses the absolute
ordinary cutoff without requiring an enumeration of reference
resultants.

Finally the curvature integral with \(h(t)=t/(e^t-1)\)
gives \(C_{j,i}\le d_i^2/m_i\). The bound
\(0\le1-h(t)\le t/2\) follows from \(-1/2\le h'(t)\le0\).
In the signed curvature derivative it gives
\[
 |B_j-B_\infty|
 \le\frac1j\sum_i\frac{|d_i|^3}{m_i^3}\le D/j,
\]
where \(|d_i|\ge1\) is used in the last step. Thus the left side
of (2) bounds every \(e_j\) for \(j\ge2\). Set
\(e_*(w)=(g/8)^{J-1}/J\). For \(j<J\), the associated
projection bound \(2(je_j)^{1/j}\) is smaller than \(g/4\).
Equation (1) now applies at every remaining terminal.

For an asymmetric family with exact limiting skew cancellation,
take blocks consisting of one increasing row \((L,L+1)\)
and nine decreasing rows \((3L+2,3L+1)\). The midpoint
of each decreasing row is three times the increasing midpoint
\(m=L+1/2\). Each block gives \(B_\infty=0\), curvature
\(4/m\), and \(D=4/(3m^3)\). For any fixed \(W_0>0\),
use \(\lfloor W_0 m/4\rfloor\) blocks as \(L\) grows.
The total curvature tends to \(W_0\), the signed skew
still cancels, and \(D=O(m^{-2})\). A bound proportional
to the largest relative gap would be of order \(m^{-1}\).
The signed test covers this family in every degree for all
sufficiently large integer \(L\).
