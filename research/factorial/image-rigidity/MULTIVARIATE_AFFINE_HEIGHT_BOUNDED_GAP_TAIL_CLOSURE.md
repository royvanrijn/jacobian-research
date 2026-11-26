# All-degree coprimality for bounded gaps at affine heights

Fix a finite ordinary two-monomial core and \(T\ge1\) possible
additional row types
\[
                    (c_i M+a_i,\ c_i M+b_i),\qquad 1\le i\le T,
\]
where the \(c_i\) are fixed positive integers, \(a_i,b_i\) are
fixed distinct integers, and the \(v_2(c_i)\) are pairwise
distinct. Each type has arbitrary nonnegative integer
multiplicity, with at least one appended row in total.

**Theorem.** For all integers \(M\) above a computable threshold
depending only on the core and the possible types, the enlarged
profile satisfies physical SFC(2) in every degree over every
characteristic-zero field. The threshold is independent of the
multiplicities.

The proof marks a first factorial factor at the even order of
each adjacent pair. Choosing a type of maximal \(v_2(c_i)\)
among the types with substantial multiplicity prevents competing
exact affine factors at the other order. A count bound controls
the types of small multiplicity.

For example, the types \((M-1,M)\) and \((2M,2M-1)\) meet
the hypotheses. Giving both types multiplicity
\(\lfloor 2\tau M/3\rfloor\), for any fixed \(\tau>0\), makes
their total signed gap zero and their total curvature tend to
\(\tau\). Their upper exponents differ by \(M\).

## A first-block mark at the even order

Take \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Consider any adjacent pair of orders with terminal in
\([4,J)\), and denote its even order by \(m\) and its odd
order by \(m'\). Put
\[
 d_i=|b_i-a_i|,\quad \ell_i=\min(a_i,b_i),\quad
 D=\max_i d_i,\quad A=\max_i|\ell_i|,\quad C=\max_i c_i .
\]
Let \(r\ge1\) be the total multiplicity. Define
\(b=4JT\), \(t_j=(r/T)b^{-j}\), and \(Q=T b^T\).
Among the \(T+1\) disjoint intervals \([t_j/b,t_j)\) for
\(0\le j\le T\), at least one contains no positive multiplicity.
For that interval set \(t=t_j\) and select the types of
multiplicity at least \(t\). At least one type is selected
because the largest multiplicity is at least \(r/T\). Each
selected type has multiplicity at least \(t\ge r/Q\), and
the total omitted multiplicity \(r_0\) satisfies
\[
                         r_0\le Tt/b=t/(4J).
\]

Choose the selected type \(s\) of largest \(v_2(c_s)\). If
necessary, reverse the two exponent columns, which reciprocates
both original polynomials and preserves coprimality. Type \(s\)
is now increasing. Strip the scalar factorial of the smaller
exponent from each row in each order. At order \(j\), the
possible remaining factors from row type \(i\) are
\[
               F=j c_i M+j\ell_i+h,\qquad 1\le h\le j d_i .
\]
At order \(m\), mark the first factor in type \(s\),
\[
                       x=m c_s M+m\ell_s+1 .
\]
It occurs at every positive index of each row of type \(s\)
and at none of their index-zero coefficients.

For any factor \(F\) in order \(j=m\) or \(m'\), set
\[
                   \Delta=(m c_s)F-(j c_i)x
                  =m c_s(j\ell_i+h)-j c_i(m\ell_s+1).
\]
If \(i\) is selected and \(j=m'\), then \(\Delta\ne0\).
Indeed the first term has \(2\)-adic valuation at least
\(v_2(m)+v_2(c_s)\), whereas the second has valuation
\(v_2(c_i)\), since \(m'\) and \(m\ell_s+1\) are odd.
A zero first term also gives nonvanishing. For \(j=m\) and
selected \(i\ne s\), divide \(\Delta\) by \(m\). The same
comparison now uses \(v_2(c_s)>v_2(c_i)\) and again gives
\(\Delta\ne0\). For \(i=s,j=m\), the only zero offset is
\(h=1\), the marked factor. An omitted row can contribute
at most one zero-offset factor in either degree because
\(\Delta\) is strictly increasing in \(h\).

## A prime-depth bound uniform in the multiplicities

Set \(A_0=2JC(J(A+D)+1)\) and \(L=\operatorname{lcm}(1,\ldots,A_0)\).
This bounds every nonzero \(|\Delta|\) and both positive
multipliers \(m c_s,j c_i\) by \(A_0\).
Let \(B_j\) be \(j!\) times the product of the full stripped
ramps of the fixed core in degree \(j\). A stripped core
coefficient, including its binomial coefficient, divides \(B_j\).

Take a prime \(p\) dividing \(x\), with \(E=v_p(x)\) and
\(L_p=v_p(L)\), and suppose \(E>L_p\). If \(\Delta\ne0\), the
identity
\[
                         (m c_s)F=(j c_i)x+\Delta
\]
gives \(v_p(F)=v_p(\Delta)-v_p(m c_s)\), between \(0\) and \(L_p\).
When \(\Delta=0\), it gives \(v_p(F)=E+v_p(j c_i)-v_p(m c_s)\),
within \(L_p\) of \(E\). Thus the coefficient ordinate at
order \(j\) is \(E\) times the number of its zero-offset
factors plus a remainder in
\[
                       [-r jD L_p,\ v_p(B_j)+r jD L_p].
\]
In order \(m'\) all zero-offset factors belong to omitted rows.
Their count lies between \(0\) and \(r_0\) at every index.
Therefore any Newton polygon slope \(q\) of that order satisfies
\[
                   |q|\le r_0 E+v_p(B_{m'})+2r m'D L_p .
\]

At marked order \(m\) every positive index has at least \(t-r_0\)
more zero-offset factors than index zero. After subtracting
\(q\) times the index, its ordinate exceeds the index-zero
ordinate by at least
\[
   \bigl(t-(m+1)r_0\bigr)E
          -v_p(B_m B_{m'}^m)-2rJ^2D L_p.
\]
Here \(m+1\le J\) and \(m+mm'<J^2\), so the count reserve
is at least \(3r/(4Q)\). This difference is strictly positive
whenever
\[
                 E>2Q\bigl(v_p(B_m B_{m'}^m)+2J^2D L_p\bigr).
\]
For such a prime no slope of the other polygon is a slope of
the marked polygon. Indeed index zero uniquely minimizes
the ordinates after subtracting any such slope, whereas a
polygon edge requires at least two minimizing indices.
By the Newton polygon product theorem a nonconstant common
rational factor would contribute a slope to both polygons.
[Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)

Consequently a common factor forces \(x\) to divide the fixed
integer
\[
                K=L^{4QJ^2D}\left(\prod_{2\le j<J}B_j^J\right)^{2Q}.
\]
For \(M\ge K+A+1\), all appended exponents are positive
and \(x>K\). This excludes every pair with terminal in
\([4,J)\). Terminals 2 and 3 use
[ordinary low-order coprimality](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
and terminals at least \(J\) use the absolute angular cutoff.
The original moment coefficients are integers, so rational
coprimality gives a Bézout identity over every
characteristic-zero field.
