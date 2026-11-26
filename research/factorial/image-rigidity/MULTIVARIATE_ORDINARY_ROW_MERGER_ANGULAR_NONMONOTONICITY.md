# Ordinary row identification can move root angles in opposite directions

Identifying integer exponent rows can increase some positive-half root
angles and decrease others in the same ordinary binomial moment
polynomial. This occurs with equal monomial total degrees and with
relative-gap curvature tending to zero. The construction works at
degree four and at every fixed \(j\ge300\) with \(j\equiv4\pmod6\).
At degree four the merger strictly increases every internal discrete
log curvature of the binomial-stripped coefficient array.

## The two integer supports

Fix a positive degree \(j\equiv4\pmod6\). For integers \(t\ge0\) put
\[
 R=3t+\frac74,\qquad s=\frac1R,\qquad
 r=\frac{8j}{3}R^2+\frac23R-\frac1{12},
 \qquad L_1=r+R-\frac12,\quad L_2=3r-R-\frac32.
\]
The split support consists of \((L_1,L_1+1)\), \((L_2,L_2+3)\),
and their two reversed rows. Merge the increasing rows and merge
the decreasing rows to obtain \((L_1+L_2,L_1+L_2+4)\) and its
reverse. Since
\[
 r=24jt^2+(28j+2)t+\frac{98j+13}{12}
       \in\mathbb Z+\frac34,
\]
both \(L_i\) are nonnegative integers for all sufficiently large
\(t\). The total degrees agree between the two monomials and are
preserved by the two mergers. For either profile,
\[
 S=\sum_{\rm rows}\frac{(b_i-a_i)^2}{\max(a_i,b_i)}
     =O_j(s^2).
\]
All assertions from the limiting calculation below hold along this
discrete sequence as \(t\to\infty\), with \(j\) fixed.

Use the ordinary factorial moment polynomial
\[
 F(z)=\sum_{k=0}^j\binom jk M(k) z^k,\qquad
 M(k)=\prod_i\frac{((j-k)a_i+kb_i)!}{(j a_i)!}.
\]
The profiles are mirrored, so \(M(k)=M(j-k)\). The
[mirrored root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md) gives
simple unit-circle roots. The roots used below also follow directly
from the real simple-root limit. We label the upper-half roots
by increasing argument in \((0,\pi)\).

## Exact lattice expansion

For a mirrored pair with lower exponent \(L\) and gap \(d\), the
factorial ratio groups into the Gauss shapes
\(\alpha=(jL+h)/d\), \(1\le h\le d\). Its contribution for each
shape is
\((\alpha)_k(\alpha)_{j-k}/(\alpha)_j\), where \((x)_k\)
denotes a rising factorial. The split and merged shape lists are
\[
\begin{array}{ll}
 \text{split:} & jr+jR-j/2+1,\quad
           jr-jR/3-j/2+h/3\quad(1\le h\le3),\\
 \text{merged:} & jr-j/2+h/4\quad(1\le h\le4).
\end{array}
\]
Write \(\alpha=A/s^2+B/s+C\), with \(A=8j^2/3\). Set
\(u=2j/3\) and \(g=-7j/12\). The merged shapes have \(B=u\)
and \(C=g+h/4\). In the split list the first shape has \(B=u+j\),
\(C=g+1\), and the three others have \(B=u-j/3\), \(C=g+h/3\).
In the following, \(\Delta\) denotes the split quantity minus
the merged quantity. Directly summing the four entries gives
\[
\begin{array}{lll}
 \Delta\sum B=0,& \Delta\sum C=1/2,&
 \Delta\sum B^2=A/2,\\
 \Delta\sum BC=u/2+j/3,&
 \Delta\sum B^3=4uj^2+8j^3/9,&
 \Delta\sum C^2=g+49/72,\\
 \Delta\sum B^2C=(4g/3+17/9)j^2,&
 \Delta\sum B^4=188j^4/27.&
\end{array}
\]
Put \(\sigma_p=\sum\alpha^{-p}\), counting all four shapes
of the indicated profile. Expansion at \(s=0\) yields
\[
\begin{split}
 \Delta\sigma_1&=A^{-3}(-g/2-67/144)s^6+O_j(s^7),\\
 \Delta\sigma_2&=(2A^3)^{-1}s^6+O_j(s^7),\\
 \Delta\sigma_3&=3(2A^4)^{-1}s^8+O_j(s^9),\\
 \Delta\sigma_p&=O_{j,p}(s^{2p+2}).
\end{split}                                                    \tag{1}
\]
For the first identity the coefficients at \(s^4,s^5\) vanish
because \(\Delta\sum B^2=A\Delta\sum C\) and
\(\Delta\sum B^3=2A\Delta\sum BC\). The coefficient at \(s^6\) is
\[
 A^{-3}\Delta\sum C^2
     -3A^{-4}\Delta\sum B^2C+A^{-5}\Delta\sum B^4.
\]
For general \(p\), equality of the counts and of the sums of \(B\)
cancels the first two terms. The coefficient at \(s^{2p+2}\) is
\[
 -pA^{-p-1}\Delta\sum C
   +\tfrac12p(p+1)A^{-p-2}\Delta\sum B^2
   =\frac{p(p-1)}{4A^{p+1}},
\]
which gives the other displayed leading coefficients.

Set \(b=k-j/2\). The exact logarithm of the multiplier is
\[
 \log M(k)=\sum_{p\ge1}\frac{(-1)^{p-1}}{p}\,\sigma_p P_p(k),
 \qquad
 P_p(k)=\sum_{h=0}^{k-1}h^p+\sum_{h=0}^{j-k-1}h^p
                      -\sum_{h=0}^{j-1}h^p .                         \tag{2}
\]
For fixed \(j\) this series converges when \(s\) is near zero,
by expanding the finitely many logarithms
\(\log(1+h/\alpha)\). Each \(P_p\) is an even polynomial in \(b\)
of degree at most \(p+1\), as follows from polynomial summation
and invariance under \(k\mapsto j-k\). In particular,
\[
 P_1=b^2-j^2/4,\qquad P_2=(j-1)P_1,\qquad
 P_3=b^4/2+\text{a polynomial of degree at most two in }b .
\]
Subtract the value of the logarithm at \(b=0\) from (2). This
multiplies \(F\) by a nonzero scalar. By (1), the leading split
minus merged quadratic and quartic terms in this logarithm are
\[
 s^6 B_2 b^2 +s^8 B_4 b^4,\qquad
 B_2=\frac{6j-31}{144A^3},\qquad B_4=\frac1{4A^4}.                    \tag{3}
\]
The \(b^2\) coefficient in (3) uses both \(p=1\) and \(p=2\).
The \(b^4\) coefficient uses \(p=3\). Terms from \(p\ge4\) have
order at least \(s^{2p+2}\) and degree in \(b\) at most \(p+1\).

## Cayley transform and the signed root comparison

Put \(w=(1+z)/(1-z)=isx\). For a degree at most \(j\) polynomial
\(G\), use the transform
\[
 {\cal T}_sG(x)=\frac{(1+isx)^j}{(2is)^j}
                    G\!\left(\frac{isx-1}{isx+1}\right).
\]
It sends \((1+z)^j\) to \(x^j\). Under this transform,
\(z\partial_z-j/2\) becomes the operator
\[
 {\cal D}_s=\frac{i}{2s}(1+s^2x^2)\partial_x
                   -\frac{ijsx}{2}
          =\frac{i}{2s}\partial_x+O_j(s)
\]
on the space of polynomials of degree at most \(j\). The leading
degree-\(j+1\) terms cancel, so this space is invariant. A term
\(s^d b^e\) in (2) therefore contributes \(O_j(s^{d-e})\) to the
operator norm in any fixed coefficient basis. In particular, the
errors beyond the \(s^6\) term for \(b^2\) and the \(s^8\) term for
\(b^4\) in (3) are \(O_j(s^5)\) after transport. For \(p\ge4\) the
transported difference contributions are \(O_{j,p}(s^{p+1})\).

To justify these estimates with the full series, truncate (2) at
\(p=j+5\). At the \(j+1\) points \(k=0,\ldots,j\) the tail,
with its midpoint value removed, is \(O_j(s^{2j+12})\).
Interpolate its values in \(b\) by a degree at most \(j\)
polynomial. Its coefficients are \(O_j(s^{2j+12})\) by the
fixed interpolation matrix, and transport costs at most \(s^{-j}\).
For the retained terms, one can substitute their full polynomial
expressions in \({\cal D}_s\), regardless of their degrees,
because the expressions agree with the original diagonal action
on each \(z^k\). These arguments also show analytic extension
of the transported operators at \(s=0\).

Let \(T_{\rm m}(s)\) and \(T_{\rm s}(s)\) be the operators obtained
from the midpoint-subtracted logarithms. They commute, since each
is a polynomial in \({\cal D}_s\) at \(s\ne0\). The normalized
Cayley polynomials are \(e^{T_{\rm m}(s)}x^j\) and
\(e^{T_{\rm s}(s)}x^j\), since \(z\partial_z-j/2\) is diagonal on
the original monomial basis. Define \(v=2/A=3/(4j^2)\).
For either operator, \(p=1\) gives the limit
\(-v\partial_x^2/2\) because \(\sigma_1=4s^2/A+O_j(s^3)\).
The \(p\ge2\) terms have transported order at least \(s^{p-1}\).
Consequently both Cayley polynomials tend to
\[
 H_j(x)=e^{-v\partial_x^2/2}x^j
       =v^{j/2}\operatorname{He}_j(x/\sqrt v).
\]
Their difference is
\[
 s^4\left(-\frac{B_2}{4}\partial_x^2+
                  \frac{B_4}{16}\partial_x^4\right)H_j(x)
            +O_j(s^5).                                               \tag{4}
\]
This follows from (3) and
\(e^{T_{\rm s}}-e^{T_{\rm m}}=(e^{T_{\rm s}-T_{\rm m}}-I)e^{T_{\rm m}}\).
All errors are coefficientwise at fixed \(j\). The transported
polynomials are real for real \(s\), since the midpoint-subtracted
logarithm is an even function on the symmetric \(b\) grid and
can be interpolated by an even real polynomial there.

For completeness, the generating function
\(\exp(Xu-u^2/2)\) of \(\operatorname{He}_j(X)\), in powers
\(u^j/j!\), gives
\(\operatorname{He}_{k+1}=X\operatorname{He}_k-k\operatorname{He}_{k-1}\)
and the differential equation
\(\operatorname{He}_j''-X\operatorname{He}_j'+j\operatorname{He}_j=0\).
The recurrence realizes the polynomial as the characteristic
polynomial of the real symmetric tridiagonal matrix with zero
diagonal and off-diagonal entries \(\sqrt{k}\), \(1\le k<j\).
Each eigenspace is one-dimensional because the first coordinate
determines the rest by the recurrence. Hence the roots are simple.
For even \(j\) they occur in opposite nonzero pairs, by the
recurrence at zero and the parity of \(\operatorname{He}_j\).

At a positive Hermite root \(X\) the two nearby Cayley roots tend
to \(\sqrt{v}X\). Subtract the root equations in (4) and divide
by the difference quotient of the merged polynomial, which
tends to \(H_j'(\sqrt v X)\ne0\). To control the error to order
\(s^5\), both analytic root branches are
\(\sqrt v X+O_j(s)\), by the simple root theorem. The differential
equation gives the derivative ratios
\(\operatorname{He}_j''/\operatorname{He}_j'=X\) and
\(\operatorname{He}_j^{(4)}/\operatorname{He}_j'
=X^3+(3-2j)X\) there. The split root minus the merged root in
the \(x\) coordinate is therefore
\[
 \frac{s^4 X B_4}{16v^{3/2}}\bigl(\Theta-X^2\bigr)+O_j(s^5),
 \qquad
 \Theta=2j-3+\frac{4v B_2}{B_4}=\frac{10j}{3}-\frac{89}{9}.
                                                                    \tag{5}
\]
On the upper half circle the coordinate is
\(s x=\cot(\arg z/2)\), a strictly decreasing function of
the angle. Thus roots with \(X^2<\Theta\) and \(X^2>\Theta\)
give opposite angular comparisons.

The smallest positive \(X^2\) is at most \(j-1\), since the
average squared eigenvalue of the tridiagonal matrix is \(j-1\),
by the trace of its square. The largest \(X^2\) is at least
\((19/10)^2(j-20)\) when \(j\ge20\): use the unit Rayleigh vector
constant on the last twenty indices and zero on the others,
where each of the nineteen included off-diagonal entries is at
least \(\sqrt{j-20}\). For \(j\ge300\),
\[
 j-1<\Theta<\frac{361}{100}(j-20).
\]
At \(j=4\) the two squared roots are \(3-\sqrt6\), \(3+\sqrt6\),
and \(\Theta=31/9\) lies between them. This proves the opposite
angular motions in both claimed ranges.

At fixed \(k\) the quadratic term of (3) is the whole change
through order \(s^6\), up to the scalar. Taking the second
difference in \(k\) gives \(2B_2s^6+O_j(s^7)\) at every internal
position. At \(j=4\), \(B_2<0\), so the merge increases all three
discrete log curvatures. At \(j\ge300\), \(B_2>0\), so the merge
decreases them. Both uniform curvature comparisons coexist with
the opposite angular motions.

## Explicit quartic and scope

At \(j=4,t=0\) the lower exponents are \(35,98\) and the merged
lower exponent is \(133\). The accompanying
[exact replay](verify_multivariate_ordinary_row_merger_angles.py)
checks this quartic with rational arithmetic, and checks the
lattice expansions and comparison signs at sample orders. The
fixed-order analytic proof gives the unbounded-order conclusion.

The comparison concerns two supports at the same order. It rules
out a universal row-merger law that orders all positive-half root
angles in one direction, even with pointwise discrete curvature
ordering. Consecutive coprimality for either individual support
is a different assertion. The calculation gives no common root
of two adjacent ordinary factorial moment polynomials.

