# The finite spectral pole of the degree-forty-two Weyl residual

**Theorem PolynomialWeylFortyTwoFinitePole.** In the hypothetical
[degree-forty-two residual](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md)
over \(\mathbb C\), translate the Fourier coordinate \(X\) to put the
leading differential coefficients of \(P,Q\) in the form
\[
       \operatorname{lc}_Y(P)=\alpha X^{12},\qquad
       \operatorname{lc}_Y(Q)=\beta X^{42}.
\]
The translation preserves the Fourier rectangle and the incoming
\((3,-1)\) face. In particular \(0\) is the only possible finite
singularity of either spectral equation.

Put \(C=XY\), and grade the polynomial Weyl algebra by \(i-j\) on
\(X^iY^j\). The lowest components are \(X^{-2k}f(C)\) and
\(X^{-7k}g(C)\), with polynomial degrees \(2\ell,7\ell\), respectively,
where
\[
       (k,\ell)=(9,15)\quad\text{or}\quad 10\le k\le\ell\le14.      \tag{1}
\]
The two finite spectral equations are regular singular at \(X=0\)
exactly for \(k=9\). For every spectral value, their formal Taylor
solution spaces there have dimensions at least \(2k,7k\), respectively.
In the regular singular case these give at least \(18,63\) independent
entire solutions, respectively.

The [regular finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_FINITE_POLE.md)
eliminates the alternatives with \(k\ge10\) in (1). The
[spectral monodromy theorem](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_SPECTRAL_MONODROMY.md)
then gives exactly \(18\) and \(63\) entire solutions, respectively.

After dividing \(f,g\) by their leading scalars and by the falling
factorials \((C)_{2k},(C)_{7k}\), the remaining roots of \(f\)
split into \(\ell-k\) pairs, each of difference \(mk\) with
\(m\in\{1,3,5,7\}\). The remaining roots of \(g\) are determined
from these pairs by (4). In the regular singular case the sum of the
paired excess roots of \(f\) has a direct formula in the
infinity coefficient \(P_4\) given in (5).

## Centering the leading coefficients

The Fourier rectangle has corners \((12,30),(42,105)\) in the
order \((X,Y)\), so the differential orders in \(Y\) are exactly
\(30,105\). The
[leading-factor residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
gives a common leading factor \(h(X)\) of degree six, with
\(\operatorname{lc}_Y(P)=\alpha h^2\) and
\(\operatorname{lc}_Y(Q)=\beta h^7\), since the order gcd is
fifteen. If \(h\) had two distinct roots, its degree would be
at least twice fifteen. Thus \(h\) is proportional to
\((X-X_0)^6\), and translating \(X\) gives the displayed form.
Indeed this translation only lowers the \(X\) exponent of each
changed monomial, so preserves the rectangle, the Bernstein top
and the \((3,-1)\) face. The maximal \(X\)-exponent faces were
already single corners. The translated maximal \(Y\)-exponent
faces are now single corners as well.

The notation \(t=X^{1/3}\), \(E=tY\) and \(P=\sum t^iP_i(E)\)
in the
[infinity coefficient constraints](POLYNOMIAL_WEYL_FORTY_TWO_SECONDARY_ROOT_CONSTRAINTS.md)
can be used with this centered coordinate. In its normal form,
\[
 P_i(E)\text{ has degree at most }29\text{ for }i\ne6,\qquad
 Q_i(E)\text{ has degree at most }104\text{ for }i\ne21.          \tag{2}
\]
Only the single maximal \(Y\)-exponent corner can give degree
\(30\) or \(105\) in \(E\), respectively. The normal-order
corrections each lower \(E\)-degree. In particular the
[lower scalar coefficient polynomials](POLYNOMIAL_WEYL_FORTY_TWO_SCALAR_FACE_COEFFICIENT_DESCENT.md)
\(C_3,C_0\) can be taken of degrees at most four and nine,
respectively, in this presentation. Their terms have \(E\)-degree
divisible by three, so their largest possible \(E\)-degree is
now \(27\).

## The lowest Euler components at zero

Use weights \((-s,1)\) on \((X,Y)\) for \(0\le s<1\). The
weights of \(P,Q\) are at least \(30-12s\) and \(105-42s\),
from their corner terms. Their sum exceeds the contraction
loss \(1-s\). Thus the leading Poisson bracket is zero, as
the actual commutator is scalar. Away from the finitely many
edge weights, the leading symbols are monomials. Their
positive weights and their zero bracket require that their
exponent vectors are proportional. The ratio of the weights
is locally constant on these open intervals and continuous at
their endpoints. At \(s=0\) it is \(2:7\), so this ratio holds
through \(s=1\).

The weight \((-1,1)\) is the negative of the exact Euler
grading. It gives maxima \(2k,7k\) with \(k\) an integer
because the ratio is \(2:7\) in integral weights. The
corner gives \(k\ge9\). At a weight sufficiently close to
\(s=1\) from below, the maximizing monomials are those
of largest \(X\)-exponent among the lowest Euler components.
They also have largest \(Y\)-exponents in those components.
Their exponent vectors are proportional in ratio \(2:7\) and
are therefore \(2(\ell-k,\ell)\) and \(7(\ell-k,\ell)\) for an
integer \(\ell\). If \(k=9\), the corner lies in the component
and has the largest possible \(Y\)-exponent, giving \(\ell=15\).
If \(k>9\), the corner lies outside the component and no other
term of maximal \(Y\)-exponent exists. Hence \(\ell\le14\).
Nonnegative \(X\)-exponents give \(\ell\ge k\), proving (1).

Write \((C)_j=C(C-1)\cdots(C-j+1)\), so
\(X^iY^j=X^{i-j}(C)_j\). In grade \(-2k\) the polynomiality
condition \(i\ge0\) gives \(j\ge2k\), which forces the factor
\((C)_{2k}\) in \(f\). Similarly \(g\) has the factor
\((C)_{7k}\). The lowest components must commute, since their
product grade \(-9k\) is negative and the bracket is scalar.
From \(a(C) X^j=X^j a(C+j)\), their equation is
\[
            f(C-7k)g(C)=f(C)g(C-2k).
\]
Write the monic remaining factors after removing the falling
factorials and leading scalars as \(F,G\).
The falling factorial products cancel since
\((C-7k)_{2k}(C)_{7k}=(C)_{2k}(C-2k)_{7k}\).
We obtain
\[
               F(C-7k)G(C)=F(C)G(C-2k).                    \tag{3}
\]
In a root orbit \(\rho+k\mathbb Z\) let \(F_\rho(z),G_\rho(z)\)
record multiplicities by the integer index in this progression.
They are Laurent polynomials with nonnegative coefficients, and
(3) gives
\[
           G_\rho(z)=\frac{1-z^7}{1-z^2}\,F_\rho(z).        \tag{4}
\]
At \(z=-1\) the divisibility gives an equal number of even and
odd indices in \(F_\rho\), counted with multiplicity.
Expanding (4) from its lowest index, the coefficient at \(z^j\)
is the number of roots of parity \(j\) at indices at most \(j\),
minus the number of roots of the other parity at indices at
most \(j-7\). Let the sorted multisets of even and odd indices
be \(e_i,o_i\). Applying these inequalities at \(j=o_i+7\) and
\(j=e_i+7\) gives \(e_i\le o_i+7\) and \(o_i\le e_i+7\).
Thus the matched differences are odd and at most seven.
Conversely, for a pair \(z^a(1+z^m)\) with \(m\) odd and
\(1\le m\le7\), the right side of (4) has only nonnegative
coefficients. After subtracting \(a\) from its indices,
its support consists of the even progression \(0,2,\ldots,m+5\)
and the odd progression \(m,m+2,\ldots,5\), with the latter
empty for \(m=7\). This also determines the remaining
\(g\) roots from the matched pairs. The degree of \(F\) is
\(2(\ell-k)\), proving the count in the theorem.

## What the finite equation shares with the infinity coefficient

For \(P-\lambda\), regular singularity at \(X=0\) is equivalent
to requiring that every coefficient \(p_j(X)\) of \(Y^j\)
has vanishing order at least \(12-(30-j)=j-18\).
Equivalently \(X^{18}(P-\lambda)\) has coefficients holomorphic
at zero in powers of \(C\), and its coefficient of \(C^{30}\)
is nonzero there. This is the \(k=9\) condition just proved.
For \(Q-\lambda\), the corresponding inequality is \(j-63\)
and is also equivalent to \(k=9\).

For the Taylor solution count write
\[
 X^{2k}(P-\lambda)=f(C)+\sum_{d\ge1} X^d f_d(C)
\]
with finitely many \(f_d\). For a formal Taylor solution
\(\sum_{n\ge0}a_nX^n\), the equation at exponent \(n\) has
diagonal coefficient \(f(n)a_n\), with the remaining terms
using strictly earlier \(a_j\). The equations below exponent
\(2k\) are identically zero because \(P-\lambda\) preserves
formal Taylor series. Each subsequent zero \(f(n)=0\)
imposes at most one constraint and adds one free \(a_n\),
giving dimension at least \(2k\). The argument for \(Q\)
is identical with \(7k\).
If \(k=9\), the degrees of \(f_d\) are bounded by the degree of
\(f\). Hence for \(n\) sufficiently large the finitely many
ratios \(f_d(n-d)/f(n)\) are bounded. The recursion then bounds
\(\lvert a_n\rvert\) by a geometric sequence, giving convergence
near zero. There are no other finite singularities. Analytic
continuation of each such germ around zero is single-valued
because the germ is holomorphic there, so the germ extends to
an entire solution.

Finally suppose \(k=9\). Let \(s_1,\ldots,s_{12}\) be the
paired excess roots of \(f\), and let
\(B(z)\) be the polynomial in the infinity coefficient
\(P_4(E)=E^2(E^3-u)^6 B(E^3)\).
The linked trace identity is
\[
               \sum_{j=1}^{12} s_j
                  =137-\frac{[z^3]B(z)}{\alpha}.          \tag{5}
\]
To see this, let \(p_{29,11}\) be the coefficient of
\(X^{11}Y^{29}\) in \(P\).
The coefficient of \(C^{29}\) in the regular indicial
polynomial is \(p_{29,11}-435\alpha\), using
\((C)_{30}=C^{30}-435C^{29}+\cdots\).
The roots of \((C)_{18}\) sum to \(153\), so the sum of
the excess roots is \(282-p_{29,11}/\alpha\).
In the conversion \(Y=t^{-1}E\), the first normal-order
correction to \(Y^n\) is
\(-n(n-1)t^{-n-2}E^{n-1}/6\), obtained by commuting
\([E,t^{-1}]=-1/(3t^3)\) at each pair of positions.
Consequently the coefficient of \(t^4 E^{29}\) in \(P\) is
\(p_{29,11}-145\alpha\), since the only other possible source
is \(\alpha X^{12}Y^{30}\). This coefficient is
\([z^3]B(z)\), proving (5). The equation uses the same
\(B\) as the simple-root condition
\(B(v)=(A_0+A_1v)^2/(4\alpha)\) at infinity.

The complete coefficient map, beyond the trace in (5), follows
from \(C=t^2E\) and \([C,t]=t/3\). Write
\([C]_{j,2/3}=\prod_{a=0}^{j-1}(C-2a/3)\).
Then \(t^iE^j=t^{i-2j}[C]_{j,2/3}\), so
\[
 \begin{aligned}
 f(C)&=\sum_{j\ge0} [E^j]P_{2j-6k}(E)\,[C]_{j,2/3},\\
 g(C)&=\sum_{j\ge0} [E^j]Q_{2j-21k}(E)\,[C]_{j,2/3}.
 \end{aligned}                                                 \tag{6}
\]
The indices select exactly the finite Euler grades, since
\(t^{i-2j}=X^{(i-2j)/3}\).
For \(k=9\) the highest four coefficients of \(f\) in this
basis are
\[
   \alpha,\quad [z^3]B(z),\quad
   [z^5]C_2(z),\quad [z^9]C_0(z).
\]
Here \(C_2,C_0\) are the polynomials from the lower scalar
coefficient descent, with the sharpened degree of \(C_0\)
from (2). Thus the entire paired-factor constraint on the
indicial polynomial applies to the infinity Laurent data
with their simple-root and fourth-multiplicity-root
conditions. If \(k\ge10\), the largest possible \(P\) index
used by (6) is \(4\ell-6k\le-4\).
