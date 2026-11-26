# Monic remainders and a simple resultant factor for increasing factorial rows

Let every unequal row of an ordinary binomial be increasing,
with entries \((a_i,a_i+d_i)\), where \(d_i\ge1\) is integral.
Strip the nonzero smaller-end factorials and write the remaining
order-\(j\) polynomial and its monic reversal as
\[
 Q_j(z)=\sum_{k=0}^j\binom jk z^k
           \prod_i\prod_{h=1}^{kd_i}(ja_i+h),
 \qquad Y_j(x)=x^j Q_j(1/x).
\]
In this note the \(a_i\) are also used as indeterminates.
For \(n\ge1\) put \(N=n+1\).

**Theorem.** If there is at least one unequal row, every
coefficient of the monic polynomial remainder
\(Y_N\bmod Y_n\) is a polynomial in the heights with nonnegative
coefficients and strictly positive constant term. The same holds
in the elementary symmetric generators within each equal-gap
group. It also holds for the continuation obtained by replacing
\(Q_N(z)\) by \(Q_N((1+t)z)\), with \(t\) an additional indeterminate.

**Proposition.** For two increasing unit-gap rows, the next
pseudo-remainder has coefficients of both signs at
\(a_1=a_2=0\) for every \(n\ge7\). A sign normalization
of the ordinary monic pseudo-Euclidean chain cannot give positive
coefficients throughout.

**Theorem.** For two increasing unit-gap rows write
\(S=a_1+a_2,\ P=a_1a_2\). At every \(n\ge1\), the resultant
\(\operatorname{Res}_z(Q_n,Q_N)\) in \(\mathbb Z[S,P]\)
has the factor \(1+S+P\) of multiplicity exactly one.

## The first subtraction

Work in the polynomial ring in the heights, or in the
elementary symmetric generators of the heights within
each group of equal gaps. Use \(\ge\) for coefficientwise
inequality. Put
\[
 A_h=\prod_i\prod_{\ell=1}^{d_i}(na_i+(h-1)d_i+\ell),
 \qquad
 B_h=(1+t)\prod_i\prod_{\ell=1}^{d_i}(Na_i+(h-1)d_i+\ell),
 \qquad \delta_h=B_h-A_h .
\]
All \(A_h,\delta_h\) have nonnegative coefficients, and
both sequences increase coefficientwise. For the symmetric
version, first group the factors with the same \((d,\ell)\).
The coefficient of the elementary generator of index \(r\) in
each group factor is \(j^r((h-1)d+\ell)^{m-r}\), with group
size \(m\). In the product, each term of fixed total height
degree \(r\) therefore has scalar \(j^r\) times a nonnegative
increasing polynomial in \(h-1\). Taking the difference of orders,
and including the extra \(t\) term, gives the same inequalities.

Write \(A^{[k]}=\prod_{h=1}^k A_h\) and \(B^{[k]}=\prod_{h=1}^k B_h\).
The coefficient of \(x^{n-k}\) in \(Y_n\) is
\(a_k=\binom nk A^{[k]}\), with \(a_k=0\) for \(k>n\).
Let \(b_k=\binom Nk B^{[k]}\) for the upper order.
The first quotient is \(x+A_1+N\delta_1\). For \(2\le k\le N\)
the remainder coefficient at \(x^{N-k}\) is exactly
\[
 e_k=
 \binom n{k-1} A^{[k-1]}(A_k-A_1)
 +\binom Nk\bigl(B^{[k]}-A^{[k]}-k\delta_1 A^{[k-1]}\bigr).     \tag{1}
\]
This uses Pascal's identity and \(k\binom Nk=N\binom n{k-1}\).
In the telescoping expansion of \(B^{[k]}-A^{[k]}\),
the term at position \(h\) is \(\delta_h\) times \(k-1\)
factors. Those factors coefficientwise dominate the \(A\)
factors with index \(h\) omitted. Since \(A_1\le A_2\le\cdots\),
the omitted product is at least \(A^{[k-1]}\).
Each term dominates \(\delta_1 A^{[k-1]}\).
The first summand of (1) has strictly positive constant
term because \(k\ge2\) and there is an unequal row.

## Where coefficient positivity breaks in the chain

Set both heights and \(t\) to zero in the two-row unit-gap
case. Then \(A_h=B_h=h^2\), \(a_k=\binom nk(k!)^2\), and
\(e_k=(k^2-1)a_{k-1}\). If \(E=Y_N\bmod Y_n\), the second
pseudo-remainder is
\[
             F=e_2^2 Y_n-(e_2 x+e_2a_1-e_3)E .
\]
Its coefficient at \(x^{n-k}\) is
\(e_2^2a_k+(e_3-e_2a_1)e_{k+1}-e_2e_{k+2}\), with missing
\(e\)'s set to zero. The coefficient at \(x^{n-2}\), divided
by \(a_2\), is \(2n(71-11n)\). In fact
\(e_2=3n,\ e_3=16n(n-1)\), and \(e_4/a_2=45(n-2)\).
For \(n\ge7\) it is negative. The constant term divided by
\(a_n\) is
\[
                  e_2^2+(e_3-e_2a_1)((n+1)^2-1)>0 ,
\]
since \(e_3-e_2a_1=n(13n-16)\). Thus even at one ordinary
profile the second step has mixed coefficient signs.

## Exact multiplicity of the forced symmetric factor

Write the two heights as \(a,b\). On \(a=-1\) the order-\(j\)
polynomial drops from degree \(j\) to degree \(j-1\), since
\((1-j)_j=0\). The monic \(Y_j\) then has a simple root at
zero over \(\mathbb Q(b)\). To show that this is the only
common root at \(a=-1\), it suffices to show that the two
degree-reduced polynomials \(Q_n,Q_N\) there are coprime at
some \(b\). As \(b\to+\infty\), \(Q_j(z/b)\) at \(a=-1\)
tends to
\[
             L_j^*(z)=\sum_{k=0}^{j-1}\binom jk(1-j)_k j^k z^k .
\]
Reversal in degree \(j-1\) gives
\[
 x^{j-1}L_j^*(1/x)=(-1)^{j-1}(j-1)!\,j^{j-1}
                                 L_{j-1}^{(1)}(x/j).             \tag{2}
\]
This follows termwise from the coefficient expansion
\(L_m^{(1)}(x)=\sum_{r=0}^m\binom{m+1}{r+1}(-x)^r/r!\).

To check coprimality in (2), let \(\lambda_{m,k}\) be the
\(k\)-th increasing zero of \(L_m^{(1)}\). Laguerre
orthogonality for weight \(x e^{-x}\) gives simple positive
zeros and strict unscaled interlacing. Also
\((m+1)\lambda_{m,k}\) is strictly decreasing in \(m\) for
each fixed \(k\). Indeed the differential equation
\(xL''+(2-x)L'+mL=0\), also immediate from the coefficients,
gives for \(w_m(y)=y e^{-y/(2(m+1))}L_m^{(1)}(y/(m+1))\)
\[
              w_m''+\left(\frac1y-\frac1{4(m+1)^2}\right)w_m=0 .
\]
The potential increases strictly with \(m\) on \(y>0\).
For completeness, take two consecutive such equations.
On the interval from zero to the first old positive root,
if the new solution has no interior zero, choose both signs
positive. The Wronskian \(w_{\rm old}w_{\rm new}'-
w_{\rm old}'w_{\rm new}\) starts at zero and strictly decreases,
while its value at the right endpoint is nonnegative.
The same argument on any interval between consecutive old roots
starts with nonpositive Wronskian and ends with nonnegative
Wronskian if the new solution has no interior zero.
This proves the strict decrease of each scaled zero.
The other side of interlacing follows from
\(\lambda_{m,k}<\lambda_{m+1,k+1}\) and positivity.
Hence
\[
           (m+2)\lambda_{m+1,k}<
           (m+1)\lambda_{m,k}<
           (m+2)\lambda_{m+1,k+1}.
\]
At \(n=1\) the lower limit is constant. Thus both limits
in (2) are coprime at every start, and continuity of the
degree-reduced resultant proves coprimality over \(\mathbb Q(b)\).

The zero branch of \(Y_j(x)\) at \(a=-1\) has derivative
\(-j(b+1)\) with respect to \(a\). Indeed the derivative of
the constant coefficient is \(j(1-j)_{j-1}(jb+1)_j\), while
the coefficient of \(x\) is \(j(1-j)_{j-1}(jb+1)_{j-1}\).
The two zero branches therefore meet to exactly first order
over \(\mathbb Q(b)\); all other cross-root differences are
units there. The product formula for the monic resultant
gives multiplicity one on \(a=-1\). By symmetry the same
holds on \(b=-1\). Reversal of both polynomials multiplies
the resultant by \((-1)^{nN}=1\). The product of the two
simple factors is \((1+a)(1+b)=1+S+P\), and division by it
has an integral symmetric quotient. This gives the stated
factor with exact multiplicity in \(\mathbb Z[S,P]\).
