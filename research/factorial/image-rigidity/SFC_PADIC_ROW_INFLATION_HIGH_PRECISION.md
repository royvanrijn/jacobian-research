# Higher-precision factorial row inflation and packing

Let \(p\ge5\) be prime. Inflating any ordinary exponent row by a
power of \(p\) preserves all moment coefficient valuations and
their normalized units modulo \(p^3\), after one positive diagonal
rescaling of the monomial coefficients. The comparison holds
simultaneously at every order and for arbitrary support size.

More generally, replacing a row already inflated by \(p^s\) with
the same row inflated by \(p^e\), \(e>s\), gives precision
\(p^{3(s+1)}\). For a pair of rows whose entries add to zero modulo
\(p^L\) in each column, simultaneous inflation from \(s\) to \(e\)
improves this precision to \(p^{3(s+1)+L}\).
These statements retain higher-order arithmetic certificates for
the corresponding original and inflated profiles.

The row packing operation in the
[first initial-form comparison](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md)
also admits this stronger precision. At the \(p\)-free factorial
level, the argument uses the block decomposition underlying
[Granville's prime-power binomial congruences](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/granville.pdf),
with an explicit odd function controlling the difference of units.

## An odd factorial block correction

For integers \(M\ge0\) put
\[
                F_p(M)=\prod_{\substack{1\le u\le M\\p\nmid u}}u .
\]
Fix \(h\ge1\), set \(P=p^h\), and define
\[
                 U_h(x)=F_p(Px)/F_p(P)^x,\qquad x\ge0
\]
for integer \(x\). Pairing each residue \(u\) with \(P-u\)
in the block from \(jP+1\) to \((j+1)P\) gives the exact
normalized block
\[
 W_h(j)=\prod_{\substack{1\le u<P/2\\p\nmid u}}
         \left(1+\frac{P^2j(j+1)}{u(P-u)}\right),
       \qquad U_h(x)=\prod_{j=0}^{x-1} W_h(j).
 \tag{1}
\]
The first reciprocal sum satisfies
\[
 \sum_{\substack{1\le u<P/2\\p\nmid u}}\frac{1}{u(P-u)}
          \equiv-\tfrac12\sum_{\substack{1\le u<P\\p\nmid u}}
                      u^{-2}\equiv0\pmod P.
\]
For the last congruence multiply the full unit group by \(2\).
This permutes its elements and multiplies the sum of reciprocal
squares by \(2^{-2}\). The difference \(1-2^{-2}\) is a
\(p\)-unit. In particular \(W_h(j)\) and \(U_h(x)\) are one
modulo \(p^{3h}\).

The normalized logarithm
\(p^{-3h}\log U_h(x)\) has an odd extension to \(\mathbb Z_p\)
that is \(1\)-Lipschitz in the \(p\)-adic norm. To prove this,
use the ordinary convergent \(p\)-adic logarithm on
\(1+p\mathbb Z_p\). Define
\[
 c_{h,r}=\frac{(-1)^{r-1}P^{2r}}{r}
           \sum_{\substack{1\le u<P/2\\p\nmid u}}
                          [u(P-u)]^{-r},\qquad r\ge1.
\]
The coefficient \(c_{h,1}\) has valuation at least \(3h\);
the other coefficients have valuation at least
\(2hr-v_p(r)\). Additivity of the logarithm and (1) give
\[
                \log U_h(x)=\sum_{r\ge1}c_{h,r}Q_r(x),
 \qquad Q_r(x)=\sum_{j=0}^{x-1}[j(j+1)]^r
                                                     \tag{2}
\]
at nonnegative integers. The polynomial \(Q_r(X)\) extends
this discrete sum to all \(X\). It is an integer linear
combination of \(\binom Xk\) for \(1\le k\le2r+1\): expand
\((X(X+1))^r\) in the binomial basis by forward differences
at zero, and sum using
\(\sum_{j=0}^{x-1}\binom jk=\binom x{k+1}\).
The polynomial \(Q_r\) is odd, because its difference
\([X(X+1)]^r\) is invariant under \(X\mapsto-1-X\), and
\(Q_r(0)=0\).

For \(x,y\in\mathbb Z_p\), Vandermonde's identity gives
\[
 v_p\!\left(\binom xk-\binom yk\right)
             \ge v_p(x-y)-\lfloor\log_p k\rfloor .
\]
Indeed each nonconstant term in the binomial sum for
\(x=y+(x-y)\) includes
\(\binom{x-y}{i}=\frac{x-y}{i}\binom{x-y-1}{i-1}\), and
binomial polynomials take \(p\)-integral values on \(\mathbb Z_p\)
by continuity from nonnegative integers.
For \(r=1\) one can use the exact
\(Q_1(X)=(X^3-X)/3\), which loses no precision for \(p\ge5\).
For \(r\ge2\), the loss satisfies
\[
             v_p(r)+\lfloor\log_p(2r+1)\rfloor\le2r-3.
\]
At \(r=2\) this is immediate. For \(r\ge3\) use
\(v_p(r)\le r-1\) and \(5^{r-1}>2r+1\) to bound the two terms.
Thus every term in (2) has differences divisible by
\(p^{3h+v_p(x-y)}\). The series converges uniformly on
\(\mathbb Z_p\), since \(Q_r\) has integral values and
\(2hr-v_p(r)\to\infty\). It is odd by the polynomial
identities term by term and proves the claimed Lipschitz
property. In particular, for nonnegative integers \(x,y\),
\[
           x+y\equiv0\pmod{p^L}
    \quad\Longrightarrow\quad U_h(x)U_h(y)\equiv1\pmod{p^{3h+L}}.
 \tag{3}
\]
To pass back from logarithms use the \(p\)-adic exponential
on \(p\mathbb Z_p\), where the valuation of
\(\exp(z)-1\) is that of \(z\). The original block units
already belong to this logarithm's domain.

## Factorial coefficients at all orders

Separate the \(p\)-divisible factors successively in a factorial.
With
\[
 \Lambda_e=p^{1+p+\cdots+p^{e-1}}\prod_{h=1}^e F_p(p^h),
                 \qquad \Lambda_0=1 ,
\]
one finds
\[
                         (p^e x)!=x!\Lambda_e^x
                                             \prod_{h=1}^e U_h(x).
 \tag{4}
\]
The empty product at \(x=0\) or \(e=0\) causes no exception.
For a fixed exponent row \(A=(A_i)\) and coefficient index
\(\nu=(\nu_i)\), its factorial argument before inflation is
\(x=\sum_i\nu_i A_i\). The factor
\((\Lambda_e/\Lambda_s)^x\) in comparing (4) for \(e\) and
\(s\) is removed by replacing each coefficient variable
\(c_i\) by \((\Lambda_e/\Lambda_s)^{-A_i}c_i\).
This is the same positive diagonal rescaling at every order.
After it, the comparison multiplies each coefficient by
\(\prod_{h=s+1}^e U_h(x)\equiv1\pmod{p^{3(s+1)}}\).
Multiple rows can be inflated independently using the product
of the indicated diagonal rescalings.

If two rows \(A,B\) are inflated together and \(A_i+B_i\)
is divisible by \(p^L\) for every \(i\), their two uninflated
factorial arguments add to a multiple of \(p^L\) at every
coefficient index. Equation (3) at each layer now proves the
stronger precision \(p^{3(s+1)+L}\).

For a binary pair at any two fixed orders, subtract the
minimum coefficient valuation of each polynomial by dividing
by a power of \(p\). After the stated diagonal rescaling,
the corresponding minima agree, and the integral polynomials
have congruent coefficients at the stated precision.
Their Sylvester resultants are therefore congruent at that
precision. If the original normalized resultant has valuation
smaller than the precision, the inflated pair is coprime.
This also preserves every higher-precision computation that
uses only the indicated normalized coefficients.

## Packing past a bounded residual row

Take existing rows \(p^s A_i\) and \(B_i\) and replace them by
one row \(p^eA_i+B_i\), where \(e>s\). At any coefficient index
write \(y=\sum_i\nu_i B_i\). After separating the inflated factor,
the remaining comparison is the exact product
\[
     \frac{(p^e x+y)!}{(p^e x)!\,y!}
                   =\prod_{i=1}^{y}(1+p^e x/i).
\]
For moment orders at most \(J\), let \(H=J\max_i B_i\).
When \(H>0\) and \(e>\lfloor\log_p H\rfloor\), the product is
one modulo \(p^{e-\lfloor\log_p H\rfloor}\), since \(y\le H\).
If \(H=0\), the product is one exactly.
Thus the packed row preserves coefficient valuations and
units, after the same diagonal rescaling, with exponent
of precision
\[
                 \min\{3(s+1),\ e-\lfloor\log_p H\rfloor\},
\]
where the second term is omitted if \(H=0\).
In particular, to carry prescribed precision \(p^r\) through
packing, it suffices to start with \(3(s+1)\ge r\) and then
separate the short row by choosing
\(e\ge r+\lfloor\log_p H\rfloor\), with \(e>s\).
The argument requires no divisibility of the original entries
\(A_i\). For two complementary high rows with separate short
rows, apply the same packing comparison to both, using
\(3(s+1)+L\) in place of \(3(s+1)\) for their joint
inflation correction.

The estimates transport established arithmetic certificates.
A residue collision at the supported precision can still occur
in either row representation.

