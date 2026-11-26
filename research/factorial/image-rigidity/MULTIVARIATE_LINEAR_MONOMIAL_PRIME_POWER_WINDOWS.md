# Prime-power consecutive windows for one shared unit gap and private prime-power exponents

Let \(p\) be prime and take the ordinary binomial with rows
\((1,2)\) and \((0,e_i)\), \(1\le i\le r\), where each
\(e_i=p^{t_i}\) for an integer \(t_i\ge0\).
Set \(M=r+1\). The two moments at each terminal \(N=p^s\),
\(s\ge1\), are coprime in the physical variable in each of
the following cases.

- There is one private row.
- There are at least two private rows, \(p\nmid M\), and
  \(\mathcal R_{p,M}(Y)=\sum_{j=0}^{p-1}(j!)^M Y^j\) is squarefree
  over \(\mathbb F_p\).

For example, the second condition holds if \(M\equiv1\pmod{p-1}\)
and \(p\nmid M\). At \(p=2\) it applies to every positive even
number of private rows. These classes allow arbitrarily many
private rows and arbitrarily large exponents. There is also a
single-private extension to shared row \((1,3)\), with
\(p=2\) and \(e_1=2^t\).

## The intermediate coefficient valuations

The [linear-monomial Euler bridge](MULTIVARIATE_LINEAR_MONOMIAL_EULER_BRIDGE.md#the-intermediate-polynomial-and-windows)
reduces coprimality to squarefreeness of the reciprocal intermediate
polynomial
\[
              J(z)=\sum_{k=0}^N c_k z^k,\qquad
              c_k=\binom Nk(N)_k\prod_{i=1}^r(e_i k)!.
 \tag{1}
\]
Put \(L=\sum_i(e_i-1)/(p-1)\) and write \(S_p(k)\) for the base-\(p\)
digit sum. Multiplication of \(k\) by any \(e_i\) preserves that
digit sum. For \(1\le k\le N\),
\[
 v_p(c_k)= Lk+2(s-v_p(k))+M v_p(k!),\qquad v_p(c_0)=0.
 \tag{2}
\]
Indeed \(v_p(N\pm j)=v_p(j)\) for \(1\le j<N\), so the product
\((N)_{\underline k}(N)_k\) contributes
\(2s+2v_p((k-1)!)\). The private factorials supply the
remaining \(k!\) powers and the linear term.

We use the following Newton polygon rule. An edge of reduced
slope denominator \(b\) has its supported indices congruent
modulo \(b\). Normalize its coefficients by their valuations
and read the resulting polynomial in the lattice variable,
whose exponents are index differences divided by \(b\).
If all such residual polynomials are squarefree, the original
polynomial is squarefree in characteristic zero. In fact the
initial polynomials of factors over \(\mathbb Q_p\) multiply
after a common fractional scaling, and a squared factor with
that slope supplies a squared factor of the lattice residual.
This uses no separability assertion after replacing the lattice
variable by its \(b\)-th power.
[Weighted initial-product argument](ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex)

## One private row

For \(M=2\), (2) has a single edge from zero to \(N\), with no
interior contacts. Subtracting the endpoint chord gives
\[
         2\left(s-v_p(k)-\frac{S_p(k)}{p-1}
                                  +\frac{k}{N(p-1)}\right)>0
                           \qquad(1\le k<N),
\]
because \(S_p(k)\le(p-1)(s-v_p(k))\).
The slope denominator is \(N\) if \(p\) is odd, since
\(2(N-1)/(p-1)\) is a \(p\)-unit. All root valuations have that
denominator, so every factor degree is a multiple of \(N\)
and \(J\) is irreducible over \(\mathbb Q_p\).

For \(p=2\) the denominator is \(N/2\). A repeated rational
factor would therefore consume the whole degree when squared.
Since \(c_0=1\), it would make the leading coefficient a
rational square. Write \(e=e_1\). The leading coefficient is
\((2N-1)!(eN)!/(N-1)!\). If \(e=1\), this is
\(N(2N-1)!\), with an odd valuation at a prime in \((N,2N)\).
If \(e\ge4\), a prime in \((eN/2,eN)\) occurs once in
\((eN)!\) and is greater than \(2N-1\). If \(e=2\), the
leading coefficient is \(2N((2N-1)!)^2/(N-1)!\), with an odd
valuation at a prime in \((N/2,N)\) when \(N\ge4\).
The required primes exist by Bertrand's postulate.
For \(e=N=2\), \(J(z)=1+8z+144z^2\), which is squarefree.

## Any number of private rows

Assume \(M\ge3\) and \(p\nmid M\). After an affine change of
ordinates, the lower hull is that of the points \((k,g(k))\),
where
\[
 g(0)=0,\qquad
 g(k)=(p-1)(v_p(c_k)-Lk)-Mk
     =2(p-1)(s-v_p(k))-M S_p(k)\quad(1\le k\le N).
 \tag{3}
\]
Let \(x_a=p^a-1\), \(1\le a\le s\). Then
\[
       g(x_a)=(p-1)(2s-Ma),\qquad g(N)=-M.
\]
The other points lie on or above the chain from the origin
through \(x_1,\ldots,x_s,N\), with the contacts described next.
It is enough to take the lower convex hull of that chain.

Between \(x_a\) and \(x_{a+1}\), write \(k=l p^a+t\) with
\(1\le l\le p-1\) and \(0\le t<p^a\). Let \(u=v_p(t)\) when
\(t>0\), and \(u=a\) when \(t=0\). The distance above the chord is
\[
 M\left((p-1)a-S_p(t)-1+\frac{t+1}{p^a}\right)
                              -2(p-1)u.                      \tag{4}
\]
If \(S_p(t)\le(p-1)(a-u)-1\), this is at least
\((M-2)(p-1)u+M(t+1)/p^a>0\). Otherwise
\(t=p^a-p^u\), and (4) equals
\[
             (M-2)(p-1)u-M\,\frac{p^u-1}{p^a}.
\]
It vanishes for \(u=0\). For \(u>0\) it is positive if
\((M-2)(p-1)\ge M\), which covers \(p\ge5\) and also
\(p=3\), since \(p\nmid M\). At \(p=2,\ M\ge5\), use
\((2^u-1)/2^a\le1/2\) when \(u=1\), and
\(2(M-2)>M\) when \(u\ge2\). So the only contacts of the layer
edge are \(k=l p^a-1\), \(l=1,\ldots,p\), including its endpoints.

The points \(0<k<x_1\) are strictly above the origin-to-\(x_1\)
chord by \(2s(p-1-k)\). The layer slopes in the transformed
ordinates are \(-M/p^a\), strictly increasing in \(a\).
The last slope is \((p-1)(M-2)s-M\), larger than the preceding
one on the relevant chain. For \(p=2,\ s=1,\ M\ge5\), this
comparison is between the last slope \(-2\) and initial
slope \(2-M\). The hull therefore starts with a chord from the
origin to a layer vertex \(x_a\), then keeps all remaining
layer edges and the last unit edge. No two layer vertices
can tie for the first chord. If they did, strict increase of
the intervening layer slopes would force adjacent vertices to
tie, which would give
\((p-1)(2s-Ma)/(p^a-1)=-M/p^a\), contrary to \(p\nmid M\).
Thus the first edge has just two supported indices, and its
width \(p^a-1\) is prime to \(p\). Its lattice residual is a
separable binomial. The last edge also has a separable residual.

For the remaining case \(p=2,\ M=3\), take \(s\ge3\) and start
the chain of layer vertices at \(x_3=7\). Formula (4) has the
same strictness for every \(a\ge3\). In the last displayed
expression, \(u=1,2\) use
\(3(2^u-1)/2^a\le3(2^u-1)/8<u\), and \(u\ge3\) uses
\(3(2^u-1)/2^a<3\le u\).
For \(1\le k\le6\) the ordinates \(g(k)\) are nonnegative
when \(s\ge4\), above the negative first chord. At \(s=3\)
they lie strictly above \(-3k/7\), the first chord.
The last slope \(s-3\ge0\) is again larger than the preceding
one. For \(s=1,2\), the ordinates \(v_p(c_k)-Lk\) are
\((0,2,3)\) and \((0,4,5,7,9)\), respectively. Both give
one primitive edge with only endpoint contacts.

It remains to find the residual on a layer edge of the hull.
The supported step is \(p^a\), which is also the slope
denominator because the valuation increment between contacts
is \(Lp^a+M(p^a-1)/(p-1)\) and \(p\nmid M\).
For \(k=l p^a-1\), the digit factorial products of
\(N!, (N+k-1)!, (N-k)!, (N-1)!, k!\) are, respectively,
proportional to
\[
                  1,\quad(l-1)!,\quad(p-l)!,\quad1,\quad(l-1)!,
\]
with signs independent of \(l\). The private powers \(e_i\)
leave these digit products for \(k!\) unchanged.
The factorial unit formula
\(p^{-v_p(t!)}t!\equiv(-1)^{v_p(t!)}\prod_j t_j!\pmod p\)
then gives the residual coefficients, up to an overall scalar
and a sign linear in \(l\), as
\[
                  \frac{((l-1)!)^{M-1}}{(p-l)!}
                         =(-1)^l((l-1)!)^M \quad\bmod p.
\]
The equality is Wilson's theorem. The valuation signs are
also linear in \(l\) along the edge, so they amount to a
unit dilation of the residual variable. Hence the lattice
residual is \(\mathcal R_{p,M}\) up to a nonzero scalar and
dilation. If \(M\equiv1\pmod{p-1}\), this polynomial equals
\(R(Y)=\sum_{j=0}^{p-1}j!\,Y^j\) modulo \(p\). It satisfies
\(R-Y(Y\,d/dY+1)R=1\). Thus it and its derivative cannot have
a common root, verifying the unconditional subcase.

## Shared gap two at power-of-two terminals

Take one private row of gap \(e=2^t\), and let the shared row
be \((1,3)\). The bridge changes \((N)_k\) in (1) to
\((N)_{2k}\). Write \(A=N/2\). After removing the linear part
\(ek\) from the valuations, the residual ordinates are
\[
             V(k)=s-v_2(k)+2v_2(k!)
                              +v_2\binom{A+k-1}{k},\qquad k\ge1 .
\]
Indeed the even and odd terms of \((N)_{2k}\) give \(k\)
plus the valuation of \((A)_k\). The last binomial valuation
equals \(s-1-v_2(k)\) for \(k<A\), zero at \(k=A,N\), and
\(s-v_2(k)\) for \(A<k<N\). This follows by counting the carries
when adding \(A-1\) to \(k\) in binary. The chord from zero to
\(N\) has slope \(2(N-1)/N\). Below \(A\), the distance above
the chord is
\(2(s-v_2(k)-S_2(k))-1+2k/N>0\), since
\(S_2(k)\le s-1-v_2(k)\).
Above \(A\) the distance is
\(2(s-v_2(k)-S_2(k))+2k/N>0\), since
\(S_2(k)\le s-v_2(k)\).
The contacts are exactly \(0,A,N\). In the lattice variable
the residual is \(1+Y+Y^2\) modulo two, which is squarefree.

The proofs give characteristic-zero coprimality of the original
moments by the exact Euler bridge and rational Bézout identity.
