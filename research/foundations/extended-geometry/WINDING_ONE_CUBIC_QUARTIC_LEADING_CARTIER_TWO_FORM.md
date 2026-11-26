# A Cartier obstruction on the double cover of a quartic-leading cubic

Let \(K\) have characteristic zero and put \(U=1+x\). Use the
regular winding-one presentation
\[
 f=\frac{U}{x}P(x,Ut),\qquad P(x,v)=v^3+d(x)v^2+b(x)v+c(x)
\]
on \(2xy+t^2=1\). Assume \(\deg d\le3\), \(\deg b\le5\), and
\(\deg c\le7\). These are the bounds for least sphere degree at most
seven in the normalized \(E=U^4\) case of the
[quartic leading-root reduction](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root).
The functional \(\Lambda\) takes the angular constant term in \(x\)
of the normalized integral in \(t\) over \([-1,1]\).

**Theorem (Double-cover Cartier obstruction).** Suppose the positive
moments \(\Lambda(f^m)\) vanish eventually. Over a number field, at all
but finitely many primes the reduction of the rational two-form
\[
 \omega=\frac{dx\wedge dv}{U Y},\qquad Y^2=xU P(x,v),
\]
has zero Cartier operator and is exact in the function field of the
double cover. In the original coordinate \(t\) this is
\(dx\wedge dt/(xZ)\), where \(Z=Y/x\) and \(Z^2=f\).
Over an arbitrary characteristic-zero field, the same obstruction
holds for every number-field specialization of the coefficient algebra.

## A marker and the selected coefficients

The fixed marker \(xU\) has infinitely many nonzero moments for any
phase in the displayed family. First work over a number field. For a
good rational prime \(p'=3k+1\), the highest height of \(xU f^k\) is
\(3k=p'-1\). Only the interval denominator at this even height is
divisible by \(p'\). Thus its scaled moment satisfies
\[
 p'\Lambda(xU f^k)\equiv [x^{k-1}]U^{4k+1}
 = [x^{k-1}]U^{p'+k}\equiv k\pmod{\mathfrak p'}.
\]
In the last step \(U^{p'}\) contributes only its constant term,
since \(k-1<p'\). This proves the claim by primes in the class one
modulo three. It also proves it over any characteristic-zero field
after one number-field specialization, since all individual moments
are rational polynomials in the coefficients.

Now use an arbitrary sufficiently large odd rational prime \(p\),
put \(m=(p-1)/2\), and reduce the coefficients in a residue field at
\(p\). Define
\[
 H_p(x)=[v^{p-1}]P(x,v)^m,\qquad
 A_0=[x^m]U^{m-1}H_p(x),\qquad
 A_1=[x^{m+p}]U^{m-1}H_p(x).
\]
A product contributing to \(H_p\) uses the four respective powers
\(v^3,v^2,v,1\) with multiplicities \(i_3,i_2,i_1,i_0\) satisfying
\(i_3=i_1+2i_0\). Its \(x\)-degree is at most
\(3i_2+5i_1+7i_0=3m-i_1-2i_0\). Hence \(\deg H_p\le3m\). The
coefficient of height \(p-1\) in \(f^m\) is
\(x^{-m}U^{m+p-1}H_p(x)\). Its angular exponents that are multiples
of \(p\), written in a variable \(X\) for \(x^p\), give
\[
                   (1+X)(A_0+A_1 X).                          \tag{1}
\]
Indeed \(U^{m+p-1}=(1+x^p)U^{m-1}\) in the residue field, and all
powers in \(x^{-m}U^{m-1}H_p\) lie between \(-m\) and \(3m-1\).
Only zero and \(p\) in that range are divisible by \(p\).

Fix \(l\ge0\). In the residue field \(f^{lp+m}=(f^l)^p f^m\).
An interval denominator divisible by \(p\) requires the height
coming from the last factor to be congruent to \(-1\) modulo \(p\).
The last factor has maximum height \(3m<2p-1\), so that height must
be \(p-1\). If the height of a term of \(f^l\) is \(j\), the total
height is \(p-1+pj\), even exactly when \(j\) is even. Multiplying
the interval moment by \(p\) then leaves denominator \(j+1\), a unit
for \(p>3l+1\). Angular matching selects just (1), so
\[
 p\Lambda(f^{lp+m})
 \equiv A_0\Lambda(U f^l)^p+A_1\Lambda(xU f^l)^p
                                      \pmod{\mathfrak p}.       \tag{2}
\]
The powers on the right are residue-field Frobenius, even when the
residue degree is greater than one. Setting \(l=0\) in (2) and using
eventual purity gives \(A_0=0\). Then use a fixed \(l\) with
\(\Lambda(xU f^l)\ne0\), whose existence was proved above. This
value is a unit at all but finitely many primes. Equation (2) gives
\(A_1=0\) as well, at all sufficiently large good primes.

## Cartier on the double cover

In characteristic \(p\), the elements \(x,v\) are a separating
transcendence basis of the function field of the double cover. In
fact the cover is separable and remains a field after reduction,
since \(xU P(x,v)\) has odd degree in \(v\) and cannot be a rational
square. Away from characteristic two we can write
\[
 \omega=(Y^{-1})^p x^m U^{m-1}P(x,v)^m\,dx\wedge dv.
\]
The Cartier formula for a top form picks, after factoring out
\(p\)th powers, the terms with both exponents congruent to \(p-1\)
modulo \(p\), then takes the \(p\)th root of the coefficient and
divides each shifted exponent by \(p\). One can see this directly
by writing the coefficient in the \(p\)-basis \(x^i v^j\),
\(0\le i,j<p\), over the subfield of \(p\)th powers. Every resulting
term with \(i<p-1\) or \(j<p-1\) can be integrated rationally in
that variable, since the other coefficient is a \(p\)th power and
the relevant integer \(i+1\) or \(j+1\) is invertible. The remaining
coefficient is the Cartier coefficient. This is the top-degree
case of the [Cartier isomorphism, stated in Katz, Theorem 7.2](https://web.math.princeton.edu/~nmk/old/nilpconn.pdf).

Here the only selected power of \(v\) is \(p-1\), since \(3m<2p-1\).
On that row the \(x\)-degree is at most \(5m-1<3p-1\). The selected
coefficients at \(x^{p-1}\) and \(x^{2p-1}\) are precisely \(A_0,A_1\)
after the factor \(x^m\) is removed. Both are zero. This proves the
asserted vanishing and rational exactness.

For specialization, take the finitely generated coefficient algebra
over \(\mathbb Q\) containing the coefficients of \(P\). Each moment
is a rational polynomial in those coefficients, and the shape and
degree bounds are preserved by every homomorphism to a number field.
The identities of a vanishing tail therefore hold in the specialized
field. The argument above applies to that field, including to primes
of arbitrary residue degree. The conclusion concerns exactness after
reduction. In particular it supplies no characteristic-zero primitive
of \(\omega\) and imposes no ordinary-reduction assertion on the
minimal surface of the double cover.
