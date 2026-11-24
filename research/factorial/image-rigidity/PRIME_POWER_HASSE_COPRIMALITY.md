# A Hasse-derivative coprimality criterion at prime-power gaps

Let \(q=p^s\), where \(p\) is prime and \(s\ge1\). Take \(n\ge1\), a positive
rational \(\gamma\) integral at \(p\), and an integer \(1\le h\le q\).
The actual pair is
\[
 T_{n,\gamma,q}(t)=\sum_{i=0}^{n}\binom ni(\gamma)_{qi}t^i,
 \qquad T_{n+1,\gamma+h,q}(t).
\]
Write \(\sigma_h(x)=[Y^2](x+Y)_h\), with the rising factorial convention.

## Constant term of the upper functional

The pair is coprime if \(v_p(\sigma_h(\gamma))<s\).
In particular a nonzero reduction of \(\sigma_h(\gamma)\) modulo \(p\) suffices.
More generally, in the notation of the
[upper-response theorem](SFC_PRIME_POWER_UPPER_RESPONSE.md), it is
coprime whenever \(T_0\ne0\) and
\[
                           v_p(T_0)<s+D,
\]
with \(D=+\infty\) allowed. The conclusions extend to every characteristic-zero field.

Use \(b=(q-1)/(p-1)\), \(B=b-s\), \(C=\gamma+q(n+1)\) and
\(\kappa=q(n+1)/(C)_h\). The closed-disk response necessary for a common
root is equation (3) in the upper-response note. Substituting its
equation (2) gives
\[
 F(v)=A_0-(-1)^q\kappa p^b T_0 v
        +p^{B+1}v A(v)
        -(-1)^q\kappa p^{b+s+D}v^2 R(v),\qquad v_p(v)\ge0 .
\]
Here \(v_p(A_0)=B\) and \(A,R\) are integral restricted series.
If \(D=+\infty\) the last term is absent by the exact functional formula.
Put \(\ell=v_p(\kappa)+b+v_p(T_0)\). The last two terms have coefficients
of valuation strictly greater than \(B\) and \(\ell\), respectively.
After dividing by \(p^{\min(B,\ell)}\), reduction therefore gives a
nonzero polynomial of degree at most one. When it is constant there is
no integral root. When it is linear the divided difference at any
two distinct integral arguments is a unit, so there is at most one
integral root. This works on the entire closed disk over a complete
algebraically closed extension.

Every common root of the original pair contributes one such \(v\).
The transformation \(v=(-1)^{q+1}/(p^b t)\) is injective and
\(T_{n,\gamma,q}(0)=1\). A unique distinct common root of two rational
polynomials is rational by conjugacy in characteristic zero. The
positive rational Gamma Gram factorization recalled in the same
upper-response theorem excludes a rational common root, because one of
the two powers is even. Rational Bézout identities give the field claim.

To compute \(T_0\), use the
[exact divided difference](SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md#1-the-functional-and-its-divided-differences).
For \(f_h(Y)=(Y)_h\), \(r=(\gamma+q-1)\bmod q\) and \(A_r=r-q+1\), it is
\[
                   T_0=(-1)^q f_h[C,C-q,A_r].
\]
All three nodes are congruent to \(\gamma\) modulo \(p^s\).
The divided difference is an integral polynomial in its nodes, including
where they coalesce. Thus
\(T_0\equiv(-1)^q\sigma_h(\gamma)\pmod{p^s}\).
If \(v_p(\sigma_h(\gamma))<s\), this is also the valuation of \(T_0\).
Every finite \(D\) is nonnegative by integrality of the actual upper
functional, so the strict inequality above applies.

## Lower-four gap nine and the offset twice the prime

At \(p=2\) or \(p=3\), the polynomial
\[
                    \sigma_4(x)=6x^2+18x+11
\]
is a nonzero constant modulo \(p\). Thus every allowed gap \(q=p^s\ge4\)
with either prime admits the offset \(h=4\) for every indicated
\(\gamma\). In particular, at \(q=9\), set \(\gamma=4n+1\).
The actual ordinary consecutive pair is coprime at every order.
The usual moment normalization and the monomial cases give all-order
Strong Factorial for the support \(\{4,13\}\) over every
characteristic-zero coefficient field.

At \(p=5\), the identity \(24\sigma_4(x)=(12x+18)^2-60\) gives
\(v_5(\sigma_4(\gamma))\in\{0,1\}\) on the entire integral disk.
The square has valuation zero or at least two, whereas \(60\) has
valuation one. Hence all allowed gaps \(5^s\) with \(s\ge2\) also
admit offset four for every positive rational \(5\)-integral shape.

## Offset twice a prime power

Take \(h=2p^w\le q\) with \(w\ge1\), so \(s\ge w+1\).
Then the mixed pair is coprime for every positive rational lower
shape integral at \(p\).

In \(f_h[C,C-q,A_r]\), expand the divided-difference product rule for
the \(h\) linear factors of \(f_h\). Each term omits two factors, with
the other factors evaluated in successive blocks at the three nodes.
All three nodes are congruent to \(\gamma\) modulo \(p^s\).
Among the \(h\) consecutive translates \(\gamma+j\), exactly two
are divisible by \(p^w\). The term omitting both of these factors has
valuation exactly
\[
                  v_p(T_0)=\sum_{\ell=1}^{w-1}(2p^{w-\ell}-2).
\]
Indeed at each \(\ell<w\) there were \(2p^{w-\ell}\) divisible translates
before the omissions, and none of the remaining factors is divisible
by \(p^w\). Every other term has valuation at least one greater.
It removes at most two divisible translates at each depth below \(w\)
and leaves at least one factor divisible by \(p^w\). The minimal term
is unique, so cancellation cannot raise its valuation.

The [third divided-difference bound](SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md)
gives \(D\ge\sum_{\ell=1}^{w-1}(2p^{w-\ell}-3)\) at this \(h\).
Thus \(v_p(T_0)\le D+w-1<s+D\) when \(D\) is finite. The same nonzero
\(T_0\) also qualifies when \(D=+\infty\).
The coprimality criterion above now applies for every order.
For \(w=1\) this recovers the unit \(T_0\) case at offset \(2p\).
For the ordinary specialization \(a=h\), the theorem includes
\[
 \{2\cdot 3^w,5\cdot 3^w\},\qquad
 \{2\cdot 5^w,7\cdot 5^w\}\qquad (w\ge1).
\]
The parameter \(\gamma=an+1\) is integral at the corresponding prime
and the moment normalization gives the claim over every
characteristic-zero coefficient field.
