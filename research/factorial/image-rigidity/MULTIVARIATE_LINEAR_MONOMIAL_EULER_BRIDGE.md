# An exact Euler bridge when one monomial is linear

Omit equal-exponent rows. When one of the two ordinary monomials
has degree one and its variable is shared, write the phase as
\(X(x+X^dV)\), with \(d\ge0\) and \(V\) a monomial in the remaining
variables. Under the factorial integral \(X\) is unit exponential
and independent of \(V\). Consecutive coprimality is exactly a
squarefreeness condition for one intermediate polynomial. For an
ordinary private monomial \(V\), this gives a uniform finite window
over all starting orders. The exact bridge also identifies the
limitation of allowing arbitrary positive laws for the private
factor.

## The intermediate polynomial and windows

Let \(G_\gamma\) denote a unit-rate Gamma random variable of shape
\(\gamma\), independent of \(V\). For terminal \(N=n+1\), put
\[
 B_N(x)=\mathbb E(x+G_N^dV)^N,\qquad E=x\,\frac{d}{dx},
 \qquad P_m(x)=\mathbb E\bigl[(X(x+X^dV))^m\bigr].
\]
Then
\[
 \frac{P_n(x)}{n!}=\frac{B_N'(x)}N,\qquad
 \frac{P_N(x)}{N!}=(1+d)B_N(x)-\frac{d}{N}E B_N(x).
 \tag{1}
\]
Indeed the factor \(X^m\) tilts the exponential to \(G_{m+1}\), and
a Gamma shape advance multiplies its \(dk\)-th moment by
\((N+dk)/N\). Therefore the consecutive gcd is
\(\gcd(B_N,B_N')\), up to scalar. The
[univariate intermediate polynomial](ORDINARY_BINOMIAL_INTERMEDIATE_SHAPE_MULTIPLICITY.md)
is the specialization \(V=1\). If the degree-one variable is
private instead, divide each order's polynomial for \(X+zV\) by
\(m!\). The result is \(\sum_{k=0}^m\mathbb E V^k\,z^k/k!\).
Successive polynomials differ by a nonzero monomial, so they
are coprime.

More generally take starting order \(n\ge1\) and terminal \(M=n+L\).
Freeze the shape at \(n+1\) in \(B(x)=\mathbb E(x+G_{n+1}^dV)^M\).
For \(0\le h\le L\), put \(m=n+h\). Up to nonzero scalars the
physical moment at degree \(m\) is
\[
    \prod_{j=0}^{h-1}\bigl(n+1+j+d(m-E)\bigr)
                     \frac{d^{\,L-h}B}{dx^{L-h}} .
 \tag{2}
\]
This follows from the fixed-shape derivative and the \(h\) unit
shape advances. At any fixed \(x\), the expression (2) is a
linear combination of \(B^{(L-h)}(x),\ldots,B^{(L)}(x)\).
The coefficient of \(B^{(L-h)}(x)\) is
\(\prod_{j=0}^{h-1}(n+1+j+dm)>0\).
Starting with \(h=0\) and going to \(L\) therefore proves that
the entire window vanishes at \(x\) precisely when \(x\) is a
root of \(B\) of multiplicity at least \(L+1\).

For \(V=\prod_i Y_i^{e_i}\) with independent unit exponentials and
integer \(e_i\ge1\), let \(D=d+\sum_i e_i\), and suppose \(D\ge1\).
Write \(J(z)=z^M B(1/z)\) and \(\mathcal E=z\,d/dz\). Comparing
consecutive coefficients gives
\[
 \mathcal E J=z(M-\mathcal E)
 \prod_{j=0}^{d-1}(n+1+j+d\mathcal E)
 \prod_i\prod_{j=1}^{e_i}(j+e_i\mathcal E)\,J .
 \tag{3}
\]
The empty product at \(d=0\) equals one. The equation is of order
\(D+1\) with nonzero leading derivative coefficient away from
\(z=0\). Since \(J(0)=1\) and \(B(0)>0\), uniqueness for the
linear differential equation implies that every root of \(B\)
has multiplicity at most \(D\).

Thus for this original factorial family, no window of \(D+1\)
consecutive moments can vanish at any nonzero projective coefficient
direction. At the first or last coordinate axis the phase is a
single monomial and all moments are nonzero. For a pair at terminal
\(N\), the degree of its gcd is at most
\(2\lfloor (N-\lceil N/D\rceil)/2\rfloor\). Indeed \(B_N\) has at
least \(\lceil N/D\rceil\) distinct roots, and the gcd in (1) has
no real roots by strict positivity of the even moment on the real
axis. All coefficients and the identities are rational for the
ordinary family, so the nonvanishing window extends to every
characteristic-zero field.

## A collision when the private factor has an arbitrary positive law

Keep the exact exponential variable \(X\) and take \(d=1\).
Let the independent positive variable \(V\) have the following
values and probabilities. Its first five moments will be used.

| Value | Probability |
| ---: | ---: |
| \(1/14\) | \(12852553/80626500\) |
| \(1/6\) | \(2643/483140\) |
| \(2/3\) | \(10617/22750\) |
| \(3/4\) | \(262336/712215\) |
| \(3\) | \(79/878220\) |
| \(5\) | \(2137/185732820\) |

The moments of degrees zero through five are
\[
               1,\quad 3/5,\quad 5/12,\quad 25/84,\quad
                                 25/112,\quad 625/3024.
\]
Multiplying the \(k\)-th moment by \((5)_k\) gives the moments of
\(G_5V\). Consequently
\[
                   B_5(x)=(x+5)(x^2+5x+25)^2 .
\]
Equation (1) gives the exact common factor \(x^2+5x+25\) of
\(P_4,P_5\). The centered repeated quadratic comes from the
[known Appell collision](GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex),
with a translation by three. The new factorization of the moments
keeps the exact Gamma shape five and both linked operators in (1).
The law in the table is the synthetic replacement for a private
factorial monomial.

~~~python
import sympy as sp

q = sp.Rational
nodes = [q(1, 14), q(1, 6), q(2, 3), q(3, 4), 3, 5]
weights = [q(12852553, 80626500), q(2643, 483140), q(10617, 22750),
           q(262336, 712215), q(79, 878220), q(2137, 185732820)]
moments = [sum(weight*node**k for weight, node in zip(weights, nodes))
           for k in range(6)]
assert all(weight > 0 for weight in weights)
assert moments == [1, q(3, 5), q(5, 12), q(25, 84), q(25, 112), q(625, 3024)]
x = sp.Symbol("x")
B = sum(sp.binomial(5, k)*sp.rf(5, k)*moments[k]*x**(5-k)
        for k in range(6))
factor = x*x+5*x+25
assert sp.expand(B-(x+5)*factor**2) == 0
lower = sum(sp.binomial(4, k)*sp.rf(5, k)*moments[k]*x**(4-k)
            for k in range(5))
upper = sum(sp.binomial(5, k)*sp.rf(6, k)*moments[k]*x**(5-k)
            for k in range(6))
assert sp.expand(lower-sp.diff(B, x)/5) == 0
assert sp.expand(upper-2*B+x*sp.diff(B, x)/5) == 0
assert sp.Poly(sp.gcd(lower, upper), x).monic().as_expr() == factor
print("Positive private law and exact linked quadratic gcd passed")
~~~
