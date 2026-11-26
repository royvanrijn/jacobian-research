# Directed factorial prefixes constrain every common factor degree

Consider two adjacent ordinary factorial moment polynomials of
distinct exponent columns. Write \(m\) for either marked order and
\(m'\) for the other. Choose one original unequal row and reverse
both columns if needed to make its entries \((a,a+d)\), \(d>0\).
Put \(Q=\max(m,m')\) and \(H=mm'd\). Strip the scalar factorial
of each row's smaller exponent at each order. All remaining
coefficient factors are positive integers.

Fix integers \(1\le A\le K\le\min(m,m')\) and \(0\le k_0<A\).
For \(x>0\) integral, write \(x_{>H}\) for its part supported on
primes greater than \(H\). Use \(C\) for the product of the
retained row factors sampled at index \(k_0\) of order \(m\),
omitting the binomial coefficient. Use \(B^+\) for the product
of the full stripped ramps at order \(m'\) of all retained
increasing rows. The chosen row is omitted from both products.

**Theorem.** A rational common factor of degree at least \(K\)
forces
\[
 \prod_{\substack{(A-1)d<h\le Kd\\m\nmid h}}(ma+h)_{>H}
               \ \mid\ \bigl(C(B^+)^{\,m-k_0}\bigr)_{>H}.       \tag{1}
\]
This condition holds for every choice of original row and marked
order. There is no assumption on upper exponent alignments, gap
ratios, dimension, or the size of the retained exponents.
Every common factor has even degree because the even factorial
moment is strictly positive on the real line. In particular one may
take \(K=2\) when testing whether a physical common root can exist.

## A bound from one anchor coefficient

Let \(p>H\) divide \(x=ma+h\) in the indicated prefix and put
\(E=v_p(x)\), \(r=\lceil h/d\rceil\). The chosen row samples
\(x\) at all indices \(k\ge r\) at the marked order.
Its other marked-order factors are units, since their nonzero
offsets from \(x\) have magnitude less than \(md\).
At the other order,
\[
                       m(m'a+t)=m'x+(mt-m'h)
                         \qquad(1\le t\le m'd).
\]
The last term is nonzero because \(m\nmid h\) and
\(\gcd(m,m')=1\). Its magnitude is at most \(H\).
Thus the chosen row contributes only units at the other order.
At both orders the binomial coefficients are units.

The marked ordinate at \(k_0<r\) is exactly \(U=v_p(C)\),
whereas every ordinate at \(k\ge r\) is at least \(E\).
The other polygon has nonnegative ordinates and leading
ordinate \(V=v_p(B^+)\); its maximum slope is at most \(V\)
by the last-edge slope formula. For any slope \(s\) there and
any \(k\ge r\), comparing the two marked ordinates after
subtracting \(s\) times their indices gives a difference
at least \(E-U-(m-k_0)V\). If this is positive, every marked
edge with a shared slope has both endpoints among indices
\(0,\ldots,r-1\). The [ordinary Newton polygon root count](https://arxiv.org/pdf/0807.2620)
then bounds common factor degree by \(r-1<K\).
Consequently \(E\le U+(m-k_0)V\) is necessary.
Since the entire chosen prefix spans less than \(H\), each
prime above \(H\) divides at most one of its terms. This proves (1).

## A quantitative necessary retained weight

Assume total absolute curvature below \(12\), so the chosen
row's upper exponent \(M=a+d\) satisfies \(M>d^2/12\).
Suppose \(d\ge24\).
Let \(w=(K-A+1)d\), and let \(k\) count the integers \(h\)
in (1). Define
\[
 \mathcal L_\pm=
   \sum_{\substack{i\ {\rm retained\ of\ orientation}\ \pm\\
                   Q M_i>H}}d_i\log(Q M_i),
       \qquad M_i=\max(a_i,b_i).
\]
Then any common factor of degree at least \(K\) forces
\[
 \begin{split}
 (k-\pi(H))\log M-\log((w-1)!)\ \le\
 &(m-k_0)\mathcal L_-\\
 &+\bigl(k_0+(m-k_0)m'\bigr)\mathcal L_+ .                      \tag{2}
 \end{split}
\]
Indeed every mark \(ma+h\) exceeds \(M\).
After removing one maximal-valuation term per prime at
most \(H\), the small-prime part of the remaining product
divides \((w-1)!\), since indices divisible by a given prime
power are separated by at least that prime power. The
large-prime part satisfies (1). In \(C\) each decreasing
retained row contributes \((m-k_0)d_i\) factors and each
increasing one \(k_0d_i\). The increasing rows contribute
\(m'd_i\) factors in \(B^+\). All are at most \(Q M_i\),
and those from the omitted rows are units at primes above
\(H\). This gives (2).

For example, \(A=1,k_0=0\) uses the entire prefix through
the \(K\)-th block. Since \(k\ge(1-1/m)Kd-1\),
the left side of (2) is at least
\[
 \bigl((1/2-1/m)Kd-1-\pi(mm'd)\bigr)\log M
                             -Kd\log(K\sqrt{12}).
\]
With the adjacent orders fixed and \(d\) tending to infinity,
the leading coefficient is \(1/2-1/m-o(1)\) times
\(Kd\log M\). Here one may use the elementary prime count
in the [retained-prime proof](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#uniform-estimates-for-the-certificate).
The dependence on the candidate common factor degree is linear.
Choosing \(k_0=A-1\) instead allows a smaller opposing-row charge
when the marked blocks are near the end of the ramp.

## Opposing large rows and the block positions

For a further necessary condition, suppose the only other row
that can contribute a nonunit at \(p\) is one decreasing row
of gap \(e\) and smaller exponent \(b\). Take
\(p>mm'\max(d,e)\) and the same marked factor as in (1).
If a common factor has degree at least \(K\ge r\), then an
opposing factor \(mb+t\) must be divisible by \(p\).
It is unique since the opposing ramp length is less than \(p\).
Put \(F=v_p(mb+t)\) and \(q=\lceil t/e\rceil\).
Then
\[
                    F\ge E,\qquad r+q\le m+1.                \tag{3}
\]
Indeed all opposite signs in the other degree give
nonpositive slopes: its coefficient valuations consist
of a nonnegative constant on an initial segment and zero
thereafter, or all zero. If \(F<E\), the marked constant
ordinate is strictly smaller than every ordinate at \(k\ge r\).
This confines shared slopes to indices below \(r\), as in
the anchor argument. If \(r+q\ge m+2\), the index \(r-1\)
has neither mark and again every index at least \(r\) lies
strictly above it for the relevant slopes. This proves (3).
When \(q\le K\) and \(m\nmid t\), reversing the two roles
gives \(E\ge F\), so the depths must be equal.

These are restrictions on every original counterexample
and its possible factor degrees. They permit mass in several
different endpoint bands and leave open whether the resulting
simultaneous directed restrictions are sufficient to exclude
every profile.
