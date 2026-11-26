# Every-degree separation at signed imbalance three for private unit pairs

Consider an ordinary two-monomial profile with exactly one private
unit row on each side and arbitrary shared rows. Write
\(\Delta=\sum_i(b_i-a_i)\) for the signed degree difference.

**Theorem.** If \(|\Delta|\ge3\), the consecutive physical
factorial moments are coprime in every degree. The
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
gives the same conclusion for \(\Delta=0\).
If a common root occurs with \(|\Delta|=2\), orient the columns
so that \(\Delta=2\). Then some decreasing shared row has
absolute gap \(d\) and smaller exponent \(\ell\) with
\(d/\ell>1/2\). Every such row has \(d\le5\) and \(\ell<2d\),
and there are at most three, counting multiplicity.

The proof compares a centered imaginary projection of the
two moments. At a common first-box angle, the exact coefficient
quotient needs one comparison between the middle and the
last lower index. A signed loss bound over that range
supplies it.

## Reduction to the middle coefficient quotient

Write \(V=\sum d_i^2/\ell_i\) and
\(S=\sum d_i^2/(\ell_i+d_i)\) over unequal shared rows,
where \(d_i=|b_i-a_i|\). We may orient the columns with
\(\Delta>0\). A reversal reciprocates the roots, and
conjugation after reversal returns a nonreal upper root to
its original argument.

The [curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
settles \(S\ge15/8\). With \(S=0\) the moments are
geometric up to their common endpoint argument scale.
The [ordinary terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
handles the first thirteen terminals. Thus take
\(0<S<15/8\), \(n\ge13\), and \(N=n+1\).
By the [higher-label separation theorem](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md),
a common nonreal upper root has first label. Its argument
satisfies
\[
                        \pi/n<\theta<2\pi/(n+2).
 \tag{1}
\]
The even moment is strictly positive on the real line,
so real roots cannot be common.

Let \(A_k\) and \(B_k\) denote the degree-\(n\) and
degree-\(N\) coefficients after dividing by their constant
coefficients, with the physical argument unchanged.
Put \(v_k=\log(B_k/A_k)\) for \(0\le k\le n\) and
\(j=\lceil n/2\rceil\).
The [exact positive Mellin representation](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md#exact-factorial-link-and-the-radius-interval)
with scale one makes \(A_k/B_k\) log-convex in \(k\).
Hence \(v_k\) is concave.

Suppose that \(v_n\ge v_j\). Concavity then gives
\(v_k\le v_j\) for \(k<j\) and \(v_k\ge v_j\) for
\(j<k\le n\). At a root \(r e^{i\theta}\), rotate the
value of each polynomial by \(e^{-ij\theta}\).
By (1) the numbers \(\sin((k-j)\theta)\) are negative
for \(k<j\) and positive for \(j<k\le N\).
Indeed \(j\le(n+1)/2\) and \(N-j\le(n+2)/2\), so both
end angles have absolute value strictly less than \(\pi\).
The lower moment's imaginary sum vanishes. Termwise,
\[
     \sum_{k=0}^n B_k r^k\sin((k-j)\theta)
       \ \ge\ e^{v_j}\sum_{k=0}^n A_k r^k\sin((k-j)\theta)=0.
\]
The additional upper term \(B_Nr^N\sin((N-j)\theta)\)
is strictly positive, contradicting an upper root.
It suffices to prove \(v_n\ge v_j\).

## The signed loss over the upper half

Define \(\eta_*=7/13\) and for \(r>0\) put
\[
        L(r)=\frac1{1-\eta_*}
                  \int_{\eta_*}^1\frac{r y}{1+r(1-y)}\,dy .
\]
For each row write \(S_i=d_i^2/(\ell_i+d_i)\) and
\(r_i=d_i/\ell_i\). Set the oriented cost
\[
              C=\sum_{\mathrm{increasing}}S_i+
                      \sum_{\mathrm{decreasing}} d_i L(r_i).
 \tag{2}
\]
We claim
\[
                 v_n-v_j\ge(n-j)(\Delta-C)\log(N/n).
 \tag{3}
\]
After cancellation of the two private unit factorials, a
constant-normalized coefficient at real order \(m\), with
integer index \(k\), has one factor
\(\Gamma(m a_i+(b_i-a_i)k+1)/\Gamma(m a_i+1)\)
per shared row. For \(n\le m\le N\), differentiate the
log of the ratio of its index \(n\) and index \(j\)
values. On an increasing row this is
\(\sum_{h=d_i j+1}^{d_i n} a_i/(m a_i+h)\).
The difference from the ideal contribution
\(d_i(n-j)/m\) is a nonnegative loss
\(\sum_h h/[m(m a_i+h)]\). Since \(h\le m d_i\),
it is at most \((n-j)S_i/m\).

On a decreasing row the derivative is
\(-\sum_{h=d_i j}^{d_i n-1}a_i/(m a_i-h)\), and the
loss relative to \(-d_i(n-j)/m\) is
\(\sum_h h/[m(m a_i-h)]\).
Here \(a_i=\ell_i+d_i\). Substituting \(y=h/(n d_i)\)
and using \(m\ge n\), the loss per term is at most
\(m^{-1}r_i y/[1+r_i(1-y)]\).
The average over the indicated \(h\) is a left Riemann
sum for the increasing function in the definition of
\(L\), on \([j/n,1]\). Since \(j/n\le7/13\) for
\(n\ge13\), this average is at most \(L(r_i)\).
The decreasing-row loss is therefore at most
\((n-j)d_i L(r_i)/m\).
Summing the signed contributions and integrating from
\(n\) to \(N\) proves (3).

## The absolute integer threshold

For \(r\le1\), the convex integrand in \(L(r)\) is bounded
on its interval by the chord through its endpoint values.
The resulting row cost divided by \(S_i\) is at most
\[
           \frac{1+r}{2}
                \left(1+\frac{\eta_*}{1+(1-\eta_*)r}\right)
                   \le\frac{26}{19}.
 \tag{4}
\]
The expression is increasing in \(r\): in the second
term the ratio \((1+r)/(1+(1-\eta_*)r)\) is increasing.
The increasing-row cost in (2) is just \(S_i\).

Under \(S<15/8\), at most one row has \(r_i>1\).
The only possible \((\ell_i,d_i)\) are \((1,2),(2,3)\):
if \(d_i\ge4\), its \(S_i=d_i r_i/(1+r_i)\) exceeds two;
for \(d_i\le3\) one checks the integral \(\ell_i<d_i\).
Both allowed rows have \(S_i\ge4/3\).
If the exceptional row is increasing, or if there is no
exceptional row, (4) gives \(C\le(26/19)S<3\).
For a decreasing \((1,2)\) row,
\[
 d_iL(r_i)=\frac{13}{2}\log(25/13)-2
                 <\frac{13}{2}\frac{131}{200}-2 .
\]
Add at most \((26/19)(13/24)\) for the other rows,
since the exceptional row costs \(4/3\) in \(S\).
The total remains less than three.
For a decreasing \((2,3)\) row,
\[
 d_iL(r_i)=\frac{65}{6}\log(22/13)-3
                 <\frac{65}{6}\frac{53}{100}-3 .
\]
The other rows cost at most \((26/19)(3/40)\) after
the exceptional \(S_i=9/5\), again giving less than three.
Thus \(\Delta\ge3\) implies \(v_n>v_j\) in (3).
The centered imaginary projection excludes the first label,
which proves the all-degree claim.

For the stated necessary condition at \(\Delta=2\),
assume all decreasing rows have \(r_i\le1/2\). Their
cost-to-\(S_i\) ratio is \((1+r_i)L(r_i)/r_i\).
It is increasing in \(r_i\) because the integrand
\((1+r)y/[1+r(1-y)]\) is increasing in \(r\) for every
\(0<y\le1\). At \(r=1/2\) this ratio is
\(3[(13/2)\log(16/13)-1]<132/125\), using
\(\log(16/13)<26/125\). Thus \(C<(132/125)(15/8)<2\)
and the first-box argument applies.
Finally a row with \(r_i>1/2\) has \(S_i>d_i/3\), so
\(d_i\le5\). It costs at least \(1/2\) in \(S\):
for \(d_i=1\) the only such smaller exponent is one,
and for \(d_i\ge2\) the cost exceeds \(2/3\).
There are at most three such rows.

The conclusions are in the original physical variable.
Ordinary coefficients are rational, so coprimality
gives a rational Bézout identity and transfers to every
characteristic-zero field.

## Exact scalar replay

```python
from fractions import Fraction as Q

def log_upper(a, b):
    x = Q(a-b, a+b)
    assert 0 < x < 1
    return 2*(x+x**3/3+x**5/5+x**7/(7*(1-x*x)))

assert log_upper(25, 13) < Q(131, 200)
assert log_upper(22, 13) < Q(53, 100)
assert log_upper(16, 13) < Q(26, 125)
eta = Q(7, 13)
assert 1+eta/(2-eta) == Q(26, 19)
assert Q(26, 19)*Q(15, 8) < 3
assert Q(13, 2)*Q(131, 200)-2+Q(26, 19)*Q(13, 24) < 3
assert Q(65, 6)*Q(53, 100)-3+Q(26, 19)*Q(3, 40) < 3
assert 3*(Q(13, 2)*Q(26, 125)-1) == Q(132, 125)
assert Q(132, 125)*Q(15, 8) < 2
```
