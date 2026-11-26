# Linked coefficient order for factorial binomials with private variables

Let \(F_m(z)=\sum_{k=0}^m c_{m,k}z^k\) be an ordinary factorial binomial
with a private variable in each monomial. Normalize its endpoint
coefficients to one by positive scaling of the polynomial and the
variable. The private exponents may be arbitrary positive integers.
An optional common Gaussian factor is allowed, with normalized
coefficient \(e^{-\tau k(m-k)/(2m)}\), \(\tau\ge0\).

Let \(N=n+1\). In either orientation of the exponent columns,
\[
 c_{n,k}^{\, (N/n)(N-k)/(n-k)}
               \ \le\ c_{N,k}\ \le\ c_{n,k}
                        \qquad(0\le k<n).                      \tag{1}
\]
The upper bound also holds at \(k=n\). In particular the upper half
coefficients at distance \(k\) from their respective top indices
satisfy the same order. In either orientation there is a positive
random variable \(W\) such that
\[
                  \frac{c_{n,k}}{c_{N,k}}=\mathbb E W^k
                             \qquad(0\le k\le n).              \tag{2}
\]
Here \(W\) is an explicit product of powers of independent Beta
variables, with a deterministic scale and an optional lognormal
factor. The quotient in (2) is nondecreasing and log-convex in the
index \(k\).

## Canceling the binomial factor

Reserve one private variable in each orientation. For an increasing
private row of exponent \(d\ge1\), the identity
\[
                   (dk)!=d^{dk} k!\prod_{h=1}^{d-1}(h/d)_k
\]
follows by grouping factors by residue modulo \(d\). In the descending
orientation use \(m-k\) in place of \(k\). The two reserved factorials
\(k!\) and \((m-k)!\) cancel the binomial up to the scalar \(m!\).
The residual Gamma factors of either reserved row have fixed
offsets \(h/d\) and affine argument \(k+h/d\) or \(m-k+h/d\).
The factors \(d^{dk}\) are exponential in the index and disappear
under endpoint normalization.

Thus the continuous normalized log coefficient, without the Gaussian,
is a sum of terms \(m A_m(k/m)\), where
\[
 A_m(s)=m^{-1}
 \left[
 \log\Gamma(m(a+(b-a)s)+\gamma)
 -(1-s)\log\Gamma(ma+\gamma)-s\log\Gamma(mb+\gamma)
 \right].                                                     \tag{3}
\]
The original retained rows have \(\gamma=1\). A reserved row can be
treated either as the sum of its fixed-offset Gamma terms, or as
the difference between its full Gamma row and its unit-exponent
Gamma factor.

## Comparing a fixed index at adjacent degrees

For an unreserved row set \(x(s)=a+(b-a)s\), and
\(h(t)=t/(e^t-1)\) for positive \(t\), with \(h(0)=1\).
The Laplace integral for the trigamma gives
\[
           A_m''(s)=(b-a)^2\int_0^\infty
                          e^{-x(s)y}h(y/m)\,dy.                \tag{4}
\]
Here \(h\) decreases and is convex. Indeed
\[
 h''(t)=\frac{e^t[(t-2)e^t+t+2]}{(e^t-1)^3},
\]
and the numerator in brackets vanishes with its first derivative
at zero and has second derivative \(t e^t\). The decreasing
property follows from the numerator \((1-t)e^t-1\) of \(h'\).

For a reserved private row of exponent \(d\), put \(s_0=s\) or
\(1-s\) according to its orientation. Subtract the reserved unit
factor from the integral for this row. Rescale the integration
variable in the full row by \(y=du\), with \(u\) the old variable.
The kernel of the difference is
\[
                   e^{-s_0 y}\{d h(y/(dm))-h(y/m)\}.            \tag{5}
\]
The function \(v\mapsto d h(v/d)-h(v)\) is nonincreasing, since
its derivative is \(h'(v/d)-h'(v)\le0\). It is nonnegative because
\(h\) is positive and decreasing. Both (4) and (5) therefore
show that \(A_m''(s)\ge0\) is nondecreasing with \(m>0\).
The integrals remain finite at a private endpoint \(x(s)=0\).

For each fixed-offset factor in (3), the exact alternative formula
is \(A_m''(s)=m(b-a)^2\psi_1(m x(s)+\gamma)\).
The decreasing trigamma and \(x(s)\ge0\) show that
\(A_N''(s)\le(N/n)A_n''(s)\) for every such factor, including
every factor in the reserved decomposition. The added Gaussian
has constant symbol second derivative \(\tau\), so we obtain the
same bounds for the entire symbol.

For \(w=k/m\), the Green formula expresses the negative log of
the coefficient as
\[
 -\log c_{m,k}
    =m\int_0^1 G(w,s)A_m''(s)\,ds,\qquad
 G(w,s)=\min(w,s)(1-\max(w,s)),
\]
with \(A_m''\) now the sum including the Gaussian. For fixed \(k\)
the factor \(mG(k/m,s)\) is the minimum of the two functions
\((m-k)s\) and \(k(1-s)\). It is nondecreasing in \(m\).
For \(k<n\) its ratio at \(m=N,n\) is bounded by
\((N-k)/(n-k)\) wherever the lower value is positive, with the
corresponding pointwise inequality also holding at zeros.
Multiplying by the two bounds on \(A_m''\) gives (1).

## A positive Mellin representation for the quotient

For a factor in (3) with \(a>0\), let \(X\) have the
Beta distribution with parameters \(na+\gamma,a\).
Its contribution to \(c_{n,k}/c_{N,k}\) is the product
of a fixed positive scale to the \(k\)-th power and
\[
 \frac{\Gamma(na+\gamma+(b-a)k)}{\Gamma(na+\gamma)}
 \frac{\Gamma(Na+\gamma)}{\Gamma(Na+\gamma+(b-a)k)}
                       =\mathbb E[X^{(b-a)k}].
\]
All these moments exist for \(0\le k\le n\) because the
first Gamma argument is at least \(\gamma>0\). A factor
with \(a=0\) contributes only its deterministic scale.
Multiply the powers of independent \(X\)'s and combine their scales.
The Gaussian quotient \(e^{\tau k^2/(2nN)}\) adds the factor
\(e^Z\), with \(Z\) a centered normal of variance \(\tau/(nN)\).
This proves (2).

The Mellin moments are log-convex by Cauchy–Schwarz and their
first ratio \(c_{n,1}/c_{N,1}\) is at least one by (1), or by
its endpoint statement when \(n=1\). Hence all consecutive
ratios are at least one. Applying the argument to the exchanged
columns proves the assertions in the other orientation. The
inequalities constrain the actual linked coefficient arrays.
Angular separation still requires control of their alternating
phase sums.
