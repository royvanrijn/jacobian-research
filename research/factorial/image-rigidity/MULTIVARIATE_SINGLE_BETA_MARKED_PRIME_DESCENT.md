# Prime-power and nearby-prime criteria for the binomial \(X+zX^2Y\)

Let \(X,Y\) be independent standard exponential variables, which
realize the ordinary factorial functional. Consider the moments
of \(X+zX^2Y\). Equal rows on other variables can be added since
they contribute only nonzero scalars.

The consecutive moments are coprime at every prime-power terminal.
More generally, at an upper order \(N=p^s m\) for an odd prime
\(p>2m-1\), the question has a residual squarefreeness test of
degree \(m\), independent of \(s\). A prime near the terminal on
either side gives another exact test of degree equal to the gap.
These conclusions use the original integer exponents.

## The single beta factor and a squarefreeness reduction

Fix \(N=n+1\) and let \(F_N(z)\) be the normalized degree-\(N\)
polynomial at the lower tilt, using \(X\) of shape \(N\) and
\(Y\) of shape one for their product \(T=XY\):
\[
 F_N(z)=\sum_{k=0}^{N}\binom Nk(N)_k k!\,z^k
       =\sum_{k=0}^{N}c_{N,k}z^k,\qquad
 c_{N,k}=\prod_{j=0}^{k-1}(N^2-j^2).
\]
The lower polynomial at degree \(n\), after a positive scalar
normalization, equals \((1-D/N)F_N\), where \(D=z\,d/dz\).
The upper actual tilt changes the shape of \(X\) from \(N\)
to \(N+1\); it gives \((1+D/N)F_N\). This is the exact
single integer beta factor since
\((N+1)_k/(N)_k=(N+k)/N\). Consequently the two actual
moments share a root exactly when \(F_N\) has a multiple root.

The recurrence for \(c_{N,k}\) also gives
\[
                 F_N=1+z(N^2-D^2)F_N.
\]
At a multiple root \(z\), this says \(z^3F_N''(z)=1\).
Thus every multiple root is double. No such root is real:
the reciprocal polynomial is the monic
\(\mathbb E(w+T)^N\), which is positive for real \(w\) in
even degree and has a strictly positive derivative in odd
degree. All roots checked here are nonzero since \(c_{N,0}=1\).

## A prime-power descent

For \(k\ge1\) the coefficient is
\[
                  c_{N,k}=N\prod_{j=N-k+1}^{N+k-1}j.
\]
Write \(N=pM\). At \(k=p\ell\), \(0\le\ell\le M\), this gives
\[
 v_p(c_{N,p\ell})=2\ell+v_p(c_{M,\ell}),\qquad
 \operatorname{unit}_p(c_{N,p\ell})
                 =\operatorname{unit}_p(c_{M,\ell})\pmod p.
                                                               \tag{1}
\]
The multiples of \(p\) are exactly the factors of the compressed
coefficient after dividing the \(2\ell\) copies of \(p\).
The \(2\ell\) full blocks of nonmultiples each have product
\(-1\) modulo \(p\) by Wilson's theorem. Their signs cancel.
The case \(\ell=0\) is the common constant coefficient one.

For \(0\le\ell<M\) and \(0<r<p\), the interval defining
\(c_{N,p\ell+r}\) has exactly the same multiples of \(p\)
as the interval defining \(c_{N,p(\ell+1)}\). In particular
\[
                  v_p(c_{N,p\ell+r})
                           =v_p(c_{N,p(\ell+1)}).
\]
The two outside ordinates increase by
\(2+v_p(M^2-\ell^2)>0\). Hence every interior point of this
block lies strictly above its chord. The lower Newton polygon
therefore retains only the \(p\)-multiple indices, with its
ordinates and units given by (1).

If \(N=p^s m\) and \(p>2m-1\), the coefficients of \(F_m\)
are \(p\)-units. Iterating (1) gives just one slope for \(F_N\),
namely
\[
                   \frac{2(1+p+\cdots+p^{s-1})}{p^s}.
\]
For \(p\) odd its reduced denominator is \(p^s\). The residual
polynomial in degree steps of \(p^s\) is exactly \(F_m\bmod p\).
If it is squarefree, so is \(F_N\). To see this, pass to a
valued extension containing a uniformizer with \(p^s\)-th
power \(p\). Scale by the stated slope, divide each factor
by its minimum coefficient valuation, and reduce. The
initial polynomials multiply. For a factor over
\(\mathbb Q_p\), the initial degrees are in one residue
class modulo the reduced slope denominator, and both endpoints
of the factor survive since all roots have the sole slope.
A repeated positive-degree factor would therefore give a
repeated positive-degree residual factor in the variable
of degree steps \(p^s\).

In particular, for \(m=1\) and odd \(p\), the slope has
denominator \(N\), forcing \(F_N\) to be irreducible over
\(\mathbb Q_p\). For \(N=2^s\), the slope denominator
instead is \(N/2\). A repeated rational factor would have
degree at least \(N/2\), so \(F_N\) would have to be a square
in \(\mathbb Q[z]\), after normalizing its constant
coefficient to one. But its leading coefficient
\(N(2N-1)!\) is divisible exactly once by a prime strictly
between \(N\) and \(2N\), by Bertrand's postulate. Hence
that leading coefficient is no square. This proves the
prime-power assertion, including \(N=2\).

## Two-sided nearby-prime descent

Let \(p\) be an odd prime. Suppose either
\(N=p+t\) with \(0<t<p/2\), or \(N=p-t\) with
\(0<t<N\). In the first case, the interval in the product
formula can cross \(p\) on the left and \(2p\) on the
right, and contains no multiple of \(p^2\). The valuations
are zero up to \(k=t\), one for \(t<k\le p-t\), and
two afterwards. The Newton polygon is flat through \(t\).
If \(4t<p\), it then has two positive slopes with rises one
and widths \(p-2t,2t\). If \(4t>p\), there is just one
positive edge, of width \(p\) and rise two. The equality
\(4t=p\) is impossible for odd \(p\).

In the second case the interval can cross only \(p\), and
does so only when \(k>t\). Thus there is a flat edge through
\(t\) followed by one positive edge of rise one.
In all cases each positive edge has reduced slope
denominator equal to its width. The Newton polygon product
rule then rules out multiple factors at those slopes:
a positive-degree factor of that slope must consume at least
the whole width of the edge.

A rational repeated factor must therefore have all its roots
\(p\)-adic units. Reduce its monic reciprocal, which is integral
because the reciprocal of \(F_N\) is monic integral. The
reduction retains its degree and has a unit constant term.
The unit-root reduction of \(F_N\) is exactly \(F_t\bmod p\),
since \(N^2\equiv t^2\bmod p\), and no later coefficient has
valuation zero. Thus any multiple rational factor has its
square dividing the corresponding unit-root reduction.
If \(F_t\bmod p\) is squarefree, this proves coprimality.
Even without that residual test, the degree of the common
factor of the two ordinary moments is at most
\(\lfloor t/2\rfloor\), using that all multiple roots are
double. In particular, prime gaps \(t\le3\) suffice by the
absence of real common roots.

The residual tests do not assert angular separation under
independent rescaling. The prime criteria leave a remaining
squarefreeness obligation at terminals not reached by them and
the known eventual theorem.
