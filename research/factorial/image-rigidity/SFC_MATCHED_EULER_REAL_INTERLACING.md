# Real angular interlacing under matched Euler factors

Write \(E=x\,d/dx\). Suppose \(f,g\) have positive leading
coefficients, degrees \(n,n+1\), parity respectively \(n,n+1\),
and simple strictly interlacing real zeros. Let \(B,C>0\).
If
\[
                               |C-B|\le1 ,
\]
then \((E+B)f\) and \((E+C)g\) also have simple strictly
interlacing real zeros. This shift interval is sharp for a
statement about all symmetric interlacing pairs. The theorem
can be iterated for any number of matched factors.

In particular, let \((X,Y)\) be a common-shape correlated
Gamma pair with fixed squared correlation \(0\le\rho<1\).
For every fixed nonnegative integer tilt \(r\), the polynomials
\[
            \mathbb E[(XY)^r(X+zY)^n],\qquad n\ge1 ,
\]
have strict consecutive angular interlacing in every degree.
One fixed positive tilt is used at all moment orders.

## A sign at every transformed lower root

Normalize \(f,g\) to be monic. Strict interlacing gives
\[
                  R(x)=\frac{g(x)}{f(x)}
                      =x-\sum_\ell\frac{w_\ell}{x-\xi_\ell},
                     \qquad w_\ell>0 ,
\]
where the \(\xi_\ell\) are the zeros of \(f\).
The constant term in the polynomial part is zero by parity.
The opposite poles have equal weights because \(R\) is odd.
A pole at zero is present just when \(f\) is odd.

The transformed lower polynomial \(F=(E+B)f\) has simple
real zeros, symmetric about zero. On the positive line its
zeros are distributed as follows. With \(f=x^\epsilon p(x^2)\)
and positive squared roots \(0<u_1<\cdots<u_k\), the ratio
\(F/f=B+\epsilon+2\sum_i x^2/(x^2-u_i)\) decreases strictly
on each positive pole-free interval, and has one zero before
each pole. This follows from the negative derivative in the
variable \(x^2\). At any positive zero \(x_0\) of \(F\),
the sign of \(F'(x_0)/f(x_0)\) is therefore negative.

For \(G=(E+C)g\) and \(h=C-B\), compute at such a zero
\[
                   \frac{G(x_0)}{f(x_0)}
                        =x_0 R'(x_0)+h R(x_0).
\]
The polynomial part contributes \((1+h)x_0\).
For each pair of poles \(\pm\xi\), \(\xi>0\), their combined
contribution is
\[
        2w x_0
             \frac{(1-h)x_0^2+(1+h)\xi^2}{(x_0^2-\xi^2)^2}.
\]
A central pole contributes \((1-h)w_0/x_0\).
All terms are nonnegative for \(-1\le h\le1\), and the
sum is strict whenever a positive zero exists. In fact
\(1+h>0\) except at \(h=-1\), where the noncentral pair
contribution is strict. Thus \(G(x_0)/F'(x_0)<0\).
The same ratio is negative at each negative zero by parity.
If \(f\) is odd, the ratio at its retained central zero is
\(C g(0)/[(B+1)f'(0)]<0\), by the original residue.
These signs, together with the positive leading coefficients,
give exactly one simple zero of \(G\) in each strict
interlacing interval of \(F\).

This interval of matched shifts cannot be enlarged. If \(C-B>1\),
use \(f=x(x^2-\eta)\) and \(g=(x^2-u)(x^2-v)\), where
\(0<u<\eta<v\). The positive zero of \((E+B)f\) lies at
\(x_0^2=(B+1)\eta/(B+3)\). At \(u=\eta\), the coefficient
of \(v\) in \((E+C)g(x_0)\) is
\(2\eta(C-B-1)/(B+3)>0\). Taking \(v\) sufficiently large
and then \(u\) just below \(\eta\) gives the wrong sign for
strict interlacing at \(x_0\).
If \(C-B<-1\), use \(f=x^2-\eta\), \(g=x(x^2-v)\), with
\(v>\eta>0\). The transformed positive roots have squares
\(B\eta/(B+2)\) and \((C+1)v/(C+3)\), respectively.
For \(v\) sufficiently close to \(\eta\), the larger degree
has the smaller outer root. Positive rescaling puts all
input roots strictly inside \((-1,1)\) in both constructions.

## A degree-linked multiplier and its physical realization

Take \(\alpha>0\) and any polynomial \(L(q)\) whose zeros are
strictly negative real numbers. For \(m\ge1\), set
\[
             H_m(x)=L\big((m+E)/2\big)\, C_m^\alpha(x),
\]
where \(C_m^\alpha\) denotes the Gegenbauer polynomial.
A nonzero constant in \(L\) only rescales each polynomial.
For every factor \(q+\beta\) of \(L\), the Euler parameters
in the two degrees are \(m+2\beta\) and \(m+1+2\beta\).
They differ by one. The original Gegenbauer pair interlaces
by its fixed positive orthogonality weight. Applying the
matched Euler theorem successively proves strict interlacing
of every \(H_m,H_{m+1}\). Each individual operator contracts
the positive roots at its own degree, so all these real roots
are strictly between \(-1\) and \(1\).

With \(\delta=1-\rho\) and \(\alpha=\gamma+r\), the
[physical Meixner factorization](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md#meixner-factors-and-contraction-of-real-roots)
gives precisely this family, with
\[
                          L(q)=M_r(-\alpha-q;\gamma,\delta).
\]
Its linear factors are \(q+\alpha+s_i\), up to a nonzero
constant, with \(s_i>0\). They are independent of the moment
order because the tilt \(r\) and the law are fixed.
The angle map \(z=\exp(2i\arccos(\sqrt{\delta}\,x))\) is
the same in both degrees and preserves strict alternation.
At \(\rho=0\) the angular polynomial is one fixed-parameter
Gegenbauer family. This proves the claimed physical
interlacing on the noncollapsed interval.

There are also other positive mixtures with this exact
degree link. Let \(\delta\in(0,1]\) and consider a positive
law with transform near the origin
\[
 F(u+v-\delta uv),\qquad
 F(w)=\sum_{k\ge0}\frac{(\alpha)_k L(k)}{k! L(0)}w^k .
\]
Its order-\(m\) angular polynomial, up to a nonzero scalar
and the same angle map, is \(H_m\). Indeed the summand of
index \(k\) in \(F(2xt-t^2)\) contributes to \([t^m]\)
only a term of \(x\)-degree \(2k-m\). Multiplying this term
by \(L((m+E)/2)\) is exactly multiplication by \(L(k)\).
For an explicit positive family, choose
\(L(q)=\prod_{i=1}^d(q+\beta_i)\) with
\(\beta_i\ge\alpha+i-1\). Then \(L(k)\) is a nonnegative
linear combination of \((\alpha+k)_j\), by iterating
\[
 (k+\beta)(\alpha+k)_j
      =(\alpha+k)_{j+1}+(\beta-\alpha-j)(\alpha+k)_j.
\]
Since \((\alpha)_k(\alpha+k)_j=(\alpha)_j(\alpha+j)_k\),
the transform is a positive mixture of common-shape correlated
Gamma transforms with shapes \(\alpha+j\). The same argument
gives all-degree interlacing for this whole family.

In the ordinary mirrored profile \(X^aY^a(X+zY)\) the tilt
is \(r=am\) at order \(m\). The Meixner polynomial then changes
degree as well as argument. Applying the paired theorem across
those different factor lists still requires an additional
comparison.
