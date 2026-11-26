# All-degree stability of an ordinary core under small shared tails and integer twists

Fix an ordinary support with two distinct exponent columns, called
the core. Append any finite collection of shared rows \((a_i,b_i)\),
with \(a_i,b_i\) positive integers, and set
\[
 \tau=\sum_{\rm tail}\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 \Delta=\sum_{\rm tail}(b_i-a_i)\in\mathbb Z.
\]
There is an explicitly computable \(\tau_C>0\) and a finite,
exactly decidable set \(E_C\subset\mathbb Z\) such that
\(\tau<\tau_C\) and \(\Delta\notin E_C\) imply consecutive
coprimality in every degree. There is no restriction on orientations
or the number of rows in the tail. The core may have zero, one or two
private sides. The computation only uses the core polynomials and
an absolute terminal cutoff. If the core has two private sides and
comparable shared rows, orient the core with \(a_i\le b_i\) on its
shared rows. In that orientation \(E_C\subset\mathbb Z_{<0}\).

## The finite set and the stability constant

Write the raw core moments as
\[
 Q_m(z)=\sum_{k=0}^{m}\binom mk
     \prod_{\rm core}((m-k)a_i+kb_i)!\,z^k.
\]
Write \(c_m=Q_m(0)\), \(l_m=[z^m]Q_m\),
\(g_m=(l_m/c_m)^{1/m}>0\), and \(F_m(u)=Q_m(u/g_m)/c_m\).
Log-convexity of the Gamma function gives
\(0<[u^k]F_m\le\binom mk\). The endpoints are one.

Take an integer \(J\ge4\) at least the
[computable absolute ordinary terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The argument uses this fixed cutoff without evaluating it.
If the core curvature \(S_C=\sum_{\rm core}(b_i-a_i)^2/\max(a_i,b_i)\)
is at least \(27\), the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
already applies to every augmentation. In this case set \(E_C=\varnothing\)
and take any positive \(\tau_C\). Equal rows, including \((0,0)\),
may be dropped throughout.

Suppose \(S_C<27\). Put \(R=729^{27}3^{392}\). For a both-private core
one may instead put \(R=2\). If a support has total curvature below
\(27\), the [uniform annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
bounds all endpoint-normalized root moduli between \(1/R\) and \(R\):
its total private degree is less than \(27\), and its
\(V_{\rm sh}\le S(S+2)<783\). For the both-private bound, cancellation of
the binomial coefficient by [two private reservations](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
makes every normalized coefficient at most one. At radius \(2\)
the leading term dominates their sum, and coefficient reversal
gives the lower bound. This applies to the core and its augmentation.

For \(2\le n<N=n+1<J\), write \(\Lambda=N/n\), \(C_n=g_N/g_n\), and set
\[
 \begin{split}
 I_n&=\{q\in\mathbb Z:(2R^2)^{-1}\le C_n\Lambda^q\le2R^2\},\\
 E_n&=\{q\in I_n:\gcd_{\mathbb Q[z]}(Q_n(z),Q_N(\Lambda^qz))\ne1\},
       \qquad E_C=\bigcup_{2\le n,\ n+1<J}E_n .
 \end{split}                                                   \tag{1}
\]
The end comparisons for \(I_n\) are rational after taking the
\(nN\)-th power, since \(C_n^{nN}=(l_N/c_N)^n(c_n/l_n)^N\).
Successive doubling and binary search find its finite endpoints.
A common core root at any integer twist already forces
\(R^{-2}\le C_n\Lambda^q\le R^2\), so the restriction to \(I_n\)
does not omit any low-order twist with zero resultant.

There are at most \(\lfloor n/2\rfloor\lfloor N/2\rfloor\)
elements of \(E_n\). A twist is positive as a rescaling of the
argument. A common root is nonreal since the even factorial moment
is strictly positive at real arguments, as the expectation of
an even power of a nonzero real polynomial on independent
exponential inputs. Upon conjugating, the two core roots are
therefore in the upper half-plane, and a pair of such roots
determines at most one \(q\). In the both-private case the
[angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
allow only equal labels \(k\) with \(4k<n+2\) and give one simple
upper root per label. Hence \(|E_n|\le\lfloor(n+1)/4\rfloor\).
In the degenerate geometric case all \(E_n\) are empty because
the arguments of the nontrivial \((n+1)\)-th and \((n+2)\)-th
roots of unity are disjoint.

For every \(q\in I_n\setminus E_n\), obtain rational Bézout
polynomials \(A_{n,q},B_{n,q}\) satisfying
\[
 A_{n,q}(z)Q_n(z)+B_{n,q}(z)Q_N(\Lambda^qz)=1.
\]
For \(H\ge0\) write \(\|\sum_j a_jz^j\|_H=\sum_j|a_j|H^j\).
Put \(H_n=R\max(1,c_n/l_n)\), \(L=2J(1+2R)^J\), and
\[
 K_C=\max\left(\{1\}\cup
 \bigl\{c_n\|A_{n,q}\|_{H_n}+c_N\|B_{n,q}\|_{H_n}
    : 2\le n<N=n+1<J,\ q\in I_n\setminus E_n\bigr\}\right),
 \qquad
 \tau_C=\min\{1/10,1/(2LK_C)\}.                               \tag{2}
\]
These quantities are rational and finite. Formula (2) requires no
sampling in the large exponents of the tail.

## Proof of stability

The eventual theorem handles \(N\ge J\) and the curvature theorem
handles total curvature at least \(27\). The first root is real
and the second factorial moment is strictly positive on the real
axis, so it suffices to consider \(2\le n<N=n+1<J\) and total
curvature below \(27\).

Let \(D_m\) be the tail's positive endpoint scale, obtained as
the \(m\)-th root of its ratio of leading and constant factorial
weights. Write \(t_{m,k}\) for its normalized factorial weights.
The full endpoint-normalized polynomial has coefficients
\(([u^k]F_m)t_{m,k}\). The
[tail bounds](MULTIVARIATE_COMPARABLE_CORE_EFFECTIVE_TAIL_STABILITY.md#one-physical-comparison-point)
hold without an orientation restriction:
\[
 0\le1-t_{m,k}\le m\tau/4,\qquad
 \log(D_N/D_n)=\Delta\log\Lambda+\epsilon,\qquad
 |\epsilon|\le\frac{\tau}{nN(1-\tau)}\le\tau .                  \tag{3}
\]
For the first bound the second derivative of the log factorial
weight, with respect to \(k\), is at most \(\tau/[m(1-\tau)]\).
Apply the endpoint chord gap bound and \(\tau\le1/10\).
The center bound follows by writing the log factorial endpoint
ratio as an integral of the digamma function over each row's
oriented exponent interval. The \(m\)-derivative of the digamma
term with its \(\log m\) subtracted has absolute value at most
\(1/(m^2\min(a_i,b_i))\). Integration from \(n\) to \(N\) and
\(\sum_i |b_i-a_i|/\min(a_i,b_i)\le\tau/(1-\tau)\) give (3).

If \(z\) is a full common root, let \(u=g_nD_nz\) and \(u'=g_ND_Nz\).
Both lie in the annulus \([1/R,R]\) in modulus. By (3),
\(u'/u=C_n\Lambda^\Delta e^\epsilon\). Since \(e^{|\epsilon|}<2\),
this forces \(\Delta\in I_n\). The coefficient discrepancy implies
\[
 |F_n(u)|,\ |F_N(u')|\le J\tau(1+R)^J/4 .
\]
Put \(v=C_n\Lambda^\Delta u=e^{-\epsilon}u'\). On the radial segment
from \(u'\) to \(v\), \(|w|\le2R\) and
\(|wF_N'(w)|\le J(1+2R)^J\). Integrating over log radius gives
\(|F_n(u)|,|F_N(v)|\le L\tau\).

In the raw argument \(\zeta=u/g_n\) we have
\(|\zeta|\le H_n\), \(Q_n(\zeta)=c_nF_n(u)\), and
\(Q_N(\Lambda^\Delta\zeta)=c_NF_N(v)\). If \(\Delta\notin E_C\),
the Bézout identity gives \(1\le K_CL\tau<1/2\). This contradiction
proves the theorem over \(\mathbb C\) and hence over every
characteristic-zero field by rational coprimality.

## The exceptional twists for comparable both-private cores

In the comparable orientation all \(q\ge0\) can be omitted.
Use the geometric representation \(Q_m(z)/c_m=\mathbb E_m G_m(zT)\),
where \(G_m(w)=1+w+\cdots+w^m\), \(T=V/U>0\), and
\(\mathbb E_m\) is the normalized tilt by \(U^m\).
Let \(\mu\) be the law of \(T\) under \(U^n\) and \(\nu\) its law
under \(U^nV\). The
[conditional comparison](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md#why-comparable-shared-ordinary-rows-satisfy-the-test)
shows \(d\nu/d\mu\) is nondecreasing. The law of \(\log T\) under
\(\mu\) has a log-concave density, as a convolution of scaled
independent log-Gamma densities after reservation. Consequently
the law \(\mu_\kappa\) of \(\kappa T\) under \(U^n\), for
\(\kappa=\Lambda^q\ge1\), is above \(\mu\) in likelihood-ratio
order. This follows directly because \(f(t-\log\kappa)/f(t)\)
is nondecreasing for a log-concave density \(f\). Scaling also
preserves the likelihood-ratio order of \(\nu\) and \(\mu\).
Thus the law \(\nu_\kappa\) of \(\kappa T\) under \(U^nV\) is
above \(\mu\). The constant-\(T\) case was handled by the arguments
of the geometric sums.

At a twisted common root \(z=r e^{i\theta}\), conjugation and
the angle boxes put both roots in one positive label. The
[geometric imaginary part](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md#the-conditional-monotonicity-test-at-the-physical-root)
\(I(s)=\operatorname{Im}G_n(s e^{i\theta})\) has a single change
from positive to negative. The lower equation gives
\(\mathbb E_\mu I(rT)=0\). Monotonicity of the density ratio
across the change of sign gives \(\mathbb E_{\nu_\kappa}I(rT)\le0\).
On the other hand, with \(X=\kappa T\) and
\(b=\mathbb E_N X>0\),
\[
 \mathbb E_N G_N(zX)=1+zb\,\mathbb E_{\nu_\kappa}G_n(zX)=0.
\]
Its consequence \(\mathbb E_{\nu_\kappa}I(rT)=
\sin\theta/(rb)>0\) is a contradiction. This proves
\(E_C\subset\mathbb Z_{<0}\) in the stated orientation.

## An arithmetic filter for the negative twists

Suppose there is exactly one private row of degree one on each side,
and the shared core rows satisfy \(b_i\ge a_i\). Dividing by the
constant coefficient after the two private cancellations leaves
integer coefficients because each shared weight ratio
\((ma_i+k(b_i-a_i))!/(ma_i)!\) is an integer. Thus
\(R_m(x)=x^m Q_m(1/x)/c_m\) is monic integral. A common root at
a twist \(q=-h<0\) implies that both \(x=1/z\) and
\((N/n)^h x\) are algebraic integers. A linear combination with
integer coefficients, using \(\gcd(n^h,N^h)=1\), makes
\(x/n^h\) an algebraic integer. Equivalently, the two polynomials
\(R_n(n^h y)\) and \(R_N(N^h y)\) share an algebraic integer root.
For every prime one may divide each polynomial by the largest
common power of that prime in its coefficients before reducing.
If the two reductions are coprime for any prime there is no
such twist, because the monic minimal polynomial of an algebraic
integer divides both primitive polynomials over the integers.

Suppose all unequal shared rows have gap one, and there are
\(D\ge1\) such rows. A negative twist \(q=-h\) with \(h\ge D\)
is impossible whenever either \(n\) or \(N\) is a prime power.
In particular for a single shared unit-gap row this excludes
all negative twists at such terminals. To see this, let
\(m=p^s\) be the prime-power order and write \(d_j=[x^j]R_m\).
The coefficient ratios are
\[
                      \frac{d_j}{d_{j+1}}
                                    =\prod_i(m(a_i+1)-j),
                      \qquad 0\le j<m,
\]
where the product runs over the unit-gap rows. Thus
\(v_p(d_0/d_1)\ge Ds\), while \(v_p(d_j/d_{j+1})=D v_p(j)<Ds\)
for \(1\le j<m\). On the Newton polygon the line through
\((0,v_p(d_0))\) and \((1,v_p(d_1))\) lies strictly below
every point to its right. Every later slope is greater than
\(-Ds\), since every slope between the remaining coefficient
points is greater than \(-Ds\). There is at most one root with
valuation at least \(Ds\), counted with multiplicity. But the
algebraic integrality proved above, at order \(m\), forces
every conjugate of the corresponding reciprocal root to have
valuation at least \(hs\ge Ds\). Its minimal polynomial must
have degree one, contradicting the nonreal common-root condition.

## Scope of the exception

For each prescribed integer \(q\), tails with this same \(\Delta=q\)
and curvature tending to zero exist. For \(q\ne0\), take \(|q|\)
oriented unit-gap rows and let their smaller exponents tend to
infinity. For \(q=0\) use an oppositely oriented pair. Equation (3)
shows that, after using the lower tail scale as the common unit,
their pairs converge coefficientwise to
\(Q_n(z)/c_n,\ Q_N(\Lambda^q z)/c_N\) at each fixed order.
Thus if \(q\in E_C\), the limiting resultant at a witnessing
order is zero. The theorem leaves those actual small-tail pairs
open. A vanishing limiting resultant alone does not give a
common root for any positive tail curvature.
