# Compactification of linked Gamma factors for angular separation

At any fixed pair of consecutive orders, angular separation for arbitrary
positive Gamma shapes reduces to a sequence of finite algebraic tests. The retained
factors keep their individual shape ratios and orientations across orders.
The discarded factors contribute two linked Gaussian parameters with an
explicit error bound. The certificate succeeds at some cutoff exactly when
the corresponding compactified family has no shared root argument.

The family below contains the Gauss-factor expansions of arbitrary ordinary
two-monomial supports. It also allows independent real bases and shape
ratios for the factors. Angular separation of this larger family implies
ordinary coprimality. A collision in the larger family leaves ordinary
coprimality undecided.

## Linked factors and comparison across orders

Fix \(n\ge2\), set \(N=n+1\), and divide the coefficient at index \(k\)
by \(\binom jk\) for \(j=n,N\). Normalize the constant and leading
coefficients to one by independent positive scalings of the argument.
An increasing row with smaller exponent \(a\), gap \(d>0\), and
Gamma shape \(\gamma>0\) splits into \(d\) factors of bases
\((ja+\gamma+i)/d\), \(0\le i<d\). A decreasing row reverses the
coefficient arrays. Thus it suffices to allow any factor with bases
\(T>0\) and \(\rho T\), where \(1\le\rho\le N/n\), and either orientation.
The two arrays in increasing orientation are
\[
 F_{j,k}(x_j)=
 \frac{\prod_{i=0}^{k-1}(1+ix_j)}
      {(\prod_{i=0}^{j-1}(1+ix_j))^{k/j}},
 \qquad x_n=1/T,\quad x_N=1/(\rho T).
\]
The arrays multiply over factors. All entries are at most one by
log convexity. Every factor with \(\rho<N/n\) is itself realized by
a real smaller exponent \(a=T(\rho-1)\), unit gap, and Gamma shape
\(\gamma=T-na>0\). The endpoint \(\rho=N/n\) is its limit with
positive shapes.

For the lower order let \(u\) be the outer geometric mean of the
normalized array. Let \(w\) denote the corresponding upper-order
mean. A factor contributes to \(-\log u\) the quantity
\[
 s=\tfrac12\log(1+(n-1)/T),
\]
since \(F_{n,1}F_{n,n-1}=(1+(n-1)/T)^{-1}\).
The upper-order factor contributes \(\tfrac12\log(1+n/(\rho T))\).
In particular the whole array has
\[
 u^{n/(n-1)}\le w\le
      \sqrt{\frac{(n^2-1)u^2}{n^2-u^2}}\le u.                 \tag{1}
\]
Indeed the ratio \(n/((n-1)\rho)\) lies in
\([n^2/(n^2-1),n/(n-1)]\).
Use \(1+c(q-1)\le q^c\) for the power comparison and
\(\prod_i(1+c(q_i-1))\ge1+c(\prod_i q_i-1)\) for the upper
bound on the mean.

A comparison also holds at every fractional index. Define the signed
log chord error for a factor at order \(j\) by
\[
 \ell_j(t)=
 \log\Gamma(T_j+jt)-(1-t)\log\Gamma(T_j)-t\log\Gamma(T_j+j),
 \qquad 0\le t\le1.
\]
In decreasing orientation use \(\ell_j(1-t)\). Sum over all factors.
Then
\[
 \ell_n(t)\ge \ell_N(t)\ge (N/n)^2\ell_n(t).                 \tag{2}
\]
More generally, every derivative of order \(p\ge2\), with its fixed
orientation sign removed, increases by a factor in
\([1,(N/n)^p]\) when passing from \(n\) to \(N\).
For increasing orientation the derivative with its sign removed is
\((p-1)!j^p\sum_{k\ge0}(T_j+jt+k)^{-p}\). Every upper-order
denominator lies between the lower-order denominator and
\(N/n\) times that denominator. Reflection has no effect on the
sign at even \(p\), so the even derivative bounds also hold for
mixtures of orientations. At \(p=2\) both endpoint chord errors
vanish. Integrating the curvature comparisons, or comparing two
convex functions with zero endpoint values, gives (2).

## The compactified coefficient family

Fix \(0<u_0\le u_1<1\), and constrain \(u\in[u_0,u_1]\).
A compactified profile consists of a finite or countable list of
factors of parameters \(s_i>0,\rho_i\in[1,N/n]\), and orientations,
and two numbers \(\tau,\nu\) with
\[
 S=\sum_i s_i+\tau,\qquad e^{-S}\in[u_0,u_1],\qquad
 \tau\ge0,\qquad n\tau/N\le\nu\le\tau.
\]
The array is obtained by multiplying the exact retained factors
and the Gaussian factors
\[
 e^{-\tau k(n-k)/(n-1)},\qquad e^{-\nu k(N-k)/(n-1)}           \tag{3}
\]
at the two orders, respectively. All infinite products converge
to positive numbers. The compactified profiles give exactly the
closure of the normalized coefficient pairs of finite profiles
in the enlarged linked-factor model.

For the estimate that proves this assertion, let a small factor
have \(s\le1/4\), and put \(x=(e^{2s}-1)/(n-1)\) and
\(x_j=x\) or \(x/\rho\). For \(0\le k\le j-2\), the second
difference in \(k\) of \(\log F_{j,k}(x_j)\) is
\(\log(1+x_j/(1+kx_j))\). Its difference from \(x_j\) is
nonpositive with absolute value at most
\((k+1/2)x_j^2\). Compare with the quadratic
\(-x_j k(j-k)/2\) using the zero endpoint values. This gives
\[
 0\le \log F_{j,k}(x_j)+\tfrac12x_j k(j-k)
                 \le \tfrac12(j-1)x_j^2 k(j-k).
\]
Since \(x\le4s/(n-1)\) and
\(0\le x-2s/(n-1)\le4s^2/(n-1)\), it follows that
\[
 \left|\log F_{j,k}(x_j)+
       \frac{s\,k(j-k)}{(n-1)\rho_j}\right|
 \le \frac{18s^2}{n-1}\,k(j-k),\qquad
 \rho_n=1,\quad\rho_N=\rho.                                 \tag{4}
\]
Reversal preserves this estimate. Consequently the error for
replacing a collection of factors of \(s_i\le\delta\le1/4\)
by (3) is at most \(18\delta\tau k(j-k)/(n-1)\) in every
normalized log coefficient, with
\(\tau=\sum_i s_i\) and \(\nu=\sum_i s_i/\rho_i\).

For necessity, sort each finite list by descending \(s_i\), pad
with zero-mass factors, and pass to a diagonal subsequence on
which every fixed index has limiting mass, shape ratio, and
orientation. Extract limits \(S_*,U_*\) of the total mass and of
\(\sum_i s_i/\rho_i\) as well. The limiting retained masses have
sum at most \(S_*\). Put \(\tau=S_*-\sum_i s_i\) and
\(\nu=U_*-\sum_i s_i/\rho_i\) for the limiting list.
The inequalities on \(\nu,\tau\) follow by applying the same
inequalities to the tails beyond each fixed index. If \(S^*\)
is an upper bound on all total masses, the largest mass in a
tail beyond index \(J\) is at most \(S^*/(J+1)\).
Equation (4) bounds its total log error by a constant depending
on \(n,S^*\) times \(1/(J+1)\). Taking the sequence limit first
and then \(J\to\infty\) gives the stated representation.

For sufficiency, retain the first \(J\) factors and add their omitted
tail budgets to \(\tau,\nu\), so that total mass is preserved exactly.
For a nonzero diffuse budget take many equal factors of mass
\(\tau_J/h\) and common ratio \(\tau_J/\nu_J\).
Letting \(h\to\infty\) gives the desired Gaussian factors by (4).
Then let \(J\to\infty\). Ratios at \(N/n\) can be approached from
below without changing the lower-order masses. These approximations
remain within the same prescribed \(u\)-interval.

## A finite algebraic test at each cutoff

Pick \(h=1/b\le1/2\) with \(b\) a positive integer. Retain precisely
the factors for which \((n-1)x\ge h\). There are at most \(M\) of
them, where \(M\) is the largest integer with
\((1+h)^M\le u_0^{-2}\). This bound can be found by exact integer
comparisons. For each possible retained length and orientation list,
use algebraic variables \(x_i>0,\rho_i\in[1,N/n]\) for their linked
parameters. Every omitted mass is at most \(h/2\). Choose a rational
\(S^*\ge-\log u_0\), for example \(u_0^{-1}-1\).
By (4), a common upper enclosure on all absolute logarithmic errors
is
\[
 \Delta_h=8(n+1)h S^* .
\]
This follows from \(k(j-k)\le(n+1)^2/4\) and
\((n+1)/(n-1)\le3\).

Use the variables \(p=e^{-\tau/(n-1)}\) and
\(q=e^{-\nu/(n-1)}\). Then \(0<p\le q\le1\) and
\(q^N\le p^n\). Equation (3) becomes the polynomial monomial
factors \(p^{k(n-k)}\) and \(q^{k(N-k)}\). Impose the budget by
\[
 u_0\le p^{n-1}\prod_i (1+(n-1)x_i)^{-1/2}\le u_1.
\]
Every retained coefficient can be expressed algebraically by its
positive \(j\)-th root, clearing the positive denominators from the
formula for \(F_{j,k}\). For cutoffs with \(\Delta_h<1/2\), admit
a multiplicative error for each nonendpoint coefficient in the
rational interval \([1-\Delta_h,1+2\Delta_h]\), and retain its
cap by the binomial coefficient. This enclosure follows from
\(e^{-v}\ge1-v\) and \(e^v\le1+2v\) for \(0\le v\le1/2\).
The endpoint coefficients stay one.

The test asks for a common root direction. In real variables
\(a,b_0,R\), put \(z=a+\mathrm i b_0\) and ask
\[
 P_n(z)=P_N(Rz)=0,\qquad R>0 .
\]
The real and imaginary equations are polynomial. Both roots have
moduli between \(2^{-j}\) and \(2^j\) at the corresponding order \(j\).
Indeed, for \(|z|\ge1\) the monic root equation and the binomial
coefficient caps give \(|z|\le2^j-1\), and reversing the polynomial
gives the lower bound. Thus one may impose
\(R\in[2^{-(2n+1)},2^{2n+1}]\). Feasibility is a finite decision
problem over the real algebraic numbers. If every test at a
cutoff is infeasible, all allowed original profiles have disjoint
root arguments. Their original centers can differ by arbitrary
positive ratios.

If a compactified profile has a common direction, retain its
factors above the given cutoff and move the rest into its diffuse
budgets. Equation (4) puts its coefficients into the admitted
error boxes, so every cutoff test is feasible. Conversely, if
arbitrarily fine cutoffs have feasible tests, their unperturbed
base arrays lie in the compactified family. They have a convergent
subsequence by the closure just proved. The perturbations tend
to zero, and the roots and relative radii stay in the stated
bounded annuli. Passing to a further subsequence gives a common
direction for the compactified pair. Hence some cutoff succeeds
exactly when this enlarged compact family has no common direction.

The [ordinary compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md)
keeps the actual ratio of the two centers and the discrete exponent
rows. The present reduction keeps positive scalings arbitrary and
permits individual shape ratios on a continuous interval. This is
suited to a proof by angular separation when the relative centers
are unavailable.

## Explicit large outer loss in every degree

The enlarged Gamma family has disjoint root arguments whenever
\[
 u^{2/N}\le\frac{1}{8nN}.                                    \tag{5}
\]
To prove this, for a single increasing factor set \(g_i=\log(1+ix_j)\)
and \(T_*=1+(j-1)x_j\). Splitting the ordered \(g_i\) below and above
\(k\) gives
\(F_{j,k}\le T_*^{-\min(k,j-k)/j}\).
For \(k\le j/2\), use
\((j-k)\sum_{i<k}g_i\le k\sum_{i=k}^{j-2}g_i\)
by comparing the weights and the largest term on the left to
the smallest on the right. For \(k>j/2\), add
\((2k-j)g_{j-1}\) to the right and use the equality of total
weights. Reversal preserves the bound. After multiplying the
factors the right side is bounded by
\(u^{2\min(k,j-k)/N}\), using \(w\le u\).
The normalized sum of the nonendpoint coefficients is therefore
at most
\[
 2\sum_{k\ge1}(j u^{2/N})^k/k!
 \le\frac{2j u^{2/N}}{1-j u^{2/N}} <\frac1{3n}
\]
under (5). For a root of modulus at most one,
\(|1+z^j|\) is at most this sum. Apply the same argument to the
reciprocal for a root of larger modulus. Each root argument lies
within \(\arcsin(1/(3n))/j\) of an odd multiple of \(\pi/j\).
The order-\(n\) and order-\(N\) reference arguments differ
by at least \(\pi/(nN)\), and the combined errors are less than
\(2(2n+1)/(3n^2N)<\pi/(nN)\). This proves (5).

Together with the
[flat-end separation](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#small-outer-coefficient-loss-for-arbitrary-positive-gamma-shapes),
this confines the angular comparison to an interior outer-budget
band at every fixed degree. The compactification provides a
criterion on that band. It does not establish that every such
test sequence succeeds.

## Euler operator relations for integer exponent rows

For exponent columns \(a,b\) of nonnegative integers and shapes
\(\gamma_i>0\), put \(d_i=b_i-a_i\) and \(\lambda_i=\gamma_i+na_i\).
In this section normalize the constant coefficient to one while
keeping the original variable scale, and call the polynomials
\(P_n,P_N\). With \(D=t\,d/dt\), define
\[
 A(D)=\prod_i\prod_{0\le j<a_i}\left(1+\frac{d_iD}{\lambda_i+j}\right),
 \qquad
 B(D)=\prod_i\prod_{0\le j<b_i}\left(1+\frac{d_iD}{\lambda_i+j}\right).
\]
Let \(c=P_N'(0)/N>0\). The coefficient recurrences give
\[
 P_N-\frac1N D P_N=A(D)P_n,\qquad
 \frac1N P_N'=cB(D)P_n,\qquad
 D A(D)P_n=c\,t(n-D)B(D)P_n.                              \tag{6}
\]
For the first identity, the reverse derivative on the left keeps
the normalized moment coefficient at index \(k\) and replaces
\(\binom Nk\) by \(\binom nk\). Increasing the lower-order shape
\(\lambda_i\) by \(a_i\) multiplies the \(k\)-th moment by
\((\lambda_i+d_ik)_{a_i}/(\lambda_i)_{a_i}\). This is the
\(i\)-th multiplier in \(A(k)\). For the ordinary derivative,
the upper-order normalized moment at \(k+1\), divided by its first
moment, has the shape shift \(b_i\) instead. This gives \(B(k)\).
Eliminate \(P_N\) using the first two identities to obtain the third.

The factors in both operators are positive on the index interval
\(0\le k\le n\), because
\(\lambda_i+j+d_i k\) is positive at both endpoints. Each inverse
factor represents a positive averaging operation on degree at
most \(n\). For \(\lambda=\lambda_i+j\) it is
\[
 (1+d_iD/\lambda)^{-1} f(t)
   =\int_0^1 \lambda x^{\lambda-1} f(t x^{d_i})\,dx .
\]
In the logarithmic variable these are convolutions with one-sided
exponential distributions. Both derivatives of the upper-order
polynomial therefore recover the same lower-order polynomial
after the two prescribed sequences of convolutions and the scalar
normalization in (6). Moreover \(B(k)/A(k)\) strictly increases
with \(k\) for distinct columns. After canceling their common
factors, a row with \(d_i>0\) contributes increasing positive
linear factors, and a row with \(d_i<0\) contributes reciprocals
of decreasing positive linear factors.
