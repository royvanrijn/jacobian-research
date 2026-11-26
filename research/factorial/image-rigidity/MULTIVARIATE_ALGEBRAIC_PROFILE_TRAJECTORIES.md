# All-degree separation along algebraic exponent trajectories

Fix finitely many exponent rows whose two coordinates are continuous
real algebraic functions of a real parameter \(t\) on an interval
\((t_0,\infty)\). Here algebraic means algebraic over
\(\overline{\mathbb Q}(t)\). Consider only parameters at which all
coordinates are nonnegative integers and the columns are distinct.
Assume these parameters are unbounded. Allowed families include
polynomial and rational profiles with algebraic coefficients as
well as continuous algebraic branches that take integer values on
an unbounded subsequence.

**Theorem.** After removing equal rows, split the family into the
fixed rows \(C\) and the remaining rows. For the latter put
\[
  \tau(t)=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
  h(t)=\sum_i(b_i-a_i).
\]
If there are infinitely many failures of ordinary SFC(2) on the
family at parameters tending to infinity, all of the following
conditions hold.

1. The fixed part \(C\) contains an unequal row, and the changing
   rows are eventually shared with \(\tau(t)\to0\).
2. Their signed gap \(h(t)\) is eventually one fixed integer \(h\)
   at the ordinary parameters.
3. At some fixed adjacent orders \(n,N=n+1\), the core polynomials
   \(Q_n(z)\) and \(Q_N((N/n)^h z)\) have a common root in
   characteristic zero.

If the fixed part has no such integer-twist contact, the ordinary
profiles on the entire trajectory satisfy SFC(2) for all sufficiently
large admissible parameters, simultaneously at every order.
At a possible contact of simple core roots one also has the
first-order arithmetic-response constraint below. It can rule
out vanishing-curvature trajectories even when their fixed cores
have exact twist contacts.

## The bounded and unbounded limits

A continuous real algebraic function has a limit in the extended
real line at \(+\infty\). Any finite limit is algebraic.
Indeed set \(u=1/t\) in a defining polynomial, clear denominators
and divide out the largest common power of \(u\). The resulting
polynomial \(H(u,y)\) has \(H(0,y)\ne0\) as a polynomial in \(y\).
Every finite accumulation value is a root of this nonzero polynomial
over the algebraic numbers. If the extended limit failed to exist,
continuity and the intermediate value theorem would give an interval
of finite accumulation values. Sums, ratios with nonzero denominator,
and maxima of finitely many continuous real algebraic functions
remain algebraic. The assertion for a maximum follows by taking
the product of the defining polynomials of its possible values.

Each bounded exponent coordinate is therefore eventually a fixed
integer on the admissible parameters. For a difference of two
coordinates, divide a polynomial relation by its maximal power
of the function variable. At every nonzero value the remaining
polynomial vanishes, and its specialization at zero is a
nonzero polynomial in the parameter. A boundary between zero
and nonzero values must therefore be a zero of that specialization,
by continuity. There are only finitely many such boundaries,
so each row is eventually equal or eventually unequal. After omitting equal rows and retaining the fixed
coordinates on the admissible parameters, every nonfixed row has
\(\max(a_i,b_i)\to\infty\). If there are no changing unequal rows,
the theorem simply reduces to the fixed pair with \(h=0\) and
an empty tail.

If \(\tau\to\infty\), the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives all-degree separation on the trajectory. If \(\tau\) has
finite limit, every changing unequal row satisfies
\[
 \frac{|b_i-a_i|}{\max(a_i,b_i)}
       \le\sqrt{\frac{\tau(t)}{\max(a_i,b_i)}}\longrightarrow0 .
\]
In particular it is shared at large parameters. If the fixed
part has no unequal rows, the
[absolute diffuse theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md)
now gives all-degree angular separation. For a nonempty fixed
core and positive finite limit of \(\tau\), that limit is algebraic.
The [arbitrary-core Gaussian transfer](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md)
has only finitely many exceptional positive strengths for each
core, all transcendental. A compact interval around the
algebraic limit therefore avoids these exceptions. The uniform
tolerance of the transfer handles the entire trajectory at large
parameters. These arguments establish the first necessary condition
when there is a nonempty changing tail.

In the remaining case \(\tau\to0\), use
[integer-twist stability of the fixed core](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md).
Its finite exact set \(E_C\subset\mathbb Z\) contains only integers
for which a core twist contact occurs. If \(h(t)\to\pm\infty\),
it eventually avoids \(E_C\). If \(h(t)\) is bounded, its limit is
finite and its values on the admissible tail are eventually a
single integer. Hence infinite failures force that integer to
belong to \(E_C\). The
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
also makes it possible to pass to a subsequence of actual failures
at one fixed terminal. The coefficient convergence and the bounded
core root neighborhoods then give a twist contact at that
same terminal, as in the proof of the stability theorem.

## Signed radial response in the vanishing tail

For any collection of weak shared rows, let
\(M_i=(a_i+b_i)/2\) and put
\[
 W=\sum_i d_i^2/M_i,\qquad \mu=\frac{\sum_i d_i/M_i}{W},
 \qquad d_i=b_i-a_i .
\]
The number \(\mu\) is a weighted average of
\(\operatorname{sgn}(d_i)/|d_i|\), with weights \(d_i^2/M_i\).
Thus \(-1\le\mu\le1\). Along the algebraic trajectory in the
theorem, a nonempty changing tail gives an algebraic limit
\(\mu_*\) of this bounded ratio. If all the changing gaps diverge
in absolute value, then \(\mu_*=0\).

Write \(D_m\) for the tail endpoint scale, so that the augmented
polynomial has roots \(w_m=D_m z\) in the tail-normalized
variable \(w_m\). At fixed orders \(n,N\), as \(W\to0\) one has
\[
 \log(D_N/D_n)-h\log(N/n)
                  =-\frac{W\mu}{2nN}+O_{n,N}(W^2).            \tag{1}
\]
Indeed, the order-\(m\) log scale is
\(\sum_i\int_{a_i}^{b_i}\psi(1+m x)\,dx\), with signed integrals.
The elementary trigamma estimate, for \(y>0\), is
\[
 0\le \psi_1(1+y)-y^{-1}+(2y^2)^{-1}\le (6y^3)^{-1}.
\]
It follows by integrating the bound
\(0\le u/(e^u-1)-1+u/2\le u^2/12\) against \(e^{-yu}\).
For the scalar bound use
\(1\le v\coth v\le1+v^2/3\), \(v=u/2\). Multiplication
by \(\sinh v\) proves both inequalities by the power
series of \(\sinh,\cosh\). Thus
\(\int_n^N[x\psi_1(1+mx)-1/m]\,dm
=-(2nN x)^{-1}+O_{n,N}(x^{-2})\).
The scale difference is consequently
\(-\sum_i\log(b_i/a_i)/(2nN)+
O_{n,N}(\sum_i |d_i|/\min(a_i,b_i)^2)\).
For small \(W\), \(|d_i|/M_i\le W\) and \(1/M_i\le W\).
The error sum is \(O(W^2)\). Likewise replacing
\(\log(b_i/a_i)=\int_{a_i}^{b_i}x^{-1}\,dx\) by \(d_i/M_i\)
costs \(O(W^2)\) in total. This gives (1) uniformly in the
number and sizes of the rows.

Let \(\zeta_n,\zeta_N=(N/n)^h\zeta_n\) be simple core roots at
a limiting contact. At a simple root define
\[
 \alpha_m(\zeta)=\frac12-\frac{1}{2m}
                           \frac{\mathcal D^2 Q_m(\zeta)}
                                {\mathcal D Q_m(\zeta)},\qquad
                              \mathcal D=w\,d/dw.
\]
The [weak-tail root response](MULTIVARIATE_MATCHED_TANH_SIMPLE_CORE_ANGULAR_RESPONSE.md#the-simple-root-response)
gives \(\log(w_m/\zeta_m)=W\alpha_m(\zeta_m)+O(W^2)\)
with the logarithm near one. A physical common root has
\(w_N/w_n=D_N/D_n\). Subtract the two local root responses,
use (1), and divide by \(W\). Along a sequence of such
contacts this yields the exact necessary equality
\[
                  2nN\bigl(\alpha_N(\zeta_N)
                        -\alpha_n(\zeta_n)\bigr)=-\mu_* .       \tag{2}
\]
This uses both the modulus and the phase. At any fixed core,
even without an algebraic trajectory, every limiting value
of \(\mu\) is real and belongs to \([-1,1]\).
If the coordinate functions are rational functions over
\(\mathbb Q\), then \(\mu_*\) is rational, since it is a
finite limit of a rational function over \(\mathbb Q\).
An irrational value on the left of (2) therefore excludes
every such rational trajectory at that contact.

**Consequence.** Consider an arbitrary fixed ordinary core, and
suppose all its integer-twist contacts below the absolute
angular cutoff are simple in both degrees and satisfy
\[
        2nN\bigl(\alpha_N(\zeta_N)-\alpha_n(\zeta_n)\bigr)
                                                \notin[-1,1].
\]
Then every sufficiently weak nonempty ordinary shared tail
gives physical SFC(2) in every degree, with no condition on the
tail's total signed gap. This follows by taking a convergent
subsequence of \(\mu\) for any proposed sequence of failures
with weak curvature tending to zero. The integer-twist
stability theorem bounds the possible twists and fixed orders.
Equation (2) excludes the limiting simple contact. Each test
here is an algebraic comparison on the exact core. A uniform
positive curvature tolerance can also be computed by root
isolation and the uniform Taylor errors in (1) and in the
root-response expansion.

For an exact local control of the radial condition, take the
positive-coefficient pair
\[
 A(w)=w^2+2w+5/4,\qquad
 B(w)=\bigl(w+\lambda/2\bigr)
            \bigl(w^2+2\lambda w+5\lambda^2/4\bigr),\qquad
 \lambda=3/2 .
\]
At orders \(2,3\) the integer twist \(h=1\) gives the contact
\(\zeta=-1+i/2,\ \lambda\zeta\).
Since \(B(\lambda w)=\lambda^3(w+1/2)A(w)\),
\[
 2\cdot2\cdot3\,(\alpha_3(\lambda\zeta)-\alpha_2(\zeta))
       =1+\zeta A''(\zeta)/A'(\zeta)-4\zeta/(\zeta+1/2)=-4 .
\]
In this control the angular responses agree and the real response
exceeds the permitted drift in magnitude. The pair consists of
rational polynomials and is a control for the method, outside
the ordinary factorial core class.
