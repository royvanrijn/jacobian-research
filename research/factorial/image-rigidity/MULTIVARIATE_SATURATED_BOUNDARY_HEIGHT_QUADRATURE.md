# Height quadrature at transverse saturated Gaussian contacts

Fix adjacent orders \(n,N=n+1\), with \(n\ge3\). Consider a
[compactified reference contact](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at strength \(T_*>0\) and saturated bias \(\lambda_*=\epsilon T_*\),
\(\epsilon\in\{1,-1\}\). Its two contacting nonreal roots are simple,
and their angular difference has nonzero derivative in strength.
The retained core and the total integer degree difference are fixed.
We use the exact log-coefficient and center conventions of the
[third-order boundary test](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md).
For a finite ordinary core the conclusions concern actual added
ordinary tails continuing the specified contact.

An added row of signed gap \(d_i\ne0\) and midpoint \(m_i\) has
\(u_i=1/m_i\), weight \(w_i=d_i^2u_i\), and signed bias
\(\ell_i=d_i u_i\). Write
\[
 W=\sum_i w_i,\quad \lambda=\sum_i\ell_i,\quad
 \delta=W-\epsilon\lambda,\qquad
 \eta=\max_i |d_i|u_i .
\]
We consider tails with \(W\to T_*,\lambda\to\epsilon T_*,\eta\to0\).
Rows with \(d_i=\epsilon\) are called aligned. Form the probability
measure \(\nu\) on \([0,\eta]\) by assigning each aligned row its
mass \(u_i/W\) at \(u_i\), and assigning the total mass of other rows
to zero.

**Theorem.** There is an analytic boundary response \(\mathcal B\) of
a finite real coefficient vector, and an analytic aligned row
curve \(K(u)\) with \(K(0)=0\), defined below. For every contact
there is an integer \(2\le s\le n+N-1\), determined by the reference,
with the following property. There is a positive probability
quadrature \(\nu^\#\) on \([0,\eta]\) with at most \(N\) nodes, sharing
the moments of \(\nu\) through degree \(s-1\). Every exact ordinary
collision continuing the contact satisfies
\[
 \mathcal B\!\left(\int K\,d\nu^\#\right)
 \ \ge\ \frac{\delta}{4nN}
        +c\!\int_0^\eta u^a(\eta-u)^b q(u)^2\,d\nu(u),           \tag{1}
\]
for some \(c>0\) independent of added rows. Here \(a,b\in\{0,1\}\)
and the quadrature's monic interior-node polynomial \(q\) have
\(2\deg q+a+b=s\) in the nondegenerate case. In the case where
\(\nu\) itself has sufficiently few nodes we take \(\nu^\#=\nu\)
and the final lower bound in (1) is interpreted as zero.

Thus nonpositivity of the response on the specified bounded-node
configurations excludes actual tails at the contact. The reduced
nodes and weights are real. For a finite core one can certify \(s\)
and the quadrature type from finitely many algebraic derivatives.

If the third-order boundary response \(A_3\) of the linked note is
strictly negative, \(s=2\) and the quadrature is the single atom
at \(\mu=\int u\,d\nu\). Write \(B(\mu)=\mathcal B(K(\mu))\).
For all sufficiently diffuse collisions,
\[
             B(\mu)\ \ge\ \frac{\delta}{4nN}
                               +c\,\operatorname{Var}_{\nu}(u).
                                                               \tag{2}
\]
In particular an eventually nonpositive \(B\), including an
identically zero \(B\), excludes the branch. If its first positive
term is of order \(k\), then a collision requires both
\(\delta=O(\mu^k)\) and \(\operatorname{Var}_\nu(u)=O(\mu^k)\).
These bounds are uniform in added row counts, gaps and heights.

For a finite ordinary core the one-height response \(B\) has
a nonzero Taylor coefficient of order at most \(n+N-1\),
regardless of the sign of \(A_3\). In the one-height case
\(A_3<0\) the first nonzero coefficient has order at most two.
In particular the equality \(A_3+T_*D=0\) is impossible there.

## The angular equation with normalized curvature

Use the retained endpoint-normalized pair \(F_j(T,z)\).
Let \(x\) have real coordinates \(x_{j,k}\) for \(1\le k<j\),
\(j=n,N\), and a center coordinate \(x_E\).
Multiply the respective coefficients by \(\exp(x_{j,k})\)
and denote the two local simple roots in their own scales by
\(z_j(T,x)\). Put
\[
 \mathcal H(T,x)=\log\frac{z_N(T,x)}{z_n(T,x)}
                       -C_R+\frac{\epsilon T}{2nN}-x_E ,
\]
where \(C_R\) is the exact retained center including the fixed
total degree difference. The logarithm is chosen real at the
contact, so \(\mathcal H(T_*,0)=0\).
Since \(\partial_T\Im\mathcal H(T_*,0)\ne0\), the real analytic
implicit function theorem gives a unique nearby \(T=T(v)\) with
\[
 \Im\mathcal H(T,Tv)=0,\qquad T(0)=T_*,
 \qquad \mathcal B(v)=\epsilon\Re\mathcal H(T(v),T(v)v).
                                                               \tag{3}
\]
Time in (3) is solved using the exact angular equation.

For a row define the vector \(X_i\) of its exact compensated log
coefficients and center. At index \(k\) of order \(j\) use
\(j(L_{j,i}+w_i\phi)(k/j)\), where \(jL_{j,i}\) is the exact
endpoint-normalized log multiplier and \(\phi(s)=s(1-s)/2\).
Use \(E_i+\ell_i/(2nN)\) for the last coordinate. Let
\(X=\sum_i X_i\). For a collision the reference at \(T=W\)
with multiplier \(X\) has exactly the two physical root values
under the original common scale. In the separate endpoint
scales their log ratio equals \(C_R+\sum_i E_i\). Therefore
\[
        \mathcal H(W,X)=\epsilon\delta/(2nN),\qquad
                  \mathcal B(X/W)=\delta/(2nN).                 \tag{4}
\]
The second equality uses uniqueness in (3).
These identities apply within the fixed neighborhood of the
specified root branches.

An aligned row of inverse midpoint \(u\) has weight \(u\).
Write \(K(u)=X_i/u\) for its normalized vector. This is analytic
at zero. For example, put \(f_{j,t}(u)=\log(1+(t/j-1/2)u)\)
for \(1\le t\le j\), and write \(k^\epsilon=k\) for
\(\epsilon=1\) and \(k^\epsilon=j-k\) for \(\epsilon=-1\).
The unnormalized coefficient coordinates are exactly
\[
 \sum_{t\le k^\epsilon} f_{j,t}(u)
    -\frac{k^\epsilon}{j}\sum_{t=1}^j f_{j,t}(u)
    +ju\phi(k/j).
\]
The center coordinate is
\[
 \epsilon\left(\frac1N\sum_{t=1}^N f_{N,t}(u)
              -\frac1n\sum_{t=1}^n f_{n,t}(u)
              +\frac{u}{2nN}\right).
\]
These follow from the finite factorial ramps for a unit gap.
Their constant and linear terms are zero, so \(K(0)=0\).
They also give an analytic continuation to small real \(u\).

For the other rows \(\sum_{d_i\ne\epsilon}w_i\le2\delta\), since
\(w_i-\epsilon\ell_i=w_i(1-\epsilon/d_i)\ge w_i/2\).
The [exact compensation estimates](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-ordinary-exclusion-on-the-aligned-cone)
give \(\|\sum_{d_i\ne\epsilon}X_i\|\le C_1\eta\delta\).
Consequently
\[
      v_0=\int K\,d\nu,\qquad
           \|v_0-X/W\|\le C_2\eta\delta .
\]
Analyticity of (3) gives a bounded Lipschitz constant in the
neighborhood. For small enough \(\eta\), (4) implies
\[
                       \mathcal B(v_0)\ge\delta/(4nN).           \tag{5}
\]
For tails with fixed signed gap and \(W\to T_*>0\), one has
\(\delta>0\) for all sufficiently small \(\eta\).
In fact \(\delta=0\) means every row is aligned, so the number
of unequal rows is the absolute signed gap. Such a bounded
number of rows has \(W\le\eta\,|\sum_i d_i|\to0\).

## A nonzero derivative of bounded order

Put \(D_0=D\mathcal B(0)\), a real linear form, and \(g(u)=D_0K(u)\).
We prove that \(g^{(s)}(0)\ne0\) for some \(2\le s\le n+N-1\).
Multiply \(g(u)\) by \(u\). The formulas above express the result
as a linear combination of the logarithms \(f_{j,t}\) plus a term
linear in \(u\). There are \(n+N-2\) distinct nonzero logarithmic
parameters \(t/j-1/2\). Coprimality of \(n,N\) gives a single
duplicate, at \(t/j=1\), and exactly one even order contributes
the zero parameter. If all the specified derivatives of \(g\)
vanish, expansion of the logarithms in powers \(u^3\) through
\(u^{n+N}\) gives a square Vandermonde system. Each coefficient
of a distinct nonzero logarithm is zero, including the combined
coefficient of the two endpoint logarithms.

Write \(C_E\) for the center entry of \(D_0\) and \(c_{j,k}\)
for its coefficient entries indexed after the reflection
\(k\mapsto k^\epsilon\). The coefficient \(\gamma_{j,t}\) of
\(f_{j,t}\) before combining orders is
\[
 \gamma_{j,t}
   =\sum_{1\le k<j}c_{j,k}\bigl(1_{t\le k}-k/j\bigr)
                   +\epsilon\,C_E\,\sigma_j/j,\qquad
               \sigma_N=1,\ \sigma_n=-1.
\]
Its sum over \(t\) is \(\epsilon C_E\sigma_j\).
In the odd order only the endpoint coefficient can remain.
The vanishing combined endpoint coefficient then fixes the
endpoint coefficient in the even order. The sum in that order
forces the coefficient at its zero logarithmic parameter to
vanish as well. Hence \(\gamma_{j,t}=0\) for both orders whenever
\(t<j\). Taking consecutive differences gives \(c_{j,k}=0\)
for \(1\le k\le j-2\).

This is impossible for \(j=N\ge4\). To see it, differentiate the
root equation at \(x=0,T=T_*\). The derivative of a log root
with respect to its \(k\)-th log coefficient is
\(-a_{j,k}z_j^k/(z_j F_j'(T_*,z_j))\), where \(a_{j,k}\) is
the reference coefficient. Eliminating the time by the angular
equation applies the same nonzero real projection
\(\epsilon(\Re w-(\Re p/\Im p)\Im w)\) to every entry, with
\(p=\partial_T\mathcal H(T_*,0)\), and multiplies by \(T_*\).
Thus the entries in an order are nonzero real coefficients
times \(\Re(C_j z_j^k)\) with \(C_j\ne0\) independent of \(k\).
Two consecutive entries cannot vanish, since their complex
arguments would then have a real ratio \(z_j\). The conclusion
is unchanged by the reflection of the index. This contradiction
proves the finite derivative assertion.

At a finite ordinary core the reference clock and the roots
in the normalized variables are algebraic after using
\(q_*=\exp(-T_* /(2nN))\), by the
[algebraic contact construction](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
In \(g^{(s)}(0)/T_*\) the rational logarithmic parameters and
the projected root derivatives are algebraic. Hence its sign
and the first suitable \(s\) can be decided for the specified
algebraic contact by exact real algebraic arithmetic.
No assumption of algebraicity is needed for existence of \(s\).

## Positive quadrature and the response comparison

Choose such an \(s\), and let \(\sigma=\operatorname{sgn}g^{(s)}(0)\).
Continuity gives a strictly uniform sign with magnitude bounded
away from zero for
\(\partial_u^s(D\mathcal B(v)[K(u)])\) over all small \(v\) and
\(0\le u\le\eta\). Choose \(a,b\in\{0,1\}\) and \(r\ge0\)
so \(s=2r+a+b\) and \(\sigma(-1)^b=-1\).
For even \(s\) this means no endpoints when \(\sigma=-1\),
and both endpoints when \(\sigma=1\). For odd \(s\) use
the endpoint zero when \(\sigma=-1\), and the upper endpoint
when \(\sigma=1\).

For completeness, take \(q\) to be the monic degree-\(r\)
orthogonal polynomial for the measure \(u^a(\eta-u)^b\,d\nu\).
If this weighted measure has at most \(r\) support points,
\(\nu\) has at most \(r+a+b\) nodes. Otherwise \(q\) exists
by the positive moment Gram matrix and has \(r\) simple
roots in \((0,\eta)\). The usual sign-change argument
proves the latter assertion by testing orthogonality against
the polynomial on its distinct odd-multiplicity zeros
in \((0,\eta)\) if their number were less than \(r\).
For \(r=0\) use \(q=1\).
Use these roots and the selected endpoints as quadrature
nodes. The weights are the integrals of the Lagrange basis
polynomials on the nodes. Polynomial division by
\(u^a(u-\eta)^b q(u)\) and orthogonality prove exactness
through degree \(s-1\). All weights are nonnegative.
For an interior node use the nonnegative polynomial
\(u^a(\eta-u)^b\) times the square of the product of the
other interior factors. Its values at all other nodes
are zero. For a chosen endpoint use \(q(u)^2\) times the
factor for the opposite endpoint when that endpoint is
selected. Each test polynomial has degree at most \(s-1\).
The weights sum to one. There are \(r+a+b\le N\) nodes.

For any \(C^s\) function whose \(s\)-th derivative has sign
\(\sigma\) and magnitude at least \(c_0>0\), interpolate at
the interior nodes with multiplicity two and at the selected
endpoints with multiplicity one. The interpolant has degree
at most \(s-1\). Rolle's theorem gives the Hermite error
as the \(s\)-th derivative at an intermediate point times
\(u^a(u-\eta)^b q(u)^2/s!\). This is at most
\(-c_0u^a(\eta-u)^bq(u)^2/s!\) on the interval by the choice
of \(b\), with the nodal values interpreted by continuity.
Integrating and using quadrature exactness bounds the
improvement of the integral below by
\(c_0\int u^a(\eta-u)^b q(u)^2\,d\nu/s!\).

Now join \(v_0=\int K\,d\nu\) to \(v_1=\int K\,d\nu^\#\) by
a segment \(v_t\). At each \(t\), apply the integral improvement
to \(u\mapsto D\mathcal B(v_t)[K(u)]\). The segment stays small
since \(K(0)=0\). Its derivative of the response is exactly
the improved integral. Integrate over \(t\) and use (5) to
obtain (1).

## An exact finite jet bound for a finite ordinary core

Write \(t=T-T_*\) and consider the angular equation with a free
compensating vector \(x\). Its solution \(t=\tau(x)\) is analytic
with \(\tau(0)=0\). Write \(b(x)\) for the signed radial response
at this solution, so \(b(0)=0\). For a finite ordinary core
both germs have algebraic real Taylor coefficients. Indeed,
after subtracting the root ratio at the reference contact,
\(\mathcal H\) uses a logarithm with value zero at the origin.
The two root equations use the algebraic reference coefficients
multiplied by exponentials of linear forms in \(t,x\). Their
simple roots have algebraic Taylor coefficients by recursive
implicit differentiation. The same holds for the subsequent
angular implicit solution because its linear denominator
\(\partial_T\Im\mathcal H(T_*,0)\) is nonzero and algebraic.
The reference center cancels before these differentiations.

Treat the positive strength as a formal independent parameter
\(A\) and solve
\[
             x(A,u)=(A+\tau(x(A,u)))K(u),\qquad
                 B_A(u)=b(x(A,u)).
                                                               \tag{6}
\]
Here \(x(A,0)=0\). Since \(K(0)=0\), recursively solving the
coefficient of each power of \(u\) shows that all coefficients
of \(x(A,u)\) and \(B_A(u)\) are polynomials in \(A\) with
algebraic coefficients. At \(A=T_*\), the series is exactly
that of \(B(u)\) by the angular uniqueness in (3).
At \(A=0\) the solution of (6) is identically zero. Differentiate
(6) formally with respect to \(A\) at zero to obtain
\[
   [A]B_A(u)=\frac{Db(0)[K(u)]}{1-D\tau(0)[K(u)]}.               \tag{7}
\]
The denominator in (7) has constant term one. The finite
derivative argument above applies to \(Db(0)K(u)\), since
\(D\mathcal B(0)=T_*Db(0)\). Thus its first nonzero term
occurs at some order \(k_0\le n+N-1\). The coefficient
of \(u^{k_0}\) in (6) is a nonzero polynomial in \(A\)
because the linear term detected in (7) is nonzero.

The algebraic contact construction gives \(q_*\in(0,1)\)
algebraic. The [Hermite–Lindemann theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)
implies that
\(T_*=-2nN\log q_*\) is transcendental. A nonzero polynomial
with algebraic coefficients cannot vanish on evaluation
at \(A=T_*\). Hence \(B(u)\) has a nonzero coefficient by
order \(k_0\). Its sign as \(u\downarrow0\) is computable.
One computes the finitely many algebraic coefficient
polynomials in (6) and takes the first one that is nonzero.
Interval evaluation at \(T_*\) determines its sign because
the value is nonzero.

This argument works as well for any fixed algebraic
probability weights and algebraic nonnegative relative
heights \(a_i\), at least one positive with positive weight.
Replace \(K(u)\) by \(\sum_i w_i K(a_i u)\) in (6).
The leading coefficient of the numerator of (7) is then
multiplied by the positive moment \(\sum_i w_i a_i^{k_0}\).
This gives the same finite jet bound for such fixed
shapes among the bounded-node configurations.

## One height and concentration at the critical boundary

Use the third-order boundary notation \(A_2,A_3,D,t_1\)
from the linked test. For the unscaled angular solution
\(\Im\mathcal H(T,x)=0\), the response has linear coefficients
\(A_2,A_3\) along the coefficients of \(u,u^2\) in \(K(u)\).
Its quadratic coefficient in the first of these directions
is \(D\), and the linear change in the solved time in the
same direction is \(t_1\). Inserting \(x=T(K(u))K(u)\)
in that unscaled solution gives
\[
 B(u)=T_* A_2 u+
       \left(T_* A_3+T_*^2D+T_*t_1 A_2\right)u^2+O(u^3).
                                                               \tag{8}
\]
Indeed the first derivative of the scaled time at zero along
\(K(u)\) is \(T_*t_1\). This accounts for the last term in (8).
In particular \(g''(0)=2T_* A_3\), so \(A_3<0\) chooses the
one-node Gaussian quadrature when \(s=2\). Its moment match
places the node at \(\mu\), and \(q(u)=u-\mu\), proving (2).
For a finite ordinary core \(A_2,A_3,D\) are algebraic.
If \(A_3<0\), either \(A_2\ne0\) or the quadratic coefficient
in (8) is nonzero by transcendence of \(T_*\).
The critical equality \(A_2=0,\ A_3+T_*D=0\) with \(A_3<0\)
can still be tested in a compactified reference with a
countable core. Its scalar response is \(O(u^3)\), with all
further terms retained in (2). If the scalar is eventually
nonpositive, (5) and \(\delta>0\) exclude collisions.

For a positive leading term of order \(k\), the bound in (2)
forces concentration near the smaller exponent \(\mu^{-1}-1/2\).
For example, at least half the curvature mass is within
\(O(\mu^{k/2-2})\) in the exponent when \(k>2\).
This follows by Chebyshev's inequality and the transformation
\(u\mapsto u^{-1}-1/2\) on \(u\asymp\mu\).
In particular \(k=3\) gives a window of width \(O(\mu^{-1/2})\).
If \(k\ge5\) the window contains at most one allowable height.
Every other aligned inverse height lies at least
\(c_1\mu^2\) from \(\mu\), since consecutive allowed values
\((L+1/2)^{-1}\) near \(\mu\) are spaced by that order.
Thus all but \(O(\mu^{k-4})\) of the curvature mass comes from
one exact aligned height. The other orientations and gaps
have total curvature \(O(\mu^k)\).

The theorem gives an exact reduction of the saturated contact
problem to bounded-node responses and a quantitative loss.
A positive value for the reduced response is a necessary
condition for a collision. Exclusion of all such positive
responses and the passage from concentration to integer
coprimality remain open.
