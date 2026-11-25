# A fixed-order coefficient-error obstruction with exact endpoint attenuation

At a sufficiently high consecutive pair of orders, independent error bars
of any fixed algebraic order in the outer deficit can admit a common root.
This remains true when their reference parameters come from the same ordinary
profile at both orders and the perturbed coefficient arrays satisfy the
[exact consecutive endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality).
The counterexample applies to a relaxation of the coefficients.

## Statement

Fix a real \(q>0\) and an even lower order \(n=2m\) with
\(m-1>\max(q,3)\). Put \(N=n+1\). Use the four ordinary rows
\[
 (\ell,\ell+1),\quad (\ell+1,\ell),\quad
 (2\ell,2\ell+1),\quad (2\ell+1,2\ell),
 \qquad \ell\in\mathbb Z_{>0}.
\]
For \(j=n,N\), let \(F_{j,k}\) be the factorial coefficient product
normalized by the endpoint chord, and put
\[
 P_j(z)=\sum_{k=0}^j\binom jk F_{j,k}z^k .
\]
In particular \(F_{j,0}=F_{j,j}=1\).
Write \(\delta=1-\sqrt{F_{n,1}F_{n,n-1}}\).

For all sufficiently large \(\ell\) there is a polynomial \(\widetilde
P_n\) such that \(\widetilde P_n,P_N\) have a common nonreal unit-circle
root. It differs from \(P_n\) only at \(k=m-1,m,m+1\), by
\(O(\delta^{m-1})=o(\delta^q)\) in both the coefficients and their
logarithms. The perturbation preserves reciprocity, positivity of the
coefficients, strict log convexity after division by \(\binom nk\), and
positivity on the negative real axis. The pair still obeys the exact
endpoint attenuation inequalities in both orientations with the actual
curvature \(S\) and the actual minimum \(M\) of the maximum exponents over
unequal rows. It also obeys the coefficient mean-defect caps using the
largest tanh parameter as in the
[quadratic model obstruction](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#a-failure-of-independent-errors-in-the-near-flat-angular-comparison).

Consequently a remainder model that allows two-sided, independent
changes of at least \(c\delta^q\) in these central coefficients, for some
fixed \(c>0\), cannot ensure angular separation at this pair of orders
using those additional constraints. The matched tanh moments and the
excess outer loss can all be held at their exact reference values in
such a model. This formulation assumes the moments parameterize the
coefficient enclosures while the stated central changes remain
admissible. An exact reconstruction of the coefficients from those
moments would impose a different constraint.

## The common root

For these rows \(\delta=\Theta(\ell^{-1})\), by the one-factor
formulas for the outer loss. The
[linked Hermite limit](MULTIVARIATE_FIXED_TERMINAL_HERMITE_SEPARATION.md#exact-cancellation-in-the-cayley-transform),
applied as in the quadratic model obstruction, gives an upper-order root
\(t=e^{\mathrm i\theta}\) with \(\pi-\theta=\Theta(\sqrt\delta)\) and
\(P_n(t)=O(\delta^m)\). Indeed a nonzero imaginary simple root of the
upper Hermite anchor gives a unique nearby root in the rescaled Cayley
coordinate. Reciprocal conjugation fixes that branch, putting its root
on the unit circle. The lower Cayley polynomial is \(O(\delta^m)\) at
this bounded rescaled argument.

Reciprocity makes \(G=(-t)^{-m}P_n(t)\) real. Add to \(P_n(z)\)
\[
 H(z)=
 \begin{cases}
   -G(-z)^m,&G\le0,\\[2mm]
   \displaystyle\frac{G}{2+2\cos\theta}(-z)^{m-1}(1+z)^2,&G>0.
 \end{cases}
\]
In the second case \((1+t)^2/(-t)=-(2+2\cos\theta)\).
In both cases \(H(t)/(-t)^m=-G\), so \(P_n+H\) vanishes at \(t\).
Both increments are reciprocal of order \(n\) and nonnegative for
negative real \(z\). The actual even moment \(P_n\) is positive there
because the two distinct monomials have a nonzero linear combination
with probability one under the independent exponential variables.
Since \(2+2\cos\theta=\Theta(\delta)\), the increment affects only the
three stated coefficients, by at most \(O(\delta^{m-1})\).
Each coefficient of \(P_n\) tends to the corresponding binomial
coefficient, so positivity persists and the same \(O(\delta^{m-1})\)
bound applies to changes of logarithms.

The second difference of a log Gauss multiplier in increasing
orientation is \(\log(1+x/(1+kx))\), which is
\(\Theta(\delta)\) for each of these four factors at the fixed orders.
Reversal just reverses the second differences. The perturbation
therefore preserves strict log convexity. It also preserves the
mean-defect caps from the quadratic model obstruction. In that proof
the defect at each perturbed index, divided by the outer loss, is a
weighted average of \(c_kT^2+O(T^4)\) at the two lower-order tanh
parameters, with \(c_k>0\). The parameter ratio tends to two and the
weights tend to \(2/3,1/3\) with the larger parameter first. The
defect is thus below the cap that uses the larger parameter by an
amount \(\Theta(\delta^3)\), since the outer loss is \(\Theta(\delta)\).
The perturbation is \(o(\delta^3)\). Since the endpoints and the outer
coefficients are unchanged, the outer losses at both orders and their
excess are unchanged. The tanh moments are those of the fixed ordinary
reference parameters.

## Slack in exact endpoint attenuation

For one row \((a,b)\) put \(d=|b-a|\), \(L=\max(a,b)\), and
\(f(s)=\psi'(1+s)\). For a real order \(j\ge k\), the derivative
with respect to \(j\) of its contribution to \(-\log F_{j,k}\) is
\[
 k\int_a^b [x f(jx)-a f(y)]\,dx,\qquad y=(j-k)a+kx.
 \tag{1}
\]
For example, the contribution itself is
\(k\int_a^b[\psi(1+jx)-\psi(1+y)]\,dx\), by the log-Gamma
endpoint chord formula. The function \(s f(s)\) increases, as seen
by writing
\[
 s f(s)=\int_0^\infty e^{-v}
             \frac{v/s}{e^{v/s}-1}\,dv\qquad(s>0)
\]
and using the monotonic decrease of \(w/(e^w-1)\). For \(x\ge a\)
this yields \(xf(jx)\ge y f(y)/j\), since \(jx\ge y\).
For \(x\le a\) this inequality reverses, as does the orientation of
the integral. Hence (1) is bounded below by
\[
 \frac{k^2}{j}\int_a^b(x-a)f(y)\,dx .
\]
The trapezoid bound for the convex function \(s\mapsto(z+s)^{-2}\)
on each unit interval gives, for \(z>0\),
\[
 \sum_{r=0}^\infty(z+r)^{-2}\ \ge\ \frac1z+\frac1{2z^2}.
\]
Indeed the integral from zero to infinity is \(1/z\) and is
bounded above by the sum with the first term halved. Use \(z=1+y\)
and \(y\le jL\). The integral of \(x-a\) with the displayed
orientation equals \(d^2/2\). We obtain the row bound
\[
 \frac{d}{dj}[-\log F_{j,k}]_{\text{row}}
 \ \ge\
 \frac{k^2d^2}{2j(1+jL)}
 +\frac{k^2d^2}{4j(1+jL)^2}.
 \tag{2}
\]
With \(S=\sum_i d_i^2/L_i\) and \(M=\min_i L_i\), the sum of
the first terms in (2) is at least \(k^2 S/[2j(j+1/M)]\).
Integrating over \([n,N]\) gives the logarithmic decrement in
the exact endpoint attenuation. Each of the four additional terms
has \(d=1\) and \(L\le2\ell+1\). Their integrated contribution is at
least
\[
 \frac{k^2}{N[1+N(2\ell+1)]^2}.
\]
For every fixed \(1\le k\le n\) this is a positive lower bound of
order \(\delta^2\) for the slack in the logarithmic inequality.
The bound is unchanged on reversing all rows. The central
logarithmic changes are \(o(\delta^2)\), and the upper coefficients
are exact. Thus every attenuation inequality for this adjacent pair
survives.

The common root belongs to the modified pair. The argument shows why
a fixed-order independent remainder bound, even with the exact
attenuation and linked moment parameters, needs further coefficient
coupling to yield an all-order angular theorem.
