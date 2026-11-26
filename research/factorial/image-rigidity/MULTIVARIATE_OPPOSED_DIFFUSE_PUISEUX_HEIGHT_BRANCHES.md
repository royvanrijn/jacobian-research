# Local height branches at opposed diffuse factorial contacts

Fix any finite ordinary two-monomial core \(C\) and adjacent orders
\(n,N=n+1\), with \(n\ge2\). Add \(r_+\) copies of the increasing row
\((L,L+d)\) and \(r_-\) copies of the decreasing row \((Y+e,Y)\),
where \(d,e\) are fixed positive integers and the nonnegative integer
heights \(L,Y\) vary independently. Fix the signed gap
\(h=d r_+-e r_-\). Use the physical Gaussian reference of the
[diffuse exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts).
Let \((T_*,\lambda_*)\) be a positive reference contact with
\(-1/e\le\lambda_*/T_*\le1/d\).

**Theorem.** In a sufficiently small neighborhood of the specified
contact, every ordinary collision with large heights lies on one of
finitely many parametrized branches. After exchanging the two
columns if necessary, let \(L\) be a height whose limiting curvature
share is positive. On each branch the count parameter \(r\) and the
height \(Y\) have the forms
\[
 r=L f(L^{-1/M}),\qquad Y=L g(L^{-1/M}),
\]
with \(M\ge1\) an integer, \(f\) real analytic at zero and
\(g\) meromorphic at zero. Here \(r\) is an integer and determines
the two multiplicities by the fixed signed-gap relation. The
other height and both counts are unique on each branch for a
given large \(L\). This holds at multiple and tangential
contacts and when the other limiting curvature share vanishes.

The number of exact collision profiles with both heights at most
\(H\) near any such positive contact is \(O(1+\log\log(3+H))\).
Any prescribed nonzero rational polynomial relation between
\(L\) and \(r\) excludes all sufficiently high collisions there.
In the interior, the same exclusion under a relation in \(L,Y\)
holds at irrational bias. It also holds on any branch with a
transcendental coefficient in the Puiseux expansion of \(Y/L\).

For fixed core and gaps, the
[absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and the [signed-gap bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
allow the conclusions to be summed over the finitely many relevant
pairs, twists, and positive contacts. In particular, fix \(s_0>0\).
Among profiles with both heights sufficiently large and added
curvature at least \(s_0\), only \(O(1+\log\log(3+H))\) profiles
of heights at most \(H\) can have a collision in any degree.
Infinitely many sufficiently sparse profiles remain possible
under the counting estimate alone.

## The exact complex one-parameter germ

Put \(g_0=\gcd(d,e)\), \(c_+=e/g_0\), and \(c_-=d/g_0\).
If \(g_0\nmid h\), there are no ordinary tuples. Otherwise fix integers
\(s_+,s_-\) with \(ds_+-es_-=h\), and write
\[
 r_\pm=c_\pm r+s_\pm,\qquad
 u=1/r,\qquad p=r/L,\qquad q=r/Y .
\]
In the regime in the theorem, \(r\to+\infty\). The limiting shares are \(c_+d^2p_*\) and \(c_-e^2q_*\). Using the smaller heights in
those shares gives the same limit as using the midpoints or the larger
heights.

At order \(j\), strip the smaller-end tail factorials and divide
by \((jY)^{j e r_-}\). Make the same variable change
\(z\mapsto L^{-d r_+}Y^{e r_-}z\) in both orders. The remaining
multiplier of the core coefficient at index \(k\) is exactly
\[
 j^{hk}\exp\left[
  (c_+/u+s_+)\sum_{a=1}^{kd}\log(1+a u p/j)
 +(c_-/u+s_-)\sum_{a=1}^{(j-k)e}\log(1+a u q/j)
 \right].
 \tag{1}
\]
This formula defines a holomorphic family in \((u,p,q)\) near
\((0,p_*,q_*)\) by using the logarithms at one. The apparent poles at
\(u=0\) are removable.

At \(u=0\), set
\[
 T=c_+d^2p+c_-e^2q,\qquad \lambda=c_+dp-c_-eq .
\]
The log multiplier in (1), apart from \(hk\log j\) and a scalar
independent of \(k\), is the sum of \(-Tk(j-k)/(2j)\) and a linear
term \(k[c_+d^2p-c_-e^2q+\lambda/j]/2\). This follows from the
sum of the first \(kd\) or \((j-k)e\) positive integers. Thus a
common argument normalization by the order-\(n\) endpoint scale
reduces the \(u=0\) pair to
\[
 H_n(w,z),\qquad H_N(w,(N/n)^h v z),\qquad
 w=e^{-T/(2nN)},\quad v=e^{-\lambda/(2nN)},
 \tag{2}
\]
where \(H_j(w,z)=\sum_k [z^k]Q_j^C\,w^{(nN/j)k(j-k)}z^k\).
Both scalar normalizations and the common argument normalization
extend holomorphically and remain nonzero near the point. The change
from \((p,q)\) to \((w,v)\) is locally invertible since the two linear
forms \(T,\lambda\) are independent and \(wv\ne0\).

Let \(A\) be the resultant of the two off-diagonal root-quotient
polynomials of the normalized family, and let \(D\) be its physical
resultant. The quotient polynomials can be taken monic and
holomorphic, since the endpoints are nonzero. They are invariant
under any nonzero variable dilation. If the common factor has degree
at least two, either two distinct common roots give a shared
quotient, or one repeated common root gives the shared off-diagonal
quotient one. An ordinary collision has degree at least two because
the even moment is positive on the real axis and a nonreal common
root comes with its conjugate. Thus \(A=D=0\) contains every ordinary
collision locally.

At \(u=0\), \(A\) is a function of \(w\) only and is a nonzero
polynomial. At \(w=0\), the reference polynomials are endpoint
binomials. Their off-diagonal quotients are precisely nontrivial
roots of unity of their respective orders, and \(n,N\) are coprime.
For fixed \(w=w_*\), \(D(0,w_*,v)\) is a nonzero polynomial in \(v\):
there are only finitely many ratios of one nonzero root from each
parent. These properties give an isolated central complex fiber
without assumptions on the roots or on the order of vanishing of
either resultant.

## Finitely many Puiseux branches

Weierstrass preparation applied to \(A\) in \(w-w_*\) gives a
monic polynomial with holomorphic coefficients in \((u,v)\).
Let its degree be \(a\). Keep its roots in a fixed small circle
about \(w_*\) and form the product of \(D(u,w,v)\) over these
roots, counted with multiplicity. This product is holomorphic
even when the roots coalesce. For instance the power sums of
the \(D\)-values are contour integrals against the logarithmic
derivative of the prepared polynomial, and Newton identities give
their product. At \(u=0\) the product equals \(D(0,w_*,v)^a\),
which is not identically zero. Weierstrass preparation in \(v-v_*\)
then gives a monic polynomial over the holomorphic germs in \(u\).

The roots of a monic polynomial over holomorphic germs in one
variable admit convergent Puiseux expansions. In this situation
one can see it by taking the squarefree part over meromorphic
germs. Its roots remain bounded, so its coefficients extend
holomorphically through zero. On a small punctured disc the
discriminant of the squarefree part is nonzero. Analytic
continuation permutes the simple roots, and a finite power
substitution makes each branch single-valued. Removable
singularities give convergence at the origin. Apply this first
to the prepared polynomial in \(v\), then to the one in \(w\)
after substituting each branch. Hence \(p,q\) for every possible
collision lie among finitely many holomorphic graphs over
\(u^{1/M}\), with \(M\) allowed to depend on the graph. A nonreal
graph cannot give accumulating real points unless its coordinates
are real analytic on the chosen positive parameter interval, by
the identity theorem for the imaginary parts.

Writing \(s=u^{1/M}\) on a real branch, choose \(p_*>0\)
by exchanging the columns if needed. The parameter
\(t=L^{-1/M}=s\,p(s)^{1/M}\) has nonzero derivative at zero and
can be inverted analytically. The required functions are
\(f(t)=p(s(t))\) and \(g(t)=p(s(t))/q(s(t))\).
A branch with \(q\) identically zero has no points with both
ordinary heights finite. Otherwise \(g\) is analytic when
\(q_*>0\) and has a pole of positive integer order when
\(q_*=0\). Recovering \(u\) from \(L\) proves uniqueness of the
two counts and the other height on the given branch.

In the boundary case \(q_*=0\), write \(g(t)\sim c t^{-\ell}\)
with \(\ell\ge1\) and \(c>0\). Then
\(Y\sim c L^{1+\ell/M}\). The
[opposed mark-product bound](MULTIVARIATE_OPPOSED_MARK_PRODUCT_HEIGHT_CONSTRAINT.md#consequences-for-arbitrarily-separated-scales)
at either marked order \(m\in\{n,N\}\), marking the row of
gap \(e\) after reversal, requires
\[
 \frac{1}{1+\ell/M}\ \ge\ \frac{2e-\lfloor2e/m\rfloor}{me}
\]
for any accumulating sequence of exact ordinary collisions.
A branch failing this inequality is excluded.

## Integral heights and counts

In (2), \(w_*\) and \(v_*\) are positive algebraic numbers with
\(0<w_*<1\). This follows from the nonzero angular eliminant and
the finite set of physical ratios at its algebraic \(w_*\).
The active count slope is
\[
 f(0)=p_*=\frac{T_*+e\lambda_*}{c_+d(d+e)}
        =-\frac{2nN}{c_+d(d+e)}\log(w_*v_*^e).
\]
It is positive. The argument of its logarithm is therefore
different from one, so \(f(0)\) is transcendental by
Hermite–Lindemann. A rational polynomial in \((L,r)\) which
vanished identically on a branch would have its leading
homogeneous part vanish at \((1,f(0))\). Thus any prescribed
polynomial relation in \((L,r)\) allows only finitely many
points on the germ.

For a height-only relation in the interior, use
\[
 g(0)=p_*/q_*=\frac{T_*+e\lambda_*}{T_*-d\lambda_*}.
\]
The slope \(\lambda_*/T_*=\log v_*/\log w_*\) is rational or
transcendental by the [Gelfond–Schneider theorem](https://www.mathnet.ru/eng/im4924).
The rational fractional transformation for \(g(0)\) is invertible,
so \(g(0)\) is transcendental when the bias slope is irrational.
More generally an algebraic branch over \(\overline{\mathbb Q}(1/L)\)
has algebraic Puiseux coefficients. Any transcendental coefficient
of \(Y/L\) therefore excludes an identically vanishing rational
relation between the two heights.

For a uniform bound on the points without a prescribed relation,
choose an integer \(D\) large enough that
\(b=(D+1)(D+2)/2\) satisfies \((b-1)/(2MD)>1\).
Use the \(b\) monomials \(L^a r^c\) of total degree at most \(D\).
At \(b\) distinct integral points in \(X\le L\le X^\gamma\),
multiply each row of the monomial evaluation matrix by \(L^{-D}\).
Its columns are now evaluations at \(t=L^{-1/M}\) of the fixed
analytic functions \(t^{M(D-a-c)}f(t)^c\). Newton divided
differences bound this determinant by a constant times
\(\prod_{i<j}|t_j-t_i|\). This is the elementary determinant
estimate of [Bombieri and Pila, Proposition 2](https://people.maths.ox.ac.uk/pila/Ovals.pdf).
Thus the original integer determinant is at most
\[
       O\big(X^{\gamma Db-b(b-1)/(2M)}\big)
\]
in absolute value. Fix \(1<\gamma<(b-1)/(2MD)\). At large \(X\),
every such determinant is zero, so if the block has at least
\(b\) points a single nonzero rational polynomial of degree at
most \(D\) vanishes on the entire block.

There is a uniform bound for the zeros near \(t=0\) of every
nonzero germ in the finite-dimensional span of the analytic
column functions. Otherwise, normalize its coordinates in a
basis to norm one and take a convergent subsequence with
arbitrarily many zeros tending to zero. The limit is a nonzero
analytic function. A small circle about zero on which it is
nonzero and uniform convergence on this circle contradict the
unbounded zero counts by Rouché's theorem. Since no nonzero
rational polynomial of degree at most \(D\) vanishes
identically on the branch, this bounds the number of integer
points in every block \([X,X^\gamma]\). Iterating
\(X\mapsto X^\gamma\) covers the heights up to \(H\) with
\(O(1+\log\log(3+H))\) blocks.

For the uniform statement above \(s_0>0\), use the
[absolute upper curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
on any collision. The signed-gap bound and the absolute terminal
cutoff leave finitely many cases. When both heights tend to infinity
with added curvature in \([s_0,12]\), the coordinate \(u\) tends to
zero and \(p,q\) lie in a bounded set with at least one bounded away
from zero. If escaping collisions were to avoid the finitely many
local neighborhoods, a subsequence would give a contact of the
central family (2) outside them, by continuity of the resultants
and compactness. The finitely many local counts therefore cover all
sufficiently high profiles at the stated curvature lower bound.

This counting scale is consistent with the existing
[Wilkie theorem on integer points of globally subanalytic curves](https://doi.org/10.2178/jsl/1096901771).
The proof here uses the finite local analytic branches of the
exact moment equations and the transcendence of their height
data. A positive contact gives a transcendental count slope
at either active height. At rational bias the height-only ratio
may still have a transcendental higher coefficient, as follows.

## A stationary contact of rational bias

For simple contacting roots let \(z_n(T),z_N(T)\) be their local
Gaussian branches before the upper physical rescaling. Define
\(\lambda_0(T)\) by the modulus comparison at those roots, so
\[
 \lambda_0(T)=2nN\left(h\log(N/n)-\Re\log(z_N(T)/z_n(T))\right).
\]
Here the raw core polynomials are used in the Gaussian deformation.
Suppose \(\sigma=\lambda_*/T_*\) is rational, the angular difference
has a zero of order \(k\ge2\), and \(\lambda_0'(T_*)\ne\sigma\).
Assume also that the coefficient of \(u\) in the angular equation
at the fixed reference \(T=T_*,\lambda=\lambda_*\) is nonzero,
with the midpoint curvature coordinates defined below. Then
every real outgoing branch continuing this specified root pair
has a transcendental height coefficient. The height-only algebraic-curve exclusion applies to those branches.

To verify the criterion, use the exact midpoint shares
\(x=r_+d^2/(L+d/2)\), \(y=r_-e^2/(Y+e/2)\), with
\(T=x+y\) and \(\lambda=x/d-y/e\). These are analytic local
coordinates together with \(u\). In these linear curvature
coordinates the [compensated factorial expansion](MULTIVARIATE_OPPOSED_INDEPENDENT_HEIGHT_RATIONAL_BIAS.md#analytic-equations-in-the-independent-curvatures)
is uniform in neighborhoods of the point. Its first non-Gaussian
coefficient vector is quadratic in \(x,y\) with rational
coefficients. In fact each per-row compensated coefficient
and center is \(a_i/M_i^2+O(M_i^{-3})\) with \(a_i\) rational
and \(M_i\) the midpoint. The factor
\(r_i/M_i^2\) is \(u\) times a rational multiple of the square
of its share, up to \(O(u^2)\). Multipliers at the reference
are algebraic by (2); thus simple-root response derivatives
and \(\lambda_0'(T_*)\) are algebraic. Since both shares at
rational bias are rational multiples of \(T_*\), the first
angular forcing is \(u\,a_0 T_*^2\) with \(a_0\) algebraic.

The modulus equation solves analytically for
\(\lambda=\lambda_0(T)+O(u)\). The angular equation after
this substitution has first terms
\(a_k(T-T_*)^k+a_0T_*^2u\), with \(a_k,a_0\) nonzero algebraic,
and remainder \(O((T-T_*)^{k+1},u(T-T_*),u^2)\) in the real
variables. At \(u=0\) the angle is independent of \(\lambda\),
so the modulus substitution leaves the coefficient \(a_0\)
unchanged. Each Puiseux solution therefore has a first change
in \(T\) of order \(u^{1/k}\) with coefficient an algebraic
nonzero multiple of \(T_*^{2/k}\). The ratio \(Y/L\) differs
from \((T+e\lambda)/(T-d\lambda)\) by \(O(u)\).
Its first change along the branch comes from the derivative
of that fraction with \(\lambda=\lambda_0(T)\), which is
\[
 \frac{(d+e)(\lambda_0'(T_*)-\sigma)}
              {T_*(1-d\sigma)^2}.
\]
Finally \(1/L=u(x_*/(c_+d^2)+o(1))\) with \(x_*/T_*\)
positive rational. In powers of \(1/L\), the first change
of \(Y/L\) is therefore at exponent \(1/k\), with coefficient
a nonzero algebraic multiple of \(T_*^{1/k-1}\).
The strength \(T_*\) is transcendental by Hermite–Lindemann,
since \(w_*\in(0,1)\) is algebraic. Hence the coefficient is
transcendental. This criterion covers the stated stationary
contacts, while the transverse rational-bias contacts are
excluded by the linked rational-bias theorem.

## Two groups of affine row types

The local branch theorem and its transcendental-slope conclusions
also hold for two groups of fixed affine height types. Within
group \(b=1,2\) allow any fixed list of smaller exponents
\(\ell_i=a_i+\alpha_i L_b\) with rational \(a_i\) and positive
rational \(\alpha_i\), fixed positive integer gaps \(d_i\), signs
\(\epsilon_i\), and multiplicities \(c_i r+s_i\), where
\(c_i\ge0\) and \(c_i,s_i\) are integers. The parameters \(L_b,r\)
are integers on the ordinary family, with all endpoints and counts
admissible. Require \(\sum_i\epsilon_i c_i d_i=0\), and set
\[
 A_b=\sum_{i\in b} c_i d_i^2/\alpha_i,\qquad
 B_b=\sum_{i\in b} \epsilon_i c_i d_i/\alpha_i .
\]
Assume \(A_1A_2>0\) and \(A_1B_2-A_2B_1\ne0\). At a
positive contact at least one group contributes positive
limiting curvature. Indeed, factor
out each \(j\alpha_i L_b\) in (1) and replace the logarithms
there by \(\log(1+(a_i+a/j)u p_b/\alpha_i)\),
where \(p_b=r/L_b\). Also use the common argument factor
\(\prod_i(\alpha_i L_b)^{-\epsilon_i d_i(c_i r+s_i)}\).
The offset \(a_i\) contributes only constants and linear
terms in the coefficient index at \(u=0\). Its linear term
does not depend on \(j\), so it cancels from the adjacent
endpoint scale comparison. The limiting clocks are therefore
still (2), with \(T=A_1p_1+A_2p_2\),
\(\lambda=B_1p_1+B_2p_2\), and \(h=\sum_i\epsilon_i s_i d_i\).
The determinant assumption gives the same invertible clock
chart. A nonzero \(p_{b,*}\) is a rational linear combination
of \(T_*,\lambda_*\). Exponentiating a nonzero integer
multiple therefore gives an algebraic number, so it is
transcendental. The same bound follows by counting \((L_b,r)\).
When both shares are positive, \(p_{1,*}/p_{2,*}\) is an
invertible rational fractional transformation of the bias slope,
which also gives the height-only exclusion at irrational bias.
