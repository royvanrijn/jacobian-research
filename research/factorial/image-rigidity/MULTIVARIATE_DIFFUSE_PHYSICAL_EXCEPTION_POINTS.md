# Finite physical exceptions for diffuse extensions of any ordinary core

Fix any finite ordinary two-monomial profile \(C\), allowing an empty
core and private rows. Append any finite multiset \(E\) of unequal
shared rows. For each added row write \(M_i=\max(a_i,b_i)\) and
\(d_i=b_i-a_i\), and set
\[
 \tau=\sum_E\frac{d_i^2}{M_i},\qquad
 \lambda=\sum_E\frac{d_i}{M_i},\qquad
 \eta=\max_E\frac{|d_i|}{M_i},\qquad h=\sum_E d_i\in\mathbb Z .
\]
Here \(|\lambda|\le\tau\). A physical common root is a root at
the same ratio of the two monomial coefficients in both orders.

**Theorem.** There is an effectively described finite set
\(\mathcal F_C\subset\{(\tau,\lambda):0<\tau\le12,\ |\lambda|\le\tau\}\).
For every compact \(K\subset\{(\tau,\lambda):\tau>0,\ |\lambda|\le\tau\}\)
disjoint from \(\mathcal F_C\), all sufficiently diffuse added ordinary
profiles with \((\tau,\lambda)\in K\) satisfy SFC(2) in every degree.
There is no bound on the number of added rows or their signed total
gap. The tolerance is computable for a finite union of closed
rational rectangles disjoint from the exceptional set and for the
compact algebraic line segments specified below.

More quantitatively, take a positive rational \(\epsilon<12\) and
\(J\ge4\) at least the
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
There are computable positive constants \(B\) and \(\eta_0\) such
that a failure anywhere in an ordinary extension with
\(\epsilon\le\tau\le12\) and \(0<\eta\le\eta_0\) forces
\[
   \operatorname{dist}\big((\tau,\lambda),\mathcal F_C\big)
      \le B\eta^{1/((J-2)(J-1))^4}.                           \tag{1}
\]
If \(\mathcal F_C\) is empty, the assertion means that no such failure
is possible. Both estimates use the existing absolute terminal cutoff,
whose value is unevaluated.

At each exceptional point there are positive algebraic numbers \(q,v\)
with \(0<q<1\) and a terminal \(N=n+1<J\) such that
\[
               \tau=-2nN\log q,\qquad \lambda=-2nN\log v.      \tag{2}
\]
In particular, an affine line \(\lambda=a\tau+b\) with real algebraic
\(a,b\) can contain an exceptional point only if \(b=0\) and \(a\) is
rational. Thus for any real algebraic irrational \(a\in(-1,1)\) and
positive rationals \(\epsilon<L\), a fixed neighborhood of the compact
segment \(\lambda=a\tau,\ \epsilon\le\tau\le L\) is covered for all
sufficiently diffuse ordinary tails. A specified compact segment of
a line with \(b\ne0\) in the same cone has the same property.

## The finite set of physical reference contacts

Write \(Q_j(z)=\sum_k b_{j,k}z^k\) for the raw core factorial
moment in order \(j\), and \(G_j=(b_{j,j}/b_{j,0})^{1/j}\). Put
\(\Lambda=N/n\). For one adjacent pair define
\[
 H_j(q,z)=\sum_{k=0}^j b_{j,k}q^{(nN/j)k(j-k)}z^k,\qquad j=n,N .
\]
The positive diffuse reference pair at integer drift \(h\) is
\[
       H_n(q,z),\qquad H_N(q,\Lambda^h v z),\qquad
       q=e^{-\tau/(2nN)},\quad v=e^{-\lambda/(2nN)}.            \tag{3}
\]
The cone condition is \(q\le v\le q^{-1}\). Both normalizations
are relative to the tail alone, so the core coefficients in (3)
are rational.

For \(\tau>0\) the even reference has no real roots. Indeed its
Gaussian multiplier is the expectation of \(e^{k(X-\tau/2)}\)
with \(X\) centered normal of variance \(\tau/j\). The even core
is nonnegative on the real axis by its moment representation.
For real nonzero \(z\), the positive random scale misses its finite
root set almost surely. At zero its constant coefficient is positive.

A contact of (3) therefore forces a common off-diagonal root
quotient of \(H_n,H_N\), by dividing each root by its complex conjugate.
For a degree-\(j\) polynomial \(F\) with nonzero endpoints write
\[
 U_j(F;u)=
       \frac{\operatorname{Res}_z(F(z),F(uz))}{(u-1)^j}.
\]
This is polynomial in the formal coefficients and \(u\). The root
product for the numerator gives the \(j\) diagonal factors and
identifies the other factors with the \(j(j-1)\) ordered quotients
of distinct root indices. Exact division by the fixed monic
polynomial preserves polynomial dependence. Multiplicities are
retained. The resultant
\[
 A_n(q)=\operatorname{Res}_u(U_n(H_n;u),U_N(H_N;u))
\]
is a nonzero rational polynomial. At \(q=0\) the references are
binomials, and the quotients are the nontrivial roots of unity
of orders dividing their respective degrees. The orders \(n,N\)
are coprime. This is the
[angular eliminant for an arbitrary core](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)
without the auxiliary quotient for a real odd root, which is
unneeded for a physical contact.

Endpoint normalizing \(H_j\) gives coefficients at most \(\binom jk\)
and endpoints equal to one by log convexity and Gaussian attenuation.
Its roots consequently satisfy \(2^{-j}\le|G_j z|\le2^j\), by the
monic root bound and reversal. At a reference contact,
\[
          \big|\log(G_N/G_n)+ h\log\Lambda\big|
              \le (n+N)\log 2+|\lambda|/(2nN).                 \tag{4}
\]
The curvature bound \(\tau\le12\) gives finitely many \(h\) uniformly
in \(q,v\). For later use allow the larger window
\[
       |h|\le N\,[n+N+24+B_n],
       \qquad B_n\ \text{a rational upper bound for }\
                                  |\log(G_N/G_n)| .
 \tag{5}
\]
This uses \(\log\Lambda\ge1/N\).

For each \(3\le N<J\), isolate the algebraic roots of \(A_n(q)\)
in \(e^{-12/(2nN)}\le q<1\). For each allowed \(h\) form
\[
 D_{n,h}(q,v)
       =\operatorname{Res}_z(H_n(q,z),H_N(q,\Lambda^h v z)).
\]
At any fixed \(q\) this is a nonzero polynomial in \(v\), of
degree at most \(nN\). At \(v=0\) its root product is the \(n\)-th
power of \(b_{N,0}\) times the \(N\)-th power of the lower leading
coefficient, which are both nonzero. Isolate its real roots
\(q\le v\le q^{-1}\), and take the union of (2) over these finite
lists. Leading coefficients are nonzero in the allowed region,
so this is exactly the set of reference contacts for the indicated
terminals. The left \(q\)-endpoint can be compared effectively
by rational interval bounds. It is transcendental, by the
one-logarithm case of Baker's theorem, and thus cannot equal an
algebraic root. The absolute cutoff handles the actual ordinary
extensions above these terminals.

Both numbers in (2) are logarithms of positive algebraic numbers
times the same rational multiplier. If \(\lambda=a\tau+b\), then
\(\log v-a\log q=-b/(2nN)\). Baker's theorem on a nonzero algebraic
linear form in logarithms forces \(b=0\). When \(a\ne0\), Baker’s
linear-independence theorem shows that
\(\log v=a\log q\) with algebraic \(a\) forces \(a\) rational.
If \(\log v=0\), then \(a=0\) since \(q<1\).
[Baker, linear forms in logarithms II](https://doi.org/10.1112/S0025579300008068) and [III](https://doi.org/10.1112/S0025579300003843).

## Uniform approximation of the physical pair

Take \(\eta\le1/2\). For each order \(j\), divide the tail contribution
by its constant coefficient and use its own positive endpoint
argument scale. The normalized tail factor at index \(k\) is
\(e^{T_{j,k}}\) with zero logarithms at the endpoints. The curvature
estimate for the
[exact row chord](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
gives
\[
           |T_{j,k}+\tau k(j-k)/(2j)|
                 \le \eta\tau k(j-k)/j.                      \tag{6}
\]
In fact the row's log-Gamma second derivative on its interval
is between \(1/(jM_i+1)\) and \(1/[jM_i(1-\eta)]\).
These bounds differ from \(1/(jM_i)\) by at most \(2\eta/(jM_i)\),
since \(M_i^{-1}\le\eta\). The normalized error at \(k\) then
follows by comparison with the quadratic chord error.

In the original common-root equation rescale the variable once
by the lower tail endpoint scale. The upper argument scale is
the ratio of the two tail endpoint scales. Its logarithm is
\(h\log\Lambda+\sum_E E(w_i)\) in the notation of the exact row
chord. For a shared row the signed center estimate there gives
\[
 \left|E(w_i)+\frac{\log(b_i/a_i)}{2nN}\right|
     \le\frac{|d_i|}{4n^3\min(a_i,b_i)^2}.
\]
The right side sums to at most \(\eta\tau/n^3\).
Also
\[
 \left|\sum_E\log(b_i/a_i)-\lambda\right|\le\eta\tau .
\]
Indeed the single-row magnitude of the difference is bounded by
\((|d_i|/M_i)^2/[2(1-\eta)]\), using the logarithm series, and
\(\sum_E |d_i|/M_i\le\tau\). Hence the actual exponent of the
upper relative scale differs from \(h\log\Lambda-\lambda/(2nN)\) by
at most \(\eta\tau(1/n^3+1/(2nN))\).

At an actual common root the same endpoint-normalized root annuli
apply by factorial log convexity. The preceding center bound shows
that (5) holds for that root, since the tail correction has magnitude
at most \(2\tau\). Thus throughout the potential common-root
calculation there are finitely many \(h\) and \(n\), and all scales
in (3) are bounded above and away from zero. Equations (6) and
the center estimate approximate the actual pair by (3) in a fixed
coefficient box, with effective error \(O_{C,J}(\eta\tau)\).
Only the tail is approximated. No coefficient of the fixed core
is replaced.

## From coefficient error to distance to the finite set

Use \(\epsilon\le\tau\le12\), take the two resultants of the actual
pair, and use the fact that they vanish at a common physical root.
For the angular one use the actual tail-normalized arrays before
the upper relative scale, since the quotients are scale invariant.
The even ordinary moment excludes real common roots, so the quotient
resultant applies. Both resultants are polynomial functions of
the coefficients being compared. They are Lipschitz on the bounded
coefficient boxes, by bounding their first derivatives there.
At an actual collision we therefore have effective estimates
\[
               |A_n(q)|\le C_1\eta,\qquad
               |D_{n,h}(q,v)|\le C_2\eta.                     \tag{7}
\]

For a nonconstant polynomial of degree \(D\), the distance
to its nearest complex root is at most
\((|P(x)|/|\operatorname{lc}P|)^{1/D}\), by its root product.
A nonzero constant polynomial instead excludes sufficiently small
values. Apply this first to \(A_n\). The degree in \(q\) of \(H_j\)
is at most \(nN j/4\), so the first resultant in \(z\) has \(q\)-degree
at most \(nN j^2/2\). The respective \(u\)-degrees are at most
\(j^2\); the second resultant gives \(\deg_q A_n\le(nN)^3\).
For sufficiently small \(\eta\) the nearby root \(q_*\) is real,
positive and lies in the closed \(q\)-range for \(\tau\in[\epsilon,12]\).
Indeed the finite other roots have positive distance from that
range. Thus \(|q-q_*|\le C_3\eta^{1/(nN)^3}\).
Lipschitz continuity of \(D_{n,h}\) in \(q\) on the compact box
gives
\[
                |D_{n,h}(q_*,v)|\le C_4\eta^{1/(nN)^3}.
\]
Apply the root product estimate to \(D_{n,h}(q_*,\cdot)\).
Its degree is at most \(nN\), giving a nearby \(v_*\) at distance
at most \(C_5\eta^{1/(nN)^4}\). There are finitely many algebraic
pairs \(q_*,v_*\) under consideration before imposing the cone
condition. Those that are nonreal or outside the closed positive
cone box for \(\tau\in[\epsilon,12]\) are at positive distance
from it and can be eliminated by reducing \(\eta_0\). All these
distances have effective positive lower bounds where required.
For the \(q\)-endpoints use the transcendence of exponentials of
nonzero rationals. The cone inequalities \(q_*\le v_*\le q_*^{-1}\)
are algebraic comparisons. The logarithms on the positive box
are uniformly Lipschitz. This proves (1) upon taking the worst
exponent and constants over the finite list of terminals and
integers in (5).

The result for \(K\) follows by its positive distance from the
finite exceptional set and from \(\tau=0\). For \(\tau\ge12\),
the [absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
already gives separation. A positive distance can be certified
for a finite union of rational rectangles, and for closed compact
algebraic segments excluded by the logarithm argument, by interval
enclosures of the finite list. The conclusion for the ordinary
augmented profiles holds over every characteristic-zero field
by the rational Bézout identity. Strict even positivity also
handles the pair of orders one and two.

## A local bound when the contacting roots are simple

Fix \(n,N,h\) and one contact \((q_*,v_*)\) of (3). Suppose every
common root is simple in both polynomials. For each corresponding
pair take the analytic root branches \(z(q),w(q)\) of \(H_n,H_N\),
respectively. They satisfy \(w(q_*)=\Lambda^h v_*z(q_*)\).
For real \(q\) near \(q_*\), let \(m\) be the vanishing order at \(q_*\)
of
\[
          \frac{w(q)}{\overline{w(q)}}
                -\frac{z(q)}{\overline{z(q)}}.
\]
It is finite, since an identically zero branch difference would
make \(A_n\) vanish on an interval. The order is effectively
computable by implicit differentiation of the root equations
at the algebraic point, using the conjugate branch for each
conjugated root. Let \(m_*\) be the maximum over the finitely
many common root pairs.

In a sufficiently small neighborhood of \((q_*,v_*)\), any ordinary
collision at the specified \(n,N,h\) and sufficiently small \(\eta\)
has
\[
                   |q-q_*|+|v-v_*|=O(\eta^{1/m_*}).            \tag{8}
\]
To prove this, isolate each of the contacting roots in a disk
where the reference derivative stays bounded away from zero.
A sufficiently small coefficient perturbation moves its simple
root by \(O(\eta)\). This follows by factoring out the reference
root and applying Rouché on a circle of radius a constant times
the coefficient error, uniformly for \(q,v\) in the neighborhood.
Root compactness excludes actual common roots outside the chosen
contact disks once the neighborhood and error tolerance are small.
Thus an actual collision gives
\[
                         |w(q)-\Lambda^h v z(q)|=O(\eta).
\]
The branches are bounded away from zero, and \(\Lambda^h v\) is
positive real. Dividing by conjugates bounds the above angular
quotient difference by \(O(\eta)\). Its nonzero leading Taylor
coefficient then gives \(|q-q_*|=O(\eta^{1/m})\) on this branch.
The same complex root equation gives
\(|v-v_*|=O(\eta+|q-q_*|)\), since \(w/z\) is analytic at \(q_*\).
This proves (8). In particular, at a contact where every branch
difference has a nonzero first derivative, the confinement is
linear in \(\eta\).

## Attainability of the limiting parameters

For any cone point with \(\tau>0\) and any prescribed integer \(h\),
take increasing and decreasing unit-gap rows of common smaller
base \(A\), in numbers
\[
   \left\lfloor (A+1)(\tau+\lambda)/2\right\rfloor,\qquad
   \left\lfloor (A+1)(\tau-\lambda)/2\right\rfloor .
\]
Their two curvature parameters tend to \((\tau,\lambda)\), with
\(\eta\) tending to zero. Correct the signed gap discrepancy \(b\)
by \(|b|\) unit-gap rows of the required orientation with larger
exponent at least \((A+1)(1+|b|)^2\). The correction to both
parameters tends to zero. This checks that an exceptional reference
point can be approached by ordinary integer tails at its required
gap. Realizing a reference contact by an actual ordinary common zero
requires an exact argument.
