# Uniform physical diffuse stability on affine bias segments

Let \(C\) be a finite ordinary two-monomial core, with equal rows
omitted, and let \(E\) be any finite added multiset of unequal shared
rows. For a row \(w=(a,b)\), put \(r_w=\max(a,b)\),
\(d_w=b-a\), and \(s_w=d_w^2/r_w\). The added parameters are
\[
        \tau=\sum_E s_w,\qquad
        \lambda=\sum_E d_w/r_w,\qquad
        \eta=\max_E |d_w|/r_w .
\]
In particular \(|\lambda|\le\tau\). Fix a rational curvature
bound \(K>0\) for the core. Take an integer \(J\ge4\) at least
the [absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Fixed segments for tight cores.** Suppose the cores have at most
\(R\) unequal rows and satisfy a common computable tail bound
\[
               \sum_{w\in C:\ r_w>B_t}s_w\le2^{-t}
                    \quad(t\ge1).
\]
Here the integers \(B_t\) increase to infinity and
\(\sum_C s_w\le K\). Fix a compact segment
\[
     \tau\in[A,B],\qquad \lambda=\alpha\tau+\beta,
               \qquad 0<A\le B ,
\]
with real algebraic data, where \(\beta\ne0\) or \(\alpha\) is
irrational. There are computable \(\delta,\eta_0>0\), uniform in
the core, such that \(\eta\le\eta_0\) and distance at most \(\delta\)
from this segment imply ordinary SFC(2) at every degree.
No core membership test for a narrower prescribed family is needed.

**Adaptive segments without tightness.** Drop the row cap and the
tail bounds on the core, retaining only \(\sum_C s_w\le K\).
Fix positive rational endpoints \(A\le B\), a real algebraic
slope \(\alpha\), and a nonempty open rational interval of intercepts in
which all segments \(\lambda=\alpha\tau+\beta\), \(\tau\in[A,B]\),
lie in the interior of the cone \(|\lambda|\le\tau\).
There is a finite list of predetermined rational intercepts in
that interval and computable \(\delta,\eta_0>0\) such that, for
each core, at least one listed segment has the conclusion above.
The selected intercept may depend on the core; it works for the
entire segment and every ordinary tail meeting its tolerance.
In particular the Gaussian strength may range over an arbitrary
compact positive interval while one chooses the signed bias.
This part requires no algebraicity of the core reference data
and also supplies reference margins for summable countably
retained cores, as described below.

Taking \(K=12\) in the adaptive result removes the curvature
restriction on ordinary cores. If an ordinary core has curvature
at least \(12\), the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
already gives angular separation after every augmentation.
Otherwise the \(K=12\) menu applies with the same uniform
tolerances. Thus no dimension, exponent, curvature, or tightness
bound is needed for the adaptive conclusion on ordinary cores.

Both results concern physical common roots and require no bound
on the number or signed total gap of the added rows. The constants
use the computable but unevaluated cutoff \(J\).

## Normalize the core and keep the integer drift

Write \(F_{j,C}\) for the core polynomial with endpoints
both one. Equivalently, for the raw core polynomial \(Q_j\),
put \(G_{j,C}=([z^j]Q_j/Q_j(0))^{1/j}\) and
\(F_{j,C}(z)=Q_j(z/G_{j,C})/Q_j(0)\).
Its coefficients lie in \([0,\binom jk]\) by log convexity
of factorial moments. An empty core has \(F_{j,C}=(1+z)^j\).
For \(N=n+1\) let \(\Lambda=N/n\) and put
\[
 E_{n,C}=\log(G_{N,C}/G_{n,C})-\Delta_C\log\Lambda,\qquad
 \Delta_C=\sum_C d_w .
\]
The center contribution of one row is
\[
 \int_a^b\left[\psi(1+Nx)-\psi(1+nx)-\log\Lambda\right]\,dx .
\]
It has magnitude at most \(2s_w\). Indeed
\(0\le x\psi_1(1+mx)\le 1/m\) gives a bound
\(|d_w|\log\Lambda\), which is at most \(2s_w\) if
\(\min(a,b)<r_w/2\). Otherwise the trigamma lower integral bound
gives \(0\le1/m-x\psi_1(1+mx)\le1/(m^2x)\).
Integration in \(m,x\) bounds the row contribution by
\(|d_w|/(nN\min(a,b))\le2s_w/(nN)\).
In particular \(|E_{n,C}|\le2K\).

At a proposed physical common root, write \(u,u'\) for
its lower and upper coordinates normalized using the entire
ordinary profile. For \(\tau\ge12\) the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
already gives separation, so restrict to \(\tau\le12\).
By the [signed center approximation](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#uniform-approximation-of-the-physical-pair)
for \(\eta\le1/2\),
\[
 \log(u'/u)=H\log\Lambda+E_{n,C}-\lambda/(2nN)+e,\quad
 H=\Delta_C+\sum_E d_w,\quad
 |e|\le\eta\tau(1/n^3+1/(2nN)).
\]
Here \(u'/u>0\). The normalized ordinary root annuli give
\(|\log|u'/u||\le(n+N)\log2\), since each degree-\(j\)
normalized coefficient is at most \(\binom jk\).
Since \(\log\Lambda\ge1/N\), a potential common root needs only
the integers \(|H|\le H_n\), where
\[
                         H_n=\left\lceil N(n+N+2K+24)\right\rceil .
\]
This bound concerns the total drift \(H\), even when the signed
gap of the added rows or that of the core is unbounded.

Write \(q=e^{-\tau/(2nN)}\) and \(v=e^{-\lambda/(2nN)}\).
Let \(F^\tau_{j,C}\) be \(F_{j,C}\) with its \(k\)-th coefficient
multiplied by \(q^{(nN/j)k(j-k)}\). The physical reference
resultant is
\[
 D_{n,H,C}(q,v)=\operatorname{Res}_z
       \left(F^\tau_{n,C}(z),
       F^\tau_{N,C}(e^{E_{n,C}}\Lambda^H v z)\right).             \tag{1}
\]
The finite-core coefficient arrays and \(e^{E_{n,C}}\) are
algebraic. Given the finite range \(2\le n<N<J\), the
\(H\) bound, and bounded positive \(\tau\) and bounded \(\lambda\),
the approximate actual and reference arrays lie in a common
computably bounded coefficient box. The
[diffuse chord and center bounds](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#uniform-approximation-of-the-physical-pair)
give coefficient error at most \(C_1\eta\) there, uniformly across
cores and row counts. Explicitly the log of each normalized
tail coefficient differs from its Gaussian value by at most
\(\eta\tau k(j-k)/j\); the scale error is the \(e\) above.
Multiplying by the bounded core coefficients and using the
exponential mean-value inequality gives the asserted box error.
Fixed-degree resultants are polynomial in those coefficients,
so at an actual common root \(|D_{n,H,C}(q,v)|\le C_2\eta\),
with computable \(C_2\). The references in (1) are also uniformly
Lipschitz in \((\tau,\lambda)\) on any fixed compact neighborhood
of the segments, by differentiating their coefficient formulas.
All constants here use only \(K,J\) and the parameter bounds.

## Uniformity on a fixed algebraic segment

At a reference contact the even Gaussian reference is strictly
positive on the real axis. Its ungauged even core polynomial
is nonnegative there, and lognormal averaging at positive
strength misses its finite real zero set almost surely.
A common nonreal root then gives a coinciding off-diagonal root
quotient. The [angular eliminant](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts)
for the normalized pair is a nonzero polynomial in \(q\) of
degree at most \((nN)^3\), by its binomial specialization at zero.
For a finite ordinary core its coefficients are algebraic.
Then \(q\) at a contact is algebraic, and so is \(v\) since
(1) is a polynomial of degree at most \(nN\) in \(v\) with
constant term one. Along the given segment these numbers would
satisfy
\(\log v-\alpha\log q=-\beta/(2nN)\).
The [logarithm argument for physical exceptions](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts)
rules this out for the stated \(\alpha,\beta\).
Thus every finite ordinary core has a positive reference
resultant margin on the entire compact segment.

To make it uniform and compute it, at step \(t\) enumerate the
finite candidate list of cores with at most \(R\) rows,
\(r_w\le B_t\), total curvature at most \(K\), and the tail
inequalities at indices \(1,\dots,t\). It contains every
truncated admissible core. The
[incomplete-test extraction](MULTIVARIATE_UNIFORM_TIGHT_CORE_DIFFUSE_STABILITY.md#a-terminating-calculation-without-testing-core-membership)
shows that any sequence of candidates at growing \(t\)
has a subsequence with each labelled row constant or
\(r_w\to\infty\) and \(s_w\to0\). For a diverging row its
normalized coefficient log and its center contribution
tend to zero, by the chord bound and \(|E_w|\le2s_w\).
Consequently a sequence with vanishing resultant margin
would limit, after also extracting the parameter and the
finite indices \(n,H\), to a reference contact of a finite
ordinary core on the segment. This proves a uniform
positive margin for the candidate lists at large \(t\).

For each candidate and interval, compute certified lower and
upper bounds for the minimum absolute resultant in (1), using
interval subdivision of an affine parameter in \([0,1]\).
All coefficients and exponentials have effective enclosures.
The minimum is positive, and refining until the lower bound
is at least half the upper bound terminates. Let \(g_t\) be
the smallest computed lower bound. The extraction argument
bounds \(g_t\) away from zero for large \(t\).
Deleting the omitted rows changes the normalized coefficient
logs by at most \(j2^{-t}/4\): such rows have total curvature
at most \(2^{-t}\le1/2\), are shared with smaller exponent at
least half the larger, and the factorial chord bound applies.
The center changes by at most \(2\cdot2^{-t}\).
Using the same fixed-box derivative bounds as above, the
change in (1) is at most a computable \(C_3 2^{-t}\).
Stop when this is less than \(g_t/2\). Then \(g_t/2\) is
a certified margin on the full cores. Taking \(\delta,\eta_0\)
so that the uniform parameter-displacement and coefficient
errors total less than this margin excludes actual collisions.
The absolute cutoff and strict even positivity handle the
remaining degrees, and rational Bézout transfers the result
to characteristic zero.

## Adaptive intercepts with no restriction on the retained rows

For this part enlarge the core data to a compact rational
semialgebraic box with constraints. For each \(2\le j<J\),
allow any real polynomial with endpoints one and coefficients
in \([0,\binom jk]\), imposing nonnegativity on the entire
real axis in even degree. Independently allow the center
scale in (1) to lie in a fixed positive rational interval
containing \([e^{-2K},e^{2K}]\).
All ordinary endpoint-normalized core data belong to this set.
Polynomial nonnegativity is a closed semialgebraic condition
by real quantifier elimination. Strict even reference
positivity at \(\tau>0\) continues to hold throughout the
enlarged set, by the same lognormal averaging.

For each fixed datum the angular quotient eliminant is
still nonzero at \(q=0\), since the endpoint coefficients
are one. Its coefficients may be arbitrary reals, but
it has at most \((nN)^3\) possible \(q\)'s for a contact.
For each \(q\), and each allowed \(H\), (1) is a polynomial
in \(v\) of degree at most \(nN\) and constant term one.
Each eligible \(q,v\) yields at most one intercept
\(\beta=\lambda-\alpha\tau\). Hence at most
\[
                  \sum_{2\le n<N<J}(2H_n+1)(nN)^4
\]
intercepts can give a contact for the datum, over the entire
positive curvature interval. Choose any larger number of
distinct rational intercepts in the specified open interval.
For every datum one of these segments has a positive margin
against all \(n,H\) and all \(\tau\in[A,B]\).

For each intercept the minimum absolute resultant on the
segment is a continuous function of the datum, by joint
continuity over the compact segment and finite indices.
The maximum of the finitely many minima is positive
everywhere in the enlarged compact set, hence has a uniform
positive lower bound. It is computable. Subdivide the rational
data box, retaining the boxes that intersect the semialgebraic
data set by quantifier elimination. Bisect the curvature
interval and refine the numerical enclosures as the data
mesh shrinks. Stop when, for each retained box, one of the
intercepts has a certified positive lower bound throughout
the segment, for all \(n,H\). Uniform positivity and uniform
convergence of the enclosures make this search terminate.
The smallest certified lower bound gives a uniform
\(\delta,\eta_0\) through the same transfer as above.

A countable multiset of ordinary core rows with
\(\sum s_w\le K\) defines coefficient logs and centers by
absolutely convergent sums of its finite-row contributions.
Indeed at most finitely many rows have \(s_w>1/2\), and
the remaining chord magnitudes and centers are bounded
by constants times \(s_w\). The limiting endpoint arrays
and centers are in the enlarged compact set. The adaptive
intercept therefore gives a reference margin also for such
a datum and nearby finite approximations. Algebraicity of
the countable products is not needed for this adaptive
statement.


## Exact ordinary placement of rational curvature and bias

Given rational numbers \(\tau>0\) and \(-\tau<\lambda<\tau\), take
a large integer \(m\) for which \(m(\tau+\lambda)/2\) and
\(m(\tau-\lambda)/2\) are positive integers. Use these many
increasing and decreasing unit-gap rows, respectively, with
larger exponent \(m\). Their parameters are exactly \(\tau,\lambda\).
Replacing one row with larger exponent \(x\) by two rows in its
orientation with larger exponent \(2x\) preserves both parameters
and changes the signed integer gap by one in the direction of
that orientation. Both orientations are available, so repeated
splitting realizes any prescribed integer signed gap. Every
split leaves the maximal relative gap at most \(1/m\), which can
be taken arbitrarily small. This is the
[exact unit-gap realization](MULTIVARIATE_UNIVERSAL_DIFFUSE_STRENGTH_MENU.md#exact-ordinary-placement-of-a-selected-strength)
with both rational masses fixed before splitting.

For a rational slope this realizes rational points of the chosen
segment exactly. For an irrational slope one can take rational
points arbitrarily close to the segment and inside the cone.
Consequently every certified segment tube that meets the cone
interior contains actual ordinary tails at arbitrarily small
relative gap, with any prescribed signed integer gap.
