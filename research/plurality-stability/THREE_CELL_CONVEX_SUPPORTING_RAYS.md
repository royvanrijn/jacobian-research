# Two supporting rays of one convex Gaussian cell

Fix \(0<\rho<1\). Suppose one cell \(A_k\) of a balanced stationary
planar Gaussian three-cell partition has a closed convex
representative \(K\) with two nonparallel half-rays in its boundary
on supporting lines. Write \(\theta\in(0,\pi)\) for the angle
between their outgoing directions. Then \(\theta\ge2\pi/3\), and
equality forces the regular simplex. For actual nonsimplex
maximizers there is a uniform positive gap above \(2\pi/3\) on
each compact positive correlation band.

If both supporting lines pass through the Gaussian origin and
have distinct opposing labels, the partition is the simplex
at any aperture. An opposing label here is one whose adjusted
score agrees with that of \(A_k\) along the supporting ray.
There is exactly one such label for each ray under these
hypotheses, as shown below.

There is also a local version at the origin. An actual maximizer
with one convex cell having two nonparallel straight boundary
segments issuing from the origin must be the simplex. These
segments need only have positive length.

## The two score limits on a supporting line

Write the two supporting halfplanes as \(n_i\cdot x\ge c_i\)
with inward unit normals, and the corresponding rays as
\(p_i+t\tau_i\) for \(t\ge0\) with unit outgoing tangents.
Each \(\tau_i\) is a recession direction of \(K\). Indeed, for
fixed \(x\in K\) and \(t\ge0\), the segment from \(x\) to
\(p_i+s\tau_i\) contains \((1-t/s)x+(t/s)(p_i+s\tau_i)\)
for \(s>t\). Its limit as \(s\to\infty\) is
\(x+t\tau_i\), which belongs to \(K\) by closedness.
Consequently \(p_i+W\subset K\) for each \(i\), where
\[
W=\{x:n_1\cdot x\ge0,\ n_2\cdot x\ge0\}.
\]
The outgoing directions span \(W\) as a cone, since
\(n_j\cdot\tau_i>0\) for \(j\ne i\).

Put \(U_j=T_\rho1_{A_j}\), and use prices \(\beta_j\) for the
stationary winner scores \(U_j-\beta_j\). On each supporting
ray label \(k\) ties with at least one other adjusted score by
continuity, since \(K\) has interior points up to the ray and
its exterior occupies the other side. The two candidate
contrast functions are real analytic on the full supporting
line. Two discrete zero sets cannot cover the ray, so one
label \(\ell_i\) satisfies
\[
U_{\ell_i}-U_k=\beta_{\ell_i}-\beta_k
\qquad\text{on }p_i+\mathbb R\tau_i.
\tag{1}
\]
This uses the [stationary score characterization and line
argument for a fixed simplex cell](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition).

Let \(\sigma=\sqrt{1-\rho^2}\). Along the forward direction
of line \(i\),
\[
U_k(p_i+t\tau_i)\ \longrightarrow\
q_i=\Phi\!\left(-\frac{(1-\rho)c_i}{\sigma}\right)
\qquad(t\to+\infty).
\tag{2}
\]
For this, sandwich \(K\) between \(p_i+W\) and its supporting
halfplane. A conditional source at observation \(p_i+t\tau_i\)
has Gaussian mean \(\rho(p_i+t\tau_i)\) and covariance
\(\sigma^2 I\). The probability of satisfying \(n_i\cdot x\ge c_i\)
is \(q_i\) for all \(t\). The probability of violating the other
inequality defining \(p_i+W\) tends to zero. In the reverse
direction \(U_k(p_i+t\tau_i)\to0\) because the probability
of satisfying \(n_j\cdot x\ge c_j\), for \(j\ne i\), tends to
zero. Since these limits differ, both opposing labels cannot
satisfy (1) on one line. If they did, the identity
\(\sum_j U_j=1\) would make \(U_k\) constant there.

## Angles and centered supports

Using \(0\le U_{\ell_i}\le1-U_k\) at the two ends of (1) gives
\[
0\le\beta_{\ell_i}-\beta_k\le1-2q_i.
\tag{3}
\]
Thus \(c_i\ge0\) for both \(i\), and \(K\subset W\).
Balance gives
\(\gamma_2(K)=1/3\le\gamma_2(W)=\theta/(2\pi)\).
At equality \(K=W\) up to null sets, and (G16) gives the
simplex. For any central \(120^\circ\) sector \(T_k\subset W\),
\[
\gamma_2(K\mathbin{\triangle}T_k)
\le2\big(\gamma_2(W)-1/3\big)=\theta/\pi-2/3.
\tag{4}
\]
For actual maximizers, the uniform one-cell comparison after
(G16) supplies \(\varepsilon_J>0\) on each compact positive
correlation band \(J\). A nonsimplex maximizer in \(J\) with
the two rays must therefore have
\(\theta\ge2\pi/3+\pi\varepsilon_J\).

Finally, suppose both supporting lines pass through the
origin and \(\ell_1\ne\ell_2\). Then \(q_1=q_2=1/2\), so (3)
makes both price differences vanish. Evaluating both line
equalities at zero gives \(U_k(0)=1/3\), or
\(\gamma_{\sigma^2 I}(K)=1/3\). Since \(K\subset W\) and every
vector of \(W\) is a recession direction, one has
\(tK\subset K\) for every \(t\ge1\). Now
\(\gamma_{\sigma^2 I}(K)=\gamma_2(\sigma^{-1}K)=\gamma_2(K)\)
implies \(\sigma^{-1}K=K\) up to null sets. The two sets are
closed convex sets with nonempty interior, so they agree as
sets. Iteration of multiplication by \(\sigma\) shows \(0\in K\).
The recession property now gives \(K=W\), which is the
regular simplex cell by balance. Apply (G16) once more.

## A convex cell at a central triple tie

Let \(K\) now be any closed convex representative of one cell
containing the origin, and suppose all three adjusted scores tie
at the origin. Write \(\bar\beta=(\beta_1+\beta_2+\beta_3)/3\).
The tie gives
\[
U_k(0)=1/3+\beta_k-\bar\beta .
\]
Convexity and \(0\in K\) give \(K\subset \sigma^{-1}K\), so
\(U_k(0)=\gamma_2(\sigma^{-1}K)\ge\gamma_2(K)=1/3\).
In particular \(\beta_k\ge\bar\beta\). Equality makes
\(\sigma^{-1}K=K\) by closed convexity with nonempty interior.
Convexity with \(0\in K\) then extends this equality to all
positive dilations. Indeed, for any \(x\in K\) and \(t>0\),
some \(\sigma^{-m}x\in K\) lies at distance at least \(t|x|\)
along the same ray. Thus \(tx\in K\) by convexity. Hence \(K\)
is a centered convex cone of Gaussian mass \(1/3\), and (G16)
again gives the simplex.

Suppose in particular that \(K\) has two nonparallel boundary
segments issuing from the origin. Their lines support \(K\).
For each segment choose its outgoing unit tangent \(\tau_i\)
and inward unit normal \(n_i\) as before. Then
\(n_j\cdot\tau_i>0\) for \(j\ne i\). The analytic continuation
in (1) holds from the finite segment, and \(U_k(t\tau_i)\to0\)
as \(t\to-\infty\) using the other supporting halfplane.
As \(U_k(0)>0\), both opposing labels cannot satisfy (1) on
one line. If the opposing labels of the two segments are
distinct, their prices are both at least \(\beta_k\) by the
reverse limits and nonnegativity of the other scores.
The two line equalities give a triple tie at the origin,
so the preceding price constraint forces \(\beta_k=\bar\beta\)
and the simplex conclusion follows.

For actual maximizers the two opposing labels must be
distinct. If one label \(\ell\) opposed both segments,
the adjusted \(k,\ell\) contrast would vanish on two
nonparallel lines, giving zero gradient at their
intersection. On either segment the third adjusted score
differs from the common \(k,\ell\) adjusted score in a
punctured neighborhood of zero, by analyticity and
uniqueness of the opposing label on the full line.
The actual \(k,\ell\) jump arc therefore
approaches the origin, contrary to the
[reciprocal-gradient bound (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions).

## Radial occupancy under exact exterior agreement

Suppose an actual maximizer \(A\) agrees outside a bounded set with
the centered simplex \(T\). The three unadjusted scores agree at the
origin by the
[exact line equalities for a compact simplex defect](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius).
Consequently
\[
\gamma_2(\sigma^{-1} A_i)=\gamma_2(A_i)=1/3
\qquad (i=1,2,3).
\tag{5}
\]
If \(A_i\) and \(\sigma^{-1}A_i\) are nested up to null sets for
any one label, (5) makes them equal up to null sets. Iterating
outward along almost every ray then identifies \(A_i=T_i\)
from the exterior agreement, and (G16) gives \(A=T\).
Thus every nontrivial compact defect has both
\(A_i\setminus\sigma^{-1}A_i\) and
\(\sigma^{-1}A_i\setminus A_i\) of positive Gaussian measure
for every label, and the two measures agree. In particular,
no cell of such a defect can be radially monotone in either
direction from the origin. A convex representative is also
impossible by (G16) and the equal-price central triple tie,
since the [central tripod germ](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md)
places the origin in its closure.

There is a related radial mass constraint. Let \(p_i(s)\)
be the fraction of the circle of radius \(s\) occupied by
\(A_i\). For each label the discrepancy \(p_i(s)-1/3\) has
integral zero against \(s e^{-s^2/2}\,ds\) and
\(\sigma^{-2}s e^{-s^2/(2\sigma^2)}\,ds\), by (5). If a
nonzero discrepancy had at most one sign change, subtract
from the ratio of the two densities its value at that
change point. The ratio strictly decreases on \(s>0\),
so its centered value times the discrepancy has one
strict sign off the zero set. The two integrals force
its integral against the first density to vanish, a
contradiction. Thus each nonconstant radial mass profile
must have at least two sign changes about \(1/3\).
At least two profiles are nonconstant in any nontrivial
compact defect, since they sum to one and
[equal circular masses force the simplex](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#circular-balance-and-compatibility-of-first-angular-moments).
