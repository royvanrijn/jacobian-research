# Junction degree and separated components of near-one Gaussian maximizers

Fix \(k\ge3\) and write \(d=k-1\), \(H=\mathbf1^\perp\subset
\mathbb R^k\) and \(e=\sqrt{2(1-\rho)}\). Work in the active
space \(H\) of an actual balanced maximizer \(A^\rho\). Align the
reference simplex \(T\) with its cell selected by the largest
coordinate, using the near-one compactness in (CT19). Write
\(L=L_\rho(A^\rho;T)\) for its score regret. The adjusted scores
\(s_i=T_\rho1_{A_i^\rho}-\beta_i\) are real analytic, with
\(A^\rho\) given almost everywhere by their winners. The results
below use the winning-label error bound (WD2) and the local
comparison after volume repair (WD5).
[Near-one localization and local comparison, (CT19) and (WD1)–(WD5)](THREE_CELL_EXTERIOR_TUBES.md#volume-density-and-outgoing-transfers-of-actual-maximizers)

## Degree on normal fibers of every reference stratum

For a fixed \(M\ge1\), let
\(\delta=C_{k,M}(e+L^{1/d})\), with the constant from (WD2)
on \(B_{M+1}\). Thus if \(j\) is an adjusted-score winner at
\(x\in B_{M+1}\), then \(\max_i x_i-x_j\le\delta\).
Let \(I\subseteq\{1,\ldots,k\}\) have size \(m\ge2\), and let
\(x_0\in B_M\) have its reference top coordinates attained
precisely at \(I\). Assume their common value exceeds every
spectator coordinate by more than \(\delta\), with this
assumption vacuous when \(m=k\).

Let
\[
 H_I=\{z\in H:z_j=0\text{ for }j\notin I\},\qquad
 F_I(z)=P_{H_I}s(x_0+z).
\]
Here the projection retains the \(I\) coordinates, subtracts
their mean on \(I\) and sets the other coordinates to zero.
For every \(\sqrt m\,\delta<r\le1\), the Brouwer degree satisfies
\[
                 \deg(F_I,B_r\cap H_I,0)=1.                  \tag{JD1}
\]
In particular that transverse disk contains an exact winning
tie of all labels in \(I\). Every zero in \(B_1\cap H_I\)
lies within distance \(\sqrt m\,\delta\) of the fiber center.
For \(I=\{1,\dots,k\}\), take \(x_0=0\) to obtain a full
\(k\)-label tie.

To prove (JD1), first note that all winners belong to \(I\)
throughout the disk. The spectator coordinates are fixed
while the largest of the coordinates in \(I\) is at least
their original common value, since \(\sum_{i\in I}z_i=0\).
Thus a spectator winner would violate (WD2). On
\(|z|=r\) in \(H_I\) the spread
\(\max_{i\in I}z_i-\min_{i\in I}z_i\) is at least \(r/\sqrt m\).
If \((1-t)z+tF_I(z)=0\) for \(0<t<1\), the score vector
\(F_I(z)\) is a negative multiple of \(z\). An index attaining
the minimum coordinate of \(z\) would then be a winner,
contrary to (WD2). At \(t=1\) a zero would make every label in
\(I\), including a minimum coordinate, a winner, giving the
same contradiction. The identity at \(t=0\) avoids zero on the
boundary. Homotopy invariance gives degree one. At an actual
zero, all label coordinates in \(I\) lie within \(\delta\) of
their maximum by (WD2); their perturbations sum to zero so
\(|z|\le\sqrt m\,\delta\).

If the zeros in the fiber are nondegenerate, the sum of their
Jacobian signs is one. For arbitrarily small regular target
values of \(F_I\), their preimages in the disk are finite and
have the same signed count: the trace on the boundary is
compact and avoids zero, so the degree is constant in a
neighborhood of zero, and Sard's theorem supplies the regular
values. The degree calculation itself allows a singular
zero set. It gives neither uniqueness of the tie nor a
comparison of its nonlocal cut with the reference junction.

## No small bounded island far from the rest of its label

Fix \(M\ge1\). There are \(c,C,m_0>0\), depending only on
\(k,M\), with the following property for \(\rho\) sufficiently
close to one. Let \(C_0\subset A_i^\rho\cap B_M\) have mass
\(m_C=\gamma_d(C_0)\) with \(0<m_C\le m_0\). Suppose \(C_0\)
is separated from \(A_i^\rho\setminus C_0\) by at least \(De\),
meaning \(|x-y|\ge De\) for almost every pair in those sets.
Then
\[
                   D\le C\sqrt{\log(2+m_C^{1/d}/e)}.         \tag{JD2}
\]
The claim applies to any selected positive-measure subset of the
label, including a connected component or group of components
when the indicated essential separation holds. Its effective
separation bound depends on the volume of the island.

**Remote volume repair.** The argument for (WD5) gives
\[
   S_\rho(U)-S_\rho(A^\rho)\le C_M e\,\gamma_d(A^\rho\ne U)
                                                               \tag{JD3}
\]
for any ordinary change in \(B_M\) of sufficiently small changed
mass. To use the argument on this ball, choose the fixed family
of spanning-tree repair patches beyond \(B_{M+2}\), where every
reference facet has a relative interior of positive area.
At the reference simplex the compactly supported smooth fluxes
on these patches give the full-rank mass matrix of tree incidence
columns. The bulk formula transfers invertibility uniformly to
\(A^\rho\) by its measure convergence, and the second mass
derivatives are uniformly bounded by smooth density pullback.
The affected patches are unchanged in \(U\), so the parameter for
exact volume repair is \(O_M(\gamma_d(A^\rho\ne U))\), provided
this mass is in a fixed small inverse-function neighborhood.

The derivative of the disagreement of \(A^\rho\) along these
flows is \(O_M(e)\) uniformly for small parameters. This is the
pullback estimate used in (WD5). The midpoint and difference
Gaussian representation absorbs the density derivative and
Jacobian in a constant times the density at noise length \(2e\).
The latter disagreement is \(O(e)\) by its Hermite comparison to
noise length \(e\). The contribution to the repair from pairs with
one endpoint in the change set differs by at most
\(C_M e^{-c_M/e^2}\gamma_d(A^\rho\ne U)\), as the changed set and
repair patches are a fixed positive distance apart. Applying
optimality to the repaired ordinary competitor gives (JD3).
[Volume repair and the pullback derivative estimate (WD5)](THREE_CELL_EXTERIOR_TUBES.md#local-comparison-after-volume-repair)

**Exit probability.** Write
\(\mu_\rho(E,F)=\mathbb P(X\in E,Y\in F)\) for the joint
Gaussian pair probability, and put
\[
     E_0=\mu_\rho(C_0,C_0^c),\qquad
     B_0=\mu_\rho(C_0,A_i^\rho\setminus C_0).
\]
Let \(v=|C_0|\) be the Lebesgue volume. Write
\(d_{C_0}(h)=\|1_{C_0}-1_{C_0}(\cdot+h)\|_{L^1(H)}\).
For positive integers \(\ell\), translation invariance and
the triangle inequality give \(d_{C_0}(\ell h)\le\ell d_{C_0}(h)\).
Set \(\ell=\max\{1,\lceil (2v/|B_1|)^{1/d}/e\rceil\}\).
The overlap \(O(h)=|C_0\cap(C_0-h)|\) has
\(\int_H O(h)\,dh=v^2\) by Fubini. The average of
\(O(\ell h)\) for \(h\) uniform on \(B_e\) is therefore at most
\(v^2/|B_{\ell e}|\le v/2\). Since \(d_{C_0}=2(v-O)\),
the average of \(d_{C_0}(h)\) over that ball is at least \(v/\ell\).

For \(z=(Y-X)/e\), the midpoint \((X+Y)/2\) and \(z\) are
independent centered Gaussian vectors with respective covariances
\(a^2 I\), \(I\), where \(a^2=(1+\rho)/2\). Hence the joint
density in \((X,z)\) is \(\phi_a(x+ez/2)\phi_1(z)\).
It is bounded below on \(x\in C_0,\ |z|\le1\) by a positive
constant depending only on \(k,M\). The one-sided Lebesgue
exit for displacement \(ez\) is half of \(d_{C_0}(ez)\),
since the two translates have the same volume. Integrating
over \(B_1\), and using the density comparison in \(B_M\),
gives
\[
       E_0\ge c_M v/\ell
          \ge c'_M m_C\min\{1,e/m_C^{1/d}\}.                 \tag{JD4}
\]
Meanwhile, conditional on \(X=x\), the displacement
\((Y-X)/e\) has covariance \(a^2 I\) and mean \(-ex/2\).
The separation gives \(B_0\le C_M e^{-cD^2}m_C\), with
\(C_M\) enlarged to cover small \(D\).

**Whole-island recoloring.** Choose \(j\ne i\) maximizing
\(\mu_\rho(C_0,A_j^\rho)\), and relabel the entirety of \(C_0\)
by \(j\), obtaining the ordinary partition \(U\). Agreement of
pairs with both endpoints in the island is unchanged. By symmetry,
the two ordered types with exactly one endpoint in the island give
\[
 \begin{split}
 S_\rho(U)-S_\rho(A^\rho)
    &=2[\mu_\rho(C_0,A_j^\rho)-B_0]\\
    &\ge \frac{2}{k-1}(E_0-kB_0).
 \end{split}                                                   \tag{JD5}
\]
The last inequality averages over the \(k-1\) available labels.
Combining (JD3)–(JD5) gives
\[
       c_M\min\{1,e/m_C^{1/d}\}\le C_M e+C_M e^{-cD^2}.
\]
For \(m_C\le m_0\) and \(e\) sufficiently small, the first term
on the right is at most half the left side. In the branch
\(e<m_C^{1/d}\) this uses \(C_M m_0^{1/d}\le c_M/2\), and in
the other branch it uses \(C_M e\le c_M/2\). Solving the
resulting exponential bound gives (JD2).

**Fixed physical separation.** For any fixed \(M\) and
\(\eta>0\), all sufficiently near-one maximizers have no
positive-measure component group \(C_0\) supported in \(B_M\)
at essential distance at least \(\eta\) from the rest of its
label. To see that the small-mass hypothesis eventually holds,
suppose such groups exist along \(\rho_j\uparrow1\) and pass to
a subsequence with fixed label \(i\). Align the maximizers to
\(T\), and take weak-* limits \(g\) and \(h\) of the group and
remainder indicators in \(L^\infty(\gamma_d)\). They are
nonnegative and satisfy \(g+h=1_{T_i}\) by measure convergence.
Their essential supports are at distance at least \(\eta\).
Indeed, sufficiently small balls centered at a pair of support
points closer than \(\eta\) would have both limiting integrals
positive, hence both integrals positive along the sequence,
contrary to separation. The two closed essential supports cover
the connected interior of \(T_i\), which cannot meet both of
them. The remainder's support meets that interior outside
\(B_M\), where \(h=1_{T_i}\), so \(g=0\) almost everywhere.
Consequently the component group masses tend to zero. Now
(JD2) rules out \(D=\eta/e\), as its right-hand side is
\(O(\sqrt{\log(1/e)})\).

The same displacement and recoloring arguments apply to a
Euclidean minimizer of the unit Gaussian cut under ordinary
compact changes. Use Lebesgue measure with density
\(\phi_1(y-x)\) for the pair interaction, take \(e=1\) in the
displacement and omit the midpoint weight and the volume-repair
term. A bounded island of volume \(v>0\) then necessarily
satisfies \(D\le C_k\sqrt{\log(2+v^{1/d})}\).
[Relative Euclidean cut and compact changes](EUCLIDEAN_SIMPLEX_JUNCTION_TRANSFER.md)

## Microscopic minimizing partitions at the ties

Let \(\rho_j\uparrow1\), \(e_j=\sqrt{2(1-\rho_j)}\), and let
\(A^{(j)}\) be the aligned actual maximizers. Take any bounded
centers \(x_j\in H\) and pass to a subsequence with \(x_j\to x_*\).
Then the rescaled partitions \(f^j=1_{A^{(j)}}(x_j+e_j\,\cdot)\)
have a further subsequence converging strongly in
\(L^1_{\mathrm{loc}}(H)\) to an ordinary partition \(u\).
Every such limit minimizes the Euclidean unit Gaussian cut
against all compact replacements, with no constraint on
the replacement volumes. In particular, if \(x_j\) are
full winning ties from (JD1), then
\[
                    (\phi_1*u_i)(0)=1/k\qquad(1\le i\le k).
                                                               \tag{JD6}
\]
Thus the limiting junction has all \(k\) phases and an exact
tie of the unit Gaussian convolution scores at zero.
[Planar antecedent of the blow-up and compactness argument](THREE_CELL_NEAR_ONE.md#blow-ups-of-moving-maximizers-minimize-the-euclidean-cut)

Start with a local weak-* limit \(u\) of the bounded functions
\(f^j\), passing to a diagonal subsequence on bounded balls.
It is simplex-valued. Prescribe an ordinary replacement inside a
fixed ball in the rescaled variables, keeping \(f^j\) outside.
The change has Gaussian mass \(O(e_j^d)\) and belongs to a
fixed physical ball, so (JD3) bounds its possible increase of
stability by \(O(e_j^{d+1})\). The pair measure, after changing
coordinates and dividing by \(e_j^d\), has density
\[
      \phi_{a_j}\!\left(x_j+\frac{e_j}{2}(p+q)\right)
                    \phi_1(q-p),\qquad a_j^2=(1+\rho_j)/2 .
\]
It tends uniformly on compact sets to
\(\phi_1(x_*)\phi_1(q-p)\). The signed agreement difference
vanishes when both endpoints are outside the replacement
ball, so the error from truncating the other endpoint is
bounded uniformly by a Gaussian tail. On the truncated
product the continuous limiting kernel is uniformly
approximated by finite sums of product functions of \(p\) and
\(q\). Local weak-* convergence then passes both pair terms
to the limit, using the dot product for soft agreement.
The positive factor \(\phi_1(x_*)\) cancels. The result is a
Euclidean comparison of \(u\) with every ordinary replacement
inside a ball.

To allow fractional replacements there, approximate any
simplex-valued profile by ordinary indicators using the same
label volumes in each cell of a refining grid. This
approximation converges weak-* on the ball. The truncated
product argument passes the comparison to that limit, keeping
\(u\) outside. If \(u\) were fractional on a set of positive
measure, two of its labels would have a common positive lower
bound on a bounded positive-measure subset. Transferring a
small amount from one such label to the other on that
subset gives a bounded, compactly supported nonzero \(h\)
such that both \(u+h\) and \(u-h\) are feasible fractional
replacements. The sum of their relative cuts is
\[
             -2\sum_i\int_H h_i(\phi_1*h_i)\,dp<0 .
\]
Strictness follows from the positive Gaussian Fourier
multiplier on \(L^2(H)\). This contradicts the two
comparisons. Thus \(u\) is ordinary, and local weak-*
convergence of the indicators to indicators gives strong
local \(L^1\) convergence. The ordinary replacement comparisons
give the asserted minimum.

At a full winning tie the raw scores are
\(T_{\rho_j}1_{A_i^{(j)}}(x_j)=1/k+\beta_i^{(j)}\), since the
adjusted scores agree, the prices are centered, and the raw
scores sum to one. The prices tend to zero by (WD5)'s
remote-flow estimate. Conditional on \(X=x_j\), the rescaled
second endpoint has law \(N(-e_j x_j/2,a_j^2 I)\). Strong
local convergence and the uniform conditional Gaussian tail
give (JD6).

For ties from a fixed \(I\)-stratum, let \(g_j\) be the
reference gap to the largest spectator coordinate at the
fiber center, and let \(\delta_j\) be the (WD2) error on the
bounded physical region. If \((g_j-\delta_j)/e_j\to\infty\),
then all spectators are absent on every fixed rescaled compact
for large \(j\). Indeed the shift from the fiber center to
the tie has zero spectator coordinates and a nonnegative
maximum on \(I\). The additional displacement \(e_j p\)
lowers the spectator gap by at most \(2e_j|p|\). On a fixed
rescaled compact that still leaves a gap greater than
\(\delta_j\), so (WD2) excludes spectator winners. The same
limit argument then gives \(u_l=0\) for spectator labels and
\((\phi_1*u_i)(0)=1/|I|\) for \(i\in I\), by equal adjusted
scores and the vanishing prices. This includes every fixed
reference stratum point with strictly positive spectator gap.

The microscopic comparison does not identify the
large-distance geometry of these Euclidean minimizers or
classify those with at least four phases.
