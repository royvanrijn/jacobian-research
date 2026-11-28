# Local radial mass caps from interface scores

For an actual balanced planar three-cell Gaussian maximizer at
\(.63\le\rho\le.68\), a bounded region surrounded by one strict winner
cannot lie in \(B_{1/2}(0)\). Consequently small circles around an
interface point in this disk have at least two jumps almost everywhere
in radius. The reciprocal-gradient bound gives a quantitative lower
bound on the largest pair-score slope at the interface point. At a
genuine triple junction this yields a cap on all circular masses near
the junction, including a formula for displacement from the polar
center. A corresponding conditional cap uses the active pair slope
at an ordinary two-label tie.

The statements concern the geometry when an interface point or junction
is present. The location of either relative to the polar center is a
separate question.

## No bounded region surrounded by one winner in the central disk

Use the score and interface conventions of the
[global interface bounds](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux).
Write \(U_i=T_\rho 1_{A_i}\), \(V_i=U_i-\beta_i\), and \(O_i\) for
the open strict winner set of \(V_i\). Use \(d=1-\rho^2\) and denote
by \(a_e=|\nabla(V_i-V_j)|\) the regular pair slope. The Gaussian
surface measure is \(d\nu_e=\phi_2\,d\mathcal H^1\), with \(d\nu\)
denoting its sum over interfaces in this note.

Let \(S\) be a bounded component of
\(\mathbb R^2\setminus\overline{O_j}\), contained in \(B_R(0)\)
for \(R\le\sqrt2\). Assume \(S\) is nonempty. Its outer boundary
is adjacent to label \(j\) in surface measure, and the interface
portions inside \(S\) separate the other two labels. These are
the bounded components described in the
[subcluster analysis](THREE_CELL_ISOLATED_SUBCLUSTERS.md#minimum-diameter-in-a-central-disk).
Write \(D\) for the union of affected outer and inner interfaces
and set \(W=\max_{i\ne j}(V_i-V_j)\) on \(S\). It is nonnegative
and has zero trace on \(\partial S\). On \(S\cap A_i\) its
gradient agrees almost everywhere with that of \(U_i-U_j\).

Apply Gaussian Gauss–Green to these score gradients in each
inside phase. On the outer boundary the outward trace is \(-a\).
On each inner interface the sum of the two outward traces is
also \(-a\), since the \(U_j\) terms cancel. The drift contribution
uses \(\int_S x\cdot\nabla W\,d\gamma_2
=-\int_S(2-|x|^2)W\,d\gamma_2\), by zero outer trace. Therefore
\[
 \int_D a\,d\nu
   =-\sum_{i\ne j}\int_{S\cap A_i}\Delta(U_i-U_j)\,d\gamma_2
       -\int_S(2-|x|^2)W\,d\gamma_2
   \le \frac{4\rho^2}{e d}\gamma_2(S).                      \tag{JC1}
\]
Here \(e=\exp(1)\). The last inequality uses \(R\le\sqrt2\),
which makes the second term nonpositive. For the first term the
conditional Gaussian derivative bounds the Laplacian of the
smoothed contrast by
\((\rho^2/d)\mathbb E\big||Z|^2-2\big|=4\rho^2/(ed)\), for
standard planar Gaussian \(Z\). The identity for the absolute
moment follows since \(|Z|^2/2\) is exponential of mean one.

Take the speed equal to one on the outer boundary in the orientation
outward from \(S\), with zero speed elsewhere. The two possible
outward label vectors have squared length two and mutual dot product
one. If their surface areas are \(\mu_1,\mu_2\), this gives
\[
 \langle f,K f\rangle\ge
  \tfrac32\min_{B_R\times B_R} M_\rho\,(\mu_1+\mu_2)^2.
\]
The nonnegative extra-coordinate form \(Q^{(1)}\) from (G1)
applies to this bounded compact speed by approximation, so its
local energy is at least \(\rho\langle f,Kf\rangle\).
Also \(a\ge0\) on the inner interfaces. The Euclidean isoperimetric
inequality and the Gaussian density on the disk give
\[
 (\mu_1+\mu_2)^2
   \ge (2\pi)^{-2}e^{-R^2}\operatorname{Per}(S)^2
   \ge (2\pi)^{-2}e^{-R^2}4\pi |S|.
\]
The Mehler density relative to Gaussian marginals satisfies
\(\min_{B_R\times B_R}M_\rho\ge d^{-1}\exp[-\rho R^2/(1-\rho)]\).
Using \(\gamma_2(S)\le |S|/(2\pi)\) in (JC1) and cancelling
the strictly positive area yields
\[
                    e^{-R^2/(1-\rho)}\le \frac{4\rho}{3e}.   \tag{JC2}
\]
In particular \(R=1/2\) is impossible for \(.63\le\rho\le.68\):
the left side is at least \(e^{-25/32}>e^{-1}\), and the right
side is less than \(e^{-1}\). This use of \(Q^{(1)}\) gives
the factor \(\rho\) in (JC2).

## A reciprocal estimate at any central jump point

Let \(\Gamma\) be the closure of the regular jump arcs. For
\(q\in\Gamma\cap B_{1/2}(0)\) put \(s_0=1/2-|q|\).
At least two distinct strict winners occur arbitrarily near \(q\):
take a sequence of regular jump points approaching \(q\) of a
fixed pair type. For almost every \(0<s<s_0\), the circle
\(\partial B(q,s)\) therefore cannot be pure. If it were pure
label \(j\), its almost-everywhere strict wins of \(j\) would
put the circle in \(\overline{O_j}\). A strict win of a different
label near \(q\) then belongs to a bounded component of
\(\mathbb R^2\setminus\overline{O_j}\) inside the circle, forbidden
by (JC2).

The analytic score contrasts and their locally finite semianalytic
stratification give ordinary transverse jump intersections at
almost every radius. To see the exceptional sets directly, the
singular zero-dimensional strata are locally finite and the
critical values of the radius function on each regular arc
have measure zero. The strict winner rule also holds for
almost every point of almost every circle by polar Fubini.
Whenever the resulting circle is nonpure, its cyclic order
of strict winners gives at least two regular jumps. The
semianalytic stratification input is Bierstone–Milman's
Corollary 2.11, used for the three analytic winner sets in the
[junction analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions).
[Primary stratification source](https://www.numdam.org/item/PMIHES_1988__67__5_0.pdf)

For a pair contrast \(T_\rho f\), \(|f|\le1\), its Hessian has
operator norm at most
\((\rho^2/d)\mathbb E|Z_1^2-1|=4\phi_1(1)\rho^2/d<.84\)
on the band. This follows by taking the second conditional
Gaussian derivative in each unit direction. Its gradient norm
is at most \((\rho/\sqrt d)\sqrt{2/\pi}<.75\) by the first
conditional derivative. Let
\[
 L_q=\max_{i,j}|\nabla(V_i-V_j)(q)|,\qquad H=.84.
\]
On regular jumps at distance \(s\) from \(q\) the slope is at
most \(L_q+Hs\). Now use the global reciprocal bound
\(\rho\sum_e\int a_e^{-1}d\nu_e<1\) of (IB1).
The coarea formula for the radius, whose tangential derivative
on each regular arc has absolute value at most one, counts
at least two jumps on almost every such circle. Since
\(c_\nu=(2\pi)^{-1}e^{-1/8}\) bounds the Gaussian density below
on \(B_{1/2}(0)\), we obtain
\[
 2c_\nu\int_0^{s_0}\frac{ds}{L_q+Hs}<\frac1\rho,\qquad
 L_q>\frac{Hs_0}{\exp[H/(2\rho c_\nu)]-1}>.007s_0.         \tag{JC3}
\]
The last estimate uses \(c_\nu>.14\), \(\rho\ge.63\), and
\(\exp(100/21)<120\). For example the latter follows from
\(\exp(5)<149\) and
\(\exp(5/21)>1+5/21+(5/21)^2/2\).

There is a tied-pair version. At a top tie \(V_i(q)=V_j(q)>V_k(q)\)
put \(c=V_i(q)-V_k(q)>0\) and
\(a_{ij}=|\nabla(V_i-V_j)(q)|\). The contrast Lipschitz bound
makes label \(k\) absent on \(B(q,c/.75)\). Repeat the same circle
count there. Every jump is now from \(i\) to \(j\) or conversely, so
use the single-pair bound \(\int_{\Sigma_{ij}} a_{ij}^{-1}d\nu_{ij}
<1/(2\rho)\), which follows by testing (IB1) against the unit
label incidence vector \((e_i-e_j)/\sqrt2\). Thus the reciprocal
integral gives an exponential denominator
\(\exp[H/(4\rho c_\nu)]-1<10\), since
\(\exp(50/21)<\sqrt{120}<11\). The bound is for the active
pair slope:
\[
               a_{ij}>.084\min(s_0,c/.75).                  \tag{JC4}
\]
At any point of \(\Gamma\) at least one top tied pair has
nonzero slope, since some regular pair arc accumulates there
and (G20) forbids such accumulation at zero slope. This
qualitative statement also covers a top tie with zero third
gap that is approached by only two strict winners.

## Local cap and displacement from the polar center

Let \(q\) be a top tie and write \(h=|q|\). For each label \(i\)
choose a comparator \(j\ne i\) such that \(V_i(q)-V_j(q)\le0\)
and \(|\nabla(V_i-V_j)(q)|\ge\ell/2\), with \(\ell>0\) common
to the three choices. This choice is available at an all-score
tie with \(\ell=L_q\), by the triangle inequality for the
three gradients. At a pair tie \(V_i=V_j\ge V_k\), take the
opposite tied label for \(i,j\) and whichever of them differs more
in gradient from \(k\) for the remaining label. This allows
\(\ell=a_{ij}\) when the active pair has nonzero slope.

For \(x=r\theta\) and label \(i\) winning at \(x\), Taylor
expansion about \(q\) with Hessian bound \(H\) implies
\[
 r\theta\cdot n\ \ge\ q\cdot n-\frac{H(r+h)^2}{\ell}
              \ \ge\ -h-\frac{H(r+h)^2}{\ell},
\]
where \(n\) is the normalized comparator gradient. Measured
about the Gaussian polar center this gives
\[
 \max_i p_i(r)\le
 \frac12+\frac1\pi\arcsin\left[
        \min\left(1,\frac h r+\frac{H(r+h)^2}{\ell r}\right)
                         \right].                           \tag{JC5}
\]
The formula applies for almost every \(r>0\). It is a cap
centered about \(q\) itself by using \(h=0\) when defining the
circles about \(q\). In particular, it gives
\(\limsup_{r\downarrow0}\max_i p_i(r)\le1/2\) for a centered
genuine triple junction. At such a junction the distinct
gradients also give exact limiting masses by linearization:
they are the angular masses of the winner regions of the
three linear gradient scores. If the gradients form a
nondegenerate triangle, the mass of a label is its exterior
angle divided by \(2\pi\). If they are collinear and distinct,
the two extreme labels have half the mass each and the
middle label has zero limiting mass.

There is a conditional improvement to radial clipping. On
radii where (JC5) gives \(\max_i p_i(r)\le1/2+\epsilon(r)\),
the centered profile obeys
\[
 |p(r)-\mathbf1/3|^2\le 1/6+2\epsilon(r)^2.
\]
The maximum on the capped simplex is attained at a permutation
of \((1/2+\epsilon,1/2-\epsilon,0)\), by convexity. Thus in a
decomposition \(p-\mathbf1/3=\alpha\phi+g\) this gives the
pointwise lower bound
\[
 |g(r)|\ge
  \left(|\alpha|\,|\phi(r)|-\sqrt{1/6+2\epsilon(r)^2}\right)_+ .
\]
The improvement only enters on the radii for which (JC5) is
effective. The universal slope lower bound in (JC3) is small,
and a displaced tie also pays the term \(h/r\).

Finally, if \(q\) is a nearest point of \(\Gamma\) to the
Gaussian origin, the disks of radii smaller than \(|q|\) are
pure up to Gaussian null sets: the three indicator
distributional derivatives vanish there, so each indicator
is constant almost everywhere on the connected disk. On
larger radii (JC5) applies with any available top-pair slope
at \(q\), or with \(L_q\) at an all-score tie. This distinction
does not place a jump or a genuine junction near the origin
for a general competitor. The argument is an internal analytic
bound for actual maximizers.
