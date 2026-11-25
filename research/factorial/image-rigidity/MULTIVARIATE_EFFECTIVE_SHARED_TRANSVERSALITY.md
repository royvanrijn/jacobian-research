# Computable transversality margins for shared stationary phases

For distinct ordinary fully shared exponent columns, remove equal
rows and put \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\),
\(r_i=d_i/(2m_i)\), and \(W=\sum_i d_i^2/m_i\). Use the
chord-normalized coefficient phase \(G\) and the upper-half-plane
stationary functions
\[
 X(q)=G'(q)+\log(1-q)-\log q,\qquad
 P(q)=G(q)-qG'(q)-\log(1-q).
\]
Thus \(P'(q)=-q X'(q)\). The centered target of a complex physical
coefficient ratio has \(-\pi<\Im w<0\) when the ratio lies in the
upper half-plane.

Define two closed rays in the target plane by
\[
 L_+=\{x-i\pi/2:x\ge 1-\log 2\},\qquad
 L_-=\{x-i\pi/2:x\le \log 2-1\}.
\]
The finite endpoints are \(w_\pm=\pm(1-\log 2)-i\pi/2\).
For \(\eta>0\), say a profile satisfies the marked-row margin
if every row with unordered exponents \((1,2)\) is accompanied
by at least \(\eta\) curvature from the other rows. Profiles
without that row satisfy the margin.

**Theorem.** Fix positive rational \(W_-\le W_+\), a finite
rational real-field bound \(M\), a rational interior strip margin
\(0<\epsilon<\pi/2\), and a rational minimum saddle separation \(d>0\).
There are the following two computable strict bounds on
\(W_-\le W\le W_+\), \(|\Re w|\le M\),
\(-\pi+\epsilon\le\Im w\le-\epsilon\).

- For two upper saddles satisfying \(X(q_1)=X(q_2)=w\),
  \(\Re P(q_1)=\Re P(q_2)\), and \(|q_1-q_2|\ge d\),
  there is a computable positive lower bound on
  \(|\Re(q_1-q_2)|\) if the profile satisfies a prescribed positive rational
  marked-row margin, or if \(w\) stays a prescribed positive rational
  distance from \(L_+\cup L_-\).
- For an upper critical preimage \(X'(q_c)=0\),
  \(X(q_c)=w\), there is a computable positive lower bound on
  \(|\Im X''(q_c)|\) if the profile satisfies a prescribed positive rational
  marked-row margin, or if \(w\) stays a prescribed positive rational
  distance from \(w_+,w_-\).

Each alternative gives a bound uniform in all exponent rows and their
number. To use it at factorial zeros, the moment contour must supply
the contributing saddles and an effective expansion error.

## The exact exceptional geometry

The [exterior-quadrant count](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md#a-bound-on-fibers-in-the-exterior-quadrants)
holds for the centered map, since centering translates the target
by a real number. It holds also for retained and diffuse limits.
A pair of preimages in an open quadrant, counted with multiplicity,
persists in interior disks under approximation by ordinary phases.
Thus any vertical two-saddle fiber in the strip, and any critical
preimage in the strip, has real part \(x\in[0,1]\).

Let \(p_i=-a_i/d_i\) be each retained pole,
\(\ell_i=|x-p_i|\), and \(\tau\ge0\) the diffuse curvature.
On a vertical fiber with \(q_j=x+iy_j,\ 0<y_1<y_2\), put
\[
   B(s)=\tfrac12\log\frac{s^2+y_2^2}{s^2+y_1^2}.
\]
The [vertical-fiber identities](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-active-arc-and-vertical-fibers)
come from comparing real parts of \(P\) and \(P+X\):
\[
 \sum_i a_i B(\ell_i)+\tfrac{\tau}{2}(y_2^2-y_1^2)=B(1-x),
 \qquad
 \sum_i b_i B(\ell_i)+\tfrac{\tau}{2}(y_2^2-y_1^2)=B(x).
\]
Since \(B\) decreases and each retained \(a_i,b_i\ge1\), all
\(\ell_i\ge\max(x,1-x)\). For a critical preimage at height \(y>0\),
differentiate the same real parts to obtain
\[
 \tau+\sum_i\frac{a_i}{\ell_i^2+y^2}=\frac{1}{(1-x)^2+y^2},
 \qquad
 \tau+\sum_i\frac{b_i}{\ell_i^2+y^2}=\frac{1}{x^2+y^2}.
\]
The same distance bound follows.

In either case let
\[
 f(y)=\Re X'(x+iy)=
 \tau+\sum_i\frac{|d_i|\ell_i}{\ell_i^2+y^2}
          -\frac{x}{x^2+y^2}-\frac{1-x}{(1-x)^2+y^2}.
\]
The ratio of the positive part to the two subtracted terms is
nondecreasing in \(y^2\). For each numerator and denominator
term, the differentiated ratio has the sign of the difference
of their squared distances. The diffuse constant gives a positive
derivative if \(\tau>0\). When \(0<x<1\) and \(\tau=0\),
the increase is also strict. If \(x\ne1/2\) the two positive
denominator distances differ, and if \(x=1/2\) every retained
pole lies strictly outside \([0,1]\). The retained sum is nonempty
in this case because \(W>0\).

For a vertical pair, \(X(q_1)=X(q_2)\) and \(P'=-qX'\) give
\(\int_{y_1}^{y_2}f(y)\,dy=\int_{y_1}^{y_2}y f(y)\,dy=0\).
A nonzero continuous function with at most one sign change
cannot have both moments zero, by subtracting the location of
the sign change times the first moment from the second. Thus
the ratio must be identically one. In the critical case,
\(X'=0\) and \(\Im X''=0\) give \(f=f'=0\) at height \(y>0\),
so the ratio has value one and derivative zero. The only
possibility in either case is an endpoint \(x=0\) or \(1\),
\(\tau=0\), and all retained distances equal to one. The
value of the ratio requires \(\sum_i|d_i|=1\), giving one row
\((1,2)\) or its reversal.

At \(x=0\) the increasing row gives
\[
 X(iy)=\log(y+1/y)+1-\log4-i\pi/2 .
\]
This is on \(L_+\), and its upper critical point at \(y=1\)
has target \(w_+\). Reversal sends the target to
\(-\overline{X}\), giving \(L_-\) and \(w_-\).

## Computing a positive margin

To use the
[effective fiber bounds](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md),
account for the change from midpoint to chord centering. Put
\(S(z)=\sum_i d_i\log(1+r_i z)\), with the linear diffuse term
in a limiting profile. Then
\[
 G'(q)=S(2q-1)-\int_0^1 S(2t-1)\,dt.
\]
For fully shared rows
\(|r_i|\le R<1\), with \(R\) computable from \(W_+\) by the
[integer row bound](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#effective-center-and-localization-on-the-real-lattice).
The absolute value of the real shift is at most
\(W_+/(2(1-R))\), by integrating the elementary bound
\(|\log(1+r_i z)|\le |r_i|/(1-R)\) for real \(|z|\le1\).
Enlarge the real field bound by this quantity before using the
effective fiber bounds. They put all relevant \(q\) in a
convex rational rectangle
\(D=\{|\Re q|\le R_0,\ h_0\le\Im q\le R_0\}\) with \(h_0>0\).
Let \(K\) be the target rectangle. Under the target-distance
alternative intersect it with the appropriate closed distance
constraint. Write \(d_K(w)\) for the distance to the original
rectangle plus, when appropriate, the shortfall in the
prescribed distance from the exceptional rays or points.
This continuous function is zero exactly on the specified
target set. It is computable and Lipschitz.

For two points in \(D\), define
\[
 {\cal E}_2=d_K(X(q_1))+(d-|q_2-q_1|)_+
   +|X(q_2)-X(q_1)|+|\Re(P(q_2)-P(q_1))|
   +|\Re(q_2-q_1)|.
\]
For the critical-point calculation use
\[
 {\cal E}_c=d_K(X(q))+|X'(q)|+|\Im X''(q)|.
\]
The preceding classification gives strict positivity of each
residual on its whole compact domain. In the marked-row
alternative the exceptional one-row profiles are excluded by
the mass margin. In the target-distance alternative they are
excluded by the location of their targets.

Here are finite interval tests that compute lower bounds for
these residuals. Use the retained and diffuse profile nets from
the effective fiber proof. At a rational cutoff \(\zeta<1/3\)
each retained row with \(|r_i|\ge\zeta\) has bounded integer
exponents and positive curvature at least \(2\zeta\), so its
list is finite. Replace the discarded rows by diffuse mass
\(\tau\) equal to their total curvature. For each retained
list mesh the interval of all allowed \(\tau\), including
its endpoints. The conditions \(W_-\le W\le W_+\) and the
optional marked-row condition are closed rational bounds on
\(\tau\) for this list. The mark has \(|r_i|=1/3\), so it is
always retained and its remaining curvature is preserved
under the approximation before meshing.

Use the profile-net error for \(S\) on a disk covering both
\(2D-1\) and \([-1,1]\). Chord centering at most doubles the
value error for \(S\), and derivatives in \(q\) introduce only
known powers of two. The error for \(X\) and its first two
derivatives tends uniformly to zero with the cutoff and mesh widths.
The same holds for differences of \(P\), since they can be
evaluated with the exact integral
\[
 P(q_2)-P(q_1)=-\int_{q_1}^{q_2}q X'(q)\,dq
\]
on the straight segment in \(D\). This also avoids evaluating
large cancelling terms in the chord constant. Uniform
derivative bounds in \(q\) follow from \(X\)'s rational
derivatives. On \(D\), \(|1+r_i(2q-1)|\) has a uniform positive
lower bound: either \(|r_i(2q-1)|\le1/2\), or the imaginary part
has magnitude at least \(h_0/\max_D |2q-1|\). The entropy
denominators have modulus at least \(h_0\).
These bounds and the total \(W\) control the necessary
derivatives and the interval quadrature error.

For each successively finer net and rectangular mesh of \(D\)
or \(D^2\), evaluate the residual at every center to increasing
precision. Use the uniform error bounds to enclose its variation
over each entire cell. The intervals therefore cover the
original family as well as its retained and diffuse limits.
Their centers belong to the compact admissible family of
retained rows and diffuse mass in the stated curvature band.
Every such profile is a limit of ordinary rows in a slightly wider
positive curvature band, by approximating the diffuse mass by
unit gaps at diverging exponents. The classification also applies
directly to the retained and diffuse family. The marked mass
inequality is closed because the marked rows remain discrete. Strict
positivity and uniform error control imply that some finite
mesh certifies a positive lower bound. Stop at the first
such mesh. At a pair or critical point satisfying the hypotheses
the other terms vanish, leaving the required quantitative margin.

For a regular pair the derivative of the stationary value
difference with respect to \(w\) is \(q_1-q_2\), up to sign.
The first bound is thus a transverse slope bound. At a
critical point let \(\Lambda\) be the Airy coordinate in a
cubic chart. Since \(\Lambda'(w)^3=2/X''(q_c)\), the second bound
also makes each rotation of the Airy direction uniformly
nonhorizontal. More explicitly, if \(|X''|\le M_2\), then
\[
 |\Im\Lambda'(w)|\ge
       \frac{2^{1/3}}{3M_2^{4/3}}|\Im X''(q_c)| .
\]
Here \(|\Im z^3|\le3|z|^2|\Im z|\) for every \(z\).

## Use in a certified regular contour chart

Suppose a regular two-saddle contour chart provides, uniformly for
large \(j\), a nonvanishing normalization of the \(j\)-th moment
of the form
\[
   A_1(w)e^{jP_1(w)}+ A_2(w)e^{jP_2(w)}+R_j(w),\qquad
   |R_j(w)|\le\frac Ej
        \bigl(|e^{jP_1(w)}|+|e^{jP_2(w)}|\bigr),
\]
with computable bounds \(0<m\le |A_i|\le B\). The leading
amplitudes \(A_i\) are the same functions at adjacent orders.
Let \(D_0=P_1-P_2\) and \(T=-A_2/A_1\). At a root let
\(z_j=e^{jD_0(w)}\). Division of the remainder by the second
exponential gives
\[
 |A_1 z_j+A_2|\le(E/j)(|z_j|+1).
\]
For \(j\ge2E/m\) the triangle inequality gives
\(m/(3B)\le |z_j|\le3B/m\) and
\[
 |z_j-T(w)|\le C/j,\qquad C=(E/m)(1+3B/m).
\]
In particular
\(|j\Re D_0-\log|T||\le C_1/j\), where \(C_1=3BC/m\).
This comparison does not differentiate a contour remainder.

For a common physical root the real displacement is
\(w_N-w_n=-\Delta\log(N/n)\), where \(N=n+1\).
Suppose first \(\Delta\ne0\) and that the joining segment of the
two fields stays in the same regular chart with
\(|\Re D_0'|\ge c>0\). Suppose also
\(|(\log|T|)'_x|\le L\) in the real direction.
Subtracting the two real logarithm estimates yields
\[
 (Nc-L)|w_N-w_n|
       \le \bigl(\log(B/m)+3C_1\bigr)/n .
\]
Indeed the left logarithms differ by
\(N(\Re D_0(w_N)-\Re D_0(w_n))+\Re D_0(w_n)\)
and \(|\Re D_0(w_n)|\le(\log(B/m)+C_1)/n\).
For \(Nc\ge2L\), the integer degree difference forces the
left side to be at least \(c/2\). This gives an explicit
terminal threshold.

In a regular saddle chart the transverse slope lower bound
can be drawn from the pair residual above. If the saddles
stay separated by \(d\), the common root makes
\(\Re(P_1-P_2)=O(1/n)\), so the residual still bounds
\(|\Re D_0'|\) at that field for large \(n\). An inverse
derivative bound on the stationary branches then gives
a computable small segment on which the slope keeps its
sign and half that magnitude. Near and far comparisons
can thus be separated in an effective contour cover.

If \(\Delta=0\), the two fields are the same. At their common
root, \(e^{D_0}=z_N/z_n\), so
\[
           |e^{D_0}-1|\le 6BC/(mn).
\]
Consequently a uniform positive lower bound on
\(|e^{D_0}-1|\) in the root-bearing part of the chart supplies
the other terminal threshold. This last bound is a phase
condition on the contour chart, separate from the
real-slope bound proved here.
