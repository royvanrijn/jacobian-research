# A sharp second-order tail interval at transverse Gaussian contacts

Fix a finite ordinary factorial core, adjacent orders \(n,N=n+1\),
and the integer signed gap \(h\) of the added rows. Use a simple
physical Gaussian reference contact of positive strength \(W_*\)
and bias \(\lambda_*\), with nonzero derivative of the angular
difference. The added shared ordinary rows can have arbitrary
counts, gaps and heights.

For a row with signed gap \(d_i\ne0\) and midpoint \(m_i\), put
\(x_i=d_i/m_i\), \(w_i=d_i^2/m_i\) and
\[
 \begin{gathered}
 W=\sum_i w_i,\qquad \lambda=\sum_i x_i,\qquad
 \eta=\max_i |x_i|,\qquad Z=\sum_i x_i^2,\\
 K=\sum_i w_i x_i,\qquad B_4=\sum_i w_i x_i^2,\qquad
 L_2=\sum_i d_i/m_i^2 .
 \end{gathered}
\]
Consider tails near the reference with \(\eta\) small.
A nonempty tail has \(Z>0\).
Write \(M_*=W_*-|\lambda_*|\).

The exact physical contact equations give an effective interval
for one measured second-order bias residual. The interval is
determined by four vertices of a joint factorial moment
polytope. Its use is uniform in all the added exponents.
In particular it gives an all-degree exclusion for reflected or
moment-balanced tails with small individual curvature whenever
the relevant positive reference contacts have independent real
clock velocities. The heights and multiplicities can vary freely
subject to the stated balance conditions. Simultaneous
simple root pairs also give an exclusion for arbitrary
tails near some higher-degree common factor selections.

## The joint skew, quartic, and center polytope

Set
\[
                    u=B_4/Z,\qquad v=K^2/Z,\qquad p=L_2/Z .
\]
Their exact limiting domain as
\((W,\lambda,\eta)\to(W_*,\lambda_*,0)\), even with arbitrary
fixed integer signed gap \(h\), is the convex hull
\[
 \mathcal T_{M_*}=
 \operatorname{conv}\bigl\{(0,0,1),(0,0,-1),
                             (M_*/2,0,0),(M_*,M_*^2,0)\bigr\}.
                                                                    \tag{1}
\]
Equivalently, for \(s=1-|p|\) it is
\[
        |p|\le1,\quad 0\le u\le M_*s,\quad
                  \max(0,2M_*u-M_*^2s)\le v\le M_*u .             \tag{2}
\]

Here is a quantitative finite-profile bound. Let \(q_i=1/d_i\),
\(z_i=x_i/\sqrt Z\). Then \(\sum z_i^2=1\) and
\(p=\sum q_i z_i^2\). Trim both weights and amplitudes:
\[
    \widetilde w_i=(1-|q_i|)w_i,\qquad
    \widetilde z_i=\sqrt{1-|q_i|}\,z_i .
\]
The trimmed mass is at most \(M=W-|\lambda|\), since
\(w_i|q_i|=|x_i|\), and
\(\sum\widetilde z_i^2\le1-|p|\).
Let \(\widetilde u=\sum\widetilde w_i\widetilde z_i^2\) and
\(\widetilde k=\sum\widetilde w_i\widetilde z_i\).
The [signed curvature triangle](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md#exact-finite-profile-constraints),
with amplitude square sum at most \(s\), gives
\((\widetilde u,\widetilde k^2,p)\in\mathcal T_M\).
Indeed for \(s>0\) divide all amplitudes by \(\sqrt s\)
and scale the two moment coordinates back by \(s\).
At \(s=0\) both trimmed moments vanish.
Furthermore
\[
 \begin{aligned}
   0\le u-\widetilde u&\le2\eta,\\
   |K/\sqrt Z-\widetilde k|&\le\tfrac32\sqrt Z,\\
   |v-\widetilde k^2|&\le3|K|+\tfrac94 Z .
 \end{aligned}                                                     \tag{3}
\]
The first follows from \(1-(1-|q_i|)^2\le2|q_i|\).
For the second use
\(1-(1-|q_i|)^{3/2}\le3|q_i|/2\) and
\(\sum |x_i||z_i|=\sqrt Z\).
The last follows by factoring the difference of squares.
As \(|K|,Z\le W\eta\), both errors in (1) are \(O_W(\eta)\).
The description by four vertices follows since the section
at fixed \(p\) is the untrimmed curvature triangle scaled
by \(1-|p|\), and its origin is already on the edge between
\((0,0,1)\) and \((0,0,-1)\).

For sharpness fix one point of (1). If \(M_*>0\) and \(s>0\),
take the at-most-two weighted amplitudes realizing
\((u/s,v/s)\) in the original triangle, by the
[ordinary curvature triangle construction](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md#sharpness-with-arbitrary-fixed-integer-signed-gap).
Multiply both amplitudes by \(\sqrt s\).
Their weights now sum to \(M_*\) and their squared
amplitudes sum to at most \(s\). For an amplitude \(y\ne0\)
and positive target weight \(w\), at scale \(\epsilon\to0\)
choose integer gap \(d\) nearest \(w/(\epsilon y)\)
and midpoint nearest \(d^2/w\) in \(\mathbb Z+d/2\).
This gives \(w_i\to w\), \(x_i/\epsilon\to y\) and
\(|d_i|\to\infty\). For an amplitude of zero and positive
weight use \(d=\lceil\epsilon^{-2}\rceil\) and the same
midpoint rounding. If \(s=0\), assign all \(M_*\) to an
amplitude of zero. If \(M_*=0\), omit these large-gap rows.

Let \(s_0\) be one minus the sum of squared amplitudes used
so far. It is at least \(|p|\).
Use two sparse unit-gap rows, with signs \(+,-\), at
midpoints giving \(x_i/\epsilon\) tending to
\(\sqrt{(s_0+p)/2}\) and \(-\sqrt{(s_0-p)/2}\), respectively.
Omit a row when its target is zero. Each has vanishing
weight and supplies the designated contribution to
\(\sum q_i x_i^2/\epsilon^2\).
Supply the remaining bias by a group of same-sign unit
rows at a midpoint \(H+1/2\) with \(H\) of order
\(\epsilon^{-4}\), of count nearest
\(|\lambda_*|(H+1/2)\). These rows contribute curvature
\(|\lambda_*|+o(1)\), bias \(\lambda_*+o(1)\), and
\(o(\epsilon^2)\) to \(Z\). Finally adjust the signed gap
to \(h\) by unit rows with total curvature \(O(\epsilon^4)\)
at still larger midpoints, as in the earlier triangle proof.
Now \(\sqrt Z/\epsilon\to1\); the three normalized moments
have the desired limits. This proves the moment-domain
sharpness and does not prescribe exact physical common roots.

## The effective bias residual

Use the response notation of the
[exact signed jet expansion](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#root-coefficients-at-an-algebraic-exceptional-clock).
In particular \(\Delta\alpha,\Delta\kappa,\Delta\gamma,\Delta\sigma,
\Delta\beta^H\) are upper-minus-lower log-root responses, and
\(c=(n^{-2}-N^{-2})/12\) is the signed center coefficient.
Write \(\Delta\alpha'\) and \(\Delta\kappa'\) for strength
derivatives. All responses below are evaluated at \(W_*\).
Set
\[
 t=\Im\Delta\alpha\ne0,\qquad
 \nu=\frac{\Im\Delta\kappa}{t},\qquad
 Q=\Delta\beta^H-\nu\Delta\kappa'
                                   +\frac{\nu^2}2\Delta\alpha' .
\]
Choose any fixed real slope \(s\), and use the real linear form
\[
       \Pi_s(\xi)=\Re\xi-
                 \frac{\Re\Delta\alpha+s/(2nN)}{t}\Im\xi .
\]
The slope lets one test deviations relative to the line
\(\lambda-\lambda_*=s(W-W_*)\). Without a specified line
one may simply take \(s=0\).
Define the measured residual and affine response by
\[
 \begin{aligned}
 R_s&=\frac{\lambda-\lambda_*-s(W-W_*)}{2nNZ}
                          +\Pi_s(\Delta\kappa)\frac KZ,\\
 A_s(u,v,p)&=\Pi_s(\Delta\gamma)+\Pi_s(\Delta\sigma)u
                                      +\Pi_s(Q)v-cp.
 \end{aligned}
\]

**Theorem.** A physical common root continuing this reference
must satisfy
\[
 R_s\in[-\max_{\mathcal T_{M_*}} A_s,\,
                  -\min_{\mathcal T_{M_*}} A_s]+[-C_0\eta,C_0\eta]
                                                                    \tag{4}
\]
with effective neighborhood and constant \(C_0\) when
\(s\) is computable.
Both ends of the reference interval are obtained by the four
evaluations (1). An original ordinary tail outside that
enlarged interval has no such common root.
There are no tail support or multiplicity restrictions.

To prove the bound, write \(x=W-W_*\).
The imaginary part of the signed jet expansion gives
\(x=-\nu K+O(Z)\) after shrinking the neighborhood,
as in the
[directional comparison](MULTIVARIATE_TRANSVERSE_QUADRATIC_DIRECTIONAL_HEIGHT.md#the-exact-first-order-direction).
In particular \(x=O(\eta)\) and \(x^2=O(Z)\).
The real equation also then gives
\(\lambda-\lambda_*=O(\eta)\).
Taylor expansion to scale \(Z\) and substitution of
\(x=-\nu K+O(Z)\) in the second-order products yield
\[
 \begin{aligned}
 0={}&\left(\Delta\alpha+\frac{s}{2nN}\right)x
    +\frac{\lambda-\lambda_*-s x}{2nN}
    +\Delta\kappa K+\Delta\gamma Z+\Delta\sigma B_4\\
    &\hspace{25pt}+Q K^2-c L_2+O(\eta Z).                         \tag{5}
 \end{aligned}
\]
For example \(Kx=-\nu K^2+O(|K|Z)\), and
\(x^2=\nu^2K^2+O(|K|Z+Z^2)\).
All these errors are \(O(\eta Z)\), as are the cubic
Taylor terms because \(x^2=O(Z)\).
Applying \(\Pi_s\) eliminates the first term.
After division by \(Z\) this gives
\(R_s=-A_s(u,v,p)+O(\eta)\).
By (3), \(A_s(u,v,p)\) is within \(O(\eta)\) of the
range on \(\mathcal T_M\).
As \(M=M_*+O(\eta)\), its four vertex values have the
same error from those at \(M_*\), proving (4).

The errors are effective by the local root isolation and
signed jet bounds of the cited comparisons. For rational
or algebraic \(s\), the response constants are algebraic
and the vertex expressions involve only those constants
and \(M_*\). The contact parameters \(W_*,\lambda_*\) are
computed from the finite algebraic clock list. In particular
strict interval margins can be certified by convergent
enclosures. When a whole compact positive strength range is
tested for a fixed ordinary core, the absolute angular and
curvature cutoffs reduce the check to finitely many reference
contacts and integer signed gaps. The physical exception
theorem gives a computable tolerance away from their
neighborhoods. The stated interval test applies to the
simple transverse branch pairs within them.

## Simultaneous pairs and common factors of degree at least four

Fix a selection of \(r\ge2\) distinct upper-half-plane root
pairs at the same physical reference \((W_*,\lambda_*,h)\)
and the same adjacent orders. Both roots in every selected
pair are simple. Regard \(\mathbb C^r\) as a real vector space,
write \(a=(\Delta\alpha_i)\), \(b=(\Delta\kappa_i)\),
\(g=(\Delta\gamma_i)\), \(f=(\Delta\sigma_i)\) for their
response vectors, and put \(e=\mathbf 1/(2nN)\). Let \(P\)
be the orthogonal projection onto the real orthogonal
complement of \(\operatorname{span}_{\mathbb R}\{a,b,e\}\).

**Corollary.** If \(a,b,e\) are real-linearly independent,
a necessary condition for arbitrarily diffuse exact tails
to have common roots continuing all the selected pairs is
\[
                 0\in\{Pg+u Pf:0\le u\le M_*/2\}.               \tag{6}
\]
If the segment has positive distance from zero, an
effective diffuseness threshold excludes this simultaneous
selection. The gate has no curvature, balance, or count
restriction on individual tail rows. For two pairs the
projection has a one-dimensional range, so strict same
signs at the two endpoints suffice.

Indeed write \(x=W-W_*\), \(\delta=\lambda-\lambda_*\).
The joint version of the signed expansion gives
\[
             ax+e\delta+bK=O(x^2+|xK|+Z).
\]
A fixed real left inverse gives
\(|x|+|\delta|+|K|=O(Z)\) after absorbing the first two
terms in a small neighborhood, using that \(K\to0\).
In particular \(v=K^2/Z=O(Z)=O(\eta)\).
The center term \(-c L_2\mathbf 1\) is in the direction \(e\)
and vanishes under \(P\). The surviving equation after
division by \(Z\) is
\[
                           0=Pg+u Pf+O(\eta).
\]
The discarded terms \(x^2,|xK|,K^2\) are \(O(Z^2)\).
By (2)–(3), \(u\) is in \([0,M_*/2]+[-C\eta,C\eta]\)
for an effective \(C\). This gives (6). The independence
test, projection and strict distance can be computed from
the algebraic response vectors and the contact parameters.
No individual angular derivative condition is needed.

There is also a rank-two gate when \(a,e\) are independent
and \(b=\nu a+\mu e\), with real \(\nu,\mu\). A left inverse
for \(a,e\) now gives \(x=-\nu K+O(Z)\), since it first
bounds \(x=O(|K|+Z)\) and \(K^2=O(Z)\). In the projected
second-order equation the skew response is
\[
       P Q,\qquad Q=(\Delta\beta_i^H)-\nu(\Delta\kappa_i')
                                  +\tfrac12\nu^2(\Delta\alpha_i').
\]
Thus the corresponding necessary condition is
\[
 0\in
 \left\{Pg+u Pf+v P Q:
     0\le u\le M_*,\
     \max(0,2M_*u-M_*^2)\le v\le M_*u\right\}.                     \tag{7}
\]
The image is a triangle or its degeneration. Here the errors
are \(O(\eta)\) after dividing by \(Z\), as in (5). Once
again the center term vanishes because \(P e=0\).

A common factor with at least four distinct roots has at
least two conjugate pairs, since one order is even and
its ordinary factorial polynomial has no real roots. In
particular this holds for an irreducible rational common
factor of degree at least four. At a physical reference
whose contacting roots are all simple, two distinct
actual upper-half-plane roots cannot converge to the
same limiting root, because its isolating disk then
contains exactly one root of either nearby polynomial.
Hence any such limiting common factor selects distinct
pairs to which (6) or (7) applies when the rank hypotheses
hold. The gates may be applied to any two-pair or larger
subset. Contact configurations involving coalescence at
a multiple limiting root are outside these tests.

## Uniform reduction to one quadratic factor

**Corollary.** Fix a finite ordinary core and a compact strength
interval \(I\subset(0,12]\). Suppose every physical Gaussian
reference contact with strength in \(I\) has simple contacting
roots, and suppose (6) or (7) excludes every selection of two
distinct conjugate pairs at each reference. For all sufficiently
diffuse ordinary tails with \(W\in I\), every adjacent rational
gcd is either constant or a single irreducible quadratic of
multiplicity one.

Indeed the absolute cutoff restricts the order to finitely many
possibilities. The [physical exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
bounds the signed gap on any sequence of common roots with bounded
curvature, and localizes their limits to the finitely many physical
references. The positive midpoint normalization used here differs
from the maximum-exponent normalization in that theorem by
\(O(\eta)\) in strength and bias. This follows by summing
\(0\le d_i^2/m_i-d_i^2/\max(a_i,b_i)\le
(|x_i|/2)d_i^2/m_i\) and its bias counterpart. If there were
counterexamples to the corollary at arbitrarily small \(\eta\),
a subsequence would have fixed order and signed gap and would
converge to one of these references. Simplicity implies that each
limiting common root has at most one nearby root counted with
multiplicity at either order. A nonconstant gcd has no real
roots by the even factorial moment representation. A gcd other
than a single simple irreducible quadratic would therefore have
two distinct conjugate pairs. Their distinct limiting pairs
contradict (6) or (7). This gives one diffuseness threshold for the
entire interval and all orders of the fixed-core extension.

## Balanced cases and the independent-clock criterion

If \(K=0\), (3) forces the limiting data into the face
\(v=0\), which is the triangle with
\((u,p)=(0,1),(0,-1),(M_*/2,0)\).
For arbitrary tails with \(K=\lambda=0\) approaching a
zero-bias reference, take \(s=0\), so \(R_s=0\) exactly.
The branch is excluded if the three numbers
\[
   \Pi_0(\Delta\gamma)-c,\qquad
   \Pi_0(\Delta\gamma)+c,\qquad
               \Pi_0(\Delta\gamma)+\tfrac{W_*}{2}\Pi_0(\Delta\sigma)
\]
are all strictly positive or all strictly negative. To justify the face restriction
effectively, \(\widetilde k^2=O(\eta)\) by (3);
the lower inequality in (2) then bounds \(u\) within
\(O(\eta)\) of the \(v=0\) triangle when \(M_*>0\).
The case \(M_*=0\) follows from the upper bound in (2).
This test needs no exact mirroring of the tail.

Suppose \(h=\lambda=K=L_2=0\). Here angular
transversality is unnecessary if
\(\Delta\alpha\ne0\).
Indeed the full complex equation at the fixed bias first
gives \(x=O(Z)\), by isolating its nonzero linear term.
After projection to any nonzero real linear form \(\Pi\)
annihilating \(\Delta\alpha\), the exact contact
equation gives
\[
                0=\Pi(\Delta\gamma)
                            +\Pi(\Delta\sigma) u+O(\eta).
\]
For an exactly reflected tail one always has \(u\le W/2\),
since the two members of each reflected pair have equal
curvature and \(B_4/Z\) is a weighted average of individual
curvatures. Thus matching strict signs at \(u=0,W_*/2\)
exclude arbitrary reflected tails near the positive contact.
If the two individual clock velocities \(\alpha_n,\alpha_N\)
are real-linearly independent, then
\(\Pi(\Delta\gamma)\ne0\). In fact
\(\Delta\gamma=-\alpha_N/(2N)+\alpha_n/(2n)\) and
\(\Pi(\alpha_N)=\Pi(\alpha_n)\ne0\).
Consequently it also suffices to require a sufficiently
small fixed bound on each actual row curvature, since
\(u\le\max_i w_i\).
It follows that sufficiently diffuse unit-gap tails with
\(h=\lambda=K=0\) are excluded at an independent-clock
contact, since \(L_2=K\) and \(\max_i w_i=\eta\) for these
tails. Such tails need not be reflected. For example use
four half-integral midpoints \(m_b=(2A+1)b/2\),
\(b=1,3,5,7\), with signed multiplicities
\(1,-27,75,-49\). For a negative multiplicity use its
absolute value of rows with gap \(-1\). These signed
multiplicities annihilate \(1,b^{-1},b^{-2}\), giving
\(h=\lambda=K=L_2=0\) exactly. The heights are distinct,
and there are no reflected pairs. Repeating this group
with multiplicity proportional to \(A\) lets \(W\) tend
to any prescribed positive strength while \(\eta\to0\).

If every branch pair at the finitely many positive zero-bias,
zero-gap contacts of a fixed core on a compact positive
strength interval passes the relevant check, the resulting
exclusion is uniform on that interval in every degree.

The interval (4) constrains arbitrary ordinary tails at
a transverse contact. It does not exclude tails that
satisfy the resulting bias tube and whose arithmetic
budgets still permit a common factor.
