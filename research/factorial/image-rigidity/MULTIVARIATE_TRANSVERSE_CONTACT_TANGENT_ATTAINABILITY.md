# Diffuse tail shapes at an interior Gaussian contact

Fix a finite ordinary two-monomial core, adjacent orders \(n,N=n+1\),
and the integer signed gap \(h\) of the added rows. In the
[physical Gaussian reference](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts),
let \(W_*>|\lambda_*|\) and let the two contacting root branches
be simple. Normalize the tail separately at each order by its
positive endpoint scale. Write \(R\) for the local logarithm of the
upper-to-lower normalized root ratio minus the exact log ratio of
the two endpoint scales. Its branch is chosen by \(R=0\) at the
reference. The equation \(R=0\) is a common root in the physical
variable.

For an added unequal shared row write \(d_i=b_i-a_i\),
\(m_i=(a_i+b_i)/2\), \(x_i=d_i/m_i\), \(w_i=d_i^2/m_i\), and
\(q_i=1/d_i\). Use the tail statistics
\[
 \begin{gathered}
 W=\sum_i w_i,\qquad \lambda=\sum_i x_i,\qquad \eta=\max_i|x_i|,
 \qquad K=\sum_iw_i x_i,\qquad Z=\sum_i x_i^2,\\
 B_4=\sum_i w_i x_i^2,\qquad L_2=\sum_i d_i/m_i^2 .
 \end{gathered}
\]
Equal rows affect neither the normalized polynomial nor its
argument scale. For nonempty unequal tails define their probability
measure
\[
      \mu_{\mathrm{tail}}=\sum_i \frac{w_i}{W}\,
                  \delta_{(q_i,x_i/\eta)} .
\]
This measure is supported on the compact set
\[
 \mathcal A=
 \{(\sigma/D,\sigma t):\sigma=\pm1,\ D\in\mathbb Z_{\ge1},\
                                      0\le t\le1\}
            \ \cup\ \{(0,y):-1\le y\le1\}.
\]
The coordinate \(q=0\) describes a limit in which the integer
gap magnitude tends to infinity. In particular,
\(\int q\,d\mu_{\mathrm{tail}}=\lambda/W\).

Suppose the reference angle difference changes sign at \(W_*\).
Then every probability measure \(\mu\) on \(\mathcal A\) with
\(\int q\,d\mu=\lambda_*/W_*\) occurs as a weak limit of tails with
\[
 (W,\lambda,\eta)\longrightarrow(W_*,\lambda_*,0),\qquad
                          \sum_i d_i=h .
\]
These limits hold both for extensions with exact contacts, in which
only two repeated unit-gap row midpoints may be real, and for
ordinary extensions with residual as small as any prescribed positive
function of \(\eta\). The ordinary approximation can simultaneously
satisfy arbitrary prescribed coefficient data at finite prime precision.
In both versions all row gaps and multiplicities are integers. The same assertion holds separately for every point of
the [second-order shape tetrahedron](MULTIVARIATE_TRANSVERSE_TAIL_RESPONSE_POLYTOPE.md#the-joint-skew-quartic-and-center-polytope)
for \((B_4/Z,K^2/Z,L_2/Z)\).

When the reference angle derivative is nonzero, there is also an
explicit necessary first-order region for the simultaneous strength
and bias deviations of exact ordinary contacts. That region is sharp
in both senses of attainability above.

## The joint first-order region

Put \(\ell_*=\lambda_*/W_*\), \(e=1/(2nN)\), and
\(c=(n^{-2}-N^{-2})/12\). Let
\(a=\Delta\alpha,\ b=\Delta\kappa,\ g=\Delta\gamma\)
be the upper-minus-lower logarithmic root responses at the
reference in the
[signed jet expansion](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#root-coefficients-at-an-algebraic-exceptional-clock).
The angular derivative is \(\Im a\).
Consider the convex hull \(\mathcal P\) of
\[
 \begin{split}
   &(q,0,0,0), &&q=0,\ \pm 1/D,\quad D\ge1,\\
   &(\sigma/D,\sigma,1/D,\sigma/D^2), &&\sigma=\pm1,\quad D\ge1,\\
   &(0,\sigma,0,0), &&\sigma=\pm1 .
 \end{split}
\]
The set of generators is compact, and so is its convex hull. For
instance linear dependence reduces every finite convex combination
in \(\mathbb R^4\) to at most five generators, making the hull the
continuous image of a compact product. Denote its bias slice by
\(\mathcal P_{\ell_*}=\{(k,z,p):(\ell_*,k,z,p)\in\mathcal P\}\).
For \(\Im a\ne0\) define \(\mathcal D_*\subset\mathbb R^2\) by the
real pair \((A,B)\) determined from each \((k,z,p)\in\mathcal
P_{\ell_*}\) through the single complex equation
\[
                  a A+ eB + W_*(b k+g z-cp)=0.                 \tag{1}
\]
The map is linear in \((k,z,p)\) because \(\Im a\ne0\) and \(e>0\).

Every exact ordinary contact sufficiently near this reference obeys
\[
 \operatorname{dist}\left(
          \left(\frac{W-W_*}{\eta},\frac{\lambda-\lambda_*}{\eta}\right),
          \mathcal D_*\right)=O(\eta).                        \tag{2}
\]
The constant is uniform in the ordinary tail for this fixed simple
reference pair. More generally, the right side is
\(O(\eta+|R|/\eta)\) for a sufficiently close pair with
\(|R|=O(\eta)\). Conversely every point of \(\mathcal D_*\) occurs as
a limit of those scaled deviations at exact real-base contacts and
at ordinary approximations with arbitrarily small residual as above.
In particular, (2) excludes exact ordinary contacts at scaled
distances greater than a constant times \(\eta\) from the region.

For the exact moment identity write \(y_i=x_i/\eta\). Then
\[
 \left(\frac{\lambda}{W},\frac{K}{W\eta},\frac{Z}{W\eta},
                \frac{L_2}{W\eta}\right)
           =\int (q,y,qy,q^2y)\,d\mu_{\mathrm{tail}}(q,y)
           \in\mathcal P .                                   \tag{3}
\]
At each \(q\), the integrand moves along the segment between the
two indicated endpoint generators. Conversely, any point in
\(\mathcal P\) is the integral of this map for a finitely
supported probability measure on \(\mathcal A\).

For the scaled response put \(X=W-W_*\), \(Y=\lambda-\lambda_*\).
Equations (4)–(5) of the signed jet expansion give
\[
                R=aX+eY+bK+gZ-cL_2+
                          O(X^2+|X|\eta+\eta^2).               \tag{4}
\]
Here \(B_4\le W\eta^2\), \(|K|,Z,|L_2|\le W\eta\), and all
simple-root responses are analytic in \(W\) near \(W_*\).
The imaginary part of (4) gives
\(|X|=O(\eta+|R|)\) by absorbing the terms quadratic in \(X\)
or with the small factor \(\eta\). Its real part then gives
\(|Y|=O(\eta+|R|)\). If \(|R|=O(\eta)\), divide (4) by \(\eta\)
and use (3) with \(W=W_*+O(\eta)\). The actual bias ratio is
\(\ell_*+O(\eta)\). It can be set to \(\ell_*\) in (3) by
mixing with either \((-1,0,0,0)\) or \((1,0,0,0)\),
as appropriate, using weight \(O(\eta)\). The strict condition
\(|\ell_*|<1\) keeps the correcting denominator bounded below.
The inverse real linear map in (1) gives (2).

The support function of \(\mathcal P\), in dual coordinates
\((A_0,B_0,C_0,D_0)\), is
\[
 \max\left(
   |A_0|,\ |B_0|,\
   \max_{\sigma=\pm1}\sup_{D\ge1}
       \left\{\sigma B_0+\frac{C_0+\sigma A_0}{D}
                              +\frac{\sigma D_0}{D^2}\right\}
   \right).                                                   \tag{5}
\]
For each sign the expression is quadratic in \(1/D\). The
endpoints \(0,1\) suffice unless its quadratic has a strict interior
maximum, in which case the two consecutive reciprocals bracketing
that point suffice together with the endpoints. One can also truncate
the generator set at \(D\le L\), keeping its limit points. This
changes its Hausdorff distance by \(O(1/L)\), because each dropped
generator differs from its limit by \(O(1/L)\). Mixing with the same
two endpoint generators shows the same order of approximation for
the interior slice, with a constant depending on \(\ell_*\).
This gives finite approximations to the necessary region with a
controlled error.

## Adding two blocks without changing the loaded rows

The following construction couples prescribed ordinary rows to a
reference contact. Fix \(0<\delta<1\). Suppose a family of primary
ordinary shared tails \(P\) has
\[
  W_P\longrightarrow(1-\delta)W_*,\quad
  \lambda_P\longrightarrow(1-\delta)\lambda_*,\quad
  \eta_P\longrightarrow0 .
\]
Write \(H=\sum_{i\in P}d_i\) for the primary signed gap.
Keep every primary exponent unchanged. To obtain exact contacts
with real bases, append \(R_+\) copies of a positive unit-gap row
and \(R_-\) copies of a negative unit-gap row, where
\[
 R_+=T+\max(h-H,0),\qquad R_-=T+\max(H-h,0)
\]
and \(T\) is a positive integer. The two midpoint controls are
\(m_\pm=R_\pm/s_\pm\), with \(s_\pm>0\). They contribute exactly
\(s_++s_-\) to \(W\) and \(s_+-s_-\) to \(\lambda\), and exactly
\(h-H\) to the signed gap. For a real midpoint the row uses the
Gamma factorial, so its unit-gap ratios are still finite products.

Here \(s_\pm\) lie near the positive numbers
\(\delta(W_*\pm\lambda_*)/2\). In \((W,\lambda)\) coordinates,
\(s_\pm=((W-W_P)\pm(\lambda-\lambda_P))/2\).
There is a fixed small rectangle about the reference \((W_*,\lambda_*)\)
on which these controls are positive for all sufficiently late
primary tails. On this rectangle the normalized log multipliers and
the reduced log relative argument scale of the two unit blocks are
\[
 (s_++s_-)F_j(k)+O(1/T),\qquad
                -e(s_+-s_-)+O(1/T),                           \tag{6}
\]
with their first partial derivatives in the same bounds.
Here \(F_j(k)=-k(j-k)/(2j)\) and the reduced scale has
\((h-H)\log(N/n)\) subtracted. In fact, the
[finite unit-gap products](MULTIVARIATE_DIFFUSE_INTERIOR_CONTACT_PERSISTENCE.md#exact-persistence-with-real-bases-and-integer-row-counts)
have analytic log terms \(f(1/m)\) at \(1/m=0\) with \(f(0)=0\).
For \(s=s_\pm\), a term \(R_\pm f(s/R_\pm)\) differs from
\(f'(0)s\), with its first derivative, by \(O(1/R_\pm)\).
Reflection accounts for the negative sign in the oriented scale.
This proves (6), uniformly in \(H\) since \(R_\pm\ge T\).

The exact ordinary primary contribution to the same data is
\(W_PF_j(k)+O(\eta_P)\) and \(-e\lambda_P+O(\eta_P)\), after
subtracting \(H\log(N/n)\) from its scale. This follows from
the signed jet expansion, and these primary terms are constant
in the two controls. The coefficients and the reduced center of
the combined extension are consequently \(C^1\)-close to their
reference values at \((W,\lambda,h)\), by \(O(\eta_P+1/T)\).
The isolated simple-root branches give the same estimate for
their complex residual. Its reference real part has derivative
\(e>0\) in \(\lambda\) and the reference imaginary part is a
function of \(W\) alone. Solve the real part on a smaller
rectangle as a function of \(W\). The imaginary part on this
slice converges uniformly to the reference angular difference.
Take a short \(W\) interval on which the angle has opposite signs
at the endpoints and its only zero is \(W_*\). For sufficiently
late primaries and sufficiently large \(T\) there is an exact
root contact between those endpoints by the intermediate value
theorem. All such selected roots tend to the reference as
\(\eta_P+1/T\to0\), since the reference angular zero is isolated
and the real part fixes its bias. When the angular zero has order
\(r\), its local Taylor bound also gives
\[
     |W-W_*|+|\lambda-\lambda_*|
                         =O((\eta_P+1/T)^{1/r}),
\]
after shrinking the rectangle if necessary. For \(r=1\) one can
equivalently solve the two real equations by their nonsingular
Jacobian. All constants here are for fixed \(\delta\) and the
chosen reference.

At any such real-base contact, round each midpoint to a
half-integer, keeping the two integer counts and their orientations.
Both controls move by \(O(1/T)\), since \(s_\pm=R_\pm/m_\pm\)
and \(m_\pm\) is of order \(R_\pm\) in the rectangle. The same
normalized logarithm estimate gives \(|R|=O(1/T)\) for the rounded
ordinary tail. The two appended blocks have maximum relative gap
\(O(1/T)\). Their contributions to \(|K|,Z,B_4,|L_2|\) are each
\(O(1/T)\) both before and after rounding. Thus for each fixed
nonempty primary tail far enough along its family, \(T\) can be
made large enough to keep the maximum relative gap exactly
\(\eta_P\), leave all selected normalized primary moments
arbitrarily close to their values with block mass supported at
\(y=0\), and make the ordinary residual arbitrarily small. This
works with no bound on \(h-H\) and keeps \(h\) exact.

## Realizing the prescribed shapes

First take a finite probability combination of atoms
\((q,y)\in\mathcal A\). At a small scale \(\epsilon>0\), assign
each atom of mass \(\beta\) primary curvature
\(\beta(1-\delta)W_*\). For \(q=\sigma/D\ne0\) use gap
\(d=\sigma D\), with midpoint in \(\mathbb Z+d/2\) nearest
\(|d|/(\epsilon |y|)\) when \(y\ne0\). For \(y=0\), take a
midpoint of order \(|d|\epsilon^{-2}\). To realize \(q=0\),
use the sign of \(y\) when it is nonzero and choose
\(|d|\) of order \(\epsilon^{-1/2}\) with the same midpoint
rule. The choice of sign for \(q=y=0\) is immaterial.
Each individual \(w_i\) tends to zero, so rounding the multiplicity
to the nearest integer realizes the assigned curvature in the limit.
All rows are shared for sufficiently small \(\epsilon\).
Add one positive unit-gap row with midpoint nearest \(\epsilon^{-1}\).
Its mass vanishes and it ensures \(\eta_P/\epsilon\to1\), even if no
other atom has \(|y|=1\).

For a general probability measure on \(\mathcal A\), approximate
it weakly by finite atomic probabilities on that compact set, for
example by partitioning it into sets of arbitrarily small diameter.
Choose scales successively small enough to approximate those finite
profiles. If the target mean is \(\ell_*\), the primary statistics
then satisfy the hypotheses of the two-block construction at fixed
\(\delta\). The two added curvature proportions tend to
\(\delta(1\pm\ell_*)/2\) and their scaled relative gaps tend to
zero. The full limiting measure at fixed \(\delta\) is therefore
\[
 (1-\delta)\mu+
       \frac{\delta(1+\ell_*)}{2}\delta_{(1,0)}+
       \frac{\delta(1-\ell_*)}{2}\delta_{(-1,0)} .
\]
Let \(\delta\downarrow0\) and take a diagonal sequence, choosing
the primary and subsequent block sizes sufficiently large for the
accuracy at each stage. This gives the promised \(\mu\) for exact
real-base contacts and for ordinary tails with
\(|R|\le b(\eta_P)\) for any positive function \(b\) on positive
scales. Necessity of the support and bias mean for any limiting
measure follows directly from (3) and compactness. Taking a
finitely supported measure representing a point of
\(\mathcal P_{\ell_*}\), with \(b(\eta)=o(\eta)\), gives
sharpness of (2) by the inverted expansion (4). That expansion
also holds for the real-base version, either by its Gamma
derivation or by (6) with the real-to-Gaussian block errors
chosen smaller than order \(\eta_P^2\).

For the second-order claim let \(M_*=W_*-|\lambda_*|\) and fix
\((u,v,p)\) in the tetrahedron with vertices
\((0,0,\pm1),(M_*/2,0,0),(M_*,M_*^2,0)\).
At fixed \(\delta\) use the
[ordinary tetrahedron attainment](MULTIVARIATE_TRANSVERSE_TAIL_RESPONSE_POLYTOPE.md#the-joint-skew-quartic-and-center-polytope)
for strength and bias \(((1-\delta)W_*,(1-\delta)\lambda_*)\),
with shape \(((1-\delta)u,(1-\delta)^2v,p)\).
Scaling the four vertices shows this target is admissible.
For each primary choose \(T\) so large that the changes in
\(K,Z,B_4,L_2\) are \(o(Z_P)\), for instance at most \(Z_P^2\)
in absolute value. Since \(K_P^2=O(Z_P)\), all three shape
ratios are then preserved asymptotically. Let \(\delta\downarrow0\)
diagonally once again. This second-order attainment is separate
from the measure assertion because \(Z_P/(W_P\eta_P)\) may tend
to zero.

## Compatibility with prescribed finite prime data

For each rounded ordinary profile apply the
[independent real and finite-prime approximation theorem](MULTIVARIATE_REAL_LOCAL_COEFFICIENT_APPROXIMATION.md).
It appends ordinary rows of signed total gap zero and arbitrarily
small total curvature while imposing any specified valuations and
unit residues of both constant-normalized factorial coefficient
arrays, at finitely many primes and any chosen finite precisions.
The resulting real coefficient arrays and hence the local residual
can be kept arbitrarily close to those before these extra rows.
Apply the rounding with residual below half the target tolerance
and then choose the arithmetic extension to stay within the other
half.

Choose the added curvature smaller than the primary maximum
relative gap and than any desired error in the shape ratios.
Then the maximum gap remains unchanged and the added
contributions to all the normalized statistics tend to zero at
any prescribed rate. For example, with added curvature \(U\)
one has \(\eta_{\mathrm{added}}\le U\) and
\(Z_{\mathrm{added}},|K_{\mathrm{added}}|,|L_{2,\mathrm{added}}|
\le U^2\), and \(B_{4,\mathrm{added}}\le U^3\).
This suffices even when the primary \(Z\) is much smaller than
\(W\eta\), since the arithmetic extension is chosen after the
primary ordinary profile. The finite primes and precisions can
vary along the approximating sequence. This allows all the
shape limits above together with prescribed finite local data,
for instance both factorial multiplier arrays congruent to one
at every tested prime and precision.

These statements are conditional on the specified finite-core
reference contact. The first-order necessary region uses a simple
transverse pair, and shape attainment also applies at a simple
odd-order angular crossing. It requires strict interior bias.
The ordinary conclusion concerns arbitrarily accurate approximate
contacts. Deciding exact ordinary coprimality inside this region
still requires an additional restriction on the integer exponent
rows or a reference contact exclusion.
