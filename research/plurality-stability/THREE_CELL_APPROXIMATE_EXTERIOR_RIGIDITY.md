# Exterior rigidity from approximate tripod scores

Let \(J\subset(0,1)\) be a compact correlation interval. There is
\(\epsilon_J>0\) with the following property. Let \(A\) be an
ordinary balanced planar adjusted-score winner at \(\rho\in J\),
with constants \(\beta_i\) and scores \(U_i=T_\rho 1_{A_i}\).
Let \(C\) be a centered three-sector partition with all three
openings strictly between \(\pi/2\) and \(\pi\), and use the same
label order for \(C\) and \(A\). If
\[
     \limsup_{R\to\infty}\ \sup_{|x|>R}\ \max_i
           |U_i(x)-T_\rho 1_{C_i}(x)| < \epsilon_J,        \tag{ER1}
\]
then \(C\) is a centered regular simplex. The prices of \(A\)
are equal and \(A=C\) outside a compact set, up to null sets.
All three adjusted score differences of \(A\) vanish on
the corresponding entire reference lines. The hypotheses
are the adjusted winner rule and the exterior score bound.
[The shifted-contour theorem for bounded isolated graph ends](THREE_CELL_TAIL_FLATNESS.md)
[Two isolated ends determine a conical exterior](THREE_CELL_TWO_END_CONICAL_EXTERIOR.md)

## A binary half-plane bounds the price difference

Use orthonormal coordinates \((s,t)\) along a reference
outward ray and its normal, with the normal pointing into
cell \(i\). In \(s>0\), \(C\) has label \(i\) above \(t=0\)
and label \(j\) below it. The third label is absent there
because the openings on the two sides of the ray are
obtuse. Write \(\sigma=\sqrt{1-\rho^2}\) and \(a=\rho/\sigma\).
For the Gaussian conditional source \((S',T')\) centered at
\((\rho s,\rho t)\), using the centered binary half-plane
formula instead of \(C\) costs at most
\(\Pr(S'\le 0)\) in each of its scores, uniformly in \(t\).
Consequently, for every sufficiently large \(s\), (ER1) gives
\[
       |U_i(s,t)-\Phi(at)|\le\delta,\qquad
       |U_j(s,t)-\Phi(-at)|\le\delta,\qquad
       U_\ell(s,t)\le\delta                             \tag{ER2}
\]
for all \(t\), where \(\ell\) is the remaining label and
\(\delta=2\epsilon_J\) suffices. The longitudinal threshold
is chosen large enough for all three reference rays.

Put \(b=\beta_i-\beta_j\) and \(\psi(u)=2\Phi(u)-1\).
The positive masses and winner equations imply \(|b|<1\)
since each raw score difference lies in \((-1,1)\).
If \(b>2\delta\), define \(t_0>0\) by
\(\psi(a t_0)=b-2\delta\). For \(s\) above the threshold,
label \(i\) can only win at \(t\ge t_0\), by (ER2) and the
pair winner inequality. Test its score at \((s,0)\) as
\(s\to\infty\). The source outside the forward half-plane
above the threshold has vanishing probability as \(s\)
increases. The normal source has law \(N(0,\sigma^2)\).
Thus
\[
      \tfrac12-\delta
       \ \le\ \limsup_{s\to\infty}U_i(s,0)
       \ \le\ \overline\Phi(t_0/\sigma).
\]
It follows that
\[
 b-2\delta=\psi(a t_0)
          \le 2\Phi\!\left(\rho\Phi^{-1}(\tfrac12+\delta)\right)-1
          \le 2\delta .
\]
Apply the same argument with the labels exchanged to
obtain \(|\beta_i-\beta_j|\le4\delta\). The three reference
rays give this bound for every price difference using
only the raw scores and stationarity.

## Stationarity reduces each end to a bounded graph

On each late forward half-plane, the larger of \(U_i,U_j\)
is at least \((1-\delta)/2\). Its adjusted score is strictly
above that of the spectator when \(\delta\) is small,
by (ER2) and the price bound. Thus only \(i,j\) can win.
Their adjusted difference is within \(6\delta\) of
\(\psi(at)\) everywhere in this half-plane. For a small
enough \(\delta\) depending on \(J\), it has the reference
sign whenever \(|t|\ge1\).

The derivative of the adjusted difference in \(t\) is
strictly positive for \(|t|\le1\). Indeed the second
transverse derivative of each raw score is uniformly
bounded on \(J\) by differentiating the conditional
Gaussian kernel. The corresponding derivative of the
binary model is also bounded. A first difference of
length \(\sqrt\delta\) in (ER2) and the second derivative
bound show that the derivative of the adjusted difference
is within \(O_J(\sqrt\delta)\) of \(2a\phi(at)\). The latter
has a positive uniform minimum for \(\rho\in J\) and
\(|t|\le1\). Choose \(\epsilon_J\) to enforce these
inequalities. Therefore the actual pair interface is one
analytic graph \(t=h(s)\), \(|h|<1\), on a whole late
forward half-plane, with the two incident labels exactly
on its two sides. The constant pair-score equation holds
on the graph by continuity.

The bounded graph straightening theorem applies to each
of these three half-plane ends. Its proof uses only the
constant pair-score equation once there is a bounded
isolated graph. Each end is consequently an exact
straight ray \(t=d_{ij}\) beyond some longitudinal
threshold. The pair-score contrast is equal to
\(\beta_i-\beta_j\) on its whole affine supporting line,
by real analyticity.

For an end along \(e\), the reverse points
\(x=-s e+d_{ij} n\) lie deep in the other two forward
half-planes as \(s\to\infty\). They are at distance bounded
below by a positive multiple of \(s\) from those
longitudinal cutoffs and the corresponding straight walls,
since all three ray angles are strictly obtuse. In that
region the actual labels are the same as the reference
labels, and the reverse direction is interior to the
third reference sector. The conditional mean \(\rho x\)
has the same properties. Hence \(U_i(x)-U_j(x)\to0\) on
the reverse supporting line, and \(\beta_i-\beta_j=0\).
On the forward line the bounded isolated-end formula
then gives
\[
     2\Phi\!\left((1-\rho)d_{ij}/\sigma\right)-1=0
\]
with a possible reversal of sign from the normal
orientation. Therefore \(d_{ij}=0\).

The three outward vectors positively span the active
plane, since their successive gaps are less than \(\pi\).
Their sufficiently late longitudinal half-planes
therefore cover the exterior of some disk. The exact
binary descriptions give \(A=C\) there. The conical
exterior rigidity theorem now applies, since there are
at least two bounded isolated ends. It makes \(C\) the
regular simplex.

## Detecting a noncompact defect in fixed-radius balls

The score hypothesis is implied by a uniform bound on
Euclidean defect volumes at infinity. More precisely,
for some \(v_J>0\), it suffices that
\[
  \limsup_{R\to\infty}\ \sup_{|y|>R}
      |\,\{A\ne C\}\cap B_1(y)\,|_{\mathbb R^2}<v_J.      \tag{ER3}
\]
To see it, choose a fixed \(L\) so the conditional Gaussian
mass outside \(B_L(\rho x)\) is smaller than
\(\epsilon_J/3\) on \(J\). Cover the ball by a fixed
finite number of unit balls and bound the conditional
density by \((2\pi\inf_{\rho\in J}(1-\rho^2))^{-1}\).
All cover centers tend to infinity uniformly as
\(|x|\to\infty\) since \(\inf J>0\). Choosing \(v_J\)
small enough gives (ER1). In particular a stationary
discrepancy of finite total Euclidean area from an
obtuse tripod forces exact regular simplex agreement
outside some disk. Conversely, a noncompact discrepancy
from any centered obtuse tripod has defect volumes
bounded below by \(v_J\) in unit balls along a sequence
of centers tending to infinity.

Two questions for a general maximizer are whether it
satisfies an exterior approximation of this strength
and whether a remaining compact defect can be excluded.
## Two observed ends with one connected common winner

The price and graph arguments also apply when only two forward
half-planes are observed. Fix two distinct pairs sharing one label
and two outward unit vectors with angle strictly between
\(\pi/2\) and \(\pi\). In coordinates \((s,t)\) along each outward
vector and its chosen normal, suppose that the two raw scores
of that pair are uniformly sufficiently close to
\(\Phi(at)\) and \(\Phi(-at)\) for all \(s>S\) and \(t\in\mathbb R\).
The spectator score is their remainder from one. The threshold
for this approximation depends only on \(J\).

With the two pairs in place of three, the derivation of (ER2)
bounds those two price differences by \(4\delta\), and the
remaining difference by \(8\delta\). Taking \(\delta\) smaller
by a uniform factor preserves the strict spectator exclusion
and the graph argument on both half-planes. The two bounded
graphs straighten. As the directions are strictly obtuse, the
two-end pressure theorem places both lines through the origin
and makes all prices equal.

If \(A\) is an actual balanced maximizer and the strict-winner
region of the common label is connected, the two-end conditional
rigidity theorem now makes \(A\) the centered regular simplex
everywhere. The hypothesis puts no condition on the third end
or on the error outside the two observed half-planes.
[Position and prices from two straightened ends](THREE_CELL_ONE_STRAIGHT_END_REFLECTION.md#two-isolated-ends-and-a-second-order-tie)
[Global rigidity with a connected common cell](THREE_CELL_TWO_END_CONDITIONAL_RIGIDITY.md)
