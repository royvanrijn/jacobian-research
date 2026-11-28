# From exclusive fibers to an affine separating line

Let \(A\) be a balanced planar Gaussian noise-stability maximizer
at \(0<r<1\). Suppose almost every vertical fiber contains at most
one of labels 1 and 2, up to null sets. Gaussian rearrangement
then makes the third cell an epigraph, after orienting the
vertical coordinate. Assume its height has a finite continuous
representative \(h\) on \(\mathbb R\).

If a height interval has bounded preimage of positive measure,
there are finitely many switches between the two base labels.
The critical-height argument then makes the base switch unique,
so the two labels are affinely separated. For an unbounded height,
variation diminution alone gives uniqueness as proved below.
[Rigidity for a finite number of switches](THREE_CELL_PROPER_SWITCH_HEIGHT_RIGIDITY.md#critical-height-accumulation-excludes-multiple-finite-switches)

Even with infinitely many switches, the absolute height is
at most linear in the horizontal coordinate.
A globally convex height forces the centered regular
simplex, and a globally concave height is impossible.

The argument applies the zero count for a Gaussian transform
to the portion of the base labels over a height band. Stationarity
gives a zero at every base switch, even if the band is supported far
from that switch.

## Sorting fibers and identifying their base labels

Condition on the two horizontal observations \(x,s\).
Set \(u_x=1_{A_3}(x,\cdot)\) and let \(\ell(x)\) be the
non-3 label of that fiber, choosing it arbitrarily when
the fiber is filled by 3. The term in the conditional agreement
depending on the vertical fiber shapes is
\[
  (1+1_{\{\ell(x)=\ell(s)\}})
             \langle u_x,T_r^{[1]}u_s\rangle_{\gamma_1}.
\]
The other terms depend only on the fiber masses.
Replacing both fibers by upper half-lines of their original
masses attains the two-set Gaussian upper bound simultaneously
for all \(x,s\). The equality case forces the nontrivial fibers
to be half-lines of one common orientation by Fubini, since
the horizontal joint density and the displayed coefficient are
strictly positive. The nontrivial fibers have positive measure
by full active rank. This is the rearrangement in (ASR1),
which does not require the base labels to occupy opposite halves.
[Mossel–Neeman, Theorems 1.1 and 1.3](https://ems.press/content/serial-article-files/32017)
[Affine separation and fiber sorting](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md#rearrangement-of-the-perpendicular-fibers)

The resulting shape is
\[
 A_3=\{y>h(x)\},\qquad
 A_i=\{\ell(x)=i,\ y<h(x)\}\quad (i=1,2)
                                                               \tag{FX1}
\]
up to null sets. Since the assumed representative \(h\) is
finite and continuous, below any fixed compact interval
there is an open vertical strip of heights at which \(A_3\)
is absent. Write \(V_i=T_r1_{A_i}-\beta_i\) for the adjusted
winner scores. On a suitable horizontal line in that strip,
the analytic contrast \(V_1-V_2\) is not identically zero
on the interval. Otherwise it would vanish also as a
two-variable analytic function. For almost every choice
of that horizontal line the strict sign agrees almost
everywhere with the base label, by the winner rule and
Fubini. The zeros on any smaller compact horizontal interval
are finite. Thus \(\ell\) is a piecewise constant labeling
whose genuine switches are locally finite.

Let \(B\) be the set of these switches. It is nonempty
since both labels have positive mass. At every \(b\in B\),
each side label is active just below \(h(b)\), so
\(V_1(b,y)=V_2(b,y)\) along an open vertical wall. Analytic
continuation gives the equality for every real \(y\).
Taking \(y\to+\infty\) and using finite heights yields
\(\beta_1=\beta_2\), by dominated convergence in the
conditional source coordinates. Set
\(\epsilon(s)=1_{\{\ell(s)=1\}}-1_{\{\ell(s)=2\}}\) and
\(\sigma=\sqrt{1-r^2}\). Deconvolution in the vertical
variable gives, for every \(b\in B\),
\[
   \int \epsilon(s)1_{\{t<h(s)\}}p_\sigma(s-rb)\,ds=0
                       \qquad\text{for almost every }t.         \tag{FX2}
\]
Thus the signed distribution of \(h(s)\) under the
measure \(\epsilon(s)p_\sigma(s-rb)\,ds\) is zero.
The half-line equality extends from almost every \(t\)
to every \(t\) by right continuity, and half-lines
determine finite signed measures. In particular,
\[
       T_r^{[1]}\bigl(\epsilon\,1_{\{h\in I\}}\bigr)(b)=0
       \qquad(b\in B)                                           \tag{FX3}
\]
for every Borel height band \(I\), including half-lines
and the entire real line.

## Counting switches with a height band

If a bounded function \(f\) on the real line is nonzero
on a set of positive measure and has at most \(N\) sign
changes after ignoring its zero set, its Gaussian
transform \(T_r^{[1]} f\) has at most \(N\) distinct real zeros.
Factoring off the positive conditional Gaussian terms
reduces the kernel to \(e^{r x s/\sigma^2}\) against
\(e^{-s^2/(2\sigma^2)}ds\). This kernel has positive
determinants on equally sized strictly increasing
real node lists. For each pair of finite lists, shift
both to positive arguments and remove the resulting
positive separable factors. The exponential series and
Cauchy–Binet then give the determinant assertion, with the
strictly positive ordinary Vandermonde term. If there
were \(N+1\) zero rows, take the determinant on those rows
with fixed columns at the sign-change locations of \(f\)
and a final variable column \(s\). Its sign alternates
exactly at those locations. Its integral against \(f\)
with Gaussian weight is therefore nonzero, while expansion
along its variable column makes that integral zero.
This is the same one-dimensional zero-count method
used for folded graph contrasts.
[Folded total positivity and zero counts](THREE_CELL_AFFINE_GRAPH_REFLECTION_TESTS.md#the-number-of-asymmetric-crossings)

If a height band \(I\) has bounded preimage of positive
measure, the contrast in (FX3) is nonzero and has only
finitely many sign changes. The zero count applied
to (FX3) bounds \(|B|\) by that finite number. Suppose
now that there are \(N\) base switches, so \(N+1\) open
base intervals have alternating labels. Each height band
with positive preimage must have positive preimage measure
in **every** base interval. A missing interval would make
its contrast have fewer than \(N\) sign changes, which
contradicts its \(N\) transform zeros. This shows that
every base interval has the same essential height infimum
and supremum, allowing infinite values.

If the global height is unbounded in value, \(N\ge2\) is
impossible. There would be a bounded base interval, and
continuity over its closed interval gives both a finite
upper and lower bound there. Thus \(N=1\), proving
affine separation. If \(|h(x)|\to\infty\) at both ends,
continuity supplies a bounded height band with a bounded
preimage of positive measure and also gives an unbounded
height range, as required.

The same identity describes a remaining infinite-switch
case. If switches escape toward \(+\infty\), the preimage
of each height band of positive measure contains arbitrarily
far right sets of positive measure of both base labels.
Otherwise the contrast in (FX3) has just one sign beyond
some \(S\). If its essential support extends arbitrarily
far right, a positive-measure subset of that sign to the
right of \(S\) dominates every opposite contribution as
\(x\to+\infty\) by exponential tilting in \(s\). If the
support ends at a finite point, local finiteness of switches
makes one sign occupy a left neighborhood of its rightmost
essential support point, and the same argument applies.
Both cases contradict the unbounded sequence of zeros.
The assertion also holds on the left with the signs of the
horizontal coordinates reversed.

## A growth bound even with infinitely many switches

Let \(\beta=\beta_3-\beta_1=\beta_3-\beta_2\).
At any base switch \(b\), the three scores tie at its graph
height. Indeed below the height the two base labels win
their respective sides, above it label 3 wins, and the
height is continuous. It follows that
\(U_3(b,h(b))=(1+2\beta)/3>0\), where \(U_i=T_r1_{A_i}\).
Thus \(\beta>-1/2\), and also \(\beta<1\) since it is a
pair price difference. At every graph point away from
the switches, label 3 ties with the base label \(i\), and
the other label \(j\) has adjusted score at most theirs.
The equality \(U_i=U_3-\beta\), together with \(0\le U_j\le U_i\)
and the score sum one, gives
\[
       0<\frac{1+2\beta}{3}
              \le U_3(x,h(x))\le\frac{1+\beta}{2}<1.          \tag{FX4}
\]
By continuity the same bounds hold at the switches.

Write \(M(t)=\max_{0\le s\le t}|h(s)|\) on the right.
At observation \((x,h(x))\), the conditional mean
is \((rx,rh(x))\). For all large \(x\), the rectangle
\(0\le X'\le rx+R,\ |Y'-rh(x)|\le C\) has probability
arbitrarily close to one for suitable fixed \(R,C\).
If \(r|h(x)|>M(rx+R)+C\), that rectangle is entirely
above or below the graph, contradicting (FX4). Thus
\(r|h(x)|\le M(rx+R)+C\) for all large \(x\), and
taking the prefix maximum gives the same estimate
for \(rM(x)\) upon increasing \(C\) for the bounded
initial interval. Put \(a=R/(1-r)\) and
\(H(t)=M(t+a)+C/(1-r)\). Then \(rH(t)\le H(rt)\) for
large \(t\). Iterating inward to a fixed compact
interval in \(t>0\) shows \(H(t)=O(t)\).
Reflecting \(x\) proves \(|h(x)|=O(1+|x|)\) globally.

## The graph at a base switch

At \(b\in B\) let \(i,j\) be the left and right base labels.
Each contrast \(V_3-V_i\) has strictly positive vertical
derivative, because its source is vertically nondecreasing
and varies on nontrivial fibers. Therefore the height
has analytic left and right extensions to \(b\). The
vertical derivatives of \(V_i,V_j\) there agree by their
equality on the axis. Also
\(\partial_x(V_j-V_i)(b,h(b))>0\) by the oriented wall
and the nonvanishing pair gradient (G20) up to the
junction. Hence
\[
 h'(b+)-h'(b-)=
       \frac{\partial_x(V_j-V_i)(b,h(b))}
            {\partial_y(V_3-V_i)(b,h(b))}>0 .                  \tag{FX5}
\]
In the bounded-height, multiple-switch case, every bounded
base interval must therefore attain the common maximum
at an interior point. This follows since the maximum
cannot occur at any base switch with this strict upward
change of slope.

## Classification for a convex height

If \(h\) is concave, its one-sided slopes have a downward
jump or agree at each base switch. Equation (FX5) rules this
out, since the switch set is nonempty.

Suppose \(h\) is convex. The strict slope jump makes it
nonconstant, and a nonconstant convex function on the
whole line is unbounded in value. Its two monotone tails
have limits in the extended real line. Choose a bounded
open interval in its range avoiding the at most two finite
tail limits, with closure still avoiding them. Its preimage
has positive measure by continuity and is bounded by the
tail limits. The height-band comparison then yields one
base switch \(b\). Equation (FX3) for the entire height
range says a conditional horizontal Gaussian assigns equal
masses to the two base half-lines at observation \(b\).
Its mean is \(rb\), so this forces \(b=0\).

The outward endpoint derivatives \(h'_+(0)\) and \(h'_-(0)\),
with \(h_\pm(s)=h(\pm s)\), sum to the strictly positive
jump in (FX5). At least one is positive. Since both outward
profiles are convex, that profile is increasing everywhere.
The common height law and the strict likelihood-ratio
comparison force the two profiles to agree.
[Reflection from one monotone profile](THREE_CELL_AFFINE_GRAPH_REFLECTION_TESTS.md#the-common-height-distribution)

The resulting reflected maximizer has a globally convex
outward profile. The reflected convex-height rigidity
theorem then makes it the centered regular simplex.
[Reflected convex-height rigidity](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#global-convexity-or-concavity-forces-the-simplex)

A general maximizer is not known to have exclusive
fibers in any direction. The finite-height and height-band
conditions in this reduction also require justification
before it can be applied to arbitrary maximizing contacts.
