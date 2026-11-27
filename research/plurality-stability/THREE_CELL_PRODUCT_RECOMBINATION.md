# Product recombination at a balanced three-cell optimum

Fix $0<\rho<1$ and let $s_*$ be the balanced three-cell
Gaussian optimum. It is independent of ambient dimension
once the dimension is at least two, by the dimension
reduction theorem used in the
[global maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
An actual nonsimplex partition tying the simplex at a
first transition is subject to all the comparisons below.

## Finite recombinations using two blocks lose strictly

Let $E,F$ be orthogonal Gaussian coordinate blocks. Suppose
an ordinary balanced partition has a label rule
$C(x,y)=b(i(x),j(y))$ for finite-valued measurable $i,j$.
If it attains $s_*$, it depends on only one of the blocks.
The statement also applies to any finite number of blocks
by grouping them in two.

First, the subspace of invariant translations splits
between $E$ and $F$. Suppose all translations by
$t(u,v)$ leave $C$ invariant up to null sets. For all
sufficiently small $|t|$, every nonnull fiber of $j$
intersects its translate by $-tv$ in positive measure.
This uses local $L^1$ continuity of translations and
the finitely many fibers. Fubini on this intersection
shows that $b(i(x+tu),j_0)=b(i(x),j_0)$ almost everywhere
for each nonnull value $j_0$. Thus the translations by
$(tu,0)$ leave $C$ invariant, first for small $t$ and
then for all $t$ by composition. Subtracting gives the
same assertion for $(0,tv)$. The orthogonal complement,
the active subspace, therefore splits as well.

If both blocks are active at a maximizer, dimension
reduction and full active rank give one active coordinate
in each. In these planar coordinates there are finitely
many distinct horizontal slice functions and finitely
many vertical ones, modulo null sets. This finiteness is
preserved on passing to the active subspace, by Fubini
and the separate invariant directions. Pick one generic
line representing each nonnull slice type. On this
line the analytic pairwise score differences are not
identically zero, since the full pairwise tie sets have
measure zero. The winning labels on the line have only
locally finitely many jumps. The two finite unions of
breakpoints thus give a locally finite rectangular grid
on which the partition has constant labels in each
open rectangle.

At a grid vertex the four incident open rectangles
have four score gradients, one for each of their labels,
evaluated at the vertex. Across a horizontal interface
the horizontal gradient component is unchanged, and across
a vertical interface the vertical component is unchanged,
by stationarity along the edge. These statements also hold
on edges without a label change. The four gradients
therefore have the form $(\alpha_L,\beta_D)$,
$(\alpha_L,\beta_U)$, $(\alpha_R,\beta_D)$,
$(\alpha_R,\beta_U)$. A genuine vertical jump implies
$\alpha_L\ne\alpha_R$: a connected analytic jump segment
cannot approach a point with vanishing pair gradient, by
(G20). Thus that jump continues vertically through the
vertex. The analogous assertion holds for a horizontal
jump. If both directions jump there, all four gradients
are distinct, which requires four labels. Hence genuine
jump segments extend along their full horizontal or
vertical lines, and no two lines of different directions
can occur. This makes the maximizer one-dimensional,
contrary to full active rank. Thus genuine dependence
on both blocks forces strict loss.

### Continuous selectors of a single source

If $i$ is the label of an actual balanced planar maximizer,
the same conclusion holds for a rule $C(x,z)=b(i(x),z)$
with an arbitrary measurable dependence on a Gaussian
auxiliary block. The active-space split still holds:
for a simultaneous invariant translation and small $t$,
intersect each source cell with its translate in positive
measure. Fubini gives invariance of every auxiliary
signature $b(i,\cdot)$ under the auxiliary translation.
Subtract this translation to get source invariance.
If the source is genuinely used, merge its input labels
only when their auxiliary signatures agree almost
everywhere. This is a nontrivial coarsening of three
labels, so it contains at least one singleton cell.
No individual cell of the planar optimizer is
cylindrical, by the [full cell rank result](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions).
Thus the source contributes two independent active
directions. Any genuine use of $z$ adds another and
precludes equality.

The same dimension obstruction holds for two actual
planar optimizer labels on transverse active planes,
even when the planes are not orthogonal. If
$V\cap W=\{0\}$, the map
$x\mapsto(P_Vx,P_Wx)$ is an isomorphism from $V+W$
onto $V\times W$. The induced Gaussian law has everywhere
positive density, so the finite-rectangle translation
argument applies in these coordinates. Genuine dependence
on each input retains a singleton cell or all three
cells of its source, as above. Both active planes are
therefore needed, giving four active coordinates after
the invertible linear change. A balanced lookup table
on these sources that uses both inputs loses strictly.

For independent optimal sources and a bounded number
of copies, the strict loss of genuine lookup tables
is uniform on compact positive correlation intervals.
Use [strong maximizer compactness](POSITIVE_CORRELATION_CONTINUATION.md)
to extract each input shape in its active plane as
the correlations converge. The transition entries
converge by continuity of the Gaussian operator and
strong convergence of the labels. For fixed arity
there are finitely many tables, and their energies
are polynomials in these entries. Hence a sequence
of gaps tending to zero would give a tied, genuine
lookup table at the limit, already excluded.

## Selectors and disagreement energy

For any vector field $f$ write
$Q(f)=\langle f,T_\rho f\rangle=\|T_{\sqrt\rho}f\|_2^2$,
summing over its labels. For arbitrary measurable
sections $F_y$ of a partition in a product of Gaussian
spaces, exact polarization gives
\[
 s_*-S_\rho(F)
   =\mathbb E_y[s_*-S_\rho(F_y)]
       +\tfrac12\mathbb E_{(Y,Y')_\rho} Q(F_Y-F_{Y'}).
                                                               \tag{PR1}
\]
If section masses are balanced, both terms are nonnegative.
The second vanishes only for a constant section almost
everywhere, by injectivity and the positive Gaussian pair
density. The Ornstein–Uhlenbeck gap in $y$ bounds the
second term below by $(1-\rho)\mathbb E_y Q(F_y-\mathbb E F_y)$.
If both directions have balanced sections, applying this
to either direction gives $S_\rho(F)-1/3\le\rho(s_*-1/3)$,
since $\mathbb E_y F_y$ is then constant. Sectionwise
balanced switching between two optimizers cannot tie them.

Different conditional masses can give a stronger test.
Let $a,b$ be the indicator vectors of two distinct
balanced global maximizers in the same Gaussian space.
Put $u=(a+b)/2$, $d=a-b$ and $c(t)=2\arcsin(t)/\pi$.
In an independent extra correlated Gaussian block $Z$,
take the indicator vector
\[
                 C(x,Z)=u(x)+\tfrac12d(x)
                                      \operatorname{sgn}(v(x)\cdot Z),
\]
where $v(x)$ is any measurable unit vector on the
disagreement set. It is ordinary and balanced because
the conditional mean in $Z$ is $u$. Orthogonality of
the mixed terms and the Gaussian sign correlation give
\[
 Q(d)\ge
 \iint M_\rho(x,y)\,d(x)\cdot d(y)\,
          c(\rho v(x)\cdot v(y))\,d\gamma(x)d\gamma(y).
                                                               \tag{PR2}
\]
Indeed $Q(u)=s_*-Q(d)/4$, and the right-hand kernel
is four times the added selector energy. If the unit
vectors take only finitely many values, (PR2) is
strict. The lifted partition then has a finite
rectangular representation and depends on $Z$ on the
disagreement set. Dependence on $x$ is also genuine,
since its conditional mean $u(x)$, supported on at most
two labels, cannot be the constant balanced vector.
Apply the strict recombination criterion.

In particular, use $v=d/\sqrt2$ on the disagreements.
Write $E_e$ for the three unordered disagreement
types, $m=1_{\{a\ne b\}}$, and $C_E=\sum_e Q(1_{E_e})$.
For the same type the signed inner products give
$2c(\rho)$ in (PR2); for different types they give
$c(\rho/2)$. Hence
\[
 Q(d)>
   2c(\rho) C_E+c(\rho/2)(Q(m)-C_E)
   \ \ge\ \tfrac23\big(c(\rho)+c(\rho/2)\big)Q(m).           \tag{PR3}
\]
The last step uses $C_E\ge Q(m)/3$ by the triangle and
Cauchy–Schwarz inequalities for $T_{\sqrt\rho}$.

Scalar signs give another exact inequality:
\[
 Q(d)>c(\rho)\max_{\sigma(x)\in\{-1,1\}}Q(\sigma d)
                         \ge \tfrac89c(\rho)\gamma(a\ne b)^2 .
                                                               \tag{PR4}
\]
The maximum exists: relax the signs to the weak-star
compact unit ball of $L^\infty$, and use compactness
of $T_\rho$ on $L^2$ in fixed dimension. A fractional
maximizer on a positive-measure part of the disagreement
set would allow two opposite small perturbations with
strictly larger average quadratic energy, by strict
positivity of $Q$. For the last bound, take the direction
of the heaviest unordered disagreement type. Select
$\sigma$ to align the projection of $d$ with this
direction. The projection of $\mathbb E(\sigma d)$ is
at least $(\gamma(a\ne b)+\gamma(E_e))/\sqrt2$.
The constant-mode contribution to $Q(\sigma d)$ gives
the claimed bound.

### Quantitative gain from moving the selector thresholds

The slack in (PR2) controls a variation of the midpoint
scores, uniformly over the selector directions. Orient each
nonnull disagreement type $e=(i,j)$, with $d_e=e_i-e_j$,
and put $q_e=T_\rho(u_i-u_j)$ on $E_e$. If $\mathcal R(v)$
denotes the right side of (PR2), then
\[
 Q(d)-\mathcal R(v)\ \ge\
 c_\rho\min_{\lambda\in\mathbb R^3}
                \sum_e\int_{E_e}(q_e-d_e\cdot\lambda)^2d\gamma
                                                               \tag{PR5}
\]
for a positive constant depending only on correlation.
The scalar-selector versions hold for each sign profile,
hence also at the maximum in (PR4).

Take the least-squares residual $f_e=q_e-d_e\cdot\lambda$.
The normal equations give $\sum_e d_e\int_{E_e}f_e=0$ and
$\sum_e\int q_e f_e=\sum_e\int f_e^2$. The fitted
differences on nonnull types are bounded by two. In each
connected component of the graph of types, sort the fitted
label values. An adjacent gap exceeding one gives a cut
on which every cross-edge oriented residual has strictly
negative flux toward the high side, since $|q_e|\le1$.
This contradicts the normal equations. There are at most
three labels, so the range within a component is at most
two and $|f_e|\le3$.

On $E_e$ move the auxiliary halfspace threshold so the
conditional probability of label $i$ becomes
$1/2+\varepsilon f_e(x)$. Orient the normal to select
label $i$ on its positive side before this move. Volumes
remain exactly balanced by the residual flux equation.
The changed strip in the scalar normal coordinate has
Gaussian mass $|\varepsilon f_e|$ and width at most
$C|\varepsilon f_e|$ for small fixed $\varepsilon$.
At the central selector section, the tangent Gaussian
average of the smoothed scores equals $T_\rho u(x)$.
Indeed the selector term of the lifted partition is
odd in the entire auxiliary block $Z$, and this
oddness is preserved by Gaussian smoothing. The
central section is symmetric. The smoothed score
derivatives in a unit auxiliary direction are bounded
uniformly in dimension and location, by differentiation
of the conditional Gaussian in that direction.
The exact energy increment is twice the linear
interaction with the original smoothed score plus
a nonnegative squared Mehler norm. The former is at
least
$2\varepsilon\sum_e\int q_e f_e-C_\rho\varepsilon^2
\sum_e\int f_e^2$. Choose $\varepsilon$ small enough
depending only on $\rho$. Comparing the improved
selector to $s_*$ proves (PR5).

In particular a nonsimplex tied maximizer and the
simplex have strictly positive residual in (PR5) in
every orientation. If the residual were zero, then
a score $q_e$ would be constant on a positive-measure
type. Analyticity, balance, and Mehler injectivity
would give $u_i=u_j$ everywhere up to null sets.
For the midpoint of two ordinary label assignments,
this implies they share the third cell. Sharing a
simplex cell at a tie forces the simplex by (G15)–(G16).

### Strict loss even for arbitrary disagreement fibers

For two distinct balanced planar maximizers, any balanced
partition on a finite auxiliary block whose conditional
essential label support at almost every $x$ is exactly
$\{A(x),B(x)\}$ has strictly smaller stability. This
allows arbitrary conditional probabilities and arbitrary
measurable sections, with both labels used on disagreement
fibers.

First their midpoint $u$ cannot be cylindrical along a
base direction. Otherwise choose coordinates $(s,t)$
with $u$ independent of $t$. For almost every $s$ on a
disagreement type $\{i,j\}$, the two partitions choose
complementary labels for almost every $t$. Their
price-adjusted score contrasts $F_A^{ij}(s,t)$ and
$F_B^{ij}(s,t)$ sum to a function $C_{ij}(s)$ independent
of $t$, since the sum of their indicators is cylindrical.
On one complementary orientation the first contrast is
above $\max(0,C_{ij}(s))$ and on the other it is below
$\min(0,C_{ij}(s))$, using strict winning almost
everywhere. If both orientations occur on positive
measure in $t$, continuity forces $C_{ij}(s)=0$:
otherwise the first contrast passes through an open
interval of forbidden values. If orientation varies on a
positive-measure set of fibers of one type, analyticity
of $C_{ij}$ forces it identically zero. Its constant price
offset then vanishes by balance, and Mehler injectivity
gives $u_i=u_j$. For ordinary midpoint values this
forces the third cell to be common to both partitions
and cylindrical along the base direction, impossible.
If there is no positive-measure set of varying fibers,
both partitions are themselves cylindrical. This
also contradicts full active rank.

Now any invariant direction of the lifted partition
has zero base component. Indeed translating a section
in the auxiliary block preserves its essential label
support, since Gaussian and Lebesgue null sets agree
there. That support determines the midpoint $u$ even
when the conditional masses are not $1/2$. A combined
invariant translation would therefore make $u$ invariant
in its base component. The dimension reduction theorem
would require every auxiliary direction to be invariant
if a lift of two active base coordinates were optimal.
This is impossible because its sections on a
positive-measure disagreement set use two labels.

## Exact optimization when one cell is held fixed

Take any balanced partition $A$. Keep its original
label $k$ and, when the original label is $i$ or $j$,
choose between the two output labels $i,j$ using an
independent Gaussian block. The selector may be
different in the two rows. Among balanced rules in
this class, the optimal value is the larger of the
original value and the independent halfspace-split
value from (G13).

Let $J_{ij}=\Pr[I_A(X)=i,I_A(Y)=j]$. In each of the
two adjustable rows take the auxiliary set whose
membership selects output $i$. Their masses sum to
one. For fixed masses the two-set Gaussian comparison
maximizes the two within-row agreements and both
cross-row agreements by same-direction nested
halfspaces: an agreement term for sets of masses
$p,q$ is $1-p-q+2\Pr(Z\in E,Z'\in F)$. Thus one may
take the two auxiliary thresholds as $t,-t$, with
$t\ge0$ after an output swap if needed. Write
$a=\Phi(-t)$ and put
\[
 A(t)=a-\Pr(Z<-t,Z'<-t),\qquad
 B(t)=a-\Pr(Z<-t,Z'>t).
\]
Relative to the original arrangement the change is
$-2(J_{ii}+J_{jj})A(t)+4J_{ij}B(t)$. The endpoint
$t=\infty$ gives the original value and $t=0$ gives
the independent balanced halfspace split.

For $\sigma_0=\sqrt{1-\rho^2}$, put
$k_-=(1-\rho)/\sigma_0$ and $k_+=(1+\rho)/\sigma_0$.
Gaussian differentiation gives
$A'(t)=-\phi(t)(2\Phi(k_-t)-1)$ and the formula for
$B'(t)$ with $k_+$ in place of $k_-$. The ratio
$B'(t)/A'(t)$ is strictly decreasing. Indeed it is
the initial-interval weighted mean on $[0,t]$ of
$(k_+/k_-)\exp(-(k_+^2-k_-^2)u^2/2)$ with weight
$\phi(k_-u)$. Since $A(t)$ decreases, $B$ is a strictly
convex function of $A$. Thus the gain is maximized at
one endpoint, and at neither interior threshold it can
equal the larger endpoint value, as $J_{ij}>0$.
For a global maximizer the independent split is strictly
worse by (G13), so every nontrivial selector in this
fixed-cell class loses. The only equality cases are
the original two deterministic row assignments, up to
swapping their output labels.

## Every selector with at most two choices per input row

For a source with $s\ge S_\rho(T)$ for which all
independent fixed-cell splits are strictly worse, a
Gaussian selector using at most two output labels in
each input row loses strictly unless it deterministically
permutes the input labels. This includes every balanced
global optimum. To prove this, the matrix of conditional output
masses is doubly stochastic. A row with a single output
fixes an exclusive output column, reducing to the
fixed-cell comparison. If every row has two nonzero
options then each column has at least two, since an
exclusive column would have an entry equal to one.
Thus the support is a six-cycle and the matrix is a
mixture of two permutations whose relative permutation
is a three-cycle, with a common mixing mass $p$.
Relabel the outputs using the more frequent permutation,
so $p\le1/2$ and the baseline is the original partition.

Write $p=\Phi(-t)$, $t\ge0$, and let $D_i$ be the
auxiliary event that row $i$ selects its cyclic successor.
The joint probability of two copies of $D_i$ is at most
$Q_s=\Pr(Z>t,Z'>t)$ by the two-set Gaussian comparison.
For two different rows the joint event has probability
at least $Q_o=\Pr(Z>t,Z'<-t)$ by applying the same
comparison to one complement. These bounds can be
used separately even when they are not jointly attained.
Within a row the agreement is $1-2p+2\Pr(D_i,D_i)$.
Across two different rows the agreement is $p-\Pr(D_i,D_j)$,
since the two rows have only one output in common.
Therefore the gain over a baseline value $s$ is at
most
\[
          -(3s-1)p+2sQ_s-(1-s)Q_o .
                                                               \tag{PR6}
\]
The right side is the integral against $\phi(z)\,dz$
on $z\ge t$ of
\[
 \psi(z)=-(3s-1)+4s\Phi(-k_- z)-2(1-s)\Phi(-k_+z),
\]
using differentiation of the Gaussian two-tail probabilities
as in the previous section. At $t=0$, this integral
divided by $p=1/2$ equals $-2s+(1+s)h_\rho$, where
$h_\rho=1-\arccos(\rho)/\pi$. It is strictly negative
for $s\ge S_\rho(T)$. For $0<\rho\le1/2$, put
$q=\arcsin(\rho)/\pi\le\rho/3$. Then
$h_\rho/(2-h_\rho)-1/3=4q/[3(3/2-q)]\le \rho/3$,
whereas $S_\rho(T)-1/3\ge(9/(8\pi))\rho>\rho/3$.
For $\rho\ge1/2$ put $b=\arccos(\rho)/\pi\le1/3$.
The angle displacement $\Delta\in[-\pi,\pi]$ of
the correlated planar Gaussian pair satisfies
$\mathbb E|\Delta|=\arccos(\rho)$, by computing a
halfspace mismatch conditional on its displacement
and a uniformly rotated starting angle. For three
equal sectors the conditional disagreement is
$\min(1,3|\Delta|/(2\pi))$, strictly below the linear
expression on a positive-measure set. Hence
$S_\rho(T)>1-3b/2\ge(1-b)/(1+b)=h_\rho/(2-h_\rho)$.

Finally $\psi(0)=0$, $\psi(\infty)=-(3s-1)<0$, and
$\psi'$ changes sign at most once, from positive to
negative. This follows by dividing its positive
Gaussian term by its negative one, since $k_+>k_-$.
Thus $\psi$ is either negative throughout $z>0$ or
nonnegative on an initial interval and negative
afterward. In either case the strict negativity of
the integral from zero implies strict negativity
of the integral from every finite $t\ge0$.
This proves strict loss for a genuine cyclic mixture.
Any improving product rule must consequently use
all three output labels in at least one row for
every choice of the input block isolated as the source.

## Transition matrix restrictions from explicit product rules

The transition matrix of a balanced Gaussian partition
is $P_{ij}=3J_{ij}$, symmetric stochastic and positive
semidefinite in uniform $L^2$ on the labels. Write
$P|_{\mathbf1^\perp}=\lambda I+E$, with
$\operatorname{tr}E=0$ and $\xi=\|E\|_{\rm op}$.
Thus $\lambda=(3S_\rho(A)-1)/2$ and the two nonconstant
eigenvalues are $\lambda\pm\xi$.

Fix a neutral input label $k$ and code the other two
as $s=1,-1$, with $s(k)=0$. Take two independent copies
and output the sign of their sum, using $k$ on a zero.
This is balanced, with three of the nine input
pairs for each output. In the uniform centered ON
basis $u=\sqrt{3/2}s$, $v=(1,1,-2)/\sqrt2$, write
$P|_{\mathbf1^\perp}=\bigl(\begin{smallmatrix}a&c\\c&b\end{smallmatrix}\bigr)$.
The output basis functions expanded on the input blocks are
\[
 U=\tfrac23(u_1+u_2)-\tfrac{u_1v_2+v_1u_2}{3\sqrt2},
 \qquad V=\tfrac{u_1u_2-v_1v_2}{\sqrt2}.
\]
Choose the second input copy with the two nonneutral
labels exchanged. Its matrix then has off-diagonal
$-c$, while $a,b$ are unchanged. The centered
combined energy $\langle U,(P_1\otimes P_2)U\rangle+
\langle V,(P_1\otimes P_2)V\rangle$ is
\[
 \tfrac89a+\tfrac19ab+\tfrac12(a^2+b^2)+\tfrac89 c^2 .
\]
This follows from the orthogonal constant, one-block
and two-block terms. The noise stability is one third
of one plus this energy. For an actual global maximizer
the centered energy is strictly below $2\lambda$ by
genuine two-block dependence. With
$x_k=(a-b)/2$ this becomes the necessary condition
\[
           x_k+\xi^2<\tfrac54\lambda(1-\lambda)
                                  \qquad(k=1,2,3).         \tag{PR7}
\]
Indeed the terms in $x_k^2$ cancel using
$(a-b)^2/4+c^2=\xi^2$. The three neutral directions
give the three samples of a cosine at equally spaced
angles after doubling the angle of each basis axis.
Thus $\max_k x_k\ge \xi/2$, yielding also
$\xi^2+\xi/2<\tfrac54\lambda(1-\lambda)$.

### A sharp Gaussian bound on the largest nonconstant eigenvalue

There is an independent cap from the two-set Gaussian
comparison. Put $t_0=\Phi^{-1}(2/3)$,
$A_0=\Pr(X>t_0,Y>t_0)$ and
$C_0=\Pr(X>t_0,Y<-t_0)$ for a scalar Gaussian pair.
Every balanced Gaussian three-label transition matrix satisfies
\[
                 \lambda+\xi\le 3(A_0-C_0).                \tag{PR8}
\]
Equality forces the balanced parallel three-slab
partition, up to a rotation and label permutation.
In particular actual global maximizers satisfy (PR8)
strictly.

For any fixed triple of scalar values, sort them and
write the resulting function of the labels as its minimum
plus a nonnegative sum of its two superlevel indicators.
After expanding the correlated product, the terms
involving its minimum are fixed by the masses.
Every superlevel pairing has a nonnegative coefficient.
The two-set Gaussian comparison bounds all of them
simultaneously by nested parallel halfspaces of
masses $1/3$ and $2/3$. Thus the scalar Rayleigh quotient
is bounded by that for the same values on the balanced
three-slab partition. The odd and centered even
eigenvalues of that transition matrix are
$3(A_0-C_0)$ and $-2+9(A_0+C_0)$, respectively, from
reflection symmetry and the row-sum identities.
Their difference is $2-6A_0-12C_0$.
Differentiating the joint Gaussian density and integrating
on the quadrants shows its derivative has the sign of
$2-\exp(2\rho t_0^2/(1-\rho^2))$. The difference
vanishes at zero and tends to zero at one, so it
is positive in between. Accordingly equality in the
largest eigenvalue bound requires the odd three-slab
values, which are distinct. Both nontrivial superlevel
pairings must attain the halfspace comparison,
and the cross pairing forces the orientation to
agree. This gives exactly the three-slab partition,
which is one-dimensional and hence is not an
actual balanced maximizer.

## Every two-copy lookup table loses unless it is a dictator

The two-copy comparison also follows from the transition
matrices alone under weaker hypotheses than global
maximality. Take any two balanced Gaussian input
partitions of common stability $s\ge S_\rho(T)$.
Assume that for each input and each fixed cell its
independent halfspace-split comparison (G13) has
strictly lower stability than $s$. Then only the
permutations of one input attain $s$ among balanced
two-copy lookup tables, with every other table strictly
below $s$. Actual global maximizers satisfy these
hypotheses, so the comparison recovers their strict
two-copy result by matrix inequalities.

Let both input channel traces correspond to $\lambda$
and set
$R=3(A_0-C_0)$ from (PR8), $\delta=R-\lambda$.
The two centered channel deviations have norms at
most $\delta$. Two bounds valid throughout the interval
are
\[
                 \delta<\lambda/2,\qquad
                 \delta<\tfrac9{20}(1-\lambda).              \tag{PR9}
\]
For the first, the odd slab Rayleigh quotient $R$ has
first-degree weight $\eta=6\phi(t_0)^2<4/5$ and
the remaining weights are of odd degrees at least three.
Thus $R/\rho\le 4/5+\rho^2/5$. The first three degrees
in the closed simplex value (D0) give
\[
 \tfrac{3\lambda}{2\rho}\ \ge\
 \tfrac94\left(\tfrac9{8\pi}+
                     \tfrac9{16\pi^2}\rho+
                     \tfrac9{64\pi}\rho^2\right)
  > .805+.128\rho+.100\rho^2 .
\]
This is larger than $4/5+\rho^2/5$.
For the second put $\theta=\arccos\rho$. Integrating one
positive summand of the exact derivative of $R$ gives
$1-R\ge (3/(2\pi))e^{-t_0^2}\theta$. The sector angle
bound used for (PR6) gives $1-\lambda\le9\theta/(4\pi)$.
Since $(2/3)e^{-t_0^2}>11/20$, the second bound follows.
For these two rational estimates on $\phi(t_0)$ and
$e^{-t_0^2}$, one can use $.430<t_0<.433$: integrate
$1-t^2/2\le e^{-t^2/2}\le1-t^2/2+t^4/8$ between zero
and the corresponding endpoints, and use
$.398<\phi(0)<.399$. The alternating Taylor bounds on
$e^{-t_0^2}$ then give the stated inequalities.

It suffices to check tables having a rainbow row and
a rainbow column, by the strict two-choice comparison.
Normalize the rainbow row and column to the first row
and column, both $(0,1,2)$. The remaining $2\times2$
word is a permutation of $(0,0,1,2)$. The groups in
the following table are under row, column and output
permutations and interchange of the two input blocks.
\[
\begin{array}{c|c|c}
\text{interior words} &\text{one representative, by rows}&
         \text{first-level weight}\\ \hline
1002 &020,\ 211,\ 012&8/27\\
1200,1020,0102,0012 &010,\ 212,\ 021&8/27\\
2100,0120,2010,0210,0201,0021 &001,\ 212,\ 120&4/27\\
2001 &012,\ 120,\ 201&0
\end{array}
\]
To check the grouping, transpose the interior word or
reverse it while swapping 1 and 2, which both preserve
the normalized first row and column. They reduce the
twelve words to $1002,1200,2100,0120,2001$.
Row, column and output permutations send the two
tables for $2100,0120$ to the third representative.
For the Latin-square representative encode the output
by $\omega^{X_1+X_2}$, $\omega=e^{2\pi i/3}$. For
any symmetric input channel the energy of $\omega^X$
is $\lambda$ by the trace identity on the centered
space. Hence this table has stability
$(1+2\lambda^2)/3<s$.

For the other three, write $W_1,W_2$ for the two
nonconstant product level weights of the indicator
vector. In the centered basis $u,v$ used for (PR7),
let $B_j^{(1)}$ be the coefficient Gram matrix of
level one on block $j$, and $B_j^{(2)}$ the Gram
of level two after summing over the other block.
If $\operatorname{spr}$ means the difference of the
two eigenvalues of a symmetric matrix, the centered
energy is at most
\[
 \lambda W_1+\lambda^2 W_2
     +\delta\sum_{j=1}^2\operatorname{spr}
                          (B_j^{(1)}+\lambda B_j^{(2)})
     +\delta^2 W_2.                                         \tag{PR10}
\]
Indeed $E_j$ is traceless with norm at most $\delta$,
so each linear term is bounded by the eigenvalue
spread. The product perturbation has norm at most
$\delta^2$.

For explicit verification the following entries are
the matrices multiplied by 27. Each is obtained by
taking dot products of
$\mu_{\alpha\beta}=\mathbb E[e_{f(X_1,X_2)}
\alpha(X_1)\beta(X_2)]$, $\alpha,\beta\in\{1,u,v\}$,
using uniform measure on the nine entries.
\[
\begin{array}{c|c|c|c|c}
\text{row of preceding table}&27B_1^{(1)}&27B_2^{(1)}&
                         27B_1^{(2)}&27B_2^{(2)}\\ \hline
1&\mathrm{diag}(4,0)&\mathrm{diag}(4,0)&5I&5I\\
2&\mathrm{diag}(4,0)&\mathrm{diag}(4,0)&
 \begin{pmatrix}2&-\sqrt3\\-\sqrt3&8\end{pmatrix}&
 \begin{pmatrix}5&-2\sqrt3\\-2\sqrt3&5\end{pmatrix}\\
3&\mathrm{diag}(4,0)&0&\mathrm{diag}(5,9)&
 \begin{pmatrix}6&-\sqrt3\\-\sqrt3&8\end{pmatrix}
\end{array}
\]
The two-coordinate weights are $10/27,10/27,14/27$.
Multiplying (PR10) by 27, the baseline loss relative
to a dictator is $10\lambda(1-\lambda)$ for the
first two tables and $14\lambda(1-\lambda)$ for the
third. The coefficient of $\delta$ is 8 for the
first, $4(\sqrt{1-3\lambda+3\lambda^2}+
\sqrt{1+3\lambda^2})$ for the second, and 4 for
the third. The second coefficient is at most 8
for $\lambda\le1/2$ and at most $16\lambda$ above
that point, by Cauchy–Schwarz. For the first two
tables, $\lambda\le12/25$ gives
$8\delta+10\delta^2<4\lambda+(5/2)\lambda^2
\le10\lambda(1-\lambda)$. For $12/25\le\lambda\le1/2$,
use the second bound in (PR9), giving
$8\delta+10\delta^2<[3.6+2.025(1-\lambda)](1-\lambda)
<10\lambda(1-\lambda)$.
For $\lambda\ge1/2$ both are bounded by
$16\lambda\delta+10\delta^2
<[7.2\lambda+2.025(1-\lambda)](1-\lambda)
<10\lambda(1-\lambda)$.
For the third table,
$4\delta+14\delta^2<14\lambda(1-\lambda)$ follows
from $\delta<\lambda/2$ for $\lambda\le1/2$ and
from $\delta<(9/20)(1-\lambda)$ for $\lambda\ge1/2$.
This proves the strict comparison for all the tables.

## A neighborhood in which all finite product tests lose

The transition restrictions from finite label
recombinations have a limitation. Fix $0<\lambda<1$.
There is $\delta>0$ depending only on $\lambda$ such that
for every finite collection of symmetric stochastic
three-label channels $P_j$ satisfying
\[
 \operatorname{tr}(P_j|_{\mathbf1^\perp})=2\lambda,\qquad
 \|P_j|_{\mathbf1^\perp}-\lambda I\|_{\rm op}<\delta ,
\]
the best balanced output rule on their independent
product is a permutation of a single input label.
The assertion is purely finite-channel and allows
different $P_j$ across the blocks and arbitrarily many
blocks. The matrices can be taken nonisotropic.

Here are the details. The product decomposition by the
set of nonconstant coordinates is orthogonal in
uniform measure. On a nonempty set of $n$
coordinates the difference from the isotropic channel
block $\lambda^n I$ has operator norm at most
$n\delta(\lambda+\delta)^{n-1}$ by telescoping.
This is at most $C_\lambda\delta$ uniformly in $n$
when $\delta$ is small. Encode a balanced three-valued
rule by $F=\omega^{f}$ with $\omega=e^{2\pi i/3}$.
It has mean zero and norm one; its correlated energy
is $(3S(f)-1)/2$. At the isotropic product, energy is
at most $\lambda-\lambda(1-\lambda)\varepsilon$,
where $\varepsilon$ is the squared weight above one
active coordinate in the product decomposition. If
the rule reaches the dictator value $\lambda$, the
operator perturbation bound gives $\varepsilon\le C_\lambda\delta$.

For clarity, small $\varepsilon$ forces proximity to
a permutation of a single input by an elementary argument.
Write $g=\sum_j g_j$, $g_j=\mathbb E[F\mid X_j]$,
and $u_j=\|g_j\|_2^2$. The $g_j$ are centered,
independent and bounded by one. The moment
$\mathbb E|g|^6$ is bounded uniformly because in its
expansion each nonzero term uses at most three indices,
each repeated at least twice, and $\sum_j u_j\le1$.
$\|g-F\|_2^2=\varepsilon$ implies
$\|g^2-F^2\|_2^2\le C\sqrt\varepsilon$ by Cauchy–Schwarz
and the sixth moment bound. Since $F^2=\bar F$, its
degree-two block weight is at most $\varepsilon$.
The two-coordinate part of $g^2$ is
$2\sum_{i<j}g_i g_j$, of squared norm
$4\sum_{i<j}u_i u_j$. Hence
$\sum_{i<j}u_i u_j\le C'\sqrt\varepsilon$.
Since $\sum_j u_j=1-\varepsilon$, some
$u_j\ge1-C''\sqrt\varepsilon$. Rounding $g_j$ to the
nearest cube root changes the original rule on mass
at most $C'''\sqrt\varepsilon$. The rounded rule must
be a permutation when this mass is less than $1/3$,
because otherwise it misses a balanced output label.

Finally, near any single-input permutation $D$ a
uniform loss estimate holds. Let $m=\Pr(f\ne D)$, and
let $m_{ij}=\Pr(D=i,f=j)$ for $i\ne j$. The outgoing
and incoming masses agree at each label. Decompose
this directed flow into cycles of lengths two and three.
The linear score cost from each transferred label is
$P_{ii}-P_{ij}$ for the input channel used by $D$
with its labels permuted. The cost per unit mass on
a two-cycle is
$\frac12(e_i-e_j)^{\mathsf T}P(e_i-e_j)\ge\lambda_-$,
the smallest nonconstant eigenvalue of that channel.
For a three-cycle the average cost is $\lambda$, by
symmetry and row sums. Thus the total linear regret
is at least $\lambda_- m$.

The quadratic transfer energy is at most $2m^{2/p}$
for some $p<2$ independent of the number of channels.
Here is a finite-space proof of that bound. Choose
$a$ with $\lambda<a<1$ and $\delta$ small. In
positive-semidefinite order every $P_j$ is bounded
by the isotropic channel $J+a(I-J)$, where $J$
is averaging. They are positive semidefinite for
$\delta<\lambda$. Tensor products preserve this
order by telescoping. The single-site isotropic
channel $J+\sqrt a(I-J)$ maps $L^p$ to $L^2$
with norm at most one for some $1<p<2$. At $p=2$
the inequality is strict off constants, so away
from constants on the unit sphere it persists for
$p<2$ sufficiently close to two. Near a positive
constant, the quadratic coefficient of the squared
$L^p$ norm in centered directions is $p-1$, versus
$a$ for the squared output $L^2$ norm. Choose $p-1>a$
and use the uniform Taylor remainder in this finite
dimension. Positivity allows reduction to nonnegative
inputs. Minkowski's inequality tensorizes the bound.
Applying the result to the three transfer indicators
and using that their total $L^1$ norm is $2m$ and
each has $L^1$ norm at most $m$ gives the asserted
energy bound.

The exact loss from $D$ is twice the regret less the
transfer energy, hence at least
$2m(\lambda_- -m^{2/p-1})$. For the rounding distance
above and $\delta$ small, this is strictly positive
if $m>0$, since $\lambda_-\ge\lambda-\delta$.
Thus every rule attaining the optimum is a dictator
or its output permutation, as claimed. This leaves
an open family of transition matrices for which
even unbounded finite product tests give no strict
improvement over a single input.
