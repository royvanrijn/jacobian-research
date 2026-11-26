# Private prime support on differences of actual shared arrays

Pairs of shared-value boxes can be matched on some or all of
their residue indices. In a large family of such pairs, each
edge of the divided difference array has prime support of
fixed positive radical weight that appears on no other edge.
A Latin family of cells transfers these disjoint packets to
primitive abc directions. The number of transferred packets
is quadratic in the box side, even though selecting an array
requires only a spanning tree.

## Pairs with prescribed residue matches

Use one three-packet tuple from (SB2), with \(c_0=1/6\),
prices \(s_F,s_G,s_M>0\), \(d=s_F+s_G+s_M\), and
\(t=c_0-d>0\). Fix positive lower bounds for the prices
and for \(t\). Write \(n=X^L\), and let \(R_i\) denote the
residue alphabet caps and \(K\) the complete-value cap
per opposite residue pair in (LB1). We may take
\(R_i=n^{s_i+O(\kappa+\eta)}\) and \(K=n^{t+O(\eta)}\),
using the explicit upper bounds there and rounding up.
A letter in an original row has radical at most \(n^C\)
for a fixed \(C\). All packet and radical assumptions
are on the actual rows supplying the arrays.
[Simultaneous packet cells (SB2)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#simultaneous-marking-theorem)
· [Actual shared arrays (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)

Let \(\mathcal G\) be the complete tripartite graph on
the \(3m\) indexed potential vertices, with \(m\ge2\)
fixed. Each vertex has an axis \(i\in\{F,G,M\}\) and
carries its corresponding residue label. Prescribe a
subset \(S\) of these vertices. Consider ordered pairs of
actual boxes whose ordered residue labels agree at every
vertex outside \(S\). Allow the two labels at a vertex
in \(S\) to differ. Put
\[
 R=R_F R_G R_M,\qquad
 R_S=\prod_{v\in S}R_{\operatorname{axis}(v)},\qquad
 r_{\rm tree}=md+(3m-1)t .
\]
Cauchy–Schwarz over the matching keys, of which there
are at most \(R^m/R_S\), gives a pair count
\[
             \gg_m R^m R_S K^{6m-2}
                           n^{-O_m(\kappa+\eta)}.              \tag{PA1}
\]
Here the matching key retains the ordered residues outside
\(S\). Indeed (LB1) gives total array count at least
\(R^mK^{3m-1}n^{-O_m(\kappa+\eta)}\) up to a constant,
and one squares that count and divides by the key count.
One can require unequal residues at every vertex in \(S\)
at a power-small loss. For each prohibited equality,
the number of complete pairs of residue lists drops by
a factor \(R_i\) for that axis, and at fixed complete
lists the two spanning trees cost at most \(K^{6m-2}\).
When \(S=\varnothing\), (PA1) is the synchronized-pair
bound in (CS3).
[Pairs of synchronized boxes (CS3)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#fresh-direction-radicals-on-synchronized-actual-boxes)

For a pair in (PA1), write \(L_e^{(1)},L_e^{(2)}\) for its
two positive labels on an edge \(e\in E(\mathcal G)\), and
put \(\Delta_e=L_e^{(2)}-L_e^{(1)}\).
The \(\Delta_e\) are signed potential differences on
\(\mathcal G\), with edge orientation fixed by the
three roles. For a prescribed initial prime pool \(P_0\),
put
\[
 P_* = P_0 \cup\bigcup_{e\in E(\mathcal G)}
                     \operatorname{supp}(L_e^{(1)}L_e^{(2)}).
\]
For each edge set
\[
 s(e)=\max\bigl(\{0\}\cup
           \{s_{\operatorname{axis}(v)}:v\in e\cap S\}\bigr).
\]
There is no extra residue exponent at an edge whose
two endpoint residues are matched.

Fix a Latin family of \(m^2\) cells that partitions the
edges into triangles. For example, number each axis
modulo \(m\) and use the cells \(k=i+j\bmod m\).
When its three differences are nonzero, a triangle gives
a primitive abc direction by dividing their common gcd,
taking absolute values and putting the largest last.
Denote the collection of these \(m^2\) outputs by
\(\mathscr D\).

## Private packets and bounded-height completions

**PA2.** Fix \(0<\sigma<2t\) by a positive margin, and
put \(b_e=2(t+s(e))-\sigma\). Outside a power-small
fraction of (PA1), all the differences are nonzero and
all \(m^2\) directions in \(\mathscr D\) are tall, with
height at least \(X^{\theta}\) for a uniform fixed
\(\theta>0\).

On the branch where the raw radicals of all directions
in \(\mathscr D\) are at most \(n^{B_0}\), for any
prescribed fixed \(B_0\), there is an additional
power-small exceptional set outside of which every edge
satisfies
\[
  \prod_{\substack{p\mid\Delta_e\\p\notin P_*\\
                   p\nmid\Delta_f\ (f\ne e)}}p>n^{b_e}.
                                                               \tag{PA2}
\]
Here \(P_0\) may have any prescribed polynomial height.
The exceptional count is relative to all of (PA1).
All packet tolerances and inner exponent allowances
can be chosen after the fixed margins and before the
height limit at small fixed budget.

### The fixed residue lists and the geometric exceptions

At fixed full residue lists for both arrays, fix the
first array and all but one of the second potentials
up to common translation. The second difference graph
after deleting the one vertex is connected. Its spanning
tree has \(3m-2\) labels, with at most \(K\) choices each.
If two second-minus-first potential values are equal,
omit one of their vertices. The equality determines
the omitted potential relative to a known potential.
The recovered potential may fail the residue constraint,
which can only reduce the count. Thus each equality
costs at most \(K^{6m-3}\) pairs per fixed full residue
lists. The test is well-defined under separate common
translations of both arrays since equality of two
potential differences means equality of their difference
across the two arrays. There are only \(O_m(1)\) tests.
Summation over the at most \(R^mR_S\) pairs of full
residue lists loses a factor \(K\) compared with (PA1),
up to tolerance losses. This proves that equal divided
potential differences are power-sparse, without any
raw direction cap. In particular every edge difference
is nonzero, and any two potential differences that
will appear as distinct neighbors in a later test are
distinct.

Also exclude coincidences between the zeros of a moving
raw difference and a moving second original label.
Write \(z_v^{(1)},z_v^{(2)}\) for the two sets of potentials.
When \(x\) is omitted and its variable is measured as a
second-minus-first potential, these zeros for neighbors
\(y,z\) are \(z_y^{(2)}-z_y^{(1)}\) and
\(z_z^{(2)}-z_x^{(1)}\), respectively. Equality requires
\[
                 z_y^{(2)}-z_z^{(2)}
                             =z_y^{(1)}-z_x^{(1)}.
\]
For \(y=z\) this is impossible because the first label
on \(xy\) is positive. For \(y\ne z\), omit \(y\) in
the second tree count, keeping \(z\). The displayed
relation determines the missing potential, so loses
the same factor \(K\) at fixed full residue lists.
A union over all these finitely many tests again costs
a power-small fraction of (PA1).

For any fixed pair of full residue lists, the two
corresponding cells lie in a union with at most \(2K\)
complete values per role and at most \(2K^2\) rows.
The short-pair bound of (CS2) in this union has at most
\(O(K^4n^{-\delta_0})\) bad pairs for some fixed
\(\delta_0>0\). One can use, for example,
\(\alpha<\frac14\liminf t\) there, then choose the
allowances sufficiently small. An identical pair is
included at cost \(O(K^2)\). Each fixed pair of
original rows fixes two spanning-tree labels in each
array, leaving at most \(K^{6m-6}\) completions per
full pair of residue lists. Sum over the finitely
many Latin cells and all lists. This excludes short
directions at a power saving. It also excludes any
zero entry of a direction, which would give height
one after dividing the common gcd.
[Short secants under the three original radical caps (CS1–CS2)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#the-secant-inequality-at-the-cell-scale)

### The omitted-vertex test

Fix a test edge \(e\), and omit an endpoint \(x\) attaining
\(s(e)\) when \(s(e)>0\), or either endpoint otherwise.
Fix the entire first array and all the other potentials
of the second array using the connected tree just
described. If \(x\notin S\), its required second
residue is known from the first array's list. If
\(x\in S\), do not specify its second residue yet.
There are at most \(R^mR_S\) pairs of lists to use in
the first case and at most \(R^mR_S/R_i\) in the
second, where \(i=\operatorname{axis}(x)\). The
omitted numerical potential has at most \(K\) or
\(K R_i\) possibilities, respectively: after its
residue is supplied, one edge to a known potential
determines it and has at most \(K\) labels in that
pair class. This is where the factor \(s(e)\) is
available, and only at an unmatched endpoint.

Use \(\xi\) for the variable second-minus-first
potential at the omitted vertex and \(\gamma_y\)
for the fixed potential differences at its neighbors.
The corresponding signed moving differences are
\(\pm(\xi-\gamma_y)\). Let \(\lambda_z=z_z^{(2)}-z_x^{(1)}\)
be the zero of a moving second original edge label.
Take a pool \(P_e\) containing \(P_0\), all first
original labels, all fixed second original labels,
and all fixed nonzero edge differences. Also include
the prime supports of \(\gamma_y-\gamma_z\) for
distinct known neighbors and of \(\gamma_y-\lambda_z\)
for every such difference and original label zero.
These integers are nonzero after the geometric
discard and have polynomial height uniformly. In
particular \(P_e\) is fixed before the omitted value
is chosen. A prime of the test form outside \(P_e\)
divides no fixed label or difference there. It
divides no other moving difference, and divides
no moving original label of the second array, by
the respective zero-difference tests. Thus it is
outside \(P_*\) and private among all the raw edge
differences.

The prescribed Latin triangle containing \(e\) joins
the omitted vertex to two known potentials. Their
fixed label difference \(H\) between the two arrays
is nonzero. The common gcd of its three differences
divides \(H\) and is supported on \(P_e\).
Consequently the outside-\(P_e\) radical of the
moving entry of its primitive direction is exactly
the outside-\(P_e\) radical of the tested \(\Delta_e\).
The direction's fixed entry has absolute value
dividing \(|H|\), at divisor cost. Fix that entry,
its position among the three absolute values, and
one dyadic height shell for the direction. A signed
direction recovers the gcd from \(H\) and its signed
fixed entry. It then determines the omitted potential
from one of the moving labels. There are only
boundedly many sign and position choices.

On the prescribed raw direction branch, (FL1) in the
derived shell bounds the completions with
outside-\(P_e\) radical at most \(n^{b_e}\) on the test
entry by \(O(n^{b_e/2+\xi_0}X^\varepsilon)\), with any
fixed small normalized allowance \(\xi_0>0\).
Indeed the shell height \(Y\) is at least \(X^\theta\),
so its raw caps are \(Y^{O_{B_0,\theta}(L)}\), the
fixed pool has polynomial height in \(Y\), and the
outside cap on the one free entry is
\(Y^{Lb_e\log X/\log Y}\). The second free cap can be
the whole raw cap. The height conversion cancels in
the main exponent. The signs, divisor records, shell
choices and inner allowances can be absorbed in
\(X^\varepsilon\). The bad completions cost a fixed
power less than \(K\), or \(KR_i\) if the omitted
residue is free, since \(b_e/2=t+s(e)-\sigma/2\).
Multiply by the corresponding pair-of-lists count
and by the \(K^{6m-3}\) choices from the known trees.
In the free-residue case the saved list factor
\(R_i^{-1}\) cancels the unrestricted completion factor
\(R_i\). Thus in both cases the total exceptions are
\(O(R^mR_SK^{6m-2}n^{-\delta_1})\) for some fixed
\(\delta_1>0\). Sum over the edges and use (PA1).
The enlarged pool test implies (PA2) because its
surviving primes are outside all the other differences.
[Prescribed-pool fixed-letter estimate (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Successive directions at seed vertices (CS8)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#a-seed-triangle-and-successive-one-vertex-directions)

## Joint output, relative budgets and common fields

Write \(B_*=\sum_{e\in E(\mathcal G)} b_e\), so
\[
              B_* = 6m^2t+2\sum_e s(e)-3m^2\sigma .
\]
The packets counted by (PA2) are pairwise disjoint.
Any prime private to a raw edge difference survives
in the corresponding entry of every primitive triangle
containing that edge, since it does not divide either
of the other two differences. In particular, the
radical of the union of all \(m^2\) Latin directions
outside \(P_*\) is greater than \(n^{B_*}\).

That joint floor is unconditional if the initial
pool has prime product at most \(n^{B_1}\) for a fixed
\(B_1\), including the empty pool. If the joint
outside radical were at most \(n^{B_*}\), each of the
Latin direction radicals would be at most
\(n^{B_*+B_1+6m^2C}\), using the original complete
radical caps in both arrays. Choose \(B_0\) above
this exponent uniformly in the prescribed parameters.
The capped-branch proof then applies to any such
violation. In exact synchronization \(S=\varnothing\),
the joint exponent approaches \(6m^2t\), whereas the
successive seed directions in (CS8) give a linear
exponent in \(m\).

There is also an individual relative-budget test.
Write \(D=FGM\), \(w=\log_X D\), and
\(w_i=\log_X D_i\) for the packet on axis \(i\).
For a Latin triangle \(T'\), each of its complete
differences is divisible by
\[
                 D\Big/\prod_{v\in T'\cap S}
                                 D_{\operatorname{axis}(v)}.
\]
Indeed an edge's own-role packet always divides both
original values. Opposite moduli are preserved at
the vertices whose residues were matched. Hence the
direction height is at most
\(X^{1-w+\sum_{v\in T'\cap S}w_{\operatorname{axis}(v)}}\).
For a fixed \(\chi>0\), the condition that a specified
Latin direction satisfies an abc budget \(\chi L\)
is power-exceptional as soon as
\[
 \chi\left(1-w+\sum_{v\in T'\cap S}w_{\operatorname{axis}(v)}\right)
                          <\sum_{e\in T'}2(t+s(e))
\]
by a fixed margin. In this test it suffices to impose
the raw cap on this one direction, because its three
edge tests use only that triangle's derived radical
cap. In exact synchronization, if \(w>6d\) by a fixed
margin, almost all pairs have no Latin direction
satisfying the original relative budget \(L\).
Thus a difference iteration preserving that budget
cannot use generic synchronized arrays in that sector.

Finally the first-tree support argument of (CS6)
applies to the matching pattern in (PA1): a given
first array has at most \(R_S K^{3m-1}\) second arrays
with the prescribed matches. Thus for every fixed
\(\tau>0\), almost all pairs have first tree radical
greater than \(n^{r_{\rm tree}-\tau}\) outside the
initial pool, simultaneously for all first trees.

For a pair of original spanning trees, one in each
array, put
\[
 r_{\rm pair}=md+\sum_{v\in S}s_{\operatorname{axis}(v)}
                     +(6m-2)t.
\]
The same support enumeration applied jointly gives
a larger generic floor \(n^{r_{\rm pair}-\tau}\)
for the radical of the product of both trees outside
\(P_0\). Indeed the outside support at cost at most
\(n^{r_{\rm pair}-\tau}\), together with \(P_0\), permits
only subpower many choices of the complete values on
both trees. The two trees reconstruct both arrays,
including their residue lists. This is power smaller
than (PA1), and there are only boundedly many tree
pairs. With either choice, one or both trees, the
Latin directions contribute radical greater than
\(n^{B_*}\) outside the original tree support since
their packets are outside both arrays.

Suppose those chosen tree labels, with height caps
\(X^{h_e}\), and all entries of the \(m^2\) Latin
directions become \(Q\)-full principal ideals in one
field \(F_{\rm lift}\), with integer \(Q\ge1\).
Every indexed vertex appears in \(m\) Latin triangles,
so the sum of the output entry height exponents is
at most \(3m^2(1-w)+3m\sum_{v\in S}w_{\operatorname{axis}(v)}\).
Equation (RL1) therefore gives, on these generic pairs,
\[
 \log_n\operatorname{rd}(F_{\rm lift})>
 r_*-\tau+B_*-
 \frac{\sum_e h_e+3m^2(1-w)+3m\sum_{v\in S}w_{\operatorname{axis}(v)}}{LQ},
                                                               \tag{PA3}
\]
where \(r_*=r_{\rm tree}\) for the first tree alone
and \(r_*=r_{\rm pair}\) for the pair of trees.
The sum of \(h_e\) is over the chosen original labels.
The field, its degree and the trees may vary between
pairs. The bound charges fresh support from the
difference outputs in addition to the original trees.
[Discriminant and normalized ideal radicals (RL1)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#the-local-cost-and-a-field-attaining-it)
· [Common-field tree count (RL8)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#a-common-field-on-an-actual-shared-value-box)


## Same-axis shifts in one actual array

A single shared-value array already controls the prime
supports of differences between two potentials on the
same axis. For every unordered pair of vertices \(x,y\)
in its potential graph put \(D_{xy}=|z_x-z_y|\).
All these values are positive. Opposite groups are
separated by the original positive edge labels, and
distinct residues within a group distinguish its
potentials. The cross-axis values are the original
complete labels.

**PA4.** Fix a sufficiently small \(\sigma>0\) and a
fixed raw cap exponent \(B_0\). Take the actual single
arrays counted by (LB1) from the surviving tuple, with
the prices and margins of (PA1). On all but a power-small
fraction of these arrays, each same-axis pair \(xy\)
in axis \(i\) for which \(\operatorname{rad}(D_{xy})\le
n^{B_0}\) satisfies
\[
         \prod_{\substack{p\mid D_{xy}\\p\notin P_0\\
                     p\nmid D_{uv}\ (\{u,v\}\ne\{x,y\})}}
                    p>n^{2(t+s_i)-\sigma}.                    \tag{PA4}
\]
Here the prescribed prime pool \(P_0\) may have polynomial
height. Only the raw radical of the tested shift is
capped. The other same-axis shifts can have arbitrary
radical. The primes counted in (PA4) are private among
all potential differences and outside all original
complete labels.

If the prime product of \(P_0\) is at most \(n^{B_1}\)
for fixed \(B_1\), the joint radical of the same-axis
shifts outside \(P_0\) and the original labels exceeds
\[
 n^{B_{\rm shift}},\qquad
 B_{\rm shift}=m(m-1)(3t+d)-3\binom m2\sigma
\]
on all but a power-small fraction of the arrays,
without any raw caps on those shifts. This unconditional
assertion concerns the union of their supports, including
primes shared by several shifts. The lower bound on the
union of their private supports uses the individual raw
caps as in (PA4).
[Forced actual-array count (LB1 and LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

### A quotient height bound on the original fibre

Fix \(\alpha>0\) and a bound \(n^C\) on each original
complete radical. There exist \(\theta_1,L_1>0\)
such that, for each fixed \(0<L<L_1\), the following
subfamily of any original fixed-complete-letter fibre
has at most \(O_{C,\alpha,L,\varepsilon}(n^\alpha
X^\varepsilon)\) rows. In that subfamily one free
letter \(U\) satisfies
\[
                U/\gcd(U,A)\le X^{\theta_1}
\]
for a fixed positive integer \(A\le X\) that may vary
with height and fibre. Use the original height shell
\(X/2<C_{\rm out}\le X\).

If this fails, take prospective cutoffs
\(\theta_j\to0\) and violating fixed budgets \(L_j\to0\).
For each \(j\), take arbitrarily large violating heights
before passing to the next budget. The selected fibres
have limiting normalized exponent at least \(\alpha\).
Use the prime support of \(A\) as the prescribed pool
on \(U\). Its entire supported part can be fixed at
height-subpower cost uniformly in \(A\); carry out the
original anatomic refinements in the inner height limit.
The outside part of \(U\) divides \(U/\gcd(U,A)\), so
its exact scaled valuation measure has first moment
at most \(\theta_j\) before taking weak limits.
For every fixed positive cutoff \(s\) its mass above
\(s\) is at most \(\theta_j/s\). The outer limiting
measure on the outside part is therefore supported
at zero. This contradicts (FT1), which gives it
nonatomic mass at least \(\alpha\) on this fixed-letter
fibre. This proves the uniform height cutoff with the
budget and height limits taken in order.
[Nonatomic mass on prescribed-pool fibres (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

### Counting a fixed label shift

For a test pair \(xy\) on axis \(i\), omit \(x\), including
its residue, and fix all the other potentials up to
translation. This costs at most \((R^m/R_i)K^{3m-2}\)
choices of residue lists and second-stage tree labels
in the single array. At the omitted vertex there
are at most \(K R_i\) possible completions before the
radical test.

Choose one known vertex \(z\) from another axis and
one known vertex \(z'\) from the remaining axis.
The current original row in cell \((x,z,z')\) has
fixed opposite complete label \(T=L_{zz'}\).
Its moving label \(U=L_{xz}\) is one of its two free
letters. Put \(A=L_{yz}\), a fixed original complete
label, with \(|U-A|=D_{xy}\), since \(x,y\) belong
to the same axis and have the same orientation to
\(z\). Each \(U\) determines \(x\) and the original
row. By the preceding height cutoff, the candidates
with \(U/g\le X^{\theta_1}\), \(g=\gcd(U,A)\), number
at most \(O(n^\alpha X^\varepsilon)\), with
\(\alpha<\frac14\liminf t\).

Put into a fixed pool \(P_x\) the primes of \(P_0\)
and all potential differences among the known vertices.
Its prime product is of polynomial height. A prime
of \(D_{xy}\) outside this pool divides no fixed
difference, and it divides no other moving difference
\(D_{xv}\), since it would divide the known
difference \(D_{yv}\). It is therefore private
in the sense of (PA4). It is enough to count
candidates with outside-\(P_x\) radical at most
\(n^{2(t+s_i)-\sigma}\) on \(D_{xy}\).

Fix \(g\mid A\) at divisor cost. The three terms
\(U/g,A/g,D_{xy}/g\) form a positive primitive
additive triple after ordering. Its fixed complete
letter is \(A/g\); its height is greater than
\(X^{\theta_1}\) on the remaining branch. The raw
radicals of the first two terms are bounded by the
original caps \(n^C\), and the third is bounded
by the tested shift's raw cap \(n^{B_0}\).
The support of \(g\) is in \(P_x\), since \(A\) is
one of the fixed labels. Thus dividing by \(g\)
leaves the test's outside radical unchanged.
Use (FL1) in each derived height shell with the
fixed complete letter \(A/g\) and the outside cap
on \(D_{xy}/g\). Because the derived height
is at least \(X^{\theta_1}\), the raw relative
caps there are uniformly \(O(L)\) and the fixed
pool has polynomial height in the shell scale.
With a small fixed normalized allowance \(\xi\),
this costs at most \(O(n^{t+s_i-\sigma/2+\xi}
X^\varepsilon)\) possibilities, including the
dyadic shells and the divisor choices after
absorbing them in the inner height allowance.
The signs have only boundedly many choices, and
\(g\) and the primitive triple recover \(U\) and
the omitted potential.
[Fixed-letter count with a prescribed pool (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Both bounds for the omitted potential save a fixed
power over \(KR_i\). Multiply by the list and
known-tree factor and sum over the \(3\binom m2\)
tests. Take the tolerances small compared with this
power in (LB1). This gives (PA4). To prove its
unconditional joint consequence, a failure of the
joint bound outside \(P_0\) and the original labels
would bound every raw shift radical by
\(n^{B_{\rm shift}+B_1+3m^2C}\). Choose \(B_0\) above
this exponent uniformly before applying (PA4).
The private prime supports then sum to more than
\(B_{\rm shift}\) in base \(n\), contradicting the
joint failure.

For the common-field application the same-axis
difference in axis \(i\) is always divisible by
\(D/D_i\), where \(D=FGM\): compare its two
complete labels against any one vertex in a
different axis. Their own packet agrees and the
opposite residue on that known vertex agrees.
The shifted values \(D_{xy}/(D/D_i)\) therefore
have height at most \(X^{1-w+w_i}\). Dividing this
fixed modulus leaves all private primes counted
in (PA4). The first-tree support floor
\(n^{r_{\rm tree}-\tau}\) holds generically for the
single arrays by the same supported-value count
as in (CS6), now comparing directly with (LB1).
If the chosen original tree labels and all these
divided same-axis shifts become \(Q\)-full principal
ideals in a common field, their root discriminant
satisfies
\[
    \log_n\operatorname{rd}(F_{\rm lift})>
       r_{\rm tree}-\tau+B_{\rm shift}
       -\frac{\sum_e h_e+\binom m2(3-2w)}{LQ}.
\]
This uses the sum of the height caps on the
original tree and the fact that
\(\sum_i(1-w+w_i)=3-2w\) for the shifts.
[Common-field counting from tree labels (RL8)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#a-common-field-on-an-actual-shared-value-box)
