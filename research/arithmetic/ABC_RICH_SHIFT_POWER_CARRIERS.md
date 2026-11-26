# Polynomial and powered descriptions of rich additive shifts

A rich shift is present in many differences between values of either
varied complete letter. These differences give an exact additive graph
with the shifts themselves as a part. Polynomial or powered
descriptions of that part can force a saving even when the shift
radicals are much larger than the original radical caps. Near-power
descriptions give a second bound, with a stronger conclusion when
the coefficient divides the shift.

## From cherries to a shift graph

Use the cherry and shift notation of (PSA1). Put \(n=X^L\),
\(c_0=1/6\), and let the complete part sizes be at most
\(n^{c_0+o(1)}\). Suppose a cohort \(\mathcal H\) of shifts
has exponent \(\tau\), with \(M_h\ge n^{1/3-o(1)}\) for
every \(h\in\mathcal H\). In either prescribed varied role,
project each cherry to its smaller and larger positive
values \(u,v\). Then \(u+h=v\), even when the chosen role
was the signed target. For fixed \(h,u,v\), choosing the
other varied role's smaller positive value determines the
remaining labels. There are at most \(n^{c_0+o(1)}\)
such choices. Consequently the simple graph of distinct
projected rows
\[
                             u+h=v                              \tag{SP1}
\]
has exponent at least \(c_0+\tau-o(1)\), and the two source
parts have sizes at most \(n^{c_0+o(1)}\). No common divisor
is removed in this projection. Neither the cap on
\(\operatorname{rad}(h)\) nor primitivity of the projected
triples is needed below. With fixed part and richness
tolerances, the same argument incurs only a fixed multiple
of those tolerances in the exponent.
[Cherries and the rich-shift count](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md#a-bound-for-rich-shifts-with-a-radical-cap)

## Exact polynomial descriptions

Suppose each \(h\in\mathcal H\) can be assigned a representation
\[
                    h=t_{s,\lambda}+sP_{s,\lambda}(z),
                                                               \tag{SP2}
\]
where \(s\) is a positive integer and \(z,t_{s,\lambda}\)
are integers.
The dictionary of \(s\)'s has size at most \(n^{f+o(1)}\).
For each \(s\) there are at most \(X^{o_X(1)}\) prescribed
types \(\lambda\), each with its own integer polynomial
\(P_{s,\lambda}\) and offset \(t_{s,\lambda}\) fixed
throughout its use. Assume
\[
              2\le\deg P_{s,\lambda}\le (\log X)^B
\]
for some fixed \(B\). All the subpower allowances here
are uniform over the types at each fixed budget. Then
\[
                    \tau\le\max\{f,\ c_0/2+f/2\}.
                                                               \tag{SP3}
\]
Thus a full-size rich shift cohort requires \(f\ge1/6\).
In particular (SP3) applies to shifted weighted exact
powers of any order \(k\ge2\) in the displayed degree
range. Positive exact powers of integers at most \(X\)
have \(k\le\log_2X\) when the base exceeds one. For a base
of one one may use one fixed allowed order. The polynomial
argument allows arbitrary coefficients and offsets, since
only differences within a type enter the proof.

Label a row of a bipartite graph by \((v,s,\lambda)\)
and a column by \(u\). Each chosen (SP1) row yields one
edge, with at most one edge at a given row and column
because they determine \(h=v-u\). There are at most
\(n^{c_0+f+o(1)}\) row labels and \(n^{c_0+o(1)}\)
columns. Fix two distinct columns \(u_1,u_2\), and
put \(D=u_2-u_1\ne0\). On a common row,
\[
       s\big(P_{s,\lambda}(z_1)-P_{s,\lambda}(z_2)\big)=D .
\]
The offset cancels and \(s\mid D\), so the choices of \(s\)
cost only \(X^{o_X(1)}\). After selecting the type, the
nonzero difference \(j=z_1-z_2\) divides \(D/s\), by
the integer polynomial difference identity. For fixed
\(j\) the polynomial \(P_{s,\lambda}(z+j)-P_{s,\lambda}(z)\)
has degree \(\deg P_{s,\lambda}-1\) and nonzero leading
coefficient. Hence it takes the prescribed value \(D/s\)
at at most \(\deg P_{s,\lambda}-1\) integers. This
determines the row label via \(v=u_1+h_1\). The common-row
count for two distinct columns is therefore \(X^{o_X(1)}\).

For \(E\) edges and \(R\) row labels, the ordered column
pairs at the rows number at least \(E^2/R-E\). The last
bound makes this at most \(n^{2c_0+o(1)}\). Solving gives
\(E\le n^{\max(c_0+f,\,3c_0/2+f/2)+o(1)}\). Compare with
(SP1) to get (SP3).
[Weighted difference count for power records](ABC_SHIFT_DICTIONARY_WEIGHTED_CODEGREES.md#a-finite-graph-bound)

A coefficient dictionary cap \(n^f\) follows, for example,
from an individual coefficient radical cap \(n^f\) and a
fixed polynomial bound on coefficient heights. The same
is true if the radical outside a prescribed polynomial-height
prime pool is so bounded. Fix the pool-supported powers
at height subpower cost, and enumerate the outside radical
and powers.
[Fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

### A small dictionary of errors

Suppose instead that \(h=t_{s,\lambda}+sP_{s,\lambda}(z)+e\)
with an error dictionary of size at most \(n^{g+o(1)}\).
Retain one error value. Apply the same proof with the
modified offset \(t_{s,\lambda}+e\). The cohort loss is
at most \(g\), so
\(\tau\le g+\max\{f,c_0/2+f/2\}\).
The numerical height of the error plays no part in this
version.

## Approximations by high powers

Suppose the chosen source values \(u,v\) have height at least
\(X^{h_*}\) and radical at most \(n^A\), where \(h_*,A>0\)
are fixed. Choose one approximant \(s x^k\) for each \(h\),
with
\[
 |h-sx^k|\le X^\theta,\qquad
 X^{h_*}\le sx^k\le 2X,\qquad x\le X^{C_0L},\qquad
 k\ge \lceil L^{-\alpha}\rceil ,
\]
where \(C_0>0,\alpha>1/2\) are fixed and \(s,x,k\) are
positive integers. In the two statements below, an arbitrary
fixed exponent allowance \(\Delta>0\) is permitted and
\(\theta>0\) is chosen sufficiently small in terms of that
allowance and the fixed data. Take the budgets sufficiently
small and fixed before taking heights. In that convention the
bounds are
\[
 \begin{array}{ll}
 \tau\le c_0/2+f+\Delta
   &\text{if }\operatorname{rad}(s)\le n^f,\\[2mm]
 \tau\le\max\{f,\ c_0/2+f/2\}+\Delta
   &\text{if }s\mid h\text{ and the }s\text{ dictionary has cap }n^f .
 \end{array}                                                    \tag{SP4}
\]
Neither statement imposes a raw radical cap on the exact
shift \(h\). The second uses a divisibility condition on its
error, since then \(h-sx^k\) is a multiple of \(s\).

For the first bound replace \(h\) in (SP1) by its chosen
approximant. The triple \((u,sx^k,v)\) determines the
exact shift \(h=v-u\), so this replacement preserves
edge cardinality. Its shift dictionary has exponent
at most \(\tau\), and its radical cap is \(n^{f+C_0}\)
by the bound on the base. The real-filtered grid bound
(OG5) applies with the two other caps \(c_0+o(1)\).
A constant enlargement of the height bound from \(X\) to
\(2X\) leaves all exponent caps unchanged. The bound gives exponent at most
\[
                \max(c_0,\ c_0+\tau/2,\ 3c_0/2+f)
\]
up to arbitrarily small fixed allowances. Comparison with
\(c_0+\tau\) proves the first bound in (SP4); the diagonal
and middle term already save if \(\tau\) is separated from
zero.
[Real-filtered grid graphs (OG5)](ABC_ONE_POWER_GRID_ADDITIVE_GRAPHS.md#real-filtered-additive-graphs-at-tall-heights)

For the second bound keep the exact graph (SP1) and label
each row by \((v,s,k)\). The \(k\)'s give only a subpower
factor, using the height bound for bases greater than one
and one fixed allowed order at base one. The number of
rows is at most \(n^{c_0+f+o(1)}\).

Fix \(u_1\ne u_2\). Their common rows obey
\[
 \begin{split}
 D&=u_2-u_1=h_1-h_2,\qquad s\mid D,\\
         |s(x_1^k-x_2^k)-D|&\le 2X^\theta .
 \end{split}
\]
Choose a small fixed \(h_0>0\) so that (SI2) applies to
the source dictionaries on intervals of width \(4X^{h_0}\),
with arbitrarily small normalized allowance. Then take
\(\theta<h_0/2\), with room for constant factors.
For \(|D|\ge X^{h_0}\), enumerate \(s\mid D\) at divisor
cost and the orders at subpower cost. Apply (WP4) to
the two powered terms and to the positive target \(|D|\),
exchanging them if necessary. All three have positive fixed
numerical heights, and the bases have the required short
height. The orders exceed the fixed auxiliary degree at
sufficiently small fixed \(L\). This gives \(n^{\eta+\varepsilon}\)
base pairs for any prescribed \(\eta>0\), with the height
allowance reduced as needed. If one base is one, it gives
a fixed term and the other base is unique per order for
large \(X\), by spacing of the tall power values. Both
bases equal to one cannot give \(|D|\ge X^{h_0}\).
For each resulting pair, \(v\) is within \(X^\theta\)
of \(u_1+s x_1^k\). Count \(v\) by (SI2) on the original
source dictionary. Thus each separated column pair has
at most \(n^{O(\eta)+o(1)}\) common rows.

For \(|D|<X^{h_0}\), every column has at most
\(n^{\eta+o(1)}\) nearby columns in its original
dictionary, by (SI2). Across all rows this accounts
for at most \(E n^{\eta+o(1)}\) ordered column pairs.
With \(R\) row labels, Cauchy--Schwarz now gives
\[
             E^2/R-E\ \ll\
                      E n^{\eta+o(1)}+n^{2c_0+O(\eta)+o(1)} .
\]
Insert the row cap and compare with (SP1), taking \(\eta\)
small in terms of \(\Delta\). This proves the second
bound in (SP4). The proof uses the size of the coefficient
dictionary. It allows the coefficient radicals to be large
when \(s\mid h\).
[Short-interval radical count (SI2)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#arbitrarily-small-normalized-counts-in-shorter-intervals)
· [Weighted-power near-equality bound (WP4)](ABC_WEIGHTED_POWER_NEAR_EQUALITY.md#both-signs-and-all-omitted-letter-orientations)

For \(\tau=c_0=1/6\), the first threshold in (SP4) is
\(f\ge1/12\) up to fixed margins and the second is
\(f\ge1/6\). A full-size cohort near tall perfect powers
of short base is excluded by the first bound, even if the
actual shift radicals are large. The exact result (SP3)
also excludes pure powers without the short-base condition,
including squares at the full numerical height.

## Tall prime shifts in a primitive graph without radical caps

The exact addition graph and the numerical sparsity of short
progressions can coexist with a full-size family of tall prime
shifts. The following countermodel omits the original radical
caps. Fix sufficiently small \(L>0\), put \(n=X^L\) and
\(N=n^{1/6}\), and take integral heights \(X\to\infty\).
There are three complete parts at heights comparable to \(X\),
each of size \(N^{1+o_X(1)}\), with primitive edges in
\(X/2<c\le X\), and \(N^{1-o_X(1)}\) distinct prime shifts
\(h\asymp X\) each with \(\gg N^2\) cherries. Each complete
part and the prime shift cohort meet any ordinary arithmetic
progression with at most \(X^{1-\sqrt L}\) terms in only
\(N^{o_X(1)}\) values.
[Earlier progression-avoidance construction](ABC_MODULAR_PROGRESSION_SLICING.md#a-primitive-additive-graph-with-uniform-progression-avoidance)

Put \(r=\lfloor\sqrt{\log N}\rfloor\),
\(D=\lfloor N^{1/r}\rfloor\) and
\(J=\lfloor X/(100rD)\rfloor\).
Take \(g_i\) independently and uniformly in \([1,J]\cap
\mathbb Z\) for \(1\le i\le r\). Write \(g\cdot u\) for
the linear map on the box \(\{0,\ldots,2D\}^r\). Its
values lie in \([0,X/50]\). Any two distinct coefficient
vectors collide with probability at most \(1/J\), by
conditioning on all but a generator whose difference
coefficient is nonzero. The map is injective with
probability \(1-o_X(1)\) since \(2\log_X N<1\).

For progression avoidance fix a positive integer \(j\).
For any \(j+1\) affinely independent coefficient vectors,
choose a nonzero \(j\)-by-\(j\) minor from their difference
matrix. Condition on the other generators. Any specified
list of the \(j\) image differences determines the remaining
generators uniquely. If the images lie in a progression
of step \(d_0\) with at most \(X^{1-\sqrt L}\) terms,
the differences have at most \(2X^{1-\sqrt L}+1\)
options each, all multiples of \(d_0\). Sum over the
\(O(X)\) possible positive steps at most \(X\) and the
affinely independent tuples. The expected number of such
tuples with a step is at most
\[
 X^{1+(j+1)L/6-j\sqrt L+o_X(1)} .
\]
For sufficiently large fixed \(j\) this tends to zero.
Steps larger than \(X\) give singletons. Thus the
coefficient vectors in each eligible intersection have
affine dimension less than \(j\), so projection onto
pivot coordinates bounds their number by
\((2D+1)^{j-1}=N^{o_X(1)}\). Fix lower generators satisfying
both conclusions. They also hold on subsets and translates
of the box image, since the tests used only differences.

Now take \(g_0\) uniformly in \([X/10,X/5]\cap\mathbb Z\).
For each vector \(\delta\) with
\(0\le\delta_i\le\lfloor D/(10r)\rfloor\), let
\(h_\delta=g_0+g\cdot\delta\). The shifts are distinct,
and there are \(N^{1-o_X(1)}\) of them. Their lower-box
offsets are at most \(X/100\). For every such offset the
range of \(h_\delta\) contains the fixed interval
\([X/9,X/6]\), up to harmless endpoints. The prime number
theorem therefore gives probability \(\gg 1/\log X\)
that a given \(h_\delta\) is prime. Fix \(g_0\) with
\(N^{1-o_X(1)}\) prime shifts by averaging.

Let \(S=g\cdot\{0,\ldots,D-1\}^r\), and choose
\(\lambda_A,\lambda_B\) independently uniformly in
\([X/4,X/3]\cap\mathbb Z\). Use the parts
\[
 \begin{split}
 A&=\lambda_A+\bigl(S\cup(g_0+S)\bigr),\\
 B&=\lambda_B+\bigl(S\cup(g_0+S)\bigr),\\
 C&=\lambda_A+\lambda_B+g_0+(S+S).
 \end{split}
\]
The displayed intervals put all targets in the stated
height shell. Each part is a union of at most two
translates of subsets of the large box image. The same
is true of the chosen prime shifts, so the established
progression bound applies. Take all positive primitive
edges on these parts.

For a chosen \(h=h_\delta\), every \(u,v\in
\{0,\ldots,D-1\}^r\) with \(u+\delta,v+\delta\) also
in that box gives a candidate cherry. Use
\(a=\lambda_A+g\cdot u\), \(b=\lambda_B+g\cdot v\)
and the two triples \((a+h,b,a+b+h)\),
\((a,b+h,a+b+h)\). There are
\(\prod_i(D-\delta_i)^2\gg N^2\) distinct candidates
uniformly in \(\delta\), since
\((1-1/(10r))^r\) is bounded below and \(D^r\sim N\).

Let \(J_0\asymp X\) be the number of possible translations
for each summand. For each of the two original pairs,
independence of the translations makes the probability of
a common prime divisor \(p\) at most \((1/p+1/J_0)^2\).
A union bound over both pairs gives probability of either
failure at most
\[
   2\sum_{p\le X}(1/p+1/J_0)^2
       \le 2\left(\frac14+\sum_{k\ge1}\frac1{(2k+1)^2}\right)
                    +O(\log X/X)
       =\frac{\pi^2}{4}-\frac32+o(1)<1 .
\]
Thus each candidate contributes a primitive cherry with
probability bounded below by an absolute constant.
Sum over all chosen prime shifts and average over the
translations. Some translations leave \(\gg N^2\) cherries
for a fixed positive fraction of those shifts, since
there are at most \(N^2+o(N^2)\) candidates per shift.
This yields the stated cohort and a near-quadratic graph,
with original radical budgets left free.
