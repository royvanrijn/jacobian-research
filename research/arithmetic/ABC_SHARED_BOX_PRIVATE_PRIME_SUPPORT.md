# Private prime support on every edge of a shared-value abc box

A complete label in a shared-value box can be counted while
leaving one of its incident potentials unrecorded. The other
potentials supply a fixed complete letter and a prescribed
pool. The only primes of the moving label outside this pool
are private to that label in the entire box.

## A bound for one position

Use the ordered \(m\)-boxes of
[the shared-value theorem (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem),
with \(m\ge2\). Denote the axis residue bounds by
\(R_1,R_2,R_3\), their product by \(R\), and the
pair-class bound by \(K\). The complete labels are the
positive integer differences of potentials in the three
axes. Every triangle with one potential in each axis is an
original primitive row in \(X/2<C\le X\). Set \(n=X^L\)
with \(L\) fixed before the height limit, and impose an
individual radical cap \(n^C\) for a fixed \(C>0\) on
every original letter. Fix a prime pool \(P\) whose prime
product is at most \(X^B\), for a fixed \(B>0\).

For a potential edge \(e\), write \(L_e\) for its complete
label. Define its private radical by
\[
 {\cal P}_{e,P}
   =\prod_{\substack{p\mid L_e,\ p\notin P\\
                     p\nmid L_f\text{ for every }f\ne e}}p .
\]
Here \(f\) ranges over all label positions of the box,
including positions of different roles. In particular, a
prime occurring in two positions with equal complete labels
is not private.

**PV1.** Fix a potential edge \(e\) joining axes \(i,j\).
For any fixed \(\xi>0\), at all sufficiently small fixed
\(L>0\), and uniformly in \(0\le b\le C\), the number of
ordered boxes with \({\cal P}_{e,P}\le n^b\) is
\[
 O_{m,B,C,\xi,L,\varepsilon}
 \left(
 X^\varepsilon R^m K^{3m-1}\,
 \frac{n^{\,b/2+\xi}}{K\max(R_i,R_j)}
 \right).
 \tag{PV1}
\]
The constants are uniform in the choice of \(P\), the
original shell family, the positions and the residue bounds.

For the proof omit one endpoint \(v\) of \(e\), say in
axis \(i\). All other ordered residue coordinates can be
recorded at cost at most \(R^m/R_i\). The induced graph of
the remaining potentials is connected, since \(m\ge2\).
Assigning the \(3m-2\) labels of a spanning tree determines
them up to common translation, at cost at most \(K^{3m-2}\).
These are the anchors. All anchor potential differences
have magnitude \(O_m(X)\).

Prescribe a pool consisting of \(P\), the supports of all
labels between two anchors, and the supports of all
nonzero differences between neighbors of \(v\). It has
uniform polynomial height at fixed \(m,B\). Its choice
depends only on the anchors. Every unknown label at \(v\)
is an affine function of its potential, up to sign,
of slope one with its zero at a known neighbor. The
neighbor potentials are distinct within each axis by
the residue tests, and across axes by positivity of the
original rows. A prime outside the prescribed pool which
divides \(L_e\) thus divides no other label incident to
\(v\), and no label between two anchors. The outside
radical of \(L_e\) in this pool is at most its private
radical relative to \(P\).

Complete \(e\) to a triangle with any known vertex on the
third axis. The edge between its two known vertices is a
fixed complete letter of an original primitive row. The
other two letters satisfy the original caps, and the
outside radical of \(L_e\) is at most \(n^b\). Apply
(FL1) in this genuine fixed-letter fibre. With the desired
normalized and inner allowances, it gives
\(O(n^{b/2+\xi}X^\varepsilon)\) choices for the row. The
complete value on \(e\) determines \(v\) because its
position and sign are known. This counts all boxes above
the anchors even though the residue of \(v\) was not
recorded. Use the endpoint in the axis with larger
\(R_i\) for (PV1). The proof also allows \(P\) to be
chosen as a function of the anchors for the specified
omitted endpoint, with the same uniform height bound.
[Prescribed-pool fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Simultaneous privacy and the outside cap

In a surviving simultaneous packet tuple, set
\(c_0=1/6\), \(d=s_1+s_2+s_3\), and \(t=c_0-d\),
where \(s_i\) is the price of the packet on the role
opposite axis \(i\). Assume fixed positive lower margins
for \(t\) and for the \(s_i\). Use the bounds
\[
 R_i=n^{s_i+O(\kappa+\eta)},\qquad
 K=n^{t+O(\eta)},\qquad
 E\ge n^{d+2t-\eta}.
\]
The alphabet bound \(R_i\) may be chosen at its indicated
scale from above; the actual alphabet also has the
corresponding lower bound by (SC1). The pair-class bound
\(K\) can be chosen at the stated scale by (SB2).
(LB4) gives at least \(n^{md+(3m-1)t-O_m(\kappa+\eta)}\)
actual ordered boxes, up to constants.
[Complete pair-class bounds (SB2)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#simultaneous-marking-theorem)
· [Actual ordered box lower bound (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

For any fixed \(\sigma>0\), take the fixed tolerances and
height allowances sufficiently small. Outside a power-small
fraction of all the boxes, every edge \(e\) joining axes
\(i,j\) satisfies
\[
              {\cal P}_{e,P}>
                 n^{\,2(t+\max(s_i,s_j))-\sigma}.
 \tag{PV2}
\]
For a positive exponent threshold at most \(C\), apply
(PV1) with that threshold and sum over at most \(3m^2\)
positions. A nonpositive threshold requires only the
nonempty-support test at \(b=0\) when needed. If a tested
threshold exceeds \(C\), apply (PV1) with \(b=C\); the
result would bound every box and contradict the lower
count. This also verifies that such a threshold cannot
arise in a retained full-exponent stratum.

The private supports in (PV2) are pairwise disjoint.
The full box radical outside \(P\) therefore exceeds
\[
 n^{\,2m^2(3t+\max(s_1,s_2)+\max(s_1,s_3)+\max(s_2,s_3))
             -3m^2\sigma}.
 \tag{PV3}
\]
For prices tending to zero, its coefficient tends to
\(m^2\), the sum of all \(3m^2\) individual radical caps
at \(2c_0\). This lower bound concerns actual integer
values with the original additive and primitive constraints.

Let \(e\) join axes \(i,j\) and write \(k\) for the
remaining axis. The fixed packet on its role has radical
at least \(n^{\beta_k}\). None of that packet support is
private, since the packet occurs on all \(m^2\) edges of
that role. With the balanced cap \(n^{2c_0+\kappa}\),
the radical on \(L_e\) left after removing its own packet
and its private primes has normalized exponent at most
\[
          2(s_k+\min(s_i,s_j))-\beta_k+\kappa+\sigma .
 \tag{PV4}
\]
In particular, it includes every prime outside its own
packet shared with another position, as well as any
prime outside its packet in the prescribed external pool.
This is an individual bound holding on almost all full
arrays. It is uniform in that external pool.

## A height-first statement for private valuation measures

Suppose a marked stratum as above has a cohort of ordered
boxes on the lower-bound exponent of (LB4) in the
height-first small-budget limit. For one fixed position
\(e\) joining axes \(i,j\), assume that the scaled measures
restricted to primes private to \(e\) and outside \(P\)
have a common weak limit \(\mu_e\), uniformly within the
retained refinement cells. This follows on a subsequence
by refining the exact layer masses at height-subpower cost
at each fixed budget before taking the small-budget limit.

**PV5.** In this setting
\[
        (\mu_e)_{\rm na}([0,\infty))
                       \ \ge\ t+\max(s_i,s_j).
 \tag{PV5}
\]
To prove it, omit the endpoint on the larger-price axis
as in (PV1). Average the full cohort over the anchor
records, at normalized cost
\(md+(3m-1)t-(s_i+t)+o(1)\) for omission of axis \(i\).
The obtained fixed-anchor family has exponent at least
\(t+s_i-o(1)\). Its rows in the known triangle are
distinct, since its missing potential is determined by
one moving complete value. Use the pool from (PV1) on
that free role and apply (FT1) to this actual fixed-letter
fibre. Its exterior exact valuation measure is dominated,
row by row, by the private measure. After the fixed-budget
refinements and weak limits its nonatomic mass is at least
\(t+s_i\), hence so is that of \(\mu_e\).
[Fixed-fibre nonatomic mass outside a prescribed pool (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

For example, if the limiting complete measure on role
\(k\) has nonatomic mass \(c_k\) and its fixed packet has
nonatomic mass \(b_k\), then necessarily
\[
     c_k-b_k\ge t+\max(s_i,s_j).
\]
The packet and private supports are disjoint and their
measures both lie in the actual complete profile. The
nonatomic part of the residual outside those two supports
has mass at most \(c_k-b_k-t-\max(s_i,s_j)\). This gives a
separate profile constraint on every role in a full array.

## Two arrays with the same residue lists

For synchronized pairs as in (CS3), first fix one
entire array. In the recovery count for any one edge of
the other array, include all labels of the fixed array
in \(P\). The missing residue is already fixed by
synchronization, so the anchor cost is \(K^{3m-2}\)
at fixed lists and first array. The same proof gives,
simultaneously on the two arrays outside a power-small
fraction of synchronized pairs, the private radical
bound \(n^{2t-\sigma}\) on each label even after excluding
every label of the other array. Indeed, each individual
exception count in a fixed pair of positions is at most
\(O(R^m K^{6m-2} n^{-\delta})\) using the remaining
choices for the fixed array. Compare with (CS3), and then
interchange the arrays. Thus the supports counted this
way are disjoint across both original arrays.
[Synchronized pair scale (CS3)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#fresh-direction-radicals-on-synchronized-actual-boxes)

A further saving requires an upper bound on the supply
of private primes in a full array, or a way to use those
primes to reconstruct the array at a smaller cost.
