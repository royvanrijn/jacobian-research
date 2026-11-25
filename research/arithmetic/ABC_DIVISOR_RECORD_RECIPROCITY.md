# Divisor record frequencies and overlaps at saturation

Individual prime frequencies at the one-third endpoint have the
reciprocal scale. For a bounded group of root coordinates, the same
scale governs arbitrary simultaneous divisor records, including
records of powers. This yields common-divisor fibres with very small
extra overlap inside the group.

Use primitive shell boxes \(X/2<c\le X\) and the anatomic coordinates
\(x_i\) with raw height caps \(X^{\alpha_i+o_X(1)}\). Choose one
representation for each row. Fix the initial coefficients, all degree
lists and all truncation parameters before taking the height limit at
each fixed \(\lambda\). Let \(\mathcal B\) be a family of \(X^{D+o_X(1)}\)
distinct rows in the box. A saturation sequence here has
\[
 \frac{D}{\lambda}\longrightarrow\frac13,\qquad
 \frac{\sum_{i\text{ on }T}\alpha_i}{\lambda}
        \longrightarrow\frac13\quad (T=A,B,C)
\]
in the outer small-budget limit. Use the anatomic truncation with
initial coefficient heights \(X^{o(\lambda)}\) and an
\(o(\lambda)\) loss in the raw radical caps.
[Bernert–Browning–Lichtman–Teräväinen, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

## Reciprocal frequencies for joint divisor records

The following partitions are permitted:

1. Partition each letter separately into the groups of (GD3).
2. Partition the two letters opposite a letter with arbitrarily cheap
   full moduli as in (CG3), using groups of (CG1), including its
   bounded shrinking-window extension.
3. Partition the tall pair when the third letter has height at most
   \(X^{1-\tau}\) for fixed \(\tau>0\), using groups of (RG1).

For a group \(J\) and a divisor vector \(g=(g_i)_{i\in J}\), write
\[
 b(g)=\sum_{i\in J}\log_Xg_i,\qquad
 \pi_J(g)=\Pr_{\mathcal B}(g_i\mid x_i\text{ for every }i\in J).
\]
Here \(\Pr_{\mathcal B}\) is uniform on the rows. For a row \(z\), set
\[
 \Delta_J(z)=\max_{g_i\mid x_i(z)\ (i\in J)}
                 \bigl(b(g)+\log_X\pi_J(g)\bigr)_+.
\]
Then for every fixed \(\eta>0\), sufficiently far along a saturation
sequence, the exceptional set
\[
 \sum_J\Delta_J(z)>\eta\lambda                                      \tag{DRR1}
\]
has proportion \(O(X^{-\delta\lambda+o_X(1)})\), for some fixed
\(\delta>0\) depending on \(\eta\). Uniformity is subject to the
group and lock conditions in the stated estimates. At fixed budget,
one may impose this simultaneously for all permitted partitions of
the original finite coordinate lists, provided the estimates hold
with uniform margins. For the cheap opposite modulus case fix the
same supplying letter and a lock for the indicated tolerance.

Indeed, write \(\alpha_J=\sum_{i\in J}\alpha_i\). At fixed truncation
each row lies in only \(X^{o_X(1)}\) divisor records for \(J\), by
the divisor bound. Hence \(\sum_g\pi_J(g)\le X^{o_X(1)}\). Bin
\(b(g)\) into intervals of mesh at most \(\xi\lambda\). In a bin
with upper endpoint \(u\), the records satisfying
\(b(g)+\log_X\pi_J(g)\ge d>0\) number at most \(X^{u-d+o_X(1)}\).
For each such record within the bin, the number of possible
root vectors from the raw ranges is at most
\[
 \prod_{i\in J}\left(\frac{X^{\alpha_i+o_X(1)}}{g_i}+1\right)
                   \le X^{\alpha_J-u+\xi\lambda+o_X(1)}.
\]
In the last estimate each \(g_i\) divides a root in the range.
Thus a group on which \(\Delta_J\ge d\) has a dictionary of exponent
at most \(\alpha_J-d+\xi\lambda+o_X(1)\), after summing the bins.
All rates in this paragraph are at fixed truncation. For \(d=0\)
use the raw dictionary.

Bin all the \(\Delta_J\) with total downward rounding error at most
\(\eta\lambda/4\), and take \(\xi\) so small that the bin errors in
the dictionary bounds sum to at most \(\eta\lambda/4\). Both choices
are made before the height limit. For each bin profile satisfying
(DRR1), the summed dictionary cost loses at least \(\eta\lambda/2\)
against the raw group caps. There are only boundedly many bin
profiles at each fixed budget and truncation. In the first case,
(GD3) bounds the count by one third of the summed cost up to the
fixed small-budget allowance. In the second case, use (CG1)
with a full divisor whose record costs less than \(\eta\lambda/16\).
In the third case, (RG1) applies directly. In the last two cases
the sum of the raw free caps is \((2/3+o(1))\lambda\), and the
dictionary cost enters with coefficient one half. In every case,
the count is power smaller than \(X^{D+o_X(1)}\). This proves
(DRR1).
[Within-letter count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [Opposite full modulus count (CG1)](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#one-divisor-and-two-free-letters)
· [Real tall pair count (RG1)](ABC_REAL_STRIP_GROUP_DICTIONARIES.md)

The other side of the frequency estimate comes from (RP2), applied
with tolerance small enough that the errors sum to at most
\(\eta\lambda\) over the partition. The group count is fixed at
each height limit, so this yields a vanishing proportion of
exceptions, also when the number of groups grows in the subsequent
budget limit. In particular, on typical rows all their divisor
tuples in any one permitted group have
\(-\log_X\pi_J(g)=b(g)+o(\lambda)\), uniformly in the record. More
generally the sum over groups of the maximal absolute error in
each group is \(o(\lambda)\) in probability in the iterated
endpoint limit. The lower frequency bound from (RP2) alone permits
arbitrary divisor records across the entire tuple. The upper
bound here requires a permitted group.
[All-coordinate lower recurrence (RP2)](ABC_RARE_PRIME_OCCURRENCE.md#joint-recurrence-of-divisor-records)

## Small extra overlap inside a common-divisor group

Fix a group \(J\) as above, \(b_0<1/3\), and \(\rho>0\). Define
\(P_J(z)=\{p:p\mid \prod_{i\in J}x_i(z)\}\). At sufficiently late
fixed budgets and then at large heights, a typical row \(z\) has
the following property for every divisor record \(g\) on \(J\)
at \(z\) with \(b(g)\le b_0\lambda\). The fibre \(\mathcal B_g\)
contains, including \(z\), a family of size at least
\[
      X^{(\min\{\rho,\,1/3-b_0\}-\sigma)\lambda}
                      \qquad(\sigma>0\text{ fixed})                 \tag{DRR2}
\]
on which distinct rows \(z_1,z_2\) satisfy
\[
 \sum_{\substack{p\in P_J(z_1)\cap P_J(z_2)\\p\nmid \prod_i g_i}}
                         \log_Xp <\rho\lambda .                    \tag{DRR3}
\]
As usual the claim of exponential size concerns positive exponents.
A stronger whole-union conclusion is available for a shorter
sequence: for every fixed \(\tau>0\) the fibre contains
\(\lfloor(\log X)^\tau\rfloor\) ordered rows starting with \(z\)
such that the sum in (DRR3), using the union of the previous
\(P_J\) in place of \(P_J(z_2)\), is still less than \(\rho\lambda\).

For the proof, first retain the typical subfamily \(\mathcal B'\)
satisfying the upper frequency bound in (DRR1) with a small tolerance
\(\eta\). Apply (RP2) again inside \(\mathcal B'\). A typical row
in \(\mathcal B'\) is in each fibre \(\mathcal B'_g\) of size at least
\(|\mathcal B'|X^{-b(g)-\eta\lambda}\). Fix an anchor \(y\) in this
fibre and put \(P=P_J(y)\setminus\{p:p\mid\prod_i g_i\}\).
A bad query \(x\) in the fibre is in the record \(h\) obtained
by multiplying \(g_i\) by the product of the distinct primes
of \(P\) occurring in position \(i\) at \(x\). These extra
prime factors are assigned to unique positions because the
positive-index roots on a row are pairwise coprime. In particular,
\[
 b(h)=b(g)+\sum_{p\in P\cap P_J(x)}\log_Xp
                             \ge b(g)+\rho\lambda .
\]
Since \(h\) occurs on \(x\in\mathcal B'\), its occurrence in the
original family is at most \(|\mathcal B|X^{-b(h)+\eta\lambda}\).
For a single anchor \(P\) has \(O(\log X/\log\log X)\) primes
since it is supported on a polynomial-height integer. Thus there
are only \(X^{o_X(1)}\) possibilities for \(h\), at fixed
truncation. At every greedy step, discarding the selected row and
the violating queries removes at most
\[
 1+|\mathcal B|X^{-b(g)-\rho\lambda+\eta\lambda+o_X(1)}
\]
rows. The fibre lower bound and the small choice of \(\eta\)
give (DRR2).

For the whole-union conclusion let \(P\) be the union of the
previous anchor supports outside \(\prod_i g_i\), of cardinality
\(O((\log X)^{1+\tau}/\log\log X)\). Choose a fixed \(a<1\) with
\(a(1+\tau)>\tau\). In the same bad query bound, sum
\(X^{-b(h)}\) over all \(h\) with extra radical of logarithmic
mass at least \(\rho\lambda\). Since each extra prime has at
most \(|J|\) possible positions, this sum is at most
\[
 X^{-b(g)-(1-a)\rho\lambda}
       \prod_{p\in P}(1+|J|p^{-a})
                     =X^{-b(g)-(1-a)\rho\lambda+o_X(1)}.
\]
The subpower factor follows by ordering the primes of \(P\) and
using \(\sum_{j\le |P|}(j+1)^{-a}=o(\log X)\). This gives a
power-smaller set to avoid at every step, provided
\(2\eta<(1-a)\rho\). Previously chosen rows are also power
smaller than the fibre. This proves the sequence assertion.

If all degrees on \(J\) have \(e_i\lambda\le v\) and its roots
are squarefree, stripping the shared support outside the mandatory
record loses at most \(v\rho\) of numerical height in base \(X\)
on each selected row. The parts left after removing the mandatory
prime support and all previous supports are pairwise coprime.
This stripping conclusion only concerns the indicated positions.

## Anisotropic intersections for an arbitrary record

For this estimate the record can extend across any number of
groups or letters. Let \(\beta_T\lambda\) denote the logarithmic
radical of its squarefree support on letter \(T\), and let \(b\lambda\)
be the total logarithmic size of the record, allowing powers in
the record itself. Thus \(\sum_T\beta_T\le b\).
For each letter \(T\) let \(P_T\) be a prescribed prime set,
and assign prices \(w_T(p)\) satisfying (EB3) uniformly,
with the ordinary price \(\log_X p\) off \(P_T\). On a row
\(z\), write \(q_T(z)\lambda\) for the total discount on the
distinct primes of the complete letter \(T\). For fixed
distinct \(U,V\), the rows with \(q_U(z)+q_V(z)\ge q\)
satisfy
\[
           \#\{z:q_U(z)+q_V(z)\ge q\}
               \le X^{(1/3-q/2+o(1))\lambda+o_X(1)}.              \tag{DRR4}
\]
Here \(o(1)\) can be taken smaller than any fixed positive
tolerance at late budgets. To see this, the coordinate
enumeration from (EB3) is valid separately on each letter,
with pairwise coprimality inside a row. The sum of the charged
caps for the coordinates on \(T\) is at most
\(\sum_{i\text{ on }T}\alpha_i-q_T(z)\lambda+o(\lambda)+o_X(1)\).
Indeed the sum of unpriced radical masses of the positive-index
roots is at most their raw cap, and any discount on a prime missed
by all those coordinates lies in the small initial coefficient. Bin the charged
caps with arbitrarily small total rounding error and apply (GD3)
to this pair at each bin. This gives (DRR4).
[Priced coordinate enumeration (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

For a prescribed anchor \(y\) belonging to the same record fibre,
take the zero-priced prime pool on every letter to be the primes
of the entire product \(abc\) of \(y\). It satisfies (EB3) uniformly
by the divisor bound for the number of distinct primes. Then
\(q_T(z)\lambda\) is the shared radical on the query letter \(T\),
and \(q_T(z)\ge\beta_T\) for every query in the fibre. By (RP2),
a typical record fibre has size at least
\(X^{(1/3-b-o(1))\lambda+o_X(1)}\) when \(b<1/3\).
In particular for each fixed \(\eta>0\), queries sharing more
than \((2b+\eta)\lambda\) radical weight with the anchor on
any specified pair of letters form a power-smaller subset of
that fibre. The same greedy selection as above gives a family
of size \(X^{(\min\{\eta/2,1/3-b\}-\sigma)\lambda}\)
such that the later row of each pair has overlap at most
\((2b+\eta)\lambda\) in every pair of its letters with the
complete prime support of the earlier row. To make this
symmetric on a saturated subfamily, first choose a prime
colouring by complete letter. A random three-colouring
retains each row with probability \(3^{-\omega(abc)}
=X^{-o_X(1)}\), so a single colouring preserves the exponent.
Apply recurrence inside the retained family. Every prime
shared by two rows lies on the same letter in both of them,
making all three overlap sums symmetric.

When the record is all on one letter, the first orientation bounds
additional overlap on that letter by roughly \(b\lambda\), because
its mandatory radical weight is already \(b\lambda\) for a squarefree
record. The local result (DRR2) makes this extra overlap arbitrarily
small inside a permitted single group.

## Reciprocal records on a large fixed-value fibre

Fix a complete letter in the shell. Its finitely many original
coordinates lie in subpower dictionaries since all divide the
fixed value. For either free letter \(U\), let \(L_U=\sum_{i\in
U}\alpha_i\) and let its complete radical cap be \(X^{s_U}\).
In the anatomic reduction \(L_U\le s_U+o(\lambda)+o_X(1)\).
Equation (GD3) with the fixed letter at cost zero gives
\[
               D\le W_U/2+o(\lambda)                       \tag{DRR5}
\]
for every allowed within-letter partition of \(U\) and its
dictionary exponent \(W_U\). In particular, the local reciprocal
record proof and both overlap selections apply to a pointwise
saturation sequence \(D/\lambda\to1/6,\ L_U/\lambda\to1/3\),
using (DRR5) in place of the original count and \(1/6\) in
place of \(1/3\) in the size of the recurrent fibre.

There are two uniform pointwise restrictions. First, suppose a
free letter \(U\) has a prime \(p\ge X^t\), with \(0\le t\le s_U\).
For a family of \(X^{D+o_X(1)}\) such rows, let \(\pi_U(p)\)
be its frequency on \(U\). Apply the same prices as in (EB7)
on \(U\) only. Every occurring prime has
\(\pi_U(p)\ge X^{-D+o_X(1)}\), so the imposed prime discounts the
charged cap by at least \((t-D)_++o_X(1)\). Binning the prices
and using (DRR5) gives
\[
        D\le \min\{s_U/2,s_U-t\}+o(\lambda).                 \tag{DRR6}
\]
If \(D\ge t\), the trivial consequence \(D\le s_U/2+o(\lambda)\)
of (DRR5) suffices whenever \(t\le s_U/2\). When \(t>s_U/2\) it
first gives \(D<t\) in the limiting exponents, so the discount
calculation applies. In particular, at a balanced fixed value,
a free prime of mass at least \((1/6+\eta)\lambda\) gives an
exponent at most \((1/6-\eta)\lambda+o(\lambda)\).
[Empirical prime prices (EB7)](ABC_EFFECTIVE_RADICAL_BUDGET.md#occurrence-frequencies-and-anatomic-entropy-at-saturation)

Second, let \(J\) be any one permitted group on \(U\). Its whole
root vector has at most \(|\mathcal B|\) possible values, giving
\[
                D\le L_U-\alpha_J+o(\lambda).               \tag{DRR7}
\]
This follows by charging \(J\) at cost at most \(D\), charging
the remaining coordinates at their raw caps, and using (DRR5).
For a fixed finite set of atoms of the limiting scaled degree
measure of \(U\), take \(J\) to be the union of shrinking windows
about those atoms. Choose their widths to decrease slowly enough
that the normalized raw mass in the windows tends to at least the
sum of the chosen atomic masses. A window about zero is permitted
by the high projection argument in (GD3). Taking the supremum over
finite sets of atoms, the total atomic mass on each free letter,
normalized by \(\lambda\), is at most \(l_U-d\), where
\(d=\lim D/\lambda\) and \(l_U=\lim L_U/\lambda\). At balance and
pointwise equality this is at most \(1/6\).

## A support-labelled obstruction across many groups

The three-residue obstruction can also be equipped with squarefree
prime supports. Use its complete graph on \(k\) vertices and
\(P=\binom{k}{2}\) edges. A row chooses one independent uniform
Latin triple \((a_e,b_e,a_e+b_e)\) in \(\mathbb F_q\) for each
edge \(e\). For each incidence of a vertex and an edge, and
for each letter, assign a list of \(q\) distinct primes of size
\(\asymp_k q\log q\), one prime per symbol in that position.
Take all lists to be disjoint. At a fixed \(k\) these lists
exist by the prime number theorem. A coordinate for a vertex
and a letter is the product of its selected incident-edge
primes. Every selected prime has logarithmic mass
\(\log q+O_k(\log\log q)\). The total raw logarithmic size is
\(6P\log q+O_k(\log\log q)\). Taking \(q\) with
\(\log_Xq=(1+o_X(1))\lambda/(6P)\), reduced by a logarithmic
factor to account for the prime sizes, gives radical at most
\(X^\lambda\) and \(q^{2P}=X^{(1/3+o_X(1))\lambda}\) rows.
[Three-residue dictionary model](ABC_THREE_RESIDUE_DICTIONARY_MOMENT_OBSTRUCTION.md)

For a divisor record inside a group using at most two letters,
distinct edge-letter symbols are independent. A single such
symbol appears twice only if both vertex positions for that
edge and letter occur in the group. For any partition into
groups of size at most \(m\), there are at most \(O(km)\)
potential repetitions in total, since each repetition uses
a distinct pair of positions in one group. The sum of group
frequency deficits from these repetitions is
\(O(km\log_X q)+o_X(1)\). Thus it is \(o(\lambda)\) when
\(m/k\to0\). Assigning evenly spaced scaled degree labels in
a compact positive interval gives the same conclusion for
groups supported in boundedly many shrinking degree windows,
which then have size \(O(1+k\rho_\lambda)\) for window width
\(\rho_\lambda\). The labels are bookkeeping for the model.

In contrast, a record on one letter choosing the primes at
both endpoints of each of \(s\) prescribed edges has logarithmic
size \(2s\log q+O_k(\log\log q)\) and frequency exactly \(q^{-s}\).
Taking \(s/P\) bounded away from zero gives the inverse square
root of the record size on the exponent scale. Thus near-reciprocal
frequencies in all the indicated local groups do not force
reciprocal frequencies for arbitrary large records. This model
has squarefree coprime supports but does not impose the additive
equation or realize the abstract residue maps by the prime
values chosen for the labels.

## The arithmetic requirement

These frequency estimates remove a combinatorial obstruction to
choosing common-modulus rows with new prime supports. They give lower
bounds on fibre sizes and restrictions on possible intersections.
An upper bound on such fibres still requires the additive equation and
the numerical weight of a common divisor. In a diffuse profile a
record with fixed positive raw mass can occupy a number of separated
windows growing in the budget limit. Reciprocal frequencies within
individual windows do not impose reciprocal joint frequencies for
such a record.
