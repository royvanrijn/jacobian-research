# Sharing opposite dictionaries across CRT packet records

Packets on the same summand use witnesses from one opposite-letter
dictionary. The witnesses can share a numerical interval and a residue
modulo a fixed full divisor of the marked summand. Counting the common
cell once improves the independent residue charges. An explicit primitive
example shows why a common dictionary without restrictions on its cell
sizes does not suffice.

## Data and shared dictionary bound

Use the complete-target bin and the admissible full packet records of
(TC1)–(TC3), with \(n=X^L\), primitive \(A+B=C\), \(X/2<C\le X\),
and raw radical caps \(n^{r_A+\tau},n^{r_B+\tau},n^{r_C+\tau}\).
The count can retain the original joint radical budget and the exact
cutoff on the least-repaired target. At fixed sufficiently small \(L\)
take the height limit first. Denote the marked summand by \(W\) and the
opposite summand by \(V\). In this note all \(\ell\) marks are on \(W\).
Use complete target multiplicities in
\([n^{u-\tau},n^{u+\tau}]\) in the chart under consideration.
Let the complete \(V\)-dictionary \(\mathcal V\) have size at most
\(n^{\alpha+\tau}\), and assume \(W\le X^{h+\tau}\), where
\(0\le h\le1\). The \(i\)-th packet \(F_i\) satisfies
\[
 n^{t_i-\tau}\le\operatorname{rad}(F_i)\le n^{t_i+\tau},
 \qquad F_i\ge X^{w_i-\tau},\qquad
 Lv_p(F_i)\le S_i\quad(p\mid F_i).
\]
The packets are full on their selected supports, and their groups and
complements satisfy the admissibility in (TC2). The different marks
can be supplied by different original rows on the same target.

Suppose a fixed divisor \(D\mid W\) is present in the chart. It is full
on its selected prime support and all marked packets are coprime to
\(D\). Assume
\[
 n^{d-\tau}\le\operatorname{rad}(D)\le n^{d+\tau},
 \qquad D\ge X^{w_D-\tau}.
\]
For the bound from a common residue below require \(w_D>0\) and
\(\lfloor\kappa/L\rfloor\)-fullness for a fixed \(\kappa>0\).
One may instead use \(D=1\), \(d=w_D=0\) and keep only the interval
information. The fixed coefficient \(D\) has zero dictionary cost
within the chart. The radical cap for \(W/D\) is at most
\(n^{r_W-d+2\tau}\), by the assumption that \(D\) is full on its
support. If one starts with a dictionary of choices for \(D\), add
its record exponent when summing the estimates over charts.
The multiplicity bin is imposed within each such chart.
[Admissible packet records and fixed-target bins (TC1–TC2)](ABC_TARGET_PACKET_CRT_MOMENTS.md#divisor-records-on-a-complete-target-bin)

Put
\[
 \beta=r_W-d-2u,\qquad
 \gamma_i=\min(t_i,\beta),\qquad
 w_{\mathrm{tot}}=w_D+\sum_{i=1}^{\ell} w_i
                        -\sum_{i=2}^{\ell} S_i\zeta_i .
\]
Use \(\beta\ge0\). If it is negative with a fixed margin,
(FL1) applied with the support of \(D\) already excludes the
multiplicity bin for small enough tolerances. Fix \(\zeta_i>\beta\)
at each nontrivial step \(i\ge2\); a trivial mark uses
\(t_i=w_i=S_i=0\). Assume the strict numerical margin
\[
                     w_{\mathrm{tot}}>h .                     \tag{SW1}
\]
In particular the total weight has a fixed positive lower bound.
One can also use the smaller margin obtained by applying (TC4)
to the final compatible residue cells. The numerical cutoff for
that version depends on the chosen counting allowance.

Set \(H=\lceil X^{h+\tau}\rceil\), and let \(I_j\) be the
successive integer intervals of length \(H\). For each unit
residue \(a\bmod D\) put
\[
   \mathcal V_{a,j}
       =\{v\in\mathcal V:\ v\equiv a\bmod D,\
                              v\in I_{j-1}\cup I_j\},
   \qquad
   \mathcal B_\ell=\sum_{a,j}|\mathcal V_{a,j}|^\ell .
\]
For \(D=1\) take one residue label. Every complete \(V\) is
in at most two of these cells and is coprime to \(D\) by
primitivity. For any fixed \(\eta>0\), by taking \(\tau,L\)
sufficiently small with the numerical margins fixed, the
number of occurring complete targets is at most
\[
 O\!\left(\mathcal B_\ell\,
        n^{\sum_i[\gamma_i-(u-t_i)_+]+\eta} X^\varepsilon\right).
                                                               \tag{SW2}
\]
In particular if every cell has at most \(n^{\alpha_{\rm loc}+o(1)}\)
values, uniformly with the usual inner height allowances, the
target exponent before arbitrary allowances is at most
\[
      \alpha+(\ell-1)\alpha_{\rm loc}
                   +\sum_i\bigl[\gamma_i-(u-t_i)_+\bigr] .
                                                               \tag{SW3}
\]
The same bound with an additional \(ku\) gives the same-target
scalar \(k\)-fold image in (SR3), retaining at most
\(\lceil n^u\rceil\) witnesses per represented radical. Thus in
the balanced bin an exponent in (SW3) strictly below \(1/3-u\),
with chart cost included when necessary, gives the corresponding
row-count and image saving.

### Extracting a packet list

Apply the proof of (TC3) inside the chart. Charge the remaining
original roots of \(W\) after removing the fixed full \(D\).
The raw complement cap used in the heavy-record count for packet
\(i\) is \(r_W-d-t_i+O(\tau)\). This is why the parameter in that
argument is now \(\beta=r_W-d-2u\). The estimate (FS1) is uniform
in the fixed coefficient \(D\). The target radical cap still
supplies its high full modulus in this estimate. In the overlap
test for the next packet take the prescribed pool to be the
support of \(D\) together with that of the earlier packets.
Since the new packet is coprime to \(D\), an overlap of normalized
radical \(\zeta_i\) with the earlier packets reduces the
outside-pool cap on \(W\) to \(r_W-d-\zeta_i+O(\tau)\).
The bound (FL1) discards these rows with exponent strictly
below \(u\).

The resulting fixed list \(F_1,\ldots,F_\ell\) is common to a
subset of complete targets of cardinality at least the total
target count times \(n^{-\sum_i\gamma_i-o(1)}\), up to the
prescribed fixed exponent allowances. Each of these targets
has at least \(K_i=n^{(u-t_i)_+-o(1)}\) distinct original
opposite witnesses from the record for \(F_i\). These lower
multiplicities refer to the original chart witnesses and retain
the original caps on all three letters. The least common
multiple of the \(F_i\) loses at most \(S_i\zeta_i\) in
numerical exponent at each overlap step. Because
\((D,F_1\cdots F_\ell)=1\), appending \(D\) supplies its full
weight. In particular
\[
            Q=\operatorname{lcm}(D,F_1,\ldots,F_\ell)
                    \ge X^{w_{\mathrm{tot}}-O_\ell(\tau)}.
\]
These are the same conditional packet and overlap estimates as
in (TC2), with the fixed support paid only once.
[Conditional packet and overlap counts (TC3)](ABC_TARGET_PACKET_CRT_MOMENTS.md#many-targets-with-one-packet-list)
· [One full divisor in a fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### Counting the joint cells

If a common target lies in \(I_j\) and has residue \(a\bmod D\),
then all of its opposite witnesses for every \(i\) lie in the
same \(\mathcal V_{a,j}\), since \(0<C-V=W\le H\) and \(D\mid W\).
Every packet \(F_i\) requires a residue for \(C\) modulo \(F_i\)
that contains at least \(K_i\) elements of this one cell.
Thus there are at most \(|\mathcal V_{a,j}|/K_i\) options
for that residue in the cell. By (SW1) and the bound on \(Q\),
a compatible residue list along with \(a,j\) determines at most
one target. If one uses (TC4), all targets are units modulo \(Q\)
by original primitivity and the count at one compatible list
is at most \(n^\eta X^\varepsilon\) for the corresponding
allowance. Summing over \(a,j\) proves (SW2), after restoring
the packet-list cost. Finally
\[
 \mathcal B_\ell
     \le \bigl(\max_{a,j}|\mathcal V_{a,j}|\bigr)^{\ell-1}
                         \sum_{a,j}|\mathcal V_{a,j}|
     \le 2|\mathcal V|\bigl(\max_{a,j}|\mathcal V_{a,j}|\bigr)^{\ell-1}.
\]
This gives (SW3). The proof shares the coarse cell index.
The witnesses for distinct packets can be distinct.
[Unit progression occupancy (TC4)](ABC_TARGET_PACKET_CRT_MOMENTS.md#unit-progression-cells-with-a-large-small-radical-step)
· [Scalar image and multiplicity bins (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)

### Unconditional bounds for the local dictionary

For \(D=1\) and fixed \(h<1\), the opposite values satisfy
\(V\asymp X\) for sufficiently small \(\tau\). The relative
interval count (SI1) gives
\[
       \alpha_{\rm loc}
          \le \min\!\left(\alpha,\frac{1+h}{2}r_V\right)
\]
up to arbitrary fixed allowances, uniformly over the height
cells. One can use the actual dictionary moment instead if it
has a better bound.
[Relative short intervals (SI1)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#an-elementary-bound-in-relative-intervals)

For the positive full lock \(D\), one also has
\(\alpha_{\rm loc}\le\min(\alpha,r_V/2)\) up to the same
allowances. Here only the original single-letter radical cap,
the fixed unit congruence and the large full step enter the
local count. To see this, truncate and fix the high valuation
tail of \(V\) at arbitrarily small normalized cost, then use
dyadic singleton root caps on the remaining exact positions.
A tail factor is a unit modulo \(D\). Fourier inversion on
the unit group modulo \(D\), followed by Hölder, gives the
single-letter estimate in the proof of (FS1), with the
sum of the singleton dictionary charges at most \(r_V+\tau\)
in base \(n\). The fixed positive numerical weight of \(D\)
and its fullness provide exactly the modulus hypotheses in
that estimate. The representative-independent tail and
dyadic losses fit the arbitrary allowance.
[Single-letter moment on a full modulus (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#proof-by-concentration-of-the-free-roots)
· [Grouped character moments and Hölder (GD1)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#diagonal-moments-for-a-joint-dictionary)

For example take \(r_W=r_V=r_C=1/3\), \(u=0.155\),
\(\alpha=0.205\), \(h=0.50\), \(\ell=2\), and a fixed
full lock with \(d=0.02,\ w_D=0.10\). Suppose the two
marks are supplied by different witnesses and have
\(t_i=0.04,\ w_i=0.22,\ S_i=5.6\). Use \(\zeta_2=0.005>\beta=1/3-0.02-0.31\).
The lower numerical weight is \(0.512>h\).
The common residue estimate gives
\[
  0.205+\tfrac16-2(0.155-0.04)
            +2(\tfrac13-0.02-0.31)=0.148333\ldots
\]
for the target exponent. Adding \(u\) and a chart price
of \(0.02\) gives a row exponent of \(0.323333\ldots\).
In this example the independent opposite charges would give
target exponent \(0.186666\ldots\), which does not give
the saving after the same costs. These parameters describe
a sector with strict margins; they make no assertion of
a family attaining the endpoint.

## A primitive example with independent residue choices

Fix integers \(\ell\ge2\) and \(e\ge2\), and real
\(0<r<1,\ 0<k<\ell+r-2\). For large \(N\) choose
distinct primes \(p_i\) with \(q_i=p_i^e\asymp N\)
and put \(Q=\prod_i q_i\). Such primes exist in disjoint
fixed-ratio intervals. Set \(R=\lfloor N^r\rfloor\)
and \(K=\lfloor N^k\rfloor\). Choose for each \(i\)
exactly \(R\) integers \(j_i\in[1,2R+O(1)]\) coprime
to \(p_i\), and take residues
\[
       a_i(j_i)=(Q/q_i)j_i\bmod q_i .
\]
All \(R^\ell\) choices give distinct complete targets
\[
       c=Q+\sum_{i=1}^\ell (Q/q_i)j_i .
\]
They are units modulo every \(q_i\) and lie in a single
interval \([Q,Q+O(QN^{r-1})]\). The numerical ratios
between \(Q\) and this interval length exceed a fixed
positive power of \(N\).

In an interval immediately preceding \(Q\) of length
\(\asymp QN^{r-1}\), place in one opposite dictionary
\(K\) successive integers from each of the \(R\) classes
\(a_i(j_i)\bmod q_i\), for every \(i\). They fit because
the length divided by \(q_i\) is of order \(N^{\ell+r-2}\).
Here successive means successive in the indicated
progression. The dictionary size is between \(RK\) and
\(\ell RK\). For a given target and its selected class
at \(i\), write \(c-b=q_i t\) for these \(K\) witnesses.
Their \(t\)'s are consecutive and positive. The condition
\((t,cp_i)=1\) makes \(a=c-b\) and \(b\) coprime and
also gives \(v_{p_i}(a)=e\), so \(q_i\) is an exact
full prime power of \(a\). The number retained is
\[
       K\,\frac{\varphi(cp_i)}{cp_i}
              +O\bigl(2^{\omega(cp_i)}\bigr)=N^{k-o(1)}
\]
uniformly. The equality in exponent follows from the
elementary uniform bounds
\(2^{\omega(m)}=m^{o(1)}\) and \(\varphi(m)/m=m^{-o(1)}\)
at polynomial height. These are actual positive primitive
sums, in one short physical target interval, and the same
opposite dictionary serves all the marks. The independent
residue exponent \(\ell r\) is attained despite these properties.

One may insert a common full lock \(D=p_0^e\) with \(p_0\) comparable in size to and
different from the \(p_i\), by applying the affine change
\(c\mapsto 1+Dc,\ b\mapsto 1+Db\) before the coprimality
exclusions. In choosing the \(j_i\), also exclude the one
class modulo \(p_i\) that would make the changed target
zero modulo \(p_i\), which still leaves \(R\) choices for
large \(N\). Excluding \(p_0\) from \(t\) makes \(D\) an
exact full power of the changed summand \(a=Dq_it\).
The same inclusion–exclusion with the changed target
and \(p_0p_i\) gives \(N^{k-o(1)}\) witnesses. All
opposite values are now in one unit class modulo \(D\),
and the numerical CRT margin is unchanged by the scale
factor \(D\).

The example imposes no small radical caps on the three
original letters. It rules out a uniform improvement
based solely on sharing a dictionary, within-target
witness multiplicity, primitivity and CRT geometry.
The improvement in (SW3) uses the small-radical control
of witness dictionaries inside shared cells.
