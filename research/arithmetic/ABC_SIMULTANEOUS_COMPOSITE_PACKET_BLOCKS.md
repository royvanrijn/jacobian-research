# Simultaneous composite packet residue blocks

A near-minimum abc graph with a full divisor assigned to each
complete letter admits simultaneous residue block partitions.
The total empirical price of the marked roles determines the
scale of every part in a block. When all three roles are marked,
each block translates to a shorter integer interval carrying
the original three outside radical constraints.

## Empirical prices and one marked role

Write \(c_0=1/6\) and \(n=X^L\), with \(L\) fixed before taking
the height limit. Use primitive positive triples \(A+B=C\),
\(X/2<C\le X\), as edges. Each of the three complete dictionaries
has at most \(n^{c_0+\kappa}\) values, and every original letter
has radical at most \(n^{2c_0+\kappa}\).

For each marked role \(i\) assign one packet \(D_i(T)\mid T\)
as a function of its complete letter \(T\). Its full prime powers
have \(Lv_p(T)\le S_0\), and its radical and numerical height satisfy
\[
 n^{\beta_i}\le\operatorname{rad}(D_i)\le n^{\beta_i+\gamma_i},
 \qquad D_i\ge X^{w_i}.
 \tag{SB1}
\]
The constants satisfy \(0<\beta_i\le\beta_i+\gamma_i<c_0\),
\(w_i>0\), and \(S_0>0\). There may be one, two or three marked
roles. Label an edge with the ordered tuple of its prescribed
packets in the marked roles.

In each role bin the packets by the number of complete original
letters assigned to them, in dyadic intervals \([H_i,2H_i)\).
Some tuple of bins retains at least \(n^{2c_0-\tau}/
O((\log X)^3)\) edges if there were at least \(n^{2c_0-\tau}\)
initially. Let \(s_i=c_0-\log_n H_i\). In this bin the role
has at most \(n^{s_i+\kappa}\) packets and each packet supports
fewer than \(2n^{c_0-s_i}\) complete letters. The same bounds
hold if the bin is supplied in advance and has enough edges.
As in (AP2), fixed-support enumeration gives
\(s_i\le\beta_i+\gamma_i+\eta\). At fixed packet in role \(i\)
and fixed complete value in a different role, (FL1) counts
edges at cost \(n^{c_0-\beta_i/2+\kappa/2+\eta/10}\) up to a
constant. This gives \(s_i\ge\beta_i/2-\eta\) for small
\(\kappa,\tau\) and large height. The allowances are normalized
in base \(n\). Both estimates use \(L\) sufficiently small in
terms of the fixed tolerances, independently of height.
[Composite packet price and one-role inverse (AP2)](ABC_AGGREGATE_DEEP_PACKET_INVERSE.md#composite-packet-inverse)

Throughout the proof the ambient classes use the full complete
dictionaries in the chosen bins. Edges may be discarded while
those ambient dictionaries remain fixed. For one marked packet
\(D_i\), denote by \(x_i\) the common residue of the other
two letters, with the sign determined by the additive equation.
For example, if \(D_i=M\) is on \(C\), take
\(x_i=A\bmod M=-B\bmod M\).

For any small fixed \(\eta>0\), take \(\kappa,\tau\le\eta/1000\)
and absorb subpower costs in \(n^{\eta/100}\) at fixed \(L\).
As in (AP2), one can remove
\(O(n^{2c_0-\eta/2})\) edges so that for every marked role,
packet and occupied \(x_i\), the two unmarked ambient classes
have sizes in
\([n^{c_0-s_i-\eta},n^{c_0-s_i+\eta}]\).
For clarity, the one-side upper test charges each oversized
class to distinct pairs in it with weight \(O(n^{-\eta})\),
using the \(2n^{c_0-s_i}\) target capacity. A nonzero pair
difference has only a divisor-many packet labels. The lower
test sums over disjoint classes on the opposite side at fixed
packet. This argument is valid on the graph already pruned
in other marked roles, with no recomputation of prices.
In particular the number of residues \(x_i\) still occupied
at any packet is at most \(n^{s_i+\kappa+\eta}\).

## Simultaneous marking theorem

**SB2.** Suppose that the bins selected in (SB1) have
\[
                  d=\sum_i s_i\le c_0-\zeta
\]
for a fixed \(\zeta>0\). Fix sufficiently small \(\eta>0\)
in terms of the packet parameters and \(\zeta\).
After the tolerance and limit choices above, there is a subgraph
with at least \(n^{2c_0-\tau-2\eta/100}\) edges with the
following properties.

For each surviving tuple of marked packets retain only the blocks
specified by the simultaneous residues \(x_i\).
Each occupied block has between \(n^{2(c_0-d)-6\eta}\)
and \(n^{2(c_0-d)+16\eta}\) edges. Each original role has
between \(n^{c_0-d-14\eta}\) and \(n^{c_0-d+8\eta}\)
different values in the block. Each surviving packet tuple
contains at least \(n^{d-17\eta}\) occupied blocks and
at most \(n^{d+3\kappa+3\eta}\). In addition, for two
marked roles the upper bounds on the three ambient classes
are \(n^{c_0-d+5\eta}\).

### Two marked roles

Use packets \(F\) on \(A\) and \(M\) on \(C\), with prices
\(s_A,s_C\), and set \(d=s_A+s_C\). The residue coordinates
can be taken as
\[
             x_F=B\bmod F=C\bmod F,\qquad
             x_M=A\bmod M=-B\bmod M .
\]
At fixed \(F,M\) an ambient \(A\)-class depends only on \(x_M\)
and a \(C\)-class only on \(x_F\). An ambient \(B\)-class
depends on both. Write \(a,b,c\) for the respective sizes.
The block edge count is at most \(\min(ab,bc,ac)\) because
two complete original letters determine each edge.

First count too large \(A\)-classes with \(a>n^{c_0-d+3\eta}\).
For each such class, sum its edges over all \(x_F\). At a
fixed \(A\) there are fewer than \(2n^{c_0-s_C}\) edges
with the packet \(M\), since each complete \(C\) determines
the edge. Thus the count in the class is at most
\(O(n^{s_A-3\eta}\binom a2)\). For each distinct pair
\(A,A'\) assigned to the same \(F\), every \(M\) counted
divides \(A-A'\), and the residue \(x_M\) is unique for this
pair and \(M\). The sum of pair counts before choosing \(M\)
is at most
\[
           \sum_F |\mathcal A_F|^2
              \le 2n^{2c_0-s_A+\kappa},
\]
where \(\mathcal A_F\) is the ambient \(A\)-dictionary assigned
to \(F\). The ordinary divisor bound therefore makes the
oversized edge count \(O(n^{2c_0+\kappa-3\eta+\eta/100})\).
The same estimate holds for too large \(C\)-classes.

After those deletions a \(B\)-class of size
\(b>n^{c_0-d+5\eta}\) contributes at most
\(O(n^{-2\eta}\binom b2)\), by the \(A\)-side bound.
For a distinct pair \(B,B'\), both \(F\) and \(M\) divide
\(B-B'\), so each pair can enter for at most
\(\tau_{\rm div}(|B-B'|)^2\) packet pairs. This again costs
\(O(n^{2c_0+2\kappa-2\eta+\eta/100})\) edges in total.

At fixed \(F,M\), if \(a<n^{c_0-d-3\eta}\), sum \(ab\) over
its blocks. The \(B\)-classes are disjoint, so this costs
at most \(n^{2c_0-d-3\eta+\kappa}\) edges. The same bound
works for small \(C\)-classes. If \(b<n^{c_0-d-3\eta}\),
sum \(ab\) using
\[
     \sum_{x_F,x_M\ {\rm occupied}} a_{x_M}
          \le n^{s_A+\kappa+\eta}|\mathcal A_F|
          \le 2n^{c_0+\kappa+\eta}.
\]
This costs \(O(n^{2c_0-d-2\eta+\kappa})\) edges per packet pair.
Sum over the at most \(n^{d+2\kappa}\) packet pairs.

The same proof works for every pair of marked roles with
the appropriate signs. When three roles are marked, perform
these deletions for all three pairs, recording the upper
bounds in each pair graph before continuing. All such upper
bounds remain valid after more edges are removed.

### Three marked roles

Use \(F,G,M\) on \(A,B,C\), respectively. At a fixed packet
tuple the residue coordinates are
\[
 \begin{split}
             x_F&=B\bmod F=C\bmod F,\\
             x_G&=A\bmod G=C\bmod G,\\
             x_M&=A\bmod M=-B\bmod M .
 \end{split}
\]
An \(A\)-class depends on \((x_G,x_M)\), a \(B\)-class on
\((x_F,x_M)\), and a \(C\)-class on \((x_F,x_G)\). Put
\(d=s_A+s_B+s_C\).

For a fixed \(A\)-class of size
\(a>n^{c_0-d+8\eta}\), sum all its edges over \(x_F\).
With \(A\) fixed, these edges lie in one of the two-role
blocks at \(G,M\), with its two residues fixed by \(A\).
The degree at this \(A\) is at most
\(n^{c_0-s_B-s_C+5\eta}\) by the two-role bounds.
The sum in the class is therefore
\(O(n^{s_A-3\eta}\binom a2)\).
An unordered pair of distinct \(A\)'s with packet \(F\)
enters for at most \(\tau_{\rm div}(|A-A'|)^2\) choices of
\(G,M\), because both packets divide their difference and
their residue pair is fixed. Sum over \(\mathcal A_F\)
and then \(F\) as before. This gives
\(O(n^{2c_0+\kappa-3\eta+\eta/100})\) exceptional edges.
Apply the same test to \(B,C\).

To bound too small \(A\)-classes, those of size
\(a<n^{c_0-d-4\eta}\), use \(e\le ac\) at a full residue
triple. At fixed \(F,G,M\), sum over all allowed residues.
There are at most \(n^{s_C+\kappa+\eta}\) settings of
\(x_M\) by the one-role cleanup, and the \(C\)-classes
are disjoint as \((x_F,x_G)\) vary. The edge loss from
small \(A\)-classes in this tuple is therefore at most
\[
        n^{c_0-d-4\eta}
        n^{s_C+\kappa+\eta}|\mathcal C_M|
           \le 2n^{2c_0-d-3\eta+\kappa}.
\]
For small \(B\)-classes use \(e\le ab\), sum over \(x_F\)
and use disjoint \(A\)-classes. For small \(C\)-classes
use \(e\le bc\), sum over \(x_G\) and use disjoint \(B\)-classes.
The same bound holds in each case. Sum over the at most
\(n^{d+3\kappa}\) tuples.

### Edge counts and actual parts

For either number of marked roles, the one-role residue
bounds give at most \(n^{d+3\kappa+3\eta}\) occupied
blocks per tuple. There are at most \(n^{d+3\kappa}\) tuples.
Discard blocks with fewer than \(n^{2(c_0-d)-6\eta}\)
edges. The loss is at most \(n^{2c_0+6\kappa-3\eta}\).
Discard also tuples with fewer than \(n^{2c_0-d-\eta}\)
surviving edges, at loss at most \(n^{2c_0+3\kappa-\eta}\).
The same conclusion for one marked role follows by applying
these block and tuple tests to the initial one-role cleanup.

All losses are power smaller than the binned edge count.
Each block has at most \(n^{2(c_0-d)+16\eta}\) edges by
the ambient class upper bounds. Every complete value in
the block has degree at most \(n^{c_0-d+8\eta}\),
again because two complete letters determine an edge.
Dividing gives the lower bound on the actual part sizes.
Finally, divide the tuple's edge lower bound by the block
upper bound to obtain \(n^{d-17\eta}\) occupied blocks
per surviving tuple. This proves (SB2).

## Simultaneous lists and the shortened abc lattice

**SB3.** On a graph as in (SB2), fix a role \(A\) for
which the opposite roles have packets \(G,M\), and an
integer \(k\ge1\). Set \(d'=s_B+s_C\).
For any fixed \(0<\xi<\min(\beta_B,\beta_C)\) and
any prescribed pools \(Q_A\le X^{B_*}\) of fixed
polynomial height, some fixed ordered list of \(k\)
packet pairs is incident to \(n^{c_0-O(k\eta)}\)
different complete \(A\)-values. At each such \(A\),
the two packets in any pair occur on a common edge,
and the union of their prime supports intersects the
pool and all earlier pairs in radical at most \(n^\xi\).
The same assertion holds with any role fixed.

To prove this, discard \(A\)-values of degree less than
\(n^{c_0-\eta/4}\). A fixed \(A\) has at most
\(n^{c_0-d'+5\eta}\) eligible edges per packet pair:
all such edges are in a single two-marked-role block
at \(G,M\). After a list prefix has been chosen, its
packets and \(Q_A\) form a pool at most \(X^{B_*+2k}\).
A pair with radical overlap \(>n^\xi\) puts mass at
least \(\xi/2\) on this pool in one of the two free
letters. Apply (FL1) with fixed complete \(A\) and this
pool on the corresponding free letter. Their total is
\(O(n^{c_0-\xi/4+\kappa/2+\eta/10})\), uniformly
in the prefix. First fix \(\xi,k,B_*\), then take
\(\eta\) small compared with \(\xi\). At least half
the incident edges are still eligible, giving
\(\gg n^{d'-6\eta}\) choices per step. Average the lists
from all high-degree vertices over the at most
\(n^{k(d'+2\kappa)}\) possibilities to get the assertion.
In particular the retained least common multiple for
each list after stripping \(Q_A\) has weight at least
\(k(w_B+w_C-S_0\xi)\), since the two same-edge
packets are coprime.

At any surviving full block of three marked roles, put
\(D=FGM\). The packets are pairwise coprime because
they occur on a common primitive edge. Thus each role
lies in one residue class modulo \(D\). Choose an
actual base triple \((A_0,B_0,C_0)\) in the block and write
\[
       A=A_0+D a,\qquad B=B_0+D b,\qquad C=C_0+D c.
\]
The three respective parameters range over intervals of
length at most \(X/D\) and satisfy \(a+b=c\).
Let \(R_U=\operatorname{rad}(D_U)\) for the packet on
each role. Because the packets use exact full powers,
\[
 \operatorname{rad}(A/F)\le n^{2c_0+\kappa}/R_A,\quad
 \operatorname{rad}(B/G)\le n^{2c_0+\kappa}/R_B,\quad
 \operatorname{rad}(C/M)\le n^{2c_0+\kappa}/R_C .
\]
For example \(A/F=A_0/F+GM a\), and the other two
quotients are corresponding linear forms of \(b,c\).
The surviving block contains a near-full additive graph
on its three short parameter dictionaries, and all
three quotient radical constraints hold on that graph.

A two-role simultaneous list carries the sum of the two
numerical weights, but its residue price is at most the
sum \(s_B+s_C\), before the small allowances. A direct
CRT height exclusion has the same weight-to-price
threshold as the one-role method (AP4).

## A sufficient count on the short cells

The total number of packet tuples and occupied residue triples
in (SB2) is at most \(n^{2d+6\kappa+3\eta}\). Consequently
a bound \(O(n^{2(c_0-d)-\epsilon_*})\), uniform in all of
these cells for a fixed \(\epsilon_*>0\), would save on the
whole selected bin. The count on each cell may use all three
linear forms on the short parameters together. It would also
suffice to bound the number of possible unit quotients in just
one of its three affine progressions by \(O(n^{c_0-d-\epsilon_*})\),
uniformly with the choice of role allowed to vary by cell.
For example, the progression on \(A\) is
\[
  0<u\le X/F,\qquad u\equiv A_0/F\pmod{GM},\qquad
  \gcd(u,FGM)=1,\qquad
  \operatorname{rad}(u)\le n^{2c_0+\kappa}/R_A .
\]
Indeed, multiply that bound by the other role's
\(n^{c_0-d+8\eta}\) ambient size and use pair determination.
The small tolerances can be chosen after the saving
\(\epsilon_*\). This isolates a uniform short-progression
or joint affine count as a sufficient additional input.

Gallagher's pair-difference argument underlies the class
pruning. Walsh's inverse
sieve instead has reductions of a common set on many
prime moduli. The bounds above concern different
edge-selected subsets in composite moduli, so that
hypothesis is a separate requirement.
[Gallagher, A larger sieve](https://doi.org/10.4064/aa-18-1-77-81)
· [Walsh, The inverse sieve problem in high dimensions](https://arxiv.org/abs/1105.1551)
