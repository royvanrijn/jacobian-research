# Incidence models for marked deep target layers

Near-minimum dictionary exponents, balanced local components, and a
bounded number of common marked powers on each pair are compatible with
exact integer additivity at the endpoint height scale. The construction
below uses target-dependent formal marks. Those marks are prescribed
incidences rather than divisors of the integers. It isolates an
arithmetic requirement in an inverse argument for deep target layers.

## Additive induced-matching layers

Let \(P\) tend to infinity and choose a three-term-progression-free
\(S\subseteq[1,P/4]\) with \(|S|=P^{1-o(1)}\) by Behrend's construction.
On the two copies of \([1,P]\), retain the pair \((i,j)\) if \(j-i\in S\),
and color it by \(i+j\). The graph has \(P^{2-o(1)}\) edges and \(O(P)\)
colors. Each color is an induced matching. Indeed, two distinct edges
of color \(l\) have differences \(l-2i\) and \(l-2i'\). Either cross pair
has difference \(l-i-i'\), which is the midpoint of these two distinct
members of \(S\). Thus neither cross pair is present.
[Behrend's progression-free construction](https://pmc.ncbi.nlm.nih.gov/articles/PMC1078964/)

Replace each base vertex by \(H\) labels with an auxiliary digit
\(0\le u<H\), and every base edge by a complete \(H\times H\) block.
The labels \(4Hi+u\) and \(4Hj+v\) have exact additive targets with coarse
digit \(i+j\) and fine digit \(u+v\). The complete dictionaries all have
size at most \(O(PH)\), and the graph has \((PH)^2 P^{-o(1)}\) edges.
For each color, its layer consists of the separate complete blocks
of that induced matching. Every actual edge lies in the product of the
two incident vertex sets of only its own color. The primitive realization
by fixed shifts described below can remove a fixed small fraction of
edges while preserving these exponent counts.

This construction has the same same-side partitions on the blocks of
all colors. Two distinct vertices over a common base vertex therefore
lie together in blocks of \(P^{1-o(1)}\) colors on the retained base
vertices of typical degree. Assigning distinct primes of exponent-scale
size \(P=n^\beta\) to all these colors, with depths at least a fixed
positive multiple of \(1/L\), would force the difference of those two
vertices to have superpolynomial height by (GL1). Thus the fixed
partitions of the induced-matching blow-up already fail a necessary
pairwise numerical constraint for actual deep congruences.
[Path congruences (GL1)](ABC_GRAPH_PATH_DIVISOR_LOCKS.md#exact-projection-of-overlapping-edge-records)

## Varying partitions with bounded pair agreement

Put \(c_0=1/6\). Fix \(0<\beta<c_0\) and \(0<\kappa<6\), then take a
sufficiently small fixed \(L>0\). Use \(q\to\infty\) through primes,
\(X=q^{1/(\beta L)}\), \(n=X^L\), and
\[
 d=\lceil c_0/\beta\rceil,\qquad
 \vartheta=c_0/\beta-(d-1)\in(0,1],\qquad
 N=q^{d-1}\lfloor q^\vartheta\rfloor,\qquad H=N/q .
\]
The two summand parts have the same coefficient box
\(\mathcal A=\prod_{i=0}^{d-1}\{0,\ldots,s_i-1\}\), with
\(s_i=q\) for \(i<d-1\) and \(s_{d-1}=\lfloor q^\vartheta\rfloor\).
Thus \(N=n^{c_0+o(1)}\) and \(H=n^{c_0-\beta+o(1)}\).
A coefficient vector \(a\) defines the polynomial
\(f_a(z)=\sum_{i=0}^{d-1}a_i z^i\) over \(\mathbb F_q\).
For each \(t\in\mathbb F_q\), mark a pair \((a,b)\) at \(t\) if
\[
             f_a(t)+f_b(t)=0,\qquad f_a+f_b\ne0.
 \tag{DI1}
\]
Use only edges with at least one mark. Before the primitive
restriction, the layer at \(t\) is the union of \(q\) blocks
\(f_a(t)=r,\ f_b(t)=-r\) for \(r\in\mathbb F_q\). Each side of
each block has exactly \(H\) vertices, since the constant digit
ranges over the whole field. The second condition in (DI1)
deletes at most one neighbor per vertex, so every block has
at least \(H(H-1)\) edges.

An edge has at most \(d-1\) marks by the root bound for a nonzero
polynomial. Any pair of distinct vertices on one side shares
a block in at most \(d-1\) layers, by the root bound for the
difference of their two polynomials. These bounds are uniform
in the shapes of the two vertex boxes.

### Exact additive values and primitivity

Encode the coefficient vectors by
\(D(a)=\sum_{i=0}^{d-1}a_i (4q)^i\).
Use actual summand labels
\[
                   A(a)=U_0+K J_X D(a),\qquad
                   B(b)=V_0+K J_X D(b)
\]
and target \(C(a,b)=A(a)+B(b)\). Here \(K\) is a fixed product of
small primes, chosen below. Take \(J_X\) to be a power of two with
\(KJ_X\asymp_d X/q^d\) and \(KJ_X\max_a D(a)\le X/100\).
This is possible at sufficiently small fixed \(L\). In particular,
every nonzero same-side difference is at least \(KJ_X\).
There are no carries in the digit
sums. In particular a complete target determines the coefficient
vector of \(f_a+f_b\), and its marks depend on the target alone.
There are at most \(\prod_i(2s_i-1)\le2^dN\) possible targets, and
at most \(2^dH\) targets with any prescribed mark \(t\). For the
latter bound fix the nonconstant sum digits. The constant sum
digit then lies in one residue class modulo \(q\), with at most
two choices. All encoded values \(D(a)\) are \(O_d(q^d)\), which
is \(o(X)\) at sufficiently small fixed \(L\).

Take \(K\) to be the product of primes at most a large fixed \(P_0\).
Choose \(U_0,V_0\) independently and uniformly from the integers
in \([X/3,3X/8]\) congruent to \(1\bmod K\). All targets lie in
\((X/2,X]\) for large \(X\). For any fixed digit pair and prime
\(p>P_0\), each divisibility probability is at most
\(1/p+O_K(1/X)\), because \(K\) is invertible modulo \(p\).
Their joint probability is the product. The smaller primes cannot
divide either summand. Summing over \(p\le X\), the probability
of losing this row to nonprimitivity is at most
\(\sum_{\ell>P_0}\ell^{-2}+o_X(1)\), uniformly in the digits.

Choose \(P_0\) large enough that this probability is at most a small
fixed \(\epsilon\). The total pre-thinning number of marked incidences
is at least \(q^2H(H-1)=(1-o(1))N^2\) and at most \(N^2\).
There is a choice of the two shifts for which deleting all
nonprimitive edges removes at most \(\epsilon N^2\) incidences.
In particular \(\gg q^2\) marked blocks still have at least
\(3H^2/4\) edges. Each such block contains a component with both
part sizes comparable to \(H\) and at least \(H^2/2\) edges. To
verify this, its largest component on the first side must have
size at least \(3H/4\), since the total number of edges is at most
that largest size times \(H\). The other components together have
at most \(H^2/4\) edges. The chosen component therefore has at
least \(H^2/2\) edges, forcing its second side to have size at least
\(H/2\). The primitive graph as a whole has \(\gg_d N^2\) edges,
since each edge has at most \(d-1\) marks. Thus each active
complete dictionary has size \(\asymp_d N\), using the \(N\)
upper bound for the degree of any one complete label.

### Prime-scale budgets satisfied by the marks

Choose \(q\) distinct rational primes \(p_t\in(q,q(\log q)^2]\),
one per marker, and put \(e=\lfloor\kappa/L\rfloor\), taking \(L\)
small enough that \(e\ge1\). Give a mark \(t\) the formal
weight \(e\log p_t\). At fixed \(L\), all the marked primes have
normalized size \(\log_n p_t=\beta+o(1)\). The marked weight
on the large components has order \(N^2\log X\) with a positive
constant depending on \(\beta,\kappa,L\). The number of distinct
marked targets per prime is \(O_d(n^{c_0-\beta+o(1)})\);
each layer has at most \(N^2/q\) edges. The product of primes
marking one target is at most
\(n^{(d-1)\beta+o(1)}\), with \((d-1)\beta<c_0\).

For two fixed distinct labels on one side, the product of all
marked powers assigned to common blocks has numerical weight at
most \((d-1)eL\beta+o(1)\) in base \(X\). On the other hand every
nonzero same-side difference has numerical weight at least
\(1-d\beta L-o(1)\). Since \((d-1)\kappa\beta<\kappa c_0<1\),
the latter bound strictly exceeds the former for sufficiently
small fixed \(L\). Thus every individual pair passes the
numerical capacity test used by (GL5).
A fixed polynomial-height prime pool has only \(O_{\beta,L,B}(1)\)
of the assigned primes if its product is at most \(X^B\).
The fraction of marked incidence weight on such a pool is
\(O_{\beta,L,B,d}(1/q)\). This also avoids a concentration
of the formal marked supports in one fixed polynomial-height pool.

These comparisons use formal assigned prime powers. The digits were
constructed to satisfy (DI1) over a common finite field.
They impose no divisibility \(p_t^e\mid C(a,b)\), and there is
no bound on the complete radicals of the three original letters.
The model therefore separates the combinatorial and numerical
tests on local layers from the actual prime-power realizability
problem. An arithmetic inverse must use the explicit congruences
or the full radical caps to rule out this remaining incidence pattern.
[Deep difference transfer (GL5)–(GL7)](ABC_GRAPH_PATH_DIVISOR_LOCKS.md#a-difference-moment-inequality-for-the-actual-edge-layers)

## No low-degree polynomial digit lift of the marked congruences

There is an absolute \(\epsilon>0\) such that the construction
in (DI1) cannot realize all but \(\epsilon N^2\) of its marked
edge incidences by bounded-height integer congruences if
the first summand labels \(A(a)\) are given by an integer
polynomial in their digits of total degree
\(D_q=o(\min_i s_i)\). This includes arbitrary affine digit
encoders. The first labels are distinct and in \([1,X]\),
and the second summand labels may be arbitrary integers.
For each \(t\) demand \(A(a)+B(b)\equiv0\bmod p_t^e\)
on its retained marked edges. The conclusion holds for the
logarithmic height normalization above, even though the
varying partitions meet each individual pair's numerical
capacity test.

Suppose that at most \(\epsilon N^2\) marked incidences are
missing. Include the identically-zero polynomial pairs in that
loss, which adds only \(q^2H=o(N^2)\) incidences.
For all but \(O(\sqrt{\epsilon}\,q)\) of the \(t\)'s, the missing
edge count over the \(q\) complete blocks at that level is
at most \(O(\sqrt{\epsilon}\,qH^2)\).
For each such level, all but \(O(\epsilon^{1/4}q)\) of its blocks
lose at most \(O(\epsilon^{1/4}H^2)\) edges. In each remaining
block the largest component on the first side has at least
\((1-O(\epsilon^{1/4}))H\) vertices, by the bound on total
edges by the largest first-side component size times \(H\).
The union of these selected first-side components therefore
omits only \(O(\epsilon^{1/4}N)\) vertices.

Choose one of these levels with \(t\in[q/4,3q/4]\),
using the integer representative of \(t\).
On the first side the map \(a\mapsto a'\) given by
\[
                 a'_1=a_1+1,\qquad
                 a'_0=a_0-t\bmod q
\]
and all other digits unchanged, preserves \(f_a(t)\).
Use the representative \(0\le a'_0<q\), and restrict to
\(a_1+1<s_1\). The map is injective on its domain.
Both wrap and nonwrap in the constant digit occur on at
least \(N/8\) pairs for large \(q\), since \(s_1\ge2\).
For sufficiently small \(\epsilon\) there remains a pair
of each type with both endpoints in the chosen union.
Such endpoints with the same evaluation lie in the same
selected component. The path congruences force equality
of their actual first-summand labels modulo \(p_t^e\).

In fact there are \(\gg N\) remaining pairs of each type.
Consequently each of the two integer difference polynomials
\[
 \begin{split}
 \Delta^-_t(x)&=A(x_0-t,x_1+1,x_2,\ldots)-A(x),\\
 \Delta^+_t(x)&=A(x_0-t+q,x_1+1,x_2,\ldots)-A(x)
 \end{split}
\]
vanishes modulo \(p_t^e\) at \(\gg N\) points of the digit box.
If its coefficients are not all divisible by \(p_t^e\),
remove their smallest common \(p_t\)-power and reduce modulo
\(p_t\). A nonzero polynomial of total degree at most \(D_q\)
vanishes on at most \(D_q N/\min_i s_i=o(N)\) points of that
Cartesian box, since the digits are distinct in \(\mathbb F_{p_t}\).
The bound follows by induction on the number of coordinates,
separating the exceptional zeros of the leading coefficient
in the final variable. This is a contradiction. Both difference
polynomials therefore vanish coefficientwise modulo \(p_t^e\).

Subtracting the polynomials and changing variables by an
integer translation gives
\(A(x_0+q,x')-A(x_0,x')\equiv0\bmod p_t^e\) coefficientwise.
Descending on the degree in \(x_0\), the top difference
coefficient is \(j q\) times the coefficient of \(x_0^j\).
Since \(p_t>q\) and \(j\le D_q<q\), it follows that all the
coefficients of positive \(x_0\)-degree are divisible by
\(p_t^e\). This holds for \(\gg q\) distinct \(t\)'s.

The grid \(\{0,\ldots,D_q\}^d\) is inside the digit box.
Interpolation there bounds every coefficient of \(A\) by
\(X\exp(O_d(D_q\log(D_q+1)))\). For instance the coefficient norm
of the numerator of a univariate Lagrange basis polynomial is
at most \((D_q+1)!\), and its denominator has absolute value
at least one. Tensor interpolation gives the displayed bound.
Its logarithm is \(o(q\log q)\) at fixed \(L\),
whereas the product of the \(\gg q\) compulsory moduli has
logarithm \(\gg q\log q\). The positive-\(x_0\) coefficients
must be zero. This contradicts distinctness of the labels and
proves the obstruction.

Pure occupancy and pair-capacity bounds thus allow the
incidence model even though every polynomial digit encoding
of the stated complexity fails the prime-power congruences.


### Exact additive realization of the whole template

For the unthinned template, exact additivity alone forces the
first digit encoder to be affine. Suppose the target label is a
function of the componentwise integer digit sum and every pair
from (DI1) is an additive edge. Fix one digit index \(i\), and
two first-side vectors \(a,a'\) at which incrementing coordinate
\(i\) stays in the box. Choose a second-side vector \(b\) with
\(b_i\ge1\) such that both \((a,b)\) and \((a',b)\) are edges.

For \(d=2\), choose \(b_1\) to avoid the two classes
\(-a_1,-a'_1\bmod q\) and choose \(b_i\ge1\); both sum
polynomials are nonconstant linear polynomials and have a root.
For \(d\ge3\), impose \(f_b(1)=-f_a(1)\) and
\(f_b(2)=-f_{a'}(2)\). Fix digits above \(b_2\) and vary \(b_2\).
The two equations fix the unrestricted digits \(b_0,b_1\),
whose slopes as functions of \(b_2\) are \(2,-3\) in
\(\mathbb F_q\). Thus for \(i=0,1\) at most one candidate
with a given varying digit is excluded by \(b_i=0\).
For \(i=2\) choose \(b_2\ge1\), and for \(i>2\) fix
\(b_i\ge1\) in advance. Exclude the at most two \(b\)'s
giving identically zero sum polynomials. Such a choice
exists for large \(q\).

The pairs \((a+{\bf e}_i,b-{\bf e}_i)\) and
\((a'+{\bf e}_i,b-{\bf e}_i)\) have the same respective
integer digit sums as \((a,b)\) and \((a',b)\), so are
also edges with the same respective targets. Exact
additivity shows that the two first-side increments both
equal \(B(b)-B(b-{\bf e}_i)\). Each coordinate increment
is constant on its whole domain, which makes \(A\)
affine in the integer digits. Hence the congruence
obstruction applies to the unthinned template with no
polynomial encoding hypothesis, provided its target
labels are assigned by those integer digit sums.

A finite-difference argument excludes near-complete CRT
realizations with arbitrary bounded-height label maps.
[Robust evaluation-layer obstruction](ABC_EVALUATION_LAYER_CRT_OBSTRUCTION.md)
