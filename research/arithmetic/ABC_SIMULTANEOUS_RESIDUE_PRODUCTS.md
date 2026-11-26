# Residue products in simultaneous abc packet blocks

When the three full divisors are fixed, near-saturation forces
the *absolute* residue triples to occupy almost the entire
Cartesian product at the level of exponents. This gives a
pair-projection test for an arithmetic saving and supplies
nondegenerate residue boxes with a precise count. Throughout,
each occupied cell consists of primitive rows satisfying all
three original radical constraints.

## Near-full products and their pair projections

Use the three-role case of
[the simultaneous marking theorem (SB2)](ABC_SIMULTANEOUS_COMPOSITE_PACKET_BLOCKS.md#simultaneous-marking-theorem).
Let \(F,G,M\) be any surviving packet tuple, let the
three corresponding prices be \(s_1,s_2,s_3\), and set
\(d=s_1+s_2+s_3<c_0=1/6\), \(t=c_0-d\).
Use the residue coordinates
\[
 x_1=B\bmod F,\qquad x_2=A\bmod G,\qquad x_3=A\bmod M .
\]
Let \(\mathcal V_i\) be their three occupied alphabets and
\(\mathcal S\subseteq\mathcal V_1\times\mathcal V_2\times\mathcal V_3\)
their set of occupied cells. Write \(U_i=n^{s_i+\kappa+\eta}\).
The one-role pruning preceding (SB2) and its final tuple
pruning give
\[
      |\mathcal V_i|\le U_i,\qquad
      |\mathcal S|\ge n^{d-17\eta}.
\]
Thus the three alphabet sizes are \(n^{s_i+O(\kappa+\eta)}\),
obtained for the lower bound by dividing the lower bound for
\(|\mathcal S|\) by the other two alphabet upper bounds.
In particular the density in the product of the three
actual alphabets is at least \(n^{-20\eta-3\kappa}\).

**SC1.** For any subfamily \(\mathcal E'\) in this fixed tuple
let \(P_{ij}\) be the set of residue pairs \((x_i,x_j)\)
it uses. With \(k\) the remaining index,
\[
       |\mathcal E'|\le |P_{ij}|\,U_k\,n^{2t+16\eta}.
                                                               \tag{SC1}
\]
This follows by summing the cell edge cap over the at most
\(U_k\) occupied completions of a pair. In particular, a
\(\sigma\)-loss from the natural pair exponent \(s_i+s_j\)
gives a \(\sigma-(17\eta+\kappa)\)-loss from the tuple
exponent \(2c_0-d\). The pair can depend on the tuple.
The same proof applies when the compression is tested on
part of a tuple, independently of the complementary rows.

For a useful algebraic special case, give each of the
alphabets for \(i,j\) any injective representation in a
field. If the two coordinates used by \(\mathcal E'\) lie
on the zero set of a nonzero polynomial of total degree
\(D\) there, then
\[
       |P_{ij}|\le D\,(|\mathcal V_i|+|\mathcal V_j|).
\]
Indeed, regard the polynomial as a polynomial in the
second variable and choose one nonzero coefficient
polynomial in the first. Its zeros among the first
coordinates contribute at most \(D|\mathcal V_j|\) pairs,
and all other first coordinates contribute at most \(D\)
pairs each. A dictionary of \(J\) such polynomials costs
a factor \(J\). Hence even a dictionary of size
\(n^{\min(s_i,s_j)-\sigma}\), with bounded \(D\),
saves on the tuple when the tolerances are small in terms
of \(\sigma\). This is a sufficient arithmetic criterion.
The three absolute residues at distinct moduli have no
automatic common-field relation of this kind.

The projection also gives a direct role-specific count.
For example \(P_{12}\) is precisely the set of occupied
complete \(C\bmod FG\) classes, under the canonical CRT
identification. The other projections give complete \(A\bmod GM\) or
\(B\bmod FM\) classes with the indicated sign convention.
In each case its cap has exponent \(s_1+s_2\).
A uniform power saving in this class count across the
surviving tuples, with the choice of projected role
allowed to vary, suffices in place of a short-cell
representation bound. Summing (SC1) over at most
\(n^{d+3\kappa}\) tuples gives the saving on the retained
family.

## Nondegenerate residue boxes

**SC2.** Let \(m\ge2\) be fixed, and assume that each
\(s_i\) is bounded below by a common positive constant.
For any subset \(\mathcal S'\subseteq\mathcal S\) of size
at least \(n^{d-\xi}\), there are
\[
       \gg_m n^{md-O_m(\xi+\eta+\kappa)}
                                                               \tag{SC2}
\]
ordered choices of distinct \(x_{i,1},\ldots,x_{i,m}\) in
each of the three alphabets such that every one of the \(m^3\)
cells \((x_{1,a},x_{2,b},x_{3,c})\) is in \(\mathcal S'\).
Here \(\xi,\eta,\kappa\) are sufficiently small in terms
of \(m\) and \(\min_i s_i\), and \(n\) is sufficiently
large. The normalized scales and limit order are those
of (SB2).

To prove it, write \(\rho=|\mathcal S'|/\prod_i|\mathcal V_i|\).
Choose \(m\) independent uniform samples from each alphabet,
with replacement. The probability of hitting only
\(\mathcal S'\) across the whole box is at least \(\rho^{m^3}\).
In detail, first sample only one point per axis. Replacing
the first point by \(m\) independent samples replaces an
inner conditional expectation by its \(m\)-th power, so
Jensen gives a lower bound \(\rho^m\). Repeat in the
second and third axes to get \(\rho^{m^3}\).
Also
\[
 \rho\ge n^{-\xi-3(\kappa+\eta)},\qquad
 |\mathcal V_i|\ge n^{s_i-\xi-2(\kappa+\eta)}
\]
by bounding the two complementary alphabets. The chance
of a repeat in some alphabet is at most
\(\binom m2\sum_i|\mathcal V_i|^{-1}\), which is \(o(\rho^{m^3})\)
whenever
\[
 (m^3+1)\xi+(3m^3+2)(\kappa+\eta)<\min_i s_i .
\]
Multiplying the remaining
probability by \(\prod_i|\mathcal V_i|^m\) and using
\(\prod_i|\mathcal V_i|\ge n^{d-\xi}\) proves (SC2).

Every cell of such a box carries actual additive rows.
For the full \(\mathcal S\) one can choose any such row in
each cell, with at least \(n^{2t-6\eta}\) choices per cell.
If \(\mathcal S'\) comes from a further edge restriction,
the same choice statement holds with whatever lower cell
size that restriction preserves. The original letter
classes within this box have a precise overlap pattern.
Each \(A\)-class depends on indices \(2,3\), each \(B\)-class
on \(1,3\), and each \(C\)-class on \(1,2\). Distinct pairs
of those indices give disjoint classes in the corresponding
letter. A box asserts simultaneous occurrence of those
classes. Choosing one *value* from each pair class that
works on every adjacent cell is a separate equation
system. In a near-full tuple the actual common-target
differences force abundant solutions to that system.
[Shared-value extraction (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)

## Arithmetic scope

(SC1) and (SC2) start with the full numerical packet
divisors and the absolute residues, so they retain the
shortened affine cells and their three outside radical
caps. Their proofs then use only the size and incidence
properties of those cells. The Cartesian product itself
attains the pair count exponent in (SC1) and has all the
boxes in (SC2). Consequently a residue graph argument
using only these properties supplies no pair-projection
saving. The arithmetic input needed for this route is a
uniform compression of the actual two-coordinate
projections, an exclusion of a positive-power proportion
of residue boxes at the scale of (SC2), or a direct
uniform bound on the jointly constrained cell.

An exclusion phrased in terms of boxes must bound their
number by \(n^{md-\sigma}\) for a fixed \(\sigma>0\) to
contradict (SC2) as \(\xi,\eta,\kappa\) go to zero.
A bound for shared-value lifts can instead be compared
with the lower bound (LB1) from the common-difference
argument. The residue supersaturation here gives the
stated number of boxes with independent row choices.

## A shared-value lift without the outside radical constraints

**SC3.** Fix \(m\ge2\) and \(\theta>0\). Let \(F,G,M\) be
pairwise coprime, \(Q=FGM\), and suppose \(Q\operatorname{rad}(Q)
\le X^{1-\theta}\). For each of the three moduli prescribe
\(m\) distinct unit residues. For all sufficiently large
\(X\) in terms of \(m,\theta\), the whole box has an exact
additive lift with one integer value per pair class. Every
row is positive and primitive, has \(X/2<C\le X\), and
contains its prescribed packet to full depth. There are no
outside radical caps in this assertion.

For each index choose the CRT lifts
\(\mu_a,\lambda_b,\nu_c\) of the corresponding prescribed
residues on \(F,G,M\), respectively, taking zero on
each of the other two moduli. Define
\[
 A_{bc}=\lambda_b+\nu_c,\qquad
 B_{ac}=\mu_a-\nu_c,\qquad
 C_{ab}=\lambda_b+\mu_a.
\]
These have the requisite absolute residues and additive
equations. To make the packet on each role full, impose
further congruences at the primes of \(Q\). At \(p\mid F\)
take \(\lambda_b/F\equiv1,\ \nu_c/F\equiv0\bmod p\).
At \(p\mid G\) take \(\mu_a/G\equiv1,\ \nu_c/G\equiv0\bmod p\).
At \(p\mid M\) take \(\lambda_b/M\equiv1,\ \mu_a/M\equiv0\bmod p\).
Each congruence is on an integral quotient. Since the old
lifts can be changed by multiples of \(Q\), each newly
specified condition fixes one residue for that change at
\(p\). The CRT thus gives progressions for the three
groups of variables modulo \(Q\operatorname{rad}(Q)\),
using arbitrary residues on any remaining primes.
They enforce \((A_{bc}/F,F)=(B_{ac}/G,G)=(C_{ab}/M,M)=1\).
At every prime of \(Q\), the two letters on the other
roles are units.

Take \(\lambda_b,\mu_a\) in \([0.30X,0.31X]\) and \(\nu_c\)
in \([0.05X,0.06X]\). Fix \(K=100m^3+1\).
At primes \(p\le K\) outside \(Q\), further impose
\(\lambda_b\equiv1,\ \mu_a\equiv\nu_c\equiv0\bmod p\).
The enlarged progression step is at most a constant
depending on \(m\) times \(X^{1-\theta}\), so each
progression in the prescribed interval has at least
\(H\gg_m X/(Q\operatorname{rad}(Q))\gg_m X^\theta\) terms,
for a common integer \(H\).
Choose uniformly and independently from \(H\) consecutive
terms for every variable. For any given row and prime
\(K<p\le H\) outside \(Q\), condition on \(\nu_c\).
The independent choices of \(\lambda_b,\mu_a\) make the
probability that \(p\) divides both \(A_{bc}\) and \(B_{ac}\)
at most \((1/p+1/H)^2\le4/p^2\). Sum this over \(p>K\).
For the primes \(p>H\) outside \(Q\), also condition on
\(\lambda_b\). There are at most \(\log X/\log H\) such
prime divisors of \(A_{bc}\). Each divides \(B_{ac}\) with
conditional probability at most \(2/H\). The total failure
probability on any row is therefore at most
\[
       \frac{4}{K-1}+\frac{2\log X}{H\log H}.
\]
A union bound over the fixed \(m^3\) rows is less than
\(1/2\) when \(X\) is large enough. There are thus at least
\(H^{3m}/2\) successful parameter choices. An ordered array
determines all parameters once \(\lambda_1\) is given, by
the formulas for \(A_{1c},C_{a1},A_{b1}\). Each array has
at most \(H\) parameter choices. This gives at least
\(H^{3m-1}/2\) distinct ordered lifts. The intervals give
the positive height shell and the residue choices
distinguish all the pair classes within a role.

This example permits arbitrary unit residue boxes and actual
shared-value lifts, even with the prescribed full divisors
and numerical height. It leaves the large primes of the three
outside cofactors unrestricted.
