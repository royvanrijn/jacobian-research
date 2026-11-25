
# Numerical order weights of secants in affordable fibres

A positive-exponent affordable fibre cannot be concentrated
in a common projective residue class whose modulus nearly
reaches the height of its cross determinants. This holds for
an arbitrary modulus. Angular thinning then gives the
cross determinants an effective numerical scale, and
counting large full divisors at that scale forces positive
numerical weight at bounded scaled orders on almost every
normalized secant of the retained subfamily.

## Setup and the projective packing gap

Use the fixed-budget height limit and then the small-budget
limit of the balanced boxes, with \(n=X^L\). On a fixed
affordable record write
\[
 U_i=F u_i,\qquad V_i=G v_i,\qquad
 T_i=\sigma F u_i+\tau G v_i=D w_i .
\]
The positive \(F,G,D\) are fixed full factors, \(D=1\) being
allowed. The signs are fixed, the complete letters are
pairwise coprime, and the targets are distinct. At least
one of \(F,G\) has numerical weight bounded away from zero
by the fixed affordable margin. This lower bound is kept
when full factors are enlarged. Set
\[
 \Delta_{ij}=u_i v_j-u_j v_i,\quad
 g_u=\gcd(u_i,u_j),\quad g_v=\gcd(v_i,v_j),\quad
 g_w=\gcd(w_i,w_j),\qquad
 H_{ij}=\frac{|\Delta_{ij}|}{Dg_u g_v g_w}.
\]
These are nonzero determinants for distinct targets and
\(H_{ij}\) is integral by (DC1). The cofactor pairs are
primitive. One may refine their numerical heights to
dyadic boxes at inner subpower cost.
[Full-support peeling and quotient determinants](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md)

Fix \(\rho>0\) and a positive lower margin for the numerical
weight of one of \(F,G\). There is a fixed \(\nu>0\) with
the following property in all sufficiently late balanced
endpoint boxes. Suppose
\[
               |\Delta_{ij}|/D\le X^{h+o_X(1)}
\]
for every pair in a family, with \(h\) bounded. For any
positive integer \(M\) with \(\log_X M\ge h-\nu\), each
projective residue class of this family modulo \(DM\)
has normalized exponent less than \(\rho\), with the
usual inner height allowances. The class is on the primitive
cofactor pairs \((u_i,v_i)\). There is no fullness or
radical cap on \(M\), and \(M\) and the class can vary with
\(X\). The assertion is uniform over the fixed factors
and the class, with their usual complete-letter height
and radical caps.

To prove this, suppose no such gap exists. Use a diagonal
sequence of violating classes with exponent at least \(\rho\)
and with \(h-\log_X M\to0\), taking every fixed-budget
height limit first. We may assume
\(DM\) has polynomial height: if it is larger than the
largest \(|\Delta_{ij}|\), its class already has at most
one point. Within each class strip and fix all full parts
of the complete letters supported on \(DM\), at subpower
cost in the inner limit. This uses the uniform fixed-support
estimate and includes the original fixed full factors if
their supports meet \(DM\). Keep the *original* cofactor
coordinates in the determinant inequalities; the newly
fixed parts serve only as records.

If any remaining degree measure is nonatomic on a
positive-mass set, choose a compact positive degree band
in it and fix a full block \(P\) of arbitrarily small
normalized root cost, less than \(\rho/4\), and fixed
positive numerical weight. This is the existing
whole-position selection: remove the large atomic
windows, then select enough small original coordinates
in continuity bands. Its prime support avoids \(DM\).
For any two points in the resulting record, \(P\mid
\Delta_{ij}\). If the block is on one of the crossed
letters, its primes are also outside the fixed full
factor of that letter and the assertion follows directly
in its residual coordinate. If the block is on the
target, it follows from the two determinant identities
with \(F,G\) in (DC1). Since \(DM\mid\Delta_{ij}\) and
\((P,DM)=1\), their product divides \(\Delta_{ij}\).
It exceeds the determinant cap for sufficiently late
budgets because the original slack tends to zero.
Two primitive positive cofactor pairs with zero
determinant are identical. This leaves at most one
point per block record, against the hypothesized exponent.

If there is no retained nonatomic measure, all residual
spectra are purely atomic. Use a positive-weight full
modulus from an original fixed \(F\) or \(G\); discarding
orders below \(\kappa_0/L\) for a small fixed
\(\kappa_0>0\) leaves positive weight because its
normalized radical cap is bounded. The two opposite
letters have a prescribed unit ratio modulo this
divisor. Their newly fixed full parts from the support
record enter the fixed unit twist. Keep finitely many
positive atom windows and a zero window, and freeze the
remaining tails at arbitrarily small root cost, including
any escaped numerical moment. Group the positive
windows and put their discrepancies and zero windows
in the wider factor of numerical height \(X^{o(1)}\).
The existing amplified full-modulus count releases the
exact ratio at arbitrarily small normalized cost.
Primitivity recovers the original triple. This is
again a contradiction and proves the projective gap.
[Amplified full-modulus count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

In particular \(M=1\) shows that an arbitrarily small
numerical diameter of the divided determinants is already
incompatible with a positive-exponent fibre. If one works
with count upper bounds with \(\epsilon\) allowances,
the same diagonal proof gives the bound \(n^{\rho+o(1)}\)
uniformly for each class after choosing its strict numerical
margin. This is the version used in the summation below.

## Fixing an effective numerical scale

Fix a lower bound \(r_0>0\) for the fibre exponent and
an arbitrarily small allowed exponent loss
\(0<\epsilon<r_0/8\). Choose the projective gap \(\nu\)
above for class exponent \(\rho=r_0/4\), and fix
\(\kappa>16/r_0\). Only after these choices take
\(0<\theta<\nu/8\) and \(0<\eta<\nu/(32(1+\kappa))\).

First apply the existing numerical and radical support
peeling with threshold \(\eta\), allowing loss less than
\(\epsilon/3\). This can enlarge \(F,G,D\) at the
same fixed record. On the resulting subfamily there is
a fixed \(\delta_1>0\) such that all but an
\(n^{-\delta_1+o(1)}\) fraction of pairs have total
residual shared numerical weight at most \(\eta\):
\[
                    \log_X(g_u g_v g_w)\le\eta .
\]
The fraction statement follows by taking each
\(u_i v_i w_i\) as the fixed prime pool for its
neighbours; the polynomial-height pool bound can be
any fixed constant greater than three. The corresponding
radical bound also holds but is not needed here.
[Support peeling (DC3)](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#removing-concentrated-support-from-a-large-fibre)

Refine the residual \(v\)'s to one dyadic interval at
subpower cost. Now angularly thin without changing the
fixed factors. Start with any numerical cap \(h\le3\)
for \(|\Delta_{ij}|/D\), for instance \(h=3\). If a
subfamily of \(N\) points has at least \(N^2 n^{-\delta_2}\)
ordered distinct pairs with
\[
                         |\Delta_{ij}|/D<X^{h-\theta},
\]
one anchor has at least \(N n^{-\delta_2}\) such
neighbours. Their slopes \(u/v\) all lie within
\(O(D X^{h-\theta}/V^2)\) of its slope, where \(V\)
is the lower endpoint of the dyadic \(v\)-range.
Thus for any pair of these neighbours the divided
determinant is \(O(X^{h-\theta})\). Keep the neighbours
and reduce the cap by \(\theta/2\), absorbing the
absolute constant in the inner limit. There are at most
\(K=\lceil 6/\theta\rceil+1\) such steps. In fact the
projective gap with \(M=1\) excludes a cap at most
\(\nu\) on any surviving fibre of exponent greater
than \(\rho\).

Choose \(\delta_2>0\) small enough that
\[
                 K\delta_2<\min(\epsilon/3,\delta_1/4).
\]
This is a power-scale thinning, with no positive-density
assumption. It retains exponent at least \(r_0-\epsilon\).
The previously exceptional support pairs are still a
vanishing fraction: their old number divided by the
square of the new count grows by at most
\(n^{2K\delta_2+o(1)}\). On the final subfamily, all
pairs have divided determinant at most \(X^{h+o_X(1)}\),
and all but a vanishing fraction have
\[
          \log_X(|\Delta_{ij}|/D)\ge h-\theta,\qquad
          \log_X H_{ij}\ge h-\theta-\eta .
\]
Here \(h\) is fixed by the final numerical cap of the
current subfamily, not assumed to equal the sum of the
full cofactor height exponents. The dyadic refinement,
support thinning and angular thinning all take place
inside the fixed-budget height limit.

## Positive weight at bounded scaled orders

Use the fixed \(\kappa>16/r_0\) above and take an
integer \(k\) with \(kL\to\kappa\) in the outer limit. For a
positive integer \(H\), write \(Q_k(H)\) for its full
part at prime powers of order at least \(k\).
The number of \(k\)-full integers at most \(Y\) is
\(O_k(Y^{1/k})\). For example, express each as
\(a^k\prod_{j=k+1}^{2k-1}b_j^j\) by separating the
nonzero valuation remainders, then sum the convergent
products \(\sum b_j^{-j/k}\) after counting \(a\).
All \(k\)-dependent constants are fixed before the height
limit.
[Full-divisor counting precedent](ABC_POWERED_PROGRESSION_OCCUPANCY.md#stopping-a-full-divisor-above-the-block-width)

Therefore there are at most \(n^{h/\kappa+o(1)}\)
possibilities for \(M=Q_k(H_{ij})\) on the final
subfamily. Suppose \(\log_X M\ge h-\nu\). Then
\(DM\mid\Delta_{ij}\), even if \(D\) and \(M\) share
primes, by the exact definition of \(H_{ij}\). Thus
\(i,j\) occupy the same projective class modulo \(DM\).
Conversely, equality of primitive projective classes
modulo a modulus is equivalent to vanishing of the
determinant modulo it. To see this over each prime
power, choose a unit coordinate of one primitive vector
and solve for the scalar. The scalar is a unit by
primitivity of the second vector, and CRT applies.

By the uniform class bound, each such \(M\) contributes
at most \(N n^{r_0/4+o(1)}\) ordered pairs on the
final \(N\)-point set, summing the squares of the
class sizes. Here the class bound may be applied to
subfamilies of this same set. Summing over \(M\)
gives at most \(N n^{h/\kappa+r_0/4+o(1)}\) pairs,
a power-small fraction because \(h\le3\) and
\(N\ge n^{r_0-\epsilon+o(1)}\).
All other pairs have \(\log_X Q_k(H_{ij})<h-\nu\).
Combining with the lower height from angular thinning
gives
\[
  \sum_{\substack{p\mid H_{ij}\\v_p(H_{ij})<k}}
                   v_p(H_{ij})\log_X p
                          \ \ge\ \nu-\theta-\eta>\nu/2
\]
for almost every ordered pair of distinct targets in the
retained family, with strict room for the inner allowances.

### Primes outside the endpoint supports

The bounded-order numerical contribution can avoid any
prescribed prime pool of polynomial height. More precisely,
one may fix a pool \(J_i\le X^C\) separately for each
anchor \(i\), for fixed \(C\). For each candidate neighbour
take \(M\) to be the full part of \(H_{ij}\) on
\(J_i FGD\) times its \(k\)-full part outside that
support. The first factor has inner subpower many
possibilities for each anchor by fixed-support
enumeration. The second has the same \(k\)-full
dictionary bound as before. Each such \(M\) with
\(\log_X M\ge h-\nu\) puts the neighbour in the
single class of its anchor modulo \(DM\). That
class has at most \(n^{r_0/4+o(1)}\) points.
Thus the same summation, now per anchor and then
over all anchors, shows bounded-order numerical
weight at least \(\nu-\theta-\eta\) outside
\(J_iFGD\) for almost all pairs.

A prime outside \(FGD\) which divides \(H_{ij}\) and
a complete endpoint letter must divide that same
letter at both endpoints. For \(U,V\) this follows
by reducing \(\Delta_{ij}\) modulo the prime, since
each cofactor pair is primitive. For \(T\), reduce
\(u_iT_j-u_jT_i=\tau G\Delta_{ij}\), using that
\(u_i\) and \(u_j\) are units if the prime divides
the corresponding target. Therefore all such primes
divide \(g_ug_vg_w\). Their normalized radical
mass is at most \(\eta\) on the good pairs from
the support peeling. At orders below \(k\) their
total numerical contribution is at most \(kL\eta\).
With the stated order of the constants,
\(\nu-\theta-(1+kL)\eta>\nu/2\) for sufficiently
small \(L\). Hence the same fixed numerical weight
can be obtained on primes outside the prescribed
pool and **all six complete endpoint letters**.

### Several secants with disjoint prime supports

For every fixed \(m\) one can choose \(m\) distinct
vertices such that all their secants have pairwise
coprime factors of order below \(k\), each of fixed
positive numerical weight and all avoiding the complete
letter supports of those \(m\) vertices.

For this choose \(\eta\) also smaller than
\(\nu/(16m(1+\kappa))\), and allow all polynomial-height
pools of degree bound larger than the combined product
height of \(m\) endpoints and their secants. Omit the
vanishing fraction of anchors for which the angular or
shared-support exceptional pairs do not have vanishing
relative degree, by Markov's inequality. The failure
count for the prime-pool conclusion above is already
uniformly power-small per anchor, since it was proved
by counting only the classes of that anchor.

Select vertices in order. For the next vertex use as
pool the product of prior raw secants and all prior
complete letters. Each old vertex has all but a
vanishing fraction of neighbours allowing a bounded-order
packet of weight greater than \(\nu/2\) on the new
secant. Since the number of old vertices is fixed, a
vertex works for all of them. All these new packets
are coprime to the earlier secants by the pool.
They are also pairwise coprime: a prime common to
the two raw secants from the new vector to two
old vectors divides the old-old determinant by
primitive projective congruence, contrary to the pool.
When choosing a new vertex, also use the peeling
bound to keep its original residual radical overlap
with the earlier packets at most \(n^\eta\).
This excludes only a vanishing fraction of vertices.
At the end strip from each packet the primes in
later endpoints. Each later endpoint costs at most
\(kL\eta\) numerical weight because all packet
orders are below \(k\), and all packets avoided
the fixed factors from the start. At least
\(\nu/4\) weight remains on every packet.

The constants \(\nu,\kappa\) depend on a fixed positive
lower bound for the fibre exponent and the original
affordable margin, and are chosen before the final
thinning parameters and the height limits. The result
forces a positive **numerical** contribution from
bounded scaled orders; the radical lower bounds alone
do not do this, since a small number of primes with
low absolute valuations can contribute normalized
radical mass but vanishing numerical weight. This
conclusion still permits large-radical, low-order
determinants, and gives no unconditional upper bound
on the fibre exponent.

