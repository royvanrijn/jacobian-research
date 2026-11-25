# A strict half-coefficient saving for two highly full summands

For fixed \(\alpha,\beta>0\) there is \(\delta=\delta(\alpha,\beta)>0\)
with the following property. For all sufficiently large fixed
integers \(k\), every fixed nonzero signed shift \(t\), and
every \(\varepsilon>0\), uniformly for \(|t|\le 2X\),
\[
 \#\{(u,v)\in[1,X]^2\cap\mathbb Z^2:
   s_Uu+s_Vv=t,\ u\text{ is }\lceil\alpha k\rceil\text{-full},
   \ v\text{ is }\lceil\beta k\rceil\text{-full}\}
 \ll_{k,\alpha,\beta,\varepsilon}
 X^{\{1/[2\max(\alpha,\beta)]-\delta\}/k+\varepsilon}.
 \tag{SHF1}
\]
Here \(s_U,s_V\in\{-1,1\}\) are prescribed. The count includes
all common divisors of \(u,v\). For positive sums take \(t=N\)
and \(X=N\). In particular (SHF1) improves the coefficient
\(1/2\) of \(1/k\) for sums of two \(k\)-full numbers at every
prescribed target.

The argument uses the limiting fixed-shift bound (FS5).
At its proposed half-coefficient endpoint, the more highly
full summand exhausts its root-height budget. All its
numerical height can then be carried by one combined root
and a cofactor of vanishing numerical height. Mutual modulus
locks force boundedly many actual high roots on both terms,
which can be released simultaneously.
[Fixed-shift half bound (FS5)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#arbitrary-shifts-with-high-free-degrees)

## A compressed term in an arbitrary fixed-shift fibre

Use the signed fixed-shift setup of (FS5) with singleton
dictionaries. In particular the two complete terms are
coprime, their fixed coefficients have polynomial height,
and the free roots have combined raw exponent
\(O(\lambda)\). At each fixed budget the finite degree
lists and truncation are fixed before the height limit.
For each original root position write \(\alpha_i\) for
its raw exponent in a dyadic box and \(W_i\) for its
dictionary exponent, at most \(\alpha_i\). Fix a positive
scaled cutoff \(a\) and write \(B\) for the summed cost of
the original positions of degree less than \(a/\lambda\)
over both terms.

Suppose that, after its fixed coefficient is removed, one
term is a product of boundedly many powered roots of
height \(X^{O(\lambda)}\) and a wide cofactor of height
\(X^{o(1)}\). Assume this description is obtained by
grouping original positions into products with
integer representative exponents and moving their signed
power corrections into the cofactor. The original
root tuples may repeat a complete value divisor-many times.
All limits of numerical exponents here and below are in the
height-first small-budget convention. Then
\[
                               D\le B+o(\lambda),
 \tag{SHF2}
\]
where \(D\) is the fibre count exponent.

Indeed, along a contrary limiting sequence write
\(d=\lim D/\lambda>\lim B/\lambda\), passing to a
subsequence to fix these limits. Name the compressed
term \(V\). A block of original positions on \(U\) whose
degrees are \(\gg1/\lambda\) and whose combined numerical
weight has a positive limit supplies a full modulus. Fix
its root tuple. The compressed \(V\) has a prescribed
unit residue, because the shift is fixed and the complete
terms are coprime. Apply the amplified wide-cofactor
bound (GW) to count \(V\) at cost \(o(\lambda)\) in the
exponent. The original position colouring recovers all
group product coordinates from the exact powered product
for amplified tuples, with divisor-bounded ordered
preimages, as in the window argument of (GD1).
[Wide-cofactor congruence (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Grouping powers and exact amplification fibres](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#diagonal-moments-for-a-joint-dictionary)

The fixed block therefore has normalized dictionary cost
at least \(d\). This condition on all positive-weight high
blocks forces the numerical weight of \(U\) into boundedly
many original positions. To see this, set aside the
positions of normalized dictionary cost at least \(d/4\).
For each fixed \(\gamma>0\), pack the remaining positions
of degree at least \(\gamma/\lambda\) into greedy blocks
of normalized dictionary cost at most \(d/2\). There are
boundedly many blocks, because every closed block except
the last costs more than \(d/4\). If their total numerical
weight had a positive limit, at least one block would
have positive weight along a subsequence, contradicting
the cost floor. The remaining degrees below
\(\gamma/\lambda\) have combined numerical weight
\(O(\gamma)\), by the combined raw root cap.
Letting \(\gamma\) decrease to zero proves the concentration.
The number of set-aside positions is bounded by the raw
cap divided by \(d\lambda/4\), independently of \(\gamma\).
After taking a subsequence discard the set-aside
positions with vanishing limiting numerical weight.

The term \(U\) now supplies its own description by boundedly
many short roots and a cofactor of vanishing numerical
height. Apply (GW) again, this time under a lock on \(V\)
made from its **original** high positions. The same cost
floor and greedy argument concentrate \(V\) into boundedly
many original positive-weight positions. All retained
positions have degree \(\gg1/\lambda\), since each has raw
root exponent \(O(\lambda)\) and positive numerical weight.
The variable remainder has vanishing numerical weight,
so the original remainder positions of degree at least
\(a/\lambda\) have total raw exponent \(o(\lambda)\).
Freeze the other remainder positions at their dictionary
cost, at most \(B\). This leaves an equation of two powered
monomials and the nonzero constant \(-t\), with boundedly
many actual short roots. They are pairwise coprime and
each original root has exponent \(\gg1/\lambda\) in its
own term and zero in the other terms. A vanishing proper
subsum would leave one monomial zero, which positivity
and \(t\ne0\) exclude. Thus (EB1) counts the remaining
roots at \(X^{o(\lambda)}\) cost, uniformly in the fixed
coefficients. If there is at most one retained variable
the count is bounded directly. This contradicts
\(d>\lim B/\lambda\).
[Exclusive-boundary release (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md#exclusive-boundary-release-for-nondegenerate-power-collisions)

## Why equality of the raw and full height budgets is impossible

First count primitive pairs with arbitrary fixed coefficients,
the respective fullness prescribed on their variable integer
parts, and these parts at most \(X\). Complete values and
coefficients may have any fixed polynomial height in \(X\).
Assume \(\beta\ge\alpha\); the roles can be exchanged.
If (SHF1) failed in this setup, the half bound (FS5) would
give a sequence \(k\to\infty\) and actual dyadic root boxes
with
\[
          \lambda=1/k,\qquad
          D/\lambda \longrightarrow 1/(2\beta).
\]
Use the anatomic factorization on the variable integer parts
as in (FS4), choosing the initial coefficient loss
\(o(\lambda)\) and the compression index at least
\(\max(\lceil\alpha k\rceil,\lceil\beta k\rceil)\).
Thus both sides have only high free degrees. At every
fixed \(k\) take the large-height limit after fixing the
coefficients and dyadic lists. The exceptional coefficients
are charged at \(o(\lambda)\), so they have no effect on
this proposed endpoint. Lemma 2.2 of
Bernert–Browning–Lichtman–Teräväinen preserves exactly the
indices up to its prescribed cutoff and moves higher ones
to the selected compression index and an exceptional factor.
[Anatomic factorization and its construction](https://arxiv.org/html/2410.12234v2#S2)

Let \(L_V\) be the total raw root exponent on the side of
fullness \(\lceil\beta k\rceil\), and set \(q=\lceil\beta k\rceil\).
The numerical height of its variable part is at most one
in base \(X\). All its free degrees \(e_i\) satisfy
\(e_i\ge q\), so \(qL_V\le1\) in the height limit.
Apply (FS5) with raw singleton dictionary costs and \(B=0\).
It implies \(\liminf L_V/\lambda\ge1/\beta\), whence
\(qL_V\to1\). Consequently the excess numerical height is
\[
       0\le\sum_{i\text{ on }V}(e_i-q)\alpha_i
                       \le 1-q L_V \longrightarrow 0 .
\]
Write the free product on \(V\) exactly as
\[
             \left(\prod_{i\text{ on }V} x_i\right)^q
                         \prod_{i\text{ on }V}x_i^{e_i-q}.
\]
The combined root has height \(X^{O(\lambda)}\) and the
correction has height \(X^{o(1)}\). This meets (SHF2) with
\(B=0\) and contradicts the positive limiting exponent.
Hence some \(\delta(\alpha,\beta)>0\) works for all
sufficiently large fixed \(k\), even with arbitrary fixed
positive coefficients on the two primitive terms. The
same proof permits the scaled floors to tend to the
specified \(\alpha,\beta\) along the budgets.

For arbitrary common divisor enumerate \(h=\gcd(u,v)\mid |t|\).
In each primitive quotient fix the parts supported on
\(h\). Each support is contained in the primes of a
polynomial-height integer and fixing both costs only
\(X^\varepsilon\) with arbitrarily small \(\varepsilon\),
uniformly in \(h\), by the fixed-support bound used in (FS2).
All remaining primes retain the original high valuations.
Apply the primitive fixed-coefficient result at shift \(t/h\)
and sum over the divisor-many \(h\). This proves (SHF1).
The bound survives a small defect budget on the two
summands, at the cost of adding that budget to the inner
exponent. Indeed the defect radicals can be enumerated
before factoring out the common divisor, and their supported
parts fixed at additional subpower cost.
[Fixed supports and common-divisor removal](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#unequal-powerfulness-at-an-arbitrary-target)

## Weighted codegrees near a high degree floor

Use the varying-shift graph of (FS2), with target values at
most \(X\), arbitrary positive polynomial-height fixed
coefficients per column, and raw root degree floors at
least \(\kappa/\lambda\) in both target columns.
Suppose \(\kappa>0\) is fixed. The coprime weighted
targets obtained by removing the gcd still have degree
floors at least \(\kappa/\lambda\) in their free factors:
first fix their parts supported on the column determinant,
then absorb those parts and the column weights into the
coefficients. The remaining free factors have numerical
height at most one because they still divide the original
target values. Their difference is a prescribed nonzero
shift of polynomial height. The preceding proof with
\(\alpha=\beta=\kappa\) and arbitrary small fixed \(\lambda\)
therefore bounds every such column-pair codegree by
\[
                 O\left(X^{\{1/(2\kappa)-\delta_\kappa\}\lambda
                                     +\varepsilon}\right),
 \tag{SHF3}
\]
for a fixed positive \(\delta_\kappa\) and sufficiently
small fixed budgets. No fullness of the column determinant
is required. The degree floors can also be
\((\kappa-o(1))/\lambda\), by the same strict endpoint
argument.
[Primitive weighted target difference](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#weighted-differences-in-a-varying-shift-graph)

For example, at row and column dictionary exponents each
at most \(1/6+o(1)\) in base \(X^\lambda\), (SHF3) gives
a strict saving from edge exponent \(1/3\) when all
target degrees are at least \((3-\tau)/\lambda\) with a
sufficiently small fixed \(\tau>0\). Use the elementary
codegree sum in (FS3) and the strict half coefficient
applied first at the limiting floor \(3/\lambda\). This
graph consequence has a degree-floor hypothesis and
leaves diffuse targets away from that floor open.

Browning and Verzobio’s Theorem 1.2 counts a coefficient-uniform
generalized Fermat surface with one powered root on each
term. Their Theorem 4.1 transfers such estimates to a
height count with all three values constrained to be
powerful. The fixed-target statement here allows all
valuation levels above the summand floors and leaves the
shift unrestricted. Its bounded-root reduction uses the
one-residue moments and the wide-cofactor release.
[Browning and Verzobio, Theorems 1.2 and 4.1](https://arxiv.org/html/2608.24512v1)


## Grouped dictionaries with a large coefficient on the short letter

In the three-letter setting of (GD3), designate one complete
letter \(L\) and assume
\[
  L\le X^{1-\nu+o_X(1)},\qquad
  \max(P,Q)\asymp X,\qquad \nu>0,
\]
where \(P,Q\) are the other two letters and \(\nu\) is fixed.
The signed relation then gives
\(P/Q=1+O(X^{-\nu+o_X(1)})\). Allow the fixed coefficient
on \(L\) to have arbitrary numerical height. Continue to
require the fixed coefficients on \(P,Q\) to have numerical
height exponent tending to zero in the outer limit. The
joint raw cap \(O(\lambda)\) concerns the original
variable roots only and does not charge this fixed
coefficient. Then (GD3) continues to hold:
\[
       D\le\tfrac12\left(W-\max_T W_T\right)+o(\lambda),
       \qquad W=\sum_{T\in\{L,P,Q\}}W_T .
 \tag{SHF4}
\]
Use the same permitted groups, truncations, and height-first
bounds as in (GD3), with the extra fixed margin \(\nu\).
The raw cap, rational height of the short cofactor in any
release, and numerical height of the short complete letter
are separate quantities.
[Grouped dictionary count](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

To prove this, use normalized limiting costs \(w_T\), their
sum \(w\), and the limiting count exponent \(d\) in units
of \(\lambda\). The bound (GD1) and its greedy concentration
argument concern the numerical weight of **variable roots**
alone. They apply with the present fixed coefficients since
they use unit congruences between complete opposite terms,
with known coefficients. Call a letter supported when
\(w_T>w-2d\). Every supported letter is a product of
boundedly many narrow roots and a wider cofactor of
vanishing variable numerical exponent, apart from its
fixed coefficient. This includes the case of entirely
vanishing variable height. Both tall letters have positive
variable numerical height, since their fixed coefficients
have vanishing numerical exponents.

First recover the quarter inequality
\(4d\le w+\min_T w_T\). Suppose \(w_C<d\) and some
other letter \(A\) is supported. Name the remaining
letter \(B\). Fix \(C\) completely. Every high
positive-weight lock on \(B\) must cost at least
\(d-w_C\), by releasing \(A\) under that lock.
Greedy concentration therefore gives the same narrow
description on \(B\), including the case of empty narrow
list. If \(C\) is tall, its positive variable height
supplies a full lock at cost at most \(w_C\). Releasing
both \(A,B\) under that modulus by (GW) contradicts
\(w_C<d\). If \(C\) is short, then \(A,B\) are tall.
Their ratio is in the stated fixed-gap strip, so (RUW1)
releases both at vanishing normalized cost, again a
contradiction. Thus \(w_C<d\) implies \(w_A,w_B\le w-2d\).
If the quarter inequality failed, a cheapest \(C\) would
have \(w_C<d\) and \(w-w_C>2(w-2d)\), a contradiction.

Now suppose (SHF4) fails and choose a most expensive
letter \(S\). It is supported, and the quarter inequality
gives \(d\le w/3\). Three supported letters would give
\(w>3(w-2d)\), impossible. With exactly two supported
letters, the third costs less than \(4d-w\le d\), contrary
to the just-proved implication. Hence only \(S\) is
supported. Write \(U,V\) for the two other letters with
\(w_V\le w_U\). Then \(w_U+w_V<2d\) and \(w_V<d\).

If \(U,V\) both have positive variable numerical heights,
use the two locks in the proof of (GD3). Fix \(V\) entirely
and release \(S\) under any full positive-weight lock
on \(U\). This gives the cost floor \(d-w_V\) for such
locks. Greedy packing gives one group on \(U\) of cost
\(\beta_U\ge d-w_V-\xi\) with positive-weight high
projection. Choose \(\xi>0\) small enough that
\(w_U-\beta_U<d\). Fix that group's complement on \(U\)
and use any positive full lock on \(V\) to release
the selected group and \(S\). The resulting positive
cost floor \(d-(w_U-\beta_U)\) gives a group on \(V\)
with positive-weight high projection, of cost \(\beta_V\).
Taking the two switches yields
\[
              d\le\beta_U+w_V-\beta_V,\qquad
              d\le\beta_V+w_U-\beta_U .
\]
This contradicts \(w_U+w_V<2d\).

If one of \(U,V\) has vanishing variable numerical
height, it is the short letter \(L\). If the other
unsupported letter costs less than \(d\), fix it and
release \(S\) by (RUW1), using the fixed relative
gap to the short letter. Otherwise that unsupported
tall letter is \(U\) since \(w_V<d\). The same first
greedy step selects a group on \(U\) with complement
cost less than \(d\). Fix just this complement and
release the whole group together with \(S\), using
(RUW1) on the ratio of the tall letters. This again
contradicts the limit exponent. Throughout, the fixed
coefficients change only the prescribed residue or
real-strip centre. Primitivity and the original
position colouring give the same exact-value
multiplicities as in (GD3).
[Unequal-height real-strip release (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

For a reduced progression with fixed step \(M\), two
values \(X<u<v\le2X\) in one unit class have
\(v-u=M j\). Put \(g=\gcd(u,v)\), so \(g\mid j\).
If \(j\le X^\tau\) and \(M X^\tau\le X^{1-\nu}\),
the normalized triple \((u/g,Mj/g,v/g)\) has the required
relative gap: the short-to-tall ratio is at most
\(X^{-\nu}\). Apply the anatomic factorization to
\(j/g\), then move \(M\) to the fixed coefficient on
the short letter. The raw cap there is charged only
the radical of \(j/g\). On the two tall quotients
apply the original anatomic factorization. If
\(\operatorname{rad}(u),\operatorname{rad}(v)\le X^\rho\)
and \(2\rho+\tau\) is sufficiently small, (SHF4) gives
the two-cheapest-letter estimate on these primitive
pairs without an upper bound on \(\operatorname{rad}(M)\).
The comparison of logs at the height \(X/g\) preserves
the fixed positive gap and the small raw-budget condition,
since \(g\le X^\tau\).


## A prescribed coefficient of arbitrary height with one low root

There is another extension of (GD3), giving an alternative
to the relative shortness assumption in (SHF4). Write the
three terms as \(P,Q,L\). Assume that \(P,Q\) have
positive limiting variable numerical height. Impose on
their free positions a fixed scaled lower degree bound
\(a>0\), so each such degree is at least \(a/\lambda\).
Both may have arbitrary fixed coefficients. On \(L\)
allow an arbitrary fixed coefficient and any positions
of degree at least \(a/\lambda\), together with at most
one additional position of arbitrary positive degree.
Use within-letter groups of bounded size, and put the
additional arbitrary-degree position in a singleton.
The original raw root cap is \(O(\lambda)\) within the
same height-first convention as (GD3). All complete values
and fixed coefficients have polynomial height in \(X\).
Then (SHF4) holds without any relative height gap. If a
preliminary record fixes other low positions at cost \(b\),
the estimate applies after fixing it and \(b\) is added
to the total exponent.

In particular the variable part of \(L\) may be any
integer \(R\le X^{C\lambda}\) for a fixed \(C\), with an
arbitrary dictionary of cost \(W_L\). Treat \(R\) itself
as one degree-one root, whose charged cap is at most
\(\min(W_L,C\lambda)\). When \(R\) is a product of
earlier group tuples, its dictionary cardinality is at
most their product of cardinalities, so the same statement
applies with \(W_L\) the sum of those costs. This permits
an arbitrary collection of originally low-degree positions
when their entire variable product has height at most
\(X^{C\lambda}\). It also permits such a cofactor alongside
high roots on \(L\), when their original prime supports
are disjoint. Primitivity keeps \(R\) coprime to both
other letters in either case.

### Releasing the high roots with a determined low root

We need the following variant of (EB1). Let \(x_1,\ldots,x_m,y\)
be pairwise coprime positive integers, each at most \(H\).
The \(x_i\) occur in two disjoint powered monomials, with all
their exponents at least \(E\). Add a third monomial \(C y^r\)
with arbitrary integer \(r\ge1\). All three coefficients
are fixed nonzero signed integers. The number of solutions
is at most
\[
 O_{m,\mathbf e,r,\varepsilon}
       \left(H^{C_m E^{-1/(2m)}+\varepsilon}\right)
 \qquad (m\ge1),
 \tag{SHF5}
\]
uniform in the coefficients. An empty one of the two
\(x_i\)-products is permitted. If \(m=0\), \(y\) is
determined directly.

Follow the dimension descent in the proof of (EB1), with
ambient dimension \(m+1\). The polynomial is nonzero
because its \(y^r\)-coefficient is nonzero. Every branch
carrying sample points is contained in the nondegenerate
equation: in a three-term relation with nonzero monomials
a proper vanishing subsum is impossible. For a geometrically
integral component of degree \(D\), use the exclusive
affine boundary test only for the high positions \(x_j\).
If some such boundary with \(x_j=0\) and all other
coordinates nonzero occurs, a general linear curve
section meets it, as in (EB1). Brownawell–Masser bounds
the degree of the ratio of the monomials which differ
in the exponent of \(x_j\) by \(O_m(D^2)\), while its
valuation at the boundary point has magnitude at least
\(E\). Thus \(D\gg_m \sqrt E\), independently of \(r\).
For components of this degree use the same finite
projection and high-degree covering, which costs
\(H^{(j+1)D^{-1/j}+\varepsilon}\) in dimension \(j\le m\)
and descends to proper subvarieties at each step.
Components without a geometrically integral rational
model are replaced by the lower-dimensional conjugate
intersections as in (EB1).

For a component below the threshold \(c_m\sqrt E\),
each high coordinate has no exclusive affine boundary.
Replace its defining equations by the bounded-height
interpolation kernel on its rational sample, exactly
as in (EB1). If the replacement drops dimension,
descend. Otherwise the component is unchanged.
For \(x_j\) use the arithmetic Nullstellensatz on the new
equations, \(x_j=0\), and \(z\,y\prod_{i\ne j}x_i=1\).
The system has no geometric zero. It gives a nonzero
integer of height \(H^{O_{m,D}(1)}\) divisible by all primes
occurring in \(x_j\) on the sample, since the other
coordinates are units modulo each such prime.
This is a fixed-support bound for every high coordinate.
Once the high coordinates are fixed, the equation gives
at most one positive \(y\), regardless of the low exponent.
There is therefore no cost for counting prime supports
of the low variable. A zero-dimensional component has
bounded cardinality directly. Summing the high-degree
covering exponents along a dimension descent of length
at most \(m\) proves (SHF5).
[Boundary and covering argument (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md#exclusive-boundary-release-for-nondegenerate-power-collisions)

### The two exceptional concentration branches

The proof of (SHF4) now applies with the following
replacements of its two real-strip terminals. In the
quarter inequality, let \(w_C<d\) and assume another
letter \(A\) is supported. If \(C\) has positive
variable numerical height, it provides the modulus
directly. Otherwise \(C=L\), by the positive variable
heights on \(P,Q\). Fixing \(C\) for the preliminary
concentration step shows that both \(P,Q\) concentrate
into boundedly many of their original groups. All their
degrees are at least \(a/\lambda\), so the complementary
variable raw caps are \(o(\lambda)\) and can be frozen.
On \(C=L\) freeze all high positions at cost \(o(\lambda)\),
since their variable numerical height vanishes. Its
remaining free part is one arbitrary powered root.
The resulting system is (SHF5), whose height is
\(H=X^{O(\lambda)}\) and whose high exponents are
\(\gg1/\lambda\), so its count costs \(o(\lambda)\).
If there is no arbitrary-degree position, insert the
dummy root \(y=1\) in (SHF5).

In the one-supported-letter branch of (GD3), only
\(L\) can be an unsupported letter of vanishing variable
height. If the other unsupported letter has cost
less than \(d\), freeze it completely. The supported
letter and \(L\) are counted at cost \(o(\lambda)\)
by the same release, after fixing the small complements.
If the other unsupported letter has cost at least
\(d\), use the first greedy lock step from (GD3) on it.
This selects one original bounded-size group of
high positions whose complement costs less than \(d\).
Fix that complement. Release the selected group and
the supported letter under (SHF5), with the other
free root originally on \(L\). This again contradicts
\(d\). All cases of positive variable height use the
same (GW) switches as before. This proves the stated
extension of (SHF4).

### Progression occupancies with highly full inputs

Fix \(\kappa>0\) and let \(\rho,\tau\ge0\) with
\(0<\lambda=2\rho+\tau\) sufficiently small. In
\((X,2X]\) consider distinct, \(k\)-full integers of
radical at most \(X^\rho\) in a block of a unit residue
class modulo an arbitrary integer \(M\). Let the
block have parameter width at most \(X^\tau\), and
let \(k\ge \lceil\kappa/\lambda\rceil\). For every
fixed \(\eta,\varepsilon>0\) and sufficiently small fixed
\(\lambda\), the number of such integers is
\[
 O_{\kappa,\rho,\tau,\eta,\varepsilon}
 \left(
 X^{\,\min\{\tau,\,\rho/2,\,(\rho+\tau)/4\}
                      +\eta\lambda+\varepsilon}
 \right).
 \tag{SHF6}
\]
There is no requirement on the radical or physical span
of the step. One can take \(k=\lceil\kappa/\lambda\rceil\)
when proving this count, since larger \(k\) only reduces
the allowed set.

For a pair of block values \(u<v\), put \(g=(u,v)\),
\(j=(v-u)/M\), and \(R=j/g\). Since the block is a unit
class, \(g\mid j\). The primitive terms \(P=u/g\),
\(Q=v/g\) have order \(X/g\) and the third term is \(M R\).
The gcd \(g\) is \(k\)-full and hence
\(\operatorname{rad}(g)\le g^{1/k}\). Any deficient prime
of \(P\) or \(Q\), meaning one with its positive valuation
below \(k\), divides \(g\). Their combined radical is at
most \(X^{\tau/k}\), so enumerate their support and exact
parts at \(X^{\tau/k+\varepsilon}\) cost by the uniform
fixed-support bound. Each prescribed deficient part
has height at most \(X^\tau\), since its valuations are
below \(k\). The original tall quotients therefore still
have positive numerical height in their high variable
parts, using \(g\le X^\tau\) and \(\tau\le\lambda\).

Factor the remaining \(k\)-full parts of \(P,Q\) with
the anatomic truncation used in the proof of (SHF1).
Their free degrees are at least \(k\), and the
exceptional coefficients are fixed at arbitrarily small
additional exponent cost. Put \(M\) in the fixed
coefficient on the third term and leave \(R\) as its
one standalone root of degree one. This is exactly the
new variant of (SHF4). Write \(s=\log_X g\) and
\(r=\max\{\log_X\operatorname{rad}(P),
\log_X\operatorname{rad}(Q)\}\), taking narrow fixed
bins for these exponents. The raw dictionary caps are
at most \(r,r,\tau-s\), up to the truncation and bin
losses. Thus the count of primitive triples is at most
\[
  X^{\,\frac12\min\{2r,r+\tau-s\}+\tau/k
                +o(\lambda)+\varepsilon}.
\]
The induction applies at the normalized shell height
\(X/g\), or at base \(X\) with positive comparable
numerical heights. The raw caps of \(P,Q\) follow from
their radicals, including for their \(k\)-full parts, by
the anatomic radical comparison with small truncation
loss.

Reconstructing \(g\) adds at most
\(\min(s,\rho-r)\) to the exponent. Indeed its size is
at most \(X^s\) in the bin, and outside the support of
the quotient with larger radical it has radical at
most \(X^{\rho-r}\). The supported part has only
subpower-many possibilities. Now
\[
  \tfrac12\min\{2r,r+\tau-s\}+\min(s,\rho-r)
              \le \min\{\rho,(\rho+\tau)/2\}.
\]
The first inequality in this minimum follows from
\(\min(s,\rho-r)\le\rho-r\); the second follows from
\(2\min(s,\rho-r)\le s+\rho-r\). The term \(\tau/k\)
is \(O_\kappa(\lambda^2)\). Taking the square root of
the pair bound and also applying the trivial index
width bound gives (SHF6).
[Scaling from primitive pairs](ABC_POWERED_PROGRESSION_OCCUPANCY.md#counting-nearby-pairs-across-all-blocks)

The same argument permits original deficient radical
caps \(X^\xi\) on each endpoint. Primes deficient after
dividing by \(g\) are then either already deficient on
an original endpoint or divide the \(k\)-full part of
\(g\) supported away from both original deficiency
supports. The combined radical of all such primes is
at most \(X^{2\xi+\tau/k}\). Thus add \(\xi\) to the
occupancy exponent, as long as the low parts fixed at
this expense still leave positive variable numerical
height on both tall quotients. A sufficient uniform
margin is \(2k\xi+2\tau\le1-\zeta\) for a fixed
\(\zeta>0\), since each combined fixed low part has
height at most \(X^{2k\xi+\tau}\). As before, the
additional coefficient truncation loss can be taken
arbitrarily small.
