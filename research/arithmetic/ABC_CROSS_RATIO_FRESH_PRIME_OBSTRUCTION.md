# Fresh prime costs of universal cross-ratio transfers

A rational identity on the directions of several primitive triples that
introduces no polynomial prime source beyond their original letters
and pairwise secants must be a cross ratio. For typical bounded
tuples from a large affordable fibre, every secant in such a cross
ratio contributes its own fresh radical and numerical weight. This
places a lower bound on the budget needed to apply the unrestricted
small-budget count to the transferred triple.

## Classification of the universal identities

For \(m\ge1\) use the ring
\[
 \mathcal R_m=\mathbb Q[z_1,\ldots,z_m,1/P],\qquad
 P=\prod_i z_i(1-z_i)\prod_{i<j}(z_i-z_j).
\]
A universal identity here is an element \(f\) such that both \(f\)
and \(1-f\) are units of \(\mathcal R_m\). Equivalently, each is a
nonzero rational constant times a Laurent monomial in the displayed
linear factors of \(P\). The localization describes exactly the
support restriction on the *polynomial divisors* of the identity,
before specialization at any arithmetic tuple.

If \(f\) is nonconstant, it is a cross ratio of four distinct
marks among \(0,1,\infty,z_1,\ldots,z_m\). In fact, on an ordered
configuration of \(m+3\) distinct points of \(\mathbb P^1(\mathbb C)\),
apply the unique Möbius transformation sending its first three
marks to \(0,1,\infty\). The remaining transformed marks give
points of the open set \(P\ne0\). Their substitution in \(f\)
defines a nonconstant holomorphic map from the spherical ordered
configuration space to \(\mathbb C\setminus\{0,1\}\). The
classification of Chen and Salter gives a four-mark cross ratio.
Restricting to configurations with the first three marks already
at \(0,1,\infty\) gives the assertion.
[Chen and Salter, Theorem 4.10](https://link.springer.com/article/10.1007/s00209-025-03893-6)

When exactly one moving mark occurs, the ratio is an anharmonic
transform of \(z_i\), meaning a Möbius transform permuting
\(0,1,\infty\). With exactly two moving marks there are three
possibilities up to those transforms and permutation of the
moving marks. They are \(z_i/z_j\),
\((1-z_i)/(1-z_j)\), and
\((z_i/(1-z_i))/(z_j/(1-z_j))\), according to the retained two
boundary marks. Substitution of \(z_i=\sigma U_i/T_i\) identifies
these with the three orientations of the normalized two-row
crossing identity. In particular, losing all moving-to-moving
secant factors at the polynomial level leaves only one-row
transforms.
[Three orientations of a quotient determinant](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#three-orientations-of-the-same-quotient-determinant)

## Exclusive prime packets on bounded tuples

Use the fixed affordable divisor setup and notation of (PS1)–(PS6).
In particular \(U_i=F u_i,\ V_i=G v_i,\ T_i=D w_i=
\sigma F u_i+\tau G v_i\), where the complete letters are positive
and pairwise coprime, with distinct targets. Write
\(d_{ij}=|u_iv_j-u_jv_i|\) and \(n=X^L\).
Fix \(m\ge2\), a positive lower bound \(r\) on the fibre exponent,
and a positive affordable margin. Work in late balanced endpoint
boxes, with the height limit at fixed \(L\) taken first. For any
fixed \(0<\epsilon<\min(1,r/10)\) there is a thinning of the
fibre with loss at most \(\epsilon\) in the normalized exponent,
and a common projective modulus \(M\) and angular scale \(h\)
with \(h\ge\nu>0\) as in (PS3). There is a cutoff
\(k=\lceil\kappa_0/L\rceil\), with \(\kappa_0\) independent of
\(\epsilon\), such that the following holds for all but a
power-smaller fraction of the ordered distinct \(m\)-tuples of
the retained fibre at each fixed \(L\).

For each edge \(ij\) of the tuple take the full prime powers
\(p^{v_p(d_{ij})}\) with \(0<v_p(d_{ij})<k\), where \(p\) divides
neither \(M\), nor any complete letter of the tuple, nor any
other raw secant of the tuple. Call their product \(Q_{ij}\).
These packets satisfy
\[
 \log_X Q_{ij}\ge (1-\epsilon)h,\qquad
 \log_n\operatorname{rad}(Q_{ij})\ge r-\epsilon .             \tag{CR1}
\]
The same packet has normalized radical mass at least
\(r/3-\epsilon\) at orders not divisible by \(q\) for every
\(q\ge3\) simultaneously, and at least \(r/4-\epsilon\) at
odd orders. All bounds are uniform in the fixed record and
the original family. The power saving in the exceptional
fraction can depend on \(L\).

For the proof run (PS2) with support bound large enough for the
fixed tuple size and the number of anchors used in (PS6), and tolerance \(\eta_L\) at most \(c_m\epsilon L\),
where \(c_m\) will be small. Choose the exponent thinning
allowance \(\gamma<\epsilon/10\) sufficiently small.
Also make \(\eta_L\) smaller than the fixed tolerances in (PS4)
and in the stated bounds below. At this fixed \(L\), the
number of steps in (PS2) can be \(O(1/\eta_L)\).
Its exceptional fraction exponent can shrink with \(L\).
The total normalized exponent loss is still at most \(\gamma\)
by the choice of the threshold in that proof, and the angular
scale still satisfies the uniform numerical gap of (PS3).

With \(\kappa_0\) chosen large in terms of \(r\), use (PS5)–(PS6)
to discard a power-smaller fraction of pairs on which orders
\(\ge k\) outside \(M\) contribute more than \(\vartheta h\)
in numerical weight. Here \(\vartheta>0\) can be arbitrarily
small and the number of anchors depends on it. The discarded
high-order part on a good pair has normalized radical mass at
most \(\vartheta h/(kL)\). For each anchor use the fixed pool of
primes of \(M U_i V_i T_i\). The fixed-anchor divisor bound
gives outside radical mass at least \(r-\gamma-o(1)\) up to
arbitrarily small fixed allowances, with a power-smaller
exceptional fraction. For precision, for any fixed \(b>0\)
the number of queries with outside radical \(\le n^{r-\gamma-b}\)
is \(O(n^{r-\gamma-b/2})\) for sufficiently large \(X\)
at the fixed budget. The counts are uniform by the
fixed-pool enumeration and dyadic cofactor intervals.
The adaptive (SD3) bounds with the same pools give both
defect thresholds with \(r-\gamma\) in place of \(r\), up to
arbitrarily small fixed allowances. Outside the pool and the
high-order part, the valuations of \(d_{ij}\) and its
quotient determinant agree. The primes also avoid the query's
complete letters because a secant prime outside \(FGD\)
identifies the projective directions of the two endpoints
modulo that prime. Angular separation and (PS2) give numerical
weight at least \((1-\vartheta)h-2\eta_L-o_X(1)\) for this
bounded-order pair packet.
[Projective peeling, packing and angular estimates](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md)
· [Fixed-anchor and adaptive-order counts](ABC_SECANT_VALUATION_DEFECT_GRAPHS.md#sparse-pair-bounds-and-adaptive-orders)

Consider now an edge \(ij\) of a bounded tuple and the support
of a letter at another vertex or of another edge. For the
letter choose one endpoint of \(ij\) as the anchor and condition
on the letter. For the other edge, if it shares one endpoint
with \(ij\), choose that endpoint as the anchor and condition
on the other edge. In the disjoint case condition on both its
endpoints and on one prescribed anchor of \(ij\). In each
case the forbidden support is fixed at polynomial height and
the other endpoint of \(ij\) remains a free query. The
simultaneous bound (PS2) applies also with this choice of
anchor. On primes outside \(M\) it bounds the discarded
numerical weight by \(\eta_L\) and hence the normalized
radical mass by \(\eta_L/L\). Endpoint letters of \(ij\)
were already excluded. Sampling across distinct vertices
only changes each conditional exceptional fraction by a
bounded factor as the retained size tends to infinity.
There are boundedly many pairs of these tests. A union bound
gives all of them for a typical tuple. Choose \(c_m\),
\(\vartheta\), \(\gamma\), and the fixed counting allowances
small enough for the accumulated losses above to fit
\(\epsilon h\) and \(\epsilon\) in the two scales. This proves
(CR1) and its simultaneous defect form.

## Primitive cross ratios pay for every moving secant

For arithmetic specialization, take the primitive homogeneous
representatives \((\sigma U_i,T_i)\) of the distinct directions
\(z_i=\sigma U_i/T_i\) and representatives \((0,1),(1,1),(1,0)\)
of the boundary marks. Write \(B_{ab}\) for the signed binary
determinant of two representatives. For an ordered quadruple
of marks \(a,b,c,d\), the three values in the Plücker identity
are, up to sign,
\[
                B_{ab}B_{cd},\quad B_{ac}B_{bd},
                                      \quad B_{ad}B_{bc}.
\]
One of their absolute values is the sum of the others.
Division by their common gcd gives the primitive triple
encoding the cross ratio. Each of the six unordered edges
occurs in exactly one product. The moving-to-moving
determinant for \(ij\) is \(\sigma\tau FG(u_i v_j-u_j v_i)\).
A moving-to-boundary determinant is, up to sign, one of the
three complete letters, and a boundary-to-boundary determinant
is \(\pm1\).

On a tuple satisfying (CR1), every \(Q_{ij}\) for two moving
marks in the chosen quadruple divides exactly its own Plücker
term and is coprime to the other two. Its valuations therefore
survive in full in the primitive triple, with the same defect
thresholds, and the distinct \(Q_{ij}\) have disjoint radical
supports. In particular, a cross ratio with \(\ell\) moving
marks contributes new radical mass at least
\[
                          \binom{\ell}{2}(r-\epsilon)         \tag{CR2}
\]
in base \(n\), on primes absent from all original letters.
For \(\ell\ge2\) the height of the output is at least
\(X^{(1-\epsilon)h}\). The output height is also polynomial
in \(X\), from the Plücker formula.

For \(\ell=2\), two original boundary letters occur at each
moving mark. These give an additional radical cost even when
one letter has many repetitions. Fix the other rows in the
tuple and put all their complete letters and the common
modulus in a prime pool of polynomial product height. For a
candidate row, the product of the two selected complete
letters has at most \(n^{r-\epsilon+o_X(1)}\) possible
values with outside radical at most \(n^{r-\epsilon}\).
This follows by recording that radical and then counting
values of polynomial height supported on the augmented
pool. Each such product can be split into the two letters
in \(X^{o_X(1)}\) ways, and any two complete letters
determine the row through the actual signed relation.
Thus each of the two moving marks contributes at least
\(r-\epsilon\) normalized radical mass outside the letters
of all other rows on a typical tuple. Every one of these
primes divides exactly one of the two raw products of
original letters in the Plücker relation, by primitivity
and the outside pool. The other raw product is a unit
there, so the prime survives primitive reduction. The
two row contributions are disjoint, and (CR2) already
excludes all original-letter primes. The total radical
lower bound is \(3(r-\epsilon)\) for \(\ell=2\).

The direct use of a single cross ratio with at least two
moving marks cannot give a uniform saving by the existing
\(1/3\) count alone. A cross ratio and any three of its
marks recover the fourth, so fixing \(\ell-1\) input rows
leaves one row to count this way. If the output has a small
radical budget, the lower radical bounds impose a floor
of \(r-\epsilon\) on the exponent obtained by applying
the \(1/3\) coefficient when \(\ell=2\) or \(3\), and of
\(2(r-\epsilon)\) when \(\ell=4\). The conversion of a
small absolute radical cap to a relative small budget
can be made in narrow output height bins because the
height has a fixed positive numerical lower bound.
[Unrestricted limiting count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

The obstruction concerns universal identities and this
single-output budget count. An identity on a proper
arithmetic subfamily can have new polynomial factors with
unexpectedly small specialized radicals. Joint restrictions
on several simultaneous cross ratios can also contain
information beyond their separate radical budgets.

## Universal arithmetic support forces a change of markings

For an ordered rational configuration \(z_1,\ldots,z_m\), use
primitive homogeneous representatives of those points and of
\(0,1,\infty\). Let its prime support be the union of the supports
of all pair determinants. This is the union of the primes of
the original three letters of each primitive direction and
the primes of all its raw moving secants. The definition
is independent of the signs of the representatives.

Consider fixed rational functions \(f_1,\ldots,f_s\), distinct
as functions and generically different from \(0,1,\infty\).
Suppose one fixed finite prime pool \(S\) contains every prime
of the output configuration \(f_1,\ldots,f_s,0,1,\infty\)
absent from the input configuration, for every rational
input at which both configurations consist of distinct marks.
It suffices that the arithmetic condition hold on the rational
points of any nonempty real open set of such inputs. Then
either the \(f_i\) are all constants, or there is a fixed
ordered selection of \(s+3\) input marks and the \(f_i\)
are the images of its last \(s\) marks under the unique
Möbius transformation taking its first three to \(0,1,\infty\).
In particular \(s\le m\) in the nonconstant case.

First consider a rational function \(f=A/C\) with coprime
integer polynomials \(A,C\), and put \(B=C-A\).
Suppose an irreducible nonconstant polynomial \(D\) dividing
one of \(A,B,C\) does not divide \(P\). Then \(D\) is coprime
to the product of \(P\) and the other two polynomials.
Choose a source variable \(t\) on which \(D\) depends.
The resultant in \(t\) of \(D\) and that product is nonzero
over the rational function field in the other variables.
Specialize all other variables to integers preserving the
nonvanishing of this resultant and the degree of \(D\).
The specialized nonconstant \(D\) has a root modulo
infinitely many primes: for any nonconstant integer polynomial
\(D_0\) and integer \(t_0\) with \(b=D_0(t_0)\ne0\),
the polynomial \(D_0(t_0+b t)/b\) is integral, nonconstant
and has constant term \(1\). If only finitely many primes
divided its nonzero values, taking \(t\) divisible by all of
them would give values of absolute value \(1\), a contradiction.
Discard the finitely many primes of \(S\), the specialized
resultant and leading coefficients. At a resulting residue
vector all input marks remain distinct, and precisely one
of the three output polynomials is zero modulo the prime.

Lift this vector to rational points in the stipulated real
open set with denominators prime to the chosen prime.
Indeed one can take an arbitrarily large common denominator
congruent to \(1\) modulo that prime and approximate each
coordinate by numerators in its designated residue class.
After evaluation and clearing these denominators, the prime
divides one of the three values and the other two are units
there. It survives primitive reduction. Thus the assumed
arithmetic containment forces \(f_i\) and \(1-f_i\) to be
units of \(\mathcal R_m\) for every \(i\).

At every prime where the input marks are distinct and
the finitely many coefficient primes are avoided, the
\(f_i\) and \(1-f_i\) now evaluate to units. If the numerator
of a reduced \(f_i-f_j\) had a factor outside \(P\), the same
resultant and lifting argument would yield two output marks
coincident modulo a new prime. Hence each \(f_i-f_j\) is
also a unit in \(\mathcal R_m\). Conversely these unit
conditions ensure the arithmetic containment at every prime
outside a fixed finite pool, because evaluation at good
input reduction makes all the units in question \(p\)-adic
units.

Apply the input normalization to the spherical ordered
configuration space as in the single-output proof. The
output \((0,1,\infty,f_1,\ldots,f_s)\) gives a holomorphic map to
the spherical ordered configuration space on \(s+3\) marks.
Chen and Salter's Theorem 4.4 says that, after reordering
input marks and applying one common holomorphic projective
transformation to the image tuple, this map is constant
or selects a subconfiguration. Normalizing its first
three output marks removes the common transformation.
The constant case stays constant after normalization.
The selected-subconfiguration case is exactly the asserted
change of markings. For \(s=1\) this also proves that the
arithmetic prime-support condition alone forces a cross
ratio in the nonconstant case.
[Chen and Salter, Theorems 4.4 and 4.10](https://link.springer.com/article/10.1007/s00209-025-03893-6)

If the selected input marks are distinct modulo a prime
\(p\), the three chosen normalizing marks determine a matrix
in \(\operatorname{PGL}_2(\mathbb Z_p)\), so the output marks
stay distinct modulo \(p\). The reverse implication can fail
when the selected marks coalesce in a common residue class
before normalization. The exclusive packets of (CR1)
exclude that degeneration at the primes used below.

Suppose the selection contains \(\ell\) moving marks from
a tuple obeying (CR1), with \(s\ge1\). For a prime of
\(Q_{ij}\) on a pair of those moving marks, every other
pair determinant among the selected marks is a unit.
If the pair meets a chosen normalizing mark, a single
output cross ratio includes both endpoints of the pair
and also sees the full power \(p^{v_p(d_{ij})}\). For a pair
of normalizing marks use any output. For a pair with
exactly one normalizing mark use its other endpoint as
the output. If both marks are outputs, the normalizing
marks have good reduction and the normalizing matrix
can be chosen in \(\operatorname{GL}_2(\mathbb Z_p)\).
Its determinant is a unit, as are the factors making
the transformed representatives primitive locally, so
the power occurs on that output secant with the same
valuation. Thus the joint output configuration pays at
least \(\binom{\ell}{2}(r-\epsilon)\) normalized fresh radical
mass, with full numerical weights and defect orders on
the respective output letters or secants. This joint
statement uses the actual simultaneous tuple condition.

## Fresh radical costs of extra polynomial factors

The preceding obstruction uses a classification of rational maps.
A separate arithmetic bound applies to a polynomial factor even
when no support-preserving identity is assumed.

Fix a number of slots \(m\), a uniform degree bound and a uniform
polynomial coefficient-height bound. Let \(\mathcal A\) be \(N\)
primitive integer vectors with pairwise distinct projective directions
and coordinate height at most \(X\), where \(N\ge n^r,\ n=X^L\).
Take \(L,r>0\) fixed before the height limit. Consider a fixed
bounded collection of nonassociated irreducible multihomogeneous
integer polynomials in the \(m\) slots. Choose one factor
\(\Phi\) that depends on at least two slots, or has degree at
least two in some slot. Then for all but a power-smaller
fraction of the ordered distinct \(m\)-tuples in \(\mathcal A\),
\(\Phi\) is nonzero and the product of primes of \(\Phi\) absent
from all other factors has normalized radical weight at least
\(r-\epsilon\), for every fixed \(0<\epsilon<r\). A fixed auxiliary
prime pool of polynomial product height can also be excluded.
The coefficients, collection and auxiliary pool can vary with \(X\)
within their uniform bounds.

For the application include the original letter and raw secant
factors among the forbidden factors, without repetitions. The
conclusion for an irreducible factor distinct from them concerns
genuinely fresh primes. If that factor occurs in one of the three
pairwise coprime homogeneous polynomials representing \(f,1-f\)
over a common denominator, its primes are present in the
arithmetic primitive output triple, after including the three
coefficient contents in the auxiliary pool. Applying the estimates to
each eligible factor in a bounded collection and taking a union
bound gives their disjoint contributions simultaneously.

### Resultant pools and fixed integer values

Choose a slot of positive degree for \(\Phi\). Fix the other
slots, called anchors, and consider \(\Phi\) as a homogeneous
binary form in the query slot. An irreducible factor of positive
query degree remains irreducible over the rational function
field of the anchor directions, as a binary form. This follows
by dehomogenizing the anchor slots and localizing, which inverts
only factors independent of the query slot. Discard the at most one direction in each anchor slot
whose second coordinate vanishes, so that the others have
finite direction coordinates. This costs a power-smaller fraction. Distinct
factors of positive query degree are coprime there for the same
localization reason.

For generic choices of anchors, the specialized \(\Phi\) is a
nonzero squarefree binary form, and its homogeneous resultant
with each other form of positive query degree is nonzero.
The values of the forms independent of the query are also
nonzero generically. The exceptional anchor choices are
\(O(N^{m-2})\) if \(m\ge2\), by counting zeros of each
nonzero polynomial of bounded degree on the Cartesian product
of the \(N\) distinct finite directions. This elementary bound
follows by induction on the number of coordinates, separating
the vanishing of a leading coefficient and then counting roots
in the last variable. When \(m=1\) the anchors are empty and
these nonvanishing conditions hold directly.

Put the supports of these resultants and fixed form values
into each anchor's exception pool, together with the auxiliary
prime pool. All have polynomial height in \(X\). At a prime
outside that pool, a primitive query gives a nonzero projective
point modulo the prime, so the resultant excludes simultaneous
vanishing of \(\Phi\) and another factor there. There are at
most a bounded number of queries per anchor where \(\Phi=0\)
because their projective directions are distinct.

For any fixed \(b>0\), the number of signed nonzero values of
\(\Phi\) of radical at most \(n^b\) outside this anchor pool
is at most \(n^{b+o_X(1)}\), with uniform inner error. For
each possible outside radical, adjoin its primes to the pool
and enumerate the integers of polynomial height supported on
the resulting fixed prime product.
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

### A slot of degree at least two

Use such a query slot when it exists. Uniformly for a squarefree
integer homogeneous binary form of fixed degree at least two,
with coefficients of polynomial height, any prescribed nonzero
integer value of polynomial height has \(X^{o_X(1)}\)
representations by pairs of height at most \(X\). Here and
below the degree and all polynomial-height bounds are fixed.
The coefficient bound also holds for the integer factors of
the binary form. To see it, separate a possible coordinate
factor, dehomogenize, use the elementary polynomial bound
on the complex roots, and note that the leading coefficient
of each integral factor divides the original leading
coefficient.

If all irreducible factors are linear, take two nonproportional
ones and enumerate their two signed divisor values. Each
pair of values fixes the two coordinates. Otherwise choose
an irreducible factor \(G(U,V)\) of degree \(d\ge2\) and
fix its signed divisor value \(b_0\ne0\).
Let \(a\ne0\) be the leading coefficient of \(G(T,1)\) and
let \(\theta\) be one of its roots. Then \(a\theta\) is an
algebraic integer of degree \(d\) with all conjugates of
polynomial height, and
\[
       \operatorname{Norm}_{\mathbb Q(\theta)/\mathbb Q}
                    (aU-a\theta V)=a^{d-1}b_0 .
\]
Each possible principal ideal generated by \(aU-a\theta V\)
has this prescribed absolute norm. The number of ideals of
norm \(B\) is at most \(\tau(B)^d\), because at most \(d\)
prime ideals lie over each rational prime and their possible
exponents for norm dividing \(B\) are bounded by the
corresponding rational prime exponent. This is subpower.

Within a principal ideal, candidate generators differ by
units. The logarithms of their conjugate absolute values
are \(O(\log X)\) in each coordinate. The upper bounds follow
from the input and coefficient heights, and the corresponding
lower bounds follow from the nonzero integral norm. In fixed
degree there is a positive uniform lower bound on the maximum
absolute conjugate logarithm of a unit with nonzero logarithm
vector. Indeed, units with all those logarithms in \([-1,1]\)
have bounded integral minimal-polynomial coefficients, so
there are only finitely many possibilities in that bounded
region across all degree-\(d\) fields. The same argument
bounds the number with zero logarithm vector. Since quotients
of units are units, packing separated logarithm vectors
gives \(O_d((1+\log X)^d)\) candidates per ideal, up to the
bounded multiplicity of each vector. Finally \(aU-a\theta V\)
determines \(U,V\) since \(\theta\) is irrational. This
proves the uniform representation estimate.

Combining this estimate with the value count at
\(b=r-\epsilon\) gives a power-smaller fraction of
small-radical choices, including the previously excluded
polynomial zeros.

### Mixed factors linear in each slot

In the remaining case choose any slot on which \(\Phi\)
depends. The specialized query form is linear. Its coefficient
direction is a nonconstant rational function of the anchor
directions. Otherwise its two coefficients would be
proportional with a constant rational ratio, making the
original irreducible form a multiple of a form in one slot.
This would contradict its dependence on at least two slots.
Consequently every coefficient direction has at most
\(O(N^{m-2})\) anchors by the same Cartesian zero bound
applied to the equation for that direction.

Put \(M_0=N^{m-1}\) and let \(E=\omega M_0N\) count the
eligible small-radical edges between anchors and queries
after the generic-anchor restrictions. For each anchor
the preceding enumeration gives at most
\(B_1=n^{r-\epsilon+o_X(1)}\) possible signed integer values.
There are at least \(\omega M_0/2\) anchors with at least
\(\omega N/2\) eligible queries, since the other anchors
account for at most \(E/2\) edges. Unless
\(\omega\ll N^{-1/2}\) already, this allows a choice of
\(D_0\asymp\omega N\) such anchors in distinct coefficient
directions. Each chosen anchor supplies at most \(B_1\)
distinct level lines in the real query plane. Across chosen
anchors the directions are distinct, so these are distinct
lines in total. All queries are among the same \(N\)
distinct real points. The point-line incidence bound gives
\[
           \omega ND_0
                \ll (N D_0 B_1)^{2/3}+N+D_0 B_1 .
\]
If the \(N\) term controls the left side, then
\(\omega\ll N^{-1/2}\). If the last term controls it,
then \(B_1\gg\omega N\). If the first term controls it,
then \(B_1\gg\omega^{3/2}(ND_0)^{1/2}
\gg\omega^2N\). Therefore
\(\omega\ll N^{-1/2}+(B_1/N)^{1/2}\), which is a
power-smaller fraction. This proves the same \(r-\epsilon\)
threshold for mixed factors linear in each slot.
[Szemerédi and Trotter, Theorem 1](https://trotter.math.gatech.edu/papers/38.pdf)

A linear factor confined to one slot has a fixed direction,
and this incidence argument supplies no such bound for its
values. For a rational transfer with at least three eligible
extra factors, the bound charges at least \(3(r-\epsilon)\)
normalized fresh radical mass before charging any original
primes or secants. This has the same budget obstruction as
(CR2) for the direct application of the limiting \(1/3\)
count, whenever the transferred height remains above a
fixed positive numerical exponent and the transferred
radical is in the small-budget range.

### One-slot linear factors in affordable records

There is an additional uniform value bound when a linear
factor depends on one slot and the given family is in an
affordable record. Take \((\sigma U,T)\) as the primitive
projective coordinates. Suppose \(V=Gv\) is the shorter
free letter and the numerical exponent of \(FG\) exceeds
the height exponent of \(V\) by a fixed positive margin.
A linear factor other than the three original letters has
the form \(c_1 F u+c_2 Gv\) in residual coordinates with
nonzero integer coefficients \(c_1,c_2\). Suppose these
coefficients are fixed. This condition holds for a fixed
rational transfer. For each fixed value of the form, the
integer solutions have \(v\) in one residue class of step
\[
              |c_1|F/\gcd(c_1F,c_2G)\gg_{c_1,c_2}F.
\]
In fact \(\gcd(F,G)=1\) implies that
\(\gcd(c_1F,c_2G)\) divides \(|c_1c_2|\). The allowed
\(v\)-interval has length \(o(F)\) by the numerical
margin. Thus a value specifies at most one possible
\(v\), and then at most one \(u\). The fixed-pool
enumeration gives the same \(r-\epsilon\) threshold
for this factor. In this affordable setting a fixed
rational transfer therefore has that threshold at
every extra irreducible factor, including one-slot
linear factors different from the original letters.

## Factors meet every pair of active slots

Write a nonconstant rational map \(f=A/C\) and \(1-f=B/C\)
with \(A,B,C\) pairwise coprime and of the same multidegree
after homogenization. A slot is active when \(f\)
depends on that direction. Every pair \(i,j\) of active
slots occurs together in at least one irreducible factor
of \(ABC\).

If no factor depends on both variables, the rational
logarithmic derivatives satisfy
\[
                 \partial_i\partial_j\log f
                    =\partial_i\partial_j\log(1-f)=0.
\]
Here \(\partial_j\log g\) means \((\partial_j g)/g\) in
the rational function field. On substituting
\(a_j=(\partial_j f)/f\), differentiation of
\(-f a_j/(1-f)\) in variable \(i\) gives
\[
              \partial_i\partial_j\log(1-f)
                  =-\frac{(\partial_i f)(\partial_j f)}
                                {f(1-f)^2},
\]
using \(\partial_i a_j=0\). Both factors in the numerator
are nonzero in characteristic zero for active slots, a
contradiction.

Thus if \(f\) uses \(\ell\ge2\) active slots and all factors
involving multiple slots involve at most \(s\) slots, it
has at least \(\lceil\binom{\ell}{2}/\binom{s}{2}\rceil\)
distinct mixed factors. The fresh radical bound applies to
each on disjoint supports, including ordinary secants if
they occur. For a chosen secant, it excludes the original
letters and the *other* factors, while retaining that
secant itself. These primes can all be charged in the
primitive output triple. This gives a cost of at least
\(\lceil\binom{\ell}{2}/\binom{s}{2}\rceil(r-\epsilon)\)
in normalized radical mass. In particular any transfer
using three or more slots and only two-slot interactions
pays at least the \(3(r-\epsilon)\) cost before counting
its one-slot factors.

## Two mixed factors and a pure free-letter term

For a fixed rational transfer in the affordable setup
above, assume \(f=A/C\) depends on at least two slots and
use a coprime multihomogeneous representation with
\(B=C-A\). There is a precise remaining factor pattern
for attempting a direct saving by the output radical.
Outside a power-smaller exceptional fraction, the
primitive output carries at least \(3(r-\epsilon)\)
normalized radical mass unless the following conditions
all hold.

- There are exactly two distinct mixed factors across
  \(A,B,C\), and they are in different terms. Here a
  mixed factor is one depending on at least two slots.
- Of the three terms, the one with no mixed factor is,
  up to its scalar coefficient, a product of powers of
  one original free letter per active slot. That term
  has no other factors.
- At least one of the two mixed factors depends on
  every active slot.

To prove the first two restrictions, recall that every
active slot has positive degree in each of \(A,B,C\).
At least one mixed factor exists by the preceding argument.
At least three already pay the asserted radical through
the mixed-factor estimates. If there are one or two mixed
factors all in the same term, there are two wholly unary
terms. In each active slot choose one factor from each of
these two terms. The factors are nonassociated because
\(A,B,C\) are pairwise coprime. Fix all other rows and use
the simultaneous resultant pool for the pair. Recording
the joint outside radical of the two chosen values costs
at most \(n^{r-\epsilon}\) when it is at most that bound,
and recording the two signed polynomial-height values
then costs subpower. Two nonproportional linear values
determine the query vector. If either factor has degree
at least two, the uniform binary-form count gives only
subpower multiplicity for its value alone. Thus the two
chosen factors supply \(r-\epsilon\) normalized radical
weight exclusive to this pair for a typical tuple. Use
two active slots to obtain \(2(r-\epsilon)\) in addition
to the mixed-factor cost. The pools exclude the other
output factors, so all this weight survives primitive
reduction on disjoint supports.

There remain two mixed factors in different terms and
a term made entirely of unary factors. Any nonboundary
factor there contributes a further \(r-\epsilon\), by
the fixed-coefficient one-slot bound or the binary-form
bound. A boundary factor corresponding to the target
also gives this cost by direct value enumeration, since
the targets of the original fibre are distinct. If the
pure unary term uses two different boundary factors in
any one slot, the same joint-value argument gives the
cost. Hence a still-unpriced pure unary term can use
only one original free-letter factor in each active
slot, with multiplicity allowed.

For the last restriction, the two mixed factor supports
together must cover every pair of active slots. If both
supports were proper subsets, a slot omitted by the
first and a different slot omitted by the second would
give an uncovered pair. If both omit the same slot, any
pair containing it is uncovered. Therefore one support
contains all active slots.

This reduction leaves possible transfers involving a
pure product of free letters whose radicals outside the
exception pools can be small, and two mixed factors with
fresh radicals. Such a
transfer requires arithmetic control of the specialized
values in the three terms to supply an upper radical
bound and reach the small-budget counting range.
