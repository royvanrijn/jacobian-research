# Global counts when one complete letter has a small feature model

A fixed-dimensional exact affine description of one complete letter
forces a uniform saving below the balanced one-third coefficient.
The integer feature vectors can vary with height. Their affine minors
control the count, with arbitrary real coefficients in the description
of the letter. This includes every fixed-rank generalized arithmetic
progression of sufficiently small positive numerical box exponent.

## A strict saving for small affine minors

Fix a positive integer \(R\). There are positive constants
\(\theta_R,\delta_R,\kappa_R,L_R\) with the following property.
Let \(0<L<L_R\) be fixed and put \(n=X^L\). Count a family of
positive primitive triples \(a+b=c\) in \(X/2<c\le X\) with
\[
                      \operatorname{rad}(W)\le n^{1/3+\kappa_R}
                      \qquad(W=a,b,c).
\]
Choose one role \(T\). Suppose each triple is assigned an integer
vector \(z\in\mathbb Z^R\), with the same real affine function
giving \(T\) on all assigned vectors. Write \(r\) for the affine
dimension of those vectors. Require that every determinant of
\(r\) feature differences on \(r\) coordinate columns, using one
common anchor and \(r\) other assigned vectors, have absolute value
at most \(X^{\theta_R}\). For \(r=0\) the determinant is \(1\).
Then the number of triples in the family is
\[
       O_{R,L,\varepsilon}
          \bigl(X^{(1/3-\delta_R)L+\varepsilon}\bigr).
                                                               \tag{CF1}
\]
The bound is uniform in the letter role, the feature vectors and the
affine function. The prescribed feature assignment need not be
injective on triples.

### Normalization and a rank-drop criterion

For an affine feature space of dimension \(r>0\), choose one
witness anchor and an independent set of \(r\) feature differences.
Project onto \(r\) coordinate columns where their determinant
\(\Delta_0\) is nonzero. Use as new coordinates \(\Delta_0\)
times the coordinates of each anchored difference in this basis.
Cramer's rule makes these \(r\) new coordinates integers. Each is
an original determinant covered by the hypothesis. The change is
an affine isomorphism on the feature space, so it preserves all
affine identities. Every minor of \(j\) differences in the new
coordinates, for \(0\le j\le r\), now has absolute value at most
\[
                         j!\,(2X^{\theta_R})^j .                \tag{CF2}
\]
This holds on every subfamily, including those of smaller affine
dimension.

The rank-drop criterion uses only these minors. Let \(W\) be an
integer-valued affine function of the features on a family.
If \(W\) varies, there are two spanning witness rows with nonzero
difference \(D\) of their \(W\)-values. Suppose a subfamily
has a common divisor \(M\) of its \(W\)-values with \((M,D)=1\).
If the subfamily has the same affine span, choose an anchored basis
from it, with a nonzero projected minor \(\Delta\). Integer Cramer's
rule expresses \(\Delta\) times each witness difference vector
from that anchor as an integer combination of the basis
differences. Apply the affine function and subtract the two
witness identities. This gives \(M\mid\Delta D\), and hence
\(M\mid\Delta\). Thus any such \(M\) exceeding the upper bound
for the corresponding minor forces an affine dimension drop.

### Limiting setup and cost of fixed prime supports

Suppose (CF1) fails uniformly. Then one can take a sequence of
increasingly small fixed budgets \(L\), with the cap tolerance
and feature-minor exponent tending to zero, and with normalized
count exponent tending to \(1/3\). More explicitly, take all four
proposed constants at most \(1/j\) and obtain a failing fixed
budget \(0<L_j<1/j\), at arbitrarily large heights for that budget.
Take the inner height exponents and the anatomic refinements at
that budget before passing to the next \(j\). The limiting
one-third bound for the original radical budget
\((1+o(1))L\) gives the upper exponent \(1/3\). The given
failure yields the matching lower exponent. The letter role can
be fixed on a subsequence. Individual radical masses tend to
\(1/3\), since each is at most \(1/3+o(1)\) and total normalized
mass less than one by a fixed amount would contradict the
one-third bound.

Use exact valuation layers in the bounded scaled-degree bands
below, with fixed truncation at each inner limit and with the
normalized truncation error tending to zero. Refine the roots and
complete heights dyadically and pass to weak limits
\(\nu_A,\nu_B,\nu_C\) of the normalized valuation measures.
Write \(d_W\) for the nonatomic mass of \(\nu_W\). A prescribed
prime pool of bounded polynomial height has
\(O(\log X/\log\log X)\) primes. In each of its exact valuation
layers its supported factor has subpower many possible values.
This uniform enumeration also applies before any reboxing of
the complementary factors.
[Limiting one-third count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Valuation measure and height conventions](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)
· [Uniform fixed-support counting](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

A compact bounded band of positive degrees carrying a positive
nonatomic mass allows whole prime-power locks of arbitrarily small
fixed radical cost. For precision, choose a compact continuity
band \([a,b]\subset(0,\infty)\). For any prescribed small cost
remove small continuity neighborhoods of the finitely many
atoms in the band whose masses exceed a sufficiently small
fraction of that cost. The nonatomic mass removed can be made
arbitrarily small. Every retained exact-degree root then has
normalized mass at most that fraction for all sufficiently late
boxes. Otherwise a sequence of large coordinates would give
a larger atom in the remaining compact set. This also bounds
the mass of each prime from these coordinates. If the remaining
mass outside a prescribed pool is bounded below, choose full
prime powers there until their combined normalized radical
mass is between \(p/2\) and \(p\), for any fixed sufficiently
small \(p>0\). Fixing their product costs at most
\(n^{p+o_X(1)}\), since its radical is at most \(n^p\) and its
value has height at most \(X\). Its numerical weight is at least
\(ap/2\). Choices of \(p,a,b\) and the atom neighborhoods are
fixed before the inner height limit.

### A target with positive nonatomic mass

Suppose \(d_T>0\). Choose a compact positive band of its
nonatomic measure, and let \(\xi>0\) be smaller than a fixed
fraction of the mass left in that band after the atom-window
removal. Choose \(p>0\) small compared with \(\xi/(R+1)\)
and \(1/(R+1)\). The atom-window removal described above can be
made for this \(p\) while preserving the chosen lower bound.

At any stage with \(T\) moving, choose a nonzero witness
difference \(D\) as in the rank-drop criterion. Its height is
at most \(X\). Discard triples with normalized radical mass at
least \(\xi\) on the prime support of \(D\). The global
effective radical bound (EB5), with zero prices on that pool,
bounds their exponent by \(1/3-\xi/4\) for sufficiently late
budgets and sufficiently small allowances. This estimate is
uniform in \(D\) and is applied to the original triples. At
most \(R\) prior locks have been paid for, each of cost \(p\),
so this discarded set has a strictly smaller exponent than
the current subfamily. On the other rows the retained band
still has positive mass outside the pool. Fix a lock coprime
to \(D\) at cost \(p\), choosing a largest lock class.
The lower bound \(ap/2\) on its numerical weight exceeds the
exponent of every feature minor for sufficiently late budgets,
by (CF2). The rank-drop criterion reduces the affine dimension.

There are at most \(R\) drops before \(T\) is constant. Its
fixed-letter fibre has exponent at most \(1/6+o(1)\) by
(FL1). The rank descent retained exponent
\(1/3-Rp-o(1)>1/6\). This is a contradiction.
[Effective count for fixed prime pools (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Fixed-letter half-radical count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### A target with zero nonatomic mass

Suppose \(d_T=0\). At a balanced limiting endpoint the other
two letters \(U,V\) have wholly nonatomic measures of mass \(1/3\).
For this pair the pooled degree inequality supplies full-divisor
records \(F\mid U,\ G\mid V\) of combined normalized cost
at most \(1/3-\sigma\) for an absolute \(\sigma>0\), and
combined numerical weight exceeding the shorter complete height
by a fixed positive margin. The prime-size hypothesis for the
selection holds because the two measures are wholly nonatomic:
within any bounded scaled-degree band the largest coordinate
mass tends to zero. The tail in the selection is accounted for
in the pooled degree bound. Choose a fixed record retaining
exponent at least \(\sigma-o(1)\). We may take \(V\) to be
the shorter letter in the box, and write \(V=Gv\).
[Zero nonatomic alternatives (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Pooled degree inequality and full-divisor selection](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-pooled-degree-inequality-and-primitive-shift-fibres)

On this record \(v\) is an affine function of the features.
To show this, choose pivot feature columns for the affine span.
Take a minor of their anchored difference columns augmented by
the column of \(v\)-differences. Multiplying its last column
by \(G\) replaces it modulo \(F\) by the signed column of
\(T\)-differences. The determinant with that column is zero
because \(T\) is affine in the features. Thus the original
integer determinant is divisible by \(F\), since \((F,G)=1\)
by primitivity. Expanding it along the \(v\) column bounds its
absolute value by a constant times \(\max v\) times the bound
(CF2) on the corresponding feature minors. This is smaller
than \(F\) for sufficiently late budgets, by the numerical
surplus of \(FG\) above the height of \(V\). So every such
augmented minor is zero. The additive equation gives the same
affine representation for \(U\). This extends the dependence
argument for small feature coordinates to a bound on the
original maximal affine minors.
[Affine lifting within one affordable record](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#the-exact-dependence-of-the-cofactors)

Choose a common compact continuity band in \((0,\infty)\)
leaving less than \(\sigma/32\) of either wholly nonatomic
measure outside it. Fix \(p<\sigma/(100(R+1))\).
While at least one of \(U,V\) moves, select and strip the
\(D_W\)-supported part on each moving letter \(W\), where
\(D_W\) is its own nonzero witness difference. Fix all the
supported parts at subpower cost by keeping a largest class.
Constant letters can be fixed outright. If both remaining
outside radicals are at most \(n^{\sigma/8}\), the count
of this part is at most \(n^{\sigma/4+o_X(1)}\): enumerate
both remaining integers by radical, and the known full
parts recover \(U,V\), hence the triple. Each current
subfamily has exponent at least \(\sigma-Rp-o(1)\), so this
part is smaller. Pigeonhole a moving letter with remaining
radical mass greater than \(\sigma/8\). The band contains
more than \(\sigma/16\) of that mass for late boxes. Wholly
nonatomicity bounds individual coordinate caps there, so a
full lock on that letter can be fixed at cost \(p\), with
positive numerical weight and coprime to \(D_W\). Apply
the rank-drop criterion. The witness differences are chosen
before the support refinements, and the criterion uses them
whenever the refined family still has the old affine span.

After at most \(R\) drops both opposite letters are constant,
and determine a single triple. The retained exponent is
still positive. This contradicts the assumed sequence and
proves (CF1).

## Progressions, polynomial descriptions and additive energy

A generalized arithmetic progression of fixed rank at most
\(R\), with coefficient lengths \(H_i\ge1\), assigns features
\(0\le z_i<H_i\). Its affine minors have size at most
\(R!\prod_i H_i\). Thus (CF1) applies when the coefficient-box
volume is at most \(X^{\theta_R/2}\), independently of the
steps and possible repetitions in the presentation. A single
ordinary progression may therefore have a fixed positive
numerical exponent for its index diameter. The theorem counts
the whole family of triples with their original three caps
and applies at arbitrary sizes of the complete-label dictionaries.

For a polynomial of fixed degree \(D\ge1\) in a fixed number
of integer parameters, take its nonconstant monomials as
features, and write \(R\) for their number. If each parameter
has absolute value at most \(X^\gamma\), all feature
differences have size at most \(2X^{D\gamma}\).
Taking \(\gamma<\theta_R/(2DR)\) gives (CF1) for these exact
polynomial descriptions, uniformly in their coefficients.

The same conclusion holds for an approximate affine description
\(|T-\ell(z)|\le X^\vartheta\), with \(\ell\) real affine,
if both the original maximal-minor exponent and \(\vartheta\)
are sufficiently small in terms of \(R\). First use the integer
normalization (CF2), with \(w\) denoting the normalized coordinates.
If adjoining the integer coordinate \(T\) preserves the affine
dimension, \(T\) is already an exact affine function of \(w\),
since projection to \(w\) is an affine isomorphism on the span.
Otherwise every nonzero maximal affine minor of the augmented
features \((w,T)\) includes the column of \(T\)-differences
and all \(r\) normalized feature columns. Subtract the
corresponding affine approximation from the \(T\) column.
Its entries are then bounded by \(2X^\vartheta\), and the
determinant is bounded by
\[
                         O_R(X^\vartheta X^{r\theta})
\]
if \(\theta\) was the original maximal-minor exponent.
Apply (CF1) in dimension at most \(R+1\) to the exact
coordinate projection for \(T\). In particular, an interval
of width \(X^{\vartheta}\) for any prescribed complete letter,
or such a neighborhood of an allowable progression or
polynomial description, gives a uniform strict saving for
sufficiently small positive \(\vartheta\).

For each fixed \(K\), a complete-letter dictionary with
doubling at most \(K\) also gives a uniform strict saving in
the balanced sector. Its size is at most
\(n^{1/3+\kappa_R+o_X(1)}\) by the radical cap.
Green--Ruzsa's theorem in the integers embeds it in a proper
progression of bounded rank \(R=R(K)\) and size
\(O_K(n^{1/3+\kappa_R+o_X(1)})\). The coefficient-box volume
then fits the numerical threshold for sufficiently small
fixed budgets. One may sum over \(n^{\delta_R/2}\)
prescribed feature descriptions of the same bounded rank
and still retain a strict saving.
[Green and Ruzsa, Theorem 1.1](https://arxiv.org/pdf/math/0505198)

In a full-count subfamily whose fibres over that letter have
comparable cardinalities, its distinct-letter dictionary
\(\mathcal S\) therefore has
\(E_+(\mathcal S)=o(|\mathcal S|^3)\) in the iterated limit.
Otherwise Balog--Szemerédi--Gowers gives a positive-density
subset of \(\mathcal S\) with bounded difference doubling,
hence bounded sum doubling by Plünnecke's inequality.
The corresponding triples retain full count exponent because
their fibres have comparable cardinalities. This contradicts
the bounded-doubling consequence. To obtain a comparable-fibre
subfamily one can pigeonhole dyadic fibre sizes at height
subpower cost. The resulting qualitative energy bound allows
doubling to grow with zero normalized exponent and does not
supply an unrestricted improvement.
[Reiher and Schoen, Theorem 1.1](https://arxiv.org/abs/2308.10245)
· [Petridis, Theorem 3.1](https://arxiv.org/abs/1101.3507)

For a uniform triple from a balanced limiting endpoint, the
maximum probability of any one interval bin of width \(X^\vartheta\)
is at most \(n^{-\delta+o(1)}\), for fixed positive absolute
\(\vartheta,\delta\) from the rank-one approximate model.
Similarly a cover of any one complete dictionary by allowable
rank-\(R\) models needs at least \(n^{\delta_R-o(1)}\) pieces
at that endpoint. Both statements follow by applying the
pointwise count to each piece.

Raghavan's quantitative extension of Freiman–Ruzsa provides
covers by convex coset progressions with dimension and logarithmic
number of pieces
\(O(\log(3K)\log\log(3K))\)
for doubling \(K\). At fixed budget, a doubling loss of a small
positive power of \(n\) can still let the rank grow with height.
The preceding fixed-rank theorem therefore does not give a
uniform saving for the covers supplied by that result.
[Raghavan, Question 1.2 and Theorem 1.4](https://arxiv.org/html/2512.11217v2#S1)
