# Dictionary counts for fixed shifts with a full divisor

A large full divisor on a fixed shift makes it possible to charge
only the cheaper of two free root dictionaries. Such a divisor
also survives the primitive normalization of the weighted
difference in a varying-shift graph.

## A fixed-shift count without a radical cap on the shift

Let \(s_U,s_V\in\{-1,1\}\), fix an integer \(t\ne0\), and count
positive pairs
\[
                  s_U U+s_V V=t,\qquad \gcd(U,V)=1.
\]
Write each letter as a fixed coefficient times a product of powered
roots, with the positive-index positions pairwise coprime in each
term. All complete values and fixed coefficients have polynomial
height in \(X\), with a fixed height exponent. Only the free roots
enter the raw caps, of sum \(O(\lambda)\). Use the allowed
within-letter groups in the
[grouped dictionary count](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count),
with summed dictionary exponents \(W_U,W_V\). Thus groups have
bounded size or boundedly many shrinking scaled-degree windows,
and nonsingletons use degrees at least \(a/\lambda\) for fixed
\(a>0\). The initial fixed coefficients can have arbitrary
prime support and positive numerical weight. All finite degree
lists and the truncation are fixed before taking the height limit
at each budget.

We assume at least one of the following, with positive fixed
margins and the position of the short letter prescribed:

- The shift \(|t|\) has a divisor \(M\ge X^\nu\) of fullness at
  least \(\lfloor\kappa/\lambda\rfloor\).
- One of \(U,V,|t|\) is smaller than their maximum by a factor at
  least \(X^\tau\) throughout the sample.

For every fixed \(\eta,\varepsilon>0\), at all sufficiently small
fixed budgets the count is
\[
               O\bigl(X^{\min(W_U,W_V)/2+\eta\lambda+\varepsilon}\bigr).
                                                        \tag{FS1}
\]
The sufficiently small budget threshold depends on the fixed
margins, the dimension and raw cap bounds, and \(\eta\).
Inner constants can also depend on \(\lambda\), the degree lists,
the truncation and \(\varepsilon\). The bound is uniform in the
coefficients and shift within their height caps.

### Proof by concentration of the free roots

Use dyadic refinements and the prime-position colouring of (GD1)
at subpower cost. Primitivity keeps the primes in either fixed
coefficient out of the opposite complete term across all rows.

Fix a block of groups on \(U\) with positive numerical weight
and minimum degree of order \(1/\lambda\). Its full product
gives a fixed unit residue on \(V\). Apply the one-letter moment
bound from the proof of (GD1) directly to this residue. The group
with raw cap \(\alpha_J\) uses the character norm of order \(2h_J\),
where \(h_J=\lceil L_V/\alpha_J\rceil\) and \(L_V\) is the raw
cap on \(V\). The norm is at most \(X^{W_J/2+o(\lambda)/(2h_J)+o_X(1)}\).
Hölder on the normalized character measure gives the single
residue count \(X^{W_V/2+o(\lambda)+o_X(1)}\). The complete \(V\)
determines the pair. Therefore the limiting count exponent \(D\)
satisfies \(D\le b+W_V/2+o(\lambda)\), where \(b\) is the cost
of fixing the block on \(U\).

Suppose (FS1) fails along a sequence with \(D/\lambda\to d>w_V/2\),
writing \(w_V=\lim W_V/\lambda\) for the smaller charged cap.
Every positive-weight full block on \(U\) then costs at least
\(d-w_V/2>0\) in normalized units. Apply the greedy concentration
argument of (GD2). With a fixed positive cost floor \(\theta\),
there are boundedly many groups of normalized cost at least
\(\theta/4\). For any fixed \(0<\gamma<a\), the other groups with
all degrees at least \(\gamma/\lambda\) fit into boundedly
many greedy blocks of cost at most \(\theta/2\). None can have
positive limiting weight. The remaining low-degree singletons
weigh \(O(\gamma)\). Letting \(\gamma\) tend to zero shows that
the boundedly many retained groups hold all the variable numerical
weight up to \(o(1)\). Compress window groups as in (GD1), so \(U\)
has boundedly many narrow roots and a wide cofactor of vanishing
numerical height. If the whole variable part has vanishing
weight, it is already such a cofactor.
[Supported letter concentration](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#supported-letters-and-the-quarter-inequality)

Fix any positive-weight full block on \(V\). The congruence on
the concentrated \(U\) now costs \(o(\lambda)\) by (GW) and its
overlap-compatible amplification. Thus every such block on \(V\)
costs at least \(d>0\) in normalized units, concentrating \(V\)
by the same argument. Each application absorbs the fixed
coefficient into the residue target. Exact values recover the
coloured products in the moment or amplification, with
divisor-many factorizations.
[Full-modulus wide-cofactor count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

Under the full-divisor hypothesis release both concentrated
letters via (GW) on the fixed unit ratio modulo \(M\). An exact
complete ratio determines the two coprime values. Colouring
keeps the numerator and denominator factors separate across
rows, so the amplified values have the same divisor-bound
property. Under the real-gap hypothesis apply (RUW1) to the
tall ratio if \(|t|\) is short. If one free letter is short,
apply (RUW1) to the other free value, which is within relative
\(O(X^{-\tau})\) of \(|t|\). The fixed coefficients enter the
strip centre in either case. First fix the chart parameter
and then let it tend to zero. The cost is \(o(\lambda)\) in all
cases, contradicting \(d>0\).
[Unequal-height real strip count](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

## Arbitrary shifts with high free degrees

There is a third terminal argument for the two-letter concentration
proof. Use groups of fixed bounded size in each letter, with the
same lower scaled-degree condition on each nonsingleton group.
Put \(B\) equal to the summed dictionary cost of the singleton
positions of degree below \(a/\lambda\). Allow an arbitrary
nonzero fixed shift \(t\), with no fullness or short-height
assumption. Then, in the same height-first sense,
\[
                  D\le\max\{\min(W_U,W_V)/2,\ B\}+o(\lambda).
                                                          \tag{FS5}
\]

Indeed if \(D\) exceeds the first expression in the maximum
by a fixed normalized margin, the one-residue argument
and then (GW) concentrate both variable numerical weights
in boundedly many of the original groups. Retain only
groups of positive limiting weight. Every retained
nonsingleton contains only degrees at least \(a/\lambda\),
and for a retained singleton its positive weight and
\(O(\lambda)\) raw height force degree \(\gg1/\lambda\).
The remainder has numerical weight \(o(1)\). All its
positions of degree at least \(a/\lambda\) therefore have
total raw cap \(o(\lambda)\). Fix the remainder at cost
at most \(B+o(\lambda)\), counting its low positions by
their dictionary costs.

There are now boundedly many unfrozen original positions.
The remaining equation has two powered monomials and the
constant \(-t\). Each unfrozen position occurs to degree
\(\gg1/\lambda\) in its own monomial and to degree zero
in the other two. The original positions are pairwise coprime
and all three monomials are nonzero. A proper vanishing subsum
would force a remaining monomial to be zero. The exclusive
boundary release (EB1) therefore counts them at \(X^{o(\lambda)}\)
cost, since each root has height \(X^{O(\lambda)}\) and the
number of roots is bounded before the budget limit. With
at most one unfrozen position the count is elementary.
This gives \(D\le B+o(\lambda)\), proving (FS5). It retains
all original positions in a selected group and never
subtracts the cost of a single marginal from the joint cost.
[Exclusive-boundary release (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md#exclusive-boundary-release-for-nondegenerate-power-collisions)

## Weighted differences in a varying-shift graph

Let the edges have distinct labels \((x,y,k)\) with
\[
                  T=Fxy+\sigma Gk>0,\qquad \sigma\in\{-1,1\},
\]
where \(F,G\) are fixed positive integers. Require the complete
terms \(Fxy,Gk,T\) to be pairwise coprime. The row label is \(x\)
and the column label is \((y,k)\). All values have polynomial
height in \(X\). In column \(i\), assign targets one allowed
grouped root representation with raw cap \(O(\lambda)\) and
summed dictionary exponent \(W_i\). Degrees and groups have
the same uniform margins as in (FS1). The fixed target
coefficient may differ between columns.

Suppose \(G\) includes a divisor of numerical weight at least
\(\nu>0\) and fullness at least \(\lfloor\kappa/\lambda\rfloor\)
for fixed \(\kappa>0\). On a common row for distinct columns
\(i,j\), the weighted difference is
\[
                 y_j T_i-y_i T_j
                    =\sigma G(y_j k_i-y_i k_j)=D_{ij}\ne0 .
\]
The nonvanishing follows because \((y_i,k_i),(y_j,k_j)\) are
distinct primitive positive pairs. The common divisor
\(d=\gcd(y_j T_i,y_i T_j)\) divides \(D_{ij}\) and is coprime
to \(G\). In fact both \(y\)'s and both targets are units
modulo \(G\) by primitivity of the incident edges.

Fix the \(D_{ij}\)-supported part of every original target root
at subpower cost. Values with prime factors in the support
of a fixed polynomial-height integer have only \(X^{o_X(1)}\)
possibilities at polynomial height, uniformly in the integer.
This follows from Rankin's bound with any fixed positive
parameter and a sublogarithmic Euler product. The support
refinement fixes \(d\), since \(d\mid D_{ij}\). Now the
normalized terms \(y_jT_i/d,y_iT_j/d\) are primitive. Their
remaining root tuples are images of the old dictionaries
under stripping, with no increase in the raw caps or costs.
The \(y\)'s, fixed supported parts and division by \(d\)
enter the two fixed coefficients. Their shift \(D_{ij}/d\)
still contains the full divisor of \(G\), independently
of its remaining prime factors. By (FS1) the common-row
count is
\[
                     O\bigl(X^{\min(W_i,W_j)/2+\eta\lambda+\varepsilon}\bigr).
                                                        \tag{FS2}
\]
A recovered \(T_i\) determines the row. The same proof
uses the real alternative in (FS1) in place of the condition
on \(G\) if the weighted difference is uniformly shorter than
the weighted terms by a fixed numerical factor \(X^\tau\).
For example this holds if the \(y\)'s lie in a common dyadic
interval and every shift \(Gk\) is smaller than every incident
target by such a factor. Division by \(d\) cancels in the
relative height comparison. A second replacement for the condition
on \(G\) is that the targets in both columns have no free degree
below \(a/\lambda\), with groups of fixed bounded size.
Apply (FS5) after the same support refinement, with \(B=0\).
If low-degree positions have summed costs \(B_i,B_j\),
the corresponding column-pair exponent is bounded by
\(\max\{\min(W_i,W_j)/2,B_i+B_j\}\), up to the usual allowance.

If \(R\) is the row count and \(E\) the edge count, summing
(FS2) over column pairs gives
\[
                 E\ll X^{\eta\lambda+\varepsilon}
                       \left(R+\sqrt R\sum_i X^{W_i/4}\right).
                                                        \tag{FS3}
\]
Indeed the row degrees contribute at least \(E^2/R-E\)
ordered column pairs. The upper bound uses
\(\min(W_i,W_j)/2\le(W_i+W_j)/4\). The result works on any
selected subgraph. It charges arbitrary grouped dictionaries
on the targets, without requiring exact target powers or a
shared power cofactor.

A shift letter with a cheap positive-weight full block can
supply \(G\) by fixing the block first, at the cost of its
prescribed root record. A positive nonatomic mass supplies
such blocks at arbitrarily small normalized cost once the
numerical margins are chosen. At the symmetric minimum-mass
complete-label profile with row and column exponents \(1/6\)
and target root cap \(1/3\), (FS3) reaches \(1/3\) in base
\(X^\lambda\). A further saving by this method needs
target dictionary compression or stronger simultaneous
control of the column-pair codegrees.
[Cheap positive-weight blocks](ABC_DIFFUSE_FACTOR_GRAPH.md#anchor-sets-for-arbitrary-nonatomic-masses)
· [Exact-power weighted differences](ABC_SHIFT_DICTIONARY_WEIGHTED_CODEGREES.md#a-finite-graph-bound)

## Unequal powerfulness at an arbitrary target

Fix \(\alpha,\beta>0\). For every \(\eta,\varepsilon>0\) and
all sufficiently large fixed integers \(k\), uniformly in \(N\),
the number of ordered representations of \(N\) as the sum
of an \(\lceil\alpha k\rceil\)-full integer and a
\(\lceil\beta k\rceil\)-full integer is
\[
       \ll_{k,\alpha,\beta,\eta,\varepsilon}
             N^{\{1/[2\max(\alpha,\beta)]+\eta\}/k+\varepsilon}.
                                                         \tag{FS4}
\]
First consider primitive sums with fixed multiplicative
coefficients allowed, and the respective fullness on the
two variable integer parts. Take \(X=N,\lambda=1/k\) and
use BBLT's anatomic factorization on these parts. Their radical
caps are at most \(X^{1/(\alpha k)}\) and \(X^{1/(\beta k)}\).
Choose its truncation loss much smaller than \(\eta/k\),
and with packing index at least both fullness thresholds.
All nontrivial positive-index positions then have degree
at least \(\min(\alpha,\beta)k\). The exceptional coefficients
are frozen inside the chosen losses. Thus (FS5) applies
with \(B=0\) and charges half of the smaller cap. This works
uniformly with fixed coefficients of arbitrary numerical height.
[Anatomic decomposition, BBLT Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

For arbitrary original summands, enumerate their common divisor
\(h\mid N\) and divide both by \(h\). Fix the full parts of the
two primitive quotients supported on the primes of \(h\).
The uniform fixed-support bound used for (FS2) makes this
a subpower enumeration. Their remaining variable parts retain
the stipulated unequal fullness, since their primes do not
divide \(h\). Apply the primitive coefficient-weighted estimate
to target \(N/h\), at scale \(X=N\). Summing over all fixed
parts and the divisor-many \(h\) proves (FS4). This gives the
half coefficient of the smaller summand radical cap even
when the fullness orders differ and \(N\) has arbitrary radical.
For equal orders, the two-root argument gives a separate proof.
[Two-root proof at equal orders](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#sums-of-two-highly-full-integers-at-every-target)
