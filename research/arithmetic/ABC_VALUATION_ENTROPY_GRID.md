# An entropy criterion for releasing actual valuation grids

A narrow union of valuation grids can be certified from
weighted distributions of the integer prime valuations,
without requiring a small set of occupied degree positions.
The condition uses an independent additive entropy estimate
on the two free valuation laws. Replacing independent draws
by even weakly dependent coupled draws fails.

## Independent entropy gives an affordable common grid

Let \(E_1,E_2\) be the laws of the actual integer valuations
on retained factors of the two free letters, obtained by
drawing primes with probabilities proportional to their
normalized radical log weights on each letter. Let \(t\)
be the sum of the two retained normalized masses.
Use \(J=1/L\) at small budget \(L\). Define
\[
 h_i=H(E_i),\qquad
 D_R=H(E_1'-E_2')-\tfrac12(h_1+h_2),
\]
where \(E_1',E_2'\) are independent with the prescribed laws.
If \(D_R\le D_0\) for an absolute fixed \(D_0\), and
\[
                      L\exp((h_1+h_2)/2)=o(1),             \tag{EV1}
\]
then for each fixed \(\epsilon>0\), a fixed-rank signed
grid of width product \(o(J)\) captures at least
\(t-\epsilon\) of the two original masses.
The little-oh bounds are uniform when the displayed
hypotheses hold uniformly. The same conclusion with a
sufficiently small fixed bound in (EV1) gives an
arbitrarily small prescribed constant in the normalized
width product.

Tao's entropy inverse sumset theorem gives from bounded
\(D_R\) both bounded entropy transport distance between the
laws and bounded entropy doubling for the first law.
It then gives one bounded-rank progression \(P\) in
\(\mathbb Z\) such that both laws have transport distance
\(O_{D_0}(1)\) to the uniform law \(U\) on \(P\).
Here one uses the triangle inequality for transport
distance, or glues the two couplings over their common
marginal and bounds the entropy of the sum of their
shifts by the sum of their entropies.
Also \(|h_1-h_2|\le 2D_R\) by the lower entropy bound on
an independent difference, and
\(|H(U)-h_i|\le O_{D_0}(1)\) from the bounded shift
entropy in either direction. Thus
\[
                 |P|\ll_{D_0}\exp((h_1+h_2)/2).
\]
If \(P\) is not initially proper, its set has doubling
bounded by \(2^{\operatorname{rank}(P)}\), so the integer
Freiman theorem covers it by a bounded-rank proper
progression of size within a fixed factor.
[Tao, Theorem 1.11 and §6](https://arxiv.org/pdf/0906.4387)

Each law can now be coupled as \(U+Z_i\) with
\(H(Z_i)=O_{D_0}(1)\). Since the expectation of
\(-\log\mathbb P(Z_i=z)\) is \(H(Z_i)\), truncating
this self-information gives a set of
\(O_{D_0,\epsilon}(1)\) shift values covering probability
at least \(1-O(\epsilon)\). Both valuation laws are
therefore captured up to the desired mass loss in
finitely many translates of the *same* progression.
Concatenate its one set of steps with the translate
steps, each of width one. This covers the degree
sample with coefficient-width product \(O_{D_0,\epsilon}
(\exp((h_1+h_2)/2))\).
The literal valuation range can then be truncated to
\(e\le B_0/L\) for a sufficiently large fixed \(B_0\).
Indeed the excluded normalized radical mass is at most
\(2/B_0\) by the two numerical height moments; entropy
is used before the truncation and need not be
recomputed afterwards.

The resulting criterion fits the rank compression and
full-modulus release bounds (VG1) and (VG2). If the two
letters are balanced and \(t>1/3\) by a fixed margin,
and if the complementary letter supplies a cheap full
divisor as below, choose the transport and height
truncation losses and the complementary divisor cost
small enough to leave a strict count saving. The
entropic hypotheses can also hold triplet by triplet:
at fixed \(L\), the clipped degree sets lie in the
fixed finite interval \([1,B_0/L]\), so one can refine
by the clipped selected set and use one admissible
grid for each such set. This costs only a fixed
constant at that budget. The small-width assumption
is uniform in the selected set. All limits take
\(X\to\infty\) at fixed \(L\) first.
[Rank compression and all-three-cap counting](ABC_SPARSE_VALUATION_GRID_COMPRESSION.md)

One may instead use the pooled valuation law on the
retained prime mass of both letters. Bounded doubling
entropy \(H(E'+E'')-H(E)=O(1)\) for two independent
copies and \(L\exp H(E)=o(1)\) directly give the same
covering by Theorem 1.11(ii).
Two separate doubling controls without bounded
common \(D_R\) also suffice with the stronger condition
\(L\exp(h_1+h_2)=o(1)\), because separate progressions
then have to be concatenated.

## A diffuse full-support model

Take two high-density progressions of fixed steps \(2g\)
and \(3g\) in an integer interval of size comparable to
\(J\), with offsets coprime to their respective steps,
so each contains \(\asymp J/g\) points for \(J\to\infty\)
at fixed large \(g\). For each law take the mixture of
weight \(1-\delta\) on its uniform progression and of
weight \(\delta\) on the uniform full interval, where
\(\delta=(\log g)^{-2}\).
Both laws have full support of size \(\asymp J\) and
\[
                 h_i=\log(J/g)+O(1+\delta\log g).
\]
The lower bound follows by conditioning on the selector,
and the upper adds at most its binary entropy.
When both draws use their progression component, their
difference occupies \(O(J/g)\) possibilities on a
single residue class modulo \(g\). Otherwise there
are \(O(J)\) possibilities. Conditioning on the two
selectors and adding their entropy gives bounded
\(D_R\), uniformly as \(g\to\infty\).
Taking \(g\) sufficiently large and then fixed meets
the fixed-small variant of (EV1). All weak scaled
degree laws are diffuse. No growing exact divisibility
order captures positive mass: on each progression
there are either no compatible positions or a
fraction \(O(1/P+g/J)\) divisible by an order \(P\);
the full interval tail has fraction \(O(1/P+1/J)\).
This is a model for the entropy criterion, not an
arithmetic existence assertion.

## Coupled difference entropy does not imply a common grid

Take integers \(g,M\to\infty\) with \(M/g\to\infty\)
and let \(J=(g+1)M\). Define two degree progressions
\[
 A=3J+g[-M,M]_{\mathbb Z},\qquad
 B=3J+(g+1)[-M,M]_{\mathbb Z}.
\]
Use the marginal law assigning equal weight to their
two uniform distributions, both for \(E_1\) and \(E_2\).
The marginal entropy is \(\log M+O(1)\), so the
entropy-effective size is \(o(J)\). The scaled laws
have a diffuse limit and exact mean \(3J\), compatible
with balanced numerical height one.

Couple the two draws by first choosing a common
uniform selector between \(A\) and \(B\), and then
drawing independently within that selector. Conditional
independence implies mutual information at most the
entropy of the selector, \(\log 2\). Their coupled
difference has entropy at most \(\log M+O(1)\):
conditional on the selector it has at most \(4M+1\)
possible values. Thus both the coupled entropy deficit
and the joint-information loss are bounded.

Nevertheless any fixed-rank progression of total
coefficient width \(o(J)\) captures at most
\(1/2+o(1)\) of this marginal law. If it met
positive proportions of both \(A\) and \(B\),
the sumset of the two intersections would have
\(\gg J\) elements. Indeed each sum has at most
\(O(1+M/g)\) presentations by the two progressions,
as equality of sums gives the index displacement
\((g+1,-g)k\). There are \(\gg M^2\) pairs across
the intersections. But the doubled box of a fixed-rank
progression of width product \(o(J)\) has at most \(o(J)\)
elements. This is a contradiction.
For two letters of mass \(1/3\), the common-grid
mass provided by these coupled data is therefore
at most \(1/3+o(1)\), which cannot furnish the
strict release margin. Independently sampled
differences instead have \(D_R\to\infty\):
half their conditional cases are cross-grid,
of entropy at least \(\log J-O(1)\) by the same
fiber bound, while the other cases have entropy
at least \(\log M-O(1)\).


## Choosing a full divisor inside exact valuation layers

The modulus assumption can be verified directly from
individual primes, even when the weak scaled-degree
limit on the complementary letter \(T\) is atomic.
Suppose its primes satisfying
\[
           v_p(T)\ge\kappa/L,\qquad p\le X^{\tau L}
\]
have combined normalized radical mass at least \(b>0\).
Here \(\kappa,b\) are fixed. For \(0<\mu<b\), greedily
choose a product of distinct eligible primes with
radical between \(X^{\mu L}\) and \(X^{(\mu+\tau)L}\).
Take their entire powers in \(T\) as the integer
modulus \(m\). Then \(m\mid T\), \(m\ge X^{\kappa\mu}\),
and \(m\) is at least \(\lfloor\kappa/L\rfloor\)-full.
The radical box bound gives at most
\(X^{(\mu+\tau)L+\varepsilon}\) records for \(m\).
The selection may depend on the triple, since the
bound sums over every admissible modulus of this
size and radical. Primitivity makes \(m\) a unit
modulus for the two free letters.

This supplies the full modulus in (VG2) at normalized
cost \(\mu+\tau\). If \(\tau=o(1)\) in the outer budget
limit, \(\mu>0\) can be any sufficiently small fixed
constant. One sufficient condition is positive mass in
a positive compact scaled-degree band with maximum
radical mass at any individual integer degree \(o(L)\),
in unnormalized units. Each prime in that band is
bounded by its entire exact-valuation root product,
so it meets the small-prime cap. More generally the
small primes can be selected from within a root
position of positive mass.
[Radical box count and full congruence release](ABC_SPARSE_VALUATION_GRID_COMPRESSION.md#a-two-letter-release-inequality-with-several-unrelated-grids)

In particular, consider a sequence of boxes in which
all but \(o(L)\) unnormalized radical mass on each of
two free letters lies in a single integer interval of
width \(o(1/L)\) and with positions \(O(1/L)\).
The two intervals can be unrelated. Represent them by
a common unit step of width equal to the larger width
and two separate offsets with width one. The combined
width product is \(o(1/L)\), and the frozen cost is
\(o(1)\) in normalized units. If the complementary
letter has fixed positive eligible mass with
\(\tau=o(1)\), (VG2) gives a count exponent at most
\(\eta L+\varepsilon\) for each fixed \(\eta>0\) and
all sufficiently small fixed budgets in the sequence.
Choose the modulus cost and the chart allowance
before taking the inner height limit. Thus the
normalized limiting exponent of this sector is zero.

This applies to the narrow prime-degree profile used
to test the moving external arithmetic rows. Its
scaled degree laws converge to atoms, while its
microscopic mass is spread over increasingly many
exact degree positions. The external rows allow that
profile at the two-thirds exponent. Actual boxes
satisfying the same narrow bands and microscopic
cap on the third letter admit the full-modulus
release just proved.
[External moving-selector limitation](ABC_EXTERNAL_MOVING_ROWS_LIMIT.md)
