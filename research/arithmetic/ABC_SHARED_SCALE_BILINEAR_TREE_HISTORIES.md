# Shared scales in exact bilinear tree histories

A common multiplicative scale can carry the original radical
supports through (BS1). On a pair of larger histories the common
scales can be replaced by a product relation in each paired block.
Even the normalized splits of the trees can retain all original
supports when two paired blocks are used and their orientation is
chosen appropriately. This gives a signed upper bound for histories
carrying original complete labels. It leaves open the lower bound
that would distinguish their actual three-letter incidence.

[Exact bilinear split (BS1)](ABC_AUXILIARY_BILINEAR_SPLIT.md#an-exact-integer-transformation-and-its-unit-laws)
[Support-product unit supply (MSU1–MSU2)](ABC_MULTIPLICATIVE_SPECTATOR_UNIT_SUPPLY.md#uniform-mixing-at-auxiliary-primes)
[Paired-product tree law (AT1–AT2)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#the-tree-estimate-and-the-law-of-shared-slots)
[Two-sided cancellation (CT1–CT2)](ABC_TWO_SIDED_AUXILIARY_TREE_CANCELLATION.md#hypotheses-and-the-two-sided-projection)

## Shared spectator scales and complete-label reconstruction

Fix \(L,H,h>0\), put \(n=X^L\), and take the height limit
at fixed \(L\). Use a fixed number \(m\) of slots, each with a
dictionary of original rows \((A_i,B_i,C_i)\) of positive integers
bounded by \(X^H\). The original caps, primitivity, and repair
conditions can be retained when present. Assume that \(B_i\)
identifies the row in its slot and
\[
                         r_i=\operatorname{rad}(B_i)\le n^h.
\]
This uniqueness applies to an induced complete-label matching.
Take an odd positive integer \(M\) fixed at the given \(L\). At a
prime \(q\), assume all original values are units, and
\(\gcd(M,q-1)=1\). Put
\[
 \mu_M(B)=H_M(B)/B,\qquad
 H_M(B)=\prod_{p\mid B}p^{\max(M,v_p(B))},\qquad
 E_i=H_M(B_i)/r_i^M.
\]
In each of a fixed set of paired blocks \(I_{1,b},I_{2,b}\),
choose spectator units \(S_i\) with
\(\prod_{I_{1,b}}S_i=\prod_{I_{2,b}}S_i\), independently of the
original rows. Here the blocks on each side partition its slots
and all blocks are nonempty. The spectators can have the uniform
law on this subgroup, or a relative density \(1+O_m(\rho)\)
against it. For the second construction, use independent
\(1+O(\rho)\) unit slots on a connecting bipartite graph for each
paired block and set \(S_i\) to the incident-slot product. The
surjective pushforward in (AT1) gives the claimed density. One
can also sample all but one \(S_i\) per component independently
and solve the product relation for the last.

Choose positive representatives for \(S_i\) when making the
integer transformation. In (BS1) set
\[
 R_i=\mu_M(B_i)S_i^M,\qquad
 B_i'=E_i(r_i S_i)^M,\qquad
 x_i=(B_i'/E_i)^{1/M}\pmod q=r_i S_i\pmod q.
                                                               \tag{SS1}
\]
The root in the field is unique. The \(T_i\) in (BS1) can
be any positive representatives of valid residues, including
choices made after the spectators. Each original additive defect
is multiplied by \(R_i\) exactly. If integer representatives are
at most \(q\le 2n^\beta\) with fixed \(\beta\), all auxiliary
integers have polynomial height at the given \(L,M\). For
\(M\ge K/L\), taking \(M\) odd and \(K>0\) fixed, the possible
\(E_i\) have a catalogue of size \(n^{H/K+o_X(1)}\) per slot.
Indeed \(E_i\le X^H\) and
\(\operatorname{rad}(E_i)^{M+1}\mid B_i\), so the support
radical is at most \(n^{H/K}\). The number of possibilities on one
such support at polynomial height is \(X^{o_X(1)}\) uniformly.
[Uniform enumeration on fixed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Conditional on the original rows, \(x_i\) has the uniform law,
or the same relative density bound, on the coset described by
\[
 \frac{\prod_{I_{1,b}}x_i}{\prod_{I_{2,b}}x_i}
 \equiv P_b\pmod q,\qquad
 P_b=\frac{\prod_{I_{1,b}}r_i}{\prod_{I_{2,b}}r_i}.
                                                               \tag{SS2}
\]
Let \(s=\max_{t,b}|I_{t,b}|\). If \(q>n^{2sh}\), the left
side uniquely determines the exact rational \(P_b\). Indeed a
reduced representative has positive numerator and denominator
at most \(n^{sh}\). Two distinct such rationals differ by a
nonzero cross determinant of absolute value less than \(q\).
The denominators on the construction are units.

Suppose the product \(J\) of primes dividing two or more of the
\(r_i\)'s is at most \(n^\sigma\). Given the \(P_b\)'s, record
\(J\) at a cost at most \(n^\sigma\). Every prime outside \(J\)
in a numerator or denominator belongs to exactly one original
slot in its corresponding block and side. Assign these primes to
their slots, then assign the primes in \(J\) to arbitrary subsets
of slots. There are \(X^{o_X(1)}\) such assignments for fixed
\(m,H,h,L\), since all the radical products have polynomial
height. Each assignment gives all the \(r_i\)'s, each \(r_i\)
gives \(X^{o_X(1)}\) choices of \(B_i\), and the injectivity
assumption gives the original row. In particular, the transformed
\(B_i'\) residues alone have at most
\[
                         n^{\sigma+mH/K+o_X(1)}
                                                               \tag{SS3}
\]
possible original labelled histories under these assumptions.
If the \(E_i\)'s are supplied, the exponent \(mH/K\) is absent.
For a single pair the shared support is \(\gcd(r_1,r_2)\).
For arbitrary fixed \(m\), separate bounds
\(\gcd(r_i,r_j)\le n^{\sigma_0}\) give
\(J\le n^{\binom m2\sigma_0}\).

In the balanced occurrence regime, \(h=d+\tau\) with \(d=1/3\)
and small fixed balance allowance \(\tau\). A matching has an
actual one-role radical dictionary of size \(n^{d-o_X(1)}\) in
the saturated limit by the same support enumeration. Consequently
for any fixed \(\beta>2sh\), one can choose \(\tau\) and the
cardinality allowance sufficiently small for (MSU1–MSU2).
There are good primes \(q\in[n^\beta,2n^\beta]\) on the fixed
progression with \(\gcd(M,q-1)=1\), at each fixed budget. All
three original letters are units there when all three radical
caps are at most \(n^h\). The relative density error is
\(O_m(n^{-\gamma})\) for a fixed \(\gamma>0\).
For primitive rows at height \(X/2<C\le X\) and total radical
at most \(n\), the excluded tuples \(J>n^\sigma\) have a power
loss below total exponent \(m/3\) for every fixed \(\sigma>0\)
and sufficiently small fixed \(L\), with the usual \(X^\varepsilon\)
allowance. This is the prescribed-pool argument in (NR2): select
a pair sharing mass \(>\sigma/\binom m2\), fix the first row
and discount its support on the second by (EB5). The budget
threshold may depend on \(\sigma,m\).
[Shared-prime exclusion in the height-first setting](ABC_PROJECTIVE_SUM_NUMERATOR_RADICALS.md#private-denominator-primes-and-support-records)

## Turning scale cosets into the signed tree law

Use \(r=2^\ell\) slots on each side and the constants and
recursion from (AT1). Use the \(x_i\)'s of (SS1) as the actual
leaf variables. For each node \(v\), let \(c_v\) be the product
of the \(r_i\)'s at its descendant leaves. Conditional on the
original rows, rewriting a diagram on the \(S_i\)'s amounts
to replacing \(C_{j,\pm}\) by \(C_{j,\pm}c_{j,\pm}\).
All \(H\)'s, \(p\)'s and arguments \(y\) are identical under
this replacement. The required child relation holds because
\(c_{j,\pm}\) is the product of the two factors at that child.
The modified coefficients are fixed before sampling spectators.
Thus (CT1–CT2) applies to the raw scale cosets, with additive
\(O_\ell(\rho)\) error by the one-tree norms in (AT2).
The constants and masks may be chosen as any functions of the
original history before sampling, provided the stated unit,
consistency, and mask bounds hold for each history.

For example, in the setting of (CT2) on both sides, with
parameters \(p_t,\nu_t\) for side \(t=1,2\), and any complex
weight \(F\) of the original histories with \(|F|\le1\), one has
\[
 \big|\mathbf E F\,W_1\overline{W_2}\big|
 \ll_\ell
      (p_1^2+\nu_1+q^{-1/4})(p_2^2+\nu_2+q^{-1/4})+\rho .
                                                               \tag{SS4}
\]
Here \(W_1,W_2\) are the diagram values from (AT1).
The parameters in this display are uniform upper bounds over
histories, and the expectation permits any distribution of
original histories chosen before the spectator units. This
allows weights obtained from the exact \(P_b\)'s, since these
ratios are fixed by the original history. The proof uses the
relative density on the raw scale coset.

## Reconstruction from normalized split data

At every valid node put \(z_j=-y_{j,-}/y_j\). Its two children
satisfy \(y_{j,+}-y_{j,-}=y_j\), so \(z_j\) is a unit
different from one. For any prescribed original genuine row
\((A,B,C)\), any fixed nonzero scale \(R\) modulo \(q\) and this
\(z_j\), there is a unique valid \(T\) modulo \(q\) in (BS1)
satisfying
\[
       B'/C'=z_j,\qquad
       (y_{j,+},y_{j,-})=(y_j A'/C',-y_j B'/C').
                                                               \tag{SS5}
\]
This follows by the explicit bijection in (BS2), and a positive
integer representative of \(T\) still preserves the exact
defect identity. At fixed \(R\) the same chosen row can be
used at every node, or a different original genuine row can be
used at each node. For (SS5), choose \(T\) after computing the
tree splits. The assertion is about the normalized split and
does not require the chosen \(T\)'s to be independent.

**SS6.** Take two paired blocks containing the two distinguished
leaves with addresses \(+ -\cdots-\) and \(- -\cdots-\) in
different blocks in each diagram, one distinguished leaf per
block per diagram. Assume each block has the same number of
leaves as its partner across the diagrams. The other paired
blocks may be arbitrary. Take the fixed tree constants
common to all original histories. On the valid part of the law,
the two arrays \((z_j)_j\), the constants and \(q\) admit lists
of at most \(n^{\sigma+o_X(1)}\) original histories with
\(J\le n^\sigma\), provided
\[
 q>n^{2h\max\{|I_{1,1}|+|I_{1,2}|,\ \max_{b\ge3}|I_{1,b}|\}}.
\]
The second maximum is zero if there are only two blocks.
This bound involves no repair tail catalogue.
In particular put each distinguished leaf in a singleton block
and partition the other leaves into pairs. The resulting
\(b=r/2+1\) blocks require just \(q>n^{4h}\) at any fixed
depth. They meet the hypothesis of (CT2) on both sides for
\(\ell\ge3\), and that of (CT1) for \(\ell=2\).

To verify the reconstruction, let \(M_v\) be the product
of the raw leaf variables below a node, and write
\[
  \frac{y_{j,+}}{y_{j,-}}
       =\frac{z_j-1}{z_j}
       =\frac{s_{j,+}H_{j,-}}{s_{j,-}H_{j,+}} .
\]
The split data therefore determine the ratio \(H_{j,-}/H_{j,+}\).
At the root, prescribe \(H_{\mathrm{root},+}=1\). The ratio
gives \(H_{\mathrm{root},-}\), hence \(p_{\mathrm{root}}\)
and both subtree products. At a child the incoming \(X\)'s
and its subtree product \(M_j\) are now known, and the ratio
gives \(M_{j,-}/M_{j,+}\). Because \(M_{j,-}M_{j,+}=M_j\),
there are at most two choices. Each choice specifies the
\(H\)'s and \(p_j\), so this continues recursively, leaving
at most \(2^{r-2}\) normalized leaf assignments per diagram.

The root prescription uses a scale freedom in the actual raw
leaf assignment. Multiply just the two leaves
\(+ -\cdots-\) and \(- -\cdots-\) by \(t\in\mathbb F_q^\times\).
Then all \(H_{j,\pm}\) and \(p_j\) are multiplied by \(t\).
To see this, the \(X\)'s at the root are unchanged. At each
child the first \(X\) is the parent's \(p\), so acquires a
factor \(t\). The second \(X\) is inherited from the selected
side of the parent. When its subtree contains a designated
leaf, that leaf always descends along its minus edges.
Thus each side has exactly one \(t\) factor, supplied by its
current \(X\) or its subtree product. The identity for the
\(H\)'s gives the identity for \(p_j\) at every node. Hence all
\(y\)'s are multiplied by \(t^{-2}\) and the \(z\)'s stay
fixed. Choosing \(t=1/H_{\mathrm{root},+}\) makes the true
leaf assignment occur on the reconstruction list.

If two such root normalizations use factors \(t_1,t_2\), the
first two block ratios in (SS2) change by \(t_1/t_2\), because
each block contains one distinguished leaf per diagram. Their
quotient is thus invariant. Each later block ratio is itself
invariant, since its leaves are unchanged. The resulting list of
vectors
\[
                   \mathcal I=(P_1/P_2,P_3,\ldots,P_b)
\]
has length at most \(2^{2r-4}\). In the first coordinate,
\[
 P_1/P_2=
 \frac{\left(\prod_{I_{1,1}}r_i\right)\left(\prod_{I_{2,2}}r_i\right)}
      {\left(\prod_{I_{2,1}}r_i\right)\left(\prod_{I_{1,2}}r_i\right)} .
\]
Each numerator and denominator has at most the number of factors
in the first two blocks of one diagram. Each later invariant has
at most the corresponding block size in factors on each side.
Thus rational reconstruction applies at the stated scale. The
same prime assignment argument as in (SS3) recovers all original
supports at price \(n^{\sigma+o_X(1)}\). Every slot lies on one
side of exactly one signed product, and \(J\) accounts for any
cancelled support. The repair tails are unused because the
inversion starts directly from the normalized tree arguments.

The valid part has relative probability \(1-O_\ell(q^{-1})\).
For the exact uniform coset the two marginal leaf arrays have
independent uniform units within each diagram. At a fixed node,
condition on outside leaves and the node's subtree product.
The inherited \(X\)'s and the status of all its ancestors are
then fixed, since all ancestor products are fixed. The product
on the plus child is uniform, so the ratio of the two child
products takes each possible value with probability at most
\(2/(q-1)\). Thus a first vanishing \(p_j\) at that node has
probability at most \(2/(q-1)\). Take a union bound and use the
relative density bound for the near-uniform law.

## Several auxiliary primes in one integer history

**SS7.** Use the singleton and pair blocks in (SS6). Take a fixed
number \(k\) of distinct good primes \(q_j\in[n^\beta,2n^\beta]\)
with \(\gcd(M,q_j-1)=1\) and put \(Q=\prod_{j=1}^k q_j\).
If \(k\beta>4h\), the two arrays of normalized split ratios at
these primes, with fixed known tree data at each prime, recover
the original histories with \(J\le n^\sigma\) at the same
\(n^{\sigma+o_X(1)}\) list size. In particular \(\beta\) may
be any fixed scale above the dictionary exponent allowed in
(MSU1–MSU2).

Choose the spectator arrays independently at the separate primes,
conditional on the original history. For each slot choose one
positive integer \(S_i\le Q\) with the prescribed residues by the
Chinese remainder theorem. The transformations (SS1) then give
the desired scale at each prime. At each internal node, a single
positive \(T\le Q\) can be chosen by the Chinese remainder
theorem from the residues prescribed by (SS5) on the valid
locus. The integer \(R\) for that node can be any positive unit
at all test primes, including an \(R\) given by (SS1). Thus all
the reductions are reductions of the same exact integer BS
identities, with polynomial auxiliary height at the given
\(L,M,k,\beta\).

At an individual prime the reconstruction up to gauge still
gives at most \(2^{2r-4}\) candidate residue vectors, even when
that prime is too small to recover the exact rational numbers.
Take their Cartesian product and combine each tuple modulo
\(Q\). The actual invariant vector is on this list of at most
\(2^{k(2r-4)}\) candidates. The original radicals are units
at every test prime. Rational reconstruction is now unique
per vector, since each numerator and denominator is at most
\(n^{2h}\) and \(Q>n^{4h}\). The prime assignment and support
bounds finish the argument. A union bound gives an exceptional
probability \(O_\ell(k n^{-\beta})\) for validity of all the
diagrams.

For \(Z_j=W_{1,q_j}\overline{W_{2,q_j}}\) let \(D_j\) denote
the right hand side bound in (SS4) at \(q_j\), including its
implied constant and the actual relative density error. The
one-prime CT1 bound with its AT2 error may also be used for
\(D_j\). Any weight \(F\) of the original histories with
\(|F|\le1\) satisfies
\[
                         \left|\mathbf E F\prod_{j=1}^k Z_j\right|
                               \le\prod_{j=1}^k D_j .
\]
Indeed the spectator arrays, and hence the diagram values,
are conditionally independent across primes before taking the
integer representatives. All tree coefficients and masks must be
fixed given the original history, before the spectators. An
individual power error \(D_j\ll n^{-\gamma_0}\) uniformly across
the chosen primes consequently gives \(O_{\ell,k}(n^{-k\gamma_0})\).
Here \(k\) is fixed before the height limit. The statement gives
an amplified signed upper bound for one actual integer history.

## The effect of discarding the paired product conditions

Consider the singleton and pair blocks, one prime, and fix
the common tree constants. Let \(E\) be the set of true
invariant vectors \(\mathcal I\) modulo \(q\) from the allowed
original histories. There are at most \(n^{2rh}\) such
vectors, since all the \(2r\) original radical supports are
integers at most \(n^h\). Under the actual law the valid
normalized split arrays admit a candidate in \(E\).
If both raw leaf arrays are instead entirely independent
uniform units, that event has probability at most
\[
                    2^{2r-4}\frac{|E|}{(q-1)^{b-1}} .
\]
In the independent law, all blocks are independent, and each
raw block ratio is a uniform unit. The \(b-1\) invariants are
therefore independent uniform units as well. Conditional on a
valid pair of normalized split arrays, every root-normalized
leaf assignment has exactly \((q-1)^2\) original gauge lifts.
All lifts have equal probability in this law. If any assignment
gives an invariant in \(E\), that conditional probability is
at least \(2^{-(2r-4)}\), by the reconstruction list size.
This proves the bound. With \(b-1=r/2\) and \(q\ge n^\beta\)
for fixed \(\beta>4h\), the probability tends to zero at a
power rate. The event has actual probability \(1-O_\ell(q^{-1})\).
Thus the arithmetic dependence persists in the normalized data
even after the common scales are suppressed.

## Coordinating the three original roles

The marked \(B\)-support can be any one of the original three
roles. For a genuine positive row \((a,b,c)\), use
\((b,a,c)\), \((a,b,c)\), or \((a,-c,-b)\) to mark \(a\),
\(b\), or \(c\), respectively. The signed version of (BS1) is
the same polynomial identity. Use \(|B|\) for the minimal
repair, and insert the known sign in the decoded residue:
\[
        x_i=(\operatorname{sgn}(B_i)B_i'/E_i)^{1/M}
              =\operatorname{rad}(|B_i|)S_i \pmod q .
\]
The auxiliary coordinates for the signed orientation may
be signed integers. The underlying original primitive row
and its three caps stay as prescribed.

An original row may occupy several slots, once for each
marked role. The prescribed repetition among slots survives
decoding, since a complete coordinate in any role identifies
its row in the original matching. In particular one can
include all three radical supports of a single genuine row
in the same statistic. Since the three original supports
are pairwise coprime, these repetitions create no within-row
shared prime in \(J\). Between different original rows \(J\)
is bounded by the shared support of those complete rows, so
the same prescribed-pool exclusion applies to independently
chosen original rows, with the multiplicity of their fixed
marked slots irrelevant to that exclusion. The local signed
upper bound still permits different masks for the marked
roles as in (CT2).

## Scope of the signed comparison

Independently chosen uniform bilinear split ratios at fixed
scales carry zero information about the original genuine rows.
In fact their law is the same for all such rows. With \(N\)
equally weighted histories, \(K\) side records and lists of
size at most \(D\), the average decoding success from such
ratios is at most \(KD/N\). To prove this, fix the ratio data.
Taking the union of the lists over the \(K\) records gives at
most \(KD\) histories and the conditional law of the original
history is uniform. This also proves the assertion for
randomized side records. If the ratio laws are \(1+O(\rho)\)
relative to a common law the bound acquires the factor
\(1+O(\rho)\).

The correlated ratios in (SS6) retain information because they
use the actual recursion on the original-root scale cosets.
Adapting the bilinear split to those ratios is consistent with
exact arithmetic. For each valid normalized split, however,
every genuine input row can realize it through a choice of
\(T\) as in (SS5). The decoration alone supplies no signed
lower bound sensitive to the genuine matching. Such a lower
bound would have to use the recovered original labels together
with the three separate residue masks, or a further arithmetic
relation coupling those labels to the spectators. The upper
bound (SS4) and the decoding theorem (SS6) supply compatible
conditions for that step, while the unrestricted abc exponent
still requires the incidence lower bound.

## A positive symmetric statistic and its Mellin overlap

**SS8.** Use two identical diagrams, with the same
paired blocks, fixed tree constants, and fixed masks.
Take the original history on each side independently
with the same law. For one diagram let \(K(v)\) be
the conditional expectation of its value \(W\) under
independent uniform leaf units with block products
fixed to \(v\). For an original history \(H\), write
\(R(H)\) for the vector of products of the marked
radicals in its blocks. With exact uniform spectator
law, the paired signed mean equals
\[
  \mathbf E_{v}\left|\mathbf E_H K(v R(H))\right|^2
    =\sum_\chi |\widehat K(\chi)|^2
                         |\mathbf E_H\chi(R(H))|^2 .
                                                               \tag{SS8}
\]
The group is the product of one unit group per block.
For the same original histories at several primes use
the product of the probe groups and the product of
their respective \(K\)'s.

In fact the common spectator block vector \(v\) is
uniform and the two arrays are independent given
that vector and the original histories. On a fixed
side, multiplication by the original radicals takes
the leaf law to the uniform law with block products
\(v R(H)\). Averaging the two independent histories
gives the squared mean. Multiplicative Parseval
gives the last expression. At several primes the
spectators are independent conditional on the same
history, so the identical argument uses its joint
radical-product vector.

This is nonnegative, but positivity at the required
power scale is a spectral-overlap question. If the
history block-product distribution is itself
\(1+O(n^{-\gamma})\) relative to the uniform law,
its contribution away from the principal character
is \(O(n^{-2\gamma})\|K\|_2^2\). The principal term
is independent of the original history law. Thus
making these original block products independently
near-uniform can erase the signal sought in the
lower bound, even when the conditional spectator
law and the label recovery both work. A transfer
needs nontrivial spectral mass aligned with the
tree kernel on the genuine three-role data.
