# Additive energy of complete occurrence dictionaries

A fixed positive mass at positive scaled valuation orders
prevents a small-radical dictionary from having cubic additive
energy. The arbitrary-coefficient half bound controls
differences with large primitive quotients. If the quotient
height is bounded on the radical scale, the common gcd
contains almost all the positive-scale mass and can be
charged directly. This yields a uniform power bound for the
complete dictionaries at the occurrence endpoint.

## Additive energy with positive-scale radical mass

**OE1.** Fix \(H,A,a>0\), put \(n=X^L\), and let
\(\mathcal D\subset[1,X^H]\) be a dictionary of \(N\) distinct
positive integers with
\[
 \operatorname{rad}(x)\le n^r,\qquad
 \prod_{\substack{p\mid x\\ L v_p(x)\ge a}}p\ge n^\sigma
 \quad (x\in\mathcal D),
 \qquad 0<\sigma\le r\le A .
\]
Write \(E_+(\mathcal D)\) for the number of ordered quadruples
\(x_1+x_2=x_3+x_4\) in \(\mathcal D^4\). For every
\(\eta,\varepsilon>0\), at sufficiently small fixed \(L\),
\[
 E_+(\mathcal D)
 \ll_{H,A,a,\eta,L,\varepsilon}
        N^2 n^{\,\max(r/2,r-\sigma)+\eta}X^\varepsilon .
                                                               \tag{OE1}
\]
The bound is uniform in \(r,\sigma\) and in the particular
dictionary under these caps.

For a prescribed nonzero difference \(t=x-y\), choose
\(g=\gcd(x,y)\mid t\) at divisor cost. The quotients
\(x/g,y/g\) have an exact nonzero fixed shift and are
coprime. Choose a fixed height threshold \(C\) in base
\(n\) sufficient for (MCB2) with the desired counting
allowance. If their maximum is at least \(n^C\), apply
(MCB2) with the shift as the third complete letter of
cost zero, permuting the signed relation into a positive
sum. Anatomic factorization gives singleton costs at
most \(r\) on each quotient after arbitrarily small
allowances. The repairless original support restrictions
can be kept by counting an enlarged set in each box.
This gives a representation bound \(n^{r/2+\eta'}X^\varepsilon\)
for these pairs, for arbitrarily small \(\eta'>0\).
[Arbitrary-coefficient dictionary bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

If the maximum quotient is less than \(n^C\), consider
the full powers in \(x\) on primes outside \(g\) whose
orders are at least \(a/L\). They all divide \(x/g\).
Their combined radical is therefore at most
\(n^{CL/a}\). At least \(\sigma-CL/a\) of the normalized
radical mass of \(x\) thus lies on the support of \(g\).
For fixed \(g\), enumerate the remaining radical at cost
\(n^{r-\sigma+CL/a}\) and then enumerate the powers on
the known support at subpower cost, uniformly at height
\(X^H\). This counts \(x\), and \(t\) determines \(y\).
Sufficiently small fixed \(L\) makes \(CL/a\) fit the
counting allowance. Both types of comparison are made
after choosing \(C\) and before the height limit.
[Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

If \(r_{\mathcal D}(t)\) counts representations of the
difference \(t\), then \(E_+(\mathcal D)=\sum_t r_{\mathcal D}(t)^2\).
The term at zero is \(N^2\) and
\(\sum_t r_{\mathcal D}(t)=N^2\). Apply the uniform
nonzero-difference bound just proved and absorb the
small allowances to obtain (OE1).

## The three actual occurrence dictionaries

**OE2.** In the balanced occurrence endpoint, work on
an actual complete-label matching of exponent
\(d=1/3\) in base \(n\), allowing the stated height-first
tolerances. There is an absolute \(\delta_0>0\) such that
each complete role dictionary \(\mathcal D_T\) satisfies
\[
           E_+(\mathcal D_T)\ll
                  n^{\,3d-\delta_0+\eta}X^\varepsilon .
                                                               \tag{OE2}
\]
At least two of these three dictionaries satisfy the
stronger bound \(n^{(5/2)d+\eta}X^\varepsilon\). In addition,
after colouring the primes by their original roles,
each pair of dictionaries from different roles has
mixed additive energy at most
\(n^{(5/2)d+\eta}X^\varepsilon\). All bounds permit
arbitrarily small fixed \(\eta,\varepsilon>0\) in the
height-first convention.

For (OE2) refine to an actual common anatomic profile.
The positive-scale assertion following (SH1) supplies
a fixed positive amount of radical mass on a compact
positive scaled-degree range on every letter. Use a
slightly smaller mass and a lower positive cutoff to
obtain the hypotheses of (OE1) on the complete
dictionaries of that box, at \(r=d+o(1)\). This gives
\(\delta_0>0\) uniformly. The pairwise nonatomic
inequalities imply that at least two roles have
nonatomic mass at least \(1/6=d/2\). On those roles
take a compact positive continuity range carrying
at least \(d/2\) minus an arbitrarily small prescribed
allowance. The weak measures give the corresponding
actual high-order masses on late boxes. Applying
(OE1) gives the stronger exponent. This step can
also be phrased by contradiction on a sequence of
failing fixed budgets: pass to one balanced limiting
profile and choose the cutoffs there before returning
to late finite boxes. A small raw cap or matching
tolerance is absorbed in \(\eta\).
[Positive scales at the balanced endpoint (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)
· [Pairwise nonatomic lower bounds](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Complete-label matching extraction (IM1)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#extracting-an-exact-induced-matching-under-all-three-caps)

For the mixed bound, the colour refinement retains
the exponent and makes every cross-role pair primitive.
For each fixed sum, the primitive terms stay at positive
numerical height on the retained boxes. Apply (MCB2)
with that sum fixed and with the two original radical
caps. Each representation count is at most
\(n^{d/2+\eta}X^\varepsilon\). Summing the squares over
the sums gives \(n^{(5/2)d+\eta}X^\varepsilon\) after
decreasing the input allowances. The mixed difference
energy is the same by rearranging the equality of
two differences. The primitive-height lower bound
follows either from the selected positive-scale mass
or from the strict short-summand sector.

## A quantitative input for centered auxiliary-prime masks

Let \(q\) be sampled among the primes of size
\(n^\beta\) for fixed \(\beta>d\), and reduce an
actual dictionary of size \(n^{d+o(1)}\) modulo \(q\).
As in (AT4), most primes retain
\(n^{d+o(1)}\) distinct residues. The average modular
additive energy of its residue image is at most the
integer energy plus \(n^{4d-\beta+o(1)}\). Indeed
choose one representative per residue so that every
modular quadruple is among those of the original
dictionary. For a nonzero integer four-term
discrepancy the probability of divisibility by a
sampled prime is \(n^{-\beta+o(1)}\).
The resulting most-prime bounds, with arbitrary
small normalized allowances, follow by Markov's
inequality. For the normalized centered transform
\(g\) of (AT3), inequality (AT5) then gives
\[
        \|g\|_1\ \ge\
          n^{-\frac12(d-\min\{\delta_0,\beta-d\})-o(1)}
\]
at the exponent level. For the two stronger
dictionaries replace \(\delta_0\) by \(d/2\).
At \(\beta>2d\), the latter masks thus have
\(\|g\|_1\ge n^{-d/4-o(1)}\), while their
local mixed-character upper bound in (AT4) is
\(n^{-d/2+o(1)}\). Unlike a one-letter mask
without the balanced arithmetic origin, all
three here have a positive power gap between
these two quantities at this prime scale.

This controls the marginal mask supply. A signed
comparison across auxiliary primes still requires
an arithmetic realization of the paired histories
and a lower bound for a positive statistic that
retains actual three-letter occurrence.
[Centered masks and modular energy (AT3–AT5)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#sparse-residue-masks-give-a-polynomial-local-input)
