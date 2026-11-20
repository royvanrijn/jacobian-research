# Scope comparison for fixed-support eventual binomial coprimality

This is a narrow source audit for the proposed ordinary result
\[
 a\geq1,\quad \sqrt{a+d}-\sqrt a\geq1
 \quad\Longrightarrow\quad
 \exists N_0(a,d)\ \forall N\geq N_0:
       \gcd(A_{N;a,d},A_{N-1;a,d})=1,
\]
where
\(A_{m;a,d}(x)=\sum_k\binom mk(am+1)_{dk}x^{m-k}\).
The cutoff is uniform in the two binomial coefficients and in
characteristic-zero coefficient fields, but is unevaluated and may depend
on the support. This audit does not establish historical priority or
assert that every support in this region was previously untreated.

## Edo--van den Essen: the exact classical reduction

The actual arXiv1304.3956v2 PDF was read at printed2--3,
Definition2.7/Conjecture2.8, and printed10--12, Section3.2.
The strong conjecture concerns every consecutive block of length equal
to the number of nonzero monomials. Theorem3.2 gives the precise
two-monomial/consecutive-polynomial equivalence. Its phrase
“have no common zero in C” is on printed10. Lemma3.3 proves the factorial
coefficient normalization; Conjecture3.4 then states the general
coprimality problem. Proposition3.5 treats two particular mixed-variable
families, not the full ordinary same-variable exponent region above.
The final paragraph of Section3.2 notes the obstruction to directly
using recurrences of order at least three.

Our ordinary family is an exact specialization of their classical
dictionary: choose exponent vectors \((a+d,0)\), \((a,0)\). Their polynomial
satisfies
\[
 P_n(t)=\frac{(an)!}{n!}\,t^nA_{n;a,d}(1/t).
\]
Thus their pair \(n,n+1\) is precisely the actual pair here, after a
nonzero scalar and reversal. The reduction itself is not a new claim.
Their listed mixed-variable cases cannot be specialized by simply
setting a variable to one while preserving the factorial functional.

Primary: [Edo--van den Essen, arXiv1304.3956v2](https://arxiv.org/pdf/1304.3956v2).
The cached17-page PDF is pinned in the adjacent source record.

## Rocks: distinguish support parameters from terminal orders

The actual2015 author preprint, dated6 September2015 and28 pages, was
read at printed13--18, Section3.4. Proposition3.14 settles several small
explicit families. Proposition3.15 concerns \(X^m(1+\lambda X)\), with
exponent gap one, and supplies the sufficient arithmetic condition
\(m\nmid((n-1)!)^n\). For each fixed integer \(m\geq3\), that condition
fails for every sufficiently large \(n\). Hence it does not provide a
fixed-support eventual-terminal theorem. Remark3.19 supplies infinitely
many irreducible orders, not an assertion about every order beyond a
cutoff. In Remark3.18, printed18, the author writes: “we haven’t succeeded
in computing this polynomial in its entirety”. The suggested Hurwitz
stability is presented there as computational evidence.

The published2017 article is correctly identified by its publisher as
*On the incompatibility of diophantine equations arising from the Strong
Factorial Conjecture*, Journal of Algebra478,341--352,
doi:10.1016/j.jalgebra.2017.01.039. Its publisher metadata and abstract
were inspected; the direct full-text request returned403. Therefore this
audit does not claim that the2015 Section3.4 numbering or content survives
unchanged in the12-page published version. No theorem from that version
is inferred from its abstract.

Primary: [Rocks2015 author preprint](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf),
[2017 publisher record](https://www.sciencedirect.com/science/article/pii/S0021869317300649).

As a nearby comparison, the actual Filaseta--Rocks2016 published PDF,
printed1, Theorem1.1, gives lower asymptotic irreducibility density
\(\log2\) for the sequence associated with \(1+\lambda X^d\), for fixed
\(d\). This is a different support class and a density statement; it does
not imply every sufficiently large consecutive pair for positive lower
exponents. Primary: [Filaseta--Rocks2016](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf).

## Subtract the repository's stronger existing cases

The proposed theorem is weaker than all-order results on overlaps. The
current repository already proves ordinary all-terminal coprimality for:

* constant-plus-power supports, and adjacent exponents \(d=1\);
* lower exponents \(a=1\) and \(a=2\), at every positive gap;
* relative gaps \(a\geq2,d\geq3a\);
* gap two at lower exponents3,4,5, with the corresponding retained
  results actually allowing every positive Gamma shape;
* additional arithmetic and specialized support families, and every
  lower exponent in the separate very-large-lower-exponent quadratic
  tail region.

These are existing repository results, not claims about what the selected
external papers prove. Their exact proof homes are pinned in the source
record. This list is a subtraction of major overlaps, not a complete
classification of previously covered individual supports.

The critical example \((a,d)=(4,5)\), with exponent support \(\{4,9\}\),
and the strict example \((3,5)\), with support \(\{3,8\}\), illustrate the
new eventual argument beyond the large all-order families just listed.
The relevant scope remains eventual, leaving finitely many terminal
orders per fixed support. The square-gap region contains infinitely many
supports with both the lower exponent and gap unbounded.

The retained `COALESCING_BINOMIAL_SADDLE_BRIDGE.md` Section4 already proves
the local simple-Airy radial transversality formula and its uniform
\(\pi/18\) angle statement in its stated complex-turning-point regime.
The new global proof should not claim a new local Airy mechanism. Its
additional content is the actual physical-pair ordering, complete cycles,
Stokes detours and real-boundary treatment, and at critical equality the
shrinking-window exact-period log-ratio argument. Classical stationary
phase, Morse topology and Gaussian estimates remain classical mechanisms.

## Conclusion and limits of this comparison

None of the precisely inspected external statements above gives the
displayed fixed-support eventual theorem throughout this square-gap
region. This is an applicability comparison, not an exhaustive novelty
or priority conclusion. It does not establish full Strong Factorial,
emptiness of any bad-terminal set, a numerical cutoff, or a cutoff
uniform in the support. No frozen mathematical package was changed.
