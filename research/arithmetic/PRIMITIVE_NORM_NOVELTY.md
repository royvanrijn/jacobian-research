# Novelty comparison and applications of classical estimates

Assessment date: 18 November 2026. The targeted review below concerns
the exact separately primitive two-block theorem. The later sections
retain the earlier comparison of quantitative corollaries and the exact
discriminant -8 classification. This does not establish priority for
the separate exact -7 or -11 classifications.

The general theorem is a structured deduction from classical estimates.
The proof-level comparison below shows that its full arbitrary-pair
power-saving asymptotic already follows by a short adaptation of
Baker 2006, Sections 4 and 6, together with elementary two-block
Möbius tails and local-series completion. This is stronger historical
overlap than a matching special case. The exact designated-block
statement was not found printed there, but its central analytic
conclusion should not carry a standalone novelty claim.
All eventual -8 existence, including both strong cores, also follows
from Plaksin by the parity transformation below. The exact all-target
classification is a separate candidate contribution.
It must be distinguished from Liouville's nineteenth-century
whole-vector primitive theorem.

## Targeted historical review of the two-block theorem

### The statement being compared

Fix arbitrary primitive positive-definite integral binary forms
\(Q_1,Q_2\), allowing different discriminants, nonprincipal form
classes and nonfundamental order discriminants. The comparison target is
\[
 R^*_{Q_1,Q_2}(N)=
 \#\{(x,y,u,v)\in\mathbb Z^4:
       Q_1(x,y)+Q_2(u,v)=N,\ (x,y)=(u,v)=1\},
\]
with the two gcds taken in their prescribed blocks. The central
quantitative statement is a local-density main term plus
\(O_{Q_1,Q_2}(N^{1-\delta})\), for some fixed \(\delta>0\),
uniformly through the positive targets. The existence consequence
must cover every sufficiently large separately primitive locally
admissible target, including arbitrary permitted target valuations.

The precise \(7/8+\epsilon\) exponent, the printed finite admissibility
modulus and the genus-invariant corollary are assessed separately.
An older theorem can precede the central result without printing
those particular consequences. A general relative asymptotic, even
without a power saving, could already precede the existence theorem.
Conversely, a theorem concerning coefficient content, gcd of the
whole four-vector, primitive lattice embeddings, target averages,
or only one quadratic field does not automatically have the full
comparison scope.

The review checked both explicit statements and short deductions from
published theorems. It did not require the older paper to use the
phrase separately primitive.

### Findings at theorem level

| Inspected primary result | Exact location and permitted scope | Consequence for the comparison |
|---|---|---|
| Brüdern--Fouvry, 1994, Lagrange's four squares theorem with almost prime variables | Theorem 3, printed pp.61--62; proof closes p.74. Uniform absolute remainder sum for four-square solutions with arbitrary overlapping squarefree coordinate divisibilities. Includes divisors 2,3 and divisors of the target. | A short deduction already gives the complete Gaussian two-block coordinate asymptotic with a positive power saving; see below. Its stated quadratic family has square diagonal coefficients. |
| Schinzel, 2013, On Sums of Four Coprime Squares | Final paragraph of printed p.111, rather than only Theorem 2. Both designated gcds are expressly imposed, with an odd prime in each Gaussian block. | Direct published eventual Gaussian existence on its full local range \(N\not\equiv0,1,5\pmod8\). |
| Brüdern, 2014, Sums of Squares Coprime in Pairs | Theorem 2, p.219, and Lemmas 2,4,5: four squares with all six pairwise gcds, a large-common-factor tail, scaled-square estimates and completed singular series. | The coprimality-sieve architecture is classical. Its main theorem alone has a different local range; its lemmas provide the appropriate Gaussian comparison. |
| Plaksin, 1984/1985, two binary forms with prime arguments | Main theorem, printed p.562, equations (17)--(18); the second form has fundamental discriminant. Only the first pair is constrained to positive primes; the free pair is unrestricted. | Genuine short primitive-existence deductions include the complete eventual discriminant -20 range. It does not state the arbitrary two-block coordinate asymptotic. |
| Bredikhin, 1965; Bredikhin--Linnik, 1966; A.I.Vinogradov, 1967 | The actual formal results concern prime plus one quadratic form or ideal norm. The 1966 and 1967 papers advertise a broader well-distributed-sequence mechanism. | Important prior genus/class-mixing methods, including a primitive-ideal specialization of the 1966 good-trajectory lemma. They do not state the required two-block main term and error; see the precise interface below. |
| Podsypanin, 1975 | Theorem 1, printed p.459, for arbitrary positive forms of rank at least four; rank-four error \(O(N^{19/20+\epsilon})\). Theorem 2, p.487, gives local positivity. | The restriction is squarefreeness of each scalar coordinate, not either designated block gcd. Its exact local range can omit infinitely many represented two-block targets. |
| Baker, 2006, The values of a quadratic form at square-free points | Definition (1.2), pp.101--102; Theorem 1, p.104, arbitrary fixed nonsingular forms in at least four variables with nonzero varying target. | The main theorem has different restrictions, but Sections 4 and 6 yield the full arbitrary-pair two-block power-saving asymptotic by the short checked adaptation below. |
| Pall, 1949, Representation by Quadratic Forms | Definition on p.346; Theorem 2 p.349; binary embedding discussion p.353. | Primitive means gcd of maximal minors. For a represented number this is whole-vector primitivity, not two prescribed block gcds. |
| Alaca--Pehlivan--Williams, 2014, sum of two binary forms | Definition on p.1396; Theorem 1.1, pp.1397--1399. | It counts unrestricted four-tuples and gives reductions/evaluations for selected coefficient families. It has no designated gcd condition or arbitrary-pair primitive asymptotic. |
| [Jancevskis, 2009](https://pcwww.liv.ac.uk/~karpenk/JournalUDT/vol04/no2/8Jancev09.pdf), Waring with convergent sieve sequences | Theorem 1, pp.133--134, under \(s>2^k\), with independent scalar sequence restrictions. | For squares this requires at least five variables and does not impose joint block gcds. |
| [Kane--Xue, 2023/2024](https://arxiv.org/pdf/2303.18071v2), primitive even-rank theta series | Theorem 1.1 and Section 3, especially preprint PDF p.8. | The primitive condition is gcd of all coordinates; inversion is over square divisors of the target. |
| [Wessel--zur Verth, 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13736505/), quadratic forms in restricted variables | Condition L1 and Remark 1.1 require an off-block rank at least five, hence at least ten variables; weights factor over scalar coordinates. | Neither the dimension nor the restriction matches the present four-variable problem. |

Retained primary editions, online comparison URLs and inspection scopes
are recorded in [the existing source ledger](primitive_norm_quaternary_source_ledger.json).
The rows do not claim a complete independent verification of the
deep analytic proofs.

### The Gaussian asymptotic is already a short 1994 deduction

Brüdern--Fouvry's Theorem 3 counts positive solutions of
\[
 x_1^2+x_2^2+x_3^2+x_4^2=N,\qquad d_i\mid x_i,
\]
and bounds the sum of absolute remainders over the entire squarefree
divisor box \(d_i\le H=N^\theta\), for \(0<\theta<1/22\).
Here squarefree means each \(d_i\) is squarefree; common primes
between different entries are allowed. There is no inherited
\(N\equiv4\pmod{24}\) restriction in this theorem: that hypothesis
enters the later almost-prime application. The source proof already
finishes before that specialization.

Use the divisor vectors \((d,d,e,e)\) and signed weights
\(\mu(d)\mu(e)\). Their absolute remainder sum is bounded by the
published full-box remainder sum, since this is an injectively indexed
subset of that box. The elementary positive-binary representation
bound gives a large-gcd tail \(O_\epsilon(N^{1+\epsilon}/H)\).
Completing the nonnegative genus coefficients with any moment
\(0<\alpha<1\) gives \(O_\epsilon(N^{1+\epsilon}H^{-\alpha})\);
the local Jacobians and inclusion-exclusion are exactly those
in the canonical two-block proof. The source's power saving,
together with sufficiently small epsilon losses, therefore gives
\[
 R^*_{x^2+y^2,x^2+y^2}(N)
       =\pi^2N\prod_p\beta_p^*(N)+O(N^{1-\delta})
       \quad\text{for some }\delta>0.
\]
The source uses positive coordinates and real factor \(\pi^2/16\).
Multiplying by 16 recovers the nonzero signed coordinates. After
block primitivity, a zero coordinate forces its mate to be \(1\)
or \(-1\), so all zero-coordinate contributions are
\(O_\epsilon(N^\epsilon)\). This proves the displayed full-coordinate
deduction without a new Gaussian analytic estimate.

Thus Gaussian two-block existence and its power-saving coordinate
asymptotic are not candidate novelties of the present work.
The deduction is recorded here; it is not represented as a theorem
printed verbatim by Brüdern--Fouvry.

This does not make the arbitrary-form theorem a variable-change
corollary. The family \(\sum d_i^2x_i^2\) has square determinant
squareclass over \(\mathbb Q\), whereas
\(x^2+y^2+u^2+2v^2\) has nonsquare determinant squareclass.
Even some self-sums differ locally: the self-sum of
\(x^2+xy+y^2\) is anisotropic over \(\mathbb Q_3\), while four
squares is isotropic there. Indeed its completed-square equation is
\(a^2+c^2+3(b^2+d^2)=0\), forcing all four primitive-scaled
coordinates divisible by three; four squares has the nonzero
\(\mathbb Q_3\)-solution \((\sqrt{-2},1,1,0)\).
Fixed congruence restrictions and finite-index sublattices do not
change these rational invariants. Extending the source's analytic
proof to other quadratic spaces is an additional argument.

### What Plaksin does and does not already imply

The primitive local count on Plaksin's p.555 is an auxiliary
quantity. Identity (9) adds common-coordinate prime powers back
to recover the unrestricted count used in his main theorem.
His local condition at p.562 requires both prime-coordinate
residues to be units and does not require the free pair primitive.

A literal second-block Möbius inversion produces
\(\phi(p,q)+d^2\psi(X,Y)=N\). For \(d>1\), the second form has
content \(d^2\) and nonfundamental discriminant. This violates
the stated hypothesis on \(\psi\). Passing to a nonmaximal-order
sublattice also needs free-coordinate progression restrictions
absent from that theorem. These are precise missing interfaces,
not assertions that his proof cannot be extended. Forgetting the
prime restriction gives a lower subset of scale \(N/\log^2N\);
it cannot supply the full order-\(N\) coordinate main term.

There is nevertheless a short classical deduction covering all
eventual primitive \(x^2+5y^2\) self-sums. Set both forms to that
fundamental-discriminant form and regroup as
\[
 p^2+5q^2+x^2+5y^2
       =Q(p,y)+Q(x,q).
\]
The bad gcds require \(p\mid y\) or \(q\mid x\).
For each fixed prime, the offending free coordinate has
\(O(\sqrt N/p+1)\) choices; the remaining fixed binary equation
has \(O_\epsilon(N^\epsilon)\) solutions. Summing gives
\(O_\epsilon(N^{1/2+\epsilon})\) discarded witnesses.

The prime-unit local range is precisely
\(N\bmod8\in\{2,3,4,6,7\}\): at two the prime block is
\(6\bmod8\), and the free form takes \(0,1,4,5,6\bmod8\).
At five, a nonzero square plus a square covers every residue.
At other odd primes fix the prime-coordinate residues to one and
use the nondegenerate free binary form. The selected roots lift.
The good local factors are \(1+O(p^{-2})\) away from target
primes, and exceed one at target primes; the finitely many bad
factors are uniformly positive. Consequently the
\(\gg N/\log^2N\) main term survives the deletion. This proves
every sufficiently large locally admissible discriminant -20
target already by a short Plaksin deduction. It supplies neither
a numerical onset nor the sharper unweighted \(NJ(N)/\log N\)
count.

### Earlier genus-to-class arguments are relevant, but not automatic

The actual Bredikhin 1965 survey, all three pages of the 1966
Doklady announcement, the 1966 Bredikhin--Linnik paper and all
nine pages of Vinogradov 1967 were checked. Their formal
equations are prime plus one norm; the survey also treats divisor
convolutions. Their general comments about replacing the prime
by a well-distributed sequence are method outlines with further
estimates required, not a quantified theorem already covering
two primitive-form sequences.

Bredikhin--Linnik's Lemma 5, printed pp.157--160, is a genuine
prior class-mixing mechanism: exact-exponent-one split primes
are independently reoriented among ideal classes. Restricting
its coprime residual ideals to primitive ideals preserves these
flips, giving an immediate primitive-ideal proof specialization.
But a primitive genus convolution main term and a bound for its
bad-trajectory contribution with the second form still have to
be proved. The paper's Theorem C concerns the prime sequence.
The 1965 norm-inflation argument additionally clears rational
denominators and does not preserve primitive coordinates.
These sources must be credited for class mixing; neither
their formal statements nor the inspected adaptations establish
the full comparison target.

### Squarefree-coordinate theorems have a different local range

This distinction can be checked without trusting a title.
Let \(Q=x^2+16y^2\). If all four coordinates are nonzero
squarefree, an odd first coordinate gives norm \(1\bmod8\).
An even first coordinate is \(2\bmod4\), giving norm
\(4\) or \(20\bmod32\). Two such norms therefore never
sum to a multiple of 32. Nevertheless
\[
 Q(8t,1)+Q(8,1)=64t^2+96\qquad(t\ge1)
\]
has both blocks primitive, all coordinates nonzero, and target
divisible by 32. Thus the published Podsypanin/Baker restrictions
do not even cover the full local target set for this allowed
nonmaximal-order pair. A new sieve adaptation could work; an
application of their stated theorems is insufficient.

Fixed block changes of basis can repair the finitely many
squarefree-coordinate local obstructions in that example. Thus the
example excludes a literal use of the fixed-form theorem, rather
than every indirect deduction. A fixed-cutoff argument alone fails
with the tail \(O_\epsilon(N^{1+\epsilon}/P)\).
The polynomial estimates in Baker's actual proof allow a growing
cutoff and repair this issue, as written out next.

### The general theorem has a short, explicit Baker adaptation

Baker's main Theorem 1 imposes a different restriction. A literal
substitution into its squared-coordinate divisor family is invalid.
The complete replacement is written in
[the explicit first-power proof](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md#explicit-first-power-adaptation-of-baker),
equations (B.1)--(B.12).

The proof reconstructs Lemmas 9--10 for \((d,d,e,e)\), including
common primes, two and all target valuations. It tracks the
first-power Jacobian, gives a full exponent budget with
\(H=P^{1/1024}\), and replaces the squarefree-point Euler product
by the exact primitive-block factors. A fixed radial smooth
weight equals one on the full shell, so the conclusion concerns
raw coordinate counts, with zero blocks separately removed.
The conservative resulting saving is \(1/16384\).

The underlying estimates and strategy are classical. This is a
short proof adaptation, not an immediate application of Baker's
printed Theorem 1 and not a claim that he printed the two-block
statement. Under the review's standard, it substantially weakens
a standalone novelty claim for generality and separate
primitiveness alone. The sharper \(7/8+\epsilon\) bound, explicit
local arithmetic and exact exceptions require their own comparisons.

### Search coverage, access limitations and verdict

The primary genealogy followed the cited theorem numbers and
definitions of primitiveness, including prime/free variables,
divisor overlaps, target valuations, order discriminants and
whether a formal generic-sequence interface was supplied.
English, French, German and Russian proper/primitive and
coprimality terminology was searched. Public OpenAlex
forward-citation queries were executed for Plaksin,
Brüdern 2014, Schinzel, and Estermann 1951/1962. These indices have duplicate
translated records and visibly sparse coverage; they are leads,
not a proof of absence.

No live MathSciNet or zbMATH search result set was obtained.
MathSciNet returned the institutional authentication page;
zbMATH web requests returned 403 and its API requests timed out.
These are access failures, not negative database searches.
The original Jones 1949, Pall/Ross 1946 and James 1938 papers remain
uninspected. A later primary Dietmann--Elsholtz account credits James
with exact three-summand Gaussian split support and a binary
approximation allowing two bad prime factors in each summand.
That confirms his relevance without excluding the uninspected original.

Actual publisher first-page scans of Estermann were obtained. The
1951 opening Theorems 1--2 require at least five individually
squarefree squares (eight for all large targets). The 1962 opening
problem is indefinite, with fixed nonzero target and growing boxes.
Their complete papers, later statements and proofs remain unread.
The precise first-page inspections and retained images are recorded
in the ledger. None of these sources is excluded merely from a title
or an inaccessible search result.

**Targeted verdict.** The full arbitrary-binary-pair power-saving
asymptotic is obtainable by the short Baker adaptation above. Gaussian
existence and its asymptotic have still earlier antecedents, and
Plaksin also gives the full eventual -20 existence result. Thus the
general separately primitive theorem alone is insufficient as the
central novelty claim of a standalone paper. The potential remaining
contribution is the sharper quantitative error, explicit ramification
and conductor criteria, and exact exception classifications.
Genus invariance is an immediate consequence, rather than a separate
analytic breakthrough. Priority for those refinements remains
unresolved; the remaining access limits do not erase the positive
classical deduction already found.

## Assessment of each claim

| Claim | Closest verified knowledge | Defensible assessment |
|---|---|---|
| Primitive norm prime-factor criterion | Classical quadratic ideal factorization; Euclidean argument in the class-one cases | Classical ingredient. |
| Full -8 shell and whole-vector primitive counts | Liouville 1860; Deutsch 2004/2005 and 2008 | Explicitly classical, including the whole-vector dyadic obstruction 16. |
| Möbius truncation, large-factor tail and singular-series completion | Brüdern--Fouvry 1994 Theorem 3; Brüdern 2014 Theorem 2 and Lemmas 2,4,5 | Classical architecture; the full Gaussian two-block power-saving asymptotic is already a short deduction. |
| Arbitrary fixed binary pair, each block primitive, every locally admissible target | Short first-power adaptation of Baker 2006, Sections 4 and 6, plus elementary two-block tails and local completion | The central power-saving theorem is a classical deduction; no standalone novelty claim for this scope alone. |
| All imaginary quadratic orders and invariance within fixed genera | The preceding arbitrary-binary-pair theorem and integral local isometries | Explicit class-group consequence, with no class-number restriction; an immediate corollary, not a new analytic input. |
| Power saving | Waibel 2021 v2, published 2022, is stronger than the previously used Rouse bound | Now improved to \(O_\epsilon(N^{7/8+\epsilon})\), by specialization and target-gcd averaging. |
| Linear same-field coordinate lower bound and mixed-field envelope | Exact local factors plus fixed-character Mertens | Immediate consequences once the asymptotic is proved. |
| Distinct-value lower and upper bounds | Corrected Henriot 2014 plus Hölder | Now \(N/\log^{\log4+\epsilon}N\) below and \(N J_{12}(N)/\log N\) above; classical inputs, new application here. |
| Fixed macroscopic balance | Harmonic theta theory plus the positive sieve tail | Classical equidistribution deduction; no shrinking-window theorem. |
| All eventual -8 targets and both strong cores | Plaksin 1984/1985, including its Gaussian specialization, an integral parity transformation and gcd deletion | Proved classical deduction below, even with prime-coordinate counts; no explicit numerical onset here. |
| Exact -8 identities, including strong cores | No matching primary antecedent identified in the searches below | An exact arithmetic candidate contribution, with historical priority unverified. |
| Exact -11 identity | A few direct norm/quadratic-residue searches and known shell/Hecke antecedents | No priority assessment as strong as the -8 search; do not extend a negative search claim to it. |
| Every-target unweighted leading asymptotic | Inspected general correlation theorems fail specific hypotheses | Still an independent obligation; no such asymptotic is established here. |
| \(C_{12}+C_{12}\) | Two independent norm conditions per summand | Outside this quaternary theorem and unchanged by these applications. |

## What the class-group conclusion establishes

The canonical proof now states AdditiveClassDisappearance explicitly:
for every fixed imaginary quadratic order of discriminant \(D<0\),
its primitive principal norm self-sum consists of all pair-primitive
locally admissible targets \(N\ge2\), apart from a finite exceptional
set. This does not require a restriction on the class group. The
sufficient local modulus is \(\prod_{p\mid2D}p^{2v_p(|D|)+1}\);
conductor primes are included.

GenusInvariantConvolution is a further immediate consequence. Replacing
either binary form by an integral form in the same genus leaves
every pair-primitive local density unchanged, hence leaves the full
leading coordinate convolution unchanged. The difference is
\(O_\epsilon(N^{7/8+\epsilon})\), and the two value sumsets differ
by finitely many targets. This identifies what disappears: finer
class obstructions within fixed genera. Local genus restrictions
remain.

These are explicit extractions of the already proved general
quaternary deduction. They do not extend the prime-factor semigroup
theorem to arbitrary class groups, and do not establish a matching
general \(N/\log N\) unweighted lower bound or an unweighted leading
asymptotic. The separate prime-discriminant and discriminant -56
arguments supply sharper unweighted results in their stated scopes.
No historical priority is asserted for these corollaries.

## Classical estimates actually applied

The [canonical analytic proof](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md)
now applies Waibel's Theorem 1 and Lemma 15 in the exact inspected
[arXiv:2012.04311v2](https://arxiv.org/pdf/2012.04311v2).
For the primitive scaled family \(u^2Q_1\oplus v^2Q_2\),
\[
 L\asymp(uv)^2,\qquad
 \Delta\asymp(uv)^4,\qquad
 \mathcal F(A,2)\asymp(uv)^4.
\]
The last comparison uses the good-prime Jordan exponents
\(0,0,2,2\); the finitely many bad-prime patterns are bounded uniformly.
The all-index estimate retains a target gcd. Averaging that gcd over
the squarefree Möbius divisors gives a truncated cusp error
\[
 O_\epsilon(N^{1/2+\epsilon}D^{3+\epsilon}
                   +N^{3/4+\epsilon}D^{1+\epsilon}).
\]
Together with the existing elementary tails, \(D=N^{1/8}\) proves
\[
 R^*(N)=JN\prod_p\beta_p^*(N)+O_\epsilon(N^{7/8+\epsilon}).
\]
Every fixed \(\eta<1/8\), in particular \(\eta=1/9\), is therefore
proved with the same local criterion. This is an application of an
existing analytic result and an elementary averaging step, not an
improvement to Waibel's theorem. Waibel's equations (1)--(2) also
independently confirm the real factor \(4\pi^2/\sqrt{\det A}\).

Henriot's **New Theorem 5**, printed p.377 of the
[2014 erratum](https://doi.org/10.1017/S0305004114000280), gives
uniform shifted fractional moments of maximal-order ideal-count
majorants for arbitrary fixed binary forms. It improves the general
distinct ordered value-pair bounds to
\[
 \frac{N}{(\log N)^{\log4+\epsilon}}
 \ll_{Q_1,Q_2,\epsilon}V(N)
 \ll_{Q_1,Q_2}\frac{N J_{12}(N)}{\log N},
\]
where
\[
 J_{12}(N)=
 \prod_{\substack{p\mid N,\ p\notin S\\
                    \chi_1(p)=\chi_2(p)=-1}}
                  \frac{p-1}{p-2}.
\]
The lower bound is for sufficiently large locally admissible targets;
the upper bound needs no asymptotic or admissibility assumption.
Both are uniform in target valuations. The canonical proof gives
the full growth-class, interval-height, Euler-factor and Hölder
checks. No lower-bound version of Henriot is used: its positive-weight
hypothesis fails for these functions at inert primes.

This also identifies the limit of this particular unweighting step.
For a fixed moment \(t>1\), the logarithmic loss is
\((2^t-2)/(t-1)>\log4>1\). Taking \(t\) toward 1 cannot reach a
matching \(N/\log N\) lower bound. A separate lower-count argument is
required for that scale.

## General-theorem antecedents and exact failed hypotheses

### 1 Brüdern: the closest explicit coprimality antecedent

Jörg Brüdern, Sums of Squares Coprime in Pairs, Bull. Polish Acad. Sci. Math. 62 (2014), 215--231, DOI 10.4064/ba62-3-3:
https://www.impan.pl/shop/en/publication/transaction/download/product/85773

Theorem 2 proves a power-saving asymptotic for four positive squares with all six pairwise coprimality restrictions. Lemmas 2,4,5 supply a large-common-factor tail, a uniform weighted four-square estimate, and absolute singular-series completion. Theorem 1 gives a linear lower bound on its admissible residue classes modulo 24.

This establishes that the sieve architecture and even stronger coprimality graphs were already used in the four-square case. It does not state the arbitrary-binary-pair theorem. Four-square all-six-gcd conditions do not imply two-block asymptotics on every two-block-admissible target: the admissible sets differ. Deleting four gcd restrictions is a fresh counting deduction, even for that special quadratic form.

Schinzel's On Sums of Four Coprime Squares (2013), cited by Brüdern, goes further than this comparison previously recorded: its final paragraph on printed p.111 explicitly imposes the two designated Gaussian block gcds on the entire eventual local range. It is a direct special-case antecedent, while remaining outside the arbitrary-form coordinate asymptotic. The targeted review above also records the complete Gaussian asymptotic deduction from Brüdern--Fouvry 1994.

### 2 Rouse and earlier quaternary theta theory

Rouse's 2018 preprint/2019 article gives uniform all-index coefficient estimates and eventual representation under whole-vector primitive local solubility. That primitivity means p does not divide all four coordinates. It does not impose both binary gcd conditions. Thus Rouse's representation theorem is not by itself an antecedent implying the desired statement. The two-block sieve is required.

Waibel improves the usable imported estimate, as proved in the canonical source. Baker already supplies a weaker power-saving route by the checked adaptation above. The present sharper exponent is a specialization and summation of classical estimates; it is not a new general theta estimate.

### 3 Browning--Heath-Brown geometric sieve

Tim Browning and Roger Heath-Brown, The geometric sieve for quadrics, arXiv:2003.09593v1 (21 March 2020), later Forum Math. 33(1) (2021), 147--165:
https://arxiv.org/pdf/2003.09593v1

Theorems 1.1--1.2 apply to homogeneous quadratic zeros of rank at least five and codimension-two subsets inside the projective quadric. They give a large-prime tail with a power saving in height. The desired positive four-variable shell is not a homogeneous zero locus of rank five. Homogenizing F(x)=N to F(x)-Nt^2=0 changes the quadratic form with N, while the theorem's constants depend on the fixed form; additionally, its count ranges over all t rather than the slice t=1. Neither issue is repaired by merely observing rank five.

This is a close conceptual antecedent, not a theorem with matching hypotheses. An N-uniform slice theorem would be an additional analytic result.

### 4 Cao--Huang affine quadrics

Yang Cao and Zhizhong Huang, Arithmetic purity of the Hardy-Littlewood property and geometric sieve for affine quadrics, arXiv:2003.07287v3 (10 December 2021), subsequently Advances in Mathematics 398 (2022):
https://arxiv.org/pdf/2003.07287v3

Their setup in Section 1.2.2 explicitly assumes indefinite q, fixed nonzero m, and a noncompact real locus. Theorem 1.4 proves the arithmetic-purity Hardy-Littlewood property for n>=4; Theorem 1.6 gives the geometric tail for such fixed affine quadrics. The codimension-two deletions corresponding to two block gcd conditions fit the geometry, but our positive-definite varying-N shells fail the indefinite/noncompact and fixed-variety hypotheses. A fixed positive shell has finitely many integral points as height grows.

The 2025 extension by Cao--Huang--Zhang to S-integral points has fixed q,m and a growing p-adic/archimedean height. It likewise cannot supply uniform asymptotics over arbitrary N by clearing denominators without a new uniformity argument.

### 5 Matthiesen correlations

Lilian Matthiesen, Linear correlations amongst numbers represented by positive definite binary quadratic forms, arXiv:1106.4690v3 (2012), Theorem 1.1:
https://arxiv.org/pdf/1106.4690v3

The theorem requires no two affine linear forms to be affinely dependent, and its error is o(N^d), without a stated power saving. The desired convolution uses forms n and N-n; these are affinely dependent. Thus this exact hypothesis fails. Making N another averaging variable can yield an averaged problem, not the requested uniform pointwise result for every target. The related 2013 paper extends the linear-correlation setting to indefinite binary forms; it does not remove this finite-complexity obstruction.

### 6 Browning--Dietmann and Jones circle-method estimates

Browning and Dietmann, On the representation of integers by quadratic forms, arXiv:math/0603358v2, Proposition 1 (PDF p.13), gives a uniform weighted asymptotic for n>=4 variables with positive target, with a specially chosen form-dependent weight:
https://arxiv.org/pdf/math/0603358v2
Its existence bounds and whole-vector local hypotheses do not already impose the two gcd conditions. It is a possible analytic input to a sieve, not an exact separately primitive theorem.

Edna Jones, The Kloosterman circle method and weighted representation numbers of quadratic forms, arXiv:2412.01780v2 (24 July 2025), Theorem 1.1 and Corollary 1.5, later Res. Math. Sci. 12,59 (2025):
https://arxiv.org/pdf/2412.01780v2
The theorem gives arbitrary smooth weights with explicit dependence on level and Hessian size. It counts unrestricted coordinates. The coprimality mentioned in its final sphere-packing section is a proposed application, not a theorem proved there.

### 7 Dietmann--Elsholtz/Greaves Gaussian special case

Dietmann--Elsholtz, Sums of Two Numbers Having Only Prime Factors Congruent to One Modulo Four, author preprint:
https://www.math.tugraz.at/~elsholtz/WWW/papers/greaves20.pdf
This explains that Greaves's 1976 weighted-quaternary result with two prime variables implies eventual Gaussian split-semigroup coverage. It concerns the particular prime semigroup 1 modulo 4, not arbitrary binary forms/ideal classes or a pair-primitive coordinate asymptotic. Its role is historical and special-case context. Its precise role is checked in the special-form comparison below.


## Classical special-form inputs and the exact -8 distinction

### 1. Liouville already gives the exact full-shell and whole-vector primitive formulas

J. Liouville, Sur la forme x^2+y^2+2(z^2+t^2), Journal de mathematiques pures et appliquees, series2,5(1860),269-272. Actual four-page article inspected, especially printedpp271-272. Primary scan:
https://www.numdam.org/item/JMPA_1860_2_5__269_0.pdf

After permuting coordinates, this is precisely Q plus Q for Q=x^2+2y^2. Liouville gives the unrestricted coefficients4,8,24 times sigma(m) at N=2^e m. He also gives exact counts for representations whose FOUR coordinates have gcd1: coefficients4,8,20,16 times m times product(1+1/p), for e=0,1,2,3, and zero for e>=4. His definition atp271 explicitly concerns a common divisor of all four coordinates. It does not concern the two prescribed binary gcds.

Thus both full-shell universality and the same dyadic obstruction16 are old facts for whole-vector primitive representation. The new statement cannot be described merely as primitive universality of this quaternary form; its content is separate primitivity of the two fixed binary blocks.

Own counterexample to the incorrect inference:13=Q(3,0)+Q(2,0) is primitive as a four-vector but neither binary block is primitive. Even re-pairing among coordinates with equal coefficients need not repair an arbitrary primitive shell point:157=2^2+3^2+2*6^2+2*6^2 has coordinate gcd1, while every pairing of a coefficient1 coordinate with a coefficient2 coordinate has gcd>1. These examples refute an implication between witnesses, not the existence theorem M8.

### 2. Deutsch is a direct quaternionic antecedent, with a different primitive definition

Jesse I. Deutsch, A Quaternionic Proof of the Representation Formula of a Quaternary Quadratic Form, J.Number Theory113(2005),149-179; actual arXiv:math/0406434v1,22June2004, inspected.
https://arxiv.org/pdf/math/0406434

Definitions27 and47 (PDFpp11,16) use the gcd of all four quaternion coefficients, with a primary normalization in the latter. Lemma28 identifies this with coordinate gcd at odd modulus. Theorems55-56 count primitive/primary quaternions; Section7 derives the same classical shell formulas and parity refinements. These theorems do not assert separate primitivity of the two norm components. They support the shell arithmetic and proof architecture, not the exact M8 conclusion.

Further primary corroboration, inspected online: Deutsch, Bumby’s Technique and a Result of Liouville on a Quadratic Form, Integers8(2)(2008),A2, equation(1.1) and introduction:
https://emis.de/ft/35897
This explicitly identifies the4,8,24 sigma formula as Liouville and develops another quaternion/ideal proof. No priority claim is justified for those shell formulas or the use of quaternion ideals to count them.

That paper's Definition 1 is a particularly close but different notion:
\(i\)-primitive means
\(\gcd_{\mathbb Z[\sqrt{-2}]}(\alpha,\beta)=1\).
It is not separate rational-coordinate primitivity of \(\alpha\)
and \(\beta\). For example \((\alpha,\beta)=(3,1)\) satisfies that
ring-gcd condition while the first rational coordinate pair \((3,0)\)
is imprimitive. Conversely
\(\alpha=\beta=1+\sqrt{-2}\) has both rational pairs primitive,
but has a common nonunit in the quadratic ring. The two conditions
are incomparable. The paper's resulting \(8\nmid N\) obstruction
for its primitive notion is also different from M8's \(16\nmid N\).
Thus its primitive representation count is not an exact antecedent
for the designated two-block condition.

### 3. Classical lifts from Plaksin cover both eventual strong cores

V.A.Plaksin, An asymptotic formula for the number of representations of a natural number by a pair of quadratic forms, the arguments of one of which are prime, Russian1984; English Math.USSR-Izvestiya25(1985),551-572, DOI10.1070/IM1985v025n03ABEH001306. Actual English PDF inspected: introductionp551 and main theoremp562, equation(18).
https://www.mathnet.ru/eng/im1517
https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1517&what=fullteng&option_lang=eng

The theorem counts N=phi(p,q)+psi(x,y), with p,q prime, phi primitive, and psi of fundamental discriminant. It gives an asymptotic under an explicit finite prime-unit congruence test. The fixed-form main term dominates N^(1/2+epsilon) on the admissible classes. These hypotheses permit phi=psi=Q=x^2+2y^2. They do not directly permit psi=2(x^2+y^2), whose content is2 and discriminant-16 is not fundamental.

Here is an independent deduction, making the gcd preservation explicit. Start from

    N=p^2+2q^2+x^2+2y^2
     = Q(p,y)+Q(x,q).

Discard solutions with p=2 or q=2; their count is O_epsilon(N^(1/2+epsilon)) by the fixed-binary representation bound. For odd p,q, failure of the first binary gcd requires p dividing y. For a fixed prime p there are O(sqrt(N)/p+1) such y; the remaining equation in x,q has O_epsilon(N^epsilon) solutions, even with primality dropped. Summing p<=sqrt(N) gives O_epsilon(N^(1/2+epsilon)) after enlarging epsilon. The same argument handles q dividing x. Therefore the asymptotic supplies representations with BOTH required gcds1.

Its exact local range is easy to check. At2, p and q odd force Q(p,q)=3 modulo8. Since Q modulo8 takes0,1,2,3,4,6, the allowable target residues are1,3,4,5,6,7 modulo8. Each has a prime-unit solution modulo8, and the derivative in p has valuation1, so univariate Hensel lifts it. At every odd prime, fix p=q=1: the nondegenerate binary Q represents every residue over the finite field, and a derivative in the prime variable is a unit, so lifting holds. This verifies the finite local test in the covered target classes.

For completeness, the singular-series main term really dominates the
discarded gcd failures in this specialization. At an odd prime
\(\ell\), put \(\chi=(-8/\ell)\) and
\[
 U_\ell(N)=\#\{u,v\in\mathbb F_\ell^\times:
                         u^2+2v^2=N\}.
\]
Plaksin's prime-unit local factor is
\[
 \theta_\ell(N)=1-\frac{\chi}{\ell}
              +\frac{\chi U_\ell(N)}{(\ell-1)^2}.
\]
If \(\ell\nmid N\), the binary finite-field count gives
\(U_\ell(N)=\ell+O(1)\), hence \(\theta_\ell=1+O(\ell^{-2})\).
If \(\ell\mid N\), it is \(1+1/\ell\) in the inert case, and
\(1-1/\ell+2/(\ell-1)>1\) in the split case.
All factors are positive and their product has a uniform positive
lower bound. The dyadic factor is uniformly positive in the listed
six residue classes by the same finite lifting. This checks the
main-term comparison without assuming it from an existence assertion.
Plaksin's actual modulus \(G\) here is also bounded (8 or 16),
since its prime product is over divisors of
\((3\Delta N,\tau)\) and \(\Delta=\tau=8\).
No coprimality of the two discriminants is required.

For odd N, the regrouped first norm is odd and the second even. Separate primitivity therefore puts them in C and2C. For even N in the allowed classes, x is odd, so both norms are in C. Consequently this classical theorem plus elementary deletion proves:

* every sufficiently large odd N belongs to C+2C;
* every sufficiently large N congruent to4 or6 modulo8 belongs to C+C.

The apparent missing class \(2\bmod8\) concerns only this first
arrangement. It can be filled without extending Plaksin to a
nonfundamental discriminant.

Apply Plaksin instead with
\(\phi=\psi=X^2+Y^2\), both primitive of fundamental discriminant \(-4\):
\[
 N=p^2+q^2+x^2+y^2.
\]
For even \(N\), its prime-unit local conditions allow exactly
\(N\equiv2,4,6\pmod8\). At every odd prime fix \(p=q=1\);
the Gaussian binary norm represents every residue, and a prime-variable
derivative is a unit. At 2, the two odd prime squares sum to 2 modulo 8,
while the free even Gaussian sum takes residues \(0,2,4\).
The derivative in \(p\) has valuation 1, so these solutions lift.
The odd-prime factors have the preceding formula with
\(\chi=(-4/\ell)\): they are \(1+O(\ell^{-2})\) at nontarget primes
and greater than 1 at target primes. The dyadic factors are 1 at
\(2,6\bmod8\), 2 at \(4\bmod8\). Thus the singular product again
has a uniform positive lower bound, independent of target valuations.

Discarding \(p=2\) or \(q=2\) costs
\(O_\epsilon(N^{1/2+\epsilon})\). In every remaining representation
of an even target, \(x\equiv y\pmod2\). Consequently the integral
bijection
\[
 u=(x+y)/2,\quad v=(x-y)/2,\qquad x=u+v,\quad y=u-v
\]
gives
\[
 x^2+y^2=2u^2+2v^2,\qquad
 N=Q(p,u)+Q(q,v).
\]
Each binary block contains a prime. Deleting \(p\mid u\) and \(q\mid v\)
costs \(O_\epsilon(N^{1/2+\epsilon})\), by exactly the fixed-binary
count and prime summation used above. Both remaining norms are odd,
positive, and primitive, hence lie in \(C\).
This proves the entire even strong-core statement eventually:
\[
 N\text{ even},\quad8\nmid N,\quad N\text{ sufficiently large}
       \ \Longrightarrow\ N\in C+C.
\]
For \(N\equiv2,6\pmod8\), respectively, \(u,v\) are both even or
both odd; for \(N\equiv4\pmod8\) they have opposite parity.
This independently confirms the residue range.

Together with the first arrangement's odd-target \(C+2C\) result,
this proves both eventual strong cores from classical inputs.
If \(v_2(N)=3\), represent \(N/2\equiv4\pmod8\) by \(C+C\)
and double both values. The transformation
\(Q(a,b)\mapsto Q(2b,a)=2Q(a,b)\) preserves primitive coordinates
when \(a\) is odd. Hence **every sufficiently large \(N\) with
\(16\nmid N\) lies in \(P_{-8}+P_{-8}\)** by this classical route.

There is a stronger coordinate consequence. Plaksin's main term
survives the deleted square-root error, giving
\(\gg N/\log^2N\) representations with a positive odd prime in each
primitive binary block: \(Q(p,u)+Q(q,v)\) on the even core, and
\(Q(p,y)+Q(x,q)\) at odd targets. Doubling gives the corresponding
count at valuation 3. These are coordinate counts, not distinct-value
asymptotics. No explicit numerical constant or onset is asserted here.

The exact M8 identities still require the explicit cutoff and complete
small-target bridge. Eventual existence, or even eventual coverage of
both strong cores, should therefore be attributed as a classical
deduction rather than presented as a new existence theorem. No
nonfundamental-discriminant extension of Plaksin is used.

### 4. Dietmann-Elsholtz, Greaves, Schinzel and Brudern: close classical architecture

R.Dietmann and C.Elsholtz, Sums of Two Numbers Having Only Prime Factors Congruent to One Modulo Four. Actual three-page author PDF inspected; theorem onp1, proof and remarkpp2-3:
https://www.math.tugraz.at/~elsholtz/WWW/papers/greaves20.pdf
Official edition:
https://www.impan.pl/shop/publication/transaction/download/product/113105

The theorem proves eventual coverage of2 modulo4 by two integers supported on primes1 modulo4. It invokes Greaves1976 prime-square-plus-prime-square-plus-two-squares, regroups with one prime coordinate in each Gaussian norm, then deletes the few witnesses where that prime divides its partner coordinate. Its remark distinguishes ineffective onset from the effective Kloosterman approach of Brudern. The Plaksin deduction above is a genuine weighted-norm analogue of this classical regrouping/deletion strategy. The Gaussian theorem itself has different prime support and target congruences, and does not imply exact M8.

Greaves primary citation supplied by that paper: On the representation of a number in the form x^2+y^2+p^2+q^2 where p,q are odd primes, Acta Arithmetica29(1976),257-274. Its theorem was not independently re-proved in this investigation. Do not silently attribute a weighted-coefficient version to it.

Brudern2014, Sums of Squares Coprime in Pairs, actual primary PDF from the repository, Theorem2p219 and Sections2-5, proves a power-saving four-square asymptotic with six pairwise gcd restrictions, using Mobius truncation, large-common-factor tails and singular-series completion:
https://www.impan.pl/shop/en/publication/transaction/download/product/85773
This is a stronger coprimality graph on a particular quadratic form, not an arbitrary-binary-pair theorem. Its effective eventual lower bound is not an exact classification of every eligible integer. The general bridge should present its architecture as classical rather than claiming a new sieve method.

### 5. Jones-Pall and ordinary ternary regularity do not supply the two gcds

B.W.Jones and G.Pall, Regular and semi-regular positive ternary quadratic forms, Acta Mathematica70(1939),165-191, DOI10.1007/BF02547347. Actual primary PDF inspected, especially introductionpp165-167 and quaternion proper-solution resultspp174-179:
https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/5599-11511_2006_Article_BF02547347.pdf

This develops ordinary regularity, proper whole-vector ternary solutions and selected congruence refinements. No theorem imposing the designated primitive binary subblock at every prime was found. Thus it cannot replace the M8 sieve merely by describing x^2+u^2+2v^2 as regular. The existing example29 makes the failure precise: that ternary value is represented, but every representation has v=0 and {|x|,|u|}={2,5}, so gcd(u,v)>1. The existing square-target equivalence with2t in C+C further prevents replacing the binary problem by an unproved exact primitive-subblock ternary classification.

### 6. Search coverage and limits

Searches used norm, semigroup, coordinate and historical formulations. Representative literal query groups included:

* x^2+2y^2 with relatively prime, primitive, coprime, sum of two;
* sum of two primitive norms, two primitive quadratic values, properly represented;
* prime factors1 or3 modulo8, integers composed of split primes, quadratic-character semigroups;
* quaternary two coprime pairs, relatively prime pairs, separately primitive, primitive representations two pairs;
* Pall/Jones/Watson with coprime, primitive, x^2+y^2+2z^2 and quaternary;
* French premiers entre eux and German teilerfremd versions;
* Greaves, Plaksin, prime-square coefficients, binary quadratic forms with two prime variables;
* quaternionic counting and parity papers of Liouville and Deutsch.

Dickson History volumesII(chapters6-8) andIII(chapter10) were downloaded as OCR and searched as bibliographic leads. They were not treated as proof of uninspected earlier primary results. References to Hermitian universality and primitive universal forms were checked for scope: ordinary representation, lattice primitivity, whole-vector primitivity and coordinatewise primitive binary blocks are distinct conditions. The Plaksin regrouping and integral parity change above do preserve both gcds after deleting a negligible set; they prove eventual existence, not the exact all-target classification.

These are substantial targeted searches, not an exhaustive bibliographic database review. Unindexed older literature, alternate terminology and inaccessible papers remain possible. The appropriate publication language is that no antecedent for the exact M8 classification has been identified in this audit; priority remains unverified. The meaningful candidate contribution is the exact separately primitive classification, supported by the explicit sieve and complete finite bridge, not the classical shell formulas, ordinary universality, norm-prime criterion or basic proof architecture.


## Unweighted correlations: close theorems that do not specialize

### Matthiesen 2012 and 2013

Lilian Matthiesen's *Linear correlations amongst numbers represented by positive definite binary quadratic forms*, Acta Arith. 154 (2012), 235--306, arXiv:1106.4690v3, counts representation functions under a pairwise affine-independence hypothesis on the linear forms. The 2013 extension, arXiv:1205.4145, preserves that condition. For the fixed-target problem, n and N-n are affinely dependent. The hypothesis fails; one cannot claim the desired pointwise theorem by specializing those results. Treating N as an extra variable produces an averaged-target theorem, whose error cannot be restricted to one target fiber.

Primary sources:
https://arxiv.org/pdf/1106.4690v3
https://www.impan.pl/shop/en/publication/transaction/download/product/83544

### Matthiesen 2020: an important closer unweighted antecedent

*Linear correlations of multiplicative functions*, Proc. London Math. Soc. 121 (2020), 372--425, DOI 10.1112/plms.12309, explicitly includes the sums-of-two-squares SUPPORT INDICATOR as Example (3) (arXiv v4 PDF p.5). It is therefore inaccurate to dismiss all Matthiesen correlation work as representation-weighted only.

However Theorems 2.1 and 2.4 retain pairwise linearly independent nonconstant linear parts in at least two variables and positive-volume averaging regions. The one-variable pair n,N-n fails that hypothesis. A line a+b=N inside a two-dimensional region has zero area, and the two-dimensional error estimate gives no relative information on one fiber. Her theorem is a close classical antecedent for unweighted averaged patterns, not the required every-target convolution.

Published primary article:
https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12309
Retained primary preprint:
https://arxiv.org/pdf/1606.04482v4

### Pretentious correlation theorems

Oleksiy Klurman, *Correlations of multiplicative functions and applications*, arXiv:1603.08453v1, advertises counts of a+b with a,b in multiplicative sets. Its useful polynomial-correlation asymptotic (Theorem 1.3, Corollary 1.4) requires finite pretentious distance to a character times n^{it}. The norm/allowed-prime support indicator vanishes at every good inert prime. For EVERY proposed twist its squared distance is at least sum_{p inert}1/p=infinity. Thus this essential hypothesis fails. An o(N) result would in any case not resolve a main term of scale N/log N.

Primary source:
https://arxiv.org/pdf/1603.08453v1

### Recent weighted-ratio work

V. V. Iudelevich, *On the sum of two squares function*, Izvestiya: Mathematics 90:3 (2026), 590--653, Theorem 1.2, proves a weighted ratio asymptotic for sum'_{n<=x}r(n)/r(n+1), of order x/(log x)^{3/4}. It does not assert an asymptotic for support indicators or for the fixed-target primitive convolution.

Primary article, exact theorem and full text:
https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=9681&wshow=paper

### Landau/Bernays and the remaining implication

One-variable Landau/Bernays type density results identify the scale N/sqrt(log N) for represented-value sets. That marginal scale does not determine the correlation along a+b=N, its target-local factors, or uniformity in N. For a fixed nonprincipal form, support is not generally multiplicative. None of the inspected theorems justifies replacing that support by the full split-prime semigroup.

The repository's subgroup semigroup theorem does give N J_12(N)/log N ORDER for its class-number-one unramified cores. It is a separate argument, not an immediate consequence of the quaternary count or the Henriot upper bound. Its coefficients are restricted to allowed semigroups; multiplying one summand by a bad prime or retaining arbitrary ramified strata needs a fresh hypothesis/application proof.

No general every-target distinct-value leading asymptotic was identified in the inspected sources. This is a precise non-applicability record, not a proof that no relevant literature exists and not a priority claim.



## Publication scope

The primitive-norm publication plan consists of the following two
papers. These are proposed article scopes; no new manuscript bundle,
submission or external review is asserted by this plan.

### Paper 1: Exact Sumsets of Primitive Quadratic Norms

The main result is the three exact primitive self-sum identities.
For \(P_D=\{Q_D(x,y)>0:\gcd(x,y)=1\}\), they are:

| \(D\) and \(Q_D\) | Exact sumset | Admissible exceptions |
|---|---|---|
| \(-7,\ x^2+xy+2y^2\) | \(\{N\ge2:49\nmid N\}\setminus\{7,77,203\}\) | \(7,77,203\) |
| \(-8,\ x^2+2y^2\) | \(\{N\ge2:16\nmid N\}\) | none |
| \(-11,\ x^2+xy+3y^2\) | \(\{N\ge2:2\mid N,\ 121\nmid N\}\) | none |

The article should explain how the primitive local obstructions,
stronger split-core identities and exceptional-quotient descent
produce these complete sets. Its central contribution is the
exact classification and the explicit proof joining finite
certificates to an infinite sieve, including the exceptional
targets in discriminant \(-7\).

Use the three current proofs as the authoritative mathematical
sources: [discriminant -7](PRIMITIVE_NORM_MINUS7_EXACT.md),
[discriminant -8](PRIMITIVE_NORM_MINUS8_EXACT.md) and
[discriminant -11](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md).
The [common classification record](PRIMITIVE_NORM_EXACT_CLASSIFICATIONS.md)
distinguishes these global identities from the bounded lists for
other discriminants.

A coherent proof order is:

1. State all three identities and prove the exact local exclusions.
2. Give the primitive-norm and split-core reductions, including the
   different roles of 2 and the ramified prime.
3. Present the shared quaternionic/spectral sieve, preserving each
   order's class-mass and local normalizations.
4. Prove that the finite certificates, intermediate rational
   intervals and infinite bounds cover every target with no gap.
5. Complete the exceptional-quotient descent and verify the actual
   exclusions \(7,77,203\).

The certificate supplement should contain the existing exact inputs,
verifiers and receipts: complete finite coverage through \(2^{24}\),
\(2^{27}\) and \(2^{25}\), respectively; 257, 107 and 35 intermediate
interval certificates; the additional 17 higher bins in the \(-7\)
proof; and the precise handoffs to the proved infinite estimates.
The mathematical proof must explain why those finite domains and
intervals exhaust what remains. A passing finite verifier does not
prove the imported spectral theorems.

The explicit \(N^{2/3}\) value-count lower bounds can be secondary
corollaries, with their recorded cutoffs and constants. They must
remain distinct from Paper 2's sharp \(NJ_D(N)/\log N\) scale.
The exact \(-19,-43,-67,-163\) exception lists remain bounded
computations and are not added to the main theorem.

### Paper 2: Primitive Norm Sums Beyond Class Number One

The centerpiece is now the
[all-order sharp value-count theorem](PRIMITIVE_NORM_ALL_ORDERS.md).
For every fixed imaginary quadratic order, with arbitrary ring class
group and conductor, every sufficiently large pair-primitive locally
admissible target has
\[
 \#\{a+b=N:a,b\in P_{\mathcal O}\}
             \asymp_{\mathcal O}
                   \frac{N J_{\mathcal O}(N)}{\log N}.
\]
These are ordered pairs of numerical values. Here \(J_{\mathcal O}\)
omits every prime dividing \(2\operatorname{disc}(\mathcal O)\);
the full conductor conditions are in an explicit finite local criterion.

The quantitative comparison with the individual local shadow has
a logarithmic saving \(\delta<1/h(\mathcal O)\), uniformly in every
target valuation, and has zero error when \(2\operatorname{Pic}(\mathcal O)=0\).
The proof uses generic good-class-prime saturation and an exact-one
prime-square sieve. Its conductor step preserves primitive invertible
local ideal types; it does not assign a unique genus to a numerical
conductor norm. Discriminant \(-100\), value \(25\), provides a useful
example of that distinction.

A coherent exposition can begin with the genuine individual class
failures at [discriminant -23](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md),
then pass through the periodic genus filters at
[discriminant -20](PRIMITIVE_NORM_CLASS_NUMBER_TWO.md) and
[discriminant -56](PRIMITIVE_NORM_CLASS_NUMBER_FOUR.md).
The [prime-discriminant odd-class theorem](PRIMITIVE_NORM_ODD_CLASS_GROUPS.md)
and [all-maximal-order theorem](PRIMITIVE_NORM_MAXIMAL_ORDERS.md)
are proved precursors. The general proof can then present the
ring-class norm invariant, local value profiles, finite-character
positive unweighting and lower-order class-error comparison.

The conductor tables include all dyadic types and the split-prime
tails. The additive lower bound requires only finitely many bounded
valuation profiles, with a proved uniform cap; it never assumes that
all individual conductor valuations are bounded. Class-prime
distribution enters only for a fixed order. The additive analytic
input remains the existing quadratic-character theorem, without
a new nonabelian distribution estimate.

The [pointwise correlation theorem](PRIMITIVE_NORM_BINARY_CORRELATION.md)
now strengthens this to an explicit leading asymptotic, using the
specified recent graph and fixed zero-free-half-plane theorems.
Keep those analytic inputs visible in any manuscript statement of its
proof method. The unconditional order-of-magnitude argument remains
an independently useful, older-input result.
Numerical general cutoffs and complete global exception lists for
\(-20,-23,-56\) also remain open. Historical priority and specialist
review of the exact value-count theorem are separate publication tasks.

### Supporting Baker/Waibel theory

The separately primitive quaternary theorem belongs in supporting
theory. Include the needed parts as an appendix or background
section of the relevant paper. A separate concise methodological
note is an option if its independent purpose is clear; it is not
the headline breakthrough of either planned paper.

Credit [Baker's convergent-sieve proof and its explicit first-power
adaptation](PRIMITIVE_NORM_ANALYTIC_BRIDGE.md#explicit-first-power-adaptation-of-baker),
the older four-square coprimality work, and Waibel's theta estimate
where used. Genus-invariant coordinate convolution is an immediate
consequence of that theory. Neither a coordinate asymptotic nor
qualitative existence establishes Paper 2's sharp value counts.

### Pointwise leading terms and the scope of the new application

The [leading-constant note](PRIMITIVE_NORM_LEADING_CONSTANT.md)
proves the conductor-prefix tail, Haar normalization and class-error
transfer. The [pointwise binary theorem](PRIMITIVE_NORM_BINARY_CORRELATION.md)
now supplies its binary correlation input, giving
\[
 R_P(N)\sim
 \frac{4c_{\mathcal C}^2}{(\varphi(L)/L)^2}
 \left(\prod_{p\mid2\Delta}d_p(N)\right)
 C_{\rm in}J_{\mathcal O}(N)\frac N{\log N}.
\]
The \(d_p(N)\) are primitive local VALUE convolution densities.
This is an asymptotic for numerical value pairs.

The same theorem treats arbitrary fixed abelian splitting subgroups
on both sides, including \(C_{12}+C_{12}\), with its explicit
\(N/(\log N)^{3/2}\) main term. It uses the weighted-dilation
root-moment theorem and the fixed Dirichlet zero-free half-plane from
Part I of the cited preprint. Their applicable written proofs and
hypotheses were inspected in
[the input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md).
Part II's stronger zero-free boundary is not needed.

The new application requires a projected indicator model, a quotient
character frame for fixed progression banks, a target-sized-shift
transfer, and a sparse-relative target-grid estimate. These steps
are supplied explicitly. The older positive-restoration sieve and
Jiang--Lau--Lü's fixed-shift prime-plus-norm theorem do not by themselves
yield the result. A targeted historical comparison of this exact
pointwise statement remains necessary; neither the underlying graph
method nor the zero-free theorem is claimed as a new contribution here.


The exact quartic-semigroup historical comparison is now available in
[the full Iwaniec--Pomykała audit](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#full-iwaniec--pomykała-theorem-audit).
Their Theorem 2 already treats coprime \(C_{12}+C_{12}\) representations
on every large admissible target, conditionally on a growing-shift
convolution estimate; it gives the expected logarithmic order.
Theorem 1 is unconditional for differences. The present leading
asymptotic also gives a [coprime corollary by Möbius inversion](SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md#coprime-summands),
with the additional allowed-prime target factor and relative error
\(O((\log N)^{-1/65536})\). Its target factor is uniformly comparable
to the weight in their theorem. The comparison now concerns the same
coprime count, its leading constant, quantitative remainder and analytic
inputs; it does not certify first-result status.

Historical priority, specialist mathematical review and complete
manuscript preparation remain separate publication obligations for
both planned papers. No first-result claim follows from the absence
of a matching source in the current searches.

## Evidence and limits

The searches covered norm and semigroup formulations, designated
binary gcds, primitive whole vectors, positive binary representation
correlations, geometric and arithmetic-purity sieves, uniform
quaternary estimates, prime-coordinate results, older ternary
regularity, and French/German terminology. Actual primary statements
were read for the sources used in the comparisons. Old Dickson OCR
was used as a bibliographic lead, not as a substitute for inspecting
a claimed primary theorem.

The actually inspected Waibel full text is arXiv v2. Publisher
metadata confirms the 2022 article, but access HTML prevented
inspection of the journal PDF; journal theorem numbering is not
claimed checked. No correction was located in targeted title/author
searches, which is not proof of absence of a correction.

The [existing source ledger](primitive_norm_quaternary_source_ledger.json)
retains exact editions and hashes, identifying mathematical inputs,
classical deductions, and comparison-only sources. The new exponent,
fractional-moment lower bound, and upper bound were separately
reread internally after integration. They are written deductions;
a finite experiment is not used to prove the analytic uniformity.

The review found a positive classical deduction for the general
two-block asymptotic. It did not identify an exact printed designated-block
statement, and it does not certify priority of the remaining quantitative
or exact arithmetic refinements. No exhaustive bibliographic database
review, external human assessment or formal verification is asserted.
