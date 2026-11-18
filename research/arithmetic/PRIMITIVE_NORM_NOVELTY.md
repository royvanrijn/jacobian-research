# Novelty comparison and applications of classical estimates

Assessment date: 18 November 2026. Scope: the general pair-primitive
binary-sum theorem, its quantitative consequences, and the exact
discriminant -8 classification. The -11 specialization is supplementary
context; its historical priority has not received the same depth of
search here.

The general proof is best presented as a structured deduction from
classical quaternary estimates and a known coprimality sieve. No exact
published antecedent for its full arbitrary-binary-pair statement was
located, but that does not make the method or analytic estimates new.
All eventual -8 existence, including both strong cores, also follows
from Plaksin by the parity transformation below. The exact all-target
classification is the stronger candidate contribution.
It must be distinguished from Liouville's nineteenth-century
whole-vector primitive theorem.

## Assessment of each claim

| Claim | Closest verified knowledge | Defensible assessment |
|---|---|---|
| Primitive norm prime-factor criterion | Classical quadratic ideal factorization; Euclidean argument in the class-one cases | Classical ingredient. |
| Full -8 shell and whole-vector primitive counts | Liouville 1860; Deutsch 2004/2005 and 2008 | Explicitly classical, including the whole-vector dyadic obstruction 16. |
| Möbius truncation, large-factor tail and singular-series completion | Brüdern 2014, Theorem 2, Lemmas 2,4,5 | Classical architecture; a new sieve method should not be claimed. |
| Arbitrary fixed binary pair, each block primitive, every locally admissible target | Rouse/Waibel coefficient estimates plus the two-block sieve; no exact statement antecedent located | A useful application/deduction, with priority unresolved; no new theta theorem. |
| All imaginary quadratic orders and invariance within fixed genera | The preceding arbitrary-binary-pair theorem and integral local isometries | Explicit class-group consequence, with no class-number restriction; an immediate corollary, not a new analytic input. |
| Power saving | Waibel 2021 v2, published 2022, is stronger than the previously used Rouse bound | Now improved to \(O_\epsilon(N^{7/8+\epsilon})\), by specialization and target-gcd averaging. |
| Linear same-field coordinate lower bound and mixed-field envelope | Exact local factors plus fixed-character Mertens | Immediate consequences once the asymptotic is proved. |
| Distinct-value lower and upper bounds | Corrected Henriot 2014 plus Hölder | Now \(N/\log^{\log4+\epsilon}N\) below and \(N J_{12}(N)/\log N\) above; classical inputs, new application here. |
| Fixed macroscopic balance | Harmonic theta theory plus the positive sieve tail | Classical equidistribution deduction; no shrinking-window theorem. |
| All eventual -8 targets and both strong cores | Plaksin 1984/1985, including its Gaussian specialization, an integral parity transformation and gcd deletion | Proved classical deduction below, even with prime-coordinate counts; no explicit numerical onset here. |
| Exact -8 identities, including strong cores | No matching primary antecedent identified in the searches below | Main candidate contribution, with historical priority unverified. |
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

Schinzel's On Sums of Four Coprime Squares (2013), cited by Brüdern, is a prior existence result in the same special setting; it is not an arbitrary-form representation theorem.

### 2 Rouse and earlier quaternary theta theory

Rouse's 2018 preprint/2019 article gives uniform all-index coefficient estimates and eventual representation under whole-vector primitive local solubility. That primitivity means p does not divide all four coordinates. It does not impose both binary gcd conditions. Thus Rouse's representation theorem is not by itself an antecedent implying the desired statement. The two-block sieve is required.

Waibel improves the usable imported estimate, as proved in the canonical source. The strongest defensible description is therefore a separately primitive deduction from classical quaternary estimates, with a new application/summation if no exact antecedent is identified. It should not be presented as a new general theta estimate.

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

A defensible presentation should lead with the exact -8 theorem and
its two strong-core identities if historical checking continues to
support their originality. Present the broad pair-primitive asymptotic
as a separately useful deduction from classical theta estimates,
with the Waibel-based error and exact local factors. The effective
sieve and complete interval/finite bridge are a separate part of the
contribution; the numerical onset alone is no longer the main
existence statement once every admissible target is classified.

Attribute the unrestricted shell formula to Liouville and its
quaternionic derivations to the inspected classical sources. Attribute
the coprimality-sieve architecture to Brüdern and distinguish the
Plaksin eventual routes for every admissible -8 target, including the
Gaussian parity transformation; Dietmann--Elsholtz supplies the close
Gaussian regrouping/deletion precedent. The local densities,
Mertens corollaries and corrected-Henriot applications should be
described as deductions, not as independent new analytic methods.

A cautious contribution sentence is:

> We deduce a pair-primitive local-global asymptotic for fixed positive
> binary quadratic forms from classical quaternary estimates, and give
> a complete effective proof of the exact primitive discriminant -8
> norm sumset, including its two strong-core identities.

This wording does not assert priority. Before an explicit first-result
claim, the remaining work is a bibliographic MathSciNet/zbMATH search
under older proper-representation terminology and, ideally, a specialist
comparison. Those steps were not performed here. No messages were sent
to authors or other researchers.

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

The negative novelty finding is limited: **no exact antecedent was
identified**, rather than a claim that none exists. No exhaustive
bibliographic database review, external human assessment, formal
verification or publication priority is asserted.
