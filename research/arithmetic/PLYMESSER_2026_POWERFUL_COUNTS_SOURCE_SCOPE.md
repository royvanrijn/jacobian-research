# Latest powerful-count literature and the prescribed-sum uniformity gate

Audit date: 24 November 2026. This note records source identities, stated estimates and the limits of their transfer to r(N), the number of ordered positive square-full pairs with sum N. It introduces no new theorem about r(N), and does not certify historical priority or a best-known bound.

## 1. The quantities behind the exponents

| Source | Actual quantity and estimate | Uniformity relevant to a prescribed sum |
|---|---|---|
| [Chan, *Twin squareful numbers* (2012), Theorems 3--4](https://doi.org/10.1017/S144678871200016X) | Pairs n,n+l with n<=x: O(d3(l)x^(2/5)log²x); for fixed l, O_l(x^(7/19)logx). | The 7/19 constant depends on l. This is not the same statement as the repository's uniform prescribed-sum 7/19 theorem. |
| [Blomer--Schoebel, *Twins of powerful numbers* (2013), Theorem 1, Corollary 2, Theorem 3](https://doi.org/10.7169/facm/2013.49.2.12) | Theorem 1 is uniform in the nonzero difference, with exponent 2/(2k+1). Corollary 2 explicitly gives the same exponent for positive prescribed sums; k=2 gives 2/5. Theorem 3 gives gamma2>61/180 for a fixed nonzero difference. | The improved fixed-difference theorem absorbs the difference into the approximation error and constants. Corollary 2 comes from the coarser sign-insensitive argument, not automatically from Theorem 3. |
| [Reuss, arXiv:1212.3150v2, Theorem 4](https://arxiv.org/abs/1212.3150v2) | Consecutive square-full pairs n,n+1: O_epsilon(x^(29/100+epsilon)). | The general determinant Theorem 1 explicitly permits constants depending on h; the application has h=1. Its exponent cannot simply be substituted for a growing prescribed target. |
| [Heath-Brown, arXiv:2601.07817v1, Theorems 1.1 and 1.3](https://arxiv.org/abs/2601.07817v1) | The main result counts ALL square-full triples u+v=w<=B, including nonprimitive ones, with exponent 3/5-3/1555+epsilon. The auxiliary theorem bounds primitive integer zeros of aX³+bY³+cZ³ by O_epsilon(H^epsilon), uniformly in all nonzero signed integer coefficients. | The global triple count requires w itself to be square-full and is not a bound for arbitrary r(N). The cubic theorem does transfer to a fixed nonzero-coefficient slice, with one coordinate 1 ensuring primitivity. It bounds completions, not the frequency of varying slices. |

The repository already retains the complete distinction between prescribed sums, fixed differences and representable-value support counts. Its current [3/10 argument](ERDOS_943_DUAL_KERNEL_GRAPHS.md) is stronger than its earlier 7/19 argument. The [full kernel-pair](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale) and [six-cycle](ERDOS_943_SIX_CYCLES.md) savings remain separate obligations.

## 2. A newly catalogued consecutive-number manuscript

Danny Plymesser's [*Counting powerful products of consecutive integers*](https://hexagonmath.org/2611.00602v1) is **hexagon:2611.00602v1**, shown by the current Hexagon version-history page as 19 November 2026. Its linked manuscript is version 1.0.1 dated 24 October 2026; the exact [GitHub v1.0.1 release](https://github.com/dannyply/math-frontier/tree/v1.0.1) and [Zenodo record](https://doi.org/10.5281/zenodo.23615335) are the author-linked evidence. This is not an arXiv identifier: [arXiv:2611.00602](https://arxiv.org/abs/2611.00602) identifies a different paper.

The manuscript claims x^(1/4+epsilon) for consecutive powerful pairs and x^(77/320+epsilon) for starts of powerful products of at least three consecutive integers, including the union over all widths. These are improvements claimed for consecutive-number counts, not r(N). It discloses that the work and its reviews were produced with AI systems, with no human mathematical review.

Its new geometric ingredient is Lemma 5.5: a primitive solution on an irreducible auxiliary curve Gamma in P1xP1, of bidegree (alpha,beta), with a bihomogeneous (3,2) equation F=h!=0, maps with fibres of size at most 2 to integral points on an affine curve of degree d=2alpha+3beta and coordinate height at most max(1,L1³T1²). The classical affine-curve estimate then gives O(d³B^(1/d)(logB+d)). This is a useful conditional method once an appropriate auxiliary-curve cover is available.

## 3. What was checked in that method, and what remains an input

The complete written component/root argument, rational-component reduction, reduced-coordinate bounds, lifted-curve proof, interval summation and final assembly were read. In the core lift, the image ring is a domain; its degree-s filtration has dimension

    (3s+1)(2s+1)-(3s-alpha+1)(2s-beta+1)
       =(2alpha+3beta)s+(alpha+beta-alpha beta).

The actual counted coordinates are integral monomials, rather than rational coordinates with an untracked denominator. Primitivity restricts the two projective scalings to signs and leaves at most two lifts. These checks support the stated geometry. They do not replace the source's optimization certificates or a formal build.

The estimate used is precisely [Castryck--Cluckers--Dittmann--Nguyen, Theorem 3](https://arxiv.org/pdf/1904.13109v1), including its extension to Q-integral curves via Corollary 4.1.4 and Remark 4.1.5. [Stewart, Theorem 1 and its consequence on printed pp. 795--796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf), allows content-one binary forms with nonzero discriminant, including reducible ones, and gives 2800 d^(1+omega(|h|)) primitive solutions. The corresponding precise source statements were checked, not their entire deep proofs.

The actual tagged Lean declarations have explicit literature hypotheses:

* fixed widths k>=3: Stewart, CCDN in ambient dimension 12, and Gelfond's coefficient-height inequality;
* consecutive pairs k=2: those three, plus Reuss at h=1 and Heath-Brown's anisotropic Theorem 15;
* the union over widths: also Tao's short-interval Lemma 3.1.

This is a 3/5/6-input distinction, rather than an unconditional kernel proof of every cited theorem. Gelfond's exact source is Bombieri--Gubler Lemma 1.6.11, printed p. 27; the factor-height consequence follows because the other nonzero integer factor has height at least 1. Tao's [Lemma 3.1](https://arxiv.org/html/2603.27990v3#S3) bounds a powerful-product interval's width by exp(log^(2/3+o(1))N).

The author's registry explicitly leaves a source-comparison question concerning a logarithmic factor in the Heath-Brown manuscript versus its proof sketch; the published LNM text was not compared there. A fixed extra power of log x would be absorbed in x^epsilon in the mathematical exponent argument, but that observation does not itself change the formal premise. This audit did not run Lean, Nanoda, the statement comparator, or the exact optimization certificates, and confers none of those assurances.

Most importantly, the determinant cover is a fixed-h construction, and the reduced coordinates require x>=3|h|M. The paper's uniform growing-shift Lemma 8.2 is only for |h|<=x^(delta/4), with a stated small fixed delta and an exponent loss. This does not include a difference comparable to x, nor the general signed shifts Nh in the balanced compatibility problem. The lifted-curve lemma itself is uniform in nonzero h; the missing step for prescribed sums is a suitable uniform auxiliary-curve cover and summation.

## 4. The 20 November quadratic companion

The same author's [*Powerful values of quadratic polynomials*](https://hexagonmath.org/2611.00699v1), **hexagon:2611.00699v1**, links [Zenodo version 1.0.1](https://doi.org/10.5281/zenodo.24348789) and [GitHub v1.1.1](https://github.com/dannyply/math-frontier/tree/v1.1.1). The manuscript is dated 18 November and was catalogued on 20 November 2026. Its statements and coefficients were inspected; its complete 34-page proof was not audited here. The companion reports a Lean formalization with three literature inputs, but this audit did not inspect or replay that companion formalization.

For a FIXED integer quadratic A n²+B n+C, A>0 and nonsquare discriminant, Theorem 1.1 claims exponent 293/1220+epsilon, with a coefficient-dependent constant. Theorem 1.2 states

    S(a,b,c;X) <= K L^nu X^(389/1600),
    a>=1, |a|,|b|,|c|<=L, b²-4ac nonsquare,

with absolute K,nu, but nu is not computed. Its coefficient-family consequence only extends to height X^eta with eta=1/[1600(nu+3)]. It does not state a coefficient-free bound.

Fixing one square-free kernel in the sum problem produces f(a)=N-b³a². This has a negative leading coefficient and N-varying coefficient height, so it does not satisfy those displayed hypotheses. The polynomial height loss cannot be discarded or relabelled as X^epsilon. The quadratic companion is therefore another method lead, not a direct prescribed-sum improvement verified by this audit.

## 5. Other current primary sources

The current arXiv records checked on the audit date still list Heath-Brown 2601.07817, Browning--Verzobio 2608.24512 and Trifonov 2609.17715 as v1. [Browning--Verzobio](https://arxiv.org/abs/2608.24512v1) counts primitive Campana triples with three fullness orders; its Theorem 1.2 is genuinely uniform in nonzero signed coefficients and rectangular side lengths >=2. These are different quantities from r(N). Its slicing comparison with Heath-Brown's Theorem 15 gives the anisotropic affine-curve bound min(X^(1/q),Y^(1/p)), already usable in unbalanced boxes but unchanged at the balanced cubic crossing.

[Trifonov's September 2026 paper](https://arxiv.org/abs/2609.17715v1) supplies a short-interval square-full result and a coefficient-uniform ternary conic estimate. Their [exact specialization and summation limits](TRIFONOV_2026_TERNARY_CONIC_SOURCE_SCOPE.md) are already retained: bounded completions for one pair/shift do not give a saving in the number of nonempty pairs.

A separate new method source is [Blomer--Pascadi, arXiv:2607.24311v1, Theorems 1.1 and 5.5](https://arxiv.org/html/2607.24311v1#S1.SS1). Theorem 1.1 controls bilinear forms in classical Kloosterman sums for ANY modulus, including cube moduli; it requires unit multiplier, interval lengths <=c, and the stated gcd condition (removable for initial intervals). At equal length sqrt(c) it saves c^(1/32); its nontrivial balanced range is c^(13/28+epsilon)<length<c^(7/12-epsilon). The existing completion phase A t³+B t^(-2) is not automatically a classical Kloosterman phase. At the current critical cube modulus c~T³, original lengths T=c^(1/3) and product lengths T²=c^(2/3) are outside that balanced range. An actual transformation, weights and length analysis are required; the unbalanced theorem is not an automatic replacement either. Only these precise statements and surrounding scope discussion were read, not the full paper proof.

The fresh search also found conditional square-full counting-function expansions and unrelated consecutive-powerful constructions. None of their inspected statements gives a uniform full prescribed-sum or long-complement frequency bound. This is a bounded literature audit, not a completeness or priority certificate.

## 6. Retained primary identities

The two new PDFs retained for this note are exact author-tagged release files, without claiming byte identity with another server's rendition:

* [Plymesser counting paper](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf), tag v1.0.1, SHA256 dbc3ac29ec50485a61c3f4119a06745c7d4475710701aeb0bb6f5597751e48a5.
* [Plymesser quadratic paper](sources/plymesser-powerful-quadratics-hexagon2611.00699-v1.0.1.pdf), repository tag v1.1.1, manuscript version 1.0.1, SHA256 1c1f06fa1f2b496582ceebbaf85fdf9ec85846a6a6363bb4a9a44dd80449178e.

Existing Chan, Reuss, Heath-Brown, Browning--Verzobio and Trifonov primary PDFs are reused at their existing paths and pins. No old source pin, principal state, gap or assurance is altered by this source-scope note.

## 7. Focused 7/19 recheck: Chan, Blomer--Schobel and Heath-Brown

Rechecked 25 November 2026 against the actual primary statements and the proof passages below. Write
\[
 S_l(X)=\#\{m\le X:m,m+l\ \text{positive square-full}\},\qquad
 r(N)=\#\{(u,v)>0:u,v\ \text{square-full},\ u+v=N\}.
\]
These are different counting problems. The first allows a constant depending on a fixed shift; the second needs a bound uniform in the prescribed target.

**Where the exponents belong.** [Chan's 2012 Theorem 4](https://doi.org/10.1017/S144678871200016X) gives S_l(X)<<_l X^(7/19)log X. Theorem 3 instead makes the shift dependence explicit in d_3(l)X^(2/5)log²X. In [Blomer--Schobel's author version](https://www.researchgate.net/profile/Anita-Schoebel-2/publication/268024970_Twins_of_powerful_numbers/links/67816b208210a977a17fb2d5/Twins-of-powerful-numbers.pdf), Theorem 1 is shift-uniform and Corollary 2 explicitly treats positive prescribed sums with exponent 2/(2k+1). Their Theorem 3 improves the fixed-difference square-full exponent to any gamma_2>61/180; it does not state that sharper exponent uniformly in a growing shift. [Reuss, exact v2, Theorems 1 and 4](https://arxiv.org/pdf/1212.3150v2) gives X^(29/100+epsilon) for the consecutive case; its general determinant estimate permits dependence on h. None of these sharper fixed-difference statements can be substituted directly for a uniform bound on r(N).

The two Chan papers must also be distinguished. [Chan, arXiv:1407.0054, Theorem 4 and Lemma 4](https://arxiv.org/pdf/1407.0054), is the coefficient-unit modular box-completion result used in the repository's own 7/19 proof. It is not the 2012 fixed-difference theorem with the same numerical exponent. Its explicit error is O_epsilon((K+L)q^(-1/2+epsilon)+q^(1/2+epsilon)).

**The exact Blomer--Schobel transfer gate.** Their Proposition 1 requires 0<Delta<1/2 and counts rational approximations with error Delta/Q². Equation (2.3), for their normalized nonzero difference d and X comparable to either powerful term, bounds the approximation error by
\[
 \frac{2|d|}{M_j^2}\frac{N_jM_j}{X}.
\]
Since Q=2M_j, the sufficient Delta supplied by this argument is proportional to |d|N_jM_j/X. The following displayed choice Delta asymp N_jM_j/X absorbs the fixed d into its constants, as their text explicitly permits. In an all-variables scale T with X asymp T^5, fixed d gives Delta asymp T^(-3). Letting |d| asymp X instead gives a supplied Delta of order T², outside Proposition 1. This explains the failure of the direct growing-shift substitution. It does not rule out a different auxiliary curve, a new parameterization, or a separate estimate for positive sums.

**The exact Heath-Brown bridge.** [Browning--Van Valckenborgh, arXiv:1106.4472v1, Theorem 2](https://arxiv.org/pdf/1106.4472v1), gives B^(3/5)log^12 B for primitive square-full triples. Its Section 4 combines diagonal-conic and diagonal-cubic slices; its discussion immediately after Theorem 2 identifies the all-six-variables scale B^(1/5) as the obstruction to improving that exponent. This is the same balanced scale that limits the repository's two-graph exponent optimization, although the two counting problems and the required frequency bounds differ. The official record still lists [arXiv:2601.07817v1](https://arxiv.org/abs/2601.07817v1). Theorem 1.1 counts all square-full triples u+v=w<=B, including nonprimitive ones, with exponent 3/5-3/1555+epsilon; it requires w itself to be square-full. Its Theorem 1.3 is the genuinely coefficient-uniform primitive diagonal-cubic estimate. Setting one coordinate equal to 1 makes it applicable to nonzero-coefficient binary cubic equations with prescribed right-hand side. The primitive hypothesis must also be checked when all three coordinates vary. Sections 4--6 use lifted bilinear congruences, rational lines and lattice determinants to retain more information than a raw modular box count. They offer a method to investigate, rather than an automatic pointwise consequence of Theorem 1.1.

[The August 2026 Browning--Verzobio preprint](https://arxiv.org/abs/2608.24512v1) still cites Heath-Brown for the square-full triple case. Its Theorem 1.2 is coefficient-uniform in rectangular generalized-Fermat boxes; its new main counting theorem concerns sufficiently large fullness orders. Neither displayed result closes the prescribed-sum cycle frequency problem here.

**Research consequence.** The repository's current unconditional pointwise claim is 3/10, established internally in [the dual-graph proof](ERDOS_943_DUAL_KERNEL_GRAPHS.md), rather than its earlier 7/19 bound. The productive transfer problem is to exploit the common target N when counting families of compatible slices, lines or cycles. Bounding each completed slice by N^epsilon does not bound the number of slices. The [critical rectangular criterion](ERDOS_943_CRITICAL_RECTANGULAR_CYCLE_CRITERION.md) records one precise sufficient cycle estimate for reaching 5/17, its equivalent residual, and the actual sectors already controlled using the uniform cubic theorem. That cycle estimate and the resulting global 5/17 bound remain open. No converse from a pointwise representation bound to that cycle estimate is asserted.

This focused search verifies the cited comparisons and identifies method gates. It does not establish a worldwide priority claim or completeness of the literature search. It adds no imported premise to the existing 7/19 or 3/10 proofs.
