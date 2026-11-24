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
