# Literature on ordinary reduction with a CM field action

Let \(A\) be an abelian variety of dimension \(g=eh\) over a number field,
with a Rosati-compatible action by a CM field \(E\) of degree \(2e\).
The comparison concerns the simultaneous bounds on the \(E\)-relative
Frobenius coefficient norms and split Newton polygons after finite
extension. It treats the action as an embedding
\(E\hookrightarrow\operatorname{End}^0_{\overline K}(A)\).
The comparison records differences in theorem scope and leaves priority open.

## Ordinary-reduction conjecture

Fité, Section 1, formulates positive-density ordinary reduction over
the original number field.
Pink, Section 7, and Cantoral Farfán, Li, Mantovan, Pries and Tang,
Section 1.1, formulate density one after finite extension. Pink also
records the sharper prediction that the connectedness field suffices.
These formulations are conventionally attributed to Serre. Here the
finite extension is allowed to contain the normal closure of the acting
field and torsion and roots of unity of an auxiliary level.
[Fité, Publ. Mat. 68 (2024), 27–40, §1](https://arxiv.org/html/2005.10185v3#S1),
[Pink, J. Reine Angew. Math. 495 (1998), 187–237, §7](https://people.math.ethz.ch/~pinkri/ftp/AMT-v3.pdf#page=45),
[Cantoral Farfán et al., arXiv:2508.11174v1, §1.1](https://arxiv.org/html/2508.11174v1#S1.SS1).

## Frobenius bounds and exterior powers

**Ogus and Sawin.** Ogus's Proposition 2.7.1, attributed there to Katz,
uses torsion to make a weight-one integral trace nonzero and rules
out divisibility at large degree-one primes by the Weil bound.
Proposition 2.7.2 adjoins roots of unity too. It makes an integral
weight-two normalized trace congruent to the representation dimension
modulo a sufficiently large level, forcing equality and hence a scalar
spectrum at each relevant prime. Corollary 2.9 applies this to abelian
surfaces. In Theorem 2.1 Sawin gives an exact density for abelian
surfaces from the normalized trace of the second exterior power of
the full cohomology and its constancy on monodromy components.
Corollary 2.2 rules out constancy on the identity component by Weil-bound
equality and the nonintegrality of the resulting weight-one trace.
The follow-up applies saturation to field norms, with the equality
contradiction using the common slope in a mixed *relative* component.
[Ogus, Hodge Cycles, Motives, and Shimura Varieties, LNM 900 (1982),
357–414, printed pp. 371–372](https://publications.ias.edu/sites/default/files/Number45.pdf#page=376),
[Sawin, C. R. Math. 354 (2016), 566–568, Theorem 2.1 and Corollary 2.2](https://www.numdam.org/item/10.1016/j.crma.2016.01.025.pdf).

**Fité.** Proposition 6 already bounds the rational \(p\)-valuation
of the coefficient-field norm of a weight-one trace by
\([E:\mathbb Q]/2\) on a density-one set after finite extension. Its
proof makes the trace nonzero using a torsion congruence and excludes
only valuations strictly above the weight threshold by magnitude.
Proposition 8 combines torsion and roots of unity with equality in
the weight-two bound when the trace itself belongs to \(p\mathcal O_E\).
Proposition 10 explicitly shows that Weil-bound equality for a
second exterior power of relative dimension at least three forces
coincident weight-one eigenvalues, and then excludes this possibility
at all but finitely many primes. Lemma 3.11 of Cantoral Farfán et al.
instead compares and bounds a selected relative exterior coefficient,
whose norm forms their invariant in equation (3.2) and
Proposition 3.12. At primes totally split in their acting field this
selected coefficient is the relative trace. These are direct precedents
for the individual ingredients of the simultaneous norm argument.
[Fité, §2, Propositions 6, 8 and 10](https://arxiv.org/html/2005.10185v3#S2),
[Cantoral Farfán et al., §3.2, Lemma 3.11 and Proposition 3.12](https://arxiv.org/html/2508.11174v1#S3.SS2).

## Exact comparisons

**Rank two and partial ordinarity.** Fité's Theorem 2 allows a
Galois number field of degree \(g\) inside the geometric endomorphism
algebra, with no CM or Rosati condition. The theorem states
\(p\)-rank at least \(\lceil g/2\rceil\) on a set of nonzero density
over the given field. The proof in Section 5 already obtains a
density-one set after finite extension with that bound. In the overlap
of Galois CM action and the present hypotheses, Proposition 6
and the local lower bound \(e-b+2r\) for the trace-norm valuation
imply \(r\leq\lfloor b/2\rfloor\), where \(b\) counts balanced
signature pairs and \(r\) counts nonordinary balanced pairs.
The strict norm estimate gives \(r\leq\lfloor(b-1)/2\rfloor\) for
\(b>0\). This improves the inferred signature-sensitive bound by
one pair when \(b\) is positive and even. Both signature-sensitive
bounds coincide for odd \(b\) in that overlap. The strict result
also covers CM fields that are not Galois over \(\mathbb Q\).
When \(b=0\), the split Hodge signatures themselves force ordinarity.
[Fité, Theorem 2 and its proof in §5](https://arxiv.org/html/2005.10185v3#S5).

**A single mixed pair of minimum multiplicity one.** Assumptions 2.1
and 2.5 of Cantoral Farfán et al. require the full geometric
endomorphism algebra to be a CM field and one conjugate signature
pair \((1,h-1)\), with all remaining pairs definite. Theorem 4.18
proves density-one ordinarity after finite extension in arbitrary
relative rank under those hypotheses. Remark 4.21 also gives the
conclusion for products of varieties satisfying their assumptions.
The split norm criterion includes the single pair \((1,h-1)\) with
an embedded Rosati-compatible action, regardless of other geometric
endomorphisms. For \(h=2\) it also permits two balanced pairs.
Theorem 4.18 determines the ordinary density over the original field
in its Galois relative-rank-at-least-two case, and Theorem 4.19 gives
density-one \(\mu\)-ordinarity over that field when the acting field
is abelian. Those conclusions address different questions from eventual
ordinarity at degree-one primes.
[Cantoral Farfán et al., Theorems 4.18 and 4.19, Remark 4.21](https://arxiv.org/html/2508.11174v1#S4.SS5).

**Relative rank four and one \((2,2)\) pair.** Fité's Theorem 4
allows geometric multiplication by an imaginary quadratic field
of signature \((2,2)\). Although the statement gives nonzero
density over the original field, Section 8 explicitly obtains
density one after finite extension by Proposition 10. This covers
the \(e=1\) case, and his statement imposes neither a full
endomorphism algebra condition nor a Rosati hypothesis on the
embedding. His Corollary 5 gives additional cases with a
noncommutative simple subalgebra and a specified Schur index.
For \(e>1\), the present criterion allows one \((2,2)\) pair and
arbitrarily many \((0,4)\) pairs for the action of \(E\), including
in dimension \(4e\). This \((2,2)\) pair does not meet the simple
signature assumption of Cantoral Farfán et al. Earlier theorems
may still apply to particular varieties or their isogeny factors.
[Fité, Theorems 1 and 4 and Corollary 5](https://arxiv.org/html/2005.10185v3#S1),
[Fité, proof of Theorem 4](https://arxiv.org/html/2005.10185v3#S8).

**Other relative-rank-two results.** Wang and Zhang's
Theorems 1.1–1.2 concern specified prime and twice-prime
dimensions over \(\mathbb Q\) and prime dimensions over quadratic
number fields. Theorem 1.3 applies in general dimension under
absolute simplicity, a degree-\(g\) endomorphism subfield, a
specified splitting behavior at an auxiliary rational prime,
and exclusion of Albert type IV. Their proofs use a normalized
square of a relative Frobenius trace norm with residue splitting
and large-image hypotheses, leading to positive density over the
original field in their setting. Suh's Theorem 4.1.1 and its
proof involve trace norms over Hecke fields for parallel
weight-two non-CM Hilbert modular forms, with arguments for
partial proximity to ordinariness. These results do not state
simultaneous arbitrary-signature bounds for split CM components.
[Wang and Zhang, Mathematika 72 (2026), e70081, §1](https://arxiv.org/html/2503.21111v2#S1),
[Suh, Compos. Math. 156 (2020), 647–678, Theorem 4.1.1](https://arxiv.org/html/2410.01182v1#S4.SS1).

**Further monodromy hypotheses.** Pink's Theorem 7.1 applies to
varieties with trivial geometric endomorphism ring under a
restriction on the root system. Thakur's Theorem 1.1 assumes a
specified noncommutative geometric division algebra, a specified
Mumford–Tate representation, and connected monodromy. The
arbitrary-rank valuation estimate here uses the CM field action
and its signatures without these group hypotheses.
[Pink, §7](https://people.math.ethz.ch/~pinkri/ftp/AMT-v3.pdf#page=45),
[Thakur, Theorem 1.1 and Definition 1.1](https://arxiv.org/html/1803.05954v2#S1).

## Example from a cyclic cover

For a cyclic cover with branch exponents \((1,1,1,1,2,2)\) modulo
8, the primitive factor has field \(\mathbb Q(\zeta_8)\), relative
rank four and dimension eight. The Hodge dimensions at characters
\(1,3,5,7\) are \(0,2,2,4\), or in the opposite order with the
inverse character convention, by the cyclic-cover character
formula. The Jacobian polarization makes the deck transformation
unitary, giving the required Rosati compatibility. Moonen's only
conductor-8 row has four branch points, and the six-branch family
is outside the examples in Cantoral Farfán et al., Corollary 5.1.
This provides a uniform application of the rank-four criterion.
The calculation uses the given cyclotomic action and makes no
claim about the generic full geometric endomorphism algebra.
[Moonen, Doc. Math. 15 (2010), 793–819, §2 and Table 1](https://pure.uva.nl/ws/files/949135/83991_Moonen_DocMath_2010.pdf),
[Cantoral Farfán et al., Corollary 5.1](https://arxiv.org/html/2508.11174v1#S5.SS1).

## Extremal polygons and the numerical \(p\)-rank formula

Oort gives the usual integral-breakpoint, convex polygon convention for
\(p\)-divisible groups and studies Newton strata. Achter and Pries,
Definition 2.4, give the most generic admissible *symmetric* Newton
polygon of given height \(2g\) and \(p\)-rank \(f\). After the \(f\) zeros
and \(f\) ones, its slopes are \(1/(g-f)\) and \(1-1/(g-f)\), with the
appropriate one-segment cases for \(g-f\le2\). This agrees with the
fixed-endpoint result in the symmetric specialization. The oriented
split-CM problem allows different endpoint counts and different
interior dimension and codimension. Its additional conclusions are
the minimum of both rounded prefix sums simultaneously, the interval
formula for their excess, and the resulting closed minimum \(p\)-rank
for a collection of conjugate pairs subject to the Frobenius norm bounds.
Those conclusions apply to numerical polygons. Realizability for a
given abelian variety involves a separate question.
[Oort, Astérisque 298 (2005), 255–269, §1.1](https://www.numdam.org/item/AST_2005__298__255_0.pdf#page=4),
[Achter and Pries, Definition 2.4](https://arxiv.org/pdf/1311.5846#page=6).

## Projective trace and the three-balanced-pair case

The function \(\operatorname{tr}(g)^2/\det(g)-1\) is the adjoint
trace of a rank-two representation. Pink's trace characterization
describes the field of a compact subgroup in an adjoint representation,
and Arias-de-Reyna, Dieulefait and Wiese use inner twists to describe
the projective field of definition for a suitable irreducible
representation. In Ribet's argument for modular forms, the proof of
Theorem 3.1 uses squared traces on a finite-index subgroup where the
determinant is rational, so these traces encode the same field.
These results give precedents for projective trace fields. The
ordinary-reduction criterion here uses equality of projective trace
functions in a different way. It uses valuation parity to keep balanced
and definite pairs in separate equivalence classes and the number
three to make the balanced failures simultaneous.
[Pink, J. Algebra 206 (1998), 438–504, Proposition 0.6](https://people.math.ethz.ch/~pinkri/ftp/LastVersion.pdf),
[Arias-de-Reyna, Dieulefait and Wiese, Theorem 2.15 and Remark 2.16](https://arxiv.org/pdf/1203.6546#page=11),
[Ribet, On \(\ell\)-adic representations attached to modular forms II,
§3, Theorem 3.1](https://math.berkeley.edu/~ribet/Articles/rankin.pdf).

For the remaining cases with independent rank-two factors, the
closest number-field results include Fité's Theorem 3, requiring that
the cubic field acting on a threefold is not contained in a field of
definition of the action, and Wang–Zhang's Theorem 1.3, requiring a
specified inert or two-prime splitting behavior and excluding Albert
type IV. Fité's §6 and Wang–Zhang's §4 exploit nonsplit coefficient
places. The criterion based on one Frobenius orbit of all three
balanced pairs has the same divisibility mechanism and its own
orbit hypothesis. Noot's Theorem 2.2 addresses an odd tensor power of
standard two-dimensional factors. The direct sum in the unlinked
case lies outside that representation hypothesis. The newer density-one
results of Can and Lin concern the transcendental Hodge structure
of K3 surfaces with an abelian endomorphism field and additional
hypotheses. Hui proves density zero for supersingular reduction of
a non-CM variety, leaving Newton polygons with just one nonordinary
balanced pair outside that conclusion. Zhu's density result applies
to a hyperelliptic family at its generic point over a rational
function field. These results leave the fixed split-prime threefold
question in the reduction of the three-balanced-pair problem open.
[Fité, Theorem 3 and §6](https://arxiv.org/html/2005.10185v3#S6),
[Wang–Zhang, Theorems 1.3 and 1.7–1.8](https://arxiv.org/html/2503.21111v2#S1),
[Noot, Compos. Math. 97 (1995), 161–171, Theorem 2.2](https://www.numdam.org/item/CM_1995__97_1-2_161_0.pdf#page=9),
[Can–Lin, arXiv:2609.34093v1, Theorem 1.1](https://arxiv.org/html/2609.34093v1#S1),
[Hui, arXiv:2504.08088v2, Theorem 1.1](https://arxiv.org/html/2504.08088v2#S1),
[Zhu, arXiv:2608.22779v1, Theorem 1.5](https://arxiv.org/html/2608.22779v1#S1).

Bhatnagar and Fu's Theorem 1.4 classifies real-multiplication isogeny
classes over finite fields with split residue characteristic and
nonzero \(p\)-rank, and includes ramified almost-ordinary factors.
This is a precedent for the type of finite-field examples in the
first residual case. The specific construction here varies the
rational prime while imposing a fixed auxiliary torsion level and
retaining exactly one failing conjugate pair.
[Bhatnagar and Fu, arXiv:2308.11105v3, Theorem 1.4](https://arxiv.org/html/2308.11105v3#S1).

## Selected Hasse–Witt entries and moment denominators

Bouw, Definition 2.1, Lemma 2.2 and §5, gives characterwise
Hasse–Witt invariants and coefficient formulas for cyclic covers.
Lin, Mantovan and Singhal, Theorem 4.6 and §5, give rank tests
along Frobenius orbits of the extended Hasse–Witt matrix and compute
the matrix for cyclic covers via the two-chart Čech complex. Their
generic ordinariness result has hypotheses on the characteristic
and the branch family. Vlasenko, Theorem 1(i), represents iterated
Hasse–Witt operations modulo \(p\) by coefficients of the
\(p^s-1\)-st power of a Laurent polynomial in the applicable
hypersurface setting. Huang, Lian, Yau and Yu, Theorem 1.2 and §5,
relate Hasse–Witt coefficients and period integrals in the
toric-hypersurface and flag-variety setting. The present moment
argument uses a prime-power denominator to isolate the coefficient
on one chosen coherent eigenline, including examples where the
full primitive factor is nonordinary. Its deduction is a
nonvanishing statement for sphere moments of the specified
initial form and its lower-weight perturbations.
[Bouw, Compos. Math. 126 (2001), 295–322, §§2 and 5](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6D1BEC9200F6849339665AC4371AB7AF/S0010437X01000343a.pdf/div-class-title-the-p-rank-of-ramified-covers-of-curves-div.pdf),
[Lin, Mantovan and Singhal, arXiv:2303.13350v3, Theorem 4.6 and §5](https://arxiv.org/html/2303.13350v3#S4),
[Vlasenko, arXiv:1605.06440v3, Theorem 1](https://arxiv.org/pdf/1605.06440v3),
[Huang, Lian, Yau and Yu, arXiv:1801.01189, Theorem 1.2 and §5](https://arxiv.org/html/1801.01189).
