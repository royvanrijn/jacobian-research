# Sources consulted

Accessed 20--22 July 2026. Because the announcement is same-day, sources are split
between exact certificates and background; none should yet be treated as settled
historical attribution.

See [PROVENANCE_AUDIT.md](PROVENANCE_AUDIT.md) for the dated announcement
trail, the timestamp discrepancy between same-day sources, and the discovery
materials that have not yet been located.

## Plane degree-frontier sources

- [T. T. Moh, *On the Jacobian conjecture and the configurations of
  roots*](https://doi.org/10.1515/crll.1983.340.140), J. Reine Angew. Math.
  340 (1983), 140--212: the historical \(\max(\deg P,\deg Q)>100\)
  result.  The publisher exposes the article's opening scan but restricts the
  remaining text; later section-level claims were cross-checked through the
  exact citations and caveats in the 2014, 2017, and 2022 GGV/GGHV sources.
- [Raymond C. Heitmann, *On the Jacobian
  conjecture*](https://www.sciencedirect.com/science/article/pii/002240499090042T),
  J. Pure Appl. Algebra 64 (1990), 35--72, MR 1055020: primary source for
  \(\gcd(\deg P,\deg Q)\ge16\), monomial-valuation restrictions, and the
  same four computational cases.  Heitmann explicitly says his paper does not
  reprove Moh's reduction-of-degree step.

- [Zenodo record 21479814](https://doi.org/10.5281/zenodo.21479814), Billel
  Helali, *Exact Computer-Assisted Exclusion of the (72,108) Frontier in the
  Two-Dimensional Jacobian Problem*, version 1.0.1, 21 July 2026: manuscript,
  exact coefficient ideals, explicit membership/unit certificates, and replay
  code.  Exact archive/PDF hashes and the version comparison are in the
  [plane-JC provenance record](../../plane-jc/PROVENANCE.md).
- [Jorge Alberto Guccione, Juan José Guccione, Rodrigo Horruitiner, and
  Christian Valqui, *Increasing the degree of a possible counterexample to the
  Jacobian Conjecture from 100 to 108*](https://arxiv.org/abs/2204.14178):
  Theorem 2.1 gives the below-125 list and Proposition 4.3 gives the two
  remaining \((8,28)\) Laurent polygons.
- [The same four authors, *Some algorithms related to the Jacobian
  Conjecture*](https://arxiv.org/abs/1708.07936): admissible-complete-chain
  algorithm, \((m,n)\)-families, and the 34 oriented cases through maximum
  degree 150.
- [Jorge Alberto Guccione, Juan José Guccione, and Christian Valqui, *On the
  shape of possible counterexamples to the Jacobian Conjecture*](https://arxiv.org/abs/1401.1784):
  minimal-pair and standard-\((m,n)\)-pair reduction.
- [The same three authors, *The two-dimensional Jacobian conjecture and the
  lower side of the Newton polygon*](https://arxiv.org/abs/1605.09430) and
  [*A system of polynomial equations related to the Jacobian Conjecture*](https://arxiv.org/abs/1406.0886):
  lower-corner restrictions and coefficient-system methods used by the 2022
  reduction.
- [Jorge Alberto Guccione, Juan José Guccione, Rodrigo Horruitiner, and
  Christian Valqui, *The Jacobian Conjecture: Approximate roots and
  intersection numbers*](https://arxiv.org/abs/1708.09367): intersection
  inequality used for the degree-84 exclusion.

- [The Jacobian counterexample, explained](https://jacobianfun.org/jacobian-explained):
  exact map, checks, weighted coordinates, inverse cubic, and proposed seed
  family.
- [Omniscience Research Agent and Jeff Pickhardt, *An Explicit Counterexample to the Dixmier Conjecture in A_3*](https://omniscienceproject.com/papers/an-explicit-counterexample-to-the-dixmier-conjecture-in-a-3-jfLENtXF):
  same-dimension inverse-Jacobian lift of the foundational map to an explicit
  injective non-surjective Weyl-algebra endomorphism, with centralizer and
  differential-order proofs of non-surjectivity.
- [V. V. Bavula, *Holonomic modules and 1-generation in the Jacobian Conjecture*](https://arxiv.org/abs/2112.03177):
  inverse-Jacobian extension of Keller maps to Weyl-algebra endomorphisms and
  centralizer statements used in the classical implication from Dixmier to
  Jacobian.
- [A. Belov-Kanel and M. Kontsevich, *The Jacobian Conjecture is stably equivalent to the Dixmier Conjecture*](https://arxiv.org/abs/math/0512171):
  the dimension-preserving implication `DC_n => JC_n`, differential-order
  pullback argument, and the deeper stable converse.
- [Juntang Zhuang, `jzkay12/jacobian_conjecture`](https://github.com/jzkay12/jacobian_conjecture),
  pinned here at commit
  [`1ff68e870f66afec8c6611f910fcc8f5522fdbce`](https://github.com/jzkay12/jacobian_conjecture/commit/1ff68e870f66afec8c6611f910fcc8f5522fdbce):
  public compilation of the expanded `F4a`, `F4b`, and `F4c` maps, their
  Island A/B/C labels, constant-Jacobian checks, and rational collision
  certificates.  The repository's resolvent and canonical-boundary analysis
  is recorded separately in the [quartic-islands audit](../../extended-geometry/EXTERNAL_QUARTIC_ISLANDS.md).
- [Dean Cureton, `deancureton/jacobian`](https://github.com/deancureton/jacobian),
  pinned here at commit
  [`0d4a9212d874226ad81ce5a926becddfa94e6a88`](https://github.com/deancureton/jacobian/commit/0d4a9212d874226ad81ce5a926becddfa94e6a88):
  external Lean 4 formalization of the determinant, collisions, unit-Jacobian
  counterexample, and complex specialization. See [Lean foundational-map audit](../../verified/LEAN_FOUNDATIONAL_MAP.md) for
  theorem scope, reproduction, attribution, and the no-license source boundary.
- [Macaulay2 `CoincidentRootLoci` package](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/CoincidentRootLoci/html/toc.html):
  classical coincident-root ideals, parameterizations, tangent spaces, and
  singular loci used by the independent degree-five-through-eight comparison.
- [Jaydeep Chipalkatti, *On equations defining coincident root
  loci*](https://arxiv.org/abs/math/0110224): equations and singular-locus
  background for coincident-root strata.
- [L. M. Feher, A. Nemethi, and R. Rimanyi, *Coincident root loci of binary
  forms*](https://arxiv.org/abs/math/0311312): classical geometry and
  equivariant-class background.
- [Simon Kurmann, *Some remarks on equations defining coincident root
  loci*](https://doi.org/10.1016/j.jalgebra.2011.10.045): normalization and
  singular-locus background cited by the Macaulay2 package.
- [Stacks Project, *Weil divisors*](https://stacks.math.columbia.edu/tag/0BE0)
  and [Picard group of projective space](https://stacks.math.columbia.edu/tag/0BXJ):
  standard divisor-class and Picard inputs used in the marked-point dimension
  barrier; the product Picard group is generated by the two pullbacks.
- [Stacks Project, Lemma 58.21.4, *Purity of branch
  locus*](https://stacks.math.columbia.edu/tag/0BMB): the precise purity input
  used to make the universal `S_n` monodromy proof algebraic.
- [Stacks Project, *Universally injective, unramified
  morphisms*](https://stacks.math.columbia.edu/tag/06ND), together with its
  sections on [universally injective morphisms](https://stacks.math.columbia.edu/tag/01S2)
  and [universal homeomorphisms](https://stacks.math.columbia.edu/tag/04DC):
  scheme-theoretic bridge criteria used in the master quotient theorem Master Quotient Theorem.
- [Stacks Project, proper-image and curve inputs](https://stacks.math.columbia.edu/tag/01W6),
  together with [proper plus locally quasi-finite implies finite](https://stacks.math.columbia.edu/tag/0F2P)
  and [the genus of a curve](https://stacks.math.columbia.edu/tag/0BY6): the
  external closure, finiteness, and arithmetic-genus inputs isolated in the
  complete Generic Discriminant Theorem proof.
- [SGA 1, *Etale covers and the fundamental
  group*](https://firmaprim.github.io/sga/sga-1/), Expose XII, Theorem 5.1 and
  Corollary 5.2: Riemann-existence comparison used after descent to `C` to
  identify finite etale covers of affine space.
- [Bass, Connell, and Wright, The Jacobian Conjecture: Reduction of Degree and
  Formal Expansion of the Inverse](https://doi.org/10.1090/S0273-0979-1982-15032-7):
  classical reduction to cubic-homogeneous form.
- [L. Andrew Campbell, Reduction Theorems for the Strong Real Jacobian
  Conjecture](https://arxiv.org/abs/1303.3853): explicit four-step stable and
  Segre implementation used for the cubic-homogeneous artifact.
- [Gorni and Zampieri, On cubic-linear polynomial
  mappings](https://arxiv.org/abs/1204.4026): explicit pairing construction
  used for the Druzkowski artifact.
- [Zihan Zhang, Direct Consequences of the Three-Dimensional Counterexample](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html):
  immediate logical consequences and same-day source trail.
- [Derksen–van den Essen–Zhao, The Gaussian Moments Conjecture and the Jacobian Conjecture](https://arxiv.org/abs/1506.05192):
  proves the implication used in the Gaussian-moments consequence.
- [Wenhua Zhao, *Images of commuting differential operators of order one with
  constant leading coefficients*](https://arxiv.org/abs/0902.0210): the
  fixed-map Special Image Conjecture input used after the cubic-homogeneous
  BCW reduction.
- [Christopher D. Long, *Small Counterexamples to the Gaussian Moments
  Conjecture*](https://arxiv.org/abs/2607.18186), arXiv:2607.18186v1,
  submitted 20 July 2026: explicit three-variable Gaussian witness and a
  separate conservative BCW route from the announced JC map to
  `not GMC(158)`.  The paper explicitly says the direct witness was not
  derived from the JC map.  The repository now reproduces the route's exact
  18-step reduction and 79-variable cubic-homogeneous collision locally.
- [I. J. Good, *Generalizations to several variables of Lagrange's expansion,
  with applications to stochastic processes*](https://doi.org/10.1017/S0305004100034666),
  Proceedings of the Cambridge Philosophical Society 56 (1960), 367--380:
  multivariable inversion formula underlying Long's Gaussian fixed-point
  method.  The repository's weighted-seed bridge derives its needed special
  case locally from the coefficientwise Abhyankar--Gurjar identity.
- [Christopher D. Long, *Counterexamples to the (xz)-Conjecture and the
  Mathieu Conjecture for (SU(2))*](https://arxiv.org/abs/2607.19012),
  arXiv:2607.19012v1, submitted 21 July 2026: direct `(xz)` witness and its
  `SU(2)` lift.
- [M. Müger and L. Tuset, *The Mathieu conjecture for `SU(2)` reduced to an
  abelian conjecture*](https://doi.org/10.1016/j.indag.2023.10.001),
  Indagationes Mathematicae 35 (2024), 114--118: Lemma 2.1 is the
  Haar-integration formula used in Long's `SU(2)` lift.  The repository now
  also proves this formula directly from `SU(2)=S^3` and Hopf coordinates.
- [Zhao, Hessian Nilpotent Polynomials and the Jacobian Conjecture](https://arxiv.org/abs/math/0409534):
  equivalence with the quartic Vanishing Conjecture and Hessian-nilpotent
  formulation.
- [van den Essen–Wright–Zhao, On the Image Conjecture](https://arxiv.org/abs/1008.3962):
  records the implication from the Image Conjecture to the Vanishing Conjecture.
- [MathOverflow: explicit cubic model with S3 monodromy](https://mathoverflow.net/questions/513387/galois-structure-of-the-new-counterexample-to-the-jacobian-conjecture-an-explic):
  primitive element `t=y+1/x`, reconstruction, discriminant, and monodromy
  argument. This is a question/discussion, not a refereed paper.
- [Andy Jiang (`@davikrehalt`), geometric interpretation](https://x.com/davikrehalt/status/2079175065695035442):
  public post adopting the binary-cubic/marked-projective-root formulation.
  The status was posted 20 July 2026 at 12:02:19 UTC and located for this
  audit 20 July 2026 at 22:43 UTC (21 July 00:43 CEST). Its affine
  cubic/reconstruction and fiber consequences overlap repository claims
  cubic marked-root and image results; the projective marked-root organization is used here with
  attribution and without a priority claim.
- [Encyclopedia of Mathematics: Jacobian conjecture](https://encyclopediaofmath.org/wiki/Jacobian_conjecture):
  established formulation, reductions, and references predating the
  announcement.
