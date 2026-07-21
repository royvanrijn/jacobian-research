# Sources consulted

Accessed 20--21 July 2026. Because the announcement is same-day, sources are split
between exact certificates and background; none should yet be treated as settled
historical attribution.

See [PROVENANCE_AUDIT.md](PROVENANCE_AUDIT.md) for the dated announcement
trail, the timestamp discrepancy between same-day sources, and the discovery
materials that have not yet been located.

- [The Jacobian counterexample, explained](https://jacobianfun.org/jacobian-explained):
  exact map, checks, weighted coordinates, inverse cubic, and proposed seed
  family.
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
