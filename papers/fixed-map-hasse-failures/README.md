# Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map

This is the focused sequel to
[`papers/common-arithmetic-fibers`](../common-arithmetic-fibers/).

The claim-by-claim assurance boundary is recorded in
[`VERIFICATION.md`](VERIFICATION.md). External mathematical inputs and the
remaining priority-review boundary are recorded in
[`LITERATURE_AUDIT.md`](LITERATURE_AUDIT.md).

The first paper lets the map vary and realizes every finite étale algebra of
rank at least three. This paper reverses the quantifiers: it fixes one
Jacobian-one map of geometric degree five and proves that a multiplicatively
large family of rational targets has full fibers

```text
Spec Q[X]/((X^3-a)(X^2+X+1))
```

that are everywhere locally soluble but have no rational point. For the
constructed family, the target count through height `B` is asymptotic to a
positive constant times `B / sqrt(log B)`. The paper also records the
arithmetic lower bound showing that geometric degree five is optimal.

The exact verification commands are:

```bash
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.FixedHassePaperCertificate
```

The Lean development proves the centered inverse factorization,
separability and rank-five finite étaleness, the exact determinant-one
normalization, the complete literal fiber representation, and the uniform
rational/real/all-`p`-adic Hasse certificate. It also proves target
primitivity, the exact height formula `H(y_a)=32a`, and target distinctness.
The multiplicative congruence/support conditions are closed under products
and powers in Lean. The noncube condition is kept separate because the
admissible set itself is not closed under multiplication. Every prime
`ℓ = 1 mod 9` is proved to supply the complete parameter certificate,
including exclusion of rational cubes.
The Selberg--Delange asymptotic and analytic input to degree minimality
remain ordinary mathematical proofs.

Build the manuscript from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The bounded enumeration is a regression check. The local theorem is proved
uniformly for every admissible parameter, and the asymptotic comes from the
Euler products and Selberg--Delange.
