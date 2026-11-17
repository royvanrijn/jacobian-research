# Gao F4: complete complex fibres and an elliptic stratum

The complete fibre classification requested for F4 in Gao's
[4 November revision](https://arxiv.org/html/2608.00222v2) is given here
by a written proof and exact arithmetic checks.

- **GAOF4FIB1:** every fibre has 0, 1, 2, 3 or 5 points; all five sizes
  occur. On the target boundary c=0 there are three points when
  `3(2v2+v4)^2-32v3 != 0`, and one otherwise. On c!=0 the full quintic
  multiplicity table determines every fibre.
- **GAOF4ELL1:** the 2+2+1 multiplicity stratum in weighted Y-space is
  a smooth elliptic curve minus thirteen points. Its target stratum
  is the product with C*. The 3+1+1 stratum is A1 minus thirteen points;
  4+1 and 3+2 each consist of five Y-points.
- **GAOF4NP1:** the nonproper-value locus has exactly two irreducible
  hypersurface components, including the entire target boundary. The
  other component is the c-saturated quintic discriminant. Its nonzero
  chart has discriminant normalization A2, explicitly parameterized.

[Proof and source-credit boundary](PROOF.md) · [checker](verify.py) ·
[receipt](receipt.json) · [source pins](sources/manifest.json).

Gao's simple-root lifting criterion and five omitted curves are prior
results. An earlier Ideosphere note proves selected exact fibres; it
does not provide the full classification in the inspected manuscript.
F6 now has a [separate complete classification](../gao-f6-fibres/README.md).
F5, F7 and JC2 remain unresolved. The proof has no external review,
formal verification or established novelty priority.

```sh
timeout 60s .venv/bin/python research/jacobian/gao-f4-fibres/verify.py
```

These are scoped, explicitly requested Jacobian follow-ups. They do not
restart the archived programme. MATH_STATUS.json remains authoritative.
