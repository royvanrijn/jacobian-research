# Non-thin rank jumps for the published R17 fibration

Date: 2026-08-31

## Theorem

Let

```text
pi_R17 : X -> P1_Q
```

be the published rootless rank-17 elliptic fibration on the certified H3 K3
surface. If `U` is the smooth-fibre locus of `pi_R17`, then

```text
{t in U(Q) : rank X_t(Q) > 17}
```

is not thin in `P1(Q)`. Equivalently, the rational fibres of rank at least 18
form a non-thin set. Adding or removing the finite singular-fibre set does not
change this conclusion.

## Hypothesis audit

We apply Pasten--Salgado, Theorem 1.1, with `K=Q` and `pi=pi_R17`.

1. **Elliptic K3 over `Q`.**
   [`EC-K3-ELKIES-2026-R17`](../MATH_STATUS.json) identifies the compact
   published model over `Q` with the exact q12/orbit5867 endpoint.
   [`EC-K3-H3-Q12O5867-ENDPOINT-QQ`](../MATH_STATUS.json) proves that this
   endpoint is the certified H3 K3 surface and that its full Mordell--Weil
   rank is exactly 17.
2. **Non-isotriviality and no non-reduced fibres.**
   [`EC-K3-H3-Q12O5867-QQ-ROOTLESS`](../MATH_STATUS.json) proves that the
   singular fibres are geometrically `24I1`. Thus every fibre is reduced. An
   `I1` fibre gives a pole of the `j`-map, so the fibration is non-isotrivial.
3. **A different elliptic fibration over `Q`.**
   [`EC-K3-H3-SOURCE`](../MATH_STATUS.json) supplies the H3
   `E7+E8/MW2` Jacobian fibration over `Q`, while the endpoint source-identity
   certificate places it on this same K3 surface. It is different from
   `pi_R17`: its reducible-fibre root lattice has rank 15, whereas `pi_R17`
   is rootless.
4. **Zariski density.**
   The published fibration has seventeen independent `Q(t)`-sections. Their
   multiples give infinitely many distinct rational sections. A proper closed
   subset of a surface contains only finitely many irreducible curves, so it
   cannot contain all these sections. Hence `X(Q)` is Zariski dense.

All hypotheses of Pasten--Salgado's theorem therefore hold. Its implication
from Zariski density to non-thin rank jumps gives

```text
{t in U(Q) : rank X_t(Q) > rank MW(X,pi_R17)}
```

non-thin in `P1(Q)`. The exact endpoint certificate gives
`rank MW(X,pi_R17)=17`, proving the claim. A finite subset of `P1(Q)` is thin,
and thin sets are closed under finite unions; consequently the choice of
convention at the 24 singular parameters is immaterial. QED.

## Scope

This is a theorem about rational specializations of the published R17
fibration. It proves neither the exact rank of an individual fibre nor an
effective height bound for finding the rank jumps. It is independent of the
field of definition of the alternate rootless rank-17 candidate: the second
fibration used here is the already exact H3 fibration over `Q`.

## External theorem

H. Pasten and C. Salgado,
[*Non-thin rank jumps for double elliptic K3 surfaces*](https://doi.org/10.1007/s00229-024-01554-2),
*Manuscripta Mathematica* **175** (2024), 771--781, Theorem 1.1.
