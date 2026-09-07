# Two-adic and real images: the generic subgroup still fills the product

Adding 2 and the real place to the [five-odd-prime comparison](BAD_PRIME_QUOTIENT_SUPPORT.md)
does not expose any exceptional quotient direction. **On all six fibres,
the generic subgroup fills the complete product of local point images at**
\[
S=\{2,3,5,7,11,13,\infty\}.
\]
The observed \(+9,+9,+7\) directions have simultaneous generic corrections
at all seven places. This result closes the declared small-place test,
not the full Selmer problem.

The [protocol](DYADIC_REAL_PROTOCOL.json),
[portable inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_dyadic_real_inputs_v1.json)
and [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_dyadic_real_support_v1.json)
retain the exact computations. Replay:

    python3 elliptic-curves/rank-jump/dyadic_real_support.py check
    sage -python elliptic-curves/rank-jump/dyadic_real_support.py verify

The second command independently checks real components with exact algebraic
roots and verifies the corrected class products with a separate local-power
calculation at every prime ideal above the six finite primes.

## What was tested

The same six curves and independent point subgroups were retained. No new
parameter, point, or field/class-group search was run. Six sequential workers,
each bounded by 30 seconds, produced separate checkpoints. The existing
local-squareclass implementation was imported without its shared cache.

At 2, an order maximal at 2 suffices. For a prime ideal \(P\) of ramification
index \(e\), the characters are valuation parity and the even cyclic
components of
\[
(\mathcal O_K/P^{2e+1})^\times.
\]
The omitted principal units are squares by Hensel's lemma: at \(z=1\),
the valuation condition for \(z^2-(1+\epsilon)\) is
\(v(\epsilon)>2v(2)=2e\). Thus these are full local squareclasses, not
low-precision proxies.

For a cubic with \(k\) local factors,
\[
\dim E(\mathbf Q_2)/2E(\mathbf Q_2)=k.
\]
Compared with an odd prime, the one-dimensional formal group contributes
one additional dimension. All generic image ranks reach this known value.

At the real place, a positive cubic discriminant gives one component bit
in \(E(\mathbf R)/2E(\mathbf R)\); a negative discriminant gives zero.
For a verified point on \(y^2=x^3+Ax+B\) with three real roots, the bounded
component is identified exactly by
\[
x<0\quad\text{or}\quad 3x^2<-A.
\]
Indeed the largest root is positive, and the positive critical point
\(\sqrt{-A/3}\) separates the bounded oval from the unbounded component.
The point-membership condition excludes the intervening interval where
the cubic is negative. The implementation checks the criterion against
exact algebraic roots, with no floating-point sign decisions.

## Paired results

The dimensions at 2 and infinity are complete local point-image dimensions.
The generic subgroup attains them and the full seven-place product
dimension \(s\) in every row.

| Fibre | Generic \(m\) | Observed \(q\) | At 2 | At infinity | Seven-place \(s\) | New quotient support |
|---|---:|---:|---:|---:|---:|---:|
| MW16-05, \(307/206\) | 16 | 9 | 2 | 1 | 9 | 0 |
| MW16-05, \(-3158/1291\) | 16 | 0 | 1 | 0 | 7 | 0 |
| MW16-04, \(-1647/91\) | 16 | 9 | 3 | 0 | 8 | 0 |
| MW16-04, \(-2177/2397\) | 16 | 0 | 2 | 0 | 6 | 0 |
| published-R17, \(-2300/843\) | 17 | 7 | 2 | 0 | 9 | 0 |
| published-R17, \(-1561/3133\) | 17 | 0 | 2 | 1 | 9 | 0 |

The original matching limitations and censored-zero interpretation remain
unchanged. These are separate curves, not a common identified global class
space on which to subtract Selmer dimensions.

The contrasts matter:

- MW16-05 has both an additional dyadic factor and an additional real
  component on the high fibre, but the generic points already account for
  both local images.
- MW16-04 has completely split two-torsion over \(\mathbf Q_2\) on the high
  fibre and only one rational nonzero two-torsion point on the control.
  The extra local dimension again lies in the generic image.
- The published-R17 control, rather than the high fibre, has the additional
  real component. Their seven-place product dimensions are equal.

Therefore neither real disconnectedness nor the number of dyadic factors
uniformly distinguishes the high fibres. More strongly, their complete
local images supply no residual class after quotienting by the generic
product image.

## Exact simultaneous normalization

The report gives new correction masks for all 25 exceptional directions
that work at all seven places together. They replace the earlier
five-prime corrections for this enlarged local set. They preserve the
same exceptional quotient directions and introduce no new rank.

The global-kernel identity from the preceding note consequently holds
with the enlarged \(S\). Writing \(V=\operatorname{Sel}_2(E)\),
\(G=\delta(M)\), and \(K_S\) for the kernel of localization at \(S\),
\[
V/G\cong K_S/(G\cap K_S).
\]
This applies to all Selmer classes because the generic subgroup alone
surjects onto the product of the allowed local images.

| Fibre | Generic kernel dimension \(m-s\) | Kernel dimension in the rational witness space \(m+q-s\) |
|---|---:|---:|
| MW16-05 high | 7 | 16 |
| MW16-05 control | 9 | 9 |
| MW16-04 high | 8 | 17 |
| MW16-04 control | 10 | 10 |
| published-R17 high | 8 | 15 |
| published-R17 control | 8 | 8 |

The full dimension of \(K_S\) is UNKNOWN. Its quotient by the generic kernel
can still contain unobserved rational directions or Sha classes. Trivial
localization here does not mean that every global or every local condition
has been solved.

## Decision for rank-jump understanding

This is an **incidence diagnostic excluding the declared local-support
mechanism**. Generic surjectivity is available without exceptional points.
It is not a positive rank predictor and must not discard any of these
high-gain fibres.

The small-place experiment is complete. Continuing to add convenient
primes would risk turning a precise hypothesis test into an open-ended
dictionary search. The remaining issue is global: a construction of several
independent classes and a reason their genus-one covers are rationally
soluble together.

The next useful bounded test should return to the retained geometric
constructions. The complete 39,120 rational-bisection census already
exists, and its failure on most rank-28 quotient directions is proved.
It should not be rebuilt or replaced by another trace shell.

Instead, examine only the branch quadratics of the **already split**
bisections on the historic controls. Distinct quadratic characters can
still have a low-genus common fibre product if their branch points overlap;
the earlier sixteen-cover test of identical quadratic characters did not
test that mechanism. Computing their character rank, branch union and
fibre-product genus would determine whether simultaneous lifting there
admits a low-genus common auxiliary curve. This needs no parameter search
or construction of a high-degree number field. Any positive block would
still need its section-independence proof; any negative result is limited
to the retained bisection construction.
