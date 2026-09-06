# The fresh103b2 jump requires at least nine strict rational directions

The fixed103b2 comparison now isolates where the additional classes must
occur. On the fresh rank-at-least27 fibre at **3726/881**, at most **one**
additional Selmer direction can come from the local boundary beyond the
17-dimensional generic image. Consequently at least **nine independent
rational directions are strict classes**, unramified everywhere in the cubic
field and split at every allowed bad prime. None belongs to the generic
subgroup, whose strict kernel is zero.

This is a rank-derived structural necessity, **not an independent numerical
class-group computation**. The boundary calculation uses only the equation
and generic sections. The existing rank lower bound enters only afterward.
No exceptional coordinate or exceptional-point-derived class is supplied.

The [arithmetic certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_class_boundary_v1.json)
and [independent verification and consequences](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_class_boundary_verification_v1.json)
cover exactly the two fixed fibres, with no point search or new parameter.

## Exact incidence accounting

Let S contain2, infinity and every bad prime, K be the cubic two-division
field, G the marked generic Kummer subgroup, and

\[
 U=\mathrm{Sel}_2^S(E),\quad
 c_S=\dim_{\mathbf F_2}\mathrm{Cl}(\mathcal O_{K,S_K})/2,\quad
 I=\operatorname{loc}_S\mathrm{Sel}_2(E).
\]

The [strict-class identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md) gives
\(\dim U=c_S\). Both fibres have \(\dim G=\dim\operatorname{loc}_S G=17\),
so \(G\cap U=0\). Therefore

\[
0\longrightarrow U\longrightarrow \mathrm{Sel}_2(E)/G
\longrightarrow I/\operatorname{loc}_S G\longrightarrow0.
\tag{1}
\]

The point-independent computation gives:

| Fixed103b2 fibre | Recorded rank ≥ | Local point product dimension ℓ | dim loc(G) | dim I | dim Sel₂/G |
|---|---:|---:|---:|---:|---:|
| 3726/881, new-20260906-71 | 27 | 19 | 17 | 17–18 | c_S+e, 0≤e≤1 |
| −1049/2296, matched observed-zero control | 17 | 20 | 17 | 17–19 | c_S+e, 0≤e≤2 |

Each row has its own c_S and e. Neither is identified across the two
nonisomorphic cubic fields. The low label is censored, not an upper bound
on its rank or its strict rational image.

### Why the boundary loses at least one dimension

For the monic cubic f with root θ and discriminant δ, define

\[
\beta=-\delta f'(\theta),\qquad N\beta=\delta^4.
\]

This is a global norm-square class defined by the equation. It is unramified
outside S: away from the polynomial discriminant it is a unit at odd primes;
at any polynomial-discriminant prime omitted from S, the certificate checks
even valuations directly. This handles nonminimal-model primes explicitly.

Both cubics have positive discriminant. At their three ordered real roots,
β has sign vector (−,+,−), while the real point image consists of
(+,+,+) and (+,−,−). The real Hilbert pairing with the nonzero point class
is nontrivial: the sign masks101 and011 have odd dot product.

Consequently

\[
 \lambda_\beta(x)=\sum_{v\in S}\langle\beta_v,x_v\rangle_v
\]

is a nonzero functional on the full local point product L_S. Global
reciprocity makes it vanish on I, so

\[
 17\le\dim I\le\ell-1.
\tag{2}
\]

This is the same [derivative reciprocity argument](DERIVATIVE_RECIPROCITY_AND_COMPLETE_BOUNDARY.md),
now applied before exceptional points are supplied to these fresh fibres.
β itself is locally inadmissible and is **not** an additional Selmer class.

The verifier also constructs abstract Lagrangian boundary completions
realizing every value in each displayed interval while retaining G and
one independent derivative constraint. This shows that these linear
conditions alone cannot select e. It does not assert that each abstract
completion occurs arithmetically on the fixed fibre.

## The necessary soluble block on the high fibre

Let W be the full rational Kummer image. There is no rational2-torsion
because the cubic has S3 Galois group. Hence \(\dim W=\operatorname{rank}E(\mathbf Q)\).
Equations(1)–(2) imply

\[
 \dim(W\cap U)\ge \operatorname{rank}E(\mathbf Q)-18.
\]

Joining the previously certified rank≥27 label now gives

\[
 \boxed{\dim(W\cap U)\ge9,\qquad c_S\ge9.}
\tag{3}
\]

Since G∩U=0, these nine directions inject into the additional quotient.
They determine at least nine independent unramified, S-split quadratic
characters of K; their compositum over K has degree at least **512**.
This is a finite field encoding their class information. It is **not** a
minimal geometric carrier whose rational points explain their simultaneous
solubility. Such a carrier remains unidentified.

This separates two claims which would otherwise look similar:

* **Proved retrospective necessity:** an independently specified arithmetic
  space U contains at least nine rational directions on this successful
  fibre. A purely local-boundary explanation cannot supply all ten.
* **Still missing for prediction:** construct enough of U independently,
  compare it with the control, and explain why a large subspace is rational
  rather than Sha. Equation(3) does not supply those representatives or a
  test of rationality before the rank is known.

Thus the small generic strict kernel in the earlier panel does not mean
that the high fibre lacks a strict block. Its successful strict block must
lie **outside the generic subgroup**. Computing the generic block alone
misses exactly the part needed to explain the jump.

## A sharper finite target than a full class-group calculation

For this high fibre the unconditional bound is

\[
 \operatorname{rank}E(\mathbf Q)\le c_S+18.
\]

A certified upper bound **c_S≤9** would therefore match the existing
rank≥27 certificate. It would force c_S=9, e=1, exact rank27,
and \(\Sha(E)[2]=0\). No such upper bound has been obtained.
A bound c_S≤8 would contradict the retained rank certificate; it would
signal an erroneous class computation or an inapplicable bound.

This gives a concrete next arithmetic target: certify an upper envelope
for the localized class group modulo2, rather than calculate every unit,
class invariant and class representative. For the control, its corresponding
rank envelope is c_S+19; observed-zero recovery supplies no upper bound.

A prospective condition of the same type is **incidence**:

\[
 \operatorname{rank}E_t\ge R\ \Longrightarrow\
 c_{S(t)}\ge R-(\ell(t)-1),
\]

whenever the same norm-square derivative and nonzero boundary-functional
hypotheses are verified. Its converse is false as a general inference:
large class incidence does not prove rational solubility. The necessary
S-class dimension is a meaningful target, but an unknown quantity is not
an operational score for Agent1.

## Bounded representation test and lesson

The previous class-group attempts exhausted a256MiB PARI stack. This turn
used one `polredbest` reduction per field with exact maps in both directions,
then transported the existing maximal-order basis into the new presentation.
Independent rational arithmetic checks the inverse maps, trace discriminants
and derivative norm identity. The field discriminants have383 and433 bits.

The high cubic's largest coefficient drops from308 to256 bits; the low's
from258 to253. Both reduced class-group attempts still exhaust the same
256MiB cap. There is no certified class-group result and no new class
representative. This bounded experiment rules out the hope that this single
reduction alone would make the existing backend sufficient; it does not
prove all class-group methods infeasible.

[PARI's primary documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#polredbest)
describes this as a smaller presentation of the same field, without a
minimality guarantee. The experiment uses exactly that role. It does not
change technical factor-base parameters or increase the memory limit.

The next useful work is a mathematically justified modulo2 S-class upper
certificate or a new independent class-generation method on this fixed
pair. Repeating the unchanged full class-group call, supplying the known
exceptional classes, or expanding a point-search budget would not close
the missing incidence-to-solubility implication.

Replay:

```sh
sage -python elliptic-curves/rank-jump/verify_matched103b2_class_boundary.py check
```

The two30-second boundary caps, two15-second reduction caps and two20-second
class caps are pinned in the protocol. The computations are sequential,
checkpointed and confined to rank-jump artifacts. No active search files,
worker settings, rank ledger or shared navigation are changed.
