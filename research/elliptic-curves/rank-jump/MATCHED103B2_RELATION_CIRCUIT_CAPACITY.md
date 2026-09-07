# A relation collector must expose non-unit parity circuits

The bounded near-root experiment produces one new principal relation on
the fixed high-gain 103b2 fibre and none on its matched control. Together
with the previously certified box relation, the high fibre has two
relations whose outside-S parity vectors are independent. Their generated
strict image is therefore **zero**, including after adjoining the generic
point subgroup. No additional CT entry becomes computable from this pilot.

The useful result is a precise gate for future relation collectors:
counting relations, smooth norms, or parity dependencies alone does not
measure the independent strict block. The difference between rational and
mod-two ranks of the **complete integer valuation matrix** bounds the
non-unit contribution obtainable from those generators.

This is a follow-up to the
[16-fibre governing/CT comparison](FRESH_RANK27_GOVERNING_AND_CT_COMPARISON.md)
and the [independent norm-relation gate](MATCHED103B2_INDEPENDENT_NORM_RELATION_GATE.md).
The panel's additional-class CT matrices remain UNKNOWN. This small test
does not replace that panel with a two-fibre rank-discrimination claim.

## Frozen experiment and exact outcome

The [protocol](MATCHED103B2_ROOT_CIRCUIT_PROTOCOL.json) uses the already
certified maximal-order binary cubic F for each masked field. For each of
its three ordered real roots alpha, it enumerates primitive distinct pairs

\[
 1\le n\le512,\qquad m=\lfloor n\alpha\rfloor+s,
 \qquad s\in\{-1,0,1\}.
\]

The smoothness bound is 400000; retained residual cofactors have at most
32 bits. The first 32 eligible relations in the declared order may be
audited. The single prior box relation is also retained. There is no
parameter expansion, elliptic point search, exceptional class input, or
adaptive enlargement. Each stage has a fixed time and memory cap.

| Quantity | High: 3726/881 | Control: −1049/2296 |
|---|---:|---:|
| Family / marked generic rank | 103b2 / 17 | 103b2 / 17 |
| Recorded rank lower bound / gain | 27 / +10 | 17 / 0 |
| Primitive near-root norm values | 2807 | 2802 |
| Eligible near-root relations | 1 | 0 |
| Audited relations including prior seed | 2 | 0 |
| Rank of complete valuation matrix over Q | 2 | 0 |
| Rank over F2 | 2 | 0 |
| Outside-S parity rank | 2 | 0 |
| Strict image of these generators together with G | 0 | 0 |

The high fibre's new pair is (m,n) = (−3164501,20), with

\[
 F(m,n)=-4484728154678871050167797299866580,
 \qquad\text{residual}=1337134327.
\]

The complete principal ideal includes the fixed a² factor in
Norm(au+vw)=a²F(m,n). Its odd outside-S support includes a prime ideal over
53, absent from the old relation; the old relation includes one over 491,
absent from the new. Thus their outside-S parity vectors are independent.
Every generic Kummer class has even valuations outside S. A product of
these two elements and generic classes can be strict only when both new
exponents are even. The previously certified G∩U=0 then proves the zero
strict image in the table. This conclusion does not require interpreting
zero search output as a class-group bound.

The control's smallest residual has 64 bits and fails the frozen threshold.
The high fibre's accepted residual has 31 bits. This difference measures
the cost of obtaining relations in these coordinates. It is neither a
Selmer-incidence statistic nor a rational-solubility statistic, and it
does not measure visibility in an elliptic point search. No candidate
score should be inferred from it.

## A capacity bound for a supplied set of principal generators

Let K be a totally real cubic field and alpha_1,...,alpha_m nonzero
elements. Let A have one row for every finite prime ideal appearing in
their divisors and entry A_(p,i)=v_p(alpha_i). This must be the **complete
integer** matrix, including bad primes and even valuations. Put

\[
 r_{\mathbf Q}=\operatorname{rank}_{\mathbf Q}A,\qquad
 r_2=\operatorname{rank}_{\mathbf F_2}(A\bmod2),\qquad
 d=r_{\mathbf Q}-r_2.
\]

Consider squareclasses of products of these alpha_i which are positive
at every real place and have even valuation at every finite place.
Their dimension is at most

\[
 \boxed{2+d.}
\]

In particular, this bounds the strict locally square subspace generated
by the supplied elements. It is an upper bound for this construction,
**not** for all of U, the Selmer group, or the ideal class group.

Proof. The integer kernel L=ker(A:Z^m→Z^n) is saturated: if kx is in L,
then kAx=0 in the torsion-free group Z^n, hence Ax=0. Thus L/2L embeds
in F2^m with dimension m−r_Q, and its image lies in ker(A mod 2), of
dimension m−r_2. The quotient has dimension d. An integer relation in L
produces an ordinary unit. If its mod-two vector represents a positive
squareclass, that unit is totally positive, since changing exponents by
even integers multiplies by a square. In a totally real cubic field,
units modulo squares have dimension three; the sign map is nonzero
because of −1, so its kernel has dimension at most two. The contribution
from L is therefore at most two, while the remaining quotient contributes
at most d. Additional local squareness conditions can only lower this
bound.

The same proof gives the sharper bound u_+ + d when the dimension u_+ of
totally positive units modulo squares is known. It also explains why the
bound need not be attained: further principal identities, squares, or
local conditions can kill candidate products. A character-rank or exact
squareclass certificate is still needed for a lower bound.

For the present high-fibre matrix, d=0, and in fact the parity kernel
itself is zero, giving the sharper dimension zero. The earlier
[nine-direction necessity](MATCHED103B2_JUMP_REQUIRES_NINE_STRICT_DIRECTIONS.md)
says that this fibre contains at least nine strict rational directions
outside G. A principal-product construction capturing nine independent
strict classes must therefore reach d≥7. This is a necessary condition
for that construction to succeed, derived using the retained rank lower
bound; it is not a condition on t or an independent rank prediction.

## What this changes in the research plan

1. **Incidence gate:** obtain an independently certified strict class
   basis or a class-group presentation with a justified generation bound.
   For relation products, retain integer valuations and check d, local
   squareness, and actual squareclass independence. A large parity kernel
   with d=0 cannot supply the required nine-dimensional block.
2. **Solubility gate:** evaluate CT on the resulting additional classes.
   A nonzero pairing obstructs simultaneous rationality; a vanishing
   pairing still leaves higher descent and rational solubility unresolved.
   Neither this pilot nor the generic governing fields supplies that gate.
3. **Weak explanations:** generic governing-field degree already failed
   the 16-fibre comparison. Easier norm smoothness fails to reach even the
   strict-class gate here. Increasing the same norm budget without a
   credible way to create non-unit circuits has little discriminating value.
4. **Next useful comparison:** compare a certified additional strict
   dimension and its CT obstruction across the frozen panel once those
   classes are available independently. Until then, report missing values
   rather than use exceptional-point representatives to fill them.

No new prospective selection feature is established. The missing chain is
still: equation-defined independent additional classes → CT information
on those classes → a sufficient rational-solubility criterion. The current
artifacts resolve neither the second arrow nor a condition on t.

## Certificates and replay

The [capture](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_root_circuits_v1.json)
binds the protocol, source forms, prior seed and computation scripts.
The [verification](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_root_circuits_verification_v1.json)
independently checks all 5609 norm values, 3072 root floors using disjoint
rational sign-changing intervals, both full principal ideals, and the
outside-S parity argument. It does not claim an independent replay of
the unused local-signature entries. All checks pass.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_matched103b2_root_circuits.py check
```

Only rank-jump-specific scripts, notes and immutable certificates change.
Active searches, worker policies and mathematical status entries are untouched.
