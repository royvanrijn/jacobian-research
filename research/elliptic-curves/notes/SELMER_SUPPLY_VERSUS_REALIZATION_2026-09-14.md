# Selmer supply versus visible Mordell–Weil gain

Interpretation of retained evidence, not a completed global descent. The
mathematical authority remains `EC-DET1092-RR-FULL-SELMER-COMPARISON-20260908`
and its [canonical comparison](DET1092_RR_FULL_SELMER_COMPARISON_2026-09-08.md).
The fixed cases09/08 have **unknown** global Selmer dimensions.

## The simultaneous phenomenon already has a witness

Write a retained prime block as `[A_p | b_p]`, where the columns of `A_p`
evaluate inherited divisors and `b_p` evaluates the marked rational class.
Each equation `A_p c_p=b_p` is soluble. Stacking the blocks gives `A c=b`,
which is insoluble: the inherited rank is16 and the augmented rank17.
The saved row functional `lambda` satisfies `lambda A=0`, `lambda b=1`.
Thus one coefficient vector cannot work at all the retained places.

The [binary audit](../../artifacts/generated-results/elliptic-curves/rr_simultaneous_compatibility_audit_v1.json)
checks these statements directly on both saved matrices:

| RR member | Rows / prime blocks | Inherited / augmented rank | Every individual block soluble |
|---|---:|---:|---|
| Historical first unlock,09 |67 /40|16 /17|Yes|
| Source-only section0 control,08 |58 /38|16 /17|Yes|

This audit checks binary algebra only, relying on the original
[divisor replay](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md) and
[control replay](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md) for arithmetic
meaning. It does not reconstruct divisors or claim new local completeness.
The separate theta argument supplies the nonzero true-Kummer kernel class
`D0`, giving true inherited and marked dimensions17 and18.

This is a failure of membership in the diagonal localization of the inherited
subspace, not failure of local solubility of the marked covering. The marked
class has an explicit global rational witness, hence zero Sha image.
The row functional is not a nonzero Cassels–Tate row.

## What the pair actually controls

Both members lie over the same elliptic fibre Curve302. The historical marking
is an extra elliptic direction; the control marking is a generic section.
Both Jacobians nevertheless have a certified18-dimensional rational Kummer
subspace. Neither full Jacobian rank is known, so these are not certified
different-rank Jacobians.

The control phenomenon has a geometric explanation already in the canonical
specificity proof: marking one point of a degree7 intersection differs from
restricting its entire divisor. Over the incidence cover `u=r0(t)`,
`Tr(xi)=Phi(S0)` and `7 xi-Phi(S0)` has trace zero. The Jacobian can acquire a
direction relative to the Picard image even when the marked elliptic point
is inherited. Global independence and simultaneous local incompatibility
therefore do not by themselves distinguish exceptional elliptic gain.

## Correct interpretation of the dimension comparison

For an abelian variety `A/Q`, the Kummer exact sequence gives

\[
0\longrightarrow A(\mathbf Q)/2A(\mathbf Q)
\longrightarrow\mathrm{Sel}_2(A)
\longrightarrow\Sha(A)[2]\longrightarrow0,
\qquad s=r+t+h,
\]

where `t=dim A(Q)[2]`, `h=dim Sha(A)[2]`. See
[Poonen–Schaefer, §13](https://math.mit.edu/~poonen/papers/descent.pdf).
Here `t=0` in both cases. Equal Selmer dimensions **together with a certified
rank difference** would force an opposite Sha[2] dimension difference.
Equality alone proves neither a rank difference nor a realization mechanism.

For this pair let `G` be the inherited true image (dimension17), `V` the
known rational image (dimension18), and `W` the full rational image. Then

\[
0\longrightarrow W/V\longrightarrow\mathrm{Sel}_2(J)/V
\longrightarrow\Sha(J)[2]\longrightarrow0,
\quad s-18=(r-18)+h.
\]

| Certified future outcome | Supported conclusion |
|---|---|
| `s09 > s08` | More global Selmer capacity on this historical RR member; no causal elliptic-gain conclusion. |
| `s09 = s08 = 18` | Both Jacobians have exact rank18 and zero Sha[2]; the abstract rational Selmer flags have the same dimensions. |
| `s09 = s08 > 18` | Equal total capacity; allocation between further rational classes and Sha remains unknown. |
| `s09 < s08` | The control has more capacity; the proposed positive supply discriminator fails on this pair. |

Within a single vector space, all flags of fixed dimensions are equivalent
under linear automorphisms. Thus “position of the Mordell–Weil subspace”
needs extra arithmetic structure: labelled localization maps, a specified
correspondence, cover equations, or a pairing. There is no default
identification of these two Jacobians' Selmer spaces.

Rational classes pair trivially against the entire Selmer group under the
pairing pulled back from Sha. Pairing only `V` yields no discriminator.
Even the full radical need not equal `W`: divisible Sha classes can survive
there. The [residual pairing theorem](../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md)
records the precise `2 Sha[4]` boundary. This is why complementary classes,
their pairings and rational witnesses matter if `s>18`.

## Next gate and hypothesis boundary

The existing two-case work item remains the correct global task: certify a
complete supported squareclass space and all required simultaneous local
conditions, retaining the true/fake distinction. The two previous bounded
order passes are unresolved; repeating those calls or the closed local
panels supplies no new upper bound. No new class-group campaign is launched
by this interpretation or binary audit.

The [4,482-record analysis](SEED_AND_AMPLIFICATION_HISTORY_2026-09-13.md)
supports failure of the tested amplification predictors to transfer. Its
records are not independent curves, and search exposure affects observed
gain. That failure does not establish that elementary features cannot work
or that global cohomology must explain amplification.

The sharpened hypothesis is that **arithmetic structure beyond marginal
local dimensions and known rational-image dimensions may discriminate
marked constructions or rational realization**. The RR pair can test the
first question. Explaining elliptic amplification additionally requires
an arithmetic bridge to the elliptic quotient direction and validation
across distinct fibres; neither follows from this pair's Selmer dimensions.

Reproduce the inexpensive binary audit from the repository root:

```sh
python3 research/elliptic-curves/cas/audit_rr_simultaneous_compatibility.py
```

## Further global-order checkpoint

The subsequent bounded order join improves the global prerequisite without
repeating either prior maximal-order call. In the common reduced field,
transport the verified binary order and retain the partial order. Their
36 pairwise basis products span an integral order: products of sums remain
in the span because each input is a commutative ring containing1. An integer
Hermite basis realizes this lattice without factoring any discriminant.

| Case | Smallest previous order discriminant, bits | Joined order, bits |
|---|---:|---:|
| Historical09 |26213|22655|
| Control08 |37559|37378|

The [retained outputs and Newton-trace replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_order_join_v1/replay.json)
certify ring closure, containment of both input orders, exact indices and
trace discriminants. The independent checker does not reconstruct the
Hermite basis. The bit sizes describe orders, not certified field
discriminants or intrinsic arithmetic complexity. Maximality, supported
squareclass completeness and both Selmer dimensions remain unknown.

Each case used a60-second alarm and an8GiB address-space cap; arithmetic
finished in less than one second per case. The command-interface failure
and initial2GiB address-space mapping failure are retained in
[execution.json](../../artifacts/generated-results/elliptic-curves/det1092_rr_order_join_v1/execution.json).
These are implementation failures, not arithmetic obstructions.

The bundled runtime replays this certificate with:

```sh
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python research/elliptic-curves/cas/verify_rr_pair_order_join.sage
```

Use these bases as improved order checkpoints. Supplying them to `nfinit`
still does not certify maximality. No new class-group run is justified on
that basis alone.

### Coprime-support successor

The [third-order containment audit](../../artifacts/generated-results/elliptic-curves/det1092_rr_order_join_v1/third-order-containment.json)
transports the remaining v2 order into the joined field. It is contained in
the join in both cases; all three previously saved orders are therefore
already represented. This closes that particular source of missing integrality.

The next attempt uses the distinct [PARI coprime-divisor input contract](https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html#nfbasis).
Gcd splitting of the polynomial discriminant, joined-order discriminant,
indices and prime divisors through1000 gives pairwise-coprime divisors with
an exact product decomposition. None of the composite divisors is assumed
prime or squarefree. One `nfinit([f,divisors])` and one `nfcertify` call per
case use a60-second whole-attempt alarm,8GiB address-space cap and512MiB PARI
stack. No supplied basis, full factorization, BNF or Selmer routine is used.

| Case | Returned order discriminant, bits | Unresolved certification cofactor, bits |
|---|---:|---:|
| Historical09 |22655|22633|
| Control08 |28881|28878|

Both attempts completed in under three seconds. Both certification calls
returned a nonempty unresolved cofactor. The
[independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_coprime_support_v1/replay.json)
checks the exact support decomposition, ring closure, Newton trace
discriminants and containment of the joined orders. It does not reproduce
the certification calls or prove the remainders squarefree. These are
successor **integral-order** checkpoints, not maximal orders. The global
Selmer dimensions and their comparison remain UNKNOWN.

```sh
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python research/elliptic-curves/cas/audit_rr_third_order.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python research/elliptic-curves/cas/verify_rr_coprime_support.sage
```

### Stopping boundary for this completion attempt

The [exact gcd check](../../artifacts/generated-results/elliptic-curves/det1092_rr_coprime_support_v1/gcd-boundary.json)
shows that neither unresolved cofactor splits against the new order index,
the complementary part of its order discriminant, or the retained support
divisors. The control's additional index has4249bits; the historical index
is1. Thus feeding back that index alone supplies no new support refinement.
This does not prove either cofactor squarefree or any general infeasibility.

The completion attempt stops blocked on certified global arithmetic. The
available order-joining and gcd-support refinements have been used, but
maximal orders, complete class/unit contributions, all required local
conditions and full true Selmer groups are still absent. Resume requires
an additional certified method or certificate addressing that global gap;
repeating the same certification inputs is not a new experiment. Neither
the equality branch nor the inequality branch has been established.
