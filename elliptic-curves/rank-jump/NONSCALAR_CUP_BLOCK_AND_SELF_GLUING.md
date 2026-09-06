# An explicit nonscalar norm-lifting block and its self-gluing explanation

The bounded nonscalar control closes a complete two-dimensional arithmetic
block. In the same cubic field, with the same strict classes and allowed
ramification, one multiplier obstructs every nonzero class and its negative
admits explicit unramified norm lifts of both generators.

The zero case has a further exact explanation: the associated genus-two
curve has two elliptic quotients isomorphic to the same curve. This is a
controlled self-gluing event, **not a new elliptic rank result**. It supplies
an independently verified nonscalar test for the unresolved production
cup calculations and prevents treating a zero cup matrix as evidence of
new Mordell–Weil directions.

## Fixed arithmetic inputs and the falsifiable test

Use the retained small control
\[
K=\mathbb Q(\theta),\quad
f(T)=T^3-11T^2-14T-1,\quad S=\{2,163,\infty\}.
\]
Its certified class group is `C2 x C2`. The previously certified strict
classes are
\[
\beta_0=\theta^2-10\theta+1,\qquad
\beta_1=\theta^2-13\theta+12.
\]
They are positive at all real places, square above `S`, and unramified
outside `S`. Both norms are `625`. Their independent strict characters
already fill the two-dimensional dual of the S-class group, so they form
the complete strict space `U` for this control. No elliptic point was used
to select these classes.

The [protocol](NONSCALAR_CUP_CONTROL_PROTOCOL.json) fixes only
\(\gamma_+=1+\theta\) and \(\gamma_-=-(1+\theta)\), with two norm equations
per multiplier and a 30-second cap per multiplier. Both multipliers are
units, so neither changes the allowed ramification set. Their real sign
patterns are `(-,+,+)` and `(+,-,-)`: neither is a rational scalar times
a square. The intervals `(-2,-1),(-1,0),(12,13)` certify the three roots
and these signs exactly.

For `F_gamma=K(sqrt(gamma))`, let
\[
B_\gamma(\beta,\psi)=\chi_\psi(\kappa_\gamma(\beta)),
\quad\kappa_\gamma(\beta)=\gamma\cup\beta\in C_S/2.
\]
The predeclared test was bilinearity against the independently known scalar
matrix: `B_gamma+ + B_gamma- = B_-1`. It did not assume either individual
matrix. Independent norm and parity-ideal calculations give
\[
\boxed{B_{\gamma_+}=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\qquad B_{\gamma_-}=0.}
\]

All four **global norm equations succeed**. The distinction lies in
ramification of the norm witnesses, not field norm solubility.

## Independently checked parity ideals

Write `p_(l,r)=(l,theta-r)` for a degree-one prime. A norm witness
`z=a+b*sqrt(gamma)` has a descended parity ideal `I`: outside `S`, its
valuation at a base prime is the common parity of the valuations of `z`
above that prime. The restricted cup is `[I]` modulo twice the S-class
group. This is the ideal formula in
[McCallum–Sharifi, Theorem 2.4](https://arxiv.org/pdf/math/0202161), applied
as in the [strict cup proof](CUP_IDEAL_AND_STRICT_LIFTING_OBSTRUCTION.md).

| Multiplier | Class | Descended parity ideal | Artin values against `(beta0,beta1)` |
|---|---|---|---|
| `gamma+` | `beta0` | `p_(5,2) p_(37,33)` | `(0,1)` |
| `gamma+` | `beta1` | `p_(5,1) p_(5,2) p_(37,33)` | `(1,0)` |
| `gamma-` | `beta0` | `p_(5,2) p_(13,2) = (theta-2)` | `(0,0)` |
| `gamma-` | `beta1` | `p_(5,1) p_(5,2) p_(5,3) = (5)` | `(0,0)` |

The raw PARI norm witnesses and valuations are retained in the
[control certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_nonscalar_cup_control_v1.json).
A [separate verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_nonscalar_cup_control_verification_v1.json)
uses rational polynomial arithmetic and Hensel lifting to precision 64,
with no number-field CAS, norm solver, class group, or ideal-valuation
routine. It checks each norm identity and every recorded valuation.
At split primes it evaluates both square-root branches; at inert primes
the unramified quadratic norm doubles the valuation.

Support is complete, not guessed: the coefficients are integral away from
their recorded denominator primes, `gamma` is a unit, and the norm target
has rational norm `5^4`. Hence no other prime can enter. The only primes
needed are 5, 13 and 37, and the cubic splits simply at each.

## Explicit simultaneous unramified lifts

For `gamma-`, divide the first norm witness by `eta0=theta-2` and the
second by `eta1=5`. The principal-ideal equalities in the table follow
from `N(theta-2)=65` and the three simple roots modulo 5. The resulting
elements `z'_i=a'_i+b'_i*sqrt(gamma-)` have
\[
N(z'_i)=\beta_i/\eta_i^2
\]
and even valuation at every prime outside `S`. Explicit coefficients are

\[
\begin{aligned}
a'_0&=(-23001-301957\theta+24991\theta^2)/169,\\
b'_0&=(80478+1053116\theta-87162\theta^2)/169,\\
a'_1&=\beta_1/25,\qquad b'_1=-\beta_1/25.
\end{aligned}
\]

Their [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_nonscalar_unramified_lifts_v1.json)
checks the identities and removal of every odd valuation. The norm images
are the original independent squareclasses, so the two lifts and their
products give the complete two-dimensional block.

More precisely, the permutation-module sequence for the quadratic extension
and Shapiro's lemma identify the norm-image obstruction with cup by `gamma`:
\[
\operatorname{im}\!\left(
H^1(G_{F_\gamma,S_F},\mu_2)\xrightarrow{N}
H^1(G_{K,S},\mu_2)\right)\cap U
=\ker(\kappa_\gamma|_U).
\]
For `gamma+` this intersection is zero; for `gamma-` it is all of `U`.
This is a statement about **unramified Kummer norm lifts**. It does not
assert the local Jacobian Kummer conditions or elliptic rational solubility.

## The entire unit orbit and the geometric reason

The cubic is cyclic, with exact automorphism
\[
\tau(\theta)=\theta^2-12\theta-2=-1/(\theta+1).
\]
For `u0=gamma-`, its conjugates are
\[
u_0=-1-\theta,\quad u_1=1+12\theta-\theta^2,
\quad u_2=\theta^2-11\theta-14,
\quad u_0u_1u_2=1.
\]
Their signatures span the even-parity plane. Together with `-1` they
span three independent unit squareclasses. A totally real cubic has unit
squareclass dimension three, so these give **all ordinary unit
squareclasses** of this field, not all S-unit squareclasses.

The strict space is Galois-stable because it is the complete S-class
character space and `S` is Galois-stable. Naturality of cup products
transports the zero map for `u0` to its conjugates. Bilinearity and the
nonzero scalar `-1` control therefore prove, for every ordinary unit `u`,
\[
\boxed{
B_u=\begin{cases}0&N_{K/\mathbb Q}(u)=+1,\\
\begin{pmatrix}0&1\\1&0\end{pmatrix}&N_{K/\mathbb Q}(u)=-1.
\end{cases}}
\]
The [orbit certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_nonscalar_cup_orbit_v1.json)
checks the automorphism, conjugate product, signature ranks and all eight
multiplier classes. Only four norm equations were computed; the remaining
statements follow from equivariance and bilinearity.

The two selected multipliers also illustrate a Galois change: the conjugate
squareclasses of `gamma+` have rank three and those of `gamma-` rank two.
Their normal closures over `Q` consequently have groups `C2 x A4` and
`A4`, of orders 24 and 12. Both individual quadratic extensions of `K`
still have degree six over `Q`. This group calculation alone is not the
solubility proof; the cup matrices and principal corrections supply it.

There is an exact geometric explanation available without exceptional
points. Put
\[
C_s:\ Z^2=f(sX^2-1),\qquad s\in\{1,-1\}.
\]
The two elliptic quotients are `E0: y^2=f(T)` and
`E_s: Y^2=s(T+1)f(T)`. For the latter, take
`v=1/(T+1)`, `w=Y/(T+1)^2`; then
\[
w^2=s(v^3+11v^2-14v+1).
\]
Its standard cubic model is
\[
E_s:\ y^2=x^3+11s x^2-14x+s.
\]
Thus **`E_- = E0` and `E_+ = E0^(-1)`**. The negative multiplier is
a self-gluing, with the nontrivial order-three relabelling of two-torsion
given by `tau`. The positive multiplier compares the curve with its scalar
twist. This explains why a nonscalar computation can reproduce the scalar
obstruction while the other sign vanishes. The quotient identities do not
produce extra rank on `E0`.

## Ranked conclusions and the next missing implication

1. **Solubility obstruction:** the labelled restricted cup map can change
   an entire two-dimensional norm-lifting block at once. Here both its
   vanishing and nonvanishing are independently certified, and the zero
   case has explicit lifts. Field norm success by itself is too weak.
2. **Incidence:** the field and complete strict class space stay fixed.
   The norm-sign criterion is proved for the ordinary units of this
   particular cyclic cubic. It is not a new Selmer-incidence predictor.
3. **Weak rank explanation:** the zero case is a self-gluing, and its
   two-congruent quotients supply the same elliptic curve. Neither the
   normal-closure reduction nor the two-dimensional lift block proves
   new Mordell–Weil directions. The production control cubics have Galois
   group `S3`; this cyclic order-three self-identification does not transfer.
4. **Missing computation:** the independent production `1+theta` norm
   witnesses, their complete parity supports, and the non-strict local
   corrections for the nine retained CT bits remain open. The present
   pipeline now has independently checked nonscalar positive and zero
   controls, including inert places and denominator-prime corrections.
5. **Next useful bridge:** identify the precise local-Jacobian obstruction
   after these explicit unramified norm lifts, before interpreting them
   as rational points. In parallel at the conceptual level, compare the
   production descent fields with this self-gluing criterion; do not assume
   an abstractly smaller Galois group preserves rational directions.

Agent 1 receives a validated arithmetic diagnostic and its scope, not a
candidate score. No production point, parameter, search policy, worker
limit, or status entry was changed. No visibility feature was measured.

Replay from the repository root:

```sh
python3 elliptic-curves/rank-jump/verify_nonscalar_cup_control.py --check
python3 elliptic-curves/rank-jump/nonscalar_unramified_lifts.py --check
python3 elliptic-curves/rank-jump/nonscalar_cup_orbit.py --check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_nonscalar_cup_control.py
```
