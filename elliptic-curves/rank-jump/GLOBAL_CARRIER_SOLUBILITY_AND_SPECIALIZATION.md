# Global carrier solubility does not select a prescribed specialization

**Follow-up:** the [labelled Cassels–Tate certificate](NATIVE_PAIR_CARRIER_HAS_A_SHA_OBSTRUCTION.md)
now proves the A,D carrier globally insoluble. Its UNKNOWN classification
below records the earlier incomplete computations, whose immutable outputs
are retained. The two positive carrier results and the prescribed-parameter
distinction remain valid.

The smallest cross-group pair from the completed bisection experiment has
an **everywhere locally and globally soluble genus-one carrier**, despite
having no simultaneous split in the 32-fibre sample. Its Jacobian has exact
rank **3** and \(\Sha[2]=0\). A pair that did split on the observed +8 fibre
also has a globally soluble carrier, with Jacobian rank **2** and
\(\Sha[2]=0\).

This closes a specific global-solubility question without searching new
original-family parameters. The missed pair is not explained by a Sha
obstruction on its fixed parameter carrier. Its rational points are infinite;
the remaining specialization condition is membership of the prescribed
parameter in their image.

A third, complementary carrier is everywhere locally soluble and has a
Jacobian with \(\Sha[2]\) of dimension two. **The class of that particular
carrier remains UNKNOWN.** Two bounded attempts to identify it failed their
completeness gates. No assertion that this carrier is Sha or has no rational
point is justified.

## Fixed comparisons and exact outcomes

Let
\[
C_{ij}:\quad u_i^2=f_i(t),\qquad u_j^2=f_j(t),
\]
using primitive quadratic forms from the existing generic R17 bisection
atlas, with only verified rational square factors removed. Distinct atlas
quadratics have disjoint geometric branch divisors, so the smooth projective
carrier has genus one. These are parameter carriers, not the elliptic
2-coverings of an already fixed fibre \(E_{t_0}\).

The [first protocol](GLOBAL_PAIR_SOLUBILITY_PROTOCOL.json) selects the
cross-group pair with smallest maximum coefficient height, then minimum
height and labels. It fixes a common cover and selects its smallest partner
from the observed +8 quartet for the positive control. This is retrospective
selection, not a validated prospective rule.

| Carrier labels | Observation in the frozen cohort | Exact Jacobian rank | Full 2-Selmer dimension | Exact Sha[2] dimension | Carrier C(Q) |
|---|---|---:|---:|---:|---|
| 1795d, 11278 | No simultaneous split in 32 fibres | 3 | 4 | 0 | Nonempty |
| 1795d, 0911e | Split together on 08234-009, observed +8 | 2 | 3 | 0 | Nonempty |
| 11278, 030cb | Complement A,D in the modulo-23 obstructed quartet | 2 | 5 | 2 | UNKNOWN |

All labels have the prefix `orbit-`. Every Jacobian in the table has one
rational 2-torsion direction. The third carrier was fixed by the
[separate complementary-pair protocol](DISJOINT_SOLUBLE_CARRIERS_PROTOCOL.json)
after the first two completed descents. Its complement B,C already splits
at the retained parameter \(-3115/2756\). The attempted hypothesis that A,D
would also be proved globally soluble was not established.

For the cross-group pair, the equations are explicitly
\[
\begin{aligned}
u^2&=409689-1439214t+328441t^2,\\
v^2&=-144492039-201200094t+18383017t^2.
\end{aligned}
\]
The exact equations, Jacobian models, local witnesses, and raw descent
outputs for all three carriers are retained in the certificates below.

## How global nonemptiness was proved

Parametrize the first conic from its already retained single-cover point
\((t_*,u_*)\). If \(f(t)=at^2+bt+c\), the line
\(u=u_*+z(t-t_*)\) gives
\[
D=z^2-a,\quad
N=t_*z^2-2u_*z+at_*+b,\quad
U=-u_*z^2+f'(t_*)z-au_*.
\]
Then \(t=N/D\), \(u=U/D\), and
\(U^2=D^2f(N/D)\). The second equation becomes a squarefree quartic
double cover of the z-line. All identities are checked exactly. This
rational degree-two divisor shows that the genus-one torsor class has
period dividing two.

The Jacobian is computed from the classical binary-quartic invariants I,J
as \(y^2=x^3-27Ix-27J\), cross-checked against Sage's Jacobian constructor,
and transported to a minimal model.

For local completion, the only potentially bad primes divide
\[
2\operatorname{disc}(f_i)\operatorname{disc}(f_j)
\operatorname{Res}(f_i,f_j).
\]
The factorizations and primality are proved, and exact square-value witnesses
are checked at every prime in their supports: 19, 18, and 21 primes for the
three rows respectively. This includes primes beyond the earlier bound 43.
A real witness is also retained. At every other prime, the smooth projective
genus-one reduction has a point by the finite-field bound, and Hensel lifts
it. Thus these are complete everywhere-local proofs, not small-prime screens.

PARI's effort-zero descents give matching unconditional rank bounds. Their
full 2-Selmer dimensions and CT quantity determine the displayed Sha[2]
dimensions. The reported upper bound is
\(C-T-s\), where C is the 2-Selmer dimension, T the rational 2-torsion
dimension, and \(s=\dim\Sha[2]/2\Sha[4]\); matching bounds identify
\(\dim\Sha[2]=s\). This uses the
[documented unconditional outputs of PARI ellrank](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellrank),
not analytic ranks or BSD.

For each of the first two rows, the returned points plus rational torsion
have independent finite Kummer images of dimensions four and three,
respectively. This separately certifies the free-rank lower bounds three
and two. Since \(\Sha(J)[2]=0\), their everywhere locally soluble
period-dividing-two torsor classes vanish. Hence each carrier has a rational
point. Its positive Jacobian rank then makes its rational points, and their
images under the finite parameter map, infinite.

The established [singleton character theorem (F4)](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md)
supplies generic rank at least 19 over the corresponding carrier function
field. This is a base-change
statement with the original marked rank-17 subgroup retained. It is not an
exact fibre rank, nor does it supply a lift over any particular t.

## Historical incomplete attempts on the third carrier

The A,D Jacobian's initial exact-rank certificate reports rank two and
Sha[2] dimension two, but returns only **one** free generator. A separate,
seeded effort-one call, fixed by the
[class-identification protocol](CARRIER_SHA_CLASS_PROTOCOL.json), also returns
one. Therefore the rational Kummer fingerprints span dimension two including
torsion; the full rational Kummer space has dimension three.

An explicit degree-four map from the quartic to its Jacobian was nevertheless
constructed and verified:
\[
X=N(z)/f(z),\qquad Y=G(z)/w^3,\qquad w^2=f(z).
\]
Its syzygy is checked as a polynomial identity. In the cubic etale algebra
\(K=\mathbf Q[\theta]/(\theta^3+A\theta+B)\), the audit also verifies
\[
N(z)-\theta f(z)=\beta h(z)^2.
\]
This gives the carrier's descent representative. Adding beta raises the
finite-character rank from two to three. **That does not prove a nonzero
Sha class**, because the missing rational direction could account for it.

The [isogeny continuation](CARRIER_ISOGENY_OBSTRUCTION_PROTOCOL.json) tests
the scalar component without increasing the earlier effort. For
\(J:y^2=x(x^2+ax+b)\) and its rational 2-isogenous curve J', the exact index
formula is
\[
\dim\alpha(J(\mathbf Q))+\dim\alpha(J'(\mathbf Q))
=\operatorname{rank}J(\mathbf Q)+2.
\]
The existing point and torsion give two scalar directions on J; beta gives
a third. Two independent scalar directions on J' would have closed the
original scalar quotient and proved an obstruction. The single effort-zero
call on J' returns no free points, however, leaving only its torsion scalar
direction. This gate also remains unclosed.

Both failed attempts are retained as reproducible UNKNOWN results. Finding
the second rational Kummer direction, or computing a class-specific
Cassels–Tate pairing, remains necessary to classify A,D. A Jacobian's
nontrivial Sha group does not by itself identify a given carrier as Sha.

## The specialization distinction

These computations separate three different statements:

1. \(C_{ij}(\mathbf Q_v)\ne\varnothing\) for each place v, possibly using a
   different local parameter at every place.
2. \(C_{ij}(\mathbf Q)\ne\varnothing\), so some rational parameter admits
   both lifts.
3. A **prescribed** rational parameter \(t_0\) lies in
   \(\pi(C_{ij}(\mathbf Q))\).

The first two rows prove statement 2, and even infinitely many rational
parameter lifts. They do not decide statement 3. This is why global
nonemptiness or the auxiliary Jacobian rank cannot explain the difference
between those specific observed fibres.

For the finite fibre of these particular parameter covers, away from branch
values, statement 3 is exactly that every \(f_i(t_0)\) is a rational square.
There is **no additional Sha gap for this finite square-root system**:
square roots at the same fixed t0 over every completion imply rational square
roots. Indeed every finite prime valuation is even and the real sign is
positive. This is different from a genus-one elliptic 2-cover attached to a
Selmer class of \(E_{t_0}\), where everywhere-local solubility can still
represent Sha.

The next global computation must therefore specify which problem it solves:
membership of t0 in a common carrier image, or triviality of a particular
Selmer torsor at t0. Replacing either by “the auxiliary curve has rational
points” loses the specialization condition.

## Mechanisms and the remaining gap

- **Solubility, proved partial construction:** simultaneous square conditions
  produce the earlier independent subblocks. The two globally soluble pair
  carriers now also have certified positive-rank Jacobians. What is missing
  for an extreme fibre is a condition forcing its prescribed parameter into
  enough compatible images, with surviving quotient independence.
- **Solubility, unresolved class-specific obstruction:** A,D has complete
  local proofs and a Jacobian with Sha[2] dimension two. Its actual torsor
  class must still be identified. This is a concrete target for a global
  obstruction computation, not an established negative control.
- **Weak explanation:** rational points somewhere on the carrier, or a high
  auxiliary Jacobian rank. The unobserved pair has rank three while the
  observed pair has rank two; this comparison supports no rank predictor.
- **Incidence:** additional product-character sections remain uncomputed.
  The known two generic singleton directions do not prove such a section.
- **Visibility:** no candidate score, chart exposure, or point-search
  performance conclusion follows from these descents.

No active selection policy or search output was changed. All new point
computations were bounded descents on fixed auxiliary Jacobians; no new
original-family parameter or high-rank candidate was generated.

## Certificates and replay

```sh
sage -python elliptic-curves/rank-jump/verify_global_carriers.py check
```

The replay checks the conic parametrizations, quartic identities, Jacobian
models, full bad supports, primality, every local witness, the two independent
finite-Kummer lower bounds, and fresh effort-zero rank/CT outputs. It also
replays both incomplete class-identification gates without promoting their
outcomes. Software versions and inputs are pinned in the artifacts.

- [Two-carrier global-solubility certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_global_pair_solubility_v1.json)
- [Complementary A,D descent](../../artifacts/generated-results/elliptic-curves/rank_jump_disjoint_soluble_carriers_v1.json)
- [Explicit cover class and incomplete rational Kummer basis](../../artifacts/generated-results/elliptic-curves/rank_jump_carrier_sha_class_v1.json)
- [Incomplete isogeny obstruction test](../../artifacts/generated-results/elliptic-curves/rank_jump_carrier_isogeny_obstruction_v1.json)
- [Verification certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_global_carrier_verification_v1.json)
