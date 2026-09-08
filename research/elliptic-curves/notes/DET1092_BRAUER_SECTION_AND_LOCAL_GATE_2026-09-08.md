# Surface Brauer classes: a rational-section tree and a local blindness theorem

## Exact conclusion

**New verified deduction.** Let `X/Q` be the completed determinant1092 K3
parent, `O` its zero section, and `alpha` any class in `Br(X)`. Subtract
its constant restriction to `O`, so `alpha|O=0`. Then its restriction to
**every generic Mordell--Weil section** is zero.

For every such `alpha` of2-power order, at all189 pairs of the unchanged
nine-address,21-place panel,

\[
\boxed{\alpha(P_v)=0\quad\text{for every }P_v\in E_t(\mathbf Q_v).}
\]

This concerns all local points, not only the known generic points or the
historical first seed. For302 the checked places include infinity,2 and
**every prime dividing its exact discriminant**.

**New reciprocity boundary.** If a normalized2-primary surface Brauer class
has a nonzero value on a rational302 point, its local evaluation must be
nonzero at at least two primes outside that footprint. Both are good primes
for the original elliptic curve. No such class has been constructed.

**Established cohomological deduction applied to the certified parent.**
The full geometric Picard group is torsion-free of rank19 with trivial
Galois action. Thus `Br1(X)/Br(Q)=0`. A useful nonconstant surface Brauer
class would have to be transcendental; algebraic Brauer classes cannot
distinguish the first seed.

These are explicit obstructions to one global-class route. They do not
exclude the original elliptic Kummer class, all Selmer-based mechanisms,
transcendental Brauer detection at other primes, or nonlinear incidence.

## A connected tree of rational sections

**Verified generic geometry.** On the proved24-I1 K3, write `G` for the
full generic height Gram. In the integral basis `O,F,phi(P_i)`, use

\[
O^2=-2,\quad O.F=1,\quad F^2=0,\quad
\phi(P_i).\phi(P_j)=-G_{ij}.
\]

The generic section `S_i` has class `O+(G_ii/2)F+phi(P_i)`. Hence for
distinct sections, including `O`,

\[
S_i.S_j=\frac{\widehat h(P_i-P_j)}2-2.
\]

**New verified construction.** Among all153 pairs of the original18
curves, the intersection-degree-one graph is connected. One deterministic
spanning tree is the following; indices are zero-based.

```text
O
├── S13
│   ├── S1 ── S15
│   ├── S7
│   ├── S8
│   ├── S9
│   ├── S10
│   └── S16
└── S14
    ├── S0
    ├── S2
    ├── S3
    ├── S4
    ├── S5
    ├── S6
    ├── S11
    └── S12
```

Its17 difference words have height6 and integer determinant1, so they
also provide an integral norm-six basis of the existing MW17 lattice.
No new section or rank direction is created.

Every edge meets at an **explicit rational point**, not merely a geometric
intersection. For example the poles joining the zero section are

\[
t(O\cap S_{13})=
-\frac{6481138169118894561298909600}{46349538549253852152496289607},
\]
\[
t(O\cap S_{14})=
-\frac{1656735636583736590771280}{52970289929165564535992671}.
\]

The [rational-intersection certificate](../../artifacts/generated-results/elliptic-curves/det1092_brauer_section_gate_v3/rational-intersections.json)
contains all17 parameter values, charts and coordinates. The independent
checker uses the19-dimensional intersection form, then directly computes
common zeros of the saved section-coordinate differences. It verifies
each affine elliptic equation and nonsingular fibre point; at the two
zero-section meetings it checks the exact pole orders2 and3. The17 derived
locations are geometry witnesses, not new control or point-search fibres.

## Why every generic section has constant normalized value

**Established literature.** `Br(P1_Q)=Br(Q)`, and Brauer evaluation on
zero-cycles factors through rational equivalence, functorially under proper
maps; see [Colliot-Thélène--Skorobogatov, Theorem5.1.3 and Proposition5.3.2](https://www.imo.universite-paris-saclay.fr/~jean-louis.colliot-thelene/BGgroup_book.pdf).

**New application.** Each section is a rational curve, so `alpha|S_i` is
constant. At a rational intersection, the constants on the two curves are
equal by evaluation at that common point. The connected tree therefore
sets all17 constants equal to `alpha|O=0`.

On the smooth generic elliptic fibre over `Q(t)`, the normalized map

\[
P\longmapsto\alpha(P)-\alpha(O)
\]

is a group homomorphism: `P -> [P]-[O]` is the identification with the
degree-zero Picard group, followed by the Brauer pairing. It kills the17
generators and hence every generic section. That section's restriction is
again a constant on `P1`, so vanishing at the generic point implies
vanishing at every specialization. This step uses no2-saturation
assumption and works for Brauer classes of **any order**.

**Established literature applied exactly.** The Hochschild--Serre sequence
injects `Br1(X)/Br(Q)` into `H1(Q,Pic(Xbar))`; see the same book,
section4.3, equation(4.9). The completed full-rational Picard19 theorem
makes this latter group zero: a continuous homomorphism from a profinite
group to the discrete torsion-free group `Z^19` has finite, hence zero,
image. This proves the algebraic Brauer assertion separately from the
section-tree argument.

## Entire local fibres, not a finite sample of points

**Previously established application, independently replayed here.** At
each of the same nine addresses and21 places, the specialized generic
sections surject onto `E_t(Q_v)/2E_t(Q_v)`. This includes wild2-adic and
real places, not just square tests at good odd primes. The
[original finite-place theorem](DET1092_SEED_KUMMER_COVER_2026-09-08.md#equation-only-finite-place-gate)
retains the exact local algebra and completeness proof.

**New deduction.** Let `H_v` be the subgroup generated by those local
generic points. From `E_t(Q_v)=H_v+2E_t(Q_v)`, induction gives

\[
E_t(\mathbf Q_v)=H_v+2^nE_t(\mathbf Q_v)\quad(n\ge1).
\]

For `alpha` of order dividing `2^n`, evaluation is a homomorphism,
vanishes on `H_v` by the section tree, and vanishes on `2^nE_t(Q_v)`.
It is therefore identically zero on the **whole** local fibre. The argument
extends the proved mod2 local completeness to all2-primary Brauer classes
without computing higher local descents.

The independent replay rebuilds the specialized equations and generic
points, checks the p-maximal cubic-order multiplication tables, and
verifies every generic relation and local-basis nonsquare using the
previous independent finite-ring helper. Its real check uses exact Sturm
signs. It never calls that helper's historical-point evaluator.

## What the first five-place character is not

**New distinction.** The earlier character

\[
\chi=\epsilon_7+\epsilon_{19}+\epsilon_{23}
       +\epsilon_{29}+\epsilon_{167}
\]

annihilates the **joint** generic code and takes value1 at the first seed.
Its separate local summands do not annihilate the generic local groups;
each such group already fills its local quotient. A normalized surface
Brauer class, in contrast, must annihilate them **at every place
separately**. Therefore this compatibility character is not the proposed
unramified surface Brauer evaluation in disguise. The elliptic Kummer
recognition certificate remains valid; no Brauer interpretation of it is
inferred.

**Established reciprocity, new boundary.** Evaluation at a rational point
is a class in `Br(Q)` whose local invariants sum to zero. The final checker
reconstructs the literal302 discriminant from its cubic and verifies the
complete saved prime factorization by multiplication and primality proofs.
All those primes,2 and infinity lie in the zero-evaluation footprint. A
nonzero global value must consequently have at least two nonzero local
invariants at other finite primes, both good for `E302`. This does not
identify the bad primes of an arithmetic model of the **surface**, or
compute any transcendental Brauer class.

## Reproducibility and retained failures

**Verified computation.** An initial4096-row prefix lookup, and then a
norm-six filter of the old bisection table, found no norm-six entries.
Both results are retained in versions1 and2. The table contains the
degree-two bisection classes, not the short section classes; these were
input-exposure failures, not nonexistence proofs. Version3 instead uses
the exact153 pairs of the original sections, with its rule frozen before
outcomes. No census was recomputed.

```sh
sage -python research/elliptic-curves/cas/construct_det1092_brauer_section_gate_v3.sage
sage -python research/elliptic-curves/cas/verify_det1092_brauer_gate_panel.sage
```

The [final replay](../../artifacts/generated-results/elliptic-curves/det1092_brauer_section_gate_v3/final-replay.json)
binds both exact checks, all inputs and source hashes. Each stage has a
25-second cap; the final independent replay finishes in about three seconds.
No exceptional point, later V3 artifact, new test fibre, point search,
Brauer/Selmer group, class group, paid backend, pilot change or detached
process was used.

**Unresolved.** The sought prospective seed and global arithmetic
characterization are not obtained. This theorem restricts a concrete
alternative to the failed Jacobian transfers; it is not a claim that every
possible global arithmetic discriminator has failed.
