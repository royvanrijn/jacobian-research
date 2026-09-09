# Switching fibrations: an exact global degree barrier

## Result and scope

**New deduction, independently verified.** Switching between the two known
elliptic fibrations does escape the previously preserved carrier parameter,
but the first such move from the cheapest constructed rational curve has
a large, exact projection-degree cost.

Use zero-based generic section indices. Let `T` be translation by the
inherited point `B` in the alternate fibration, and let `t_q` be translation
by any original generic section `q in M17`. The preceding construction
selected the rational curve

\[
 C=T(S_1),\qquad \deg(\pi|_C)=13.
\]

It is the unique smallest-degree member of the six previously frozen
generic constructions. For the two families

\[
 C_{\epsilon,q}=T^\epsilon t_q(C),\qquad \epsilon\in\{-1,1\},
\]

the **global**, not box-bounded, minima are:

| Final alternate move | Degree when `q=0` | Minimum over every `q != 0` | Number of minimizers |
|---|---:|---:|---:|
| `T^-1` | 1: the original section `S1` | 442 | 4 |
| `T` | 63: `T^2(S1)` | 272 | 1 |

The forward minimum is attained by `q=S15`. One inverse minimizer is
`q=S10-S1`; all four integral words are recorded in the
[certificate packet](../../artifacts/generated-results/elliptic-curves/det1092_two_fibration_action_v1/).
An arbitrary further original translation does not change these degrees
or the rational image of the curve on the original parameter line.

This rules out a low-degree shortcut in these two infinite translation
families. It does not exclude other starting curves, other alternate
translating sections, longer alternating words, or rational incidence
on the high-degree curves. It neither constructs nor excludes a seed on302.

## Reconstructing the actual automorphisms

**Verified application of established elliptic-surface theory.** Work in
the full rational generic Neron--Severi lattice, with coordinates
`(a,b,v)=aO+bF+phi(v)` and pairing

\[
 (a,b,v)\cdot(c,e,r)=-2ac+ae+bc-v^tGr.
\]

The rootless original fibration has
`S_q=(1,q^tGq/2,q)`. Its translation action is therefore

\[
 t_q(a,b,v)=
 \left(a,\ b+v^tGq+\frac a2q^tGq,\ v+aq\right).
\]

This is an integral isometry and takes `S_v` to `S_(v+q)`.
The underlying section/height and trivial-lattice decomposition are
established theory; see [Schuett--Shioda, sections6 and11](https://arxiv.org/pdf/0907.0298).

**New verified application to the alternate fibration.** Its fibre is
`D=2O+4F+phi(w)`, with `w=S14-S15`, and its origin is `S14`.
The ten previously exhibited vertical section curves form exactly five
pairs. Each pair sums to `D`, its members have intersection two, and
different pairs are disjoint. These give five distinct two-component
fibres. Their component groups have exponent two; elliptic inversion
fixes each component class. The full inherited rank12 theorem leaves
no room for further independent reducible-fibre roots.

Let `V` be the span of `D`, the alternate origin, and one nonidentity
component from each pair. Its rank is seven and its Gram determinant
is32. Alternate inversion `iota` is identity on `V` and minus identity
on its orthogonal complement. Thus

\[
 \iota=2\operatorname{pr}_V-1.
\]

The exact matrix is integral. An independent eigenspace construction
checks it against the projector construction, verifies all component
actions, and checks the intersection form.

**New deduction identifying the translation.** On the pointed quartic,
the deck involution is the original-fibre operation
`j(P)=S_w-P`. On the alternate elliptic curve it is `j(A)=B-A`.
Consequently

\[
 T=j\circ\iota.
\]

These are actual automorphisms, not merely abstract lattice isometries.
Extension of elliptic inversion and section translation to the minimal
surface is established; compare
[Garbagnati--Salgado, section3.2.1](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/rank-jumps-and-multisections-of-elliptic-fibrations-on-k3-surfaces/602807FB00055D106E3CEA5418DE08F7).
The resulting action reproduces all six previously equation-verified
degrees `27,13,26,14,25,15`. It is not calibrated to an exceptional point.

## A quadratic formula for every original translation

**New deduction.** Write `C=(a,b,v)` and
`K=T^(-epsilon)F=(c,e,r)`. Then

\[
 \deg(\pi|_{C_{\epsilon,q}})
 =C\cdot K+(cv-ar)^tGq+\frac{ac}{2}q^tGq.
\]

Here `a=13` and `c=19` for both signs. In particular the exact formulas are

\[
 d_-(q)=1+\frac{247}{2}q^tGq-r_-^tGq,\qquad
 d_+(q)=63+\frac{247}{2}q^tGq-r_+^tGq,
\]

with

```text
r_minus = (6,-19,-6,-6,0,-12,0,0,6,-6,12,-6,0,12,9,3,0)
r_plus  = (32,-19,-32,-32,0,-64,0,0,32,-32,64,-32,0,64,48,16,0).
```

The coefficient247 is the product of the two relevant projection
degrees, `13*19`. This supplies a geometric explanation for the cost:
adjoining an original translation changes the second projection by a
large positive quadratic term, with only a linear cancellation term.
It does not merely change a chart presentation.

## Proving the minima without an unbounded calculation

**Verified application and new deduction.** The original generic lattice
is even and rootless, so a nonzero vector has norm four or at least six.
The exact closed norm-four ellipsoid contains2,436 nonzero vectors.
An independent rational LDL enumeration visits21,876 nodes, proves
completeness, and agrees with the separate PARI candidate producer.
No floating-point completeness assertion is used.

On this entire shell,

\[
 \max r_-^tGq=53,\qquad \max r_+^tGq=285.
\]

Therefore its degree minima are `1+494-53=442` and
`63+494-285=272`, with the multiplicities in the first table.

For any vector of norm `n>=6`, Cauchy--Schwarz gives

\[
 d_\epsilon(q)\ge d_\epsilon(0)+\frac{247}{2}n
                   -\sqrt{n\,r_\epsilon^tGr_\epsilon}.
\]

The right-hand side increases for `n>=6`. The exact norms are
`r_minus^t G r_minus=1216` and `r_plus^t G r_plus=31844`.
At six, its strict comparison with the claimed minima follows by squaring
positive quantities:

```text
(1 + 741 - 442)^2 - 6*1216  = 82704 > 0
(63 + 741 - 272)^2 - 6*31844 = 91960 > 0.
```

Thus no higher-norm translation improves either minimum. This proves
the global theorem over all of `M17`, including translations with
arbitrarily large coordinates. The statement and proof are intrinsic
under unimodular coordinate changes; no coordinate-prefix restriction
enters the final minimization.

## Evidence history, limits and constructive boundary

**Verified checkpoints.** The preliminary70-word calculation used zero
and both signs of the17 displayed sections. It motivated, but does not
prove, the global result. Two subsequent exact ellipsoids through degree63
each visited18 nodes and retained only `q=0`. The final short-shell plus
Cauchy argument supersedes the degree63 bound with global nonzero minima.
All three layers remain in the packet; no failed or partial history is
cleaned away.

The first minimizer failed before enumeration because this runtime's
`Fraction(Sage rational)` conversion retained a method-valued denominator.
The original source and protocol are retained. Version2 uses explicit
integer numerators and denominators with unchanged geometry and limits.

Every process completed within its25-second cap. There was no target
parameter evaluation, point search, new high-degree rational-map expansion,
class-group work, production modification or detached process. Only generic
section/lattice data and previous generic map certificates enter the theorem.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_two_fibration_gap.sage
```

**Open goal, not a conjectural conclusion.** We still do not have a
prospectively selected seed through302. The old carrier-preserving
obstruction and this degree barrier concern different operations and
neither is a universal arithmetic obstruction. They now justify avoiding
a blind expansion of this particular alternating-map route. A different
source/fibration operation or a direct admissible nongeneric Kummer class
remains necessary for the requested low-complexity seed construction.
