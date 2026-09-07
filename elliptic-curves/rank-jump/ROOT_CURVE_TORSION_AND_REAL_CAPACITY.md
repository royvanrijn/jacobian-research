# The cubic root curve carries global incidence; real topology cuts its capacity

The original-parameter global cover pool has room for at most **two**
dimensions beyond the marked generic subgroup on every family in the
sixteen-fibre panel. This improves the previous bound of three. The
improvement comes from exact real topology of the cubic root curve, without
using exceptional points. A preceding finite-field parity test did not
improve the bound.

The root curve is an **incidence carrier** for a restricted set of global
classes. Rational 2-torsion on its Jacobian supplies elliptic 2-cover
classes, not rational points on those covers. It is not the simultaneous
rational-solubility carrier requested for the large exceptional quotient.
That distinction is essential here.

## A concrete description of the restricted global pool

Keep the notation of the
[generic Selmer-capacity theorem](LARGE_JUMPS_EXCEED_THE_GENERIC_SELMER_POOL.md).
For E/Q(t), let L⊂H1(Q(t),E[2]) consist of classes whose restrictions to
Qbar(t) satisfy the local Kummer condition at every geometric base place.
Write C for the smooth projective normalization of

\[
 x^3+A(t)x+B(t)=0,
\]

and J=Jac(C). This is a connected degree-three cover C→P1, not an elliptic
curve. Its genus is 10 for the 24-I1 families and 9 for the 22-I1-plus-I2
families. The latter plane/ruled-surface model has an ordinary node over
the I2 base place. Those genera and fibre types are independently verified
in the preceding work.

For the R17 families one has the natural identification

\[
\boxed{L\simeq J[2](\mathbb Q).}
\tag{1}
\]

For the MW16 families there is an exact sequence at its first two terms

\[
\boxed{0\longrightarrow J[2](\mathbb Q)\longrightarrow L
       \stackrel{\epsilon}{\longrightarrow}\mathbb F_2.}
\tag{2}
\]

Surjectivity onto the last F2 is not asserted. It records the single
allowed odd-valuation pattern at the two normalization branches over
the I2 node. Equations (1)-(2) specify a useful arithmetic object even
when the cubic field varies with t.

Here is the divisor proof, including the descent issue over Q. In the
standard cubic Kummer description, H1(Q(t),E[2]) is the norm-square
kernel in Q(C)^*/Q(C)^*2. At a good geometric base place, the local Kummer
image is zero. It is also zero at I1: Tate uniformization over Qbar((z))
gives E/2=0 when the discriminant valuation is odd. Thus these conditions
force even valuations at all corresponding places of C.

At I2 the local E/2 has dimension one. Its nonzero cubic valuation
pattern is (1,1,0): the two roots meeting at the node carry odd valuation,
and the third carries even valuation. Both patterns (0,0,0) and (1,1,0)
are invariant under exchanging the two branches, so epsilon takes values
in one F2 even if the branches are conjugate over Q. Its kernel is exactly
the norm-square classes with even valuation everywhere on C.

For such a class alpha, write div(alpha)=2D. The map alpha↦[D] lands in
J[2](Q). If [D]=0, then alpha=c*h^2 for c∈Q^*. The norm condition makes
c^3 a square in Q(t), hence c a square in Q. This proves injectivity.

Conversely, a Q-rational 2-torsion point of J has a rational divisor
representative. Indeed, its Picard descent obstruction in Br(Q) is
killed by two, and it is killed by the degree-three rational zero-cycle
given by a fibre of C→P1. It therefore vanishes. A rational line bundle
has a rational divisor representative, for example after adding a large
multiple of that degree-three divisor and applying Riemann-Roch. Choose
alpha with div(alpha)=2D. Its norm has even divisor on P1, so
N(alpha)=c*h(t)^2 for c∈Q^*. Replacing alpha by c*alpha makes its norm
c^4*h(t)^2, a square. All geometric local square conditions now hold.
This proves (1) and the kernel assertion in (2). No rational point on C
or on an exceptional elliptic cover was used.

## Bounded finite-field test: no improvement

The [protocol](ROOT_CURVE_FROBENIUS_CAPACITY_PROTOCOL.json) selects the first
three good odd primes at most 499 for each of the six panel family
presentations. The first three are 131,137,151 in every case. Selection
checks smoothness of the normalized root curve through the preserved
discriminant factors, c4 coprimality, good infinity and ordinary I2 nodes.
It does not depend on the point counts.

| Family | #C(F131) | #C(F137) | #C(F151) |
|---|---:|---:|---:|
| 074d9 | 124 | 148 | 140 |
| 103b2 | 128 | 138 | 154 |
| 11952 | 110 | 150 | 150 |
| a1-fibration-01 | 130 | 150 | 160 |
| historic R17 lineage, 074d9 chart | 124 | 148 | 140 |
| curve398-p16875 | 130 | 150 | 160 |

These are counts on the normalized cubic root curves, not on elliptic
fibres. The raw nodal count is corrected by +1 for a split ordinary node
and -1 for a nonsplit node. At infinity, x is replaced by x/t^4.

Why could parity have helped? The generic sections and (1)-(2) force
dim J[2](Q)>=2g-3. At good odd p this subgroup injects into J[2](Fp).
The characteristic polynomial of Frobenius on J[2] is reciprocal and has
a factor (T+1)^(2g-3). Its remaining reciprocal cubic is determined by
the trace parity. If #C(Fp) is odd, it is T^3+1, giving total multiplicity
2g-2 for eigenvalue one and hence dim J[2](Q)<=2g-2. With the I2 allowance,
that would leave at most one extra global dimension beyond G.

Every observed count is even. The mod-two characteristic polynomial at
these primes is therefore (T+1)^(2g); it yields no improvement. This is
not evidence that all 2-torsion is rational. The verifier constructs
symplectic matrices with that same characteristic polynomial and each
fixed-space dimension 2g,2g-1,2g-2,2g-3, by putting zero through three
transvections on separate hyperbolic planes. Thus the ambiguity is exact,
not merely a limitation of numerical precision. These three samples do
not determine the full Galois image, either.

## Exact real topology: one global dimension is excluded

The second [protocol](ROOT_CURVE_REAL_COMPONENTS_PROTOCOL.json) computes
components of C(R) from the ordered real cubic roots. It isolates all
real zeros of rad(discriminant)*B in rational intervals. B has constant
nonzero sign around every discriminant root. At a simple branch value
the double root is -3B/(2A), whose sign is the sign of B because A<0.
Positive B merges the upper two roots; negative B merges the lower two.

The intervals between branch values have one or three real roots according
to the discriminant sign. Their strands are glued through the folds.
The transition x↦x/t^4 preserves their order at infinity. At an I2 node,
normalization uses crossing branches if they are real, and removes the
isolated real node if its branches are conjugate. For a node (t0,x0),
the tangent discriminant is

\[
 A'(t_0)^2-6x_0\bigl(A''(t_0)x_0+B''(t_0)\bigr).
\]

Both MW16 nodes are nonsplit over R. They occur at t=-2 in the compact
a1-fibration-01 model and u=0 in curve398's inverted presentation. Counting
the isolated singular point as a real component would give a wrong bound.

| Family type | Genus g | Real I1 branch values | Real I2 nodes | Components s of the normalization | dim J[2](R) |
|---|---:|---:|---:|---:|---:|
| Each R17 presentation | 10 | 18 | 0 | 10 | 19 |
| Each MW16 presentation | 9 | 16 | 1 nonsplit | 9 | 17 |

For a smooth real curve with s>0 components, dim J[2](R)=g+s-1.
See Kummer, [A signed count of 2-torsion points on real abelian varieties,
§5](https://arxiv.org/pdf/2301.10621), where the components are indexed
0 through s and the notation consequently has one fewer than our s.
One can also see the dimension from J(R)^0≅(R/Z)^g and its component
group (Z/2)^(s-1): every component contains a 2-torsion point.

Rational 2-torsion is real. Applying (1)-(2) therefore proves

\[
\boxed{\dim L\leq19\quad(\mathrm{R17}),\qquad
       \dim L\leq18\quad(\mathrm{MW16}),\qquad
       \dim(L/G)\leq2.}
\tag{3}
\]

The bounds are unconditional. Equality, a rational 2-torsion basis and
the extra classes themselves remain UNKNOWN.

## Consequence for the large-jump panel

Specialization of L at a smooth rational parameter is defined as in the
preceding theorem. Joining the retained rank lower bounds only after the
equation-only computations gives:

| Fibres | Retained rank >= | Forced rational Kummer dimensions outside the global pool >= |
|---|---:|---:|
| Seven fresh R17 highs | 27 | 8 |
| Fresh MW16 high | 27 | 9 |
| Historic R17 356 and 385 | 29 | 10 |
| Historic MW16 398 | 30 | 12 |
| Same-family censored low controls | 17 | 0 |

These necessities are recorded for all sixteen fibres in the
[comparison certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_root_curve_capacity_comparison_v1.json).
The zero lower bound for a low control is not an exact low-rank theorem.
The high-case numbers depend on known rank labels and must not become
prospective selectors. The real component count is constant within each
family; it cannot discriminate its high and low specializations.

For a Q(t) block W with enough classes to specialize to the full retained
rank, the kernel of its geometric local-obstruction map is W∩L. Thus the
necessary obstruction-rank bounds from the preceding note improve from
7/8/9/11 to **8/9/10/12**, respectively. This does not identify those
obstructions or supply a sufficient condition on t.

## Implications for the next mechanism test

The class carrier is now explicit: the cubic root Jacobian's rational
2-torsion, with at most one extra node-parity bit in MW16. The full jump
still exceeds that pool. The useful remaining avenues are, in order:

1. Specialization-dependent classes in the varying cubic number field,
   with an independently computed simultaneous-solubility obstruction.
   This remains the direct missing input to the fresh/historic CT panel.
2. Ramified family-level cover constructions or auxiliary base changes
   with enough independent local-obstruction dimensions. Their rational
   solubility must be proved separately.
3. Computing the actual Galois fixed space of J[2], if a concrete model
   makes it cheap. Ordinary trace parity did not determine that space in
   the frozen test. Enlarging the same prime sample is not the next task.

The weak explanations remain inherited generic CT switches, universal
governing-field degree and visibility. This note adds no candidate score
for Agent1. It narrows the class-construction problem without resolving
the missing implication from additional incidence to rational solubility.

## Verification and non-interference

The [independent verifier](verify_root_curve_capacity.py) counts finite
roots using gcd(f,x^p-x), verifies the real isolation intervals and their
completeness with PARI Sturm counts, and recounts topology components by
graph traversal. The workers use direct finite-field root enumeration,
Sage rational root isolation and a disjoint-set graph calculation.
All twelve workers completed within their separate 30-second caps; the
independent replay completes within one second.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_root_curve_capacity.py check
```

The [Frobenius artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_root_curve_frobenius_capacity_v1.json),
[real-topology artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_root_curve_real_components_v1.json)
and [verification](../../artifacts/generated-results/elliptic-curves/rank_jump_root_curve_capacity_verification_v1.json)
are immutable. The Picard-descent and real/Frobenius implications above
are mathematical proofs, not formal proof-assistant certificates.
No rational elliptic points, prospective parameters, class groups, active
search files or mathematical-status entries were changed or searched.
