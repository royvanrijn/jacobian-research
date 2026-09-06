# A direct native-equation test gives a local solubility obstruction

The matched relation construction can now be rejected locally **without
computing its characteristic-zero intersection polynomial**. For control 2,
the native equations and generic section word give no relation point over
F131, with every surface and cover reduction condition verified. Proper
reduction therefore proves that the relation has no Q131 point.

The successful control survives. The other globally insoluble control has
points over both Q131 and Q137, showing explicitly why the new test is
necessary rather than sufficient. All conclusions concern the specified
relation constructions, not ranks or solubility of the original fibres.

## Inputs and fixed tests

Keep the three lattice-selected words S0,S1,S2 and the native triple
A=`01333`, B=`0b2d0`, D=`19e45` from the
[matched-translate experiment](MATCHED_TRANSLATES_SEPARATE_RATIONAL_SPLITTING_FROM_LATTICE_CAPACITY.md).
The relation is

\[
P_A-P_B+P_D=S_i.
\]

The detector receives only the original generic Weierstrass coefficients,
the seventeen published generic sections, the three primitive quadratic
forms, their scaled polynomial lift maps, and the selected generic words.
It receives no characteristic-zero intersection polynomial, exceptional
point, original rational parameter, or observed rank.

The [first protocol](DIRECT_RELATION_LOCAL_GATE_PROTOCOL.json) fixes primes
53, 61 and 73 from the preceding finite-scheme analysis. All fail the
surface-reduction gate and correctly return UNKNOWN. Their immutable
results are retained. The [second protocol](DIRECT_RELATION_GOOD_PRIMES_PROTOCOL.json)
then fixes the first four primes at least 131: 131,137,139,149. No failed
prime is replaced adaptively.

## Exact local obstruction criterion

Let p>3 and suppose the following conditions hold for the native equations:

1. All coefficients are p-integral, with the stated weighted degree bounds.
2. The original discriminant remains squarefree of degree 24 modulo p,
   and the fibre at infinity remains smooth.
3. The three primitive quadratic branch polynomials retain degree two,
   remain squarefree, and have pairwise disjoint branch divisors disjoint
   from the original discriminant.
4. All generic section and native lift identities hold modulo p.

The simultaneous carrier is the proper smooth genus-five curve

\[
C_3:\quad U_A^2=q_A(T,Z),\quad
U_B^2=q_B(T,Z),\quad U_D^2=q_D(T,Z)
\subset\mathbf P^4_{\mathbf Z_p},
\]

where all coordinates have weight one and the q's are homogeneous
quadratics. Every Qp point extends to a Zp point by properness and reduces
to a point of C3(Fp).

The polynomial native maps extend over both base charts. On the chart
T=1, their coordinates are obtained from the leading weighted coefficients:
x0 has weight four, x1 weight three, y0 weight six, y1 weight five, and
each cover root has weight one. Thus infinity is explicitly included.

Every singular original fibre modulo p is nodal. The generic sections and
native maps avoid its singular point. Indeed, the covers are unramified
above these fibres. At a hypothetical nodal image, differentiation of the
Weierstrass identity with respect to the local base parameter would force
the base partial derivative to vanish. That contradicts the simple zero
of the discriminant. The same argument applies to the published sections.

Consequently all maps land in the smooth group of the Weierstrass
fibration, including at the nodal fibres. The generic word S_i can be
specialized by adding the reduced published sections in that group. A
Qp solution of the relation must therefore reduce to a solution in C3(Fp).

Define the finite set

\[
R_{i,p}=\{c\in C_3(\mathbf F_p):
P_A(c)-P_B(c)+P_D(c)=S_i(t(c))\}.
\]

The proved criterion is

\[
\boxed{R_{i,p}=\varnothing\quad\Longrightarrow\quad
\text{the specified relation has no Qp point}.}
\]

The calculation enumerates every t in P1(Fp), every square-root choice,
and the exact group equation. It performs no rational-function elimination
or characteristic-zero factorization. Nonempty R is only a necessary
condition for local or global solubility; the detector does not promote
such points to Qp or rational points.

## Results, including the failed primes

| Prime | Surface/cover gate | #R0,p | #R1,p | #R2,p |
|---:|---|---:|---:|---:|
| 53 | UNKNOWN: discriminant not squarefree; branch A meets it | — | — | — |
| 61 | UNKNOWN: discriminant not squarefree; branch D meets it | — | — | — |
| 73 | UNKNOWN: discriminant not squarefree; branch D meets it | — | — | — |
| 131 | Pass | 1 | 1 | **0** |
| 137 | Pass | 1 | 1 | 1 |
| 139 | UNKNOWN: discriminant not squarefree | — | — | — |
| 149 | UNKNOWN: discriminant not squarefree | — | — | — |

The preceding finite intersection scheme gave valid local obstructions at
53 and 61. That does not imply that the ambient K3 and its native maps have
good reduction there. The failed geometry gates identify a limitation of
this direct implementation, not a contradiction or an exclusion.

At 131 the full genus-five carrier has 116 F131 points. Exactly eight lie
above the nodal original fibre t=42. Those eight were included in all
relation checks. At 137 the carrier has 120 F137 points and there is no
rational nodal original fibre. Infinity is checked at both primes and
has no carrier points. In total the successful panels cover all 270
projective base fibres and all 236 carrier points.

The known successful control at 137 occurs above t=4 with its A-root
equal to zero. This is a branch point in the reduction, despite the
nonzero rational roots at the characteristic-zero specialization. A test
that insisted on nonzero square values modulo p would incorrectly reject
this successful construction. Branch values must be included.

These finite-field counts are coverage certificates, not proposed rank
features. The carrier is identical for all three controls; the generic
relation word accounts for the different outcomes.

## Independent group and chart verification

The producer uses explicit modular chord-and-tangent formulas. The
independent verifier uses Sage elliptic-curve groups on smooth fibres.
On a nodal fibre it uses the normalization instead of repeating those
formulas. Writing the cubic as

\[
y^2=(x-r)^2(x+2r),\qquad k^2=3r,
\]

the smooth group is represented by

\[
\lambda(x,y)=\frac{y/(x-r)-k}{y/(x-r)+k},\qquad \lambda(O)=1.
\]

Addition becomes multiplication. For a nonsplit node, this calculation is
performed in Fp²; the rational smooth group is the corresponding norm-one
subgroup. At the observed nodal fibre modulo 131, r=20 and the node is
split. The independent normalization rechecks all eight carrier points
above it.

Every finite root tuple and every infinity tuple agrees between the two
implementations. Weighted degree bounds and polynomial identities are
checked explicitly, so the infinity computation does not rely on omitted
higher coefficients or a heuristic affine search.

Only after that independent equation-only replay does the verifier read
the previous characteristic-zero intersection certificates. At 131 and
137 their degree-twelve polynomials have squarefree reductions preserving
degree. Their simple root counts match the direct relation counts. Hensel
lifting then gives the separate conclusions:

| Relation | Q131 solubility | Q137 solubility | Previously certified Q solubility |
|---|---|---|---|
| S0 | Yes | Yes | No |
| S1 | Yes | Yes | Yes |
| S2 | No | Yes | No |

The first row is a concrete demonstration that passing this finite local
panel does not close global solubility. It remains obstructed at 61 by
the earlier finite-scheme certificate.

## Mathematical use and remaining gap

This is a **solubility** gate for a specified simultaneous relation, derived
from generic equations before characteristic-zero elimination. The absence
of residue solutions is an exact exclusion of that construction. It is
not an exclusion of exceptional points on the original fibre, nor of other
relations on the same carrier.

The negative controls show that the relevant arithmetic is more than the
trace norms, carrier genus, or generic section capacity. The positive
control still requires the rational Galois-fixed component established in
the preceding characteristic-zero calculation, followed by its independent
two-direction quotient certificate. Local survival alone supplies neither.

For eventual use by Agent 1, the detector could cheaply discard proposed
relation constructions before compiling their rational intersection maps.
It must retain the bad-reduction UNKNOWN outcome and include infinity,
branch roots, and smooth points on nodal fibres. It has not been installed
in a selector or search pipeline. The cover triple remains retrospective,
and the full +7 jump is still unexplained beyond its previously certified
rank-three quartet.

The next substantive gap is the global lift from the minimal genus-one
pair carrier to the genus-five three-cover carrier, beyond a particular
chosen relation word. Characterizing that lift's Jacobian/cover obstruction
would address simultaneous solubility itself; a larger residue-prime panel
would only improve this necessary filter.

## Replay and immutable evidence

Both workers had a 60-second bound. Sage 10.9 and PARI 2.17.3 were used;
the detector itself uses finite-field polynomial identities and modular
group arithmetic, without a descent or point search.

```sh
sage -python elliptic-curves/rank-jump/direct_relation_local_gate.py check
sage -python elliptic-curves/rank-jump/direct_relation_good_primes.py check
sage -python elliptic-curves/rank-jump/verify_direct_relation_local_gate.py check
```

Immutable [generic-only input](../../artifacts/generated-results/elliptic-curves/rank_jump_direct_relation_local_gate_inputs_v1.json),
[initial UNKNOWN panel](../../artifacts/generated-results/elliptic-curves/rank_jump_direct_relation_local_gate_v1.json),
[fixed larger panel](../../artifacts/generated-results/elliptic-curves/rank_jump_direct_relation_good_primes_v1.json),
and [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_direct_relation_local_gate_verification_v1.json)
bind the scripts, protocols and source equations. Active search files,
worker policies, and mathematical status entries remain untouched.
