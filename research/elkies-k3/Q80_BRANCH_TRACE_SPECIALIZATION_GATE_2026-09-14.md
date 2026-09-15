# A common halving defect is necessary at every branch fibre

The trace-norm formulation gives a concrete arithmetic obstruction on the
branch fibres. Let M be the full generic Mordell--Weil group, let a quadratic
cover have branch points b, and write

```
K_b = ker(M/2M -> E_b(k(b))/2E_b(k(b))).
```

Every realized trace parity belongs to the **same intersection** of all K_b.
For Q80 under the good-branch-reduction hypotheses of the
[trace-norm theorem](Q80_GOOD_BRANCH_CODE_AND_TRACE_NORMS_2026-09-14.md),

```
genus0: rank gain <= dim(intersection_b K_b),
genus1: rank gain <= dim(intersection_b K_b)+1.          (1)
```

Thus correlated gain requires a shared loss of generic parity under
specialization, through an additional half or a specialization relation:
at least two common parity directions in genus0, and at least one in
genus1 in this reduction scope. A class specializing to O is included.
This is necessary, not a substitute for global trace-torsor solubility.

On the literal Q80 parent, **K_0=0**. A small exact finite-group certificate
proves it. Consequently any good-reduction quartic cover ramified at t=0
has gain at most one; a good-reduction quadratic branch divisor containing
zero has gain zero. No height or coefficient bound is imposed.

There is also an explicit genus-one control with infinitely many rational
points, an actual MW17 parent, and a compatible nonzero131 branch code,
but **exactly zero generic gain**. It illustrates the global obstruction
that a single residue-code calculation misses. The
[positive two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open.

## 1. Branch specialization of a trace

Let K=Q(t), L=Q(C), G=E(L), and let sigma be the deck involution. Assume
every characteristic-zero branch fibre is smooth. At the unique point
above a closed branch value b, a section Q specializes to Q_b in E_b(k(b)).
Since sigma fixes that point and its residue field,

```
Tr(Q)_b = (Q+sigma(Q))_b = 2Q_b.                       (2)
```

Specialization is a group homomorphism on the smooth fibre, including O.
It follows immediately that

```
W := Tr(G)/2M  subset  intersection_b K_b.              (3)
```

This is stronger than asking independently for some halving defect at
each fibre: the generic parity class must be the same at every branch.
Equivalently, the entire branch divisor must lift, over its residue fields,
to the division curve 2R=T for each realized trace parity T. Two independent
trace parities require simultaneous lifts to both division curves.
These branch lifts are only necessary; they need not extend to points of
the function-field torsors H_T(m)=d*s^2.

For rootless Q80, the retained unramified Kummer theorem identifies the
branch-zero anti-invariants with actual anti-traces. The resulting exact
rank formula is

```
rank gain = dim W + dim(actual branch image).          (4)
```

Use the previous good-reduction code bound to obtain (1). Without that
local reduction hypothesis, (3) and (4) still hold in the smooth-branch
setting, but the numerical bound on the branch image must be supplied
separately. In particular, if every branch fibre has no nonzero k(b)-rational
2-torsion, the branch image is zero at every height and every base genus.
Then one branch with K_b=0 forces zero gain, even without good branch
reduction at131.

## 2. Exact specialization certificate at zero

The [producer](scripts/audit_q80_branch_trace_specialization.py) uses only
the retained direct11952 equation and its full saturated17-section basis.
The parameter t=0 and the ascending prime pool5 through997 are fixed in
the script before arithmetic. No exceptional point, desired specialized
rank, or target j-invariant is used. The only early stop is exact row rank17.

For a good odd prime p and a rational cubic root e, the standard2-Kummer
character is x-e modulo squares, with value3e^2+A at (e,0), and the square
class at O. It is a homomorphism killing doubles. Applying it to the17
specialized generic sections gives a binary row. The retained matrix has
24 rows and rank17;17 independent rows use these14 primes:

```
17,43,53,59,71,73,83,89,101,109,113,127,139,157.
```

The [independent checker](scripts/verify_q80_branch_trace_specialization.py)
does not use those characters to establish rank. At each witness prime it
enumerates the entire finite elliptic group, constructs its subgroup of
doubles, and labels the quotient by explicit cosets. The specialized17
points give quotient matrices of combined rank17. The quotient row spaces
also agree with the producer's character row spaces. Exact rational point
identities are rechecked at t=0. Parent rank and generic saturation are
inherited from the retained parent theorem.

If an odd generic parity combination specialized to twice a rational
point, it would vanish in every one of these finite quotients, contradicting
rank17. Therefore K_0=0. No upper bound for the rational rank of E_0 is
asserted; the result concerns the specialized generic subgroup modulo two.

The [result](../artifacts/generated-results/elkies-k3-q80-branch-trace-specialization-v1/result.json)
records the exact specialized points, all character rows, selected pivots,
source hashes, skipped primes and declared15 CPU-second/1GiB limit.
Producer time was about0.005 seconds. The
[independent replay](../artifacts/generated-results/elkies-k3-q80-branch-trace-specialization-v1/independent-replay.json)
records every finite-group order and quotient, in about0.041 seconds under
20 CPU seconds/1GiB.

## 3. A whole rational neighbourhood has the same obstruction

The certificate extends beyond the single parameter. Normalize every
rational-function denominator of the generic sections by its constant
coefficient. The checker verifies that at each of the14 witness primes
the constant terms are integral and every positive-degree coefficient c_i
satisfies

```
v_p(c_i)+i >= 1.
```

The same check applies to A and B. Thus for every rational b in p*Z_p,
the literal curve and all17 generic points reduce to their t=0 values.
Denominators remain units. All14 quotient witnesses apply simultaneously
whenever

```
b in intersection_(p in witness list) p*Z_p.            (5)
```

In particular all integer multiples of

```
M0=301744286059269190181732689
```

satisfy K_b=0. The statement also includes rational b satisfying (5),
not only integers. Such fibres are smooth because their reductions are
smooth at the witness primes. This is a certified arithmetic neighbourhood,
not a bounded parameter sweep or a global exclusion of all rational branches.

## 4. Infinite elliptic base, admissible residue code, but no gain

Take the prospective literal quartic

```
d=t*(t-2)*(t-76)*(t-103),
C: w^2=d.
```

It is squarefree and its four branch fibres on Q80 are smooth. Modulo131,
the assignment O at0 and roots54,3,26 at2,76,103 has the inherited codes

```
0,71746,50167,121781,
71746 XOR 50167 XOR 121781 = 0.
```

This is one of the four compatible three-site patterns from the previous
certificate. The global obstruction is visible at other primes:

| Rational branch | Root-free cubic witness prime | A mod p | B mod p |
|---:|---:|---:|---:|
|0|23|3|22|
|2|29|2|13|
|76|37|25|17|
|103|31|0|24|

All four reductions are smooth and their monic cubics have no root. A
rational root would be integral at that prime and reduce to a root, so
none of these branch fibres has rational nonzero2-torsion. Thus the actual
branch image is zero. Since K_0=0, (3)--(4) prove

```
rank E(Q(C)) = rank E(Q(t)) = 17.                     (6)
```

The base really has infinitely many rational points. Put k=-15656. The
birational maps

```
X=k/t,       Y=k*w/t^2,
t=k/X,       w=k*Y/X^2
```

identify its smooth projective model with

```
Y^2=X^3+8186*X^2+2833736*X+245110336.
```

The point P=(0,-15656) has

```
2P=(17/4,128325/8).
```

The integral short model obtained by x=9X+3*8186, y=27Y has integral
coefficients (81*a4-27*a2^2,729*a6-243*a2*a4+54*a2^3). On it,
x(2P)=98385/4 is not integral. The Lutz--Nagell integrality theorem therefore
proves that P is nontorsion. The maps and group-law calculations use exact
rational arithmetic; no numerical rank or BSD assumption is involved.

Consequently this is a control with an actual arithmetic MW17 parent and
an infinite genus-one rational base, but no new generic section direction.
It does not bound the ranks of its individual specialized rational fibres.
The control was selected from generic branch-code data and the fixed zero
branch; it is not a proposed rank discovery.

## Next positive gate and proof boundary

For a positive construction, first seek a branch divisor whose fibres have
a **common** nonzero halving defect in M/2M, or two independent common defects
when the branch contribution is zero. Then solve the corresponding global
trace torsors on the same cover. Alternatively, a branch2-torsion contribution
requires actual residue-field2-torsion, not merely roots modulo one prime.
Neither finite local compatibility nor an infinite rational base closes
these lifting obligations.

```
python3 research/elkies-k3/scripts/verify_q80_branch_trace_specialization.py
```

The finite quotient calculation and rational identities replay independently.
The smooth-branch trace specialization, unramified descent and rank accounting
are written proofs using the retained hypotheses. No formal verification,
external review, large search or new positive MW17 gain is claimed.
