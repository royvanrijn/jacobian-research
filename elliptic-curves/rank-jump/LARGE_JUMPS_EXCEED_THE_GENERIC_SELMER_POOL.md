# Large jumps exceed the generic geometric Selmer pool

For the families in the frozen fresh27/historic comparison, classes defined
over `Q(t)` and satisfying the **geometric local Kummer conditions at every
base place** have room for at most **three dimensions beyond the marked
generic subgroup**. This remains true when their cubic fields vary with t.
Consequently they cannot, on their own, explain the recorded +10, +11, +12
or +14 jumps.

This is an incidence-capacity theorem for a precisely restricted source of
classes. It neither computes the additional arithmetic Selmer quotient nor
proves a criterion for its rational solubility. The
[point-masked governing/CT comparison](FRESH_RANK27_GOVERNING_AND_CT_COMPARISON.md)
still has UNKNOWN additional-quotient CT entries.

## Exact family and fibre checks

The equation/generic-only [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_selmer_capacity_verification_v1.json)
checks seven family presentations and all sixteen specialized models. Each
family is a minimal elliptic K3 surface with section, and its geometric
two-division monodromy is S3. Each model transport is checked by exact
`A_target = u^4 A_family(t)`, `B_target = u^6 B_family(t)` with rational u.
The marked generic points' mod-two independence is replayed using their
retained good-prime Kummer characters. Exceptional coordinates never enter.

| Family presentation | Geometric singular fibres | Root-lattice rank c | Geometric Selmer dimension d |
|---|---|---:|---:|
| 074d9, 103b2, 11952 | 24 I1 | 0 | 20 |
| a1-fibration-01, a1-fibration-05 | 22 I1 + I2 | 1 | 19 |
| historic published-R17 lineage, 074d9 chart | 24 I1 | 0 | 20 |
| curve398-p16875 | 22 I1 + I2 | 1 | 19 |

The last two configurations are checked from their own equations, rather
than inferred from their labels. Curve398's original presentation has its
I2 fibre at infinity. The initial worker, which required good infinity,
returned UNKNOWN. A separate immutable completion reverses the coefficient
arrays under `t=1/u`, verifies that inversion exactly, and passes the same
geometric tests. This changes coordinates on the same family and retained
specialization. The unsupported original row remains preserved.

The K3 assertion follows from a minimal degree-(8,12) Weierstrass model
with discriminant degree 24 after this coordinate change: the fundamental
line bundle has degree two, the canonical bundle is trivial and q=0.
Coprimality with c4 makes every finite bad fibre multiplicative; good
infinity and the discriminant multiplicities account for all Euler number
24. A retained irreducible cubic specialization and nontrivial odd inertia
give geometric S3 monodromy, as in the
[earlier geometric argument](FIXED_CUBIC_TRANSFER_REQUIRES_HIGH_GENUS.md).

## The restricted pool and its capacity

Put K=Q(t), F=Qbar(t), and let E/K be the generic elliptic curve. For every
geometric base place v, write F_v=Qbar((z_v)). Define

\[
 S_{\rm geom}=\{\alpha\in H^1(F,E[2]):
 \operatorname{res}_v\alpha\in\delta_v(E(F_v)/2E(F_v))\ \text{for every }v\},
 \qquad
 L=\operatorname{res}^{-1}(S_{\rm geom})\subset H^1(K,E[2]).
\]

These local fields are completions **in the parameter direction**. They
are not Q_p. In particular L is not being called the arithmetic Selmer
group over a number field. Its specializations need not be Q_p-soluble.

Let rho be the geometric Picard rank and let
`c = sum_v(number of geometric fibre components - 1)`. Geometric S3
monodromy implies E(F)[2]=0. The Kummer exact sequence gives

\[
0\longrightarrow E(F)/2E(F)\longrightarrow S_{\rm geom}
\longrightarrow\Sha(E/F)[2]\longrightarrow0.
\]

Here Sha uses the local conditions at **all** geometric base places. For
an elliptic surface with section, geometric Sha identifies with the
surface's Brauer group. For a characteristic-zero K3 surface its 2-torsion
dimension is 22-rho; Shioda-Tate gives geometric MW rank rho-2-c. These
are standard results; see Huybrechts,
[Lectures on K3 surfaces](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf),
Corollary 11.3.4, Proposition 11.5.6 and §18.1, especially (1.14).
Only the locally trivial Sha group is used, not the unrestricted
Weil-Chatelet group. Therefore

\[
\boxed{\dim S_{\rm geom}=(\rho-2-c)+(22-\rho)=20-c=d.}
\tag{1}
\]

No exact Picard rank is assumed. The certificate checks the cancellation
for every possible characteristic-zero K3 Picard rank compatible with c.
The geometric MW group is finitely generated here, and absence of
2-torsion makes its mod-two dimension equal its rank.

Inflation-restriction for K⊂F has kernel
`H^1(Q,E[2](F))`. This is zero because E[2](F)=0. Consequently restriction
on H^1 is injective, giving

\[
\boxed{\dim L\leq d.}
\tag{2}
\]

Let G be the Kummer span of the marked generic sections. The retained
generic-family certificates license this identification; the finite-fibre
mod-two checks give dimension m=17 or 16, respectively. Their Kummer
classes belong to L, so `dim(L/G) <= d-m = 3` in every family above.
This does not assert that three additional arithmetic classes exist.

## Specialization and the panel

At a smooth rational base value t0, E(Qbar((t-t0)))/2 is zero. Indeed both
the good-reduction special-fibre group and the formal group are
2-divisible. Thus every element of L is geometrically unramified at t0.
Since E[2] extends as a finite etale group scheme there, unramified descent
gives a canonical linear specialization

\[
 s_{t_0}:L\longrightarrow H^1(\mathbb Q,E_{t_0}[2]).
\]

Equivalently, restrict to the henselian parameter disc and use
inflation-restriction from its residue field Q. Inertia acts trivially
on E[2], and the geometric local restriction vanishes, so the class is
the inflation of a unique residue-field class. This also shows why an
arbitrary ramified H^1 class need not have such a specialization.

Write P_t for the rational Kummer image E_t(Q)/2E_t(Q), and V_t=s_t(L).
If a retained certificate proves rank(E_t)>=R, then dim P_t>=R, whether
or not a particular displayed subgroup is 2-saturated. Hence

\[
\boxed{\dim\bigl(P_t/(P_t\cap V_t)\bigr)\geq\max(0,R-d).}
\tag{3}
\]

| Retained group of fibres | R | m | d | Forced rational Kummer dimensions outside V_t |
|---|---:|---:|---:|---:|
| Seven fresh R17 highs, cases 00/02/04/06/08/09/10 | 27 | 17 | 20 | at least 7 |
| Fresh MW16 high, case 11 | 27 | 16 | 19 | at least 8 |
| Four R17 low controls, cases 01/03/05/07 | 17 | 17 | 20 | lower bound 0 |
| MW16 low control, case 12 | 17 | 16 | 19 | lower bound 0 |
| Historic 356 and 385, cases 13/14 | 29 | 17 | 20 | at least 9 |
| Historic 398, case 15 | 30 | 16 | 19 | at least 11 |

The [comparison artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_selmer_capacity_comparison_v1.json)
joins these retained rank labels only after the point-masked computation.
The last column is a consequence of labels, **not a prospective feature**.
A zero lower bound for a censored low control does not prove absence of
extra classes. Exact full ranks remain UNKNOWN in this panel.

## A falsifiable obstruction for larger global blocks

A construction can evade (2) by allowing geometric local obstructions.
For a finite-dimensional proposed block W⊂H^1(K,E[2]), define the linear
obstruction map

\[
 o:W\longrightarrow\bigoplus_v
 H^1(F_v,E[2])/\delta_v(E(F_v)/2E(F_v)).
\]

Its kernel is W∩L. If W specializes to at least R independent Kummer
classes at a smooth fibre where it is unramified, then dim W>=R, and

\[
\boxed{\operatorname{rank}(o|_W)\geq R-d.}
\tag{4}
\]

Thus a candidate global explanation for a full fresh R17 rank-27 block
needs obstruction rank at least seven; the MW16 block needs at least
eight. This is a different invariant from the pair masks or collision
defects already measured at rational primes.

There is also a support budget. At a good geometric base place the local
quotient has dimension two: inertia is procyclic, E[2] is constant, and
the local Kummer image is zero. At an I_n place its dimension is one.
To see the latter, Tate uniformization gives local Kummer dimension zero
for odd n and one for even n. Tame monodromy on E[2] is respectively a
transvection or the identity, so H^1 has dimension one or two. Subtracting
again gives one. Consequently if all obstructions of W are supported
on g good and b bad **geometric base places**, then

\[
\dim W\leq d+2g+b,\qquad 2g+b\geq R-d.
\tag{5}
\]

Closed points over Q are counted with their geometric degree. If no bad
base place is used, the fresh highs require at least four good geometric
places, the historic +12 controls at least five, and +14 at least six.
These are necessary bounds, not measured support sizes or sufficient
solubility conditions. A common low-genus auxiliary curve may still
exist: the result neither bounds the genus of all auxiliary constructions
nor applies unchanged after a base extension.

The small falsifiable next test for any proposed family-level cover block
is therefore to compute this obstruction matrix before searching for
rational points. If its rank is too small, that block cannot explain the
whole jump. A successful matrix passes only an incidence-capacity gate;
the specialization still needs arithmetic local tests, CT and a rational
solubility argument. No such additional block is supplied by this note.

## What this changes in the mechanism search

The strongest surviving route is to produce **specialization-dependent
arithmetic classes in the varying cubic field**, then explain why several
become rational together. The earlier
[fixed-incidence six-direction twist switch](FIXED_INCIDENCE_SIX_DIRECTION_SOLUBILITY_SWITCH.md)
proves that a large solubility change is possible once a suitable class
block exists. The missing implication on the fresh high/low panel remains
an equation-only construction of its additional classes and CT structure.

A ramified Q(t) construction or an auxiliary base change is another
possible route, now subject to (4)-(5). A fixed, everywhere geometrically
locally soluble pool over the original Q(t) is insufficient. A fixed
cubic field obtained by low-genus base change is also excluded by the
previous genus bound. Inherited generic CT-switch rank, the universal
governing-field degree and search-coordinate visibility remain weak or
irrelevant explanations of the additional rank.

Agent1 can eventually use a certified additional-class incidence bound
and a solubility obstruction on that block. The numbers in this note
should **not** be added to a candidate score: they are constant within a
family, or are derived retrospectively from the known rank label. No
active search policy, worker, candidate or mathematical-status entry was
changed.

## Reproduction and boundaries

The [protocol](GENERIC_SELMER_CAPACITY_PROTOCOL.json),
[worker/exporter](generic_selmer_capacity.py),
[base-inversion completion](complete_generic_selmer_geometry.py) and
[portable verifier](verify_generic_selmer_capacity.py) are narrowly scoped.
The two new historical geometry workers and the one coordinate completion
each had a 30-second cap. There were no new parameter or elliptic point
searches and no class-group computations.

Rational Euclid in the first verification attempt exhausted its 60-second
cap on large historical coefficients, before producing a certificate.
The final verifier uses degree-preserving coprime reductions modulo fixed
small primes to certify the same rational coprimality claims; polynomial
product and transport identities still use exact rational arithmetic.
It finishes within one second. The input whitelist and mathematical test
are unchanged.

```sh
timeout 60 python3 elliptic-curves/rank-jump/verify_generic_selmer_capacity.py check
```

The verifier establishes the polynomial hypotheses, all sixteen model
transports, generic mod-two independence and the numerical consequences
of the stated standard theorems. The cohomological proof is mathematical
prose, not a formally checked theorem. The exact arithmetic dimension of
L, its image at each fibre, and additional-quotient CT are all UNKNOWN.
