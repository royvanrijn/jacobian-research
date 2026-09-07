# Compact S-unit corrections: evaluation works, research bases remain unknown

The fixed 103b2 high/low pair still has no independently constructed
S-unit basis. A bounded alternative to the previous default class-group
calls avoided their immediate stack exhaustion and reached relation
collection, but both calls timed out without returning a BNF or unit
basis. No new strict class or rank comparison follows.

The positive control did complete. Its compact S-unit factorizations
were certified and evaluated in local squareclasses without expanding
them in the main evaluation path. This prepares the correction step of
[the generic-half-ideal incidence theorem](GENERIC_IDEALS_AND_SUNIT_CORRECTIONS_CARRY_INCIDENCE.md);
it does not supply its missing research-field inputs.

## Fixed experiment and provenance

The [protocol](COMPACT_SUNIT_BASIS_PROTOCOL.json) retains precisely the
already reduced cubic fields of 103b2 at 3726/881 and -1049/2296. Its
worker input projects only the defining cubic, certified integral basis,
field discriminant and complete S. No generic or exceptional point,
Artin character, rank label, or active-search output is projected into it.

The calls use PARI 2.17.3, bnfinit flag 1, tech=[0.03,0.03,-1], one
worker, thirty seconds per field, and the previous 256 MiB stack limit.
Flag 1 retains compact unit information; the technical settings request
a smaller factor base and random relations. These interfaces and their
roles are documented in
[Allombert's PARI tutorial, slides 2–3 and 9–10](https://pari.math.u-bordeaux.fr/Events/PARI2026/talks/bnfinit.pdf).

The explicit small control is K=Q(cuberoot(11)), with S above 2,3,11.
It ran through bnfinit, bnfunits, and bnfcertify before either research
field. Source and protocol hashes, complete logs and terminal states
are retained in the
[artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_compact_sunit_basis_v1.json).

| Input | Initial factor-base columns in log | Relation request during retained random phase | Returned BNF | Returned S-unit basis |
|---|---:|---:|---|---|
| Small cubic control | recorded in artifact | completed | certified | certified |
| 103b2, 3726/881 | 350 | 273, unchanged through timeout | no | no |
| 103b2, -1049/2296 | 421 | 328, unchanged through timeout | no | no |

PARI later reduces the working ideal counts to 268 and 323. Its internal
secondary bounds are larger than the initial factor-base cutoffs; the
logs retain these details rather than interpreting the technical input
as a fixed mathematical generating-set bound. Neither research run
returned an object before its thirty-second timeout. The processes are
terminal and are not left running in the background.

The unchanged logged relation requests locate the observed bottleneck
in this particular run. They are not a proof that the required relations
do not exist, or a reliable estimate of the time needed to obtain them.
No larger resource allocation or parameter retry is part of this run.

## Exact compact evaluation control

For a compact product u=product(a_j^e_j), write pi(a)=N(a)*a. Norm
projection is multiplicative, and squareclass characters are linear:

\[
 \lambda_v(\pi(u))=
 \sum_j(e_j\bmod2)\lambda_v(\pi(a_j)).
\tag{1}
\]

Thus large exponents need not be expanded to evaluate the projected
S-unit's local squareclass. The individual factors need not themselves
be S-units: the complete compact expression certifies that property,
and (1) evaluates its squareclass.

On the small control the verifier checks all six supplied S-unit
generators by exact ideal valuations. It independently expands these
small examples to verify (1), and checks local zero/nonzero status with
PARI's separate local-power interface. The ordinary class group is C2;
inverting the four prime ideals above 2,3,11 kills its two-primary part.

Dirichlet's theorem and the odd-degree norm map give

    dim E_S = (1+1-1)+4-3 = 2.

The compactly evaluated norm projections have local rank **two**, so
they span the full norm-square S-unit space and its strict kernel is
**zero**. This matches the independently computed S-class two-rank zero.
It is a calibration of the local correction machinery, not an additional
high-rank or class-creation example.

The [verification artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_compact_sunit_basis_verification_v1.json)
contains the factor expressions, S valuations, projected local bits,
dimension checks, and the audited research terminal states. The first
verifier attempt used an unavailable PARI accessor; switching to the
documented BNF nf component fixed that interface before a verifier
artifact was produced.

## What remains necessary

1. Obtain a sufficiently large **equation-derived** S-unit squareclass
   span on a research fibre, with exact support witnesses. A timeout,
   a partial class-group rank, or an assumed unit rank does not supply it.
2. Use (1) and the generic local columns to solve the joint strict
   correction equations. Remove global squareclass dependencies before
   asserting a dimension gain; a coefficient-kernel dimension alone
   can overcount dependent classes.
3. Certify independence of the resulting strict classes modulo G at
   auxiliary good primes, then study CT and rational solubility. This
   is where the prospective arithmetic must remain independent of
   exceptional-point characters.

There is no new incidence or visibility feature for Agent 1. The
generic-ideal criterion remains valid; its explicit correction input is
still missing. The existing isolated-prime tests and the present general
class-engine attempt have not supplied the mixed principal relations
needed to construct it.

Only the small control and saved terminal evidence are replayed by:

```sh
timeout 30 sage -python elliptic-curves/rank-jump/verify_compact_sunit_basis.py check
```

Active searches, their worker limits, candidate populations, and
mathematical-status entries remain unchanged.
