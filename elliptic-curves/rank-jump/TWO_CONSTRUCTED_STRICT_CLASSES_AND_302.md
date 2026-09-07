# Two constructed strict classes; 302 remains the primary target

The constructive calibration has succeeded. From the MW16-05 equation at
`t=3/17`, its sixteen generic Kummer classes, and equation-only principal
relations, we have constructed **two independent strict Selmer classes
outside the generic subgroup**. The independent certificate proves combined
mod-2 rank **18**. Exceptional points and historical point-derived class
anchors were not construction or verification inputs.

This supersedes the unmet calibration endpoint in the
[previous construction report](POSITIVE_CLASS_CONSTRUCTION_AND_THE_302_PIVOT.md).
It does **not** yet construct a new class on302, prove rational solubility
of the two classes, or identify a specialization condition creating them.

## Exact objects and certificate

Let `K=Q[theta]/(f)`, where

```
f(z) = z^3 + z^2
       - 2919231625641258502793755607986240*z
       + 45440201616242830029801770634418828098464819545088.
```

Write `pi(alpha)=Norm(alpha)*alpha`. Its norm is a fourth power, and
changing alpha by a rational scalar changes pi(alpha) by a fourth power.
The [independent certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_reference_additional_strict_verification_v1.json)
contains both exact factored representatives

\[
 \beta_j=\prod_{i\in I_j}\pi(\alpha_i)
          \prod_{k\in J_j}\gamma_k,\qquad j=1,2,
\]

where the gamma are the frozen generic representatives. No unspecified
square root or class-group generator occurs in this definition. The
certificate lists every alpha, factor index, generic correction and prime
ideal. The two products use1676 and1572 principal atoms, respectively, with
generic corrections

```
J1 = [3,5,6,7,9,10]
J2 = [0,1,4,5,6,7,10]          (zero-based generic indices).
```

In the elimination ordering, the closing atoms are

```
alpha_4717 = -45363542651566580064 - 741*theta
alpha_4718 = -155295826912877673888 - 2547*theta.
```

These small closing atoms are not themselves the new strict classes. Their
large accompanying principal-ideal dependencies and local generic
corrections are essential.

The [independent verifier](verify_reference_additional_strict.py) checks:

- Exact norms and42786 prime-ideal valuations on the2411 distinct used atoms.
- Even projected valuations at every prime outside the full bad set S.
- At every prime above S, zero Hilbert pairings against a full local
  squareclass basis. The Hilbert Gram ranks certify completeness of each
  local test basis. This is independent of the producer's ideal-log tests.
- Positive signs at all three real embeddings.
- An18-by-18 good-prime character minor, recomputed from the factors, proving
  independence of the sixteen generic classes and both additions.

Thus this is an unconditional certificate of two additional strict Selmer
classes. It is not a class-group or Selmer upper bound, and it does not use
the reference's GRH-conditional complete descent.

## Which construction step worked

The first adaptive relation run completed501 targets in1073.61seconds.
Its formal supported quotient fell from350 to73; its final strict extractor
still detected only the six inherited strict classes.

We then constructed ten ordinary unramified characters from the generic
sections alone. At odd primes these impose even valuations; at two they
impose Hilbert orthogonality to local units; at real places they impose
positivity. Independently detected small-prime ideal classes provide
protected inherited anchors. The same procedure gives six such generic
characters on302. These character counts are **inherited controls**, not
new directions or estimates of class-group rank.

The new selector targets small-prime representatives of unresolved quotient
combinations, including columns that were already pivots. It also adds only
the sixteen generic principal parity rows. No historical exceptional-point
anchor is read. In352 completed targets, with strip bound1024, the formal
quotient fell from70 to18. No target class-group dimension is a stopping
rule: the endpoint is an additional strict independence witness.

At that checkpoint the outside-S parity kernel acquired two large new
dependencies. The old verification window below1009 then failed to detect
even generic rank16: requiring every factor in the enlarged dictionary to
be a unit at a test prime had eliminated too many primes. The run stopped
fail-closed. The preserved snapshot was re-extracted with fixed primes
50001 through55000, above the relation smoothness bound50000. No new norm
relations were searched. This recovered detected strict rank8, and the
separate Hilbert/ideal verifier certified both additions.

The methodological lesson is constructive: use small representatives of
unresolved ideal combinations, with independently justified generic anchors,
and evaluate large factored classes away from their factor support. A raw
free-column count was a poor guide to the cost of finding the missing
principal relations. None of these runtime diagnostics is a rank predictor.

## Solubility is now a concrete second-stage question

Each explicit beta determines a2-cover without supplying a point. Write

\[
 \beta(u+v\theta+w\theta^2)^2
   =Q_0(u,v,w)+Q_1(u,v,w)\theta+Q_2(u,v,w)\theta^2.
\]

Its standard genus-one degree-four carrier is the intersection in P3

\[
 C_\beta:\quad Q_2=0,\qquad Q_1+s^2=0.
\]

A rational point with s nonzero maps to the cubic model by
`x=Q0/s^2` and
`y=sqrt(Norm(beta))*Norm(u+v*theta+w*theta^2)/s^3`.
The square root of the norm is supplied by the factored construction.
These equations specify the covers exactly from the certificate; large
expanded quadratic coefficient tables have not been generated.

Strictness proves local solubility. Global rational solubility is **UNKNOWN**.
A nonzero Cassels--Tate pairing against a certified Selmer class would prove
that a candidate is not rational. A zero pairing does not prove a rational
point. No common smaller carrier for both classes has been established.

## What this transfers to the new 302 parent

The [302 arithmetic certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_strict_constructor_arithmetic_v1.json)
already supplies the complete bad set, maximal order, ramification/index
data, seventeen generic Kummer classes and complete local boundary.
Its generic strict dimension is zero. Thus a nonzero strict class on302
would automatically lie outside the marked generic subgroup. The total
local point-Kummer dimension is22; the joint generic local image has
dimension17, even though each individual local projection is full.

The saturated new parent gives
`D/sp(M17)=Z^14` in the certified independent rank31 subgroup. The whole
rational group and exact rank of302 remain unknown. The parent's
retrospective construction provenance is retained: freezing its generic
sections does not make it a prospectively selected family.

The calibration now supplies a verified algorithmic path to an individual
class without requiring full class-group completion. The remaining302 work
is to obtain an appropriate fibre-specific principal dependency in its
larger cubic field, then apply the same norm projection, local correction
and independent-character certificate. Copying the smaller field's
parameter or smoothness budgets would not itself justify that step.

The exact specialization event is still unidentified. We have proved

\[
 \text{explicit principal-ideal dependency at }t_0
 \Longrightarrow \text{two extra strict classes on the reference},
\]

but not either of the missing implications

\[
 \text{interpretable condition on }t
 \Longrightarrow \text{that dependency},\qquad
 \text{strict classes}\Longrightarrow\text{rational soluble block}.
\]

In particular, discovery of a principal relation is not evidence that a
new ramified prime caused it. The next mechanism investigation should track
these explicit ideal/cover identities in the marked family, rather than
returning to generic capacity exclusions or unlabeled score correlations.

## Replay and ownership

The [checkpoint archive](../../artifacts/generated-results/elliptic-curves/rank_jump_positive_class_evidence_v1.zip)
and [manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_positive_class_evidence_v1.json)
preserve the completed runs, the failed small-prime verification and its
corrected extraction. They supplement the previous constructive archive.
Use an empty replay checkout when restoring them. Producers refuse to
overwrite evidence. The immutable extractor is
[reference_large_support_extraction.py](reference_large_support_extraction.py);
the independent verifier is linked above. Both have explicit bounded
protocols beside them.

Only rank-jump scripts, notes and artifacts were added or edited. No live
search protocol, candidate output, worker setting or mathematical-status
entry was changed by this work.
