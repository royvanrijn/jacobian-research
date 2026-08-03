# Conditional simple-spectator gluing for the F2 `A_6` packet

## Result and claim boundary

The terminal F2 residue cover has the certified branch-cycle passport

\[
(5,1)\mid(3,3)\mid(3,1,1,1),
\qquad
\sigma_0\sigma_\infty\sigma_*=1,
\]

and geometric monodromy `A_6`.  On the squarefree cofactor stratum there are
also two simple `R` Kummer orbits.  The existing Kummer theorem proves only
that these are not further above-bisectrix F2 continuations.  It does **not**
give their residue degrees, target components, or branch cycles.

This note makes one narrow bridge assumption and exhausts the resulting
permutation problem:

> each simple `R` orbit gives one separate simple branch value, hence one
> transposition; the terminal triple fixes all sheets outside its six-sheet
> packet; and these are the only five branch cycles on one compact target
> `P^1`.

If a spectator sheet is required, the abstract cover is then forced to have
degree seven.  After fixing the six core letters and the terminal five-cycle,
there are exactly `30` normalized factorizations and `6`
simultaneous-conjugacy classes.  Equivalently, there are `30,240` fully
sheet-labelled tuples.  Every one has genus zero and generates `S_7`.  Thus
this conditional model gives six tiny residual permutation gluings, not a
contradiction.  Retaining the certified endpoint/interior markings refines
the census: five classes can specialize the connector anchor to the
five-cycle point `s=-1`, all six can specialize it to one of the two
denominator-root attachments, and the strongest naive requirement that the
anchor avoid both toric source endpoints leaves exactly three classes of
signature `(5,3,1)`.  Marked incidence therefore prunes but still does not
contradict the conditional model.

The exact checker and pinned output are
[`cas/verify_f2_a6_simple_spectator_gluing.py`](cas/verify_f2_a6_simple_spectator_gluing.py)
and
[`../artifacts/generated-results/jc2_f2_a6_simple_spectator_gluing.json`](../artifacts/generated-results/jc2_f2_a6_simple_spectator_gluing.json).
The artifact SHA-256 is
`34ed6f4711bef592c87f86e7e4e19f70984a35259fab68dec6c35c71cd9b80a2`;
the checker uses only the Python standard library.  Intentional regeneration
uses the reproduction command below with `--refresh`.

The result is a **conditional finite enumeration**.  It is not a
classification of the actual F2 spectator orbits and does not exclude
`(75,125)`.

## 1. Certified input

The extracted terminal row has source ray `(12,-17)`, target ray `(5,2)`,
transverse index `e=1`, and residue degree `f=6`.  It is centered at target
infinity, has zero transverse different, and has residue-different packet
`(4,2,2,2)`.  The permutation audit uses the residue cover only; the companion
terminal checker verifies these toroidal statements directly.

Fix one of the five terminal triples on six letters.  Its ramification
contributions are

\[
6-2=4,
\qquad
6-2=4,
\qquad
6-4=2.
\]

Their sum is

\[
4+4+2=10=2\cdot6-2.
\]

The checker independently reproduces the five triples with a fixed
five-cycle, verifies their product is one, and verifies that each generates
the `360` elements of `A_6`.  This part replays the certified theorem in
[`F2_TERMINAL_RESIDUE_COVER.md`](F2_TERMINAL_RESIDUE_COVER.md).

The squarefree cofactor row has two simple Kummer orbits by
[`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md).  No statement
in either source identifies one such orbit with a transposition.  That
identification is an assumption from this point onward.

## 2. Two simple cycles force the candidate degree

Embed the terminal triple in `S_N` by fixing the other `N-6` letters, and
append transpositions `tau_1,tau_2`.  Because the terminal product is already
one, the compact meridian relation gives

\[
\tau_1\tau_2=1.
\]

Transpositions are involutions, so

\[
\boxed{\tau_1=\tau_2.} \tag{1}
\]

The total ramification is `10+1+1=12`.  Riemann--Hurwitz for a connected
cover of `P^1` gives

\[
2g-2=-2N+12,
\qquad
g=7-N. \tag{2}
\]

If there is a genuine sheet outside the terminal packet, then `N>=7`.
Nonnegative genus in (2) gives `N<=7`.  Hence

\[
\boxed{N=7,\qquad g=0.} \tag{3}
\]

Equation (1) also proves that both simple branch values use the same
connector.  Transitivity forces it to be

\[
\tau_1=\tau_2=(a\ 7)
\]

for some terminal sheet `a` in `{1,...,6}`.  A transposition supported
strictly away from the terminal block leaves the `A_6` orbit disconnected
and cannot be a global connector.

For comparison, the exact degree audit is:

| degree | product-one normalized tuples | transitive normalized tuples | RH genus |
| ---: | ---: | ---: | ---: |
| 6 | 75 | 75 | 1 |
| 7 | 105 | 30 | 0 |
| 8 | 140 | 0 | impossible |

The degree-six row adds two internal branch points but no spectator sheet.
For degree at least eight, equation (2) already gives negative genus; the
repeated transposition also cannot connect more than one new sheet.

## 3. The six exact residual gluings

For each of the five terminal triples with the core and five-cycle fixed, and
each of the six choices of `a`, the checker gets one normalized
factorization, for `5*6=30` total.  It then quotients by simultaneous
conjugacy in `S_7`, while retaining the five branch-value labels.  Exactly
six classes remain.  Each tuple generates `S_7`, whose center is trivial, so
each simultaneous-conjugacy orbit has `7!=5,040` sheet labellings.  The six
classes therefore give `6*5,040=30,240` fully labelled tuples.

The terminal cover has trivial target-fixed deck group.  Correspondingly,
its six possible connector sheets are not identified after the terminal
triple is fixed.  The cycle lengths containing the connector sheet at the
three terminal branch values give the following coarse census:

| connector signature | conjugacy classes |
| --- | ---: |
| `(1,3,3)` | 1 |
| `(5,3,1)` | 3 |
| `(5,3,3)` | 2 |

The signature is not complete: the three middle classes and the two last
classes remain distinct.  The pinned JSON records an exact branch-cycle
representative for each of the six classes.

Every completed tuple generates `S_7`.  Indeed `A_6` is transitive on the
six terminal letters.  Conjugating `(a 7)` by its elements gives every
`(i 7)`, and those transpositions generate `S_7`.

Consequently product one, transitivity, and Riemann--Hurwitz do not exclude
the conditional squarefree row.  They reduce it to six marked permutation
gluings.

### 3.1 Endpoint/interior marking filter

The terminal cover identifies more than three unmarked cycle types:

- over `0`, the fixed point is the source endpoint `s=0`, while the
  five-cycle is the interior attachment `s=-1`;
- over `infinity`, the two three-cycles are the two interior
  denominator-root attachments; and
- over `125/729`, the three-cycle is the source endpoint `s=infinity`,
  while the three fixed points are the simple interior roots of the cubic.

For a connector anchor `a`, its cycle lengths in the terminal triple record
which of these points it specializes to.  Applying these labels to the six
classes gives:

| conditional marked requirement | surviving classes |
| --- | ---: |
| anchor specializes to the interior point over `0` | 5 |
| anchor specializes to an interior point over `infinity` | 6 |
| anchor avoids the source endpoint over `125/729` | 3 |
| anchor is interior over `0` and avoids both source endpoints | 3 |

The last three classes all have signature `(5,3,1)`.  This is the strongest
filter available from the marked residue triple alone, and it still
survives.

There is a geometric caveat.  Specialization of a generic sheet label to a
marked point of the residue cover is not itself a proof that a spectator
boundary component glues at that source node.  Promoting the last table to
an F2 statement requires a toroidal node-compatibility theorem identifying
the spectator valuation with one of the three forced boundary attachments.
The checker therefore records this as a conditional marked filter, not as a
new bridge theorem.

## 4. Why a larger remaining degree also survives coarse constraints

The two-cycle hypothesis is essential.  For every `k>=1`, append a repeated
star connector for every new sheet:

\[
(1\ 7),(1\ 7),
(1\ 8),(1\ 8),
\ldots,
(1\ 6{+}k),(1\ 6{+}k). \tag{4}
\]

Together with the terminal triple, (4) has:

\[
N=6+k,
\qquad
\prod\sigma_i=1,
\qquad
R=10+2k=2N-2. \tag{5}
\]

It is transitive, has genus zero, and generates `S_(6+k)`.  The checker
constructs and verifies these witnesses for `k=1,...,19`; the displayed
argument proves them for every `k`.

In particular, if one *incorrectly* identified the common-edge polynomial
degree `25` with global geometric degree, a genus-zero simple completion
would require `38` transpositions.  Two assumed orbit-level transpositions
would leave `36` additional cycles to locate.  That is missing data, not a
contradiction.  The canonical handoff explicitly does not identify common-
edge degree with field degree.

There is an analogous abstract witness for two disjoint terminal packets on
twelve letters.  Their core ramification is `20`; one repeated cross-packet
transposition adds `2`, gives `22=2*12-2`, makes the action transitive, and
generates `S_12`.  The double-root F2 stratum itself supplies no certified
simple cycles, so this only shows that two connectors from some other global
row would evade the coarse permutation filters.

## 5. Missing bridge to F2

The enumeration cannot be promoted to an F2 theorem until one proves all of
the following geometric statements:

1. a simple `R` Kummer orbit determines a branch value at all;
2. its inertia is one transposition rather than an unramified residue row or
   a different cycle profile;
3. it lies over the same target component as the terminal packet;
4. the terminal and spectator data assemble into the connected compact
   curve cover used above;
5. the purity-forced affine ramification row and every other boundary row
   contribute no further cycles; and
6. the remaining sheet degree is known independently.

The actual surface valuation equality may have several source divisors over
one target divisor, and its decomposition action can have several orbits.
Global monodromy transitivity does not by itself turn their separate residue
covers into the single `P^1` cover assumed here.  This is the main reason the
permutation audit is conditional rather than an exclusion.

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py
# Intentional refresh only after reviewing a changed enumeration:
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py --refresh
```

Expected final marker:

```text
F2_A6_SIMPLE_SPECTATOR_GLUING_CONDITIONAL_PASS
```
