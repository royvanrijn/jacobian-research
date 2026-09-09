# No hidden unit relation in the sealed302 multiplier bank

## Result and boundary

**New deduction, independently checked application.** Let `a_1,...,a_567`
be exactly the multipliers saved by the old generic-only302 ideal-reduction
experiment, in the certified maximal order of its cubic field `K`. Then

\[
 \boxed{\langle\mathbf Q^*,a_1,\ldots,a_{567}\rangle
          \cap\mathcal O_K^*=\{1,-1\}.}
\]

Exponents are arbitrary integers. This excludes longer multiplicative
relations, not just the81 repeated-ideal pairs tested previously. Neither
the bank nor its reduction directions was enlarged or rerun.

**Verified application.** The567 saved identities reduce to486 distinct
primitive ideals. Of these,479 have successive private unramified rational
prime witnesses; their multipliers are multiplicatively independent modulo
rational scalars and units. The other88 trials have multiplier exactly1.
There are no primitive-ideal repetitions across distinct starting ideals.

**Boundary.** This explains why further relation processing of this bank
cannot produce a unit. It does not exclude other unit algorithms, prove
either seed-relative ideal nonprincipal, or decide seed incidence. The
[two explicit integral norm equations](DET1092_SEED_VIRTUAL_UNIT_AND_IDEAL_PARITY_2026-09-09.md)
remain locally soluble everywhere and globally unresolved. No exceptional
point is an input to this computation.

## Private-prime lemma

**New deduction from elementary ideal factorization.** Let `J` be an
integral ideal primitive with respect to rational scaling: its HNF entries
in an integral basis have gcd1. Suppose an unramified rational prime `p`
divides `N(J)`, but divides neither any starting-ideal norm nor the norms
of the other currently retained primitive ideals.

Some prime above `p` has positive `J`-valuation because `p | N(J)`.
Some prime above `p` has valuation zero: otherwise `J` would be contained
in `p O_K`, contradicting primitiveness. Unramifiedness is essential;
the valuation vector of a rational principal ideal is constant at the
primes above `p`.

Consequently a relation in which `J^n` times the other retained ideals is
rational principal forces `n=0`, by comparing its positive and zero
components above `p`. Remove this ideal and repeat. This argument also
works for negative exponents and has no search-height bound. **Without the
unramified condition, a totally ramified prime would invalidate this
argument.**

## Exact application to the bank

**Verified application.** Each retained record has

\[
 J_i(a_i)=I_{s(i)},\qquad J_i=c_i J_i^{\rm prim},\quad c_i\in\mathbf Z_{>0}.
\]

For repeated primitive ideals from the same source the replay checks
`c_i a_i/(c_j a_j) = +/-1`. Thus repetitions introduce only rational
scalars and signs. Their exponents can be grouped before applying the lemma.

To certify a private prime without integer factorization, put

\[
 B=|\operatorname{disc}K|\prod_s N(I_s)
                  \prod_{j\ne i\text{ still retained}}N(J_j^{\rm prim}).
\]

Starting from `N(J_i^prim)`, repeatedly divide by its gcd with `B`.
The final remainder `d_i>1` is coprime to `B`, so any prime dividing it
is a private unramified prime. The certificate retains every gcd removal
and every remainder. The independent checker verifies the divisions and
coprimality directly; it does not run the selection algorithm.

There are479 successful eliminations. The six remaining nonunit ideals
are starting reductions, and the remaining unit ideal is the identity
source. All88 trials belonging to these seven groups have `a_i=1`.
Applying the lemma successively therefore makes any surviving product a
rational scalar. A rational scalar is an algebraic unit exactly when it
is `+/-1`. Conversely those two units are already in `Q*`. This proves
the boxed equality.

**Interpretation, not a new incidence claim.** The old weighted reductions
generated independent ideal-divisor changes, not unit cycles. Their small
reduced ideals did not supply the multiplicative cancellations needed for
a unit. This is an exact reason for that particular negative experiment,
not an explanation of the absence or existence of rational points on a
control fibre.

## Reproduction and retained evidence

**Verified application.** The support construction completed in0.57s and
the independent replay in0.63s, each capped at25s. An initial API failure
was retained; it produced no mathematical conclusion. No factorization,
new reduction, class group, unit group, point search, parameter sweep,
production change, or detached job occurred.

```bash
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/audit_det1092_saved_unit_relations.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_saved_unit_relations.sage
```

- [Frozen protocol and input hashes](../../artifacts/generated-results/elliptic-curves/det1092_saved_unit_relations_v1/protocol.json).
- [479-step factor-free support certificate](../../artifacts/generated-results/elliptic-curves/det1092_saved_unit_relations_v1/support.json).
- [Independent replay](../../artifacts/generated-results/elliptic-curves/det1092_saved_unit_relations_v1/independent-replay.json).
- [Original experiment](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_arithmetic_unit_class_v1.json) and [sealed raw-evidence archive](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_transfer_attempts_evidence_v1.zip).

The local source files referenced by the old experiment are hash-bound in
the new protocol; the sealed archive preserves their historical evidence.
The status authority records this as a theorem about one fixed input bank,
not the solution of the seed-construction objective.
