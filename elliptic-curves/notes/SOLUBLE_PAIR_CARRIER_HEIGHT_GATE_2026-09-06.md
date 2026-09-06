# Height gate for the newly soluble pair carrier

The recent [global-carrier calculation](../rank-jump/GLOBAL_CARRIER_SOLUBILITY_AND_SPECIALIZATION.md)
proves that the pair `orbit-1795d`, `orbit-0911e` has a soluble actual
genus-one carrier and a Jacobian of exact rank2. That gives a concrete source
of parameters with two rational bisection lifts. It does not ensure those
parameters give practical search models or that the two specialized points
remain independent modulo the generic subgroup.

A separate fixed audit now maps twelve specified combinations of the two
recorded Jacobian generators onto that actual carrier. Reversing its quartic
makes the retained rational point at infinity an explicit affine origin. The
pointed Weierstrass model, rational inverse, two square conditions and both
generic lift polynomial identities are checked exactly. Each finite image is
transported to compact08234 and its17 generic sections plus the two constructed
points receive an exact finite-reduction independence check through997.
No original-family point search or score calculation occurs.

The [result](../../artifacts/generated-results/elliptic-curves/soluble_pair_carrier_height_v1.json)
and its replay pass in6.385 and6.295 seconds. Of twelve fixed words:

- Ten give certified lower bounds19. Their recorded homogeneous short models
  have640–5044 coefficient bits. These are seven distinct parameter images;
  three repeats are retained explicitly.
- The word `[1,1]` returns the existing anchor774/149, with a256-bit model.
  The19 supplied points certify only18 modulo2; a separate
  [modulo3/5 audit](../../artifacts/generated-results/elliptic-curves/soluble_pair_carrier_anchor_modl_v1.json)
  also gives18. Those finite checks alone did not determine the exact span;
  the subsequent rational relation below closes that question.
- `[0,-2]` exceeds the declared512-bit compact-parameter cap and is censored
  before elliptic specialization. Its parameter has596 bits.

Thus none of these fixed images supplies a new model under the declared360-bit
search-model gate. Coefficient sizes refer to the recorded homogeneous models,
not globally minimal models. The bound19 is an independent-point certificate,
not a whole-curve rank. These are a deterministic subgroup sample and a bounded
cost gate, not an exclusion of small parameters elsewhere on the carrier.

## The height failure survives integral normalization

A subsequent [factor-free invariant certificate](../../artifacts/generated-results/elliptic-curves/soluble_pair_carrier_invariant_height_v1.json)
closes the normalization gap for these fixed images. All ten nonanchor rows,
representing seven distinct parameters, admit **no normalized integral model
within360 coefficient bits**, including a normalized global minimal model.
The smallest certified necessary coefficient size is610 bits. The
[independent Sage replay](../../artifacts/generated-results/elliptic-curves/soluble_pair_carrier_invariant_height_replay_v1.json)
checks all eleven available invariant pairs and exact inequalities.

For integral source invariants letG=gcd(|c4|,|c6|). Under a rational
isomorphism of scaleu=a/b in lowest terms, the target invariants are
c4/u^4 andc6/u^6. Their integrality impliesa^4 dividesG, so

    |c6_target|^2 >= c6^2/G^3.

For a normalized integral equation witha1,a3 in{0,1}, a2 in{-1,0,1},
and|a4|,|a6|<=M, the invariant formulas give

    |b2|<=5, |b4|<=2M+1, |b6|<=4M+1,
    |c6_target|<=1224M+521.

Every nonanchor row satisfiesc6^2>G^3(1224M+521)^2 forM=2^360-1.
This excludes all such normalized equations under arbitrary rational
scaling; integral short models within the limit are excluded too.
It does not assert a bound for arbitrary unnormalized equations with all
five coefficients bounded. The anchor remains unexcluded, and the twelfth
word remains censored before specialization. No factorization, explicit
global minimal model or new point search is needed for this proof.

Replay `certify_carrier_invariant_height_gate.py --check` and
`verify_carrier_invariant_height_gate.sage --check` under `../cas/`.

The newer [quartet geometry](../rank-jump/SOLUBLE_QUARTETS_REQUIRE_HIGHER_GENUS_LIFTS.md)
also distinguishes these two supplied directions from the whole observed
jump: two further fixed square conditions require a degree-four lift from
this elliptic carrier to a genus17 curve. The carrier's infinitely many
rational points do not supply infinitely many such simultaneous lifts.

## Generator coverage

The anchor's [exact relation certificate](../../artifacts/generated-results/elliptic-curves/soluble_pair_carrier_anchor_relation_v1.json)
now proves that its supplied19 points span exactly18 dimensions. Writing
P1,...,P17 for the generic sections in the certificate's order, andQ1,Q2 for
the two supplied carrier lifts, the exact rational group law gives

    Q1+Q2 = P1-4P2-6P4+4P5+4P7+P8-2P10+2P12+P13-3P15+P17.

The17 generic points andQ1 are independently certified. Hence the pair's
quotient over the generic span has **exact rank1** on this anchor. A
384-bit height calculation proposed this bounded word; the Sage-free
verifier checks the rational group identity and both independence proofs.
The whole curve's rank and other carrier fibres are unchanged.
Replay `../cas/certify_carrier_anchor_relation.py --check`.

A subsequent fixed [generator audit](../../artifacts/generated-results/elliptic-curves/soluble_pair_generator_saturation_v1.json)
proves that the recorded rank2 subgroup together with rational2-torsion is
saturated at2,3,5. At2 the two free columns and torsion column have finite rank3;
at3 and5 the free columns have finite rank2 and separate reductions exclude
rational ell-torsion. The argument includes the rational2-torsion explicitly.
Build and replay pass in0.640 and0.744 seconds. Thus missing2,3,5 divisions do
not enlarge this particular carrier subgroup. Larger-prime saturation and a
full integral basis remain open.

The degree-four parameter map also identifies some distinct Jacobian words,
so counting words as distinct candidate curves would overstate coverage.
Further carrier sampling or group completion would need its own finite
protocol; the current audit supplies no reason for automatic expansion.

Source: `../cas/audit_soluble_pair_carrier_height.sage`. The frozen protocol,
raw logs, source hashes and the separate19-point anchor audit are retained in
`artifacts/local/elliptic-curves/soluble-pair-carrier-height-v1`.
This gate does not change the full11952 cohort or its point budgets.
