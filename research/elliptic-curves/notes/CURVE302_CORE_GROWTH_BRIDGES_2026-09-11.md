# Curve302 decisive common-core bridge analysis

This is a positive-evidence-only post-processing experiment over the sealed
Curve302 chart replay.

It restricts attention to historical gain stages where at least one explicitly
recorded positive quotient direction gives a strictly positive exact one-step
intersection-rank gain with the next observed common integral core.  In the
current sealed replay this is expected to be 27 stages.

For each decisive stage it records:

- current and next common-core intersection ranks;
- actual gain word(s) and every explicitly positive alternative;
- exact one-step common-core gain after `Sat_Z(S + Z v)`;
- whether the historical gain is uniquely best or tied by positive alternatives;
- a canonical quotient line for `v mod span_Q(S)`;
- which of the 14 named quotient axes become newly contained after saturation;
- when the one-step core gain is one, the canonical new line inside the target
  common core modulo its current intersection.

The optional raw-schema enrichment reads the original chart definitions only for
charts exposing actual or tied-best bridge directions in these decisive stages.
It retains exact `centre`, `mapping`, and `search` payloads and reports both
schema-only and exact-value hashes.  It performs no point search and no point
recognition.

Missing chart hits are never treated as negative evidence.  The analysis does
not depend on chart-completeness flags.
