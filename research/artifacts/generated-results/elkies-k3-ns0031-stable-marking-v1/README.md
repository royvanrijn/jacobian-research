# NS0031 full stable arithmetic obstruction

Canonical proof: [full stable marking obstruction](../../../elkies-k3/NS0031_FULL_STABLE_MARKING_OBSTRUCTION_2026-09-14.md).

- [certificate.json](certificate.json): exact stable group, arithmetic level,
  genus-10 curve, degree-four genus-2 quotient, and its twelve rational points.
- [magma-output.xml](magma-output.xml): unchanged response from Magma 2.29-10.
- [magma-output.txt](magma-output.txt): readable extraction of that response.
- [magma-replay.json](magma-replay.json): input/output hashes and execution receipt.

The [Sage producer/checker](../../../elkies-k3/scripts/certify_ns0031_stable_marking.sage)
recomputes the finite group, modular symbols and Sturm-bound calculations;
`--check` compares the certificate and audits the retained Magma proof.
The [Magma input](../../../elkies-k3/scripts/certify_ns0031_genus2_rational_points.m)
separately replays the unconditional full Mordell–Weil and elliptic-Chabauty
proof. The raw run used the public Magma calculator, with a 60-second CPU
limit, and completed in under six CPU seconds.

The old norm-one obstruction and missing-reflection counter-witness remain
separate retained artifacts. This packet closes the corrected problem by a
new quotient argument; it does not rehabilitate the old group containment.
