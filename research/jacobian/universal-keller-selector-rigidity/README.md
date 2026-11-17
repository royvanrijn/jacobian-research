# Arbitrary source selectors still preserve one-direction factors

The common first-integral gate is now computed in every polynomial degree. A Keller selector can use an active feedback block only through its retained first output; it may use all source coordinates of an inactive block. The selector is independent of the auxiliary source.

- [General proofs](PROOF.md): exact rational/polynomial common kernels, arbitrary-selector classification, and polynomial processing of active blocks, auxiliary map, then inactive blocks.
- For the auxiliary cubic in any constant nonzero direction, EVERY allowed source selector still factors. A displayed atomic chain has length N+1; no classification of all chains is claimed.
- [Literature](LITERATURE_AND_GAPS.md), [checker](verify.py), [receipt](receipt.json), [sources](sources.json), [integration](integration.json), [focused validation](validation.json).

Replay `python verify.py` with [requirements](requirements.txt). Canonical absolute atomicity remains open outside this architecture. No independent whole-proof replay, formal verification, external review or priority certificate.

The [joint-selector extension](../universal-keller-joint-feedback/PROOF.md) now classifies ALL polynomial beta(P,U) via a Bertini field argument. It factors every cubic member and every auxiliary with polynomial or stably polynomial full first-integral ring. Joint dependence alone is no longer an unclassified atomicity gate.
