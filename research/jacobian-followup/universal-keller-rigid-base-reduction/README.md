# Rigid bases and the quartic canonical obstruction

The fixed universally quartic map remains strongly absolutely stably atomic, but its quotient W has metabelian automorphism group and cannot support the required canonical S4 action. This obstruction survives every finite coupled affine cylinder.

- [Written structural proofs](PROOF.md): cylinder rigidity, automatic finiteness of dominant torus-quotient maps, the exact quartic automorphism group, classification-free Sylow degree budgets, and reduction of nonfinite attempts to affine fibres.
- [Literature and remaining gates](LITERATURE_AND_GAPS.md).
- [Exact controls](verify.py), [receipt](receipt.json), [source pins](sources.json), [ledger integration](integration.json), [focused validation](validation.json).

Replay: `python verify.py` with [requirements](requirements.txt). The controls check matrices and finite groups; they are not an independent proof replay. The original affine canonical descent problem remains open.

Further extension: [all-degree holonomy families](../universal-keller-holonomy-families/PROOF.md) supply prime-square strong atoms, the exact cubic norm-torsor obstruction, complete Z_m factor lengths and a canonical quartic exclusion even when the enlarged base contains S4 automorphisms.

The [atomic-cylinder extension](../universal-keller-atomic-cylinder/PROOF.md) removes the metabelian hypothesis for nonfinite self-atomic data on every torus-quotient base. Its relative restriction yields affine-space forms as a pair, without chosen coordinates or an identified self-map.
