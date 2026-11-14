# Complete quintic cusp locus and the excess-one singleton exclusion

14 November 2026. PF2R7B5CUSP2 classifies all nonimmersive high-contact
degree(3,5) candidates:161 normalized(C,B) pairs over160 first-image
coefficient points. The159 nonspecial pairs have one ordinary cusp and
three nodes; the two special pairs have one cusp and one ordinary triple.
Every cusp is clean, lies off C, and has infinity contact38.

PF2R7B5B1 therefore excludes singleton fibres throughout the immersed
excess-one row B in this chart. Together with the previous nonimmersed
row-N exclusion, only immersed rows A or C can explain a singleton
under these degree/contact hypotheses.

Actual B is not forced to have degrees(3,5) or contact at least38.
The global Keller covering and source realization remain open.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/jc2-r7-quintic-cusps/verify.py
```

[Proof](PROOF.md) · [classification and modular certificates](classification.json) ·
[replay](receipt.json) · [validation](validation.json) ·
[author-side audit](REVIEW.md) · [discovery](discovery.json).

Authored finite-algebra and fibre interpretation with exact certificates;
no formalization, independent external review or JC2 closure. Having no
singleton is not itself a contradiction to a Keller map.
