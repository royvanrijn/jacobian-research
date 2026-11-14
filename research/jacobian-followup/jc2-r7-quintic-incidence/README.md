# Quintic incidence: a finite cusp gate and a nonimmersed singleton exclusion

14 November 2026. Three conditional results refine the seven high-contact
degree(3,5) candidate images:

* Generically all seven have four nodes and twelve transverse intersections
  with C, with mutual singular-locus avoidance. Local branch counts allow
  row C and both shared-image row-A groupings; they do not exclude them.
* The excess-one row B must lie over160 normalized first-image coefficient
  points, at most54 cube-root orbits, given by two explicit polynomials.
* At both nonimmersed C classes, every candidate avoids the known cusp.
  The oriented row N therefore has no singleton fibre in this chart.

Actual B is not forced to have degrees(3,5) or infinity contact at least38.
The global Keller covering and source realization remain open.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/jc2-r7-quintic-incidence/verify.py
```

[Proof](PROOF.md) · [matrices, norms and witnesses](incidence.json) ·
[replay](receipt.json) · [validation](validation.json) ·
[author-side audit](REVIEW.md) · [discovery and failed runs](discovery.json).

Written finite-algebra and geometric interpretation with exact modular
ranks and rational determinants; no formalization, independent external
review, global Keller construction or JC2 closure.

**Finite cusp gate completed:** [PF2R7B5CUSP2/PF2R7B5B1](../jc2-r7-quintic-cusps/README.md)
classify161 nonimmersive image pairs over the160 coefficient points.
Every cusp is ordinary, clean and off C. The immersed excess-one rowB
therefore has no singleton fibre under the same degree/contact hypothesis.
Together with rowN, this leaves immersed rows A/C in that chart.
