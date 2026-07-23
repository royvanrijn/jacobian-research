# Hurwitz--LL compactification of rerooting

This paper contains the admissible-cover closure, rerooting quotient stack,
marked-zero-fiber LL degree, formal collision comparison, and the boundary
Picard and caustic/Maxwell calculations.

The missing global comparison is resolved negatively as stated by the
[global comparison obstruction](global-comparison-obstruction.tex).  In
degree five, the stable zero-root compactification forgets the target
cross-ratio with leading term `x^3/y^2`; the admissible-cover model resolves
it by the normalized blowup of `(x^3,y^2)`.  Thus there is no global morphism
from the former compactification to the latter.  The correct package is their
normalized graph correspondence followed by the normalization/Stein
contraction to the finite incidence.  The standalone
[DVR marking audit](dvr-marking-audit.tex) remains the canonical source for
the abstract finite-cover valuative theorem and quotient collision models.

ACV and Deopurkar provide the compactification technology cited in the paper.
They do not externally review the repository's specialized LL-degree or
boundary-class calculations; no such external review is recorded.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```
