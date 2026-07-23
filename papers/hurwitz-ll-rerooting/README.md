# Hurwitz--LL compactification of rerooting

This paper contains the admissible-cover closure, rerooting quotient stack,
marked-zero-fiber LL degree, formal collision comparison, and the boundary
Picard and caustic/Maxwell calculations.

The H1--H3 dependency chain is now assessed as partial/high-risk.  The
standalone [DVR marking audit](dvr-marking-audit.tex) proves the abstract
finite-cover valuative theorem and the quotient/coarse collision models, and
isolates the missing global comparison with the repository-specific
admissible-cover contraction.  It is the canonical status source for H1--H3.

ACV and Deopurkar provide the compactification technology cited in the paper.
They do not externally review the repository's specialized LL-degree or
boundary-class calculations; no such external review is recorded.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```
