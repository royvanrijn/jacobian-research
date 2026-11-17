# A Sharp Coefficient-Ring Obstruction for Gaussian Marker Generation

**Unpublished · partial draft.**

For $P_b(Y)=Y^2+bY$ over a commutative Q-algebra, the manuscript computes the exact quotient of all marked periods by the module generated from the unit period using polynomial time coefficients and time derivatives. It identifies the invertible/nilpotent criterion and proves that generation at every field specialization or finite coefficient jet need not imply generation over the coefficient ring. This family is not pure, so it is not a counterexample to a Gaussian moment theorem.

[PDF](gaussian_marker_generation_ring_defect.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Verification](verification.json)

This note supplies the coefficient obstruction for the [base-marker reduction](../gaussian-constant-fibre-base-marker-certificates/). Its parameter-ring issue is separate from the field-valued classification in the [dimension-degree boundary paper](../gaussian-boundary/) and the [one-base Jordan-cancellation paper](../gaussian-quadratic-fibre-jordan-cancellation/).

Build from this directory with:

```sh
latexmk -pdf -jobname=gaussian_marker_generation_ring_defect -interaction=nonstopmode -halt-on-error main.tex
```

This writes `gaussian_marker_generation_ring_defect.pdf`. Only `main.tex` is required to compile; use a TeX installation with `latexmk`, `pdflatex`, and the packages listed in its preamble. The verification record identifies the current source and PDF and describes the checks performed. Earlier audit receipts in a source bundle, where present, retain their original scope and were not replayed during this revision.

Independent mathematical review and formal verification remain pending. Next: Review the source-kernel calculation over arbitrary Q-algebras and the idempotent criterion. Compare the specialization defect with the coefficient conditions in base-marker generation.
