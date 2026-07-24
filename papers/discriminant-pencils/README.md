# Generic discriminants of polynomial tangent pencils

This directory contains the short replacement paper about the reduced
discriminant curve of

\[
H(W)-sW+t.
\]

The paper proves only the generic cusp--node theorem. The former real-chamber
theorem and repository-specific endpoint slice have deliberately been left
out so that the manuscript has one question and one main result.

The accompanying Lean development in
`../../formalization/discriminant-pencils/` is currently partial.  It checks
the polynomial identities, contact calculations, local analysis at infinity,
the complete affine normalization-fiber exhaustion under the three explicit
contact exclusions, and the distinct cusp count without axioms.  It does not
yet contain a complete machine-checked proof of the main theorem.  Neither
this paper nor the repository should claim that the formalization is
complete.

Build with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The internal adversarial proof audit remains in `AUDIT.md`.
