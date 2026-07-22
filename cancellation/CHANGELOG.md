# Cancellation programme changelog

Closed milestones moved out of the active research roadmap.

## Parameter equivalence

For fixed `(m,r)`, the common finite normalization together with its marked
affine reconstruction open recovers the selected parameter root.  The visible
target quotient is the weighted scaling torus

\[
 (Q,R)\mapsto(uQ,u^{-m}R),
\]

and it fixes every parameter root.  Any residual target congruence-kernel
element that lifts and preserves the affine reconstruction open must preserve
the two affine factors `P=AB`; the source identities force

\[
 A\mapsto A,\qquad y\mapsto uy,\qquad B\mapsto u^{m+1}B.
\]

It therefore fixes the reconstruction-pole residue

\[
 q=\left.\frac{B}{y^{m+1}}\right|_{A=0}.
\]

Hence the `mr` distinct geometric roots of `M_(m,r)` give `mr` distinct
stable left--right classes.  The proof, including stabilization and the
`m=1` factor-swap exclusion, is the
[cancellation-parameter faithfulness theorem](../papers/marked-root-multiplicity/cancellation-parameter-faithfulness.tex).
Computing the target kernel itself remains an optional automorphism-group
problem, not an input to parameter faithfulness.
