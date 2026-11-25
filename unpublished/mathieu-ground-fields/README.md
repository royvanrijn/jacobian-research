# Mathieu's Conjectures over Fields of Characteristic Zero

*Split Homogeneous Spaces, Compact Groups, and Ground-Field Obstructions*

Companion manuscript to the [released complex classification](https://hexagonmath.org/2611.00556v1).

The main theorem classifies the Reynolds kernel for split identity component over every characteristic-zero field and proves the obstruction from an active split rank-one subgroup in nonsplit groups. It also treats compact Hausdorff homogeneous spaces, connected reductive groups over real closed fields, the universal rational-algebra formulation, and the ground-field toral nullcone.

The connected real closed homogeneous-space classification reduces to the two compact exceptional products. The manuscript records an explicit obstruction to the unmodified discriminant-density seed. For anisotropic nonsplit groups it records the limits of the split-root argument and a low-degree obstruction on anisotropic norm-one quadrics.

[PDF](mathieu-ground-fields.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) · [Verification](verification.json)

## Publication boundary

The complex counterexamples and their exceptional moment calculations are in the deposited paper. The present proof starts from those rank-one moment identities and adds principal-orbit transfer and arithmetic specialization. The [quadratic complement manuscript](../mathieu-quadratic-complements/) concerns a different trace and its sharp dimension classification. The [finite mixed-moment manuscript](../finite-tests-for-mixed-moments/) develops an effective general test rather than a homogeneous-space classification.

Build with `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`. Mathematical review of the new principal-orbit and arithmetic specialization arguments remains open. The proof does not claim a full classification of nonsplit homogeneous spaces.
