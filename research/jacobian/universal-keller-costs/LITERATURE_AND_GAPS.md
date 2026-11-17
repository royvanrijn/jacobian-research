# Literature transfer and exact boundaries of the cost results

14 November 2026. The decisive recent input is the planar torus classification. The other ingredients are classical torsor and covering theory, together with the retained compiler-independent fixed-space argument.

| Primary source | Use here |
|---|---|
| [T. Shaska, Graded Keller maps and the Jacobian Conjecture, arXiv:2607.20210v2, Theorem3.4](https://arxiv.org/html/2607.20210v2) | Every characteristic-zero planar Keller map with nontrivial G_m equivariance is an automorphism. The proof here only needs linear actions and reproduces their elementary all-signature argument. The theorem is credited to the source, not claimed as new. |
| [Santibanez-Leal, foundational study, latest linked record v0.13](https://zenodo.org/records/22834207) | Its primary record description independently discusses opposite-sign and zero-weight planar rigidity and reconciles the all-signature result with Shaska. We used Shaska's accessible full theorem as the proof dependency; the Zenodo PDF could not be retrieved through the web tool. |
| [Stacks Project, Kummer theory, Section59.28](https://stacks.math.columbia.edu/tag/03PK) | The exact units--torsors--Picard sequence. On affine space it makes tame cyclic torsors constant; on a stack, torsion character lines can retain the missing data. |
| [J.-P. Serre, Revetements de courbes algebriques, Bourbaki749,1992, Sections2--3](https://www.numdam.org/article/SB_1991-1992__34__167_0.pdf) | Classical quasi-p obstruction for covers of the affine line. We use and prove only the elementary necessary direction, not a sufficiency claim for polynomial Keller maps. |
| [Stacks Project, Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) | A tame connected cover of P1 branched only at infinity has degree one. The prime-to-p condition is essential. |
| [Jensen--Ledet--Yui, Generic Polynomials, AppendixA.1](https://library.slmath.org/books/Book45/files/book45.pdf) | Additive Hilbert90 for the one-variable translation-torsor family. The earlier UAFROB1 proof already supplies the affine-group field-versality input. |

## What is new in the programme

UACANCOST1 combines the actual natural N-point orbit with the characteristic-free fixed-space test. It proves the N+1 ambient bound without a catalogue of small S_N representations, and handles the modular scalar-line quotient and the N=2, characteristic-two exception explicitly.

UATORUSCOST1 applies the published planar theorem to the fixed plane at that bound. If the torus is trivial there, irreducibility and weighted homogeneity force the moving component linear, and UANORMAL1 closes the other case. Thus a commuting nontrivial torus raises the canonical linear bound to N+2 without assuming JC(2).

UAFINITE1 moves the constant-discriminant obstruction from univariate maps to every globally finite Keller map in odd characteristic. The general Kummer criterion also identifies exactly why a character line on BH_q avoids that affine-target obstruction.

UAFROBCOST1 simplifies canonical marked descent to a line bundle. The same symmetry which provides this simplification forces geometric translation monodromy in one variable; composite degree therefore loses absolute atomicity. The plane construction remains optimal when that stronger property or one fixed framed universal map is required.

## What these inputs do not justify

- The canonical linear bound is not a lower bound for an unmarked universal map. Its hypotheses include a natural marked fibre and linear ambient representations.
- The torus bound is characteristic zero and requires compatible commuting linear actions. A trivial action gives no information about JC(2).
- The fixed-map affine-group minimum includes a specified G torsor family on a nonempty open target chart. It is not silently a theorem about arbitrarily selected unstructured subgroup fibres.
- The quasi-p criterion is necessary, not sufficient for a source A1 Keller map or for field-versality. A general covering curve is not automatically affine space.
- A canonical family on bundles can vary its coefficients with the input algebra. It need not be one fixed universal polynomial.
- Arithmetic primitive monodromy can prove atomicity over the original field while absolute atomicity fails after extending constants. Coefficient promotion changes this situation; identity stabilization does not.
- The modular omitted-point plane control is an exact written countercontrol, not a new priority claim or a construction at the canonical ambient lower bound.

The exact controls do not replace mathematical review of the all-degree proofs.
