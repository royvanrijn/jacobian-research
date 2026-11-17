# Literature transfer and gap map

14 November 2026. Primary-source review for the quadratic-trace extension. Search results were used to locate the sources; the mathematical imports below were checked in the linked papers or author-hosted notes. The general generator-scheme and quadratic-form methods are established literature.

| Source | Exact input | Consequence here and limit |
|---|---|---|
| [Arpin--Bozlee--Herr--Smith I, arXiv:2108.07185](https://arxiv.org/abs/2108.07185) and [II, arXiv:2205.04620v2](https://arxiv.org/abs/2205.04620) | Schemes of monogenerators; configuration spaces, twists and scaling quotients, especially II Sections3--4. | Treat choices as a functorial space, not a canonical generator. Representability of generator spaces is not new here. |
| [Edens--Reichstein, arXiv:2308.10096v1](https://arxiv.org/html/2308.10096v1) | Lemma3.1 treats the integral root quadric; Lemma4.2 and Remarks4.3--4.4 use binary blocks and odd-degree descent to produce trace-isotropic generators. | The quadratic-trace method can be transferred to the explicit Keller compiler. Their positive-characteristic hypotheses yield field-versality in N-3 coordinates here. Their N-4 essential-dimension compression does not transport the ambient Keller map automatically. |
| [Brassil--Reichstein, The Hermite--Joubert problem over p-fields, Theorem1.3 and Corollary8.1](https://www.math.uni-bielefeld.de/LAG/man/560.pdf) | At prime two, a binary trace-form criterion characterizes simultaneous first and second trace elimination over all extension fields; primitive generators are supplied when N is invertible. | This gives necessity as well as sufficiency for field-versality of D_N in that range. The underlying criterion is external, not a new theorem claimed here. |
| [Duncan, arXiv:0908.3220, Theorem1](https://arxiv.org/html/0908.3220v1) | Essential dimension of S7 is four over every characteristic-zero field. | Combined with our field-versal A4 construction over bases containing sqrt(-14), gives ktdim(S7)=4 on those bases. This lower-bound theorem is external. |
| [Casselman, Quadratic forms over local fields, Proposition9.2](https://personal.math.ubc.ca/~cass/research/pdf/GammaQ.pdf) | Quadratic forms of dimension greater than four over nonarchimedean local fields are isotropic, including dyadic completions. | The trace criterion has no finite-place obstruction over number fields for N>=7. |
| [Milne, Elliptic Curves, second edition, Aside2.9](https://www.jmilne.org/math/Books/EC2.pdf) and [Bayer-Fluckiger--Parimala--Serre, Hasse principle for G-quadratic forms, introduction](https://arxiv.org/abs/1305.3161) | Classical Hasse--Minkowski for quadratic forms over number fields. | Real signatures give the exact number-field realizability criterion for D_N. The trace form is an ordinary quadratic form; no Hasse principle for more general equivariant forms is assumed. |
| [Stacks Project, quotient stacks](https://stacks.math.columbia.edu/tag/04UI) and [inertia](https://stacks.math.columbia.edu/tag/06PA) | Quotient objects retain stabilizers, and sections induce compatible homomorphisms on automorphism groups. | Needed for the next receiver gate; no new global stable Keller moduli stack is supplied by this note. |

## What changed in the research picture

The earlier real-split failure of double depression was correct, but did not mean that the N-3 map was useless. Its obstruction is exactly the trace quadratic form. This now gives one all-degree theorem covering existence, arithmetic failures, good-generator denominators, generic monodromy and stable atomicity.

A crucial quantifier distinction emerged: realizing every algebra over a number field k does not imply a field-versal map over k. Universality for ktdim must include all extensions L/k. The binary-block test is designed to preserve that quantifier. In rank seven the single constant delta^2=-14 makes the test work over every extension, allowing an exact dimension result.

## Rules and remaining gates

| Question | Proved rule | Remaining gate |
|---|---|---|
| Can scaling be removed without roots? | The retained positive weight vector must be primitive; a polynomial monomial section then exists. A coprime pair is sufficient but unnecessary. | Source/target normalization using target coordinates needs a separate affine polynomial landing proof. |
| Can one add another natural linear constraint? | Only the scalar line and trace hyperplane are proper S_N-invariant linear subspaces. | Use nonlinear constraints or a different construction. |
| Does double trace zero realize A? | Exactly when the trace quadratic form has a zero outside its radical. | For N invertible the binary test is exact, so its negative cases need another compiler; a smooth-locus converse in modular ranks remains separate. |
| Does a field-size or u-invariant condition prove versality? | Only if its conclusion holds over every extension field; the constant binary test does. | Do not promote fixed-field arithmetic to an essential-dimension upper bound. |
| Is the compiler canonical after optimal compression? | Its double-trace family still has N-6 stable-class parameters for each realizable A. | Another compiler, or a correspondence retaining the choices. |
| Does a quotient stack solve descent? | It retains the automorphisms actually present in its decorations. | A section over BS_N must support the full S_N automorphism group; compute receiver inertia before attempting gluing. |

## Retained next receiver idea

For every N>=5, use all primitive trace-zero generator lines with only d1*d3 nonzero, rather than imposing d4*d5 nonzero. If H_N is that open in the split trace projective space, the quotient [H_N/S_N] should identify with the full trace-zero incidence quotient by G_m. The cubic and linear coefficients force any geometric homothety stabilizer to satisfy lambda^2=1. A proof of this identification and its exact inertia would test whether the receiver can admit a section from BS_N.

This is an explicitly recorded next gate, not a proved result of the present note. Likewise the singularity of affine toric coarse quotients and the nonconstant units in C=1 source levels need a separate written scope before being used as obstructions.

## A recent AI-assisted lead, and its exact compatibility limit

[Chebolu--Minac--Nikzad--Ure, A septic covariant and the Hermite--Joubert problem in degree seven, arXiv:2609.15786v1](https://arxiv.org/html/2609.15786v1) proves the characteristic-zero field-extension case using a covariant divided by the polynomial derivative (Theorems1.2--1.3). Section8 describes AI participation in discovery and verification. Its wider all-characteristic finite-etale statement is announced for future work; we do not import that announcement as a proved input.

This is a useful alternative to the quadratic-trace constraint: it kills the first and third Newton traces rather than the first and second. Our compatibility check shows why it does not immediately yield a rank-seven universal A4 map over Q. In the compiler's ascending coefficient notation it kills a6 and a4, leaving seed weights six and eight. Their gcd is two, so the remaining seed-scaling normalization has a mu2 obstruction. Its denominator and good-coefficient conditions also require a separate audit on a covariant image, unlike the rational-quadric density argument above. A different polynomial landing or a method retaining that Kummer datum is the concrete transfer gate.

This distinguishes the new published generator breakthrough from a completed Keller compression; no large covariant expansion or unsupported universality claim was needed for the comparison.
