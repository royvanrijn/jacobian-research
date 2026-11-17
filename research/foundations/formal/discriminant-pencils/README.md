# Polynomial tangent pencils in Lean

This is a **partial formalization** of [Generic Discriminants of Polynomial Tangent Pencils](../../../../unpublished/discriminant-pencils/README.md). The project is axiom-free, but it does not prove the manuscript's main theorem in full. Toolchain versions are pinned locally.

```sh
lake build
```

The modules verify the tangent parametrization and jets, contact-factor algebra, projective parametrization, formal graph at infinity, cusp counting and collision exclusions, and the stated affine-exhaustion interfaces.

The complete generic geometric theorem remains a separate written argument. The formal source files specify precisely which geometric hypotheses are assumed and which algebraic consequences are proved; they should be consulted before citing a machine-checked result.
