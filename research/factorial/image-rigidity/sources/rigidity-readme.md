# Polynomial composition rigidity via critical values

Proof of Furter's rigidity conjecture R(m,n) for all m, n.

**Theorem.** Let a, b be complex polynomials with a∘b(0) = 0, (a∘b)'(0) = 1 and a∘b ≠ z. Then

    ord₀(a∘b − z) ≤ deg a + deg b − 1.

Consequences: the composition-coefficient map A^(m+n) → A^(m+n) is finite flat of degree C(m+n, m);
the closure of every length-two multidegree stratum of Aut(A²) is the one conjectured by Furter
(the length-two Polydegree Conjecture); the Strong Factorial Conjecture holds for
X₁⋯X_m(μ₁X₁ + ⋯ + μ_mX_m).

- `proof.pdf` — the paper
- `rigidity.typ`, `refs.yml` — Typst source; build with `typst compile rigidity.typ`
