# Labelled node saturation and the corrected H2 quotient

## Outcome

Two parts of the former H1/H2 gap are formal once the base is corrected.

1. The normalization and inertia of an admissible-cover node are determined
   by one finite phase quotient together with the character of the
   **label-preserving** cover automorphism group.
2. If \(\Gamma_N^{\mathrm{lab}}\) is the labelled normalized graph of the
   polynomial-cover map, then the compactified rerooting morphism is
   canonically
   \[
     [\Gamma_N^{\mathrm{lab}}/S_{n-1}]
       \longrightarrow
     [\Gamma_N^{\mathrm{lab}}/S_n],
     \qquad n=N-2.                                  \tag{0.1}
   \]
   It is finite, representable, etale, and of degree \(n\).

Thus the quotient-stack and finite-incidence parts of H2 do not require the
false identification of the admissible graph with
\(\overline M_{0,N}\).  They work unchanged over the corrected graph.
Moreover, full root-label equivariance is not a separate chart-gluing
problem: it is built into the global normalized graph closure.

The [recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) now identifies
\(\Gamma_N^{\mathrm{lab}}\) explicitly with the recursive logarithmic
branch-scale construction for arbitrary tame collision types.

## 1. One target node

Consider the Harris--Mumford completed node ring

\[
 R=
 k[[\tau,s_1,\ldots,s_r]]/
 (s_1^{e_1}-\tau,\ldots,s_r^{e_r}-\tau).            \tag{1.1}
\]

Put

\[
 L=\operatorname{lcm}(e_1,\ldots,e_r),\qquad
 P=\prod_{j=1}^r\mu_{e_j}.                          \tag{1.2}
\]

On a normalized branch one may write

\[
 \tau=q^L,\qquad
 s_j=\omega_jq^{L/e_j},\qquad \omega_j\in\mu_{e_j}. \tag{1.3}
\]

Replacing \(q\) by \(\zeta q\), with \(\zeta\in\mu_L\), changes the phase
tuple by

\[
 \Delta_L(\zeta)
   =(\zeta^{L/e_1},\ldots,\zeta^{L/e_r})\in P.       \tag{1.4}
\]

The map \(\Delta_L:\mu_L\to P\) is injective.  Consequently:

> **Labelled node-normalization theorem.**  
> The normalized branches of (1.1) form the phase quotient
> \[
>   P/\Delta_L(\mu_L),
> \]
> and their number is
> \[
>   \frac{\prod_j e_j}{L}.                          \tag{1.5}
> \]

This is the standard normalization multiplicity in admissible-cover
deformation theory; see Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626), Section 4.2.  It is not an inertia
formula.  The normalized stack has the twisted-cover interpretation of
Abramovich--Corti--Vistoli,
[Twisted bundles and admissible
covers](https://arxiv.org/abs/math/0106211).

## 2. Actual inertia

Assume the source nodes above the target node are labelled, so cover
automorphisms do not permute them.  Every label-preserving automorphism acts
on the source smoothing coordinates through a character

\[
 \chi:
 \operatorname{Aut}_{\mathrm{lab}}(f)
   \longrightarrow P.                              \tag{2.1}
\]

An automorphism preserves one normalized lift precisely when its phase
change can be absorbed by the reparametrization (1.4).  Hence:

> **Labelled inertia theorem.**
> \[
>   I_{\widetilde f}
>   =\chi^{-1}\bigl(\Delta_L(\mu_L)\bigr).           \tag{2.2}
> \]

This formula cleanly separates three objects that were previously
conflated:

\[
\begin{array}{c|c}
\text{object}&\text{finite group or set}\\ \hline
\text{raw node phases}&P\\
\text{normalized branches}&P/\Delta_L(\mu_L)\\
\text{inertia of one labelled lift}&
  \chi^{-1}(\Delta_L(\mu_L)).
\end{array}                                         \tag{2.3}
\]

For \(r\) unanchored cyclic tails of common index \(e\), the character is
surjective, so there are \(e^{r-1}\) normalized branches and diagonal
inertia \(\mu_e\).  More generally, if some component phases are forced to
be \(1\), put

\[
 M=\operatorname{lcm}\{e_j:\text{the \(j\)-th phase is forced to be }1\}.
\]

The surviving diagonal subgroup has order \(L/M\), with \(M=1\) when no
phase is constrained.  Thus one anchored tail kills the diagonal for equal
indices, but need not do so for unequal indices.  For example, at indices
\((2,3)\), anchoring the index-two phase leaves order \(3\).

The sextic calculation is exactly the case \(e=2\), after matching all
vertex actions along the radial chain:

- every normalized radial lift has trivial inertia, although its raw nodal
  cover can have connector automorphisms;
- a pairwise Maxwell node has two normalization branches and inertia
  \(\mu_2\);
- a triple Maxwell node has four normalization branches and inertia
  \(\mu_2\).

For several target nodes, intersect the conditions (2.2) node by node.
Shared connector components merely make the characters dependent; they do
not introduce a new kind of gluing datum.

## 3. The global labelled graph

Let \(\mathcal R_N^{\mathrm{lab}}\) be the fully labelled root-stable
compactification and let \(\mathcal H_N^{\mathrm{lab}}\) be the corresponding
labelled normalized admissible-cover stack.  On their common polynomial
open there is an \(S_n\)-equivariant map, where \(S_n\) permutes the
\(n=N-2\) simple zero-fiber points.  Define

\[
 \Gamma_N^{\mathrm{lab}}
 =
 \operatorname{Norm}
 \overline{\Gamma\bigl(
   \mathcal R_N^{\mathrm{lab},\circ}
     \dashrightarrow
   \mathcal H_N^{\mathrm{lab}}
 \bigr)}.                                           \tag{3.1}
\]

Graph closure is functorial for the \(S_n\)-action, and normalization is
functorial for automorphisms.  Therefore \(S_n\) acts canonically on
\(\Gamma_N^{\mathrm{lab}}\).

This observation removes one previously stated gap:

> **Equivariant-gluing proposition.**  
> The labelled collision charts of the corrected graph have canonical
> \(S_n\)-equivariant gluing because they are restrictions of (3.1).
> No additional gerbe cocycle or independently chosen permutation descent
> datum is required.

The recursive resonance theorem gives the explicit nested-screen
presentation of (3.1).  Degree five identifies its first weighted blowup,
and the degree-six `(2,2,2)` calculation remains the smallest complete
closed-form chart.

## 4. H2 over the corrected graph

Let \(S_{n-1}\subset S_n\) fix the selected regular zero-fiber label.  Define

\[
 \Gamma_N^{\mathrm{mark}}
   =[\Gamma_N^{\mathrm{lab}}/S_{n-1}],
 \qquad
 \Gamma_N^{\mathrm{unmark}}
   =[\Gamma_N^{\mathrm{lab}}/S_n].                  \tag{4.1}
\]

For any algebraic stack with an \(S_n\)-action, subgroup inclusion gives a
cartesian-locally constant fiber \(S_n/S_{n-1}\).  Therefore:

> **Corrected stacky-mark theorem.**  
> The forgetful map
> \[
>   \Gamma_N^{\mathrm{mark}}
>      \longrightarrow
>   \Gamma_N^{\mathrm{unmark}}                      \tag{4.2}
> \]
> is finite, representable, etale, and has degree
> \([S_n:S_{n-1}]=n=N-2\).

The theorem is independent of whether \(\Gamma_N^{\mathrm{lab}}\) is a
blowup, a root stack, or a more general logarithmic modification.  It is the
correct replacement for the old formula using
\([\overline M_{0,N}/S_\bullet]\).

There is also a useful normalization reformulation.  Over the generic point
of \(\Gamma_N^{\mathrm{unmark}}\), selecting one regular root gives a
degree-\(n\) etale algebra \(K_{\mathrm{sel}}\).  The source of (4.2) is
normal and finite with this generic algebra, hence

\[
 \Gamma_N^{\mathrm{mark}}
 \cong
 \operatorname{Norm}_{\Gamma_N^{\mathrm{unmark}}}
   (K_{\mathrm{sel}}).                              \tag{4.3}
\]

Thus the selected-root extension and the normalized finite incidence are
the same object by uniqueness of normalization.  The completed
root-collision rings and multicluster conductor calculations already proved
in the repository apply to this finite factor; additional branch-scale
coordinates live on the base \(\Gamma_N\) and are not discarded.

## 5. Simplified dependency structure

The corrected `H1-COARSE`/H2 logic is:

\[
\boxed{
\begin{array}{c}
\text{construct or identify }\Gamma_N^{\mathrm{lab}}
\\ \Downarrow\\
\text{take the formal subgroup quotient (4.2)}
\\ \Downarrow\\
\text{apply uniqueness of finite normalization (4.3)}
\\ \Downarrow\\
\text{transport the intrinsic conductor and reconstruction marking}.
\end{array}}
\]

Accordingly:

- the abstract corrected graph exists canonically by (3.1);
- \(S_n\)-equivariant label gluing is automatic;
- the compactified degree-\(N-2\) rerooting map is formal over that graph;
- the finite normalization and conductor package needs no new
  branch-scale blowups;
- the stable-target half is the normalized wonderful pullback, and
  source-vertex rigidity makes the component maps unique once the weighted
  target-flag fibers are known;
- the general radial source theorem constructs those fibers, maps, and node
  partitions on every ordered first-scale stratum;
- polynomial monodromy subforests determine every nested simple-branch
  resonance source dual graph and node partition;
- finite normalization identifies the coarse polynomial source graph with
  the normal wonderful target graph and makes the corrected H2 quotient
  unconditional;
- monodromy centralizers compute every radial and generic simple-resonance
  label-preserving character;
- recursive normalized initial forms give the algebraic flag positions,
  residue compatibility, contraction formulas, and simultaneous inertia
  characters at higher-codimension nested resonances; the selected-factor
  log-étale comparison supplies the global coverage and descent for
  `H1-STACK`.

This closes the former separate “H2 gluing” problem.  It does not require
`H1-STACK`, which is a separate theorem with stronger proof obligations.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_labelled_node_saturation.py
.venv/bin/python scripts/verify_degree_six_stack_inertia.py
.venv/bin/python scripts/verify_degree_six_stacky_fan_descent.py
.venv/bin/python scripts/verify_multicluster_ll_comparison.py
```

The first command exhausts 1,554 node-index profiles with one to four source
nodes and indices at most six, checks every permutation of their labels,
tests the anchored and unanchored inertia families, and verifies the
subgroup quotient degree through eight labels.  The remaining commands
specialize the theorem to the sextic chart and replay the intrinsic
multicluster conductor comparison.
