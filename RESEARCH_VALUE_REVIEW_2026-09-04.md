# Final research-value review — 4 September 2026

This is a dated assessment and preservation plan, not a second mathematical
status ledger. `MATH_STATUS.json` and its canonical sources remain authoritative.
The 5 September follow-up extracted a GVC corollary and corrected the HC4
master status; both changes are recorded in that registry. The [maintenance retrospective](DISCOVERY_RETROSPECTIVE_AUDIT_2026-09-04.md)
has a separate purpose: repository integrity and computational failure modes.

## Judgment

GVC is the right headline. The strongest additional results to preserve are
the prescribed finite étale fiber construction and the two-pair Special Image
counterexample. Both have mathematical value independent of whether HC4 or
JC2 is ever settled. The fixed-map Hasse sequel is also a coherent result,
with a separate analytic review requirement.

The follow-up audit found a concrete normalization gap in the claimed HC4
master reduction. Its full equivalence with JC2 is now partial. A new exact
prolongation closes one motion sign, while the other remains open; see the
[correction and surviving branch](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md).
The remaining full HC4 and JC2 problems are not close merely because many
special cases have been eliminated. Their surviving global compatibility
questions are substantive.

The main missed opportunity is extraction: state finished results without
making them wait for a larger conjecture, and turn successful proof mechanisms
into compact, checkable interfaces. One concrete extraction from binary GVC
is described below.

## Scope and evidence

The reviewed working snapshot has 1,202 registry entries. Of these, 898 have
canonical sources outside `elkies-k3/` and `elliptic-curves/`: 827 proved,
36 partial, 15 open, 16 parked, and 4 falsified. These are registry labels,
not 827 independently reverified theorems. This review mapped that ledger
and read selected load-bearing proof chains; it did not reprove every entry
or audit the separate active K3 and elliptic-rank campaigns.

Snapshot identifiers:

- Git HEAD: `8f4ea112757000d9d89d7fd0843c0e77ef20e702`.
- Working `MATH_STATUS.json` SHA-256:
  `f4bfb93cd063a6c29acb92801ca043f2b9688f72e62fe2b6afd1307d932b9d62`.
- The worktree already contained substantial unrelated changes. No research
  certificate, status row, or active proof was rewritten in the initial
  review phase. The dated follow-up below makes the stated GVC and HC4 changes.

The original determinant and collision were replayed independently. The
[AFP entry](https://isa-afp.org/entries/Jacobian_Counterexample.html) also
records an external Isabelle verification of that foundational map. Its
discovery is external to this repository; the value assessed here is the
subsequent mathematics.

## Results to extract and preserve

| Priority | Result | Why it warrants its own readable record | Remaining assurance or presentation work |
|---|---|---|---|
| 1 | Binary GVC and the ternary failure | Exact dimensional boundary; one positive proof and a concrete negative mechanism | Independent review and completion of the binary formal bridges |
| 2 | Every characteristic-zero finite étale algebra of rank at least three is a full Keller fiber in dimension three | A universal construction, not another isolated collision; actual algebra and geometric degree are controlled | External mathematical review; preserve the small paper-facing Lean certificate |
| 3 | SIC fails in exactly the pair dimensions at least two | A finished dimensional result with a short all-order coefficient identity | A focused two-pair manuscript or deliberately versioned sequel to the frozen three-pair paper |
| 4 | Quantitative Hasse failures for one fixed degree-five Keller map | Changes the quantifier from a map for each algebra to many failures on one map, with optimal geometric degree | Independent review of the Euler product, asymptotic constant, and minimality inputs |
| 5 | HC4 reductions and the reopened negative motion sign | The audit corrects an unsupported closure and supplies a valid positive-sign obstruction | Close the remaining negative sign; the full pencil equivalence is partial |

The priority order concerns preservation and review effort, not a prediction
of journal acceptance or a claim that a literature search proves novelty.

### GVC: protect the actual proof boundary

The [manuscript](papers/generalized-vanishing-two-variables/main.tex) has a
plausible and economical binary mechanism. In shifted-ray separation, the
least-ordinary-degree endpoint contributes valuation `s`. A non-Frobenius
term gets two coefficient divisibilities, while its factorial contribution
is at least `s-1`; every other integral endpoint has larger degree. This
gives a unique least-valuation term. The moving envelopes then convert that
local statement into global support separation. I found no contradiction
in this chain on close reading; this is not an independent referee report.

The [Lean audit](formal/gvc/README.md) now covers the concrete ternary
counterexample, its characteristic-zero and dimension extensions, and the
paper's arbitrary-profile failure theorem. The binary theorem still needs
the translated constant-term/no-matching bridge, field and coordinate
packaging, number-field shifted-ray transfer, and global envelope closure.
The finite Hall core and common-threshold support argument are already
formalized. The paper registry's older description understated the negative
side's formal coverage; that navigation paragraph was corrected in this review.

The right review assignment is therefore precise: audit Proposition 4.1 and
its use at zero output shift, then the Hall-to-envelope bridge. Another
bounded moment table would add little assurance to those steps.

The comparison with [Wilson's Gaussian theorem](https://arxiv.org/abs/2607.23887)
and the origin in [Long's Gaussian example](https://arxiv.org/abs/2607.18186)
are substantive and already acknowledged. A current presentation should also
position the ternary construction against
[Dvorsky's five-variable GVC counterexample](https://arxiv.org/abs/2608.07338).
These are related constructions, not interchangeable fixed-dimensional
statements. In particular, `Lambda = Delta^6` does not settle the ordinary
Laplacian problem.

### Finite étale fibers: the strongest independent sequel

The [paper](papers/common-arithmetic-fibers/README.md) prescribes
`Spec K[T]/(P)` as the whole fiber, with geometric degree `deg P` and
coordinate degree at most `6 deg P + 2`. The
[formal package](formal/finite-etale-keller/README.md) reaches the actual
characteristic-zero map, determinant, function-field degree, fiber
representation, and abstract finite étale corollary.

That is a clean theorem to explain to someone who has never seen the
repository. Keep it focused. The supplied positive-characteristic
presentation, noncanonical primitive-element choice, universal promoted map,
monodromy, and stable moduli have different hypotheses or assurance levels.
They should not obscure the already complete characteristic-zero statement.

The [fixed-map Hasse sequel](papers/fixed-map-hasse-failures/README.md) is
worth retaining separately. Its constructed targets have count asymptotic
to a positive constant times `B / sqrt(log B)` and its geometric degree is
five. The formal local-solubility and height results do not formalize the
analytic counting argument. The existing verification/literature notes
already give a sensible review boundary.

### Two-pair SIC: do not wait for minimum degree

The [canonical counterexample](extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
has bidegree `(4,4)` and an exact nonzero mixed contraction

\[
\mathcal E_2(QF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}
\qquad(m\ge1),
\]

while every pure contraction vanishes. The known one-pair theorem makes
the minimum failing pair dimension exactly two. The
[frozen three-pair paper](papers/three-pair-image-counterexample/README.md)
correctly preserves its older bound, but the stronger finished result needs
equally prominent publication navigation.

There is no conflict with binary GVC. The SIC witness has coefficient-matrix
rank five; it is not a separable symbol-times-polynomial input of the form
used by GVC. This distinction itself is a useful explanation of why the
two dimensional classifications differ.

The unresolved bidegree `(3,3)` problem asks for minimum balanced degree.
It does not weaken the two-pair counterexample. Its semistable Rodrigues
survivor is already known to be SIC-safe, so restoring moment–nullcone
equality is not a viable target. The most informative continuation would
be a relative-period/recurrence certificate on an actual moment-zero
component, including exceptional parameter fibers and endpoint terms.
The [holonomic probe](extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md)
explicitly warns that its generic benchmark fibers are not moment-zero
fibers. More fitted moments there would not close the classification.

## GVC consequence extracted on 5 September

The proposed support-separation corollary has become `GVC2SC`, recorded in
[MATH_STATUS.json](MATH_STATUS.json) and proved in the
[finite-certificate note](extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md).
The additional deductions give a unique Hall direction over the original
field, a finite rational decision procedure, and the explicit uniform mixed
cutoff `m > (deg P + deg Lambda) deg Q`. Section 6.1 of the active GVC
manuscript now includes the result.

The universal proof depends on the existing Hall/shifted-ray/envelope
argument; the exact checker verifies positive certificates and finite
regressions. Neither the full binary proof nor this corollary has received
external review or a complete formal verification. The
[original proposal](archive/RESEARCH_VALUE_REVIEW_GVC_PROPOSAL_2026-09-04.md)
is retained as history; its weaker field-extension and cutoff statements
are superseded by the canonical note.

## HC4: preserve the reduction, attack the correct boundary

The [relative-nilpotent master reduction](HC4_RELATIVE_NILPOTENT_MASTER_REDUCTION.md)
is now partial. The [5 September audit](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md)
computes the frozen-to-adapted transition explicitly: its normalized motion
determinant is `pq/a^2`, so `d(pq)=0` does not follow. The old unit ideal
remains correct for the augmented system assuming that extra equation.

Differentiating the actual branch relations closes `p=q=a` independently.
The sign `p=q=-a` has a compatible exact finite jet with nonzero `d(pq)`;
its all-order integrability and polynomial realization remain unknown.
This is a concrete reopened opportunity, not a counterexample. The
nonzero-Hessian-direction hypothesis has also been added explicitly to
the master statement. The lower-rank reductions and `PHC4 => JC2` survive;
`JC2 => PHC4` is not currently proved by this route.

For full HC4 the two useful continuations remain distinct:

- **Direct quintic classification:** nonreduced Hessian–Schur modules,
  nonlinear denominator components and positive-defect strata, together
  with the separate top-rank-at-most-two synchronization problem. The
  [direct frontier](HC4_PROJECTIVE_POLAR_GEOMETRY.md) and the scope of
  `OP-HC4-D5` must be read together; many generic and split-linear rows are
  already closed.
- **Descent from Meng–Yang:** polynomial termination of the uniquely
  determined formal graph, or a genuinely moving/nonlinear Schur pivot.
  The [graph theorem](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md) supplies formal
  normal solvability in every order. More successful formal jets alone
  therefore cannot distinguish a polynomial counterexample.

A concrete literature update matters here. Zixiang Ni's
[*The Quartic Hessian Conjecture in Dimension Four*](https://arxiv.org/abs/2608.14217)
was submitted on 14 August 2026 and proves the degree-at-most-four case,
overlapping `HC4CQ1`. The oldest local commit touching that proof note is
`6de7fd1`, dated 28 July; a local commit date alone does not establish public
priority. The checked HC4 notes do not cite Ni. Compare the proofs and
document their relationship before presenting the quartic theorem as a new
standalone contribution. The partial pencil reduction and direct quintic
work have different scopes.

## JC2: two degree notions and a real global obstruction

The externally reduced larger-coordinate-degree frontier `125` is not
geometric degree `125`. Likewise the geometric-degree-four classification
is not a search through quartic coordinate polynomials.

The strongest conceptual JC2 package is the finite free normalization,
its freely generated missing-boundary class group, and the contrast between
locally valid covers and a global affine-plane source. The
[quartic programme](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md) has explicit
spectator countermodels with the desired local packets whose deleted open
is `A1 x Gm`. Thus a local packet contradiction cannot prove what those
countermodels already realize.

The unfinished geometric-degree-four gate is a **finite-cover-derived
bound** on connector number, minimized marked height, and contact, followed
by the actual braid/meridian constraints and boundary deletion. The
[endpoint-semigroup result](plane-jc/QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md)
shows that the abstract conductor data alone do not give such finiteness.
This is a missing theorem, not a missing enumeration run.

For the `(75,125)` F2 branch, the
[degree-six Stein reduction](plane-jc/F2_GEOMETRIC_DEGREE_SIX_STEIN_REDUCTION.md)
is a more focused continuation than indiscriminate Laurent expansion.
It localizes all affine branching to one terminal cubic packet and closes
the normal even rows `6,8,10`. It leaves normal odd rows and nonnormal
terminal slices; smoothness of the ambient Stein surface does not make
that slice normal. A conductor/attachment theorem that controls these
remaining slices would be a meaningful advance. Degree six itself would
still not settle all F2 degrees or all JC2.

The [residue-degree budget](plane-jc/OREVKOV_RESIDUE_DEGREE_BUDGET.md),
cyclic logarithmic-cokernel positivity, and arbitrary-degree
[support-at-most-six classification](plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md)
are useful standalone results. Their value need not depend on a JC2 proof.
The budget is a consequence/refinement of Orevkov's inputs, so publication
should identify the added content carefully. Generic affine conductor
length must not be charged as a relative logarithmic defect: strict étale
base change is already a countertest to that inference.

## Other work worth preserving, with explicit limits

- **Stable atomic maps and boundary moduli.** The
  [atomic spectrum](verified/ETALE_MONOID_ATOMIC_SPECTRUM.md), explicit
  families of stably inequivalent maps, and
  [Tschirnhaus non-descent](verified/GENERIC_TSCHIRNHAUS_NON_DESCENT.md)
  describe a substantial post-counterexample geometry. A focused theorem
  about the displayed families is preferable to claiming a global moduli
  stack that has not been constructed.
- **Characteristic-two lifting obstruction.** The
  [Cartier/de Rham obstruction modulo four](verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md)
  is an all-degree obstruction for Mondello's plane map, invariant under
  polynomial left–right equivalence, with a contrasting stabilized Witt
  tower. This is a coherent short topic. The tower is not a polynomial
  characteristic-zero plane lift. Optimization at the next Witt digit is
  less important than preserving that structural contrast. The original
  plane construction is [Mondello's external theorem](https://arxiv.org/abs/2608.02634).
- **LND ideal images.** The
  [complete principal-ideal result](extended-geometry/LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md)
  for `D=u*d_x+v*d_y`, and the nonprincipal plinth-power class, deserve
  comparison with the existing LNED literature. They do not solve LNED for
  arbitrary derivations and ideals. Preserve the positive theorem rather
  than treating every unsuccessful counterexample search as the result.
- **Dixmier rank two.** The
  [rank-two programme](extended-geometry/FIXED_RANK_DIXMIER_REDUCTION.md)
  has a genuine classical compression, but polynomial quantization and
  descent remain essential. Formal-local quantization does not close them.
  The intrinsic obstruction cocycles may be the publishable result; this
  review found no short route to a Weyl-algebra counterexample in rank two.
- **Ordinary Laplacian and dimension compression.** The
  [GVC power frontier](extended-geometry/GVC3_POWER_TAIL_AND_MINIMUM_FRONTIER.md)
  proves minimum power six within its one-profile architecture. A smaller
  witness must change that architecture. The
  [dimension-twelve cubic reduction](verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md)
  is useful evidence and infrastructure, but a new dimension record needs
  comparison with current external constructions. Neither bounded graph
  failures nor a template minimum prove a global minimum.

## Practical preservation order

1. Give the current GVC paper an immutable, synchronized source/PDF/formal
   snapshot and a precise independent-review brief. Consider the finite
   support criterion and explicit cutoff above as a corollary to review.
2. Finalize the finite étale paper around its existing paper certificate.
3. Give the two-pair SIC theorem a short current manuscript, preserving the
   frozen three-pair deposit as history. Do not wait for bidegree `(3,3)`.
4. Review the fixed-map Hasse paper's analytic layer separately.
5. Extract the HC4 restricted-equivalence proof only after auditing its
   global bridges; retain direct quintic and polynomial-termination work
   as clearly scoped continuations.
6. Preserve one concise JC2 theorem package and one explicit global restart
   gate. Avoid making its useful local/module results depend editorially
   on solving the conjecture.

These are recommendations. No external message, submission, or publication
was made during this review.

## Checks performed in this review

- `python3 scripts/verify_counterexample_independent.py`: passed the exact
  determinant, three-point collision, and coordinate degrees.
- `scripts/verify_gvc3_homogeneous_counterexample.py`: replayed through
  `m=6`, redirecting its output to a temporary directory; the generated
  bytes matched the preserved artifact exactly.
- `.venv/bin/python scripts/verify_finite_etale_keller_fibers.py`: passed
  degrees three through five, both reconstruction directions, normalization,
  collision algebra decompositions, and the fixed-family regression.
- `.venv/bin/python scripts/verify_hc4_affine_plane_prolongation.py
  --audit-existing-only`: passed artifact integrity and its explicit local
  proof boundary. This was not a fresh symbolic prolongation replay.
- `python3 scripts/check_lean_placeholders.py`: passed for 126 Lean files;
  its two permitted explicit axioms belong to the documented GMC(2)
  allowlist. Absence of placeholders is not a kernel rebuild.
- `python3 scripts/check_paper_certificate_imports.py`: passed both paper
  certificate import boundaries.
- The weighted GVC example above passed an independent exact SymPy spot
  check. No new large search or full Lean rebuild was run.
- All 26 local links in this review resolve, and the edited navigation passes
  whitespace checks. The whole-repository Markdown check found one unrelated
  unresolved link: `archive/elliptic-curves/external-audit-2026-09-04/README.md`
  points to the absent `elliptic-curves/notes/EXTERNAL_AUDIT_2026-09-04.md`.
  That separate archive was not modified.

The review's main output is a separation of finished contributions,
unverified proof interfaces, and actual open gates. It offers no new HC4,
JC2, ordinary-Laplacian, or global-minimality claim.
