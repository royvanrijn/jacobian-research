# AGENTS.md

This directory is a separate research programme on elliptic curves over
\(\mathbb Q\).  It inherits the repository-level instructions.  In particular,
`../MATH_STATUS.json` remains the sole mathematical-status authority.

## Targets and terminology

- The first search target is a curve with at least 21 rigorously independent
  rational points and natural-log conductor `log(N) < 182.72`.  If exact rank
  21 is claimed, a matching upper bound is additionally required.
- The second search target is a curve with at least 30 rigorously independent
  rational points.  This is the operational meaning of beating the public
  rank-at-least-29 record.
- A list of points is only a rank lower bound after independence is certified.
  A search score, analytic-rank estimate, or Selmer upper bound is not a rank
  lower bound.
- “Exact rank” requires lower and upper bounds that agree.  Conditional upper
  bounds must name their hypotheses.
- Conductor always means the exact conductor of a global minimal model.  A raw
  discriminant radical is only a heuristic proxy.

## Computation and evidence

- Keep family equations and their provenance in `families/`.
- Keep deterministic code in `ecsearch/` and command-line entry points in
  `scripts/`; tests belong in `tests/`.
- Compact pinned manifests belong in
  `../artifacts/generated-results/elliptic-curves/`.  Raw searches and
  checkpoints belong in the ignored `../artifacts/local/elliptic-curves/`.
- Record equations, parameter normalization, bounds, prime sets, random seeds,
  software versions, exact commands, and whole-file hashes.
- A bounded Hensel/CRT/lattice scan is an experiment.  Do not turn survival,
  smoothness, a Nagao score, or a PARI rank guess into a theorem.
- At a shaping prime forced into the discriminant, score the bad-reduction
  data.  Ordinary good-reduction `a_p` scores belong to different primes.
- Re-minimize every specialization before using discriminant valuations to
  infer anything about conductor.  Treat 2 and 3 separately from clean
  semistable primes.

## Structural-search discipline

- Once a certified subgroup is available, score newly found points by exact
  finite-quotient escape before spending more search height. Escape is strong
  positive evidence; non-escape at a finite set of quotients is not a
  dependence proof.
- A residual Selmer class can represent a new Mordell--Weil direction or a
  Tate--Shafarevich class. Do not promote Selmer dimension to rank.
- Néron--Severi vector enumeration produces divisor-class candidates only.
  Effectiveness, nefness, irreducibility and field of definition remain
  separate gates.
- Search local conditions on both charts of `P^1(Q_p)`. A homogeneous
  congruence near infinity still needs chart-unit checks and global
  minimalization after specialization.
- Every external CAS task must pin exact inputs, required outputs, resource
  limits and the interpretation of timeout/failure.

## Elkies--K3 bridge

- Start from `../elkies-k3/ELKIES_K3_PROCESS_ATLAS.md` and the nested
  `../elkies-k3/AGENTS.md`; do not infer the construction from a single local
  artifact.
- The K3 equation route is certified through `D12/MW5`, while its marked
  lattice route reaches rootless `MW17`. The orbit42 A11 output is provisional
  until `MATH_STATUS.json` admits it.
- A K3 MW rank does not by itself certify the rank of a rational
  specialization. Specialize, minimize, verify every point, and prove
  independence again.

## Candidate promotion gate

A target candidate is not promoted until a replay checks nonsingularity, a
minimal integral model, exact conductor, every displayed point, point
independence, and the claimed rank status.  A record candidate also needs an
independent implementation before any public record claim.
