# Sources and provenance

Primary sources are preferred.  This file also states exactly which claim is
being imported, so a later computation is not mistaken for a published
theorem.

## Construction and search methods

- Jean-François Mestre, [*Courbes elliptiques de rang >= 12 sur
  Q(T)*](https://gallica.bnf.fr/ark:/12148/bpt6k57325582/f175), C. R. Acad.
  Sci. Paris Sér. I Math. 313 (1991), 171--174.  Primary source for the
  generic-rank-12 construction underlying the six-root search space.
  Mestre's author page also hosts his ICM Zürich proceedings paper
  [*Constructions polynomiales et théorie de
  Galois*](https://www.imj-prg.fr/wp-content/uploads/2020/prix/mestre1994.pdf),
  useful background for the polynomial-construction viewpoint but not used
  here as a rank theorem.
- Stéphane Fermigier, [*Une courbe elliptique définie sur Q de rang >=
  22*](https://doi.org/10.4064/aa-82-4-359-363), Acta Arithmetica
  82.4 (1997), 359--363.  Source for the six-root Mestre construction, generic
  rank at least 12, the parameter `19754/39`, the published `E22` model, its 22
  independent points, and the historical score table.  The normalized
  parameter and conductor replay in this repository are independent
  computations.
- Koh-ichi Nagao, [*Construction of high-rank elliptic
  curves*](https://doi.org/10.24546/E0003610), Kobe Journal of Mathematics 11
  (1994), 211--219; [repository
  PDF](https://da.lib.kobe-u.ac.jp/da/kernel/E0003610/E0003610.pdf).  Primary
  source for the two six-root tuples, the quadratic base change with a
  thirteenth section, the printed rank-13 specialization, and the printed
  rank-21 curve and independent points.  This repository checks formulas,
  models, points, conductors, and local data but cites Nagao for independence.
- Jasper Scholten, [*Elliptic curves of high rank over function
  fields*](https://arxiv.org/abs/math/9709235), 1997.  Independent geometric
  analysis of the Mestre--Nagao examples.  It contextualizes the geometric
  rank and the role of base change; it is not evidence for any new
  specialization in this repository.
- Matthias Schuett and Tetsuji Shioda, [*Elliptic
  Surfaces*](https://arxiv.org/abs/0907.0298), 2009.  Reference for elliptic
  surface fiber geometry, Neron--Severi groups, Mordell--Weil lattices and the
  Shioda--Tate formula used in the exact section-7 generic-rank argument.
- Ronald van Luijk, [*K3 surfaces with Picard number one and infinitely many
  rational points*](https://arxiv.org/abs/math/0506416), Algebra & Number
  Theory 1 (2007), 1--15.  Reference for bounding characteristic-zero Picard
  rank using good-reduction Frobenius eigenvalues and exact point counts.  The
  present certificate needs no Tate-conjecture equality: specialization and
  the root-of-unity eigenvalue condition supply only an upper bound.
- Noam D. Elkies and Zev Klagsbrun, [*New Rank Records for Elliptic Curves
  Having Rational Torsion*](https://arxiv.org/abs/2003.00077), ANTS XIV (2020).
  Source for the rank-nine `Z/2Z` K3 family, its sections, modular
  specialization sieve, staged cutoffs, skew search regions, and discussion of
  Bayesian scoring.  It does not describe this repository's combined
  prime-power/CRT/Gauss objective.
- Noam D. Elkies and Mark Watkins, [*Elliptic curves of large rank and small
  conductor*](https://arxiv.org/abs/math/0403374), ANTS VI (2004).  Background
  for joint rank/conductor searches and for distinguishing conductor from
  discriminant height.
- Sang Yook An, Seog Young Kim, David C. Marshall, Susan H. Marshall, William
  G. McCallum and Alexander R. Perlis, [*Jacobians of Genus One
  Curves*](https://doi.org/10.1006/jnth.2000.2632), Journal of Number Theory 90
  (2001), 304--315.  Primary journal source for obtaining the Jacobian and
  covering map of a double cover of `P^1` from the classical covariant
  syzygy.
- Tom Fisher, [*The invariants of a genus one
  curve*](https://doi.org/10.1112/plms/pdn021), Proceedings of the London
  Mathematical Society 97 (2008), 753--782.  Authoritative invariant-theory
  framework for binary quartics as degree-two genus-one models and their
  Jacobians.  The exact normalization used here is additionally checked
  directly: source membership, target membership, and the covariant syzygy are
  evaluated in rational arithmetic rather than trusted from an unchecked
  formula transcription.

- John Cremona, Tom Fisher and Michael Stoll,
  [*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic
  curves*](https://arxiv.org/abs/0908.1741), Algebra & Number Theory 4 (2010),
  763--820. Algorithmic basis for converting residual Selmer classes into
  smaller explicit covering models before rational-point search.
- Nils Bruin and Michael Stoll, [*The Mordell--Weil sieve: proving
  non-existence of rational points on curves*](https://arxiv.org/abs/0906.1934),
  LMS Journal of Computation and Mathematics 13 (2010), 272--306. Reference
  for replacing unbounded height searches on higher-genus transport covers by
  exact local sieving against a known Jacobian subgroup.
- Noam D. Elkies, [*Three lectures on elliptic surfaces and curves of high
  rank*](https://arxiv.org/abs/0709.2908), 2007. Reference for the
  Néron--Severi/Mordell--Weil lattice viewpoint, alternate elliptic fibrations
  and K3 constructions used in the structural-search queue.

## Rank records and calibration data

- The public ICARM database entries for
  [curve 281](https://elliptic-rank.icarm.cloud/curve/281),
  [curve 282](https://elliptic-rank.icarm.cloud/curve/282),
  [curve 285](https://elliptic-rank.icarm.cloud/curve/285), and
  [curve 286](https://elliptic-rank.icarm.cloud/curve/286) provide the four
  `7fff-zip` equations, points, conductors and timestamps analyzed in the
  pinned 2026-08-22 snapshot.  The repository independently replays point
  membership and independence; conductor minimality is still source-reported.

- Andrej Dujella's current [history of elliptic-curve rank
  records](https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html) lists the
  Nagao--Kouya rank-at-least-21 curve in 1994, the Elkies--Klagsbrun
  rank-at-least-29 curve in 2024, and the Alpöge--Howell rank-at-least-30
  curve in 2026. His [rank-at-least-21
  page](https://web.math.pmf.unizg.hr/~duje/tors/rk21.html) provides a public
  transcription of the former model and its independent points.  The primary
  independence source used here remains Nagao's paper above.
- Dujella's [rank-at-least-30 record
  page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html) gives the 2026
  Alpöge--Howell curve and lists 30 rational points.  This is the same model
  and point set previously published on the ICARM leaderboard by the submitter
  `ranksunbounded`.  The repository independently proves the points
  independent; no exact-rank statement is imported.
- Andrej Dujella's [rank-at-least-29 record
  page](https://web.math.pmf.unizg.hr/~duje/tors/rk29.html) gives the 2024
  Elkies--Klagsbrun model and 29 independent points.  This supports an
  unconditional rank lower bound of 29.
- Noam Elkies and Zev Klagsbrun's public account, reproduced in
  [MathOverflow's background thread](https://mathoverflow.net/questions/477849/background-for-the-elkies-klagsbrun-curve-of-rank-29/478050),
  describes the rank-17 K3 sieve, 12 extra points, conditional exact-rank
  statement, height-matrix computation, discriminant factorization, and
  conductor support.  Status: public computational account, not yet a complete
  paper for the rank-29 construction.  The rank lower bound comes from the
  independent points; exact rank 29 is stated conditional on GRH.
- Zev Klagsbrun, Travis Sherman and James Weigandt, [*The Elkies Curve Has Rank
  28 Subject Only to GRH*](https://arxiv.org/abs/1606.07178), Mathematics of
  Computation 88 (2019), 837--846.  Source for the analytic/class-group method
  later applied conditionally to the rank-29 curve; it is about the earlier
  rank-28 curve, not a proof about the new curve.
- The generated [rank-21/conductor public-data
  audit](../artifacts/generated-results/elliptic_rank21_conductor_public_data_audit.json)
  and [rank-30 public-source
  audit](../artifacts/generated-results/elliptic_rank30_public_source_audit.json)
  pin and hash the exact source snapshots inspected through 2026-08-14.  They
  replay all 225 printed points on Dujella's rank-21 through rank-29 pages and
  find no public curve meeting either target in their declared scope. This
  historical audit predates the 2026 rank-30 publication and is not evidence
  against it.

## Software semantics

- [PARI/GP elliptic-curve function
  reference](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html),
  especially `ellminimalmodel`, `ellglobalred`, `elllocalred`, `ellrank`,
  `ellheightmatrix`, and `ellsaturation`.  The documentation states that
  `ellrank(E)` returns unconditional bounds `r1 <= rank(E(Q)) <= r2` and a list
  of independent non-torsion points.  A nonzero effort invokes a randomized
  point search.
- The Magma handbook chapters on [elliptic-curve rank and
  descent](https://magma.maths.usyd.edu.au/magma/handbook/text/1570),
  [elliptic-curve models and
  maps](https://magma.maths.usyd.edu.au/magma/handbook/text/1544), and
  [rigorous class-group
  computation](https://magma.maths.usyd.edu.au/magma/handbook/text/415)
  document the routines used in the recorded `u=42` diagnostic.  The exact
  17-point verification is a software computation.  The anonymous public
  calculator's observed 311.34 MB memory limit on the descent probes produced
  no rank upper bound; a resource failure is not mathematical evidence.
- John Cremona's [eclib source repository](https://github.com/JohnCremona/eclib)
  is the software provenance for the bounded `mwrank` diagnostic.  The pinned
  replay used source tag `20260707`, commit
  `538457c5fba5c89b5040e24c0c1e51bab0ad5997`, without installing it systemwide.
  Its machine-integer range failure, and the infeasible range exposed by a
  temporary arbitrary-integer patch, supplied no rank bound.

## Repository-local provenance

The formulas implemented here were transcribed into exact rational arithmetic
and then checked through internal identities and benchmark specialization.
Those checks are **verified computations**, not replacements for the cited
generic-rank proofs.  Conversely, a published generic-rank theorem does not
certify that its sections remain independent at every specialization.
