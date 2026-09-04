# Higher $2$-power Selmer programme for the R17 record pair

Date: 2026-09-04  
Status: **protocol fixed; complete $2$-Selmer groups and all higher-filtration dimensions remain `UNKNOWN`**

## Question

For ICARM curves 356 and 385, let $E_t$ be the record fibre, let $G_t$
be the specialized generic MW17 subgroup, and write

\[
H_t=\delta_2(G_t)\subseteq \operatorname{Sel}_2(E_t/\mathbf Q).
\]

The exact record certificates prove that $H_t$ has dimension 17 and that
the twelve displayed exceptional points span a subspace

\[
W_t\simeq \mathbf F_2^{12}
\quad\text{in}\quad
\operatorname{Sel}_2(E_t)/H_t.
\]

The full $2$-Selmer groups are not yet known.  In particular, no current
certificate says that $W_t$ is the whole residual Selmer group.

The higher-$2$-power object to compute is not a literal inclusion of groups
called \(S_2,S_4,S_8\).  Define instead the image filtration

\[
F_j(t)=
\frac{\operatorname{im}\!\left(
  \operatorname{Sel}_{2^j}(E_t/\mathbf Q)
  \longrightarrow \operatorname{Sel}_2(E_t/\mathbf Q)
\right)}{H_t},
\qquad j=1,2,3.
\]

Then

\[
F_1(t)\supseteq F_2(t)\supseteq F_3(t)\supseteq W_t.
\]

This is the precise meaning here of a residual $2/4/8$-Selmer filtration.

## Two necessary corrections

Every rational point class lifts through every $2^j$-Selmer level.  Hence
the twelve known exceptional directions automatically lie in all three
spaces.  Higher descent does not discover or explain that fact; it measures
which *additional* $2$-Selmer classes remain compatible with deeper
divisibility.

Also, the first drop has even codimension.  The Cassels--Tate pairing on the
$2$-Selmer group is alternating and its kernel is the image of the
$4$-Selmer group.  Consequently

\[
\dim F_1(t)-\dim F_2(t)\equiv0\pmod 2.
\]

Thus the illustrative profile $16\to13\to12$ cannot occur.  A compatible
example is

\[
16\to14\to12.
\]

For these two curves the known rational $2$-torsion is trivial.  If
\(\Sha_t=\Sha(E_t/\mathbf Q)\), the filtration separates the stable
Mordell--Weil contribution from the successive images

\[
\Sha_t[2],\qquad 2\Sha_t[4],\qquad 4\Sha_t[8]
\]

inside the residual $2$-Selmer quotient.  It is therefore a diagnostic for
extra Mordell--Weil directions versus increasingly divisible
Tate--Shafarevich classes, not a causal theorem for why the twelve displayed
points occur simultaneously.

## Decisive profiles

The following interpretations are unconditional once the corresponding
complete descents are certified.

| profile $(\dim F_1,\dim F_2,\dim F_3)$ | conclusion |
|---|---|
| `(12,12,12)` | The known residual $\mathbf F_2^{12}$ is the full residual $2$-Selmer group. Both the rank lower and upper bounds are 29. |
| `(16,12,12)` | Four extra $2$-Selmer dimensions are removed by the first Cassels--Tate pairing; the stable image is exactly the known rational block. The rank is 29. |
| `(16,14,12)` | Two dimensions die at each higher stage; the image of $8$-Selmer in $2$-Selmer is exactly the known rational block. The rank is 29. |
| `(*,*,>12)` | Higher descent still leaves room beyond the displayed rank-29 subgroup. No extra rational direction follows without a point or another exact argument. |

The first row needs no higher descent and must stop after the complete
$2$-Selmer computation.  More generally, once a certified stage equals
$W_t$, later stages are unnecessary for the rank bound.

## Computation order

1. Complete the unconditional global $2$-Selmer computation for 356 and
   385.  The active quotient-native class/unit route remains the current
   bottleneck.  No GRH class-group bound may be mixed into the first
   certificate.
2. Freeze a basis of $F_1(t)$, the complete all-place condition matrix, and
   the local contribution ledger before using the twelve exceptional points
   as labels.
3. Load the held-out exceptional points, verify that their quotient image is
   exactly $W_t\simeq\mathbf F_2^{12}$, and extend those rows to a basis of
   $F_1(t)$.
4. If $\dim F_1>12$, extend $W_t$ to a basis of $F_1$, construct two-covers
   only for the complementary basis, and compute the Cassels--Tate matrix on
   $F_1/W_t$.  Pairings with $W_t$ are already zero because these are
   rational-point classes.  Adjoining $W_t$ to the matrix kernel gives
   $F_2$.  Do not enumerate all $2^{\dim F_1}-1$ nonzero residual classes.
5. Quotient the kernel by the known stable floor $W_t$.  Apply the
   four-cover/two-cover pairing or explicit `EightDescent` only to that
   complement.  Preserve its full linear structure: individual basis-cover
   outcomes do not by themselves determine the image dimension because a
   liftable linear combination need not be a chosen basis element.  Running
   higher descent on the twelve known rational classes is tautological and
   wastes the dominant computation.
6. Record the dimensions, pairing ranks, local pairing contributions, exact
   cover models, software version, class-group assumptions, wall/RSS limits,
   and every incomplete terminal.  An empty or capped higher-descent call is
   `UNKNOWN` until its semantics are independently checked.
7. Compare 356 and 385 through basis-independent data: the dimension profile,
   pairing ranks, and sorted local suppression/contribution profiles.  Their
   raw Kummer coordinates live in different cubic fields and have no direct
   coordinatewise comparison.

The existing `build_fermigier_rank20_residual_selmer.py` is useful API
precedent, but its all-nonzero-cover pairing mode is not suitable at residual
dimension 12 or larger.  The record implementation must stay basis-level.

## Prospective search consequence

For a candidate fibre, the conceptual score is the dimension of a global
intersection of local Selmer conditions.  It cannot be estimated from the
known MW17 local image alone.  If $V_t$ is a certified global squareclass
envelope and $A_t$ is the stacked all-place condition map, the relevant
quantity is

\[
s_2^{\mathrm{res}}(t)=\dim\ker A_t-17,
\]

followed, for survivors, by the higher image dimensions.  Both the ambient
dimension and the condition rank may vary with $t$.  Therefore a raw
codimension or a matched list of local MW17 fingerprints is not yet a Selmer
score.

The current quotient-rank escape detector already has the correct
post-descent interface: it reports complete local-condition incidence and
leave-one-place-out suppression only after the global domain and all relevant
places are certified.  Prospective ranking should be reopened when that
measurement is available on both positive controls.  The frozen CRT
experiment's local MW17 matches and detector-insensitive point-search misses
do not supply this calibration.

## Literature boundary

Pasten--Salgado prove a non-thin set of strict rank jumps for suitable double
elliptic K3 surfaces.  The repository has already applied that theorem to the
published R17 fibration.  It proves abundance of fibres of rank above 17, not
a twelve-dimensional jump at a specified fibre.

Garbagnati--Salgado explain how special multisections produce new sections
after base change.  This validates the carrier mechanism but does not identify
twelve simultaneous independent directions on 356 or 385.

Mazur--Rubin and Klagsbrun--Mazur--Rubin show that quadratic twisting changes
Selmer ranks through local conditions and can produce prescribed or broadly
distributed $2$-Selmer ranks under stated hypotheses.  Those results supply
the right language and comparison model; they are not a specialization
theorem for this non-isotrivial R17 surface and do not turn Selmer classes into
rational points.

Watkins--Donnelly--Elkies--Fisher--Granville--Rogers used successive
$2$-, $4$-, and $8$-Selmer tests and analytic bounds in a large
quadratic-twist computation.  That is the relevant filtering precedent.  Its
fast higher-pairing layers exploit the rational $2$-torsion/isogenies of the
congruent-number family, whereas curves 356 and 385 have no rational
$2$-torsion.  Magma does implement general four- and eight-descent over
\(\mathbf Q\), but the precedent is an architecture, not a drop-in algorithm
for the record pair.

No theorem identified in this literature pass packages the specialization
$17\to29$ as one coherent twelve-dimensional event.  The resulting open
problem is to relate large specialization jumps on a fixed non-isotrivial
elliptic surface to the incidence and higher-divisibility filtration of its
specialization-dependent Selmer conditions.

## References

- H. Pasten and C. Salgado,
  [*Non-thin rank jumps for double elliptic K3 surfaces*](https://doi.org/10.1007/s00229-024-01554-2).
- A. Garbagnati and C. Salgado,
  [*Rank jumps and multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).
- B. Mazur and K. Rubin,
  [*Ranks of twists of elliptic curves and Hilbert's tenth problem*](https://doi.org/10.1007/s00222-010-0252-0).
- Z. Klagsbrun, B. Mazur, and K. Rubin,
  [*Disparity in Selmer ranks of quadratic twists of elliptic curves*](https://doi.org/10.4007/annals.2013.178.1.5).
- M. Watkins, S. Donnelly, N. D. Elkies, T. Fisher, A. Granville, and
  N. F. Rogers,
  [*Ranks of quadratic twists of elliptic curves*](https://doi.org/10.5802/pmb.9).
- T. Fisher, E. F. Schaefer, and M. Stoll,
  [*The yoga of the Cassels--Tate pairing*](https://arxiv.org/abs/0710.2079).
- [Magma V2.29 handbook: Mordell--Weil groups and descent methods](https://magma.maths.usyd.edu.au/magma/handbook/text/1570).

## Claim boundary

- `dim F_1`, `dim F_2`, and `dim F_3` are all `UNKNOWN` for both records.
- The only certified residual dimension is the lower bound
  `dim W_t = 12` inside the rational-point image.
- No incomplete BNF, Cassels--Tate, four-descent, or eight-descent run is a
  rank or Selmer upper bound.
- No bounded point-search miss is evidence that a Selmer class is not
  rational.
- The literature comparison motivates the programme but is not a theorem
  explaining the record fibres.
