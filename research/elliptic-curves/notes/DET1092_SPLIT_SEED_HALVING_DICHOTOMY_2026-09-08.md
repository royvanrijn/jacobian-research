# Productive versus inherited splits: an exact halving-or-cycle criterion

## Result and classification

**New deduction/application of standard Mordell--Weil descent.** Given a
certified mod2-injective generic subgroup `M` on a rational elliptic fibre
without rational2-torsion, there is a terminating, equation-only procedure
which decides whether any constructed rational point is outside `M tensor Q`.
It requires only one fixed generic finite-quotient footprint, rational group
operations, and rational roots of degree-four duplication polynomials.
Neither a full Mordell--Weil basis nor a rank upper bound is needed.

For a split degree-two multisection with inherited trace `Z`, this gives
an exact criterion for whether its branches add a direction. In particular,
the inherited2-Kummer class of the anti-trace `2Q-Z` need not be a dead end:
one extra halving layer reveals precisely the same residual direction as `Q`.

**Verified bounded application.** All38 prescribed decisions on five old
fibres finish within the frozen eight-halving-step limit. A separate
implementation checks every rank/dependence conclusion and exact group
identity; another checker reconstructs the generic17 prefixes and RR incidence.

| Constructed points / reference subgroup | Outcome | Successful halving layers |
|---|---|---:|
| Both branches of the calibrated302 genus-one carrier, over M17 | Independent | 0 |
| Their anti-traces, over M17 | Independent | 1 |
| Both branches at each of `1926/2699`, `2953/1671`, `5193/35630`, over M17 | Independent | 0 |
| Their anti-traces, over M17 | Independent | 1 |
| Both branches and anti-traces at `-528/3635`, over M17 | Inherited; exact cycles | 2–3 |
| Nine old blinded RR recoveries, over their respective M16 cores | Independent over the core | 0 |
| The same nine points, over full M17 | Inherited; exact cycles | 2–3 |

The two branches and their anti-traces supply only **one** quotient direction
per productive fibre. The nine M16-to17 results are generic recoveries, not
specialization gains above17. No new rank or new fibre is claimed. The
`1926/2699` comparison uses only its seed, not its later rank21 discoveries.

**Calibration boundary.** The302 carrier member was previously fitted to the
first historical unlock. Its equations, not exceptional point coordinates,
are inputs here, but this does not make its selection prospective. The
generic conic and the nine generic RR reconstruction controls use generic
data. The theorem decides whether a *specified constructed split* is
productive; choosing a useful302 member from generic data remains open.

## 1. What was already known

**Verified prior applications.** The
[conic progression](DET1092_CONIC_SEED_PROGRESSION_2026-09-08.md) already supplies
an infinite equation-only seed family. The
[trace/twist identity](DET1092_TRACE_TWIST_KUMMER_OBSTRUCTION_2026-09-08.md) proves
that the anti-trace transfers to an inherited2-Kummer class. The
[dependent split](DET1092_FUNNEL_FIRST_SEEDS_2026-09-08.md) already has exact
generic words, found after an inconclusive one-layer halving attempt. The
[first-seed covering note](DET1092_SEED_KUMMER_COVER_2026-09-08.md) already uses
mod2-injectivity to identify the generic rational-span image modulo2.

What is added is a complete halving-or-cycle decision theorem, its bounded
implementation, and source-only reconstruction of the old control points
followed by an exact uniform classification. The classifier is never given
the controls' known dependence words. This is an application of standard
descent/saturation ideas; literature novelty of the general algorithm is
not claimed. The standard exact rational-division operation is documented in
[Sage's point API](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_point.html#sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_field.division_points).

## 2. Hypotheses and generic-only footprint

Let `E/Q` be smooth with `E(Q)[2]=0`, and let `M=<P1,...,Pr>`.
Suppose a fixed homomorphism

\[
f:E(\mathbf Q)\longrightarrow V,\qquad 2V=0,
\]

sends the `Pi` to independent vectors. In the application, `f` is the product
of complete finite quotients `E(Fp)/2E(Fp)` at good odd primes. Odd-order
good reduction at one prime excludes rational2-torsion. The footprint is
constructed using **only** the generic sections, stopping when its column
rank is17 and that torsion exclusion has been found, before constructing
any candidate. The fixed prime cap is1009; cap failure is unresolved, not
a theorem about the fibre. All five frames pass with final primes at most179.

**Elementary consequence.** The `Pi` are independent and

\[
M\cap2E(\mathbf Q)=2M.
\tag{1}
\]

Indeed applying `f` to a relation or a doubled generic word forces all its
coefficients even. No rational2-torsion permits repeated division of a
relation, proving independence; the same argument proves (1).
If `M_sat` is the preimage of the rational span of `M` in `E(Q) tensor Q`,
then `M_sat/M` is finite of odd order. Thus rational-span membership is not
the same as integral membership in the displayed basis.

## 3. The exact algorithm

Start with the constructed point `P_0=P`. At step `n`:

1. Solve `f(P_n)=sum e_i f(P_i)` over `F2`.
   If insoluble, return **independent direction**.
2. There is at most one solution. Put `T_n=sum e_i P_i`, `e_i in {0,1}`.
   Solve the exact rational equation `2P_{n+1}=P_n-T_n`.
   If it has no rational solution, return **independent direction**.
3. If zero is reached or an exact point repeats, return the resulting
   **rational dependence certificate**. Otherwise continue.

Every successful half is unique because `E(Q)[2]=0`. The implementation caps
successful halving layers at8 and reports `UNKNOWN_STEP_CAP` if necessary.
There is no finite universal step bound asserted by the theorem.

**Proof of correctness.** Any generic word which could be subtracted to
make `P_n` divisible by2 has the unique parity in step1. Subtracting a
different word of that same parity changes the target by twice a generic
point, and therefore does not change rational divisibility. Thus either
escape proves that the residual class of `P_n` in `E(Q)/M` is not divisible
by2. It cannot be torsion of odd order. By (1) that quotient has no
2-torsion, so this proves that `P_n`, and hence the original point, is
outside the generic rational span.

For an explicit dependence witness, retain integer words

\[
W_n=\sum_{j<n}2^jT_j,\qquad P=2^nP_n+W_n.
\]

Here `W_n` denotes either the word or its generic group value. If
`P_a=P_b` with `a<b`, then

\[
\boxed{(2^{b-a}-1)P=2^{b-a}W_a-W_b\in M.}
\tag{2}
\]

The multiplier is odd. Zero gives an integral word directly. A cycle is
therefore a proof of rational dependence, not a heuristic stall.

**Proof of termination.** The Mordell--Weil group is finitely generated.
In its real height space, if `L` bounds the norms of the finite set of
binary words `T_n`, then

\[
\|P_{n+1}\|\le\tfrac12(\|P_n\|+L).
\]

Consequently every point in an indefinitely continuing chain has bounded
canonical height. The Mordell--Weil lattice is discrete and its torsion is
finite, so there are only finitely many such rational points. An exact point
must repeat. This argument proves termination without computing a height,
regulator, saturation index or an effective height bound in the algorithm.

## 4. The precise arithmetic meaning of the halving depth

**New deduction.** Put `A=E(Q)/M` and let `A_free=A/A_tors`. Its torsion has
odd order by (1). For an independent candidate, the number of successful
halves before escape is exactly

\[
\nu_2(\overline P)=\max\{n:\overline P\in2^n A_{\rm free}\}.
\tag{3}
\]

Indeed the step succeeds if and only if the class in `A` is divisible by2:
one direction follows from the doubling identity, and the other from the
unique compatible generic parity. Odd torsion imposes no extra obstruction.

For **any** split quadratic multisection with trace `Z in M`,

\[
Q+Q'=Z,\quad R=Q-Q'=2Q-Z,
\qquad\overline R=2\overline Q.
\]

Thus

\[
Q\notin M\otimes\mathbf Q\iff R\notin M\otimes\mathbf Q,
\qquad\nu_2(\overline R)=\nu_2(\overline Q)+1.
\tag{4}
\]

This predicts the extra halving layer on every productive control. The
inherited2-Kummer image of `R` is caused by multiplication by2 in the
residual quotient; it does not imply that `R` is inherited. On a dependent
split the residual class is odd torsion or zero, and the process cycles
instead. Changing from an M16 core to full M17 changes the reference quotient:
that is precisely why the nine generic recoveries switch classifications.

## 5. A small, independently checkable302 obstruction

On a short model `y^2=x^3+Ax+B`, for a nonzero target `D=(a,b)`, the abscissa
of a rational half must solve

\[
g_D(X)=X^4-4aX^3-2AX^2-(8B+4Aa)X+A^2-4Ba=0.
\tag{5}
\]

This follows by clearing denominators in the duplication formula. Factor
over `Q`, test each rational root for a rational ordinate, and verify the
sign by doubling; this is a finite degree-four computation, not point search.

**Verified application.** The generic-only302 footprint has exactly17
rows and rank17. Therefore it cannot distinguish *any* candidate from a
generic parity using those rows alone. For branch0 of the already calibrated
carrier, its unique generic parity word is

```
(0,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0).
```

For `D=Q-T`, the primitive integer polynomial (5) reduces modulo191 to

\[
\boxed{15X^4+56X^3+181X^2+107X+109.}
\]

It has no root in `F191`, and its leading coefficient is nonzero. Thus its
homogenization has no projective root over `F191`, excluding any rational
half of `D`. Together with the generic17 injection this is an exact rank18
certificate. The full integer polynomial, all191 nonzero evaluations and
the exact target are retained in the replay package.

The rank21 seed and small-conic seed have analogous root-free certificates
at157 and137. Their anti-traces yield such a certificate after one half.
The `2953/1671` frame instead already has an extra quotient row and detects
the branch by a finite escape. No new point or full-rank computation is used.

**Important arithmetic distinction.** The obstruction above concerns the
finite four-point torsor of halves of a **given rational point**. It is not
a nonsolubility claim for a genus-one2-covering, and not a nonzero Sha class.
The Kummer classes of the actual rational points remain rational Kummer
images, hence have zero image in Sha. The decision concerns their position
relative to the generic rational-span image and their residual divisibility.

## 6. Prospective use and unresolved selection

For a specified degree-two RR construction, at a smooth parameter with valid
maps and a passed generic footprint, the complete pointwise criterion is:

\[
\boxed{\text{residual quadratic splits rationally}
\quad\text{and its halving chain escapes rather than cycles}.}
\]

This is algorithmic, not an existence test referring to an unknown full
Mordell--Weil basis. The quadratic roots and their elliptic maps are explicit;
each escape or cycle has a finite exact certificate. The bounded version
preserves unresolved outcomes. Nonsplitting only rejects this multisection's
branches, and failure of the generic-footprint gate is not seed absence.
Degenerate coefficients, singular fibres and map poles require separate
charts and are not silently assigned a rank.

This can classify future conic-produced candidates without reading known
exceptional points, catalogue ranks, extra-prime rank searches or search
landscape scores. **No running factory or V3 policy has been modified.**

The remaining problem is constructive member selection: the conic misses302,
while the useful302 genus-one member was retrospectively calibrated. This
theorem removes the ambiguity between productive and inherited **specified
splits**; it does not yet supply the generic-only rule that selects a useful
cover through302. The broader seed-incidence goal therefore remains open.

## Checkpoints and reproduction

The [protocol](../../artifacts/generated-results/elliptic-curves/det1092_split_descent_v1/protocol.json)
freezes the five old parameters, nine old successful arms, generic-only prime
rule, eight-layer cap and source hashes. Each generic frame is saved before
candidate construction. Every completed classification is immutable.

The [independent descent replay](../../artifacts/generated-results/elliptic-curves/det1092_split_descent_v1/independent-replay.json)
uses Sage's finite elliptic groups rather than the producer's elementary
integer kernel. It checks exact doubling and odd-cycle relations. It replaces
the nonhalving factorization claims by elementary root-free modular proofs,
with a separately fixed proof-prime cap257. The
[independent incidence replay](../../artifacts/generated-results/elliptic-curves/det1092_split_descent_v1/incidence-replay.json)
reconstructs every generic17 prefix and checks every RR branch/anti-trace and
old generic recovery against source equations. Neither imports the classifier.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_split_descent.sage
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_split_incidence.sage
timeout 25s sage -python -m unittest discover \
  -s research/elliptic-curves/tests -p test_split_seed_descent.py
```

Eight focused tests include a genuine odd-index cycle, an independent point
hidden behind two successful halves, exact half reconstruction, cap semantics,
frame tampering, absence of artifact reads during classification, and basis
sign changes. The synthetic curve is a unit-test fixture, not a new research
fibre. All research invocations complete in seconds under25-second caps.
No point search, parameter sweep, class group, Selmer census, pilot mutation
or detached process was started.
