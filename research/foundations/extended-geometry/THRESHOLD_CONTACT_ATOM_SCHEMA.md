# Threshold contact atoms: an abstract component theorem

This note isolates the contact-atom mechanism from the weighted Keller map
that first produced it.  The theorem is conditional only in the useful sense:
it lists the geometric inputs a reconstruction problem must supply, and then
deduces its components, dominant dimensions, quasipolynomial behavior, and
abc separation formally.  No Keller map occurs in the statement.

Fix an algebraically closed field `k` of characteristic zero and an integer
`r>=2`.  Put

\[
 S_r=\{r,r+1,r+2,\ldots\},\qquad
 A_r=\{r,r+1,\ldots,2r-1\}.                         \tag{1}
\]

For partitions, write `alpha <=_c lambda` when `lambda` is obtained from
`alpha` by merging parts.  Thus `alpha` is the separated contact type and
`lambda` is one of its collision types.

## 1. The threshold reconstruction package

A **threshold-`r` reconstruction package** in degree `n` consists of a
parameter space `X_n`, a failure locus `F_n`, and locally closed strata
`E_lambda` indexed by partitions `lambda` of `n`, subject to the following
hypotheses.

1. **Exact threshold.** A marked reconstruction chart is singular at a root
   precisely when that root has multiplicity at least `r`.  Complete
   reconstruction failure therefore has the disjoint stratification

   \[
   F_n=\bigsqcup_{\substack{\lambda\vdash n\\\lambda_i\ge r}}E_\lambda.
                                                               \tag{2}
   \]

2. **Collision law.** Every merger of marked roots occurs in a one-parameter
   collision, and no other specialization changes the multiplicity type:

   \[
   E_\lambda\subseteq\overline{E_\alpha}
   \quad\Longleftrightarrow\quad
   \alpha\le_c\lambda.                              \tag{3}
   \]

3. **Geometric recovery.** Every `E_lambda` in (2) is nonempty and
   irreducible, its root-incidence model is finite onto its image, and no
   component is added at infinity in the chosen parameter space.

The last hypothesis may be replaced by any problem-specific result that gives
irreducibility of the stratum closures and rules out hidden boundary
components.  It is stated in the form used by coincident-root incidences.

## 2. Threshold contact-atom theorem

> **Theorem (threshold contact atoms).** In a threshold-`r` reconstruction
> package, the irreducible components of `closure(F_n)` are exactly
>
> \[
> C_\alpha=\overline{E_\alpha},\qquad
> \alpha\vdash n,\quad \alpha_i\in A_r.             \tag{4}
> \]
>
> Every other failure stratum is a collision stratum in at least one of these
> components.  Consequently the number of primitive components is
>
> \[
> c_r(n)=[t^n]\prod_{m=r}^{2r-1}\frac1{1-t^m}.       \tag{5}
> \]

**Proof.** No integer in `A_r` is a sum of two elements of `S_r`.  Conversely,
every `m>=2r` splits as `r+(m-r)`, with both summands in `S_r`.  Hence `A_r`
is exactly the set of atoms of the additive contact semigroup.

Apply this splitting repeatedly to every part of `lambda`.  It produces an
atomic partition `alpha` with `alpha <=_c lambda`; the collision law puts
`E_lambda` in `C_alpha`.  An atomic partition admits no proper threshold
refinement, so its closure is maximal among the stratum closures.  Geometric
recovery makes these maximal closures precisely the irreducible components.
Finally, choosing the number `a_m` of parts of each atomic size gives
`sum(m*a_m)=n`, whose generating function is (5).  QED.

This separates the universal implication from its hypotheses.  Merely finding
a locus of polynomials with all root multiplicities at least `r` does **not**
establish a component theorem; (3) and geometric recovery still have to be
proved.

## 3. Dominant components are an integer program

Suppose the incidence calculation gives an affine-linear dimension formula

\[
 \dim E_\alpha=\delta(n)+
 \sum_{m=r}^{2r-1}w_m a_m,
 \qquad \sum_{m=r}^{2r-1}m a_m=n,\quad a_m\in\mathbb Z_{\ge0}.  \tag{6}
\]

Then the dominant components are exactly the optimizers of the one-row
integer program

\[
 \max\left\{\sum_{m=r}^{2r-1}w_m a_m:
       \sum_{m=r}^{2r-1}m a_m=n, a_m\ge0\right\}.               \tag{7}
\]

For the common root-incidence formula `dim E_alpha=length(alpha)-delta`, all
weights equal one.  If `n=qr+s`, with `0<=s<r`, the partition

\[
 (r+s,r,\ldots,r)                                   \tag{8}
\]

has `q=floor(n/r)` parts and is atomic.  No threshold partition can have more
than `q` parts.  Thus

\[
 \max\ell(\alpha)=\left\lfloor\frac nr\right\rfloor,\qquad
 \dim F_n=\left\lfloor\frac nr\right\rfloor-\delta.             \tag{9}
\]

If the ambient dimension `d(n)` is a quasipolynomial, the codimension

\[
 d(n)-\left\lfloor\frac nr\right\rfloor+\delta                  \tag{10}
\]

is a quasipolynomial whose period divides the least common multiple of `r`
and a period of `d(n)`.  More generally, a fixed finite-weight program (7) has
an eventually quasi-linear optimum on residue classes.  Equation (5) also
makes `c_r(n)` a quasipolynomial of degree `r-1`, with period dividing
`lcm(r,r+1,...,2r-1)` and leading term

\[
 \frac{n^{r-1}}{(r-1)!\,r(r+1)\cdots(2r-1)}.                     \tag{11}
\]

This is the one-dimensional vector-partition phenomenon; see
[Milev's exposition and algorithm](https://arxiv.org/abs/2302.06894).

## 4. Mason separation above the critical threshold

Let `P,Q` have degree `n`, with every root of each polynomial having
multiplicity at least `r`.  Suppose

\[
 P-Q=R\ne0,\qquad \deg R\le d,                                \tag{12}
\]

and `P,Q,R` are pairwise coprime.  Write `ell(P)` for the number of distinct
roots.  Mason--Stothers gives

\[
 n\le\deg\operatorname{rad}(PQR)-1
 \le\ell(P)+\ell(Q)+d-1
 \le2\left\lfloor\frac nr\right\rfloor+d-1.                   \tag{13}
\]

Therefore separation is automatic whenever

\[
 n-2\left\lfloor\frac nr\right\rfloor-d+1>0.                  \tag{14}
\]

For affine pencils, `d=1`.  Every `r>=3` then gives a strict contradiction for
every admissible `n`.  At `r=2`, (13) is borderline: all-double types attain
the support bound, and the excess/factorization argument from the weighted
theory is genuinely necessary.  For constrained pencils of difference degree
`d>1`, (14) gives the exact eventual range rather than an informal claim that
larger thresholds are always enough.

## 5. Where to look next

The following audit distinguishes a natural threshold detector from a proved
threshold reconstruction package.

| Setting | Threshold detector | What is already available | Missing package input | Assessment |
|---|---|---|---|---|
| Higher derivative reconstruction | `p(a)=...=p^(r-1)(a)=0` iff `mult_a(p)>=r` | An exact local threshold and an additive zero-divisor under collision | A global reconstruction map whose chart is regular whenever one of these jets is nonzero | Strong local model |
| Marked jets / Hermite data | A marked `(r-1)`-jet vanishes exactly at an `r`-fold root | Confluent Vandermonde matrices organize arbitrary derivative data; see [Niu--Zhang--Zhou](https://arxiv.org/abs/2207.01852) | Ordinary confluent determinants detect collisions of nodes, not automatically `r`-full root divisors; the incidence must impose jet vanishing | Best explicit laboratory, not yet an application |
| Higher ramification incidence / Lyashko--Looijenga maps | The ramification divisor records `e_x-1`, and these orders add when ramification points collide | LL maps already stratify Hurwitz spaces by ramification data; see [Lando--Zvonkine](https://arxiv.org/abs/math/0303218) | Prove that failure of the chosen inverse/marking map is exactly the `r`-full ramification locus and that its compactification has no extra components | Strongest moduli candidate |
| Iterated tangent maps | An osculating construction could force simultaneous vanishing through order `r-1` | The repository's first tangent map realizes `r=2` | Naive composition gives a union of earlier critical divisors, not their simultaneous intersection; an osculating rather than iterated map is needed | Weak without redesign |
| Constrained polynomial pencils | If two exceptional members differ in a degree-`d` constraint space, (13) applies verbatim | Multiple-root loci are partition-indexed ([Lee--Sturmfels](https://arxiv.org/abs/1508.00202)); polynomial pencils meet these strata ([Borcea--Shapiro](https://arxiv.org/abs/math/0404215)) | Construct a finite marked-jet recovery map and prove exact threshold, collision, and irreducibility | Best route to the Mason/quasipolynomial package |

The most economical first target is therefore a constrained osculating pencil:
choose a finite-dimensional space `V` of low-degree differences, mark a root
together with its first `r-1` derivatives, and seek a reconstruction formula
regular exactly off their simultaneous vanishing.  If its incidence map is
finite and its collision strata are nonempty and irreducible, the threshold
contact-atom theorem immediately supplies a non-Keller family of component
and codimension results.

## 6. Application checklist

Before declaring a new instance, verify all six items.

1. Identify the reconstruction denominator or Fitting ideal.
2. Prove its marked-root zero set is exactly `mult>=r`, in both directions.
3. Prove complete failure means **every** root meets the threshold.
4. Construct every merger by a flat one-parameter collision and exclude other
   closure relations.
5. Prove nonemptiness, irreducibility, finite recovery, and boundary
   exhaustion for every relevant stratum.
6. Record the dimension weights `w_m` and the degree `d` of differences before
   invoking the optimization and Mason corollaries.

## 7. Executable certificate

Run:

```bash
python scripts/verify_contact_atom_principle.py
```

Besides the original `r=2` checks, the script now verifies atomic component
counts, the dominant-length formula, the bounded-difference Mason margin, and
the quasipolynomial finite-difference identities for thresholds through five.
