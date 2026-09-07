# ADE mass as an asymptotic neighbour distribution

<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS e7589727ca8f7e50 -->

## 1. The invariant

For a positive even genus `G`, put

```text
Phi(L) = {x in L : x^2=2},
m(G)   = min_[L in G] #Phi(L).
```

Thus `m(G)=0` is exactly the existence of a rootless isometry class in the
genus.  This is finer than determinant.  In the rank-17 controls relevant to
the foundry,

- determinant at most 28 is excluded by the Hermite bound;
- the determinant-78 genus has 1,549 classes and all are rootful;
- the determinant-948 and determinant-950 genera contain explicit rootless
  classes.

The useful local surrogate is not the integer `m(G)` itself but the
**rootless mass**.  Write

```text
mass(G) = sum_[L in G] 1/|O(L)|,
mu_0(G) = sum_[L in G, Phi(L)=empty] 1/|O(L)|.
```

All summands are positive, so

```text
m(G)=0  if and only if  mu_0(G)>0.
```

Siegel representation averages determine `mu_0(G)` from local genus data by
a finite triangular calculation.  The first row of that calculation is a
particularly cheap sufficient test.

## 2. The first-root-moment gate

Normalize

```text
theta_L(q) = sum_(x in L) q^(x^2/2),
Theta_G(q) = mass(G)^(-1) sum_[L in G] theta_L(q)/|O(L)|.
```

Then

```text
a_1(G) := [q]Theta_G
        = mass(G)^(-1) sum_[L in G] #Phi(L)/|O(L)|.
```

### Theorem 2.1: cheap rootless-mass certificate

For every positive even genus,

```text
mu_0(G)/mass(G) >= 1-a_1(G)/2.
```

In particular,

```text
a_1(G)<2  implies  m(G)=0.
```

#### Proof

Every root occurs with its negative.  A nonempty `Phi(L)` therefore has at
least two elements.  If `p_0=mu_0(G)/mass(G)`, the weighted mean satisfies
`a_1(G)>=2(1-p_0)`.  Rearrangement gives the inequality. QED.

This is one-sided.  A mean at least two says nothing about existence: a small
positive rootless mass can be hidden among rootful classes having many roots.

### Exact rank-17 local formula

Let `L` have rank `2k+1=17`, determinant `d`, and let

```text
D = fundamental_discriminant((-1)^k 2d),
chi_D = Kronecker character of D.
```

For `k=8`, the Siegel product used by the checker is

```text
a_1(G) = 2^(2k) |B_(k,chi_D)|
         -------------------------------- sqrt(2|D|/d)
                    |D|^k |B_(2k)|

         * product_(p divides 2d)
             alpha_p(1) /
             ((1-p^(-2k))/(1-chi_D(p)p^(-k))).
```

Here `alpha_p(1)` is the exact local representation density of `1` for
`Q(x)=x^2/2`.  Hence every quantity in the formula is determined by the local
genus.  The script
[`certify_rootless_genus_first_moment.sage`](scripts/certify_rootless_genus_first_moment.sage)
evaluates it exactly.

The three controls give:

| determinant | exact `a_1(G)` | decimal | known `m(G)=0`? | `<2` gate |
| ---: | ---: | ---: | :---: | :---: |
| 78 | `2913380886349/59299224796` | 49.13016816613631 | no | inconclusive |
| 948 | `7957563723128755857618/562456712956783562285` | 14.14786869783556 | yes | inconclusive |
| 950 | `4967763637986279936/352882035745379473` | 14.07768924108891 | yes | inconclusive |

For determinant 78, the complete 1,549-class census independently gives

```text
mass(G) = 1463420154787/4131952105881600,
sum #Phi(L)/|O(L)|
        = 2013146192467159/115694658964684800.
```

Their ratio is exactly the local-density value above.  This is both a check
of conventions and a useful negative control.  The determinant-948 and
determinant-950 rows show that the first moment is too coarse for the current
positive controls.

## 3. The exact local decision: ADE mass inversion

The correct completion of the first-moment idea is the prescribed-root-system
mass method used by King for unimodular lattices.  The argument itself is not
restricted to the unimodular genus.

Let `R_0=0,R_1,...,R_s` be all ADE root lattices of rank at most `rank(G)`, or
just a locally certified superset of the root systems that can occur in `G`.
Order them by increasing number of roots.  Define

```text
mu_j = sum_[L in G, <Phi(L)> isometric to R_j] 1/|O(L)|,
A_i  = sum_[L in G] r(L,R_i)/|O(L)|,
```

where `r(L,R_i)` is the number of isometric embeddings `R_i -> L`.  An
embedding sends roots to roots, so it factors through the full root lattice
of `L`.  Consequently

```text
A_i = sum_j r(R_j,R_i) mu_j.                         (3.1)
```

Put `U_(i,j)=r(R_j,R_i)`.  If `U_(i,j)` is nonzero, then
`#Phi(R_i)<=#Phi(R_j)`.  Equality forces the embedded roots to be the entire
root system and hence `R_i` and `R_j` to be isometric.  Thus `U` is upper
triangular in the chosen order, with

```text
U_(i,i)=|O(R_i)|.
```

It is therefore invertible over the rationals.  Siegel's degree-`rank(R_i)`
weighted representation theorem expresses each normalized `A_i/mass(G)` as
a product of local representation densities.  Hence

```text
(mu_0,...,mu_s)^t = U^(-1) (A_0,...,A_s)^t            (3.2)
```

is determined by local genus data, and

```text
m(G)=0  if and only if  (U^(-1)A)_0>0.                (3.3)
```

This is the exact genus-theoretic decision procedure.  The `A1` row is the
first-root moment from Section 2.  Higher rows are not merely heuristic
moments: they count compatible ADE configurations and separate the masses of
the possible full root systems.  King's paper, [*A mass formula for
unimodular lattices with no roots*](https://arxiv.org/abs/math/0012231),
implements precisely this triangular strategy using Katsurada's Fourier
coefficient formula.

Two cheaper intermediate gates are available before full inversion.

1. Local discriminant and rank restrictions can delete impossible `R_j`.
   If every remaining nonempty root system has at least `b_G` roots, then
   `a_1(G)<b_G` already proves `mu_0(G)>0` and
   `mu_0(G)/mass(G)>=1-a_1(G)/b_G`.
2. A selected set of low-rank representation averages gives a rational
   linear program in the nonnegative unknown masses `mu_j`.  A positive
   lower bound for `mu_0` is again a certificate, even when the full matrix
   is not yet computed.

## 4. The mass-to-neighbour theorem

The ADE inversion and large-prime neighbour statistics are two halves of one
statement.  The marked version needs slightly more notation because a good
prime fixes a target spinor genus, not merely a connected component after all
good-prime displacements are allowed.

Let `V` be positive definite over `QQ`, of dimension `n>2`, let

```text
K subset O_V(A_f)
```

be compact open, and let `X(K)` be Chenevier's finite marked class set.  For
`y in X(K)`, write `Gamma_y` for its marked stabilizer and `s(y) in S(K)` for
its spinor genus.  The mass of every spinor genus is

```text
m_sp(K) = m_K/|S(K)|.
```

Assume that the marked class `y` functorially determines a positive even
lattice `W_y`; in the unmarked application, `W_y` is just the lattice `y`.
Let `R(y)` be the complete ADE root lattice of `W_y`.  For a spinor genus
`tau` define its marked root-system mass by

```text
mu_R(tau) = sum_[y in X(K), s(y)=tau, R(y)=R] 1/|Gamma_y|.
```

### Theorem 4.1: ADE masses are large-prime neighbour frequencies

Fix `x in X(K)` and `a in S_1(K)`, and put `tau=a*s(x)`.  Let
`N_(p,R)(x)` be the number of marked `p`-neighbours of `x` whose completed
root lattice is `R`.  As `p` tends to infinity through good primes satisfying

```text
delta_p=a,
```

one has

```text
N_(p,R)(x)/c_V(p) = mu_R(tau)/m_sp(K) + O(p^(-1/2)).       (4.1)
```

Here `c_V(p)` is the total number of `p`-neighbour lines.  If `n>4`, and in
particular for the rank-17 foundry genera, the error is `O(p^(-1))`.
The compatible primes form a nonempty union of arithmetic progressions and
have Dirichlet density `1/|S_1(K)|`.

#### Proof

Chenevier's Theorem 5.9 gives, for every fixed target class `y` with
`s(y)=tau`,

```text
N_p(x,y)/c_V(p)
  = (1/|Gamma_y|)/m_sp(K) + O(p^(-1/2)),
```

with `O(p^(-1))` for `n>4`.  Sum over the finite set of target classes with
`R(y)=R`.  The main terms sum to `mu_R(tau)/m_sp(K)`, and a finite sum does
not change the order of the error.  The assertion about compatible primes is
Chenevier's Remark 5.11. QED.

For an unmarked genus `G` consisting of one spinor genus,

```text
Gamma_y = O(y),       m_sp(K)=mass(G),       mu_R(tau)=mu_R(G),
```

so (4.1) is exactly

```text
# {p-neighbour lines producing R}
---------------------------------- = mu_R(G)/mass(G) + O(p^(-1))
       # {p-neighbour lines}
```

in rank 17.  Combining this with `(3.2)` gives the promised local-to-search
formula

```text
N_(p,R)(x)/c_V(p)
  = (U^(-1) A)_R/mass(G) + O(p^(-1)).                  (4.2)
```

For `R=0`, the limiting density is the rootless mass fraction.

Two qualifications are essential.

1. If the level class set has several spinor genera, the limit must be taken
   in a fixed compatible prime progression.  Mixing progressions can
   oscillate when their root-system masses differ.
2. The ordinary Siegel ADE inversion computes genus masses.  It supplies the
   marked numerator in (4.1) only when the relevant marked/spinor refinement
   is also computed, or when the one-spinor-genus unmarked specialization
   applies.

Thus the reciprocal score

```text
construction_score_R(G) = mass(G)/mu_R(G)
```

is the asymptotic mean number of uniformly sampled compatible neighbour
lines per `R`-hit.  It is a baseline for a foundry queue, not a claim that
guided lines are independent or that the asymptotic is accurate at a
particular small prime.

The external input is Gaetan Chenevier,
[*Statistics for Kneser p-neighbors*](https://doi.org/10.24033/bsmf.2852),
Theorem 5.9 and Remarks 5.10--5.11.

## 5. Exact three-control calibration

The bounded checker
[`certify_ade_neighbor_mass_score.sage`](scripts/certify_ade_neighbor_mass_score.sage)
now joins the mass and neighbour sides without claiming the still-missing
local inversion.

### Determinant 78: the low-rank LP does not close

There are `3,768` abstract ADE types of rank at most 17, including zero.  The
complete determinant-78 census realizes `621` distinct nonzero types among
its `1,549` classes.  The checker aggregates the exact reciprocal-
automorphism mass of every one of those strata and computes all embedding
averages with source root lattice of rank at most four:

```text
A1, 2A1, A2, 3A1, A1+A2, A3,
4A1, 2A1+A2, 2A2, A1+A3, A4, D4.
```

It then maximizes the possible rootless mass over nonnegative rational
stratum masses.  To make the test as favorable as possible, the candidate
columns are only zero and the 621 root systems known from the complete
census; every genuine locally admissible superset gives a weakly larger
feasible region.  Exact primal and dual LP certificates give:

| representation rows used | maximum allowed `mu_0/mass(G)` |
| --- | ---: |
| `A1` | `22110891977563/25024272863912 = 0.8835778...` |
| all rank at most 2 | `166140250586/4818062014675 = 0.0344828...` |
| all rank at most 3 | `16449900414854628311/581488895115467353000 = 0.0282893...` |
| all rank at most 4 | `11817053639/4286486820968 = 0.00275682...` |

The true census value is `mu_0=0`, but low-rank data do **not** force it.
Consequently the proposed low-rank LP is a useful gate, but it cannot be
required to reproduce the determinant-78 zero before higher rows are added.
This is a stronger negative calibration than the first moment alone.

The checker emits the complete abstract rank-at-most-17 ADE list and the
complete census-realized determinant-78 list.  It does **not** relabel either
as the complete locally admissible list: integral local embedding filters at
`2,3,13` remain to be implemented.

### Rootless construction scores

All three control genera have one proper spinor genus according to the exact
Sage spinor-kernel quotient.  The unmarked specialization of Theorem 4.1
therefore applies.

| determinant | rootless mass information | mass of genus | asymptotic rootless fraction | reciprocal score |
| ---: | ---: | ---: | ---: | ---: |
| 78 | `mu_0=0` exactly | `1463420154787/4131952105881600` | `0` | infinity |
| 948 | `mu_0=3/4` exactly | `77731517730627488307787/925557271717478400` | `694167953788108800/77731517730627488307787 = 8.930328...e-6` | `111977.969...` |
| 950 | `mu_0>=7/4` from four distinct classes | `152400929187535759901875/888534980848779264` | at least `1554936216485363712/152400929187535759901875 = 1.020293...e-5` | at most `98011.049...` |

The determinant-948 value uses the complete two-class rootless Niemeier
classification: the automorphism orders are two and four.  It corrects the
older statement that its exact rootless mass was unknown.  For determinant
950, the three foundry classes have automorphism orders `4,2,2`, and the new
class has order two and a distinct norm-four count, giving the stated lower
bound.  Completeness of that rootless stratum remains open.

### Finite-prime boundary

For odd rank 17,

```text
c_V(p)=1+p+...+p^15.
```

At the first good prime this is `38,146,972,656` lines for determinants 78
and 948 (`p=5`), and `21,523,360` lines for determinant 950 (`p=3`).  The
determinant-78 finite frequency is exactly zero at every good prime because
the entire genus is rootful.  Exact determinant-948 and 950 finite-prime
frequencies have not been enumerated.  They remain a separate computational
benchmark; the asymptotic theorem must not be presented as their value.

The replay artifact is
[`elkies-k3-ade-neighbor-mass-score-v1.json`](../artifacts/generated-results/elkies-k3-ade-neighbor-mass-score-v1.json).

## 6. Foundry integration

The foundry should rank discriminant forms before it generates frame classes:

```text
discriminant form / local genus
        |
        +-- Hermite obstruction                         reject
        |
        +-- exact a_1(G), with local ADE exclusions
        |       |
        |       +-- certified positive rootless mass    prioritize
        |       +-- inconclusive
        |
        +-- selected higher ADE averages / LP
        |       |
        |       +-- certified positive rootless mass    prioritize
        |       +-- inconclusive
        |
        +-- full triangular ADE mass inversion          decide m(G)=0
                |
                +-- mu_0=0                              reject
                +-- mu_0>0                              construct classes
```

Class generation is then reserved for genera already proved
rootless-capable, or for calibration controls.  Once `mu_0>0`, a neighbour
enumerator need only find and mass-close the rootless stratum if the project
requires explicit Gram matrices; it need not discover rootless capability by
blindly traversing the whole genus.

Chenevier--Taibi's 2026 inductive classification machinery is relevant at
that last stage.  It combines efficient inductive neighbour/orbit generation,
strong isometry invariants, automorphism data, and independent mass closure,
and it reaches genera orders of magnitude larger than prior complete
classifications.  Its published scope is specific: rank-29 unimodular
lattices and even ranks at most 28 with prime (half-)determinant at most 7.
The available code is therefore a blueprint, not an off-the-shelf classifier
for the composite determinant-948 and determinant-950 genera.  See
[*Unimodular lattices of rank 29 and related even genera of small
determinant*](https://arxiv.org/abs/2601.19780) and its
[code/data companion](https://olitb.net/pro/uni29/).

## 7. Proof boundary and next calculation

Proved here:

- the `<2` sufficient criterion and its quantitative rootless-mass bound;
- the finite ADE mass-inversion criterion `(3.1)--(3.3)`;
- the marked spinor/level ADE-frequency formula `(4.1)` and its unmarked
  one-spinor-genus specialization `(4.2)`;
- exact first moments for determinants 78, 948, and 950;
- exact agreement of the determinant-78 local formula with its complete
  class census;
- all determinant-78 ADE stratum masses from the census, exact rank-at-most-
  four embedding averages, and exact primal/dual LP optima;
- exact determinant-948 rootless mass `3/4` and the determinant-950 lower
  bound `7/4`, converted into asymptotic neighbour scores.

Not yet computed:

- the complete locally admissible ADE lists for any of the three controls;
- any higher-degree representation average for those genera;
- the local-density ADE inversion for any of the three controls;
- the exact determinant-950 rootless mass `mu_0`;
- exact finite-good-prime rootless-line frequencies for determinants 948 and
  950;
- a generic implementation of the Chenevier--Taibi classifier for these
  composite discriminants.

The next bounded task is to add exact integral local-embedding filters to the
3,768-type list, starting with determinant 78 at `2,3,13`, and then compute
the next ADE rows selected by the LP dual.  Determinant 78 remains the
calibration target: a full local inversion must return `mu_0=0` before it is
trusted on 948 or 950.  The rank-at-most-four LP is now known not to suffice.

The original first-moment replay artifact is
[`elkies-k3-rootless-genus-first-moment-v1.json`](../artifacts/generated-results/elkies-k3-rootless-genus-first-moment-v1.json).
