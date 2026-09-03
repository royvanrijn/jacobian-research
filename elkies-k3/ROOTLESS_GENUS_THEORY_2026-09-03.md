# Rootless existence as genus theory

<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS 2f5b874c0c22133b -->

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

## 4. Foundry integration

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

## 5. Proof boundary and next calculation

Proved here:

- the `<2` sufficient criterion and its quantitative rootless-mass bound;
- the finite ADE mass-inversion criterion `(3.1)--(3.3)`;
- exact first moments for determinants 78, 948, and 950;
- exact agreement of the determinant-78 local formula with its complete
  class census.

Not yet computed:

- the locally admissible ADE lists for the determinant-948 and 950 genera;
- any higher-degree representation average for those genera;
- their exact rootless masses `mu_0`;
- a generic implementation of the Chenevier--Taibi classifier for these
  composite discriminants.

The next bounded task is to build the locally admissible root-system list and
the smallest low-rank representation matrix that separates `mu_0` for the
three controls.  Determinant 78 is the calibration target: every partial or
full inversion must return `mu_0=0` before it is trusted on 948 or 950.

The exact replay artifact is
[`elkies-k3-rootless-genus-first-moment-v1.json`](../artifacts/generated-results/elkies-k3-rootless-genus-first-moment-v1.json).
