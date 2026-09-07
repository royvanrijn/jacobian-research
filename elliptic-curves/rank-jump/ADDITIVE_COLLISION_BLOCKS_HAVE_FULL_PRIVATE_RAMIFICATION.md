# The tested additive collision blocks have full private ramification

The first additive constructor beyond the exhausted multiplicative radical
pool produces **680 independent ramified directions on every frozen MW17
fibre and 560 on every MW16 fibre**. Its entire span has zero intersection
with Selmer. Adding the generic subgroup still produces no new Selmer class.

This holds on all sixteen retained high/low fibres, including the fresh
rank-27 examples and historic +12/+14 controls. It is a block exclusion:
**every nonzero product of the tested classes remains ramified at an odd
prime of good reduction**. The conclusion does not merely exclude each
generator separately.

The experiment introduces addition among generic Kummer radicals, which
is not covered by the preceding
[multiplicative radical-closure theorem](ALL_IN_FIELD_ROOTS_OF_THE_RETAINED_POOL_ARE_EXHAUSTED.md).
It therefore tests a different possible source of classes. It supplies
no new unramified class or sufficient specialization condition.

## A symmetric three-radical constructor

For a frozen equation E: y²=f(x)=x³+Ax+B, let K=Q(θ), f(θ)=0. Its marked
generic sections give α_i=x_i−θ. For every unordered triple i,j,k, put

\[
 q_{ijk}(\theta)=(\alpha_i+\alpha_j-\alpha_k)^2-4\alpha_i\alpha_j
 =-3\theta^2+2s_1\theta+s_1^2-4s_2,
\tag{1}
\]

where s1=x_i+x_j+x_k and s2=x_ix_j+x_ix_k+x_jx_k. This expression is
symmetric in the triple. It is the Heron polynomial in the three α's:

\[
 \prod_{\epsilon\in\{\pm1\}^3}
  (\epsilon_1\sqrt{\alpha_i}+\epsilon_2\sqrt{\alpha_j}
                         +\epsilon_3\sqrt{\alpha_k})=q_{ijk}^2.
\tag{2}
\]

Thus q measures an additive collision of three radicals. Its zero locus
is where some signed sum vanishes. Equation (2) is verified as a formal
polynomial identity, not inferred from numerical values.

Clear rational scalar content to make q a primitive integral polynomial
of degree two, and form

\[
       a_q=N_{K/\mathbb Q}(q)q,
       \qquad N(a_q)=N(q)^4.
\tag{3}
\]

This is a cubic norm-square Kummer class before local conditions. Rational
rescaling of q changes a_q by a fourth power, so the primitive normalization
does not change its class. Translation of all x_i and θ leaves (1)
unchanged; rational Weierstrass scaling changes it by a rational square.
The construction still depends on the specified marked generic sections,
not just the unmarked elliptic curve.

There are C(17,3)=680 triples in a MW17 row and C(16,3)=560 in a MW16
row. The [first protocol](THREE_RADICAL_INCIDENCE_PROTOCOL.json) includes
**every** such triple on **every** frozen row. No exceptional point,
outcome-dependent triple choice, new parameter, or point search is used.

## A factorization-free ramification obstruction

Let Δ=disc(f), and let q∈Z[θ] be primitive of degree below three. At an
odd prime p not dividing Δ, the cubic order is maximal and étale over Zp.
Write v_l for the valuations of q at its factors, and

    n = v_p(N(q)) = sum_l f_l*v_l.

The valuation of a_q at the l-th factor is v_l+n. If n is odd and all
these valuations were even, every v_l would be odd. Then q would belong
to p O_K at every factor over p. Because the power basis is a Zp-basis,
all coefficients of q would be divisible by p, contradicting primitivity.
Therefore

\[
 p\nmid2\Delta,\quad v_p(N(q))\text{ odd}
 \quad\Longrightarrow\quad [a_q]\text{ ramified at a prime above }p.
\tag{4}
\]

At an odd good prime every elliptic Selmer class is unramified. Hence
(4) is an **incidence exclusion before Selmer membership**, not a CT
obstruction distinguishing rational classes from Sha.

To certify the existence of such a prime, no factorization is needed.
Repeatedly divide |N(q)| by its gcd with Δ. A nonsquare remainder has an
odd prime exponent outside Δ. For the Selmer conclusion one must also
exclude p=2. The independent verification explicitly removes two from
every final private factor below and checks that the odd part remains
nonsquare. This matters on case05, whose integral cubic discriminant is
odd. All 10,520 witnesses pass this stronger odd-prime test.

## Why no combination cancels the ramification

An individual exclusion would leave the main block question unanswered.
The [second protocol](THREE_RADICAL_PRIVATE_RAMIFICATION_PROTOCOL.json)
therefore tests cancellation across each fibre's **entire** triple space.

Let n_i be the part of |N(q_i)| left after removing Δ-supported factors.
Starting with n_i, repeatedly divide by its gcd with the product of all
other n_j. Call the result P_i, and remove its power of two. The
certificate proves

    P_i divides |N(q_i)|,
    gcd(P_i, 2*Δ*product_(j != i) N(q_j)) = 1,
    P_i is not a square.

There is consequently an odd prime p_i at an odd norm exponent belonging
to q_i and to no other q_j. By (4), a_qi has a nonzero valuation-parity
coordinate at some prime ideal above p_i. Every other a_qj is a unit
there. These coordinates give a diagonal identity matrix on the candidate
classes, even though the primes need not be factored or explicitly named.

If H is their Kummer span and n the number of triples, this proves

\[
 \boxed{\dim H=n,\qquad
 \operatorname{rank}(\text{outside ramification}|_H)=n,\qquad
 H\cap\operatorname{Sel}_2(E)=0.}
\tag{5}
\]

For the generic subgroup G, whose classes are unramified at these odd
good primes, it also proves

\[
 (G+H)\cap\operatorname{Sel}_2(E)=G.
\tag{6}
\]

In fact (6) holds with G replaced by any subspace of Selmer. Thus this
is not an artefact of missing generic corrections or a chosen basis.
Each fibre has at least n distinct private odd rational primes supporting
these obstructions. The actual prime factorization and full support size
are not computed.

## Paired results

| Frozen fibre | Generic rank | Retained gain | Candidate dimension | Proved outside ramification rank | Added Selmer dimension |
|---|---:|---:|---:|---:|---:|
| 074d9, 2818/1535 | 17 | +10 | 680 | 680 | 0 |
| 074d9, 2824/885 | 17 | 0 | 680 | 680 | 0 |
| 103b2, 3726/881 | 17 | +10 | 680 | 680 | 0 |
| 103b2, -1049/2296 | 17 | 0 | 680 | 680 | 0 |
| 11952, -2448/11 | 17 | +10 | 680 | 680 | 0 |
| 11952, -1171/1683 | 17 | 0 | 680 | 680 | 0 |
| 11952, 110314/102227 | 17 | +10 | 680 | 680 | 0 |
| 11952, 130349/28916 | 17 | 0 | 680 | 680 | 0 |
| 11952, 2012/211 | 17 | +10 | 680 | 680 | 0 |
| 11952, 2828/2015 | 17 | +10 | 680 | 680 | 0 |
| 11952, 4286/1881 | 17 | +10 | 680 | 680 | 0 |
| MW16 a1-01, -1867/270 | 16 | +11 | 560 | 560 | 0 |
| MW16 a1-01, -3187/3697 | 16 | +1 | 560 | 560 | 0 |
| ICARM356 | 17 | +12 | 680 | 680 | 0 |
| ICARM385 | 17 | +12 | 680 | 680 | 0 |
| ICARM398 | 16 | +14 | 560 | 560 | 0 |

The labels come from the existing
[frozen manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_governing_panel_manifest_v1.json)
after arithmetic. Low gains remain censored outcomes, not exact ranks;
the later public rank-28 bound is not substituted for case06's retained
rank-27 bound. The same-family matching and exposure limitations are
unchanged.

The 103b2 pair, both 11952 pairs, the 074d9 pair and the MW16 pair all
give the same exact failure. The historic controls do too. For this
constructor, a large jump does not coincide with a small collision-defect
rank: the rank is maximal everywhere.

The private odd cofactors range from 398 to 1942 bits across the panel.
These are composite support witnesses, **not asserted prime values**.
The proof uses their coprimality and nonsquareness, which are exact.

## Controls and interpretation of the specialization condition

All 269 inherited linear Kummer representatives pass the necessary
outside-norm square test. Their actual Selmer membership is already
known from the generic points. All sixteen derivative representatives
f'(θ) pass too. Their norm projection is the
[four-division derivative class](FOUR_DIVISION_RAMIFICATION_CANNOT_SUPPLY_THE_JUMP.md),
which fails Selmer locally on every row. These deliberate false positives
show why passing the norm test is not sufficient.

For a family with generic coordinates x_i(t), the construction provides
explicit functions

\[
 R_I(t)=\operatorname{Res}_X\bigl(f_t(X),q_I(t,X)\bigr).
\]

At t0 first replace q_I(t0,X) by its primitive integral representative
and let N_I(t0) be its integer norm. These normalized norms, rather than
the raw rational values R_I(t0), are essential: removing rational content
changes the norm by a cube. Define the private part P_I(t0) of |N_I(t0)|
by removing primes supported on 2Δ(t0) or on any other N_J(t0). A necessary
condition for a word of these classes to become unramified is

\[
     c_I=1\quad\Longrightarrow\quad P_I(t_0)\in\mathbb Q^{\times2}
     \quad\text{for each participating }I.
\tag{7}
\]

Equivalently, any surviving odd norm valuation must be shared with
another participating class or supported at the excluded bad places.
This is an explicit necessary condition on specialization for this
constructor. It is not a sufficient condition, an independently solved
low-degree carrier, or a rank predictor. Even removing all private
obstructions leaves factor-by-factor parity and elliptic local conditions
to check. Global rational-versus-Sha solubility comes later.

On every frozen t0, (7) forces all c_I=0. These particular additive
collision covers therefore cannot carry the successful jump block.

## Consequence for the next construction

The useful lesson is more specific than another failed correlation:
arbitrary additive norms generate many new standard S3 directions in
ambient Kummer cohomology, but they bring private good-prime ramification.
Neither taking products nor correcting by rational generic directions
removes it here. A useful class-creation construction needs arithmetic
that forces those new divisors into even valuations or shared cancellable
support at the specialization.

The next high-value gate is therefore a **symbolic factorization or
shared-divisor identity in a proposed family of class generators**, before
another specialization panel. A promising identity would remove the
private branch factors in advance and leave a small explicit parity
kernel. Its specialization would then need an exact S-split unramified
certificate and independence from the old global pool. Nothing in this
result licenses a larger search over radical coefficients or parameters.

Ranked conclusions are: (1) the new S-split unramified standard class
module remains the necessary incidence mechanism; (2) unweighted
three-radical additive collisions are excluded as its source throughout
this panel; (3) a shared-divisor construction, independent class
representatives and their rational-solubility theorem are still missing.
Agent 1 gets a capacity gate for evaluating future mechanisms, not a new
candidate score. No visibility feature is introduced.

## Evidence and replay

The [compact manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_three_radical_evidence_v1.json)
summarizes all sixteen rows. The
[compressed evidence](../../artifacts/generated-results/elliptic-curves/rank_jump_three_radical_evidence_v1.zip)
preserves the two full JSON certificates byte for byte: every primitive
quadratic, exact norm, discriminant gcd sequence and private cofactor.
The [independent verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_three_radical_incidence_verification_v1.json)
checks the ZIP hashes, reconstructs the generic-only equations and
coefficients, and recomputes all norms using multiplication-matrix
determinants rather than resultants. It checks every private odd witness
against the product of the other norms and verifies the formal identity (2).

All **10,520** triple norms, **285** controls and **10,520** private odd
ramification witnesses passed. No norm was factored; no class group,
new number field or rational point was computed. Each construction worker
has a 30-second cap and a checkpoint. The full independent replay takes
seconds locally and reads the compressed bundle directly:

```sh
timeout 60 python3 elliptic-curves/rank-jump/verify_three_radical_incidence.py check
```

`package_three_radical_evidence.py unpack` restores the original JSONs
without overwriting changed files. The construction scripts' `capture`
modes regenerate them in a clean output location; the package tool's
`package` mode writes the deterministic ZIP and manifest. All outputs
are immutable. Only new rank-jump-specific files are committed; Agent 1's
active search and mathematical-status entries remain untouched.
