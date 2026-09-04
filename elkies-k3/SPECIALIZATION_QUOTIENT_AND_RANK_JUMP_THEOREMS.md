# Specialization quotients and rank-jump theorems

This note is the canonical arithmetic companion to
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
The fibration-hopping note controls changes of a marked elliptic fibration on
one K3 surface.  The statements below instead control the arithmetic quotient
between a generic Mordell--Weil group and one specialized elliptic fibre.

They separate three logically different objects:

```text
certified points  ->  specialization quotient  ->  residual Selmer space.
```

The first gives a rank lower bound, the middle is the actual rank jump when
the generic rank and its specialization are certified, and only a complete
Selmer calculation gives the last, upper-bound side.  The height statements
describe how to search the quotient.  The Kummer statement explains a
class-group obstruction forced by already known points; it is not a Selmer
upper bound or a prospective rank predictor.

## 1. Setup

Let `K` be a number field, let `E_eta/K(T)` be an elliptic curve, and let
`E_t/K` be a smooth specialization.  Let `M_eta` be the certified full free
part of the generic Mordell--Weil group and let

```text
M = sp_t(M_eta) subset E_t(K).
```

Throughout the rank-jump statements, assume that `rank(M)=rank(M_eta)=r`.
This is an input: it can follow from an injective specialization certificate
or from an independent rank calculation on the displayed specialized
generators.  Let `L` be any certified finitely generated subgroup with

```text
M subset L subset E_t(K).
```

Ranks mean free abelian ranks.  Torsion is retained in the group quotients
but does not contribute to rank.

## 2. The rank-jump sandwich

### Theorem S1: certified quotient below, residual Selmer above

For every prime `ell`, let

```text
delta_ell:E_t(K)/ell*E_t(K) -> Sel_ell(E_t/K)
```

be the Kummer injection, and let `bar(M)_ell` be the image of `M` under this
map.  Then

\[
 \boxed{
 \operatorname{rank}(L/M)
 \;\leq\;
 \Delta_t
 \;\leq\;
 \dim_{\mathbf F_\ell}
 \bigl(\operatorname{Sel}_\ell(E_t/K)/\overline M_\ell\bigr),
 }
 \tag{S1}
\]

where

\[
 \Delta_t=
 \operatorname{rank}E_t(K)-\operatorname{rank}E_\eta(K(T))
 =\operatorname{rank}(E_t(K)/M).
\]

#### Proof

The left inequality follows from `L subset E_t(K)` and the equality of the
generic and specialized ranks of `M`.  For the right inequality, write
`A=E_t(K)`.  The image of `A/ell*A` modulo the image of `M` is

\[
 A/(M+\ell A).
\]

For a finitely generated abelian group this vector space has dimension at
least `rank(A)-rank(M)`; finite `ell`-primary saturation defects can only make
it larger.  The Kummer map injects it into
`Sel_ell(E_t/K)/bar(M)_ell`, proving the upper bound. QED.

Thus an exact Smith computation of `L/M` certifies the lower side even when
`L` is not the full Mordell--Weil group.  Conversely, a local-signature span,
an incomplete descent, or a bounded point search does not certify the upper
side.  If the rank of `M` after specialization has not been proved to remain
`r`, then `rank(L/M)` is a **quotient gain**, not a certified jump from the
generic rank.

The theorem concerns arithmetic ranks of fibres `E_t/K`.  It is distinct
from the Shioda--Tate identity

```text
Delta(MW rank) = Delta(Picard rank) - Delta(fibre-root rank)
```

for specialization in a family of elliptic surfaces, and from the fixed-K3
rank balance under a change of fibration.

## 3. Quotient height

### Definition S2: the real Mordell--Weil quotient metric

Put

\[
 V=(E_t(K)/E_t(K)_{\rm tors})\otimes_{\mathbf Z}\mathbf R,
 \qquad U=M\otimes_{\mathbf Z}\mathbf R,
\]

and equip `V` with the Neron--Tate pairing.  If `pr_U` is orthogonal
projection to `U`, define

\[
 \langle \bar P,\bar R\rangle_{/M}
 =\langle P-\operatorname{pr}_U(P),
          R-\operatorname{pr}_U(R)\rangle
\]

on `V/U`, and

\[
 \boxed{
 \widehat h_{/M}(\bar P)
 =\inf_{u\in U}\widehat h(P+u)
 =\widehat h(P-\operatorname{pr}_U(P)).
 }
 \tag{S2}
\]

If `m_1,...,m_r` is a basis of `M`, `G=(<m_i,m_j>)`, and
`b=(<m_i,P>)`, then

\[
 \widehat h_{/M}(\bar P)=\widehat h(P)-b^{\mathsf T}G^{-1}b.
\]

For several quotient generators this is exactly the Schur-complement Gram
matrix.  The definition is independent of bases and of the chosen lift in
`V`.  It is not the generally nonquadratic minimum among integral
representatives of a fixed coset of `M` in `E_t(K)`.

## 4. Midpoints and orthogonal decomposition

### Theorem S3: midpoint decomposition for a pointed degree-two chart

In a short Weierstrass model

\[
 E_t:y^2=x^3+Ax+B,
 \qquad Q=(x_Q,y_Q),
\]

the pointed coordinate

\[
 t_Q(R)=\frac{y_R+y_Q}{x_R-x_Q}
\]

satisfies

\[
 t_Q(R)^2=x(R)+x(Q-R)+x(Q),
\]

and gives the quartic

\[
 w^2=t^4-6x_Qt^2-8y_Qt-3x_Q^2-4A.
\]

Let `Q in M`.  This pointed quartic chart has fibre involution

```text
R -> Q-R
```

and is canonically centered at `Q/2 in U`.  For `R in E_t(K)`, put
`R_U=pr_U(R)`.  Then

\[
 \boxed{
 \frac14\widehat h(2R-Q)
 =\widehat h(R-Q/2)
 =\widehat h_{/M}(\bar R)
  +\widehat h(R_U-Q/2).
 }
 \tag{S3.1}
\]

Equivalently, the two points in the chart fibre satisfy

\[
 \boxed{
 \widehat h(R)+\widehat h(Q-R)
 =\frac12\widehat h(Q)
  +2\widehat h_{/M}(\bar R)
  +2\widehat h(R_U-Q/2).
 }
 \tag{S3.2}
\]

#### Proof

The line of slope `t_Q(R)` through `-Q` and `R` has third intersection
`Q-R`.  Comparing the `x^2` coefficient after substituting that line in the
Weierstrass equation gives the displayed identity for `t_Q(R)^2`; eliminating
the two `x`-coordinates gives the quartic, and changing its ordinate sign
exchanges `R` with `Q-R`.

The vectors `R-pr_U(R)` and `pr_U(R)-Q/2` are orthogonal, which gives
(S3.1) by the quadraticity of the canonical height.  Applying the
parallelogram law to

```text
R = (R-Q/2)+Q/2,
Q-R = -(R-Q/2)+Q/2
```

and then substituting (S3.1) gives (S3.2). QED.

This is the exact separation used by half-lattice searches.  The quotient
term measures genuinely new real directions.  The second term measures
misalignment inside the already known real span.  A deep midpoint can exclude
old lattice points from a search region, but it does not create a new
arithmetic direction there.

## 5. Half-lattice depth and covering radius

### Corollary S4: midpoint holes are a discrete covering-radius certificate

For a positive-definite Mordell--Weil lattice `M` and `c in M/2M`, define

\[
 \mu_2(c)=\min_{Q\in c}\widehat h(Q).
\]

For any representative `Q` of `c`, the squared distance from its midpoint to
the old lattice is

\[
 \boxed{
 D_2(c)=\min_{P\in M}\widehat h(P-Q/2)
       =\frac14\mu_2(c).
 }
 \tag{S4.1}
\]

Consequently the `2`-torsion-restricted covering radius is

\[
 \boxed{
 \rho_2(M)=\max_{c\in M/2M}\sqrt{D_2(c)}
 =\frac12\sqrt{\max_c\mu_2(c)}.
 }
 \tag{S4.2}
\]

If `rho(M)` is the full Euclidean covering radius of `M`, then

\[
 \rho(M)\geq\rho_2(M).
 \tag{S4.3}
\]

#### Proof

Since `2P-Q` runs through the class `c` as `P` runs through `M`, quadraticity
gives (S4.1).  Taking the maximum over the finite group `M/2M` gives (S4.2).
The points `Q/2 mod M` form only the two-torsion subset of the real torus
`(M tensor R)/M`, so their maximum distance is bounded above by the supremum
over the whole torus, proving (S4.3). QED.

Thus a complete parity-coset census is an exact discrete covering-radius
calculation and a lower bound for the full covering radius.  On the published
R17 lattice, the complete census has `max mu_2=12`, hence

```text
rho_2(R17)=sqrt(3),  and therefore rho(R17) >= sqrt(3).
```

This is an old-point exclusion statement, not a success probability for a
rank-jump search.

## 6. Everywhere-even Kummer classes

The companion
[`RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md`](RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md)
develops the global subspace inside the residual Selmer quotient, its
Cassels--Tate obstruction, and the distinction between chart labels and
cover classes. It also gives the exact eleven-fibre soluble-cover replay.

### Lemma S5: the ideal square-root map and its unit loss

Let `F` be a number field and let `W` be a finite-dimensional
`F_2`-subspace of `F^*/F^{*2}`.  For a class represented by `alpha`, let

\[
 v([\alpha])=(\operatorname{ord}_{\mathfrak p}(\alpha)\bmod 2)_{\mathfrak p}
\]

be its finite-prime valuation-parity vector.  On the everywhere-even kernel
there is a well-defined homomorphism

\[
 \boxed{
 c:\ker(v|_W)\longrightarrow\operatorname{Cl}(F)[2],
 \qquad (\alpha)=\mathfrak a^2\longmapsto[\mathfrak a].
 }
 \tag{S5.1}
\]

Its kernel consists of unit squareclasses that lie in `W`.  If all classes
under consideration have positive square norm, the kernel is contained in
the norm-positive unit squareclasses.  Write their dimension as `u_F^+`.

Now let `G subset B subset F^*/F^{*2}` be finite-dimensional subspaces of
positive-square-norm classes, all supported modulo squares on a declared
finite set `S` of prime ideals.  Let `v_S` record valuation parity on `S`, and
put

\[
 \Pi(B,G)=
 \frac{c(\ker(v_S|_B))}{c(\ker(v_S|_G))}.
\]

Then

\[
 \boxed{
 \dim_{\mathbf F_2}\Pi(B,G)\geq
 \dim(B/G)
 -\bigl(\operatorname{rank}v_S(B)-\operatorname{rank}v_S(G)\bigr)
 -u_F^+.
 }
 \tag{S5.2}
\]

For an odd-degree field, `Norm(-1)=-1`, so the norm-sign map on unit
squareclasses is onto and Dirichlet's unit theorem gives

\[
 u_F^+=r_1+r_2-1.
 \tag{S5.3}
\]

#### Proof

Even valuations give a fractional ideal `a` with `(alpha)=a^2`; changing
`alpha` by a square changes `a` by a principal ideal.  Hence (S5.1) is
well-defined and lands in the two-torsion of the ideal class group.  Its
kernel is represented by units.  Rank--nullity gives

\[
 \dim\ker(v_S|_B)-\dim\ker(v_S|_G)
 =\dim(B/G)-\bigl(\operatorname{rank}v_S(B)
                         -\operatorname{rank}v_S(G)\bigr).
\]

Passing through `c` can lose at most `u_F^+` further dimensions, which proves
(S5.2).  Finally, the full unit squareclass space has dimension `r_1+r_2`;
in odd degree `-1` maps to the nontrivial norm sign, leaving the kernel
dimension in (S5.3). QED.

The lemma gives a lower bound for a known-point-forced quotient of the full
ideal class-group `2`-torsion.  It computes neither the full class group nor
an `S`-class group: after the primes in `S` are inverted, their ideal classes
may kill some or all of these directions.  It also gives no Selmer upper
bound.  When `B` contains exceptional point classes, the exceptional points
and their quotient gain are inputs, so the result is explanatory rather than
prospective.

## 7. Certified applications and boundaries

The general statements above organize the existing exact applications:

- [`ELKIES_RANK_JUMP_FINGERPRINTS.md`](../elliptic-curves/notes/ELKIES_RANK_JUMP_FINGERPRINTS.md)
  computes displayed lower quotients and their Schur-complement height
  lattices without promoting them to full specialized Mordell--Weil groups;
- [`HALF_LATTICE_HEIGHT_COMPRESSION_MECHANISM_2026-09-04.md`](../elliptic-curves/notes/HALF_LATTICE_HEIGHT_COMPRESSION_MECHANISM_2026-09-04.md)
  proves the pointed-quartic involution and midpoint identities and audits the
  corresponding search mechanism;
- [`QUOTIENT_GEOMETRY_TABLE_2026-09-04.md`](../elliptic-curves/notes/QUOTIENT_GEOMETRY_TABLE_2026-09-04.md)
  tabulates the displayed quotient geometry and every independent recovery
  event for the R17 controls, refreshed ladder, and A1/MW16 hits; its exact
  subspace comparison rules out intrinsic quotient height as a sufficient
  scalar detector;
- [`R17_MULTISECTION_DIVERSITY_2026-09-02.md`](R17_MULTISECTION_DIVERSITY_2026-09-02.md)
  supplies the complete R17 parity-coset census;
- [`R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md`](R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md)
  applies Lemma S5 to six fibres, with the specialized generic `MW17` as `G`.

None of these statements turns a bounded search miss into point absence,
identifies a displayed lower quotient with the full specialization quotient,
or replaces an unconditional complete Selmer calculation on the upper side
of (S1).

<!-- status-consumer: EC-SPECIALIZATION-QUOTIENT-RANK-JUMP-THEOREMS e68a8f4b00720de3 -->
