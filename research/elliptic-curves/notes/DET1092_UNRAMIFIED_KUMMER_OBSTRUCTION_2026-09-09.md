# Every globally unramified cubic Kummer lift is inherited

The [subsequent singular-fibre-only theorem](DET1092_BAD_FIBRE_KUMMER_OBSTRUCTION_2026-09-09.md)
uses this result to exclude lifts of the primitive302 seed class ramified
only over singular fibres. The broader unramified theorem below remains
the proof input; its historical scope statements are not a claim that the
follow-up route is still open on the fixed panel.

## Result and exact boundary

**New deduction from established cohomology and verified parent arithmetic.**
Let `K=Q(t)`, let `E/K` be the determinant1092 parent, and use its completed-square
model `Y^2=f_t(X)`. Put `L=K(theta)`, where `f_t(theta)=0`, and let
`g:Gamma -> P1_Q` be the smooth projective normalization of this cubic cover.
Write `delta` for ordinary cubic 2-descent and `M17=E(K)` for the already
certified full generic group. Then

\[
\boxed{
 \{[\alpha]\in L^*/L^{*2}:N_{L/K}(\alpha)\in K^{*2},\quad
       v_z(\alpha)\equiv0\pmod2\text{ for every }z\in\Gamma^{(1)}\}
       =\delta(M_{17}).}
\tag{1}
\]

Every closed point is required, including those over the bad parameter
fibres and infinity. There is **no degree bound**: arbitrary rational
coefficients in `1,theta,theta^2` are allowed. This is an obstruction to
globally unramified lifts of a nongeneric **mod-2 class**, not an obstruction
to all nongeneric rational points or all multisections.

In particular, at a smooth rational parameter `t0`, a class
`delta(P) outside delta(M17(t0))` cannot be the specialization of such an
unramified lift. Specialization means restriction of the extended torsor;
it does not require a chosen rational-function representative to be nonzero
and finite at `t0`.

This is stronger in degree than the earlier
[polynomial-lift obstruction](DET1092_SEED_KUMMER_COVER_2026-09-08.md#polynomial-lift-obstruction-an-entire-family-not-one-copied-abscissa),
but has a stronger unramifiedness hypothesis. That earlier theorem allowed
ramification over bad parameter fibres. Its degree26 opening for that
different class of lifts is **not** closed by (1).

## 1. Inputs already established

**Verified prior applications.** The full parent geometry and the
[surface-Brauer theorem](DET1092_SURFACE_BRAUER_TRIVIALITY_2026-09-09.md)
give the following inputs, without using an exceptional point:

- `X -> B=P1_Q` is a smooth K3 surface with zero section and exactly24
  irreducible nodal fibres; the fibre at infinity is smooth.
- The full integral Picard group is generated over Q by `O`, a fibre `F`,
  and the17 generic sections. Geometric and rational MW groups agree and
  have no torsion.
- Pullback identifies `Br(X)=Br(Q)`. In particular this holds for2-torsion.
  This was proved from the fixed149/151 reductions, not conjectured from
  a bounded absence of Brauer representatives.

The cubic cover is geometrically connected: otherwise a cubic would have
a root over `Qbar(t)`, giving geometric rational2-torsion. It is simply
ramified at the24 nodal parameters and unramified at infinity. Consequently
its genus is10 by Riemann--Hurwitz. No Jacobian group computation is needed.

## 2. Identify the unramified norm kernel

**Established ingredients, explicit application.** Let `j:U -> B` be the
smooth-fibre open. The Weil pairing gives the generic split exact sequence

\[
0\longrightarrow E[2]\longrightarrow
   \operatorname{Res}_{L/K}\mu_2\xrightarrow{N}\mu_2\longrightarrow0.
\]

The diagonal splits the norm because the cover has odd degree. It extends
over the entire base as

\[
0\longrightarrow j_*E[2]\longrightarrow g_*\mu_2
   \xrightarrow{N}\mu_2\longrightarrow0.
\tag{2}
\]

Here is the bad-fibre check that cannot be omitted. At an `I1` parameter,
inertia interchanges two cubic roots. In additive `F2` notation, the stalk
of `g_*mu2` consists of triples `(a,a,b)`. The norm is `b`, and its kernel
is `(a,a,0)`, exactly the inertia-invariant part of the even-parity
permutation module `E[2]`. At a good parameter the kernel consists of all
even-parity triples. These checks also show that the diagonal splitting
persists at the nodal fibres. Infinity is a good-fibre stalk in minimal
coordinates.

Finite proper `g` has no higher direct images of `mu2`, by proper base
change and the cohomology of its finite geometric fibres. Thus (2) gives

\[
 H^1(B,j_*E[2])=
 \ker\bigl(H^1(\Gamma,\mu_2)\xrightarrow{N}H^1(B,\mu_2)\bigr).
\tag{3}
\]

The Kummer sequence on a regular curve identifies `H1(Gamma,mu2)` with
the squareclasses in `L` whose valuations are all even. The field norm
agrees with the sheaf norm. Also `H1(B,mu2)=Q*/Q*2` injects into `K*/K*2`.
Therefore the right side of (3) is exactly the left side of (1), denoted `A`.

## 3. Why every such class comes from a generic section

**Established literature.** For this smooth minimal Weierstrass surface,

\[
 R^1\pi_*\mu_2\simeq j_*E[2].
\tag{4}
\]

One reference is the proof of Proposition2.7 in
[Feng--Landesman--Rains](https://link.springer.com/article/10.1007/s00208-022-02429-1).
Their squarefree-discriminant setting applies here in characteristic zero:
both sheaves identify with the torsion Neron sheaf. In particular (4) includes
the nodal stalks; it is not merely a statement on the smooth open.

**New deduction.** Consider the Leray spectral sequence for `pi` and `mu2`.
The term `H1(B,R1 pi_*mu2)=A` has no incoming differential. Its only possible
outgoing differential is

\[
 d_2^{1,1}:A\longrightarrow H^3(B,\mu_2).
\]

The zero section splits pullback `H3(B,mu2) -> H3(X,mu2)`, so no class in
that base term can be killed. The differential is therefore zero. Hence
every `a in A` is the `(1,1)` component of a class `b in F1 H2(X,mu2)`,
where `F1` is the kernel of restriction to geometric fibre cohomology.

The Kummer sequence and `Br(X)=Br(Q)` now write

\[
 b=c_1(D)+\pi^*c,\qquad D\in\operatorname{Pic}(X),\quad
 c\in H^2(B,\mu_2).
\tag{5}
\]

To see that constants really can be removed at this level, lift the
constant Brauer image of `b` through `H2(Q,mu2) -> Br(Q)[2]`, and subtract
its pullback. The remainder is a divisor class by Kummer exactness.
The constant contribution is in `F2` and has zero `(1,1)` component.
Because `b` has zero geometric fibre degree, `D.F` is even. Subtracting
`(D.F)O` makes this degree zero without changing its mod-2 class.

Since `Pic(X)` is generated by `O,F,S1,...,S17`, the resulting degree-zero
divisor is an integral combination of `F` and `Si-O`. The fibre contributes
zero to `(1,1)`. The contribution of `Si-O` is precisely `delta(Si)`.
For this last compatibility, on the smooth base take local halves of the
relative degree-zero line bundle `O(Si-O)`; their differences are the
`E[2]` cocycle of the torsor `2Q=Si`. This is both the Kummer--Leray
component of its divisor class and the usual elliptic Kummer class.
Under the Weil-pairing identification (2), its generic representative is
`X(Si)-theta`.

It follows that `A` is contained in `delta(M17)`. Conversely the same
divisor construction puts the Kummer class of every generic section in
`H1(B,j_*E[2])=A`. This proves (1).

The cohomological tools used here are the Kummer sequence, proper base
change and Leray; see also
[Colliot-Thelene--Skorobogatov, sections2.2 and3.1](https://www.imo.universite-paris-saclay.fr/~jean-louis.colliot-thelene/BGgroup_book.pdf).
No identification of a rank-two Picard surface's Brauer group with a
spectral Jacobian is used. No generic Selmer-dimension calculation or
specialized Sha vanishing is asserted.

## 4. A sharper interpretation: the missing information is ramification

**New deduction.** Define `Nker=ker(L*/L*2 -> K*/K*2)`. The parity divisor
induces an injection

\[
 \boxed{N_{\rm ker}/\delta(M_{17})\ \hookrightarrow\
      \bigoplus_{z\in\Gamma^{(1)}}\mathbf F_2,\qquad
 [\alpha]\longmapsto(v_z(\alpha)\bmod2)_z.}
\tag{6}
\]

Indeed two norm-square classes with identical residue patterns have an
everywhere-even quotient; (1) says that quotient is inherited. Conversely
generic Kummer multiplication never changes the pattern. We do not claim
every finite pattern occurs. Principal-divisor, norm and rationality
conditions constrain the image of (6).

This gives a useful exact formulation for future construction: a primitive
seed class needs a nonzero ramification pattern on the generic spectral
curve, followed by a specialization at which the class is rational Kummer.
Nonzero generic ramification is necessary for such a lift, **not sufficient**
for a seed, a soluble specialized covering, or independence.

## 5. Why the completed controls behave differently

**Verified application, retrospective only.** In the previously frozen
first-seed continuation set

\[
 k=\frac{3192654717569013009686535431767390554676}{1234321},\quad
 D(t)=f_t(k),\quad \alpha=D(t)(k-\theta).
\]

The generic cubic is

\[
 f_t(X)=X^3+5X^2+(16a_4(t)+8)X+64a_6(t)+16.
\]

Exact polynomial replay gives `N(alpha)=D^4`; `D` has degree12, is
squarefree and is coprime to the discriminant. Above each geometric root
of `D`, the three valuations of `alpha` are `(2,1,1)`. Thus its generic
ramification is real, at24 points over12 good parameter values. At infinity
the minimal coordinate is `X'=t^-4 X`; `alpha` has pole order16 on each
branch, so it adds no ramification. This is consistent with (1), not an
exception to it.

At302 (`t=0`), `D(0)` is a rational square and the class is exactly the
certified first-seed Kummer class. On each of the eight unchanged controls,
the old factor-free certificates find a nonsquare remainder coprime to all
bad-prime support. Therefore an odd valuation at a good arithmetic prime
obstructs Selmer membership of this **particular continuation**. It is not
a Sha class there and cannot have a rational covering point. This explains
the observed split for the calibrated family; copying `k` from302 is still
an oracle, and no prospective discriminator has been obtained.

**Verified control distinctions reused, not new point computations.**

| Old control | Class relative to full generic MW17 | Consequence of (1) |
|---|---|---|
| First302 branches; first rank21 seed; two other productive conic seeds | Outside generic mod2 image, zero halving layers | No globally unramified lift of their branch class |
| Their anti-traces `R=2P-Z` | Inherited mod2 class; one halving layer to escape | An inherited lift exists despite independence of `R` |
| Nine blinded M16 recoveries, viewed over M17 | Inherited rational span; exact halving cycles | Inherited mod2 lift; missing only from the M16 reference |
| Dependent conic split | Inherited rational span; exact cycles | Inherited mod2 lift; no extra direction |

The anti-trace row is essential: `delta(R)=delta(Z)` does **not** imply
`R` is generic. The completed
[halving-or-cycle theorem](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md)
supplies the exact primitive-class distinction, under its verified mod2
injection hypotheses. Different covers or squareclasses alone remain
insufficient for a rank claim.

## 6. Checkpoint, open goal and limits

**Verified computation.** The accompanying checker verifies the pinned
prior inputs,24I1 hypotheses, good infinity, the two stalk calculations,
the calibrated norm and residue regression, all eight factor-free remainder
obstructions, and the exact classification fields of the38 old controls.
The prior point-independence and halving certificates are reused with their
source hashes; they are not presented as newly recomputed points. The
cohomological proof above is a written proof, not a formal CAS verification.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_unramified_kummer_obstruction.sage
```

The [checkpoint packet](../../artifacts/generated-results/elliptic-curves/det1092_unramified_kummer_v1/)
has an explicit25-second cap, no point searches, no new parameter values,
no integer factorization and no group census.

**Open constructive goal.** Selecting a useful ramified class or base-change
member from generic equations alone, and proving that its specialization
produces a seed outside generic MW17, remains open for302. The theorem
closes the *globally unramified* route, not ramified lifts, classes on open
surfaces, new base changes, bad-fibre-only ramification, or individual
rational fibre points. No production search or later amplification artifact
was used or modified, and no detached job was started.
