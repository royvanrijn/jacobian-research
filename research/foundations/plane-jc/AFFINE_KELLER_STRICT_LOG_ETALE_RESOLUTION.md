# Affine Keller strict-log-étale resolution theorem

> **Status.**  Let `Phi:S->T` be a plane Keller map in characteristic zero
> and let `C` be any reduced affine target curve.  Every embedded resolution
> of `(T,C)` pulls back to an embedded resolution of `(S,Phi^(-1)C)`, and the
> resolved map of logarithmic pairs is strict étale.  Its relative
> logarithmic cotangent cokernel, all its Fitting ideals, and every localized
> Chern correction are therefore zero.  Ordinary normalization and conductor
> lengths still pull back fiberwise, but they are invariants of the curve,
> not finite-length pieces of the relative logarithmic module.  In
> particular, the `k=1` cusp-plus-three-node carrier discriminant cannot by
> itself supply the missing boundary `Fitt_1/ch_2` class.

Exact chain-rule and blowup-chart regressions for a node, cusp, tacnode, and
ordinary triple point are checked by
[`verify_affine_keller_strict_log_etale.py`](../scripts/verify_affine_keller_strict_log_etale.py).

## 1. Statement

Let `k` be an algebraically closed characteristic-zero field, let

\[
 \Phi:S=\mathbb A_k^2\longrightarrow T=\mathbb A_k^2
\]

have nonzero constant Jacobian, and let `C` be a reduced curve in `T`.
Choose an embedded resolution by point blowups

\[
 \rho:(\widetilde T,D_T)\longrightarrow(T,C),      \tag{1.1}
\]

where `D_T` is the reduced total transform and is SNC.  Form the Cartesian
square

\[
\begin{array}{ccc}
 \widetilde S=S\mathbin\times_T\widetilde T&\xrightarrow{\widetilde\Phi}
   &\widetilde T\\
 \downarrow&&\downarrow\rho\\
 S&\xrightarrow{\Phi}&T.
\end{array}                                      \tag{1.2}
\]

Put `D_S=widetildePhi^(-1)D_T` with its reduced structure.

### Theorem 1.1 -- affine strict-log-étale resolution

The following hold.

1. `widetilde S` is smooth, and its projection to `S` is the sequence of
   point blowups obtained by pulling back the centers in (1.1).  It is an
   embedded resolution of `Phi^(-1)C`, with reduced SNC total transform
   `D_S`.
2. `widetildePhi` is étale and the map of divisorial log pairs
   `(widetilde S,D_S)->(widetilde T,D_T)` is strict.
3. Pullback of logarithmic differentials is an isomorphism:

   \[
   \widetilde\Phi^*\Omega^1_{\widetilde T}(\log D_T)
   \mathrel{\cong}
   \Omega^1_{\widetilde S}(\log D_S).             \tag{1.3}
   \]

   Consequently

   \[
   \boxed{
   \mathcal T^{\log}_{\widetilde\Phi}=0,
   \qquad \operatorname{Fitt}_i
      (\mathcal T^{\log}_{\widetilde\Phi})=\mathcal O,
   \qquad ch_2^{\rm loc}=0.}                     \tag{1.4}
   \]

This holds for every reduced affine plane-curve singularity and every
number of affine points above it.

## 2. Proof

A Keller map is étale on the affine plane.  In particular it is flat, and
étaleness is preserved by base change.  Blowup commutes with flat base
change.  Applying these facts at every center of (1.1) shows inductively
that `widetilde S->S` is the corresponding sequence of pulled-back point
blowups and that `widetildePhi` is étale.  Each pulled center is a finite
étale collection of smooth points.  Thus `widetilde S` is smooth and `D_S`
is the étale pullback of an SNC divisor, hence is SNC.

The divisorial log structure of `D_S` is the pullback of that of `D_T`, so
the resolved morphism is strict.  A strict étale morphism is log-étale.  In
local resolved coordinates `(r,s)` with `D_T` equal to `r=0` or `rs=0`, the
pulled coordinates `(r circ widetildePhi,s circ widetildePhi)` are étale
coordinates on `widetilde S`; in the corresponding logarithmic bases the
matrix in (1.3) is the identity.  This proves (1.3)--(1.4).

The same argument also shows directly why no singularity type is special.
A target cusp may require more blowups than a node, but every blowup chart
and every exceptional divisor is simply reproduced by étale base change.
There is no relative logarithmic matrix left on which a cusp correction
could live.

## 3. Ordinary conductor survives, but in a different module

Let `z` be a singular point of `C`, with delta invariant `delta_z`, and let
`N_z` be the number of affine source points over it.  Normalization and the
conductor commute with étale base change.  Therefore the pulled curve still
has ordinary normalization-quotient length

\[
 N_z\delta_z                                      \tag{3.1}
\]

and conductor-divisor degree `2*N_z*delta_z`.  These numbers are useful:
they constrain finite fibers and record how much of a degree-`d` fiber must
escape to the boundary.  But (1.4) proves that they are not lengths of
\(\mathcal T_{\widetilde\Phi}^{\log}\) and cannot be inserted as the finite sheaf in a Chern
filtration of that relative module.

For the F2 `k=1` curve, the all-stratum affine delta remains four and the
conductor degree remains eight.  Those statements are unchanged.  Their
correct role is ordinary curve/fiber bookkeeping, not an automatic
localized logarithmic correction.

## 4. Separation from carrier-parameter ramification

The raw carrier-jet map is a morphism from a coefficient-parameter space to
a jet space.  Its Jacobian determinant, including the factor
`Res(p',q')`, lives in that coefficient ring.  The relative logarithmic
matrix in (1.3) lives on the resolved source surface.  These objects are on
different spaces and there is no canonical Fitting-ideal identification
between them.

On the generic `k=1` nonimmersion divisor the target curve has one ordinary
cusp and three nodes, and the raw seven-jet map has corank one.  Theorem 1.1
shows what happens after an actual Keller pullback: every affine preimage of
that cusp-plus-node packet is strict log-étale on the pulled embedded
resolution.  Hence the carrier rank loss is not the missing source-surface
point class.  It can affect the target stratum and conductor distribution,
but additional boundary geometry is required before it can enter the
global obstruction module.

## 5. Remaining calculation

All possible nonzero terms are now confined to the compactification
boundary, where the extended map is not étale.  At each unresolved boundary
attachment one must:

1. extract the actual completed pullbacks of two target SNC parameters;
2. compute their exponent matrix and logarithmic unit first jets;
3. apply the rank and first-jet gates of the
   [`tame-node theorem`](F2_AFFINE_K1_TAME_NODE_PACKET.md); and
4. only for a surviving rank-one packet, compute the higher determinant,
   normalized branch module, and signed finite quotient.

For a unibranch value there is already one exact conditional endpoint: the
[`unibranch attachment theorem`](LOG_UNIBRANCH_ATTACHMENT_FITTING.md) gives
local point length `q_p*m_C` when the boundary attachment is minimally ramified,
transversely unimodular, reduced, and SNC.  What remains is to prove that an
actual F2 attachment has those properties—or to compile the higher packet
identified by their failure.  Over a complete residue-degree-`f` fiber of
minimal attachments the total is `m_C*f`; `q_p` is a local residue
ramification index and is not the transverse divisorial index `e`.

Thus neither the four-node target packet nor its cusp degeneration needs a
separate affine logarithmic calculation.  The missing fixed-coordinate F2
Laurent pair and its unresolved boundary factorization remain the upstream
input.  This theorem does not exclude `(75,125)` or prove `JC(2)`.

## Sources

- The Stacks Project, [étale morphisms](https://stacks.math.columbia.edu/tag/02GH),
  for stability under base change.
- The Stacks Project, [blowing up](https://stacks.math.columbia.edu/tag/085P),
  for compatibility of blowups with flat base change.
- K. Kato, [*Logarithmic structures of Fontaine--Illusie*](https://www.math.brown.edu/dabramov/LOGGEOM/Kato-log.pdf),
  especially the chart criterion for log-étaleness.

## Reproduction

```bash
.venv/bin/python scripts/verify_affine_keller_strict_log_etale.py
```
