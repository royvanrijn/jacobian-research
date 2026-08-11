# F2 affine-purity puncture attachment theorem

> **Status.**  Every F2 nonproperness component with normalization degrees
> `(3k,5k)` meets the already extracted target divisor `(5,2)`.  In the
> target coordinates `pi=b^3/a` and `h=b^5/a^2=P^5/(-Q)^3`, its puncture has
> contact order `k` and residue
> `lambda=p_lead^5/(-q_lead)^3`.  Hence the `k=1` target is smooth and
> transverse there.  A formal direct-SNC comparison with a terminal divisor
> would force `lambda=125/729` and `e=3`.  The certified terminal
> neighborhood, however, is already a regular resolved morphism; every
> divisor extracted above it maps to one target point.  Thus it has no
> available direct slot for a divisor dominating the affine target curve.
> The value `125/729` is a target-incidence condition, not an actual source
> attachment or a component-count improvement.  The unconditional source
> floors remain `28/49`.

The toric orders, terminal branch-value comparison, and direct-SNC local
matrix are replayed by
[`verify_f2_affine_purity_puncture_attachment.py`](../scripts/verify_f2_affine_purity_puncture_attachment.py).

## 1. The affine target puncture lies on `(5,2)`

Let

\[
 \nu(t)=(p(t),q(t)),\qquad
 (\deg p,\deg q)=(3k,5k)                         \tag{1.1}
\]

be the normalization of an irreducible nonproperness component `C`, and
write

\[
 p(t)=A t^{3k}+\cdots,\qquad
 q(t)=B t^{5k}+\cdots,qquad AB\ne0.             \tag{1.2}
\]

At the `Q`-dominant target-infinity point use the certified coordinates

\[
 a=(-Q)^{-1},\qquad b=P/(-Q).                    \tag{1.3}
\]

With `u=1/t`, the curve has

\[
 \operatorname{ord}_u(a,b)=(5k,2k).             \tag{1.4}
\]

Its primitive ray is therefore `(5,2)`, independently of `k`.  On the
extracted divisor use

\[
 \pi=\frac{b^3}{a},\qquad
 h=\frac{b^5}{a^2}=\frac{P^5}{(-Q)^3}.           \tag{1.5}
\]

Substitution gives

\[
 \operatorname{ord}_u\pi=k,qquad
 h\longrightarrow
 \boxed{\lambda=\frac{A^5}{(-B)^3}\in k^\times}. \tag{1.6}
\]

Thus the strict transform of `C` meets the smooth interior of the target
divisor `(5,2)`, never either toric node `h=0,infinity`.  Its intersection
multiplicity with the divisor is `k`.  In particular, on `k=1`, `pi` is a
uniformizer on the normalized curve and the intersection is transverse.

For `k>1`, (1.6) identifies the first unresolved target datum precisely: the
order of `h-lambda` and its later characteristic exponents determine the
additional target proximity chain.

## 2. Direct attachment transfers the purity index to the residue map

The following local lemma is independent of F2.  Let a smooth target have
transverse curves

\[
 L=(\pi=0),\qquad C=(z=0),                       \tag{2.1}
\]

and let an SNC source node have components

\[
 T=(x=0),\qquad E=(y=0),                         \tag{2.2}
\]

with `T` dominating `L` with normal index one and `E` dominating `C` with
normal index `e`.  Since the generic point of `E` does not lie over `L`, and
the generic point of `T` does not lie over `C`, divisor pullback gives

\[
 f^*\pi=x\cdot\text{unit},\qquad
 f^*z=y^e\cdot\text{unit}.                       \tag{2.3}
\]

Restricting to `T` shows that its residue map to `L` has ramification index
exactly `e` at `T intersect E`.  Equivalently, the local logarithmic exponent
matrix is

\[
 \begin{pmatrix}1&0\\0&e\end{pmatrix}.           \tag{2.4}
\]

Therefore a purity component with `e>1` cannot attach directly to an
unramified point of a boundary component that dominates `L`.

## 3. Comparison with the terminal `A_6` cover

Every certified principal terminal component dominates `(5,2)` with residue
map

\[
 h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}.        \tag{3.1}
\]

Its branch values and partitions are

\[
 \begin{array}{c|c}
 h=0&(5,1)\\
 h=\infty&(3,3)\\
 h=125/729&(3,1,1,1).
 \end{array}                                     \tag{3.2}
\]

The affine target puncture has `lambda` finite and nonzero, so the first two
values in (3.2) are impossible.  The direct-SNC lemma now gives

\[
 \boxed{
 \text{direct ramified attachment}
 \Longrightarrow
 \lambda=125/729,\quad s=\infty,\quad e=3.}       \tag{3.3}
\]

The three other points over `125/729` are simple and cannot support `e>1`
under a direct SNC attachment.

## 4. The certified terminal neighborhood has no affine extraction

The implication (3.3) is only a compatibility test.  It does not construct
a direct slot.  On the certified source model the point `s=infinity` is
already the node between the terminal divisor and its `(3,7)`-arm neighbor.
In regular parameters `(tau,w)` the smooth-target endpoint has leading form

\[
 \pi=\tau w\cdot\text{unit},\qquad
 h-125/729=w^3\cdot\text{unit}+\tau(\cdots),       \tag{4.1}
\]

and the exact logarithmic cokernel is `R/(w^3)`.  In particular, the map is
regular on a neighborhood of this source node.

The following elementary birational observation closes the apparent slot.
Let `f:X->Y` be a morphism and `g:X'->X` a proper birational morphism.  If a
curve `F` is contracted by `g` to `x`, then

\[
 (f\circ g)(F)=f(x).                              \tag{4.2}
\]

Thus every exceptional divisor of every further blowup centered in the
resolved terminal neighborhood maps to a point.  The same conclusion holds
for an arbitrary competing resolution after passing to a common
resolution.  No such divisor can dominate the one-dimensional affine target
curve `C`.  This applies at `s=infinity` and, more generally, at every
already resolved terminal point.

Consequently the purity-forced divisor is **not** obtained by direct
extraction over the terminal residue cover.  Its attachment must be found at
another, globally unresolved source-boundary locus.  The terminal branch
passport does not determine its transverse index `e`.

This also corrects a tempting but invalid component count: from
`lambda!=125/729` one cannot infer the conditional floors `29/50`.  A global
attachment outside the terminal neighborhood may add only the one component
already required by purity.  The proved unconditional floors are still

\[
 N_{\rm source}^{\rm squarefree}\ge28,
 \qquad
 N_{\rm source}^{\rm double}\ge49.               \tag{4.3}
\]

## 5. Remaining local data

This theorem closes both the target-infinity location and one false source
attachment route.  The remaining `k=1` calculation is now:

1. decide whether the leading-coefficient residue equals `125/729`;
2. evaluate the finite carrier-normalized jet test
   `b=min(ord_u(w|_C),8)` when it does;
3. locate the affine divisor at an unresolved nonterminal source locus and
   compute its `(e,f,E^2)` data and point-supported `Fitt_1` corrections; and
4. factor the implicit quintic pullback to determine which local branch is
   boundary and which branches are affine companions.

The theorem does not decide the leading residue, construct the source
attachment, or prove that its Chern term cannot cancel the extraction-root
contribution `27`.  The subsequent
[`generic k=1 affine-row Chern theorem`](F2_AFFINE_K1_LOG_CH2.md) computes
the divisorial conormal term for arbitrary `(e,f,E^2,b)`, leaving precisely
the source attachment and point corrections.  This note therefore does not
exclude `(75,125)` or prove `JC(2)`.

On the generic cusp face, the later affine-row and cusp-dichotomy theorems
show that exact minimal attachments have incidence-sensitive total
`2f-h+c`, ranging from `f` for smooth unramified folds to `2f` for an
all-node fiber.  This tightens the doubled conditional residual numerators
by `2(2f-h+c)` without selecting a terminal index.

<!-- status-consumer: PF2K1L1 5221f5659fc19729 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_purity_puncture_attachment.py
```
