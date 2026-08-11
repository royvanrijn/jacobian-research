# F2 outgoing terminal-tail theorem

> **Status.**  This note compiles the entire boundary tail leaving the
> terminal divisor through its `s=0` endpoint.  In carrier coordinates
> `q=y_old`, `r=x*y_old^5-1`, the forced terminal endpoint terms give the
> unimodular exponent map `(a,b)=(q,q^-2*r)` up to units and a target
> translation.  It sends the complete source fan
> `(5,12),(2,5),(1,3),(0,1)` to the already extracted target fan
> `(5,2),(2,1),(1,1),(0,1)`.  Every outgoing boundary node and the final
> boundary/nonboundary endpoint is log-etale; no further source or target
> blowup, logarithmic cokernel, normalization defect, or localized-Chern
> correction occurs on this tail.  The proof is stable under every unexposed
> lower Laurent term allowed by the certified terminal support halfspaces.
> It does not compile the affine purity row or close the global Chern ledger.

The support inequalities, fan map, intermediate node coordinates, and final
endpoint Jacobian are replayed by
[`verify_f2_outgoing_terminal_tail.py`](../scripts/verify_f2_outgoing_terminal_tail.py).

## 1. Coordinates at the outgoing side

At the principal carrier point put

\[
 q=y_{\rm old},\qquad r=xy_{\rm old}^5-1.         \tag{1.1}
\]

Then

\[
 x=(1+r)q^{-5},\qquad
 X=q^{-1}(1+r)^{1/5},qquad
 y_{\rm tr}=q-X^{-1}.                            \tag{1.2}
\]

Since `y_tr=q*r*unit`, the terminal coordinate satisfies

\[
 s=X^{17}y_{\rm tr}^5=q^{-12}r^5\cdot\text{unit}. \tag{1.3}
\]

The `s=0` endpoint terms of the certified terminal block are nonzero and give

\[
 P=q^{-3}r\cdot\text{unit},\qquad
 R=-Q=q^{-1}\cdot\text{unit}                    \tag{1.4}
\]

on the open outgoing cone.  Near the `R`-dominant target point use, after a
constant target translation of `P` when necessary,

\[
 a=R^{-1},\qquad b=\frac{P-c_0}{R}.              \tag{1.5}
\]

Equations (1.4)--(1.5) have exponent map

\[
 M_{\rm out}=
 \begin{pmatrix}1&0\\-2&1\end{pmatrix},
 \qquad \det M_{\rm out}=1.                     \tag{1.6}
\]

The translation replaces the old target coordinate by `b-c_0*a`, a regular
unipotent shear.  It does not change the target boundary fan.

## 2. Stability under all permitted lower terms

After grouping the full polynomial map as a Laurent polynomial in `q` and a
polynomial in `r`, write a monomial as `q^A r^B`, with `B>=0`.  The terminal
support halfspaces give

\[
 5A+12B\ge-3\quad(P),\qquad
 5A+12B\ge-5\quad(R).                            \tag{2.1}
\]

For the ray `(2,5)`, the `P` order relative to the endpoint `(-3,1)` obeys

\[
5(2A+5B+1)=2(5A+12B+3)+B-1.                    \tag{2.2}
\]

For `B>=1`, (2.2) is nonnegative and equality forces `(A,B)=(-3,1)`;
for `B=0`, (2.1) forces `A>=0`, giving a strictly larger order.  On the next
ray `(1,3)`,

\[
5(A+3B)=(5A+12B+3)+3B-3.                       \tag{2.3}
\]

The only order-zero terms are the forced endpoint `(-3,1)` and a possible
constant `(0,0)`.  The latter is precisely the removable target value `c_0`.

For `R`, the corresponding identities are

\[
\begin{aligned}
5(2A+5B+2)&=2(5A+12B+5)+B,\\
5(A+3B+1)&=(5A+12B+5)+3B.
\end{aligned}                                    \tag{2.4}
\]

They show that `(-1,0)` is the unique leading `R` term on both outgoing
boundary rays.  Thus no unexposed lower band can change (1.6).  At the last
endpoint, higher `B=0` terms vanish with positive `q`-order after subtracting
`c_0`; the coefficient of `q^-3*r` remains nonzero and supplies a simple
tangential parameter.

## 3. Complete fan map

The source rays from the terminal divisor to the strict nonboundary curve are

\[
 (5,12),(2,5),(1,3),(0,1).                      \tag{3.1}
\]

Applying (1.6) gives

\[
 (5,2),(2,1),(1,1),(0,1),                       \tag{3.2}
\]

which are already consecutive rays of the minimal target extraction.  Every
adjacent determinant in (3.1) and (3.2) equals one.  Hence no common fan
refinement is missing.

The first cone `(5,12)|(2,5)` is the previously certified unramified
`s=0` terminal node.  For the middle cone `(2,5)|(1,3)`, regular source
parameters are

\[
 c=q^3/r,\qquad d=r^2/q^5,
 \quad q=c^2d,\quad r=c^5d^3.                    \tag{3.3}
\]

The target cone `(2,1)|(1,1)` has parameters `a/b` and `b^2/a`.  Equations
(1.4) give

\[
 f^*(a/b)=c\cdot\text{unit},\qquad
 f^*(b^2/a)=d\cdot\text{unit}.                  \tag{3.4}
\]

Its logarithmic exponent matrix is therefore the identity modulo the maximal
ideal.  The middle node is log-etale.

## 4. The final boundary/nonboundary endpoint

For the last cone `(1,3)|(0,1)`, take

\[
 \alpha=q,qquad \beta=r/q^3.                    \tag{4.1}
\]

Here `alpha=0` is the last source-boundary component and `beta=0` is the
strict transform of the affine curve `xy_old^5=1`.  The endpoint coefficient
in (1.4) is nonzero, so after setting `c_0=P(0,0)` the full functions have

\[
 a=\alpha\cdot\text{unit},\qquad
 P-c_0=\beta\cdot\text{unit}+\alpha\cdot H.      \tag{4.2}
\]

The implicit-function theorem replaces `beta` by `P-c_0`.  Relative to the
target boundary coordinate `a` and tangential coordinate
`b/a=P-c_0`, the full differential is invertible.  Thus the endpoint is
log-etale even if lower bands change `c_0` or the higher function `H`.

## 5. Consequence and claim boundary

The outgoing source components `(2,5)` and `(1,3)` map to the target
components `(2,1)` and `(1,1)` with transverse index one.  Their common node,
the `s=0` node, and the final point where the boundary meets the affine curve
all have zero logarithmic cotangent cokernel.  The existing `27/48` source
component bounds and the target extraction are unchanged.

This closes the outgoing-tail entry in the F2 boundary ledger.  It does not
construct the separate purity-forced ramification row, compile an unknown
noncyclic global attachment, determine the global degree, or prove that the
exact root contribution `27` cannot be cancelled.  Therefore it does not
exclude `(75,125)` or prove `JC(2)`.

The subsequent affine-purity frontier proves that this separate row requires
at least one new component, so the global source floors become `28/49`.  Its
coarse ledger survives at every remaining degree, leaving the target curve,
pullback factorization, and localized-Chern term open.

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_outgoing_terminal_tail.py
```
