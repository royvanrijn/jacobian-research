# A cup-product ideal class controls strict Jacobian lifting

The fixed-cubic solubility switch is governed by more than the ordinary
half ideals. On the strict class space, the Jacobian lifting obstruction
has the form
\[
\boxed{\kappa_u(\beta)=(1-u\theta)\cup\beta
 \quad\in\operatorname{Cl}(\mathcal O_{K,S_K})/2.}
\]
This cup product has zero image in the field Brauer group. Its surviving
ideal-class component measures failure to lift **without introducing
ramification outside \(S\)**.

For two strict classes, evaluation of this ideal class by the second
class's unramified character is the CT difference:
\[
\boxed{\Delta_u(\beta,\psi)
 =\chi_\psi\bigl(\kappa_u(\beta)\bigr).}
\]
On the retained anchor-rational space, \(\Delta_u=\operatorname{CT}_{E_u}\).

The existing strict CT matrix therefore detects at least four independent
cup-obstruction values at \(u=-1\). Two detected dimensions project to the
previously certified elementary class factor, and at least two lie in its
complement. These coordinates are inferred from the retained CT values
using the identity above. Explicit norm witnesses, and independent
recovery of the nine target CT bits, remain uncomputed.

This is a **solubility-obstruction mechanism** for the labelled fixed-cubic
comparison. It is not a rational-point construction or a production-family
rank predictor.

## The first cohomological lifting obstruction vanishes

Write \(T=E_0[2]\simeq E_u[2]\) and \(M=J_u[2]\). The glued Jacobian gives
\[
0\longrightarrow T\longrightarrow M\longrightarrow T\longrightarrow0.
\]
The [local Kummer diagram](JACOBIAN_LOCAL_CONDITIONS_AND_CT.md) identifies
the right local condition with \(L_{0,v}\cap L_{u,v}\).

For the cubic algebra \(K\), the odd-degree norm sequence splits:
\[
\operatorname{Res}_{K/\mathbf Q}\mu_2\simeq\mu_2\oplus T.
\]
Thus \(H^2(\mathbf Q,T)\) is a direct summand of
\(\operatorname{Br}(K)[2]\). The Brauer local-global injection implies
\[
\ker\!\left(H^2(\mathbf Q,T)\to
            \prod_v H^2(\mathbf Q_v,T)\right)=0.
\]
Every class in the common Selmer group \(S_C\) lifts locally through \(M\);
its global connecting obstruction is therefore locally zero and hence
zero globally. Consequently
\[
S_C\subset
\operatorname{im}\!\left(H^1(\mathbf Q,M)\to H^1(\mathbf Q,T)\right).
\]

This proves existence of an unrestricted global cohomology lift. It does
not put that lift in \(\operatorname{Sel}_2(J_u)\). The remaining issue is
whether a global correction can impose the required local conditions.
The distinction is the one in
[Morgan–Smith, Section 3.1](https://arxiv.org/pdf/2103.08530).

## The extension cocycle as a projected cup product

Let \(P=\mathbf F_2^3\) be the permutation module on the three cubic roots.
Identify
\[
T=\{(v_1,v_2,v_3):v_1+v_2+v_3=0\},\qquad
\operatorname{pr}(z)=z+\operatorname{Tr}(z)(1,1,1).
\]
In the coordinates used by the
[six-branch-point model](TORSION_DIFFERENCE_AND_CT.md),
\(v_i=e_2(T_i,v)\). Dot product on this even-parity subspace is the
elliptic Weil pairing.

The branch pairs are \(\pm\sqrt{\gamma_i}\), where
\(\gamma_i=1-u\theta_i\). If a Galois element has root permutation \(g\)
and sign vector \(\epsilon\), its extension term is
\[
c(\epsilon)v=\sum_i\epsilon_i e_2(T_i,v)T_i.
\]
In root coordinates this is exactly
\[
c(\epsilon)v=\operatorname{pr}(\epsilon_i v_i)_i.
\]
The new finite certificate checks this identity for every root permutation,
sign vector and element of \(T\): 192 cases.

The connecting map is therefore the projection of the ordinary cup
product over \(K\):
\[
\delta(\beta)=
\operatorname{pr}_*\bigl(\gamma_u\cup\beta\bigr),\qquad
\gamma_u=1-u\theta.
\]
This keeps the actual labelled extension. Replacing it by its abstract
degree-48 Galois group would discard \(\gamma_u\).

## Restrict ramification before taking the obstruction

Fix \(S\) containing 2, infinity and the bad primes of both curves and the
Jacobian. The retained sets have this property: outside primes dividing
\(2uD(u)\operatorname{disc}(f)\), the displayed genus-two model has good
reduction; all those primes are included.

Let
\[
U=\operatorname{Sel}_2^S(E_0)=\operatorname{Sel}_2^S(E_u),\qquad
C=\operatorname{Cl}(\mathcal O_{K,S_K}).
\]
Here strict means zero localization at \(S\). For \(\beta\in U\),
\(\gamma_u\cup\beta\) is zero in every field completion:

- above \(S\), \(\beta\) is square;
- outside \(S\), \(\gamma_u\) is a unit and \(\beta\) is unramified;
- at real places, \(\beta\) is positive.

Hence the field Brauer class is zero. In cohomology with ramification
restricted to \(S_K\), it can still be nonzero:
\[
\gamma_u\cup\beta\in
\ker\!\left(H^2(G_{K,S_K},\mu_2)\to\operatorname{Br}(K)[2]\right)
 \simeq C/2.
\]
Define this class to be \(\kappa_u(\beta)\). The Kummer description of this
kernel and its ideal formula are given by
[McCallum–Sharifi, Section 2 and Theorem 2.4](https://arxiv.org/pdf/math/0202161).

Under the norm projection from \(K\) to \(\mathbf Q\), its ideal-class
component is unchanged, since the rational \(S\)-ideal class group is
zero. The preceding extension calculation thus gives
\[
\boxed{\beta\text{ lifts to }H^1(G_{\mathbf Q,S},M)
       \quad\Longleftrightarrow\quad\kappa_u(\beta)=0.}
\]
The unrestricted lift always exists, but an \(S\)-unramified lift need not.

This obstruction is different from
\(\Phi(\beta)=[\mathfrak J_\beta]\in\operatorname{Cl}(K)[2]\), where
\((\beta)=\mathfrak J_\beta^2\). The old \(\Phi\) depends on \(\beta\);
the new \(\kappa_u\) also depends on the deformation through \(\gamma_u\).

## A concrete norm-and-ideal construction

Put \(F=K(\sqrt{\gamma_u})\) with involution \(\sigma\).
The vanishing field Hilbert symbols imply that \(\beta\) is a global norm:
choose \(z\in F^\times\) with
\[
N_{F/K}(z)=\beta.
\]
This is a norm witness, not a rational point on the elliptic cover.

Outside \(S_K\), the extension \(F/K\) is unramified. At a split prime,
the two valuations of \(z\) have equal parity because their sum is
\(v_{\mathfrak p}(\beta)\), which is even. At an inert prime there is
only one valuation. Therefore there is an \(S_K\)-ideal \(\mathfrak I\)
whose extension to \(F\) records the parity of \((z)\):
\[
(z)=\mathfrak I\mathcal O_{F,S_F}\,\mathfrak A^2
\]
for a fractional ideal \(\mathfrak A\). Then
\[
\boxed{\kappa_u(\beta)=[\mathfrak I]\pmod{2C}.}
\]

For comparison with the cited ideal formula, write
\[
(\beta)=\mathfrak b^2,\qquad
(z)=\mathfrak c^{\,1-\sigma}\mathfrak b\mathcal O_{F,S_F}.
\]
Modulo square ideals, \(1-\sigma=1+\sigma\), so
\[
[\mathfrak I]=[N_{F/K}\mathfrak c]+[\mathfrak b]\quad\text{in }C/2,
\]
which is that formula for \(n=2\).

The choice of norm witness does not change this class. Replacing \(z\)
by \(z\,t/\sigma(t)\) changes the descended parity ideal by the principal
norm ideal of \(t\), modulo squares. Replacing \(\beta\) by a square
multiple likewise changes it by a principal ideal.

This gives an executable target: an exact norm witness, its valuation
parity ideal, then Artin evaluations. It does not require mistaking a
successful norm equation for rational solubility.

## Why strict CT is Artin evaluation of this ideal

The norm witness gives an unrestricted \(M\)-valued cohomology lift:
use its Kummer class in the six-point permutation module, take the
even-parity submodule, and quotient by the all-ones vector.

In the CT formula for a global cocycle lift, compare it with local
Jacobian Kummer lifts and pair the difference with \(\psi\in U\).
At places in \(S\), \(\psi\) is locally zero, so the contributions vanish.
Outside \(S\), local Jacobian Kummer classes are unramified. The ramified
part of the global lift is the projected valuation parity of \((z)\),
namely the parity ideal \(\mathfrak I\).

Pairing this inertia class with the unramified class \(\psi\) gives its
Frobenius value at that prime. Summing these contributions yields
\[
\Delta_u(\beta,\psi)
 =\sum_{\mathfrak p\notin S_K}
      v_{\mathfrak p}(\mathfrak I)\chi_\psi(\mathfrak p)
 =\chi_\psi([\mathfrak I]).
\]
The norm projection introduces no extra term: \(\psi\) has norm-square
class, so its root-coordinate trace is zero. This proves the boxed
cup/CT identity at the beginning.

Both arguments must be strict for this simplified formula. The original
three partner masks \(1,6,128\) need not be strict. The full nine-bit target
therefore also requires the local correction terms at \(S\); this note
does not replace them with the simpler strict formula.

## Detected block dimensions

Restricting the retained CT forms to the already certified strict spaces
gives the following consequences for \(\kappa_u\):

| \(u\) | Retained strict dimension | Detected cup-image dimension, at least | Cup-kernel dimension in that space, at most | Full retained CT cross-rank |
|---:|---:|---:|---:|---:|
| \(-3\) | 2 | 0 | 2 | 2 |
| \(-2\) | 1 | 0 | 1 | 1 |
| \(-1\) | 5 | 4 | 1 | 5 |
| \(1\) | 3 | 2 | 1 | 3 |
| \(2\) | 2 | 2 | 0 | 2 |
| \(3\) | 4 | 2 | 2 | 4 |

These are bounds using the retained strict characters, not full
computations of \(\kappa_u\) or \(C\). A zero displayed rank does not
make the cup map zero. Uncomputed strict characters could detect more.

At \(u=-1\), in mask order
\(17108,34628,65575,404296,528076\), the strict CT matrix is
\[
\begin{pmatrix}
0&0&0&0&0\\
0&0&1&0&1\\
0&1&0&0&0\\
0&0&0&0&1\\
0&1&0&1&0
\end{pmatrix}.
\]
Its rank is four. Thus at least four independent directions already
fail the \(S\)-unramified cohomology lift, before imposing the full
Jacobian Kummer conditions at \(S\).

The class with mask 17108 is invisible to these five strict test
characters but has a nonzero CT cross-row against the larger retained
Selmer space. We cannot yet decide whether an uncomputed strict character
detects its cup class, or whether its cup class is zero and the remaining
obstruction lies in the local Kummer conditions at \(S\).

## Projection onto the known elementary factor

For the common \(u=-1\) class group, retain the
[three-dimensional direct factor](SAME_CLASS_FACTOR_DIFFERENT_SOLUBILITY.md)
with dual ideal words
\[
D_0=\mathfrak J_0\mathfrak J_2,\quad
D_1=\mathfrak J_4,\quad D_2=\mathfrak J_2.
\]
Evaluation by strict characters \(0,1,2\) defines its retraction.
The cup/CT identity gives the projections
\[
\begin{array}{c|ccccc}
\beta\text{ index}&0&1&2&3&4\\ \hline
\operatorname{pr}_H\kappa_{-1}(\beta)&
0&[\mathfrak J_2]&[\mathfrak J_4]&0&[\mathfrak J_4].
\end{array}
\]
They span dimension two. Since the detected total image has dimension
at least four, its intersection with the complementary class factor
has dimension at least two.

These projections are inferred from existing CT data. They are not
independent norm computations, and they do not identify complete ideal
representatives of the cup values.

## Mechanism and the remaining computation

The strongest refined candidate is now
\[
\text{specialization of }\gamma_u
\ \Longrightarrow\
\text{change in the cup map }\kappa_u
\ \Longrightarrow\
\text{change in simultaneous higher lifts}.
\]
Its kernel is a necessary solubility structure. It is still followed
by local Kummer compatibility and the distinction between a Selmer lift
and a rational point.

At the anchor \(u=0\), \(\gamma_0=1\), so the cup obstruction vanishes.
The separate split-Jacobian construction proves rational lifts for the
known anchor space there. Vanishing of the cup map alone does not
prove rationality of arbitrary unknown Selmer classes.

For Agent 1, the distinction is concrete:

- **Incidence:** cubic \(S\)-class characters supply candidate directions.
- **Solubility obstruction:** the cup map uses the actual labelled
  deformation, and can exclude several directions together even when
  every original local class is zero.
- **Visibility:** neither map measures chart exposure or point height.

This is still a relative fixed-cubic construction. A production-family
selector would need a compatible reference curve and class space
constructed independently of exceptional points. No such selector has
been supplied or installed.

The nine target CT bits remain open as an independent calculation.
The next arithmetic input is an exact norm witness in
\(K(\sqrt{1+\theta})/K\) for each of the three fixed strict classes,
followed by its parity ideal and the local Jacobian corrections.
No large class-group or norm-search campaign was launched to obtain it.

## Certificate and replay

The [protocol](STRICT_CUP_OBSTRUCTION_PROTOCOL.json) records the
retrospective inspection, finite limits and uncomputed endpoints.
The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_cup_obstruction_v1.json)
contains all 192 finite cup-action checks, the six \(S\)-unit norm-support
identities, strict CT restrictions, rank bounds and the projected ideal
words.

The new checks are exact finite arithmetic and binary linear algebra.
Strictness and support rely on the preceding independently replayed
local-square certificates; the CT values remain the retained exact
elliptic certificates.

```sh
python3 elliptic-curves/rank-jump/strict_cup_obstruction.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_strict_cup_obstruction.py
```

No new curve, point, CT value, norm witness or full class group was
computed. No active-search file or mathematical-status entry was changed.
