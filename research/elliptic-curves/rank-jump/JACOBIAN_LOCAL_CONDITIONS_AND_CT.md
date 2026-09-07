# Ordinary Jacobian local conditions and the CT comparison

The local condition left open in
[the torsion-difference note](TORSION_DIFFERENCE_AND_CT.md) is now
identified: it is the ordinary 2-Kummer image of the genus-two Jacobian.
Consequently the CT difference is the isogeny-descent pairing attached
to the factorization of multiplication by two on that Jacobian.

A bounded calculation at \(u=2\) certifies the connecting map on local
rational 2-torsion at all 17 retained places. It includes 2 and infinity,
with exact Hensel error bounds and 55 independent local-power comparisons.
This checks a necessary part of the local lifting data. It does not yet
supply complete middle Kummer bases or an independently computed CT entry.

## The compatible Kummer diagram

Let \(A=E_0\times E_u\), and keep the labelled dual isogenies
\[
\Phi:A\longrightarrow J_u,\qquad
\Psi:J_u\longrightarrow A,\qquad
\Phi\Psi=[2]_{J_u},\quad\Psi\Phi=[2]_A.
\]
Both kernels identify with \(V=E_0[2]\simeq E_u[2]\).
The exact torsion sequence is the sequence of kernels of this
composition:
\[
0\longrightarrow J_u[\Psi]
\longrightarrow J_u[2]\xrightarrow{\Psi}A[\Phi]\longrightarrow0.
\]
More simply, with both labelled kernels written as \(V\), it is
\(0\to V\to J_u[2]\to V\to0\). The left-hand copy embeds through
\(\Phi(V\times0)=\Phi(0\times V)\).

Fix a completion \(k=\mathbb Q_v\). Write \(L_0,L_u\) for the two
elliptic 2-Kummer images, \(C=L_0\cap L_u\), \(D=L_0+L_u\), and
\(\mathcal K_J\) for the usual image of \(J_u(k)/2J_u(k)\).
The isogeny Kummer maps satisfy
\[
\delta_\Psi(A(k))=D,\qquad \delta_\Phi(J_u(k))=C.
\]
For the first identity, take a half of \((P_0,P_u)\) in \(A(\bar k)\)
and apply \(\Phi\); its cocycle is the sum of the two elliptic Kummer
classes. For the second, \(\Psi(j)=(P_0,P_u)\) has equal elliptic
Kummer classes. Conversely, any pair with equal classes has zero
\(\Psi\)-Kummer obstruction, and therefore comes from \(j\in J_u(k)\).
The diagonal inclusion \(V\to V\oplus V\) is split, so it introduces
no ambiguity in the common class.

The Kummer diagram for \([2]=\Phi\Psi\) now proves
\[
i^{-1}(\mathcal K_J)=D,\qquad \Psi_*(\mathcal K_J)=C.
\]
It remains to identify this middle condition with the previously
constructed local Baer sum, rather than merely match its endpoints.

For that, let
\[
N=\{(R_0,R_u)\in E_0[4]\times E_u[4]:2R_0=2R_u\},
\]
with the stated torsion identification. The map \(\Phi:N\to J_u[2]\)
is the pushout map. A cocycle in \(N\) whose two components are
4-Kummer classes can be written as the cocycle of \(R\in A(\bar k)\)
with \(4R\in A(k)\). The \(N\)-condition says that \(2\Phi(R)\)
is rational. Its image is therefore the ordinary 2-Kummer class of
\(2\Phi(R)\).

Conversely, given \(j\in J_u(k)\), choose \(R\in A(\bar k)\) with
\(\Phi(R)=j/2\). Then \(4R=\Psi(j)\in A(k)\), and the cocycle of \(R\)
takes values in \(N\). Its image is the 2-Kummer class of \(j\).
Thus the image of these local pullback cocycles is exactly
\(\mathcal K_J\). The additional left-endpoint pushout contributes
\(i(D)\), already contained in \(\mathcal K_J\) by the preceding Kummer
diagram. This proves
\[
\boxed{\mathcal W_{\rm sum}=\mathcal K_J.}
\]
The argument works over every completion, including the real place.

## The ordinary isogeny-descent pairing

We therefore have an exact sequence of Selmer objects
\[
0\longrightarrow(V,D)\longrightarrow(J_u[2],\mathcal K_J)
\longrightarrow(V,C)\longrightarrow0.
\]
Its pairing on \(S_C=S_0\cap S_u\) is
\[
\Delta_u(x,y)=\operatorname{CT}_{E_0}(x,y)+
             \operatorname{CT}_{E_u}(x,y).
\]
The pairing identity and its kernel follow from
[Morgan–Smith, Theorem 1.3 and Proposition 4.4](https://arxiv.org/pdf/2103.08530).
The new step here is the compatible Kummer diagram identifying the
middle object with the ordinary Jacobian condition.

Because the fixed cubic is irreducible, \(V^{G_{\mathbb Q}}=0\).
The resulting global sequence is consequently
\[
\boxed{0\longrightarrow S_D\longrightarrow
\operatorname{Sel}_2(J_u/\mathbb Q)
\longrightarrow\operatorname{rad}(\Delta_u|_{S_C})\longrightarrow0.}
\]
In particular a class in the common elliptic Selmer space survives
the difference pairing exactly when it lifts through the Jacobian's
2-Selmer group. This is still a Selmer lift, not a rational point.

## The local connecting map must not be omitted

Locally \(V^{G_k}\) can be nonzero. Let
\[
t_v=\operatorname{rank}\bigl(V^{G_k}\longrightarrow H^1(k,V)\bigr)
\]
for the connecting map of the torsion sequence. Exactness gives
\[
\dim\mathcal K_J=\dim D_v-t_v+\dim C_v.
\]
Simply adding the endpoint dimensions overcounts when \(t_v>0\).

For a rational root \(a\in k\) of \(f\), the corresponding nonzero
point \(T_a\in E_0[2](k)\) maps to the sum of its two elliptic Kummer
classes. At the three cubic embeddings indexed by \(a,\theta_j,\theta_k\),
the ratio of their torsion Kummer representatives is
\[
(\gamma_j\gamma_k,\ \gamma_k,\ \gamma_j),\qquad
\gamma_i=1-u\theta_i.
\]
The distinguished component uses the derivative value at the torsion
point, rather than the zero value of \(x-\theta\).

There is a convenient expression in the local cubic algebra:
\[
e_a=\frac{\theta^2+a\theta+A+a^2}{3a^2+A},\qquad
b_a=1+u(a+\theta)+e_a\bigl(-ua+u^2(A+a^2)\bigr).
\]
Here \(e_a\) is the idempotent of the distinguished factor.
The local squareclass of \(b_a\) is exactly the connecting class.
The new script computes these classes for every rational local root.
Where all three roots are rational, their three classes sum to zero,
as required by the 2-torsion relation.

## A frozen \(u=2\) experiment with certified root errors

The [protocol](JACOBIAN_LOCAL_PROTOCOL.json) uses the existing complete
support and local coordinate maps, with one worker capped at 60 seconds.
It performs no point search or class-group computation. Each place has
an immutable checkpoint.

If \(a_0\) is a rational root approximation, put
\[
h=v_p(f(a_0)),\quad d=v_p(f'(a_0)),\quad s=h-d.
\]
The checker requires \(h>2d\). Hensel's lemma then supplies a unique
root with \(v_p(a-a_0)\ge s\). Since \(A,B,u,\theta,a_0\) are integral
at the places used, the formula above gives
\(v_{\mathfrak p}(b_a-b_{a_0})\ge e_{\mathfrak p}(s-2d)\).
For every prime of the cubic algebra it verifies
\[
e_{\mathfrak p}(s-2d)-v_{\mathfrak p}(b_{a_0})
>2v_{\mathfrak p}(2).
\]
Thus \(b_a/b_{a_0}\) is a square by Hensel's lemma. This is an exact
stability certificate, not agreement between two numerical precisions.
The root count agrees with the known local cubic decomposition, and
the Hensel balls for distinct roots are certified disjoint.

All finite-place stability margins are positive; the smallest is 15.
The real calculation uses exact algebraic signs. Every connecting class
lies in the already certified \(D_v\). An independent PARI local-power
test checks all 55 subset products at finite places and confirms exactly
the same kernels as the squareclass-coordinate calculation.

| Place | \(\dim D_v\) | \(\dim C_v\) | \(t_v\) | \(\dim\mathcal K_J\) |
|---|---:|---:|---:|---:|
| 2 | 2 | 2 | 1 | 3 |
| 3 | 1 | 1 | 0 | 2 |
| 5 | 1 | 1 | 1 | 1 |
| 7 | 2 | 0 | 0 | 2 |
| 13 | 3 | 1 | 2 | 2 |
| 17 | 0 | 0 | 0 | 0 |
| 31 | 3 | 1 | 2 | 2 |
| 79 | 1 | 1 | 0 | 2 |
| 233 | 3 | 1 | 1 | 3 |
| 647 | 3 | 1 | 2 | 2 |
| 1049 | 1 | 1 | 1 | 1 |
| 71889448247 | 1 | 1 | 1 | 1 |
| 40200713707633 | 1 | 1 | 1 | 1 |
| 26558335042564254253 | 2 | 0 | 0 | 2 |
| 47769497729370851316266401 | 2 | 0 | 0 | 2 |
| 491007790268548705232623905732119 | 1 | 1 | 1 | 1 |
| infinity | 2 | 0 | 2 | 0 |

These are dimensions and connecting classes, not a complete explicit
basis of \(\mathcal K_J\subset H^1(k,J[2])\).
For example the endpoint dimensions at 2 sum to four, but the actual
Jacobian Kummer dimension is three. At infinity it is zero despite a
two-dimensional left endpoint.

## Consequences for the common unknown Selmer quotient

The [relative full-Selmer theorem](RELATIVE_FULL_SELMER_THEOREM.md)
already proves \(\dim S_D=s_0+1\), where \(s_0=\dim S_0\).
Write \(\epsilon=s_0-20\ge0\), still UNKNOWN. Let \(k_u\) be the
radical dimension of the retained CT form on \(W_u\): it is two at
\(u=-1\), and one at the other five controls.

Since \(S_C/W_u\) has dimension \(\epsilon\), extension of an alternating
form gives
\[
\max(0,k_u-\epsilon)\le
\dim\operatorname{rad}(\Delta_u|_{S_C})\le k_u+\epsilon.
\]
The exact sequence above therefore proves
\[
21+\epsilon+\max(0,k_u-\epsilon)
\le\dim\operatorname{Sel}_2(J_u)
\le21+k_u+2\epsilon.
\]
In particular, if the single anchor completeness question were resolved
as \(\epsilon=0\), all six Jacobian Selmer dimensions would follow:
23 at \(u=-1\), and 22 at each other control. These numbers are
**conditional on \(\epsilon=0\)**; no such completeness certificate is
currently available. They are not Mordell–Weil rank claims.

This ties the earlier large Sha blocks to a standard isogeny-descent
filtration on one Jacobian. Killing an elliptic Sha class in \(H^1(J)\),
lifting a common class through \(\operatorname{Sel}_2(J)\), and
producing a rational Jacobian point are distinct statements.

## Next target and use for rank-jump analysis

The local-condition identification is closed. The next missing computation
is to build explicit middle local Kummer representatives, then calculate
a small retained CT submatrix through this isogeny sequence. The current
17-place calculation supplies the dimensions and torsion corrections
against which such representatives must be checked.

This is a **solubility** mechanism: the obstruction is failure to lift
through a specific isogeny-descent sequence. Its vanishing still does
not prove a rational point. For the recent MW17/MW16 fibres, a useful
prospective feature would additionally require a point-independent
candidate class space and compatible auxiliary construction. Neither
has been supplied by this fixed-field control experiment.

## Replay

The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_u2_isogeny_local_conditions_v1.json)
contains all rational root approximations, cubic representatives,
valuation/error bounds, squareclass signatures and source hashes.

```sh
python3 elliptic-curves/rank-jump/jacobian_local_conditions.py check
sage -python elliptic-curves/rank-jump/jacobian_local_conditions.py verify
sage -python elliptic-curves/rank-jump/verify_jacobian_local_power.py
```

No new CT value, rational point, parameter or active-search output was
computed or changed.
