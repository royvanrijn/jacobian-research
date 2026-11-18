# Prime augmentation gaps at unbounded bridge number

Date: 18 November 2026.

For `b>=5`, let `d=b-2` and use the braids from [UNIFORM_PROTECTED_AUGMENTATION_GAPS.md](UNIFORM_PROTECTED_AUGMENTATION_GAPS.md):

\[
\delta_b=\sigma_2\sigma_3\sigma_4\sigma_3^{-1}
                  \sigma_5\cdots\sigma_{b-1},\qquad
K_{b,k}=\widehat{U\delta_b^{(d+1)(1+d^2k)}}\quad(k\ge0).
\]

**Theorem.** There is an integer `B` such that, for every `b>=B`, infinitely many distinct prime knots among the displayed `K_(b,k)` satisfy

\[
\operatorname{ar}_{\mathbb C}(K_{b,k})\le b-1<
\operatorname{mr}(K_{b,k})=\operatorname{br}(K_{b,k})=b.
\]

More precisely, every `b>=5,k>=0` has the displayed augmentation and bridge bounds, and for every sufficiently large `b` there is a threshold `k_0(b)` beyond which all these knots are hyperbolic. The bounds `B,k_0(b)` are not made effective here. The previously certified bridge-five and bridge-six prime tails remain valid independently of those thresholds.

The new step is a uniform exterior homeomorphism to fillings of one fixed hyperbolic three-cusped manifold. It is proved by two finite Reidemeister tangle cells with fixed boundary ports and an explicit induction. Bounded checks supplement that proof; they do not supply its quantifiers.

## 1. The uniform transfer to one fixed manifold

Set

\[
\rho_b=\sigma_{b-1}\cdots\sigma_4,\quad
\theta_d=\sigma_2\cdots\sigma_{b-2},\quad
Q'_b=\rho_bU\rho_b^{-1}\theta_d,
\]

and

\[
P_b=\sigma_b\sigma_{b-1}\cdots\sigma_3\sigma_2^2
                    \sigma_3\cdots\sigma_{b-1}\sigma_b^{-1}.
\]

Let `L'_b=closure(Q'_b P_b)`. The last strand is an unknot `C` bounding the disk through positions two through `b-1`. Let `M'_b` be this two-component link exterior. The earlier parent `M_b`, with base exponent `d+1` instead of one, differs by one full disk twist around `C`, so `M_b` and `M'_b` are homeomorphic as unmarked exteriors. This change can alter the circle peripheral marking; the knot meridian is preserved.

Use the fixed six-braid

\[
R=\sigma_4\sigma_3\sigma_2^2\sigma_3\sigma_4\sigma_5^{-1}
               \sigma_4U\sigma_4^{-1}\sigma_2\sigma_3\sigma_5,
\]

and add a seventh strand with word

\[
D=\sigma_6\sigma_5\sigma_4^2\sigma_5\sigma_6^{-1}.
\]

Write `N` for the exterior of `closure(RD)`. The new component `D_0` bounds the disk meeting one knot arc and one old-circle arc at positions four and five of `R`. Its preferred disk slope `(1,-q)` inserts `q` positive full twists between those arcs, giving

\[
N\bigl(D_0:(1,-q)\bigr)
\cong S^3\setminus\widehat{R\sigma_4^{2q}}.
\tag{1}
\]

The fixed link, component marking, disk-surgery sign and rigorous hyperbolicity of `N` are retained in [the universal-parent geometry note](PROTECTED_UNIVERSAL_PARENT_CANDIDATE.md) and its exact inputs and interval receipt. The following diagram lemma supplies the uniform identification.

**Uniform diagram lemma.** For every `b>=5`, there is an oriented link isotopy

\[
L'_b\cong\widehat{R\sigma_4^{2(b-5)}}
\tag{2}
\]

which preserves the knot and old-circle component identities. Consequently

\[
M_b\cong N\bigl(D_0:(1,-(b-5))\bigr)
\tag{3}
\]

as unmarked exteriors.

We prove the diagram lemma below with all outside ports fixed. Equation (3) does not specify the old-circle longitude in the original three-cusped peripheral basis. That marking is unnecessary for parent hyperbolicity. The second filling argument uses the geometrically marked circle of the original `M_b` separately.

## 2. Port convention and the source ladder

A crossing has cyclic ports `0,1,2,3`; the even branch is under and the odd branch over. At every positive crossing, ports `0,3` are incoming and `1,2` outgoing. A positive braid letter uses incoming ports `(0,3)` and outgoing ports `(1,2)` in left-to-right strand order. A negative letter uses incoming `(1,0)` and outgoing `(2,3)`. This fixes the literal spherical braid projection, including the closing arcs.

Let `r=b-4`. Number the crossings in the written source word `Q'_b P_b` from zero. First cancel the outer-circle detour crossings numbered `3b+3` and `5b` by a Reidemeister II move. Then for `i=0,...,r-1` cancel the pair

\[
i,\quad 2r+13-i.
\]

Here is the universal port contract for those cancellations. In every pair `(X,Z)`, including the outer-circle pair, the bigon edges are `X.2-Z.0` and `X.3-Z.3`. Removing the outer pair joins `rho_0.3` to `(2b+5).3`, creating the first rho bigon, and joins `u_0.3` to `v_0.2`. Before rho step `i`, write `X=i`, `Z=2b+5-i`. Their four outside neighbors are

\[
X.0-v_i.1,\quad
X.1-\begin{cases}(i+1).3&i<r-1,\\F_2.0&i=r-1,\end{cases}
\quad Z.1-(Z-1).3,\quad Z.2-Y_{r-i}.3.
\]

RII therefore joins `v_i.1` to `Y_(r-i).3` and creates the next `X.3-Z.3` bigon. At the last step the second join is `F_2.0-F_13.3`. This proves availability of the whole indexed cancellation by induction from the literal braid wiring. Opposite signs and the displayed cyclic ports give genuine cancelling RII disks. None meets a crossing of `U`.

The resulting graph is the following closed-form ladder. This reduction is uniform, rather than a choice of moves found separately for each degree.

Name the fourteen unchanged crossings `F_0,...,F_13`. Name the `d-1` Coxeter crossings `Y_0,...,Y_(d-2)`. In the circle suffix name the descending crossings `u_0,...,u_(d-1)` and the ascending crossings in reverse order `v_0,...,v_(d-1)`. Their original source indices are

\[
F_j=r+j,\quad Y_i=2b+6+i,\quad
u_j=3b+4+j,\quad v_j=5b-1-j,
\]

with the same crossing names as in the table below.

For readability write `X.i` for port `i` of crossing `X`, and put

\[
A=F_{12}.1,\quad B=F_{13}.2,\quad
C=F_0.3,\quad E=F_1.3.
\]

These are the four attachment ports to the fixed `U` block. The following rows list the neighbors of ports `0,1,2,3`, respectively. Every edge appears reciprocally; the table determines the entire outside graph.

| Crossing | Port 0 | Port 1 | Port 2 | Port 3 |
|---|---|---|---|---|
| `Y_0` | `A` | `u_(d-1).0` | `Y_1.0` | `B` |
| `Y_i`, `1<=i<=d-2` | `Y_(i-1).2` | `u_(d-1-i).0` | `Y_(i+1).0`, or `u_0.0` if `i=d-2` | `v_(d-2-i).1` |
| `u_0` | `Y_(d-2).2` | `u_1.3` | `v_0.3` | `v_0.2` |
| `u_j`, `1<=j<=d-2` | `Y_(d-1-j).1` | `u_(j+1).3` | `v_j.3` | `u_(j-1).1` |
| `u_(d-1)` | `Y_0.1` | `v_(d-1).0` | `v_(d-1).3` | `u_(d-2).1` |
| `v_0` | `v_1.2` | `Y_(d-2).3` | `u_0.3` | `u_0.2` |
| `v_j`, `1<=j<=d-3` | `v_(j+1).2` | `Y_(d-2-j).3` | `v_(j-1).0` | `u_j.2` |
| `v_(d-2)` | `v_(d-1).2` | `E` | `v_(d-3).0` | `u_(d-2).2` |
| `v_(d-1)` | `u_(d-1).1` | `C` | `v_(d-2).0` | `u_(d-1).2` |

The table follows directly by tracing the descending and ascending suffix and the consecutive Coxeter crossings; its two end cases correspond to the doubled `sigma_2` and the two fixed `U` attachments. All its crossings are positive. The rest of the graph is the same fixed `F` block for every `d`.

## 3. Two finite embedded-disk cells

The contract of each cell includes the cyclic order of the boundary endpoints, their incoming/outgoing directions, the crossing cyclic orders and over/under branches. Thus it is a disk-tangle rewrite, not a Gauss-word equivalence.

For a triple `(X,a),(Y,b),(Z,c)`, the face edges are

\[
X.(a+1)-Y.b,\quad Y.(b+1)-Z.c,\quad Z.(c+1)-X.a.
\]

The indicated RIII moves below have three distinct crossings and the required over/under pattern. These edges form a complete triangular face in the spherical rotation system. Since the projection is connected and planar, the face is an empty disk. The moves preserve that embedding. The RI monogon below is likewise a one-corner face, so is empty. These conditions exclude a strand or a disconnected component hidden in a purported move region.

**Clasp-slide cell.** Use crossings `u,t,v,w` with triangle edges

\[
u.0-t.2,\quad t.3-v.1,\quad v.2-u.3,
\]

and an existing clasp

\[
u.1-w.0,\qquad u.2-w.3.
\]

The six external ports attach as

\[
t.0-e_0,\ t.1-e_1,\ v.0-e_2,\ v.3-e_3,
\ w.1-e_4,\ w.2-e_5.
\]

Their cyclic boundary order is `(e_0,e_1,e_4,e_5,e_3,e_2)`; the inputs are `e_0,e_2,e_3`. Apply the two RIII faces

\[
(u,3),(t,2),(v,1);\qquad
(w,3),(v,2),(t,1).
\]

The old clasp edges `u.1-w.0,u.2-w.3` remain, while the other internal edges become

\[
w.1-t.0,\quad w.2-v.0,\quad t.3-v.1.
\]

The new boundary attachments are

\[
u.0-e_0,\ t.1-e_1,\ u.3-e_2,\ v.3-e_3,
\ t.2-e_4,\ v.2-e_5.
\]

Thus a moving triangle has passed one clasp, with every outside endpoint fixed and its orientation preserved.

**Cap cell.** Use the same three triangle edges on `u,t,v`, add the cap edge `u.2-v.3`, and attach

\[
t.0-e_0,\quad t.1-e_1,\quad u.1-e_2,\quad v.0-e_3.
\]

The boundary order is `(e_0,e_1,e_2,e_3)` and the inputs are `e_0,e_3`. Apply RIII at `(u,3),(t,2),(v,1)`. It creates the monogon edge `t.2-t.3`. Remove that positive curl by RI. The remaining two crossings form a positive clasp with edges

\[
u.1-v.0,\quad u.2-v.3,
\]

and external attachments

\[
u.0-e_0,\quad v.1-e_1,\quad v.2-e_2,\quad u.3-e_3.
\]

The [finite input](threaded_curl_uniform_transfer_inputs.json) retains both complete before/after pairings and the moves. The standard-library checker verifies the disk boundary order, Euler characteristic, oriented strands, triangular faces, over/under admissibility and RI monogon. It rejects an altered RIII port and a switched external endpoint.

## 4. The induction and its endpoint

Let `q=d-3`. For stage `s=0,...,q-1`, put

\[
t=Y_{d-2-s}.
\]

For `j=0,...,s`, perform

\[
(u_j,3),(t,2),(v_s,1),
\]

and, if `j<s`, then perform

\[
(v_j,3),(v_s,2),(t,1).
\]

After `j=s`, remove `t` by the RI monogon at its ports two and three.

Here is the invariant that proves every step is available for arbitrary `d`. At the beginning of stage `s`, there are `s` already extracted clasp pairs `(u_j,v_j)`, `j<s`, with edges

\[
u_j.1-v_j.0,\quad u_j.2-v_j.3.
\]

The moving triangle enters at `u_0`; its edges are

\[
u_0.0-t.2,\quad t.3-v_s.1,\quad v_s.2-u_0.3.
\]

For each existing clasp, the following pair is connected by the two rails

\[
v_j.1-u_{j+1}.0,\qquad v_j.2-u_{j+1}.3.
\]

Thus the clasp-slide cell moves the triangle from `u_j` to `u_(j+1)` while preserving all previous clasp pairs. At `u_s`, the cap edge `u_s.2-v_s.3` is present. The cap cell deletes exactly `t` and creates the next clasp `(u_s,v_s)`.

The undeleted ladder supplies the remaining attachments: before these slides, `t.0` leads to the preceding `Y` crossing, `t.1` to `u_(s+1).0`, `u_s.1` to `u_(s+1).3`, and `v_s.0` to `v_(s+1).2`. Tracking those four ends through the two cells gives the same moving-triangle and rail conditions for stage `s+1`. The next cap edge was untouched. This proves the invariant by induction from the source table at `s=0`. All four attachment ports `A,B,C,E` and all ports of every `F_j` remain fixed.

At the end, `Y_2,...,Y_(d-2)` have disappeared. The remaining base ladder is the `d=3` ladder on

\[
Y_0,Y_1;\quad u_q,u_{q+1},u_{q+2};\quad
v_q,v_{q+1},v_{q+2},
\]

with `q` positive clasp pairs inserted on the two arcs formerly joining

\[
Y_1.2-u_q.0,\qquad v_q.2-u_q.3.
\]

More explicitly, these arcs now start at `Y_1.2-u_0.0` and `v_q.2-u_0.3`; each pair `(u_j,v_j)`, `j<q`, has the two clasp edges and then the rail connections to `u_(j+1)` displayed above. All other base-ladder edges are unchanged. This describes the endpoint graph for every `q`.

The raw target `R sigma_4^(2q)` has `26+2q` crossings. Its two fixed RII cancellations normalize the base `R` to 22 crossings, giving `22+2q`, as does the source count `14+(d-1)+2d-(d-3)`. The final graph is identified with this normalized target, not with its unreduced raw projection.

Here is the finite target-cut binding, with both selected arcs held open. Number the 26 crossings of `R` from zero; `U` occupies crossings eight through 21. Before reduction, the fourth closure arc is `(24,2)-(0,0)` and the fifth is `(25,1)-(0,3)`. Open them and attach the four distinct named boundary endpoints shown below. Cancel RII pairs `(6,25)` and `(7,22)`. Both remain valid with the two arcs open, and both miss crossing zero and the entire `U` block.

| Named boundary endpoint | Before the two RII moves | After the two RII moves |
|---|---|---|
| `Kin` | `(0,0)` | `(0,0)` |
| `Kout` | `(24,2)` | `(24,2)` |
| `Cin` | `(0,3)` | `(0,3)` |
| `Cout` | `(25,1)` | `(5,2)` |

The `d=3` source ladder maps `Y_1,u_0,v_0` to target crossings `24,0,5`, respectively, with no port rotation, and fixes every `U` port. Thus reclosing the named pairs identifies the two reduced target cut arcs exactly as `Y_1.2-u_0.0` and `v_0.2-u_0.3`. The open four-port contract retains the insertion band: it is the small standard braid slab immediately before the first positive `sigma_4` crossing zero. The two RII disks are disjoint from that slab. Consequently this identifies the actual fourth/fifth strand band, not merely two abstract arcs of the same components.

The cap cell's two-crossing output has internal edges `u.1-v.0,u.2-v.3`, incoming boundary ports `u.0,u.3` and outgoing boundary ports `v.1,v.2`. It is exactly the positive two-strand braid `sigma_4^2` at that marked band. Each repetition inserts the same two-rail cell, so the indexed endpoint description agrees with the normalized `R sigma_4^(2q)` for every `q`. Crossing cyclic orders, over/under branches, named cut endpoints and `U` attachments agree. This proves the unrestricted endpoint identification rather than only a bounded comparison of complete diagrams.

There are `1+3+...+(2q-1)=q^2` RIII moves and `q` RI moves, after the initial `r+1` RII cancellations. The fixed `U` block lies in the knot component, so preserving all its ports identifies that component; the other component is the old circle. All moves preserve the oriented link. This proves (2), and (1) gives (3). ∎

## 5. Two Dehn filling arguments and the quantifiers

The fixed `N` has a rigorous complete finite-volume hyperbolicity certificate. By Thurston's hyperbolic Dehn filling theorem, all but finitely many slopes `(1,-q)` on `D_0` give hyperbolic two-cusped manifolds while the two old cusps stay complete. Since `q=b-5` tends to infinity, equation (3) proves that `M_b` is hyperbolic for every sufficiently large `b`.

Fix any such `b`. The original geometrically marked parent, independently of the marking change in (3), has

\[
M_b\bigl(C:(1,-d(d+1)k)\bigr)
\cong S^3\setminus K_{b,k}.
\]

Those slopes tend to infinity as `k` grows, so a second application of hyperbolic Dehn filling makes every sufficiently large `K_(b,k)` hyperbolic and prime. The thresholds may depend on `b`. No uniform threshold in `k`, effective list of exceptional bridge numbers or assertion of hyperbolicity for every `b>=5` is needed.

For each fixed `b`, infinitely many distinct knot types occur. In an embedded circle horocusp, the Euclidean lengths of the displayed slopes tend to infinity. Futer–Kalfagianni–Purcell's Theorem 1.1 bounds the filled volumes below by a factor tending to one times `vol(M_b)`; Thurston's Theorem 6.5.6 makes every proper hyperbolic filled volume strictly smaller. The volumes therefore approach `vol(M_b)` from below. A finite collection of hyperbolic exterior types would have a finite maximum volume strictly below it, contradicting convergence and rigidity. Exact primary editions and statements are retained in the source ledger. This is the fixed-parent argument already given in the bridge-five proof, now applied separately for each sufficiently large `b`.

Finally, the uniform alternating quotient and protected paired-span theorem give `ar_C<=b-1<mr=br=b` for every parameter. Combining these statements proves the theorem.

The circle framing deserves separation. The diagram isotopy (2) preserves the link components and their meridians. The disk twist relating `M'_b` and `M_b`, and the filling of `D_0` which links `C` once, can change the expression of the preferred `C` longitude in a pre-existing basis. We use (3) only to transfer hyperbolicity. The second filling uses the original source `M_b` disk marking; no unproved simultaneous slope formula is invoked.

## 6. Exact evidence

The new [checker](verify_threaded_curl_uniform_transfer.py), [input](threaded_curl_uniform_transfer_inputs.json), [receipt](threaded_curl_uniform_transfer_verification.json) and [source ledger](threaded_curl_uniform_transfer_source_ledger.json) retain the finite port proof. The checker verifies both local disk contracts exactly. As supplementary controls, it constructs the literal source and target braid rotation systems for `5<=b<=64`, performs the indexed RII/RIII/RI sequence, and verifies a spherical diagram isomorphism fixing all 56 ports of the fourteen `U` crossings. It uses no floating point arithmetic, Sage topology recognition or sampled hyperbolic isometries.

```bash
.venv/bin/python research/meridional-rank/representations/verify_threaded_curl_uniform_transfer.py \
  --write-receipt research/meridional-rank/representations/threaded_curl_uniform_transfer_verification.json
```

The separate fixed-parent geometry input and Sage replay remain indispensable imported evidence. No formal verification, independent replay or external human review is claimed by this checker. An internal independent-agent audit is tracked separately if completed. No novelty or priority claim is made here.
