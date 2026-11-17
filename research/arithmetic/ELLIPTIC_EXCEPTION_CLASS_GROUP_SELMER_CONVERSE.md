# The level-one class-group condition is a strict Selmer condition

Under explicit good-prime and image hypotheses, all equivariant
class-group maps in EWC have a precise description by a strict
Selmer group. Its rational-point part measures the formal depth;
the remaining part comes from the Tate--Shafarevich group. No
finiteness or vanishing of that group is assumed.

## 1. Hypotheses and the global restriction map

Retain the fixed rank-two curve and rational group of
[EWC](ELLIPTIC_EXCEPTION_WIEFERICH_CLASS_GROUP_BRIDGE.md):

\[
 E:y^2=X^3-10758382X+6600759744,\qquad
 \Gamma=E(\mathbb Q),\qquad \Gamma_{\rm tors}=\langle T\rangle\simeq C_2.
\]

For each bad rational prime `ell`, retain
`m_ell=|E(Q_ell)/E_1(Q_ell)|`. In this note fix a prime `p` satisfying

\[
 \begin{gathered}
 p\ge7,\quad E\text{ has good reduction at }p,\quad
 V=E[p],\quad K=\mathbb Q(E[p]),\\
 G=\operatorname{Gal}(K/\mathbb Q)=\operatorname{GL}_2(\mathbb F_p),
 \qquad p\nmid m_\ell\quad\text{for every bad }\ell.              \tag{1.1}
 \end{gathered}
\]

The image condition is a stated hypothesis, not an inference from
full translation in a Kummer extension. The bad-local exclusions
form a fixed finite set. No computed free basis or displayed
saturation index is required.

Let

\[
 \mathcal C_p=\operatorname{Hom}_G(\operatorname{Cl}(K)/p,V),
 \qquad S_p=\operatorname{Sel}_p(E/\mathbb Q),
\]

where `Cl(K)/p` means `Cl(K)/p Cl(K)`. Write

\[
 S_p^0=\ker\!\left(S_p\longrightarrow H^1(\mathbb Q_p,V)\right).
                                                               \tag{1.2}
\]

This imposes zero localization at `p` and the usual Kummer local
condition at every other place. Under (1.1) it is the usual fine
Selmer condition as well, since the cohomology at each bad prime
will be zero. We use the explicit definition (1.2).

**Theorem ECS1 (class-group maps equal strict Selmer classes).**
Restriction to `K`, followed by the unramified Artin identification,
gives a canonical isomorphism

\[
                         \boxed{S_p^0\simeq\mathcal C_p.}        \tag{1.3}
\]

**Global cohomology.** The central subgroup `Z={I,-I}` has order
two and acts on `V` with zero invariants. Since two is invertible
in `F_p`, `H^j(Z,V)=0` for every `j>=0`. The Hochschild--Serre
spectral sequence therefore gives `H^i(G,V)=0` for every `i>=0`.
Inflation--restriction consequently gives the full isomorphism

\[
 H^1(\mathbb Q,V)\xrightarrow{\ \sim\ }
 H^1(K,V)^G
       =\operatorname{Hom}_G(G_K,V).                            \tag{1.4}
\]

The cohomology is continuous. The primary central-scalar statement
is Prasad--Shekhar, arXiv:1912.12928v1, Lemma 3.2, printed/PDF p.6;
(1.4) is also the isomorphism used in their section 4 diagram on
p.11. The spectral-sequence interface is Milne, *Arithmetic Duality
Theorems*, second edition, I Theorem 0.3 and Example 0.7, printed
pp.5,7/PDF pp.13,15:
[retained Prasad--Shekhar v1](sources/prasad-shekhar-2019-class-group-selmer-v1.pdf),
[version metadata](sources/prasad-shekhar-2019-class-group-selmer-v1.json),
[retained Milne](sources/milne-arithmetic-duality-2006.pdf),
[version metadata](sources/milne-arithmetic-duality-2006.json).

Under (1.4), `C_p` consists exactly of classes whose restrictions
vanish on every inertia group `I_w` over `K`. The point of the
following local argument is that this condition alone would not,
for a general curve and prime, be a classical or a strict Selmer
condition.

## 2. The local conditions in the converse

**At a bad prime `ell`.** The local argument in EWC1 shows that
multiplication by `p` is an automorphism of `E(Q_ell)` under
`p not dividing m_ell`. Hence `E(Q_ell)/p=0`. Local Tate duality,
with elliptic selfduality, identifies `H^1(Q_ell,E)[p]` with the
Pontryagin dual of `E(Q_ell)/p`, so it is zero too. The local Kummer
exact sequence gives

\[
                         H^1(\mathbb Q_\ell,V)=0.              \tag{2.1}
\]

The exact duality interface is Milne ADT, I Theorem 3.2 and
Corollary 3.4, printed pp.41,43/PDF pp.49,51. Thus no additional
bad-prime local classes survive in the present range.

**At a good prime `ell!=p`.** The torsion representation is
unramified, so `I_w=I_ell` for `w|ell`. The condition obtained from
`Cl(K)` is therefore ordinary unramified cohomology over `Q_ell`.
It equals the Kummer image `E(Q_ell)/p`: reduction identifies this
quotient with `E(F_ell)/p`, because the formal subgroup is
`p`-divisible, and finite-field Kummer theory identifies the latter
with `H^1(F_ell,V)`. Here `H^1(F_ell,E)=0` follows from
surjectivity of `Frob_ell-1` on `E(Fbar_ell)`: this nonconstant
separable endomorphism is surjective over the algebraic closure.
Finally `H^1(F_ell,V)=H^1_ur(Q_ell,V)`.

**At the residue prime `p`.** First,

\[
                             E(\mathbb Q_p)[p]=0.               \tag{2.2}
\]

Indeed rational two-torsion makes `N_p=#E(F_p)` even, and Hasse's
bound gives `N_p<2p` for `p>=7`; hence `p` does not divide `N_p`.
The formal group is torsionfree by the graded multiplication proof
in EWC1, and its finite torsion complement has order `N_p`.

Let `H=I_p/I_w` be the actual inertia image on `V`. The following
local calculation is the one explicitly written in the proof on
p.11 of the retained Prasad--Shekhar v1:

\[
                                H^1(H,V)=0.                    \tag{2.3}
\]

For completeness, in the ordinary case inertia acts in a basis by
matrices `((chi,*),(0,1))`, where `chi` is the mod-`p` cyclotomic
character, surjective onto `F_p^*` on inertia. If `H` has no
`p`-part, averaging proves (2.3). Otherwise its normal Sylow
`p`-subgroup is the full unipotent group `P` of order `p`. If
`A` is its invariant line and `B=V/A`, then

\[
 H^1(P,V)\simeq V/(g-1)V\simeq B,
 \qquad H^1(P,V)\simeq H^1(P,B).
\]

The last identification is induced by `V->B`. The action of `H/P`
on this cohomology is `chi^(-1)`, since it acts trivially on `B`
and by `chi` on `P`; it has no invariants. Averaging over `H/P`
therefore gives (2.3). In the supersingular case the good-reduction
inertia image is cyclic of order prime to `p`, so averaging again
applies. The ordinary representation and this supersingular inertia
input are the precisely scoped local imports from Prasad--Shekhar,
equation (3.1), p.6, and the p.11 calculation. We do not import
that paper's full converse theorem with its separate hypotheses.

If a local class becomes unramified over `K_w`, its restriction
to `I_p` lies in `H^1(H,V)` by inflation--restriction. Thus (2.3)
makes it unramified already over `Q_p`. But (2.2) gives

\[
 H^1_{\rm ur}(\mathbb Q_p,V)
    =V^{I_p}/(\operatorname{Frob}_p-1)V^{I_p}=0:                 \tag{2.4}
\]

on this finite vector space, the kernel of `Frob_p-1` is
`V^(G_Qp)=E(Q_p)[p]=0`, so its cokernel is zero as well.
Consequently the class-group local condition at `p` is exactly
zero localization. This rules out a residue-prime contribution
outside the strict Selmer group in the stated range.

At the infinite place, cohomology with these odd-order coefficients
is zero. Combining all places with (1.4) proves ECS1 in both
directions. The Artin map is canonical and equivariant as in EWC1.
\(\square\)

## 3. The rational contribution and the Tate--Shafarevich obstruction

Put

\[
 M_p=\Gamma/p\Gamma,\qquad L_p=E(\mathbb Q_p)/pE(\mathbb Q_p),
 \qquad X_p=\Sha(E/\mathbb Q)[p].
\]

Then `dim M_p=2`, `dim L_p=1`. Localization on `S_p` takes values
in the local Kummer line `L_p`; write it as `loc:S_p->L_p`.
The standard exact sequence is

\[
                       0\to M_p\to S_p\to X_p\to0.            \tag{3.1}
\]

It is equation (23) in the retained Milne *Elliptic Curves*,
IV section 2, printed p.110/PDF p.118; its Kummer construction is
equation (22), printed p.109/PDF p.117:
[retained book](sources/milne-elliptic-curves-2006.pdf).

For `xi in X_p`, choose a Selmer lift `c` in (3.1) and define

\[
 \delta_p(\xi)=\operatorname{loc}(c)
                   \bmod\operatorname{loc}(M_p)
       \quad\text{in }L_p/\operatorname{loc}(M_p).              \tag{3.2}
\]

Changing the lift changes its localization by `loc(M_p)`, so the
map is well-defined and linear.

**Theorem ECS2 (exact rational and Sha contributions).** There is
an exact sequence

\[
 0\longrightarrow M_p^0\longrightarrow\mathcal C_p
              \longrightarrow\ker\delta_p\longrightarrow0,
 \qquad M_p^0=\ker(M_p\to L_p).                                \tag{3.3}
\]

The rational term has dimension

\[
 \boxed{\dim M_p^0=
       \begin{cases}1,&e_p=1,\\2,&e_p>1.\end{cases}}           \tag{3.4}
\]

Here `e_p` is the depth of the actual full rational formal image,
as in EWC. In particular, if `e_p=1` then `delta_p=0` and

\[
                \boxed{\dim\mathcal C_p=1+\dim X_p.}           \tag{3.5}
\]

If `e_p>1`, then `delta_p:X_p->L_p` and

\[
                \boxed{\dim\mathcal C_p=2+\dim\ker\delta_p.}   \tag{3.6}
\]

**Proof.** The kernel-cokernel sequence applied to (3.1) and
localization gives

\[
 0\to M_p^0\to S_p^0\to X_p
    \xrightarrow{\delta_p}L_p/\operatorname{loc}(M_p)
    \to L_p/\operatorname{loc}(S_p)\to0.
\]

Use ECS1 and take the image in `X_p` to obtain (3.3). The formal
coordinate constructed in EWC4 identifies `E_1(Q_p)` with `Z_p`
and the full rational formal closure with `p^(e_p-1)Z_p`.
Thus the map `M_p->L_p` is onto if `e_p=1`, and zero if
`e_p>1`. Its kernel dimensions give (3.4); the target of
`delta_p` is respectively zero or `L_p`, proving (3.5)--(3.6).
This uses the full abstract group `Gamma`; no displayed-point
index is involved. \(\square\)

Thus the exact extra source is the subgroup of `Sha[p]` whose
classes admit a Selmer lift with zero localization at `p`.
At a nondeep prime every `Sha[p]` class admits such a lift, because
the rational points already fill the local Kummer line. At a deep
prime the possible obstruction to such a lift is the one-dimensional
map (3.2). Its value is not asserted to vanish.

## 4. An exact interpretation of the second module copy

**Corollary ECS3 (two-copy criterion).** Under (1.1), the following
are equivalent:

\[
 \boxed{\operatorname{Cl}(K)/p\twoheadrightarrow_G E[p]^2}
 \quad\Longleftrightarrow\quad
 \boxed{e_p>1\ \text{or}\ \Sha(E/\mathbb Q)[p]\ne0.}             \tag{4.1}
\]

If `Sha(E/Q)[p]=0`, more precisely,

\[
 \dim\mathcal C_p=\begin{cases}1,&e_p=1,\\2,&e_p>1.\end{cases}   \tag{4.2}
\]

The minimal value `dim C_p=1` occurs exactly when both
`e_p=1` and `Sha(E/Q)[p]=0`.

**Proof.** Since `G=GL_2(F_p)`, `V` is absolutely irreducible and
`End_G(V)=F_p`. Two independent homomorphisms from `Cl(K)/p`
to `V` give a surjection to `V^2`: a proper image inside the
semisimple module `V^2` would be one copy of `V`, making its
coordinate maps scalar multiples, contrary to independence.
Conversely a surjection gives two independent maps. Thus the
left side of (4.1) is equivalent to `dim C_p>=2`, as an actual
quotient, not merely an occurrence in a semisimplification.

If `e_p>1`, the rational term in (3.3) already has dimension two.
If `e_p=1`, equation (3.5) says that dimension is at least two
exactly when `X_p!=0`. This proves (4.1); taking `X_p=0`
proves (4.2). \(\square\)

At primes covered by (1.1), a second module copy at a nondeep
prime is therefore explained exactly by `Sha[p]`. The local
calculation of ECS1 leaves no other non-Selmer residue-prime source
in this range. Ordinary `p`-regularity is unavailable because the
one-copy baseline is forced, as already shown in EWC4.

This gives a pointwise interpretation of the second-copy density
hypothesis in EWC3. It does not show that the primes with
`Sha[p]!=0` have density zero, that `Sha` is finite, or that
`e_p>1` has density zero. If one separately assumes `Sha` finite,
its `p`-torsion is zero outside finitely many primes, and (4.1)
then agrees with the depth condition outside those primes among
the primes satisfying (1.1). No new Galois-image theorem is
asserted here.

The group in (1.2) is a finite-level strict Selmer group. No
control theorem in a cyclotomic or division-field tower, Iwasawa
invariant, generic `Sha` vanishing, BSD statement, weighted
recurrence or all-primes Hasse conclusion is inferred.

## 5. Scope of the primary inputs

The retained Prasad--Shekhar source is explicitly the 2019 arXiv
v1, not an unexamined later published version. We use its written
local representation and inertia-cohomology calculation, together
with the direct arguments above. We do not import an upper bound
by a guessed subtraction from the ordinary Selmer dimension.
The exact global kernel and its rational/Sha sequence are proved
here using the stated local hypotheses.

The standard class-field and Kummer/Selmer interfaces are the
retained Milne sources specified above. The level-one argument does
not require the uniform full-translation theorem of EWC, or a
computed Mordell--Weil basis. No numerical class-group computation,
artificial finite checker, formal verification or external review
is claimed.
