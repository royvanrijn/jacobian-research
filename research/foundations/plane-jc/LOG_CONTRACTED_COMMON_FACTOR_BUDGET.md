# Contracted common factors, boundary minimality, and the cubic packet

> **Status.** Exact global `K`-theory identity and boundary-blowup lemma.
> If the logarithmic differential has a common rank-two factor along an
> effective contracted cycle `H`, its positive-dimensional quotient is
> `Omega_X^1(log D_X)/Omega_X^1(log D_X)(-H)` and has
> `ch_2=L_X.H-H^2`.  On a relative SNC-minimal compactification of `A^2`, a
> contracted boundary component has self-intersection at most `-2`; hence a
> common factor `hT` has strictly positive charge
> `h(v-2)+h^2(-T^2)`.
>
> For the isolated cubic E8 packet of the contracted-divisor Smith theorem,
> this closes the `I=3` attachment gate: its required remaining point class
> is negative on every relative minimal model.  The connected-cycle theorem
> below then saturates the **entire** contracted fiber.  Anti-nefness absorbs
> the fiber geometry; the generated-ratio Wronskian shows that any residual
> vertical determinant lies only on `H.T=0` components and has trivial kernel
> degree, while the parameter-ideal intersection number gives `-H^2>=2`.
> Consequently the
> complete saturated budget is at most three, while the forced local quotient
> has length four.  This excludes the one-component degree-six cubic E8
> equality row even when the contracted fiber is a larger connected cycle.
> A further horizontal determinant incidence is necessarily another affine
> ramification row and lies outside that equality row.  The theorem does not
> exclude higher geometric degree, multiple affine components, or the
> `k=2,...,24` target charts.

The formulas and inequalities are replayed by
[`verify_log_contracted_common_factor_budget.py`](../scripts/verify_log_contracted_common_factor_budget.py).

## 1. The common-factor exact sequence

Let

\[
 \theta:E\longrightarrow F,
 \qquad F=\Omega_X^1(\log D_X)
\]

be the logarithmic differential between rank-two bundles.  Suppose every
entry of `theta` is divisible by an effective Cartier divisor `H`.  Then
`theta` factors through `F(-H)` and gives the exact sequence

\[
 0\longrightarrow\operatorname{coker}(E\to F(-H))
 \longrightarrow\operatorname{coker}\theta
 \longrightarrow F/F(-H)\longrightarrow0.      \tag{1.1}
\]

The final term is the full positive-dimensional rank-two packet.  From

\[
 \operatorname{ch}(F/F(-H))
 =\operatorname{ch}(F)(1-e^{-H})
\]

and `rank(F)=2`, `c_1(F)=L_X`, one obtains

\[
 \boxed{\operatorname{ch}_2(F/F(-H))=L_X.H-H^2.} \tag{1.2}
\]

For `H=hT`, where `T` is a rational boundary component of valency `v` and
`T^2=-n`, logarithmic adjunction gives `L_X.T=v-2`, hence

\[
 \boxed{Q_{hT}=h(v-2)+h^2n.}                   \tag{1.3}
\]

This is the correct global replacement for the local notation
`(R/(t^h))^2`.  It is a curve class, not a point length.

## 2. Why relative minimality forces `n>=2`

Start from `(P^2,L_infinity)` and obtain a smooth SNC compactification of
`A^2` by blowing up boundary points.  Every boundary component `C` obeys

\[
 \boxed{C^2+v(C)\le1.}                         \tag{2.1}
\]

Indeed, the original line starts with `(C^2,v)=(1,0)`.  A blowup at a
smooth point of `C` changes this pair by `(-1,+1)`, while a blowup at a
boundary node on `C` changes it by `(-1,0)`.  A new exceptional component
starts as `(-1,1)` or `(-1,2)`, respectively.  These are all transitions,
so (2.1) follows inductively.

In particular, every boundary `(-1)` curve has valency at most two.  Its
blowdown preserves the SNC property and does not change the affine open.
If a completed morphism contracts that curve to a point, the morphism
factors through the blowdown.  Repeating this operation gives a relative
SNC-minimal model in which every contracted boundary component satisfies

\[
 \boxed{T^2\le-2.}                              \tag{2.2}

\]

Equations (1.3) and (2.2) imply

\[
 Q_{hT}\ge-h+2h^2>0.                           \tag{2.3}

Thus the common split quotient itself never supplies a negative correction
on the relative minimal model.  This statement is deliberately narrower
than positivity of the whole logarithmic cokernel: after saturation, an
unequal generic Smith row can leave a cyclic contracted determinant whose
kernel-line and self-intersection terms still need separate control.

## 3. Saturating the complete Chern budget

Write the complete logarithmic determinant as

\[
 D_{\log}=D'+2H.                                \tag{3.1}

\]

Suppose `U` is the point budget obtained before removing the common factor,
after the affine conormal and the complete half-square
`D_log^2/2` have been subtracted.  Replacing that square by `D'^2/2` and
subtracting (1.2) gives the exact saturated budget

\[
 \boxed{
 U_{\rm sat}
 =U+2D'.H+3H^2-L_X.H.}                         \tag{3.2}

\]

For `H=hT`, `I=D'.T`, this is

\[
 \boxed{
 U_{\rm sat}
 =U+2hI-3h^2n-h(v-2).}                         \tag{3.3}

\]

This identity is additive in `K_0`; it makes no generic-cyclicity
assumption about the unsaturated matrix.

## 4. The isolated cubic E8 gate closes

The minimal cubic E8 normal form has:

- common order `h=1` on the contracted component `T`;
- unsaturated complete-chain budget `U=2`;
- saturated determinant incidence `I=3` when the cusp row is the only
  residual determinant meeting `T`; and
- a forced saturated cyclic-submodule quotient of length four.

Equation (3.3) gives

\[
 U_{\rm sat}=4-v+2I-3n.
\]

After the forced length-four quotient, every other isolated point term has
total

\[
 \boxed{P_{\rm other}=2I-v-3n.}                \tag{4.1}

\]

For the isolated incidence `I=3`, relative minimality `n>=2` yields

\[
 P_{\rm other}=6-v-3n<0,
\]

which contradicts isolated-`Fitt_1` positivity.  Therefore

\[
 \boxed{\text{the one-component, one-contracted-packet cubic `I=3` row is impossible.}}
\]

This removes the three apparent `(-1)` survivors in the earlier attachment
table: boundary blowup geometry allows a `(-1)` component only at valency
one or two, and those curves are not relative minimal.

## 5. The complete connected contracted cycle

Let `H=sum_i h_iT_i` be the full divisorial common factor over one affine
target point `p`.  Resolve the boundary base ideal generated by local target
parameters `(x,y)` so that

\[
 (x,y)O_X=O_X(-H).                              \tag{5.1}
\]

The two generators make `O_X(-H)` relatively globally generated.  Therefore

\[
 \boxed{H\mathbin\cdot T_i\le0\quad\hbox{for every }T_i\subset H.} \tag{5.2}
\]

Thus `H` is anti-nef.  Contract every removable boundary `(-1)` curve on
which the map is constant.  On the resulting relative SNC-minimal model
every `T_i^2=-n_i` has `n_i>=2`.

Assume the connected support of `H` meets the rest of the boundary in one
edge, at `T_0`, and `h_0=1`.  This is the incidence of the minimal cubic
normal form.  Let `d_i` be the valency of `T_i` **inside** the support of
`H`, and put

\[
 a_i=-H\mathbin\cdot T_i
     =n_ih_i-\sum_{T_j\sim T_i}h_j\ge0,
 \qquad q=-H^2=\sum_i h_i a_i.                 \tag{5.3}
\]

The full boundary valency is `d_i` except at `T_0`, where it is `d_0+1`.
Logarithmic adjunction consequently gives

\[
 L_X\mathbin\cdot H
 =\sum_i h_i(d_i-2)+1.                          \tag{5.4}
\]

Summing (5.3) without the factors `h_i` yields the exact tree identity

\[
 \sum_i a_i+L_X\mathbin\cdot H
 =1+\sum_i(n_i-2)h_i\ge1.                      \tag{5.5}
\]

Since `h_i>=1` and `a_i>=0`, one has `q>=sum_i a_i`; hence

\[
 \boxed{L_X\mathbin\cdot H\ge1-q.}             \tag{5.6}
\]

There is also an intrinsic interpretation of `q`.  Take the Stein
factorization `X->Z->Y` near the connected fiber.  Equation (5.1) and the
local intersection formula for the parameter ideal give

\[
 \boxed{q=-H^2=\ell\bigl(O_{Z,z}/(x,y)\bigr).}  \tag{5.7}
\]

Indeed, two general generators have divisors `H+A_x` and `H+A_y`; relative
principal-divisor intersection gives `A_x.H=A_y.H=-H^2`, while basepoint
freeness makes `A_x.A_y=0` over `z`, leaving local intersection `-H^2`.
The normal surface local ring `O_(Z,z)` is Cohen--Macaulay, so `(x,y)` is a
parameter ideal.  A nonempty fiber on the relative minimal model lies over a
singular point of `Z`.  If the length in (5.7) were one, `(x,y)` would be the
maximal ideal, making this two-dimensional local ring regular.  Therefore

\[
 \boxed{q\ge2.}                                 \tag{5.8}
\]

This is the elementary parameter-ideal form of the maximal-ideal-cycle
multiplicity bound.  It agrees with Wagreich's intersection theorem for
maximal ideal cycles.

There is one more consequence of (5.1) that controls the possible saturated
vertical determinant.  After dividing `(x,y)` by the local equation of `H`,
the two leading coefficients generate

\[
 O_X(-H)|_{T_i},\qquad
 \deg O_X(-H)|_{T_i}=-H.T_i=a_i.                \tag{5.9}
\]

Their ratio defines a morphism `phi_i:T_i->P^1` of degree `a_i`.  In
characteristic zero, if `a_i>0`, its Wronskian is generically nonzero.  The
logarithmic differential after removal of the common factor then has full
generic rank on `T_i`, so `T_i` is not a component of the residual
determinant `D'`.  If `a_i=0`, the generated line bundle is trivial and the
ratio is constant.  A fixed target covector annihilates that leading
direction, so the kernel line of any residual cyclic determinant packet on
`T_i` has degree zero.  Consequently, if `R` is the residual determinant
cycle supported on `H`,

\[
 \boxed{R.H=0,\qquad \deg K_R=0.}               \tag{5.10}
\]

On a connected zero-degree subgraph the constant directions agree at the
nodes because the ratio is the restriction of the single morphism to the
exceptional line of the blowup of `p`.  Any remaining failure of cyclicity is
isolated, and the cyclic-submodule theorem makes its quotient effective.
Thus a residual vertical determinant cannot hide a negative kernel-degree
correction in the argument below.

## 6. Connected-cycle cubic exclusion

Return to the degree-six cubic E8 equality row.  Its unsaturated
complete-chain budget is `U=2`.  Write

\[
 D_{\log}=D'+2H.                                \tag{6.1}
\]

The unique ramified horizontal row is the cubic divisor `E` with coefficient
three, and `E.H=h_0=1`.  Every other residual determinant component meeting
this connected fiber is supported on `H`: another horizontal component would
be another ramified affine row and would change the equality passport.
Write `D'=3E+R+D_{\rm disj}`, where `R` is effective and supported on `H`,
and `D_disj.H=0`.  Equation (5.10) gives

\[
 \boxed{D'.H=3+R.H=3,}                         \tag{6.2}
\]

and also shows that `R` contributes no negative kernel-line degree.

Apply the general saturation identity (3.2) to the whole cycle.  Using
(5.6)--(5.8),

\[
\begin{aligned}
 U_{\rm sat}
 &=2+2D'.H+3H^2-L_X.H\\
 &\le8-3q-(1-q)\\
 &=7-2q\le3.                                   \tag{6.3}
\end{aligned}
\]

At `E cap T_0`, the saturated cubic Smith normal form remains
`[[r,0],[t^2,r^2]]`.  Its cyclic-submodule quotient has forced length four.
Equation (6.3) gives fewer than four total units for that quotient before any
other effective isolated correction is added.  The residual vertical kernel
degrees vanish by (5.10), so there is no signed divisorial term left that can
offset the deficit.  Hence

\[
 \boxed{
 \text{the one-component degree-six cubic E8 equality row is impossible.}}
                                                               \tag{6.4}
\]

This closes both former internal escapes: enlarging the connected contracted
fiber and adding residual **vertical** determinant multiplicity cannot help.
An additional horizontal determinant incidence is a different affine
ramification row, not a completion of the equality row.  The affine strict
log-etale theorem removes the four target-node resolution packets, and the
terminal, carrier, root, and outgoing-tail packets are already compiled, so
there is no unclassified packet inside the stated one-component row that can
restore the missing unit.

## 7. Remaining scope

The connected-cycle theorem is local-to-global but its F2 application is
specific to the degree-six equality passport.  Full `(75,125)` closure still
requires exclusion of:

1. higher geometric degrees through `28`;
2. two or more affine nonproperness components (including a second row
   through the cubic target value);
3. non-E8 `k=1` strata outside the already completed degree-six atlas; and
4. the normalization-contact charts `k=2,...,24`.

Those are now the genuine escapes; a larger connected contracted cusp fiber
is no longer one of them.

For the classical multiplicity comparison, see P. Wagreich,
[*Elliptic singularities of surfaces*](https://doi.org/10.2307/2373333),
Theorem 2.7.  The argument above only uses its elementary two-generated
parameter-ideal special case.
