# Explicit even alternating Keller maps from cyclic cones

For every even `N>=6`, the formulas below give a polynomial map of affine three-space with constant Jacobian determinant `-1` and generic geometric monodromy `A_N`. For every odd `m>=3`, the member of degree `N=m^2+1` is defined over the **same** field `E=Q(sqrt(-3))` and has regular arithmetic monodromy `A_N` over `E`.

These are explicit two-term members of the cyclic-cone construction of Kistner–Shaska, [arXiv:2608.02863v1](https://arxiv.org/pdf/2608.02863v1), Section 3. That paper already establishes the cyclic-cone family and its monodromy. The purpose here is an elementary closed formula, its exact arithmetic constants, and its application to the alternating-Keller realization question. Generic alternating polynomial Keller realization in dimension three is therefore supplied by that existing literature. No novelty or rational descent claim is made.

## Theorem 1. A closed formula in every even degree

Let `d>=3`, `N=2d`, and work over

\[
 K_d=\mathbb Q(D),\qquad D^2=-3(2d-1).
\]

Set

\[
 a=\frac{-3-D}{2(d-2)},\qquad
 b=\frac{2d-1+D}{2(d-2)}.
\tag{1}
\]

Both coefficients are nonzero. In independent variables `Lambda,w`, define

\[
\begin{split}
A={}&w+a\Lambda w^4+b\Lambda^{d-1}w^{2d},\\
B={}&1+4a\Lambda w^3+\frac{2db}{d-1}\Lambda^{d-1}w^{2d-1}
 +\frac{8a^2}{3}\Lambda^2w^6\\
 &+\frac{4d^2ab}{(d-1)(d+1)}\Lambda^dw^{2d+2}
 +\frac{2d^2b^2}{(d-1)(2d-1)}\Lambda^{2d-2}w^{4d-2}.
\end{split}
\tag{2}
\]

Put `Lambda=1+xy`, `w=xz-1`, and

\[
 F_d(x,y,z)=\left(\frac{A(1+xy,xz-1)}x,
                       \frac{B(1+xy,xz-1)}x,
                       x(1+xy)\right).
\tag{3}
\]

Then (3) is polynomial, with component degrees `(6d-3,12d-9,3)`, determinant `-1`, and generic degree `2d`. Its generic geometric group is `A_(2d)`. If

\[
 C_d=(-1)^{d+1}(2d-1),
\tag{4}
\]

is a square in `K_d`, its arithmetic group is also `A_(2d)` and the normal closure is regular over `K_d`. Otherwise the arithmetic group is `S_(2d)` and its full constant field is exactly `K_d(sqrt(C_d))`.

In this particular quadratic coefficient field, `C_d` is square precisely when `d` is odd and `2d-1` is a rational square. Multiplying the first coordinate by `-1` gives determinant `1` without changing any field extension or group.

**Proof.** We prove the lift identities, the field identification, and the monodromy separately.

### The seed and the polynomial lift

Consider

\[
 \varphi(W)=aW^2+bW^d.
\]

Direct substitution in (1), using `D^2=-3(2d-1)`, gives

\[
 a+b=1,\qquad
 I:=\int_0^1\frac{\varphi(W)^2}{W^2}\,dW
   =\frac{a^2}{3}+\frac{2ab}{d+1}+\frac{b^2}{2d-1}=0.
\tag{5}
\]

For clarity, eliminating `b=1-a` gives the exact identity

\[
 3(d+1)(2d-1)I
 =2(d-2)^2a^2+6(d-2)a+3(d+1).
\tag{6}
\]

The quadratic on the right has discriminant `-12(d-2)^2(2d-1)`. Thus (1) supplies its two conjugate roots directly.

Define the canonical chain

\[
\begin{split}
A_1(W)&=4aW+\frac{2db}{d-1}W^{d-1},\\
A_0(W)&=\frac{8a^2}{3}W^3
 +\frac{4d^2ab}{(d-1)(d+1)}W^{d+1}
 +\frac{2d^2b^2}{(d-1)(2d-1)}W^{2d-1}.
\end{split}
\tag{7}
\]

They satisfy `A_1'=2 phi'/W`, `A_0'=phi' A_1`, and both have zero constant term. Set

\[
 T=\Lambda w,\quad W=\Lambda w^2,\quad
 Q=T+\varphi(W),\quad
 P=\frac{T^2}{W}+A_1(W)T+A_0(W).
\tag{8}
\]

Here `T^2/W=Lambda`, so (8) is polynomial in `Lambda,w`. Substitution gives exactly `Q=Lambda A`, `P=Lambda B`, with (2). Also

\[
 \Pi(W,Y)=A_0(W)+A_1(W)(Y-\varphi(W))
                    +\frac{(Y-\varphi(W))^2}{W}
\tag{9}
\]

satisfies `P=Pi(W,Q)` and

\[
 \Pi_W=-\frac{(Y-\varphi(W))^2}{W^2},\qquad
 \Pi(W,0)=-\int_0^W\frac{\varphi(s)^2}{s^2}\,ds.
\tag{10}
\]

At `(Lambda,w)=(1,-1)` one has `T=-1`, `W=1`. Equations (5) and (10) give `Q=P=0`, hence `A=B=0`. On setting `x=0` in their substitutions into (3), these are exactly the same constant values, so both numerators are divisible by `x`. Thus (3) is polynomial on all affine three-space. Its highest terms in (2) give the stated component degrees.

### The Jacobian and the embedded extension

The canonical-chain identities imply

\[
 \operatorname{Jac}_{T,W}(P,Q)=\frac{T^2}{W^2},\qquad
 \operatorname{Jac}_{\Lambda,w}(T,W)=W.
\]

Consequently `Jac_(Lambda,w)(P,Q)=Lambda`. Equivalently,

\[
 \Lambda(B_\Lambda A_w-B_wA_\Lambda)+BA_w-AB_w=1.
\tag{11}
\]

Writing `S=x Lambda`, the coordinate change `(x,y,z)->(Lambda,w,x)` has determinant `x^2`. The target change from `(F_1,F_2,S)` to `(Q=SF_1,P=SF_2,S)` has determinant `S^2`. The determinant of `(Q,P,S)` with respect to `(Lambda,w,x)` is `-Lambda^2`. These identities give `det DF_d=-1` on a dense open set and hence everywhere.

The substitution (8) is birational: `w=W/T` and `Lambda=T^2/W`. The element `S` is independently transcendental, with

\[
 x=S/\Lambda,\quad y=(\Lambda-1)/x,\quad z=(w+1)/x.
\]

Thus the source and target fields of (3) identify with

\[
 K_d(W,Y,S)\supset K_d(P,Y,S),\qquad P=\Pi(W,Y).
\tag{12}
\]

For generic `Y!=0`, (9) has poles of orders `2d-1` at infinity and `1` at zero. Its leading coefficient is `-b^2/(2d-1)` and its coefficient of `W^(-1)` is `Y^2`. This proves the generic degree `2d` in (12).

### Full geometric monodromy

Over an algebraic closure of `K_d(Y)`, the `d` roots of `phi(W)=Y` are nonzero and distinct. They are all the finite critical points of (9), each with derivative of order two and local degree three. Their critical values are `A_0(W)`.

These values are distinct for generic `Y`. Indeed the degree of the map
`W -> (phi(W),A_0(W))` divides both `deg phi=d` and `deg A_0=2d-1`, so it is one. Thus it is birational onto its image; outside a finite exceptional set of `Y`, no two critical points have the same value.

The resulting connected cover has `d` finite branch cycles, each an individual three-cycle, and infinity profile `(2d-1,1)`. Its finite branch cycles generate its transitive monodromy group. A transitive permutation group generated by three-cycles is the full alternating group: their supports form a connected three-uniform hypergraph, and adjoining a support to an already obtained alternating group enlarges that group to the alternating group on their union. For an overlap of one symbol, conjugating `(a u v)` to `(b u v)` first supplies a three-cycle on two old symbols and one new symbol; the remaining new symbol is then added. An overlap of two symbols directly adds the new symbol. This proves the claim by induction from `A_3`.

This establishes `A_(2d)` over the algebraic closure of the slice field `K_d(Y)`. The ambient geometric group and arithmetic constants require the next calculation.

### The exact quadratic constant

Use an independent target coordinate `Z`, and put

\[
 f(W)=W(\Pi(W,Y)-Z),\qquad
 g(W)=Y-\varphi(W),\qquad R=\operatorname{Res}_W(f,g).
\]

The degree of `f` is `N=2d`, its leading coefficient is `L=-b^2/(2d-1)`, and its constant coefficient is `Y^2`. Equation (10) gives

\[
 Wf'-f=-g^2.
\tag{13}
\]

At each root `alpha` of `f`, therefore, `alpha f'(alpha)=-g(alpha)^2`. Using the product of the roots, `Y^2/L`, and `deg g=d=N/2`, the discriminant-resultant formula gives the exact identity in `K_d(Y,Z)`

\[
 \boxed{\operatorname{Disc}_W(f)
 =(-1)^{d+1}(2d-1)
       \left(\frac{R}{bY}\right)^2.}
\tag{14}
\]

In particular there is no nonconstant quadratic twist in `Y`. The ambient geometric group over `K_dbar(Y,Z)` is contained in `A_N` by (14), and contains the slice group already proved to be `A_N`. It equals `A_N`; adjoining the independent `S` in (12) preserves this assertion. The arithmetic group contains the geometric group, so it is `A_N` or `S_N` exactly as (4) prescribes. In the second case the square root in (14) supplies `sqrt(C_d)`, and the index two leaves no larger constant field.

For a rational number to become square in `Q(D)`, either it or its quotient by `D^2` must be square over `Q`. Here

\[
 C_d/D^2=(-1)^d/3
\]

is never a rational square. Thus `C_d` is square in `K_d` precisely when `d` is odd and `2d-1` is a rational square. This completes the proof. ∎

## Corollary 2. One quadratic field for all admissible degrees

For every odd `m>=3`, set `d=(m^2+1)/2` and choose `D=m sqrt(-3)`. Formula (1) becomes

\[
 a=\frac{-3-m\sqrt{-3}}{m^2-3},\qquad
 b=\frac{m^2+m\sqrt{-3}}{m^2-3}.
\]

The polynomial map (3) is defined over `E=Q(sqrt(-3))`, has generic degree `N=m^2+1`, and has regular arithmetic and geometric monodromy `A_N`. Its component degrees are `(3N-3,6N-9,3)`.

**Proof.** Such `d` is odd and `2d-1=m^2`, so `C_d=m^2` is square. Apply Theorem 1. ∎

For `N=10` the two seed coefficients are `a=(-1-sqrt(-3))/2`, `b=(3+sqrt(-3))/2`, and the map degrees are `(27,51,3)`. For `N=26` they are `a=(-3-5sqrt(-3))/22`, `b=(25+5sqrt(-3))/22`, and the map degrees are `(75,147,3)`. These are specializations of one formula, not isolated coefficient searches.

## Corollary 3. Alternating whole fibres on one affine target line

For each admissible `N=m^2+1` in Corollary 2, one fixed map `F_d` over `E` has infinitely many complete degree-`N` field fibres with normal-closure group `A_N` along one affine target line. Its ambient arithmetic and geometric group is also `A_N`.

More precisely, choose a nonzero `Y_0` in `Q` outside the finite exceptional set where the critical points or critical values of `Pi(W,Y_0)` collide. On the literal line

\[
                  \Gamma(Z)=(Y_0,Z,1)
\]

and away from its finite branch set, the full scheme-theoretic fibre algebra is

\[
 E[W]/\bigl(W(\Pi(W,Y_0)-Z)\bigr).
\tag{15}
\]

**Proof.** The finite exceptional set exists by the birational critical-image argument in Theorem 1 and separability of `phi(W)-Y`. The slice has distinct three-cycle finite branches, hence geometric `A_N`; (14), with `C_d=m^2`, makes its arithmetic group `A_N` too and proves regularity over `E`.

We check the whole fibre algebra. In an actual fibre above `(Y_0,Z,1)`, the identity `x Lambda=1` makes both `x` and `Lambda` units. The identity `Q=Y_0` makes `w` a unit, because `Q=Lambda w+phi(Lambda w^2)` is divisible by `w`. Thus `T=Lambda w` and `W=Lambda w^2` are units and satisfy the polynomial equation in (15).

Conversely, in the algebra (15), `W` is a unit because the polynomial's constant coefficient is `Y_0^2`. Put `T=Y_0-phi(W)`. Off the branch set, the polynomial is separable, so (13) makes `T^2=-W f'(W)` a unit. Then set

\[
 \Lambda=T^2/W,\quad w=W/T,\quad x=W/T^2,
 \quad y=(\Lambda-1)/x,\quad z=(w+1)/x.
\tag{16}
\]

All denominators are units. These formulas recover `1+xy=Lambda`, `xz-1=w`, `x Lambda=1`, `Q=Y_0` and `P=Z`, hence every polynomial coordinate equation of the fibre. The two constructions are inverse as algebra maps, proving (15) for the entire scheme rather than selected points.

Hilbert irreducibility for the regular `A_N` extension over the number field `E` now gives infinitely many `Z` in `E` for which (15) is a field with full normal-closure group `A_N`. The Hilbert theorem is the same retained number-field input used for [alternating pencil specializations](../../arithmetic/ALTERNATING_PENCIL_SPECIALIZATIONS_AND_KELLER_FIBRES.md). ∎

This uses a new three-cycle pencil over `E`. The [fixed rational degree-ten compiler](FIXED_KELLER_MAP_ALTERNATING_FAMILY.md) carries the original rational double-transposition pencil and has ambient `S_10`; the two fibre-family statements concern different pencils and coefficient fields.

## Relation to the paired-pencil and plane questions

The quotient covers here have individual **three-cycle** finite inertia. The knot-derived pencils have **double-transposition** finite inertia. Both have alternating monodromy and poles `(N-1,1)`, but their different passports make these different constructions. Corollary 2 does not give rational degree-26 paired-critical pencils or their rational descent.

There is an elementary real obstruction in this particular cyclic-cone chart: for a nonzero real seed, the integral in (5) is strictly positive. Thus this level-two seed construction has no real-coefficient member. This is the real lift-form obstruction already present in the cited primary paper, not an obstruction to arbitrary coordinate forms of polynomial Keller maps.

The earlier [fixed-Laurent normalization-unit obstruction](ALTERNATING_LAURENT_FEASIBILITY.md) concerns a fixed nonzero pole coefficient over a retained affine base. Here the pole coefficient is `Y^2` and vanishes on `Y=0`; the target invariant products are `P=F_2F_3`, `Q=F_1F_3`. Inverting `W` is valid after inverting `Y`, but that gives no nonconstant unit over the whole polynomial target. The hypotheses of the earlier obstruction do not apply to (12).

The theorem answers generic alternating polynomial realization in dimension three over number fields. A plane polynomial realization or universal plane obstruction remains open here. No rational-coefficient descent, formal verification, independent external replay, or external human review is asserted.

## Exact evidence and attribution

The retained [primary PDF](sources/kistner-shaska-graded-keller-2608.02863v1.pdf) has SHA-256 `3cd7e239d8b01c79239789a97c26e466208d9cb216d9773574562c468c915a2e`. Its Section 3, Propositions 3.2, 3.5, 3.6 and 3.9 and Theorem 3.11 give the canonical chain, lift conditions and existence; Section 6 computes monodromy. Equations (1)–(14) above supply a direct specialization and proof of every assertion needed here rather than importing an uninspected general claim.

The [checker](verify_uniform_even_alternating_cone_keller.py) checks (5) and (6) symbolically in `d`, the canonical chain and exact polynomial identities in bounded degrees, and the actual expanded three-variable polynomial lift and Jacobian in degrees six and ten. It also checks the full discriminant identity (14) in degree six, the algebraic reconstruction identities (16), and rejects corrupted seed coefficients. These are finite controls for the written uniform proof. The [receipt](uniform_even_alternating_cone_keller_verification.json) and [source ledger](uniform_even_alternating_cone_keller_source_ledger.json) pin their scope and source bytes.

```sh
.venv/bin/python research/jacobian/monodromy-theory/verify_uniform_even_alternating_cone_keller.py \
  --write-receipt research/jacobian/monodromy-theory/uniform_even_alternating_cone_keller_verification.json
```
