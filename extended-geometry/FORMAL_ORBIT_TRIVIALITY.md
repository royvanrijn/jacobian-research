# Formal orbit triviality and the algebraization gap

The cubic dual-number calculation is an instance of a general statement, but
the general statement does not stop at first order.  A Keller map is formally
invertible for composition near the identity even when it is not globally
invertible as a polynomial map.  Consequently the quotient by the full
polynomial source-automorphism functor has no finite-order deformation data.

The operational consequence is immediate: a dual-number or higher Artin
direction found inside a bounded-support ansatz is a feature of that slice,
not genuine global moduli.  Further searches must impose bounded-degree or
factorization complexity, study algebraization over a reduced base, or use
global normalization-boundary invariants.  The
[complexity-filtered follow-up](COMPLEXITY_FILTERED_CONTACT.md) is the direct
continuation adopted in this repository.

Throughout, `k` is a field.  Characteristic zero is needed only for the
reduced-curve interpolation theorem in Section 4.  Write

\[
 F:\mathbb A_k^n\longrightarrow\mathbb A_k^n,
 \qquad \det DF=c\in k^\times.
\]

We distinguish
`Aut_1(A^n)={A:det DA=1}` from
`SAut(A^n)=<Ga-subgroups>`.  No equality of these groups is assumed.

## 1. First-order theorem

For every polynomial deformation vector `H`, put

\[
 V=(DF)^{-1}H=c^{-1}\operatorname{adj}(DF)H.
\]

Then `V` is polynomial and `H=DF V`.  Differentiating
`F\circ(\mathrm{id}+\epsilon V)` gives

\[
 \left.\frac d{d\epsilon}\right|_{0}
 \det D(F+\epsilon H)=c\,\operatorname{div}V.       \tag{1.1}
\]

Hence, for fixed determinant `c`,

\[
 \boxed{
 T_F\operatorname{Kell}_c
 =DF\bigl(\operatorname{Der}_{\operatorname{div}=0}
 k[x_1,\ldots,x_n]\bigr).}                         \tag{1.2}
\]

In dimension three and characteristic zero, the divergence-free shear lemma
from [the cubic orbit audit](CUBIC_DUAL_NUMBER_ORBIT_TANGENCY.md) writes every
field on the right of (1.2) as a finite sum of locally nilpotent derivations.
Composing their exponential flows gives a reduced polynomial `A^1`-curve in
`SAut(A^3)` with the prescribed tangent.  Thus every fixed-Jacobian
first-order deformation of a Keller map on `A^3` is source-orbit tangent.

## 2. Formal source-triviality theorem

Let `(R,m)` be a local Artin `k`-algebra with residue field `k`.  Let

\[
 \mathcal F:\mathbb A_R^n\longrightarrow\mathbb A_R^n
\]

be a polynomial map whose reduction modulo `m` is `F`.

### Theorem 2.1

There is a unique polynomial automorphism

\[
 \alpha\in\operatorname{Aut}_R(\mathbb A_R^n),
 \qquad \alpha\equiv\mathrm{id}\pmod m,
\]

such that

\[
 \boxed{\mathcal F=F\circ\alpha.}                    \tag{2.1}
\]

If

\[
 \det D\mathcal F=c,
\]

then `det D alpha=1`.

### Proof

Use the finite filtration by powers of `m`.  Suppose `alpha_r` solves (2.1)
modulo `m^r`.  The error modulo `m^(r+1)` is a vector with coefficients in
`m^r/m^(r+1)`.  Replacing `alpha_r` by `alpha_r+V` changes `F\circ alpha_r`
by

\[
 DF\,V\pmod {m^{r+1}}.
\]

Here `DF` may be evaluated at the identity: evaluating at `alpha_r` changes
it by an `m`-valued matrix, which vanishes after multiplication by the
`m^r`-valued correction.  The polynomial matrix `(DF)^(-1)` therefore gives
a unique `V`.
Induction terminates because `m` is nilpotent.  A polynomial endomorphism
congruent to the identity modulo a nilpotent ideal is an automorphism: its
inverse is obtained by the same finite induction.

Finally, the chain rule in (2.1) gives

\[
 \det D\mathcal F
 =c\,(\det D\alpha).
\]

This proves the last assertion.  Notice that neither this proof nor its
uniqueness assertion assumes that `F` is injective.  QED

For a one-parameter deformation

\[
 F_t=F+tH_1+t^2H_2+\cdots,
 \qquad
 \widehat\alpha_t=\mathrm{id}+tV_1+t^2V_2+\cdots,
\]

the recursion begins

\[
 V_1=(DF)^{-1}H_1,
\]

\[
 V_2=(DF)^{-1}
 \left(H_2-\frac12D^2F[V_1,V_1]\right).             \tag{2.2}
\]

It continues uniquely to every order.  Therefore

\[
 \widehat{\operatorname{Kell}}_{c,F}/
 \widehat{\operatorname{Aut}}_{1,\mathrm{id}}
\]

is the one-point functor.  In particular, ordinary second-order lifting in
the full polynomial source quotient cannot detect moduli either.

## 3. Compatible formal trivialization

Applying Theorem 2.1 to every `k[t]/(t^(m+1))` produces compatible jets and
hence a unique formal automorphism

\[
 \widehat\alpha(t)\in
 \operatorname{Aut}_{k[[t]]}(\mathbb A^n),
 \qquad F_t=F\circ\widehat\alpha(t).                 \tag{3.1}
\]

For a fixed-Jacobian family it lies in formal `Aut_1`.  Equation (3.1) is a
coefficientwise statement: every coefficient of `t` is polynomial in the
source variables.  It does not assert that the series terminates in `t`, is
rational in `t`, or defines an automorphism over a reduced finite-type base.

This is the formal-to-algebraic gap in which the stable moduli can live.

## 4. Every finite special jet has a reduced representative in dimension three

Assume now that `char(k)=0`.  Let

\[
 \bar\alpha\in
 \operatorname{Aut}_{1,k[t]/(t^{m+1})}(\mathbb A^3),
 \qquad \bar\alpha_0=\mathrm{id}.
\]

### Theorem 4.1

There is

\[
 A(t)\in\operatorname{SAut}_{k[t]}(\mathbb A^3)
\]

whose reduction modulo `t^(m+1)` is `bar alpha`.

### Proof

Suppose a reduced polynomial curve `A_(r-1)(t)` matches through order `r-1`.
After composing with its inverse, the order-`r` discrepancy has the form

\[
 \mathrm{id}+t^rV_r\pmod {t^{r+1}}.
\]

The determinant-one condition gives `div(V_r)=0`.  By the divergence-free
shear lemma, write `V_r=D_1+...+D_s` with every `D_i` locally nilpotent.  Then

\[
 C_r(t)=\exp(t^rD_1)\circ\cdots\circ\exp(t^rD_s)
\]

is a polynomial `SAut` curve and

\[
 C_r(t)=\mathrm{id}+t^rV_r\pmod {t^{r+1}}.
\]

If the discrepancy was defined as `A_(r-1)^(-1) bar alpha`, composing
`A_(r-1)` on the right with `C_r` supplies the next coefficient without
altering the earlier ones.  Induction proves the claim.
QED

It follows that the unfiltered invariant

\[
 \operatorname{contactord}(F_t,\operatorname{Orb}(F))
\]

is infinite for every fixed-Jacobian family if order `m` merely asks for
some polynomial source--target automorphism curve matching modulo
`t^(m+1)`.  Source automorphisms already suffice.  The curve is allowed to
depend on `m`; Theorem 2.1 shows that even the compatible formal system is
always present.

## 5. What remains nontrivial

The meaningful question is whether the canonical formal trivializer (3.1)
algebraizes over the reduced parameter base.  Exact source-orbit membership
requires a single regular automorphism family, not a different increasingly
complicated curve for every finite jet.

This suggests replacing raw contact order by a filtered contact profile.  A
prototype is

\[
 \kappa_m(F_t,F)=\min_A\max_i\deg_{x,t}A_i,           \tag{5.1}
\]

where `A` ranges over polynomial special-automorphism curves satisfying
`F_t=F\circ A mod t^(m+1)`.  More intrinsic versions should record the
ind-group degree filtration, stabilization dimension, or the number and
degree of shear factors.  The sequence is finite termwise by Theorem 4.1;
its growth measures the cost of approximating the formal orbit.

This proposal is developed with a two-sided degree filtration and tested on
the degree-five stable-moduli arc in
[complexity-filtered orbit contact](COMPLEXITY_FILTERED_CONTACT.md).
The finite-type automorphism strata, constructible bounded-orbit relation,
full parameter-height/pole-order spectrum, and the precise limitation of the
pointwise-to-family bridge are developed in
[ind-group constructibility](IND_GROUP_CONSTRUCTIBILITY.md).

There are three natural levels:

1. **Artin quotient:** trivial by Theorem 2.1.
2. **Finite reduced contact:** unfiltered contact is infinite by Theorem 4.1.
3. **Exact reduced orbit:** algebraization may fail and global stable
   normalization-boundary invariants can vary.

The decorated-normalization construction belongs to the third level.
Normalization is formed from the reduction and is known to commute with
smooth base change; it is not automatically a nilpotent-base deformation
functor.  Consequently the generically unramified decorated seed map does not
contradict formal source triviality.  Rather, the two results together show
that the positive-dimensional stable invariant is genuinely global and
algebraic, not a finite-order quotient invariant.

The bounded exact regression
[`verify_formal_orbit_triviality.py`](../scripts/verify_formal_orbit_triviality.py)
reconstructs the first three coefficients of the canonical source
trivializer for the foundational noninjective Keller map and checks its
determinant-one identity.
