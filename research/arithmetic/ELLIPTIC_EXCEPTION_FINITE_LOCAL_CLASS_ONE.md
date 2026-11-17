# Finite-place positivity for the cyclotomic quadratic class

The quadratic field `Q(sqrt(-3))` has a positive finite-place target
coefficient on the existing elliptic slice for every fixed finite set
of places. This follows from the formal neighbourhood of the identity;
it does not require an everywhere locally soluble seed in this class.

Retain `E`, `Gamma=E(Q)`, the four-class map `alpha`, the actual target
height and coefficient `C_E` from
[EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md). The class `a=1`
corresponds to `K_1=Q(sqrt(-3))`. Let `N_{S,1}(B)` count the admissible
full quadratic-times-cubic field-fibre targets of this class, of height
at most `B`, locally soluble at the prescribed finite set of places
`S`. Their rational source-point set is empty because both field
components have degree greater than one.

**Theorem ECO1 (positivity in class one for every finite set).** For
every fixed finite `S`, the Haar coefficient of EFA2 satisfies

\[
 \boxed{\mu_{S,1}\in\mathbb Q,\qquad
                         0<\mu_{S,1}\le\frac14.}    \tag{1}
\]

Moreover

\[
 \boxed{N_{S,1}(B)=\mu_{S,1}C_E\log B
                            +O_{E,S}(\sqrt{\log B}).}          \tag{2}
\]

In particular there are infinitely many such full field-fibre targets
for every fixed finite `S`.

**Proof.** Write `S_f` for the finite primes. EQA1 in
[the quantitative finite-place note](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md)
supplies an open product inside the actual joint rational closure,

\[
 \{1_A\}\times\prod_{p\in S_f}E_{d_p}(p)
                                      \subseteq H_S.          \tag{3}
\]

For each `p`, choose a point of sufficiently large formal depth `3k_p`
whose parameter `t=-X/y` is exactly `p^{3k_p}`. The formal parameter
identifies a sufficiently small neighbourhood of the identity with
a `p`-adic ball, so these choices exist in `E_{d_p}(p)`. They may be
chosen arbitrarily close to the identity and away from the finite
inadmissible set.

Each chosen point supplies a local source point. At two, the formal
criterion of ETF2 is divisibility of its depth by three. At three,
ETC1 requires the parameter itself to be a cube, as it is here. At
five, ETF3 again requires divisibility of formal depth by three.
At any other prime `p=2 mod3`, ELV1 gives
`v_p(128x-25865)=-2(3k_p)`, divisible by three. Primes `p=1 mod3`
and the real place are automatic. These exact tests are retained in
[ETF2--ETF3](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md),
[ETC1](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md), and
[ELV1](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md).

Thus the tuple in (3) lies in `U_{S,1}`. That set is open in `H_S`
by EFA1, so its nonemptiness gives positive Haar measure. Projection
to the four-element class group makes its measure at most `1/4`.
EQA2 makes the measure rational. For `S_f` empty these assertions
reduce directly to `mu_{S,1}=1/4`.

The dense rational image also explains the field-fibre assertion
directly. A small basic open coset about the chosen tuple pulls back
to a finite-index coset of `Gamma`, entirely in class one and satisfying
the prescribed local conditions. Removing the finite global
inadmissible and reducible-cubic exceptions leaves infinitely many
full quadratic-times-cubic field fibres. No claim that the exhibited
rational points form a Mordell--Weil basis is used.

Finally apply Theorem ESQ3 of
[the sharp finite-place error note](ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md)
to obtain (2). Its elliptic-logarithm input and assurance scope are
retained there; positivity here uses only the open product and the
explicit local tests. \(\square\)

Together with the certified class `10p_0`, this proves finite-place
positivity for two of the four quadratic fields for every finite
`S`. Whether the other two classes have this property remains
undecided here. The coefficient can shrink as `S` grows; this proof
gives neither an all-places target nor an all-primes limiting measure
or tail estimate for class one.

This is a written deduction from the cited canonical results. No new
finite experiment, formal verification or external review is claimed.
