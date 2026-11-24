# Cubic moment-kernel transfer and an exact resultant valuation

Internal written proof, 21 November 2026. Repository base:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Result and source boundary

The main result is a uniform determinant theorem. For every order `n`,
arbitrary monic coefficient perturbations of two explicit consecutive model
polynomials, with three-adic errors of valuation **strictly greater than one**,
have a nonzero resultant whose valuation is exactly the degree of their common
residual factor. That residual degree need not be small. The theorem allows
ramified perturbation coefficients; it does not incorrectly declare their
actual roots unramified.

Applications include a cubic unit-deformation theorem for entire moment
kernels, and ordinary factorial families obtained by composing with linear
forms in arbitrarily many variables. The linear forms satisfy an explicit
arithmetic coefficient condition. Integer Gamma shapes and low-degree markers
give additional applications. This does not close arbitrary supports or
unrestricted Strong Factorial.

The existing
[adjacent powers-of-three proof](STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md)
supplies the cubic model, its common-root simplicity identity, and its trace
identity. Their precise use is reconstructed below. The new statements here
are the sharp perturbation threshold, exact resultant valuation, kernel
transfer, and its stated applications. The retained
[dyadic theorem](DYADIC_BINOMIAL_SFC.tex)
already proves a different, arbitrary-two-adic-unit transfer; its application
below is explicitly a consequence. The falling-factorial deformation theorem
in [PRIME_POWER_TWO_TERM_SFC.tex](PRIME_POWER_TWO_TERM_SFC.tex)
concerns squarefreeness for constant-plus-power phases and is not substituted
for the both-nonconstant kernel used here.

## 1. The model and its simple common residue disks

Define integer polynomials

\[
 Q_n(x)=\sum_{k=0}^n(-1)^k\binom{n+2k}{3k}x^{n-k},\qquad
 \sum_{n\ge0}Q_n(x)t^n=\frac{(1-xt)^2}{(1-xt)^3+t}.
 \tag{1}
\]

Write `R_n=Q_n mod3` and `g_n=gcd(R_n,R_(n+1))`, made monic.
Every common residual root is nonzero, since `R_n(0)=(-1)^n`.

For completeness, let

\[
 \sum_{m\ge0}A_m(y)t^m=\frac1{1+t-yt^3},\quad
 A_0=1,\quad A_{-1}=A_{-2}=0.
\]

In characteristic three,

\[
 A_m=-A_{m-1}+yA_{m-3},\quad
 R_n=A_n(x^3)+xA_{n-1}(x^3)+x^2A_{n-2}(x^3),
\]
\[
 R_n'=A_{n-1}(x^3)-xA_{n-2}(x^3),
\]
\[
 x^2A_{n-2}(x^3)=(x-1)R_n+xR_n'-R_{n+1}.
 \tag{2}
\]

If a common residual root also killed `R_n'`, (2) would force three
consecutive `A` values to vanish. The recurrence, run backwards using
`x^3!=0`, would then force `A_0=0`. At `n=1`, the contradiction is directly
`R_1'=1`. Thus every root of `g_n` is simple in `R_n`, and `g_n` is squarefree.
This asserts neither coprimality of the reductions nor squarefreeness of
each entire `R_n`.

Let `tau_i` be the three roots of `tau^3+tau^2-x`, and set

\[
 S_j(x)=\sum_i\tau_i^j,\qquad
 V_j(x)=\sum_{i<k}(\tau_i\tau_k)^j.
\]

These are integer polynomials. Their defining recurrences are

\[
 S_0=3,\ S_1=-1,\ S_2=1,\quad S_j=-S_{j-1}+xS_{j-3},
\]
\[
 V_0=3,\ V_1=0,\ V_2=2x,\quad V_j=xV_{j-2}+x^2V_{j-3}.
\]

For every `N>=1`,

\[
 Q_N+2xQ_{N-1}=S_{3N}
   =S_N^3-3S_NV_N+3x^N.
 \tag{3}
\]

The first equality follows by taking the logarithmic derivative of
`D=(1-xt)^3+t=product_i(1-tau_i^3 t)`:
`-tD'/D=(1+2xt)sum Q_n t^n-1`. The second is the third-power Newton
identity. These are polynomial identities, including at repeated cubic roots.

## 2. Sharp local separation and the resultant theorem

Normalize `v(3)=1`, and work in a complete algebraically closed valued
extension of `Q_3`. Coefficients and points may be ramified. A unit means
valuation zero.

**Lemma 1 (unit-disk separation).** For every `n>=1` and every unit `alpha`,

\[
 \min\{v(Q_n(alpha)),v(Q_{n+1}(alpha))\}\le1.
 \tag{4}
\]

**Proof.** If the residue of `alpha` is not common to both reductions,
one valuation is zero. Otherwise (2) makes it a simple lower root. Lift it
to a root `beta` of the integer polynomial `Q_n` in the finite unramified
extension of `Q_3` containing that residue field. This follows by ordinary
Newton iteration, since the residual derivative is a unit.

Put `N=n+1`. The value `Q_N(beta)` has positive valuation. By (3),
`S_(3N)(beta)=Q_N(beta)`. Reducing the second equality of (3) modulo three
shows `S_N(beta)=0 mod3`. This is a statement in the unramified valuation
ring, whose uniformizer is three. Consequently

\[
 Q_N(beta)\equiv3beta^N\pmod9,
 \qquad v(Q_N(beta))=1.
 \tag{5}
\]

Polynomial division gives `Q_n(alpha)=(alpha-beta)U(alpha)`, where `U`
has integral coefficients and `U(alpha)` is a unit: its residue is the
nonzero lower derivative. If `v(alpha-beta)<=1`, this proves (4). If
`v(alpha-beta)>1`, the integral divided difference of `Q_N` gives
`Q_N(alpha)-Q_N(beta)` valuation greater than one, so (5) gives
`v(Q_N(alpha))=1`. This proves the other case. The actual point `alpha`
was never asserted to be unramified. QED.

**Theorem 2 (exact resultant under three-adic perturbations).** Let `n>=1`.
Let `P` and `T` be monic polynomials of degrees `n` and `n+1`, respectively,
with integral coefficients in the valued extension. Suppose every coefficient
of `P-Q_n` and `T-Q_(n+1)` has valuation strictly greater than one. Then

\[
 \boxed{\ v\operatorname{Res}(P,T)=\deg g_n.\ }
 \tag{6}
\]

In particular the resultant is nonzero. For rational or unramified integral
coefficients, coefficientwise congruence modulo nine is sufficient.

**Proof.** The constant coefficients remain units. All roots of `P` are
therefore units, by monicity and the valuation of its constant coefficient.
For a lower root whose residue is not common to the two reductions,
`T(alpha)` is a unit and contributes zero to the resultant valuation.

For each common residual root, lower residual simplicity gives exactly one
root `alpha` of `P` in that residue disk, of multiplicity one. Let `beta`
be the unramified root of the model used in Lemma 1. Since
`P(alpha)=0`, coefficient closeness gives `v(Q_n(alpha))>1`. The same
unit divided difference gives `v(alpha-beta)>1`. Equation (5) and the
coefficient closeness of `T` then give `v(T(alpha))=1`.

There are exactly `deg g_n` such disks over the algebraic closure of the
residue field. The product formula for the resultant of monic polynomials
sums these valuations over all lower roots, with multiplicity, giving (6).
QED.

The common residual degree can grow. For example it is `10` at `n=39`
and `13` at `n=52`. Formula (6) still applies with the same coefficient
precision. No bounded residual rank, bounded common-factor degree, or
uniformly simple full reduction is assumed.

## 3. A cubic unit-deformation theorem for moment kernels

Let `(u_L)_(L>=0)` be rational three-adic units, and suppose there is a
rational three-adic unit `zeta` such that

\[
 \frac{u_L}{u_0}\equiv\zeta^L\pmod9\qquad(L\ge0).
 \tag{7}
\]

Define a weighted functional `Lambda_u(X^L)=L!u_L` and

\[
 M_n(t)=\sum_{k=0}^n\binom nk(n+2k)!u_{n+2k}t^k.
 \tag{8}
\]

**Theorem 3 (cubic kernel transfer).** The polynomials `M_n,M_(n+1)`
are coprime over `Q` for every `n>=1`. Thus every nonzero
`X(A+BX^2)` has a nonzero weighted moment in every consecutive positive
pair, over every characteristic-zero coefficient field. In addition,

\[
\boxed{
 v_3\operatorname{Res}(M_n,M_{n+1})
 =n(n+1)+(n+1)v_3(n!)+n v_3((n+1)!)+\deg g_n.
}
 \tag{9}
\]

**Proof.** Put `eta=-2zeta^2`, a unit, and use the same change of variable
at both orders:

\[
 G_n(x)=\frac{x^n}{n!u_n}\,
                 M_n\!\left(\frac1{3eta x}\right).
 \tag{10}
\]

Let `U(k)=product_(j=1)^k(3j-2)(3j-1)`. Its exact factorial identity is
`(3k)!=3^k k!U(k)`, and `U(k)=2^k mod9`. The coefficient of `x^(n-k)`
in (10) is exactly

\[
 \binom{n+2k}{3k}U(k)\frac{u_{n+2k}}{u_n}\,eta^{-k}
 \equiv(-1)^k\binom{n+2k}{3k}\pmod9.
 \tag{11}
\]

All denominators are units. Thus `G_n` is integral and monic, with unit
constant coefficient, and `G_n=Q_n mod9`. Apply Theorem 2.

The reversal/scaling law, with `c=1/(3eta)` and `N=n+1`, is

\[
 \operatorname{Res}(G_n,G_N)
 =\frac{c^{nN}\operatorname{Res}(M_n,M_N)}
        {(n!u_n)^N(N!u_N)^n}.
\]

There is no reversal sign because `nN` is even. Taking valuations gives
(9). Rational coprimality gives a rational Bezout identity and hence
noncollision over every characteristic-zero field. If either phase
coefficient vanishes, the weighted monomial moments are nonzero because
each `u_L` is nonzero. QED.

**The extra unit scaling matters.** At `u_L=1`, (10) uses `eta=-2`,
and its first polynomials are

\[
 G_1=x-1,\qquad G_2=x^2-4x+10\equiv Q_2\pmod9.
\]

With the older, fixed normalization `eta=1`, the lower-exponent-one
coefficient differs from the model. The additional common unit scaling
is precisely what repairs that discrepancy. Thus including exponent
`3^0=1` in the applications below is justified, not inferred from the
earlier fixed-normalization statement.

The local version of Theorem 3 permits ramified unit sequences whose
errors in (7) have valuation greater than one. Theorem 2 applies to their
normalized polynomials without requiring the perturbed roots to be
unramified. Rationality is used only for transfer to an arbitrary
characteristic-zero coefficient field.

## 4. An ordinary multivariate composition family

Let `S=sum_(i=1)^m c_i X_i`, where the `c_i` are rational and three-adically
integral. Suppose

\[
 \prod_{i=1}^m(1-c_i t)\equiv1-c t\pmod9,
 \qquad c\in\mathbb Z_{(3)}^\times.
 \tag{12}
\]

This is a finite coefficient test. Equivalently, exactly one of the
`c_i` is a three-adic unit, all the others are divisible by three,
and the sum of those other coefficients is divisible by nine.

To see the equivalence, first reduce (12) modulo three. Degree over the
field `F_3` forces exactly one nonzero coefficient. The quadratic
elementary symmetric function is then the unit coefficient times the
sum of the others, modulo nine. All higher elementary symmetric functions
already contain at least two coefficients divisible by three. The converse
is the same calculation. One may take `c=sum c_i` in (12).

**Corollary 4 (ordinary many-variable transfer).** For every `r>=0`,
`q=3^r`, and every `A,B` in a characteristic-zero field, every nonzero phase

\[
 F(X)=S^q(A+B S^{2q})
 \tag{13}
\]

has a nonzero **ordinary factorial** moment in every consecutive positive
pair. This gives SFC for these phases, with no bound on dimension or actual
monomial count. The coefficient restriction (12) is part of the theorem.

**Proof.** Write `h_L(c_1,...,c_m)` for the complete homogeneous polynomial.
Direct multinomial expansion gives

\[
 \Lambda(S^L)=L!h_L,\qquad
 \sum_{L\ge0}h_Lt^L=\prod_i(1-c_i t)^{-1}.
 \tag{14}
\]

Thus `h_L=c^L mod9` and every `h_L` is a unit. Let

\[
 V_q(L)=\prod_{j=0}^{r-1}U(3^jL),\qquad
 C_q=3^{(q-1)/2},
\]

with the empty product equal to one at `r=0`. Exact factorial splitting
gives

\[
 (qL)!=C_q^L L!V_q(L),\qquad
 V_q(L)\equiv2^{(q-1)L/2}\pmod9.
\]

Consequently

\[
 \Lambda(S^{qL})=C_q^L L!u_L,
 \quad u_L=V_q(L)h_{qL}
 \equiv\left(2^{(q-1)/2}c^q\right)^L\pmod9.
 \tag{15}
\]

Expanding (13), the harmless factor `C_q^n` and common parameter scaling
`B/A -> C_q^2 B/A` reduce its moment pair to Theorem 3. Degenerate
coefficients give nonzero moments directly from (15). When both are
nonzero the two distinct homogeneous degrees ensure at least two actual
monomials; otherwise every moment is nonzero. Hence the actual support-size
window also follows. QED.

For example `S=X+3Y-3Z` satisfies (12). The conclusion holds for
`S^q(A+B S^(2q))` for every `q=1,3,9,...`, and one may add arbitrarily
many coefficients divisible by three whose sum is divisible by nine.
The proof does not assert a general moment-preserving linear change of
variables. The arithmetic identity (14) and condition (12) are essential.

### Dyadic counterpart, explicitly using the retained theorem

Suppose the rational coefficients of `S` are two-adically integral and
exactly one is odd. Then `h_L` in (14) is odd for every `L`. For `q=2^r`,
the retained factorial splitting gives `(qL)!=2^((q-1)L)L!W_q(L)` with
`W_q(L)` odd. Therefore the effective units `W_q(L)h_(qL)` satisfy the
arbitrary-unit hypothesis of `lem:units` in
[DYADIC_BINOMIAL_SFC.tex](DYADIC_BINOMIAL_SFC.tex).
It follows that every nonzero `S^q(A+B S^q)` has a two-moment window.
This is a transfer of that existing dyadic theorem, not a new scalar
dyadic result. Both multivariate conclusions retain their stated
arithmetic restrictions on the linear form.

## 5. Integer Gamma shapes and marked kernels

For `q=3^r`, `r>=1`, and an integer `1<=beta<=q`, put `delta=beta-1`.
The exact Gamma kernel is

\[
 (beta)_{qL}=C_q^L L! V_q(L)
              \prod_{j=1}^{delta}\left(1+\frac{qL}{j}\right).
 \tag{16}
\]

Every factor on the right after `C_q^L L!` is a three-adic unit.
Modulo nine, only `j=q/3` and `j=2q/3` can affect the last product:
all other `q/j` have valuation at least two. Thus

\[
 \prod_{j=1}^{delta}\left(1+\frac{qL}{j}\right)
 \equiv(1+3epsilon_delta)^L\pmod9,
\]

where `epsilon_delta=1` when `q/3<=delta<2q/3`, and zero otherwise.
If both exceptional factors occur, their linear terms add to
`3L(1+1/2)=9L/2`, which vanishes modulo nine. Theorem 3 therefore
proves the two-moment assertion for `X^q(A+B X^(2q))` for all these
integer Gamma shapes, uniformly in the terminal order.

Using `S=X_1+...+X_beta`, the exact identity
`Lambda(S^L)=(beta)_L` gives another ordinary many-variable family,
`S^q(A+B S^(2q))`, now with `1<=beta<=q`. This is a separate route
from (12); equal unit coefficients need not satisfy (12).

More generally, let a rational divided-power marker be

\[
 M(X)=\sum_{e=0}^{q-1}a_eX^e/e!,\quad
 C=\sum_ea_e\ne0,\quad a_e/C\in\mathbb Z_{(3)}.
\]

The normalized marked kernel is

\[
 \frac{\Lambda(M(X)X^{qL})}{C(qL)!}
 =\sum_e\frac{a_e}{C}\prod_{j=1}^e(1+qL/j)
 \equiv1+3hL\equiv(1+3h)^L\pmod9,
\]

where `h` is the sum of `a_e/C` over `q/3<=e<2q/3`, reduced modulo
three. The same transfer proves the marked two-moment assertion. Neither
positive marker coefficients nor positivity of the functional is assumed.
The mass and integrality conditions are indispensable hypotheses here.

## 6. Sharp controls and the obstruction to an all-prime shortcut

### Modulo three is insufficient, even for the cubic kernel

Set `u_L=1` except

\[
 u_3=220,\qquad u_4=-12860.
 \tag{17}
\]

All these rational units are congruent to one modulo three. Nevertheless
the weighted moments of `X(1-X^2/6)` at orders three and four both vanish.
With the baseline `u_L=1`, their normalized values `M_3/3!` and
`M_4/4!` are exactly `-219` and `12861`. Changing `u_3` and `u_4`
cancels these constants. The two moments use disjoint odd/even index
sets, so the changes do not interfere.

After the normalization (10), this gives monic polynomial perturbations
with errors of valuation at least one and a common root. Thus the strict
threshold in Theorem 2 cannot be replaced by valuation at least one.
This is a counterexample for modified moment kernels, not for the ordinary
factorial functional.

### The direct model for general odd primes

For an odd prime `p` define

\[
 Q_n^{(p)}(x)=\sum_{k=0}^n(-1)^k
      \binom{n+(p-1)k}{pk}x^{n-k},\quad
 \sum Q_n^{(p)}t^n=\frac{(1-xt)^{p-1}}{(1-xt)^p+t}.
 \tag{18}
\]

Here is the exact factorial normalization underlying this model. Put
`a=p^r`, `r>=1`, and

\[
 T_n(t)=\sum_{k=0}^n\binom nk(an+1)_{(p-1)ak}t^k,
 \qquad P_n(x)=x^nT_n(1/(p^a x)).
\]

Define `U_p(m)=product_(1<=j<=pm, p not dividing j)j`. Iterating
`(pm)!=p^m m!U_p(m)` gives the exact coefficient formula

\[
 [x^{n-k}]P_n
 =\binom{n+(p-1)k}{pk}U_p(k)
   \prod_{j=0}^{r-1}
       \frac{U_p(p^j(n+(p-1)k))}{U_p(p^j n)}.
 \tag{18a}
\]

The denominators are units. Since the original normalization divides only
by a power of `p`, the coefficients are integers. The leading coefficient
is one and the constant coefficient is the unit `U_p(an)`.
Writing `w=(p-1)!`, each nonmultiple block is `w mod p^2`: expanding
`product_(j=1)^(p-1)(pL+j)` to first order gives
`w[1+pL sum j^(-1)]`, and the reciprocal sum is zero modulo `p`.
The total unit exponent in (18a) is `ak`. Since `w=-1 mod p`,
`w^(p^r)=-1 mod p^2` for `r>=1`. Thus `P_n=Q_n^(p) mod p^2`.
The model is therefore mathematically natural, but the cubic separation
theorem does not extend to it merely by changing the prime.

Uniformly for every odd prime,

\[
 R_{2p-2}^{(p)}=(x+1)^{p-2}(x^p+x+1),
\]
\[
 R_{2p-1}^{(p)}=(x+1)^{p-2}(x^{p+1}-x-1),
\]
\[
 \gcd(R_{2p-2}^{(p)},R_{2p-1}^{(p)})=(x+1)^{p-2}.
 \tag{19}
\]

For a derivation, reduce (18) to
`(sum_(j=0)^(p-1)(xt)^j)/(1+t-x^p t^p)`. The coefficients of the
reciprocal denominator for `0<=m<2p` are `(-1)^m` when `m<p`, and
`(-1)^m[1-(m-p+1)x^p]` when `m>=p`. Summing gives the first two
formulas. The cofactor difference
`(x^(p+1)-x-1)-x(x^p+x+1)=-(x+1)^2` and the nonzero cofactor values
at `x=-1` prove the exact gcd. For every `p>=5`, the common root is
multiple in the lower polynomial. The previous unramified-lift argument
cannot be reused at that disk.

The analogous arbitrary-unit deformation assertion modulo `p^2` is in
fact false already at `p=5`. Direct exact evaluation gives

\[
 Q_{156}^{(5)}(8)\equiv25\pmod{125},\qquad
 Q_{157}^{(5)}(8)\equiv100\pmod{125}.
 \tag{20}
\]

Thus subtracting these constant values produces monic integral polynomials
in exactly the prescribed modulo-25 coefficient classes, with common root
eight. This is not just a failure of residual simplicity.

There is also a single coherent unit-kernel realization. Put `t=-1/960` and

\[
 E_m=\frac1{m!}\sum_{k=0}^m\binom mk(m+4k)!t^k.
\]

Take `u_L=1` except `u_156=1-E_156` and `u_157=1-E_157`. All `u_L`
are rational five-adic units congruent to one modulo 25. To verify this,
use `(5k)!=5^k k!U_5(k)` and `U_5(k)=24^k mod25`; at the chosen
phase, `E_m=8^(-m)Q_m^(5)(8) mod25`. The normalized values are
indeed zero modulo 25 by (20). The modified moments at orders 156 and
157 vanish exactly, since their coefficient index sets are disjoint
modulo four. Full rational values of the two modified units are retained
in the checker receipt.

These controls rule out the broad unit-only all-prime shortcut. They do not
disprove ordinary SFC, and do not exclude an argument that uses more than
the first two factorial-unit digits or additional arithmetic relations.

## 7. Evidence and remaining general problem

The standard-library [main checker](verify_sfc_cubic_unit_deformation.py) checks the model
identities, exact determinant valuations on finite representatives,
coefficient normalizations including `r=0`, the ordinary linear-form
kernel, Gamma and marker congruences, the all-prime residual formulas at
recorded prime guards, and both explicit deformation counterexamples.
The all-order and all-parameter conclusions use the written arguments;
finite examples do not supply those quantifiers.

Its [receipt](sfc_cubic_unit_deformation_verification.json) records eight
exact model resultants and eight coefficient perturbations, sixteen raw
unit-kernel resultants, nine ordinary cubic-composition cases, and the
explicit deformation counterexamples, together with the other stated finite
guards. The checker rejects optimized Python execution, which would disable
its assertions, and verifies fixed reviewed SHA-256 pins for this note and
the imported adjacent-powers and dyadic proof sources before checking the
mathematics. A source mismatch requires a reviewed update of the pins.

A separately implemented
[resultant audit](verify_sfc_cubic_resultant_independent.py) uses
`python-flint` and FLINT's exact integer resultant algorithm. It checks
34 model pairs (`n=1,...,32,39,52`), 48 coefficient perturbations
(three variants at each `n=1,...,16`), and 12 raw weighted-kernel pairs
(`n=1,...,12`), as well as the modulo-three counterexample. Its
[receipt](sfc_cubic_resultant_independent_verification.json) records the
dependency version, source hashes, and complete finite cases. This is an
independent implementation and internal replay, not external refereeing
or a proof-assistant verification of the unrestricted theorem.

Replay both checkers from the repository root. The main checker requires only
the standard library. The independent audit requires `python-flint==0.9.0`
(FLINT `3.6.0` in this replay), available to the selected Python interpreter
either as an installed package or through `PYTHONPATH`:

```sh
python3 -B research/factorial/image-rigidity/verify_sfc_cubic_unit_deformation.py \
  --repo-root . \
  --output research/factorial/image-rigidity/sfc_cubic_unit_deformation_verification.json
python3 -B research/factorial/image-rigidity/verify_sfc_cubic_resultant_independent.py \
  --repo-root . \
  --output research/factorial/image-rigidity/sfc_cubic_resultant_independent_verification.json
```

The remaining unrestricted problem is not removed: arbitrary support
moments need not lie in this model's coefficient neighborhood, and arbitrary
linear-form coefficients need not satisfy (12). The modulo-25 counterexample
shows why geometric unit residues alone cannot provide an all-prime replacement.
No generic rank, complex positivity, Riemann hypothesis, external refereeing,
or formal-verification premise is used or asserted.
