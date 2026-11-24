# Ordinary binomial Strong Factorial for adjacent positive powers of three

Internal proof, 21 November 2026. Base repository revision
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## 1. Statement and scope

**Theorem.** Let `K` be any field of characteristic zero, let `r>=1`,
and put `a=3^r`. With the ordinary factorial functional
`Lambda(X^j)=j!`, every nonzero polynomial

\[
 f(X)=X^a(A+BX^{2a}),\qquad A,B\in K,
\]

satisfies

\[
 \bigl(\Lambda(f^n),\Lambda(f^{n+1})\bigr)\ne(0,0)
 \qquad\text{for every integer }n\ge1.
 \tag{1}
\]

Thus the ordinary univariate binomial supports
`{3,9}`, `{9,27}`, `{27,81}`, and all subsequent adjacent positive
powers of three satisfy consecutive-moment nonvanishing at every
terminal order `N=n+1>=2`. A coefficient is allowed to vanish; only
the zero polynomial is excluded.

This proves a family of ordinary **univariate binomial** SFC(2)
statements. It does not prove unrestricted Strong Factorial, the
multivariate two-term conjecture, or arbitrary Gamma-shape variants.
No ordering or positivity assumption on `K` is made. The proof uses
neither analytic positivity nor a Riemann-hypothesis input.

The main mechanism is an exact congruence modulo nine, shared by the
whole family. Its model has many common roots modulo three, so plain
residue coprimality fails. We prove that every such common root is
simple in the lower polynomial. A hypothetical characteristic-zero
common root must consequently be unramified at three. A cubic trace
identity then contradicts the congruence modulo nine.

## 2. Actual moment normalization

For fixed `a=3^r`, define

\[
 T_n(t)=\sum_{k=0}^n\binom nk(an+1)_{2ak}t^k,
 \qquad H_n(u)=u^nT_n(u^{-1}),
 \tag{2}
\]

where `(b)_m=b(b+1)\cdots(b+m-1)` and `(b)_0=1`.
The reciprocal in (2) is the **positive** reciprocal, with no sign
change. Put

\[
 P_n(x)=3^{-an}H_n(3^a x)
       =\sum_{k=0}^n b_{n,k}x^{n-k},\qquad
 b_{n,k}=3^{-ak}\binom nk(an+1)_{2ak}.
 \tag{3}
\]

The same scale `u=3^a x` is used at both consecutive orders. For
`A,B` nonzero,

\[
 \Lambda(f^m)=A^m(am)!T_m(B/A),\qquad
 u=A/B,\quad x=\frac{A}{3^aB}.
 \tag{4}
\]

Thus a consecutive zero pair would give a common root of `P_n` and
`P_(n+1)`. The cases `A=0` or `B=0` are monomials: each positive
factorial moment is a nonzero coefficient power times a factorial,
and is nonzero in characteristic zero.

### Lemma 1. The common integral model modulo nine

For `m>=0`, set

\[
 U(m)=\prod_{j=1}^m(3j-2)(3j-1),\qquad U(0)=1.
\]

The actual coefficients in (3) satisfy the exact rational identity

\[
 b_{n,k}=\binom{n+2k}{3k}\,U(k)
    \prod_{j=0}^{r-1}
      \frac{U(3^j(n+2k))}{U(3^j n)}.
 \tag{5}
\]

In particular, `P_n` is a monic integer polynomial, its constant
coefficient is `U(an)`, a three-adic unit, and

\[
 P_n(x)\equiv Q_n(x)\pmod9,\qquad
 Q_n(x)=\sum_{k=0}^n(-1)^k\binom{n+2k}{3k}x^{n-k}.
 \tag{6}
\]

**Proof.** Removing the multiples of three gives

\[
 (3m)!=3^m m!U(m),
\]

and iterating this identity gives

\[
 (3^r m)!=3^{m(3^r-1)/2}m!
             \prod_{j=0}^{r-1}U(3^j m).
 \tag{7}
\]

Apply (7) to the factorial quotient
`(a(n+2k))!/(an)!` in (3), then use
`(3k)!=3^k k!U(k)`. The power of three left over is exactly zero,
and the remaining factors are (5).

Every denominator in (5) is a three-adic unit, so `b_(n,k)` is
three-adically integral. Its definition in (3) puts it in
`Z[1/3]`; these two facts give `b_(n,k) in Z`. At `k=0` the
coefficient is one. At `k=n` the unit product telescopes:

\[
 U(n)\prod_{j=0}^{r-1}\frac{U(3^{j+1}n)}{U(3^j n)}=U(an).
\]

Finally,

\[
 (3j-2)(3j-1)=9j(j-1)+2\equiv2\pmod9,
 \qquad U(m)\equiv2^m\pmod9.
\]

All denominators remain invertible modulo nine. The unit multiplier
in (5) is congruent to

\[
 2^{k+2k\sum_{j=0}^{r-1}3^j}=2^{ak}\equiv(-1)^k\pmod9,
\]

because `a=3^r` is congruent to three modulo six for `r>=1`.
This proves (6). \(\square\)

Every root of `P_n` in an algebraic closure of `Q_3` is a unit.
Indeed, monicity and integral coefficients give nonnegative root
valuations; the product of the roots is the unit constant
coefficient, so all valuations are zero. Equivalently, a root of
positive valuation would leave the constant term as the unique
term of least valuation, and a root of negative valuation would
leave the leading term as the unique term of least valuation.

For `r=1`, (5) also gives the exact valuation identity

\[
 v_3\!\left(\binom nk(3n+1)_{6k}\right)-3k
   =v_3\binom{n+2k}{3k}.
\]

More generally, the left subtraction is `ak` for `a=3^r`. Thus
every root of the positive reciprocal `H_n` has valuation exactly
`a`; the passage to `P_n` removes that common valuation.

## 3. Generating function and residual simplicity

The model in (6) has the exact generating function

\[
 G(t,x)=\sum_{n\ge0}Q_n(x)t^n
       =\frac{(1-xt)^2}{(1-xt)^3+t}.
 \tag{8}
\]

To prove (8), put `j=n-k` and sum in the formal power-series ring:

\[
 \sum_{j,k\ge0}(-1)^k\binom{j+3k}{3k}x^jt^{j+k}
 =\frac1{1-xt}\sum_{k\ge0}
        \left(\frac{-t}{(1-xt)^3}\right)^k.
\]

Each coefficient receives only finitely many contributions, so this
rearrangement is formal and requires no convergence hypothesis.
The initial polynomials and exact recurrence are

\[
 Q_0=1,\quad Q_1=x-1,\quad Q_2=x^2-4x+1,
\]
\[
 Q_n=(3x-1)Q_{n-1}-3x^2Q_{n-2}+x^3Q_{n-3}
 \quad(n\ge3).
 \tag{9}
\]

Let `R_n` be the reduction of `Q_n` modulo three. Define
`A_m(y) in F_3[y]` by

\[
 \sum_{m\ge0}A_m(y)t^m=\frac1{1+t-yt^3},\qquad
 A_0=1,\quad A_{-1}=A_{-2}=0.
\]

For `m>=1`,

\[
 A_m=-A_{m-1}+yA_{m-3},
 \tag{10}
\]

with the given negative-index values. Reducing (8) gives

\[
 R_n(x)=A_n(x^3)+xA_{n-1}(x^3)+x^2A_{n-2}(x^3),
\]
\[
 R_n'(x)=A_{n-1}(x^3)-xA_{n-2}(x^3).
 \tag{11}
\]

The derivative of each `A_m(x^3)` is zero in characteristic three.

### Lemma 2. Every common residual root is simple

For each `n>=1`, if `xi` in an algebraic closure of `F_3` satisfies
`R_n(xi)=R_(n+1)(xi)=0`, then `xi!=0` and `R_n'(xi)!=0`.

**Proof.** Formula (6) gives `R_n(0)=(-1)^n`, so `xi!=0`.
Using (10) in (11) gives the polynomial identity in `F_3[x]`

\[
 x^2 A_{n-2}(x^3)
   =(x-1)R_n(x)+xR_n'(x)-R_{n+1}(x).
 \tag{12}
\]

For clarity, if `B=A_n(x^3)`, `C=A_(n-1)(x^3)`, and
`a_0=A_(n-2)(x^3)`, then

\[
 R_n=B+xC+x^2a_0,\qquad R_n'=C-xa_0,
\]
\[
 R_{n+1}=(x-1)B+x^2C+x^3a_0,
\]

which verifies (12) directly.

Suppose also `R_n'(xi)=0`. Identity (12) gives
`A_(n-2)(xi^3)=0`, and then (11) gives
`A_(n-1)(xi^3)=A_n(xi^3)=0`. For `n=1` this already contradicts
`R_1'=1`. For `n>=2`, the recurrence (10) can be run backwards
because `xi^3!=0`. Three consecutive zero values force all preceding
values to be zero, including `A_0=1`, a contradiction.
\(\square\)

The lemma concerns common roots, not the squarefreeness of each
entire `R_n`. It is stronger than the finite observation that sampled
adjacent gcds are squarefree. Direct coprimality modulo three is
false: `R_3` and `R_4` have gcd `x-1`. At larger indices the common
gcd can have higher degree; no bounded root-count argument is used.

## 4. The indispensable unramified step

### Lemma 3. An actual common root is unramified

If `P_n` and `P_(n+1)` have a common root `alpha` over an algebraic
closure of `Q_3`, then `alpha` belongs to a finite unramified
extension of `Q_3` and is a unit in its valuation ring.

**Proof.** Lemma 1 gives `v_3(alpha)=0`. Its residue `bar alpha`
belongs to a finite field `F_(3^d)`. By (6) and Lemma 2,
`bar alpha` is a simple root of the reduction of `P_n`.

Take the finite unramified extension `E/Q_3` with residue field
`F_(3^d)` and choose a lift of `bar alpha` in its valuation ring
`O_E`. The derivative of `P_n` at that lift is a unit. Newton
iteration

\[
 b_{j+1}=b_j-\frac{P_n(b_j)}{P_n'(b_j)}
\]

stays in `O_E`, keeps the derivative a unit, and at least doubles the
valuation of the error at each step. Completeness therefore gives
a root `beta in O_E` of `P_n` with residue `bar alpha`.

This root is unique even among roots in ramified extensions with
the same residue. Indeed, if `delta=alpha-beta`, then
`v_3(delta)>0` and polynomial division gives

\[
 P_n(alpha)-P_n(beta)
   =(alpha-beta)\bigl(P_n'(beta)+(alpha-beta)J\bigr),
\]

where `J` is integral because `alpha`, `beta`, and the coefficients
of `P_n` are integral. The expression in parentheses is a unit.
Both polynomial values are zero, so `alpha=beta in E`.
\(\square\)

The unramified conclusion cannot be omitted. In `O_E`, reduction
modulo three is a field and every valuation is an integer. Thus
`z^3 in 3O_E` implies `z in 3O_E`. For a general ramified extension
this implication fails: an element satisfying `z^3=3` has valuation
`1/3`. Merely knowing that the actual root is integral would not
justify the next argument.

## 5. Cubic trace identity and contradiction modulo nine

Let `tau_1,tau_2,tau_3` be the roots, counted with multiplicity, of

\[
 \tau^3+\tau^2-x=0.
\]

For `m>=0` define symmetric power sums

\[
 S_m(x)=\sum_{i=1}^3\tau_i^m,\qquad
 V_m(x)=\sum_{1\le i<j\le3}(\tau_i\tau_j)^m.
 \tag{13}
\]

These are integer polynomials, not chosen algebraic branches. For
example they can be defined independently by

\[
 S_0=3,\ S_1=-1,\ S_2=1,\qquad
 S_m=-S_{m-1}+xS_{m-3}\quad(m\ge3),
\]
\[
 V_0=3,\ V_1=0,\ V_2=2x,\qquad
 V_m=xV_{m-2}+x^2V_{m-3}\quad(m\ge3).
 \tag{14}
\]

The second recurrence follows because the three pairwise products
have polynomial `v^3-xv-x^2`.

### Lemma 4. The two exact trace formulas

For every `N>=1`,

\[
 Q_N(x)+2xQ_{N-1}(x)=S_{3N}(x),
 \tag{15}
\]
\[
 S_{3N}(x)=S_N(x)^3-3S_N(x)V_N(x)+3x^N.
 \tag{16}
\]

**Proof.** The elementary symmetric functions of the `tau_i` are
`-1,0,x`. Consequently the sum, pair sum, and product of their
cubes are `3x-1,3x^2,x^3`. Hence

\[
 D(t,x):=\prod_{i=1}^3(1-\tau_i^3t)=(1-xt)^3+t.
\]

Taking the formal logarithmic derivative and using (8) gives

\[
 \sum_{N\ge1}S_{3N}(x)t^N
   =-\frac{t\partial_tD}{D}
   =(1+2xt)G(t,x)-1.
 \tag{17}
\]

The last equality is the polynomial identity
`(1+2xt)(1-xt)^2-D=-t partial_t D`, so involves no division by
anything with a nonunit constant coefficient. Extracting positive
coefficients proves (15). The identity is deliberately stated for
`N>=1`; the constant terms of `S_0` and `Q_0` differ by two.

Apply the elementary identity

\[
 y_1^3+y_2^3+y_3^3
 =(y_1+y_2+y_3)^3
  -3(y_1+y_2+y_3)(y_1y_2+y_1y_3+y_2y_3)
  +3y_1y_2y_3
\]

to `y_i=tau_i^N`. Since their product is `x^N`, this proves (16).
Both formulas are identities of integer polynomials, so remain
valid after every specialization, including repeated cubic roots.
\(\square\)

**Proof of the theorem over `Q`.** Suppose `P_n` and `P_(n+1)` had a
common root `alpha` over an algebraic closure of `Q_3`, and put
`N=n+1`. Lemma 3 puts `alpha` in the valuation ring of a finite
unramified extension `E/Q_3`, with `alpha` a unit. By (6),

\[
 Q_n(alpha),\ Q_N(alpha)\in9O_E.
\]

Formula (15) therefore gives `S_(3N)(alpha) in 9O_E`. Reducing
(16) modulo three, in the field `O_E/3O_E`, gives

\[
 S_N(alpha)^3=0\pmod3,
\]

so `S_N(alpha) in 3O_E`. As `V_N(alpha)` is integral, (16) now gives

\[
 S_{3N}(alpha)\equiv3alpha^N\pmod {9O_E}.
\]

The right side has valuation one because `alpha` is a unit. It
cannot belong to `9O_E`, which is the desired contradiction. Thus
`P_n` and `P_(n+1)` have no common root over the algebraic closure
of `Q_3`, and their gcd over `Q` is one.

**Transfer to `K`.** Reversal and the common nonzero scaling in
(2)–(3) preserve the existence of a nonzero common root. All
constant coefficients of `T_n` are one. Thus `T_n,T_(n+1)` are
coprime in `Q[t]`. A Bezout identity with rational coefficients
remains an identity over every characteristic-zero field `K`.
Together with (4) and the monomial cases, this proves (1) for all
the stated coefficient fields. \(\square\)

## 6. Exact checks, proof coverage, and limitations

The companion
[checker](verify_stepanov_adjacent_powers_of_three.py) uses Python's
standard library and exact integer polynomial arithmetic. It checks
the universal low-degree polynomial identities used in the proof,
the actual factorial normalization at explicitly recorded finite
indices and exponents, the integral recurrences and trace formulas,
and finite residual-gcd/Hensel-obstruction controls. It retains the
noncoprime reduction at `n=3` as a negative control; the lower lift
in that disk is `x=7 mod9`, where the upper model is `3 mod9`.

The [receipt](stepanov_adjacent_powers_of_three_verification.json)
records the exact finite domains, current source hashes, and the
outcome of an internal replay. The arbitrary-`r`, arbitrary-`n`
theorem follows from the algebraic arguments in Sections 2–5;
neither the finite samples nor the receipt establish those
quantifiers independently.

The proof was reconstructed by the root reviewer and the parallel
moment/determinant reviewer. The latter explicitly checked the
positive reciprocal convention, the factorial unit splitting, the
residual derivative identity, the unramified lifting step, and both
trace identities. These are internal AI mathematical reviews, not
external refereeing or formal proof-assistant verification.

The result closes `{3,9}` without requiring a small bound for the
degree of a residual gcd. The same proof simultaneously covers
every `{3^r,3^(r+1)}` with `r>=1`. No assertion about arbitrary
lower exponents, arbitrary Gamma shapes, or unrestricted Strong
Factorial follows from the congruence alone.
