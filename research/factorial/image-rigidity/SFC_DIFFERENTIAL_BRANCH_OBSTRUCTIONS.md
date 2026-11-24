# Differential scalarization: exact obstructions and the branch information still needed

Status: proved operator obstructions and exact finite controls, 2026-11-21.
This note does not establish a new SFC closure or an SFC counterexample.
The computations below concern actual ordinary factorial frames, not
arbitrary observables substituted into their function fields.

## 1. Verdict

The proposed universal route

> find a homogeneous scalar equation of Euler width at most the support
> size, with no sufficiently late positive-integer indicial roots,

is false in its simplest precise form. The retained positive-exponent
binary frame already supplies an obstruction at index 3. A second exact
specialization of the same frame supplies an obstruction at index 4 and
also disproves the stronger-looking repair that the narrow scalar
operator should have only constants in its polynomial kernel.

This is an obstruction to a statement about the whole scalar solution
space. It is not a factorial counterexample. The actual factorial series
is one distinguished algebraic branch. Linear differential annihilation
necessarily retains other conjugate branches and their constant linear
combinations; some of those combinations have forbidden zero windows.

The useful next differential target therefore needs either a genuine
exceptional-parameter analysis or an algebraic condition selecting the
distinguished branch. Merely increasing the differential order or
removing apparent singularities does not provide that condition.

## 2. Actual series and the exact propagation lemma

For positive integers `r_i` and distinct nonzero complex coefficients
`c_i`, set

\[
 P(z)=\prod_i(1-c_i z),\qquad
 A(z)=\prod_i(1-c_i z)^{r_i},\qquad a(z)=zA(z),
\]
\[
 H(0)=0,\quad H'=A/P,\qquad
 D=P+zP A'/A,\qquad a'=H'D.
\]

The actual normalized factorial series is

\[
 \rho(t)=\sum_{m\ge0}g_m t^m=\frac1{D(a^{-1}(t))},\qquad
 g_m=\frac{\Lambda(F^m)}{m!\prod_i(mr_i)!},
 \quad F=\Bigl(\prod_iX_i^{r_i}\Bigr)\sum_i c_iX_i.
\]

These identities and their mixed-zero-exponent extension are proved in
[the contact note](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md)
and [the uniform-window note](SFC_ALGEBRAIC_UNIFORM_WINDOWS.md).

Let `theta=t d/dt` and suppose

\[
 L=\sum_{j=0}^{w}t^jP_j(\theta),\qquad L\rho=0.
\]

The coefficient relation is

\[
 \sum_{j=0}^{\min(w,m)} P_j(m-j)g_{m-j}=0.
 \tag{1}
\]

For width `w<=s`, a zero block `g_n,...,g_(n+s-1)` propagates forward if
`P_0(m)` never vanishes for `m>=n+s`. Thus the originally proposed
uniform sufficient head condition was `P_0(k)!=0` for `k>=s+1`.
Requiring it for every positive integer was unnecessarily strong even
before the counterexamples below.

**Truncation lemma.** If any series `y` annihilated by `L` has
`y_n=...=y_(n+s-1)=0` and `w<=s`, then

\[
 Q(t)=\sum_{m<n}y_m t^m\quad\hbox{satisfies}\quad LQ=0.
 \tag{2}
\]

Indeed `ord L(y-Q)>=n+s`, whereas `deg LQ<=n+s-1`. Since their sum is
zero, both vanish. If the first nonzero coefficient after the gap is at
index `k`, then also `P_0(k)=0`.

There is a valid conditional repair: if an actual frame has a width-`s`
annihilator whose polynomial kernel contains only constants, the
[initial resultant theorem](SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md)
proves its SFC conclusion. Equation (2) would force the actual truncation
to be `1`, and the alleged block would then give `g_1=...=g_s=0`, contrary
to that theorem. Section 4 shows that this operator hypothesis itself
does not hold universally.

## 3. Why the retained residue countercontrol constrains every annihilator

Write `u` for the coefficient parameter, reserving `t` for the series
variable. Consider the actual profile

\[
 r=(1,2),\quad c=(1,u),\quad
 a=z(1-z)(1-uz)^2,\quad H'=1-uz,
\]
\[
 D=1-(2+3u)z+4uz^2,\qquad
 \omega=\frac{H'}a\,dz=\frac{dz}{z(1-z)(1-uz)}.
\]

For a zero `alpha` of `a`, let `C_alpha` enclose that zero and no other
one. The residue-generating series is exactly

\[
 \sum_{m\ge0}\operatorname{Res}_{\alpha}
       \frac{H'}{a^{m+1}}\,dz\ t^m
 =\frac1{2\pi i}\int_{C_\alpha}\frac{H'(z)}{a(z)-t}\,dz
 =\sum_{z_i(t)\to\alpha}\frac{H'(z_i(t))}{a'(z_i(t))}
 =\sum_{z_i(t)\to\alpha}\frac1{D(z_i(t))}.
 \tag{3}
\]

The equality first holds for small nonzero `t` and then as a convergent
power series at zero. At `alpha=0,1`, the fibers are simple when `u!=1`,
so these two sums are individual analytic conjugate branches.

Every rational differential annihilator of the branch at zero also
annihilates every algebraic conjugate: the unique extension of
`d/dt` to a separable algebraic extension commutes with conjugation.
It therefore annihilates every constant linear combination of (3).
This fact does not depend on using a minimal-order operator.

Put `v_m=(x_m,y_m)` for the residues at 0 and 1. Exact formulas valid
for every nonnegative integer `m` are

\[
 x_m=\sum_{k=0}^{m}\binom{2m-k}{m-k}\binom{2m+k}{k}u^k,
\]
\[
 y_m=\sum_{k=0}^{m}
 \frac{(-1)^{k+1}\binom{2m-k}{m-k}\binom{2m+k}{k}u^k}
 {(1-u)^{2m+1+k}}.
 \tag{4}
\]

The second formula is the coefficient of `(z-1)^m` in the local
binomial expansion, including the factor `(-1)^(m+1)` from `1-z`.
For orders 0 through 3, direct repeated application of
`(1/a')d/dz` to `1/D`, evaluated at 0 and 1 and divided by `m!`, gives
the same residues. This independently checks the series/conjugate
identification on the finite obstruction.

### 3.1. A width-two equation must have the resonance 3

Take

\[
 5u^2-2u+2=0,\qquad u=(1+3i)/5,
 \qquad \ell_m=y_1x_m-x_1y_m.
\]

The retained exact determinant identities give

\[
 \ell_1=\ell_2=0,\qquad
 \ell_0=\frac{158}{625}-\frac{14u}{125}\ne0,
\]
\[
 \frac{\ell_3}{\ell_0}
 =\frac{-3984+9360u}{625}\ne0.
 \tag{5}
\]

Consequently every width-at-most-two rational differential annihilator
of the actual `rho` annihilates a series with the gap in (5). It must
annihilate the constant function and satisfy `P_0(3)=0`.

The actual factorial coefficient `g_1=2+3u` is nonzero. Thus (5) is
specifically a counterexample to universal scalar-space propagation,
not to SFC.

### 3.2. The order/width/resonance tradeoff occurs exactly

An exact function-field calculation for the same actual frame gives a
minimal operator of differential order 3 and Euler width 3. Its Euler
polynomials are

\[
\begin{aligned}
 P_0(k)&=(11/9600+39i/12800)k(k-1)(k-1/2),\\
 P_1(k)&=(3871/192000-553i/8000)
                 k(k^2-71k/79+37/79),\\
 P_2(k)&=(-919/2400+919i/3200)
                 k(k^2+560k/919+129/919),\\
 P_3(k)&=(k+1/2)(k+3/4)(k+5/4).
\end{aligned}
 \tag{6}
\]

The common factor removed from the raw cofactor operator has degree 13
and origin order 3. The subsequent left-multiple cancellation below
removes the monic factor `t+(u-1)/12`, whose origin order is zero. Thus
neither the final width nor its head is inferred from uncancelled
cofactor degrees.

In ordinary derivative form, its zeroth coefficient is `p_0=15t/32`.
The left multiple `(d/dt-1/t)L`, after polynomial denominator clearing
and exact common-factor removal, has order 4 and **width 2**. One
normalization gives

\[
\begin{aligned}
 \widetilde P_0(k)&=(21/4096-9i/512)k(k-3)(k-1)(k-1/2),\\
 \widetilde P_1(k)&=(-759/5120+2277i/20480)
       k(k^3-k^2-247k/253-68/253),\\
 \widetilde P_2(k)&=(15/32)k(k+1/2)(k+3/4)(k+5/4).
\end{aligned}
 \tag{7}
\]

Thus increased order really attains support-size width here, but the
head acquires the necessary root 3. This is not just a failure to find a
better scalarization. Equation (5) proves that every scalarization of
that width has the resonance.

For this particular operator the polynomial kernel is exactly the
constants: on a polynomial of positive degree `d`, the coefficient of
`t^(d+2)` is its leading coefficient times
`(15/32)d(d+1/2)(d+3/4)(d+5/4)`, which is nonzero. The valid conditional
repair from Section 2 therefore applies to this one specialization.
The next specialization prevents extending that repair universally.

## 4. A second actual specialization forces a nonconstant polynomial kernel

Now choose either root of

\[
 12u^2-5u+5=0
 \quad\bigl(u=(5\mathbin\pm i\sqrt{215})/24\bigr),
 \qquad \ell_m=y_2x_m-x_2y_m.
\]

The canonical determinant formulas immediately give

\[
 \ell_2=\ell_3=0,
 \qquad
 \ell_1=-\frac{5u^3(u-2)}{4(u-1)^6}\ne0,
 \tag{8}
\]

because `5u^2-2u+2=(u-1)/12` on this quadratic. Formula (4) gives

\[
 \ell_0=\frac{953743u}{995328}+\frac{14108185}{11943936}\ne0.
\]

After normalization by `ell_0`, the first five coefficients are exactly

\[
 \left(1,\ \frac{300u}{707}+\frac{25}{101},\ 0,\ 0,
       \frac{2253875u}{203616}-\frac{5000125}{349056}\right).
 \tag{9}
\]

The degree-one and degree-four entries are nonzero since `u` is
nonreal. For an additional direct check, the fourth residue pair is

\[
 x_4=5(99u^4+165u^3+135u^2+63u+14),
\]
\[
 y_4=\frac{5(476u^4-680u^3+408u^2-119u+14)}{(u-1)^{13}}.
\]

By (2), **every** width-at-most-two rational annihilator of the actual
factorial branch at this parameter must annihilate the nonconstant
linear polynomial

\[
 1+\left(\frac{300u}{707}+\frac{25}{101}\right)t.
 \tag{10}
\]

It also satisfies `P_0(4)=0`. Thus allowing the resonance at 3, or
requiring the polynomial kernel to consist only of constants, does not
repair the universal scalar-operator proposal. Again the actual
factorial branch is not the period (9); no factorial zero window is
claimed.

The exact operator construction realizes this obstruction too. Its
minimal order-three operator sends the polynomial (10) to
`(4725/3232+2025u/808)t^2`. The left multiple
`(d/dt-2/t)L`, after denominator clearing and common-factor removal,
has order 4 and width 2. With monic highest ordinary-derivative
coefficient its Euler head and highest Euler coefficient are

\[
 P_0(k)=\frac{7-19u}{320}k(k-1)(k-1/2)(k-4),
\]
\[
 P_2(k)=(k-1)(k+1/2)(k+3/4)(k+5/4).
\]

Thus any nonzero polynomial solution has degree 1; constants are not solutions,
and the polynomial kernel is exactly the span of (10). The raw scalar
gcd again has degree 13 and origin order 3. The final cancellation
removes `t+3(u-1)/26`, of origin order zero.

## 5. Primitive and mixed-zero exact operator controls

The same retained exact connection/cofactor algorithm was run on two
other actual frames. It checks a nonzero first-`q` derivative minor and
an all-coordinate relation for derivative `q`, so the stated orders
are minimal over the rational function field. It also directly checks
normalized factorial coefficients through order 8.

| Actual profile and coefficients | Minimal order | Reduced Euler width | Integer roots of its Euler head |
|---|---:|---:|---|
| `r=(1,1,2), c=(1,2,-2)` | 4 | 6 | `0,1,2` |
| `r=(0,1), c=(1,2)` | 2 | 2 | `0,1` |

### 5.1. The primitive positive profile

Here

\[
 a=8z^5-4z^4-6z^3+z^2+z,\qquad D=20z^3-18z^2+1.
\]

The raw cofactor operator has coefficient degrees `(36,37,38,39,40)`.
Its exact common factor is a nonzero rational constant times

\[
 t^6\Delta_0(t)^9,\qquad
 \Delta_0=12500t^3+7056t^2-1107t-216.
\]

Removing it changes the origin valuation by 6. The reduced coefficient
degrees are `(3,4,5,6,7)` and their origin orders are `(0,0,0,0,1)`.
Consequently multiplying the reduced ordinary operator by `t^3`
gives Euler width 6. Its head is

\[
 -146124k(k-2)(k-1)(2k-1).
\]

The leading ordinary-derivative coefficient factors as

\[
 t\Delta_0(t)\Psi(t),\qquad
 \Psi=3906t^3-19740t^2+3092t+1353.
\]

The factor `Psi` is apparent. Adjoining the constant solution by
`(d/dt-p_0'/p_0)L`, followed by exact scalar-gcd removal, increases
order to 5 and removes `Psi`, but inserts

\[
 \Phi=273420t^3-2443231t^2+840714t+669834
\]

in its place. The width remains **6**, and the head acquires the extra
integer root 3. Removing an apparent factor therefore did not reduce
the width in this control. No claim is made that width 3 is impossible
for this primitive coefficient specialization.

For reproducibility, the primitive order-four operator's Euler
coefficients, in addition to the head above, are

\[
\begin{aligned}
 P_1&=-9k(k-1)(240627k^2-93309k+279494),\\
 P_2&=12k(865647k^3+1115451k^2+6635450k-884648),\\
 P_3&=12(4978178k^4+11225837k^3+27819389k^2+25762488k+8038008),\\
 P_4&=-2(52479691k^4-17971833k^3-226267222k^2-228092208k-60531408),\\
 P_5&=-168(1304698k^4+7261724k^3+12576851k^2+8709829k+2094198),\\
 P_6&=78120(5k+3)(5k+4)(5k+6)(5k+7).
\end{aligned}
\]

### 5.2. A mixed-zero frame where the narrow recurrence works

For `r=(0,1),c=(1,2)`,

\[
 a=z-2z^2,\qquad D=1-5z+4z^2.
\]

The exact common scalar factor of the raw operator is
`(t+1)^2(8t-1)^2/64`, of origin order zero. The reduced equation is

\[
 [12+2(14t+5)\partial_t+(t+1)(8t-1)\partial_t^2]\rho=0.
\]

Its Euler coefficients are

\[
 P_0(k)=-k(k-1),\quad
 P_1(k)=k(7k+3),\quad
 P_2(k)=4(k+1)(2k+3).
\]

There are no head resonances at integers at least 3, so this control
has the desired propagation property. The obstruction is to a
universal theorem, not to using scalar recurrences on individual
structural families.

## 6. Why a rational inhomogeneous source cannot select the branch here

For positive exponents, `H` is polynomial with `deg H<deg a=d`. The
trace of `H` over all inverse branches is constant, by the Newton
identities for `a(z)-t`. Therefore

\[
 \operatorname{Tr}\rho
 =\partial_t\operatorname{Tr}H(z_i(t))=0.
\]

If an operator `L in C(t)[d/dt]` satisfied `Lrho=R(t)` with rational
`R`, taking the trace would give

\[
 dR=L(\operatorname{Tr}\rho)=0.
\]

Thus `R=0`. A nonzero rational inhomogeneous right-hand side cannot
repair the branch-selection failure in this class. This argument does
not exclude nonrational auxiliary data, nonlinear equations, or
retaining the marked inverse branch in a larger algebraic system.

## 7. What the primary differential-algebra literature actually supplies

* Bostan--Chyzak--Lecerf--Salvy--Schost,
  [Differential Equations for Algebraic Functions](https://specfun.inria.fr/chyzak/Publications/BostanChyzakLecerfSalvySchost-2007-DEA.pdf),
  ISSAC 2007, Theorems 1--3. For a separable algebraic equation in
  characteristic zero they bound differential order, coefficient
  degrees, and recurrence size, including improved bounds after
  increasing order. Their equations annihilate all conjugate roots.
  The bounds depend on the algebraic equation's degrees; they do not
  supply support-size width with the required indicial control.

* Chen--Kauers--Singer,
  [Desingularization of Ore Operators](https://singer.math.ncsu.edu/papers/desing.pdf),
  Section 2 and Theorem 6. Their construction removes removable
  factors using higher-order left multiples. The solution space is
  retained or enlarged. The introduction explicitly distinguishes
  removal with a supplied order budget from proving a bound on the
  necessary budget. Thus it cannot remove the conjugate-period
  obstruction above.

* Chen--Jaroschek--Kauers--Singer,
  [Desingularization Explains Order-Degree Curves for Ore Operators](https://arxiv.org/pdf/1301.0917),
  Section 2 and Theorem 9. Order can be traded for coefficient degree,
  but removing a leading factor alone does not control all other
  coefficients. This is precisely why recurrence width and indicial
  roots must be audited after every cancellation, as in Sections 3
  and 5.

* Van der Hoeven,
  [Creative telescoping using reductions](https://www.texmacs.org/joris/redtel/redtel.html),
  introduction and reduction framework. This supplies a practical
  way to construct differential/recurrence annihilators for the
  relevant parameterized integrals. It does not identify one marked
  cycle from the full period space or imply the sharp zero-window
  theorem.

These are algorithmic and structural inputs. None is cited as proving
the missing SFC assertion.

## 8. Ranked next targets, with precise failure tests

**1. A branch-sensitive Padé or jet certificate.** Keep the actual
algebraic constraint

\[
 a(z(t))=t,\quad z(0)=0,\qquad D(z(t))\rho(t)=1
 \tag{11}
\]

through the elimination. The sharp missing lemma is a bound
`ord(H-b(a))<=deg b+s` for this special root profile, or an equivalent
saturated jet-ideal certificate. A computational first target is an
identity in the jet ideal of (11), using the nonzero Jacobian
`a'(0)=1`, that excludes a block of `s` actual moment coefficients.
The exact initial resultant gives the base certificate. Any candidate
general identity must use (11) in an essential way: it must reject
the spurious period (9), whose constant term is also normalized to 1.
Replacing (11) by a scalar annihilator loses the necessary information.

**2. A generic narrow-operator theorem with an explicit exceptional
divisor and a separate certificate on that divisor.** This is weaker
than the refuted universal constants-only-kernel theorem and could
still be useful. The exact target is a parameter-stratified theorem
giving width, polynomial kernel, and resonance set simultaneously.
It must then prove actual moment regularity on the exceptional
strata, rather than delete them. Equations (5) and (9) are mandatory
tests. A generic nonvanishing computation alone cannot settle all
complex coefficients.

**3. Improved uniform windows, rather than claiming sharp closure.**
The exact connection and cofactor machinery already proves finite
windows in the canonical note. Better choices of a cyclic observable,
controlled Ore left multiples, and exact apparent-factor removal may
improve those bounds. The primitive profile `(1,1,2)` is a concrete
benchmark: its minimal actual operator above has width 6, while its
requested SFC window is 3. A rigorous reduction to a smaller finite
list of shifted windows would be useful even if the sharp universal
bound remains out of reach.

**Stop conditions.** Discard any universal proof that only uses linear
span of conjugate periods, support-sized module rank, or cancellation
of apparent factors without auditing the resulting width and
polynomial kernel. Sections 3 and 4 give exact counterexamples to
those stronger hypotheses. Do not interpret failure of these
operator hypotheses as evidence against the actual SFC statement.

## 9. Finite audit provenance

Replay from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_differential_branch_obstructions.py
```

The standard-library checker pins this note, the retained contact note,
and both imported exact helper sources. It writes
`sfc_differential_branch_obstructions_verification.json`. Its complete scope is
two exact quadratic-field cases, 20 residue/conjugate coefficient
comparisons, 18 direct factorial coefficient comparisons, two minimal
cofactor operators, and two exact narrow left multiples with polynomial
kernel witnesses. It rejects optimized execution, changed note bytes,
and missing or changed sources.

The rational controls use the exact function-field construction in
`verify_sfc_uniform_rational_window.py`, whose `check_case` checks the
actual factorial expansion and a complete derivative-vector relation.
The quadratic-field control uses the same arithmetic over exact
Gaussian rationals, not floating point. The new residue control (9)
uses (4) and exact reduction modulo `12u^2-5u+5`.

The optional rational controls in Section 5 give their polynomial
operators explicitly. They use the same retained `check_case` construction
on `r=(1,1,2),c=(1,2,-2)` and `r=(0,1),c=(1,2)`; they are not additional
cases counted in the canonical two-field replay. These finite controls
are not exhaustive parameter checks. All-order claims in this note are
the explicitly proved propagation, conjugation, truncation, and trace
arguments.
