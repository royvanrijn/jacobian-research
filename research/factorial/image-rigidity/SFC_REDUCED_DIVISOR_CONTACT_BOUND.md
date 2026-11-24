# A reduced logarithmic-divisor contact bound

Internal partial theorem, 21 November 2026. This improves the retained
degree-based contact estimate when coefficient groups have neutral
exponent totals or the critical quotient has repeated roots. It does not
prove the desired sharp `n+S` bound or
unrestricted Strong Factorial. The original contact proof is retained unchanged.

## 1. Exact statement

Work over any characteristic-zero field. Consider the ordinary phase

\[
 F=\left(\prod_{j=1}^t\prod_{l=1}^{p_j}Y_{j,l}^{r_{j,l}}\right)
                        \sum_{j=1}^t c_j\sum_{l=1}^{p_j}Y_{j,l},
\]

where `t>=1`, the `c_j` are distinct and nonzero, `p_j>=1`, and
`r_(j,l)>=0`.
Set `R_j=sum_l r_(j,l)`, `S=sum_j p_j`, and

\[
 P=\prod_j(1-c_jz)^{p_j},\qquad
 A=\prod_j(1-c_jz)^{R_j},\qquad
 a=zA,\qquad W=\frac AP.
 \tag{1}
\]

Let `H(0)=0,H'=W` be the formal primitive and put `B=H composed with
a^{-1}`. Then `a=z+O(z^2)` and `H=z+O(z^2)`. Define

\[
 D=\frac{a'}W=P+zP\frac{A'}A,
 \qquad
 J=\prod_{j:R_j\ne p_j}(1-c_jz),\qquad T=\deg J.
 \tag{2}
\]

The polynomial `D` has degree exactly `S` and constant term one: its
leading coefficient is `(1+sum_j R_j)` times that of `P`. The
polynomial `J` records the distinct finite zeros and poles of the
rational function `W`; its constant term is one. In the positive
polynomial-observable case `R_j>=p_j`, `T` is the number of groups with
strict excess `R_j>p_j`. For distinct coefficients and positive individual
exponents it is precisely the number of exponents greater than one.
Let `K=rad(D)` be the product of the distinct factors of `D`, normalized
by `K(0)=1`, and set `d=deg K`, `G=D/K`. Thus `1<=d<=S` and
`G(0)=1`.
Finally put `L=gcd(J,K)`, normalized by `L(0)=1`, and set
`J_0=J/L`, `T_0=deg J_0`. Thus `T_0` counts the logarithmic-divisor
points not already roots of the critical quotient, and `0<=T_0<=T`.

**Theorem.** For every integer `n>=1` and every polynomial `b` of degree
at most `n`,

\[
 \boxed{\operatorname{ord}_0(H-b\circ a)
       =\operatorname{ord}_0(B-b)
       \le nd+(n-1)T_0+1
       \le nd+(n-1)T+1
       \le nS+(n-1)T+1.}
 \tag{3}
\]

Consequently every consecutive moment window beginning at `n>=1` with
length

\[
                   \boxed{(d+T_0-1)n-T_0+1}
 \tag{4}
\]

contains a nonzero ordinary factorial moment.

Repeated coefficients and zero common exponents are included with their
actual multiplicities. For positive exponents and distinct coefficients, if `u` of the
`S` exponents equal one, (3) improves the retained bound by exactly
`(n-1)u` at the `S`-level, with a further improvement `n(S-d)`
when `D` has repeated roots:

\[
 nS+(n-1)T+1=(2n-1)S+1-(n-1)u.
 \tag{5}
\]

In this positive, distinct-coefficient case `J` and `K` are coprime:
at `z=1/c_j`, the term `zP A'/A` is nonzero. Thus `T_0=T` there.
For grouped coefficients, the additional overlap removal can improve
the bound by a further `(n-1)(T-T_0)`.

This is a support- and order-uniform bound, not an extrapolation from
finitely many starting indices.

## 2. Remove the factors that the logarithmic derivative does not need

The exact logarithmic derivative is

\[
 \frac{W'}W=-\sum_{j:R_j\ne p_j}
                  \frac{(R_j-p_j)c_j}{1-c_jz}.
\]

Therefore

\[
                        E=J\frac{W'}W
 \tag{6}
\]

is polynomial of degree at most `T-1`. If `T=0`, set `J=1,E=0`.
Using the full coefficient polynomial `P` as a denominator introduces
unnecessary factors at every neutral group and at repeated copies of
the same group. Formula (6) retains only the actual logarithmic divisor.

The chain rule gives `B'(a)=1/D`. For integers `k>=2`, define
polynomials by

\[
\begin{aligned}
 N_2&=-D',\\
 N_{k+1}&=JD N_k'
  -\bigl((k-1)ED+(2k-1)JD'+(k-2)J'D\bigr)N_k.
\end{aligned}
 \tag{7}
\]

**Exact derivative identity.** For every `k>=2`,

\[
             B^{(k)}(a(z))
       =\frac{N_k(z)}{W(z)^{k-1}D(z)^{2k-1}J(z)^{k-2}}.
 \tag{8}
\]

For `k=2`, differentiate `B'(a)=1/D` and divide by `a'=WD`.
For the induction step differentiate the right side of (8) with respect
to `z`, divide by `WD`, and use `W'/W=E/J`. Multiplying the bracket by
`JD` gives exactly (7). This proves the identity as a rational-function
identity after composition; all its denominators are units at zero.

Each multiplier in (7) increases degree by at most `S+T-1`, so

\[
                 \deg N_k\le S-1+(k-2)(S+T-1).
 \tag{9}
\]

This includes `T=0`. No `N_k` is identically zero: otherwise (8) would
make `B` a polynomial in characteristic zero. The polynomial identity
`D(z)B'(a(z))=1` would then be impossible because `deg D=S>=1`.

### Remove repeated roots of the critical quotient

Define

\[
 U=\frac{D'}G,\qquad C_2=-U,
\]

\[
 C_{k+1}=JKC_k'
 -\bigl((k-1)EK+kJU+(k-1)JK'+(k-2)J'K\bigr)C_k.
 \tag{9a}
\]

The polynomial `U` has degree exactly `d-1`: in characteristic zero the
polynomial `G=D/rad(D)` divides `D'`. Direct differentiation gives

\[
 B^{(k)}(a(z))=
 \frac{C_k(z)}{W(z)^{k-1}D(z)^kK(z)^{k-1}J(z)^{k-2}}.
 \tag{9b}
\]

Indeed the logarithmic derivatives needed in the induction are
`W'/W=E/J` and `D'/D=U/K`. Differentiate (9b), divide by `WD`, and
multiply its remaining bracket by `JK`; the result is exactly (9a).
No coprimality between `J` and `K` is assumed or needed. Comparison
with (8) proves the full polynomial divisibility identity

\[
                         N_k=G^{k-1}C_k.
 \tag{9c}
\]

Every multiplier in (9a) increases degree by at most `d+T-1`, including
when `T=0`. Consequently

\[
                  \deg C_k\le d-1+(k-2)(d+T-1).
 \tag{9d}
\]

These polynomials cannot vanish identically, by (9c) and the preceding
nonpolynomiality argument. This also covers `T=0`: `W=1`, but
`B'(a)=1/D` still cannot be the composition of a polynomial derivative
with `a` when the polynomial `D` has positive degree.

### Remove a logarithmic pole already counted by the critical quotient

There is one further exact cancellation when `J` and `K` overlap.
With `L,J_0,T_0` defined above, put `C_hat_2=-U` and

\[
\begin{split}
 \widehat C_{k+1}={}&J_0K\widehat C_k'\\
 &-\left((k-1)\frac KL E+kJ_0U
        +(k-1)J_0K'+(k-2)J_0'K\right)\widehat C_k.
\end{split}
 \tag{9e}
\]

Every term is polynomial. The same differentiation gives

\[
 B^{(k)}(a)=
 \frac{\widehat C_k}{W^{k-1}D^kK^{k-1}J_0^{k-2}},
 \qquad C_k=L^{k-2}\widehat C_k.
 \tag{9f}
\]

The only changed logarithmic term is
`K J_0 W'/W=(K/L)E`. Its degree is at most `d+T_0-1`, as are the
other multipliers. Hence the final nonzero numerator satisfies

\[
             \deg\widehat C_k\le d-1+(k-2)(d+T_0-1).
 \tag{9g}
\]

## 3. Contact and moment windows

Composition with the normalized formal automorphism `a` preserves
vanishing order. Let `M=ord_0(B-b)`. If `M<=n`, (3) already holds.
Otherwise differentiating `n+1` times, using `b^(n+1)=0` and (9f), gives

\[
 M-n-1=\operatorname{ord}_0\widehat C_{n+1}
 \le d-1+(n-1)(d+T_0-1).
\]

Rearrangement is (3). The argument also proves that the contact order is
finite, since `B` cannot be polynomial.

The ordinary normalized moments are exactly

\[
 g_m=\frac{\Lambda(F^m)}{m!\prod_{j,l}(mr_{j,l})!}
     =[z^m]A^{-m}P^{-1}
     =(m+1)[t^{m+1}]B(t).
 \tag{10}
\]

If a block `g_n,...,g_(n+ell-1)` vanished, the degree-`n` truncation of
`B` would have contact at least `n+ell+1`. Taking the value of `ell` in
(4) exceeds (3) by one. Every removed factorial is a nonzero integer,
so the same nonvanishing conclusion holds for ordinary moments over
every characteristic-zero field.

At `n=1`, the guaranteed window is `d<=S`, retaining and sometimes
improving the known initial-window result. For the primitive positive profile `(1,1,2)`,
`S=3,T=1`, so the contact bound is `4n` and the window is `3n`, improving
the previous `6n-2` contact and `5n-2` window. The desired contact is
`n+3`; that remains strictly stronger for later starting orders.

More generally, if only `k` distinct coefficient groups have nonneutral
totals, then `T<=k` and (3)--(4) are uniform in the size of their
exponents. No upper bound on those exponents is needed.

## 4. Relation to the previous numerator and to pullback

Let `M_k` be the numerator in the retained derivative formula using
`P` instead of `J`. Comparing its exact rational identity with (8)
gives the polynomial identity

\[
                         M_k=(P/J)^{k-2}N_k.
 \tag{11}
\]

The quotient `P/J` is polynomial: each factor of `J` corresponds to a
coefficient group of positive multiplicity in `P`. Thus the improvement
removes a proven factor; it does not assume generic cancellation or a
nonvanishing determinant.

If the frame first descends by an admissible normalized polynomial
pullback, `a=a_0(h), H=H_0(h)`, then
`H-b(a)=(H_0-b(a_0)) composed with h`. Since `h=z+O(z^2)`, the contact
order is unchanged. The bound can therefore be applied to a descended
frame with its values of `d,T_0`.

## 5. Limits and exact-control scope

The final degree bound (9g) is actually an equality for these factorial
frames. Write `R=sum_j R_j`. At infinity `B'(a)=1/D` has degree
`-S`. Each further operation `(1/a')d/dz` lowers that degree by
`R+1`; its leading derivative multiplier is
`-S-(k-2)(R+1)`, which never vanishes. Thus `B^(k)(a)` has degree
`-S-(k-1)(R+1)`. Comparing with the denominator in (9f) gives exactly
`deg C_hat_k=d-1+(k-2)(d+T_0-1)`. There is no further universal
cancellation of the leading degree in this numerator. A sharper
conclusion would need a bound on its order at zero, or another method.

The proof uses the actual logarithmic divisor but does not yet exploit
that every zero of `H'` is a zero of `a` strongly enough to prove
`n+S`. Indeed the argument works more generally whenever `a` is a
normalized polynomial, `W` is a rational unit at zero,
`D=a'/W` is polynomial and nonconstant, and `J` clears the finite
logarithmic divisor of `W`. The retained counterexample with only
`a'=H'D` therefore remains relevant. This partial theorem does not
repair that failed shortcut.

The companion standard-library checker reconstructs (7), (9a), and (9e), compares
the resulting rational derivatives with independently computed inverse
series, verifies (11), (9c), and the overlap factor in (9f), and checks the claimed degrees on its stated
finite inputs. Those exact controls test the formulas; Sections 1--3
prove the unrestricted bound. No nonsingular scalar recurrence,
positivity of complex moments, or correspondence Fatou theorem is used.

The controls include an actual positive frame with repeated `D` roots
(`r=(1,1,1),c=(1,1,1)`, so `D=(1-z)^2(1-4z)`) and an actual positive
frame in which `J` and `K` overlap (`r=(1,2),c=(1,1)`, so `J=1-z`
and `D=(1-z)(1-4z)`). Neutral `T=0` frames, cyclic coefficient frames,
and formal frames with zero common exponents are tested separately.

The retained proof being sharpened is
[SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md](SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md),
Sections 1--3. The exact pullback relation used in Section 4 is in
[SFC_RAMIFIED_PULLBACK_CLOSURE.md](SFC_RAMIFIED_PULLBACK_CLOSURE.md).

The proof and scope received a separate internal mathematical review;
this is not an external review or a formal verification claim. Replay
the finite exact controls from the repository root with:

```sh
python3 research/factorial/image-rigidity/verify_sfc_reduced_divisor_contact.py
```
