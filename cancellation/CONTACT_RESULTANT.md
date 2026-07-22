# Endpoint reduction for the cancellation contact resultant

For positive integers `m,r`, put

\[
 p(v)=v(1-v)^m,
 \quad
 K_{m,r}(w)=w^{-r-1}\int_0^w p(v)^r\,dv,
\]

and

\[
 L_{m,r}(w)=w^{-r-1}\int_0^w\{p(w)-p(v)\}^r\,dv.
\]

The full cancellation prime-intersection diagram requires
`Res_w(K_{m,r},L_{m,r}) != 0`.  This note records a uniform endpoint
reduction and proves the two complete columns `r=1,2`.  It does **not** claim
the remaining all-parameter theorem for `r>=3`.

## 1. Triangular beta-primitive identity

Define `M_k(w)=K_{m,k}(w)`, including `M_0=1`.  Scaling `v=wt` and expanding
the binomial gives

\[
 \boxed{
 L_{m,r}(w)=\sum_{k=0}^r(-1)^k\binom rk
 (1-w)^{m(r-k)}M_k(w).}                               \tag{1}
\]

The last summand is `(-1)^r K_{m,r}`.  Thus modulo `K_{m,r}`, the contact
coefficient is a triangular combination of lower incomplete-beta
primitives.  Also

\[
 K_{m,r}(w)=\frac1{r+1}
 {}_2F_1(-mr,r+1;r+2;w),                              \tag{2}
\]

so (1) is the useful contiguous structure behind the proposed
hypergeometric search.

## 2. Fixed-`r` endpoint equations

Put

\[
 y=1-w,\qquad z=y^m.
\]

After `u=1-wt`,

\[
 M_k(w)=\frac{1}{(1-y)^{k+1}}
 \left(\beta_k-yz^kT_k(y)\right),                    \tag{3}
\]

where

\[
 \beta_k=\frac{k!}{\prod_{j=1}^{k+1}(mk+j)},\qquad
 T_k(y)=\sum_{j=0}^k\frac{(-1)^j\binom kj}{mk+j+1}y^j.
                                                                    \tag{4}
\]

Equations (1), (3), and (4) turn a degree-`mr` resultant in `w` into two
degree-`r` equations in the auxiliary endpoint variable `z`, together with
the single exponential relation `z=y^m`.  For each fixed `r`, eliminating
`z` has degree independent of the growing degree `mr`.  This is a sharper
symbolic-search formulation than interpolating the raw Sylvester resultant.

The points `y=1` and `y=0` never give a common root: they are respectively
`w=0`, where `K_{m,r}=1/(r+1)`, and `w=1`, where
`K_{m,r}=Beta(r+1,mr+1)`.

## 3. The complete column `r=1`

Equation (1) becomes

\[
 L_{m,1}+K_{m,1}=(1-w)^m.                             \tag{5}
\]

A common root would therefore be `w=1`, but

\[
 K_{m,1}(1)=B(2,m+1)=\frac1{(m+1)(m+2)}\ne0.
\]

More precisely,

\[
 \boxed{
 \operatorname{Res}_w(K_{m,1},L_{m,1})
 =\big((m+1)(m+2)\big)^{-m}.}                         \tag{6}
\]

This also gives a literal Bézout certificate.  Let `U_m(w)` be the truncation
through degree `m-1` in `1-w` of the Taylor series `1/K_{m,1}(w)` at `w=1`.
Then `U_mK_{m,1}-1` is divisible by `(1-w)^m=K_{m,1}+L_{m,1}`.  Dividing it
produces explicit polynomials `A_m,B_m` with

\[
 A_mK_{m,1}+B_mL_{m,1}=1.                             \tag{7}
\]

## 4. The complete column `r=2`

At a common root, `M_2=K_{m,2}=0`.  Equation (1) reads

\[
 L_{m,2}=z^2-2zM_1+M_2.
\]

Neither `y` nor `z=y^m` is zero, so `L_{m,2}=0` forces `2M_1=z`.
Using (3)--(4), this is

\[
 zA_m(y)=2,                                           \tag{8}
\]

with

\[
 A_m(y)=m(m+1)y^2-2m(m+2)y+(m+1)(m+2).
\]

Substitute `z=2/A_m(y)` into `M_2=0` and clear the displayed nonzero
denominators.  Direct simplification gives

\[
 (m+1)^2(y-1)^3\bigl(m^2y-(m+2)^2\bigr)=0.            \tag{9}
\]

The first candidate was excluded above.  The second is

\[
 y=\frac{(m+2)^2}{m^2}>1.
\]

At this value

\[
 A_m(y)=\frac{2(m+2)(5m^2+10m+4)}{m^3}>2,
\]

while `z=y^m>1`; hence `zA_m(y)>2`, contradicting (8).  Therefore

\[
 \boxed{\operatorname{Res}_w(K_{m,2},L_{m,2})\ne0
 \quad\text{for every }m\ge1.}                       \tag{10}
\]

## 5. Remaining step

For `r=3`, eliminating `z` factors off the excluded endpoint `y=1` but leaves
a genuine degree-six polynomial in `y` with coefficients in `Q[m]`.  It does
not reduce to the linear factor seen in (9).  The all-parameter problem is
therefore narrowed, not closed: one must show that the fixed-`r` eliminant
has no root compatible with `z=y^m`, uniformly for every `r>=3`.

Promising next steps are creative telescoping applied to (1)--(4), or a
subresultant recurrence in `r` for the endpoint equations.  Raw interpolation
of `Res_w(K,L)` obscures this structure because its degree grows with `mr`.

The exact symbolic certificate and bounded general-identity regression are
[`verify_contact_resultant_endpoint_reduction.py`](../scripts/verify_contact_resultant_endpoint_reduction.py).
