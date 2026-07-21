# Classification with one additional resolvent factor

This note closes the second minimal relaxation of the C24 cancellation
skeleton.  The inverse derivative is allowed a second normalized factor
defined by an arbitrary polynomial.  The constant-Jacobian equation forces
that factor either to coincide with the original factor or to be trivial.
The same proof classifies an arbitrary finite product of normalized factors.
Consequently no polynomial cancellation branch beyond C24 occurs anywhere
in this factorized extension.

Throughout, `k` is a characteristic-zero field, `f_1,f_2 in k[y]`, with
`f_1` nonconstant, and `a,b>=1`.  To avoid confusing the second factor with
the correction term, write

\[
 A=1+xf_1(y),\qquad B=A^bz+\gamma(y,A),                       \tag{1}
\]

\[
 P=A^aB,\qquad Q=y+xA^{a-1}B,\qquad s={x\over A}.             \tag{2}
\]

Fix integers `u>=0` and `v>=1`, and set

\[
 D_i(t)=1-tf_i(Q-Pt)\quad(i=1,2),                             \tag{3}
\]

\[
 R=C\int_0^sD_1(t)^uD_2(t)^vdt,qquad C\in k^*.              \tag{4}
\]

These are the most general two normalized factors of this form; a nonzero
scalar multiplying `t f_i` is absorbed into `f_i`.  The expressions are
initially interpreted in `k[x,y,z,A^(-1)]`.

## 1. The two-factor Jacobian

Put

\[
 n=a+b-2.
\]

The two-weight determinant calculation gives

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}=-A^n.             \tag{5}
\]

On the source, `Q-Ps=y` and

\[
 D_1(s)=1-{x\over A}f_1(y)=A^{-1}.                            \tag{6}
\]

In the rational function field `k(y)`, put `rho=f_2/f_1`.  Since
`x=(A-1)/f_1`, the second factor is

\[
 D_2(s)
 ={A-xf_2(y)\over A}
 ={\rho+(1-\rho)A\over A}.                                  \tag{7}
\]

Holding `(P,Q)` fixed in (4) and combining (5)--(7) therefore gives

\[
 \boxed{\quad
 \det {\partial(P,Q,R)\over\partial(x,y,z)}
 =-C A^{n-u-v}\{\rho+(1-\rho)A\}^v.
 \quad}                                                       \tag{8}
\]

The correction `gamma(y,A)` again disappears completely.

## 2. Complete factor classification

### Theorem 2.1

The expression in (8) is a nonzero constant if and only if exactly one of
the following reductions occurs:

1. `f_2=f_1` and `n=u+v`; the two factors coincide and (4) has the single
   derivative `D_1^(u+v)`;
2. `f_2=0` and `n=u`; the second factor is identically one and (4) has the
   single derivative `D_1^u`.

In particular, no genuinely independent second normalized resolvent factor
is compatible with the Keller condition.

**Proof.**  Work in the Laurent polynomial ring
`k(y)[A,A^(-1)]`.  If (8) is a nonzero constant, then

\[
 L(A)=\rho+(1-\rho)A
\]

must be a unit of that ring, because `v>=1`.  Its units are precisely
`lambda A^j` with `lambda in k(y)^*` and `j in Z`.  Since `L` has degree at
most one and its two displayed coefficients sum to one, it is a unit only in
the following cases:

- `rho=1`, when `L=1`; or
- `rho=0`, when `L=A`.

In the first case (8) is `-C A^(n-u-v)`, which is constant exactly when
`n=u+v`.  In the second it is `-C A^(n-u)`, which is constant exactly when
`n=u`.  Conversely, direct substitution proves that both listed reductions
give the constant `-C`.  QED

### Corollary 2.2 (no polynomial branch beyond C24)

If `(P,Q,R)` in (1)--(4) is polynomial and has nonzero constant Jacobian,
then, up to translations, nonzero scalings, and polynomial left--right
equivalence, it is a C24 map.

**Proof.**  In both cases of Theorem 2.1, the inverse derivative reduces to

\[
 \{1-tf_1(Q-Pt)\}^n.                                       \tag{9}
\]

If `n=0`, then `R=Cs=Cx/A`, which is not polynomial.  Thus every polynomial
member has `n>=1`.  The
[generalized cancellation theorem](GENERALIZED_CANCELLATION_MECHANISM.md)
then applies with `e=n`: its weight rigidity forces `a=1`, `b=n+1`, and its
spectral classification proves that every polynomial member is
left--right equivalent to C24.  QED

Together with the
[one-additional-weight theorem](THREE_WEIGHT_CANCELLATION_CLASSIFICATION.md),
this settles both minimal relaxations proposed in Option B: neither an
independent third monomial weight nor one additional normalized resolvent
factor produces anything beyond C24.

## 3. Arbitrarily many normalized factors

The two-factor argument does not stop at two factors.  Let
`q_1,...,q_r in k[y]`, let `v_i>=1`, and replace (4) by

\[
 R=C\int_0^s\prod_{i=1}^r
 \{1-tq_i(Q-Pt)\}^{v_i}dt.                                \tag{10}
\]

The polynomial `f_1` defining `A=1+xf_1(y)` remains fixed.  Put
`rho_i=q_i/f_1`.  Exactly as in (7)--(8), the Jacobian is

\[
 -C A^{n-\sum_i v_i}
 \prod_{i=1}^r\{\rho_i+(1-\rho_i)A\}^{v_i}.                \tag{11}
\]

### Theorem 3.1 (finite factorized classification)

The expression in (11) is a nonzero constant if and only if every `q_i` is
either `f_1` or zero and

\[
 n=\sum_{q_i=f_1}v_i.                                      \tag{12}
\]

All factors with `q_i=0` are identically one, while all remaining factors
coalesce to

\[
 \{1-tf_1(Q-Pt)\}^n.                                      \tag{13}
\]

If the resulting coordinates are polynomial, the map is left--right
equivalent to C24.

**Proof.**  In `k(y)[A,A^(-1)]`, the product in (11) can be a unit only if
each factor `rho_i+(1-rho_i)A` is a unit.  The proof of Theorem 2.1 shows this
is equivalent to `rho_i in {0,1}`.  A factor with `rho_i=0` contributes
`A^(v_i)` and cancels its own `A^(-v_i)` contribution; a factor with
`rho_i=1` contributes no numerator power.  The remaining power of `A` is
therefore `n-sum_(q_i=f_1)v_i`; its vanishing is exactly (12).  This proves
the constant-Jacobian classification and (13).  If `n=0`, the primitive is
again `Cx/A` and is not polynomial.
For `n>=1`, the generalized cancellation theorem gives the C24 conclusion as
in Corollary 2.2.  QED

The result does not classify derivatives with factors depending essentially
on both `t` and the target parameters `(P,Q)`, nor coordinate skeletons with
additional source functions or variables.

## 4. Exact regression

Run

```bash
.venv/bin/python scripts/verify_two_factor_resolvent.py
```

The script independently differentiates representative nonlinear two- and
multi-factor instances and exhausts bounded weights, powers, and
second-factor polynomials.  The all-parameter result is the Laurent-unit
proof above.
