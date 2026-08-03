# Audit of the normalized-boundary \(A_4\) assembly

## 1. Outcome

The normalized-boundary incidence proposed at the
[two-mask frontier](A4_TWO_MASK_FACTORIZATION_FRONTIER.md) can be completed
to an explicit polynomial automorphism of affine five-space.  It therefore
gives a canonical candidate for the desired factorization.

The candidate is not a Keller map:

1. the two inverse-mask numerators are not divisible by the target cubic;
2. its rational Jacobian is a nonconstant quotient of two target-cubic
   pullbacks.

More generally, no polynomial automorphism of the target can both move the
generic mask-zero section into the target cubic and turn this factorization
into a Keller map.  A surviving construction must use a nonautomorphic
log-crepant incidence rechart satisfying a new relative Jacobian identity.
The first such incidence exists explicitly, but an order argument excludes
the entire radial-base family containing it, even with arbitrary polynomial
mask outputs.

Thus the ordinary polynomial Keller theorem has **not** been proved.
What is proved is an exact obstruction to the first complete assembly.

## 2. Unimodular completion of the normalized boundary

Let

\[
\begin{aligned}
f_1(r)&=2r^3-3r^2+3r-1,\\
f_2(r)&=-r^3+3r-1,\\
f_3(r)&=r(r-1).
\end{aligned}
\]

Then

\[
\mathcal B\bigl(\lambda f_1(r),
                 \lambda f_2(r),
                 \lambda f_3(r)\bigr)=0.            \tag{2.1}
\]

The column \(f=(f_1,f_2,f_3)^t\) is unimodular over \(\mathbb Q[r]\).
Indeed,

\[
(2r-1)f_2+(2r^2+r-5)f_3=1.                          \tag{2.2}
\]

Put

\[
a=2r-1,\qquad b=2r^2+r-5
\]

and

\[
\boxed{
M(r)=
\begin{pmatrix}
f_1&1&0\\
f_2&0&b\\
f_3&0&-a
\end{pmatrix}.
}                                                    \tag{2.3}
\]

Equation (2.2) gives

\[
\boxed{\det M(r)=1.}                                \tag{2.4}
\]

Hence

\[
\begin{aligned}
\alpha_{\rm norm}:
(\lambda,r,c,S,T)\longmapsto
\bigl(
&p=f_1(r)\lambda+S,\\
&q=f_2(r)\lambda+b(r)T,\\
&\rho=f_3(r)\lambda-a(r)T,\\
&A=c,\ B=r
\bigr)                                              \tag{2.5}
\end{aligned}
\]

is a polynomial automorphism of \(\mathbb A^5\) with Jacobian one.  On
\(S=T=0\), its first three coordinates are the normalized target-boundary
parametrization (2.1).

This proves that ambient polynomial completion of the normalized incidence
is not the obstruction.

## 3. The assembled rational candidate

Let

\[
\mathcal Y
=\widehat\Phi\times\operatorname{id}
=\left(
WN_1,\ WN_2,\ WH,\ \frac{WL}{4}z_1,\ z_2
\right).
                                                               \tag{3.1}
\]

The proposed assembly is

\[
F_{\rm rat}
=\beta_2^{-1}\circ\alpha_{\rm norm}\circ\mathcal Y, \tag{3.2}
\]

where \(\beta_2\) is the two-mask blowdown with

\[
\det D\beta_2=\mathcal B(p,q,\rho).
\]

In abstract coordinates \((\lambda,r,c,S,T)\), put

\[
D_{\rm new}
=\mathcal B\bigl(
f_1\lambda+S,\
f_2\lambda+bT,\
f_3\lambda-aT
\bigr).                                             \tag{3.3}
\]

Let

\[
X=p-q,\qquad
\mathcal U=27\rho^2-X^2-3(X-3\rho)q.
\]

The two inverse-mask numerators are

\[
\begin{aligned}
n_1&=c-Xr,\\
n_2&=-\mathcal Uc+27\rho^3r.                        \tag{3.4}
\end{aligned}
\]

Exact division gives

\[
\boxed{
D_{\rm new}\nmid n_1,\qquad
D_{\rm new}\nmid n_2.
}                                                    \tag{3.5}
\]

Therefore \(F_{\rm rat}\) has poles and is not a polynomial map.

## 4. The determinant mismatch

Since

\[
\det D\mathcal Y
=\mathcal B(\lambda,r,c),
\qquad
\det D\alpha_{\rm norm}=1,
\]

the chain rule gives

\[
\boxed{
\det DF_{\rm rat}
=\frac{\mathcal B(\lambda,r,c)}
       {D_{\rm new}}.
}                                                    \tag{4.1}
\]

At \(S=T=0\),

\[
D_{\rm new}=0
\]

by (2.1), while \(\mathcal B(\lambda,r,c)\) is generically nonzero.
Thus (4.1) is not constant even as a rational function.

The proposed assembly therefore fails both necessary tests independently:
polynomiality and constant Jacobian.

## 5. Automorphic incidence obstruction

The mismatch is structural.  Let

\[
\alpha\in\operatorname{Aut}(\mathbb A^5)
\]

have constant Jacobian \(u\ne0\), and consider

\[
F_{\alpha,\rm rat}
=\beta_2^{-1}\circ\alpha\circ\mathcal Y.
\]

The chain rule gives

\[
\det DF_{\alpha,\rm rat}
=u\,
\frac{\mathcal B(P,Q,R)}
{\mathcal B\bigl(
\operatorname{pr}_{1,2,3}(\alpha(P,Q,R,S,T))
\bigr)}.                                            \tag{5.1}
\]

The map \(\mathcal Y\) is dominant.  Hence, if (5.1) were a nonzero
constant, then already in the target polynomial ring

\[
\boxed{
\mathcal B\circ\operatorname{pr}_{1,2,3}\circ\alpha
=v\,\mathcal B(P,Q,R)
}                                                    \tag{5.2}
\]

for some \(v\in\mathbb Q^*\).

Equation (5.2) says that \(\alpha\) preserves the divisor

\[
V(\mathcal B)\times\mathbb A^2.
\]

In particular, it cannot send the generic point of the old mask-zero
section

\[
\{S=T=0\}
\]

into \(V(\mathcal B)\), because \(\mathcal B(P,Q,R)\) is nonzero at that
generic point.

Therefore:

\[
\boxed{
\begin{gathered}
\text{No automorphic incidence rechart can both move the generic}\\
\text{mask-zero section into \(\mathcal B=0\) and produce a Keller
factorization.}
\end{gathered}
}                                                    \tag{5.3}
\]

This closes the normalized-boundary strategy in its automorphic form.

## 6. Correct surviving equation

A nonautomorphic incidence map

\[
\alpha:\mathbb A^5\longrightarrow\mathbb A^5
\]

could still work, but it must satisfy the log-crepant identity

\[
\boxed{
\mathcal B\bigl(
\operatorname{pr}_{1,2,3}\alpha
\bigr)
=u\,\mathcal B(P,Q,R)\det D\alpha
}                                                    \tag{6.1}
\]

for some \(u\in\mathbb Q^*\), together with polynomial divisibility of both
inverse-mask numerators.

Indeed, if

\[
F=\beta_2^{-1}\circ\alpha\circ\mathcal Y
\]

is polynomial, (6.1) gives \(\det DF=u^{-1}\).

The next search must therefore solve three coupled requirements:

1. the log-Jacobian equation (6.1);
2. both adjugate divisibilities from (3.4);
3. generic reconstruction of the original \(A_4\) target field.

Dropping (6.1) and solving only the boundary embedding problem cannot
produce a Keller map.  This is the correction supplied by the present
audit.

## 7. First nonautomorphic solution and its residue obstruction

The log-crepant equation (6.1) is not empty.  Because \(\mathcal B\) is
homogeneous of degree three, the polynomial map

\[
\boxed{
\alpha_{\rm sc}(P,Q,R,S,T)
=(SP,SQ,SR,S,T)
}                                                    \tag{7.1}
\]

is generically birational and satisfies

\[
\det D\alpha_{\rm sc}=S^3,
\qquad
\mathcal B(SP,SQ,SR)=S^3\mathcal B(P,Q,R).
                                                               \tag{7.2}
\]

Thus (6.1) holds exactly with \(u=1\).  This is the first explicit
nonautomorphic incidence rechart passing the determinant gate.

It nevertheless fails the adjugate gate.  Write

\[
X=P-Q,
\qquad
\mathcal U=27R^2-X^2-3(X-3R)Q.
\]

For the scaled base one has \(X_{\rm sc}=SX\) and
\(\mathcal U_{\rm sc}=S^2\mathcal U\).  The two inverse-mask numerators
become

\[
\begin{aligned}
n_{1,{\rm sc}}
 &=S-SXT
   =S(1-XT),\\
n_{2,{\rm sc}}
 &=-S^2\mathcal U\,S+27(SR)^3T
   =S^3(-\mathcal U+27R^3T),                       \tag{7.3}
\end{aligned}
\]

whereas the denominator is

\[
\mathcal B(SP,SQ,SR)=S^3\mathcal B(P,Q,R).         \tag{7.4}
\]

After cancelling the one visible copy of \(S\), the first numerator is
congruent to \(1\) modulo \((S,T)\).  Consequently (7.4) cannot divide
it.  In particular this failure persists after the pure-target lift,
whose first mask is \(S=(WL/4)z_1\): setting \(z_1=T=0\) leaves the same
unit residue.

Therefore

\[
\boxed{
\text{homogeneous scaling solves the log-Jacobian equation but cannot
solve the first inverse-mask divisibility.}
}                                                    \tag{7.5}
\]

The remaining search must alter the mask outputs together with the radial
base scaling; a base-only log-crepant contraction is insufficient.  In fact
the next argument shows that arbitrary mask alterations do not rescue this
radial base either.

### 7.1 All polynomial masks over the radial base are impossible

Let \(A,C\in\mathbb Q[P,Q,R,S,T]\) be arbitrary and consider

\[
\alpha_{A,C}=(SP,SQ,SR,A,C).                        \tag{7.6}
\]

Put

\[
\mathfrak m=(P,Q,R,S),
\qquad
E=P\partial_P+Q\partial_Q+R\partial_R,
\qquad
\mathcal D=S\partial_S-E.
\]

An exact block-determinant calculation gives

\[
\boxed{
\det D\alpha_{A,C}
=S^2\bigl((\mathcal D A)C_T-A_T(\mathcal D C)\bigr).
}                                                    \tag{7.7}
\]

If (6.1) holds, homogeneity of \(\mathcal B\) and cancellation in the
polynomial domain give

\[
(\mathcal D A)C_T-A_T(\mathcal D C)=u^{-1}S.        \tag{7.8}
\]

Write

\[
a=\operatorname{ord}_{\mathfrak m}A,
\qquad
c=\operatorname{ord}_{\mathfrak m}C.
\]

The operators \(\mathcal D\) and \(\partial_T\) do not decrease
\(\mathfrak m\)-adic order.  Therefore the left side of (7.8) has order at
least \(a+c\), while the right side has order one.  Hence

\[
\boxed{a+c\leq1.}                                   \tag{7.9}
\]

Now pull back by the pure-target fivefold map

\[
\mathcal Y=\left(WN_1,WN_2,WH,\frac{WL}{4}z_1,z_2\right).
\]

For every nonzero polynomial \(f(P,Q,R,S,T)\),

\[
\operatorname{ord}_W\mathcal Y^*f
=\operatorname{ord}_{\mathfrak m}f.                \tag{7.10}
\]

Indeed, the leading form of \(f\) is evaluated on

\[
[N_1:N_2:H:(L/4)z_1].
\]

This map is dominant onto \(\mathbb P^3\): the rational map
\((U,V)\mapsto[N_1:N_2:H]\) is dominant onto \(\mathbb P^2\), as witnessed
by its nonzero Jacobian \(4K^3L/H^3\), and \(z_1\) supplies the fourth
projective coordinate.  Thus a nonzero leading form cannot vanish
identically.

For (7.6), the first inverse numerator is

\[
n_1=A-S(P-Q)C.                                      \tag{7.11}
\]

Its denominator is \(S^3\mathcal B(P,Q,R)\).  After \(\mathcal Y\), that
denominator is

\[
\left(\frac{WLz_1}{4}\right)^3
W^3K^3L^2
=\frac{W^6K^3L^5z_1^3}{64}.                        \tag{7.12}
\]

Polynomial divisibility would force
\(\operatorname{ord}_W\mathcal Y^*n_1\geq6\), hence by (7.10)

\[
\operatorname{ord}_{\mathfrak m}n_1\geq6.          \tag{7.13}
\]

But (7.11) then gives

\[
a=\operatorname{ord}_{\mathfrak m}A
\geq\min(6,c+2)\geq2,                              \tag{7.14}
\]

contradicting (7.9).  Therefore

\[
\boxed{
\begin{gathered}
\text{No radial-base incidence }(SP,SQ,SR,A,C)\text{ can satisfy both}\\
\text{the log-crepant identity and polynomial factorization after
\(\mathcal Y\).}
\end{gathered}
}                                                    \tag{7.15}
\]

This is an all-degree theorem and uses only the first adjugate numerator.
It closes radial contraction of the cubic base, not general
nonautomorphic incidence maps whose first three coordinates are coupled
nonradially.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_normalized_boundary_assembly.py
```

The checker verifies the normalized boundary, its unimodular completion,
the two failed divisions, the explicit nonconstant determinant ratio, and
the homogeneous nonautomorphic log-crepant incidence together with its
unit inverse-mask residue.
