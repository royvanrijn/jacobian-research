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

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_normalized_boundary_assembly.py
```

The checker verifies the normalized boundary, its unimodular completion,
the two failed divisions, and the explicit nonconstant determinant ratio.
