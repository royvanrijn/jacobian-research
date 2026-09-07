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

## 8. Where the nonradial coupling lives

The two scalar divisibility tests in (3.4) are the rows of one matrix
factorization.  For arbitrary cubic-base coordinates \((p,q,\rho)\), put

\[
 X=p-q,\qquad
 \mathcal U=27\rho^2-X^2-3(X-3\rho)q
\]

and

\[
 \mathsf N=
 \begin{pmatrix}
  1&-X\\
  -\mathcal U&27\rho^3
 \end{pmatrix},qquad
 \mathsf M=\operatorname {adj}(\mathsf N)=
 \begin{pmatrix}
  27\rho^3&X\\
  \mathcal U&1
 \end{pmatrix}.                                    \tag{8.1}
\]

The defining cubic identity is exactly

\[
 \boxed{
 \det\mathsf N=\det\mathsf M=\mathcal B(p,q,\rho),
 \qquad
 \mathsf N\mathsf M=\mathcal B(p,q,\rho)I_2.}      \tag{8.2}
\]

The kernel after extension to the full source ring is not yet the desired
object: because the upper-left entry of \(\mathsf N\) is one, that kernel is
tautologically generated.  The nontrivial question is whether a target
polynomial acquires a new boundary factor only after the pure \(A_4\) lift.
Let

\[
 \mathcal A=\mathbb Q[P,Q,R,S,T],\qquad
 \Gamma=\mathbb Q[U,V,W,z_1,z_2],\qquad
 \iota=\mathcal Y^*:\mathcal A\hookrightarrow\Gamma,
\]

and for a chosen nonradial base triple put

\[
 B_\pi=\mathcal B(p,q,\rho),\qquad d=\iota(B_\pi).
\]

For a target mask pair \(v=(A,C)^t\), set \(g=A-XC\).  The two numerator
rows reduce exactly to

\[
 \mathsf Nv=
 \binom{g}{B_\pi C-\mathcal U g}.                  \tag{8.3}
\]

Consequently both rows become divisible by \(d\) after pullback if and only
if \(\iota(g)\in d\Gamma\).  Define the contraction ideal and the descended
coupling module by

\[
 \boxed{
 \begin{aligned}
 \mathcal J_\pi&=\iota^{-1}(d\Gamma)
                 =\mathcal A\cap d\Gamma,\\
 \mathcal C_\pi&=
 \frac{\{(A,C)\in\mathcal A^2:A-XC\in\mathcal J_\pi\}}
 {\mathcal A(X,1)+\mathcal A(B_\pi,0)}
 \cong\frac{\mathcal J_\pi}{(B_\pi)}.
 \end{aligned}}                                    \tag{8.4}
\]

This is the exact location of a different selector.  If \(\iota\) were flat
along \(d=0\), then \(\mathcal J_\pi=(B_\pi)\) and the module would vanish.
A surviving class records the nonflat boundary contraction exploited by the
root-chart selectors.  It should be computed by elimination and saturation
in the target ring, not as an unrestricted kernel over \(\Gamma\).  Once a
class \(g\) and a mask \(C\) are chosen, write \(\iota(g)=dk\); the two
inverse masks are then exactly

\[
 k,\qquad \iota(C)-\iota(\mathcal U)k.
\]

Concretely, in
\(\mathcal A\otimes_{\mathbb Q}\Gamma\), let \(G_\iota\) be the five graph
relations \(Y_i-\mathcal Y_i(U,V,W,z_1,z_2)\), where
\(Y_i=P,Q,R,S,T\).  Then

\[
 \mathcal J_\pi=(G_\iota+(d))\cap\mathcal A.
\]

This is a direct elimination computation.  Saturating before contraction
separates genuine horizontal classes from the fixed \(W,K,L,z_1\)
components already present in the ledger.

There is a decisive support test before that elimination.  Suppose
\(B_\pi\) is irreducible and an irreducible factor \(h\) of \(d\) does not
divide \(\det D\mathcal Y\).  At the generic point of \(h=0\), the map
\(\mathcal Y\) is etale.  Its restriction to \(h=0\) consequently has rank
four and maps dominantly to the irreducible fourfold \(B_\pi=0\).  If
\(g\in\mathcal J_\pi\), then \(h\mid\iota(g)\), so \(g\) vanishes on a
dense subset of \(B_\pi=0\).  Therefore \(B_\pi\mid g\), and

\[
 \boxed{
 h\nmid\det D\mathcal Y
 \quad\Longrightarrow\quad
 \mathcal J_\pi=(B_\pi),\qquad \mathcal C_\pi=0.}
\]

Thus a nonzero coupling class requires every component of the pulled
boundary to lie in the critical divisor.  For the pure \(A_4\) lift,

\[
 \det D\mathcal Y=W^3K^3L^2,
\]

so, for an irreducible prospective boundary, the necessary condition is

\[
 \boxed{
 \mathcal C_\pi\ne0
 \quad\Longrightarrow\quad
 \operatorname {rad}\!\left(\iota(B_\pi)\right)\mid WKL.}
\]

This criterion closes the normalized nonradial triple (2.5) without a
large graph elimination.  Its boundary \(B_\pi=D_{\rm new}\) is
irreducible because \(\alpha_{\rm norm}\) is an automorphism and
\(\mathcal B\) is irreducible.  On \(W=0\), the pure lift is

\[
 (\lambda,r,c,S,T)=(0,0,0,0,z_2),
\]

and the normalized base triple becomes

\[
 (p,q,\rho)=(0,-5z_2,z_2).
\]

Hence the exact boundary specialization is

\[
 d\big|_{W=0}=\mathcal B(0,-5z_2,z_2)=-133z_2^3.
\]

It follows at once that none of \(W,K,L\) divides \(d\): if \(K\) or \(L\)
did, the same factor, being independent of \(W\), would divide the displayed
nonzero constant in \(U,V\).  Thus \(d\) is coprime to the complete critical
divisor, and

\[
 \boxed{
 \mathcal J_{\rm norm}=(D_{\rm new}),\qquad
 \mathcal C_{\rm norm}=0.}
\]

Operationally, saturation away from \(WKL\) is therefore a rejection test,
not merely cleanup: any surviving horizontal factor already forces the
contraction module to vanish.  The expensive contraction and syzygy
calculation should be run only after the pulled boundary passes the
exceptional-support condition above.

The old cubic boundary does pass that test, and its contraction module is
nonzero.  In fact its first class is unexpectedly small.  Put

\[
 C_K=Q^2+3QR+9R^2.
\]

The cubic and pure-lift identities are

\[
\begin{aligned}
 \mathcal B(P,Q,R)
   &=P^3+(2Q+3R-3P)C_K,\\
 \iota(C_K)&=W^2\rho_VK^3,\qquad
   \rho_V=V^2+3V+9,\\
 \iota(\mathcal B)&=W^3K^3L^2.
\end{aligned}
\]

Consequently

\[
 \boxed{
 g_{\min}=S^2C_K\in\mathcal J_{\mathcal B},\qquad
 \frac{\iota(g_{\min})}{\iota(\mathcal B)}
   =\frac{W\rho_Vz_1^2}{16}.}
\]

The class is nonzero modulo \((\mathcal B)\).  Moreover, it is the first
one.  Decomposition by the independent \(z_1,z_2,W\) powers reduces target
degree at most four to homogeneous divisibility kernels in
\(\mathbb Q[P,Q,R]\).  Exact coefficient matrices give

\[
\begin{aligned}
 K^3L^2\mid F(N_1,N_2,H),\ \deg F\le4
   &\Longleftrightarrow F\in(\mathcal B),\\
 K^3L\mid F(N_1,N_2,H),\ \deg F\le3
   &\Longleftrightarrow F\in(\mathcal B),\\
 K^3\mid F(N_1,N_2,H),\ \deg F\le2
   &\Longleftrightarrow F\in(C_K).
\end{aligned}
\]

It follows that the contraction quotient has no nonzero class below degree
four and

\[
 \boxed{
 (\mathcal J_{\mathcal B})_4
 =\mathcal B\langle P,Q,R,S,T\rangle
   \oplus\mathbb Q\,S^2C_K.}
\]

The contraction can in fact be completed in every degree.  Exact
elimination of the nonreduced root-chart ideal gives

\[
 \iota_0^{-1}((K^3))=(C_K,P^3)
 \quad\text{in }\mathbb Q[P,Q,R],
\]

where \(\iota_0(P,Q,R)=(N_1,N_2,H)\).  On the other hand, \(L=0\) maps
nonconstantly to the irreducible projective cubic \(\mathcal B=0\): the
tangential derivative of \(N_1/H\) is nonzero modulo \(L\).  Its image is
therefore dense in that cubic.  Every homogeneous base form whose pullback
has an \(L\)-factor is consequently divisible by \(\mathcal B\).
Coefficientwise decomposition in \(W,z_1,z_2\) now gives the full formula

\[
 \boxed{
 \mathcal J_{\mathcal B}
 =(\mathcal B,\ S^2C_K,\ S^2P^3),\qquad
 \mathcal C_{\mathcal B}
 \cong
 \frac{S^2(C_K,P^3)}
      {(\mathcal B)\cap S^2(C_K,P^3)}.}
\]

In particular every reduced representative of a nonzero contraction class
is double along \(S=0\).  This also closes every base-fixed use of the
module.  If \(g=\mathcal Bh+S^2f\) and \(A-(P-Q)C=g\), then

\[
 \det\frac{\partial(A,C)}{\partial(S,T)}
 =g_SC_T-g_TC_S\in(\mathcal B,S).
\]

It cannot be the nonzero constant required by the log-crepant identity
when the three base outputs are fixed.  A surviving use of the module must
therefore combine a \(\mathcal B\)-multiple representative, a nonradial
change of all three base outputs, and the log equation in one incidence.

Even that does not suffice if the new base preserves the old boundary.
Suppose

\[
 \mathcal B(p,q,\rho)=v\mathcal B(P,Q,R),\qquad v\in\mathbb Q^*,
\]

and both inverse masks divide after the pure lift.  Then
\(g=A-(p-q)C\) lies in \(\mathcal J_{\mathcal B}\), so

\[
 g=\mathcal Bh+S^2f.
\]

On the smooth locus \(\mathcal B=S=0\), one has \(dg=h\,d\mathcal B\).
Differentiating the displayed boundary identity puts \(d\mathcal B\) in
the span of \(dp,dq,d\rho\).  The identity

\[
 dA-(p-q)dC-C(dp-dq)=dg
\]

then makes \(dp,dq,d\rho,dA,dC\) linearly dependent.  Hence
\(\det D\alpha=0\) there.  But the log-crepant equation (6.1) would force
\(\det D\alpha=(uv)^{-1}\), a nonzero constant.  This contradiction proves

\[
 \boxed{
 \text{No boundary-preserving, possibly nonradial incidence can use }
 \mathcal J_{\mathcal B}\text{ to produce the Keller factorization}.}
\]

Thus the surviving search must genuinely replace the base boundary while
keeping its entire pullback supported on \(WKL\).  It must also reconstruct
the old \(A_4\) function field; replacing \(\mathcal B\) by a high power
typically introduces an additional radical extension and is not harmless.

The first degrees of that replacement search are now closed.  If
\(\iota(B_\pi)\) is supported on \(WKL\), independence of \(z_1,z_2\)
forces \(B_\pi\) to be independent of \(S,T\).  Distinct homogeneous
degrees in \(P,Q,R\) have distinct \(W\)-orders, so the same support
condition forces \(B_\pi\) to be homogeneous and gives an identity

\[
 B_\pi(N_1,N_2,H)=cK^aL^b,\qquad c\in\mathbb Q^*.
\]

For every degree \(1\le d\le6\), all pairs with
\(2a+3b\le4d\) were tested by exact finite-field linear algebra.  Modulo
101 the homogeneous substitution columns remain independent, and the only
members \(K^aL^b\) in their span are

\[
\begin{array}{c|c}
d&(a,b)\\ \hline
3&(3,2)\\
6&(6,4).
\end{array}
\]

The exact characteristic-zero identities are \(\mathcal B\) and
\(\mathcal B^2\).  Because a primitive rational identity would reduce to a
nonzero modular identity, this proves

\[
 \boxed{
 \text{through base degree six, every exceptional-support boundary is a
 power of }\mathcal B.}
\]

Thus a genuinely different reduced boundary starts in degree at least
seven.  The degree-six power \(\mathcal B^2\) is a separate nonreduced
possibility, but no polynomial triple \((p,q,\rho)\) with
\(\mathcal B(p,q,\rho)=\mathcal B^2\), the required log Jacobian, and
field reconstruction is supplied by this sieve.

The lowest-degree version of that nonreduced possibility can in fact be
closed.  There is no triple

\[
 (p,q,\rho)\in\mathbb Q[P,Q,R,S,T]^3,
 \qquad \max(\deg p,\deg q,\deg\rho)\le2,
\]

and no \(\nu\in\mathbb Q^*\) such that

\[
 \boxed{\mathcal B(p,q,\rho)=\nu\mathcal B(P,Q,R)^2.}       \tag{8.4a}
\]

Here is the proof.  First set \(S=T=0\) and take degree-six homogeneous
parts.  A putative solution of (8.4a) gives homogeneous quadrics
\(f_0,f_1,f_2\) on \(\mathbb P^2\) satisfying

\[
 \mathcal B(f_0,f_1,f_2)=\nu\mathcal B^2.                 \tag{8.4b}
\]

They have no common polynomial factor: the cube of such a factor would
divide the square of the irreducible cubic \(\mathcal B\), while its degree
is at most two.  The rational map
\(f=[f_0:f_1:f_2]\) is dominant.  Indeed, if its image were a curve not
contained in \(\mathcal B=0\), the generic point of the irreducible source
cubic would map to one point of the finite intersection of those two
curves.  Two independent target lines through that point would then pull
back to quadrics vanishing on an irreducible cubic, hence to zero.  This
would make \(f\) constant, contrary to (8.4b).

If there is no basepoint, \(f\) is a holomorphic endomorphism of
\(\mathbb P^2_{\mathbb C}\) of degree two for which the irreducible cubic
is totally invariant.  This is impossible: a totally invariant curve of a
holomorphic endomorphism of \(\mathbb P^2\) of degree at least two is a
union of at most three lines
([Favre--Jonsson, Proposition 1.1](https://aif.centre-mersenne.org/articles/10.5802/aif.1985/)).

Suppose instead that basepoints occur.  A proper basepoint is simple: if
all three quadrics vanished to order two there, they would be binary
quadrics in two linear coordinates and the map would not be dominant.  The
cubic is irreducible and has its unique ordinary node at
\([P:Q:R]=[1:1:0]\); in the chart \(P=1,Q=1+x,R=y\), its tangent cone is

\[
 3x(x+3y).
\]

At a proper simple basepoint the left side of (8.4b) has order at least
three, whereas the right side has order twice the multiplicity of the
cubic.  Hence the proper basepoint must be the node.  After blowing it up,
the three linear initial forms map the exceptional line into
\(\mathcal B=0\): the degree-three initial term on the left of (8.4b) must
vanish because the right side has order four.  A line cannot lie in the
irreducible cubic, so this exceptional map is constant.  The three linear
forms consequently share one linear factor, giving exactly one
infinitely-near basepoint.  It is simple, since multiplicity two there would
give self-intersection at most \(4-1-4<0\) for the pullback of a line under
a dominant map.  Dividing its cubic common exceptional
factor from (8.4b), the boundary divisor after the first blowup is

\[
 E_1+2\widetilde{\mathcal B}.                            \tag{8.4c}
\]

That second basepoint must therefore be one of the two transverse
points \(E_1\cap\widetilde{\mathcal B}\); away from them the displayed
divisor has multiplicity only one, not the required three.  After blowing
up that point, the new exceptional component has coefficient zero in the
pulled boundary; above it the remaining coefficients are only one and two,
so no further infinitely-near basepoint can occur.  The other branch point
was not the unique common zero of the initial forms.  Thus the two
basepoints just found resolve the quadratic map.  On the resulting surface,
in the total-transform basis \((H,E_1,E_2)\),

\[
\begin{aligned}
 D&=2H-E_1-E_2,\\
 K_X-f^*K_{\mathbb P^2}&=3H-2E_1-2E_2,\\
 [\widetilde{\mathcal B}]&=3H-2E_1-E_2.
\end{aligned}
\]

Normal ramification index two says that
\(\widetilde{\mathcal B}\) is a component of the effective ramification
divisor.  The residual class would be \(-E_2\), but

\[
 D\mathbin\cdot(-E_2)=-1,
\]

contradicting nefness of \(D=f^*H\).  This proves (8.4a).  In particular,
the \(\mathcal B^2\) hit is not a quadratic boundary replacement, even if
the quadratic base outputs use both mask variables.  Any attempt to use
that hit must start in degree at least three and arrange cancellation of
all terms above degree six.

This distinction also explains why the birational invariant-cubic
literature does not directly provide the missing coupling.  Identity
(8.4b) is *strong total invariance*, involving the complete pullback.
Cantat explicitly separates it from strict-transform invariance and notes
that cubic-preserving birational maps acquire exceptional components in
their total transforms
([Section 3](https://archive.intlpress.com/site/pub/files/_fulltext/journals/mrl/2010/0017/0005/MRL-2010-0017-0005-a003.pdf)).
The many invariant cubics available for birational surface maps therefore
do not satisfy the square-pullback identity needed here.

There is nevertheless a useful genuinely different boundary rechart hiding
in the same nodal geometry.  Move the normalization surface in the node
direction and define

\[
 \boxed{
 \chi(\mu,\lambda,t)=
 \bigl(
  \mu+\lambda f_1(t),
  \mu+\lambda f_2(t),
  \lambda f_3(t)
 \bigr).}                                             \tag{8.4d}
\]

Direct calculation gives

\[
\begin{aligned}
 \mathcal B\circ\chi
   &=27\mu\lambda^2t^3(t-1)^3,\\
 \det D\chi
   &=3\lambda t^2(t-1)^2.                             \tag{8.4e}
\end{aligned}
\]

Thus \(\chi\) is dominant and monomializes the cubic boundary.  This is
exactly the sort of different base map suggested by the normalization and
birational literature; it is not found by the homogeneous
\((P,Q,R)\)-identity sieve because of the node translation.

It still does not provide the missing coupling.  Put
\((\mu,\lambda,t)=(c,\lambda,r)\) in the old target coordinates.  Under the
pure lift,

\[
 \lambda=WMK,\qquad \mu=WH,\qquad t=WN_2,
 \qquad t-1=WN_2-1,
\]

and hence

\[
 \iota(\mathcal B\circ\chi)
 =27(WMK)^2(WH)(WN_2)^3(WN_2-1)^3.                   \tag{8.4f}
\]

Each target factor retains an irreducible etale witness: respectively

\[
 M,\qquad H,\qquad N_2,\qquad WN_2-1.
\]

All four are coprime to \(WKL\).  Applying the etale-component valuation
argument with the displayed multiplicities shows that every
\(g\in\mathcal J_\chi\) is divisible by
\(\lambda^2\mu t^3(t-1)^3\).  Therefore

\[
 \boxed{
 \mathcal J_\chi=(\mathcal B\circ\chi),\qquad
 \mathcal C_\chi=0.}                                 \tag{8.4g}
\]

Permuting the node-chord parameters among the available target coordinates
does not help.  For the coordinate hyperplanes the etale witnesses are

\[
\begin{array}{c|ccccc}
Y&P&Q&R&S&T\\ \hline
\text{etale factor of }\iota(Y)&M&N_2&H&z_1&z_2.
\end{array}
\]

For each of the five choices of \(Y\), the shifted pullback
\(\iota(Y)-1\) is also irreducible and coprime to \(WKL\).  Hence the same
valuation proof rejects all
\(5\cdot4\cdot3=60\) injective coordinate placements of
\((\mu,\lambda,t)\).

The old boundary itself can replace the failed normal coordinate.  Put

\[
 \chi_{\mathcal B}(\lambda,t)
 =\chi(\mathcal B,\lambda,t),
 \qquad
 h=\lambda^2t^3(t-1)^3.
\]

Then

\[
 \boxed{
 \mathcal B\circ\chi_{\mathcal B}=27\mathcal B h.}       \tag{8.4h}
\]

This boundary is genuinely different and, unlike (8.4d) with a coordinate
normal parameter, it has a nonzero coupling module.  For coordinate choices
of \(\lambda,t\), the etale witnesses force the exact factor \(h\) in every
contracted numerator.  Cancellation in the source domain then gives

\[
 \boxed{
 \mathcal J_{\chi_{\mathcal B}}=h\mathcal J_{\mathcal B},
 \qquad
 \mathcal C_{\chi_{\mathcal B}}
 \simeq h\mathcal J_{\mathcal B}/(h\mathcal B).}          \tag{8.4i}
\]

For example, take \((\lambda,t)=(P,Q)\) and

\[
\begin{aligned}
 (p_B,q_B,\rho_B)
   &=\bigl(\mathcal B+Pf_1(Q),\mathcal B+Pf_2(Q),Pf_3(Q)\bigr),\\
 h&=P^2Q^3(Q-1)^3,\\
 g_B&=27hS^2C_K,\\
 C_B&=T,\qquad A_B=(p_B-q_B)T+g_B.
                                                               \tag{8.4j}
\end{aligned}
\]

Both inverse-mask quotients are polynomial and are again

\[
 k_{\min}=\frac{W\rho_Vz_1^2}{16},\qquad
 z_2-\iota(\mathcal U_B)k_{\min},
\]

where \(\mathcal U_B\) is formed from \((p_B,q_B,\rho_B)\).  Thus this gives
an explicit polynomial different-boundary coupled map

\[
 \boxed{
 F_B=\bigl(
   \iota(p_B),\iota(q_B),\iota(\rho_B),
   k_{\min},z_2-\iota(\mathcal U_B)k_{\min}
 \bigr).}                                               \tag{8.4k}
\]

It is still not Keller.  Its exact Jacobian is

\[
 \boxed{
 \det DF_B=
 \iota\!\left(
  6PQ^2(Q-1)^2S C_K\,\partial_R\mathcal B
 \right),}                                              \tag{8.4l}
\]

which is nonconstant.

Thus (8.4k) is a verified polynomial coupling, not a Keller map; no
reconstruction of the original (A_4) function field is asserted for it.

In fact every coordinate version of this \(\mathcal B\)-twist is closed,
even with arbitrary mask feedback.  Let
\(g=A-(p_B-q_B)C\).  Since (8.4i) forces \(g=hg_0\), while

\[
 dp_B\wedge dq_B\wedge d\rho_B
 =3\lambda t^2(t-1)^2
   d\mathcal B\wedge d\lambda\wedge dt,
\]

the log equation \(\det D\alpha=27h\) would require

\[
 d\mathcal B\wedge d\lambda\wedge dt\wedge dg\wedge dC
 =9\lambda t(t-1)
  \,dP\wedge dQ\wedge dR\wedge dS\wedge dT.             \tag{8.4m}
\]

The left side is divisible by
\(h=\lambda^2t^3(t-1)^3\), because the terms containing \(dh\) disappear
in the wedge; the right side is not.  This contradiction rejects all
\(5\cdot4=20\) ordered coordinate choices of \((\lambda,t)\), for every
polynomial pair \((A,C)\), not just the pair in (8.4j).

This near miss identifies the selector more precisely.  A successful
node-chord deformation must replace a coordinate parameter by a descended
function whose pullback has a smaller exceptional multiplicity than the
full target factor forced by etale valuation.  Equivalently, it must evade
the factor \(h\) in (8.4i), not merely make the coupling module nonzero.
The exceptional divisors isolated in the root-chart selector calculations
are therefore the right raw material, but they must be inserted into the
base rechart before contraction rather than appended as an independent
mask.

Returning to the old-boundary calibration, the class (g_{min}) makes an
explicit polynomial coupled map.  Take

\[
 (A_{\min},C_{\min})
 =\bigl((P-Q)T+S^2C_K,\ T\bigr).
\]

Then \(A_{\min}-(P-Q)C_{\min}=g_{\min}\), and after the pure lift the two
inverse-mask quotients are

\[
 k_{\min}=\frac{W\rho_Vz_1^2}{16},\qquad
 z_2-\iota(\mathcal U)k_{\min}.
\]

Thus

\[
 F_{\min}=
 \left(
 WN_1,\ WN_2,\ WH,\
 \frac{W\rho_Vz_1^2}{16},\
 z_2-\iota(\mathcal U)\frac{W\rho_Vz_1^2}{16}
 \right)
\]

is polynomial.  It is not the desired map.  The incidence mask block has

\[
 \det\frac{\partial(A_{\min},C_{\min})}{\partial(S,T)}
 =2SC_K,
\]

and the composite Jacobian is

\[
 \boxed{
 \det DF_{\min}
 =\frac{W^3K^3L\rho_Vz_1}{2},}
\]

not a unit.  The computation therefore finds the different coupling but
also proves why a block-triangular use of it fails.  Any further deformation
must genuinely change the base boundary and its Jacobian simultaneously while
preserving the two quotient identities; merely adding feedback to this
reduced representative cannot turn the displayed block into a Keller map.

The columns of \(\mathsf M\) give only the zero classes in (8.4).  A
base-independent block-triangular use of them has determinant divisible by
\(\mathcal B\), so feedback from the two mask coordinates into the three
base coordinates is essential.  This is the precise algebraic meaning of
*distributed cancellation* in this problem.

Indeed, the direct adjugate construction

\[
 (P,Q,R,S,T)\longmapsto
 \left(P,Q,R,27R^3S+XT,\mathcal U S+T\right)
\]

makes both numerator rows divisible because
\(\mathsf N\mathsf M=\mathcal B I_2\).  But its full Jacobian is
\(\det\mathsf M=\mathcal B\).  With the base fixed, (6.1) would therefore
read \(\mathcal B=\mathcal B^2\).  Thus even the tautological coupled pair
fails: a solution must change the mask-zero base map and the mask block in
one incidence, rather than append this matrix factorization
block-triangularly.

For a representative target pair \(v=(A,C)\), put

\[
 g=A-XC,\qquad
 I_\pi=(B_\pi,g,B_\pi C-\mathcal U g)\subset\mathcal A.
\]

The natural target-side space on which the two rational quotients are
regular is the affine modification

\[
 X_\pi=\operatorname {Spec}\mathcal A[I_\pi/B_\pi]. \tag{8.5}
\]

The condition \(g\in\mathcal J_\pi\) is exactly what makes \(\mathcal Y\)
lift polynomially to \(X_\pi\), with quotient coordinates
\(k\) and \(\iota(C)-\iota(\mathcal U)k\).  The recognition problem is
whether a suitable lifted chart, possibly after the authorized
stabilization, is actually affine five-space.

This matches the repository's earlier
[log-suspension model](../cancellation/LOG_GEOMETRY_OF_SUSPENSIONS.md),
[controlled-boundary conclusion](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md),
and [conductor-three obstruction](CONDUCTOR_THREE_BOUNDARY_COX_FILL_OBSTRUCTION.md):
separated masks, symmetric Cox fills, and tensor products are too rigid;
an asymmetric mixed Rees algebra or distributed cancellation is the
surviving locus.  The affine-modification language is standard in
[Kaliman--Zaidenberg](https://arxiv.org/abs/math/9801076) and its global
Rees-algebra form is developed by
[Dubouloz](https://arxiv.org/abs/math/0503142).

Two tempting shortcuts can now be removed from the search.  First, the
quadratic polar map \(\nabla\mathcal B\) is not log-crepant.  At the smooth
boundary point

\[
 (P,Q,R)=(9,-3,2)
\]

one has

\[
 \mathcal B=0,\quad
 \nabla\mathcal B=(162,54,-648),\quad
 \mathcal B(\nabla\mathcal B)=-8533918944,\quad
 \det\operatorname {Hess}(\mathcal B)=1259712.      \tag{8.6}
\]

Second, multiplying \(\mathcal B\) by one arbitrary hyperplane does not
produce a Saito-free quartic.  For
\(\ell=e_0P+e_1Q+e_2R\), the coefficient matrix for a linear Jacobian
syzygy of \(\mathcal B\ell\) is \(15\)-by-\(9\), and the ideal of its
maximal minors is exactly

\[
 (e_0,e_1,e_2)^9.                                  \tag{8.7}
\]

Hence it has full column rank for every nonzero \(\ell\).  A free reduced
plane quartic would have exponents summing to three and therefore a
non-Euler syzygy of degree at most one, so (8.7) excludes every single
adjoint hyperplane.  This uses Saito's logarithmic criterion
([original paper](https://repository.dl.itc.u-tokyo.ac.jp/records/39646)).
Results in which an adjoint can complete a discriminant to a free divisor
require additional normalization and stability hypotheses; see
[Mond--Schulze](https://arxiv.org/abs/1001.1095).  Finite-flat discriminants
and their normalization matrices remain a plausible source of structured
logarithmic frames
([Buchweitz--Ebeling--von Bothmer](https://arxiv.org/abs/math/0612119)), but
none of these results supplies the required affine-five-space recognition
or the \(A_4\) field reconstruction.

The next bounded computation is therefore concrete:

1. choose a degree and monomial support for a nonradial triple
   \((p,q,\rho)\) with source-mask feedback, and first solve the exceptional
   support condition
   \(\operatorname {rad}(\iota(\mathcal B(p,q,\rho)))\mid WKL\), rejecting
   boundaries proportional to \(\mathcal B\); the reduced-boundary search
   now starts in degree seven; the degree-six nonreduced hit
   \(\mathcal B^2\) has no realization by base outputs of degree at most
   two, so its first possible realization has degree at least three and
   must cancel every term above degree six; use the node-chord map (8.4d)
   as the first calibration, replacing one of its four parameter factors
   by a non-coordinate selector with only \(W,K,L\)-supported pullback
   rather than one of the etale witnesses in (8.4f), since all sixty
   coordinate placements are closed; the first positive calibration is the
   \(\mathcal B\)-twist (8.4h), and its next search should choose
   non-coordinate \((\lambda,t)\) for which contraction does not force the
   full factor \(h\), since (8.4m) closes all twenty coordinate pairs;
2. for survivors, compute (8.4) and its syzygies; use
   \((\mathcal B,S^2C_K,S^2P^3)\) as the exact calibration showing why the
   old boundary cannot work;
3. impose (6.1), dominance, and both exact quotient identities;
4. verify reconstruction of the original \(A_4\) function field;
5. for (8.5), test smoothness, units, factoriality/class group, conductor,
   and the available motivic and topological obstructions before claiming
   affine space.

The first four steps are polynomial equations for every fixed support.  The
fifth is a recognition gate, not something implied by the preceding
identities.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_normalized_boundary_assembly.py
```

The checker verifies the normalized boundary, its unimodular completion,
the two failed divisions, the explicit nonconstant determinant ratio, and
the homogeneous nonautomorphic log-crepant incidence together with its
unit inverse-mask residue.  It also verifies the determinant-\(\mathcal B\)
matrix in (8.2), the direct-adjugate Jacobian obstruction, the polar
counterexample (8.6), and the maximal-minor ideal (8.7).  Finally, it checks
the pure-lift critical divisor and the specialization
\(d|_{W=0}=-133z_2^3\) used to prove that the normalized contraction module
vanishes.  For the old boundary it verifies the complete contraction
kernels through target degree four, eliminates the full \(K^3\)-contraction
as \((C_K,P^3)\), checks dominance of \(L=0\) onto the projective cubic, and
thereby proves
\(\mathcal J_{\mathcal B}=(\mathcal B,S^2C_K,S^2P^3)\).  It also checks
both polynomial inverse-mask quotients and the nonconstant Jacobian of the
resulting explicit minimal coupled map.  The accompanying differential
argument uses that full module formula to close every incidence whose new
base boundary is proportional to \(\mathcal B\).  Finally, an exhaustive
mod-101 linear sieve proves that through base degree six the only
exceptional-support identities are \(\mathcal B\) and \(\mathcal B^2\).
For the latter, the checker verifies the ordinary-node tangent cone and the
negative residual ramification intersection used in the written proof that
no degree-at-most-two base triple can pull the cubic back to its square.
The basepoint-free case of that proof uses the cited holomorphic
totally-invariant-curve theorem and is not a computer-algebra assertion.
It also verifies the dominant node-chord rechart (8.4d), its monomial
boundary and Jacobian, and the four irreducible etale factors proving that
its descended coupling module is zero.  The five coordinate-hyperplane
witnesses and five shifted-coordinate witnesses then reject all sixty
injective coordinate placements of its three parameters.  Finally, it
verifies the \(\mathcal B\)-twisted chord, the nonzero class (8.4j), both
polynomial quotients, and the nonconstant Jacobian (8.4l).  The
exterior-form divisibility argument (8.4m), together with the same etale
witnesses, closes all twenty coordinate pairs even for arbitrary masks.
