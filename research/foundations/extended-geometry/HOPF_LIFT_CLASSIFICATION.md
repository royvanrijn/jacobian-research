# Hopf-lift classification and minimal \(t^2\)-linear rigidity

## 1. Scope

This note isolates the mechanism behind the adjacent Taylor coefficients in
the two-pair Image-Mathieu counterexample.  It first classifies a parametric
one-profile class, then proves uniqueness in the complete minimum-degree
\(t^2\)-linear numerator class.  Neither result is a search in the full space
\(V_d=\operatorname{End}(\operatorname{Sym}^d)\).

On the Hopf sphere use

\[
 t^2+2xy=1,
 \tag{1.1}
\]

where the phase of \(x\) is uniform and \(t\) is uniform on \([-1,1]\).
Fix

\[
 X=1+ax,\qquad a\ne0.
 \tag{1.2}
\]

The **one-profile Hopf class** consists of Laurent presentations

\[
 p_{r,q,R}
 =c\,x^{-r}X^qR(t^2X^2),
 \qquad c\ne0,\quad r\ge1,\quad q\ge0,
 \tag{1.3}
\]

where \(R\in k[z]\) has a zero of order \(s\) at \(z=1\), with \(s\ge r\).
Here \(k\) has characteristic zero; for the nonvanishing statements below
we use the real integral and assume its displayed values are nonzero.

The restrictions in (1.3) have geometric content.  There is one affine
phase profile \(X\), the \(t\)-dependence sees only the Hopf-invariant
combination \(t^2X^2\), and \(x^{-r}\) records the phase winding.  Arbitrary
linear combinations \(\sum_jB_j(x)t^{2j}\) are not claimed to reduce to
this class.  If \(R(z)=\sum_j\rho_jz^j\), then (1.3) is the explicit
subclass of the general ansatz
\[
 \frac{A(x)}{x^r}\sum_jB_j(x)t^{2j}
 \quad\text{with}\quad
 A(x)=cX^q,\qquad B_j(x)=\rho_jX^{2j}.
\tag{1.4}
\]

## 2. Exact polynomiality criterion

Write

\[
 R(z)=(1-z)^sS(z).
 \tag{2.1}
\]

In the sphere coordinate ring,

\[
\begin{aligned}
 1-t^2X^2
 &=1-t^2(1+ax)^2\\
 &=x\left(2y-t^2(2a+a^2x)\right).
\end{aligned}
\tag{2.2}
\]

Thus \(R(t^2X^2)\) is divisible by \(x^s\) modulo (1.1).  The condition
\(s\ge r\) is therefore a transparent sufficient condition for (1.3) to
represent a polynomial on the sphere.  It is also necessary.

> **Proposition 2.1 (polynomiality equals endpoint contact).**  Suppose
> \(a,c\ne0\), \(q\ge0\), and \(R\ne0\).  Then the Laurent presentation
> \(p_{r,q,R}\) in (1.3), temporarily without assuming \(s\ge r\),
> represents an element of the sphere coordinate
> ring if and only if
> \[
> \operatorname{ord}_{z=1}R(z)\ge r.
> \tag{2.3}
> \]

Indeed, eliminate \(t^2\) by (1.1) and work in \(k[x,y]\).  If
\[
 z=(1-2xy)(1+ax)^2,
\]
then
\[
 1-z=xD(x,y),\qquad D(0,y)=2(y-a)\ne0.
\tag{2.4}
\]
Also \(X^q\equiv1\pmod x\) and \(S(z)\equiv S(1)\ne0\pmod x\).  Hence the
numerator in (1.3) has exact \(x\)-adic order \(s\).  Division by \(x^r\)
is polynomial exactly when \(s\ge r\).

Thus polynomiality and the Taylor cancellation are controlled by the same
integer: contact order with the Hopf endpoint \(t^2X^2=1\).

## 3. The classification theorem

Define the phase-integrated polynomial

\[
 H_{m,q,R}(X)
 =X^{qm}\int_0^1R(v^2X^2)^m\,dv
 \tag{3.1}
\]

and put

\[
 C_m=\int_0^1R(v^2)^m\,dv.
 \tag{3.2}
\]

The letter \(v\) in (3.1) is the uniform Hopf-height variable.  It is
renamed \(w=vX\) in the proof.

> **Theorem 3.1 (one-profile classification).**  Let \(R\) have order
> \(s\ge r\) at \(1\), and suppose \(C_m\ne0\) for every \(m\ge1\).  Then
> the following are equivalent:
>
> 1. for every \(m\ge1\),
>    \[
>    [u^{rm}]H_{m,q,R}(1+u)=0,\qquad
>    [u^{rm-1}]H_{m,q,R}(1+u)\ne0;
>    \tag{3.3}
>    \]
> 2. \(q=r\).
>
> When these conditions hold, the adjacent coefficient is exactly
> \[
> [u^{rm-1}]H_{m,r,R}(1+u)=C_m.
> \tag{3.4}
> \]

For the primitive winding \(r=1\), (3.3) is precisely

\[
 [u^m]H_m(1+u)=0,\qquad
 [u^{m-1}]H_m(1+u)=C_m\ne0.
 \tag{3.5}
\]

### Proof

Changing variables \(w=vX\) in (3.1) gives

\[
 H_{m,q,R}(X)
 =X^{qm-1}J_m(X),\qquad
 J_m(X)=\int_0^X R(w^2)^m\,dw.
 \tag{3.6}
\]

Because \(R(w^2)\) has a zero of order \(s\) at \(w=1\),
\(J_m'(w)=R(w^2)^m\) has a zero of order \(sm\) there.  Hence

\[
 J_m(1+u)=C_m+O(u^{sm+1})
          =C_m+O(u^{rm+1}).
 \tag{3.7}
\]

Through degree \(rm\), the Taylor polynomial of (3.6) is therefore the
Taylor polynomial of

\[
 C_m(1+u)^{qm-1}.
\tag{3.8}
\]

If \(q=0\), its coefficient at degree \(rm\) is
\((-1)^{rm}C_m\ne0\), so (3.3) fails.  We may therefore assume \(q\ge1\).
The adjacent coefficient in (3.3) can be nonzero only if
\(qm-1\ge rm-1\), so \(q\ge r\).  The coefficient at degree \(rm\) can
vanish only if \(qm-1<rm\), so \(q\le r\).  Thus \(q=r\).  Conversely,
\((1+u)^{rm-1}\) has zero coefficient at degree \(rm\) and coefficient one
at degree \(rm-1\), which proves (3.3)--(3.4).  \(\square\)

The proof also shows why the two coefficients are genuinely adjacent:
polynomiality supplies at least \(r\) endpoint zeros, and the two coefficient
conditions pin the profile exponent to the same integer \(r\).

The same proof gives more than the adjacent pair.

> **Corollary 3.2 (full lower-jet ladder).**  Under the hypotheses of
> Theorem 3.1, with \(q=r\), one has for
> \(1\le\ell\le rm\)
> \[
> \boxed{
> [u^{rm-\ell}]H_{m,r,R}(1+u)
> =\binom{rm-1}{\ell-1}C_m.}
> \tag{3.9}
> \]
> The coefficient at degree \(rm\) is zero.

Indeed, all these degrees are at most \(rm\), where (3.7) replaces \(J_m\)
by \(C_m\), and
\[
 [u^{rm-\ell}](1+u)^{rm-1}
 =\binom{rm-1}{\ell-1}.
\tag{3.10}
\]
Thus the adjacent phenomenon is the first step of a complete binomial jet
ladder.

## 4. Phase moments and a parametric family

Before imposing (1.1), phase averaging is constant-term extraction in \(x\).
For \(q=r\) and \(1\le\ell\le rm\),

\[
\begin{aligned}
 \int_{\mathrm{phase}}p_{r,r,R}^m
 &=c^ma^{rm}[u^{rm}]H_{m,r,R}(1+u)=0,\\
 \int_{\mathrm{phase}}x^\ell p_{r,r,R}^m
 &=c^ma^{rm-\ell}[u^{rm-\ell}]H_{m,r,R}(1+u)\\
 &=c^ma^{rm-\ell}\binom{rm-1}{\ell-1}C_m.
\end{aligned}
\tag{4.1}
\]

Thus every \(R\) satisfying the hypotheses gives an all-order parametric
family, and every allowed positive phase multiplier detects it.  A simple
real sufficient condition for \(C_m\ne0\) is

\[
 R(z)=(1-z)^sS(z),\qquad S(z)>0\quad(0\le z\le1).
\tag{4.2}
\]

In particular \(R(z)=(1-z)^s\) gives

\[
 C_m=\int_0^1(1-v^2)^{sm}\,dv
 =\frac{2^{sm}(sm)!}{(2sm+1)!!}\ne0.
\tag{4.3}
\]

This family is larger than taking powers of the primitive example: an
arbitrary factor \(S(t^2X^2)\) satisfying (4.2) may be inserted without
destroying the adjacent jet.

## 5. Degree-minimal uniqueness and Long's polynomial

Now take primitive winding \(r=1\).  The theorem forces \(q=1\), while
polynomiality forces \(R(1)=0\).  Among nonzero \(R\), the minimum possible
degree is one, and then necessarily

\[
 R(z)=\lambda(1-z).
\tag{5.1}
\]

Absorbing \(\lambda\) into \(c\) gives the unique degree-minimal
one-profile lift

\[
 p_{a,c}
 =c\,\frac{1+ax}{x}
   \left(1-t^2(1+ax)^2\right).
\tag{5.2}
\]

Using (2.2), this is the five-term sphere polynomial

\[
\boxed{
 p_{a,c}
 =c(1+ax)\left(2y-(2a+a^2x)t^2\right).
 }
\tag{5.3}
\]

The torus change \(x'=ax,\ y'=a^{-1}y\) preserves (1.1), so all
\(a\ne0\) give one orbit, up to the overall scalar \(c\).  At
\(a=1,c=\tfrac12\),

\[
 p=y+xy-t^2-\frac32xt^2-\frac12x^2t^2,
\tag{5.4}
\]

which is Long's five-term polynomial in Hopf coordinates.  Its
phase-integrated function is

\[
 H_m(X)=X^m\int_0^1(1-v^2X^2)^m\,dv,
\tag{5.5}
\]

and (3.5) recovers the adjacent-jet argument used by the two-pair lift.

> **Corollary 5.1 (minimal uniqueness).**  Up to the sphere torus
> \(x\mapsto ax,\ y\mapsto a^{-1}y\) and nonzero scaling, Long's polynomial
> is the unique member of the one-profile Hopf class with primitive winding
> and linear \(R\).  Equivalently, it is the unique member of minimum
> \(R\)-degree in that class having the adjacent-jet phenomenon.

This is a uniqueness theorem for the Long/Hopf mechanism inside (1.3).
It is not global uniqueness among all low-degree sphere polynomials.  In
particular, it makes no claim about general \(V_4\), about ansatzes with
several independent phase profiles, or about cancellations not explained
by endpoint multiplicity.

## 6. Minimal-degree \(t^2\)-linear classification

The factorization in the preliminary ansatz
\[
 \frac{A(x)}x\left(B_0(x)+B_1(x)t^2\right)
\]
is redundant: only the products
\[
 C(x)=A(x)B_0(x),\qquad D(x)=A(x)B_1(x)
\tag{6.1}
\]
enter the sphere polynomial and its moments.  The invariant low-degree
question should therefore bound \(C,D\), not a chosen factorization.

Consider the minimum numerator degrees containing Long's polynomial:

\[
 p=\frac{C(x)+D(x)t^2}{x},
 \qquad \deg C\le1,\quad\deg D\le3.
\tag{6.2}
\]

Put
\[
 K_m(x)=\int_0^1(C(x)+v^2D(x))^m\,dv.
\tag{6.3}
\]

> **Theorem 6.1 (three-jet rigidity in the minimal linear class).**  For
> (6.2), the following are equivalent, up to multiplying \(p\) by a nonzero
> scalar:
>
> 1. \(p\) is polynomial on the sphere and, for every \(m\ge1\),
>    \[
>    [x^m]K_m(x)=0,\qquad [x^{m-1}]K_m(x)\ne0;
>    \tag{6.4}
>    \]
> 2. \(p\) is polynomial,
>    \[
>    [x]K_1=[x^2]K_2=[x^3]K_3=0,
>    \qquad [x]K_2\ne0;
>    \tag{6.5}
>    \]
> 3. for some \(a\ne0\),
>    \[
>    C(x)=1+ax,\qquad D(x)=-(1+ax)^3.
>    \tag{6.6}
>    \]
>
> Consequently Long's five-term polynomial is the unique polynomial with
> the all-order adjacent-jet phenomenon in the full \(t^2\)-linear class
> (6.2), up to the sphere torus and nonzero scaling.

### Proof

Eliminating \(t^2\) gives
\[
 C(x)+D(x)(1-2xy).
\tag{6.7}
\]
Thus (6.2) is polynomial exactly when \(C(0)+D(0)=0\).  Moreover,
the adjacent condition \([x]K_2\ne0\) forces \(C(0)\ne0\), since otherwise
both \(C,D\) are divisible by \(x\) and \(K_2\) is divisible by \(x^2\).
After scaling, write
\[
 C=1+ax,\qquad
 D=-1+b_1x+b_2x^2+b_3x^3.
\tag{6.8}
\]

Direct integration gives
\[
\begin{aligned}
 [x]K_1
 &=\frac{3a+b_1}{3},\\
 [x^2]K_2
 &=\frac{
 15a^2+10ab_1+3b_1^2+4b_2}{15},\\
 [x^3]K_3
 &=\frac1{35}\bigl(
 35a^3+35a^2b_1+21ab_1^2+28ab_2\\
 &\hspace{31mm}
 {}+5b_1^3+12b_1b_2+8b_3\bigr).
\end{aligned}
\tag{6.9}
\]
They form a triangular system.  Successive vanishing yields
\[
 b_1=-3a,\qquad
 [x^2]K_2=\frac4{15}(b_2+3a^2),\qquad
 [x^3]K_3=\frac8{35}(b_3+a^3),
\tag{6.10}
\]
so
\[
 b_1=-3a,\qquad b_2=-3a^2,\qquad b_3=-a^3.
\tag{6.11}
\]
Therefore \(D=-(1+ax)^3=-C^3\).  Under the first pure equation,
\[
 [x]K_2=\frac{8a}{15},
\tag{6.12}
\]
so the adjacent condition in (6.5) is exactly \(a\ne0\).

Now (6.6) is the primitive profile (5.2), and Theorem 3.1 proves all the
identities in (6.4).  Finally \(x\mapsto ax,\ y\mapsto a^{-1}y\) removes
the remaining nonzero parameter.  \(\square\)

Theorem 6.1 closes the first multi-coefficient island beyond the composite
ansatz (1.3).  Its finite-jet certificate is conceptual and triangular; it
does not enumerate supports or points of \(V_4\).

The triangularity is uniform in the degree.

> **Proposition 6.2 (universal triangular reconstruction).**  Normalize
> \[
> C(x)=1+\sum_{j\ge1}a_jx^j,\qquad
> D(x)=-1+\sum_{j\ge1}b_jx^j.
> \tag{6.13}
> \]
> For every \(m\ge1\), the pure jet \([x^m]K_m\) depends only on
> \(a_1,\ldots,a_m,b_1,\ldots,b_m\), and its coefficient of the new
> variable \(b_m\) is
> \[
> \tau_m
> =m\int_0^1v^2(1-v^2)^{m-1}\,dv
> =\frac{2^{m-1}m!}{(2m+1)!!}\ne0.
> \tag{6.14}
> \]
> Consequently the equations
> \([x^m]K_m=0\), \(1\le m\le N\), solve
> \(b_1,\ldots,b_N\) uniquely and successively in terms of
> \(a_1,\ldots,a_N\).

Indeed, terms of degree greater than \(m\) cannot contribute to
\([x^m]\).  Perturbing \(D\) by \(\epsilon x^m\), the coefficient linear in
\(\epsilon x^m\) is
\[
 m\int_0^1v^2(C(0)+v^2D(0))^{m-1}\,dv
 =m\int_0^1v^2(1-v^2)^{m-1}\,dv,
\tag{6.15}
\]
which evaluates to (6.14).  No term containing \(b_m^2\) can have degree
\(m\).  This proves the claim.

For the rectangle \((\deg C,\deg D)\le(d,d+2)\), the first \(d+2\) pure
jets therefore reconstruct all of \(D\).  After normalizing \(a_1=1\), the
entire remaining problem is an elimination problem in the \(d-1\) weighted
parameters \(a_2,\ldots,a_d\). Proposition 6.2 gives this reduction; a
classification of those higher-degree residual ideals is outside the
retained proof above.

## Verification scope

The written classification above is the all-order argument. The [rank-one Gaussian checker](../scripts/verify_three_real_gmc_rank_one_classification.py) verifies the associated exact coefficient identities. Finite checks do not replace the universal proof.
