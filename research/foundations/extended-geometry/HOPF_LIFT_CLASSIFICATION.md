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

## 7. Real profiles with mixed height powers

**Status and scope.** The results in this section are internal deductions,
with no independent mathematical review, formalization or novelty claim.
They extend the endpoint argument to real polynomials in the signed height.
The real-coefficient hypothesis is essential to the nonvanishing argument
used here; no classification for complex profile coefficients is asserted.
The only proof inputs are the sphere relation (1.1), constant-term phase
averaging, polynomial integration and the radial contraction identity in
[the two-pair note, Sections 2 and 4](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
Neither a bounded search nor a new external theorem is used.

Let
\[
 A=\mathbb R[x,y,t]/(t^2+2xy-1),\qquad X=1+x,
\]
and write \(\Lambda\) for the normalized angular functional: substitute
\(y=(1-t^2)/(2x)\), take the constant term in \(x\), and then integrate
\(t\) against \(dt/2\) on \([-1,1]\). Fix
\[
 r\geq1,\quad q\geq1,\quad \gamma\in\mathbb R\setminus\{0\},\quad
 0\ne P\in\mathbb R[w],\qquad (1-w^2)^r\mid P(w),
\]
where \(r,q\) are integers, and set
\[
 p=\gamma x^{-r}X^qP(tX),\qquad
 C_m=\frac12\int_{-1}^1P(w)^m\,dw.
 \tag{7.1}
\]
This is a single affine phase profile, allowing arbitrary mixtures of
height powers inside \(P\). It represents a polynomial in \(A\): if
\(P(w)=(1-w^2)^rS(w)\), then
\[
 p=\gamma X^q\bigl(2y-(2+x)t^2\bigr)^r S(tX).
 \tag{7.2}
\]
In particular \(C_{2j}>0\) for every \(j\geq1\).

> **Theorem 7.1 (real signed-height profile classification).** For (7.1),
> \[
>  \Lambda(p^m)=\gamma^m\binom{qm-1}{rm}C_m
>  \qquad(m\geq1).
>  \tag{7.3}
> \]
> Thus all pure moments vanish exactly when \(q\leq r\). If \(q=r\),
> \[
>  \Lambda(xp^m)=\gamma^m C_m,
>  \tag{7.4}
> \]
> which is nonzero at every even \(m\). If \(q<r\), then for every fixed
> polynomial marker \(B\in A\),
> \(\Lambda(Bp^m)=0\) for all sufficiently large \(m\).
> The cutoff may depend on \(B\), as specified in Lemma 7.2.

To prove (7.3)--(7.4), use the polynomial
\[
 H_m(X)=X^{qm-1}J_m(X),\qquad
 J_m(X)=\frac12\int_{-X}^{X}P(w)^m\,dw.
 \tag{7.5}
\]
Substitution \(w=tX\) and phase extraction give
\[
 \Lambda(p^m)=\gamma^m[x^{rm}]H_m(1+x),\qquad
 \Lambda(xp^m)=\gamma^m[x^{rm-1}]H_m(1+x).
 \tag{7.6}
\]
Because \(P\) has zeros of order at least \(r\) at both endpoints,
\[
 J_m'(X)=\frac{P(X)^m+P(-X)^m}{2}
\]
is divisible by \((X-1)^{rm}\). Consequently
\[
 H_m(X)\equiv C_mX^{qm-1}\pmod{(X-1)^{rm+1}}.
 \tag{7.7}
\]
This proves (7.3). When \(q\leq r\), the binomial coefficient there is
zero. When \(q>r\), it is positive, and the even moments are nonzero
because \(\gamma^{2j}C_{2j}>0\). For \(q=r\), the coefficient of degree
\(rm-1\) in (7.7) proves (7.4). No assertion that odd marked moments
are nonzero is needed.

> **Lemma 7.2 (all fixed polynomial markers when \(q<r\)).** Suppose
> \(1\leq q<r\). For nonnegative integers \(a,b,h\),
> \[
>  \Lambda(x^ay^bt^h p^m)=0
>  \quad\text{if}\quad
>  m\geq\max\left\{1,
>       \left\lceil\frac{h+2b+1}{q}\right\rceil,
>       \left\lceil\frac{a}{r-q}\right\rceil\right\}.
>  \tag{7.8}
> \]
> Therefore \(m\geq2D+1\) suffices for every marker having a polynomial
> representative of total degree at most \(D\) in \(x,y,t\).

Here is the endpoint estimate needed for arbitrary markers. Put
\[
\begin{aligned}
 N&=rm-a+b,& e&=qm-h-2b-1,\\
 K_m(X)&=\frac12\int_{-X}^{X}
       w^h(X^2-w^2)^bP(w)^m\,dw,\\
 K_m^0(X)&=\frac12\int_{-1}^{1}
       w^h(X^2-w^2)^bP(w)^m\,dw.
\end{aligned}
 \tag{7.9}
\]
The superscript zero means fixed integration endpoints; it is not a
constant in \(X\). It is a polynomial of degree at most \(2b\).
After substituting for \(y\) and changing height variable,
\[
 \Lambda(x^ay^bt^h p^m)
   =\gamma^m2^{-b}[x^N](1+x)^eK_m(1+x).
 \tag{7.10}
\]
The endpoint difference has the stronger divisibility
\[
 K_m(X)-K_m^0(X)\quad\text{is divisible by}\quad
 (X-1)^{rm+b+1}.
 \tag{7.11}
\]
Indeed, split the difference into the integrals from \(1\) to \(X\)
and from \(-X\) to \(-1\). In the first, put
\(X=1+u\), \(w=1+v\). The factor \(P(1+v)^m\) is divisible by
\(v^{rm}\), while
\(X^2-w^2=(u-v)(2+u+v)\). Each term integrated from \(v=0\) to
\(v=u\) therefore has a factor \(u^{rm+b+1}\). The second integral
has the same estimate after putting \(w=-1-v\), using the zero of
\(P\) at \(-1\). These are polynomial integrations, so the estimate
is an exact divisibility statement.

Under (7.8), \(e\geq0\), \(N\geq0\), and \(N\leq rm+b\).
Thus (7.11) permits replacing \(K_m\) by \(K_m^0\) in (7.10).
The resulting polynomial has degree at most
\[
 e+2b=qm-h-1,
\]
whereas
\[
 N-(qm-h-1)=(r-q)m-a+b+h+1>0
\]
by (7.8). Its coefficient of degree \(N\) is zero, proving the lemma.
For \(a+b+h\leq D\), the simpler bound \(m\geq2D+1\) implies
both inequalities in (7.8). A polynomial marker is a finite linear
combination of these monomials. This completes Theorem 7.1.

### Minimal balanced degree and the original Image condition

Let \(n=\deg P\), \(s=n-2r=\deg S\), and write
\(S(w)=\sum_{j=0}^s s_jw^j\), where \(s_s\ne0\). Use the original
bilinears \(R,Z,W,T\), put \(L=R+Z\), and use
\(D=2WR^2-2RT^2-ZT^2\), as in the two-pair note.

> **Proposition 7.3 (primitive homogenization and degree parity).**
> The least balanced degree of a homogeneous lift of (7.1) is
> \[
>  d=2n+q-r.
>  \tag{7.12}
> \]
> Its lift in this degree is
> \[
>  F=\gamma L^qD^r
>       \sum_{j=0}^s s_jR^{2(s-j)}T^jL^j.
>  \tag{7.13}
> \]
> It is \(R\)-primitive. For \(q=r\), this is an Image failure
> witness of even balanced degree \(d=2n\), with marker \(Z\) and
> nonzero marked contractions at every even power. For \(q<r\), it
> obeys the eventual Image conclusion for every fixed original
> polynomial multiplier, including unbalanced multipliers. For a
> multiplier of bidegree \(u,v\), \(m\geq2v+1\) suffices.

Each term of (7.13) has balanced degree
\(q+3r+2s=2n+q-r\), and setting \(R=1\) recovers (7.2).
The highest total-degree term of (7.2), before sphere reduction, is
\[
 \gamma(-1)^rs_s x^{q+n-r}t^n.
\]
Its image is nonzero in the associated graded sphere ring
\(\mathbb R[x,y,t]/(t^2+2xy)\), which is a domain. This proves the
minimal degree assertion: the highest term cannot be removed by the
sphere relation. Equivalently, reduction using \(t^2=1-2xy\) gives
a nonzero term of that same degree in the unique representative of
height degree at most one.

At \((\xi_1,\xi_2,z_1,z_2)=(1,1,1,-1)\), one has
\(R=0\), \(T=2\), \(L=-1\), \(D=4\). Only the \(j=s\) term
in (7.13) remains, and its value is
\(\gamma s_s(-1)^{q+s}2^{s+2r}\ne0\). Thus \(R\nmid F\).
For \(q=r\), radial homogeneity and (7.4) give
\[
 \mathcal E_2(ZF^m)=(dm+2)!\,\gamma^m C_m,
 \tag{7.14}
\]
which is nonzero at every even \(m\), while all pure contractions
vanish. For \(q<r\), Lemma 7.2 first gives eventual zero contraction
against every balanced marker: such a marker is a polynomial in the
four bilinears, and its sphere restriction belongs to \(A\).

To check arbitrary original multipliers, split a multiplier \(Q\) into
its finitely many bihomogeneous parts, of bidegrees \((u,v)\). If
\(u>v\), then \(\mathcal E_2(QF^m)=0\) automatically, since the
number of differentiations exceeds the coordinate degree. If \(u\leq v\),
this contraction is a homogeneous polynomial in \(z\) of degree
\(v-u\). For each multi-index \(\alpha\) with \(|\alpha|=v-u\),
\[
 \partial_z^\alpha\mathcal E_2(QF^m)
     =\mathcal E_2(\xi^\alpha QF^m).
 \tag{7.15}
\]
The marker on the right is balanced. It eventually gives zero by the
preceding paragraph. There are only finitely many such scalar tests;
once all vanish, the homogeneous output polynomial is zero. Taking
the maximum cutoff over the finitely many parts of \(Q\) proves the
full eventual Image conclusion. More explicitly, the balanced markers
in (7.15) have balanced degree \(v\), so they are sums of products of
\(v\) bilinears. Their sphere restrictions have representatives of
total degree at most \(v\). Lemma 7.2 therefore gives the sufficient
cutoff \(m\geq2v+1\).

The classification therefore rules out primitive odd-degree Image
failures within (7.1): the only failing case is \(q=r\), whose least
balanced degree is even. Homogenization in a higher degree introduces
\(R\) padding and loses primitivity. This obstruction allows arbitrary
mixed height powers in a single real profile, but makes no claim about
multiple phase profiles or the full primitive degree-five space.


## Verification scope

The written classifications are the all-order arguments. The preexisting
[rank-one Gaussian checker](../scripts/verify_three_real_gmc_rank_one_classification.py)
is associated with the earlier classification material; it is not claimed
to verify Section 7. That section rests on its displayed elementary proof.
Its endpoint estimates, degree calculation and arbitrary-marker bridge
received focused internal AI proof checks. No retained computational
checker, independent replay, formal verification or external mathematical
review is claimed for Section 7. Finite checks do not replace the
universal proof.
