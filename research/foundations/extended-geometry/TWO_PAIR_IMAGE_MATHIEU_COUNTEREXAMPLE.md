# A bidegree-\((4,4)\) counterexample in two contraction pairs

## 1. Statement

Work over a characteristic-zero field and use the contraction pairs

\[
 (\xi _1,z_1),\qquad(\xi _2,z_2).
\]

Put

\[
 \mathcal M_2
 =(\partial_{z_1}-\xi _1)k[\xi _1,\xi _2,z_1,z_2]
  +(\partial_{z_2}-\xi _2)k[\xi _1,\xi _2,z_1,z_2],
 \tag{1.1}
\]

and define the four bilinears

\[
\begin{aligned}
 R&=\xi _1z_1+\xi _2z_2,&
 Z&=\xi _1z_2,\\
 W&=2\xi _2z_1,&
 T&=\xi _1z_1-\xi _2z_2.
\end{aligned}
\tag{1.2}
\]

They satisfy the rank-one-quadric identity

\[
 T^2=R^2-2ZW.
\tag{1.3}
\]

Set

\[
 \boxed{
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),
 \qquad Q=Z.}
\tag{1.4}
\]

> **Theorem 1.1.** For every \(m\geq1\),
> \[
>  \boxed{
>  \mathcal E_2(F^m)=0,\qquad
>  \mathcal E_2(QF^m)
>   =\frac{(4m+2)!\,m!}{(2m+1)!!}\ne0.}
> \tag{1.5}
> \]
> Consequently \(F^m\in\mathcal M_2\) but
> \(QF^m\notin\mathcal M_2\) for every \(m\geq1\).
> Thus \(\operatorname{SIC}(2)\) is false.

The one-pair Special Image Conjecture is known, so the minimum failing pair
dimension is exactly two.  The witness \(F\) is bihomogeneous of bidegree
\((4,4)\), has ordinary total degree eight, and has sixteen expanded terms.
The multiplier \(Q\) has bidegree \((1,1)\).

In the monomial bases
\(\xi _1^i\xi _2^{4-i}\) and \(z_1^jz_2^{4-j}\), the coefficient matrix of
\(F\), with \(i,j=0,\ldots,4\), is
\[
 \begin{pmatrix}
 -1&2&0&0&0\\
 -3/2&2&6&0&0\\
 -1/2&3/2&6&6&0\\
 0&1&3/2&2&2\\
 0&0&-1/2&-3/2&-1
 \end{pmatrix},
\qquad \det=48.
\tag{1.6}
\]
In particular \(F\) has tensor rank five, so this is a genuinely
nonseparable two-pair witness rather than a constant-coefficient GVC point
on the rank-one Segre cone.

## 2. Contraction as circular Gaussian expectation

For a monomial, let

\[
 \mathcal E_2(\xi^\alpha z^\beta)
 =\partial_z^\alpha z^\beta.
\tag{2.1}
\]

Zhao's image-kernel identity gives

\[
 \mathcal M_2=\ker\mathcal E_2.
\tag{2.2}
\]

Both \(F^m\) and \(QF^m\) have equal dual and coordinate degrees.  Their
contractions are therefore scalars.  If \(G_1,G_2\) are independent
standard circular complex Gaussians, the Wick rule gives, for every
balanced polynomial \(H\),

\[
 \mathcal E_2(H)
 =\mathbb E\,
 H(\overline G_1,\overline G_2,G_1,G_2).
\tag{2.3}
\]

All coefficients in (1.4) are rational, and both sides of (1.5) are
rational scalars obtained by formal differentiation.  It therefore
suffices to prove the identities over \(\mathbb C\) using (2.3); the
resulting rational identities then extend to every characteristic-zero
field.

Write

\[
 S=|G_1|^2+|G_2|^2,\qquad
 (G_1,G_2)=\sqrt S\,(U_1,U_2),
\tag{2.4}
\]

where \(U\) is uniform on the unit sphere in \(\mathbb C^2\) and is
independent of \(S\).  The variable \(S\) has the gamma distribution of
shape two, hence

\[
 \mathbb E(S^n)=(n+1)!.
\tag{2.5}
\]

On the sphere put

\[
 x=\overline U_1U_2,\qquad
 y=2\overline U_2U_1,\qquad
 t=|U_1|^2-|U_2|^2.
\tag{2.6}
\]

Then \(t^2+2xy=1\), and radial homogeneity gives

\[
 F=S^4p,\qquad Q=Sx,
\quad
 p=(1+x)\left(y-\frac12(2+x)t^2\right).
\tag{2.7}
\]

It remains to calculate two angular moments.

## 3. The all-order angular identity

The Hopf coordinates may be chosen so that \(t\) is uniform on
\([-1,1]\), the phase of \(x\) is uniform, and

\[
 xy=\frac{1-t^2}{2}.
\tag{3.1}
\]

Phase averaging is therefore constant-term extraction after substituting
\(y=(1-t^2)/(2x)\).  Equation (2.7) becomes

\[
 p=\frac{1+x}{2x}
 \left(1-t^2(1+x)^2\right).
\tag{3.2}
\]

Since the integrand is even in \(t\), define

\[
 H_m(X)
 =X^m\int_0^1(1-s^2X^2)^m\,ds.
\tag{3.3}
\]

The phase coefficient and the \(t\)-integral give

\[
\begin{aligned}
 \mathbb E_U(p^m)
 &=2^{-m}[u^m]H_m(1+u),\\
 \mathbb E_U(xp^m)
 &=2^{-m}[u^{m-1}]H_m(1+u).
\end{aligned}
\tag{3.4}
\]

After the change of variables \(v=sX\),

\[
 H_m(X)
 =X^{m-1}J_m(X),
\qquad
 J_m(X)=\int_0^X(1-v^2)^m\,dv.
\tag{3.5}
\]

Now \(J_m'(X)=(1-X^2)^m\) has a zero of order \(m\) at \(X=1\).
Consequently the Taylor coefficients through order \(m\) of
\(H_m(X)\) at \(X=1\) are the same as those of
\(J_m(1)X^{m-1}\).  Equation (3.4) immediately yields

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=2^{-m}J_m(1).
\tag{3.6}
\]

The remaining beta integral is

\[
\begin{aligned}
 J_m(1)
 &=\int_0^1(1-v^2)^m\,dv\\
 &=\frac12B\left(\frac12,m+1\right)
 =\frac{2^m m!}{(2m+1)!!}.
\end{aligned}
\tag{3.7}
\]

Thus

\[
 \mathbb E_U(p^m)=0,\qquad
 \mathbb E_U(xp^m)=\frac{m!}{(2m+1)!!}.
\tag{3.8}
\]

The endpoint-zero mechanism in this calculation is classified in
[`HOPF_LIFT_CLASSIFICATION.md`](HOPF_LIFT_CLASSIFICATION.md).  Inside the
one-profile class treated there, the phase winding forces the profile
exponent, and the present five-term polynomial is the unique member of
minimum profile degree, up to the sphere torus and overall scaling.

Combining (2.3), (2.5), (2.7), and (3.8) gives

\[
\begin{aligned}
 \mathcal E_2(F^m)
 &=(4m+1)!\,\mathbb E_U(p^m)=0,\\
 \mathcal E_2(QF^m)
 &=(4m+2)!\,\mathbb E_U(xp^m)
 =\frac{(4m+2)!\,m!}{(2m+1)!!},
\end{aligned}
\]

which proves Theorem 1.1.

## 4. A second proof by direct coefficient extraction

This section gives an algebraic proof independent of the circular-Gaussian
and Hopf-coordinate argument above.  In particular, it uses neither a
Gaussian measure nor radial decomposition nor the Taylor-jet argument in
(3.5)--(3.6).

Let \(I:\mathbb Q[t]\to\mathbb Q\) be the linear functional
\[
 I(t^j)=\frac1{j+1}.
\]
It satisfies the formal fundamental theorem
\(I(P')=P(1)-P(0)\), first on monomials and hence on every polynomial.
For \(a\geq1\), applying this to
\(P=t^a(1-t)^{b+1}\), whose two endpoint values vanish, gives
\[
 (b+1)I\!\left(t^a(1-t)^b\right)
 =aI\!\left(t^{a-1}(1-t)^{b+1}\right).
\]
Induction on \(a\), starting from
\(I((1-t)^b)=1/(b+1)\), proves the entirely algebraic beta identity
\[
 I\!\left(t^a(1-t)^b\right)=\frac{a!\,b!}{(a+b+1)!}.
 \tag{4.1}
\]
We continue to write \(I(P)=\int_0^1P(t)\,dt\), but no analytic integration
is involved.  For a polynomial \(H\) bihomogeneous of bidegree \((n,n)\),
direct monomial inspection now gives the coefficient-extraction identity
\[
 \boxed{\displaystyle
 \mathcal E_2(H)
  =(n+1)!\operatorname {CT}_u
    \int_0^1H\left(1,u,t,\frac{1-t}{u}\right)\,dt.}
 \tag{4.2}
\]
Indeed, for
\[
 H=\xi _1^a\xi _2^{n-a}z_1^bz_2^{n-b},
\]
the constant term in \(u\) is zero unless \(a=b\); when \(a=b\), (4.1)
turns the right-hand side of (4.2) into \(a!(n-a)!\), exactly its
contraction.

Under the substitution in (4.2), put
\[
 v=2t-1,\qquad x=Z=\frac{1-v}{2u}.
\]
Then
\[
 R=1,\qquad T=v,\qquad ZW=\frac{1-v^2}{2},
\]
and hence
\[
 F=\frac{1+x}{2x}\left(1-v^2(1+x)^2\right),\qquad Q=x.
 \tag{4.3}
\]
Here is the precise constant-term change.  Before applying \(I\), work over
the coefficient field \(\mathbb Q(v)\).  If \(c=(1-v)/2\), the substitution
\(x=c/u\) sends \(u^j\) to \(c^jx^{-j}\).  Thus exponent zero is preserved
and no nonzero exponent can contribute to it:
\(\operatorname {CT}_u=\operatorname {CT}_x\).  After simplification,
(4.3) lies in \(\mathbb Q[v][x,x^{-1}]\), so the formal functional \(I\)
applies coefficientwise.  The expression is even in \(v\).  Expanding it
in (4.2) is especially short:
\[
\begin{aligned}
 \operatorname {CT}_x F^m
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \operatorname {CT}_x\!
    \left(x^{-m}(1+x)^{m+2k}\right)v^{2k}\\
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \binom{m+2k}{m}v^{2k},\\
 \operatorname {CT}_x(xF^m)
 &=2^{-m}\sum_{k=0}^m(-1)^k\binom mk
   \binom{m+2k}{m-1}v^{2k}.
\end{aligned}
\]
Moreover the formal fundamental theorem already used above gives
\[
 I\!\left((2t-1)^{2k}\right)
 =I\!\left(\left(\frac{(2t-1)^{2k+1}}
 {2(2k+1)}\right)'\right)
 =\frac1{2k+1}.
\]
Since \(F^m\) and \(QF^m\) have bidegrees \((4m,4m)\) and
\((4m+1,4m+1)\), respectively, (4.2) now gives, for every \(m\geq1\),
\[
\begin{aligned}
 \frac{\mathcal E_2(F^m)}{(4m+1)!}
  &=2^{-m}\sum_{k=0}^m
    \frac{(-1)^k\binom mk\binom{m+2k}{m}}{2k+1},
    \tag{4.4}\\
 \frac{\mathcal E_2(QF^m)}{(4m+2)!}
  &=2^{-m}\sum_{k=0}^m
    \frac{(-1)^k\binom mk\binom{m+2k}{m-1}}{2k+1}.
    \tag{4.5}
\end{aligned}
\]
These formulas came directly from formal contraction, not from Wick
expectation or phase averaging.

We now evaluate both sums by finite differences.  If
\(\Delta A(X)=A(X+1)-A(X)\), induction on \(m\) gives
\[
 \sum_{k=0}^m(-1)^k\binom mk A(k)=(-1)^m\Delta^mA(0).
\]
Each application of \(\Delta\) lowers the degree of a nonconstant
polynomial by one.  Consequently, for every polynomial \(A(X)\) of degree
less than \(m\),
\[
 \sum_{k=0}^m(-1)^k\binom mk A(k)=0.
 \tag{4.6}
\]
For the pure sum,
\[
 \frac{\binom{m+2X}{m}}{2X+1}
   =\frac1{m!}\prod_{j=2}^{m}(2X+j)
 \tag{4.7}
\]
is a polynomial of degree \(m-1\).  Equations (4.4), (4.6), and (4.7)
immediately prove \(\mathcal E_2(F^m)=0\).

The same observation gives a useful general residue principle.  If
\(\deg A<m\), polynomial division at \(X=-1/2\) gives
\[
 A(X)=A(-1/2)+(2X+1)D(X),\qquad \deg D<m.
\]
Consequently
\[
 \boxed{\displaystyle
 \sum_{k=0}^m(-1)^k\binom mk\frac{A(k)}{2k+1}
 =A(-1/2)
   \sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1}.}
 \tag{4.8}
\]
Thus the alternating quotient transform on all polynomials of degree
less than \(m\) has rank one: it is evaluation at \(-1/2\), multiplied by
one universal scalar.

There is a useful higher-degree-denominator form of the same argument.  Fix
a polynomial \(L\) such that \(L(k)\ne0\) for \(0\leq k\leq m\), and write
\[
 {\cal T}_{m,L}(A)
 =\sum_{k=0}^m(-1)^k\binom mk\frac{A(k)}{L(k)}.
\]
If
\[
 A=R+LD,\qquad \deg D<m,
\]
then (4.6), applied to \(D\), gives the exact invariance
\[
 \boxed{{\cal T}_{m,L}(A)={\cal T}_{m,L}(R).}
 \tag{4.8a}
\]
In particular, if \(L\) is normalized to be monic of degree \(r\geq1\)
and \(\deg A<m+r\), polynomial division gives \(\deg D<m\) and
\[
 {\cal T}_{m,L}(A)={\cal T}_{m,L}(A\bmod L).
 \tag{4.8b}
\]
Thus, on this degree range, the transform factors through the
\(r\)-dimensional quotient \(\mathbb Q[X]/(L)\): it sees only the remainder
modulo \(L\).  (As a single scalar-valued functional its linear-map rank is
of course at most one.)  Formula (4.8) is the case \(r=1\).  For a repeated
denominator
\(L=(X+\tfrac12)^r\), the remainder is equivalently determined by the
first \(r\) Taylor coefficients of \(A\) at \(-1/2\).  This jet
interpretation is a structural extension of the finite-sum certificate;
by itself it does not assert the existence of further SIC witnesses.

The repeated-pole case admits an explicit all-order evaluation.  For
\(s\geq0\), put
\[
 B_{m,s}=\sum_{k=0}^m
 \frac{(-1)^k\binom mk}{(2k+1)^s}.
\]
Thus \(B_{m,0}=0\) for \(m\geq1\), and the sum \(B_m\) in (4.10) is
\(B_{m,1}\).  For \(m,s\geq1\), the termwise identity
\[
 \frac{2m+1}{(2k+1)^s}
 =\frac1{(2k+1)^{s-1}}
  +\frac{2(m-k)}{(2k+1)^s}
\]
and \((m-k)\binom mk=m\binom{m-1}k\) give the triangular recurrence
\[
 \boxed{(2m+1)B_{m,s}=B_{m,s-1}+2mB_{m-1,s}.}
 \tag{4.8c}
\]
The boundary values are \(B_{0,s}=1\) and \(B_{m,0}=0\) for \(m\geq1\).
Hence (4.8c) evaluates every repeated-pole residue using only rational
arithmetic.

More explicitly, suppose \(r\geq1\), \(\deg A<m+r\), and expand the
remainder in the local coordinate \(\lambda=2X+1\):
\[
 A(X)=\sum_{j=0}^{r-1}c_j\lambda^j+\lambda^rD(X),
 \qquad \deg D<m.
\]
Then (4.6) gives the exact jet formula
\[
 \boxed{\displaystyle
 \sum_{k=0}^m(-1)^k\binom mk
 \frac{A(k)}{(2k+1)^r}
 =\sum_{j=0}^{r-1}c_jB_{m,r-j}
 =\sum_{j=0}^{r-1}
 \frac{A^{(j)}(-1/2)}{2^j j!}\,B_{m,r-j}.}
 \tag{4.8d}
\]
There is also a compact closed description.  With an auxiliary variable
\(a\), the same first-order recurrence as (4.11) gives
\[
 G_m(a):=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+a}
 =\frac{2^m m!}{\prod_{\ell=0}^m(a+2\ell)}.
\]
Formal differentiation in \(\mathbb Q(a)\) yields
\[
 B_{m,s}
 =\frac{(-1)^{s-1}}{(s-1)!}
 \left.\frac{d^{\,s-1}}{da^{s-1}}G_m(a)\right|_{a=1}
 \qquad(s\geq1).
 \tag{4.8e}
\]
If \(h_q\) denotes the complete homogeneous symmetric polynomial of degree
\(q\), expansion at \(a=1\) makes this still more concrete:
\[
 \boxed{\displaystyle
 B_{m,s}
 =\frac{2^m m!}{(2m+1)!!}\,
 h_{s-1}\left(1,\frac13,\ldots,\frac1{2m+1}\right)>0.}
 \tag{4.8f}
\]
Indeed,
\[
 \frac{G_m(1+y)}{G_m(1)}
 =\prod_{\ell=0}^m
  \left(1+\frac{y}{2\ell+1}\right)^{-1},
\]
whose coefficient of \(y^q\) is
\((-1)^qh_q(1,1/3,\ldots,1/(2m+1))\).  Consequently the jet in (4.8d)
cannot vanish when its local coefficients \(c_j\) are all nonnegative and
at least one is positive.  Equations (4.8c)--(4.8f) make the surviving
\(r\)-jet effective; realizing such a jet by an admissible SIC multiplier
remains a separate problem.

For the mixed sum, take
\[
 A_m(X)=\binom{m+2X}{m-1}
       =\frac1{(m-1)!}\prod_{j=2}^{m}(2X+j).
 \tag{4.9}
\]
It has degree \(m-1\) and \(A_m(-1/2)=1\).  Formula (4.8) therefore
reduces the mixed sum to
\[
 B_m=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1}.
 \tag{4.10}
\]
This last sum has a purely finite recurrence.  Since
\[
 m\binom{m-1}{k}=(m-k)\binom mk,
\]
termwise subtraction gives
\[
\begin{aligned}
 (2m+1)B_m-2mB_{m-1}
 &=\sum_{k=0}^m(-1)^k\binom mk\\
 &=0.
\end{aligned}
\tag{4.11}
\]
Together with \(B_0=1\), this proves
\[
 B_m=\prod_{j=1}^m\frac{2j}{2j+1}
     =\frac{2^m m!}{(2m+1)!!}.
 \tag{4.12}
\]
Substitution in (4.5) yields
\[
 \mathcal E_2(QF^m)
  =(4m+2)!\,2^{-m}B_m
  =\frac{(4m+2)!\,m!}{(2m+1)!!}.
\]
This proves both identities in Theorem 1.1 for all \(m\) a second time.
The only general inputs are the monomial contraction rule, the formal beta
identity (4.1), and the elementary finite-difference identity (4.6).

## 5. Propagation to every balanced degree \(d\geq4\)

The quartic witness propagates without changing its angular cancellation.
The mechanism is multiplication by powers of the invariant \(R\).

> **Lemma 5.1 (radial shift).** If \(H\) is bihomogeneous of bidegree
> \((n,n)\), then for every \(a\geq0\),
> \[
>  \mathcal E_2(R^aH)
>  =\frac{(n+a+1)!}{(n+1)!}\mathcal E_2(H).
>  \tag{5.1}
> \]

There are again two proofs.  On circular Gaussians, \(R=S\), so radial
independence gives (5.1).  Directly and algebraically, the substitution in
(4.2) sends \(R\) to \(1\).  Formula (4.2), applied in bidegrees \(n\) and
\(n+a\), therefore gives
\[
 \mathcal E_2(R^aH)
 =(n+a+1)!\operatorname {CT}_u\int_0^1
 H\left(1,u,t,\frac{1-t}{u}\right)\,dt
 =\frac{(n+a+1)!}{(n+1)!}\mathcal E_2(H).
\]
Thus the propagation also has a proof over an arbitrary
characteristic-zero field with no Gaussian input.

For \(d\geq4\), put
\[
 \boxed{F_d=R^{d-4}F,\qquad Q=Z.}
 \tag{5.2}
\]
Then \(F_d\in V_d\).  Apply Lemma 5.1 first to \(H=F^m\), which has
bidegree \((4m,4m)\), with \(a=(d-4)m\).  Apply it again to \(H=QF^m\),
which has bidegree \((4m+1,4m+1)\).  Theorem 1.1 gives
\[
\begin{aligned}
 \mathcal E_2(F_d^m)&=0,\\
 \mathcal E_2(QF_d^m)
 &=\frac{(dm+2)!}{(4m+2)!}\mathcal E_2(QF^m)\\
 &=\boxed{\frac{(dm+2)!\,m!}{(2m+1)!!}}\ne0
 \qquad(m\geq1).
\end{aligned}
\tag{5.3}
\]

> **Theorem 5.2.** The moment--nullcone assertion
> \(\mathrm{MN}_d\) fails for every \(d\geq4\).

The fixed multiplier \(Q\) is independent of \(d\), and every mixed moment
in (5.3) is nonzero.  Since a one-sided nullcone point has eventual mixed
vanishing for every fixed multiplier, \(F_d\) is outside the nullcone.
Thus (5.3) supplies a moment-zero semistable point in every balanced degree
\(d\geq4\).

The radial order can be bounded independently of \(d\) by combining
ordinary powers with invariant multiplication.  More generally, let
\[
 d=4r+k,\qquad r\geq1,\quad k\geq0,
\]
and put
\[
 G_{r,k}=R^kF^r.
 \tag{5.4}
\]
The same calculation, now starting from \(F^{rm}\) and \(QF^{rm}\), gives
\[
\boxed{
 \mathcal E_2(G_{r,k}^m)=0,\qquad
 \mathcal E_2(QG_{r,k}^m)
 =\frac{((4r+k)m+2)!\,(rm)!}{(2rm+1)!!}\ne0.}
 \tag{5.5}
\]
Thus for \(r=\lfloor d/4\rfloor\) and \(k=d-4r\), every degree \(d\geq4\)
has a witness whose \(R\)-adic order is at most three.  Moreover \(F\) is
not divisible by \(R\): at
\[
 (\xi _1,\xi _2,z_1,z_2)=(1,1,1,-1)
\]
one has \(R=0\) and \(F=-2\).  Since \(R\) is irreducible, the exact
\(R\)-adic order of \(G_{r,k}\) is \(k\).  In particular:

> **Corollary 5.3.** Every degree divisible by four has an
> \(R\)-primitive counterexample \(F^{d/4}\), and every other degree
> \(d\geq4\) has one divisible by no more than \(R^3\).

The primitive family in degrees divisible by four need not consist of
proper powers.  Put
\[
 D=2WR^2-2T^2R-ZT^2
 \tag{5.6}
\]
and, for \(h\geq1\), define
\[
 \Phi_h=(R+Z)D
 \sum_{j=0}^{h-1}\binom{h-1}{j}
 R^{4(h-1-j)}T^{2j}(R+Z)^{2j}.
 \tag{5.7}
\]
Every summand has balanced degree \(4h\).  On the Hopf sphere this is the
one-profile polynomial from the
[Hopf-lift classification](HOPF_LIFT_CLASSIFICATION.md)
\[
 x^{-1}(1+x)\rho_h\!\left(t^2(1+x)^2\right),
 \qquad
 \rho_h(z)=(1-z)(1+z)^{h-1}.
 \tag{5.8}
\]
The endpoint-contact argument gives
\[
\begin{aligned}
 \mathcal E_2(\Phi_h^m)&=0,\\
 \mathcal E_2(Q\Phi_h^m)
 &=(4hm+2)!\,C_{h,m},\\
 C_{h,m}
 &=\int_0^1(1-v^2)^m(1+v^2)^{(h-1)m}\,dv>0.
\end{aligned}
\tag{5.9}
\]
The integral is rational by direct polynomial integration.  At the point
used above, \(R=0\), \(Z=-1\), and \(T=2\); only the \(j=h-1\) summand
survives and is nonzero.  Hence \(\Phi_h\) is \(R\)-primitive.  The factor
\(R+Z\) has multiplicity exactly one: neither \(D\) nor the sum in (5.7)
vanishes identically on \(R+Z=0\).  Thus \(\Phi_h\) is not a proper power.
Notice that \(\Phi_1=2F\).

> **Theorem 5.4.** Every balanced degree divisible by four admits an
> explicit \(R\)-primitive, non-proper-power failure of
> \(\mathrm{MN}_d\).

The congruence restriction is intrinsic to the entire one-profile class,
not just to (5.8).  If the profile polynomial has degree \(h\), its leading
term in
\[
 c\,x^{-r}(1+ax)^r\rho\!\left(t^2(1+ax)^2\right)
\]
is a nonzero scalar multiple of \(x^{2h}t^{2h}\), of total angular degree
\(4h\).  Endpoint contact makes the full Laurent expression polynomial and
introduces no term of larger degree.  Hence its minimal balanced
homogenization has degree \(4h\); placing it in any larger \(V_d\) adds
the radial factor \(R^{d-4h}\).  Therefore no one-profile Hopf lift can
produce an \(R\)-primitive witness in degrees \(1,2,\) or \(3\) modulo
four.  Reaching those classes requires several phase profiles, odd height
dependence, or a different cancellation mechanism.

There is a sharp obstruction internal to the invariant-multiplier
strategy.  The first fundamental theorem for one vector and one covector
gives
\[
 k[\xi _1,\xi _2,z_1,z_2]^{\mathrm{SL}_2}=k[R].
 \tag{5.10}
\]
Consequently every homogeneous invariant multiplier of balanced degree
\((k,k)\) is a scalar multiple of \(R^k\).  No alternative invariant
multiplier can remove the residual factors \(R,R^2,R^3\) in the three
nonzero congruence classes.  Equivalently, invariant-power propagation is
zero in the primitive quotient \(V_d/RV_{d-1}\) whenever it uses radial
padding.

This also records exactly what the two successful operations contribute.
Ordinary powers, the simplest plethystic operation here, produce the
\(R\)-primitive congruence class \(d\equiv0\pmod4\).  Tensoring with the
invariant line spanned by \(R^k\) fills the other degrees.  General Cartan
projections, plethysms other than ordinary powers, polarizations, block
embeddings, and differential intertwiners have no automatic scalar
radial-shift identity, so no moment-preservation claim is made for them.
The existence of \(R\)-primitive witnesses in degrees not divisible by
four remains open.  A Cartan or trace-free projection would have to solve
new moment identities rather than inherit them from the radial-shift
lemma.

## 6. Frontier consequences

The family begins in balanced bidegree \((4,4)\).  It is therefore
compatible with the complete positive theorem in bidegree \((2,2)\) and
does not decide the still-open full bidegree-\((3,3)\) stratum.  It
falsifies the proposed two-pair moment--nullcone equality in every degree
\(d\geq4\): the pure moments of \(F_d\) all vanish, while the fixed
multiplier \(Q\) proves that \(F_d\) cannot be one-sided.

Because the coefficient matrix has full rank, the witness does not
contradict the split-symbol GVC theorem.  That theorem concerns rank-one
forms \(A(\xi)P(z)\), whereas (1.4) is nonseparable.  The two-variable GVC
frontier and the ordinary-Laplacian polarization problem therefore remain
open in exactly their previously stated nonhomogeneous or degree-raising
forms.

## Reproduction

Run

```bash
python3 scripts/verify_two_pair_image_mathieu_counterexample.py
python3 scripts/audit_two_pair_image_mathieu_coefficient_extraction.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.SIC2C4FiniteSum
```

The dependency-free checker builds \(F\) from (1.2), verifies (1.3), checks
the sixteen-term expansion and coefficient-matrix determinant, and performs
exact sparse contraction through \(m=8\).  It also checks (5.3) for
\(4\leq d\leq10\) and \(1\leq m\leq4\), checks the bounded-radial-order
family (5.4)--(5.5) for \(4\leq d\leq15\) and \(1\leq m\leq3\), checks
the non-power profiles (5.7)--(5.9) for \(1\leq h\leq4\) and
\(1\leq m\leq3\), and replays through \(m=99\) the two finite sums
obtained by expanding (3.3).  The arguments in Sections 2--5, rather than
any finite cutoff, prove (1.5), (5.3), (5.5), and (5.9) for every allowed
parameter.

The second dependency-free audit is deliberately separate.  It checks the
formal coefficient-extraction substitution, the polynomial-division
residue identity (4.8), the finite-difference cancellations, and the recurrence
(4.11) using exact rational arithmetic.  Its finite cutoff is only a
regression; the displayed degree bounds, divisibility at \(X=-1/2\), and
termwise recurrence are the independently checkable all-order certificate.
The Lean module formalizes the general finite-difference cancellation, the
rank-one residue identity (4.8), the specialized normalized products
(4.7)--(4.9), and the recurrence and explicit double-factorial evaluation
of \(B_m\).  It also formalizes the coefficient functional and beta
recurrence, the linear assembly of (4.2), the literal Laurent witness (4.3),
its all-order binomial expansion, both constant-term extractions, and the
closed values in (4.4)--(4.5).  It now also defines the original
four-variable \(F,Q\) as `MvPolynomial` objects and identifies their
displayed substitution with the Laurent witness, including all pure and
mixed powers.  The only remaining Lean integration seam is a wrapper
identifying the generic balanced coefficient array in the contraction
theorem with `MvPolynomial.coeff`; the algebraic proof itself is formalized
on both sides of that representation boundary.

No term-count minimality or literature-priority claim is made here.

## Source for the Image framework

A. van den Essen, D. Wright, and W. Zhao,
[*On the Image Conjecture*](https://arxiv.org/abs/1008.3962),
J. Algebra 340 (2011), 211--224, supplies the contraction
image-kernel identity and the one-pair positive theorem used for the sharp
dimension conclusion.
