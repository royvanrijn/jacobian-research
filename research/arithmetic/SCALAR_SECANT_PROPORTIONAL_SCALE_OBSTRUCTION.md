# Uniform scalar-secant obstructions for paired Laurent pencils

For an odd Laurent pencil with even numerator degree, two natural
simplifications of the quadratic critical-algebra chart are impossible:
a constant secant and proportional reciprocal scale over R, or a
constant secant and constant reciprocal scale even over C. These exclude
both affine-scale endpoints in every degree 4h, without assuming real
critical points.

The conclusion concerns this sublocus of the
[paired critical trace chart](PAIRED_CRITICAL_TRACE_CHART.md). The general
nonreciprocal locus and a nonconstant secant element remain outside the
obstructions. For a scalar secant and affine scale rho=r_0+cM, the case
r_0*c!=0 remains outside these uniform arguments. In degree twelve,
the separate [rational affine-scale obstruction](DEGREE_TWELVE_SCALAR_AFFINE_OBSTRUCTION.md)
now excludes that case over Q as well. Higher degrees and irrational
real affine parameters are not covered by that rational result.

## Statement

Let h>=1, L=4h-1, and let N in R[u] be monic of degree 2h.
Assume R_N=2uN'-N has 2h distinct nonzero roots. At these roots assume
M=N(u)^2/u takes h distinct nonzero values, each exactly twice.
Thus f(z)=N(z^2)/z has the good paired-critical passport in degree 4h.

In the critical algebra A=R[u]/(R_N/L), write nu=N(u) and M=nu^2/u.
The trace-chart theorem supplies the paired-value algebra B=R[M] and
unique units a,rho in B satisfying nu=a(u-rho).

**Theorem.** It is impossible that a belongs to R* and rho=cM for
some c in R*. Equivalently, there are no real constants a,c != 0 with

\[
                    acN(u)^2+uN(u)-au^2=0
                    \quad\bmod R_N.
\tag{1}
\]

The same exclusion holds over every subfield of R, including Q.
There is no assumption on the locations of the critical points.

## The two constant slopes

Assume such a model exists and set Z=nu/u. All u,nu,M are units.
Dividing nu=a(u-cM) by u and using M/u=Z^2 gives

\[
                         acZ^2+Z-a=0.                       \tag{2}
\]

Let z_1,z_2 be its complex roots. Every equal-M pair has one of each:
if Z had the same value on the two members, M=Z^2u would make their
distinct u coordinates equal. Thus the roots z_1,z_2 are distinct and
nonzero, and each occurs at exactly h critical roots.

Put t=z_1/z_2. If alpha is a critical root with slope z_1 and beta
its equal-M partner, then

\[
 M=z_1^2\alpha=z_2^2\beta,\qquad
 \beta=t^2\alpha,\qquad
 N(\beta)=tN(\alpha).                                     \tag{3}
\]

In particular t != 0 and t^2 != 1, since the paired roots are distinct.

Let S_1 be the monic polynomial on the h slope-z_1 critical roots,
and let S_2 be its partner polynomial. They have nonzero constants and

\[
 S_2(u)=t^{2h}S_1(u/t^2),\qquad
                     R_N=L S_1S_2.                         \tag{4}
\]

At every root alpha of S_1, the critical equation gives
N'(alpha)=z_1/2 and N'(t^2alpha)=z_2/2. Therefore

\[
                         H(u)=N(t^2u)-tN(u)
\]

has both value and derivative zero at alpha. Its h distinct double
zeros show that S_1^2 divides H.

The coefficient of u^(2h) in H is t^(4h)-t. It cannot vanish:
otherwise H would have degree below 2h, hence be zero, whereas
H(0)=(1-t)N(0) != 0. This also handles every possible leading-degree
drop, including nontrivial roots of t^(4h-1)=1. Consequently

\[
                 H(u)=(t^{4h}-t)S_1(u)^2.                   \tag{5}
\]

Evaluating R_N at zero and using (4) gives

\[
                  N(0)=-L t^{2h}S_1(0)^2.
\]

Compare constants in (5) and cancel the nonzero S_1(0)^2. Since
t != 0,1, the result is

\[
 \boxed{\sum_{j=0}^{4h-2}t^j=(4h-1)t^{2h-1}.}               \tag{6}
\]

## Why real scalar parameters cannot satisfy the identity

The two roots of the real quadratic (2) are either both real or
complex conjugates. Hence t is real, or |t|=1.

If t>0, the arithmetic-geometric mean inequality applied to the
4h-1 positive terms 1,t,...,t^(4h-2) makes the left side of (6)
at least (4h-1)t^(2h-1). Equality requires all terms equal,
hence t=1, already excluded.

If t<0, the left side is

\[
                       \frac{t^{4h-1}-1}{t-1}>0,
\]

whereas the right side of (6) is negative. This is impossible.

Finally, if |t|=1 and t!=1, divide (6) by t^(2h-1). It becomes

\[
                 1+2\sum_{j=1}^{2h-1}\operatorname{Re}(t^j)
                                      =4h-1.
\]

Every real part is at most one and the first is strictly below one.
The left side is strictly smaller than the right side. This last
contradiction proves the theorem.

The argument uses real scalar a,c; it does not assert the same
obstruction for arbitrary complex scalar parameters.

## A polynomial trace witness

The ratio also has the real scalar trace

\[
                    x=t+t^{-1}=-2-\frac1{a^2c}.             \tag{7}
\]

Define T_0=2, T_1=x, T_(j+1)=xT_j-T_(j-1), and

\[
 B_h(x)=\frac{1+\sum_{j=1}^{2h-1}T_j(x)-(4h-1)}{x-2}.
\tag{8}
\]

This is a monic integer polynomial of degree 2h-2. Equation (6)
would require B_h(x)=0. In fact B_h(x)>0 for every real x:
for x>2 use the strict arithmetic-geometric mean argument; for
x<-2 use the opposite signs for negative t; for -2<=x<2 use the
strict unit-circle bound. At x=2,

\[
                         B_h(2)=\sum_{j=1}^{2h-1}j^2>0,
\]

as follows by differentiating T_j at 2.

For degree twelve the witness is

\[
                         B_3(x)=x^4+3x^3+2x^2+x+5.
\]

This witness makes the obstruction independent of any numerical
critical-root search or rational-denominator assumption.

## Proposition 2. Both scalar elements are impossible over C

Keep the good critical-pair hypotheses, now allowing N to have complex
coefficients, and put m=2h. It is impossible that both a and rho are
nonzero base scalars. In particular this excludes c=0 in the affine
scale rho=r_0+cM, including nonreal critical points and coefficients.

**Proof.** Write S=R_N/(2m-1). The monic polynomial
N-a(u-rho) vanishes at all m distinct roots of S, so equality of
degree and leading coefficient gives

\[
                 S=N-au+a\rho,\qquad R_N=(2m-1)S.
\tag{9}
\]

Euler coefficient comparison in (9) forces

\[
\begin{split}
 N(u)&=u^m+\frac{(2m-1)a}{2m-2}u
                      -\frac{(2m-1)a\rho}{2m},\\
 S(u)&=u^m+\frac{a}{2m-2}u+\frac{a\rho}{2m}.
\end{split}                                                   \tag{10}
\]

The canonical trace-chart norm relation is u*tau(u)=rho^2.
Consequently the roots of S are stable under u->rho^2/u,
and their product, obtained by multiplying the h paired products,
is rho^m. Since m is even,

\[
 S(0)=\rho^m,\qquad
 S(u)=\frac{u^m}{\rho^m}S(\rho^2/u),\qquad
                         a=2m\rho^{m-1}.                    \tag{11}
\]

If m>=4, the actual coefficient of u^(m-1) in (10) is zero.
The reciprocal expression in (11) instead gives this coefficient as

\[
 \frac{\rho^2}{\rho^m}\frac{a}{2m-2}
                            =\frac{m\rho}{m-1}\ne0.
\]

This is a contradiction. If m=2, (11) gives rho=a/4 and

\[
                         S(u)=(u+\rho)^2,
\]

contrary to the simple critical-root hypothesis. This proves the
proposition over C, and hence over every characteristic-zero
subfield. No positivity or real-structure argument is used. The
original theorem and this proposition cover r_0=0 and c=0,
respectively, for a scalar-secant affine scale. They do not exclude
the case r_0*c!=0.

## Evidence and limits

The [checker](verify_scalar_secant_proportional_scale.py) verifies
the constant-slope equation, the double-zero jets, the constant-term
identity, and the polynomial trace witnesses in bounded degrees.
It uses exact Sturm counts to confirm their absence of real roots
in those controls. It also verifies the Euler coefficient identities
and reciprocal coefficient contradiction for Proposition 2.
The uniform statements rest on the written proofs, not on the bounded
checks or a Groebner experiment.

The [receipt](scalar_secant_proportional_scale_verification.json)
and [source ledger](scalar_secant_proportional_scale_source_ledger.json)
record exact source bytes and internal checking. No general rational
point theorem, complete degree-twelve obstruction, novelty claim,
formal verification, independent external replay or external human
review is asserted.
